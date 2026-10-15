# Positivity limiter

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| PV-1 to PV-5 | Ideal MHD | M3 ([#47](https://github.com/szaghi/adam/issues/47), M3-P5); PV-5 [#50](https://github.com/szaghi/adam/issues/50) | CPU, FNL (PV-1 as a sphere: FNL only) | yes | `blast-limiter`, `blast3d-limiter`, `blast-amr-limiter` |

## Case and scope

The limiter (`positivity_limiter = cell`) blends every face flux with the first-order Lax–Friedrichs backbone so that
each stage keeps the density and the pressure positive. It must make the hard problems run without floors and leave the
smooth ones untouched. All legs run with the floors disabled, on the splitting scheme and on `weno-riemann` (HLLD for
MHD, HLLC for Euler, characteristic interpolation).

A leg passes when the run reaches its final time and no stage meets an inadmissible backbone (the log line
`positivity limiter: N faces limited, M inadmissible backbones` has $M = 0$ at every stage), and no cell has a
non-finite high-order flux. The Euler legs, which have no floors to stop a run, also assert a positive density and
pressure in the last checkpoint.

What the case catches:
- **loss of positivity** on strongly magnetised blasts (PV-1, PV-2) and near vacuum (PV-3);
- **a limiter that acts on smooth flow** (PV-4): with the limiter on, smooth runs must be bitwise equal to the runs
  without it;
- **positivity lost at 2:1 seams** (PV-5): the end-of-step reflux and the seam ghost fill must not undo the limiter;
  the mass must stay constant and the coarse seam flux must equal the mean of the fine ones at every stage.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Limiter | cell-based parametrised flux limiter: each face flux blended with the first-order backbone; the relative floor keeps $\rho$ and the internal energy $\ge \kappa$ times their backbone values, $\kappa = 0.1$ | `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_positivity_factors`, `compute_backbone_fluxes`, `blend_positivity_fluxes`, [`adam_flume_cpu_positivity_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_positivity_kernels_agnostic.INC); `POSITIVITY_LIMITER_KAPPA`, [`adam_flume_parameters.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_parameters.F90); [numerics](../../numerics#positivity-limiter) |
| MHD backbone | Lax–Friedrichs flux with the Wu speed (no-cleaning and EGLM variants) | `mhd_backbone_flux`, `mhd_eglm_backbone_flux`, [`adam_flume_mhd_riemann_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_riemann_library.F90) |
| EGLM sources under the limiter | sources taken with the limited arithmetic; unlimited where $\Lambda = 1$ | `add_eglm_sources_limited`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90) |
| Seams (PV-5) | per-stage seam flux synchronisation; ghost positivity blend of the tricubic seam ghosts | `seam_sync_fluxes`, `seam_sync_factors`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `blend_inadmissible_ghosts`, [`adam_flume_cpu_ghost_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_ghost_kernels_agnostic.INC) |
| Equations | ideal MHD with EGLM (PV-1, PV-2, PV-5), $\gamma = 1.4$; Euler (PV-3), $\gamma = 5/3$ (LeBlanc) or 1.4 | [models](../../models#divergence-control-eglm) |
| Space | WENO5 flux splitting; `weno-riemann` HLLD (MHD) or HLLC (Euler), characteristic, `6th`, `weno` sensor | `compute_residuals_weno`, `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | MHD: SSP Runge–Kutta 5 stages, CFL 0.4 (the Orszag–Tang setup); Euler: SSP-33, CFL 0.5 (the Sod setup) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Floors | disabled | [models](../../models#positivity-floors) |
| FNL backend | the same limiter on the device | `limit_positivity_dev`, `seam_sync_fluxes_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Leg | Problem | Grid | Final time |
|---|---|---|---|
| PV-1 | Balsara–Spicer blast: $\rho = 1$, $p = 0.1$ outside, $p = 1000$ for $r < 0.1$ around $(0.5, 0.5)$ (cell-centre sampling), $\lvert\mathbf{B}\rvert = 100/\sqrt{4\pi}$ at 45°, ambient $\beta = 2.5\cdot10^{-4}$; EGLM, $c_h = 60$ | $[0, 1]^2$ periodic, the Orszag–Tang grid, $128^2$ | 0.01 |
| PV-1, sphere | the same blast as a sphere | $[0, 1]^3$ periodic, octree $4 \times 4 \times 4$ blocks, $64^3$ | 0.01 |
| PV-2 | Wu & Shu (2018, Example 4.4): $p = 10^4$ in the disc, $\mathbf{B} = (1000/\sqrt{4\pi}, 0, 0)$, $\beta = 2.51\cdot10^{-6}$; EGLM, $c_h = 400$ | $128^2$ | 0.001 |
| PV-3, LeBlanc | $\gamma = 5/3$, $[0, 9]$, $(\rho, u, e) = (1, 0, 0.1)$ for $x < 3$, $(10^{-3}, 0, 10^{-7})$ beyond | 400 cells (4 blocks of 100), the Sod grid, extrapolation | 6 |
| PV-3, double rarefaction | $\gamma = 1.4$, $[-1, 1]$, $(\rho, u, p) = (7, \mp1, 0.2)$ | 400 cells ($h = 0.005$) | 0.6 |
| PV-3, planar Sedov | $\gamma = 1.4$, $[0, 4]$ shifted by $h/2$ (a cell centred at $x = 2$), $(\rho, u, p) = (1, 0, 4\cdot10^{-13})$, $p = 2.56\cdot10^8$ in the centre cell | 800 cells ($h = 0.005$) | 0.001 |
| PV-4 | the isentropic vortex ([V2](../euler/isentropic-vortex), Euler, $64^2$, $128^2$) and the 1-D fast wave (EGLM, 16, 32 cells, amplitude $10^{-7}$) | their verification grids | theirs |
| PV-5 | PV-1 at $64^2$ base with $[0.25, 0.75]^2$ refined one level (octree): the outer shock crosses the 2:1 seams | `amr_box.py` | 0.01 |

All legs: 2 ranks. Initial conditions: the blasts and the near-vacuum problems use
[`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions) regions (one
strip region per cell row for the disc). PV-3 setups as Zhang & Shu (2010) and Guo et al. (arXiv:1402.5618, Examples
8–10).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/positivity-probe/make_blast.py` | writes the blasts (`--cells`, `--eglm`, `--3d`, `--b-axis x`, `--p-in`, `--b0`, `--time-max`, `--glm-ch`, `--refine-box`, `--ratio`) |
| `src/tests/flume/verification/mhd/positivity/make_vacuum.py` | writes the PV-3 inputs from `sod/sod-x.ini` (`--problem leblanc\|double-rarefaction\|sedov`) |
| `src/tests/flume/verification/mhd/positivity/check.sh` | runs the legs (`--legs pv1,pv1-3d,pv2,pv3,pv4,pv5`; default `pv1,pv2,pv3,pv4`) and the verdicts |
| `src/tests/flume/verification/make_doc_runs.sh` groups `blast`, `vacuum` | re-run PV-1, PV-2 (4 runs) and PV-3 (6 runs) with checkpoints kept |
| `src/tests/flume/verification/make_doc_figures.py` `fig_blast`, `fig_vacuum` | draw `blast.png`, `blast-wushu.png` and `near-vacuum.png` (`--only blast near-vacuum`) |

```bash
cd src/tests/flume/verification/mhd/positivity
./check.sh                       # PV-1, PV-2, PV-3, PV-4 (about 70 min on the CPU, PV-3 takes 55)
./check.sh --legs pv1-3d         # the blast as a sphere on 64^3 (about 1 h per scheme on the CPU)
./check.sh --legs pv5            # the blast across 2:1 AMR seams (about 14 min per scheme on the CPU, 19 on FNL)
```

## Results

![Balsara-Spicer blast with the limiter](/flume/blast.png)

*Balsara–Spicer blast, $\beta = 2.5\cdot10^{-4}$, $128^2$, $t = 0.01$, EGLM, limiter on: density and magnetic pressure,
splitting scheme (top) and `weno-riemann` HLLD (bottom).*

![Wu-Shu blast with the limiter](/flume/blast-wushu.png)

*Wu–Shu blast, $\beta = 2.51\cdot10^{-6}$, $128^2$, $t = 0.001$, the same layout.*

![Euler near vacuum with the limiter](/flume/near-vacuum.png)

*LeBlanc, the double rarefaction and the planar Sedov blast with the limiter, splitting scheme and `weno-riemann`
HLLC, with the exact solutions of LeBlanc and of the double rarefaction.*

| Problem | Without the limiter | Splitting + limiter | Riemann + limiter |
|---|---|---|---|
| PV-1, Balsara–Spicer blast, $\beta = 2.5\cdot10^{-4}$, $128^2$, EGLM | splitting fails at step 3; HLLD runs | passes, up to 1540 faces limited per stage | passes, up to 280 |
| PV-1 as a sphere, $64^3$ (FNL) | not run | passes, up to 51566 | passes, up to 2976 |
| PV-2, Wu–Shu blast, $\beta = 2.51\cdot10^{-6}$, $128^2$, EGLM | splitting fails at step 1, HLLD at step 41 | passes, up to 10333 | passes, up to 352 |
| PV-3, LeBlanc shock tube, 400 cells | both run | passes, never limited | passes, up to 512 |
| PV-3, double rarefaction, 400 cells | splitting runs; HLLC fails at step 79 | passes, never limited | passes, up to 1024 |
| PV-3, planar Sedov blast, 800 cells | splitting runs; HLLC fails at step 69 | passes, never limited | passes, up to 768 |
| PV-4, isentropic vortex ($64^2$, $128^2$) and fast wave (16, 32 cells, EGLM) | — | bitwise equal to the run without the limiter | — |
| PV-5, PV-1 at $64^2$ with $[0.25, 0.75]^2$ refined 2:1 (issue #50) | not run (uniform: splitting fails at step 3) | passes; mass to round-off, seam flux mismatch $3.6\cdot10^{-12}$ | passes; the same |

The face counts include the copies along the null directions (256 per 1-D face). They are those of the CPU; the device
ones differ in the splitting blasts (1468 and 10955), whose limited runs are not bitwise reproducible across backends,
and agree elsewhere. The 3-D blast was run with the relative floor on the device only.

| Assertion | Bound | CPU | FNL |
|---|---|---|---|
| final time reached, inadmissible backbones, non-finite fluxes (every leg) | $t_\text{end}$, 0, 0 | pass | pass |
| PV-3: final $\rho$, $p$ positive | $> 0$ | $\rho \ge 9.9\cdot10^{-4}$, $p \ge 4\cdot10^{-13}$ | pass |
| PV-4: interior of the last checkpoint, limiter on against off | bitwise | bitwise | bitwise |
| PV-5: relative mass drift | $10^{-12}$ | round-off | round-off |
| PV-5: forest $\max\lvert F_\text{coarse} - F_\text{fine}\rvert$ at the end-of-step reflux | $10^{-9}$ | $3.6\cdot10^{-12}$ | the same |

Related cases: [EGLM cleaning](./eglm) shows the blast failing without the limiter; [quadtree AMR](./quadtree-amr) runs
the limited blast on a quadtree; [the MHD blast symmetry](../runtime-amr/mhd-blast-symmetry) studies the blast under
runtime AMR; [Sod](../euler/sod) and [Lax](../euler/lax) are the Euler Riemann problems the HLLC path is checked on.

## Findings

- Both schemes reach the final time, but not with the same quality: the density of the splitting scheme carries
  grid-scale noise around the blast, mild at $\beta = 2.5\cdot10^{-4}$ and plain at $\beta = 2.51\cdot10^{-6}$, where
  HLLD stays clean. The limiter guarantees positivity, not accuracy: on strongly magnetised blasts `weno-riemann` with
  HLLD is the scheme to use.
- The Euler splitting scheme needs no limiting on any of the three near-vacuum problems at these resolutions: the
  limiter matters for `weno-riemann`.
- On LeBlanc both schemes put the shock ahead of its exact position $x = 7.96$: by 0.33 at 400 cells and by 0.19 at
  800 (measured on the splitting scheme). The error converges slowly under refinement, as reported for this problem in
  the literature, and the limiter does not cause it (the splitting scheme never limits here). On the double rarefaction
  the velocity is ill-defined where the density nearly vanishes, hence the kink at $x = 0$. On Sedov HLLC reaches a lower
  centre density than the splitting scheme ($4\cdot10^{-3}$ against $7\cdot10^{-3}$), and the ambient pressure stays at
  its initial $4\cdot10^{-13}$.
- Sedov on HLLC set the limiter's relative floor. With the absolute floor alone the run stopped at step 10, at every
  CFL from 0.2 to 0.5: one stage drained the centre cell to $\rho = 10^{-13}$ with its energy kept (the sound speed grew
  by $10^6$), and the next stage, whose time step is that of the step's first state, had no admissible backbone. With
  the relative floor ($\kappa = 0.1$) the limiter acts in 15 stages and the run completes; the smooth cases of PV-4
  stay bitwise.
- The splitting scheme takes 614 steps on PV-2 against 284 for HLLD (779 with the absolute floor alone): its limited
  states keep low-density cells with large Alfvén speeds, which cost time steps but not admissibility.
- PV-4 is bitwise because a face whose two cells keep $\Lambda = 1$ is not touched, and the EGLM sources with
  $\Lambda = 1$ take the unlimited arithmetic.
- Before issue [#50](https://github.com/szaghi/adam/issues/50) neither scheme reached the end of PV-5. The splitting
  scheme stopped at step 184: the end-of-step reflux drove coarse cells beside the seam to $\rho e < 0$ (it replaced the
  limited coarse flux with the restricted fine one, a correction no cell factor bounds). HLLD stopped at step 261, after
  the tricubic seam ghost fill had produced inadmissible fine ghosts. With the ghost positivity blend and the per-stage
  seam flux synchronisation both schemes reach $t = 0.01$ on the CPU and on FNL; the coarse seam flux is the mean of the
  fine ones at every stage, so the end-of-step reflux corrects round-off only.
