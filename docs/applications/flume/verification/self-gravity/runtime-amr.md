# Self-gravity under runtime AMR

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-11 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5c | CPU, FNL | yes | — |

## Case and scope

Under runtime AMR the grid, and with it the solver, change during the run. Every regrid must:
- set up the right solver: the FFT on a uniform grid, the multigrid on an AMR one;
- solve $\Phi$ from the regridded state;
- rebuild the flux register the seam gravity work rides on.

**What changes and what does not.** The regrid itself conserves $\rho$, $\rho\mathbf u$ and $E$ (conservative
prolongation and restriction). $W = \tfrac12\sum\rho\Phi V$, though, is the energy of the grid's own discrete
Laplacian, so a regrid that changes the grid changes $W$, and $E + W$ jumps by a truncation-sized amount. Between
regrids the scheme stays exact ([SV-10](./amr-jeans)). The leg separates the two.

What the case catches:
- a regrid that keeps the old solver or the old potential: the step changes of $E + W$ between regrids leave the
  time-integration level;
- a jump of $E + W$ off a regrid step (a leak that is not the change of the discrete $W$);
- a momentum or mass leak through the regrid.

Isolated boundaries under runtime AMR, refused here until [#87](https://github.com/szaghi/adam/issues/87), are
verified with the polytrope ([SV-14](./polytrope)), with the refusal of a regrid that refines a boundary block.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, $\gamma = 5/3$, periodic self-gravity | [models](../../models#self-gravity) |
| Regrid | every 5 steps (`[amr] frequency = 5`), one level above the base, a Löhner estimator on the density (mode 4): refine above 0.3, derefine below 0.1 | `regrid_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Solver after a regrid | `set_potential` sets the solver up on the new grid (FFT if the leaves are at one level, multigrid otherwise) and solves $\Phi$ from the regridded state; the FNL backend copies it to the device | `set_potential`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `refresh_potential`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Seam gravity work | in the per-stage reflux | `add_seam_gravity_work`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space, time | WENO5 flux splitting (`weno`), SSP-RK3, CFL 0.4 | [numerics](../../numerics#scheme-space-weno-flux-splitting) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Case | `jeansc2d`: the collapse of [SV-4](./jeans), $G = 9\pi/2$, amplitude $10^{-2}$, the box $[1/8, 9/8]^2$ ([`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), `wave_mode = jeans`) |
| Grid | base $32^2$ (4×4 blocks of $8^2$, two uniform levels), at most one level more by the regrid; 4 ghosts |
| Boundaries | every active face `periodic` |
| Reflux | per stage (`reflux_cadence = stage`) |
| Ranks | 2 |
| End | $t = 0.4$ (the figure's run: a checkpoint every 10 steps) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/make_jeans.py` | `--regrid N` writes the Löhner marker and `[amr] frequency` |
| `src/tests/flume/verification/selfgravity/jeans.sh` | leg `sv11` (`FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/selfgravity/jeans_oracle.py` | `regrid`: the grid must change; every jump of $E + W$ (a step change above 100 times the median) on a regrid step; the other step changes, the momentum and the mass bounded |
| `src/tests/flume/verification/make_doc_runs.sh`, group `amr` | `sv11-regrid`, a checkpoint every 10 steps |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_amr_regrid` | draws the figure (`--only selfgravity-amr-regrid`) |

```bash
cd src/tests/flume/verification/selfgravity && ./jeans.sh --legs sv11                       # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./jeans.sh --legs sv11                   # FNL
```

## Results

![SV-11: the change of E + W through the regrid, the blocks, and the density on the final composite grid](/flume/selfgravity-amr-regrid.png)

*Left: $E + W$ drifts at the time integration's level until the regrid at $t = 0.35$ refines the ridge (16 → 52
blocks), then jumps once and stays flat. Right: the final density, the refined blocks along the ridge.*

| Asserted | Bound | CPU | FNL (WSL GPU) |
|---|---|---|---|
| the grid changes (blocks of the first and the last checkpoint) | must differ | 16 → 52 | 16 → 52 |
| every jump of $E + W$ on a regrid step | yes | one jump, at step 60: $-2.3\cdot10^{-4}$ | the same, $-2.3\cdot10^{-4}$ |
| the largest other step change of $E + W$, relative | $10^{-6}$ | $7.6\cdot10^{-8}$ | $7.6\cdot10^{-8}$ |
| max $\lvert$momentum$\rvert$; mass drift | $10^{-14}$; $10^{-14}$ | $2.9\cdot10^{-18}$; 0 | $4.7\cdot10^{-18}$; $1.1\cdot10^{-16}$ |

Related cases: the same collapse on a static AMR grid is [SV-10](./amr-jeans), on a uniform grid [SV-4](./jeans);
runtime regridding without gravity is [runtime regrids](../runtime-amr/runtime-regrids).

## Findings

- **The regrid jump is the change of the discrete $W$.** At the regrid $E$ is unchanged and $W$ moves by
  $1.4\cdot10^{-3}$ of itself, the size of the $O(h^2)$ truncation at $32^2$. It does not converge in a resolution
  ladder, because the refinement happens at different stages of the collapse (at base $64^2$ the ridge crosses the
  threshold later, when it is sharper). The jumps are reported, not bounded.
- **Correcting $E$ at each regrid** would make the regrid conserve $E + W$; it is not done.
- **The solver follows the grid.** With refine and derefine thresholds of 0.1 and 0.03 the same collapse refines
  everywhere by $t = 0.27$, so the grid becomes uniform again and the FFT solver takes over. That run passes the same
  oracle, with jumps at steps 40 and 50.
- No regression case freezes self-gravity under runtime AMR; `selfgravity-amr-stage` covers the per-stage reflux on a
  static grid.
