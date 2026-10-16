# MHD blast symmetry on the flux-splitting path

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| #77 study, bounded by av3c | Runtime AMR (M5) | M5 ([#74](https://github.com/szaghi/adam/issues/74)), [#77](https://github.com/szaghi/adam/issues/77) | CPU, FNL | yes | `blast-limiter` |

## Case and scope

The Balsara–Spicer blast of `blast-limiter` and [av3b](./accuracy) starts exactly point-symmetric. Ideal MHD maps it
onto itself under the rotation by π about the centre composed with $\mathbf B \to -\mathbf B$: ρ, ρE, B even, ρu, ψ
odd. On a uniform grid the run broke the symmetry at O(0.1) and the two backends disagreed at O(0.1), with no AMR
involved. #77 asks where it comes from, and either fixes it or explains and bounds it.

What the study would catch: an asymmetric defect of the scheme (an operator that is not invariant under the
rotation), a grid-scale instability, or a backend divergence hidden behind a symmetric initial condition.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 1.4$, EGLM ($c_h = 60$, $\alpha = 0.18$), floors disabled | [models](../../models#ideal-mhd), [EGLM](../../models#divergence-control-eglm) |
| Space, `weno` | WENO5 flux splitting, characteristic reconstruction; one Lax–Friedrichs speed per MHD wave | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_eglm_face_split_fluxes`, `mhd_eglm_eigenvectors`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Space, `weno-riemann` | HLLD with and without the 6th-order correction, LLF, HLLD with primitive interpolation | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_mhd_eglm_hlld_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_hlld_kernels.F90); [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Positivity limiter | `positivity_limiter = cell` | `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Time | `runge-kutta-ssp-54`, CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| FNL backend | the same scheme on the device | `compute_residuals_weno_dev`, `limit_positivity_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null, `periodic` on the $x$ and $y$ faces |
| Grid | the uniform blast of av3b: quadtree, $N = 128$ ($8 \times 8$ blocks of $16^2$ cells); $N = 64$ and $256$ in the resolution probe |
| Ranks | 2 |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions) regions: $\rho = 1$; $p = 1000$ in the disc $r < 0.1$ around $(0.5, 0.5)$ (sampled at the cell centres of $N = 64$), $p = 0.1$ outside; $\mathbf B = 100/\sqrt{4\pi}$ at 45°, $\beta = 2.5\cdot10^{-4}$ outside the hot disc |
| Final time | $t = 0.01$ (428 steps at $N = 128$); a checkpoint every 10 or 20 steps |
| Asymmetry of a field | $\max\lvert f(\mathbf x) - s f(\mathbf 1 - \mathbf x)\rvert / \max\lvert f\rvert$, $s = \pm1$ the parity of the field |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/regression/blast-limiter/input.ini` | the blast (written by `verification/mhd/positivity-probe/make_blast.py --cells 64 --eglm`); its golden is per backend |
| `src/tests/flume/verification/mhd/positivity-probe/make_blast.py` | the blast generator (`--b0` for the β probe) |
| `src/tests/flume/verification/regrid/make_regrid.py` (`blast-av`) | the av3b and av3c inputs |
| `src/tests/flume/verification/regrid/accuracy.sh --leg av3c` | the symmetry bound on the Riemann-flux path |
| `src/tests/flume/verification/regrid/av_oracle.py` (`blast`) | the point-symmetry measure, bound $10^{-10}$ |
| `src/tests/flume/verification/make_doc_runs.sh <runs> blast77` (FNL), `make_doc_figures.py --only mhd-blast-symmetry` | re-runs the uniform blast on the splitting scheme and on HLLD with a checkpoint every 20 steps |

```bash
cd src/tests/flume/verification/regrid && ./accuracy.sh --leg av3c          # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./accuracy.sh --leg av3c   # FNL
```

The probes of the table below are investigation runs, not legs of a script.

## Results

![The MHD blast point-symmetry defect](/flume/mhd-blast-symmetry.png)

*The uniform $N = 128$ blast on FNL (WSL GPU, 2 ranks), a checkpoint every 20 steps. Left: the splitting scheme's density at $t = 0.01$. Centre: its point-symmetry defect $|\rho(\mathbf x) - \rho(\mathbf 1 - \mathbf x)| / \max\rho$, spread over the whole blast, not along the grid axes. Right: the largest defect of $\rho$ and $B_y$ against the step: on the splitting scheme it leaves round-off by step 20, grows exponentially and saturates at 0.18 (this measure) by step 260; on `weno-riemann` HLLD it stays near $10^{-12}$, under the av3c bound.*

| Probe | Result |
|---|---|
| HEAD, after the #79 fix (CPU) | ρ 0.16 pointwise; in L1 only 2.3·10⁻²: the pointwise figure is a shock shifted by a few cells. CPU against FNL 0.16 |
| first steps | 1.2·10⁻¹⁰ after one step, at the front cells in the low-β gas; exponential growth from step 60 (e-folding in 10–15 steps), saturation at O(0.1) by step 220 |
| N = 64 / 128 / 256 | the asymmetry reaches 10⁻² at t = 0.0094 / 0.0047 / 0.0024: the growth rate scales as 1/Δx, a grid-scale instability, not a physical one (whose rate would converge) |
| CFL 0.4 → 0.2 | 10⁻² at t = 0.0047 / 0.0041: the spatial operator, not the time integrator |
| β × 10 / × 100 (B / √10, B / 10) | 6.8·10⁻⁸ at t = 0.01, still growing / 2·10⁻¹⁰, flat: a low-β effect |
| conservative reconstruction | the same growth, about 100 steps later |
| limiter off (floors 0.01), EGLM or GLM | NaN at step 8: the limiter keeps the run alive |
| `weno-riemann`: HLLD with and without the 6th-order correction, LLF, HLLD with primitive interpolation | symmetric to 10⁻¹² – 10⁻¹⁰ through t = 0.01 |
| positivity limiter activity | flux splitting: faces limited at every step (428 of 428; β × 100: 46 of 384); `weno-riemann`: only in the first 17–27 steps. The backbone flux the limiter falls back to is the same on both paths |

The bound that closes the study (av3c, `weno-riemann` with `riemann_solver = hlld`):

| Check | Bound | CPU | FNL |
|---|---|---|---|
| point symmetry, uniform $N = 128$ | $10^{-10}$ | $6.1\cdot10^{-12}$ | $2.6\cdot10^{-12}$ |
| point symmetry, tracked (3 regrids) | $10^{-10}$ | $3.3\cdot10^{-12}$ | $2.6\cdot10^{-12}$ |
| shock radius on the four half-axes | — | the same on both runs | — |

Related cases: the regridded blast is av3b and av3c of [AV](./accuracy); the limiter itself is the
[positivity limiter](../mhd/positivity-limiter) case; the Riemann-flux path is [MHD on `weno-riemann`](../mhd/weno-riemann).

## Findings

- **Diagnosis.** At β ≈ 2.5·10⁻⁴ the `weno` flux splitting keeps producing nearly inadmissible states at the blast
  front, the positivity limiter clips them at every step, and the round-off asymmetry of the arithmetic grows at the
  grid scale from 10⁻¹⁰ to O(0.1). The Riemann-flux path does not produce those states past the first steps and stays
  symmetric. The limiter's part rests on correlation and elimination (it acts at every step exactly where the scheme
  is unstable, and the scheme cannot run without it), not on a direct measurement of its amplification.
- **A fix tried and not kept.** One Lax–Friedrichs speed for the seven MHD waves (the stencil's fast bound) instead of
  one per wave, the cleaning pair unchanged, delays the growth (N = 128: 5.8·10⁻⁵ at t = 0.01; N = 256 saturates at
  4·10⁻³ instead of 0.17) but does not stop it, keeps the limiter acting at every step, and raises the MV-5 errors 2
  to 4 times on smooth waves (orders kept, every fixed bound failed).
- **Resolution (accept and bound).** On the flux-splitting path at β ≲ 10⁻³ the point symmetry, and the agreement of
  the two backends, hold only to the grid-scale instability above: the shock positions agree within a few cells, the
  pointwise fields do not. That path keeps its golden (`blast-limiter`: each backend is compared with its own golden,
  which a deterministic run reproduces; symmetry is not asserted). For low-β problems where symmetry or backend
  agreement matter, use `scheme_space = weno-riemann` with `riemann_solver = hlld`; av3c bounds that path.
