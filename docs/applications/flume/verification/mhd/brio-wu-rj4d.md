# Brio–Wu and Ryu–Jones 4d

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| — (`mhd/riemann`, M2-P3 exit) | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes | `brio-wu` |

## Case and scope

The Brio–Wu shock tube ($\gamma = 2$, $B_x = 0.75$, $B_y = \pm 1$) develops a compound wave; Ryu–Jones 4d a
switch-on/switch-off structure. Both stress the MHD eigensystem near its degeneracies, where a scheme that is not clean
needs the positivity floors to survive. The leg also exercises the floors themselves on a state no clean solver
produces: a negative pressure from the first stage.

What the case catches:
- **an unclean scheme:** with the floors armed at $10^{-10}$ the runs must reach the final time with zero floored cells
  in every stage, and a positive, finite final density and pressure;
- **floors that do not engage:** with `p_floor = 1e-3` the negative-pressure state must be floored and the run must
  complete;
- **floors that hide a failure:** with the floors disabled the same state must stop the run with its message.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD; Brio–Wu and the negative state $\gamma = 2$, Ryu–Jones 4d $\gamma = 5/3$ ($c_v = 1$, $c_p = \gamma$) | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Divergence control | `none`; Brio–Wu also with GLM (`glm_ch = 4.0`, `glm_alpha = 0.18`, `glm_damping_length = 1.0`, `glm_ch_check = warning`) | `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#divergence-control-mixed-glm) |
| Floors | armed at $10^{-10}$ (clean legs); `p_floor = 1e-3` (floored leg); disabled (stop leg) | `apply_floors`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_cpu_mhd_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_kernels_agnostic.INC); [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | zeroth-order extrapolation | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant) |
| FNL backend | the same kernels on the device | `compute_residuals_weno_dev`, `apply_floors`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]$ along $x$; $y$ and $z$ null |
| Grid | the Sod grid of [V1](../euler/sod), $N = 512$: 64 blocks of $128 \times 4 \times 4$ cells, 4 along $x$ |
| Ranks | 2 |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), discontinuity at $x = 0.5$ |
| Boundaries | `extrapolation` on all six faces |

States $(\rho, u, v, w, p, B_x, B_y, B_z)$ (code units; Brio–Wu from Athena++ `athinput.bw`, Ryu–Jones 4d from the
Athena C test `athinput.rj4d`):

| Leg | $\gamma$ | Left | Right | $t$ |
|---|---|---|---|---|
| Brio–Wu, `none` and GLM | 2 | $(1, 0, 0, 0, 1, 0.75, 1, 0)$ | $(0.125, 0, 0, 0, 0.1, 0.75, -1, 0)$ | 0.1 |
| Ryu–Jones 4d, `none` | 5/3 | $(1, 0, 0, 0, 1, 0.7, 0, 0)$ | $(0.3, 0, 0, 1, 0.2, 0.7, 1, 0)$ | 0.16 |
| negative pressure (floored; stop) | 2 | the Brio–Wu left state | $(0.125, 0, 0, 0, -0.1, 0.75, -1, 0)$ | 0.02 |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/riemann/make_riemann.py` | writes the input from `sod/sod-x.ini` (`--case bw\|rj4d\|negative`, `--cells`, `--divergence-control`, `--rho-floor`, `--p-floor`) |
| `src/tests/flume/verification/mhd/riemann/check.sh` | runs the legs and checks the logs and the final checkpoint (`--np N`, `--numerics SPEC`) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_mhd_riemann` | draws `brio-wu.png` and `rj4d.png` from `<runs>/brio-wu` and `<runs>/rj4d` (`--only mhd-riemann`) |

```bash
cd src/tests/flume/verification/mhd
riemann/check.sh                                                     # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl riemann/check.sh    # FNL
riemann/check.sh --numerics hlld                                     # weno-riemann, HLLD
```

`make_doc_runs.sh <runs> mhd` writes the figure runs (`<runs>/brio-wu`: 512 cells, GLM; `<runs>/rj4d`: 512 cells, no cleaning); `make_doc_figures.py --runs <runs>` draws it.

## Results

![Brio-Wu](/flume/brio-wu.png)

*Brio–Wu shock tube with GLM: density, normal velocity, pressure and the tangential fields at $t = 0.1$.*

![Ryu-Jones 4d](/flume/rj4d.png)

*Ryu–Jones 4d at $t = 0.16$: the switch-on and switch-off structure.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| Brio–Wu (`none`, GLM), Ryu–Jones 4d: stages with floored cells (no `MHD floors:` log line) | 0 | 0 | 0 |
| the same: final $\rho$, $p$ positive and finite | $> 0$ | pass | pass |
| negative pressure, `p_floor = 1e-3`: floors engaged, run completes | engaged | pass | pass |
| negative pressure, floors disabled: the run stops with "cells with a non-positive density or pressure" | stop | pass | pass |

On `weno-riemann` (HLLD, either interpolation) the floored-cell count and the HLLD fallbacks to HLL are both 0
([MHD on `weno-riemann`](./weno-riemann)); the fallbacks are reported, not asserted.

Related cases: [Ryu–Jones 2a](./ryu-jones-2a) has an exact solution on the same grid; [positivity limiter](./positivity-limiter)
removes the floors from the hard problems.

## Findings

- Both problems run clean on the splitting scheme: the eigensystem holds near its degeneracies without floors.
- A negative-pressure state is floored when the floors are armed and stops the run when they are not: the floors never
  act silently on a run that does not ask for them.
- The figures have no exact solution overlaid; the leg asserts cleanliness and positivity, not accuracy.
