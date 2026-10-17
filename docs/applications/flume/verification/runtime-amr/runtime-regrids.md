# Runtime regrids in FLUME

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| RG (rg3 to rg6) | Runtime AMR (M5) | M5 ([#74](https://github.com/szaghi/adam/issues/74)) | CPU, FNL | yes | `vortex-regrid` |

## Case and scope

FLUME regrids for real during the run, on cases where the grid keeps changing: a vortex travelling across a periodic
quadtree, and a shock over an immersed cylinder.

What the case catches:
- **a regrid that changes what the scheme conserves:** a drift of a volume integral through a prolongation, a
  restriction, the redistribution or the flux register rebuilt at every regrid, with seams that appear and disappear;
- **a run that is not reproducible across a restart:** the restart must hold the grid the next step runs on;
- **a grid-dependent piece out of step with the grid:** the flux register, the solids' distance function;
- **a device state that does not survive the host round trip** of the FNL regrid.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space | WENO5 flux splitting, characteristic variables; reflux on | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | vortex: `runge-kutta-ssp-54`, CFL 0.4; cylinder: `runge-kutta-ssp-33`, CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Regrid step | every `[amr] frequency` steps, up to `[amr] iters` sweeps: ghosts refreshed, markers combined, the library regrid, the distance function recomputed on a changed grid; then the new state checked (density and pressure positive) | `regrid_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mark_runtime`, `check_regridded_state`, `is_regrid_step`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `amr_update`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90) |
| Markers | Löhner estimator; solid marker (blocks crossed by the surface) | `mark_by_lohner`, `mark_by_solid`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Prolongation | `conservative` (default) or `linear` | `refine_conservative`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); [round trip](./regrid-round-trip) |
| Forest | regrids after the step closes, before the per-step output; rebuilds the 2:1 flux register on a change | `regrid`, `share_flux_register`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90) |
| FNL regrid | a host round trip: the state copied to the host, marked and regridded there, the blocks in use copied back packed and transposed | `regrid_forest`, `copy_blocks_transposed`, `copy_cpu_gpu`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Leg | Case | Grid and AMR | Run |
|---|---|---|---|
| rg3 | the [isentropic vortex](../../initial-conditions#isentropic-vortex) of V2 (strength 5, radius 0.07, centre $(0.5, 0.5)$, $\rho = p = 1$, $(u, v) = (1, 1)$), `periodic` in $x$ and $y$ | quadtree on $[0, 1]^2$, $8 \times 8$ cells per block, base level 1 ($2 \times 2$ blocks), `max_level` 3, `iters` 2; Löhner marker on the density (`refine_tol` 0.3, `derefine_tol` 0.1, buffer 2), a regrid every 5 steps | 200 steps, 2 ranks, `conservative` and `linear` |
| rg4 | the rg3 conservative case | as rg3 | 30 steps, against 20 steps, a restart saved at step 20, and 10 more |
| rg5 | the [shock over the cylinder](../euler/shock-cylinder) (V6): inflow $(\rho, u, p) = (3.7333, 1.25, 4.5)$ into $(1.4, 0, 1)$, the shock at $x = 0.2$, the cylinder of radius 0.13 at $(0.5, 0.5)$; `inflow` at $x = 0$, `extrapolation` at $x = 1$, `wall-inviscid` on the $y$ faces | octree on $[0, 1]^2$ ($z$ null), $16 \times 16 \times 4$ cells per block, base level 2, `max_level` 3; the solid marker and a Löhner marker on the density (`refine_tol` 0.6, `derefine_tol` 0.2), a regrid every 5 steps | 60 steps, 2 ranks |
| rg6 | the rg3 conservative runs of both backends | as rg3 | — |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/regrid/make_regrid.py` | `vortex` and `cylinder` write the inputs |
| `src/tests/flume/verification/regrid/check.sh --legs rg3,rg4,rg5,rg6` | runs the legs (`FLUME_EXE` selects the backend; rg6, not a default leg, compares the work directories of both) |
| `src/tests/flume/verification/make_doc_runs.sh <runs> regrid`, `make_doc_figures.py --only regrid-runtime` | re-runs rg3, conservative and linear, with a checkpoint every 5 steps (2 runs, about 2 min) |
| `src/tests/flume/verification/conservation/conservation_oracle.py` | the drift bounds (rg3), the bitwise and tolerance comparisons (rg4, rg6) |
| `src/tests/flume/regression/vortex-regrid/input.ini` | the frozen rg3 conservative case to 100 steps: six regrids that refine and coarsen |

```bash
cd src/tests/flume/verification/regrid && ./check.sh --legs rg3,rg4,rg5                         # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs rg3,rg4,rg5                 # FNL
./check.sh --legs rg6                                                                           # CPU against FNL
```

## Results

![Runtime regridding of the travelling vortex](/flume/regrid-runtime.png)

*The rg3 vortex on a periodic quadtree (levels 1–3, a Löhner marker on the density, a regrid every 5 steps), CPU: the density and the blocks at steps 0, 100 and 200 (4, 43 and 52 blocks; the grid refines ahead of the vortex and coarsens behind it), and the relative drift of the mass and energy integrals per step. With the conservative prolongation they stay at round-off through 40 regrids; with the linear one each regrid that moves data shifts them (to $3\cdot10^{-5}$), the negative control.*

| Leg | Asserted |
|---|---|
| rg3 | the grid refines and coarsens during the run; with `conservative` the five volume integrals constant within $10^{-13}$; with `linear` a drift above $10^{-10}$ (the negative control) |
| rg4 | a regrid at step 20 that changes the grid; the last fields bitwise and the histories identical |
| rg5 | the run completes with a regrid that changes the grid, each regrid checking the new state (density and pressure positive) |
| rg6 | the same regrids at the same steps into the same grids (the regrid log lines equal) and the same fields within $10^{-10}$, relative to each variable's largest magnitude |

rg3, rg4 and rg5 run on both backends (np 2); the measured values are the same on CPU and FNL except where noted:

| Leg | Regrids that changed the grid | Blocks refined / families coarsened | Measured |
|---|---|---|---|
| rg3 `conservative` | 10 | 20 / 4 | drift $2.2\cdot10^{-16}$ (ρ), $1.3\cdot10^{-16}$ (ρE); FNL at most $2.6\cdot10^{-16}$ |
| rg3 `linear` | 7 | 14 / 6 | drift $2.6\cdot10^{-5}$ (ρ), $2.2\cdot10^{-4}$ (ρu) |
| rg4 | 3 (steps 5, 20, 25) | — | 40 blocks bitwise, histories identical |
| rg5 | 1 (176 → 344 blocks) | 24 / 0 | completes, admissible |
| rg6 | 10, the same on both | 20 / 4 | fields within $7.0\cdot10^{-14}$; FNL round trips 0.015 to 0.029 s each (2.3 to 3.1 s before #75) |

Related cases: the transfers alone are the [round trip](./regrid-round-trip); the accuracy of the regridded runs is
[AV](./accuracy) (rg5's accuracy against the initial-AMR and uniform runs is AV-7, #74 P4); the cases without runtime
AMR are the [isentropic vortex](../euler/isentropic-vortex) and the [shock over the cylinder](../euler/shock-cylinder).

## Findings

- Nothing crosses the periodic boundary of rg3, so a drift can only come from a regrid or a seam: the leg covers the
  prolongation, the restriction, the redistribution and the flux register rebuilt at every regrid, with seams that
  appear and disappear.
- The two rg3 runs regrid differently: the linear prolongation changes the data the Löhner estimator reads, so the
  grids part after the first regrid. The integrals still separate the two by eleven orders of magnitude.
- The forest regrids before the per-step output, so the restart holds the grid the next step runs on; saved before the
  regrid, it would not.
- The FNL backend regrids on the host, so the grid decisions are the same code on both backends; rg6 checks that the
  device state reaches the host and comes back intact.
- On FNL the round trips of the conservative rg3 run cost 0.015–0.031 s each on a grid of at most 52 blocks: the
  copies move the blocks in use, packed on the device. They moved the whole device state, sized by the block capacity
  (17195 blocks per rank here), and cost 2.3–3.1 s each until [#75](https://github.com/szaghi/adam/issues/75).
- P0 also fixed the capacity check of a regrid: the new blocks take indices on their parent's rank before the
  redistribution, so the bound is per rank (`nb`), not the former `procs_number · nb` over every node of the
  replicated tree, which let one rank overflow its arrays undetected (`check_blocks_number`).
