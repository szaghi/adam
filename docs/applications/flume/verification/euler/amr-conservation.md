# Conservation across AMR seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| V3, RV-4 | Euler | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes | `amr-periodic-reflux`, `amr-periodic-reflux-stage` |

## Case and scope

A finite-volume update conserves only if every interior face flux leaves one cell and enters its neighbour. At a 2:1
coarse–fine face the two sides compute different fluxes (different stencils, different ghost fills), so the grid leaks
mass, momentum and energy unless the Berger–Colella reflux replaces the coarse flux with the restricted fine one. The
case is a triply periodic box with one refined octant: nothing crosses the boundary, so the five volume integrals must
stay constant to round-off.

What the case catches:
- **a missing or wrong reflux** at any of the 6 coarse–fine faces, periodic ones included;
- **a reflux register that does not hold the flux the step used:** FLUME accumulates the seam fluxes of every
  Runge–Kutta stage weighted by its SSP coefficient;
- **a vacuous test:** the run without reflux is the negative control and must drift well above round-off, otherwise the
  case does not exercise the seams;
- **block layouts the 2:1 seams cannot handle**, refused at initialization: an odd cell count along a refined axis
  ([#39](https://github.com/szaghi/adam/issues/39); a 2:1 child boundary in the middle of a parent cell gave NaN
  silently), and a count below $2\,n_{gc}$ ([#66](https://github.com/szaghi/adam/issues/66); the outer coarse ghost
  layers would be restricted from the fine block's own, stale ghosts).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space, `weno-riemann` (RV-4) | WENO5 interpolation (characteristic), HLLC, the 6th-order flux correction with the WENO-weight sensor | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR | geometric box marker, one more level on the octant; 2:1 ghost restriction and interpolation | `mark_by_geometry`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Reflux | Berger–Colella, stage fluxes weighted by the SSP coefficients | `accumulate_seam_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `accumulate_coarse_flux`, `accumulate_fine_flux`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90); `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Block-cell checks | odd counts and counts below $2\,n_{gc}$ along a refined axis refused | `check_amr_block_cells`, [`adam_realm_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_realm_object.F90) |
| Boundaries | periodic on all six faces | [boundary conditions](../../boundary-conditions#periodic) |
| Integrals | volume integrals of the conservative variables, every step | `compute_conservation`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^3$ |
| Grid | octree, one uniform level (8 blocks of $8^3$ cells), then the octant $[0, 0.5]^3$ refined once more: 15 blocks, 6 coarse–fine faces (the log must report `registered intra-realm AMR seam faces: +6`), `ngc = 3` |
| Ranks | 2 |
| State | uniform flow $(\rho, u, v, w, p) = (1, 1, 1, 1, 1)$ with a 1% seeded density and pressure perturbation (`s = 0.01`), so every seam carries a flux |
| Initial condition | [`uniform`](../../initial-conditions#uniform-seeded-uniform-state) |
| Boundaries | `periodic` on all six faces |
| Final time | 20 steps |

| Leg | Change | Expected |
|---|---|---|
| reflux on | — | the five integrals constant within $10^{-13}$ |
| reflux off | `reflux = .false.` | drift of at least $10^{-10}$ (negative control) |
| odd cells | `ni = 7` | refused at initialization (`[grid].(ni)=+7 is odd: 2:1 refinement`) |
| thin blocks | `ni = 4` ($< 2\,n_{gc}$) | refused at initialization (`[grid].(ni)=+4 is less than 2*ngc`) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/conservation/amr-periodic.ini` | the case |
| `src/tests/flume/verification/conservation/check.sh` | V3: the four legs (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/conservation/conservation_oracle.py` | the drift of the integrals (`--conserved`, `--leaky`); `--compare` for two runs cell by cell |
| `src/tests/flume/verification/riemann-flux/check.sh --leg conservation` | RV-4: the reflux on and off legs on `weno-riemann` |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_conservation`) | the figure, from the conservation histories of `conservation/` |

```bash
cd src/tests/flume/verification/conservation && ./check.sh               # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh              # FNL
cd ../riemann-flux && ./check.sh --leg conservation                      # RV-4, HLLC
```

## Results

![Conservation across AMR seams](/flume/conservation.png)

*Relative drift of the five volume integrals against time, reflux on (solid) and off (dashed).*

| Run | Bound | CPU | FNL |
|---|---|---|---|
| reflux on, `weno` | $10^{-13}$ | $\le 6\cdot10^{-15}$ | $2.2\cdot10^{-16}$ |
| reflux on, `weno-riemann` HLLC | $10^{-13}$ | 0 | $2.2\cdot10^{-16}$ |
| reflux off | $\ge 10^{-10}$ | $1$–$2\cdot10^{-5}$ | $1$–$2\cdot10^{-5}$ |
| `ni = 7` | refused | refused | refused |
| `ni = 4` | refused | refused | not run |

Related cases: [Shock over a cylinder](./shock-cylinder) and [Multi-realm runs](./multi-realm) rely on the same reflux;
the dissipative fluxes across 2:1 seams are checked in [VV-7](../dissipative/seam-conservation).

## Findings

- With the stage fluxes weighted by the SSP coefficients the register holds the flux the committed step used, and the
  integrals stay at round-off on 1 and 2 ranks (measured at M1-P5).
- Without the reflux the integrals drift by $1$–$2\cdot10^{-5}$ in 20 steps: the seams carry flux, and the reflux leg is
  not vacuous.
- The MHD version (MV-11, `mhd/conservation`) adds the three components of $\mathbf{B}$ and $\psi$: drift
  $\le 4.7\cdot10^{-14}$.
- The two refused layouts were silent before their checks: NaN for the odd count, rank-dependent results for the thin
  block (measured 0.2 in density between 1 and 2 ranks on the step forest, see
  [Forward-facing step](./forward-facing-step)).
- The `ni = 4` leg has no FNL run on record (no FNL work directory).
