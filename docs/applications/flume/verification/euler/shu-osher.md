# Shu–Osher shock / entropy-wave interaction

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| RV-3 (`shu-osher` leg) | Euler | M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno-riemann` HLLC) | yes | — |

## Case and scope

A Mach 3 shock runs into a sinusoidal density field (Shu & Osher 1989, *J. Comput. Phys.* 83, example 8). Behind the
shock the density wave is compressed into high-frequency entropy waves, which a dissipative scheme smears. The problem
has no closed-form solution; the reference is a much finer run.

What the case catches:
- **excess dissipation** of the `weno-riemann` scheme at the entropy waves: $L_1(\rho)$ against the reference must stay
  within $1.05\times$ the splitting scheme's (owner decision, 2026-09-30);
- **a flux-correction sensor that trips on smooth waves:** a jump sensor switches the correction off in the entropy
  waves and does worse than splitting;
- **loss of positivity** behind the strong shock: density and pressure must stay positive;
- **direction-dependent defects:** the $x$, $y$ and $z$ runs must be bitwise equal after the permutation of the axes,
  and every transverse copy bitwise identical.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` (baseline) | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Space, `weno-riemann` | WENO5 interpolation of the face states (characteristic), LLF, HLL or HLLC, the 6th-order flux correction with the WENO-weight sensor | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_riemann_llf`, `compute_riemann_hll`, `compute_riemann_hllc`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Initial condition | the `shu-osher` kind | `set_initial_conditions`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| Boundaries | zeroth-order extrapolation | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Reference | the splitting scheme in NumPy (`hybrid_proto.py`, FLUME's `weno` scheme) on 16× finer cells, averaged to the grid | `reference_density`, `src/tests/flume/verification/riemann-flux/shu_osher_oracle.py` |

## Initial and boundary conditions

The case is the [Sod](./sod) input rewritten by `riemann-flux/check.sh` (`set_shu_osher`).

| Item | Value |
|---|---|
| Domain | $[-5, 5]$ along the active axis; the other two axes null |
| Grid | 400 cells: 64 blocks of $100 \times 4 \times 4$ cells (4 along the active axis), one uniform level, `ngc = 3` |
| Ranks | 2 |
| Left state, $x \le -4$ | $(\rho, u_n, p) = (3.857143, 2.629369, 10.33333)$ |
| Right state, $x > -4$ | $\rho = 1 + 0.2 \sin 5x$, $u = 0$, $p = 1$ |
| Initial condition | [`shu-osher`](../../initial-conditions#shu-osher-shock-and-density-wave-interaction) (`interface = -4`, `rho_amplitude = 0.2`, `rho_wavenumber = 5`) |
| Boundaries | `extrapolation` on all six faces |
| Final time | $t = 1.8$ |
| Reference | 6400 cells, averaged to the 400 cells of the run; max $\rho = 4.672$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/sod/sod-{x,y,z}.ini` | the base inputs, rewritten into the Shu–Osher problem |
| `src/tests/flume/verification/riemann-flux/check.sh --leg shu-osher` | rewrites the inputs (`set_shu_osher`, `set_numerics`), runs the three directions; `--solver`, `--so-l1-max` |
| `src/tests/flume/verification/riemann-flux/shu_osher_oracle.py` | the reference (computed once, a few minutes, cached in the temporary directory), the $L_1$ bound, positivity, the bitwise checks |
| `src/tests/flume/verification/riemann-flux/hybrid_proto.py` | the NumPy splitting scheme that computes the reference |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_shu_osher`) | the figure, from the `split` and `hllc` work directories and the cached reference |

```bash
cd src/tests/flume/verification/riemann-flux
./check.sh --leg shu-osher                                         # HLLC (default), bound 0.274
./check.sh --leg shu-osher --solver split                          # the splitting baseline
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg shu-osher   # FNL
```

## Results

![Shu-Osher](/flume/shu-osher.png)

*Density at $t = 1.8$, 400 cells: the splitting scheme and `weno-riemann` HLLC against the 6400-cell reference; right,
the entropy waves behind the shock.*

| Scheme | Bound on $L_1(\rho)$ | $L_1(\rho)$, CPU | max $\rho$ (reference 4.672), CPU | FNL |
|---|---|---|---|---|
| `weno` (split) | — | 0.2605 | 4.563 | not run |
| `weno-riemann` LLF | explicit (exceeds the default) | 0.3009 | 4.550 | not run |
| `weno-riemann` HLL | — | 0.2490 | 4.590 | not run |
| `weno-riemann` HLLC (default) | 0.274 ($1.05\times$ split) | **0.2425** | 4.592 | the same |
| $x$, $y$, $z$ runs | bitwise | bitwise | | bitwise |

Related cases: [Sod](./sod) and [Lax](./lax) compare the same solvers on simpler Riemann problems.

## Findings

- HLLC resolves the entropy waves best: 7% below the splitting scheme's $L_1(\rho)$, and the highest density peak.
  LLF is the most dissipative and exceeds the default bound.
- The reference is the splitting scheme itself, so the comparison is biased toward it by the reference's own error,
  about 1/16 of the coarse error at the shocks. HLLC wins despite that bias.
- In the NumPy prototype study (`riemann-flux/README.md`) the jump sensor switched the correction off in the entropy
  waves and did worse than splitting; the WENO-weight sensor is the default for that reason.
- Only the default HLLC run exists on FNL.
