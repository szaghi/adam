# Lax problem

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| RV-3 (`lax` leg) | Euler | M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno-riemann` HLLC) | yes | — |

## Case and scope

The Riemann problem of Lax (1954): $(\rho, u, p) = (0.445, 0.698, 3.528)$ on the left, $(0.5, 0, 0.571)$ on the right.
Its contact is stronger than Sod's, so the dissipation of the face flux shows there. The exact solution is known, and the
error is measured against it.

What the case catches:
- **a dissipative Riemann solver or a degraded flux correction** of `scheme_space = weno-riemann`: $L_1(\rho)$ must stay
  within $1.05\times$ the splitting scheme's (owner decision, 2026-09-30);
- **direction-dependent defects:** the problem runs along $x$, $y$ and $z$, and the three runs must be equal **bit for
  bit** after the permutation of the axes;
- **broken 1-D consistency:** every transverse copy of the solution, in every block, must be bitwise identical;
- **an oracle that ignores the velocity:** unlike Sod, the left state moves, so the exact solution must read the normal
  velocity of the input.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = c_p/c_v = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` (baseline) | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Space, `weno-riemann` | WENO5 interpolation of the face states (characteristic), an LLF, HLL or HLLC Riemann flux, the 6th-order flux correction with the WENO-weight sensor | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_riemann_face_fluxes`, [`adam_flume_cpu_riemann_face_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_riemann_face_kernels_agnostic.INC); `compute_riemann_llf`, `compute_riemann_hll`, `compute_riemann_hllc`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | zeroth-order extrapolation | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

The two schemes are described in [numerics](../../numerics#scheme-space-weno-flux-splitting) and
[numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction).

## Initial and boundary conditions

The case is the [Sod](./sod) input with the two region states and the final time rewritten by `riemann-flux/check.sh`
(`set_lax`).

| Item | Value |
|---|---|
| Domain | $[0, 1]$ along the active axis; the other two axes null |
| Grid | 200 cells: 64 blocks of $50 \times 4 \times 4$ cells (4 along the active axis), one uniform level, `ngc = 3` |
| Ranks | 2 |
| Left state, $x \le 0.5$ | $(\rho, u_n, p) = (0.445, 0.698, 3.528)$ |
| Right state, $x > 0.5$ | $(\rho, u_n, p) = (0.5, 0, 0.571)$ |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), two regions |
| Boundaries | `extrapolation` on all six faces |
| Final time | $t = 0.13$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/sod/sod-{x,y,z}.ini` | the base inputs, rewritten into the Lax problem |
| `src/tests/flume/verification/riemann-flux/check.sh --leg lax` | rewrites the inputs (`set_lax`, `set_numerics`), runs the three directions; `--solver split\|llf\|hll\|hllc`, `--correction`, `--sensor`, `--lax-l1-max` |
| `src/tests/flume/verification/sod/sod_oracle.py` | the exact solution (normal velocity read from the input), the $L_1$ bound, the bitwise checks |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_lax`) | the figure, from the `split`, `llf` and `hllc` work directories of `riemann-flux/` |

```bash
cd src/tests/flume/verification/riemann-flux
./check.sh --leg lax                                    # HLLC (default), bound 1.12e-2
./check.sh --leg lax --solver split                     # the splitting baseline
./check.sh --leg lax --solver llf --lax-l1-max X        # LLF needs an explicit bound X
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg lax   # FNL
```

## Results

![Lax problem](/flume/lax.png)

*Density at $t = 0.13$, 200 cells: the splitting scheme, `weno-riemann` LLF and HLLC against the exact solution.*

| Scheme | Bound on $L_1(\rho)$ | CPU | FNL |
|---|---|---|---|
| `weno` (split) | — | $1.064\cdot10^{-2}$ | not run |
| `weno-riemann` LLF | explicit (exceeds the default) | $1.214\cdot10^{-2}$ | not run |
| `weno-riemann` HLL | — | $1.125\cdot10^{-2}$ | not run |
| `weno-riemann` HLLC (default) | $1.12\cdot10^{-2}$ ($1.05\times$ split) | $1.081\cdot10^{-2}$ | the same |
| $y$ and $z$ runs against $x$ | bitwise | bitwise | bitwise |

HLLC by flux correction (CPU):

| Correction | $L_1(\rho)$ |
|---|---|
| 6th, WENO-weight sensor (default) | $1.081\cdot10^{-2}$ |
| 6th, no sensor | $1.146\cdot10^{-2}$ |
| 4th | $1.077\cdot10^{-2}$ |
| none | $1.084\cdot10^{-2}$ |

Related cases: [Sod](./sod) and [Shu–Osher](./shu-osher) compare the same solvers; the scheme study behind the
defaults is in `src/tests/flume/verification/riemann-flux/README.md`.

## Findings

- Lax is the one problem where the splitting scheme beats every Riemann solver: HLLC is 1.6% worse. The correction is
  not the cause (1.8% worse without it). The gap belongs to interpolating the face states and solving the Riemann
  problem, on a problem measured almost entirely at its discontinuities.
- The WENO-weight sensor is needed: without it HLLC is 6% worse.
- LLF exceeds the default bound and is run with an explicit one.
- The first Lax measurement was wrong: `sod_oracle.py` took both velocities as zero (Sod only). It now reads the normal
  velocity of the input; the Sod results are unchanged.
- Only the default HLLC run exists on FNL; the other schemes were run on the CPU only.
