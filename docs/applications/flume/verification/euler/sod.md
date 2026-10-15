# Sod shock tube

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| V1, RV-3 | Euler | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes | `sod-x`, `sod-y`, `sod-z`, `sod-hllc` |

## Case and scope

The shock tube of Sod (1978): two gases at rest, separated by a diaphragm at $x = 0.5$, high pressure on the left. At
$t > 0$ the solution is a left-running rarefaction, a contact discontinuity and a right-running shock. The exact
solution is known (Toro 2009, §4.2–4.5), so the error of the scheme is measured, not judged by eye.

What the case catches:
- **direction-dependent defects:** the problem is run along $x$, $y$ and $z$ (the other two axes null). CHASE, the
  application FLUME replaced, only ever ran x-aligned problems, and its transposed characteristic projection and
  singular y/z eigenvectors went unnoticed. Here the three runs must be equal **bit for bit** after the permutation of
  the axes;
- **loss of accuracy** of the shock-capturing scheme: $L_1(\rho)$ against the exact solution under a recorded bound;
- **broken 1-D consistency:** every transverse copy of the solution, in every block, must be bitwise identical;
- **wall boundaries** (the `--wall` variant): a double Sod in a closed box must stay mirror symmetric and keep its mass
  and energy.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = c_p/c_v = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`reconstruction_variables = characteristic`) | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Characteristic projection | left and right eigenvectors of the flux Jacobian in direction $d$ | `compute_eigenvectors`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Space, `weno-riemann` (RV-3) | WENO interpolation of the face states, then an LLF, HLL or HLLC Riemann solver | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_euler_hllc_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_euler_hllc_kernels.F90) (and the `_llf`, `_hll` twins) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | zeroth-order extrapolation; reflecting inviscid walls in the `--wall` variant | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

The scheme options are described in [numerics](../../numerics).

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]$ along the active axis; the other two axes null |
| Grid | 200 cells: 64 blocks of $50 \times 4 \times 4$ cells (4 along the active axis), one uniform level, `ngc = 3` |
| Ranks | 2 (the 8 sibling groups split 32/32) |
| Left state, $x \le 0.5$ | $(\rho, u, p) = (1, 0, 1)$ |
| Right state, $x > 0.5$ | $(\rho, u, p) = (0.125, 0, 0.1)$ |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), two regions |
| Boundaries | `extrapolation` on all six faces |
| Final time | $t = 0.2$ |
| `--wall` variant | walls at both ends, the high state in the middle, $t = 0.4$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/sod/sod-{x,y,z}.ini` | the three directions |
| `src/tests/flume/verification/sod/sod-wall-{x,y,z}.ini` | the reflecting-wall double Sod |
| `src/tests/flume/verification/sod/check.sh` | runs the cases (`--wall`, `--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/sod/sod_oracle.py` | the exact solution, the $L_1$ bound, the bitwise and mirror checks |
| `src/tests/flume/verification/riemann-flux/check.sh --leg sod` | the `weno-riemann` runs (RV-3) |

```bash
cd src/tests/flume/verification/sod && ./check.sh             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh   # FNL
./check.sh --wall                                             # the reflecting-wall variant
```

The figure is drawn by `make_doc_figures.py` (`fig_sod`) from the work directories of the two scripts.

## Results

![Sod shock tube, 200 cells, t = 0.2](/flume/sod.png)

*Density, velocity and pressure at $t = 0.2$: the splitting scheme and `weno-riemann` HLLC against the exact solution.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| $L_1(\rho)$, `weno` (split, characteristic) | $3.31\cdot10^{-3}$ (measured + 2%) | $3.244\cdot10^{-3}$ | $3.244\cdot10^{-3}$ (CPU − FNL: $7.2\cdot10^{-14}$) |
| $L_1(\rho)$, `weno-riemann` LLF / HLL / HLLC | — | $3.426$ / $3.013$ / $2.942\cdot10^{-3}$ | the same |
| $y$ and $z$ runs against $x$ | bitwise | bitwise | bitwise |
| transverse copies of the 1-D solution | bitwise | bitwise | bitwise |
| `--wall`: mirror symmetry about the centre | $10^{-12}$ | pass | pass |
| `--wall`: relative drift of mass and energy | $10^{-11}$ | pass | pass |

Related cases: [Lax](./lax) and [Shu–Osher](./shu-osher) compare the same Riemann solvers on harder problems; the MHD
analogue is [Ryu–Jones 2a](../mhd/ryu-jones-2a).

## Findings

- The three directions are the same arithmetic on permuted data (cyclic tangents), so they agree to the last bit: the
  direction-dependent eigenvector defects of CHASE cannot recur unnoticed.
- The CPU and FNL backends agree to $7.2\cdot10^{-14}$, not bitwise: the compilers contract differently into FMAs. They
  agree that closely because the host and device WENO weights share one exponent (`weno_weights_exponent`).
- HLLC is the most accurate of the Riemann solvers on the contact, and slightly more accurate than the splitting
  scheme.
- Symmetry alone does not prove that the walls reflect: extrapolation at both ends is symmetric too. The `--wall`
  variant therefore also requires a closed box (extrapolated ends would leak about 10% of the mass).
