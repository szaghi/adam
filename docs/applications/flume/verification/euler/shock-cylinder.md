# Shock over a cylinder

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| V6 | Euler | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU, FNL | yes | `shock-cylinder-ib` |

## Case and scope

A Mach 2 planar shock runs into a cylinder imposed by the immersed boundary, in a channel with slip walls. The flow has
no closed form, but the setup is mirror symmetric about $y = 0.5$, and so is the exact solution. The oracle asserts what
the setup guarantees: the refinement, the symmetry and the positivity of the fluid.

The case exercises every immersed-boundary piece together: the distance function, the eikonal extrapolation into the
solid (Jacobi sweeps) and the Euler-wall inversion, the cut spacing of the fluid cells next to the surface, the solid
masks of the Runge–Kutta stages, and the solid AMR marker with the coarse–fine reflux of the refined ring.

What the case catches:
- **any defect that breaks the symmetry:** indexing and ordering errors in the distance function, the eikonal
  extrapolation, the cut spacing or the solid marker, stencils reading the wrong ghost, data races. A wrong wall law
  that is itself symmetric is not caught;
- **a solid marker that misses blocks:** every block crossed by the surface must be refined;
- **loss of positivity** next to the surface.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space, `weno-riemann` | WENO5 interpolation (characteristic), HLLC, the 6th-order flux correction with the WENO-weight sensor | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_riemann_hllc`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Immersed boundary | analytical circle, Euler wall (`bc_type = 2`), `n_eikonal = 2` Jacobi sweeps | `compute_phi_analytical_circle`, `evolve_eikonal`, `invert_eikonal`, [`adam_ib_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_ib_object.F90); [numerics](../../numerics#immersed-boundary-euler) |
| AMR | init-time solid marker (blocks crossed by the surface refined once), reflux at the 2:1 faces | `mark_by_solid`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90); `accumulate_coarse_flux`, `accumulate_fine_flux`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Boundaries | inflow at $x_{\min}$, extrapolation at $x_{\max}$, inviscid walls in $y$ | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#wall-inviscid-slip-wall-by-mirroring-with-a-sign-vector) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | octree (the 2:1 reflux restriction needs two tangential axes), 64 blocks of $16 \times 16 \times 4$ cells (base spacing $1/64$), the blocks crossed by the cylinder surface refined once (`delta_fine = 0.0078125`); measured: 176 blocks, 96 crossed by the surface, 32 coarse–fine faces |
| Ranks | 2 |
| Cylinder | centre $(0.5, 0.5)$, radius 0.13, axis $z$; its surface crosses the block faces $x = 0.375$ and $0.625$ at a shallow angle |
| Post-shock state, $x \le 0.2$ and inflow | $(\rho, u, p) = (3.7333, 1.25, 4.5)$ (Mach 2 shock) |
| Pre-shock state, $x > 0.2$ | $(\rho, u, p) = (1.4, 0, 1)$ (sound speed 1) |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), two regions, `amr_iterations = 1` |
| Boundaries | $x_{\min}$ `inflow` (post-shock state), $x_{\max}$ `extrapolation`, $y$ `wall-inviscid`, $z$ `extrapolation` |
| Final time | $t = 0.25$ (the shock has passed the cylinder; 270 steps) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/shock-cylinder/shock-cylinder.ini` | the case (`[solids]`, `[solid_1] definition = analytical_circle`, an `[amr]` solid marker) |
| `src/tests/flume/verification/shock-cylinder/check.sh` | V6 (`--np N`; `FLUME_EXE` selects the backend), mirror tolerance $10^{-10}$ |
| `src/tests/flume/verification/shock-cylinder/shock_cylinder_oracle.py` | the refinement, mirror-symmetry and fluid-sanity checks; `--compare` for two runs |
| `src/tests/flume/verification/riemann-flux/check.sh --leg cylinder` | the same case on `weno-riemann` |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_cylinder`) | the figure, from the HLLC work directory of `riemann-flux/` |

```bash
cd src/tests/flume/verification/shock-cylinder && ./check.sh                  # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                   # FNL
cd ../riemann-flux && ./check.sh --leg cylinder                               # weno-riemann HLLC
```

## Results

![Shock over a cylinder](/flume/shock-cylinder.png)

*Density with the block outlines (the refined ring on the surface) and Mach number at $t = 0.25$, `weno-riemann` HLLC.*

| Check | Bound | `weno`, CPU | `weno`, FNL | HLLC, CPU | HLLC, FNL |
|---|---|---|---|---|---|
| crossed blocks refined | all | pass | pass | pass | pass |
| relative asymmetry about $y = 0.5$ | $10^{-10}$ | $1.2\cdot10^{-11}$ | $4.8\cdot10^{-12}$ | $4.9\cdot10^{-12}$ | $1.9\cdot10^{-12}$ |
| min $\rho$, min $p$ in the fluid | positive, finite | positive | positive | 1.365, 0.968 | positive |
| CPU against FNL, relative, every field | — | $1.1\cdot10^{-11}$ | | | |

Related cases: [Conservation across AMR seams](./amr-conservation) checks the reflux this case relies on; the
regression case `shock-cylinder-ib` freezes the input at 120 steps.

## Findings

- The radius 0.13 makes the surface cross block faces at a shallow angle, so fluid and solid cells face each other
  across block faces. The case therefore depends on the all-solids $\phi$ summary in the ghost cells: with the former
  interior-only summary the solution moves by $3.9\cdot10^{-2}$.
- Both flux paths keep the symmetry to round-off and the fluid positive; HLLC is the more symmetric of the two.
- The symmetry check is blind to a wrong but symmetric wall law; the case verifies the immersed-boundary plumbing, not
  its accuracy.
