# Accuracy across 2:1 seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-8 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | — |

## Case and scope

The diagonal shear wave of [VV-1](./viscous-waves) and the diagonal Ohmic wave of [VV-5](./ohmic-alfven) on a
quadtree with its centre refined 2:1, against the exact linear solution.

What the case catches: conservation says nothing about accuracy, and a seam can conserve ([VV-7](./seam-conservation))
and still degrade the order. The case measures the order of the composite grid, and, with the reflux off, whether the
dissipative flux crossing the seam is the one the register replaces.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Euler with $\mu = 0.01$ (shear); MHD without divergence control (`mhd-none`) with $\eta = 0.01$, $\mu = 0$ (Ohmic); $\gamma = 1.4$ | [models](../../models#compressible-euler), [ideal MHD](../../models#ideal-mhd) |
| Dissipative fluxes | `dissipative_order = 4` | `add_dissipative_fluxes_o4`, `add_resistive_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Coarse ghosts | mean restriction (a coarse ghost is the mean of the fine point values under it) | [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Fine ghosts | the tricubic coarse-to-fine fill | [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Reflux | Berger–Colella reflux of point-value fluxes; off on the negative control | `accumulate_seam_fluxes`, `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree, $4 \times 4$ blocks of $(N_0/4)^2$ cells, $N_0$ = 24 / 48 / 96 at the base level; the blocks whose centroid lies in $[0.25, 0.75]^2$ refined 2:1 (48 / 96 / 192 in the refined box) |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification) at 45°: `shear` mode, $A = 10^{-4}$; `magnetic` mode, $A = 10^{-5}$; $\rho_0 = p_0 = 1$, wavelength $1/\sqrt2$ |
| Boundaries | `periodic` on the $x$ and $y$ faces, `extrapolation` on the $z$ faces |
| Final time | $t = 0.1$ |
| Error measure | RMS against the exact linear solution, weighted by the cell area |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/waves.py` | `make --refine 0.25 0.25 0.75 0.75 --ratio 4 [--reflux false]` writes each input; `oracle` judges the ladder |
| `src/tests/flume/verification/viscous/check.sh --leg vv8` | runs the three ladders (`--np N`; `FLUME_EXE` selects the backend) |
| `make_doc_figures.py` (`fig_viscous_fields`) | the figure, from the CPU work directories of `viscous/check.sh` |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --leg vv8             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg vv8       # FNL
```

## Results

![VV-1 and VV-8 in 2-D: the shear wave and its error on a uniform grid; the shear and Ohmic errors with the centre refined 2:1](/flume/viscous-fields.png)

*Top: the diagonal shear wave and its pointwise error on the uniform $96^2$ grid; bottom: the shear and Ohmic errors on
the $48^2$ grid with its centre refined 2:1 ($96^2$ cells in the box).*

Oracle: order $\ge 1.8$ on the finest pair; the Ohmic ladder without reflux, the negative control, $\le 1.3$.

| Leg | N (finest) | RMS error | Order |
|---|---|---|---|
| shear, $\mu$ | 48 / 96 / 192 | $4.26\cdot10^{-8}$ / $1.09\cdot10^{-8}$ / $2.78\cdot10^{-9}$ | +1.96, +1.98 |
| shear, $\mu = 0$ (probe) | 48 / 96 / 192 | $5.39\cdot10^{-8}$ / $1.36\cdot10^{-8}$ / $3.42\cdot10^{-9}$ | +1.99, +1.99 |
| Ohmic, $\eta$ | 48 / 96 / 192 | $1.15\cdot10^{-9}$ / $2.94\cdot10^{-10}$ / $8.12\cdot10^{-11}$ | +1.97, +1.86 |
| Ohmic, $\eta = 0$ (probe) | 48 / 96 | $7.9\cdot10^{-15}$ / $2.6\cdot10^{-15}$ | round-off |
| Ohmic, no reflux | 48 / 96 / 192 | $2.37\cdot10^{-9}$ / $1.18\cdot10^{-9}$ / $5.82\cdot10^{-10}$ | +1.00, +1.02 |

CPU and FNL identical (2 ranks). The rows without coefficients are probes, not legs of `check.sh`.

Related cases: conservation across the same seams is [VV-7](./seam-conservation); the inviscid order across static
2:1 seams, measured and then lifted in M6, is [order across seams](../seams/order-across-seams).

## Findings

- The composite grid is second order, against fourth on a uniform grid. The rows without coefficients say where that
  comes from. The shear wave is second order without viscosity too: its seam error is the inviscid one. With the mean
  restriction and the Berger–Colella reflux of point-value fluxes (a conservative but $O(H^2)$-inconsistent coarse
  flux), the composite error is second order ([issue #68](https://github.com/szaghi/adam/issues/68) F1: each of the two
  caps the order at 2 on its own; the #21 bound on restriction-compatible prolongations does not apply to the tricubic
  fill).
- A static field has no inviscid flux at first order, so the Ohmic composite at $\eta = 0$ stays at round-off, and the
  Ohmic composite error is the dissipative flux crossing the seam: second order with reflux, first without (the
  gradient of mean-restricted ghosts is an $O(1)$ error, which the reflux replaces by the fine flux).
- The plan in #65 asked for the seam error "within a stated factor of the uniform fine run"; between a second-order
  seam and a fourth-order interior that factor grows with $N$, so the leg asserts the order instead.
- The error maps show where the seam error lives. On the uniform grid the error of the shear wave is smooth and spread
  over the domain, near $10^{-13}$. On the refined grid the Ohmic error is a ring on the 2:1 faces, two orders above
  the rest of the domain: the dissipative flux at the face. The shear error fills the refined box and its surroundings
  instead, being the inviscid seam error that the flow carries and the viscosity spreads.
