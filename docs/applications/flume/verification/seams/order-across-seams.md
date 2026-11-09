# Order across static 2:1 seams

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| AO (A1, C1, D1, D2) | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) · A1, C1, D1, D2 | P | CPU, FNL | yes | — |

## Case and scope

A smooth vortex starts centred on the corner of a refined quadrant, so its core straddles an $x$ seam, a $y$ seam and
their corner, and the free stream $(1, 1)$ carries it into the fine quadrant. The exact solution is the convected
vortex, so the error is measured cell by cell, split into the seam band and the interior.

Every smooth order test ([V2](../euler/isentropic-vortex), MV-5 to MV-7 in
[order of accuracy](../mhd/order-of-accuracy)) runs on a uniform grid; [VV-8](../dissipative/seam-accuracy) measures
a seam with dissipative waves. No test measured the inviscid order of a smooth feature crossing a seam, where the #68
audit (F1) predicts that the mean restriction of the coarse ghosts and the Berger–Colella reflux of point-value fluxes
each cap the order at 2.

What the case catches:
- **a loss of order at the seams**, in the band and in the regions the flow reaches after crossing a seam;
- **a degradation of the seam operators**: the finest refined errors are bounded by the CPU baseline plus 2 %, so an
  improvement of the transfers or of the reflux passes and a degradation fails;
- **a broken instrument**: the uniform control (no refined quadrant) must reach design order.

Phases C (high-order restriction), D (6-point coarse-to-fine fills) and E (reflux policy) of #68 are measured against
these baselines; the legs `euler-l6`, `euler-w6`, `euler-q6`, `euler-ww`, `euler-mp` and `euler-tq` are their entries in the
same ladder.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations, `euler` | compressible Euler, ideal gas, $\gamma = c_p/c_v = 1.4$ | [models](../../models#compressible-euler) |
| Equations, `mhd` | ideal MHD with mixed GLM cleaning, $\gamma = 5/3$, `glm_ch` 3.0, `glm_alpha` 0.18, `glm_damping_length` 1.0 | [models](../../models#divergence-control-mixed-glm); `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Fine-to-coarse ghosts | `[amr] seam_restriction`: `mean` (default, the 2:1 mean), `lagrange6` (C1), `weno6` (C3) | `update_ghost_local`, `restrict_seam_high_order`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `seam_restriction_code`, `SEAM_W_LAGRANGE6`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); key read in `initialize`, [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90) |
| Coarse-to-fine ghosts | `[amr] seam_ghost_fill`: `tricubic` (default), `quintic` (D1), `weno6`, `quintic-mp` (D2) | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `seam_interpolate_tricubic`, `seam_interpolate_quintic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| Reflux | Berger–Colella, `[numerics] reflux = .true.`; `reflux_policy = always` (default) or `troubled` (`euler-tq`) | `accumulate_seam_fluxes`, `apply_reflux_to_stage_forest`, `flag_troubled_seam_cells`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `reduce_fine_sums`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90) |
| Boundaries | periodic in $x$ and $y$; `extrapolation` on the null $z$ faces | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device; the restriction rows in the local ghost kernel (C2) | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); `update_ghost_local_gpu_dev`, [`adam_fnl_field_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_field_kernels.F90) |

The seam operators are described in [numerics](../../numerics#adaptive-mesh-refinement).

**C1, the 6-point restriction** (`[amr] seam_restriction = lagrange6`, #68 C1, C2). The coarse seam ghosts take the
tensor 6-point Lagrange interpolation of the fine cells to the coarse centre instead of their mean: exact 1-D weights
$k/256$, the footprint shifted inside the donor block at its edge, as [B3](./design-record#footprint-of-the-6-point-restriction)
recommends. A row whose fine pair does not lie inside its donor block falls back to the mean. On FNL the restriction
runs in the device twin of C2.

**D1, the quintic fill** (`[amr] seam_ghost_fill = quintic`, #68 D1). A fine ghost takes the tensor 6-point Lagrange
interpolation of the coarse cells at its centre: 1-D weights $k/8192$ (exact in binary), the footprint shifted inside
the coarse donor block at its edge as the tricubic one is (blocks of at least 6 cells along the refined axes); along an
unrefined axis only the anchor is read. The anchor shifts travel in the fill meta of every map row (bits 18–26,
`seam_meta_pack`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90)).
FNL sums straight from the device field; the inter-realm rows use the same fill.

**D2, the limited fills** (#68 D2), both order 6 at the centred anchor positions:
- `weno6`: per direction the WENO-Z combination of the cubics of the 4-point sub-stencils that contain the ghost's
  interval (the kernel of the C3 restriction), with positive linear weights that reproduce the Lagrange interpolant of
  their union; a sub-stencil holding a jump is dropped, and the clamped linear interpolation is used when none is left.
  Order 5 one cell from the donor block's edge, 4 at it (`seam_interpolate_weno6`);
- `quintic-mp`: the D1 quintic clamped per direction to a Suresh–Huynh monotonicity-preserving bound, the interval's
  values and the curvature-corrected linear value $L + \theta(\theta-1)/2\, d^{M4}$, widened by the curvature term; at
  the block edge the extrapolation is allowed only when the edge-side intervals hold no jump
  (`seam_interpolate_quintic_mp`).

The other fills compared on the shock crossings are `injection` (order 1) and `restriction-compatible` (order 2, exact
$R \circ P = I$ with the `mean` restriction). Every fill but `tricubic` is opt-in.

## Initial and boundary conditions

| Item | `euler` | `mhd` |
|---|---|---|
| Domain | $[0, 1]^2$, null $z$ | $[-7, 7]^2$, null $z$ |
| Initial condition | [V2 isentropic vortex](../../initial-conditions#isentropic-vortex): radius 0.07, strength 5, centre $(0.5, 0.5)$; free stream $(\rho, u, v, p) = (1, 1, 1, 1)$ | [MV-7 magnetised vortex](../../initial-conditions#mhd-vortex-magnetised-vortex): radius 1, $\kappa = 1$, $\mu = 1$, centre $(0, 0)$; background $(\rho, u, v, p) = (1, 1, 1, 1)$, $\mathbf{B} = 0$ |
| Base grid | quadtree (`ratio = 4`), $4 \times 4$ base blocks of $N/4 \times N/4 \times 1$ cells, `ngc = 3`; $N$ = 64 / 128 / 256 | the same |
| Refined region | quadrant $[0.5, 1]^2$: four base blocks refined one 2:1 level at initialisation (`primitive-box` marker, `max_level = 3`, no regridding) | quadrant $[0, 7]^2$, the same |
| Uniform control | the same case without the box (`max_level = 2`), the coarse spacing | the same |
| Boundaries | periodic in $x$, $y$; `extrapolation` in $z$ | the same |
| Final time | $t = 0.1$ (CFL 0.4; V2's 0.1 gives the same errors to five digits: uniform $N = 64$ L1 $8.2184\cdot10^{-5}$ both) | $t = 0.5$ |
| Ranks | 2 | 2 |
| `euler-octree` (on request) | the euler leg on an octree (`ratio = 8`, `nk = 4`), $N$ = 64 / 128; it refines along the null $z$ too, about 16 times the cost of the quadtree | — |

Legs of the euler ladder (each: refined and uniform runs, $N$ = 64 / 128 / 256):

| Leg | `seam_restriction` | `seam_ghost_fill` | Reflux | Phase |
|---|---|---|---|---|
| `euler` | `mean` | `tricubic` | on, `always` | A1 |
| `euler-l6` | `lagrange6` | `tricubic` | off | C1 |
| `euler-w6` | `weno6` | `tricubic` | off | C3, [C3c](./shock-crossing#the-weno-restriction) |
| `euler-q6` | `lagrange6` | `quintic` | off | [D1](#the-quintic-fill) |
| `euler-ww` | `weno6` | `weno6` | off | [D2](#the-limited-fills) |
| `euler-mp` | `lagrange6` | `quintic-mp` | off | [D2](#the-limited-fills) |
| `euler-tq` | `lagrange6` | `quintic` | on, `troubled` | [E2](./reflux-policy) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/make_amr_order.py` | writes the inputs from the V2 input `vortex/vortex-n{064,128,256}.ini` (euler) or `mhd/vortex/make_mhd_vortex.py` (mhd); `--ratio`, `--uniform`, `--no-reflux`, `--seam-restriction`, `--seam-ghost-fill`, `--reflux-policy`, `--time-max`, `--it-save` |
| `src/tests/flume/verification/mhd/amr_box.py` | `refine_box`: the refined quadrant, one 2:1 level at initialisation |
| `src/tests/flume/verification/amr-order/amr_order_oracle.py` | the per-region errors against the exact vortex, the orders, the bounds |
| `src/tests/flume/verification/make_doc_figures.py --only order-across-seams` | the figure, from the kept work directories of `amr-order/check.sh --keep` (CPU) |
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | the D1 and D2 shock crossings of [AS](./shock-crossing) (`--seam-restriction`, `--seam-ghost-fill`) |
| `src/tests/flume/verification/ghosts/check.sh --seam-restriction lagrange6 --seam-ghost-fill quintic --degree 5` | the [ghost probe](../infrastructure/ghost-probe) on a degree-5 field (D1); `--seam-ghost-fill F` on a cubic field (D2) |
| `src/tests/flume/verification/make_doc_figures.py --only quintic-fill` | the D1 figure, from `ghosts/check.sh --degree 5 --keep` and the kept AO work directories (CPU) |
| `src/tests/flume/verification/make_doc_runs.sh <runs> seams`, `make_doc_figures.py --only limited-fills` | the D2 figure: re-runs the crossings `d-{quintic,weno6,quintic-mp}-{c2f,f2c}` (CPU) |
| `src/tests/flume/verification/amr-order/check.sh` | runs the legs (`--legs euler,mhd,euler-octree,euler-l6,euler-w6,euler-q6,euler-ww,euler-mp,euler-tq`, default all but `euler-octree`; `--np N`, `--keep`; `FLUME_EXE` selects the backend) |

```bash
cd src/tests/flume/verification/amr-order && ./check.sh             # CPU, default legs (all but euler-octree)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh         # FNL
./check.sh --legs euler-octree                                      # the octree leg, on request
```

**Oracle** (`amr_order_oracle.py`). Each cell of the last checkpoint is compared with the exact convected vortex: the
density for Euler; for MHD the eight conservative variables, combined as MV-7 does ($\sqrt{\sum_v \bar e_v^2}$ of the
per-variable means). The error is split into the **seam band**, the cells within `ngc` of their own spacing from the
boundary of the refined quadrant (periodic images included), and the **interior**, the rest; a uniform run is all
interior. For each region: L1 (the volume-weighted mean of $|e|$), L∞ and the orders of successive runs.

| Bound | Leg | Value |
|---|---|---|
| finest refined L1 (all, band) and L∞ (all) | every leg but `euler-octree` | at most the CPU baseline plus 2 % |
| uniform control, interior L1 order | every leg | $\ge 4.7$ |
| refined L1 / L∞ order, finest pair | `euler-l6`, `euler-w6`, `euler-ww` | $\ge 4.8$ / $\ge 4.5$ |
| refined L1 / L∞ order, finest pair | `euler-q6`, `euler-mp`, `euler-tq` | $\ge 5.5$ / $\ge 5.5$ |

Cost of `euler-l6`: about 10 min on the CPU, 2 on FNL.

## Results

![AO ladders and error maps](/flume/order-across-seams.png)

*CPU, from the work directories `amr-order/check.sh --legs euler,euler-l6,euler-q6 --keep` leaves: $L_1(\rho)$ over all cells and $L_\infty(\rho)$ in the seam band against $N$ (the 128 → 256 orders in the legend), the default operators (mean + tricubic, reflux on) capped at 2.01 and 1.46, lagrange6 without reflux at the uniform control, lagrange6 + quintic below it; and the density error maps at $N = 128$ (the refined quadrant in white): the default operators' error concentrates on the seams, lagrange6 + quintic leaves the vortex's own.*

Baselines (A1). CPU and FNL identical to every printed digit, 2 ranks; orders of the finest pair, $N$ = 128 → 256:

| Leg | Refined L1 at N = 256 | Uniform L1 at N = 256 | L1 order: all / band / interior | L∞ order | Uniform L1 order |
|---|---|---|---|---|---|
| euler | $4.05\cdot10^{-6}$ | $8.38\cdot10^{-8}$ | 2.01 / 1.71 / 1.91 | 1.46 | 5.79 |
| mhd | $1.83\cdot10^{-6}$ | $2.62\cdot10^{-8}$ | 2.07 / 1.85 / 1.92 | 1.89 | 5.34 |

C1, the 6-point restriction (`euler-l6`). CPU and FNL (the device twin of C2) agree to every printed digit, and on the
same run to $3\cdot10^{-14}$:

| N = 256 | L1 all | L∞ all | Orders 128 → 256, L1 / L∞ |
|---|---|---|---|
| euler (mean, reflux on, default) | $4.05\cdot10^{-6}$ | – | 2.01 / 1.46 |
| euler-l6 (lagrange6, reflux off) | $8.55\cdot10^{-8}$ | $6.36\cdot10^{-6}$ | **5.33 / 5.53** |
| uniform at the coarse spacing | $8.38\cdot10^{-8}$ | $5.34\cdot10^{-6}$ | 5.79 / 7.38 |

`euler-w6` reports in [C3c](./shock-crossing#the-weno-restriction) and `euler-tq` in [E2](./reflux-policy).

### The quintic fill

Phase D1: the `quintic` coarse-to-fine fill.

![D1: the quintic fill in situ and on the AO vortex](/flume/quintic-fill.png)

*CPU. Left and centre: the ghost probe on a degree-5 field, box3d (octree) and channel2d (quadtree), 2 ranks, the largest error per ghost group: mean + tricubic misses the field by up to $5\cdot10^{-5}$ in every group a seam operator fills (the expected failure: exact to degree 1 and 3), lagrange6 + quintic matches it to $10^{-14}$. Right: the AO ladder without reflux, lagrange6 with the tricubic and the quintic fill.*

FNL, 2 ranks, `lagrange6` restriction in every run. The shock crossings are the [AS](./shock-crossing) `c2f` and
`f2c` runs (Sod $N = 200$, RJ2a $N = 256$, reflux on, the null axes widened to 8 cells); plateau values in units of
$10^{-2}$, band values absolute:

| Check | lagrange6 + tricubic | lagrange6 + quintic |
|---|---|---|
| AO euler, reflux off, N = 256: L1 all / band, L∞ | $8.55\cdot10^{-8}$ / $2.61\cdot10^{-7}$, $6.36\cdot10^{-6}$ | $6.03\cdot10^{-8}$ / $1.74\cdot10^{-7}$, $4.78\cdot10^{-6}$ |
| AO orders L1 / L∞ (128 → 256) | 5.33 / 5.53 | 6.11 / 6.35 |
| Sod c2f plateau $\rho$ / $p$, band $\rho$ | $2.457$ / $0.691$, $0.633\cdot10^{-2}$ | $2.631$ / $0.858$, $0.795\cdot10^{-2}$ |
| Sod f2c plateau $\rho$ / $p$, band $\rho$ | $1.575$ / $1.015$, $0.343\cdot10^{-2}$ | $1.573$ / $0.985$, $0.338\cdot10^{-2}$ |
| RJ2a c2f plateau $\rho$ / $p$, band $\rho$ | $1.410$ / $0.622$, $1.03\cdot10^{-4}$ | $2.297$ / $0.986$, $0.68\cdot10^{-4}$ |
| PV-5 split / HLLD: most faces limited, inadmissible stages | 46336 / 17920, 0 | 45644 / 17920, 0 |

Backend agreement: np 1 / 2 / 3 bitwise on the CPU, np 1 / 2 on FNL, CPU against FNL $2.7\cdot10^{-13}$; the
`euler-q6` ladder agrees on CPU and FNL to every printed digit. [Multirealm](../euler/multi-realm) leg 6 runs the fill
across a refined inter-realm seam, bitwise.

### The limited fills

Phase D2: the limited fills `weno6` and `quintic-mp`.

![D2: the plateau with the 6-point fills](/flume/limited-fills.png)

*Sod at $N = 200$, lagrange6, reflux on, CPU: the plateau behind the shock with the quintic, weno6 and quintic-mp fills (and the tricubic, c2f), coarse to fine and fine to coarse. The fills differ by at most about $10^{-3}$ in density, against a defect of $2.5\cdot10^{-2}$ of the jump: the defect is not the fill.*

Every fill on the `c2f` crossings (FNL, 2 ranks, `lagrange6` restriction, reflux on; plateau and band defects as in
[AS](./shock-crossing)):

| Fill (order) | Sod plateau $\rho$ | Sod band $\rho$ | Sod shock (cells) | RJ2a plateau $\rho$ | RJ2a band $\rho$ |
|---|---|---|---|---|---|
| `injection` (1) | $2.98\cdot10^{-2}$ | $7.91\cdot10^{-3}$ | 0.173 | $2.49\cdot10^{-2}$ | $2.16\cdot10^{-4}$ |
| `restriction-compatible` (2) | $2.32\cdot10^{-2}$ | $5.31\cdot10^{-3}$ | 0.164 | $1.98\cdot10^{-2}$ | $1.28\cdot10^{-4}$ |
| `tricubic` (4) | $2.46\cdot10^{-2}$ | $6.33\cdot10^{-3}$ | 0.158 | $1.41\cdot10^{-2}$ | $1.03\cdot10^{-4}$ |
| `quintic` (6) | $2.63\cdot10^{-2}$ | $7.95\cdot10^{-3}$ | 0.159 | $2.30\cdot10^{-2}$ | $0.68\cdot10^{-4}$ |
| `weno6` (6, limited) | $3.01\cdot10^{-2}$ | $7.75\cdot10^{-3}$ | 0.178 | $1.48\cdot10^{-2}$ | $1.52\cdot10^{-4}$ |
| `quintic-mp` (6, limited) | $3.14\cdot10^{-2}$ | $8.82\cdot10^{-3}$ | 0.173 | $2.47\cdot10^{-2}$ | $6.25\cdot10^{-4}$ |

Smooth flow (AO, $N$ = 128 → 256):

| Leg | L1 order | L∞ order | L1 at N = 256 | L∞ at N = 256 |
|---|---|---|---|---|
| `euler-ww` (`weno6` + `weno6`) | 5.32 | 5.50 | $8.60\cdot10^{-8}$ | $6.61\cdot10^{-6}$ |
| `euler-mp` (`lagrange6` + `quintic-mp`) | 6.10 | 6.33 | $6.04\cdot10^{-8}$ | $4.77\cdot10^{-6}$ |
| `euler-q6` (`lagrange6` + `quintic`) | 6.11 | 6.35 | $6.03\cdot10^{-8}$ | $4.78\cdot10^{-6}$ |

PV-5 passes with no inadmissible stage under both limited fills; the ghost probe holds a cubic field to round-off with
either fill; multirealm leg 6 runs both across a refined inter-realm seam, bitwise.

Related cases: [AS](./shock-crossing) is the shock baseline of the same seam operators; [AL](./integral-drift) logs
the drift of this vortex without the reflux; [B1](./design-record#the-caps-in-a-2-d-prototype) reproduces the caps
in a 2-D prototype; [E0](./design-record#the-shock-crossing-defect-is-a-capturing-error) explains the crossing defect
D2 could not reduce; the uniform-grid orders are [V2](../euler/isentropic-vortex) and
[MV-5 to MV-7](../mhd/order-of-accuracy).

## Findings

- **The seams cap the whole solution at second order, not only the band.** The error made at a seam is carried into
  the interior by the flow, so every region a feature reaches after crossing a seam is second order. The plan's
  expectation of "design order in the interior" holds only for regions the flow does not reach from a seam.
- **On this smooth flow the refinement makes the solution worse.** At $N = 256$ the refined composite error is 48
  times (Euler) and 70 times (MHD) the error of the uniform run at the coarse spacing, which reaches design order.
- **Euler's L∞ order is 1.46,** below the 2 of the 1-D prototype in #68 (Appendix A). The vortex crosses the seam
  corner, where the two seam directions meet; the 2-D prototypes of Phase B show that the corner does not explain it
  ([B1, B2](./design-record#the-caps-in-a-2-d-prototype)).
- **C1: the restriction was the seam error.** With the 6-point restriction and without the reflux the refined
  composite error equals the uniform coarse run's, 47 times below the default operators': at these resolutions the restriction was the
  seam error, and the tricubic prolongation's $O(h^4)$ term has not shown yet.
- **C1: the reflux caps the order by itself.** With the reflux on, the restriction changes nothing (band L∞ order 1.43
  against 1.46 at $N$ = 64 → 128), as [B1](./design-record#the-caps-in-a-2-d-prototype) predicts. The high order
  with the reflux kept at the shocks is [E2](./reflux-policy).
- No restriction row fell back to the mean: every fine pair lies inside its donor block.
- **C1 cost on FNL:** the restriction rows read 36 cells instead of 4 on this quadtree, in the same kernels as the
  copies.
- **D1: the smooth order reaches the scheme's.** With the reflux off, the 6-point restriction and the quintic fill,
  AO euler is sixth order (6.11 / 6.35) and its error at $N = 256$, $6.03\cdot10^{-8}$, is below the uniform coarse
  run's ($8.38\cdot10^{-8}$): the 2:1 seams no longer cost accuracy on this flow.
- **D1: unlimited, the quintic fill worsens a shock entering a finer region.** On the `c2f` crossings the plateau
  defect grows (Sod +7 % in $\rho$, +24 % in $p$; RJ2a +63 % / +59 %): the shock's coarse values reach the fine ghosts
  through a 6-point stencil with larger overshoots than the 4-point one. `f2c` is unchanged, the shock positions are
  the same, mass is exact, and PV-5 passes with no inadmissible stage.
- **D2: the crossing defect has no trend with the fill's order, and limiting makes it worse.** The best fill differs by
  case and by metric (`restriction-compatible` on Sod, `tricubic` on the RJ2a plateau, `quintic` on its band); every
  unlimited fill sits within about $\pm 30\,\%$ of `tricubic`, and the two limited fills are the worst on Sod, their
  shock furthest off. The plateau defect is not an overshoot of the fill, so no fill can be held to "no worse than
  `tricubic`"; [E0](./design-record#the-shock-crossing-defect-is-a-capturing-error) shows it is a capturing error.
- **D2: on smooth flows both limited fills keep their order.** `euler-ww` has the `weno6` restriction's accuracy;
  `euler-mp` equals `euler-q6`: the MP bound never clips the vortex.
- `quintic`, `weno6` and `quintic-mp` stay opt-in; the default fill is `tricubic`.
