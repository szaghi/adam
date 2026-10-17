# Order across static 2:1 seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| AO (A1, C1) | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL | yes | — |

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

Phases C (high-order restriction), D (quintic prolongation) and E (reflux policy) of #68 are measured against these
baselines; the legs `euler-l6`, `euler-w6`, `euler-q6`, `euler-ww`, `euler-mp` and `euler-tq` are their entries in the
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
$k/256$, the footprint shifted inside the donor block at its edge, as [B3](../seams-studies/high-order-restriction)
recommends. A row whose fine pair does not lie inside its donor block falls back to the mean. On FNL the restriction
runs in the device twin of C2.

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
| `euler-w6` | `weno6` | `tricubic` | off | C3, [C3c](./shock-crossing#c3c-the-weno-restriction) |
| `euler-q6` | `lagrange6` | `quintic` | off | [D1](../seams-studies/quintic-fill) |
| `euler-ww` | `weno6` | `weno6` | off | [D2](../seams-studies/limited-fills) |
| `euler-mp` | `lagrange6` | `quintic-mp` | off | [D2](../seams-studies/limited-fills) |
| `euler-tq` | `lagrange6` | `quintic` | on, `troubled` | [E2](../seams-studies/reflux-policy) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/make_amr_order.py` | writes the inputs from the V2 input `vortex/vortex-n{064,128,256}.ini` (euler) or `mhd/vortex/make_mhd_vortex.py` (mhd); `--ratio`, `--uniform`, `--no-reflux`, `--seam-restriction`, `--seam-ghost-fill`, `--reflux-policy`, `--time-max`, `--it-save` |
| `src/tests/flume/verification/mhd/amr_box.py` | `refine_box`: the refined quadrant, one 2:1 level at initialisation |
| `src/tests/flume/verification/amr-order/amr_order_oracle.py` | the per-region errors against the exact vortex, the orders, the bounds |
| `src/tests/flume/verification/make_doc_figures.py --only order-across-seams` | the figure, from the kept work directories of `amr-order/check.sh --keep` (CPU) |
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

*CPU, from the work directories `amr-order/check.sh --legs euler,euler-l6,euler-q6 --keep` leaves: $L_1(\rho)$ over all cells and $L_\infty(\rho)$ in the seam band against $N$ (the 128 → 256 orders in the legend), today's operators (mean + tricubic, reflux on) capped at 2.01 and 1.46, lagrange6 without reflux at the uniform control, lagrange6 + quintic below it; and the density error maps at $N = 128$ (the refined quadrant in white): today's error concentrates on the seams, lagrange6 + quintic leaves the vortex's own.*

Baselines (A1). CPU and FNL identical to every printed digit, 2 ranks; orders of the finest pair, $N$ = 128 → 256:

| Leg | Refined L1 at N = 256 | Uniform L1 at N = 256 | L1 order: all / band / interior | L∞ order | Uniform L1 order |
|---|---|---|---|---|---|
| euler | $4.05\cdot10^{-6}$ | $8.38\cdot10^{-8}$ | 2.01 / 1.71 / 1.91 | 1.46 | 5.79 |
| mhd | $1.83\cdot10^{-6}$ | $2.62\cdot10^{-8}$ | 2.07 / 1.85 / 1.92 | 1.89 | 5.34 |

C1, the 6-point restriction (`euler-l6`). CPU and FNL (the device twin of C2) agree to every printed digit, and on the
same run to $3\cdot10^{-14}$:

| N = 256 | L1 all | L∞ all | Orders 128 → 256, L1 / L∞ |
|---|---|---|---|
| euler (mean, reflux on, today) | $4.05\cdot10^{-6}$ | – | 2.01 / 1.46 |
| euler-l6 (lagrange6, reflux off) | $8.55\cdot10^{-8}$ | $6.36\cdot10^{-6}$ | **5.33 / 5.53** |
| uniform at the coarse spacing | $8.38\cdot10^{-8}$ | $5.34\cdot10^{-6}$ | 5.79 / 7.38 |

The other legs of the ladder report on their study pages: `euler-w6` in [C3c](./shock-crossing#c3c-the-weno-restriction),
`euler-q6` in [D1](../seams-studies/quintic-fill), `euler-ww` and `euler-mp` in [D2](../seams-studies/limited-fills),
`euler-tq` in [E2](../seams-studies/reflux-policy).

Related cases: [AS](./shock-crossing) is the shock baseline of the same seam operators; [AL](./integral-drift) logs
the drift of this vortex without the reflux; [B1](../seams-studies/caps-prototype) reproduces the caps in a 2-D
prototype; the uniform-grid orders are [V2](../euler/isentropic-vortex) and
[MV-5 to MV-7](../mhd/order-of-accuracy).

## Findings

- **The seams cap the whole solution at second order, not only the band.** The error made at a seam is carried into
  the interior by the flow, so every region a feature reaches after crossing a seam is second order. The plan's
  expectation of "design order in the interior" holds only for regions the flow does not reach from a seam.
- **On this smooth flow the refinement makes the solution worse.** At $N = 256$ the refined composite error is 48
  times (Euler) and 70 times (MHD) the error of the uniform run at the coarse spacing, which reaches design order.
- **Euler's L∞ order is 1.46,** below the 2 of the 1-D prototype in #68 (Appendix A). The vortex crosses the seam
  corner, where the two seam directions meet; the 2-D prototype of Phase B was to check whether the corner explains
  it ([B1](../seams-studies/caps-prototype) and [B2](../seams-studies/sensor-gated-reflux) found it does not).
- **C1: the restriction was the seam error.** With the 6-point restriction and without the reflux the refined
  composite error equals the uniform coarse run's, 47 times below today's: at these resolutions the restriction was the
  seam error, and the tricubic prolongation's $O(h^4)$ term has not shown yet.
- **C1: the reflux caps the order by itself.** With the reflux on, the restriction changes nothing (band L∞ order 1.43
  against 1.46 at $N$ = 64 → 128), as [B1](../seams-studies/caps-prototype) predicts. Reaching the high order with the
  reflux kept at the shocks is Phase E's job ([E2](../seams-studies/reflux-policy)).
- No restriction row fell back to the mean: every fine pair lies inside its donor block.
- **C1 cost on FNL:** about 16 % of the wall time of this run ($N = 256$, reflux off: 56 s against 48 s with the mean;
  WSL, indicative only): the restriction rows read 36 cells instead of 4 on this quadtree, in the same kernels as the
  copies.
