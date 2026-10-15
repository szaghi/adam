# Limited coarse-to-fine fills

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| D2 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL (crossings); CPU, FNL (AO legs) | no ([#85](https://github.com/szaghi/adam/issues/85)) | — |

## Case and scope

The unlimited quintic fill of [D1](./quintic-fill) reaches sixth order on smooth flows but worsens a shock entering a
finer region. D2 asks whether a limited 6-point fill keeps the order and brings the shock crossings of A3 back to
`tricubic`'s level (the Phase D exit). Neither did, so D2 then ran every fill on the `c2f` crossings to find out what
the shock-crossing defect is not.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| `seam_ghost_fill = weno6` | per direction the WENO-Z combination of the cubics of the 4-point sub-stencils that contain the ghost's interval (the C3 restriction's kernel, now shared); linear weights reproducing the Lagrange interpolant of their union (all positive); sub-stencils holding a jump dropped; the clamped linear interpolation of the interval when none is left. Order 6 at the centred anchor positions, 5 one cell from the donor block's edge, 4 at it | `seam_interpolate_weno6`, `seam_fill_weno6_1d`, `seam_weno6_core`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| `seam_ghost_fill = quintic-mp` | the D1 quintic clamped per direction to a Suresh–Huynh monotonicity-preserving bound: the interval's values and the curvature-corrected linear value $L + \theta(\theta-1)/2\, d^{M4}$, widened by the curvature term; at the block edge, where the target lies a quarter cell outside the footprint, the extrapolation is allowed only when the edge-side intervals hold no jump. Order 6 at every target | `seam_interpolate_quintic_mp`, `seam_fill_quintic_mp_1d`, `seam_minmod4`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| The other fills compared | `injection` (order 1), `restriction-compatible` (order 2, exact $R \circ P = I$ with the `mean` restriction), `tricubic` (4, default), `quintic` (6) | `seam_interpolate_compatible`, `seam_interpolate_tricubic`, `seam_interpolate_quintic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Key | `[amr] seam_ghost_fill`, read in `initialize`; both limited fills opt-in, the default stays `tricubic` | [`adam_adam_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_adam_object.F90) |
| Restriction | `lagrange6` on the crossings; `weno6` (`euler-ww`) or `lagrange6` (`euler-mp`) on AO | [C1, C3c](../seams/shock-crossing) |
| Scheme, time, reflux | as in [AO](../seams/order-across-seams) (reflux off) and [AS](../seams/shock-crossing) (reflux on) | [numerics](../../numerics#adaptive-mesh-refinement) |

## Initial and boundary conditions

| Run | Setup |
|---|---|
| Crossings | the [AS](../seams/shock-crossing) `c2f` crossings, Sod $N = 200$ and RJ2a $N = 256$, `lagrange6` restriction, reflux on, every fill; FNL, 2 ranks |
| AO `euler-ww` | the [AO](../seams/order-across-seams) Euler vortex, `weno6` restriction and fill, reflux off, $N$ = 64 / 128 / 256 |
| AO `euler-mp` | the AO Euler vortex, `lagrange6` + `quintic-mp`, reflux off, $N$ = 64 / 128 / 256 |
| PV-5 | the limited blast across 2:1 seams of the [positivity limiter](../mhd/positivity-limiter), under both limited fills |
| GP | the [ghost probe](../infrastructure/ghost-probe) on a cubic field, under both limited fills |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/check.sh --leg euler-ww`, `--leg euler-mp` | the AO legs: orders $\ge 4.8$ / 4.5 (`euler-ww`) and $\ge 5.5$ (`euler-mp`), the errors at most the baseline plus 2 % |
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | the crossings (the fill key is added to the input by hand: the script has no `--seam-ghost-fill` option) |
| `src/tests/flume/verification/ghosts/check.sh --seam-ghost-fill F` | GP |
| `src/tests/flume/verification/multirealm/check.sh` | leg 6 runs `weno6` and `quintic-mp` across a refined inter-realm seam |

```bash
cd src/tests/flume/verification/amr-order && ./check.sh --leg euler-ww --leg euler-mp
```

## Results

No figure yet ([#85](https://github.com/szaghi/adam/issues/85): plateau profiles for the limited fills).

Every fill on the `c2f` crossings (FNL, 2 ranks, `lagrange6` restriction; plateau and band defects as in
[AS](../seams/shock-crossing)):

| Fill (order) | Sod plateau $\rho$ | Sod band $\rho$ | Sod shock (cells) | RJ2a plateau $\rho$ | RJ2a band $\rho$ |
|---|---|---|---|---|---|
| `injection` (1) | $2.98\cdot10^{-2}$ | $7.91\cdot10^{-3}$ | 0.173 | $2.49\cdot10^{-2}$ | $2.16\cdot10^{-4}$ |
| `restriction-compatible` (2) | $2.32\cdot10^{-2}$ | $5.31\cdot10^{-3}$ | 0.164 | $1.98\cdot10^{-2}$ | $1.28\cdot10^{-4}$ |
| `tricubic` (4) | $2.46\cdot10^{-2}$ | $6.33\cdot10^{-3}$ | 0.158 | $1.41\cdot10^{-2}$ | $1.03\cdot10^{-4}$ |
| `quintic` (6) | $2.63\cdot10^{-2}$ | $7.95\cdot10^{-3}$ | 0.159 | $2.30\cdot10^{-2}$ | $0.68\cdot10^{-4}$ |
| `weno6` (6, limited) | $3.01\cdot10^{-2}$ | $7.75\cdot10^{-3}$ | 0.178 | $1.48\cdot10^{-2}$ | $1.52\cdot10^{-4}$ |
| `quintic-mp` (6, limited) | $3.14\cdot10^{-2}$ | $8.82\cdot10^{-3}$ | 0.173 | $2.47\cdot10^{-2}$ | $6.25\cdot10^{-4}$ |

Smooth flows (AO, $N$ = 128 → 256):

| Leg | L1 order | L∞ order | L1 at N = 256 | L∞ at N = 256 |
|---|---|---|---|---|
| `euler-ww` (`weno6` + `weno6`) | 5.32 | 5.50 | $8.60\cdot10^{-8}$ | $6.61\cdot10^{-6}$ |
| `euler-mp` (`lagrange6` + `quintic-mp`) | 6.10 | 6.33 | $6.04\cdot10^{-8}$ | $4.77\cdot10^{-6}$ |
| `euler-q6` (`lagrange6` + `quintic`, [D1](./quintic-fill)) | 6.11 | 6.35 | $6.03\cdot10^{-8}$ | $4.78\cdot10^{-6}$ |

PV-5 passes with no inadmissible stage under both limited fills; GP holds a cubic field to round-off with either fill.

Related cases: [D1](./quintic-fill) is the unlimited fill; [E0](./capturing-error) isolates the cause of the
crossing defect; the baselines are [AO](../seams/order-across-seams) and [AS](../seams/shock-crossing).

## Findings

- **The defect has no trend with the fill's order, and limiting makes it worse.** The best fill differs by case and by
  metric (`restriction-compatible` on Sod, `tricubic` on the RJ2a plateau, `quintic` on its band); every fill but the
  limited ones sits within roughly $\pm 30\,\%$ of `tricubic`, and the two limited fills are the worst on Sod, their
  shock furthest off. The plateau defect behind a crossing shock is not an overshoot of the fill, so no fill can be
  held to "no worse than `tricubic`".
- **The best Sod result belongs to the one fill with exact $R \circ P = I$:** the defect points to the consistency of
  the fill with the restriction at the crossing, and to the reflux (A3: reflux on/off moves it by up to 30 %), which is
  Phase E's ground, not the fill's order. [E0](./capturing-error) tested both and ruled both out.
- **On smooth flows both limited fills keep their order.** `euler-ww` has the `weno6` restriction's accuracy (the
  fill's 6-point WENO does not lose to the tricubic); `euler-mp` equals `euler-q6`: the MP bound never clips the
  vortex.
