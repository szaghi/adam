# Seam design record

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| B1, B2, B3, E0, G2 | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) · B1–B3, E0, G2 | P | NumPy prototypes (B); FNL, CPU re-runs (E0, G2) | yes | — |

This appendix keeps the studies of M6 that chose a design or ruled one out without adding a key of their own. The
verification cases that freeze the outcome are [AO](./order-across-seams), [AD](./divergence-across-seams),
[AS](./shock-crossing), [AL](./integral-drift), [reflux policy](./reflux-policy) and [EGLM reflux](./eglm-reflux).

| Study | Question | Conclusion |
|---|---|---|
| [B1](#b1-the-caps-in-a-2-d-prototype) | do the 1-D order caps of a 2:1 seam hold in 2-D, and do the corners explain AO's L∞ order? | yes; no |
| [B2](#b2-a-sensor-gated-reflux-prototype) | does a sensor-gated reflux recover the order and keep the shocks? | yes, near-conservatively: built as `reflux_policy = troubled` |
| [B3](#b3-footprint-of-the-6-point-restriction) | shift-inward or two-pass footprint for the 6-point restriction? | shift-inward |
| [E0](#e0-the-shock-crossing-defect-is-a-capturing-error) | does a seam operator cause the defect behind a crossing shock? | no: a capturing error |
| [G2](#g2-the-glm-damping-length-on-a-seam) | would a per-level GLM damping length matter on a seam? | no: not added |

The prototypes run from `src/tests/flume/verification/amr-proto/check.sh` (`--legs b1,b2,b3`, a leg runs only when the
ones before it pass; `--long` for the finer ladders). They are NumPy models in `proto2d.py`, not FLUME runs: a 2:1 box
$[0.25, 0.75]^2$ on the periodic unit square, the coarse cells under the box holding the restriction of the fine ones,
and the FLUME seam operators (`mean` or 6-point `lag6` restriction, `cubic` or `quintic` prolongation, reflux on or
off) reimplemented on it.

## B1: the caps in a 2-D prototype

**Tried.** $u_t + u_x + u_y = 0$, $\sin 2\pi x \sin 2\pi y$ to $t = 1$, conservative finite differences with the
fifth-order linear upwind flux, RK4 with $\Delta t \propto h^{5/4}$. Every arm of restriction × prolongation × reflux
against the uniform coarse control, $N$ = 128 → 256 (`--long`); the L∞ error split into corner, edge and interior
regions (within 3 coarse cells of a box corner, of a box edge, the rest). `check.sh` asserts lower bounds only: uniform
L∞ order $\ge 4.8$, `lag6`-`quintic`-off $\ge 4.7$, `lag6`-`cubic`-off $\ge 3.7$.

![B1: order per arm and where the error lives](/flume/seam-proto-b1.png)

*The `check.sh` ladder, $N$ = 32 / 64 / 128: $L_1$ and $L_\infty$ per arm against the order-2 and order-5 slopes, and
the $L_\infty$ error at $N = 128$ by region. Every arm with the reflux on, or with the mean restriction, is capped at
order 2; only lag6 + quintic without reflux reaches 5, its error spread evenly over corners, edges and interior.*

| Arm | L1 order | L∞ order | Mass drift at N = 256 |
|---|---|---|---|
| uniform (coarse) | 5.01 | 5.01 | $5\cdot10^{-18}$ |
| mean, cubic, reflux on (the default) | 2.85 | 1.93 | $1\cdot10^{-18}$ |
| mean, cubic, off | 2.01 | 1.95 | $5\cdot10^{-8}$ |
| lag6, cubic, on | 2.85 | 1.93 | $2\cdot10^{-18}$ |
| lag6, cubic, off | 3.98 | 4.01 | $1.3\cdot10^{-9}$ |
| lag6, quintic, on | 2.85 | 1.93 | $6\cdot10^{-19}$ |
| lag6, quintic, off | **4.99** | **4.99** | $9.1\cdot10^{-11}$ |

**Conclusion.** The mean restriction caps the order at 2 whatever the prolongation; the reflux caps it at 2 whatever
the transfers. The 6-point restriction with the quintic prolongation and no reflux recovers order 5 with the corners
present. The corners hold the maximum error in every arm but converge at the same rate as edges and interior, so they
do not explain the L∞ order 1.46 of [AO](./order-across-seams)'s Euler vortex.

## B2: a sensor-gated reflux prototype

**Tried.** The B1 grid with WENO5-JS and Lax–Friedrichs splitting (the weights of `adam_weno_object`). Three reflux
policies on `lag6`/`quintic`: `always`, `never`, and `troubled`, a coarse box face refluxed only where the sensor
$r = \min_k w_k/d_k$ of the coarse face or of either fine face is below $\tau = 0.2$. Legs: the smooth B1 advection,
and Burgers pulses crossing the box coarse to fine (`c2f`), fine to coarse (`f2c`), and a weak $10^{-3}$ pulse that
the JS sensor misses (`weak`, the missed-detection control). `--check` asserts the smooth orders of `troubled`
$\ge 4.7$, `troubled` within 0.05 cells of `always` on the shocks with a mass defect $\le 10^{-3}$, and `never`'s mass
defect $\ge 10^{-3}$ (the negative control).

![B2: the smooth order and the shock crossings per reflux policy](/flume/seam-proto-b2.png)

*Left: the smooth leg, $L_1$ against $N$ per arm; `always` is capped near order 3, `troubled` keeps the uniform order.
Centre and right: the Burgers pulses at $N = 64$, the shock position error and the mass defect per policy: `troubled`
places the strong shocks where `always` does with a mass defect below $10^{-3}$ of the pulse; `never` misplaces the c2f
shock by a cell and loses $1.8\cdot10^{-2}$; on the weak pulse `troubled` equals `never`.*

| Leg, N | always | troubled | never | Faces flagged |
|---|---|---|---|---|
| smooth, L1 / L∞ order 64 → 128 | 2.96 / 2.10 | **5.01 / 4.95** | 5.01 / 4.95 | 0 % |
| c2f, 128: shock (local cells), mass defect | −0.50, $2\cdot10^{-16}$ | −0.50, $3.6\cdot10^{-5}$ | −1.07, $8.6\cdot10^{-3}$ | 8 % |
| f2c, 128: shock, mass defect | −0.25, $4\cdot10^{-16}$ | −0.25, $6.9\cdot10^{-6}$ | +0.06, $8.4\cdot10^{-3}$ | 7 % |
| weak, 128: shock, mass defect | −0.015, $10^{-12}$ | −0.022, $3.2\cdot10^{-5}$ | −0.022, $3.2\cdot10^{-5}$ | 0 % |

**Conclusion.** On smooth data the sensor flags no seam face, so `troubled` recovers order 5. On shocks it places them
within 0.003 cells of `always`; the flagged faces (7 to 15 % of the seam-face evaluations) are what holds the shock.
It is near-conservative, not conservative: a mass defect of $7\cdot10^{-6}$ to $6.5\cdot10^{-5}$ of the pulse from the
unflagged smooth faces. A missed weak shock costs 0.007 to 0.011 cells and $3\cdot10^{-5}$ to $4\cdot10^{-4}$ of its
mass. The order is recovered only together with the high-order transfers: with the mean restriction it stays at 2
(B1). Built in FLUME as [`reflux_policy = troubled`](./reflux-policy).

## B3: footprint of the 6-point restriction

**Tried.** The 6-point restriction reaches 2 to 3 fine cells beyond the coarse cell, so near the edge of a fine block it
needs data the block does not hold. Two footprints, on the B1 problem without reflux, with the fine box split into
blocks of $B$ = 4, 6, 8, 16 cells or left whole, $N$ = 48 / 96 / 192 (`--long`):
- **shift-inward:** the stencil shifts inside the donor block (one-sided, needs $B \ge 6$);
- **two-pass:** the fine ghosts are filled first and the stencil stays centred (any $B$, a second synchronisation).

`check.sh` asserts the L∞ orders of shift-inward ($B = 6$) and two-pass $\ge 3.7$ with the cubic prolongation and
$\ge 4.7$ with the quintic.

![B3: the restriction footprint against the fine block size](/flume/seam-proto-b3.png)

*The `check.sh` ladder, $N$ = 24 / 48 / 96: $L_\infty$ per footprint, cubic (left) and quintic (right) prolongation.
Shift-inward at 6 or 8 fine cells per block matches the whole box; two-pass is slightly more accurate with the cubic
fill and identical with the quintic one.*

| Footprint, prolongation | Block size | L1 at N = 192 | Order L1 / L∞ |
|---|---|---|---|
| shift-inward, cubic | 6 / 8 / 16 / box | $9.2$ to $9.4\cdot10^{-8}$ | 3.98 / 4.00 |
| two-pass, cubic | any | $2.4\cdot10^{-8}$ | 4.14 / 4.03 |
| shift-inward, quintic | 6 / 8 / 16 / box | $2.38\cdot10^{-9}$ | 4.98 / 4.99 |
| two-pass, quintic | any | $2.42\cdot10^{-9}$ | 5.00 / 5.00 |
| shift-inward | 4 | infeasible | |

**Conclusion.** Shift-inward is insensitive to the block size from 6 cells up (2.5 % in L1) and, with the quintic
prolongation, equivalent to two-pass. It keeps `update_ghost_local` a single pass and needs 6 fine cells per refined
axis, which the [#66](https://github.com/szaghi/adam/issues/66) guard ($n \ge 2\,n_{gc}$) already requires. Two-pass
buys a factor 4 in L1 only with the cubic prolongation, at the cost of a second local pass and a second exchange round.
C1 implements shift-inward (`seam_restriction_code`,
[`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90)).

## E0: the shock-crossing defect is a capturing error

**Tried.** [AS](./shock-crossing) finds an $O(1)$, non-converging defect behind a shock that crosses a 2:1 face, and
[D2](./order-across-seams#d2-the-limited-fills) shows that no fill reduces it. E0 ran the AS `c2f` crossings (Sod
$N = 200$ and 400, RJ2a $N = 256$; FNL, 2 ranks, CPU re-runs agree to the printed digits) for every combination of fill
(`restriction-compatible`, the one fill with $R \circ P = I$ under `mean`, and `tricubic`), restriction (`mean`,
`lagrange6`) and reflux (on, off). No `check.sh` leg runs E0; the inputs come from
`amr-shock/make_amr_shock.py` (`--seam-ghost-fill`, `--seam-restriction`, `--no-reflux`).

![E0: the plateau defect against fill, restriction and reflux](/flume/capturing-error.png)

*The table below as bars: no combination of fill, restriction and reflux removes the defect.*

![E0: the Sod plateau and the defect per arm, CPU](/flume/capturing-error-sod.png)

*CPU re-runs (`make_doc_runs.sh <runs> seams`): left, the plateau after the c2f crossing with the tricubic and the
compatible fills at $N$ = 200 and 400, the bump behind the face the same at both resolutions; right, the defect of every
arm.*

Plateau density defect behind the shock:

| Fill + restriction | Sod, reflux on | Sod, reflux off | RJ2a, reflux on | RJ2a, reflux off |
|---|---|---|---|---|
| compatible + mean ($R \circ P = I$) | $2.32\cdot10^{-2}$ | $1.82\cdot10^{-2}$ | $1.97\cdot10^{-2}$ | $1.83\cdot10^{-2}$ |
| compatible + lagrange6 | $2.32\cdot10^{-2}$ | $1.85\cdot10^{-2}$ | $1.98\cdot10^{-2}$ | $1.84\cdot10^{-2}$ |
| tricubic + mean | $2.46\cdot10^{-2}$ | $3.10\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |
| tricubic + lagrange6 | $2.46\cdot10^{-2}$ | $3.17\cdot10^{-2}$ | $1.41\cdot10^{-2}$ | $1.68\cdot10^{-2}$ |

Sod with the reflux on at $N = 400$: $2.35\cdot10^{-2}$ (tricubic), $2.14\cdot10^{-2}$ (compatible), while the band
defect next to the seam converges ($6.35\cdot10^{-3} \to 1.71\cdot10^{-3}$, order 1.9).

**Conclusion.** The restriction is irrelevant (within 0.5 %); $R \circ P = I$ is best on Sod and worse on RJ2a; the
reflux moves the defect by up to 30 % either way; doubling $N$ lowers it by 5 to 8 %. No seam operator controls it. The
interpretation (not measured directly) is the start-up error of a captured shock whose numerical structure, a few cells
wide, is squeezed into finer cells: the transition emits a spurious wave whose amplitude scales with the shock's
strength, not with the spacing. The remedy is to keep shocks off the seams by refining ahead of them
([runtime AMR](../runtime-amr/), [numerics](../../numerics#adaptive-mesh-refinement)).

## G2: the GLM damping length on a seam

**Tried.** The GLM damping rate is $c_h^2/c_p^2 = \alpha\, c_h / L$ with `[mhd] glm_damping_length` $L$ a fixed length
or `min-cell` (the finest spacing; under runtime AMR frozen at the spacing of `max_level`). G2 brackets a per-level
length between the coarse and the finest spacings on the source-free CPAW of [AD](./divergence-across-seams) (`cpaw`
leg, `glm_ch` 2.0, `glm_alpha` 0.18, refined and uniform $N$ = 64 / 128, 10 periods, FNL, 2 ranks), measuring the seam
floor of D_flux (`make_amr_divb.py cpaw --glm-damping-length`, `amr_divb_oracle.py growth`).

![G2: the seam D_flux per damping length](/flume/glm-damping.png)

*FNL: the seam maximum of D_flux on the source-free CPAW for the three damping lengths, refined $N = 64$ and 128: every
length bounded, `min-cell` doubling the floor at $N = 128$.*

| `glm_damping_length` | rate $c_h^2/c_p^2$ (N = 128) | seam D_flux floor N = 64 / 128 | × uniform (N = 128) | growth N = 128 |
|---|---|---|---|---|
| 1.0 (A2) | 0.36 | $5.39\cdot10^{-4}$ / $5.43\cdot10^{-5}$ | 119 | 0.72 |
| coarse spacing $L_x/N$ | 20.6 | $6.96\cdot10^{-4}$ / $5.32\cdot10^{-5}$ | 116 | 1.02 |
| `min-cell` (finest spacing) | 41.2 | $8.30\cdot10^{-4}$ / $1.07\cdot10^{-4}$ | 233 | 1.01 |

**Conclusion.** Every length keeps the seam divergence bounded (growth 0.72 to 1.02 over 10 periods). Strong damping
raises the seam floor: `min-cell` doubles it at $N = 128$. A per-level length would sit between the two spacings, so
the frozen `min-cell` length costs at most that factor 2 and keeps the damping rate fixed through regrids. No
`per-level` option exists.
