# Footprint of a high-order restriction

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| B3 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | yes | — |

## Case and scope

Phase C replaces the mean restriction with a 6-point one. Its stencil reaches 2 to 3 fine cells beyond the two cells of
the coarse cell, so near the edge of a fine block it needs data the block does not hold. Two footprints are possible:
- **shift-inward:** the stencil shifts inside the donor block (one-sided, at least 6 fine cells per block);
- **two-pass:** the fine block's ghosts are filled first and the stencil stays centred (any block size, at the price of
  an ordering constraint).

B3 asks which one Phase C should implement: how accurate each is, how the block size matters, and what each costs in
FLUME's ghost update.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | the [B1](./caps-prototype) problem ($u_t + u_x + u_y = 0$, 5th-order linear upwind flux, RK4) with the reflux off, where the restriction's order shows | `Composite`, `run`, `proto2d.py` |
| Fine blocks | the fine box split into blocks of $B$ fine cells per axis | `Composite`, `proto2d.py` |
| `lag6` (shift-inward) | the 6-point restriction shifted inside its donor fine block (needs $B \ge 6$) | `Composite._restriction_stencils`, `proto2d.py` |
| `lag6c` (two-pass) | the centred 6-point restriction reading the fine ghosts (exact fine neighbours inside the box, the prolongation at its edge) | `Composite._restriction_stencils`, `proto2d.py` |
| Prolongation | `cubic` or `quintic` | `Composite.pad_fine`, `proto2d.py` |
| FLUME counterpart | shift-inward, adopted in C1: `seam_restriction_code` keeps the footprint inside the donor block | `seam_restriction_code`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); `update_ghost_local`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain, initial condition | the periodic unit square, $\sin 2\pi x \sin 2\pi y$, one period |
| Block sizes | $B$ = 4, 6, 8, 16 and the whole box (the `proto2d.py` default; `check.sh` runs 4, 6, 8 and the box) |
| Grids | $N$ = 48 / 96 / 192 (so that 6, 8 and 16 divide the fine cells; `check.sh --long`); $N$ = 24 / 48 / 96 (`check.sh`) |
| Reflux | off |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-proto/proto2d.py` | `b3 [--cells ...] [--blocks 4 6 8 16 0] [--assert ARM:NORM:MIN]` (block 0: the whole box) |
| `src/tests/flume/verification/amr-proto/check.sh` | runs B1, B2 and B3 (`--legs b3` for B3 alone) |
| `src/tests/flume/verification/make_doc_figures.py --only seam-proto-b3` | draws the figure (`fig_seam_proto_b3`; NumPy only, about 20 s) |

```bash
cd src/tests/flume/verification/amr-proto && ./check.sh --legs b3  # B3 at N = 24 / 48 / 96, about 20 s
./check.sh --legs b3 --long                                        # B3 at N = 48 / 96 / 192
```

Bounds asserted by `check.sh`: the L∞ orders of shift-inward ($B = 6$) and two-pass at least 3.7 with the cubic
prolongation and 4.7 with the quintic; shift-inward at $B = 4$ is reported infeasible.

## Results

![B3: the restriction footprint against the fine block size](/flume/seam-proto-b3.png)

*The `check.sh` ladder, $N$ = 24 / 48 / 96, reflux off: $L_\infty$ per footprint, with the cubic (left) and quintic
(right) prolongation (48 → 96 orders in the legend). Shift-inward at 6 or 8 fine cells per block matches the whole
box; the centred two-pass lag6c is slightly more accurate with the cubic fill and identical with the quintic one.*

$N = 192$; orders of 96 → 192:

| Footprint, prolongation | Block size | L1 | L∞ | Order L1 / L∞ |
|---|---|---|---|---|
| shift-inward, cubic | 6 / 8 / 16 / box | $9.2$ to $9.4\cdot10^{-8}$ | $3.7$ to $3.9\cdot10^{-7}$ | 3.98 / 4.00 |
| two-pass, cubic | any | $2.4\cdot10^{-8}$ | $1.9\cdot10^{-7}$ | 4.14 / 4.03 |
| shift-inward, quintic | 6 / 8 / 16 / box | $2.38\cdot10^{-9}$ | $7.83\cdot10^{-9}$ | 4.98 / 4.99 |
| two-pass, quintic | any | $2.42\cdot10^{-9}$ | $7.85\cdot10^{-9}$ | 5.00 / 5.00 |
| shift-inward | 4 | infeasible | | |

Related cases: [C1](../seams/order-across-seams) implements the recommendation; [B1](./caps-prototype) is the base
model; [D1](./quintic-fill) adds the quintic prolongation in FLUME.

## Findings

- **The block size does not matter for shift-inward.** From 6 fine cells per block to the whole box the error changes
  by at most 2.5 % in L1 and 5 % in L∞, at the same order.
- **With the quintic prolongation the two footprints are equivalent** (within 1.5 %), both at order 5.
- **With today's cubic prolongation, two-pass is 4 times more accurate in L1** (2 times in L∞), at the same order 4:
  the one-sided shift has a larger error constant than the centred stencil reading 4th-order ghosts.
- **What each costs in FLUME.** `update_ghost_local` fills copies, coarse-to-fine interpolations and fine-to-coarse
  restrictions in one pass, in map order, the restriction reading only the donor block's interior cells. Shift-inward
  keeps that single pass and needs 6 fine cells per refined axis, which the [#66](https://github.com/szaghi/adam/issues/66)
  guard already requires ($n \ge 2\,n_{gc} = 6$). Two-pass needs the fine ghosts complete before any restriction: a
  second local pass and, on several ranks, a second exchange round, because the sender must update its fine ghosts
  before it packs.
- **Recommendation for Phase C: shift-inward.** It fits the current ghost update and block guard, and with the quintic
  prolongation of Phase D it is as accurate as two-pass. Two-pass would buy a constant factor only while the cubic
  prolongation stays, at the cost of a second synchronisation per ghost update.
