# 2:1 caps in a 2-D prototype

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| B1 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | no ([#85](https://github.com/szaghi/adam/issues/85)) | — |

## Case and scope

#68 Appendix A measured in 1-D that the mean restriction and the reflux of point-value fluxes each cap the order at 2,
and that a 6-point restriction with a quintic prolongation and no reflux recovers order 5. A 1-D model has no
tangential direction: no tangential averaging in the register (the coarse face takes the mean of two fine point fluxes,
one more $O(H^2)$ term, #68 F1), no tangential footprint of the prolongation, no seam corners. Phase B repeats the
measurement in 2-D, in NumPy, without ADAM runs, and locates the maximum error (corner, edge or interior) to test
whether the seam corners explain the L∞ order 1.46 of A1's Euler vortex ([AO](../seams/order-across-seams)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equation | $u_t + u_x + u_y = 0$ (velocity $(1, 1)$) | `Law`, `proto2d.py` |
| Space | conservative finite differences of point values with the 5th-order linear upwind flux (the WENO5 optimal weights), dimension by dimension, 3 ghost layers | `face_flux`, `proto2d.py` |
| Time | classical RK4, one global step $\Delta t = 0.4\,h\,(64h)^{1/4}/\max(a, b)$, so $\Delta t \propto h^{5/4}$ and the spatial error dominates | `advance`, `proto2d.py` |
| Composite grid | the box $[0.25, 0.75]^2$ refined 2:1 ($h = H/2$); the coarse cells it covers hold the restriction of the fine ones before every evaluation, as FLUME's coarse ghosts do | `Composite`, `proto2d.py` |
| Restriction | `mean` (the $2 \times 2$ mean, FLUME's) or `lag6` ($6 \times 6$ tensor Lagrange at the coarse centre, shifted inside the box at its edge) | `Composite.restrict`, `proto2d.py` |
| Prolongation | `cubic` (FLUME's tricubic) or `quintic`, from coarse cells outside the box only (shifted away from the box along the normal axis, as FLUME's donor-block footprint; centred along the other axis) | `Composite.pad_fine`, `proto2d.py` |
| Reflux | on: a coarse box face takes the mean of its two fine fluxes (the 2-D register); off: each grid keeps its own flux | `Composite.rhs`, `proto2d.py` |
| FLUME counterparts | the seam operators of [AO](../seams/order-across-seams) | [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | the periodic unit square |
| Initial condition | $\sin 2\pi x \sin 2\pi y$ |
| Final time | $t = 1$, one period |
| Grid | coarse $N$ cells per side; $N$ = 32 / 64 / 128 (`check.sh`), 64 / 128 / 256 (`--long`) |
| Regions | corner: within 3 coarse cells of a box corner; edge: within 3 coarse cells of a box edge otherwise; interior: the rest |
| Arms | uniform coarse control; restriction (`mean`, `lag6`) × prolongation (`cubic`, `quintic`) × reflux (on, off) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-proto/proto2d.py` | the prototype; `b1 [--cells ...] [--velocity a b] [--assert ARM:NORM:MIN]` |
| `src/tests/flume/verification/amr-proto/check.sh` | runs B1, B2 and B3 (`--long` for the finer ladders) |

```bash
cd src/tests/flume/verification/amr-proto && ./check.sh          # about 1 min for B1
./check.sh --long                                                # B1 at N = 64 / 128 / 256, about 4 min
```

Bounds asserted by `check.sh`, lower bounds only (a transfer that recovers order passes): the uniform control L∞ order
at least 4.8; `lag6`-`quintic`-off at least 4.7 in L1 and L∞; `lag6`-`cubic`-off at least 3.7 in L∞. The caps of the
mean restriction and of the reflux are reported.

## Results

No figure yet ([#85](https://github.com/szaghi/adam/issues/85): order against $N$ per arm, and the corner and edge
errors).

Orders of $N$ = 128 → 256; mass drift after one period at $N = 256$:

| Arm | L1 | L∞ | L∞ corner / edge / interior | Mass drift |
|---|---|---|---|---|
| uniform (coarse) | 5.01 | 5.01 | – / – / 5.01 | $5\cdot10^{-18}$ |
| mean, cubic, on (today) | 2.85 | 1.93 | 1.93 / 1.97 / 1.97 | $1\cdot10^{-18}$ |
| mean, cubic, off | 2.01 | 1.95 | 1.95 / 1.97 / 1.97 | $5\cdot10^{-8}$ |
| lag6, cubic, on | 2.85 | 1.93 | 1.93 / 1.97 / 1.95 | $2\cdot10^{-18}$ |
| lag6, cubic, off | 3.98 | 4.01 | 4.01 / 3.96 / 3.95 | $1.3\cdot10^{-9}$ |
| lag6, quintic, on | 2.85 | 1.93 | 1.93 / 1.95 / 1.96 | $6\cdot10^{-19}$ |
| lag6, quintic, off | **4.99** | **4.99** | 4.99 / 4.98 / 4.98 | $9.1\cdot10^{-11}$ |
| mean, quintic, off | 2.02 | 1.93 | 1.93 / 1.98 / 1.97 | $7\cdot10^{-8}$ |

Related cases: [AO](../seams/order-across-seams) measures the same caps in FLUME; [B2](./sensor-gated-reflux) adds the
nonlinear weights and a sensor-gated reflux; [B3](./high-order-restriction) studies the restriction's footprint.

## Findings

- **The 1-D findings hold in 2-D.** The mean restriction caps the order at 2 whatever the prolongation; the reflux caps
  it at 2 whatever the transfers (its tangential mean of point fluxes included). The 6-point restriction with the
  quintic prolongation and no reflux recovers order 5 in L1 and L∞, with the tangential footprint and the corners
  present; with the cubic prolongation, order 4.
- **Without the reflux, the mass drift converges at the scheme's order** when the transfers are high-order
  ($9\cdot10^{-11}$ at $N = 256$, against $5\cdot10^{-8}$ with the mean restriction).
- **The corners hold the maximum error, not a lower order.** In every arm the maximum sits within 3 cells of a box
  corner, but the corner, edge and interior maxima converge at the same rate. In this linear model the corners do not
  explain the L∞ order 1.46 of A1's Euler vortex; what the model lacks is the nonlinear WENO weights and the
  characteristic projection (#68 F6), which B2's nonlinear prototype brings in.
