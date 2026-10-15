# Sensor-gated reflux prototype

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| B2 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | NumPy prototype (no ADAM run) | yes | — |

## Case and scope

[B1](./caps-prototype) recovers order 5 only without the reflux, and [A3](../seams/shock-crossing) showed that a shock
needs it. Phase E would apply the reflux only where a WENO sensor flags the face (`reflux_policy = troubled`); D-M6-4
builds it only if a prototype shows the order recovered with the shock positions intact. B2 is that go/no-go: does a
sensor-gated reflux recover the smooth order, keep the shocks where the full reflux puts them, and how much
conservation does it give up?

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Composite grid | the [B1](./caps-prototype) grid: the box $[0.25, 0.75]^2$ refined 2:1 on the periodic unit square | `Composite`, `proto2d.py` |
| Space | WENO5-JS with Lax–Friedrichs splitting, the weights of `adam_weno_object` ($\epsilon = 10^{-6}$, exponent 3) | `weno5`, `weno_flux`, `proto2d.py`; [`adam_weno_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_weno_object.F90) |
| Time | classical RK4, $\Delta t = 0.4\,h\,(64h)^{1/4}/\text{speed}$ | `advance`, `proto2d.py` |
| Reflux policies | `always` (today), `never`, `troubled`: a coarse box face takes the mean of its fine fluxes only where the sensor $r = \min_k w_k/d_k$ of the coarse face or of either fine face is below $\tau = 0.2$ (the weno-riemann threshold) | `Composite.rhs`, `proto2d.py` |
| Arms | today's `mean`/`cubic`/`always`, and `lag6`/`quintic` under each policy | `b2`, `proto2d.py` |
| FLUME counterpart | `reflux_policy = troubled` ([E2](./reflux-policy)) | `flag_troubled_seam_cells`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |

## Initial and boundary conditions

| Leg | Law | Initial condition | Crossing | Final time |
|---|---|---|---|---|
| `smooth` | $u_t + u_x + u_y = 0$ | the B1 advection, $\sin 2\pi x \sin 2\pi y$ | — | $t = 1$ |
| `c2f` | Burgers along $x$ | $u = 1$ on $[0, 0.15]$, 0 elsewhere | $x = 0.25$, coarse to fine, at $t = 0.2$ | $t = 0.5$ |
| `f2c` | Burgers along $x$ | $u = 1$ on $[0.45, 0.6]$, 0 elsewhere | $x = 0.75$, fine to coarse, at $t = 0.3$ | $t = 0.6$ |
| `weak` | Burgers along $x$ | a $10^{-3}$ pulse on $[0, 0.15]$ riding on $u = 1$ | $x = 0.25$ | $t = 0.5$ |

The rarefaction catches the shock on the way. The `weak` pulse's smoothness indicators sit near $\epsilon$, so the JS
sensor misses it: the missed-detection control. The pulses start from exact cell averages, so the discrete mass is the
exact one (sampling the step at the centres offsets it by up to a cell, which made the first version of the shock
metric erratic). Grids: smooth leg $N$ = 32 / 64 / 128; shock legs $N = 64$ (`--long`: 64 and 128).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-proto/proto2d.py` | `b2 [--cells ...] [--shock-cells ...] [--leg L ...] [--tau 0.2] [--check]`, L one of `smooth`, `c2f`, `f2c`, `weak` |
| `src/tests/flume/verification/amr-proto/check.sh` | runs B1, B2 and B3 |
| `src/tests/flume/verification/make_doc_figures.py --only seam-proto-b2` | draws the figure (`fig_seam_proto_b2`; NumPy only, about 3 min) |

```bash
cd src/tests/flume/verification/amr-proto && ./check.sh     # B2: about 3 min
./check.sh --long                                           # shock legs at N = 64 and 128, about 5 min
```

**Metrics** against the exact solution: the mass defect over the pulse's own mass, and the mass-equivalent shock
position of every row (the cells within $6H$ of the exact shock), in local cells. **Bounds** (`--check`): the smooth L1
and L∞ orders of `troubled` at least 4.7; on `c2f` and `f2c`, `troubled` within 0.05 cells of `always` and its mass
defect at most $10^{-3}$, while `never`'s is at least $10^{-3}$ (the negative control: the metrics see a missing
reflux).

## Results

![B2: the smooth order and the shock crossings per reflux policy](/flume/seam-proto-b2.png)

*Left: the smooth leg, $L_1$ against $N$ per arm (64 → 128 orders in the legend); `always` is capped near order 3,
`troubled` keeps the uniform order. Centre and right: the Burgers pulses at $N = 64$, the mean shock position error
and the mass defect per policy. `troubled` places the strong shocks exactly where `always` does and keeps the mass
defect below $10^{-3}$ of the pulse; `never` misplaces the c2f shock by a cell and loses $1.8\cdot10^{-2}$; on the
weak pulse the sensor flags nothing, so `troubled` equals `never` (the missed-detection control).*

Smooth leg, orders of $N$ = 64 → 128:

| Arm | L1 | L∞ | L∞ corner / edge | Faces flagged |
|---|---|---|---|---|
| uniform | 5.01 | 4.95 | – | – |
| mean, cubic, always (today) | 3.03 | 2.22 | 2.22 / 2.05 | – |
| lag6, quintic, always | 2.96 | 2.10 | 2.10 / 1.96 | – |
| lag6, quintic, troubled | **5.01** | **4.95** | 5.10 / 5.02 | 0 % |
| lag6, quintic, never | 5.01 | 4.95 | 5.10 / 5.02 | – |

Shock legs (`lag6`/`quintic`; mean shock position in local cells, mass defect over the pulse's mass, the uniform
control in brackets):

| Case, N | always | troubled | never | Faces flagged |
|---|---|---|---|---|
| c2f, 64 | −0.53, $0$ | −0.53, $5.0\cdot10^{-5}$ | −1.12, $1.8\cdot10^{-2}$ | 15 % |
| c2f, 128 (−0.30) | −0.50, $2\cdot10^{-16}$ | −0.50, $3.6\cdot10^{-5}$ | −1.07, $8.6\cdot10^{-3}$ | 8 % |
| f2c, 64 | −0.25, $2\cdot10^{-16}$ | −0.25, $6.5\cdot10^{-5}$ | +0.05, $1.7\cdot10^{-2}$ | 15 % |
| f2c, 128 (−0.37) | −0.25, $4\cdot10^{-16}$ | −0.25, $6.9\cdot10^{-6}$ | +0.06, $8.4\cdot10^{-3}$ | 7 % |
| weak, 64 (+0.02) | −0.003, $10^{-12}$ | −0.014, $3.9\cdot10^{-4}$ | −0.014, $3.9\cdot10^{-4}$ | 0 % |
| weak, 128 (+0.01) | −0.015, $10^{-12}$ | −0.022, $3.2\cdot10^{-5}$ | −0.022, $3.2\cdot10^{-5}$ | 0 % |

Related cases: [B1](./caps-prototype) is the linear prototype; [E2](./reflux-policy) measures the policy in FLUME;
[AL](../seams/integral-drift) measured the sensor's false positives on the Euler vortex.

## Findings

- **`troubled` recovers the order.** On the smooth leg the sensor flags no seam face, so `troubled` is `never`: order 5
  in L1 and L∞, as the uniform control, against about 2 with the reflux always on.
- **`troubled` keeps the shock where `always` puts it.** The two differ by at most 0.003 cells; `never` moves the shock
  by 0.3 to 0.6 cells and loses $10^{-2}$ of the pulse's mass. The reflux at the flagged faces (7 to 15 % of the
  seam-face evaluations) is what holds the shock.
- **`troubled` is near-conservative, not conservative.** Its mass defect, $7\cdot10^{-6}$ to $6.5\cdot10^{-5}$ of the
  pulse, comes from the unflagged smooth faces (the rarefaction crossing the seam) and from the steps before the shock
  reaches the sensor's stencil; it is two to three orders below `never`'s.
- **A missed detection is a small error here.** The weak shock is never flagged, so `troubled` runs as `never`: its
  shock moves by 0.007 to 0.011 cells against `always` and it loses $3\cdot10^{-5}$ to $4\cdot10^{-4}$ of its own mass.
  What the JS sensor misses is what the reflux matters least for, but not nothing.
- **The corners are again not the L∞ cause.** Under the nonlinear weights corner and edge maxima still converge at the
  same rate; A1's Euler L∞ order 1.46 is not reproduced by a scalar model, and stays open for the system (the
  characteristic projection).
- **Verdict (D-M6-4): go for Phase E**, with three conditions on record. The order is recovered only together with the
  high-order transfers of Phases C and D: with today's mean restriction the order stays at 2 whatever the reflux does
  (B1). The policy is near-conservative ($10^{-5}$ of the shocked mass), not exact. And A4
  ([AL](../seams/integral-drift)) measured 5 % false positives on the Euler vortex with the sensor read on primitive
  variables, where this scalar model flags none; E must measure its own false-positive rate, since each flagged smooth
  face brings the order cap back.
