# Reflux policy: troubled against always

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| E2 (E1) | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL (flags checked equal to the CPU's); AO `euler-tq`: CPU, FNL | yes | — |

## Case and scope

The Berger–Colella reflux keeps a 2:1 seam conservative, and caps the composite at second order
([AO](../seams/order-across-seams), [B1](./caps-prototype)): it replaces the coarse face flux by the restricted fine one
everywhere. [B2](./sensor-gated-reflux) showed in a prototype that a sensor-gated reflux recovers the order and keeps
the shocks in place (D-M6-4: go). E1 builds `[numerics] reflux_policy = troubled`; E2 measures it against `always` on
the four checks of the plan: the smooth order (A1), the shock crossings (A3), the integral drift (A4) and the seam
divergence (A2).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| `reflux_policy = always` | default: every seam face cell takes the restricted fine flux | `apply_reflux_to_stage_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| `reflux_policy = troubled` | the reflux applies only at the seam face cells flagged in a stage of the step: those the positivity limiter blended, and those where the WENO-5 sensor $\min_k w_k/d_k$ of the density, the normal velocity or the pressure across the face is below 0.2; a cell flagged at any stage is refluxed | `flag_troubled_seam_cells`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); `seam_weno5_sensor`, [`adam_flume_seam_sync_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_seam_sync_object.F90); `FLUX_CORRECTION_SENSOR_TAU`, [`adam_flume_parameters.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_parameters.F90) |
| Key | `[numerics] reflux_policy`, with `reflux = .true.` | [`adam_flume_numerics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_numerics_object.F90) |
| Seam operators | default (`mean` + `tricubic`) or `lagrange6` + `quintic` | [C1](../seams/order-across-seams), [D1](./quintic-fill) |
| Scheme, time | as in each baseline case | [numerics](../../numerics#adaptive-mesh-refinement) |

## Initial and boundary conditions

| Check | Setup |
|---|---|
| A1 | the [AO](../seams/order-across-seams) euler and mhd ladders, $N$ = 64 / 128 / 256; euler also with `lagrange6` + `quintic` (`euler-tq`) |
| A3 | the [AS](../seams/shock-crossing) crossings, Sod and RJ2a, `c2f` and `f2c` |
| A4 | the [AL](../seams/integral-drift) vortex ($N = 64$ to $t = 1$) and Sod `c2f` |
| A2 | the [AD](../seams/divergence-across-seams) `cpaw`, $N$ = 64 / 128 |

FNL, 2 ranks.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/check.sh --legs euler-tq` | the AO leg under `troubled` with `lagrange6` + `quintic`: orders $\ge 5.5$, errors at most the `euler-q6` baseline plus 2 % |
| `src/tests/flume/verification/amr-order/make_amr_order.py --reflux-policy troubled` | the other AO runs |
| `src/tests/flume/verification/amr-shock/make_amr_shock.py --reflux-policy troubled` | the crossings |
| `src/tests/flume/verification/make_doc_runs.sh <runs> seams`, `make_doc_figures.py --only reflux-policy` | re-runs the crossings `e2-troubled-{c2f,f2c}` and `d-quintic-{c2f,f2c}` (CPU); the AO ladders from the kept work directories |

```bash
cd src/tests/flume/verification/amr-order && ./check.sh --legs euler-tq
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs euler-tq
```

The A2 and A4 runs under `troubled`: `amr-divb/make_amr_divb.py --reflux-policy troubled` and
`amr-log/check.sh --reflux-policy troubled` (issue #85 R2b; before, the policy was added to the inputs by hand).

## Results

![E2: troubled against always](/flume/reflux-policy.png)

*CPU. Left: the AO ladder, `troubled` with lagrange6 + quintic at the order of the run without reflux (6.11), today's operators at 2. Centre and right: the Sod shock after the crossing, lagrange6 + quintic under `always` and `troubled`: the same profile (the shock faces are flagged, so they are refluxed).*

| Check | `always` | `troubled` |
|---|---|---|
| AO euler, default operators, N = 256: L1, order | $4.05\cdot10^{-6}$, 2.01 | $1.80\cdot10^{-6}$, 2.34 |
| AO euler, `lagrange6` + `quintic`, N = 256: L1, order | $4.05\cdot10^{-6}$, 2.01 | $6.03\cdot10^{-8}$, **6.11** |
| AO mhd, default operators, N = 256: L1, order | $1.83\cdot10^{-6}$, 2.07 | $1.10\cdot10^{-6}$, 2.04 |
| AO face cells flagged, N = 64 / 128 / 256 | — | 3–5 % / 0–0.6 % / 0 % |
| A3 shock position (Sod, RJ2a, c2f and f2c) | reference | equal within $2\cdot10^{-4}$ cells |
| A3 plateau and band defects | reference | within 1 % (RJ2a c2f band +44 %) |
| A3 mass drift | $\le 3\cdot10^{-10}$ | $10^{-11}$ to $1.5\cdot10^{-7}$ |
| A4 drift, vortex N = 64 to $t = 1$ (r, ru, rv, rE) | $\le 3\cdot10^{-16}$ | $1.8$–$4.4\cdot10^{-5}$, growth 1.3–1.6 (rv 8.8) |
| A4 drift, Sod c2f | $3\cdot10^{-16}$ | $1.9\cdot10^{-7}$ |
| CPAW D_flux seam floor, N = 64 / 128 (growth) | $5.39\cdot10^{-4}$ / $5.43\cdot10^{-5}$ (0.95 / 0.72) | $5.52\cdot10^{-4}$ / $3.80\cdot10^{-5}$ (0.76 / 0.68) |

On FNL the flags equal the CPU's to the cell.

Related cases: [B2](./sensor-gated-reflux) is the prototype; [AL](../seams/integral-drift) measured the sensor's false
positives; [G3a](./reflux-not-ghosts) uses `troubled` as the #78 workaround.

## Findings

- **The reflux was what capped the high-order seam operators.** With `lagrange6` and `quintic` the reflux on
  (`always`) holds the vortex at second order, $4.05\cdot10^{-6}$; under `troubled` no smooth face is flagged once
  resolved and the run is sixth order, $6.03\cdot10^{-8}$, the error of the uniform coarse run (67 times smaller). AO
  leg `euler-tq` asserts it.
- **Shocks are unchanged:** positions equal to $2\cdot10^{-4}$ cells, defects within 1 % but one band value; the flags
  cover 35 to 42 % of the face-cell steps of a crossing.
- **Conservation is the price:** an unflagged face keeps its flux mismatch, so the integrals drift at the reflux-off
  rate on smooth flows ($10^{-5}$ per period at $N = 64$, bounded but one component accumulating) and by $10^{-7}$
  across a shock.
- **MHD:** the D_flux seam floor is bounded and 1.4 times lower at $N = 128$ (83 against 119 times the uniform floor).
- **Decision (D-M6-2, owner):** the default stays `always`, exactly conservative; `troubled` is the opt-in that makes
  the high-order seam operators worth their cost with the reflux on.
