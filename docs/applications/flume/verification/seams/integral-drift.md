# Integral drift and the WENO sensor near seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| AL (A4) | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (the two backends log the same values) | yes | — |

## Case and scope

The reflux policy of Phase E (`reflux_policy = troubled`) turns the reflux off at the seam faces that a WENO sensor
judges smooth. Two quantities decide whether that is safe and where its threshold sits, and no check logged either:
- how the conserved integrals drift when the reflux is off (bounded, or accumulating?);
- what the sensor reads on the seam faces of a smooth flow (false positives) and of a shock crossing (detections).

AL logs both on two legs, each run with the reflux on and off. It asserts only the conservation of the reflux; the rest
is a log.

What the case catches:
- **a reflux that stops conserving:** with the reflux on, every logged integral must drift at most $10^{-13}$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `vortex`: SSP Runge–Kutta 5 stages (`runge-kutta-ssp-54`), CFL 0.4; `sod`: SSP Runge–Kutta 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Seam operators | `mean` restriction, `tricubic` fill | [numerics](../../numerics#adaptive-mesh-refinement) |
| Reflux | Berger–Colella on, or off (`--no-reflux`, a measurement, not a supported configuration) | `apply_reflux_to_stage_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Conservation history | `[diagnostics] conservation_history_save` | `compute_conservation`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Sensor (post-processing) | $r = \min_k w_k/d_k$ of the WENO-5 weights of `adam_weno_object` (Jiang–Shu, $\epsilon = 10^{-6}$, exponent 3), the weno-riemann sensor of [#47](https://github.com/szaghi/adam/issues/47), threshold $\tau = 0.2$ | `amr_log.py`; [`adam_weno_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_weno_object.F90); `FLUX_CORRECTION_SENSOR_TAU`, [`adam_flume_parameters.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_parameters.F90) |
| Boundaries | `vortex`: periodic in $x$, $y$; `sod`: `extrapolation` | [boundary conditions](../../boundary-conditions) |

## Initial and boundary conditions

| Leg | Setup | Conserved integrals logged |
|---|---|---|
| `vortex` | the [AO](./order-across-seams) Euler vortex: quadtree, quadrant $[0.5, 1]^2$ refined, $N = 64$, periodic, run to $t = 1$, one period of the free stream $(1, 1)$, so the vortex crosses every seam | $\rho$, $\rho u$, $\rho v$, $E$ (every integral is conserved) |
| `sod` | the [AS](./shock-crossing) Sod `c2f`: octree, $x > 0.75$ refined, $N = 200$, $t = 0.2$; the shock crosses the 2:1 face at $t = 0.143$ | $\rho$, $E$ (no flux through the ends before $t = 0.2$) |

Both legs run with the reflux on and off, a checkpoint every 10 steps, 2 ranks.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-order/make_amr_order.py` | the `vortex` inputs (`--cells 64 --time-max 1.0 --it-save 10 [--no-reflux]`) |
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | the `sod` inputs (`--cells 200 --side c2f --it-save 10 [--no-reflux]`) |
| `src/tests/flume/verification/amr-log/amr_log.py` | `drift` and `sensor` logs |
| `src/tests/flume/verification/make_doc_figures.py --only integral-drift` | the figure, from the kept work directories of `amr-log/check.sh` (CPU) |
| `src/tests/flume/verification/amr-log/check.sh` | runs the legs (`--legs vortex,sod`, `--np N`, `--keep`) |

```bash
cd src/tests/flume/verification/amr-log && ./check.sh             # CPU: about 11 min
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh       # FNL: about 5 min
```

**Logs** (`amr_log.py`):
- **drift:** the relative drift $(I(t) - I(0))/|I(0)|$ of the integrals of the conservation history, its maximum, and
  its growth (the mean $|$drift$|$ over the last quarter of the run over that of the second quarter: about 1 when
  bounded, $7/3$ for a linear drift; "a jump" when the second quarter holds less than 1 % of the maximum, as before a
  shock first crosses a face);
- **sensor:** on every checkpoint, $r = \min_k w_k/d_k$ of every face along every non-null axis, the minimum over
  $\rho$, $u_n$, $p$ (primitive variables) and over the two reconstructions beside the face, recomputed from the
  checkpoints, whose ghosts are filled right before the save as the stencils read them. Per region (the seam band, the
  faces within `ngc` of their own spacing from the 2:1 faces; the interior): the minimum $r$ and the fraction of faces
  below $\tau = 0.2$.

The default scheme weights its characteristic split fluxes, so the map is the sensor a seam policy would read, not a
dump of the scheme's own weights.

## Results

![AL: the integrals across static seams](/flume/integral-drift.png)

*CPU, from the work directories `amr-log/check.sh --keep` and `--reflux-policy troubled --keep` leave: the relative drift of mass and energy per step. The reflux keeps both at round-off; without it the vortex drifts to $2\cdot10^{-5}$ and the Sod crossing to $2\cdot10^{-4}$; under `troubled` the smooth vortex drifts as without reflux, so the sensor leaves its seam faces unrefluxed, the Sod crossing only to $2\cdot10^{-7}$ (the shock faces are flagged).*

Drift (2 ranks):

| Leg | Bound, reflux on | Reflux on: max drift | Reflux off: max drift (mass, energy) | Reflux off: growth |
|---|---|---|---|---|
| vortex | $10^{-13}$ | $\le 4.5\cdot10^{-16}$ | $2.2\cdot10^{-5}$, $2.5\cdot10^{-5}$ (momentum $3.9$, $2.7\cdot10^{-5}$) | 1.4 to 7.7 |
| sod | $10^{-13}$ | $\le 3.3\cdot10^{-16}$ | $1.2\cdot10^{-4}$, $2.4\cdot10^{-4}$ | a jump at the crossing |

Sensor, faces with $r < \tau$:

| Leg | Seam band: largest (mean) fraction | Interior: largest (mean) fraction |
|---|---|---|
| vortex | 0.110 (0.053) | 0.051 (0.020) |
| sod | 1.000 (0.172) | 0.128 (0.078) |

Related cases: [E2](../seams-studies/reflux-policy) builds the `troubled` policy on these logs;
[B2](../seams-studies/sensor-gated-reflux) is its scalar prototype; the reflux conservation on uniform-level seams is
[conservation across AMR seams](../euler/amr-conservation).

## Findings

- **Without the reflux, a smooth flow's drift accumulates.** Every seam crossing adds to it: by one period the vortex
  has lost about $2\cdot10^{-5}$ of its mass and energy, and the growth (1.4 to 7.7) says the drift is not bounded. A
  shock adds $10^{-4}$ in one crossing. The reflux keeps both at round-off.
- **The sensor detects the crossing shock.** In the Sod run the seam band reads $r = 1$ until the shock reaches it,
  67 to 100 % of its faces are below $\tau$ while the shock crosses ($t = 0.126$ to 0.152), and 0 to 13 % after,
  where the A3 defect sits.
- **It also flags smooth seam faces.** On the vortex, up to 11 % of the band faces (mean 5 %) read $r < \tau$, 2.6
  times the interior's mean fraction: the JS weights drop at the extrema of a smooth profile, and more so where the
  seam's second-order error adds kinks. For `reflux_policy = troubled` a false positive keeps the reflux, the
  conservative choice: the cost is the order on those faces, which Phase B2 measures. The threshold is not tuned here.
