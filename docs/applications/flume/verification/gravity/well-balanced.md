# The well-balanced mode

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GV-2b, GV-5, GV-8, GV-10, GV-7 (G2b) | External gravity | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes | — |

## Case and scope

Without balancing, the truncation errors of the flux difference and of the gravity source do not cancel on a
stratified atmosphere at rest: they drive a stationary flow (GV-2: Mach $6\cdot10^{-4}$ on ten scale heights at
$N = 32$, $10^{-2}$ across a 2:1 seam) that swamps any perturbation smaller than it. G2b keeps a known equilibrium as
an exact discrete steady state by residual subtraction, the deviation method (Edelmann et al. 2021, §4.3).

What the case catches:
- **an equilibrium that is not an exact discrete steady state** (GV-2b, GV-5): isothermal, polytropic and magnetised
  atmospheres at rest must stay at rest, $\max|\mathbf u|/c < 10^{-12}$, on both face families;
- **a seam that breaks the balance** (GV-8): across 2:1 seams the residual and the reflux must both cancel the
  equilibrium's share;
- **a balance that hides the physics** (GV-10): the response to a small pressure pulse must be linear in its
  amplitude, the purpose of the mode;
- **backend divergence** (GV-7) on the balanced pulse and the balanced AMR box.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$; GV-5: ideal MHD with EGLM (`glm_ch = 3`, `glm_alpha = 0.18`, `glm_damping_length = 1`) | [models](../../models#external-gravity), [EGLM](../../models#divergence-control-eglm) |
| Well-balanced mode | `[gravity] equilibrium = atmosphere`: the initial atmosphere at rest $q_e$ (zero velocity, region 1's field, no pulse); its residual $R(q_e)$, computed with the run's own space operator (ghosts, boundary conditions, seam fill, limiter, every source) once per grid (at the first step, after a restart, after every regrid), is subtracted from every stage residual | `refresh_equilibrium_residual`, `subtract_equilibrium_residual`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `set_equilibrium`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90); `subtract_equilibrium_residual_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90); `refresh_equilibrium_residual`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Reflux at 2:1 seams | the equilibrium's coarse–fine flux mismatch is stored per register face (`seam_e`) and the end-of-step reflux subtracts it | `apply_reflux_to_stage_forest`, `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Refusals | multi-realm runs; `reflux_policy = troubled` | `initialize_flume`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `check_gravity`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Gravity sources | momentum $\rho\mathbf g$, fourth-order energy source | `add_gravity_sources`, `compute_gravity_correction`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Space, `weno-riemann` (GV-5) | WENO5 interpolation, HLLD, sixth-order flux correction | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_mhd_eglm_hlld_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_hlld_kernels.F90); [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | x `periodic`, y `hydrostatic` | `hydrostatic_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90) |

## Initial and boundary conditions

Common: $[0, 1]^2$, null z, quadtree 4×4 blocks of $8^2$ cells ($N = 32$), two uniform levels, `ngc = 4`, np 2,
CFL 0.4, region 1 $(\rho_1, p_1) = (1, 1)$ at rest at the origin, x `periodic`, y `hydrostatic`; initial condition
[`isothermal-atmosphere` or `polytropic-atmosphere`](../../initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity).

| Leg | Case | Atmosphere | Gravity | Extra | Final time |
|---|---|---|---|---|---|
| GV-2b | `strat2d` | isothermal, ten scale heights (density ratio $e^{10} = 2.2\cdot10^4$) | $(0, -10)$ | — | 0.5 |
| GV-2b | `polystrat2d` | polytropic, $G = 1.2$ (density ratio 243) | $(0, -4)$ | — | 0.5 |
| GV-5 | `strat2d-mhdy` | as `strat2d`, EGLM | $(0, -10)$ | $\mathbf B = (0, 0.05, 0)$, along $\mathbf g$ | 0.2 |
| GV-5 | `strat2d-mhdx` | as `strat2d`, EGLM | $(0, -10)$ | $\mathbf B = (0.05, 0, 0)$, across $\mathbf g$ | 0.2 |
| GV-8 | `strat2d-amr` | as `strat2d` | $(0, -10)$ | the box $[0.25, 0.75]^2$ refined one level (2:1 seams across the stratification; 28 blocks) | 0.5 |
| GV-10 | `strat2d` + pulse | as `strat2d` | $(0, -10)$ | relative pressure pulse $p\,(1 + A e^{-|\mathbf x - \mathbf x_0|^2/w^2})$, $A = 10^{-6}$ and $10^{-4}$, $\mathbf x_0 = (0.5, 0.5)$, $w = 0.05$ | 0.1 |

With $|\mathbf B| = 0.05$ the plasma $\beta = 2p/|\mathbf B|^2$ is 800 at the floor and $0.036$ at the top
($p = e^{-10}$). The balanced legs run with `--balanced`; the unbalanced runs of `strat2d` and `strat2d-amr` are
reported for comparison.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes the cases (`--balanced`, `--pulse A`, `--time-max`, `--riemann hlld`) |
| `src/tests/flume/verification/gravity/check.sh` | runs the G2b legs with the force legs; `--agree TAG` adds GV-7 |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | `rest` ($\max|\mathbf u|/c$, `--max 1e-12`), `linear` (GV-10, $10^{-2}$), `agree` (GV-7, $10^{-10}$) |
| `src/tests/flume/verification/make_doc_runs.sh` (group `gravity`) | re-runs `{strat2d,polystrat2d,strat2d-amr}-{unb,bal}` saved every 5 steps and `lin-{unb,bal}-{6,4}` |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_gravity_balance`) | draws the figure |

```bash
cd src/tests/flume/verification/gravity
./check.sh --np 2 --cases iso2d --keep                                   # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np 2 --cases iso2d --keep \
   --agree adam_flume_cpu                                                # FNL, and GV-7 against the CPU runs
```

## Results

![GV (G2b): spurious flows unbalanced and balanced, the refined box, the pulse responses](/flume/gravity-balance.png)

*Top: $\max|\mathbf u|/c$ over the run, and the refined box at $t = 0.5$ unbalanced and balanced; bottom: the pulse
responses at $t = 0.1$, balanced and unbalanced, and the balanced $10^{-6}$ response field.*

Top left, $\max|\mathbf u|/c$ over the run: the unbalanced atmospheres reach their stationary flow within a few steps
(`strat2d` levels at $6\cdot10^{-4}$, `polystrat2d` at $1.5\cdot10^{-3}$), the refined box keeps growing to
$1.2\cdot10^{-2}$; balanced, `strat2d` and `polystrat2d` stay exactly at rest (drawn at $10^{-18}$) and the refined box
at round-off. Top centre, the unbalanced refined box at $t = 0.5$: the flow peaks at $(0.52, 0.80)$, in the second
coarse row above the top seam ($y = 0.75$), where the coarse and fine truncation errors meet. Top right, the same run
balanced: round-off, largest in the fine cells ($3.3\cdot10^{-15}$), exactly zero in parts of the coarse rows. Bottom,
the pulse at $t = 0.1$ along $x = 0.484$: balanced, the responses to $10^{-6}$ and $10^{-4}$ coincide (peak
$(p/p_e - 1)/A = 0.40$); unbalanced, the curve is the stationary truncation flow divided by $A$ (peak 333 for
$10^{-6}$, 3.4 for $10^{-4}$) and hides the pulse entirely. Bottom right, the balanced $10^{-6}$ response: the acoustic
front, stronger upward, where the density falls.

| Leg | Asserted | Bound | Unbalanced | Balanced, CPU | Balanced, FNL |
|---|---|---|---|---|---|
| GV-2b `strat2d`, $t = 0.5$ | $\max|\mathbf u|/c$ | $< 10^{-12}$ | $5.9\cdot10^{-4}$ | 0 | 0 |
| GV-2b `polystrat2d`, $t = 0.5$ | as above | $< 10^{-12}$ | $1.5\cdot10^{-3}$ | 0 | 0 |
| GV-5 `strat2d-mhdy`, `weno` / `weno-riemann` HLLD | as above | $< 10^{-12}$ | $5.9\cdot10^{-4}$ | 0 / 0 | 0 / 0 |
| GV-5 `strat2d-mhdx`, `weno` / `weno-riemann` HLLD | as above | $< 10^{-12}$ | $2.3\cdot10^{-4}$ | 0 / 0 | 0 / 0 |
| GV-8 `strat2d-amr`, $t = 0.5$ | as above | $< 10^{-12}$ | $1.2\cdot10^{-2}$ | $3.3\cdot10^{-15}$ | $2.8\cdot10^{-15}$ |
| GV-10 pulses $10^{-6}$, $10^{-4}$ on `strat2d`, $t = 0.1$ | responses per unit amplitude agree | $10^{-2}$ | differ by 97 (reported) | $5.7\cdot10^{-5}$ | the same |
| GV-7 FNL against CPU, the $10^{-4}$ pulse and `strat2d-amr` | per field | $\le 10^{-10}$ | — | — | $\le 7.8\cdot10^{-16}$ |

The unbalanced MHD values are at $t = 0.5$, from the trial; `check.sh` runs the MHD legs to $t = 0.2$.

Related cases: [stratified atmospheres](./atmospheres) (the equilibria and their ghosts); [the gravity force](./force)
(GV-2, the unbalanced baseline); [Rayleigh–Taylor](./rayleigh-taylor) keeps its unseeded layers with the same mode
(GV-4b); [conservation across AMR seams](../euler/amr-conservation) for the reflux without gravity.

## Findings

- The balanced rest is exactly zero where the state equals the stored equilibrium bit for bit; on the AMR box the
  reflux arithmetic leaves round-off.
- GV-10 shows what the mode is for: the balanced response is linear in the amplitude to its second order
  ($5.7\cdot10^{-5}$ against $A = 10^{-4}$), while the unbalanced one is dominated by its own stationary flow.
- **The seam mismatch belongs to the reflux, not to the residual.** At $q = q_e$ the residual and the subtracted
  $R(q_e)$ are the same operator on the same state, at 2:1 seams too. The end-of-step reflux would still add the
  equilibrium's coarse–fine flux mismatch, so it is stored per register face and the reflux subtracts it. Folding it
  into $R(q_e)$ instead was tried and fails: the coarse seam cells then carry a stage residual, the stages leave the
  equilibrium and the fine ghosts follow (Mach $6\cdot10^{-5}$ on the AMR box; exact with the reflux off). The residual
  and the reflux each conserve for the state and for the equilibrium, so the difference conserves.
- **GV-7 momentum scale.** GV-7 scales the momentum by at least $\max\rho c$: a balanced run's momentum is the
  perturbation, round-off at rest, and the backends' round-off of the $O(1)$ residuals that the subtraction cancels
  ($\le 6\cdot10^{-15}\,\rho c$) would read $O(1)$ against its own maximum.
- **The alternative.** Xing & Shu 2013 and Li & Xing 2018 (a modified Lax–Friedrichs split with the source
  reconstructed by the same WENO operator) was prototyped in NumPy and matched the deviation method to three digits,
  but covers the flux-splitting path only.
- **Limits.** The subtracted residual is outside the positivity limiter's factor: at the equilibrium the factor is 1
  and the balance exact, but the limiter's guarantee does not cover that term (the equilibrium's truncation error,
  small where the atmosphere is resolved). Refused: multi-realm runs (the peer realms' equilibrium would be needed at
  the seams) and `reflux_policy = troubled` (its reflux cells change per step). Cost: one field and one residual
  evaluation per grid.
- **Found on the way** (FNL, nvfortran 26.1 with FiNeR): a character `get` right after the failed `get` of an absent
  numeric key returns an empty value with no error; `equilibrium = atmosphere` read as "" after the absent `froude`.
  The gravity section now reads its string keys first (`load_from_file`,
  [`adam_flume_gravity_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_gravity_object.F90)).
