# 2:1 seam div(B): truncation order

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-amr-fd`, `rmf-amr-fd-pulse` | [#21](https://github.com/szaghi/adam/issues/21), [#22](https://github.com/szaghi/adam/issues/22), [#25](https://github.com/szaghi/adam/issues/25), [#26](https://github.com/szaghi/adam/issues/26), [#29](https://github.com/szaghi/adam/issues/29), [#145](https://github.com/szaghi/adam/issues/145) ([#99](https://github.com/szaghi/adam/issues/99) open) | P | CPU, FNL | no | `rmf-amr-fd`, `rmf-amr-fd-pulse` (digest + residuals, CPU and FNL) |

## Case and scope

On a uniform grid $\nabla_h\cdot\mathbf B$ stays at round-off ([uniform grid](./uniform-grid)). At a 2:1 face the
coarse and fine sides difference with stencils of different resolution, the curl and the divergence no longer commute,
and the seam injects a $\nabla\cdot\mathbf B$ source. The source is truncation-order, $O(h^p)$: it converges under
refinement and does not decay at fixed $h$ (hyperbolic cleaning with damping keeps it low, below). Five local fixes were ruled out and the resolution is accept-truncation
([hazards](../hazards#seam-div-b), [#29](https://github.com/szaghi/adam/issues/29)).

The two cases pin that state. The historical pointwise criterion "seam div(B) ≤ 1e-13" is retired
([#21](https://github.com/szaghi/adam/issues/21)). The acceptance is now:
- the identity holds wherever there is no resolution jump (marker-off control, and a grid refined uniformly so that
  every block face is 1:1);
- the seam $\nabla\cdot\mathbf B$ sits on a pinned baseline, within 5 %: a rise is a regression of the seam fill or the
  exchange, a drop is an improvement to rebaseline deliberately;
- the seam source converges under refinement (opt-in leg).

`rmf-amr-fd` carries the four coils. Its coil filament makes $\nabla\cdot\mathbf J$ scale as $h^{-2}$, so its seam
metrics are contaminated by the source and are no order measurement. `rmf-amr-fd-pulse` is the source-free twin: with
$\mathbf J = 0$ both divergences are conserved, and any signal comes from the seam alone. Its pulse travels along $+y$,
transverse to the seam normal: an $x$-travelling pulse varies only in $x$ and is blind to the seam.

What the cases catch:
- a seam fill or a ghost exchange that changes the seam source (baseline legs);
- a 1:1 block face that leaks (the uniform-refined leg isolates the jump as the only source);
- a ghost fill that no longer reproduces constants: $\nabla\cdot\mathbf D$ of the $+y$ pulse is structurally zero
  across the seam;
- an untruthful $\nabla\cdot\mathbf J$ diagnostic ([#26](https://github.com/szaghi/adam/issues/26): unstamped J
  ghosts and a stamp lagged by one step put it six orders off);
- an RK scheme the staged forest path cannot split, running instead of refusing
  ([#25](https://github.com/szaghi/adam/issues/25)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI; 4 AC coils (`rmf-amr-fd`) or no source (`rmf-amr-fd-pulse`) | [models](../models#models) |
| Space | centred finite differences of order 6 | `compute_curl_fd_centered`, `compute_divergence_fd_centered`, [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5; SSP (1,1) and `runge-kutta-1` in the RK legs | [`adam_stepper_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_stepper_object.F90) (SELAH) |
| Seam ghosts | coarse→fine tricubic interpolation over the full ghost slab (`[amr] seam_ghost_fill`, default `tricubic`), fine→coarse restriction | `seam_interpolate_tricubic`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90); [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| RK contract | the staged forest path accepts SSP schemes only | `stages_per_step_forest`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90) |
| Diagnostic | max $\lvert\nabla\cdot\mathbf D\rvert$, $\lvert\nabla\cdot\mathbf B\rvert$, $\lvert\nabla\cdot\mathbf J\rvert$ per step, `-divergence_history.dat` | `save_divergence_history`, [`adam_prism_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_common_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-0.06, 0.06]^3$ m |
| Grid | octree, 8 blocks of $16^3$ at level 1 ($\Delta x = 3.75\cdot10^{-3}$); the $x < 0$ blocks refined to level 2 ($\Delta x = 1.875\cdot10^{-3}$): one 2:1 plane at $x = 0$; `ngc = 3` |
| Initial condition | `rmf-amr-fd`: `vacuum` with 4 AC coils. `rmf-amr-fd-pulse`: [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `r0 = 0`, `sigma = 0.01`, `k0 = 0`, `B0 = 10` ($D_z$, $B_x$) |
| Boundaries | `Neumann` on all six faces (the pulse moves $2.4\cdot10^{-3}$ m in 5 steps and never reaches them) |
| Steps, ranks | 5, 1 rank |
| Control | `markers_number = 0` (no seam) |
| Uniform leg | `box_xmax = 0.06`: the whole domain at level 2, every face 1:1 |
| Convergence leg | `ni = 16`, 5 steps against `ni = 32`, 10 steps: the same final time |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-amr-fd/input.ini`, `check.sh` | the coil case: control, seam $\nabla\cdot\mathbf B$, seam $\nabla\cdot\mathbf J$ |
| `src/tests/prism/regression/rmf-amr-fd-pulse/input.ini`, `check.sh` | the pulse: control, uniform, seam, RK legs, quadtree leg ([quadtree](./quadtree)), `--convergence`, `--cleaning` |

```bash
cd src/tests/prism/regression/rmf-amr-fd && ./check.sh
cd ../rmf-amr-fd-pulse && ./check.sh --convergence
./check.sh --cleaning                                               # the seam with hyperbolic cleaning of B
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

Bounds as coded; the baselines are CPU-pinned and the 5 % band is the cross-backend and cross-compiler allowance.

| Case, leg | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| `rmf-amr-fd`, control | max $\lvert\nabla\cdot\mathbf B\rvert$ | $\le 10^{-13}$ | $1.96\cdot10^{-19}$ | $1.03\cdot10^{-19}$ |
| `rmf-amr-fd`, seam | max $\lvert\nabla\cdot\mathbf B\rvert$ | $1.555802\cdot10^{-7} \pm 5\%$ | $1.555802\cdot10^{-7}$ (baseline) | $1.555801\cdot10^{-7}$ |
| `rmf-amr-fd`, seam | max $\lvert\nabla\cdot\mathbf J\rvert$ | $3.217958\cdot10^{-5}$, factor 3 | $3.22\cdot10^{-5}$ | $2.62\cdot10^{-5}$ |
| pulse, control | max $\lvert\nabla\cdot\mathbf D\rvert$, max $\lvert\nabla\cdot\mathbf B\rvert$ | $\le 10^{-13}$ | 0, 0 | 0, 0 |
| pulse, uniform level 2 (1:1 faces) | max $\lvert\nabla\cdot\mathbf D\rvert$, max $\lvert\nabla\cdot\mathbf B\rvert$ | $\le 10^{-13}$ | 0, 0 | 0, 0 |
| pulse, seam | max $\lvert\nabla\cdot\mathbf D\rvert$ | $\le 10^{-13}$ | $3.3\cdot10^{-15}$ | $3.3\cdot10^{-15}$ |
| pulse, seam | max $\lvert\nabla\cdot\mathbf B\rvert$ (step 5) | $12.65004 \pm 5\%$ | 12.65004 | 12.65004 |
| pulse, `runge-kutta-1` | non-zero exit with `not stage-splittable` in the log | refusal | refused, exit 1 | refused, exit 1 |
| pulse, SSP (1,1) | max $\lvert\nabla\cdot\mathbf D\rvert$; max $\lvert\nabla\cdot\mathbf B\rvert$ | $\le 10^{-13}$; $11.81251 \pm 5\%$ | $2.3\cdot10^{-15}$; 11.81251 | $3.3\cdot10^{-15}$; 11.81251 |
| pulse, `--convergence` | $p_{obs} = \log_2(E_{16}/E_{32})$ of the seam $\nabla\cdot\mathbf B$ at the same time | $\ge 0.8$ | $E_{16} = 12.65$, $E_{32} = 5.700$: 1.150 | $E_{16} = 12.65$, $E_{32} = 5.700$: 1.150 |
| pulse, `--cleaning` (hyperbolic cleaning of $\mathbf B$, $\chi = 2$, $\alpha = 6.4\cdot10^{10}$ s$^{-1}$, 10 steps) | max $\lvert\nabla\cdot\mathbf D\rvert$; seam $\lvert\nabla\cdot\mathbf B\rvert$ at the final time of the seam leg | $\le 10^{-13}$; $6.637193 \pm 5\%$ and below the seam leg | $4.3\cdot10^{-15}$; 6.637193 | $3.0\cdot10^{-15}$; 6.637193 |

The second-rung ratio of the refinement ladder, at `ni = 64`, is 1.62, toward the predicted 2.

### Long horizon and hyperbolic cleaning

Measurement, not a leg ([#145](https://github.com/szaghi/adam/issues/145)). The pulse case runs to 100 and 400 steps of the seam leg (final times $1.56\cdot10^{-10}$ and
$6.25\cdot10^{-10}$ s, the times matched when $\chi > 1$ shortens the step), with and without hyperbolic cleaning of
$\mathbf B$ (`constrained_transport = B`, `divergence_correction = hyperbolic`: a tenth variable $\psi$, cleaning speed
$\chi c$, damping rate $\alpha$, [hazards](../hazards#divergence-correction)). The divergence history gives the maximum
over the domain at every step; checkpoints every 5 steps locate it. Two boundary layouts: the case as committed
(Neumann on all faces, so the seam plane $x = 0$ meets the four $y$ and $z$ faces), and the same case with periodic $y$
and $z$ faces, where the seam meets no face.

| Faces | Cleaning | Step 5 | $t = 1.56\cdot10^{-10}$ s (100 steps) | $t = 6.25\cdot10^{-10}$ s (400 steps) | max over $t \le 6.25\cdot10^{-10}$ s |
|---|---|---|---|---|---|
| Neumann (the case) | off | 12.65 | 1002 | 6985 | 6985 |
| Neumann | $\chi = 1.01$, $\alpha = 0$ | 7.716 | 173.0 | $7.77\cdot10^{4}$ | $7.77\cdot10^{4}$ |
| Neumann | $\chi = 1.01$, $\alpha = 6.4\cdot10^{10}$ | 8.214 | 191.6 | — | — |
| Neumann | $\chi = 1.01$, $\alpha = 6.4\cdot10^{11}$ ($\alpha\Delta t = 0.99$) | 10.35 | 553.7 | 3866 | 3866 |
| Neumann | $\chi = 2$, $\alpha = 0$ | 7.464 | 129.1 | — | — |
| Neumann | $\chi = 2$, $\alpha = 6.4\cdot10^{10}$ | 6.637 | 140.1 | $1.68\cdot10^{4}$ | $1.68\cdot10^{4}$ |
| Neumann | $\chi = 4$, $\alpha = 0$ | 7.027 | 1177 | — | — |
| periodic $y$, $z$ | off | 12.65 | 34.86 | 44.64 | 92.46 |
| periodic | $\chi = 1.01$, $\alpha = 0$ | 7.716 | 12.66 | 53.18 | 53.18 |
| periodic | $\chi = 1.01$, $\alpha = 6.4\cdot10^{10}$ | 8.214 | 3.626 | 4.465 | 9.943 |
| periodic | $\chi = 1.01$, $\alpha = 6.4\cdot10^{11}$ | 10.35 | 7.885 | 27.92 | 28.80 |
| periodic | $\chi = 2$, $\alpha = 0$ | 7.464 | 13.58 | $1.65\cdot10^{5}$ | $1.65\cdot10^{5}$ |
| periodic | $\chi = 2$, $\alpha = 6.4\cdot10^{10}$ | 6.637 | 1.202 | 1.983 | 6.766 (the first step) |
| periodic | $\chi = 4$, $\alpha = 0$ | 7.027 | 94.45 | — | — |

*Step 5 is the value at $t = 7.82\cdot10^{-12}$ s (step 10 for $\chi = 2$, step 20 for $\chi = 4$). CPU and FNL agree to
every printed digit (the Neumann rows on both backends; the periodic rows on CPU, and on FNL for off and $\chi = 1.01$,
$\alpha = 0$ to 400 steps and $\chi = 2$, $\alpha = 6.4\cdot10^{10}$ to $1.56\cdot10^{-10}$ s). Step cost of the
cleaning scalar, the case at 1 rank (500 steps against 5 on CPU, 2000 against 5 on FNL, twice each): CPU 0.3270 s
without, 0.3871 s with (+18 %); FNL 9.87 ms without, 11.43 ms with (+16 %). $\chi = 2$ halves the step, so per unit of
simulated time the cleaning costs 2.4 times (CPU) and 2.3 times (FNL) the run without it. Measured: ADAM `0bcdf870` +
the night-2026-11-22 tree (HEAD `95ffd728` with the working tree of 2026-11-23), CPU gfortran 16.2.0 + Open MPI 5.0.11
(`iac_gnu`) on medooza (gnode25), FNL nvfortran 26.9, CUDA 12.9 (`iac_nvf`) on NVIDIA A30 (gnode23), 2026-11-23 (jobs
1369621, 1369624, 1369625, 1369629 CPU; 1369620, 1369626 FNL; `check.sh --cleaning`, every leg: 1369635 CPU, 1369636 FNL).*

What the table shows:
- **The growth "6.9 to 1002 in 100 steps" is not the seam alone.** From step 70, when the pulse ($\sigma = 0.01$ m,
  cut at $5\sigma$, moving in $+y$) reaches the $y$ faces, the maximum sits where the seam plane meets the Neumann $y$
  and $z$ faces (at step 100: $x = -0.9$ mm, $y = 57$ mm, $z = 53$ mm). Away from the faces ($|y|, |z| < 0.04$ m) the
  same run has 89.9 at step 100. With periodic $y$ and $z$ faces the seam source stays between 28 and 92 from step 10
  to step 400: bounded at fixed $h$ on this horizon, with no growth trend.
- **Cleaning with damping bounds it.** With $\chi = 2$ and $\alpha = 6.4\cdot10^{10}$ s$^{-1}$ ($\alpha\Delta t = 0.05$)
  the maximum stays between 0.37 and 2.2 from the time of step 50 to step 400 (beside the seam, $|x| < 7.5$ mm: 0.10
  to 2.0): 22 to 29 times below the run without cleaning at 100 and 400 steps. $\chi = 1.01$ with the same $\alpha$
  does nearly as well (1.0 to 7.3 over the same interval, 3.6 at 100 steps).
- **Cleaning without damping does not.** It moves the divergence away from the seam at $\chi c$ instead of removing it:
  with $\chi = 1.01$ the maximum leaves the seam for the Neumann $x$ faces and grows to 53 by step 400 (40 beside the
  seam); with $\chi = 2$ and
  $\alpha = 0$ it runs away at the Neumann $x$ faces ($1.65\cdot10^{5}$ at 400 steps). Too strong a damping is no better:
  $\alpha\Delta t = 0.99$ gives 28 at 400 steps.
- **Where the seam meets a Neumann face, no setting helps.** On the case as committed every variant grows at the face
  edges, the damped $\chi = 2$ one to $1.68\cdot10^{4}$ in 400 steps, against 6985 without cleaning: the $\psi$
  boundary condition at those faces, not the seam, limits the cleaning there.


Related pages: [quadtree](./quadtree) runs the same seam on a quadtree; [inter-realm 2:1 seam](./inter-realm-refined)
reproduces the seam source across a realm boundary; [absorbing layers on AMR](./absorbing-layers-amr) and the
[`[IO] seam_divB_tol` monitor](../hazards#seam-div-b) are the guard-rails.

## Findings

- The jump is the only source: on the uniformly refined grid, with the same resolution, the same block count and the
  same ghost exchange, both divergences stay at round-off.
- $\nabla\cdot\mathbf D$ of the $+y$ pulse stays at round-off across the seam: it is $\partial D_z/\partial z$ with no
  $z$ variation, and every partition-of-unity ghost fill preserves it.
- The seam source converges under refinement ($p_{obs} = 1.15$ on the first rung, 1.62 on the second). At fixed $h$ it
  stays between 28 and 92 over 400 steps where the seam meets no face; the growth to $10^3$ in 100 steps of the case
  sits where the seam meets its Neumann $y$ and $z$ faces ([#145](https://github.com/szaghi/adam/issues/145)). It is a
  constraint violation, not an energy instability.
- Hyperbolic cleaning of $\mathbf B$ with damping ($\chi = 2$, $\alpha = 6.4\cdot10^{10}$ s$^{-1}$) keeps the seam
  source below 2.2, 22 to 29 times below the uncleaned run, for 18 % more per step and twice the steps; without damping it
  only moves the divergence, and with $\chi = 2$ it runs away at the Neumann faces.
- The $\nabla\cdot\mathbf J$ band is a factor 3, not 5 %: $\nabla\cdot\mathbf J$ of the coil filament cancels gradients
  of order $10^{6}$ down to $3\cdot10^{-5}$, and the cross-compiler noise of the coil geometry ($10^{-11}$ relative on
  $\mathbf J$) leaves a $10^{-5}$ absolute floor.
- A low-storage RK scheme on the seam case stops with the contract message; before #25 it crashed.
- The control leg of `rmf-amr-fd-pulse/check.sh` stopped at initialization since `[IO] unused_keys = error` became the
  default (#80 P5): `markers_number = 0` leaves `[amr_marker_1]` unread. The leg now sets `unused_keys = warn` for
  that run ([#145](https://github.com/szaghi/adam/issues/145)); every leg passes again on both backends.
