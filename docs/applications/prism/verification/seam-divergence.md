# 2:1 seam div(B): truncation order

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-amr-fd`, `rmf-amr-fd-pulse` | [#21](https://github.com/szaghi/adam/issues/21), [#22](https://github.com/szaghi/adam/issues/22), [#25](https://github.com/szaghi/adam/issues/25), [#26](https://github.com/szaghi/adam/issues/26), [#29](https://github.com/szaghi/adam/issues/29) ([#99](https://github.com/szaghi/adam/issues/99) open) | P | CPU, FNL | no | `rmf-amr-fd`, `rmf-amr-fd-pulse` (digest + residuals, CPU and FNL) |

## Case and scope

On a uniform grid $\nabla_h\cdot\mathbf B$ stays at round-off ([uniform grid](./uniform-grid)). At a 2:1 face the
coarse and fine sides difference with stencils of different resolution, the curl and the divergence no longer commute,
and the seam injects a $\nabla\cdot\mathbf B$ source. The source is truncation-order, $O(h^p)$: it converges under
refinement and it grows in time at fixed $h$. Five local fixes were ruled out and the resolution is accept-truncation
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
| Time | SSP Runge–Kutta (5,4), CFL 0.5; SSP (1,1) and `runge-kutta-1` in the RK legs | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
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
| `src/tests/prism/regression/rmf-amr-fd-pulse/input.ini`, `check.sh` | the pulse: control, uniform, seam, RK legs, quadtree leg ([quadtree](./quadtree)), `--convergence` |

```bash
cd src/tests/prism/regression/rmf-amr-fd && ./check.sh
cd ../rmf-amr-fd-pulse && ./check.sh --convergence
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

Long-horizon measurement, not a leg: on the pulse at fixed $h$ the seam max $\lvert\nabla\cdot\mathbf B\rvert$ grows
from 6.9 to 1002 in 100 steps while the control stays at 0 and $\lVert\mathbf B\rVert$ changes by −3 %. The
second-rung ratio of the refinement ladder, at `ni = 64`, is 1.62, toward the predicted 2.

Related pages: [quadtree](./quadtree) runs the same seam on a quadtree; [inter-realm 2:1 seam](./inter-realm-refined)
reproduces the seam source across a realm boundary; [absorbing layers on AMR](./absorbing-layers-amr) and the
[`[IO] seam_divB_tol` monitor](../hazards#seam-div-b) are the guard-rails.

## Findings

- The jump is the only source: on the uniformly refined grid, with the same resolution, the same block count and the
  same ghost exchange, both divergences stay at round-off.
- $\nabla\cdot\mathbf D$ of the $+y$ pulse stays at round-off across the seam: it is $\partial D_z/\partial z$ with no
  $z$ variation, and every partition-of-unity ghost fill preserves it.
- The seam source converges under refinement ($p_{obs} = 1.15$ on the first rung, 1.62 on the second) and grows without
  bound in time at fixed $h$. It is a constraint violation, not an energy instability.
- The $\nabla\cdot\mathbf J$ band is a factor 3, not 5 %: $\nabla\cdot\mathbf J$ of the coil filament cancels gradients
  of order $10^{6}$ down to $3\cdot10^{-5}$, and the cross-compiler noise of the coil geometry ($10^{-11}$ relative on
  $\mathbf J$) leaves a $10^{-5}$ absolute floor.
- A low-storage RK scheme on the seam case stops with the contract message; before #25 it crashed.
