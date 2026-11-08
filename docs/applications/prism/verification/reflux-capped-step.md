# FV reflux on a capped last step

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-amr-fv-pulse-timed` | [#38](https://github.com/szaghi/adam/issues/38) | P | CPU, FNL | no | — |

## Case and scope

On a time-driven run (`it_max <= 0`) PRISM caps its last step so that it lands on `time_max`. The Berger–Colella
reflux correction must be scaled by that capped step, the one the update used. Before `4c1dd1a0` it was scaled by the
forest's uncapped $\Delta t$, and the coarse–fine faces were not conservative on the last step.

With $\mathbf J = 0$ the finite-volume update with reflux keeps the domain integrals of $\mathbf D$ and $\mathbf B$
constant, up to the fluxes through the domain boundary. The case measures that conservation step by step.

What the case catches: a reflux correction scaled by any step other than the one the update used, as a jump of
$\int D_z$ on the capped step. $\int B_x$ cannot see the defect: $B_x$ has no $x$-flux, so the reflux at the
$x$-normal seams never touches it. $\int D_z$ (its $x$-flux is $H_y$) can.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, no source | [models](../models#models) |
| Space | centred finite volumes of order 6 (`fv_centered`) | [numerics](../numerics#centred-finite-volumes-fv-centered) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5, time-driven (`it_max = -1`): $\Delta t = \mathrm{CFL}\,\Delta x_{\min}/(2c) = 1.5635\cdot10^{-12}$ s, `time_max` $= 4.5\,\Delta t$ | [numerics](../numerics#time-integration) |
| Reflux | correction scaled by the step of the update | `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); `apply_reflux_to_stage_forest`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $x, z \in [-0.06, 0.06]$, $y \in [-0.12, 0.12]$ m |
| Grid | 64 blocks of $8\times16\times8$ at level 2 (4 per axis, cubic cells, $\Delta x = 3.75\cdot10^{-3}$); the band $x \in [-0.03, 0.03]$ refined to level 3: two 2:1 planes, both $x$ boundaries coarse; `ngc = 3` |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, transverse to the seams, `sigma = 0.01`, `B0 = 10` |
| Boundaries | `Neumann` on all six faces |
| Steps | 4 full steps and one of half a step (`time_max = 7.0361171362237545e-12` s); a checkpoint every step |
| Ranks | 2 (`--np N`) |

The $y$ extent is twice that of `rmf-amr-fd-pulse`, 12σ on each side: with the narrower box the Gaussian tails leak
through the $y$ faces and even the grid without seams drifts by $10^{-8}$ per step.

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-amr-fv-pulse-timed/input.ini` | the case |
| `src/tests/prism/regression/rmf-amr-fv-pulse-timed/check.sh` | runs it and integrates $\mathbf D$ and $\mathbf B$ over every checkpoint (`--np N`; `PRISM_EXE`) |

```bash
cd src/tests/prism/regression/rmf-amr-fv-pulse-timed && ./check.sh
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

| Leg | Asserts | Bound | CPU | FNL |
|---|---|---|---|---|
| capped step | the last time step shorter than 0.99 of the first | — | pass | pass |
| seams active | the `forest: reflux max\|F_coarse-F_fine_sum\|` line is printed | — | pass | pass |
| conservation | at every checkpoint, the drift of each $D$ integral relative to the largest $D$ integral, and of each $B$ integral relative to the largest $B$ integral | $\le 10^{-9}$ | $2.2\cdot10^{-11}$ | $2.2\cdot10^{-11}$ |

$\int D_z$ drift per step (CPU, 2 ranks):

| Step | 1 | 2 | 3 | 4 | 5 (capped) |
|---|---|---|---|---|---|
| fixed run | $1.8\cdot10^{-15}$ | $1.4\cdot10^{-13}$ | $1.8\cdot10^{-12}$ | $1.0\cdot10^{-11}$ | $2.2\cdot10^{-11}$ |
| fix reverted (negative control) | the same | the same | the same | the same | $1.4\cdot10^{-7}$ |

$\int B_x$ stays at $10^{-16}$.

Related pages: [static AMR and the FV reflux register](./amr-reflux) for the register itself; FLUME's
[conservation across AMR seams](/applications/flume/verification/euler/amr-conservation) measures the same property on
the Euler equations.

## Findings

- With the correction scaled by the capped step, the drift grows smoothly to $2.2\cdot10^{-11}$: the seams'
  disturbance reaching the $x$ boundaries, not a reflux defect.
- With the fix reverted the drift jumps to $1.4\cdot10^{-7}$ on the capped step, 6500 times larger. The bound,
  $10^{-9}$, sits 45 times above the fixed run and 140 times below the reverted one.
