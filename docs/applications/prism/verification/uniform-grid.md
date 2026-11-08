# Uniform grid: energy and divergence

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `periodic-pulse-blocks` | — (base) | P | CPU, FNL | no | `periodic-pulse-blocks` (digest + residuals, CPU and FNL) |

## Case and scope

A Gaussian pulse travels once around a box periodic on every face. On a uniform grid the source-free Maxwell system
keeps two discrete invariants: the electromagnetic energy, up to the time integrator, and the divergences
$\nabla_h\cdot\mathbf D$ and $\nabla_h\cdot\mathbf B$, exactly. The centred curl and divergence are built from the same
antisymmetric first-derivative stencil, so $\nabla_h\cdot(\nabla_h\times\cdot) = 0$ holds to round-off.

What the case catches:
- **a dissipative or unstable space operator:** the energy must stay flat over one crossing;
- **a broken discrete identity:** the divergence histories must stay at round-off;
- **broken periodicity across blocks and ranks:** the 64 blocks are split 32/32 over two ranks, so the pulse crosses
  block faces owned by different ranks and the periodic faces themselves. The former same-block wrap fails the
  digest on all 192 residual rows.

Every other page builds on these identities: the AMR and inter-realm pages measure what a seam does to them.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell in $\mathbf D$, $\mathbf B$, SI units (`physical_model = EM`), no sources | [models](../models#models) |
| Space | centred finite differences of order 6 (`fd_centered`, `fdv_order = 6`) | `compute_curl_fd_centered`, `compute_divergence_fd_centered`, [`adam_fdv_operators_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_fdv_operators_library.F90); [numerics](../numerics#the-discrete-curl-and-divergence) |
| Time | SSP Runge–Kutta, 5 stages, 4th order (`runge-kutta-ssp-54`), CFL 0.5: $\Delta t = \Delta x/(4c)$ | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90); [numerics](../numerics#time-integration) |
| Periodic faces | the library periodicity (`BC_PERIODIC`): periodic tree neighbours and the ghost exchange across blocks and ranks | [boundary conditions](../boundary-conditions#face-kinds) |
| FNL backend | the same kernels on the device | [`adam_prism_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/fnl/adam_prism_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-0.1, 0.1]^3$ m |
| Grid | $48^3$ cells: 64 blocks of $12^3$ (`iu_ref_levels = 2`), $\Delta x = 0.2/48$, `ngc = 3` |
| Ranks | 2 (32 blocks each) |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+x$, `r0 = 0.07`, `sigma = 0.015`, `k0 = 0`, `B0 = 10` |
| Boundaries | `periodic` on all six faces |
| Steps | 192 = one crossing of the box |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/periodic-pulse-blocks/input.ini` | the case |
| `src/tests/prism/regression/periodic-pulse-blocks/golden/{cpu,fnl}/` | the digest and the residuals log `run.sh` compares |
| `src/tests/tutorials/compare_saves.py` | compares two saves of a run block by block (the return error) |

```bash
REGRESSION_CASES=periodic-pulse-blocks src/tests/prism/regression/run.sh cpu --no-build   # the digest
```

The case is also [tutorial 1](/tutorials/prism/pulse), which walks through the run and its outputs.

## Results

| Check | Bound | CPU | FNL |
|---|---|---|---|
| digest and residuals against the golden | `rtol = 1e-6`, `atol = 1e-3` | pass | pass |
| total energy $W_D + W_B$, relative drift over 192 steps | — | $1.2\cdot10^{-6}$ | ⟨run⟩ |
| $W_D$ against $W_B$ | equal (travelling wave) | equal to 15 digits | ⟨run⟩ |
| $\lvert\nabla\cdot\mathbf D\rvert$, $\lvert\nabla\cdot\mathbf B\rvert$, last row of the divergence history | round-off | 0, 0 | ⟨run⟩ |
| 64 blocks on 2 ranks against one $48^3$ block | — | $\approx 2\cdot10^{-15}$ | — |
| return after one crossing, $\max\lvert q(192) - q(0)\rvert / \max\lvert q(0)\rvert$, $B_z$ and $D_y$ | — | $3.961\cdot10^{-2}$ | — |

## Findings

- The energy drift, $1.2\cdot10^{-6}$ over one crossing, is the Runge–Kutta's: the centred scheme is non-dissipative.
- The divergences stay at zero: on a uniform grid the identity $\nabla_h\cdot\nabla_h\times = 0$ is exact.
- The 4 % return error is not the scheme's. The initial Gaussian is not wrapped across the periodic faces, so the field
  jumps by 18 % of the peak at $x = 0.1$; a 1-D model of the run gives 4.0 % with the cut and 0.14 % with a wrapped
  pulse ([tutorial 1](/tutorials/prism/pulse#what-to-look-at)).
- The block layout does not change the solution beyond round-off: the 64-block run matches the single-block run to
  $2\cdot10^{-15}$.
