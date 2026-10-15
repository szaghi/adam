# EGLM reflux leaves B alone

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| G3b | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | CPU, FNL | no ([#85](https://github.com/szaghi/adam/issues/85)) | `blast-amr-limiter` |

## Case and scope

[G3a](./reflux-not-ghosts) showed that [#78](https://github.com/szaghi/adam/issues/78), EGLM's
$\langle|B_z|\rangle$ 70 times its uniform level on a static 2:1 grid, is caused by the reflux. G3b looked for the fix,
and checked that it does not break the regressions, the MV-5 orders and the multirealm suite.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD with EGLM (Derigs sources), $\gamma = 5/3$; GLM for comparison | [EGLM](../../models#divergence-control-eglm); `add_eglm_sources`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90) |
| Reflux under EGLM (the fix) | the flux mismatch of $\mathbf{B}$ is zeroed: the EGLM reflux corrects $\rho$, $\rho\mathbf u$, $E$ and $\psi$ only; GLM unchanged (refluxes $\mathbf{B}$); no input key (owner decision) | `apply_reflux_to_stage_forest`, `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| First attempt (not kept) | the Derigs source's $\nabla\cdot\mathbf B$ from the $\psi$ face fluxes, $\sum_d (F_\psi(i+\tfrac12) - F_\psi(i-\tfrac12)) / (c_h \Delta x_d)$, the divergence of the interface field, with the matching source correction at the reflux | — |
| Space, time | WENO5 flux splitting on characteristic variables; SSP Runge–Kutta 5 stages, CFL 0.4 (field loop) | as in [AD](../seams/divergence-across-seams) |
| Seam operators | default (`mean` + `tricubic`), reflux on | [numerics](../../numerics#adaptive-mesh-refinement) |

## Initial and boundary conditions

| Run | Setup |
|---|---|
| Field loop | the A2 field loop of [AD](../seams/divergence-across-seams): refined $N = 64$ (box $[0, 1] \times [-0.5, 0.5]$) against uniform $N = 128$, $t = 1$, 2 ranks, default seam operators, reflux on; CPU and FNL |
| MV-5 | the linear-wave orders of [order of accuracy](../mhd/order-of-accuracy) under EGLM |
| Regressions | the FLUME regression suite on CPU and FNL; the one case with EGLM across 2:1 seams is `blast-amr-limiter` |
| Multirealm | the [multi-realm](../euler/multi-realm) suite |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-divb/check.sh --leg loop` | the field loop, GLM and EGLM, refined against uniform |
| `src/tests/flume/regression/blast-amr-limiter/input.ini` | the regression case with EGLM across 2:1 seams (golden digest, rtol $10^{-6}$) |
| `src/tests/flume/regression/run.sh` | the regression suite |

```bash
cd src/tests/flume/verification/amr-divb && ./check.sh --leg loop
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --leg loop
```

## Results

No figure yet ([#85](https://github.com/szaghi/adam/issues/85): the field-loop magnetic energy against time, GLM and
EGLM, reflux on and off).

**The first attempt, the source divergence from the face fluxes,** keeps the MV-5 orders (4.96 to 4.97 on every wave)
but raises EGLM's uniform $\langle|B_z|\rangle$ from $2.72\cdot10^{-6}$ to $1.28\cdot10^{-4}$ (47 times; refined
$N = 64$: $2.64\cdot10^{-4}$). The code was not kept.

**The fix.** The field loop (refined $N = 64$, default seam operators, reflux on) with the reflux of $\mathbf{B}$
removed, the rest refluxed:

| EGLM, reflux | $\langle\lvert B_z\rvert\rangle/A_0$ | × uniform |
|---|---|---|
| all components (before) | $1.91\cdot10^{-4}$ | 70 |
| all but B (G3b) | $9.11\cdot10^{-6}$ | **3.3** |
| all but B and $\psi$ | $9.31\cdot10^{-6}$ | 3.4 |
| none | $9.31\cdot10^{-6}$ | 3.4 |

The field loop on FNL confirms the table: EGLM $9.11\cdot10^{-6}$, GLM $1.10\cdot10^{-3}$ (unchanged). MV-5 under EGLM
keeps its orders (4.96 to 4.97 on every wave), and the multirealm suite passes.

**Gate.** The FLUME regressions pass unchanged on CPU and FNL (35 of 35 each at the time). `blast-amr-limiter` moves
but stays inside its digest tolerance (rtol $10^{-6}$), so its golden is kept:
- step 0 is identical, and $B_z$ (zero in 2D) does not move;
- at step 100 seven of 26 digest rows differ (`r`, `ru`, `rv`, `rE`, `bx`, `by`, `psi`), the largest by
  $2.4\cdot10^{-7}$ relative on CPU and $9.8\cdot10^{-8}$ on FNL;
- the integrals drift exactly as in the golden: $\rho$ to round-off ($3\cdot10^{-16}$), the rest at the golden's levels
  ($\rho\mathbf u$ $2.5\cdot10^{-4}$, $E$ $1.9\cdot10^{-2}$, B $4.5\cdot10^{-4}$, $\psi$ $1.3\cdot10^{-4}$, all set by
  the blast's boundaries and EGLM's non-conservative sources, not by the reflux).

Related cases: [G3a](./reflux-not-ghosts) located the cause; [AD](../seams/divergence-across-seams) is the A2 case;
the uniform-grid EGLM is [EGLM cleaning](../mhd/eglm).

## Findings

- **EGLM's cancellation needs the centred derivative.** EGLM's cancellation of the $B_z$ the flux splitting generates
  needs the centred derivative of the cell values, not the face divergence: the source divergence from the face fluxes
  makes the uniform $B_z$ 47 times worse.
- **Refluxing B is the whole of #78.** Under EGLM the induction equation carries the Derigs source
  $-\mathbf u\,\nabla\cdot\mathbf B$, so $\int \mathbf B$ is not a conserved integral: refluxing $\mathbf{B}$ buys a
  conservation the model does not have and breaks the source/flux cancellation it is built on.
- **Since G3b** (owner decision, no key) the EGLM reflux corrects $\rho$, $\rho\mathbf u$, $E$ and $\psi$ and leaves
  $\mathbf{B}$ alone, on CPU and FNL; GLM is unchanged (it refluxes $\mathbf{B}$, and its $\langle|B_z|\rangle$ does not
  depend on it, [G3a](./reflux-not-ghosts)). #78 drops from 70 to 3.3 times the uniform level.
