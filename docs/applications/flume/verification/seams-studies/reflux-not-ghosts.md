# #78 is the reflux, not the restricted ghosts

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| G3a | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | FNL only | no ([#85](https://github.com/szaghi/adam/issues/85)) | — |

## Case and scope

[A2](../seams/divergence-across-seams) traced EGLM's 70-fold $\langle|B_z|\rangle$ on a static 2:1 grid
([#78](https://github.com/szaghi/adam/issues/78)) to the centred divergence of its Derigs source reading the
second-order restricted ghosts. Phases C to E gave the seams sixth-order ghosts ([C1](../seams/order-across-seams),
[D1](./quintic-fill)) and a reflux policy ([E2](./reflux-policy)), so G3a asks first whether #78 survives them, and so
which of the two, the ghosts or the reflux, causes it.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$, GLM or EGLM, `glm_ch` 4.0, `glm_alpha` 0.18, `glm_damping_length` 1.0 | [GLM](../../models#divergence-control-mixed-glm), [EGLM](../../models#divergence-control-eglm); `add_eglm_sources`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90); `add_eglm_sources_dev`, [`adam_flume_fnl_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_mhd_eglm_kernels.F90) |
| Space, time | WENO5 flux splitting on characteristic variables; SSP Runge–Kutta 5 stages, CFL 0.4 | as in [AD](../seams/divergence-across-seams) |
| Seam operators | default (`mean` + `tricubic`) or `lagrange6` + `quintic` | [C1](../seams/order-across-seams), [D1](./quintic-fill) |
| Reflux | on (`always`), `troubled` (flags no face of this smooth loop, so it equals the reflux off bit for bit), or off | `apply_reflux_to_stage_forest`, `flag_troubled_seam_cells`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Case | the A2 field loop of [AD](../seams/divergence-across-seams): $[-1, 1] \times [-0.5, 0.5]$, quadtree, the box $[0, 1] \times [-0.5, 0.5]$ refined one 2:1 level |
| Initial condition | [field loop](../../initial-conditions#field-loop-advected-field-loop), radius 0.3, $A_0 = 10^{-3}$; $(\rho, u, v, w, p) = (1, 2, 1, 1, 1)$ |
| Grid | refined $N = 64$ against uniform $N = 128$ |
| Boundaries | periodic in $x$, $y$; `extrapolation` in $z$ |
| Final time | $t = 1$ |
| Ranks | 2, FNL |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-divb/make_amr_divb.py loop` | the field loop inputs (`--cleaning glm` or `eglm`, `--seam-ghost-fill`, `--seam-restriction`, `--reflux-policy`, `--no-reflux`) |
| `src/tests/flume/verification/amr-divb/amr_divb_oracle.py loop` | the $\langle\lvert B_z\rvert\rangle/A_0$ report |

## Results

No figure yet ([#85](https://github.com/szaghi/adam/issues/85): the field-loop magnetic energy against time, GLM and
EGLM, reflux on and off).

| Seam operators, reflux | EGLM $\langle\lvert B_z\rvert\rangle/A_0$ (× uniform) | GLM (× uniform) |
|---|---|---|
| uniform N = 128 | $2.72\cdot10^{-6}$ (1) | $6.86\cdot10^{-4}$ (1) |
| default, `always` | $1.91\cdot10^{-4}$ (70) | $1.10\cdot10^{-3}$ (1.6) |
| default, off (= `troubled`) | $9.31\cdot10^{-6}$ (**3.4**) | $1.18\cdot10^{-3}$ (1.7) |
| `lagrange6` + `quintic`, `always` | $2.89\cdot10^{-4}$ (106) | $1.12\cdot10^{-3}$ (1.6) |
| `lagrange6` + `quintic`, off (= `troubled`) | $3.14\cdot10^{-5}$ (11.5) | $1.56\cdot10^{-3}$ (2.3, run with `troubled`) |

Related cases: [AD](../seams/divergence-across-seams) is the A2 attribution this corrects; [G3b](./eglm-reflux) is the
fix; the uniform-grid EGLM is [EGLM cleaning](../mhd/eglm); the loop itself is [field loop](../mhd/field-loop).

## Findings

- **The reflux causes #78.** With the default operators and no reflux, EGLM's $B_z$ is 3.4 times its uniform level
  instead of 70; higher-order ghosts do not help (106 with the reflux, 11.5 without). GLM does not care (1.6 to 2.3
  times in every row). The A2 attribution to the restricted ghosts is wrong.
- **Mechanism (inferred, not measured directly):** EGLM's Derigs source
  $-\nabla\cdot\mathbf B\,(0, \mathbf B, \mathbf u\cdot\mathbf B, \mathbf u, \ldots)$ cancels, on a uniform grid, the
  $B_z$ that the flux divergence generates. The reflux corrects a coarse skin cell's state by the flux mismatch
  $F_c - \sum F_f$ after the step, but the source that step used saw the coarse fluxes' divergence: the correction
  changes the discrete divergence the fluxes imply without the matching source term, so the cancellation fails exactly
  at the refluxed cells.
- **Workaround before G3b:** EGLM with `reflux_policy = troubled` keeps #78 at 3.4 times on smooth fields (the reflux
  is applied only where a shock or the limiter flags a face), at the conservation price measured in
  [E2](./reflux-policy).
- **The fix:** [G3b](./eglm-reflux): the EGLM reflux leaves $\mathbf{B}$ alone (3.3 times).
