# EGLM reflux leaves B alone

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| G3a, G3b | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) · G3a, G3b; [#78](https://github.com/szaghi/adam/issues/78) (follow-up) | P | CPU, FNL (G3a table: FNL) | yes | `blast-amr-limiter` |

## Case and scope

On a static 2:1 grid, EGLM's $\langle|B_z|\rangle$ on the advected field loop rises to 70 times its uniform-grid level
([#78](https://github.com/szaghi/adam/issues/78), reproduced by [AD](./divergence-across-seams)). The divergence the
cleaning acts on, D_flux, is the same under GLM and EGLM, so the excess is not a divergence difference. Two candidate
causes remain: the second-order restricted ghosts that EGLM's centred source divergence reads at the seam, and the
Berger–Colella reflux.

G3a separates them by running the field loop with sixth-order ghosts ([C1, D1](./order-across-seams)) and with the
reflux on, off and under `reflux_policy = troubled` ([E2](./reflux-policy)). G3b fixes the cause and checks that the fix
keeps the regressions, the MV-5 orders and the multirealm suite.

What the case catches:
- **a reflux that breaks EGLM's source/flux cancellation:** $\langle|B_z|\rangle$ on the refined grid must stay within a
  small factor of the uniform level;
- **a fix that changes GLM or the conserved variables:** GLM must be unchanged, and the regression goldens must hold.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$, GLM or EGLM (Derigs sources), `glm_ch` 4.0, `glm_alpha` 0.18, `glm_damping_length` 1.0 | [GLM](../../models#divergence-control-mixed-glm), [EGLM](../../models#divergence-control-eglm); `add_eglm_sources`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90); `add_eglm_sources_dev`, [`adam_flume_fnl_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_mhd_eglm_kernels.F90) |
| Space, time | WENO5 flux splitting on characteristic variables; SSP Runge–Kutta 5 stages, CFL 0.4 | as in [AD](./divergence-across-seams) |
| Seam operators | default (`mean` + `tricubic`) or `lagrange6` + `quintic` | [C1, D1](./order-across-seams) |
| Reflux (G3a) | on (`always`), off, or `troubled` (flags no face of this smooth loop, so it equals the reflux off bit for bit) | `apply_reflux_to_stage_forest`, `flag_troubled_seam_cells`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Reflux under EGLM (G3b, the fix) | the flux mismatch of $\mathbf{B}$ is zeroed: the EGLM reflux corrects $\rho$, $\rho\mathbf u$, $E$ and $\psi$ only; GLM unchanged (refluxes $\mathbf{B}$); no input key | `apply_reflux_to_stage_forest`, `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Case | the A2 field loop of [AD](./divergence-across-seams): $[-1, 1] \times [-0.5, 0.5]$, quadtree, the box $[0, 1] \times [-0.5, 0.5]$ refined one 2:1 level |
| Initial condition | [field loop](../../initial-conditions#field-loop-advected-field-loop), radius 0.3, $A_0 = 10^{-3}$; $(\rho, u, v, w, p) = (1, 2, 1, 1, 1)$ |
| Grid | refined $N = 64$ against uniform $N = 128$ |
| Boundaries | periodic in $x$, $y$; `extrapolation` in $z$ |
| Final time | $t = 1$ |
| Ranks | 2 |
| MV-5 | the linear-wave orders of [order of accuracy](../mhd/order-of-accuracy) under EGLM |
| Regressions | the FLUME regression suite on CPU and FNL; the one case with EGLM across 2:1 seams is `blast-amr-limiter` |
| Multirealm | the [forest seams](../euler/multi-realm) suite |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-divb/make_amr_divb.py loop` | the field loop inputs (`--cleaning glm` or `eglm`, `--seam-ghost-fill`, `--seam-restriction`, `--reflux-policy`, `--no-reflux`) |
| `src/tests/flume/verification/amr-divb/amr_divb_oracle.py loop` | the $\langle\lvert B_z\rvert\rangle/A_0$ report |
| `src/tests/flume/verification/amr-divb/check.sh --legs loop` | the field loop, GLM and EGLM, refined against uniform |
| `src/tests/flume/regression/blast-amr-limiter/input.ini` | the regression case with EGLM across 2:1 seams (golden digest, rtol $10^{-6}$) |
| `src/tests/flume/regression/run.sh` | the regression suite |
| `src/tests/flume/verification/make_doc_runs.sh <runs> divb` (FNL), `make_doc_figures.py --only field-loop-seams` | re-runs the six field loops with a checkpoint every 20 steps |

```bash
cd src/tests/flume/verification/amr-divb && ./check.sh --legs loop
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs loop
```

## Results

![The field loop across the 2:1 box under GLM and EGLM](/flume/field-loop-seams.png)

*FNL, with the G3b reflux: the A2 field loop under GLM and EGLM, refined $N = 64$ with and without reflux, and uniform $N = 128$, a checkpoint every 20 steps. Left: the magnetic energy; right: $\langle|B_z|\rangle/A_0$. The re-runs end at the values of the tables (EGLM $9.11\cdot10^{-6}$ with the reflux, $9.31\cdot10^{-6}$ without, $2.72\cdot10^{-6}$ uniform; GLM $1.10\cdot10^{-3}$, $1.18\cdot10^{-3}$, $6.86\cdot10^{-4}$). The reflux of B (70 times the uniform level) is no longer in the code, so it is not drawn.*

**G3a, the cause.** $\langle\lvert B_z\rvert\rangle/A_0$ at $t = 1$, with B refluxed as before G3b (FNL, 2 ranks):

| Seam operators, reflux | EGLM (× uniform) | GLM (× uniform) |
|---|---|---|
| uniform N = 128 | $2.72\cdot10^{-6}$ (1) | $6.86\cdot10^{-4}$ (1) |
| default, `always` | $1.91\cdot10^{-4}$ (70) | $1.10\cdot10^{-3}$ (1.6) |
| default, off (= `troubled`) | $9.31\cdot10^{-6}$ (**3.4**) | $1.18\cdot10^{-3}$ (1.7) |
| `lagrange6` + `quintic`, `always` | $2.89\cdot10^{-4}$ (106) | $1.12\cdot10^{-3}$ (1.6) |
| `lagrange6` + `quintic`, off (= `troubled`) | $3.14\cdot10^{-5}$ (11.5) | $1.56\cdot10^{-3}$ (2.3, run with `troubled`) |

**G3b, the fix.** EGLM, default seam operators, refined $N = 64$, by the components the reflux corrects:

| EGLM, reflux | $\langle\lvert B_z\rvert\rangle/A_0$ | × uniform |
|---|---|---|
| all components | $1.91\cdot10^{-4}$ | 70 |
| all but B (G3b) | $9.11\cdot10^{-6}$ | **3.3** |
| all but B and $\psi$ | $9.31\cdot10^{-6}$ | 3.4 |
| none | $9.31\cdot10^{-6}$ | 3.4 |

FNL reproduces the table: EGLM $9.11\cdot10^{-6}$, GLM $1.10\cdot10^{-3}$ (unchanged). MV-5 under EGLM keeps its
orders (4.96 to 4.97 on every wave), and the multirealm suite passes.

**Regressions.** The FLUME regression suite passes on CPU and FNL. `blast-amr-limiter` moves inside its digest
tolerance (rtol $10^{-6}$), so its golden is kept:
- step 0 is identical, and $B_z$ (zero in 2-D) does not move;
- at step 100 seven of 26 digest rows differ (`r`, `ru`, `rv`, `rE`, `bx`, `by`, `psi`), the largest by
  $2.4\cdot10^{-7}$ relative on CPU and $9.8\cdot10^{-8}$ on FNL;
- the integrals drift as in the golden: $\rho$ to round-off ($3\cdot10^{-16}$), the rest at the golden's levels
  ($\rho\mathbf u$ $2.5\cdot10^{-4}$, $E$ $1.9\cdot10^{-2}$, B $4.5\cdot10^{-4}$, $\psi$ $1.3\cdot10^{-4}$), set by the
  blast's boundaries and EGLM's non-conservative sources, not by the reflux.

Related cases: [AD](./divergence-across-seams) is the A2 case; [E2](./reflux-policy) measures the conservation cost of
`troubled`; the uniform-grid EGLM is [EGLM cleaning](../mhd/eglm); the loop itself is [field loop](../mhd/field-loop).

## Findings

- **The reflux causes #78, not the ghosts.** With the default operators and no reflux, EGLM's $B_z$ is 3.4 times its
  uniform level instead of 70; sixth-order ghosts do not help (106 with the reflux, 11.5 without). GLM is insensitive
  (1.6 to 2.3 times in every row).
- **Mechanism (inferred from G3a and G3b, not measured directly).** EGLM's Derigs source
  $-\nabla\cdot\mathbf B\,(0, \mathbf B, \mathbf u\cdot\mathbf B, \mathbf u, \ldots)$ cancels, on a uniform grid, the
  $B_z$ that the flux divergence generates. The reflux corrects a coarse skin cell by the flux mismatch
  $F_c - \sum F_f$ after the step, but the source of that step saw the coarse fluxes' divergence: the correction changes
  the discrete divergence without the matching source term, and the cancellation fails at the refluxed cells.
- **Refluxing B buys a conservation EGLM does not have.** Under EGLM the induction equation carries the source
  $-\mathbf u\,\nabla\cdot\mathbf B$, so $\int \mathbf B$ is not a conserved integral. The EGLM reflux therefore corrects
  $\rho$, $\rho\mathbf u$, $E$ and $\psi$ and leaves $\mathbf{B}$ alone, on CPU and FNL; GLM still refluxes
  $\mathbf{B}$. #78 drops from 70 to 3.3 times the uniform level.
- **The cancellation needs the centred derivative of the cell values.** Taking the source divergence from the $\psi$
  face fluxes instead, $\sum_d (F_\psi(i+\tfrac12) - F_\psi(i-\tfrac12)) / (c_h \Delta x_d)$, keeps the MV-5 orders but
  raises EGLM's uniform $\langle|B_z|\rangle$ 47 times ($2.72\cdot10^{-6} \to 1.28\cdot10^{-4}$).
