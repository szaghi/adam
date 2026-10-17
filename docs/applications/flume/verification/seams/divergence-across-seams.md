# Divergence across seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| AD (A2) | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)), [#78](https://github.com/szaghi/adam/issues/78) | FNL (tables); CPU for the loop leg and the D_flux cross-check | yes | — |

## Case and scope

Every div(B) column of the MHD checks is the centred divergence of the cell values, which is not the divergence GLM
and EGLM act on (#68 F4). `[diagnostics] divb_flux` adds **D_flux**, the divergence of the interface normal field read
from the scheme's own $\psi$ face fluxes, the coarse 2:1 faces taking the fine fluxes (the reflux), so that it sums
exactly (Gauss) over a composite grid. AD measures D_flux where the seams matter:
- `cpaw`: a source-free, exact nonlinear solution with $\nabla\cdot\mathbf{B} = 0$ that crosses two seams all the time:
  is the seam D_flux bounded or does it grow, and does it converge?
- `loop`: the field loop of [MV-9](../mhd/field-loop) under GLM and EGLM, the case of
  [#78](https://github.com/szaghi/adam/issues/78) (EGLM's $\langle|B_z|\rangle$ 70 times its uniform level on a static
  2:1 grid): is #78 a divergence difference?

What the case catches:
- **a secular growth of the seam divergence** under GLM (PRISM's centred GLM lets it run away,
  [#29](https://github.com/szaghi/adam/issues/29));
- **a seam divergence that does not converge** under refinement;
- **a broken D_flux instrument**: its signed integral over a periodic composite domain must be round-off.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$; mixed GLM (`cpaw`, `loop`) or EGLM (`loop`); `glm_alpha` 0.18, `glm_damping_length` 1.0; `glm_ch` 2.0 (`cpaw`, above the fast speed, about 1.08) or 4.0 (`loop`) | [GLM](../../models#divergence-control-mixed-glm), [EGLM](../../models#divergence-control-eglm); `add_eglm_sources`, [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_glm_face_split_fluxes`, `mhd_eglm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Seam operators | `mean` restriction, `tricubic` fill, reflux on (`always`) | [numerics](../../numerics#adaptive-mesh-refinement) |
| D_flux | $\sum_d (F_\psi(i+\tfrac12) - F_\psi(i-\tfrac12))/\Delta x_d$ over $c_h^2$ (GLM) or $c_h$ (EGLM), the coarse 2:1 faces corrected to the fine fluxes; the seam column over the cells within `ngc` of a seam; key `[diagnostics] divb_flux` | `compute_dflux_norms`, `compute_divb_history`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) and [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); [`adam_flume_diagnostics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_diagnostics_object.F90) |
| Boundaries | periodic in $x$ and $y$; `extrapolation` on the null $z$ faces | [boundary conditions](../../boundary-conditions#periodic) |

## Initial and boundary conditions

| Item | `cpaw` | `loop` |
|---|---|---|
| Domain | $[0, \sqrt5] \times [0, \sqrt5/2]$, null $z$ | $[-1, 1] \times [-0.5, 0.5]$, null $z$ |
| Initial condition | [MV-6 CPAW](../../initial-conditions#mhd-cpaw-circularly-polarised-alfven-wave), right polarised, inclined at $\tan\alpha = 2$ (63.43°), $B_\parallel = 1$, amplitude 0.1, wavelength 1, $\rho = 1$, $p = 0.1$ | [MV-9 field loop](../../initial-conditions#field-loop-advected-field-loop), radius 0.3, $A_0 = 10^{-3}$, centre $(0, 0)$; $(\rho, u, v, w, p) = (1, 2, 1, 1, 1)$ |
| Grid | quadtree (`ratio = 4`), $4 \times 4$ base blocks of $N/4 \times N/8 \times 1$ cells, `ngc = 3` | the same |
| Refined region | the two middle block columns, $x \in [0.6, 1.6]$: two $x$ seams the wave crosses all the time | the box $[0, 1] \times [-0.5, 0.5]$ |
| Ladder | refined and uniform $N$ = 64 / 128 (`--long` adds 256) | refined $N = 64$ against uniform $N = 128$ |
| Boundaries | periodic in $x$, $y$; `extrapolation` in $z$ | the same |
| Final time | 10 periods ($t = 10$) | $t = 1$ |
| History cadence | every 10 steps | every step |
| Ranks | 2 | 2 |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-divb/make_amr_divb.py` | writes the inputs from `mhd/cpaw/make_cpaw.py` and `mhd/field-loop/make_field_loop.py`; `--uniform`, `--periods`, `--cleaning`, `--time-max`, `--history` |
| `src/tests/flume/verification/mhd/amr_box.py` | `refine_box`: the refined region, one 2:1 level at initialisation |
| `src/tests/flume/verification/amr-divb/amr_divb_oracle.py` | `growth` (the `cpaw` floors, growth and orders) and `loop` (a report) |
| `src/tests/flume/verification/amr-divb/check.sh` | runs the legs (`--legs cpaw,loop`, `--long`, `--np N`, `--keep`) |
| `src/tests/flume/verification/make_doc_runs.sh <runs> divb` (FNL), `make_doc_figures.py --only divergence-across-seams` | re-runs the field loops (GLM, EGLM; refined, without reflux, uniform) with the history every step; the CPAW panel reads the histories of `check.sh` |

```bash
cd src/tests/flume/verification/amr-divb
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh          # FNL: about 12 min
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --long   # adds N = 256 (30 more minutes)
./check.sh --legs loop                                               # the loop leg on the CPU
```

The `cpaw` leg costs about 2 h on the CPU ($N = 128$ dominates): run it on FNL.

**Oracle** (`amr_divb_oracle.py`, reading the `*-divb_history.dat` columns). `cpaw`: per run the seam D_flux **floor**
(the median over the second half of the run, the step-0 row excluded), its **growth** (the mean over the last quarter
over the mean over the second quarter: about 1 when bounded) and the floor's order over the refined ladder; a uniform
run has no seam, its floor is that of the global D_flux. Bounds: growth at most 1.5 on every refined run, floor order
at least 1 on the finest pair. `loop`: the final seam and global D_flux and $\langle|B_z|\rangle/A_0$, a report.

## Results

![AD: the seam D_flux against time](/flume/divergence-across-seams.png)

*FNL, 2 ranks. Left: the seam maximum of D_flux on the source-free CPAW over 10 periods, refined $N = 64$ and 128 (the histories `check.sh` keeps), against the global maximum of the uniform $N = 128$ run: the seam floor is bounded and falls 10 times from 64 to 128. Right: the A2 field loop under GLM and EGLM, the refined seam maximum against the uniform global one. No divergence map: FLUME writes no divergence field, and a centred difference of the checkpoints would read unfilled ghosts exactly at the seams.*

The instrument: the signed integral of D_flux over the periodic composite domain is round-off, $\le 10^{-17}$ every
step on the refined MHD vortex, against $4\cdot10^{-7}$ without the fine fluxes; CPU and FNL agree to $10^{-11}$
relative.

`cpaw` (FNL, 2 ranks):

| N | Seam D_flux floor | Growth | Uniform floor | Seam / uniform |
|---|---|---|---|---|
| 64 | $5.39\cdot10^{-4}$ | 0.95 | $9.97\cdot10^{-6}$ | 54 |
| 128 | $5.43\cdot10^{-5}$ | 0.72 | $4.57\cdot10^{-7}$ | 119 |
| 256 | $2.12\cdot10^{-5}$ | 0.80 | $6.74\cdot10^{-9}$ | 3146 |

`loop` (FNL, 2 ranks), $t = 1$:

| Field loop, $t = 1$ | D_flux (seam, refined) | D_flux (max, uniform) | ⟨\|B_z\|⟩/A0 refined | uniform |
|---|---|---|---|---|
| GLM | $2.35\cdot10^{-4}$ | $2.03\cdot10^{-4}$ | $1.10\cdot10^{-3}$ | $6.86\cdot10^{-4}$ |
| EGLM | $2.45\cdot10^{-4}$ | $2.05\cdot10^{-4}$ | $1.91\cdot10^{-4}$ | $2.72\cdot10^{-6}$ |

Related cases: [G2](../seams-studies/glm-damping) varies the GLM damping length on the same `cpaw` leg;
[G3a](../seams-studies/reflux-not-ghosts) and [G3b](../seams-studies/eglm-reflux) rerun the `loop` leg and fix #78;
the uniform-grid cleaning is [GLM divergence cleaning](../mhd/glm-cleaning) and [EGLM cleaning](../mhd/eglm).

## Findings

- **GLM keeps the seam divergence bounded.** Over 10 periods of a source-free flow the seam D_flux has no secular
  growth at any resolution: the upwinded $(B_n, \psi)$ pair of FLUME bounds what PRISM's centred GLM lets run away
  ([#29](https://github.com/szaghi/adam/issues/29); #68 F3).
- **The seam floor converges at first order** (orders +3.31 then +1.36: the first pair is pre-asymptotic), the $O(H)$
  divergence of the $O(H^2)$ flux inconsistency of a 2:1 face (#68 F1), while the uniform floor converges at 4 to 6:
  the seams are the divergence error of a refined run, 3146 times the uniform level at $N = 256$.
- **#78 is not a divergence difference.** GLM and EGLM leave the same D_flux, refined and uniform; EGLM's 250-fold
  smaller $B_z$ on a uniform grid comes from its source term, whose centred divergence cancels the $B_z$ the flux
  generates. The $B_z$ rises 70 times on the refined grid, reproduced here to the printed digits of #78.
- **Attribution, corrected.** A2 attributed the rise to the centred source divergence reading the second-order
  restricted ghosts at the seam, where the cancellation fails, and placed the fix in Phase G3 (a source divergence
  consistent with the face fluxes at the seams). [G3a](../seams-studies/reflux-not-ghosts) showed that the ghosts are
  not the cause: the reflux is. [G3b](../seams-studies/eglm-reflux) fixed it by not refluxing $\mathbf{B}$ under EGLM.
