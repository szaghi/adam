# GLM damping length on a seam

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| G2 | Seam studies | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | FNL only | yes | — |

## Case and scope

The GLM damping rate is $c_h^2/c_p^2 = \alpha\, c_h / L$ with `[mhd] glm_damping_length` $L$ a fixed length or
`min-cell` (the finest spacing; with runtime AMR it is frozen at the spacing of `max_level`, D-M5-9). A `per-level`
length, each level damping at its own spacing, does not exist. G2 asks whether it would matter, bracketing it between
the coarse and the finest spacings on A2's source-free CPAW.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD with mixed GLM, $\gamma = 5/3$, `glm_ch` 2.0, `glm_alpha` 0.18 | [models](../../models#divergence-control-mixed-glm) |
| Damping length | `glm_damping_length` = 1.0, the coarse spacing $L_x/N$, or `min-cell` | `set_glm_damping`, [`adam_flume_mhd_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_object.F90) |
| Space, time | WENO5 flux splitting on characteristic variables; SSP Runge–Kutta 5 stages, CFL 0.4 | as in [AD](../seams/divergence-across-seams) |
| Seam operators | `mean` + `tricubic`, reflux on | [numerics](../../numerics#adaptive-mesh-refinement) |
| D_flux | the seam floor of `[diagnostics] divb_flux` | `compute_dflux_norms`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Case | the `cpaw` leg of [AD](../seams/divergence-across-seams): source-free CPAW, $[0, \sqrt5] \times [0, \sqrt5/2]$, quadtree, two $x$ seams ($x \in [0.6, 1.6]$ refined) |
| Grid | refined and uniform $N$ = 64 / 128 |
| Boundaries | periodic in $x$, $y$; `extrapolation` in $z$ |
| Final time | 10 periods |
| Ranks | 2, FNL |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-divb/make_amr_divb.py cpaw` | the CPAW inputs (`--glm-damping-length`: 1.0, the coarse spacing, `min-cell`) |
| `src/tests/flume/verification/amr-divb/amr_divb_oracle.py growth` | the seam floor and its growth |
| `src/tests/flume/verification/make_doc_runs.sh <runs> divb` (FNL), `make_doc_figures.py --only glm-damping` | re-runs the CPAW with the coarse spacing and `min-cell` (`--glm-damping-length`), $N = 64$ and 128; the 1.0 runs are those of AD |

## Results

![G2: the seam D_flux per damping length](/flume/glm-damping.png)

*FNL: the seam maximum of D_flux on the source-free CPAW for the three damping lengths, refined $N = 64$ and 128. The re-runs reproduce the table below (the medians of the second half): every length bounded, `min-cell` doubling the floor at $N = 128$.*

| `glm_damping_length` | rate $c_h^2/c_p^2$ (N = 128) | seam D_flux floor N = 64 / 128 | × uniform (N = 128) | growth N = 128 |
|---|---|---|---|---|
| 1.0 (A2) | 0.36 | $5.39\cdot10^{-4}$ / $5.43\cdot10^{-5}$ | 119 | 0.72 |
| coarse spacing $L_x/N$ | 20.6 | $6.96\cdot10^{-4}$ / $5.32\cdot10^{-5}$ | 116 | 1.02 |
| `min-cell` (finest spacing) | 41.2 | $8.30\cdot10^{-4}$ / $1.07\cdot10^{-4}$ | 233 | 1.01 |

Related cases: [AD](../seams/divergence-across-seams) is the A2 baseline (first row); the uniform-grid cleaning is
[GLM divergence cleaning](../mhd/glm-cleaning).

## Findings

- **Every length keeps the seam divergence bounded** (growth 0.72 to 1.02 over 10 periods).
- **Strong damping raises the seam floor:** `min-cell` doubles it at $N = 128$ against the coarse spacing or the
  domain-scale length; a per-level length (between the coarse and the finest spacings) would sit between them.
- **The D-M5-9 freeze stays** (the plan's choice): its cost on this test is at most that factor 2 on the floor, with the
  damping rate fixed through regrids. No `per-level` option is added.
