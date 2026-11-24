# Orszag–Tang vortex

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| MV-12 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) · P6a | P | CPU, FNL | yes | `orszag-tang` |

## Case and scope

The standard 2-D test of MHD turbulence and shock interactions (Stone et al. 2008, section 8.4): interacting shocks, a
current sheet and low-$\beta$ regions from a smooth, divergence-free, doubly periodic start. It has no exact solution,
so the oracle checks what the equations guarantee:

- **symmetry:** the 180° rotation about $(1/2, 1/2)$ maps the initial state to itself (scalars even, vectors odd), so
  the discrete solution must keep that symmetry. A wrong ghost fill or a sign error gives a truncation-level defect
  ($\sim10^{-3}$); a correct scheme leaves a round-off seed amplified by the flow;
- **conservation:** the periodic box conserves the integrals of $\rho$, $\rho\mathbf{u}$, $E$ and $\mathbf{B}$;
- **positivity:** the floors are disabled, so a non-positive density or pressure stops the run; a completed run has
  zero floored cells.

With `--amr` the centre box is refined once, so 2:1 seams are crossed by the shocks: symmetry and conservation must
hold there too, and the seam must not be a larger div(B) source than the flow features themselves.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$ | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Divergence control | GLM, $c_h = 4$ (above the fastest signal, `glm_ch_check = error`), `glm_alpha = 0.18`, `glm_damping_length = 1.0` | [models](../../models#divergence-control-mixed-glm) |
| Floors | disabled | [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_stepper_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_stepper_object.F90) (SELAH) |
| AMR (`--amr`) | geometric `primitive-box` marker at initialisation, no regrid; 2:1 restriction, tricubic ghost fill, reflux | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Diagnostics | conservation history every step; div(B) history (seam-local column on AMR) | `compute_conservation`, `compute_divb_history`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Boundaries | periodic | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid, uniform | the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $128^2$ cells ($32^2$ per block) |
| Grid, `--amr` | $32^2$ base ($8^2$ per block), the blocks whose centroid lies in $[0.25, 0.75]^2$ refined once (2:1); quadtree with `nk = 1` (`mhd/amr_box.py`; `--ratio 8` gives the octree with null $z$ and `nk = 4` of the measurements before #142); plus the uniform $32^2$ run as reference |
| Ranks | 2 |
| Initial condition | [`orszag-tang`](../../initial-conditions#orszag-tang), evaluated as odd functions of $2x - 1$, $2y - 1$ (bitwise symmetric) |
| Boundaries | periodic on $x$ and $y$ |
| Final time | $t = 0.5$ (1281 steps at $128^2$) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/orszag-tang/make_orszag_tang.py` | writes the input from `vortex/vortex-n064.ini` (`--cells`, `--refine-box`, `--ratio`, `--nk`, `--divergence-control`, `--it-save`) |
| `src/tests/flume/verification/mhd/amr_box.py` | refines a box by one 2:1 level at initialisation |
| `src/tests/flume/verification/mhd/orszag-tang/check.sh` | runs the uniform leg, or `--amr`, or `--divergence-control eglm` (`--numerics SPEC`, `--legs mv12`) |
| `src/tests/flume/verification/mhd/orszag-tang/orszag_tang_oracle.py` | the 180° symmetry, the conservation, the positivity |
| `src/tests/flume/verification/mhd/seam_divb_oracle.py` | the seam div(B) against the uniform reference (`--amr`) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_orszag_tang` | draws `orszag-tang.png` and `orszag-tang-amr.png` from `<runs>/orszag-tang` and `<runs>/orszag-tang-amr` (`--only orszag-tang`) |

```bash
cd src/tests/flume/verification/mhd
orszag-tang/check.sh                                                     # uniform 128^2, CPU (about 5 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl orszag-tang/check.sh    # FNL
orszag-tang/check.sh --amr                                               # 2:1 seams (about 30 s on the quadtree)
```

`make_doc_runs.sh <runs> mhd` writes the figure runs (`<runs>/orszag-tang`: $128^2$; `<runs>/orszag-tang-amr`: $32^2$ with $[0.25, 0.75]^2$ refined); `make_doc_figures.py --runs <runs>` draws it.

## Results

![Orszag-Tang](/flume/orszag-tang.png)

*Density and pressure at $t = 0.5$, $128^2$.*

![Orszag-Tang on AMR](/flume/orszag-tang-amr.png)

*Density on the $32^2$ base with the centre box refined; the block outlines show the 2:1 seams the shocks cross.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| 180° rotational symmetry at $t = 0.5$, relative defect | $10^{-8}$ | $2.7\cdot10^{-9}$ | $6.4\cdot10^{-9}$ |
| integrals of $\rho$, $\rho\mathbf{u}$, $E$, $\mathbf{B}$, relative drift | $10^{-13}$ | $1.4\cdot10^{-17}$ | $1.1\cdot10^{-16}$ |
| positivity (floors disabled) | run completes | pass | pass |
| `--amr`: symmetry | $10^{-8}$ | $4.6\cdot10^{-13}$ | $6.9\cdot10^{-13}$ |
| `--amr`: conservation | $10^{-13}$ | $5.6\cdot10^{-17}$ | $1.1\cdot10^{-16}$ |
| `--amr`: peak seam $\max\lvert\nabla\cdot\mathbf{B}\rvert$ / peak $\max\lvert\nabla\cdot\mathbf{B}\rvert$ of the uniform $32^2$ run | 1.0 | 0.76 ($0.884/1.170$) | the same |

*Measured: ADAM `0bcdf870` + the night-2026-11-22 tree (the paper commit, tag pending), CPU gfortran 16.2.0 + Open
MPI 5.0.11 (`iac_gnu`) on medooza (gnode22), FNL nvfortran 26.9, CUDA 12.9 (`iac_nvf`) on NVIDIA A30, 2026-11-23
(`orszag-tang/check.sh`, `orszag-tang/check.sh --amr`, jobs 1369376 CPU, 1369378 FNL).*

The symmetry defect of the uniform run grows from $1.7\cdot10^{-13}$ at step 100 to $3\cdot10^{-11}$ at step 300, when
the shocks form, and $1.3\cdot10^{-9}$ at $t = 0.5$ (CPU, M2-P6a); at the paper commit it ends at $2.7\cdot10^{-9}$ (CPU)
and $6.4\cdot10^{-9}$ (FNL), a round-off seed amplified by the flow on both backends (FMA contraction on FNL), below the
bound. $w$ and $B_z$ stay exactly 0 on both backends in these runs. The `--amr` rows run on the quadtree with
`nk = 1`; the octree values of M2 (null $z$, `nk = 4`) were $9.1\cdot10^{-13}$, $1.1\cdot10^{-16}$ and 0.76: the
trees agree to round-off, as MV-15 measures. The bounds are fixed tolerances, not octree baselines.

**Symmetry defect against the grid** ([#145](https://github.com/szaghi/adam/issues/145)). The uniform leg at $32^2$ and
$64^2$ (`check.sh --cells N`), and the defect at every 10th step of the uniform runs, the `--amr` run and two more
layouts (checkpoints every 10 steps; `orszag_tang_oracle.py` at each):

| Run | Cells, finest spacing | Steps to $t = 0.5$ | Defect at $t = 0.07$ | Defect at $t = 0.1$ | Defect at $t = 0.5$ |
|---|---|---|---|---|---|
| uniform $32^2$ (`--cells 32`) | 1024, 1/32 | 320 | $3.8\cdot10^{-14}$ | $2.5\cdot10^{-13}$ | $1.5\cdot10^{-13}$ |
| `--amr`: $32^2$, $[0.25, 0.75]^2$ at $64^2$ | 1792, 1/64 | 640 | $3.0\cdot10^{-14}$ | $8.8\cdot10^{-14}$ | $4.6\cdot10^{-13}$ |
| uniform $64^2$ (`--cells 64`) | 4096, 1/64 | 640 | $1.6\cdot10^{-13}$ | $1.2\cdot10^{-10}$ | $2.7\cdot10^{-10}$ |
| every block refined through the AMR path ($32^2$ base, box $[0, 1]^2$) | 4096, 1/64 | 640 | $1.6\cdot10^{-13}$ | $1.2\cdot10^{-10}$ | $2.7\cdot10^{-10}$ (bitwise the uniform $64^2$) |
| $64^2$, $[0.25, 0.75]^2$ at $128^2$ | 7168, 1/128 | 1281 | $2.3\cdot10^{-13}$ | $2.0\cdot10^{-11}$ | $1.8\cdot10^{-10}$ |
| uniform $128^2$ | 16384, 1/128 | 1281 | $9.2\cdot10^{-15}$ | $2.6\cdot10^{-11}$ | $2.7\cdot10^{-9}$ |

*The columns at $t = 0.07$ and $0.1$ are the nearest checkpoints ($32^2$: $t = 0.0625$ and $0.094$). Measured: ADAM
`0bcdf870` + the night-2026-11-22 tree (HEAD `95ffd728` with the working tree of 2026-11-23), CPU gfortran 16.2.0 +
Open MPI 5.0.11 (`iac_gnu`) on medooza (gnode25), 2 ranks, 2026-11-23 (`orszag-tang/check.sh --cells
32`, `--cells 64`; the series with `make_orszag_tang.py --it-save 10`; jobs 1369617, 1369627).*

The defect does not grow with the number of steps. It stays at $10^{-14}$–$10^{-13}$ until $t \approx 0.07$ on every
grid, then, on the grids that are fine outside the centre box, grows by three decades within $0.03$ time units
($64^2$: $4.7\cdot10^{-14}$ at $t = 0.0625$, $9.3\cdot10^{-11}$ at $t = 0.094$; $128^2$: $9.2\cdot10^{-15}$ at
$t = 0.070$, $2.3\cdot10^{-11}$ at $t = 0.094$), at cells near $(0.48, 0.96)$ and their images, outside
$[0.25, 0.75]^2$; afterwards it changes by less than a factor 4 on $64^2$, while $128^2$ keeps growing to
$2.7\cdot10^{-9}$. On
$32^2$ cells that growth does not happen. The `--amr` run has $32^2$ cells there: at the same step count and the same
finest spacing as the uniform $64^2$ run its defect is 600 times smaller. Refining every block through the AMR path
reproduces the uniform $64^2$ run bit for bit, so the AMR machinery (seam fill, reflux) does not make the run more
symmetric: the resolution outside the box does.

Related cases: [rotor](./rotor) uses the same symmetry oracle; [EGLM cleaning](./eglm) runs this problem under EGLM
(EV-4); [quadtree AMR](./quadtree-amr) runs the `--amr` problem on three trees; [MHD on `weno-riemann`](./weno-riemann)
runs it on HLLD.

## Findings

- The initial state is bitwise symmetric and the scheme is symmetric in exact arithmetic, but its left- and
  right-biased reconstructions round differently: the defect is a round-off seed amplified by the flow, four to five
  decades below what a wrong fill or a sign error produces.
- On AMR the refined box is symmetric about the centre and the reflux keeps the integrals: symmetry ($4.6\cdot10^{-13}$)
  and conservation ($5.6\cdot10^{-17}$) hold across the seams.
- The adaptive run is more symmetric than the uniform one because it is coarse where the defect grows, not because
  it takes fewer steps ([#145](https://github.com/szaghi/adam/issues/145)): the seed is amplified by three decades
  between $t = 0.07$ and $0.1$ near $(0.48, 0.96)$ and its image on $64^2$ and finer cells, and not on $32^2$ cells.
  With every block refined, the AMR path gives the uniform $64^2$ run bit for bit.
- The div(B) beside the seams stays below the uniform-grid peak: the shocks, not the seams, set the div(B) maximum.
- The AMR variant runs on a quadtree with `nk = 1`, true 2-D, since [#142](https://github.com/szaghi/adam/issues/142)
  (A8). It ran on an octree with null $z$ and `nk = 4` (four times the cells) because, before
  [#46](https://github.com/szaghi/adam/issues/46), a quadtree seam produced a spurious $z$ dependence. On this problem
  the two trees agree to round-off at $t = 0.2$ ($4.1\cdot10^{-12}$); [MV-15](./quadtree-amr) keeps the octree as its
  reference.
