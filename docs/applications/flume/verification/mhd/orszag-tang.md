# Orszag–Tang vortex

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-12 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes | `orszag-tang` |

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
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR (`--amr`) | geometric `primitive-box` marker at initialisation, no regrid; 2:1 restriction, tricubic ghost fill, reflux | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `apply_reflux_corrections`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Diagnostics | conservation history every step; div(B) history (seam-local column on AMR) | `compute_conservation`, `compute_divb_history`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Boundaries | periodic | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid, uniform | the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $128^2$ cells ($32^2$ per block) |
| Grid, `--amr` | $32^2$ base ($8^2$ per block), the blocks whose centroid lies in $[0.25, 0.75]^2$ refined once (2:1); octree with null $z$ and `nk = 4` (`mhd/amr_box.py`); plus the uniform $32^2$ run as reference |
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
orszag-tang/check.sh                                                     # uniform 128^2, CPU (about 4.5 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl orszag-tang/check.sh    # FNL
orszag-tang/check.sh --amr                                               # 2:1 seams (about 6 min)
```

`make_doc_runs.sh <runs> mhd` writes the figure runs (`<runs>/orszag-tang`: $128^2$; `<runs>/orszag-tang-amr`: $32^2$ with $[0.25, 0.75]^2$ refined); `make_doc_figures.py --runs <runs>` draws it.

## Results

![Orszag-Tang](/flume/orszag-tang.png)

*Density and pressure at $t = 0.5$, $128^2$.*

![Orszag-Tang on AMR](/flume/orszag-tang-amr.png)

*Density on the $32^2$ base with the centre box refined; the block outlines show the 2:1 seams the shocks cross.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| 180° rotational symmetry at $t = 0.5$, relative defect | $10^{-8}$ | $1.3\cdot10^{-9}$ | pass |
| integrals of $\rho$, $\rho\mathbf{u}$, $E$, $\mathbf{B}$, relative drift | $10^{-13}$ | $1.2\cdot10^{-14}$ | pass |
| positivity (floors disabled) | run completes | pass | pass |
| `--amr`: symmetry | $10^{-8}$ | $9.1\cdot10^{-13}$ | pass |
| `--amr`: conservation | $10^{-13}$ | $1.1\cdot10^{-16}$ | pass |
| `--amr`: peak seam $\max\lvert\nabla\cdot\mathbf{B}\rvert$ / peak $\max\lvert\nabla\cdot\mathbf{B}\rvert$ of the uniform $32^2$ run | 1.0 | 0.76 ($0.887/1.170$) | pass |

The symmetry defect of the uniform run grows from $1.7\cdot10^{-13}$ at step 100 to $3\cdot10^{-11}$ at step 300, when
the shocks form, and $1.3\cdot10^{-9}$ at $t = 0.5$ (CPU, M2-P6a). $w$ and $B_z$ stay exactly 0 on the CPU; on FNL they
carry round-off noise ($\sim10^{-12}$ of the in-plane field), which the oracle scales by the vector magnitude. FNL ran
every leg; the tables record the CPU values.

Related cases: [rotor](./rotor) uses the same symmetry oracle; [EGLM cleaning](./eglm) runs this problem under EGLM
(EV-4); [quadtree AMR](./quadtree-amr) runs the `--amr` problem on three trees; [MHD on `weno-riemann`](./weno-riemann)
runs it on HLLD.

## Findings

- The initial state is bitwise symmetric and the scheme is symmetric in exact arithmetic, but its left- and
  right-biased reconstructions round differently: the defect is a round-off seed amplified by the flow, four to five
  decades below what a wrong fill or a sign error produces.
- On AMR the refined box is symmetric about the centre and the reflux keeps the integrals: symmetry ($9.1\cdot10^{-13}$)
  and conservation ($1.1\cdot10^{-16}$) hold across the seams.
- The div(B) beside the seams stays below the uniform-grid peak: the shocks, not the seams, set the div(B) maximum.
- The AMR variant runs on an octree with null $z$ and `nk = 4` because, before
  [#46](https://github.com/szaghi/adam/issues/46), a quadtree seam produced a spurious $z$ dependence; the comparison
  with a quadtree is [MV-15](./quadtree-amr).
