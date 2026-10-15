# Field loop advection

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-9 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)); `--agree` [#79](https://github.com/szaghi/adam/issues/79) | CPU, FNL | yes | `field-loop` |

## Case and scope

A weak magnetic loop, $A_z = A_0 (R - r)$ with $A_0 = 10^{-3}$, $R = 0.3$ (plasma $\beta \approx 10^6$), is advected
obliquely across the periodic box (Gardiner & Stone 2005; Mignone & Tzeferacos 2010, section 4.4.1). The field is
divergence-free pointwise but jumps at $r = R$, and a collocated scheme cannot keep the discrete div(B) zero there.
With the out-of-plane velocity $w = 1$ the equations give

$$\partial_t B_z = w\,\nabla\cdot\mathbf{B},$$

so the generated $\langle\lvert B_z\rvert\rangle$ measures the divergence error that GLM leaves, and the decay of the
magnetic energy measures the numerical dissipation of the loop.

What the case catches:
- **a divergence error above the GLM level:** $\langle\lvert B_z\rvert\rangle/A_0$ under a bound and decreasing under
  refinement;
- **excess dissipation:** $E_B(T)/E_B(0)$ above a bound;
- **a seam div(B) source that does not decay** (`--amr`): the seam div(B) must fall from its peak;
- **backend disagreement on $\mathbf{B}$** (`--agree`): the CPU and FNL runs compared field by field.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$ | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables; `weno-riemann` under `--numerics` | `mhd_glm_face_split_fluxes`, `mhd_eigenvectors_core`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Divergence control | GLM, $c_h = 4$ (`glm_ch_check = error`), `glm_alpha = 0.18`, `glm_damping_length = 1.0` | [models](../../models#divergence-control-mixed-glm) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR (`--amr`) | `primitive-box` marker at initialisation; 2:1 restriction, tricubic ghost fill, reflux; seam div(B) column of the history | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `compute_divb_history`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Boundaries | periodic | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[-1, 1] \times [-0.5, 0.5]$, $z$ null |
| Grid | the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks) stretched to the box, $N \times N/2$ square cells, ladder $N = 64, 128$ |
| Grid, `--amr` | ladder $N = 32, 64$, $[0, 1] \times [-0.5, 0.5]$ refined once (2:1 seams at $x = 0$ and $x = 1$, crossed twice by the loop); octree with null $z$, `nk = 4` |
| Ranks | 2 |
| Initial condition | [`field-loop`](../../initial-conditions#field-loop-advected-field-loop): $\rho = 1$, $p = 1$, $\mathbf{v} = (2, 1, 1)$, loop at the origin, $R = 0.3$, $A_0 = 10^{-3}$ |
| Boundaries | periodic on $x$ and $y$ |
| Final time | $t = 2$ (two crossings along each axis); `--amr` $t = 1$ (one crossing) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/field-loop/make_field_loop.py` | writes the input from `vortex/vortex-n064.ini` (`--cells`, `--time-max`, `--glm-ch`, `--vz`, `--refine-box`) |
| `src/tests/flume/verification/mhd/field-loop/check.sh` | runs the ladder (`--amr`, `--agree`, `--numerics SPEC`) |
| `src/tests/flume/verification/mhd/field-loop/field_loop_oracle.py` | $\langle\lvert B_z\rvert\rangle/A_0$, its rate, $E_B(T)/E_B(0)$, the CPU–FNL agreement |
| `src/tests/flume/verification/mhd/seam_divb_oracle.py` | the seam div(B) decay (`--amr`) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_field_loop` | draws `field-loop.png` from `<runs>/field-loop` (`--only field-loop`) |

```bash
cd src/tests/flume/verification/mhd
field-loop/check.sh                                                     # CPU (about 10 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl field-loop/check.sh    # FNL
field-loop/check.sh --agree                                             # CPU against FNL, after both runs
field-loop/check.sh --amr                                               # 2:1 seams (about 28 min)
```

`make_doc_runs.sh` has no group that writes `<runs>/field-loop`: the figure run is not reproducible from the scripts
alone ([#85](https://github.com/szaghi/adam/issues/85)).

## Results

![Field loop](/flume/field-loop.png)

*Magnetic pressure $\lvert B\rvert^2/2$ of the loop at the final time.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| $\langle\lvert B_z\rvert\rangle/A_0$, $N = 64 / 128$ | $1.040\cdot10^{-3}$ / $6.003\cdot10^{-4}$ (baseline + 2%) | $1.02\cdot10^{-3}$ / $5.89\cdot10^{-4}$ | pass |
| its rate under refinement | $> 0$ | 0.79 (0.20 from 128 to 256) | pass |
| $E_B(T)/E_B(0)$, $N = 64 / 128$ | $\ge 0.87011$ / $0.93544$ (baseline − 0.2%) | 0.872 / 0.937 | pass |
| `--amr`: $\langle\lvert B_z\rvert\rangle/A_0$, $N = 32 / 64$ | $2.277\cdot10^{-3}$ / $1.118\cdot10^{-3}$ | $2.23\cdot10^{-3}$ / $1.10\cdot10^{-3}$, rate +1.03 | pass |
| `--amr`: $E_B(T)/E_B(0)$, $N = 32 / 64$ | $\ge 0.80066$ / $0.89888$ | 0.802 / 0.901 | pass |
| `--amr`: final seam div(B) / its peak | 0.1 | 0.033 / 0.011 | pass |
| `--agree`: CPU against FNL, field by field ($\mathbf{B}$ relative to $A_0$) | $10^{-10}$ | — | $\sim10^{-14}$ at step 4 |

The uniform $N = 256$ run (26 more minutes, not in the ladder) gives $5.13\cdot10^{-4}$ and $E_B$ 0.967. The uniform
$32/64$ runs that the `--amr` ladder is compared with give $1.99\cdot10^{-3}$, $1.11\cdot10^{-3}$ and $E_B$ 0.759,
0.884. The `weno-riemann` energies are in [MHD on `weno-riemann`](./weno-riemann).

Related cases: [glm cleaning](./glm-cleaning) measures the GLM transport and damping on an exact solution;
[runtime regrids](../runtime-amr/runtime-regrids) and [accuracy of runtime AMR](../runtime-amr/accuracy) track this
loop with a moving refinement (regression `field-loop-regrid`); the M6 seam studies
[G3a](../seams-studies/reflux-not-ghosts) and [G3b](../seams-studies/eglm-reflux) run it across a 2:1 seam under GLM
and EGLM.

## Findings

- The loop field jumps at $r = R$, so the discrete div(B) there is $\sim A_0/h$ over a ring $\sim h$ wide: the $B_z$
  source is $O(1)$ in $h$ and the decrease is slow. No convergence order is expected; the measured level matches the
  GLM level of Mignone & Tzeferacos ($\sim10^{-3}$).
- On AMR the refined half dissipates less ($E_B$ 0.802, 0.901 against 0.759, 0.884 uniform); the seams add $B_z$ at
  $N = 32$ and none at $N = 64$. The seam div(B) peaks at the first step (the initial-data error of the loop edge) and
  GLM removes it; the collocated PRISM seam runs away instead ([#29](https://github.com/szaghi/adam/issues/29)).
- The backends agree to round-off on $\mathbf{B}$ as on the fluid. Until [#79](https://github.com/szaghi/adam/issues/79)
  they differed by $\sim10^{-6}$ on $\mathbf{B}$: the slow-wave normalisation of the eigenvectors cancelled
  catastrophically; it is now cancellation-free ([numerics](../../numerics#scheme-space-weno-flux-splitting)).
