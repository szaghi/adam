# Rotated shock tube

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-8 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | no ([#85](https://github.com/szaghi/adam/issues/85)) | `rotated-shock-tube` |

## Case and scope

Ryu–Jones 1a rotated in the $x$–$y$ plane (Tóth 2000, section 6.3.2) makes the dimension-by-dimension scheme solve a
1-D problem on truly 2-D data. The exact normal field stays uniform, $B_\parallel = 5/\sqrt{4\pi}$, so its deviation
measures the divergence error of the base scheme with GLM, and the whole solution must match the same problem solved in
1-D.

What the case catches:
- **a divergence error of the multi-dimensional scheme:** $\delta B_\parallel$ (Tóth eq. 45, against the analytic value)
  under a bound, decreasing with $N$;
- **a 2-D solution that departs from the 1-D one:** the mean relative $L_1$ error of $\rho$, $v_\parallel$, $v_\perp$,
  $p$, $B_\parallel$, $B_\perp$ against a 1-D reference run, under a bound, decreasing with $N$.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$ | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Divergence control | GLM, $c_h = 24$ (the fastest signal, $\sim$20–24 in the ramp; `glm_ch_check = error`), `glm_alpha = 0.18`, `glm_damping_length = 1.0` | [models](../../models#divergence-control-mixed-glm) |
| Floors | disabled | [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | periodic (FLUME has no shifted-periodic boundary) | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Problem | Ryu & Jones (1995) 1a in the frame of the normal: $(\rho, v_\parallel, v_\perp, p, B_\parallel, B_\perp) = (1, 10, 0, 20, 5/\sqrt{4\pi}, 5/\sqrt{4\pi})$ left, $(1, -10, 0, 1, 5/\sqrt{4\pi}, 5/\sqrt{4\pi})$ right |
| Initial condition | [`rotated-riemann`](../../initial-conditions#rotated-riemann-1-d-riemann-problem-rotated-in-x-y-periodic-along-its-normal), $s = n_x x + n_y y$, period 1: the jump at $s = 1/4 + 1/(4N)$ (off every cell centre), a linear ramp of the primitive variables back to the left state over $s \in (1/2, 9/10]$ |
| Normals | $(1, 2)$, Tóth's $\tan^{-1} 2 = 63.4°$; $(1, 1)$, 45° |
| 2-D grid | $[0, 1]^2$, the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $N^2$ cells, $N = 128, 256$, periodic |
| 1-D reference | the same problem along $x$ with normal $(\lvert n\rvert, 0)$ on $[0, 1/\lvert n\rvert]$, the Sod grid of [V1](../euler/sod) (4 blocks along $x$), $4 \times 256 = 1024$ cells per period (4 times the normal resolution of $N = 256$), periodic |
| Ranks | 2 |
| Final time | $t = 0.04/\lvert n\rvert$: the fan of the jump (fast speeds $\sim5$) stays within $0.2/\lvert n\rvert$ of it, as in Tóth ($t = 0.08\cos a$ on his strip; fast shocks at 0.1 and 0.9 at $t = 0.08$) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/rotated-shock-tube/make_rotated_shock_tube.py` | writes the 2-D input from `vortex/vortex-n064.ini` or the 1-D reference from `sod/sod-x.ini` (`--angle 45\|63`, `--cells`, `--reference`, `--cfl`, `--glm-ch`) |
| `src/tests/flume/verification/mhd/rotated-shock-tube/check.sh` | runs the reference and the ladder for both normals (`--np N`; about 30 min) |
| `src/tests/flume/verification/mhd/rotated-shock-tube/rotated_shock_tube_oracle.py` | $\delta B_\parallel$ and the mean relative error against the reference |

```bash
cd src/tests/flume/verification/mhd
rotated-shock-tube/check.sh                                                     # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl rotated-shock-tube/check.sh    # FNL
```

No figure yet: [#85](https://github.com/szaghi/adam/issues/85) lists the $B_\parallel$ deviation map at 63.4° and 45°
and the 1-D cut against the exact profile.

## Results

| Check | Bound (CPU baseline M2-P6d + 2%) | CPU | FNL |
|---|---|---|---|
| $\delta B_\parallel$, 63.4°, $N = 128 / 256$ | $6.146\cdot10^{-3}$ / $2.842\cdot10^{-3}$ | below / $2.79\cdot10^{-3}$ | pass |
| mean relative error, 63.4°, $N = 128 / 256$ | $8.989\cdot10^{-2}$ / $4.095\cdot10^{-2}$ | below | pass |
| $\delta B_\parallel$, 45°, $N = 128 / 256$ | $6.424\cdot10^{-3}$ / $4.252\cdot10^{-3}$ | below / $4.17\cdot10^{-3}$ | pass |
| mean relative error, 45°, $N = 128 / 256$ | $5.440\cdot10^{-2}$ / $2.600\cdot10^{-2}$ | below | pass |
| both measures decreasing with $N$ | strict | pass | pass |

Reference values of Tóth's second-order base scheme at 63.4°, $N = 256$ on his strip (twice the cells across the fan
of $N = 256$ here): $\delta B_\parallel = 3.7\cdot10^{-3}$, mean error 0.0238.

FNL ran both normals (work directories of both backends); the tables record the CPU values.

Related cases: [Ryu–Jones 2a](./ryu-jones-2a) is the axis-aligned MHD Riemann problem; [field loop](./field-loop)
measures the divergence error of a smooth field; [GLM cleaning](./glm-cleaning) the cleaning itself.

## Findings

- $\delta B_\parallel$ is below Tóth's base-scheme value at 63.4° ($2.79\cdot10^{-3}$ against $3.7\cdot10^{-3}$) with
  half his cells across the fan.
- The 45° normal is not degenerate for this scheme, unlike Tóth's CT and CD schemes: $B_\parallel$ is conserved to
  truncation, not to round-off.
- The convergence of $\delta B_\parallel$ is about first order at 63.4° and slower at 45° (rates about 1.1 and 0.6
  from 128 to 256, computed from the baselines behind the bounds).
- A periodic line forces the velocity jump to be undone. A second jump would be the mirror problem, an expansion close
  to vacuum (density $\sim2\cdot10^{-4}$, Alfvén speed $\sim100$) where $\max(\lvert u\rvert + c_f)$ grows with the
  resolution of the fan and no $c_h$ holds; the resolved ramp (density $\sim1/3$ at the final time) avoids it, and the
  1-D reference solves the same problem, ramp and wave interactions included.
