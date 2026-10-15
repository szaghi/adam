# GLM divergence cleaning

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-3, MV-10 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)); the MV-10 monitor M6 ([#68](https://github.com/szaghi/adam/issues/68), G1) | CPU, FNL | yes | `glm-pulse` |

## Case and scope

**MV-3, the GLM pulse.** A Gaussian pulse in $B_x$ along $x$ on a fluid at rest is a pure divergence error. In 1-D the
pair $(B_x, \psi)$ of mixed GLM then obeys a (damped) telegraph equation decoupled from the fluid,

$$\partial_t B_x + \partial_x \psi = 0, \qquad \partial_t \psi + c_h^2\,\partial_x B_x = -\frac{c_h^2}{c_p^2}\,\psi,$$

with an exact solution: the error splits into two pulses travelling at $\pm c_h$ and is damped at the rate
$c_h^2/c_p^2 = \alpha c_h / L$ ($\alpha$ = `glm_alpha`, $L$ = `glm_damping_length`). The oracle computes it mode by
mode (d'Alembert without damping).

What MV-3 catches:
- **a wrong GLM block in the face flux:** $L_1$ of $B_x$ and $\psi/c_h$ against the exact solution under a bound, and
  the observed order;
- **a wrong damping source or rate:** a wrong rate is an $O(1)$ relative error (1.23 with no damping at $N = 128$,
  against $9.2\cdot10^{-5}$), and the `min-cell` damping length is checked on its own;
- **$c_h$ missing from the time step:** $c_h = 2$ is above the fastest wave (1.29), so it sets $\Delta t$, and
  `glm_ch_check = error` must not trip;
- **a wrong wall parity of $\psi$:** reflecting walls ($B_x$ odd, $\psi$ even);
- **a leaking `none` path:** without cleaning $B_x$ must stay frozen, bit for bit.

**MV-10, the Dedner peak.** A peak in $B_x$ on a uniform magnetised flow (Dedner et al. 2002, section 5.1) has a
non-zero $\nabla\cdot\mathbf{B}$ from $t = 0$. Without cleaning the error is only advected; GLM transports it away at
$c_h$ and damps it, so $\|\nabla\cdot\mathbf{B}\|_1$ must decay, on a uniform grid and across 2:1 AMR seams. The leg
also checks the div(B) history, the derived output fields and the `divb_tol` monitor.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD with mixed GLM, $\gamma = 5/3$ | [models](../../models#divergence-control-mixed-glm); `mhd_glm_flux`, `mhd_glm_eigenvectors`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Space | WENO5 flux splitting on characteristic variables; the $(B_n, \psi)$ pair is the 2×2 block of speeds $\mp c_h$ | `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Damping | source $-(c_h^2/c_p^2)\,\psi$, $c_h^2/c_p^2 = \alpha c_h/L$; $L$ a length or `min-cell` (the minimum cell size) | `add_glm_damping`, [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90); `set_glm_damping`, [`adam_flume_mhd_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_object.F90) |
| $c_h$ check | `glm_ch_check = error` (MV-3, MV-10) | `check_glm_ch`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| div(B) history and monitor | max, $L_1$ and seam-local max of $\nabla\cdot\mathbf{B}$ (host FD, null directions weighted zero); `divb_tol` on the normalised $D_\text{flux}$, $\max h\lvert D_\text{flux}\rvert/\max\lvert B\rvert$ | `compute_divb_history`, `compute_dflux_norms`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `save_divb_row`, [`adam_flume_diagnostics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_diagnostics_object.F90) |
| Time | SSP Runge–Kutta, 3 stages, CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR (MV-10) | 2:1 seams: restriction, tricubic ghost fill, reflux | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Boundaries | periodic; reflecting walls (`wall-inviscid`: $B_x$ odd, $\psi$ even) in the wall leg | [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); [boundary conditions](../../boundary-conditions#wall-inviscid-slip-wall-by-mirroring-with-a-sign-vector) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90), [`adam_flume_fnl_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_mhd_glm_kernels.F90) |

## Initial and boundary conditions

**MV-3** (`mhd/glm-pulse`):

| Item | Value |
|---|---|
| Domain | $[0, 1]$ along $x$; $y$ and $z$ null |
| Grid | the Sod grid of [V1](../euler/sod), 4 blocks along $x$, blocks of $N/4 \times 2 \times 2$ cells |
| Ranks | 2 |
| Initial condition | [`glm-pulse`](../../initial-conditions#glm-pulse-gaussian-pulse-in-one-field-component): background $(\rho, u, v, w, p) = (1, 0, 0, 0, 1)$, $\mathbf{B} = 0$; pulse in $B_x$ at $x = 0.5$, width 0.1, amplitude 0.1 |
| Cleaning | $c_h = 2$ |

| Leg | $N$ | `glm_alpha` | `glm_damping_length` | $x$ boundaries | $t$ |
|---|---|---|---|---|---|
| d'Alembert | 64, 128, 256 | 0 | 1.0 | periodic | 0.3 (the pulses cross the boundary) |
| damped | 64, 128, 256 | 0.5 | 0.1 ($c_h^2/c_p^2 = 10$) | periodic | 0.3 |
| `min-cell` | 128 | 0.1 (rate 25.6) | `min-cell` | periodic | 0.3 |
| walls | 128 | 0 | 1.0 | `wall-inviscid` (both pulses reflected once) | 0.4 |
| `none` | 128 | — | — | periodic | 0.3 |

**MV-10** (`mhd/divb-peak`):

| Item | Value |
|---|---|
| Domain | $[-0.5, 1.5]^2$, $z$ null, periodic |
| Grid | the periodic grid of [V3](../euler/amr-conservation) (octree) at level 2: $64^2$ cells, $16^2$ per block, `nk = 4` (the 2:1 seam fill needs 4 cells on every axis, null ones included) |
| AMR leg | $8^2$ cells per base block ($32^2$), the blocks of $[-0.5, 0.5]^2$ refined once more (a 2:1 seam the peak crosses near $t = 0.375$) |
| Ranks | 2 |
| Initial condition | [`divb-peak`](../../initial-conditions#divb-peak-dedner-et-al-peak-in-b-x): $\rho = 1$, $u = v = 1$, $w = 0$, $p = 6$, $\mathbf{B} = (B_x, 0, 1/\sqrt{4\pi})$, $B_x \mathrel{+}= (1 - (r/R)^2)^2/\sqrt{4\pi}$, $R = 1/8$, at the origin |
| Cleaning | GLM, $c_h = 5$, `glm_alpha = 0.4`, damping length `min-cell`; and `none` |
| Final time | $t = 0.5$; the monitor legs 5 steps |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/glm-pulse/make_glm_pulse.py` | writes the MV-3 input from `sod/sod-x.ini` |
| `src/tests/flume/verification/mhd/glm-pulse/check.sh` | runs the MV-3 legs (`--np N`, `--divergence-control glm\|eglm`) |
| `src/tests/flume/verification/mhd/glm-pulse/glm_pulse_oracle.py` | the exact telegraph solution, the $L_1$ bounds, the order, $\int B_x$, the frozen `none` run |
| `src/tests/flume/verification/mhd/divb-peak/make_divb_peak.py` | writes the MV-10 input from `conservation/amr-periodic.ini` |
| `src/tests/flume/verification/mhd/divb-peak/check.sh` | runs the MV-10 legs and the monitor legs |
| `src/tests/flume/verification/mhd/divb-peak/divb_peak_oracle.py` | the decay, $\int\mathbf{B}$, $\max\lvert\psi\rvert$, the derived fields, the seam column |
| `src/tests/flume/verification/make_doc_figures.py` `fig_glm_pulse` | draws `glm-pulse.png` from `<runs>/glm-pulse` (`--only glm-pulse`) |

```bash
cd src/tests/flume/verification/mhd
glm-pulse/check.sh                                                     # MV-3, CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl glm-pulse/check.sh    # MV-3, FNL
divb-peak/check.sh                                                     # MV-10
```

`make_doc_runs.sh` has no group that writes `<runs>/glm-pulse`: the figure run is not reproducible from the scripts
alone ([#85](https://github.com/szaghi/adam/issues/85)).

## Results

![GLM pulse](/flume/glm-pulse.png)

*The undamped case (`glm_alpha = 0`, periodic, 128 cells, $c_h = 2$, $t = 0.3$): $B_x$ and $\psi$ of FLUME on the
exact solution, at $t = 0$ and at the final time.*

MV-3 (bounds: CPU baseline M2-P4a + 2%):

| Check | Bound | CPU | FNL |
|---|---|---|---|
| d'Alembert, $L_1$ of $B_x$ and $\psi/c_h$, $N = 64 / 128 / 256$ | $1.584\cdot10^{-4}$ / $9.342\cdot10^{-6}$ / $7.487\cdot10^{-7}$ | $1.553\cdot10^{-4}$ / $9.159\cdot10^{-6}$ / $7.340\cdot10^{-7}$ | pass |
| damped, $L_1$, $N = 64 / 128 / 256$ | $3.393\cdot10^{-5}$ / $1.962\cdot10^{-6}$ / $1.572\cdot10^{-7}$ | $3.327\cdot10^{-5}$ / $1.923\cdot10^{-6}$ / $1.541\cdot10^{-7}$ | pass |
| observed order, both ladders | $\ge 3.0$ (SSP-33 in time) | 4.1 and 3.6 (WENO5 in space still visible) | pass |
| `min-cell`, $L_1$ | $1.487\cdot10^{-7}$ | $1.458\cdot10^{-7}$ | pass |
| walls, $L_1$ | $1.364\cdot10^{-5}$ | $1.337\cdot10^{-5}$ | pass |
| $\int B_x$ drift (d'Alembert) | $10^{-14}$ ($\int B_x \approx 1.8\cdot10^{-2}$) | $\le 5.8\cdot10^{-16}$ | pass |
| `none`: $B_x$ last against first checkpoint | bitwise | bitwise | bitwise |

MV-10 (CPU baselines M2-P4b):

| Check | Bound | CPU | FNL |
|---|---|---|---|
| GLM final / initial $\|\nabla\cdot\mathbf{B}\|_1$, uniform $64^2$ (never above the initial value) | 0.1 | 0.057 | pass |
| GLM / `none` final $\|\nabla\cdot\mathbf{B}\|_1$ | 0.1 | 0.027 | pass |
| AMR: GLM final / initial | 0.2 | 0.091 | pass |
| $\int B_x$, $\int B_y$, $\int B_z$ drift ($\int B_z = 1.13$) | $10^{-11}$ | $\le 8.3\cdot10^{-13}$ | pass |
| $\max\lvert\psi\rvert$ | $0.141$ ($c_h A/10$) | $3.6\cdot10^{-3}$ uniform, $1.3\cdot10^{-2}$ AMR | pass |
| saved `pt`, `beta`, `bmag` against their definitions; saved $\max\lvert\nabla\cdot\mathbf{B}\rvert$ against the history | $10^{-12}$ | exact | pass |
| AMR seam column: zero at step 0, never above the global maximum, positive at the end | — | pass | pass |
| monitor, `divb_tol = 0.01`, `divb_error = .false.` (the peak starts near 0.2) | warnings, run completes | pass | pass |
| monitor, `divb_error = .true.` | run stops with its message | pass | pass |
| `divb_tol > 0` without `[diagnostics] divb_flux` | refused at initialisation | pass | pass |

FNL ran every leg (work directories of both backends); the tables record the CPU values, and the backends agree to
the printed digits.

Related cases: [EGLM cleaning](./eglm) runs MV-3 under EGLM (EV-2); [Orszag–Tang](./orszag-tang) and
[field loop](./field-loop) measure the seam div(B) of shocked and smooth flows; the seam div(B) studies of M6 are in
[divergence across seams](../seams/divergence-across-seams).

## Findings

- The observed order is that of SSP-33 in time with the WENO5 spatial error still visible (4.1 and 3.6), above the
  asserted 3.0.
- The $\|\nabla\cdot\mathbf{B}\|_1$ of the Dedner peak is not monotone under GLM: a few-percent transient while the peak
  spreads into outgoing waves. The leg therefore asserts "never above the initial value" and the final ratio, not a
  monotone decay.
- On AMR the final ratio doubles (0.091 against 0.057): the 2:1 seam injects div(B) (the seam source of
  [#29](https://github.com/szaghi/adam/issues/29)); GLM keeps it bounded.
- Since #68 G1 the `divb_tol` monitor acts on the normalised $D_\text{flux}$, which needs `[diagnostics] divb_flux` and
  GLM; a `divb_tol` without `divb_flux` is refused at initialisation instead of monitoring nothing.
