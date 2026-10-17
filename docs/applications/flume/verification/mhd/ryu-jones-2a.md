# Ryu–Jones 2a: all seven waves

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-4 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes | `rj2a-x`, `rj2a-y`, `rj2a-z` |

## Case and scope

The 1-D Riemann problem of Ryu & Jones (1995), case 2a, produces the full MHD wave fan from one discontinuity: two
fast shocks, two rotational discontinuities, two slow shocks and a contact. Its exact solution is known (7 states,
Dai & Woodward 1994, Tables Ia and Ib, as used by the Athena++ problem generator `shock_tube.cpp`), so every
conservative variable is compared with it, in the frame (normal, tangent 1, tangent 2).

What the case catches:
- **defects anywhere on the MHD path:** the eigensystem at the face average, the characteristic splitting, the fluxes
  and the rotation of the transverse components;
- **direction-dependent defects:** the problem runs along $x$, $y$ and $z$; the $y$ and $z$ runs must equal the $x$ run
  **bitwise** in the rotated frame;
- **loss of accuracy:** the $L_1$ sum over the 8 conservative variables under a recorded bound, and decreasing from 256
  to 512 cells;
- **a leaking cleaning block:** in 1-D $B_n$ is uniform, so with GLM the $(B_n, \psi)$ block is inert; the GLM run must
  equal the run without cleaning **bitwise**, with $\psi$ exactly zero.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$ ($c_p = 2.5$, $c_v = 1.5$), code units ($B = B_\text{Gauss}/\sqrt{4\pi}$) | [models](../../models#ideal-mhd); [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `mhd_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Characteristic projection | the 7×7 Roe–Balsara core of Stone et al. (2008) at the arithmetic face average | `mhd_face_average`, `mhd_eigenvectors`, `mhd_eigenvectors_core`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Order-independent sums | $\lvert u\rvert^2$, $\lvert B\rvert^2$, $u\cdot B$ summed independently of the term order (the bitwise x/y/z property) | `mhd_sum3`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Divergence control | `none` (main legs); GLM leg: `glm_ch = 2.0`, `glm_alpha = 0.18`, `glm_damping_length = 1.0`, `glm_ch_check = warning` | `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#divergence-control-mixed-glm) |
| Floors | disabled (`rho_floor = p_floor = 0`) | [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | zeroth-order extrapolation | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant) |
| FNL backend | the same kernels on the device | `compute_residuals_weno_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

The GLM leg sets $c_h = 2$ below the fastest wave (the left state alone has $\lvert u\rvert + c_f = 2.89$), so the fluid
sets the time step and the GLM run takes the same steps as the run without cleaning; the $c_h$ check warns, as it should.

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]$ along the active axis; the other two axes null |
| Grid | the Sod grid of [V1](../euler/sod): 64 blocks (octree), 4 along the active axis, refined to $N$ cells: $N = 256$ (blocks of $64 \times 4 \times 4$) and $N = 512$ ($128 \times 4 \times 4$), `ngc = 3` |
| Ranks | 2 |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), two regions, discontinuity at $0.5$ |
| Left state $(\rho, u_n, u_{t1}, u_{t2}, p)$ | $(1.08, 1.2, 0.01, 0.5, 0.95)$, $\sqrt{4\pi}\,(B_n, B_{t1}, B_{t2}) = (2, 3.6, 2)$ |
| Right state | $(1, 0, 0, 0, 1)$, $\sqrt{4\pi}\,(B_n, B_{t1}, B_{t2}) = (2, 4, 2)$ |
| Frame | tangents cyclic as in FLUME and Athena++: axis $y \to (y, z, x)$, axis $z \to (z, x, y)$ |
| Boundaries | `extrapolation` on all six faces |
| Final time | $t = 0.2$ |
| Legs | $N = 256$ along $x$, $y$, $z$ (no cleaning); $N = 512$ along $x$; $N = 256$ along $x$ with GLM; $N = 256$ along $x$ with EGLM ([EV-3](./eglm)) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/rj2a/make_rj2a.py` | writes the input from `sod/sod-{x,y,z}.ini` (`--axis`, `--cells`, `--divergence-control none\|glm\|eglm`) |
| `src/tests/flume/verification/mhd/rj2a/check.sh` | runs the legs (`--np N`, `--numerics SPEC`, `--legs LIST`: `mv4`, `rv6` (with `--numerics hlld[:...]` only), `ev3`; `rv6` and `ev3` need `mv4`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/mhd/rj2a/rj2a_oracle.py` | the exact 7-state solution, the $L_1$ bound, the 1-D consistency, the bitwise direction and pair checks |
| `src/tests/flume/verification/mhd/numerics.sh` | the shared `--numerics` option ([MHD on `weno-riemann`](./weno-riemann)) |
| `src/tests/flume/verification/make_doc_runs.sh` group `rj2a` | re-runs `rj2a` (splitting) and `rj2a-hlld` at 256 cells with checkpoints kept |
| `src/tests/flume/verification/make_doc_figures.py` `fig_mhd_riemann` | draws `rj2a.png` (`--only mhd-riemann`) |

```bash
cd src/tests/flume/verification/mhd
rj2a/check.sh                                                     # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl rj2a/check.sh    # FNL
rj2a/check.sh --numerics hlld                                     # weno-riemann, HLLD, characteristic
```

## Results

![Ryu-Jones 2a](/flume/rj2a.png)

*Density, normal velocity, pressure and the two tangential fields at $t = 0.2$, 256 cells: the splitting scheme and
`weno-riemann` HLLD (characteristic interpolation) against the exact solution; at this scale the two are
indistinguishable.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| $L_1$ sum, $N = 256$, $x$ | $3.919\cdot10^{-2}$ (baseline + 2%) | $3.842\cdot10^{-2}$ | $3.842\cdot10^{-2}$ |
| $L_1$ sum, $N = 512$, $x$ | $2.159\cdot10^{-2}$ (baseline + 2%) | $2.117\cdot10^{-2}$ | $2.117\cdot10^{-2}$ |
| $L_1$ sum decreasing from 256 to 512 | strict | pass | pass |
| $y$ and $z$ runs against $x$, rotated frame | bitwise (`DIR_TOL = 0`) | bitwise | bitwise |
| transverse copies and shared block columns | bitwise | bitwise | bitwise |
| GLM against no cleaning, 8 shared variables | bitwise, $\psi$ exactly 0 | bitwise | bitwise |

The baselines are the M2-P3b measurements, the same on both backends (`check.sh`). The $L_1$ decrease is first order:
the discontinuities dominate the error.

Related cases: [Sod](../euler/sod) is the Euler problem on the same grid; [Brio–Wu and Ryu–Jones 4d](./brio-wu-rj4d)
stress the eigensystem near its degeneracies; [MHD on `weno-riemann`](./weno-riemann) holds the HLL/HLLD measurements
of this problem; [EGLM cleaning](./eglm) runs the EV-3 leg of this script; [multi-realm runs](../euler/multi-realm)
split this problem in two realms.

## Findings

- The $y$ and $z$ runs equal the $x$ run bit for bit. Before M2-P3d they differed by $1.4\cdot10^{-12}$; the
  order-independent three-term sums (`mhd_sum3`) and the projections summed in the frame order of the direction made a
  cyclic rotation reproduce the $x$ arithmetic exactly.
- In 1-D the $(B_n, \psi)$ block of the splitting scheme is inert and block-diagonal, so GLM changes nothing. On
  `weno-riemann` the interpolated uniform $B_n$ carries a round-off seed into $\psi$, so `check.sh` relaxes the pair
  tolerance to $10^{-11}$ under `--numerics` (measured $\psi$ $9.2\cdot10^{-16}$ and shared variables
  $1.3\cdot10^{-12}$, HLLD primitive). The same tolerance applies to the `si` WENO weights (`WEIGHTS_EXE`, issue #49),
  whose smoothness indicator of a uniform field is round-off, not zero.
- The splitting scheme and HLLD with characteristic interpolation are indistinguishable on the figure; their $L_1$ sums
  are in [MHD on `weno-riemann`](./weno-riemann).
