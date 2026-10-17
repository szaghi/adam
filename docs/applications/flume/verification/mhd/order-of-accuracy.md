# Order of accuracy: linear waves, CPAW, magnetised vortex

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| MV-5, MV-6, MV-7 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) | CPU, FNL | yes | — |

## Case and scope

Three smooth problems whose exact solution after one period (or one crossing) is the initial state, so every
smooth-flow ingredient (eigen-decomposition, WENO5 characteristic reconstruction, fluxes, periodic ghost fill, time
integration) is measured against an exact answer:

- **MV-5, linear waves** (Stone et al. 2008, section 8.2): a small-amplitude eigenmode per family (fast, Alfvén, slow,
  entropy), 1-D along $x$ and 2-D inclined. It catches a wrong eigenvector or a family-specific loss of order.
- **MV-6, circularly polarised Alfvén wave** (Tóth 2000): an exact nonlinear solution at any amplitude. At
  $B_\perp = 0.1$ it exercises the WENO nonlinear weights, which the $10^{-7}$ linear waves leave at their linear
  values. Left and right polarisations are mirror images under $z \to -z$, so their errors must be equal **bitwise**
  (rounding is symmetric under a sign change).
- **MV-7, magnetised vortex** (Balsara 2004): a steady 2-D equilibrium of pressure, centrifugal force, magnetic tension
  and magnetic pressure convected by the free stream, with every variable varying in both directions and a
  divergence-free field. It checks the $L_1$ and $L_\infty$ order and that $\max\lvert\nabla\cdot\mathbf{B}\rvert$ stays
  at truncation level.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 5/3$ | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables (all three); `weno-riemann` HLLD under `--numerics` (MV-5, MV-6) | `mhd_face_split_fluxes`, `mhd_eigenvectors`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [numerics](../../numerics#weno-reconstruction) |
| Divergence control | MV-5, MV-6: `none` (MV-5 also GLM or EGLM with `--divergence-control`, $c_h = 1.1(\lvert u_n\rvert + c_f)$, no damping); MV-7: GLM, $c_h = 3$, `glm_alpha = 0.18`, `glm_damping_length = 1.0` | [models](../../models#divergence-control-mixed-glm) |
| div(B) operator (MV-7) | second-order centred difference (`fdv_order = 2`) | `compute_divb_history`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | SSP Runge–Kutta, 5 stages order 4 (`runge-kutta-ssp-54`); CFL 0.1 (MV-5, MV-6: the spatial error dominates), 0.4 (MV-7: errors equal to 4–5 digits at 0.1) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90); [numerics](../../numerics#time-integration) |
| Boundaries | periodic | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Leg | Domain and grid | Initial condition | Final time |
|---|---|---|---|
| MV-5, 1-D | $[0, 1]$ along $x$, the Sod grid of [V1](../euler/sod) (4 blocks along $x$), $N = 16, 32$ cells per wavelength (the figure adds 64) | [`mhd-linear-wave`](../../initial-conditions#mhd-linear-wave-small-amplitude-mhd-eigenmode): $\rho = 1$, $p = 3/5$, $\mathbf{B} = (1, \sqrt2, 1/2)$ in the wave frame ($c_f = 2$, $c_a = 1$, $c_s = 1/2$), $u_n = 1$ for the entropy wave, amplitude $10^{-7}$, wavelength 1 | one period, $1/\lvert u_n + c\rvert$ |
| MV-5, 2-D | $[0, \sqrt5] \times [0, \sqrt5/2]$, the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $N \times N/2$ square cells, $N = 32, 64$; wave at $\tan a = 2$ ($a = 63.43°$), one wavelength along each axis | the same | one period |
| MV-6 | the MV-5 2-D grid, $N = 32, 64, 128$ | [`mhd-cpaw`](../../initial-conditions#mhd-cpaw-circularly-polarised-alfven-wave): $\rho = 1$, $p = 0.1$, $B_\parallel = 1$, $v_\parallel = 0$, $B_\perp = 0.1$, wavelength 1 (Alfvén speed and period 1); right and left polarisations | $t = 1$ |
| MV-7 | $[-7, 7]^2$ (tails $e^{-24}$ at the box edge), the V2 grid, $N^2$ cells, $N = 64, 128, 256$ | [`mhd-vortex`](../../initial-conditions#mhd-vortex-magnetised-vortex): $\rho = 1$, $p = 1$, $\mathbf{v} = (1, 1, 0)$, $\mathbf{B} = 0$ background, $\kappa = \mu = 1$, radius 1 | one crossing, $t = 14$ |

All legs: 2 ranks, periodic on the active axes.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/linear-wave/make_linear_wave.py` | writes the MV-5 input (`--wave`, `--geometry 1d\|2d`, `--cells`, `--amplitude`, `--divergence-control`) |
| `src/tests/flume/verification/mhd/linear-wave/check.sh`, `linear_wave_oracle.py` | MV-5: the error $\epsilon$ of Stone et al. (last against first checkpoint), the bounds and the order (`--numerics`, `--divergence-control`, `--legs LIST`: `mv5-1d`, `mv5-2d`) |
| `src/tests/flume/verification/mhd/cpaw/make_cpaw.py` | writes the MV-6 input (`--polarisation`, `--geometry`, `--cells`) |
| `src/tests/flume/verification/mhd/cpaw/check.sh`, `cpaw_oracle.py` | MV-6: bounds, order, left = right (`--numerics`, `--legs mv6`) |
| `src/tests/flume/verification/mhd/vortex/make_mhd_vortex.py` | writes the MV-7 input (`--cells`, `--half-width`, `--cfl`) |
| `src/tests/flume/verification/mhd/vortex/check.sh`, `mhd_vortex_oracle.py` | MV-7: $L_1$, $L_\infty$, orders, div(B), initial checkpoint against the exact vortex (`--legs mv7`) |
| `src/tests/flume/verification/make_doc_runs.sh` group `order` | the 1-D ladders 16/32/64 on the two schemes (24 runs) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_order` | draws `linear-wave-order.png` (`--only order`) |

```bash
cd src/tests/flume/verification/mhd
linear-wave/check.sh                       # MV-5
cpaw/check.sh                              # MV-6
vortex/check.sh                            # MV-7 (about 40 min on the CPU, 31 of them at N = 256)
linear-wave/check.sh --numerics hlld       # MV-5 on weno-riemann
```

## Results

![MHD linear waves, order of accuracy](/flume/linear-wave-order.png)

*The 1-D linear waves on 16, 32 and 64 cells per wavelength, splitting scheme and `weno-riemann` HLLD; the legend
gives the orders of the two pairs.*

| Test | Bound | CPU | FNL |
|---|---|---|---|
| MV-5, every family, observed order | $\ge 4.7$ (design 5 minus 0.3); 4.8 under `hlld:primitive`/`hlld:characteristic` | about 5: 4.96–4.97 from 16 to 32 cells (1-D); 4.96–5.02 on the asserted ladders | pass |
| MV-5, error of every run | per run, CPU baseline + 2% (e.g. 1-D fast $1.512\cdot10^{-10}$ / $4.843\cdot10^{-12}$ at 16 / 32) | below | pass |
| MV-6, observed order on 32/64/128 | $\ge 4.7$ (4.8 under the HLLD specs) | 4.89 and 4.99 | pass |
| MV-6, $\epsilon$ at $N = 32 / 64 / 128$ | $1.049\cdot10^{-3}$ / $3.537\cdot10^{-5}$ / $1.112\cdot10^{-6}$ | below | pass |
| MV-6, left against right polarisation | bitwise (`LR_TOL = 0`) | bitwise | bitwise |
| MV-7, $L_1$ at $N = 64 / 128 / 256$ | $4.327\cdot10^{-4}$ / $2.224\cdot10^{-5}$ / $6.900\cdot10^{-7}$ | below | pass |
| MV-7, $L_\infty$ at $N = 64 / 128 / 256$ | $7.812\cdot10^{-3}$ / $5.194\cdot10^{-4}$ / $2.439\cdot10^{-5}$ | below | pass |
| MV-7, $L_1$ order (finest pair) | $\ge 4.7$ | 5.01 | pass |
| MV-7, $L_\infty$ order (finest pair) | $\ge 4.0$ | 4.41 (pre-asymptotic: 3.41, 3.91, 4.41 on 32/64, 64/128, 128/256) | pass |
| MV-7, $\max\lvert\nabla\cdot\mathbf{B}\rvert$ convergence | order $\ge 1.7$ (second-order operator) | truncation level | pass |
| MV-7, initial checkpoint against the exact vortex | $10^{-14}$ | pass | pass |

FNL ran every ladder (work directories of both backends); the tables record the CPU values, and the backends agree to
the printed digits.

Related cases: [isentropic vortex](../euler/isentropic-vortex) is the Euler order test on the same grid;
[MHD on `weno-riemann`](./weno-riemann) holds the orders of the HLLD runs; [EGLM cleaning](./eglm) repeats MV-5 under
EGLM; [order across seams](../seams/order-across-seams) measures the order across 2:1 seams.

## Findings

- Every linear-wave family converges at 4.96–4.97 from 16 to 32 cells. At 64 cells the error is $10^{-13}$, six
  decades below the amplitude, and the round-off floor lowers the slow and entropy orders to 4.85: MV-5 asserts on the
  16/32 pair (1-D) and on 32/64 (2-D) for that reason.
- The amplitude is $10^{-7}$, not the $10^{-6}$ of Stone et al.: at $10^{-6}$ the $O(A^2)$ nonlinear error of the mode
  ($\epsilon$ scales as $A^2$ at $N = 128$: $3.26\cdot10^{-10}$ at $A = 10^{-5}$, $3.26\cdot10^{-12}$ at $10^{-6}$)
  floors a fifth-order scheme already at $N = 64$–128, and below about $10^{-8}$ the accumulated round-off floors it
  ($1.5\cdot10^{-14}$ at $A = 10^{-8}$). At 2-D $N = 128$ the slow and Alfvén waves reach the floor (orders 2.6 and 4.3
  on 64/128).
- With GLM or EGLM on MV-5, $c_h$ sets the time step, so the errors differ from the `none` baseline (2-D slow wave
  2.4×).
- MV-7 runs GLM because with `none` the discrete div(B) of the vortex grows exponentially at a rate proportional to
  $1/h$: the run stops on a negative pressure at $t = 13.4$, 8.2, 4.3 for $N = 64$, 128, 256, whatever the CFL.
- The MV-7 $L_\infty$ order is still pre-asymptotic at $N = 256$ ($N = 512$ costs 2.5 h on the CPU); 32 cells is
  pre-asymptotic in $L_1$ too (order 3.65), so the ladder starts at 64.
