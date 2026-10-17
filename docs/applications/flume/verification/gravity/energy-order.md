# The order of the energy source

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GV-9 (G2a) | External gravity | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes | — |

## Case and scope

The fluxes of a finite-difference scheme are point values of the primitive function $h$ of $\rho\mathbf u$
($\rho u = \frac1{\Delta x}\int h$ over the cell), not of $\rho\mathbf u$. G1's energy source
$\tfrac12 g(F^\rho_{i-1/2} + F^\rho_{i+1/2})$ is then $g\rho u_i + O(\Delta x^2)$, and it took the whole scheme to
second order on any smooth flow under gravity. None of the G1 legs measured an order: free fall is uniform, and the
conservation legs are exact at any order.

G2a adds the correction $\tfrac1{24} g(-F^\rho_{i-3/2} + F^\rho_{i-1/2} + F^\rho_{i+1/2} - F^\rho_{i+3/2})$, which
makes the source the cell value of $h$ to fourth order. The case that measures it is a density wave falling freely in
a periodic box: the momentum source is exact on it, the energy source is not (its $(\rho u)''$ is not zero), and $p$
is uniform, so any pressure error is the energy source's. A second-order energy source fails the asserted order of
3.5.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#external-gravity) |
| Energy source | $\sum_d g_d(-F^\rho_{i-3/2} + 13F^\rho_{i-1/2} + 13F^\rho_{i+1/2} - F^\rho_{i+3/2})/24$: the face form plus the fourth-order correction, from the high-order face mass fluxes of the stage, before any blending | `add_gravity_sources`, `compute_gravity_correction`, `compute_gravity_correction_faces`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_gravity_correction_faces_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |
| Block-end faces | the faces $i \pm 3/2$ of the block-end cells lie beyond the block: two edge kernels, one per face family, compute them from the ghosts | `add_gravity_edge_correction`, [`adam_flume_cpu_face_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_face_kernels_agnostic.INC) (`weno`), [`adam_flume_cpu_riemann_face_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_riemann_face_kernels_agnostic.INC) (`weno-riemann`); `add_gravity_edge_correction_dev`, [`adam_flume_fnl_face_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_face_kernels_agnostic.INC) |
| Ghost cells | one more layer than the WENO stencil, `ngc >= S + 1` (with `weno-riemann` `max(S, 2) + 1`), fatal otherwise | `check_ngc_number`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Limiter | the correction is a high-order cell source, scaled by $\Lambda$ with the rest of the face form | `add_gravity_limited_correction`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Space, `weno-riemann` | WENO5 interpolation, HLLC, sixth-order flux correction (`flux_correction = 6th`, sensor `weno`) | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_euler_hllc_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_euler_hllc_kernels.F90); [numerics](../../numerics#scheme-space-weno-riemann-interpolation-riemann-flux-and-correction) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL $0.4\,(32/N)^{2/3}$ so the RK3 error stays below the spatial one | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | every face `periodic` | [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Case | `wave2d` |
| Domain | $[0, 1]^2$, null z |
| Grid | quadtree, 4×4 blocks of $n^2$ cells, $n = 8, 16, 32$ ($N = 32, 64, 128$), two uniform levels, `ngc = 4` |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), `wave_mode = entropy`, `wave_angle = 45`, `wave_amplitude = 0.2`, `wavelength` $= 1/\sqrt2$: $\rho = 1 + 0.2\sin(2\pi(x + y))$ |
| State | $\mathbf u_0 = (1, 0.5)$, $p = 1$ |
| Gravity | $\mathbf g = (-0.6, -0.8)$ |
| Exact solution | the wave advected by $\mathbf u_0 t + \mathbf g t^2/2$, $\mathbf u = \mathbf u_0 + \mathbf g t$, $p$ uniform |
| Boundaries | every face `periodic` |
| Final time | $t = 0.2$ |
| Asserted | the $L_1$ orders of $\rho$ and $p$ on the finest pair at least 3.5, on `weno` and on `weno-riemann` |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes `wave2d` (`--cells`, `--cfl`, `--time-max`, `--riemann hllc`) |
| `src/tests/flume/verification/gravity/check.sh` | runs the ladder on both schemes with the force legs (GV-9, `--legs gv9`) |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | `wave`: the $L_1$ errors against the exact falling wave and their orders (`--min-order 3.5`) |
| `src/tests/flume/verification/make_doc_runs.sh` (group `gravity`) | re-runs `wave-{weno,riemann}-{8,16,32}` |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_gravity_order`) | draws the figure |

```bash
cd src/tests/flume/verification/gravity
./check.sh --np 2 --cases iso2d                                          # CPU (GV-9 runs with the force legs)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np 2 --cases iso2d   # FNL
```

## Results

![GV-9: L1 errors of the falling density wave against N, weno and weno-riemann](/flume/gravity-order.png)

*$L_1$ errors of $\rho$ and $p$ of the falling density wave against $N$, `weno` and `weno-riemann` (HLLC).*

The two schemes lie on top of each other, $\rho$ at slope 5 and $p$ at 4.6 on the finest pair. The pressure error is
the energy source's: on this wave the momentum source is exact and $p$ is uniform, so any pressure error comes from the
energy equation. The control is not drawn (it needs a special build); its errors are in the table.

CPU and FNL agree to every digit shown:

| Scheme | $L_1$ $p/p_1$ at $N = 32, 64, 128$ | Orders ($\rho$, $p$), finest pair: bound | CPU | FNL |
|---|---|---|---|---|
| G1 face form alone (`weno`, negative control) | $1.1\cdot10^{-5}$, $2.8\cdot10^{-6}$, $7.0\cdot10^{-7}$ | — | 2.12, 2.00 | the same |
| `weno` | $3.3\cdot10^{-7}$, $1.2\cdot10^{-8}$, $5.2\cdot10^{-10}$ | $\ge 3.5$ | 5.01, 4.57 | the same |
| `weno-riemann` (HLLC) | $3.4\cdot10^{-7}$, $1.3\cdot10^{-8}$, $5.3\cdot10^{-10}$ | $\ge 3.5$ | 5.00, 4.59 | the same |

At $N = 128$ the pressure error falls 1300-fold. The negative control is a build with the correction switched off,
not a run-time option; it is not part of `check.sh`.

The other GV legs re-run with the correction (CPU and FNL) are tabulated in [the gravity force](./force): GV-1
unchanged, GV-3 conserving $W$ to $\le 1.2\cdot10^{-14}$, `box2d-lim` with 107 984 limited faces, GV-6 and NV-5
bitwise, GV-7 $\le 7.9\cdot10^{-14}$.

Related cases: [the gravity force](./force) (the face form and the closed boxes); the second-order potential path of
self-gravity runs the same wave in [SV-1](../self-gravity/potential).

## Findings

- The face form alone is second order on a finite-difference scheme; the fourth-order correction restores the order
  of the scheme (5 in $\rho$, 4.6 in $p$) on both face families.
- **Conservation.** Each face's weights still sum to one over the cells, so $E + \rho\Phi$ stays conserved exactly
  with periodic boundaries and across block interfaces.
- **The invariant at reflecting walls (owner decision: $W$ is the invariant).** At a wall the correction of the wall
  cell reads the face beyond it: the mirrored ghosts give $F_{-1} = -F_1$ and $F_0 = 0$, and the wall row's mass
  balance $F_1 = -\Delta x\,d\rho_1/dt$ turns the left-over $g\Delta x(F_1 - F_{-1})/24$ into
  $d/dt(-g\Delta x\rho_1 V/12)$. The exactly conserved total energy is then
  $W = \sum(E + \rho\Phi)V + \sum_{\text{lower wall}} g_n\Delta x_n\rho V/12 - \sum_{\text{upper wall}} g_n\Delta x_n
  \rho V/12$, an $O(\Delta x^2)$ wall term away from the plain sum. Dropping the correction flux at walls would keep
  the plain sum but cost the wall cells an $O(\Delta x)$ energy source error. GV-3 asserts $W$ and reports the plain
  sum ([the gravity force](./force)).
- **One more ghost layer.** The faces $i \pm 3/2$ of the block-end cells lie beyond the block, so gravity takes
  `ngc >= S + 1`, fatal otherwise; every gravity input has `ngc = 4`.
- **With the positivity limiter** the correction is a high-order cell source, scaled by $\Lambda$ with the rest of the
  face form; `box2d-lim` limits more faces than at G1 (107 984 against 74 304) for that reason.
