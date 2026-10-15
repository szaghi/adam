# Isentropic vortex: order of accuracy

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| V2, RV-2 | Euler | M1 ([#35](https://github.com/szaghi/adam/issues/35)), M3 ([#47](https://github.com/szaghi/adam/issues/47)) | CPU (all schemes), FNL (`weno`, `weno-riemann` LLF and HLLC) | yes | `vortex-periodic` |

## Case and scope

An exact smooth solution of the Euler equations: an isentropic vortex (Shu 1998, ICASE 97-65, §5.1) convected
diagonally by the free stream through a doubly periodic box. Every run is compared cell by cell with the exact solution,
and the observed order between successive resolutions, $p = \log_2(E_N / E_{2N})$, must reach the design order.

What the case catches: a loss of design order on smooth flows, which [Sod](./sod) cannot see (WENO is first order at
discontinuities):
- a wrong stencil weight or reconstruction table;
- a wrong ghost fill across periodic block faces, across blocks and across ranks;
- a first-order boundary treatment.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Space, `weno` | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90) |
| Space, `weno-riemann` (RV-2) | WENO5 interpolation (characteristic), LLF, HLL or HLLC, the 6th-order flux correction with the WENO-weight sensor | `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `initialize_interpolation`, [`adam_weno_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_weno_object.F90) |
| Time | SSP Runge–Kutta, 5 stages, 4th order (`runge-kutta-ssp-54`), CFL 0.1 (small, so the space error dominates the time error) | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Initial condition | the `isentropic-vortex` kind | `isentropic_vortex`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| Boundaries | library periodicity in $x$ and $y$ (across blocks and ranks) | `update_ghost`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#periodic) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree (`ratio = 4`, `nk = 1`), 16 blocks (4 per side) of $16^2$, $32^2$, $64^2$ cells: $64^2$, $128^2$, $256^2$ cells, `ngc = 3` |
| Ranks | 2 (8/8 blocks) |
| Free stream | $(\rho, u, v, w, p) = (1, 1, 1, 0, 1)$ |
| Vortex | centre $(0.5, 0.5)$, radius 0.07 (tails $e^{-25}$ at the box edge), strength 5 |
| Initial condition | [`isentropic-vortex`](../../initial-conditions#isentropic-vortex) |
| Boundaries | `periodic` on $x$ and $y$; `extrapolation` on the null $z$ faces |
| Final time | $t = 0.1$ (the vortex centre moves by $(0.1, 0.1)$) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/vortex/vortex-n{064,128,256}.ini` | the ladder |
| `src/tests/flume/verification/vortex/check.sh` | V2: runs the ladder (`--np N`; `FLUME_EXE` selects the backend), order bound 4.5 |
| `src/tests/flume/verification/vortex/vortex_oracle.py` | the exact density, $L_1$ and $L_\infty$ errors, the finest-pair order, the null-$z$ copies identical |
| `src/tests/flume/verification/riemann-flux/check.sh --leg vortex` | RV-2: the ladder on `weno-riemann`, order bound 4.8 |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_vortex`) | the figure, from the work directories of the two scripts |

```bash
cd src/tests/flume/verification/vortex && ./check.sh                    # V2, CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh             # V2, FNL
cd ../riemann-flux && ./check.sh --leg vortex --solver hllc             # RV-2 (also llf, hll)
```

## Results

![Isentropic vortex](/flume/vortex.png)

*Left: density at $t = 0.1$, $128^2$ cells, `weno-riemann` HLLC. Right: $L_1(\rho)$ against the cells per side for the
four schemes, with the order-5 slope.*

| Scheme | Bound (finest-pair order) | $L_1(\rho)$, $256^2$ | Order 128 → 256 | FNL |
|---|---|---|---|---|
| `weno` (split) | 4.5 (V2) | $8.38\cdot10^{-8}$ | 5.79 | the same |
| `weno-riemann` LLF | 4.8 (RV-2) | $1.42\cdot10^{-7}$ | 5.64 | the same |
| `weno-riemann` HLL | 4.8 (RV-2) | $8.98\cdot10^{-8}$ | 6.26 | not run |
| `weno-riemann` HLLC | 4.8 (RV-2) | $8.04\cdot10^{-8}$ | 6.49 | the same |

The full ladders on the CPU:

| Scheme | $L_1(\rho)$, $64^2$ / $128^2$ / $256^2$ |
|---|---|
| `weno` (split) | $8.218\cdot10^{-5}$ / $4.621\cdot10^{-6}$ / $8.379\cdot10^{-8}$ |
| `weno-riemann` LLF, WENO-weight sensor | $1.182\cdot10^{-4}$ / $7.073\cdot10^{-6}$ / $1.415\cdot10^{-7}$ |
| `weno-riemann` LLF, sensor `none` | $9.113\cdot10^{-5}$ / $4.748\cdot10^{-6}$ / $1.415\cdot10^{-7}$ |

Related cases: [Sod](./sod) checks the same schemes at discontinuities; the MHD orders are in
[Order of accuracy](../mhd/order-of-accuracy).

## Findings

- The order is asserted on the finest pair only. From $64^2$ to $128^2$ it is lower (4.15 for the splitting scheme):
  at $64^2$ the vortex radius spans 4.5 cells and the WENO weights are still far from the linear ones.
- At $256^2$ the WENO-weight sensor of `weno-riemann` is inactive (the same digits with and without it); at $128^2$ it
  trips at the vortex extrema and costs a factor 1.5. The library weights use the exponent of $S = 3$, the NumPy
  prototype 2.
- LLF is 1.7× more dissipative than the splitting scheme on the vortex; HLL and HLLC close the gap, and HLLC is the
  most accurate of the four.
- The quadtree layout is bitwise identical to the octree one with a null $z$ axis, and 16 times cheaper.
- The earlier verification page and the figure's panel title called the run "one period"; the input runs to $t = 0.1$,
  one tenth of the period of the diagonal convection through the unit box (the title is corrected since #85 R1).
