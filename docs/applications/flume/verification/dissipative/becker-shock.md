# Becker's viscous shock

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-4 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | `becker-shock` |

## Case and scope

Becker's (1922) exact profile of a shock resolved by the viscosity, with constant $\mu$ and $Pr = 3/4$, in the shock
frame. It is an exact steady solution of the compressible Navier–Stokes equations, so, like [Couette](./couette), the
run starts from the exact solution and the departure at the final time is the error of the scheme, with no fit of a
shock position.

What the case catches: the strongest nonlinear viscous flux of the suite, with the energy flux carrying both the
viscous work and the heat flux. A wrong coefficient, a wrong viscous work or a wrong heat flux moves the profile.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Navier–Stokes, $\gamma = 1.4$, $R = 1$, constant $\mu$, $Pr = 0.75$ | [models](../../models#compressible-euler) |
| Dissipative fluxes | `dissipative_order = 4` | `add_dissipative_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); `compute_dissipative_flux`, [`adam_flume_dissipation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_library.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Initial profile | Becker's implicit solution, solved cell by cell | `becker_shock`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| Boundaries | supersonic `inflow` upstream, `extrapolation` elsewhere | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#inflow-prescribed-state-dirichlet-on-the-full-conservative-vector) |
| FNL backend | the same kernels on the device | `add_dissipative_fluxes_o4_dev`, [`adam_flume_fnl_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_dissipation_kernels.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $x \in [-0.25, 0.25]$, the shock at $x = 0$; $y$ and $z$ null |
| Grid | quadtree, $4 \times 4$ blocks of $N/4 \times 1$ cells, $N$ = 64 / 128 / 256 along $x$ |
| Ranks | 2 |
| Initial condition | [`becker-shock`](../../initial-conditions#becker-shock-becker-s-exact-viscous-shock-verification): upstream state $\rho_1 = 1$, $p_1 = 1/\gamma$ (so $a_1 = 1$), $u_1 = M$ |
| `bc_x_min` | `inflow` with the upstream state |
| Other faces | `extrapolation` |
| Final time | $t = 0.05$ |

| Leg | Mach | $\mu$ | Cells across the shock at $N = 256$ |
|---|---|---|---|
| `becker-m2-o4` | 2 | 0.005 | about 8 |
| `becker-m3-o4` | 3 | 0.01 | about 8 |

At Mach 3 with $\mu = 0.005$ the shock spans 4 cells, WENO's nonlinear weights still act, and the measured order was
+2.25; the leg therefore takes $\mu = 0.01$ at Mach 3.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/profiles.py` | `make --case becker` writes each input; `oracle` measures $(u - u_\mathrm{exact})/(u_1 - u_2)$ |
| `src/tests/flume/verification/viscous/check.sh --legs vv4` | runs the two ladders (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/regression/becker-shock/input.ini` | the frozen case: Mach 2, 64 cells, $t = 0.01$ |
| `make_doc_figures.py` (`fig_viscous`) | the figure, from the CPU work directories of `viscous/check.sh` |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --legs vv4             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs vv4       # FNL
```

## Results

![VV-1 to VV-5: convergence, and Becker's shock against the exact profile](/flume/viscous.png)

*Right panel: $(u - u_2)/(u_1 - u_2)$ at 64 and 256 cells against Becker's profile, Mach 2 and 3; the left panel
holds the two Becker ladders.*

Oracle: order $\ge 3$ on the RMS and on the maximum of $(u - u_\mathrm{exact})/(u_1 - u_2)$ (finest pair).

| Mach | 64 | 128 | 256 | Orders (RMS) |
|---|---|---|---|---|
| 2 | RMS $1.05\cdot10^{-2}$, max $6.2\cdot10^{-2}$ | $1.92\cdot10^{-3}$, $1.15\cdot10^{-2}$ | $1.30\cdot10^{-4}$, $8.3\cdot10^{-4}$ | +2.45, +3.88 |
| 3 | RMS $1.47\cdot10^{-2}$, max $8.9\cdot10^{-2}$ | $3.09\cdot10^{-3}$, $1.89\cdot10^{-2}$ | $2.77\cdot10^{-4}$, $1.7\cdot10^{-3}$ | +2.25, +3.48 |

CPU, 2 ranks. Every leg passes on FNL (2 ranks) with the same numbers to the printed digits.

Related cases: the other exact steady viscous profile is [Couette](./couette); the backend agreement of these runs,
the largest of the suite, is in [CPU and FNL agree](./cpu-fnl-agreement).

## Findings

- The coarse pair converges at about second order and the finest pair above third: the shock width is fixed by $\mu$,
  so it spans about 2 cells at 64 cells and about 8 at 256. With 4 cells across it the order was still pre-asymptotic
  (+2.25).
- The Becker runs hold the largest CPU–FNL differences of the suite ($4.3\cdot10^{-13}$ of the momentum scale), being
  its steepest gradients.
