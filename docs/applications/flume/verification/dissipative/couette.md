# Compressible Couette flow

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-3 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | `couette` |

## Case and scope

A channel between a fixed isothermal wall at $y = 0$ and an adiabatic no-slip wall moving at $U$ along $x$ at
$y = H$. The exact steady solution has uniform $p$, the linear velocity $u = U y/H$ and the viscous-heating parabola
$T = T_w + \mu U^2/(kH)\,(y - y^2/(2H))$. It is an exact steady solution of the compressible Navier–Stokes equations
with constant $\mu$ and $k$, so the run starts from it and the departure at the final time is the error of the scheme,
with no relaxation transient.

What the case catches: the first test of the no-slip walls with the viscous fluxes, and of a nonlinear energy flux
(the viscous heating $u\,\tau$): a wrong wall ghost, a wrong wall velocity or temperature, or a wrong viscous work term
shows as a departure from the exact profile.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Navier–Stokes, $\gamma = 1.4$, $R = 1$, $\mu = 0.01$, $Pr = 0.75$ ($k = \mu c_p/Pr$) | [models](../../models#compressible-euler); `load_from_file`, [`adam_flume_dissipation_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_object.F90) |
| Dissipative fluxes | `dissipative_order = 4` | `add_dissipative_fluxes_o4`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Wall ghosts | velocity mirrored about the wall's (exact for the linear $u$); the moving adiabatic wall mirrors $T$ symmetrically (exact for the symmetric $T$); the isothermal wall sets the geometric mirror $T_g = T_w^2/T$ and $\rho_g = p/(R T_g)$ (second order) | `wall_noslip_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls) |
| FNL backend | the same fluxes and wall ghosts on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1] \times [0, H]$, $H = 1$, $z$ null |
| Grid | quadtree, $4 \times 4$ blocks of $6 \times N/4$ cells: $x$ active and thin (a null $x$ would freeze $\rho u$), $N$ = 32 / 64 / 128 cells across the channel |
| Ranks | 2 |
| Initial condition | [`couette`](../../initial-conditions#couette-exact-compressible-couette-flow-verification): the exact profile, $p = 1$, $U = 1$, $T_w = 1$ |
| `bc_x_min`, `bc_x_max` | `periodic` |
| `bc_y_min` ($y = 0$) | `wall-isothermal`, `wall_temperature = 1` (at rest) |
| `bc_y_max` ($y = 1$) | `wall-noslip`, `wall_u = 1` (adiabatic, moving) |
| `bc_z_min`, `bc_z_max` | `extrapolation` |
| Final time | $t = 0.2$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/profiles.py` | `make --case couette` writes each input; `oracle` measures the RMS of $u/U$ and $T/T_w$ against the exact profile |
| `src/tests/flume/verification/viscous/check.sh --legs vv3` | runs the ladder (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/regression/couette/input.ini` | the frozen case, $t = 0.05$ |
| `make_doc_figures.py` (`fig_couette`, `fig_viscous`) | the figures, from the CPU work directories of `viscous/check.sh` |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --legs vv3             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs vv3       # FNL
```

## Results

![VV-3: the Couette channel, its fields, the profiles against the exact ones and where the error lives](/flume/couette.png)

*Top: $u/U$ and $T/T_w$ at 64 cells and $t = 0.2$, and the 32-cell profiles on the exact ones; bottom: the pointwise
errors of $u$ and $T$ at 64 cells, and the temperature error averaged along $x$ at every resolution.*

![VV-1 to VV-5: convergence, and Becker's shock against the exact profile](/flume/viscous.png)

*The Couette temperature ladder is the "VV-3 Couette (T)" curve of the left panel.*

Oracle: order $\ge 1.9$ on $u/U$ and $T/T_w$ (finest pair).

| Quantity | N | RMS error (CPU) | Order | FNL |
|---|---|---|---|---|
| $u/U$ | 32 / 64 / 128 | $1.71\cdot10^{-7}$ / $3.84\cdot10^{-8}$ / $8.72\cdot10^{-9}$ | +2.15, +2.14 | the same; the last digit at 128 cells differs |
| $T/T_w$ | 32 / 64 / 128 | $3.01\cdot10^{-6}$ / $6.62\cdot10^{-7}$ / $1.49\cdot10^{-7}$ | +2.19, +2.15 | the same |

Both backends ran on 2 ranks.

Related cases: the wall-free viscous waves are [VV-1 and VV-2](./viscous-waves); the other exact steady profile is
[Becker's shock](./becker-shock).

## Findings

- The run is second order overall: the moving-wall mirror is exact for the linear $u$ and the symmetric $T$, the
  isothermal geometric mirror $T_w^2/T$ is second order.
- The fields show where that second order comes from. Velocity and temperature are uniform along the periodic $x$, and
  the profiles lie on the exact ones (the linear $u$, the viscous-heating parabola of $T$). The error enters at the
  isothermal wall and diffuses into the channel: at $t = 0.2$ its front has reached about $y = 0.5$, and beyond it the
  solution is exact to round-off at 128 cells. The moving adiabatic wall adds nothing.
- With the linear isothermal mirror of P1 ($2T_w - T$) the errors were 18% lower, at the same orders. The geometric
  mirror was adopted because the linear one made negative ghost densities beside hot gas; the two differ by
  $O(\Delta x^2)$.
