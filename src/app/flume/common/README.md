# FLUME Common

The `common/` directory contains the backend-independent modules shared by the CPU and FNL backends of FLUME. Following
the PRISM design, every type defined here is aggregated into `flume_common_object`, the base type both backends extend.
The pointwise libraries (`*_library.F90`) are `pure`, take explicit-size arguments and are tagged `!$acc routine seq` +
`!$omp declare target`, so the CPU loops and the GPU kernels call the same source.

> **Status: development** (milestones M1 Euler and M2 ideal MHD complete; M3, issue #47, in progress).

## Modules

| File | Content |
|------|---------|
| `adam_flume_parameters.F90` | Variable indexes (conservative `IQ_*`, auxiliary `IA_*`), model ids, accepted option values, `strip_control` |
| `adam_flume_physics_object.F90` | `[physics]`: physical model (`euler`, `mhd-ideal`), `cp`, `cv`, variable counts; primitive-to-conservative dispatch |
| `adam_flume_dissipation_object.F90` | `[physics]` dissipative coefficients (issue #65): viscosity, conductivity, resistivity, each as a coefficient or its number, and the power-law viscosity |
| `adam_flume_dissipation_library.F90` | Pointwise dissipative physics shared by both backends: the viscous and heat-conduction flux of a direction, the Ohmic flux (MHD), the diffusivity of the time-step limit |
| `adam_flume_mhd_object.F90` | `[mhd]`: divergence control (`glm`, `none`), GLM speed, damping and check, div(B) monitor, positivity floors |
| `adam_flume_numerics_object.F90` | `[numerics]`: `scheme_space` (`weno`, `weno-riemann`), `reconstruction_variables`, Riemann solver, flux correction and sensor, `reflux`, `positivity_limiter`, `dissipative_order` |
| `adam_flume_euler_library.F90` | Pointwise Euler physics: conversions, fluxes, Roe average and eigenvectors, flux splitting, face states, LLF/HLL/HLLC Riemann solvers |
| `adam_flume_mhd_library.F90` | Pointwise ideal MHD physics: conversions, fluxes, fast speed, Roe–Balsara eigensystem (with and without GLM), flux splitting |
| `adam_flume_mhd_riemann_library.F90` | MHD face states and Riemann solvers (LLF, HLL, HLLD; exact GLM subsystem) of `weno-riemann` |
| `adam_flume_bc_object.F90` | `[bc_*]`: boundary condition types (extrapolation, inflow, inviscid wall, no-slip and isothermal walls, periodic), inflow states, wall velocities and temperatures; `wall_noslip_ghost` (the no-slip ghost state); `realm_edge_face`, `realm_edge_donor` (the face condition and the donor of a ghost at a realm edge or corner, shared by the CPU and FNL backends) |
| `adam_flume_ic_object.F90` | `[initial_conditions]`: 17 initial conditions (uniform, isentropic vortex, Riemann regions, Shu–Osher, rotated Riemann, GLM pulse, div(B) peak, MHD linear wave, CPAW, magnetised vortex, Orszag–Tang, rotor, field loop, linear, sine wave, Couette, Becker shock); init-time AMR passes |
| `adam_flume_time_object.F90` | `[time]`: CFL, iteration and time limits |
| `adam_flume_diagnostics_object.F90` | `[diagnostics]`: conservation (and div(B)) history cadence |
| `adam_flume_common_object.F90` | `flume_common_object`: initialization and cross-section checks, AMR markers (box, gradient, solid), seam flux accumulation for reflux, immersed-boundary spacing, restart, slices, auxiliary and MHD derived fields |
| `adam_flume_common_library.F90` | Barrel re-export of all common modules |

The user documentation (models, numerics, input reference, verification gallery) is in `docs/applications/flume/`.

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
