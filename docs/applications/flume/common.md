# FLUME Common

The `common/` directory contains the backend-independent modules shared by the CPU and FNL backends of FLUME. Following
the PRISM design, every type defined here is aggregated into `flume_common_object`, the base type both backends extend.
The pointwise libraries (`*_library.F90`) are `pure`, take explicit-size arguments and are tagged `!$acc routine seq` +
`!$omp declare target`, so the CPU loops and the GPU kernels call the same source.

> **Status: development** (milestones M1 to M13 complete; the list is in `docs/applications/flume/index.md`).

## Modules

| File | Content |
|------|---------|
| `adam_flume_parameters.F90` | Variable indexes (conservative `IQ_*`, auxiliary `IA_*`), model ids, accepted option values, `strip_control` |
| `adam_flume_physics_object.F90` | `[physics]` and `[eos]`: physical model (`euler`, `euler-de`, `mhd-ideal`), passive scalars, `gamma` or `cp`, `cv`, the equation of state (ideal, isothermal, barotropic), variable counts; primitive-to-conservative dispatch |
| `adam_flume_reference_object.F90` | `[reference]`: dimensional input converted to code units on the host (issue #49) |
| `adam_flume_gravity_object.F90` | `[gravity]`: external gravity (issue #73) and the self-gravity keys, the Poisson solver and its boundaries (issue #81) |
| `adam_flume_eos_tables_library.F90` | Barotropic equation-of-state tables, continuous piecewise polytropes (issue #100) |
| `adam_flume_dissipation_object.F90` | `[physics]` dissipative coefficients (issue #65): viscosity, conductivity, resistivity, each as a coefficient or its number, and the power-law viscosity |
| `adam_flume_cooling_object.F90` | `[cooling]` radiative cooling (issue #93): the node table, the particle mass, the time-step limiter, the precomputed exponents and evolution function, the radiated energy |
| `adam_flume_cooling_library.F90` | Townsend's exact integration shared by both backends: the cooling function, the evolution function and its inverse, the cooled internal energy of a cell |
| `adam_flume_sinks_object.F90` | `[sinks]` sink particles (issue #94): the replicated list (masses, positions, velocities, accelerations), the restart file and the sinks history; the deposit and the gather live in the common object |
| `adam_flume_dissipation_library.F90` | Pointwise dissipative physics shared by both backends: the viscous and heat-conduction flux of a direction, the Ohmic flux (MHD), the diffusivity of the time-step limit |
| `adam_flume_mhd_object.F90` | `[mhd]`: divergence control (`glm`, `eglm`, `none`), GLM speed, damping and check, div(B) monitor, positivity floors |
| `adam_flume_numerics_object.F90` | `[numerics]`: `scheme_space` (`weno`, `weno-riemann`), `reconstruction_variables`, Riemann solver, flux correction and sensor, `reflux` and its policy and cadence, `positivity_limiter`, `dissipative_order` |
| `adam_flume_euler_library.F90` | Pointwise Euler physics: conversions, fluxes, Roe average and eigenvectors, flux splitting, face states, LLF/HLL/HLLC Riemann solvers; its body `adam_flume_euler_library_body.INC` is shared with the dual-energy library |
| `adam_flume_euler_de_library.F90`, `adam_flume_euler_ps_library.F90`, `adam_flume_euler_de_ps_library.F90` | The Euler physics with dual energy (issue #91), passive scalars (issue #92), or both |
| `adam_flume_euler_iso_library.F90`, `adam_flume_euler_baro_library.F90` | Isothermal and barotropic Euler physics (issue #100), one body `adam_flume_euler_iso_library_body.INC` behind the EOS hooks |
| `adam_flume_mhd_iso_library.F90`, `adam_flume_mhd_baro_library.F90` | Isothermal and barotropic ideal MHD physics and Riemann solvers (issue #100), one body `adam_flume_mhd_iso_library_body.INC` |
| `adam_flume_mhd_library.F90` | Pointwise ideal MHD physics: conversions, fluxes, fast speed, Roe–Balsara eigensystem (with and without GLM), flux splitting |
| `adam_flume_mhd_riemann_library.F90` | MHD face states and Riemann solvers (LLF, HLL, HLLD; exact GLM subsystem) of `weno-riemann` |
| `adam_flume_bc_object.F90` | `[bc_*]`: boundary condition types (extrapolation, inflow, inviscid wall, no-slip and isothermal walls, periodic, hydrostatic), inflow states, wall velocities and temperatures; `wall_noslip_ghost` (the no-slip ghost state); `realm_edge_face`, `realm_edge_donor` (the face condition and the donor of a ghost at a realm edge or corner, shared by the CPU and FNL backends) |
| `adam_flume_ic_object.F90` | `[initial_conditions]`: 25 initial conditions (uniform, isentropic vortex, Riemann regions, Shu–Osher, rotated Riemann, GLM pulse, div(B) peak, MHD linear wave, CPAW, magnetised vortex, Orszag–Tang, rotor, field loop, linear, polynomial, sine wave, Couette, Becker shock, isothermal and polytropic atmospheres, Rayleigh–Taylor, polytrope, singular isothermal sphere, Bonnor–Ebert, Boss–Bodenheimer); init-time AMR passes |
| `adam_flume_seam_sync_object.F90` | Per-stage seam flux synchronisation data of the positivity limiter at 2:1 seams (issue #50) |
| `adam_flume_time_object.F90` | `[time]`: CFL, iteration and time limits |
| `adam_flume_diagnostics_object.F90` | `[diagnostics]`: conservation (and div(B)) history cadence |
| `adam_flume_common_object.F90` | `flume_common_object`: initialization and cross-section checks, AMR markers (box, gradient, Löhner, solid, Jeans) and runtime regrids, gravity and self-gravity sources (the solve by the library's FLAIL), sink particles' deposit and gather, seam flux accumulation for reflux, immersed-boundary spacing, restart, slices, auxiliary and MHD derived fields |
| `adam_flume_common_library.F90` | Barrel re-export of all common modules |

The user documentation (models, numerics, input reference, verification gallery) is in `docs/applications/flume/`.
