# FLUME Common

The `common/` directory contains the backend-independent modules shared by the CPU and FNL backends of FLUME. Following
the PRISM design, every type defined here is aggregated into `flume_common_object`, the base type both backends extend.

> **Status: development** (milestone M1, compressible Euler).

## Modules

| File | Content |
|------|---------|
| `adam_flume_parameters.F90` | Variable indexes (conservative `IQ_*`, auxiliary `IA_*`), accepted option values, `strip_control` |
| `adam_flume_physics_object.F90` | `[physics]`: physical model (`euler`), `cp`, `cv`, variable counts |
| `adam_flume_numerics_object.F90` | `[numerics]`: `scheme_space`, `reconstruction_variables`, `reflux` |
| `adam_flume_euler_library.F90` | Pointwise Euler physics shared by host and device: conservative/primitive conversions, fluxes, Roe eigenvectors, flux splitting |
| `adam_flume_bc_object.F90` | `[bc_*]`: boundary condition types and inflow states |
| `adam_flume_ic_object.F90` | `[initial_conditions]`: uniform (seeded perturbation), isentropic vortex, Riemann regions; init-time AMR passes |
| `adam_flume_time_object.F90` | `[time]`: CFL, iteration and time limits |
| `adam_flume_diagnostics_object.F90` | `[diagnostics]`: conservation history |
| `adam_flume_common_object.F90` | `flume_common_object`: initialization, AMR markers (box, gradient, solid), seam flux accumulation for reflux, immersed-boundary spacing, restart, slices, auxiliary fields |
| `adam_flume_common_library.F90` | Barrel re-export of all common modules |

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
