# PRISM Common

The `common/` directory contains all backend-independent modules shared by the CPU and FNL backends of PRISM. Every type described here is aggregated into `prism_common_object`, which serves as the base type for both backends.

**Documentation:** <https://szaghi.github.io/adam/applications/prism/>

## Source Layout

| File | Module | Description |
|------|--------|-------------|
| `adam_prism_parameters.F90` | `adam_prism_parameters` | Physical constants, eigenvalue and eigenvector matrices |
| `adam_prism_physics_object.F90` | `adam_prism_physics_object` | Physical model selection, variable counts, eigenvector dispatch |
| `adam_prism_numerics_object.F90` | `adam_prism_numerics_object` | Temporal/spatial scheme selection, divergence correction config |
| `adam_prism_common_object.F90` | `adam_prism_common_object` | Base aggregate type for all backends |
| `adam_prism_bc_object.F90` | `adam_prism_bc_object` | Boundary conditions handler |
| `adam_prism_ic_object.F90` | `adam_prism_ic_object` | Initial conditions handler |
| `adam_prism_io_object.F90` | `adam_prism_io_object` | I/O configuration and output routines |
| `adam_prism_time_object.F90` | `adam_prism_time_object` | Time integration control |
| `adam_prism_coil_object.F90` | `adam_prism_coil_object` | Current-carrying coil current density sources |
| `adam_prism_fWLayer_object.F90` | `adam_prism_fWLayer_object` | Far-wave absorbing boundary layer |
| `adam_prism_external_fields_object.F90` | `adam_prism_external_fields_object` | Externally prescribed electromagnetic fields |
| `adam_prism_rk_bc_object.F90` | `adam_prism_rk_bc_object` | RK sub-step integrator for boundary conditions |
| `adam_prism_pic_object.F90` | `adam_prism_pic_object` | Particle-in-Cell configuration and particle-grid weighting |
| `adam_prism_particle_injection_object.F90` | `adam_prism_particle_injection_object` | Initial particle space and velocity distribution |
| `adam_prism_leapfrog_pic_object.F90` | `adam_prism_leapfrog_pic_object` | Leapfrog (Boris) PIC particle time integrator |
| `adam_prism_rk_pic_object.F90` | `adam_prism_rk_pic_object` | Runge-Kutta PIC particle time integrator |
| `adam_prism_pic_bc_object.F90` | `adam_prism_pic_bc_object` | Particle closure at the domain faces (periodic or legacy) |
| `adam_prism_pml_object.F90` | `adam_prism_pml_object` | Perfectly matched layer (`[PML]`) |
| `adam_prism_rk_pml_object.F90` | `adam_prism_rk_pml_object` | RK support for the PML auxiliary variables |
| `adam_prism_absorbing_layer_geometry.F90` | `adam_prism_absorbing_layer_geometry` | Layer geometry shared by fWLayer and PML |
| `adam_prism_probes_object.F90` | `adam_prism_probes_object` | Field probes (`[probes]`): time histories at fixed points |
| `adam_prism_grms_object.F90` | `adam_prism_grms_object` | Grms diagnostic (`[grms]`) |
| `adam_prism_magnetic_field_at_center_domain_object.F90` | `adam_prism_magnetic_field_at_center_domain_object` | Field-at-domain-centre diagnostic |
| `adam_prism_riemann_library.F90` | `adam_prism_riemann_library` | Maxwell Riemann solvers (LLF, HLL) and convective flux routines |
| `adam_prism_common_library.F90` | `adam_prism_common_library` | Barrel re-export of all common modules |

---

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
