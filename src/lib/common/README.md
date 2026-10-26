# ADAM Common Library — `src/lib/common`

The common library is the **portable, physics-agnostic core** of the ADAM SDK. The GPU objects (`fnl`; the unbuilt `nvf`, `gmp`) are separate types, initialised from these objects passed as arguments; every solver application (PRISM and FLUME are the maintained ones) composes them through `realm_object`. No GPU-specific code lives here — the entire directory compiles with any standard Fortran 2003+ compiler.

**Documentation:** <https://szaghi.github.io/adam/concepts/library>

## Module summary

| Module | Category | File |
|--------|----------|------|
| `adam_parameters` | Infrastructure | `adam_parameters.f90` |
| `adam_tree_node_object` | AMR | `adam_tree_node_object.f90` |
| `adam_tree_bucket_object` | AMR | `adam_tree_bucket_object.f90` |
| `adam_tree_object` | AMR | `adam_tree_object.F90` |
| `adam_grid_object` | Grid/Field | `adam_grid_object.F90` |
| `adam_field_object` | Grid/Field | `adam_field_object.F90` |
| `adam_maps_object` | Communication | `adam_maps_object.F90` |
| `adam_weno_object` | Numerics — spatial | `adam_weno_object.F90` |
| `adam_fdv_operators_library` | Numerics — spatial | `adam_fdv_operators_library.F90` |
| `adam_riemann_euler_library` | Numerics — spatial | `adam_riemann_euler_library.F90` |
| `adam_rk_object` | Numerics — temporal | `adam_rk_object.F90` |
| `adam_leapfrog_object` | Numerics — temporal | `adam_leapfrog_object.F90` |
| `adam_blanes_moan_object` | Numerics — temporal | `adam_blanes_moan_object.F90` |
| `adam_cfm_object` | Numerics — temporal | `adam_cfm_object.F90` |
| `adam_eos_ic_object` | Physics | `adam_eos_ic_object.F90` |
| `adam_ib_object` | Physics | `adam_ib_object.F90` |
| `adam_amr_object` | Physics | `adam_amr_object.F90` |
| `adam_flail_object` | Physics | `adam_flail_object.F90` |
| `adam_flail_parameters` | Physics | `adam_flail_parameters.F90` |
| `adam_poisson_fft_object` | Physics | `adam_poisson_fft_object.F90` |
| `adam_poisson_mg_object` | Physics | `adam_poisson_mg_object.F90` |
| `adam_lattice_green_library` | Physics | `adam_lattice_green_library.F90` |
| `adam_mpih_object` | Infrastructure | `adam_mpih_object.F90` |
| `adam_io_object` | Infrastructure | `adam_io_object.F90` |
| `adam_slices_object` | Infrastructure | `adam_slices_object.F90` |
| `adam_adam_object` | Infrastructure | `adam_adam_object.F90` |
| `adam_realm_object` | Forest | `adam_realm_object.F90` |
| `adam_forest_object` | Forest | `adam_forest_object.F90` |
| `adam_forest_manifest` | Forest | `adam_forest_manifest.F90` |
| `adam_flux_register_object` | Forest | `adam_flux_register_object.F90` |
| `adam_common_library` | Entry point | `adam_common_library.F90` |
| `adam_mpih_global` | Singleton | `adam_mpih_global.F90` |
| `adam_globals` | Singleton aggregator | `adam_globals.F90` |

---

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
