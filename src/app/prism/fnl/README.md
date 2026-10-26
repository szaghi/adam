# PRISM — FNL (OpenACC) Backend

The FNL backend (`prism-fnl`) is the GPU-accelerated implementation of PRISM using the
**FNL library** (Fortran aNd LAPACK, OpenACC-based GPU framework). It inherits the full PRISM
common infrastructure and extends it with GPU-resident arrays, OpenACC-parallelized FDV operators,
and GPU-aware MPI communication via the FUNDAL library.

**Documentation:** <https://szaghi.github.io/adam/applications/prism/backends>

## Source files

| File | Module | Description |
|------|--------|-------------|
| `adam_prism_fnl.F90` | entry point | Parses CLI argument, delegates to `prism%simulate` |
| `adam_prism_fnl_object.F90` | `adam_prism_fnl_object` | Main type, all operators and integrators |
| `adam_prism_fnl_library.F90` | `adam_prism_fnl_library` | Barrel re-export of all FNL PRISM modules |
| `adam_prism_fnl_coil_object.F90` | `adam_prism_fnl_coil_object` | GPU coil current source handler |
| `adam_prism_fnl_fWLayer_object.F90` | `adam_prism_fnl_fWLayer_object` | GPU far-wave layer handler |
| `adam_prism_fnl_external_fields_kernels.F90` | `adam_prism_fnl_external_fields_kernels` | GPU external field sub/add kernels |
| `adam_prism_fnl_pml_object.F90` | `adam_prism_fnl_pml_object` | GPU perfectly matched layer |
| `adam_prism_fnl_rk_pml_object.F90` | `adam_prism_fnl_rk_pml_object` | GPU RK support for the PML auxiliary variables |
| `adam_prism_fnl_pic_object.F90` | `adam_prism_fnl_pic_object` | GPU Particle-in-Cell (weightings, deposits, gathers) |
| `adam_prism_fnl_pic_conserving_object.F90` | `adam_prism_fnl_pic_conserving_object` | GPU charge-conserving currents (londrillo, totorica, esirkepov_Q<N>) |
| `adam_prism_fnl_leapfrog_pic_object.F90` | `adam_prism_fnl_leapfrog_pic_object` | GPU PIC leapfrog integrator |
| `adam_prism_fnl_rk_pic_object.F90` | `adam_prism_fnl_rk_pic_object` | GPU PIC Runge-Kutta integrator |
| `adam_prism_fnl_probes_object.F90` | `adam_prism_fnl_probes_object` | GPU field probes |

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
