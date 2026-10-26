# PRISM

**PRISM** — **P**lasma **R**esearch us**I**ng **S**imulation **M**ethods — is an electromagnetic and plasma simulation application built on the ADAM framework. It solves the Maxwell equations for computational electromagnetics and, optionally, couples them to charged-particle dynamics via a Particle-in-Cell (PIC) method.

**Documentation:** <https://szaghi.github.io/adam/applications/prism/>

## Source Layout

```
src/app/prism/
├── common/                                      # Shared across all backends
│   ├── adam_prism_parameters.F90                # Physical constants and eigenvector matrices
│   ├── adam_prism_physics_object.F90            # Physical model, variable counts, eigenvectors
│   ├── adam_prism_numerics_object.F90           # Temporal and spatial scheme selection
│   ├── adam_prism_common_object.F90             # prism_common_object base type (t0 elliptic solves on FLAIL)
│   ├── adam_prism_bc_object.F90                 # Boundary conditions handler
│   ├── adam_prism_ic_object.F90                 # Initial conditions handler
│   ├── adam_prism_rk_bc_object.F90              # RK sub-stepping for boundary conditions
│   ├── adam_prism_io_object.F90                 # I/O handler
│   ├── adam_prism_time_object.F90               # Time integration handler
│   ├── adam_prism_coil_object.F90               # Current-carrying coil sources
│   ├── adam_prism_fWLayer_object.F90            # Far-wave layer (absorbing BC layer)
│   ├── adam_prism_pml_object.F90                # Perfectly matched layer
│   ├── adam_prism_rk_pml_object.F90             # RK support for the PML auxiliary variables
│   ├── adam_prism_absorbing_layer_geometry.F90  # Layer geometry shared by fWLayer and PML
│   ├── adam_prism_external_fields_object.F90    # Externally prescribed fields
│   ├── adam_prism_pic_object.F90                # Particle-in-Cell configuration
│   ├── adam_prism_pic_bc_object.F90             # Particle closure at the domain faces
│   ├── adam_prism_particle_injection_object.F90 # Particle injection handler
│   ├── adam_prism_leapfrog_pic_object.F90       # Leapfrog PIC time integrator
│   ├── adam_prism_rk_pic_object.F90             # Runge-Kutta PIC time integrator
│   ├── adam_prism_probes_object.F90             # Field probes: time histories at fixed points
│   ├── adam_prism_grms_object.F90               # Grms diagnostic
│   ├── adam_prism_magnetic_field_at_center_domain_object.F90 # Field-at-domain-centre diagnostic
│   ├── adam_prism_riemann_library.F90           # Maxwell Riemann solvers (LLF, HLL) and flux routines
│   └── adam_prism_common_library.F90            # Barrel re-export of all common modules
├── cpu/                                         # CPU backend
│   ├── adam_prism_cpu.F90                       # Entry point (program)
│   └── adam_prism_cpu_object.F90                # prism_cpu_object implementation
└── fnl/                                         # OpenACC GPU backend
    ├── adam_prism_fnl.F90                       # Entry point (program)
    ├── adam_prism_fnl_object.F90                # prism_fnl_object implementation
    ├── adam_prism_fnl_external_fields_kernels.F90 # OpenACC external-fields kernels
    ├── adam_prism_fnl_coil_object.F90           # GPU coil handler
    ├── adam_prism_fnl_fWLayer_object.F90        # GPU far-wave layer handler
    ├── adam_prism_fnl_pml_object.F90            # GPU PML
    ├── adam_prism_fnl_rk_pml_object.F90         # GPU RK support for the PML auxiliary variables
    ├── adam_prism_fnl_pic_object.F90            # GPU Particle-in-Cell
    ├── adam_prism_fnl_pic_conserving_object.F90 # GPU charge-conserving currents
    ├── adam_prism_fnl_leapfrog_pic_object.F90   # GPU PIC leapfrog integrator
    ├── adam_prism_fnl_rk_pic_object.F90         # GPU PIC Runge-Kutta integrator
    ├── adam_prism_fnl_probes_object.F90         # GPU field probes
    └── adam_prism_fnl_library.F90               # Barrel re-export of all FNL modules
```

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
