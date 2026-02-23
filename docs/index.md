---
layout: home

hero:
  name: ADAM
  text: Adaptive Mesh Refinement with Immersed Boundary
  tagline: High-performance fluid dynamics solver for GPU computing — from laptop to exascale.
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: Applications
      link: /applications/nasto/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/adam

features:
  - icon: 🌲
    title: Adaptive Mesh Refinement
    details: Morton-order linearized octree/quadtree for efficient AMR with automatic refinement and coarsening driven by user-defined markers.
  - icon: 🕸️
    title: Immersed Boundary Method
    details: Handle complex and moving geometries without body-fitted meshes using accurate IB techniques with eikonal-based distance fields.
  - icon: ⚡
    title: Multi-backend GPU Computing
    details: CUDA Fortran (NVF), OpenACC (FNL), and OpenMP offloading (GMP) — a single source tree targeting exascale HPC from one to thousands of GPUs.
  - icon:
      src: /lgplv3-logo.png
    title: Free and Open Source
    details: Released under the GNU Lesser General Public License v3.0 (LGPLv3) — free to use, study, modify, and redistribute, including in proprietary software linking against the library.
    link: /license
    linkText: View license
  - icon: ⚛️
    title: Multiple Physics
    details: Compressible Navier-Stokes (NASTO) and Maxwell equations (PRISM) share the same AMR, IB, and I/O infrastructure.
  - icon: 🔧
    title: INI-driven Configuration
    details: Simulations are fully configured through human-readable INI files — no recompilation needed to change physics, numerics, or grid parameters.
---

## A modular SDK for CFD

### Agnostic physics equations

ADAM is designed as a **physics-agnostic SDK**: its core infrastructure — AMR, IB, WENO numerics, multi-backend GPU acceleration, and parallel I/O — is entirely decoupled from any specific set of governing equations. Solvers are built on top by composing these building blocks, so adding a new physics application requires only implementing the problem-specific terms while inheriting the full HPC stack for free. Applications already built on ADAM span a wide spectrum: from compressible Navier-Stokes flows and shock-driven phenomena to electromagnetic plasma simulations governed by Maxwell equations, with further models in development. The same INI-driven configuration system, the same restart and output formats, and the same GPU backends apply uniformly across all of them.

### GPU-accelerated enabled and backward support for traditional HPC CPU-based architectures

ADAM runs equally well on **GPU-accelerated clusters** — exploiting CUDA Fortran, OpenACC, or OpenMP offloading across thousands of GPUs interconnected via NVLink or InfiniBand — and on **traditional CPU-based HPC clusters**, where the CPU backend with MPI and OpenMP delivers full functionality without requiring any GPU hardware. The choice of backend is a compile-time switch; the application source and input files are identical in both environments.

## Authors

- Andrea di Mascio — [andrea.dimascio@univaq.it](mailto:andrea.dimascio@univaq.it)
- Federico Negro — [federico.negro.01@gmail.com](mailto:federico.negro.01@gmail.com)
- Giacomo Rossi — [giacomo.rossi@amd.it](mailto:giacomo.rossi@amd.it)
- Francesco Salvadore — [f.salvadore@cineca.it](mailto:f.salvadore@cineca.it)
- Stefano Zaghi — [stefano.zaghi@cnr.it](mailto:stefano.zaghi@cnr.it)

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

## Citing ADAM

If you use ADAM in work that leads to a scientific publication, please cite the following paper:

> S. Zaghi, F. Salvadore, A. Di Mascio, G. Rossi —
> **Efficient GPU parallelization of adaptive mesh refinement technique for high-order compressible solver with immersed boundary** —
> *Computers and Fluids*, 266 (2023) 106040.
> DOI: [10.1016/j.compfluid.2023.106040](https://doi.org/10.1016/j.compfluid.2023.106040)

The paper describes the ADAM framework architecture, the AMR/IB coupling strategy, the GPU parallelization approach (CUDA Fortran), and demonstrates strong scaling on a shock–sphere interaction benchmark. A preprint is available in [docs/papers/zaghi-2023-computer_fluids.pdf](../papers/zaghi-2023-computer_fluids.pdf).

BibTeX entry:

```bibtex
@article{zaghi2023adam,
  author  = {Zaghi, S. and Salvadore, F. and {Di Mascio}, A. and Rossi, G.},
  title   = {Efficient {GPU} parallelization of adaptive mesh refinement technique
             for high-order compressible solver with immersed boundary},
  journal = {Computers \& Fluids},
  volume  = {266},
  pages   = {106040},
  year    = {2023},
  doi     = {10.1016/j.compfluid.2023.106040},
}
```

## Quick start

### The core library is the shared foundation

`src/lib` provides the physics-agnostic building blocks. Every application in `src/app` composes the same set of core objects and adds only the physics-specific layer on top:

```fortran
! src/lib/common — shared by ALL applications
use adam_mpih_object   ! MPI handler
use adam_grid_object   ! block-structured grid + AMR geometry
use adam_field_object  ! 5D field arrays  (nv, ni, nj, nk, nb)
use adam_amr_object    ! AMR refinement / coarsening
use adam_ib_object     ! Immersed Boundary
use adam_weno_object   ! WENO reconstructor (orders 3–11)
use adam_rk_object     ! Runge-Kutta integrator
```

### NASTO — compressible Navier-Stokes (`src/app/nasto`)

NASTO composes the core objects and adds the NS physics layer:

```fortran
type :: nasto_common_object
   !--- reused from src/lib (identical in every application) ---
   type(mpih_object)  :: mpih    ! MPI handler
   type(grid_object)  :: grid    ! AMR grid
   type(field_object) :: field   ! conservative variables q(nv,ni,nj,nk,nb)
   type(amr_object)   :: amr     ! refinement markers
   type(ib_object)    :: ib      ! solid bodies
   type(weno_object)  :: weno    ! spatial reconstruction
   type(rk_object)    :: rk      ! time integration
   !--- NS-specific layer (src/app/nasto/common) ---------------
   type(nasto_physics_object) :: physics  ! ideal-gas EOS, fluxes
   type(nasto_bc_object)      :: bc       ! inflow / wall / periodic BCs
   type(nasto_ic_object)      :: ic       ! shock-tube, vortex ICs
   type(nasto_io_object)      :: io       ! HDF5 output
   type(nasto_time_object)    :: time     ! CFL control
end type
```

### PRISM — Maxwell equations + PIC (`src/app/prism`)

PRISM reuses the **same** core objects and replaces the physics layer with electromagnetics:

```fortran
type :: prism_common_object
   !--- reused from src/lib (identical to NASTO) ---------------
   type(mpih_object)  :: mpih
   type(grid_object)  :: grid
   type(field_object) :: field   ! EM field  E, B, J  (nv,ni,nj,nk,nb)
   type(amr_object)   :: amr
   type(ib_object)    :: ib
   type(weno_object)  :: weno
   type(rk_object)    :: rk
   !--- EM-specific layer (src/app/prism/common) ---------------
   type(prism_physics_object)            :: physics   ! Maxwell curl operators
   type(prism_coil_object)               :: coil      ! current sources
   type(prism_pic_object)                :: pic       ! Particle-in-Cell
   type(prism_external_fields_object)    :: external_fields
   type(prism_bc_object)                 :: bc
   type(prism_ic_object)                 :: ic
   type(prism_io_object)                 :: io
   type(prism_time_object)               :: time
end type
```

### Adding a new solver follows the same pattern

A new application needs only a physics-specific layer; the entire HPC stack is inherited:

```fortran
type :: my_solver_object
   !--- reused from src/lib — zero extra work ---
   type(mpih_object)  :: mpih
   type(grid_object)  :: grid
   type(field_object) :: field
   type(amr_object)   :: amr
   type(ib_object)    :: ib
   type(weno_object)  :: weno
   type(rk_object)    :: rk
   !--- only this part is new ---
   type(my_physics_object) :: physics
   type(my_bc_object)      :: bc
   type(my_ic_object)      :: ic
   type(my_io_object)      :: io
end type
```

### Run any application with an INI file

All solvers share the same INI-driven configuration and build workflow:

```ini
[numerics]
scheme_space = weno
fdv_order    = 5

[runge_kutta]
scheme = runge-kutta-ssp-54

[time]
time_max = 0.2
CFL      = 0.4

[grid]
ni = 128 ; nj = 1 ; nk = 1 ; ngc = 3
emin_x = 0.0 ; emax_x = 1.0

[initial_conditions]
type = sod-x
```

```bash
# Navier-Stokes on GPU
FoBiS.py build -mode nasto-nvf-cuda && mpirun -np 4 exe/adam_nasto_nvf

# Maxwell / plasma on GPU
FoBiS.py build -mode prism-fnl-nvf-oac && mpirun -np 4 exe/adam_prism_fnl

# Any solver on CPU-only cluster — same source, same input file
FoBiS.py build -mode nasto-cpu-gnu && mpirun -np 64 exe/adam_nasto_cpu
```

Results are written as HDF5 files, readable with ParaView or any HDF5 tool. For complete test cases see the [Tests](/tests/) section.
