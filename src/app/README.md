# Applications

ADAM ships five solver applications, each built on the same SDK layer of core objects
(`adam_grid_object`, `adam_field_object`, `adam_weno_object`, etc.).
Applications share physics-agnostic infrastructure — AMR, ghost-cell exchange, I/O, IB —
and specialise only in the equations being solved and the numerical methods they require.

## At a glance

| Application | Equations solved | Physics domain | Backends | Status |
|-------------|-----------------|----------------|----------|--------|
| [NASTO](#nasto) | Compressible Navier-Stokes | Turbulent compressible CFD | CPU · NVF · FNL · GMP | Production |
| [PRISM](#prism) | Maxwell + PIC | Plasma/electromagnetics | CPU · FNL | Development |
| [CHASE](#chase) | Euler (inviscid) | Inviscid compressible flow | CPU | Experimental |
| [PATCH](#patch) | Poisson (elliptic) | Potential / pressure fields | CPU | Research |
| [ASCOT](#ascot) | — (utility) | Post-processing | — | Complete |

## Common design pattern

Every application follows the same directory layout:

```
src/app/<name>/
├── common/        # Shared physics, BC, IC, IO — backend-independent
├── cpu/           # CPU-only: MPI + OpenMP
├── nvf/           # CUDA Fortran (NVIDIA HPC SDK)
├── fnl/           # OpenACC via FUNDAL (device-agnostic GPU)
└── gmp/           # OpenMP target offloading (experimental)
```

Each `<name>_<backend>_object` extends the `<name>_common_object`, specialising
only the compute kernels while inheriting configuration, I/O, and physics logic.
Backends that are not yet developed for a given application simply do not have the
corresponding subdirectory.

### Configuration

All applications are configured through a single INI file (parsed by
[FiNeR](https://github.com/szaghi/FiNeR)).  Typical sections:

| Section | Purpose |
|---------|---------|
| `[IO]` | Output basename, save frequency, restart options |
| `[time]` | Maximum simulated time/iterations, CFL number |
| `[schemes]` | Temporal and spatial numerical schemes |
| `[grid]` | Domain bounds, cell counts `(ni, nj, nk)`, ghost-cell width `ngc` |
| `[physics]` | Species count, thermodynamic/electromagnetic constants |
| `[amr]` | Refinement levels, pruning thresholds, marker definitions |
| `[bc_*]` | Boundary conditions for each domain face |
| `[solids]` | Immersed boundary geometry (spheres, STL meshes) |
| `[slices]` | Optional slice-output sampling locations |

---

## NASTO

> ADAM for compressible Navier-Stokes equations — turbulent and viscous flows.

NASTO is the flagship ADAM application for compressible, viscous CFD.
It solves the three-dimensional compressible Navier-Stokes conservation equations
for ideal gases with full support for AMR, IB, and multi-GPU parallelism.

### Equations

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla \cdot \mathbf{F}^c(\mathbf{U})
= \nabla \cdot \mathbf{F}^d(\mathbf{U}) + \mathbf{S}
$$

where $\mathbf{U} = (\rho, \rho\mathbf{v}, E)^T$ are the conserved variables
(density, momentum, total energy), $\mathbf{F}^c$ the convective flux,
$\mathbf{F}^d$ the diffusive flux (viscosity + heat conduction), and $\mathbf{S}$
body-force source terms.

### Key features

- **High-order WENO**: upwind reconstructions of orders 3, 5, 7, 9, 11 (configurable per run)
- **Immersed boundary**: ghost-cell IB for arbitrary solid geometries
- **Runge-Kutta**: multi-stage explicit schemes (SSP-RK2, SSP-RK3, low-storage RK4)
- **Boundary conditions**: supersonic inflow · extrapolation outflow · solid wall · periodic
- **Initial conditions**: two-region Riemann · vortex advection · user-defined
- **GPU acceleration**: production-ready NVF (CUDA Fortran) and FNL (OpenACC) backends
- **Equation of state**: calorically perfect ideal gas; multi-species support in development

### Backends

| Backend | Subdirectory | Parallelism | Status |
|---------|-------------|-------------|--------|
| CPU | `cpu/` | MPI + OpenMP | Production |
| NVF | `nvf/` | MPI + CUDA Fortran | Production |
| FNL | `fnl/` | MPI + OpenACC (FUNDAL) | Production |
| GMP | `gmp/` | MPI + OpenMP target | Development |

### Source layout

```
src/app/nasto/
├── common/
│   ├── adam_nasto_common_object.F90   # Base class; reads INI, owns grid/field
│   ├── adam_nasto_physics_object.F90  # EOS, viscosity, heat-conduction closures
│   ├── adam_nasto_bc_object.F90       # Boundary conditions
│   ├── adam_nasto_ic_object.F90       # Initial conditions
│   ├── adam_nasto_io_object.F90       # HDF5 output, restart, slices
│   ├── adam_nasto_time_object.F90     # CFL-based time-step control
│   ├── adam_nasto_eos_object.F90      # Ideal-gas EOS
│   ├── adam_nasto_schemes_object.F90  # Scheme selection (WENO order, RK stages)
│   └── adam_nasto_common_library.F90  # Utility procedures
├── cpu/
│   ├── adam_nasto_cpu.F90             # Main program (CPU)
│   └── adam_nasto_cpu_object.F90      # CPU backend object
├── nvf/
│   ├── adam_nasto_nvf.F90             # Main program (CUDA Fortran)
│   ├── adam_nasto_nvf_object.F90      # NVF backend object
│   ├── adam_nasto_nvf_kernels.F90     # Device kernels (RK, BC, IB)
│   └── adam_nasto_nvf_cns_kernels.F90 # CNS flux kernels (CUDA device)
├── fnl/
│   ├── adam_nasto_fnl.F90             # Main program (OpenACC)
│   ├── adam_nasto_fnl_object.F90      # FNL backend object
│   ├── adam_nasto_fnl_kernels.F90     # OpenACC kernels (RK, BC, IB)
│   ├── adam_nasto_fnl_cns_kernels.F90 # CNS flux kernels (OpenACC)
│   └── adam_nasto_fnl_library.F90     # FNL-specific utilities
└── gmp/
    ├── adam_nasto_gmp.F90             # Main program (OpenMP target)
    ├── adam_nasto_gmp_object.F90      # GMP backend object
    ├── adam_nasto_gmp_kernels.F90     # OpenMP target kernels
    └── adam_nasto_gmp_cns_kernels.F90 # CNS flux kernels (OpenMP target)
```

### Build

```bash
# CPU (GNU compiler)
FoBiS.py build -mode nasto-cpu-gnu

# CUDA Fortran (NVIDIA HPC SDK)
FoBiS.py build -mode nasto-nvf-cuda

# OpenACC (NVIDIA HPC SDK)
FoBiS.py build -mode nasto-fnl-nvf-oac

# Debug build
FoBiS.py build -mode nasto-cpu-gnu-debug
```

---

## PRISM

> ADAM for Maxwell equations with Particle-In-Cell plasma dynamics.

PRISM solves the three-dimensional Maxwell equations coupled with a
Particle-In-Cell (PIC) framework for kinetic plasma simulations.
It is the most physics-rich ADAM application, featuring multiple temporal
integrators tailored for electromagnetic and plasma problems.

### Equations

$$
\frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \mathbf{E}, \qquad
\frac{\partial \mathbf{E}}{\partial t} = \frac{1}{\varepsilon_0}\left(
  \nabla \times \mathbf{B}/\mu_0 - \mathbf{J}
\right)
$$

where $\mathbf{E}$ is the electric field, $\mathbf{B}$ the magnetic field,
and $\mathbf{J}$ the current density computed from particle trajectories.

### Key features

- **Maxwell equations**: full electromagnetic field evolution on AMR structured grids
- **Particle-In-Cell**: kinetic plasma dynamics with particle injection and tracking
- **Coil modelling**: user-defined electromagnetic coil sources
- **External fields**: prescribed background electromagnetic fields
- **Multiple integrators**:
  - Runge-Kutta (explicit, multi-stage)
  - Leapfrog (symplectic, time-reversible)
  - Blanes-Moan (high-order symplectic)
  - Commutator-Free Magnus (CFM) for time-varying Hamiltonians
- **fWLayer**: forward-backward layer boundary treatment
- **FLAIL**: linear algebra solver for implicit field equations

### Physical constants (hard-coded)

| Constant | Symbol | Value |
|----------|--------|-------|
| Vacuum permeability | $\mu_0$ | $1.256637 \times 10^{-6}$ H/m |
| Vacuum permittivity | $\varepsilon_0$ | $8.854187 \times 10^{-12}$ F/m |
| Speed of light | $c$ | $1/\sqrt{\mu_0 \varepsilon_0}$ |
| Electron charge | $q_e$ | $1.602176 \times 10^{-19}$ C |
| Electron mass | $m_e$ | $9.109383 \times 10^{-31}$ kg |

### Backends

| Backend | Subdirectory | Parallelism | Status |
|---------|-------------|-------------|--------|
| CPU | `cpu/` | MPI + OpenMP | Production |
| FNL | `fnl/` | MPI + OpenACC (FUNDAL) | Development |

### Source layout

```
src/app/prism/
├── common/
│   ├── adam_prism_common_object.F90            # Base class
│   ├── adam_prism_physics_object.F90           # EM field properties
│   ├── adam_prism_pic_object.F90               # PIC handler
│   ├── adam_prism_leapfrog_pic_object.F90      # Leapfrog PIC integrator
│   ├── adam_prism_rk_pic_object.F90            # RK PIC integrator
│   ├── adam_prism_coil_object.F90              # Electromagnetic coil sources
│   ├── adam_prism_external_fields_object.F90   # Background field prescription
│   ├── adam_prism_fWLayer_object.F90           # Forward-backward layer BC
│   ├── adam_prism_particle_injection_object.F90 # Particle generation
│   ├── adam_prism_numerics_object.F90          # Scheme parameters
│   ├── adam_prism_rk_bc_object.F90             # RK boundary conditions
│   ├── adam_prism_bc_object.F90                # Boundary conditions
│   ├── adam_prism_ic_object.F90                # Initial conditions
│   ├── adam_prism_io_object.F90                # HDF5 output, restart
│   ├── adam_prism_time_object.F90              # Time-step control
│   ├── adam_prism_riemann_library.F90          # Riemann fluxes (EM)
│   └── adam_prism_common_library.F90           # Utility procedures
├── cpu/
│   ├── adam_prism_cpu.F90                      # Main program (CPU)
│   └── adam_prism_cpu_object.F90               # CPU backend object
└── fnl/
    ├── adam_prism_fnl.F90                      # Main program (OpenACC)
    ├── adam_prism_fnl_object.F90               # FNL backend object
    ├── adam_prism_fnl_kernels.F90              # OpenACC kernels
    ├── adam_prism_fnl_coil_object.F90          # GPU coil sources
    ├── adam_prism_fnl_external_fields_kernels.F90 # GPU external fields
    ├── adam_prism_fnl_fWLayer_object.F90       # GPU fWLayer BC
    └── adam_prism_fnl_library.F90              # FNL utilities
```

### Build

```bash
# CPU (GNU compiler)
FoBiS.py build -mode prism-gnu

# OpenACC (NVIDIA HPC SDK)
FoBiS.py build -mode prism-fnl-nvf-oac
```

---

## CHASE

> ADAM for Euler equations — inviscid compressible flow.

CHASE solves the three-dimensional compressible Euler equations (Navier-Stokes
without viscosity or heat conduction).  It is a simplified solver sharing the
same AMR, IB, and WENO infrastructure as NASTO, making it useful as a
low-cost testbed for new numerical schemes and AMR strategies.

### Equations

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla \cdot \mathbf{F}^c(\mathbf{U}) = 0
$$

Same conserved variables $\mathbf{U} = (\rho, \rho\mathbf{v}, E)^T$ as NASTO,
but the diffusive flux $\mathbf{F}^d \equiv 0$ (inviscid limit).

### Key features

- Compressible Euler equations (inviscid, ideal gas)
- WENO high-order spatial reconstructions
- Riemann solver-based convective fluxes
- Immersed boundary method
- Adaptive mesh refinement
- Same INI configuration format as NASTO

### Backends

| Backend | Subdirectory | Parallelism | Status |
|---------|-------------|-------------|--------|
| CPU | `cpu/` | MPI + OpenMP | Experimental |

### Source layout

```
src/app/chase/
├── common/
│   ├── adam_chase_common_object.F90   # Base class
│   ├── adam_chase_physics_object.F90  # Fluid thermodynamics
│   ├── adam_chase_bc_object.F90       # Boundary conditions
│   ├── adam_chase_ic_object.F90       # Initial conditions
│   ├── adam_chase_io_object.F90       # I/O handling
│   ├── adam_chase_time_object.F90     # Time-step control
│   ├── adam_chase_riemann_library.F90 # Euler Riemann solvers
│   └── adam_chase_common_library.F90  # Utility procedures
└── cpu/
    ├── adam_chase_cpu.F90             # Main program (CPU)
    └── adam_chase_cpu_object.F90      # CPU backend object
```

### Build

```bash
FoBiS.py build -mode chase-gnu
```

---

## PATCH

> ADAM for the Poisson equation — elliptic potential-field solver.

PATCH solves the three-dimensional Poisson (elliptic) equation on AMR structured
grids.  It is the smallest ADAM application and serves primarily as a research
vehicle for elliptic solvers and as the pressure-Poisson step in future
incompressible flow extensions.

### Equations

$$
\nabla^2 \phi = \rho
$$

where $\phi$ is the potential field and $\rho$ the source term (charge density,
forcing, etc.).

### Key features

- Elliptic Poisson solver on AMR structured grids
- FLAIL linear algebra solver integration
- Single scalar field storage (minimal memory footprint)
- Immersed boundary support
- Same INI configuration format as other ADAM applications

### Backends

| Backend | Subdirectory | Parallelism | Status |
|---------|-------------|-------------|--------|
| CPU | `cpu/` | MPI + OpenMP | Research |

### Source layout

```
src/app/patch/
├── common/
│   ├── adam_patch_common_object.F90   # Base class
│   ├── adam_patch_bc_object.F90       # Boundary conditions
│   ├── adam_patch_ic_object.F90       # Initial conditions
│   ├── adam_patch_io_object.F90       # I/O handling
│   ├── adam_patch_time_object.F90     # Iteration control
│   └── adam_patch_common_library.F90  # Utility procedures
└── cpu/
    ├── adam_patch_cpu.F90             # Main program (CPU)
    └── adam_patch_cpu_object.F90      # CPU backend object
```

### Build

```bash
FoBiS.py build -mode patch-gnu
```

---

## ASCOT

> ADAM Slices Converter — binary-to-ASCII post-processing utility.

ASCOT is a standalone utility that converts ADAM slice binary output files
into human-readable ASCII (Tecplot-compatible) format.  It has no dependency
on the ADAM core libraries and requires no parallel runtime.

### Usage

```bash
ascot -i <input_binary_slice.bin> [-o <output.dat>] [-v "VARIABLES = rho u v w p"]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Input binary slice file (required) | — |
| `-o` | Output ASCII file | `slice.dat` |
| `-v` | Tecplot variable header string | (none) |

### Source layout

```
src/app/ascot/
└── ascot.F90    # Standalone main program; no ADAM objects
```

### Build

```bash
FoBiS.py build -mode ascot-nvf-cuda
```
