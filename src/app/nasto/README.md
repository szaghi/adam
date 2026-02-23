# NASTO

> ADAM for compressible Navier-Stokes equations.

NASTO is the flagship application of the ADAM framework, solving the
three-dimensional compressible Navier-Stokes conservation equations for ideal
gases. It combines high-order WENO finite-difference schemes, Immersed Boundary
method, and Adaptive Mesh Refinement with full multi-GPU acceleration.

## Governing equations

$$
\frac{\partial \mathbf{U}}{\partial t}
+ \nabla \cdot \mathbf{F}^c(\mathbf{U})
= \nabla \cdot \mathbf{F}^d(\mathbf{U}) + \mathbf{S}
$$

| Symbol | Meaning |
|--------|---------|
| $\mathbf{U} = (\rho,\, \rho\mathbf{v},\, E)^T$ | Conserved variables: density, momentum, total energy |
| $\mathbf{F}^c$ | Convective (inviscid) flux |
| $\mathbf{F}^d$ | Diffusive flux — viscosity + heat conduction |
| $\mathbf{S}$ | Body-force source terms |

The equation of state closes the system for a calorically perfect ideal gas:
$p = (\gamma - 1)\left(E - \tfrac{1}{2}\rho |\mathbf{v}|^2\right)$.

## Key features

- **High-order WENO** — upwind reconstructions of orders 3, 5, 7, 9, 11 (configurable per run via INI)
- **Immersed boundary** — ghost-cell IB method for arbitrary solid geometries embedded in the Cartesian grid
- **Runge-Kutta** — explicit multi-stage schemes: SSP-RK2, SSP-RK3, low-storage RK4
- **Adaptive mesh refinement** — Morton-order octree AMR with configurable refinement levels and gradient-based markers
- **Boundary conditions**: supersonic inflow · extrapolation outflow · solid wall · periodic
- **Initial conditions**: two-region Riemann (shock-sphere interaction) · vortex advection · user-extensible
- **Parallel I/O** — HDF5 restart and field output; slice sampling for in-situ analysis
- **Multi-GPU** — production-ready NVF (CUDA Fortran) and FNL (OpenACC) backends; GMP (OpenMP target) in development

## Backends

| Backend | Subdirectory | Parallelism | Status |
|---------|-------------|-------------|--------|
| CPU | `cpu/` | MPI + OpenMP | Production |
| NVF | `nvf/` | MPI + CUDA Fortran | Production |
| FNL | `fnl/` | MPI + OpenACC (FUNDAL) | Production |
| GMP | `gmp/` | MPI + OpenMP target | Development |

Each backend object extends `adam_nasto_common_object`, specialising only the
compute kernels (flux evaluation, time-stepping, IB masking) while inheriting
configuration, physics, I/O, and AMR logic from `common/`.

## Source layout

```
src/app/nasto/
├── common/
│   ├── adam_nasto_common_object.F90    # Base class: reads INI, owns grid/field/AMR
│   ├── adam_nasto_physics_object.F90   # EOS, viscosity, heat-conduction closures
│   ├── adam_nasto_bc_object.F90        # Boundary condition dispatch
│   ├── adam_nasto_ic_object.F90        # Initial condition setup
│   ├── adam_nasto_io_object.F90        # HDF5 output, restart, slice sampling
│   ├── adam_nasto_time_object.F90      # CFL-based time-step control
│   ├── adam_nasto_eos_object.F90       # Ideal-gas equation of state
│   ├── adam_nasto_schemes_object.F90   # WENO order and RK stage selection
│   └── adam_nasto_common_library.F90   # Shared utility procedures
├── cpu/
│   ├── adam_nasto_cpu.F90              # Main program (CPU)
│   └── adam_nasto_cpu_object.F90       # CPU backend — MPI + OpenMP loops
├── nvf/
│   ├── adam_nasto_nvf.F90              # Main program (CUDA Fortran)
│   ├── adam_nasto_nvf_object.F90       # NVF backend object
│   ├── adam_nasto_nvf_kernels.F90      # Device kernels: RK stages, BC, IB
│   └── adam_nasto_nvf_cns_kernels.F90  # CNS flux kernels (`attributes(device)`)
├── fnl/
│   ├── adam_nasto_fnl.F90              # Main program (OpenACC)
│   ├── adam_nasto_fnl_object.F90       # FNL backend object
│   ├── adam_nasto_fnl_kernels.F90      # OpenACC kernels: RK stages, BC, IB
│   ├── adam_nasto_fnl_cns_kernels.F90  # CNS flux kernels (`!$acc routine seq`)
│   └── adam_nasto_fnl_library.F90      # FNL-specific utilities
└── gmp/
    ├── adam_nasto_gmp.F90              # Main program (OpenMP target)
    ├── adam_nasto_gmp_object.F90       # GMP backend object
    ├── adam_nasto_gmp_kernels.F90      # OpenMP target kernels: RK stages, BC, IB
    └── adam_nasto_gmp_cns_kernels.F90  # CNS flux kernels (`!$omp target`)
```

## Build

NASTO is built with [FoBiS](https://github.com/szaghi/FoBiS). All build modes
are defined in `fobos` at the repository root.

```bash
# List available NASTO build modes
FoBiS.py build -lmodes | grep nasto

# CPU — GNU compiler
FoBiS.py build -mode nasto-cpu-gnu

# CUDA Fortran — NVIDIA HPC SDK
FoBiS.py build -mode nasto-nvf-cuda

# OpenACC — NVIDIA HPC SDK
FoBiS.py build -mode nasto-fnl-nvf-oac

# Debug (bounds checking, traceback)
FoBiS.py build -mode nasto-cpu-gnu-debug
```

Executables are placed in `exe/` at the repository root.

## Configuration

A simulation is fully defined by a single INI file. The default search path is
`adam_nasto.ini` in the working directory; pass an alternative with
`-i <file>`.

Representative sections:

```ini
[IO]
  output_basename = nasto_run
  save_frequency  = 100        ! steps between HDF5 snapshots
  restart         = .false.

[time]
  t_max     = 10.0             ! maximum simulated time
  iter_max  = 100000           ! maximum iteration count
  cfl       = 0.8              ! Courant–Friedrichs–Lewy number

[schemes]
  temporal  = rk-ssp-3        ! SSP-RK3 time integration
  convective = weno-upwind     ! WENO convective flux
  diffusive  = central-4       ! 4th-order central diffusive flux

[schemes_weno_upwind]
  order = 3                    ! WENO order (3 = 5th-order stencil)

[grid]
  ni = 8  nj = 8  nk = 8      ! cells per block
  ngc = 4                      ! ghost-cell layers (≥ WENO half-stencil)
  xmin = 0.0  xmax = 20.0
  ymin = 0.0  ymax = 20.0
  zmin = 0.0  zmax = 20.0

[physics]
  n_species = 1
  gamma = 1.4
  mu    = 1.8e-5               ! dynamic viscosity (Pa·s)
  Pr    = 0.72                 ! Prandtl number

[amr]
  levels_max = 5
  marker     = gradient-density

[bc_xmin]
  type = supersonic-inflow
[bc_xmax]
  type = extrapolation

[solids]
  n_solids = 1
  solid_1_type   = sphere
  solid_1_centre = 10.0 10.0 10.0
  solid_1_radius = 1.0
```

## Running

```bash
# Single MPI rank (CPU)
mpirun -np 1 exe/adam_nasto_cpu

# Multi-GPU (one MPI rank per GPU)
mpirun -np 4 exe/adam_nasto_nvf

# With custom INI file
mpirun -np 4 exe/adam_nasto_fnl -i my_case.ini
```

## Tests

Integration test cases are in `src/tests/nasto/`:

| Test | Description | Reference |
|------|-------------|-----------|
| Sod-X | Sod shock tube along x-axis | Sod (1978) |
| Sod-Y | Sod shock tube along y-axis | Sod (1978) |
| Sod-Z | Sod shock tube along z-axis | Sod (1978) |
| Shock-sphere | Normal shock interacting with a rigid sphere | Zaghi et al. (2023) |

Build and run a test:

```bash
FoBiS.py build -mode test-nasto-sod-x-gnu
mpirun -np 1 exe/test_nasto_sod_x
```

## Copyrights

NASTO is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

## Citing ADAM

If you use NASTO in work that leads to a scientific publication, please cite:

> S. Zaghi, F. Salvadore, A. Di Mascio, G. Rossi —
> **Efficient GPU parallelization of adaptive mesh refinement technique for
> high-order compressible solver with immersed boundary** —
> *Computers and Fluids*, 266 (2023) 106040.
> DOI: [10.1016/j.compfluid.2023.106040](https://doi.org/10.1016/j.compfluid.2023.106040)

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
