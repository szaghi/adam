# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

See also:
- [Architecture Reference](.claude/docs/architecture.md) - Class hierarchy, singletons, key type members, library aggregation
- [Fortran Style Guide](.claude/docs/fortran-style.md) - Zen of Fortran, coding style, modern patterns and examples
- [HPC Optimization Guide](.claude/docs/hpc-optimization.md) - GPU, MPI, memory optimization, benchmarking, pitfalls
- [Agent Workflows](.claude/docs/agent-workflows.md) - Task protocols, debugging checklists

## Development Philosophy

**Target**: Exascale HPC applications with strong scaling to thousands of GPU nodes
**Priorities**: Correctness first, then performance, maintainability, and portability
**Evidence-based optimization**: Always profile before optimizing, benchmark after changes
**Documentation standard**: Well-referenced explanations citing authoritative sources over speculation

## Project Overview

ADAM (Adaptive Mesh Refinement with Immersed Boundary) is a high-performance fluid dynamics solver for GPU computing. It solves compressible Navier-Stokes and Maxwell equations using:
- Adaptive Mesh Refinement (AMR) with Morton-order linearized octree/quadtree
- Immersed Boundary (IB) method for complex geometries
- High-order WENO finite difference schemes
- Multiple GPU acceleration backends (CUDA Fortran, OpenACC, OpenMP offloading)

**Performance Characteristics**:
- Target architecture: Multi-GPU clusters with InfiniBand/NVLink interconnects
- Parallel paradigms: MPI (distributed memory) + OpenMP/OpenACC/CUDA (node-level parallelism)
- Memory model: 5D field arrays `(nv, ni, nj, nk, nb)` - variables, i/j/k cells, blocks
- Communication pattern: Nearest-neighbor ghost cell exchange (AMR complicates stencil)
- Scalability goal: Strong scaling to O(1000) GPUs with >70% parallel efficiency

## Code Architecture

### Directory Structure

```
src/
├── app/              # Applications built on ADAM framework
│   ├── nasto/        # Navier-Stokes solver (CPU, NVF, FNL, GMP backends)
│   ├── prism/        # Maxwell equations solver
│   ├── chase/        # CFD application
│   ├── patch/        # Patch-based application
│   └── ascot/        # Binary-to-ASCII converter utility
├── lib/              # Core libraries
│   ├── common/       # Base ADAM objects (CPU, portable)
│   ├── fnl/          # OpenACC GPU implementations
│   ├── nvf/          # CUDA Fortran GPU implementations
│   └── gmp/          # OpenMP offloading (in development)
├── tests/            # Unit and integration tests
└── third_party/      # Git submodules (PENF, StringiFor, FiNeR, VTKFortran, etc.)
```

### Core Objects (src/lib/common/)

| Object | Purpose |
|--------|---------|
| `adam_grid_object` | Structured block grid management with AMR |
| `adam_tree_object` | Octree/quadtree with Morton ordering |
| `adam_field_object` | Field variable storage and interpolation |
| `adam_weno_object` | High-order WENO reconstruction |
| `adam_rk_object` | Runge-Kutta temporal integration |
| `adam_ib_object` | Immersed boundary method |
| `adam_io_object` | HDF5 and restart file I/O |
| `adam_mpih_object` | MPI wrapper and communication |
| `adam_mpih_global` | Program-scope MPI handler singleton (module variable, not embedded in types) |
| `adam_grid_global` | Program-scope grid singleton (module variable, not embedded in types) |
| `adam_field_global` | Program-scope field singleton (module variable, not embedded in types) |
| `adam_maps_global` | Program-scope maps singleton (module variable, not embedded in types) |
| `adam_tree_global` | Program-scope tree singleton (module variable, not embedded in types) |
| `adam_weno_global` | Program-scope WENO singleton (module variable, not embedded in types) |
| `adam_ib_global` | Program-scope IB singleton (module variable, not embedded in types) |
| `adam_rk_global` | Program-scope RK singleton (module variable, not embedded in types) |
| `adam_globals` | Convenience aggregator re-exporting all 8 CPU singletons (`mpih, grid, field, maps, tree, weno, ib, rk`) |
| `adam_fdv_operators_library` | Gradient, divergence, curl, Laplacian operators |
| `adam_riemann_euler_library` | Euler equation Riemann solvers |

### FNL Backend Singletons (src/lib/fnl/ and src/app/prism/fnl/)

FNL GPU objects are also exposed as program-scope singletons, eliminating composition-by-pointer in solver types:

| Singleton | Module | Purpose |
|-----------|--------|---------|
| `mpih_fnl` | `adam_fnl_mpih_global` | FNL MPI handler |
| `field_fnl` | `adam_fnl_field_global` | FNL field (GPU arrays, ghost maps) |
| `ib_fnl` | `adam_fnl_ib_global` | FNL immersed boundary |
| `rk_fnl` | `adam_fnl_rk_global` | FNL Runge-Kutta integrator |
| `weno_fnl` | `adam_fnl_weno_global` | FNL WENO reconstructor |
| `coil_fnl` | `adam_prism_fnl_coil_global` | FNL coil source (PRISM only) |
| `fwlayer_fnl` | `adam_prism_fnl_fwlayer_global` | FNL fWLayer (PRISM only) |

**Initialization pattern**: Solver `initialize` must copy CPU value singletons before calling FNL inits:
```fortran
ib = self%ib ; rk = self%rk ; weno = self%weno  ! populate CPU singletons
call field_fnl%initialize(verbose=.true.)
call ib_fnl%initialize()
call rk_fnl%initialize()
call weno_fnl%initialize()
```

### Backend Pattern

Each application has multiple backends sharing common code:
```
app/<name>/cpu/     -> CPU-only (MPI/OpenMP)
app/<name>/nvf/     -> CUDA Fortran
app/<name>/fnl/     -> OpenACC
app/<name>/gmp/     -> OpenMP target offloading
app/<name>/common/  -> Shared across all backends
```

GPU backends extend base objects with GPU-accelerated implementations in `lib/fnl/`, `lib/nvf/`, or `lib/gmp/`.

### Preprocessor Macros

**Active Macros** (defined via FoBiS compiler flags):
- `_MPI_` - Enable MPI (always defined for this project)
- `_NVF` - CUDA Fortran backend (mutually exclusive with `_FNL`, `_GMP`)
- `_FNL` - OpenACC backend using FNL library
- `_GMP` - OpenMP target offloading (experimental)
- `DEV_OAC` - Marks device-side code for OpenACC compiler
- `COMPILER_NVF` / `COMPILER_GNU` / `COMPILER_INTEL` - Compiler-specific workarounds

**Usage Patterns**:
```fortran
#ifdef _NVF
  ! CUDA Fortran specific code
  attributes(device) :: array
#elif defined(_FNL)
  !$acc declare device_resident(array)
#else
  ! CPU code
#endif
```

**Guidelines**:
- Minimize preprocessor use - prefer Fortran `select case` or procedure pointers for backend dispatch
- Use macros only for: backend selection, compiler bugs/workarounds, conditional compilation of entire routines
- Document WHY macro is needed (reference compiler version, bug report, or feature unavailability)
- Never use macros for math expressions or loop bounds (hinders debugging and readability)

## Source File Conventions

### File Extensions and Structure
- All Fortran source files use `.F90` extension (uppercase F = preprocessor-enabled)
- One module per file (exception: closely related types in single file)
- Filename matches module name: `adam_grid_object.F90` contains `module adam_grid_object`

### Language Features and Style
- Object-oriented Fortran 2003+ with derived types and type extension
- Use `submodule` for implementation when module interface is stable (reduces recompilation)
- Always `implicit none` at module and procedure level
- Explicit `intent(in/out/inout)` for all dummy arguments

### Array Memory Layout (Critical for Performance)

**Field Storage Convention**: 5D arrays with shape `(nv, ni, nj, nk, nb)`
- `nv`: Number of variables (e.g., density, velocity components, energy)
- `ni, nj, nk`: Number of cells in i, j, k directions (including ghost cells)
- `nb`: Number of blocks (for AMR and domain decomposition)

**Fortran Column-Major Layout**:
```fortran
! GOOD: Stride-1 access on leftmost index (nv)
do k = 1, nk
  do j = 1, nj
    do i = 1, ni
      do iv = 1, nv  ! Innermost loop - contiguous in memory
        field(iv, i, j, k, ib) = ...
      end do
    end do
  end do
end do

! BAD: Strided access (cache-inefficient)
do iv = 1, nv
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni  ! Non-contiguous access
        field(iv, i, j, k, ib) = ...
      end do
    end do
  end do
end do
```

**GPU Considerations**:
- For GPU kernels, collapse outer loops: `!$acc parallel loop collapse(3) gang vector`
- Ensure coalesced access: adjacent threads should access adjacent `iv` values
- Use `contiguous` attribute to guarantee no hidden copies: `real(rp), contiguous :: field(:,:,:,:,:)`

### Ghost Cells (`ngc`)
- Configurable ghost cell width for high-order stencil operations
- Typical values: `ngc = 3` (WENO5), `ngc = 4` (WENO7)
- Array bounds: `field(nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb)`
- Ghost cell exchange via MPI before stencil operations

### Derived Type Best Practices
```fortran
type :: my_object
  integer :: n
  real(rp), allocatable :: data(:,:,:)  ! Prefer allocatable over pointer
contains
  procedure :: init
  procedure :: destroy
  procedure :: compute
  final :: cleanup  ! Automatic cleanup on deallocation
end type

! Use type-bound procedures for encapsulation
call obj%init(params)
call obj%compute()
call obj%destroy()
```

**Allocatable vs Pointer**:
- Prefer `allocatable`: compiler optimizes better, automatic deallocation
- Use `pointer` only when: aliasing required, linked data structures, or C interop

## Build Commands

The build system uses [FoBiS](https://github.com/szaghi/FoBiS) — CLI binary `fobis` (3.8+). All flags are double-dash long form; the legacy `FoBiS.py build -mode X` short-form is no longer accepted.

```bash
# List all available build modes
fobis build --lmodes

# Build NASTO (Navier-Stokes solver)
fobis build --mode nasto-cpu-gnu          # CPU with GNU compiler
fobis build --mode nasto-nvf-cuda         # NVIDIA CUDA Fortran
fobis build --mode nasto-fnl-nvf-oac      # OpenACC with NVF

# Build PRISM (Maxwell equations solver)
fobis build --mode prism-gnu              # CPU with GNU compiler
fobis build --mode prism-fnl-nvf-oac      # OpenACC with NVF

# Build other applications
fobis build --mode chase-gnu              # CHASE app
fobis build --mode patch-gnu              # PATCH app
fobis build --mode ascot-nvf-cuda         # ASCOT converter utility

# Build ADAM core libraries
fobis build --mode adam-com-gnu           # Common library (GNU)
fobis build --mode adam-com-nvf           # Common library (NVF)
fobis build --mode adam-nvf-cuda          # NVF CUDA library

# Debug builds (append -debug to mode name)
fobis build --mode nasto-cpu-gnu-debug
```

## Running Tests

Build and run FDV (finite difference/volume) operator tests:

```bash
# Build a specific test
fobis build --mode test-fdv-gradient-trigonometric-gnu
fobis build --mode test-fdv-divergence-trigonometric-gnu
fobis build --mode test-fdv-curl-trigonometric-gnu
fobis build --mode test-fdv-laplacian-trigonometric-gnu
fobis build --mode test-fdv-operators-trigonometric-gnu
fobis build --mode test-fdv-operators-step-gnu

# Run with MPI (executables are in exe/)
mpirun -np <N> exe/<test_executable>
```

Integration test cases are in `src/tests/nasto/` (Riemann problems, shock-sphere).

## Testing Strategy

### Multi-Level Testing Pyramid

1. **Unit Tests** (`src/tests/`):
   - Finite difference operators: gradient, divergence, curl, Laplacian (trigonometric analytical solutions)
   - WENO reconstruction accuracy (known smooth and discontinuous test functions)
   - Riemann solvers (exact Riemann problem solutions: Sod, Lax, Shu-Osher)
   - Build and run: `fobis build --mode test-fdv-<operator>-<test>-<compiler>`

2. **Integration Tests** (`src/tests/nasto/`):
   - Riemann problems (1D/2D/3D): Sod shock tube, shock-sphere interaction
   - Known CFD benchmarks: Taylor-Green vortex, Kelvin-Helmholtz instability
   - AMR correctness: verify refinement/coarsening conserves quantities

3. **Regression Tests**:
   - Maintain reference solutions for key test cases
   - Check-in expected output (or checksums) to detect unintended changes
   - Automate comparison of field L2/L-infinity norms

4. **Scaling Tests**:
   - Weak scaling: verify efficiency >80% to target scale (e.g., 1024 GPUs)
   - Strong scaling: identify parallel efficiency knee (optimal process count)
   - Communication profiling: check MPI time percentage <20% for good scaling

### Compiler Testing Matrix

All changes affecting core numerics or parallel logic must pass:
- **GNU gfortran** (latest stable): baseline CPU, strict bounds checking
- **Intel ifort/ifx**: aggressive optimization, vectorization reports
- **NVIDIA nvfortran**: GPU targets (CUDA Fortran, OpenACC)

Use `-fbounds-check -fcheck=all` (GNU) or `-check all -traceback` (Intel) during development.

## External Dependencies

### Compilers and Runtime
- **NVIDIA HPC SDK** (nvfortran): Required for GPU builds (OpenACC, CUDA Fortran)
  - Version sensitivity: OpenACC 3.x features require SDK >=22.x
  - CUDA compatibility: Verify nvfortran CUDA version matches driver
  - Environment: Set `HDF5_nvf` for NVIDIA-compiled HDF5 libraries

- **GNU Compiler Collection** (gfortran): CPU builds and baseline testing
  - Minimum version: 9.x for Fortran 2008 submodules
  - Recommended: Latest stable (14.x) for improved diagnostics and optimization

- **Intel oneAPI** (ifort/ifx): Optional, for additional optimization and validation
  - ifx (LLVM-based) recommended for new projects
  - Excellent vectorization reports: use `-qopt-report=5`

### Parallel Libraries
- **MPI**: OpenMPI >=4.0 or MPICH >=3.3 for distributed computing
  - OpenMPI recommended for GPU-aware MPI (UCX transport)
  - Verify CUDA-aware MPI if using GPU-direct communication
  - Set `OMPI_MCA_btl_openib_allow_ib=1` for InfiniBand clusters

### I/O and Serialization
- **HDF5**: Parallel I/O library (version >=1.10 required)
  - **Critical**: Must be compiled with same MPI library as application
  - Requires ZLIB and SZIP for compression support
  - Library paths configured in fobos `[common-variables]`:
    - GNU: `$HDF5_gnu = lib/hdf5/develop/gnu/14.2.0`
    - NVF: `$HDF5_nvf = lib/hdf5/develop/nvf/25.11`
  - Collective I/O optimization: Set `H5Pset_dxpl_mpio` for MPI-IO hints

### Build System
- **FoBiS** (Fortran Building System): Manages complex multi-backend builds
  - Install: `pip install FoBiS.py`
  - Centralizes compiler flags, preprocessor macros, dependency tracking

### Cluster-Specific Considerations
- **Job Scheduler**: Slurm (assumed)
  - Use `#SBATCH --gpus-per-node=<N>` for GPU allocation
  - MPI binding: `srun --cpu-bind=cores --gpu-bind=closest` for affinity

- **Module System**: Load correct versions before building
  ```bash
  module load nvhpc/25.11  # NVIDIA compilers
  module load openmpi/4.1.6-cuda  # GPU-aware MPI
  module load hdf5-parallel/1.14.0
  ```

- **Network Topology**: High-performance interconnect (InfiniBand, Slingshot)
  - Verify RDMA support: `ibv_devinfo` or `fi_info`
  - UCX tuning for OpenMPI: `export UCX_TLS=rc_x,sm,cuda_copy,cuda_ipc`

## Configuration Files

Applications use INI files (parsed by FiNeR library) for simulation parameters:
- Grid parameters (domain bounds, resolution)
- Physics parameters (gas properties, Reynolds numbers)
- Numerical parameters (WENO order, RK stages)
- Boundary/initial conditions
- I/O options and AMR refinement criteria

## Authoritative References

### Standards and Specifications
- **MPI Standard**: [MPI Forum](https://www.mpi-forum.org/docs/) - cite specific MPI version and section
- **OpenACC Specification**: [OpenACC.org](https://www.openacc.org/specification) - current version 3.3
- **Fortran Standards**: ISO/IEC 1539-1 (Fortran 2018), J3 papers for upcoming features
- **OpenMP Specification**: [OpenMP.org](https://www.openmp.org/specifications/) - version 5.2 for GPU offloading

### Performance and Algorithms
- **HPC Texts**:
  - "Introduction to High Performance Computing for Scientists and Engineers" (Hager & Wellein)
  - "Parallel Programming for Science and Engineering" (Eijkhout)
  - "Using MPI" 3rd ed. (Gropp, Lusk, Thakur)

- **GPU Programming**:
  - NVIDIA CUDA C++ Best Practices Guide
  - "Programming Massively Parallel Processors" (Kirk & Hwu)
  - OpenACC Best Practices Guide (NVIDIA)

- **Numerical Methods**:
  - "Riemann Solvers and Numerical Methods for Fluid Dynamics" (Toro)
  - "Computational Fluid Dynamics" (Anderson) - for CFD context
  - "Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes" (Shu, 2009) - WENO reference

- **Fortran Style Guides**:
  - [Zen of Fortran](https://github.com/szaghi/zen-of-fortran) - opinionated coding guidelines
  - "Modern Fortran: Style and Usage" (Clerman & Spector)
  - [Fortran Best Practices](http://www.fortran90.org/src/best-practices.html)
  - [European Standards For Writing and Documenting Exchangeable Fortran 90 Code](http://research.metoffice.gov.uk/research/nwp/numerical/fortran90/f90_standards.html)

### Tools and Profiling
- **NVIDIA Nsight Compute**: [Documentation](https://docs.nvidia.com/nsight-compute/)
- **Intel VTune Profiler**: Performance analysis for CPU/GPU
- **Arm Forge (DDT/MAP)**: Parallel debugger and profiler

**Citation Practice**: Always cite specific sections/versions when referencing standards. Example: "According to MPI-3.1 section 3.2.5, non-blocking collectives..."
