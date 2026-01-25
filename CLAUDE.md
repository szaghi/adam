# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Philosophy

**Target**: Exascale HPC applications with strong scaling to thousands of GPU nodes
**Priorities**: Correctness first, then performance, maintainability, and portability
**Evidence-based optimization**: Always profile before optimizing, benchmark after changes
**Documentation standard**: Well-referenced explanations citing authoritative sources over speculation

## Agent Mode: Performance-Driven Development Workflow

When analyzing, modifying, or optimizing code in this repository, follow this systematic approach:

### 1. **Profile Before Action**
- **Always ask first**: "Do you have profiling data?" (Nsight Compute, nvprof, gprof, Intel VTune)
- Identify bottlenecks: memory bandwidth, compute intensity, communication overhead, load imbalance
- Never optimize without quantitative evidence of the problem

### 2. **Correctness Verification Strategy**
- Propose changes with explicit verification plan
- Specify test cases: unit tests, integration tests (Riemann problems), regression tests
- Multi-compiler validation: gfortran, Intel ifort/ifx, NVIDIA nvfortran
- Enable runtime checks during development: `-fbounds-check`, `-fcheck=all` (GNU), `-check all` (Intel)
- Memory debugging: Valgrind (CPU), cuda-memcheck (GPU)

### 3. **Performance Optimization Patterns**

#### GPU-Specific (OpenACC, CUDA Fortran)
- **Data movement**: Minimize host-device transfers, maximize data reuse on device
- **Loop directives**: Explicitly specify `gang`, `vector`, `seq` collapse depth
- **Atomic operations**: Red flag for serialization - consider alternatives (graph coloring, privatization, reduction)
- **Memory coalescing**: Ensure contiguous access patterns in innermost loops
- **Kernel fusion**: Combine multiple kernels to reduce launch overhead and improve data locality

#### MPI Best Practices
- **Communication/computation overlap**: Use non-blocking MPI (Isend/Irecv, Iallreduce)
- **Domain decomposition**: Check load balancing across ranks - propose block-splitting if needed
- **Collective optimization**: Use MPI-3 neighborhood collectives for AMR stencil operations
- **One-sided communication**: Consider RMA (Put/Get) for ghost cell exchange

#### Memory Layout (Fortran-specific)
- **Column-major awareness**: Inner loops on first index for cache efficiency
- **Array contiguity**: Use `contiguous` attribute, avoid unnecessary copies
- **Derived types**: Watch for alignment, padding, and GPU transfer efficiency
- **Avoid temporary arrays**: Use in-place operations or explicit buffers

### 4. **Benchmarking Requirements**
For any performance-related change, specify:
- **Metrics**: Wall time, speedup, scaling efficiency, memory bandwidth utilization
- **Weak scaling**: Fixed problem size per process/GPU
- **Strong scaling**: Fixed total problem size, varying parallelism
- **Target**: Justify "acceptable performance" (e.g., >80% efficiency to 1024 GPUs)

### 5. **Documentation of Changes**
- Comment optimizations with: rationale, expected impact, profiling evidence
- Reference specific sections of standards/papers when using advanced features
- Note compiler-specific behavior (especially GCC vs Intel vs NVIDIA differences)

## Code Modification Guidelines

### The Zen of Fortran

*Opinionated coding guidelines for this codebase, inspired by the [Zen of Fortran](https://github.com/szaghi/zen-of-fortran):*

1. **Standard compliance** is better than *custom or extended*
2. **Beautiful** is better than *ugly*:
   - **Readability** counts
   - **Explicit** is better than *impl.*
   - **Simple** is better than *CoMpleX*
   - **CoMpleX** is better than *c0mp1|c@ted*
   - **Flat** is better than *nested*
   - **S p a r s e** is better than *dense*
3. **Fast** is better than *slow*:
   - **Vector** is better than *loop*
   - **Matrix** is better than *vector*
   - **Strided** is better than *scattered*
   - **Contiguous** is better than *strided*
   - **Broadcasting** is a great idea, use where possible
4. **Slow** is better than *unmaintainable*
5. Make it look like the **math**
6. **Special cases** aren't special enough to break rules...
7. Although **practicality** beats *purity*
8. **Pure** procedure is better than *impure*...
9. Although **practicality** beats *purity* again
10. **Private** is better than *public*
11. **Errors** should never pass *silently*...
12. Unless **errors** are explicitly *silenced*

#### Standard Compliance

Standard-compliant code is **portable**:
- Free from particular OS/hardware architecture (ensures long-time code life, easy enabling on new hardware)
- Free from particular compiler vendor (can use different compilers for cross-checking/debugging)
- Use compiler flags to check standard adherence: `-std=f2018` (GNU), `-std18` (Intel)

#### Readability Guidelines

**Use free source form**:
- 132 characters are better than 72 (use modern wide screens)
- For mathematical formulas, avoid splitting on multiple lines when possible
- Balance: don't fill each line to 132 chars, but don't limit yourself to 72

**Indentation**:
- Use consistent number of spaces (project standard: 3 spaces)
- Increase indentation when data scope changes
- Indent blocks within all control constructs
- Indent code after named constructs so names stand out
- Break lines in logical places; indent continued lines double the block indent
- Use blank lines to separate related parts of a program

**White spaces**:
- Avoid hard tabs
- Use spaces around operators: `if (foo > bar) then`
- Align similar code statements (declarations, comments)
- Place a space after all commas: `x(i, j, k) = foo(i, j, k:k+1)`

**Example of good style**:
```fortran
module well_formatted_module
implicit none
private
public :: foo_type

type :: foo_type
   private
   logical                       :: is_good = .false.
   character(len=:), allocatable :: name
   contains
      private
      procedure, public, pass(self) :: init
endtype foo_type

contains
   pure subroutine init(self, name)
   class(foo_type),        intent(inout) :: self
   character(*), optional, intent(in)    :: name

   if (allocated(self%name)) deallocate(self%name)
   check_name: if (present(name)) then
                  if (len_trim(adjustl(name)) == 0) exit check_name
                  self%name = trim(adjustl(name))
               else
                  self%name = 'J. Doe'
   endif check_name
   if (allocated(self%name)) self%is_good = .true.
   endsubroutine init
endmodule well_formatted_module
```

#### Private by Default

- All module entities should be `private` by default
- Explicitly declare `public` only what is part of the API
- Type components should be `private` unless external access is necessary

#### Error Handling

- Never ignore error codes from allocate, MPI, I/O operations
- Use `stat=` and `errmsg=` with allocate/deallocate
- Check MPI return codes explicitly
- Provide meaningful error messages with context
- Use `error stop` with informative messages for unrecoverable errors

### Fortran Style and Standards

**Language Features**:
- Use modern Fortran (2003/2008/2018): `submodule`, `associate`, `block`, `error stop`
- Explicit interfaces via `module` and `submodule` (mandatory for all public procedures)
- `implicit none` everywhere (no exceptions)
- Use `iso_fortran_env` and `iso_c_binding` for portability

**Array Programming**:
- Prefer whole-array syntax when vectorizable: `a(:,:,:) = b(:,:,:) + c(:,:,:)`
- Use assumed-shape arrays for flexibility: `real(rp), intent(in) :: a(:,:,:)`
- Add `contiguous` attribute for performance-critical paths: `real(rp), intent(in), contiguous :: a(:,:,:)`
- Avoid assumed-size arrays: `real(rp) :: a(*)` (use assumed-shape or explicit-shape)

**Type-Bound Procedures**:
- Use type-bound procedures for polymorphism and encapsulation
- Mark `procedure :: method => implementation` for clarity
- Use `nopass` only when justified (static-like methods)

**Preprocessor Usage**:
- Minimize preprocessor macros - prefer Fortran abstractions when possible
- When unavoidable, use for: compiler detection, backend selection, MPI/GPU conditionals
- Never use macros for code that can be written in standard Fortran

### Common Pitfalls to Actively Check

1. **Atomic Operations on GPUs**: 
   - Cause severe serialization (observed: 10-100x slowdown in vertex-based operations)
   - Solution: Graph coloring, data reordering, or reduction strategies
   
2. **Array Reallocation with Derived Types**:
   - GCC gfortran bounds checking differs from Intel/NVIDIA (observed: silent corruption)
   - Always verify bounds after `allocate`/`deallocate` with multiple compilers
   
3. **OpenACC Data Directives**:
   - Unstructured data regions (`enter data`/`exit data`) require careful lifecycle management
   - Missing `update device/host` causes stale data bugs (hard to debug)
   - Use `present` clause to verify data residency assumptions

4. **MPI Datatypes**:
   - Derived type MPI commits can fail silently - always check error codes
   - Padding/alignment issues when transferring structs across heterogeneous nodes

5. **Floating-Point Reproducibility**:
   - Non-associative reductions in MPI may cause answer changes with different process counts
   - GPU atomics introduce non-determinism - document or eliminate

### Performance Anti-Patterns

**Do NOT**:
- Optimize without profiling data
- Use `!$acc kernels` without loop directives (compiler may not parallelize optimally)
- Ignore compiler warnings (especially `-Wunused-variable`, `-Wimplicit-interface`)
- Mix MPI library versions or HDF5 versions across build
- Assume performance portability (GPU code often needs per-architecture tuning)

**DO**:
- Start with `!$acc parallel loop` with explicit clauses
- Use compiler feedback: `-Minfo=accel` (NVIDIA), `-qopt-report` (Intel)
- Validate correctness with sanitizers: AddressSanitizer (CPU), cuda-memcheck (GPU)
- Profile with multiple problem sizes (small, medium, large)
- Test on target architecture before claiming "optimization complete"

## Testing Strategy

### Multi-Level Testing Pyramid

1. **Unit Tests** (`src/tests/`):
   - Finite difference operators: gradient, divergence, curl, Laplacian (trigonometric analytical solutions)
   - WENO reconstruction accuracy (known smooth and discontinuous test functions)
   - Riemann solvers (exact Riemann problem solutions: Sod, Lax, Shu-Osher)
   - Build and run: `FoBiS.py build -mode test-fdv-<operator>-<test>-<compiler>`

2. **Integration Tests** (`src/tests/nasto/`):
   - Riemann problems (1D/2D/3D): Sod shock tube, shock-sphere interaction
   - Known CFD benchmarks: Taylor-Green vortex, Kelvin-Helmholtz instability
   - AMR correctness: verify refinement/coarsening conserves quantities

3. **Regression Tests**:
   - Maintain reference solutions for key test cases
   - Check-in expected output (or checksums) to detect unintended changes
   - Automate comparison of field L2/L∞ norms

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

## Authoritative References

When providing suggestions or explanations, prioritize these sources:

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

**Citation Practice**: Always cite specific sections/versions when referencing standards. Example: "According to MPI-3.1 §3.2.5, non-blocking collectives..."

## Common Agent Tasks

### Task: GPU Kernel Optimization

**Input**: User reports kernel is underperforming (low occupancy, low utilization)

**Response Protocol**:
1. Request profiling data: occupancy, achieved memory bandwidth, warp execution efficiency
2. Analyze memory access pattern: coalesced vs. strided, transactions per request
3. Check for: bank conflicts (shared memory), atomic contention, excessive register usage
4. Propose optimization with quantitative prediction: "Expected to improve bandwidth from X to Y GB/s"
5. Suggest verification: rerun profiler, compare before/after metrics

**Example**: If atomic operations detected → suggest graph coloring (as used for vertex-based AMR operations in this codebase)

### Task: MPI Load Balancing

**Input**: User observes poor strong scaling or high MPI wait time

**Response Protocol**:
1. Request timing breakdown: computation vs. communication time per rank
2. Check domain decomposition: are blocks evenly distributed?
3. For AMR: analyze refinement pattern - clustered refinement causes imbalance
4. Propose solutions:
   - Dynamic load balancing (redistribute blocks between MPI ranks)
   - Space-filling curve (Morton/Hilbert) reordering for locality
   - Overdecomposition with task-based parallelism (MPI+X)

### Task: Fortran Modernization

**Input**: User has legacy Fortran 77/90 code to modernize

**Response Protocol**:
1. Identify obsolescent features: `common` blocks, `equivalence`, assumed-size arrays, `goto`
2. Propose replacements: `module` data, derived types, assumed-shape arrays, structured control flow
3. Introduce `submodule` for large modules (reduces compilation dependencies)
4. Convert to assumed-shape arrays: `real :: a(:,:,:)` instead of `real :: a(nx,ny,nz)`
5. Use `iso_fortran_env` for `real64`, `int32` instead of `selected_real_kind`

### Task: Debugging Numerical Instability

**Input**: Simulation crashes or produces non-physical results

**Response Protocol**:
1. Check CFL condition: `dt < CFL * min(dx, dy, dz) / max_wavespeed`
2. Verify boundary conditions: wall, symmetry, inflow/outflow implementations
3. Check for NaN/Inf: enable floating-point exception trapping (`-ffpe-trap=invalid,zero,overflow`)
4. Review discretization: upwind/central bias, dissipation terms, limiter functions
5. Suggest diagnostic output: write fields at crash point, check conservation errors

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

## Build Commands

The build system uses [FoBiS](https://github.com/szaghi/FoBiS). Common build commands:

```bash
# List all available build modes
FoBiS.py build -lmodes

# Build NASTO (Navier-Stokes solver)
FoBiS.py build -mode nasto-cpu-gnu          # CPU with GNU compiler
FoBiS.py build -mode nasto-nvf-cuda         # NVIDIA CUDA Fortran
FoBiS.py build -mode nasto-fnl-nvf-oac      # OpenACC with NVF

# Build PRISM (Maxwell equations solver)
FoBiS.py build -mode prism-gnu              # CPU with GNU compiler
FoBiS.py build -mode prism-fnl-nvf-oac      # OpenACC with NVF

# Build other applications
FoBiS.py build -mode chase-gnu              # CHASE app
FoBiS.py build -mode patch-gnu              # PATCH app
FoBiS.py build -mode ascot-nvf-cuda         # ASCOT converter utility

# Build ADAM core libraries
FoBiS.py build -mode adam-com-gnu           # Common library (GNU)
FoBiS.py build -mode adam-com-nvf           # Common library (NVF)
FoBiS.py build -mode adam-nvf-cuda          # NVF CUDA library

# Debug builds (append -debug to mode name)
FoBiS.py build -mode nasto-cpu-gnu-debug
```

## Running Tests

Build and run FDV (finite difference/volume) operator tests:

```bash
# Build a specific test
FoBiS.py build -mode test-fdv-gradient-trigonometric-gnu
FoBiS.py build -mode test-fdv-divergence-trigonometric-gnu
FoBiS.py build -mode test-fdv-curl-trigonometric-gnu
FoBiS.py build -mode test-fdv-laplacian-trigonometric-gnu
FoBiS.py build -mode test-fdv-operators-trigonometric-gnu
FoBiS.py build -mode test-fdv-operators-step-gnu

# Run with MPI (executables are in exe/)
mpirun -np <N> exe/<test_executable>
```

Integration test cases are in `src/tests/nasto/` (Riemann problems, shock-sphere).

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
| `adam_fdv_operators_library` | Gradient, divergence, curl, Laplacian operators |
| `adam_riemann_euler_library` | Euler equation Riemann solvers |

### Backend Pattern

Each application has multiple backends sharing common code:
```
app/<name>/cpu/     → CPU-only (MPI/OpenMP)
app/<name>/nvf/     → CUDA Fortran
app/<name>/fnl/     → OpenACC
app/<name>/gmp/     → OpenMP target offloading
app/<name>/common/  → Shared across all backends
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

## External Dependencies

### Compilers and Runtime
- **NVIDIA HPC SDK** (nvfortran): Required for GPU builds (OpenACC, CUDA Fortran)
  - Version sensitivity: OpenACC 3.x features require SDK ≥22.x
  - CUDA compatibility: Verify nvfortran CUDA version matches driver
  - Environment: Set `HDF5_nvf` for NVIDIA-compiled HDF5 libraries
  
- **GNU Compiler Collection** (gfortran): CPU builds and baseline testing
  - Minimum version: 9.x for Fortran 2008 submodules
  - Recommended: Latest stable (14.x) for improved diagnostics and optimization
  
- **Intel oneAPI** (ifort/ifx): Optional, for additional optimization and validation
  - ifx (LLVM-based) recommended for new projects
  - Excellent vectorization reports: use `-qopt-report=5`

### Parallel Libraries
- **MPI**: OpenMPI ≥4.0 or MPICH ≥3.3 for distributed computing
  - OpenMPI recommended for GPU-aware MPI (UCX transport)
  - Verify CUDA-aware MPI if using GPU-direct communication
  - Set `OMPI_MCA_btl_openib_allow_ib=1` for InfiniBand clusters
  
### I/O and Serialization
- **HDF5**: Parallel I/O library (version ≥1.10 required)
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

---

## Quick Reference: Common Agent Workflows

### "I want to optimize this GPU kernel"
1. ✅ Ask for profiling data (Nsight Compute report)
2. ✅ Identify bottleneck: compute, memory bandwidth, or atomic serialization
3. ✅ Propose specific optimization (e.g., loop fusion, memory layout change)
4. ✅ Predict quantitative improvement ("expect 2x speedup based on...")
5. ✅ Suggest verification: re-profile, compare metrics

### "This MPI code doesn't scale well"
1. ✅ Request timing breakdown (computation vs. communication)
2. ✅ Check load balance (are some ranks waiting?)
3. ✅ Analyze communication pattern (point-to-point vs. collective)
4. ✅ Propose: non-blocking MPI, load rebalancing, or communication reduction
5. ✅ Suggest scaling test: weak/strong scaling plots

### "Help me modernize this Fortran code"
1. ✅ Identify obsolescent features (`common`, `equivalence`, `goto`, assumed-size arrays)
2. ✅ Propose modern replacements (`module`, derived types, assumed-shape arrays)
3. ✅ Introduce `submodule` for large modules
4. ✅ Use `iso_fortran_env` for portable kinds
5. ✅ Test with multiple compilers (GNU, Intel, NVIDIA)

### "The simulation produces non-physical results"
1. ✅ Check CFL condition (`dt` vs. grid spacing and wave speed)
2. ✅ Verify boundary conditions (implementation and values)
3. ✅ Enable FP exception trapping (`-ffpe-trap=invalid,zero,overflow`)
4. ✅ Inspect discretization (upwinding, limiters, dissipation)
5. ✅ Request diagnostic output (field dumps, conservation checks)

### "I need to add a new feature"
1. ✅ Understand existing architecture (which backend? CPU/GPU/both?)
2. ✅ Propose design: module structure, type extensions, interfaces
3. ✅ Identify test cases for verification
4. ✅ Consider multi-compiler compatibility
5. ✅ Document with references (papers, standards)

---

## Emergency Debugging Checklist

**Compilation Fails**:
- [ ] Check module dependencies (`FoBiS.py build -lmodes` to verify mode exists)
- [ ] Verify HDF5/MPI library paths in fobos configuration
- [ ] Load correct module environment (compiler, MPI, HDF5 versions must match)
- [ ] Check preprocessor macros (only one of `_NVF`, `_FNL`, `_GMP` should be defined)

**Runtime Crash**:
- [ ] Enable bounds checking: `-fbounds-check` (GNU), `-check bounds` (Intel)
- [ ] Enable FP exceptions: `-ffpe-trap=invalid,zero,overflow` (GNU)
- [ ] Run with debugger: `gdb --args mpirun -np 4 ./exe/nasto`
- [ ] Check stack size: `ulimit -s unlimited` (Fortran recursive calls)
- [ ] Verify MPI library: `ldd exe/nasto | grep mpi` (should match loaded module)

**GPU Errors**:
- [ ] CUDA error codes: check return values, run with `cuda-memcheck`
- [ ] OpenACC: set `ACC_SYNCHRONOUS=1` for immediate error reporting
- [ ] Device memory: verify sufficient GPU memory (`nvidia-smi`)
- [ ] Data residency: check `!$acc update device/host` for stale data

**Performance Regression**:
- [ ] Profile before/after: compare wall time, memory bandwidth, kernel launch overhead
- [ ] Check compiler optimization flags: should be `-O3` or equivalent
- [ ] Verify no debug flags in production build (`-g`, `-fbounds-check`)
- [ ] Test on same hardware (CPU model, GPU model, memory size)
- [ ] Review recent changes: `git diff <last-good-commit>`

**MPI Hangs**:
- [ ] Deadlock detection: check for unmatched Send/Recv pairs
- [ ] Barrier analysis: ensure all ranks reach collective operations
- [ ] Set `MPICH_ASYNC_PROGRESS=1` or `OMPI_MCA_opal_progress_threads=1`
- [ ] Reduce process count to isolate issue (e.g., test with 2 ranks)
- [ ] Enable MPI debug: `export MPICH_DEBUG=1` or `--mca btl_base_verbose 30`

