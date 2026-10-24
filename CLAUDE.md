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
│   ├── chase/        # Euler solver (deprecated, superseded by FLUME)
│   ├── patch/        # Patch-based application
│   ├── flume/        # Compressible Euler (M1), MHD in development; supersedes CHASE
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
| `adam_realm_object` | Abstract app base: **owns every library object as a component** (`io, amr, slices, weno, ib, rk, flail, flail_legacy, adam{grid,tree,field,maps}`; `flail` is the elliptic solver library, issue #84) + the `_forest` TBP contract |
| `adam_forest_object` | Time-loop orchestrator (`simulate`, `simulate_from_manifest`): single-realm fast path, AMR-staged path, multi-realm seams |
| `adam_mpih_global` | The **only** CPU program-scope singleton (`mpih`); `adam_globals` re-exports just this |
| `adam_fdv_operators_library` | Gradient, divergence, curl, Laplacian operators |
| `adam_riemann_euler_library` | Euler equation Riemann solvers |

### Object ownership: realm composition, not singletons

The only program-scope singletons are **`mpih`** (`adam_mpih_global`) and **`mpih_fnl`** (`adam_fnl_mpih_global`). The former `adam_{grid,field,maps,tree,weno,ib,rk}_global` and `adam_fnl_{field,ib,rk,weno}_global` / `adam_prism_fnl_{coil,fwlayer}_global` modules **no longer exist**. Every other object is a per-realm component:

- CPU library objects live in `realm_object` (`self%adam%grid`, `self%adam%field`, `self%adam%maps`, `self%rk`, `self%weno`, `self%ib`, …) and are passed **as arguments** to library calls (e.g. `self%adam%field%update_ghost_local(q=q, grid=self%adam%grid, maps=self%adam%maps)`).
- FNL helpers are components of the backend type (`prism_fnl_object%field_fnl/ib_fnl/rk_fnl/weno_fnl/coil_fnl`, `adam_prism_fnl_object.F90:47-56`), initialised **after** the common init by argument passing:

```fortran
call self%prism_common_object%initialize(...)                                        ! builds self%adam, self%rk, ...
call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps)
call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
call self%weno_fnl%initialize(weno=self%weno)
```

Do **not** create new `*_global` singletons. NASTO/CHASE/PATCH predate this refactor (they still reference the deleted singletons); **only PRISM is maintained against the current library** (CHASE verified not to build, 2026-09-23; its successor is FLUME, [issue #35](https://github.com/szaghi/adam/issues/35)).

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

**NVF/OpenACC modes require a varset.** `$NVF_CC` (the `-gpu=ccXX` compute
capability) is defined only in the NVF varsets, not the `local_gnu` default.
Build any `*-fnl-nvf` / `*-nvf*` mode with `--varset local_nvf` (sets
`$NVF_CC=cc89` for this workstation's Ada GPU); omitting it leaves `-gpu=$NVF_CC`
unexpanded and the build dies on the first nvfortran invocation. Other varsets:
`leonardo`, `iac` (both `cc80`). nvfortran must be on PATH first: `module load nvhpc` on WSL (26.9, the default,
matching `local_nvf`); on medooza `source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh` (see the varset table below).

```bash
# List all available build modes
fobis build --lmodes

# Build NASTO (Navier-Stokes solver)
fobis build --mode nasto-cpu-gnu                          # CPU with GNU compiler
fobis build --mode nasto-fnl-nvf --varset local_nvf       # OpenACC with NVF

# Build PRISM (Maxwell equations solver)
fobis build --mode prism-cpu-gnu                          # CPU with GNU compiler
fobis build --mode prism-fnl-nvf --varset local_nvf       # OpenACC with NVF

# Build other applications (NASTO, CHASE, PATCH are NOT maintained against the current library)
fobis build --mode chase-cpu-gnu                          # CHASE app (CPU/GNU) — does not build (drifted)
fobis build --mode patch-cpu-gnu                          # PATCH app (CPU/GNU)
fobis build --mode patch-fnl-nvf --varset local_nvf       # PATCH app (OpenACC/NVF)
fobis build --mode ascot-nvf-cuda --varset local_nvf      # ASCOT converter (CUDA Fortran)

# Build ADAM core libraries
fobis build --mode adam-com-gnu                           # Common library (GNU)
fobis build --mode adam-com-nvf --varset local_nvf        # Common library (NVF)

# Debug builds (append -debug to mode name)
fobis build --mode nasto-cpu-gnu-debug
fobis build --mode prism-fnl-nvf-debug --varset local_nvf
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
  - Build NVF modes with `--varset local_nvf` (or another NVF varset), which
    resolves both `$HDF5_PREFIX` and `$NVF_CC` — see "I/O and Serialization"

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
  - The library path is the fobos variable **`$HDF5_PREFIX`**, defined per
    `[varset:*]` in the repo `fobos` (there are no per-compiler
    `$HDF5_gnu`/`$HDF5_nvf` variables; the only `[common-variables]` entry is
    `$GPU_ARCH` for the AMD modes). Selecting a varset selects the HDF5 build:

    | Varset | Machine and toolchain | `$HDF5_PREFIX` | `$NVF_CC` |
    |---|---|---|---|
    | `local_gnu` (**default**) | WSL: gfortran 16.2.0 + Open MPI 5.0.11 | `lib/hdf5/develop/gnu/16.2.0` | — |
    | `local_nvf` | WSL: nvfortran 26.9 | `lib/hdf5/develop/nvf/26.9` | `cc89` |
    | `ci_gnu` | GitHub Actions: Ubuntu 24.04 apt gfortran-14 + Open MPI | `lib/hdf5/develop/gnu/14.2.0` | — |
    | `leonardo` | Spack path (NVHPC 24.5) | Spack | `cc80` |
    | `iac_gnu` | medooza: gfortran 16.2.0 + Open MPI 5.0.11 (`source /ifs/hpc/home/szaghi/gnu/gnu.sh`) | `lib/hdf5/develop/iac/gnu/16.2.0/openmpi-5.0.11` | — |
    | `iac_nvf` | medooza: nvfortran 26.9, CUDA 12.9 (`source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh`) | `lib/hdf5/develop/iac/nvf/26.9` | `cc80,cuda12.9` |
    | `spacehpc` | Cray PE | `/opt/cray/pe/hdf5/1.14.3.3/nvidia/23.3` | `cc90` |

  - **The toolchain must match the varset exactly** (since 2026-10-16, HDF5 2.3.0 develop `ab5df552` for the
    `local_*`/`iac_*` rows): Fortran `.mod` files work only with the compiler version that wrote them. A mismatch
    stops the build at `use hdf5` ("created by a different version of GNU Fortran") or, under nvfortran, as an
    internal compiler error ("mpi_f08_types ... member symbol not found"). Build clean after switching compilers.
    On medooza never use the cluster modules (nvhpc 25.5, gcc 12). CI runs `--varset ci_gnu` (jobs pinned to
    `ubuntu-24.04`) because its apt gfortran-14 cannot read the 16.2 HDF5.
  - The `local_*` HDF5 builds were made on Ubuntu 24.04 and need glibc ≥ 2.38: an `undefined reference to
    __isoc23_strtol` means the HDF5 was built on a newer glibc than the machine that links, so rebuild it there
    (the `iac_*` ones are built on medooza). With the WSL `cuda_13.3` NVHPC bundle, an older Windows driver needs the
    `cuda_multi` bundle and `-gpu=cuda12.9`; medooza's nodes run CUDA 12.9 only.

  - These are fobos *variables*, **not** shell environment variables —
    exporting them has no effect. Build a new machine by adding a `[varset:*]`
    block; see the [bring-up tutorial](docs/tests/prism-regression-tutorial.md).
  - `scripts/hdf5_build.sh -build` fetches and builds szip + zlib + parallel
    HDF5 with `CC=mpicc FC=mpif90` — load the intended compiler/MPI modules
    *before* running it. Full form: `scripts/hdf5_build.sh -build -sdk NVF|GNU -abs-paths -lsrc /tmp/hdf5-src
    -hdf5 $PWD/lib/hdf5/develop/<sdk>/<compiler-version>`, then point the varset's `$HDF5_PREFIX` at it. It builds
    the HDF5 `develop` branch (the repo's builds are commit `ab5df552`, 2.3.0); on Rocky/RHEL HDF5 2.3 installs into
    `lib64/`, so add `-DCMAKE_INSTALL_LIBDIR=lib` to the script's cmake line (ADAM expects `lib/`).
  - Collective I/O optimization: Set `H5Pset_dxpl_mpio` for MPI-IO hints

### Build System
- **FoBiS** (Fortran Building System): Manages complex multi-backend builds
  - Install: `pip install FoBiS.py`
  - Centralizes compiler flags, preprocessor macros, dependency tracking

### Cluster-Specific Considerations
- **medooza, CPU (gfortran 16.2 + Open MPI 5.0.11):** hcoll segfaults in `MPI_Finalize` on the compute nodes, so a
  correct run exits 139; `gnu.sh` exports `OMPI_MCA_coll_hcoll_enable=0` since 2026-10-17. The regression harnesses
  fail a case whose run exits non-zero, so a crash after the last checkpoint can no longer pass unseen.

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

## PRISM-Specific Gotchas

### Divergence correction: `hyperbolic` needs CT, `poisson` is refused, `physics` owns the `phi`/`psi` slots

The `phi`/`psi` GLM scalars used to have two sources of truth: `physics%nv` added their slots only under `divergence_correction = hyperbolic`, while `io_initialize` named them from the `numerics%constrained_transport_*` flags alone, so CT + `poisson` wrote 10–11 names into 9 slots (heap corruption, [issue #11](https://github.com/szaghi/adam/issues/11) → [#61](https://github.com/szaghi/adam/issues/61)). Current state:

- **`poisson` stops at initialization** (`adam_prism_numerics_object.F90`, since `68457200`): no backend implements it. Its code is dead, not removed: `impose_div_free` (CPU and FNL) dispatches only on `DIV_CORR_VAR_POISS`, which nothing assigns, and the FNL `impose_ct_correction_dev` is a stub (solve and gradient commented out, `dq_max` read uninitialised). Do not re-enable `poisson` without implementing both.
- **`hyperbolic` without CT stops at initialization** (#61): there is no scalar to carry, and the `physics` layout chain had no `else`, so `var_Jx/Jy/Jz` stayed undefined; `has_current_density()` then read `q(var_Jx:var_Jz)` at garbage indices and sent the run into the t0 elliptic B solve. `physics%initialize` also `error_stop`s on that combination, for any caller that bypasses `numerics`.
- **One source for the slots:** `physics%var_phi` / `physics%var_psi` (0 when absent) are set in the same branch that sets `nv_cl` and `var_J*`; `io_initialize` names φ/ψ from them and `error_stop`s if the names do not number `physics%nv`. The CPU flux setup still derives `fv_ivar_phi/psi` from `nv_c` in its own branch chain; it agrees by construction today.
- **CT with `divergence_correction = no` is ignored** with a warning (`numerics` clears both flags); 16 research inputs rely on that.

| `constrained_transport` | `divergence_correction` | Result |
|---|---|---|
| `no` | `no` | ✅ `nv = 9` |
| `D` / `B` / `DB` | `no` | ✅ CT ignored (warning), `nv = 9` |
| `D` / `B` | `hyperbolic` | ✅ `nv = 10` (φ or ψ in slot 7) |
| `DB` | `hyperbolic` | ✅ `nv = 11` (φ 7, ψ 8); regression case `rmf-hyperbolic` |
| `no` | `hyperbolic` | ⛔ refused at init |
| any | `poisson` | ⛔ refused at init |

When adding a state variable, follow the `rho` pattern below (one predicate drives width and names) or set an index in `physics` and derive the names from it.

### Forest substage state is per-realm (`stage_active`), NOT a shared module scalar

**History note:** an earlier design used a single program-scope `adam_forest_global%forest_active_substage` shared by every realm. That module scalar **no longer exists** (zero hits in `src/`); the P2 hazard it carried — "silently wrong under realm-local subcycling" — was resolved by promoting the substage marker to **per-realm state**: `realm_object%stage_active` (`adam_realm_object.F90:135`). `begin_stage_forest(k)` sets `self%stage_active = k`; `close_step_forest` clears it to 0.

The seam-fill peer-buffer selection reads `self%stage_active` **and** `peer%stage_active`: `> 0` selects that realm's substage buffer `rk%q_rk(:,...,stage_active)`, `0` selects committed `q` (the N=1 / end-of-step path). Because each realm carries its own marker, the machinery is already correct for realms at different substages — the constraint that survives is β **admissibility** (both endpoint realms must agree on `scheme_time`/`rk_scheme`/`nv`/K), enforced at init by `check_beta_admissibility`, not the old lockstep assumption. Cadence itself is now **per-seam** (`maps%seam_local_cadence(:)`, `CADENCE_END_OF_STEP`/`CADENCE_STAGE_COINCIDENT`), gated per-seam in `evolve_one_step`. When touching the seam-fill buffer selection, preserve the `stage_active` read/clear discipline (a mid-step write to `q` outside `begin/close_step_forest` corrupts peer reads). Related history: [issue #13](https://github.com/szaghi/adam/issues/13) P2.

### 2:1 AMR seam injects an O(h^p) div(B) source — no local collocated-FD fix exists (accept-truncation)

PRISM keeps `div(B) = 0` to machine precision on **uniform** grids: `compute_curl_fd_centered` and `compute_divergence_fd_centered` are both built from the same antisymmetric `FD1_CC` primitive, so `div_h(curl_h) ≡ 0` by operator commutation (a genuine mimetic / discrete-de-Rham identity, Ranocha 2019 Remark 2.8). **At a 2:1 AMR seam this breaks.** The coarse (Δx) and fine (Δx/2) sides difference with different-resolution stencils, the composite is no longer a tensor product, commutation fails, and the seam injects a `div(B)` source that lives in the grid-oscillation subspace `ker D*`. Measured on `rmf-amr-fd-pulse` (source-free, `-np1`, `ni=16`): seam `max|div(B)|` runs away **6.9 → 1002 over 100 steps at fixed h** (control = 0.0), while `‖B‖` stays flat (−3%) — it is a **divergence-constraint violation, not an energy instability**. The source is **O(h^p) truncation-order**: it *converges under refinement* (committed `--convergence` leg, `p_obs ≈ +1.15 → +2`) but is **unbounded in t at fixed h**, dt-independent, and scheme-universal (FD ≡ FV to 4-5 digits).

**No local collocated-FD fix exists — five fix classes ruled out with evidence (issue #29):**

| Fix class | What it tried | Verdict |
|---|---|---|
| **E1** | matched ghost VALUES (single-valued seam fill) | RED — 1D-normal-only, behaved as q=2 |
| **E2** | single-valued flux VALUES (auxiliary EMF) | NULL — fluxes already agree 0.24%; reframed the source to OPERATOR-not-value |
| **E3** | matched-difference seam OPERATOR (telescoping-exact) | CLOSED — eigenvalue-**unstable** (both-operator horizon → 1e18); the non-SBP construction Mattsson-Carpenter 2010 / Ranocha 2019 predict; telescoping-exactness ≠ time-stability |
| **E** | complete the Dedner GLM (add parabolic damping) | CLOSED for the **centred** GLM only — implemented + tunable `[physics].c_r` (since `68457200` replaced by the rate `[physics].alpha`, 1/s), c_r sweep 0.05→5.0 all still run away. **Reopened as E′ ([#99](https://github.com/szaghi/adam/issues/99)):** centred `grad_h` annihilates the π mode, so a seam error in the grid-oscillation subspace drives ψ into a mode that exerts no force on B (#68 F3); an upwinded or dissipative (B_n, ψ) coupling is untested, and FLUME's upwinded GLM keeps its seam div(B) bounded (#68 G3) |
| **A** | SBP-norm-compatible seam + weak SAT | CLOSED (analysis) — SBP bounds `‖B‖` (already fine), NOT `div(B)`; a value-transfer cannot restore operator commutation across the resolution jump (the `R∘P≠I` wall, #21) |

**RESOLUTION = accept-truncation, pending E′ (#99).** The O(h^p) refinement-convergent floor is accepted; its growth law versus h at a fixed long physical time is not yet measured (#68 G7, part of #99), so "refine h" is not a quantified mitigation. **Route 2** ([#30](https://github.com/szaghi/adam/issues/30): Stokes-form CT, B staged on faces while the FD cell update stays; #68 C7) is **parked** behind #99, activating only if E′ fails or a production workload exceeds tolerance. Guard-rail: `[IO].seam_divB_tol` (default off) + `[IO].seam_divB_error` — a run-time monitor in `save_divergence_history` (`adam_io_object.F90` / `adam_prism_common_object.F90`) that warns (or `error_stop`s) when seam `max|div(B)|` exceeds the tolerance on AMR runs, turning the otherwise-**silent** growth into a visible signal. **Mitigation to extend AMR runs:** enable `divergence_correction = hyperbolic` (needs `constrained_transport = D/B/DB`, refused without it) and tune the damping `[physics].alpha` — this *delays* the runaway (~5.5× margin, measured with the former `c_r`, not re-measured with `alpha`) but does **not** cure it. **Safe usage:** short/moderate-horizon AMR runs, or refine h, or accept the monitored floor. Tracked in [issue #29](https://github.com/szaghi/adam/issues/29).

### Inter-realm 1:1 mirror seam needs β cadence to stay div-free (issue #31, CLOSED)

Distinct from #29: a **1:1 same-resolution INTER-realm mirror seam** (two realms glued via a `[forest]` manifest, `coupling = mirror`) must be div-free to machine precision — it is numerically identical to a 1:1 intra-realm block interface, which is proven div-free. It is NOT a fundamental obstruction. Two defects made it leak (`~93 → 1600` div(B) over 100 steps on a source-free pulse), both fixed:

1. **α cadence (`coupling_cadence = end_of_step`) leaves the seam ghosts unfilled during RK substages.** `fill_seam_from_peer_forest` runs once per step (Phase 5) under α, so substages 2..N read *zero* seam ghosts in the stage buffer `q_rk(:,...,k)` → the curl evolves the seam skin from garbage. **For a 1:1 same-dt seam use `coupling_cadence = stage_coincident` (β)** — it fills every substage and stays div-free. α's end-of-step lag is a Berger-Oliger AMR-subcycling convention, not div-preserving for a 1:1 mirror seam. (α on a div-constrained 1:1 seam still leaks — a separate, harder item if α-div-preservation is ever required.)
2. **The div(B) diagnostic did not refill the inter-realm seam.** `post_step_forest` refreshed ghosts via `self%update_ghost`, which fills intra-realm ghosts + physical BCs but NOT the inter-realm seam (the forest owns it) → the divergence stencil read stale seam ghosts → a spurious div(B) even on a div-free field. Fixed by re-filling the seam from peers on committed `q` before `compute_divergence`, `present(realm)`-guarded (both CPU and FNL).

Reproducer / regression anchor: `src/tests/prism/regression/rmf-2realm-fd-pulse/` (β cadence, holds div(B)=div(D)=0). When touching `post_step_forest`, `fill_seam_from_peer_forest`, or any forest seam path, remember the diagnostic must refill the inter-realm seam and 1:1 seams want β. Tracked in [issue #31](https://github.com/szaghi/adam/issues/31).

### FNL buffered transposed copies: staging buffer extent must equal the device destination (issue #31)

History: the fWLayer field used to be copied host→device through FUNDAL's **buffered** transposed HtoD (`dev_memcpy_to_device(..., buf=self%buf_5D_R8P)`) with the **q-shaped** staging buffer (last dim `nv=9`, not `3`). The device destination then reaches an assumed-shape, **lbound-remapped dummy** `dst(bb(1,1):,…)`; the non-identity remap makes nvfortran materialise a **host copy-in temporary**, `c_loc(dst)` yields a host address and `cuMemcpyHtoDAsync` fails with `CUDA_ERROR_INVALID_VALUE` (at `np>1`). **Not a VRAM issue** (the WSL `free/total memory` print is `/dev/dxg` garbage; the fault is only visible under `compute-sanitizer`). Since `7e39318b` the fWLayer holds **no device state** and `copy_cpu_gpu` copies only `q`, whose buffer extent equals `q_gpu`'s (identity remap, no temporary) — the specific crash cannot recur, but **the rule stands for every buffered transposed copy**: pass `buf=` only when the staging buffer's shape equals the destination's; otherwise use whole-array `dev_assign_to_device(..., ij=[1,5])`. Device analogue of the gfortran assumed-shape pointer-section + explicit-lbound copy-temp trap.

`[fWLayer]` takes a **physical** `width` (real, required even when no layer is used); the per-face cell count is derived, `C_face = min(ni, ceiling(width/ds))` (`adam_prism_fWLayer_object.F90`). All regression cases are migrated, and so are the 48 research inputs that kept `C = 0` (#58, which also added the coil counts and `[fdv]` section they lacked); `src/tests/lint/check-removed-input-keys.sh`, run at the top of the PRISM sweep, rejects `[fWLayer] C` in any committed input; `rmf-fwl` (`width = 0.045`, active layer) was migrated and its FNL golden refreshed in `43981538`.

**Resolved (`bad09954`, 2026-07-28):** `rmf-amr` — the only `fv_centered` regression case — briefly failed its digest golden after `28625dbf` made coil current-density stamping **scheme-aware** (four call sites that unconditionally used `compute_curl_fd_centered` now dispatch on `fdv_scheme`, so FV cases stamp with `compute_curl_fv_centered`; before the fix an FV solve evolved an FD-built source). The change was adjudicated, not rubber-stamped: step 0 is bit-identical including all twelve `coil_*_j_vec_*` rows, divergence appears only at step 5 confined to evolved fields, magnitudes are coherent (max ratios 1.2–1.5), and CPU and FNL show the same `By` shift to 8 digits. Goldens refreshed on both backends. **Caveat on record:** the bump rests on the new path being consistent by construction, not on a measured accuracy win — confirming that needs a refinement study against a coil case with a known analytic vector potential.

### PML: `[PML]` input contract, and `CLASSIC_DIRECT` bypasses the RK-PML path

The perfectly-matched-layer absorbing boundary lives in `prism_pml_object` (+ `prism_rk_pml_object` for the auxiliary-variable time integration), with FNL device twins. Storage is **face-local and reduced** — only boundary blocks carry it, compacted face by face, four split variables per face ordered *electric first, magnetic second* (e.g. `q_pml_x_* = [psi_Ey_x, psi_Ez_x, psi_Hy_x, psi_Hz_x]`).

`[PML] PML_type` accepts five kinds (`adam_prism_pml_object.F90:229-241`), each with its **own required parameter set** — supplying the wrong set is a hard `error_stop`, not a default:

| `PML_type` | Accepted spellings | Required params (beyond `width` + face flags) |
|---|---|---|
| `NONE` | `NO`/`no`/`None`/`NONE`… | — (PML disabled) |
| `CLASSIC` | `PML`, `classic`, `STANDARD` | `gamma_max` (≥0), `gamma_exponent` (>0) |
| `CLASSIC_DIRECT` | `direct`, `CLASSIC-DIRECT` | `gamma_max`, `gamma_exponent` |
| `CFS` (CPML) | `cfs`, `CPML` | `gamma_max`, **`alfa_max`** (note the spelling), `k_max` (>0) |
| `BERMUDEZ` | `bermudez` | `beta` (≥0), `gamma_exponent` (>0) |

All kinds also require a **physical** `width` (>0, an `error_stop` otherwise) and the six face flags `{x,y,z}_{minus,plus}_layer`. Enabling the section with **all six flags false** is an `error_stop` ("enable at least one").

Two traps:

1. **An unrecognised `PML_type` is fatal** (`error_stop` listing the accepted spellings, `0259860f`); before that commit it silently disabled the PML. Values are passed through `strip_control` first, so CRLF input parses.
2. **`CLASSIC_DIRECT` is not an RK-integrated PML.** Every `rk_pml` call site is guarded `/= 'CLASSIC_DIRECT'` — `initialize`, `initialize_stages`, `compute_stage`, `assign_stage`, `update_q_pml` are all skipped (`adam_prism_cpu_object.F90:4096-4195`, `adam_prism_fnl_object.F90:657`). It instead applies direct per-face damping kernels (`apply_{x,y,z}_face_direct_damping_dev`). It exists **for testing** (`d64369b5`); treat it as an instrument, not the production path, and do not assume PML auxiliary state advances under it.

CPU time integration is implemented only for `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ` and `CFS` — anything else `error_stop`s at `adam_prism_cpu_object.F90:342-346`.

### PIC: `physical_model = PIC` widens `nv` by one (`rho`) — gated on ONE flag, unlike the `phi`/`psi` hazard

`[physics] physical_model` selects `EM` (default), `ADIM_EM` (non-dimensional), or `PIC`; an unknown string **does** `error_stop` here (`adam_prism_physics_object.F90:447`). Under `PIC`, `nv_pic = 1` and the state vector widens — `nv = nv_c + nv_s + nv_cl + nv_pic`, so a plain `PIC` run is `nv = 10`, carrying charge density `rho` as the extra slot (`:410-414`).

**This was never the issue #11 hazard**: `io_initialize`'s `add_rho` is gated on `physics%physical_model == PIC_PHYSICAL_MODEL`, the **same** predicate that set `nv_pic`, so the names array and the allocation cannot disagree. The `phi`/`psi` slots were gated on `numerics%constrained_transport_*` while `nv_cl` came from `divergence_correction`, two sources that disagreed under `poisson`; since #61 they follow `physics%var_phi/var_psi`. When adding a new state variable, keep one source for its width and its name.

PIC configuration lives in its own `[PIC]` section (`adam_prism_pic_object.F90`, 2114 lines — the largest PRISM common file):

- `problem_type` — `plasma` or `single_particle`; `plasma_domain` ∈ `Uniform_domain` / `Uniform_cilinder` / `Uniform_cell` (note the spelling **`cilinder`**), with `cilinder_{radius,length,x_center,y_center,z_center,axis}` when cylindrical, plus `plasma_density` and `neutral_fraction`.
- `scheme_time` — `LEAPFROG` or Runge-Kutta; each has its own object pair (`prism_leapfrog_pic_object` / `prism_rk_pic_object`) on both backends.
- `particle_weighting_model` / `current_weighting_model` / `field_weighting_model` — the full ladder (NGP, CIC, TSC, cubic, quartic, quintic, Gaussian for charge and current; 0D–5D, Gaussian and bspline for fields), every one mirrored as a `_dev` device routine in `prism_fnl_pic_object`. Gaussian kinds additionally read `sigma` and `cutoff_sigma`.
- `initialization` (`coherent` → elliptic-solver field init) ties PIC into the FLAIL elliptic solver (the t0 solves,
  `*_field_elliptic_initialization`). `elliptic_correction` is removed (#84 E3d): `.true.` stops the run, and the
  removed-keys lint rejects it in committed inputs.

Particle state is `q_pic(8, particle_number)` with `pic_fields(6, particle_number)`, allocated only under PIC (`adam_prism_common_object.F90:275-285`).

**Particle closure at the domain faces** (`adam_prism_pic_bc_object.F90`, component `pic%bc`, FNL mirror in `prism_fnl_pic_object`): derived from the field `[bc_*] type` by a total map (no input key of its own). `periodic` → `PIC_BC_PERIODIC`: a particle beyond a periodic face is located in the ghost halo of the block on that face (periodic ghosts are images of the opposite interior, so standard deposits fold back through `reduce_ghost_*` and gathers read the images) and its position is wrapped by ±L at the **end of the step only**, never inside RK stages (their positions combine linearly), and **after the end-of-step deposit, just before the migration**: wrapped earlier, a particle jumps by L to the blocks of another rank and leaves every local halo. Every other kind → `PIC_BC_LEGACY` + an init warning (particle outside the domain stops the run).

**Conserving currents with periodic faces** (totorica, esirkepov_Q<N>; londrillo and `esirkepov_tail_cleanup` stop): the matched current is built **particle by particle on the trajectory box and folded** (CPU `deposit_matched_current`, FNL `deposit_matched_current_box` on the virtual block), exploiting separability: `J_d = Q_d P_e1 P_e2 F_d^E` with every factor 1-D (prefix sum and Q on the longitudinal weight change, P on the transverse weights), zero outside the widened box by construction, folded modulo n on writing — the current through the periodic face (the k=0 mode a fold-then-prefix-sum would lose) is included. The stage currents are combined as the stage sources were (`combine_stage_source`, linear). The matched charge is built the same way (`deposit_matched_charge_box`). The coherent matched gather applies `P^T` circularly along periodic axes, and the ghosts of the filtered field beyond a periodic face are the periodic images (CPU `coherent_field_weighting` + `transposed_filter_axis(periodic=)`, FNL scatter of `gq_gpu`), so the gather stays the adjoint of the deposit there and the self-force stays cancelled. **Without periodic faces the previous line solves and block deposit run unchanged, bitwise** (the box path reproduces them to roundoff only; that was verified). Single-block PEC/Neumann results are therefore unchanged. The non-conserving charge-conservation diagnostic zeroes the normal current only in the physical boundary ghosts (`local_map_bc_crown`); it used to zero every block ghost, which broke its `R_cont` at block interfaces and periodic faces.

### `[bc_*] type`: accepted spellings, and PRISM `periodic` is library periodicity

The six face sections `[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]` each take a `type`. The parser (`adam_prism_bc_object.F90`) accepts exactly six spellings:

`extrapolation` · `Neumann` · `Dirichlet` · `Silver_Muller` · `periodic` · `PEC` (or `pec`)

Matching is case-sensitive; an unrecognised value is **fatal** (`case default` → `error_stop` naming section, value and accepted spellings, `0259860f`), values pass through `strip_control` (CRLF-safe; one public copy in `adam_prism_parameters` since #57), and `bc_type(6)` is default-initialised to `0_I4P`. `radiative` is parsed and **refused** (#57): the CPU stopped only at the first BC application, and the FNL kernel had no branch for it, so the face got no BC.

**PRISM `periodic` is the library's `BC_PERIODIC = -1`** (`adam_parameters.f90`), so `grid%set_bc_type` sets `is_ijk_periodic`, the tree builds periodic neighbours and the ghost exchange (`update_ghost_local/mpi`) fills periodic ghosts across blocks and ranks; periodic faces are not boundary rows of `local_map_bc_crown`, and `set_boundary_conditions` has no periodic branch. The two faces of an axis must be both periodic or neither (`error_stop`). Two traps fixed with it: `realm_object%initialize` builds the tree and maps before the BC types are known, and `make_neighborhood` never refreshes an old boundary neighbour, so `prism_common_object%initialize` calls `tree%rebuild_neighborhood` + `make_comm_local_maps_ghost_bc` right after `set_bc_type` (without it an unrefined single block kept its periodic faces as unfilled boundary rows); on an all-periodic domain `local_map_bc_crown` stays **unallocated** (the FNL maps read that as "no rows"), so the FLAIL callers pass a zero-row copy. FLUME calls `set_bc_type` at the same point and has the same unrefined single-block trap. Before this, PRISM used an app constant (`BC_PERIOD = 5`) and a same-block wrap, correct only with one block along the periodic axis. **`Dirichlet` writes zero**, not a prescribed value, and edge/corner crown rows (`fec > 6`) receive no BC at all.

Note the underscore spelling `Silver_Muller` — the hyphenated form appears only in directory names (`BC_tests/Silver_Muller/`), never as an accepted input value. `PEC` is the only kind with a lowercase alias.

There is a **second, derived** taxonomy for the elliptic solver: `build_elliptic_bc_types` maps each face BC into `ELL_BC_{DIRICHLET,PERIODIC,EXACT_OPEN,PEC}` (`Neumann` and `Silver_Muller` collapse to `EXACT_OPEN`; so would `radiative`, which the parser now refuses). Unlike the EM parser, **this mapping is total** — `BC_EXTRAPOLATION` and any unmapped value both `error_stop`. So a misconfigured face can pass EM setup silently and only fail later, when an elliptic solve (PIC initialisation, divergence correction) first touches it.

The t0 elliptic solves run on FLAIL since [#84](https://github.com/szaghi/adam/issues/84) E3b (`solve_time_zero` in
`impose_pic_fields_time_zero`): `ELL_BC_PERIODIC` → `FLAIL_BC_PERIODIC`, `ELL_BC_DIRICHLET` → `FLAIL_BC_ZERO_GHOSTS`,
`ELL_BC_PEC` → `FLAIL_BC_DIRICHLET` (φ, tangential A) or `FLAIL_BC_NEUMANN` (normal A), `ELL_BC_EXACT_OPEN` →
`FLAIL_BC_FREE_GHOSTS` (E3c: the free-space potential of the source in the ghosts, summed over the whole domain on rank 0;
the legacy sum was rank-local); a null axis takes the periodic kind, and a refined grid is refused (E7). Every PRISM build mode links an FFT
library for it (FFTW, or cuFFT on the NVIDIA FNL modes).

When adding a BC kind, update **all three** `select case`s (the EM parser, the elliptic map, `solve_time_zero`).

## Suspected, unverified defects (do not fix before verifying)

Recorded in [issue #36](https://github.com/szaghi/adam/issues/36); each needs a numerical confirmation test first:

- **S-1 — confirmed and fixed (2026-10-17).** The library LLF and HLLC-LM solvers added `f(1)*(v²+w²)/0.5` to an energy flux `ρ u_n H` that already holds the tangential kinetic energy: with `qL = qR` both returned `F(q)·n + 2ρu_n(v²+w²)` in energy (measured, x and y). The four lines are removed; both now return `F(q)·n` to round-off. Consumers: NASTO CPU (LLF; unmaintained) only — FLUME and PRISM call neither.
- **S-2** NASTO FNL eigenvectors: `uvw_r2` evaluates to `-u` for the y direction (`adam_nasto_fnl_cns_kernels.F90:103`); reported `|L·R − I| ≈ 1` in y.

Confirmed by code reading (not planned for a fix — CHASE is superseded by FLUME, #35): CHASE's characteristic projection is applied transposed and its y/z right-eigenvector matrix is singular; only 1-D x Sod inputs ever exercised it.

## Development Environment: WSL2 GPU+MPI Caveats

The primary GPU development box is **WSL2**. Its GPU+MPI stack is broken in several places a real cluster is not, so the FNL (OpenACC/nvfortran) backend runs there for **correctness checks only — never performance**. Do not quote WSL timings as benchmarks, and **never copy WSL MPI/UCX settings into a Slurm/cluster job script** (they are perf regressions or no-ops on real InfiniBand + GPUDirect RDMA). The workarounds live in the `nvhpc` modulefile and `src/tests/prism/regression/run-fnl-local.sh`:

- `UCX_RNDV_THRESH=inf` — **issue #12**: the FNL halo exchange (`update_ghost_mpi_gpu`) posts MPI directly on device-resident ghost buffers (GPU-direct, by design — there is **no host-staging fallback in the code**). With ≥2 ranks UCX moves them via the **rendezvous** protocol, whose device-memory transports (`cuda_copy`/`gdr`) cannot retrieve the GPU primary context through WSL's `/dev/dxg` shim → SIGABRT in `ucp_proto_rndv_send_start` at the first exchange. Forcing every message **eager** sidesteps rendezvous. On real GPUDirect RDMA this is the *opposite* of what you want (rendezvous is the fast path).
- `UCX_TLS=^cma` — CMA same-node shared-memory transport is broken/restricted under WSL (`ptrace_scope`).
- `OMPI_MCA_coll_hcoll_enable=0` — hcoll HW-offloaded collectives need Mellanox IB this box lacks.
- `UCX_MEMTYPE_CACHE=n` — WSL memtype detection is unreliable.

Full rationale and the WSL-vs-cluster contrast table: [HPC Optimization Guide](.claude/docs/hpc-optimization.md) → "Development Environment: WSL2 GPU+MPI Caveats". Cluster GPU-aware tuning is the inverse shape (see "Cluster-Specific Considerations" above, e.g. `UCX_TLS=rc_x,sm,cuda_copy,cuda_ipc`).
