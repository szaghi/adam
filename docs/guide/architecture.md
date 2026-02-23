# Architecture

## Overview

ADAM is a **physics-agnostic SDK** for building high-performance CFD solvers. Its core infrastructure — block-structured AMR, immersed boundary method, high-order WENO numerics, Runge-Kutta time integration, and parallel I/O — is fully decoupled from any specific set of governing equations. Solvers are assembled by composing these building blocks and adding only the physics-specific layer on top.

ADAM targets the full spectrum of modern HPC hardware without changing application source code or input files:

- **CPU-based clusters** — MPI distributed-memory parallelism with shared-memory OpenMP threading
- **CPU+GPU accelerated clusters** — node-level GPU parallelism via CUDA Fortran (NVIDIA), OpenACC, or OpenMP offloading; multi-node scaling via MPI with GPU-aware communication over NVLink or InfiniBand

The choice of hardware backend is a compile-time switch; everything above the backend layer — physics, numerics, I/O, configuration — is identical.

```mermaid
graph TD
    A[Solver applications<br/>nasto · prism · chase · patch] --> B[ADAM SDK — src/lib/common<br/>AMR · IB · WENO · RK · I/O · MPI]
    A --> C{Hardware backend}
    C -->|CPU| CPU[CPU backend<br/>MPI + OpenMP]
    C -->|_NVF| D[NVF backend<br/>CUDA Fortran]
    C -->|_FNL| E[FNL backend<br/>OpenACC]
    C -->|_GMP| F[GMP backend<br/>OpenMP offloading]
    B --> G[Third-party libraries<br/>PENF · StringiFor · FiNeR · HDF5]
    CPU --> G
    D --> G
    E --> G
    F --> G
```

## SDK Layer — `src/lib`

The SDK provides physics-agnostic building blocks reused identically by every application.

### Core objects (`src/lib/common`)

| Object | Purpose |
|--------|---------|
| `adam_grid_object` | Block-structured grid management with AMR geometry |
| `adam_tree_object` | Octree/quadtree with Morton-order linearization |
| `adam_field_object` | 5D field arrays `(nv, ni, nj, nk, nb)` — storage, interpolation, ghost exchange |
| `adam_weno_object` | High-order WENO reconstruction (orders 3–11) |
| `adam_rk_object` | Runge-Kutta temporal integration (SSP schemes) |
| `adam_ib_object` | Immersed boundary method with eikonal distance fields |
| `adam_io_object` | Parallel HDF5 output and restart files |
| `adam_mpih_object` | MPI wrapper and nearest-neighbor ghost cell communication |
| `adam_fdv_operators_library` | Gradient, divergence, curl, Laplacian finite difference operators |
| `adam_riemann_euler_library` | Riemann solvers for the Euler equations |

### Backend libraries

Each backend extends the common objects with hardware-specific implementations:

| Directory | Backend | Parallelism model |
|-----------|---------|-------------------|
| `src/lib/common` | CPU | MPI + OpenMP |
| `src/lib/nvf` | NVF | CUDA Fortran (NVIDIA GPUs) |
| `src/lib/fnl` | FNL | OpenACC (NVIDIA/AMD GPUs) |
| `src/lib/gmp` | GMP | OpenMP target offloading (experimental) |

## Application Layer — `src/app`

Applications sit on top of the SDK and contribute only the physics-specific layer. The full HPC stack is inherited for free.

### Directory structure

```
src/
├── lib/                  # ADAM SDK
│   ├── common/           # Physics-agnostic core objects (portable, CPU)
│   ├── nvf/              # CUDA Fortran GPU backend
│   ├── fnl/              # OpenACC GPU backend
│   └── gmp/              # OpenMP offloading backend (in development)
├── app/                  # Solvers built on the SDK
│   ├── nasto/            # Compressible Navier-Stokes solver
│   ├── prism/            # Maxwell equations / plasma solver
│   ├── chase/            # CFD application
│   ├── patch/            # Patch-based application
│   └── ascot/            # Binary-to-ASCII output converter
├── tests/                # Unit and integration tests
└── third_party/          # Git submodules (PENF, StringiFor, FiNeR, VTKFortran, …)
```

### Backend pattern

Every application exposes the same set of backends via a parallel subdirectory layout:

```
app/<name>/common/    # Physics layer shared across all backends
app/<name>/cpu/       # CPU-only entry point (MPI + OpenMP)
app/<name>/nvf/       # CUDA Fortran entry point
app/<name>/fnl/       # OpenACC entry point
app/<name>/gmp/       # OpenMP offloading entry point
```

### Adding a new solver

A new physics application requires only implementing the problem-specific layer; the entire SDK is reused unchanged:

```fortran
type :: my_solver_object
   ! ---- reused from src/lib — zero extra work ----
   type(mpih_object)  :: mpih     ! MPI handler
   type(grid_object)  :: grid     ! AMR block grid
   type(field_object) :: field    ! conservative variables
   type(amr_object)   :: amr      ! refinement markers
   type(ib_object)    :: ib       ! solid bodies
   type(weno_object)  :: weno     ! spatial reconstruction
   type(rk_object)    :: rk       ! time integration
   ! ---- only this part is new ----
   type(my_physics_object) :: physics
   type(my_bc_object)      :: bc
   type(my_ic_object)      :: ic
   type(my_io_object)      :: io
end type
```

## AMR data design: inverse indexing

High-performance AMR requires resolving a fundamental tension: the grid topology changes dynamically at runtime (refinement, coarsening, load rebalancing), yet the numerical kernels must operate on dense, contiguous memory with no indirection overhead. ADAM resolves this by splitting the AMR data into two structurally opposite objects with complementary roles.

### Tree — flexible topology on the CPU

`adam_tree_object` is a **hash-map living entirely in CPU memory**. Its keys are 64-bit Morton indices that linearise the four octree coordinates `(level, bx, by, bz)` into a single integer; its values are lightweight block descriptor objects. This structure provides:

- **O(1) insertion and deletion** of blocks during refinement and coarsening steps
- **O(1) neighbour lookup** — the Morton key of any face/edge/corner neighbour can be computed arithmetically, with no pointer chasing
- **Spatial locality** — Morton ordering clusters geometrically adjacent blocks in index space, minimising ghost-cell communication volume across MPI ranks

The tree is never touched by numerical kernels. It is only consulted during the AMR update phase (marking → refinement/coarsening → load rebalancing) and to regenerate the index maps needed by the field object.

### Field — contiguous arrays for parallel computing

`adam_field_object` is a **dense, contiguous 5D array** allocated once at the beginning of the simulation:

```
field%q(nv, ni, nj, nk, nb)
```

| Dimension | Meaning |
|-----------|---------|
| `nv` | Physical variables (density, momenta, energy, …) |
| `ni, nj, nk` | Cell indices within a block, including ghost cells |
| `nb` | Block index — a compact integer from 1 to the current block count |

The block index `nb` is **not** the Morton key. It is a compact sequential integer assigned so that all blocks owned by an MPI rank occupy a contiguous slice of the array. This is the inverse of the tree's hash-map addressing — hence *inverse indexing*: the tree maps Morton keys → block descriptors, while the field maps compact block indices → raw data.

This layout guarantees:
- **Stride-1 access on the innermost dimension** (`nv`) in Fortran column-major order, enabling coalesced reads across CUDA threads or OpenACC/OpenMP SIMD lanes
- **No dynamic allocation during time integration** — the array is pre-allocated to the maximum block count and reused across all Runge-Kutta stages
- **Direct offload to device memory** — a single `!$acc data` or `cudaMemcpy` transfers the entire field; no pointer-based scatter/gather is needed

### The mapping layer

A lightweight mapping array (`maps`) bridges the two worlds. It is regenerated only when the AMR topology changes (a rare, synchronised event) and is otherwise invisible to numerical kernels:

```
maps%b2m(nb)   ! block index → Morton key  (field → tree lookup)
maps%m2b(key)  ! Morton key  → block index (tree  → field lookup)
```

During the computation phase — which accounts for the overwhelming majority of runtime — kernels iterate over the compact block index `nb` with no hash-map access, achieving the same memory access pattern as a static structured grid.

```mermaid
graph LR
    T["adam_tree_object<br/>(CPU hash-map)<br/>Morton key → block descriptor<br/>AMR topology, O(1) insert/delete"]
    M["maps object<br/>b2m / m2b arrays<br/>regenerated on AMR update only"]
    F["adam_field_object<br/>(contiguous 5D array)<br/>q(nv, ni, nj, nk, nb)<br/>MPI · OpenACC · CUDA · OpenMP"]
    T <-->|"AMR update phase"| M
    M <-->|"index translation"| F
    F -->|"computation phase<br/>(no tree access)"| F
```

### Memory and parallelism summary

- **Ghost cells**: configurable width `ngc` (typically 3 for WENO5, 4 for WENO7); exchanged via MPI before each stencil sweep using the compact block index, not Morton keys.
- **Load balancing**: blocks are redistributed across MPI ranks by reordering the compact index using Morton space-filling curves, keeping the field array layout optimal after each rebalancing step.
- **Scalability target**: strong scaling to O(1000) GPUs with >70% parallel efficiency.

## Configuration

All applications are configured through human-readable INI files (parsed by FiNeR). No recompilation is needed to change physics, numerics, grid, or I/O settings:

- Grid parameters — domain bounds, resolution, ghost cell width
- Physics parameters — gas properties, Reynolds/Mach numbers
- Numerical parameters — WENO order, Runge-Kutta scheme
- Boundary and initial conditions
- I/O options and AMR refinement criteria
