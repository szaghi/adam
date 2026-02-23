# ADAM Library — `src/lib`

The `src/lib` directory contains the **ADAM SDK**: the physics-agnostic, hardware-portable core upon which all ADAM solver applications are built. It is organised into four subdirectories — one portable common layer and three hardware backends — described in detail below.

---

## Design principles

The library is structured around two orthogonal axes of variability:

- **Physics axis** — handled entirely by the application layer (`src/app`). The SDK is equation-agnostic: it knows nothing about Navier-Stokes, Maxwell, or any other PDE.
- **Hardware axis** — handled by backend selection at compile time. A single preprocessor flag (`_NVF`, `_FNL`, or `_GMP`) switches between CUDA Fortran, OpenACC, and OpenMP offloading without touching application code.

The result is a matrix of `N_equations × N_backends` combinations served by a single shared source tree.

---

## Directory overview

| Directory | Role | Parallelism |
|-----------|------|-------------|
| [`common/`](#common--src-lib-common) | Portable core objects — grid, tree, field, numerics, I/O | CPU · MPI · OpenMP |
| [`fnl/`](#fnl--src-lib-fnl) | OpenACC GPU backend (FUNDAL-based) | OpenACC · MPI |
| [`nvf/`](#nvf--src-lib-nvf) | CUDA Fortran GPU backend | CUDA Fortran · MPI |
| [`gmp/`](#gmp--src-lib-gmp) | OpenMP target offloading backend *(experimental)* | OpenMP target · MPI |

---

## `common/` — `src/lib/common`

The common library provides all CPU-portable, physics-agnostic building blocks. Every GPU backend extends these classes rather than replacing them. All modules are re-exported by the aggregate entry point `adam_common_library`.

### AMR data structures

| Module | Purpose |
|--------|---------|
| `adam_tree_object` | Octree/quadtree with Morton-order linearization. CPU-resident hash-map keyed on 64-bit Morton indices; provides O(1) block insertion, deletion, and neighbour lookup during AMR update phases. |
| `adam_tree_node_object` | Individual tree node: stores the Morton key, refinement level, spatial coordinates `(bx, by, bz)`, and pointers to parent/children/neighbours. |
| `adam_tree_bucket_object` | Bucket (linked-list) storage for the hash-map underlying `adam_tree_object`. |
| `adam_field_object` | Dense 5D field array `q(nv, ni, nj, nk, nb)` allocated once at simulation start. Uses a compact sequential block index (`nb`) rather than Morton keys, enabling stride-1 GPU-coalesced access and zero-scatter device offload. |
| `adam_maps_object` | Mapping layer between tree and field: `b2m(nb)` translates block index → Morton key and `m2b(key)` does the inverse. Regenerated only on AMR topology changes; invisible to numerical kernels during time integration. Drives ghost-cell exchange and field interpolation during refinement/coarsening. |
| `adam_grid_object` | Block-structured Cartesian grid geometry: cell coordinates, mesh spacing, domain bounds, and ghost-cell layout for a given AMR configuration. |
| `adam_amr_object` | AMR workflow driver: applies user-defined refinement markers, enforces 2:1 level balance across block faces, and orchestrates the mark → refine/coarsen → rebalance cycle. |
| `adam_adam_object` | Top-level composer that binds `adam_tree_object` and `adam_field_object` into a single object, coordinating AMR updates and CPU↔GPU synchronisation. |

### Numerical methods

| Module | Purpose |
|--------|---------|
| `adam_weno_object` | Weighted Essentially Non-Oscillatory (WENO) finite difference reconstruction, orders 3 through 11. Provides smooth and shock-capturing reconstruction of interface values from cell-centred data. |
| `adam_rk_object` | Explicit Runge-Kutta time integration. Supports classical RK1–3, strong-stability-preserving SSP22/SSP33/SSP54, and Yoshida symplectic schemes. |
| `adam_leapfrog_object` | Leapfrog time integrator with optional Robert–Asselin–Williams filter for long-time stability. |
| `adam_blanes_moan_object` | Blanes–Moan 4th- and 6th-order symplectic splitting integrator for Hamiltonian systems. |
| `adam_cfm_object` | Commutator-Free Magnus (CFM) 4th-, 6th-, and 8th-order exponential integrators for non-autonomous ODEs. |
| `adam_fdv_operators_library` | Finite difference/volume spatial operators: gradient, divergence, curl, and Laplacian at configurable order. |
| `adam_riemann_euler_library` | Riemann solvers for the Euler equations: exact solver, HLLC, HLL, Local Lax–Friedrichs, Roe, and first-order Godunov flux. |

### Physics support

| Module | Purpose |
|--------|---------|
| `adam_eos_ic_object` | Ideal compressible gas equation of state (EOS): thermodynamic relations, speed of sound, primitive↔conservative variable conversions. |
| `adam_ib_object` | Immersed Boundary (IB) method for complex and moving geometries. Computes signed-distance/eikonal fields and modifies fluxes at solid boundaries. Supports spheres, cylinders, and user-defined shapes. |
| `adam_flail_object` | Linear algebra interface: smoothing operators, successive over-relaxation (SOR), and multigrid building blocks for implicit solvers. |

### Infrastructure

| Module | Purpose |
|--------|---------|
| `adam_mpih_object` | MPI wrapper: process topology, nearest-neighbour ghost-cell communication using compact block indices, non-blocking send/receive staging. |
| `adam_io_object` | Parallel HDF5 output (field snapshots) and restart files via MPI-IO. Handles dataset chunking and collective I/O hints. |
| `adam_slices_object` | Extracts 2D planar slices from the 3D field for lightweight diagnostic output. |
| `adam_parameters` | Named integer constants for boundary condition types, AMR flags, and field-variable index conventions shared across all modules. |
| `adam_common_library` | Aggregate re-export module: `use adam_common_library` exposes the entire common API in a single `use` statement. |

---

## `fnl/` — `src/lib/fnl`

The FNL backend accelerates compute-intensive kernels with **OpenACC**, relying on the [FUNDAL](https://github.com/szaghi/FUNDAL) library for GPU memory management and device-side data movement. It extends the common objects via Fortran type extension; all non-kernel logic is inherited unchanged.

The aggregate entry point is `adam_fnl_library`.

### Extended objects

| Module | Extends | Purpose |
|--------|---------|---------|
| `adam_fnl_field_object` | `adam_field_object` | Adds `!$acc declare device_resident` annotations and GPU-aware ghost-cell exchange. |
| `adam_fnl_weno_object` | `adam_weno_object` | Dispatches WENO reconstruction to OpenACC device kernels. |
| `adam_fnl_rk_object` | `adam_rk_object` | Drives Runge-Kutta stages entirely on device, avoiding CPU round-trips. |
| `adam_fnl_ib_object` | `adam_ib_object` | GPU-accelerated eikonal solve and IB flux correction. |
| `adam_fnl_maps_object` | `adam_maps_object` | Ghost-cell packing/unpacking and field interpolation on device. |
| `adam_fnl_mpih_object` | `adam_mpih_object` | Thin alias to FUNDAL's MPI handler; enables GPU-aware MPI transfers. |

### Device kernel libraries

| Module | Purpose |
|--------|---------|
| `adam_fnl_field_kernels` | OpenACC kernels for field operations: transposition, gradient stencils, variable-index reductions. |
| `adam_fnl_weno_kernels` | OpenACC kernels for WENO interface reconstruction across all spatial directions. |
| `adam_fnl_rk_kernels` | OpenACC kernels for Runge-Kutta stage updates and residual accumulation. |
| `adam_fnl_ib_kernels` | OpenACC kernels for eikonal distance solve and IB flux modification. |
| `adam_fnl_fdv_operators_library` | OpenACC device kernels for gradient, divergence, curl, and Laplacian operators. |

---

## `nvf/` — `src/lib/nvf`

The NVF backend uses **CUDA Fortran** (NVIDIA HPC SDK `nvfortran`) for maximum low-level control over thread blocks, shared memory, and warp-level primitives on NVIDIA GPUs. It follows the same extension pattern as FNL.

The aggregate entry point is `adam_nvf_library`.

### Extended objects

| Module | Extends | Purpose |
|--------|---------|---------|
| `adam_field_nvf_object` | `adam_field_object` | Manages `device` attribute arrays and CPU↔GPU copy scheduling via `cudafor`. |
| `adam_weno_nvf_object` | `adam_weno_object` | Launches WENO CUDA kernels with tuned grid/block dimensions. |
| `adam_rk_nvf_object` | `adam_rk_object` | Drives Runge-Kutta stages on device with stream-based kernel overlap. |
| `adam_ib_nvf_object` | `adam_ib_object` | CUDA-accelerated eikonal solver and IB flux correction. |
| `adam_maps_nvf_object` | `adam_maps_object` | Ghost-cell exchange and field interpolation via CUDA kernels. |
| `adam_mpih_nvf_object` | `adam_mpih_object` | CUDA-aware MPI: direct GPU-buffer send/receive without staging through host memory. |

### Device kernel libraries and utilities

| Module | Purpose |
|--------|---------|
| `adam_field_nvf_kernels` | CUDA Fortran kernels for field operations (transposition, stencil helpers). |
| `adam_ib_nvf_kernels` | CUDA kernels for eikonal distance computation and IB forcing. |
| `adam_memory_nvf_library` | GPU memory allocation/deallocation wrappers (`cudaMalloc`, `cudaFree`, `cudaMemcpy`) presented as Fortran interfaces. |

---

## `gmp/` — `src/lib/gmp`

The GMP backend targets **OpenMP target offloading** (OpenMP 5.x `!$omp target` directives), enabling portability to non-NVIDIA accelerators (AMD, Intel GPUs) without vendor-specific extensions. It is currently experimental.

The aggregate entry point is `adam_gmp_library`.

### Extended objects

| Module | Extends | Purpose |
|--------|---------|---------|
| `adam_field_gmp_object` | `adam_field_object` | Manages `omp_target_alloc` device arrays and host↔device transfers. |
| `adam_weno_gmp_object` | `adam_weno_object` | Dispatches WENO reconstruction to OpenMP target regions. |
| `adam_rk_gmp_object` | `adam_rk_object` | Runge-Kutta stages executed inside `!$omp target` regions. |
| `adam_ib_gmp_object` | `adam_ib_object` | OpenMP target kernels for eikonal solve and IB forcing. |
| `adam_maps_gmp_object` | `adam_maps_object` | Ghost-cell communication and interpolation via OpenMP target. |
| `adam_mpih_gmp_object` | `adam_mpih_object` | Extended MPI handler with OpenMP target device management. |

### Device kernel libraries and utilities

| Module | Purpose |
|--------|---------|
| `adam_field_gmp_kernels` | OpenMP target kernels for field operations. |
| `adam_rk_gmp_kernels` | OpenMP target kernels for Runge-Kutta stage updates. |
| `adam_ib_gmp_kernels` | OpenMP target kernels for eikonal and IB forcing. |
| `adam_memory_gmp_library` | Device memory management wrappers over `omp_target_alloc` / `omp_target_free` / `omp_target_memcpy`. |
| `adam_gmp_utils` | Low-level OpenMP target utilities: allocation guards, memcpy helpers, and device-presence queries presented as clean Fortran interfaces. |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
