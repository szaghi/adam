# ADAM Library — `src/lib`

The `src/lib` directory contains the **ADAM SDK**: the physics-agnostic, hardware-portable core upon which all ADAM solver applications are built. It is organised into four subdirectories — one portable common layer and three hardware backends — described in detail below. PRISM and FLUME build the common layer and the FNL backend; no maintained application builds NVF or GMP.

---

## Design principles

The library is structured around two orthogonal axes of variability:

- **Physics axis** — handled entirely by the application layer (`src/app`). The SDK is equation-agnostic: it knows nothing about Navier-Stokes, Maxwell, or any other PDE.
- **Hardware axis** — handled by backend selection at compile time. `_FNL` selects the GPU backend: OpenACC with nvfortran (`*-fnl-nvf` modes), or the same kernels as OpenMP offloading with amdflang (`*-fnl-omp-amd` modes, `DEV_OMP`). `_NVF` and `_GMP` select the CUDA Fortran and OpenMP-offloading libraries, which no PRISM or FLUME mode builds.

The result is a matrix of `N_equations × N_backends` combinations served by a single shared source tree.

---

## Directory overview

| Directory | Role | Parallelism |
|-----------|------|-------------|
| [`common/`](#common) | Portable core objects — grid, tree, field, numerics, I/O | CPU · MPI · OpenMP |
| [`fnl/`](#fnl) | GPU backend (FUNDAL-based), built by PRISM and FLUME | OpenACC or OpenMP target · MPI |
| [`nvf/`](#nvf) | CUDA Fortran GPU backend; no build mode builds it | CUDA Fortran · MPI |
| [`gmp/`](#gmp) | OpenMP target offloading backend; no build mode | OpenMP target · MPI |

---

## `common/`

Source: `src/lib/common`.

The common library provides all CPU-portable, physics-agnostic building blocks. The GPU objects do not extend them: they are separate types, initialised from the common objects passed as arguments. All modules are re-exported by the aggregate entry point `adam_common_library`.

### AMR data structures

| Module | Purpose |
|--------|---------|
| `adam_tree_object` | Octree/quadtree with Morton-order linearization. CPU-resident hash-map keyed on 64-bit Morton indices; provides O(1) block insertion, deletion, and neighbour lookup during AMR update phases. |
| `adam_tree_node_object` | Individual tree node: the Morton code, the owning rank, the block index in the field arrays (`block_index`), the codes of its 26 neighbours, and the bucket links (`next`, `previous`). |
| `adam_tree_bucket_object` | Bucket (linked-list) storage for the hash-map underlying `adam_tree_object`. |
| `adam_field_object` | Per-block geometry and metrics, ghost exchange (`update_ghost_local`, `update_ghost_mpi`), refinement/derefinement of block data, restart blocks. The state array itself, `q(nv, 1-ngc:ni+ngc, …, nb)`, is owned by the application; blocks are addressed by a compact index, not by Morton keys. |
| `adam_maps_object` | Index maps built from the tree: local and MPI ghost-cell maps, boundary-condition maps (`local_map_bc_crown`), block reordering (`local_map`) and the inter-realm seam maps. Regenerated only on topology changes. |
| `adam_grid_object` | Block-structured Cartesian grid geometry: cell coordinates, mesh spacing, domain bounds, and ghost-cell layout for a given AMR configuration. |
| `adam_amr_object` | AMR workflow driver: applies user-defined refinement markers, enforces 2:1 level balance across block faces, and orchestrates the mark → refine/coarsen → rebalance cycle. |
| `adam_adam_object` | `adam_object`, the composer of `grid`, `tree`, `field` and `maps`: AMR update (`amr_update`, `adapt`), load redistribution, restart files (`save_restart_files`, `load_restart_files`), slices and VTK output. |

### Numerical methods

| Module | Purpose |
|--------|---------|
| `adam_weno_object` | Weighted Essentially Non-Oscillatory (WENO) reconstruction: upwind orders 1, 3, 5, 7, 9 (`weno-u-N`) and centred orders 2, 4, 6, 8 (`weno-c-N`), Jiang–Shu or scale-invariant weights, ROR order reduction. |
| `adam_rk_object` | Explicit Runge-Kutta time integration: classical RK1–3, strong-stability-preserving SSP11/SSP22/SSP33/SSP54, and the 4-stage Yoshida symplectic scheme. |
| `adam_leapfrog_object` | Leapfrog time integrator with optional Robert–Asselin–Williams filter for long-time stability. |
| `adam_blanes_moan_object` | Blanes–Moan 4th- and 6th-order symplectic splitting integrator for Hamiltonian systems. |
| `adam_cfm_object` | Commutator-Free Magnus (CFM) 4th-, 6th-, and 8th-order exponential integrators for non-autonomous ODEs. |
| `adam_fdv_operators_library` | Finite difference/volume spatial operators: gradient, divergence, curl, and Laplacian at configurable order. |
| `adam_riemann_euler_library` | Riemann solvers for the Euler equations: exact, HLLC, HLLC-LM, HLLEM, local Lax–Friedrichs and two-shock. Used by NASTO only; FLUME has its own Riemann libraries. |

### Physics support

| Module | Purpose |
|--------|---------|
| `adam_eos_ic_object` | Ideal compressible gas equation of state (EOS): thermodynamic relations, speed of sound, primitive↔conservative variable conversions. |
| `adam_ib_object` | Immersed Boundary (IB) method. Computes signed-distance/eikonal fields and modifies fluxes at solid boundaries. Solids: `analytical_sphere`, `analytical_circle`, `analytical_rectangle`; a `file.off` definition is refused at load. |
| `adam_flail_object` | FLAIL, the elliptic solver library: a facade (`setup`, `solve`, `solve_packed`) over the FFT Poisson solver `adam_poisson_fft_object` (uniform grids; periodic, Dirichlet, Neumann, isolated faces) and the AMR multigrid `adam_poisson_mg_object` (FCG with a V-cycle, symmetric 2:1 composite operator). |

### Infrastructure

| Module | Purpose |
|--------|---------|
| `adam_mpih_object` | MPI wrapper: initialize/finalize, `error_stop`, barrier, rank-prefixed messages, timers. The ghost exchange lives in `adam_field_object`. |
| `adam_io_object` | `[IO]` configuration, history files (residuals, energy, divergence, charge conservation) and the XH5F field writer `save_field`: one HDF5 file per rank plus one XDMF index. Restart files are written by `adam_object`. |
| `adam_realm_object` | `realm_object`, the base of every application type: owns every library object as a component and declares the `_forest` contract. |
| `adam_forest_object` | `forest_object`, the time loop: single-realm fast path, staged path for 2:1 faces and multi-realm forests, runtime regridding (`regrid_forest`). |
| `adam_slices_object` | Extracts 2D planar slices from the 3D field for lightweight diagnostic output. |
| `adam_parameters` | Named integer constants for boundary condition types, AMR flags, and field-variable index conventions shared across all modules. |
| `adam_common_library` | Aggregate re-export module: `use adam_common_library` exposes the entire common API in a single `use` statement. |

---

## `fnl/`

Source: `src/lib/fnl`.

The FNL backend accelerates compute-intensive kernels with **OpenACC** (nvfortran) or, from the same sources, **OpenMP target offloading** (amdflang, `DEV_OMP`), relying on the [FUNDAL](https://github.com/szaghi/FUNDAL) library for GPU memory management and device-side data movement. Its objects are separate types, components of the application's FNL type, initialised from the realm's common objects after the common initialisation; only the Poisson solvers extend their common types.

The aggregate entry point is `adam_fnl_library`.

### Objects

| Module | Initialised from | Purpose |
|--------|------------------|---------|
| `adam_fnl_field_object` | `grid`, `field`, `maps` | Device copies of the field data (FUNDAL `dev_alloc`) and GPU-aware ghost-cell exchange. |
| `adam_fnl_weno_object` | `weno` | Device WENO coefficients and reconstruction. |
| `adam_fnl_rk_object` | `grid`, `field`, `rk` | Runge-Kutta stages on the device. |
| `adam_fnl_ib_object` | `grid`, `field`, `ib` | Device eikonal solve and IB flux correction. |
| `adam_fnl_maps_object` | `maps` (a component of `field_fnl_object`) | Device copies of the ghost-cell and BC maps. |
| `adam_fnl_mpih_object` | — | `mpih_fnl_object`, an alias of FUNDAL's MPI handler; the program-scope `mpih_fnl` also owns the device context. |
| `adam_fnl_poisson_fft_object` | extends `poisson_fft_object` | cuFFT solves, device-resident slab transposes. |
| `adam_fnl_poisson_mg_object` | extends `poisson_mg_object` | The multigrid's FCG iteration and V-cycle on the device. |

### Device kernel libraries

| Module | Purpose |
|--------|---------|
| `adam_fnl_field_kernels` | Device kernels (`!$acc` with `!$omp` twins) for field operations: transposition, gradient stencils, variable-index reductions. |
| `adam_fnl_weno_kernels` | OpenACC kernels for WENO interface reconstruction across all spatial directions. |
| `adam_fnl_rk_kernels` | OpenACC kernels for Runge-Kutta stage updates and residual accumulation. |
| `adam_fnl_ib_kernels` | OpenACC kernels for eikonal distance solve and IB flux modification. |
| `adam_fnl_fdv_operators_library` | OpenACC device kernels for gradient, divergence, curl, and Laplacian operators. |

---

## `nvf/`

Source: `src/lib/nvf`.

The NVF backend uses **CUDA Fortran** (NVIDIA HPC SDK `nvfortran`). No build mode builds it: its only consumer is NASTO NVF (`adam_nasto_nvf_object.F90`), which has no mode, and `ascot-nvf-cuda` compiles ASCOT, which uses PENF only. Its objects are separate types, except `mpih_nvf_object`, which extends `mpih_object`.

The aggregate entry point is `adam_nvf_library`.

### Objects

| Module | Counterpart | Purpose |
|--------|-------------|---------|
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

## `gmp/`

Source: `src/lib/gmp`.

The GMP backend targets **OpenMP target offloading** (OpenMP 5.x `!$omp target` directives). No build mode uses it: AMD GPUs are served by the FNL sources compiled as OpenMP offloading (`*-fnl-omp-amd`). Its objects are separate types, except `mpih_gmp_object`, which extends `mpih_object`.

The aggregate entry point is `adam_gmp_library`.

### Objects

| Module | Counterpart | Purpose |
|--------|-------------|---------|
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

ADAM is dual-licensed under the MIT License and the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0) (LGPLv3) (`LICENSE.mit.md`, `LICENSE.lgpl3.md`); you may choose either.

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
