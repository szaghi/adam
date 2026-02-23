# ADAM NVF Library — `src/lib/nvf`

The NVF backend is ADAM's **CUDA Fortran GPU acceleration layer** targeting NVIDIA GPUs via the `nvfortran` compiler (NVIDIA HPC SDK). It follows the same two-tier pattern as the FNL backend — a **wrapper object** extending the common CPU class plus a **kernel module** containing the device-side computation — but replaces OpenACC directives with explicit CUDA Fortran constructs: `allocatable, device` arrays, `!$cuf kernel do()` launches, `dim3` grid/block dimensions, and direct `cudafor` runtime calls.

No physics or algorithmic logic is duplicated from `src/lib/common`. All equations, coefficients, and data structures are defined once in the common layer and mirrored to the GPU by the NVF layer.

The aggregate entry point `adam_nvf_library` re-exports the entire NVF API together with `adam_common_library`; a single `use adam_nvf_library` statement in application code exposes both layers.

---

## Contents

- [NVF vs FNL at a glance](#nvf-vs-fnl-at-a-glance)
- [Memory management library](#memory-management-library)
- [MPI handler and CUDA device management](#mpi-handler-and-cuda-device-management)
- [Field](#field)
- [WENO reconstruction](#weno-reconstruction)
- [Runge-Kutta integration](#runge-kutta-integration)
- [Immersed boundary method](#immersed-boundary-method)
- [Communication maps](#communication-maps)
- [Module summary](#module-summary)

---

## NVF vs FNL at a glance

| Aspect | FNL (OpenACC) | NVF (CUDA Fortran) |
|--------|--------------|-------------------|
| Device arrays | OpenACC `declare device_resident` | `allocatable, device` |
| Kernel launch | `!$acc parallel loop` | `!$cuf kernel do(N) <<<*,*>>>` |
| Grid/block dims | Compiler-chosen | Manual `dim3` via `set_cuda_dimensions` |
| Synchronisation | Implicit at data regions | Explicit `cudaDeviceSynchronize()` |
| Memory alloc | `dev_alloc` (FUNDAL) | `alloc_var_gpu` (local memory library) |
| Host-device copy | `dev_assign_to_device` (FUNDAL) | `assign_allocatable_gpu` (local memory library) |
| Device procedures | `!$acc routine seq` | `attributes(device)` |
| Error checking | Implicit | Explicit `check_cuda_error` |
| CUDA device binding | Not applicable | `CudaSetDevice` per MPI rank |
| Vendor portability | OpenACC (multi-vendor) | NVIDIA only |

---

## Memory management library

### `adam_memory_nvf_library` — GPU allocation and transfer utilities

Centralises all GPU memory operations with optional verbose reporting and `cudaMemGetInfo` bookkeeping. All wrapper objects call these routines rather than issuing raw CUDA allocations.

**Generic interface `alloc_var_gpu`** — allocate a device array:

| Specialisation | Type | Rank |
|----------------|------|------|
| `alloc_var_gpu_R8P_1D` … `_6D` | `real(R8P)` | 1 – 6 |
| `alloc_var_gpu_I4P_1D`, `_5D` | `integer(I4P)` | 1, 5 |
| `alloc_var_gpu_I8P_1D` … `_3D` | `integer(I8P)` | 1 – 3 |

All specialisations share the same signature pattern:
```fortran
subroutine alloc_var_gpu_*(var, ulb, msg, verbose)
  ! var     — allocatable, device: deallocated then reallocated to requested bounds
  ! ulb     — lower/upper bounds (reshaped for multi-dimensional arrays)
  ! msg     — optional label for memory-usage log lines
  ! verbose — if .true., report free/total GPU memory via cudaMemGetInfo before and after
```

**Generic interface `assign_allocatable_gpu`** — copy host array to device:

| Specialisation | Notes |
|----------------|-------|
| `assign_allocatable_gpu_R8P_1D` … `_4D` | Standard copy |
| `assign_allocatable_gpu_R8P_2D` (transposed variant) | Supports in-place 2D transposition during copy |
| `assign_allocatable_gpu_I4P_1D`, `_1D_rhs_allocated`, `_5D` | Integer variants |
| `assign_allocatable_gpu_I8P_2D`, `_3D` | 64-bit index maps |

Only copies if the source is allocated and non-empty; no-op otherwise.

**Utility `save_memory_gpu_status(file_name, tag)`** — appends a free/total GPU memory snapshot to a text file; useful for tracking device memory usage across initialisation stages.

**Interface `transpose_a`** — standalone 2D array transposition helper (`transpose_a_R8P_2D`).

---

## MPI handler and CUDA device management

### `adam_mpih_nvf_object` — extends `mpih_object` with GPU device binding

This is the only NVF object that uses Fortran type extension (`type, extends(mpih_object)`). It adds GPU device selection, CUDA error propagation, and kernel launch dimension management on top of the common MPI wrapper.

**Additional fields:**

| Field | Type | Purpose |
|-------|------|---------|
| `mydev` | `I4P` | GPU device index assigned to this MPI rank |
| `local_comm` | `I4P` | Intra-node communicator used for device assignment |
| `iercuda` | `I4P` | Most recent CUDA error code |
| `cblk(3)` | `I4P` | Default CUDA thread-block dimensions `(x=32, y=8, z=1)` |
| `grid`, `tBlock` | `dim3` | Current kernel launch configuration |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(do_mpi_init, do_device_init, …)` | Calls parent `%initialize`, then uses `MPI_COMM_SPLIT_TYPE(MPI_COMM_TYPE_SHARED)` to assign a unique GPU per intra-node rank via `CudaSetDevice(mydev)`. Queries `cudaGetDeviceProperties` and stores total device memory in `memory_avail`. |
| `set_cuda_dimensions(cgrd, cblk)` | Computes `tBlock = dim3(cblk)` and `grid = dim3(⌈cgrd/tBlock⌉)` for a subsequent kernel launch. |
| `check_cuda_error(error_code, msg)` | Calls `cudaDeviceSynchronize` then `cudaGetLastError`; aborts with a formatted message on any CUDA error. |
| `print_device_properties(device_properties)` | Logs total memory, shared memory, warp size, compute capability, SM count, L2 cache size, and clock rate. |
| `load_from_file(file_parameters, go_on_fail)` | Reads `[cuda]` section from INI file: `block_x`, `block_y`, `block_z` → `cblk(1:3)`. |

INI configuration section: `[cuda]`.

---

## Field

### `adam_field_nvf_object` — GPU field wrapper

Holds all device-side field and geometry arrays. Unlike FNL, the CPU and GPU layouts use the **same index order** `(nv, ni, nj, nk, nb)` — no permanent layout transposition. A separate transposed scratch array `q_t_gpu` is used only during host↔device copies.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_gpu` | `(nv, 1-ngc:ni+ngc, …, nb)` | Primary conservative field |
| `q_t_gpu` | same shape | Scratch for transposed GPU↔CPU copy |
| `x_cell_gpu`, `y_cell_gpu`, `z_cell_gpu` | `(nb, ni)` / `(nb, nj)` / `(nb, nk)` | Cell centroid coordinates |
| `dxyz_gpu` | `(3, nb)` | Mesh spacing `(dx, dy, dz)` per block |
| `fec_1_6_array_gpu` | `(nb)` | Face enumeration code for IB ghost-cell lookup |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(field, nv_aux, verbose)` | Allocate all device arrays via `alloc_var_gpu`; allocate CPU transposition buffer `q_t`; call `copy_cpu_gpu`. |
| `copy_cpu_gpu(verbose)` | Transfer coordinate and spacing arrays and communication maps to device via `assign_allocatable_gpu`. |
| `copy_transpose_cpu_gpu(nv, q_cpu, q_gpu)` | Transpose and upload: `q_cpu(nv,…,nb)` → `q_gpu(nv,…,nb)` via loop-based transposition on the host. |
| `copy_transpose_gpu_cpu(nv, q_gpu, q_cpu)` | Download and transpose: calls `copy_transpose_gpu_cpu_cuf` kernel, then transfers `q_t_gpu` to host. |
| `update_ghost_local_gpu(q_gpu)` | Intra-rank ghost update; calls `update_ghost_local_gpu_cuf`. |
| `update_ghost_mpi_gpu(q_gpu, step)` | Three-step asynchronous MPI ghost exchange: (1) pack send buffer on GPU, (2) post `MPI_IRECV`/`MPI_ISEND`, (3) `MPI_WAITALL` and unpack. The optional `step` argument enables pipelining with computation. |
| `compute_q_gradient(b, ivar, q_gpu, gradient)` | AMR criterion: max `‖∇q‖/‖q‖` over block `b`; calls `compute_q_gradient_cuf`. |

### `adam_field_nvf_kernels` — field CUDA kernels

All kernels use `!$cuf kernel do(N) <<<*,*>>>` with dynamic grid/block sizing and are followed by `!@cuf iercuda=cudaDeviceSynchronize()`.

| Kernel | Directive | Purpose |
|--------|-----------|---------|
| `compute_q_gradient_cuf` | `do(3)` + `reduce(max:)` | Centred-difference gradient magnitude for AMR marking |
| `compute_normL2_residuals_cuf` | `do(4)` + `reduce(+:)` | Per-variable L2 norm `√(Σ dq²)` for convergence monitoring |
| `copy_transpose_gpu_cpu_cuf` | `do(4)` | `q_t_gpu(v,i,j,k,b) ← q_gpu(v,i,j,k,b)` (index reorder on device) |
| `populate_send_buffer_ghost_gpu_cuf` | `do(1)` | Pack ghost cells into MPI send buffer; `mode=1` direct copy, `mode=8` 8-cell AMR average |
| `receive_recv_buffer_ghost_gpu_cuf` | `do(1)` | Unpack MPI receive buffer into ghost cells |
| `update_ghost_local_gpu_cuf` | `do(2)` | Intra-rank block-to-block ghost update with AMR coarse↔fine averaging |

---

## WENO reconstruction

### `adam_weno_nvf_object` — GPU WENO coefficient wrapper

Extends no common type directly (holds a pointer to `weno_object`). Mirrors all coefficient arrays computed once on the CPU to device memory.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `a_gpu` | `(2, 0:S-1, S)` | Optimal WENO weights per sub-stencil and interface |
| `p_gpu` | `(2, 0:S-1, 0:S-1, S)` | Polynomial reconstruction coefficients |
| `d_gpu` | `(0:S-1, 0:S-1, 0:S-1, S)` | Smoothness indicator coefficients |
| `ror_schemes_gpu` | `(:)` | ROR fallback scheme orders near solid walls |
| `ror_ivar_gpu` | `(:)` | Variable indices checked by ROR |
| `ror_stats_gpu` | `(:,:,:,:,:)` | ROR statistics counters |
| `cell_scheme_gpu` | `(nb, ni, nj, nk)` | Per-cell effective reconstruction order |

`initialize(weno)` copies all arrays via `assign_allocatable_gpu`.

### `adam_weno_nvf_kernels` — WENO device procedures

All procedures carry `attributes(device)` — they are callable only from within a running CUDA kernel, not from the host.

| Procedure | Visibility | Purpose |
|-----------|-----------|---------|
| `weno_reconstruct_upwind_device(S, a, p, d, zeps, V, VR)` | public | Top-level reconstruction: given stencil `V(2, 1-S:S-1)`, return left/right interface values `VR(2)` |
| `weno_compute_polynomials_device(S, p, V, VP)` | private | `VP(f,s1) = Σ_{s2} p(f,s2,s1,S) · V(f, s1-s2)` |
| `weno_compute_weights_device(S, a, d, zeps, V, w)` | private | Smoothness indicators → normalised weights via `w(f,k) = a(f,k,S)/(zeps+IS)² / Σ` |
| `weno_compute_convolution_device(S, VP, w, VR)` | private | `VR(f) = Σ_k w(f,k) · VP(f,k)` |

The three private procedures are the standard WENO algorithm decomposition; `weno_reconstruct_upwind_device` calls them in sequence. Equation solvers invoke this from inside their spatial flux loops.

---

## Runge-Kutta integration

### `adam_rk_nvf_object` — GPU RK stage manager

Holds a pointer to the common `rk_object` for scheme metadata and allocates all stage storage on the device.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_rk_gpu` | `(nv, ni, nj, nk, nb, nrk_stage)` | Stage values; `nrk_stage=1` for low-storage, `nrk_stage=nrk` for SSP |
| `alph_gpu` | `(nrk, nrk-1)` | SSP alpha coefficients |
| `beta_gpu` | `(nrk)` | SSP beta coefficients |
| `gamm_gpu` | `(nrk)` | SSP gamma coefficients |

**Scheme allocation strategy:**

| Scheme | `nrk_stage` | Mode |
|--------|------------|------|
| `RK_1`, `RK_2`, `RK_3` | 1 | Low-storage: overwrites `q_rk_gpu` in place |
| `RK_SSP_22`, `RK_SSP_33` | 2 / 3 | Multi-stage |
| `RK_SSP_54` | 5 | Multi-stage |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(rk, nb, ngc, ni, nj, nk, nv)` | Allocate coefficient and stage arrays sized to scheme requirements |
| `initialize_stages(q_gpu)` | Broadcast `q_gpu` into all stage slots via `rk_initialize_stages_cuf` |
| `assign_stage(s, q_gpu, phi_gpu)` | Copy `q_gpu` into stage `s`, skipping solid cells (if `phi_gpu` present) |
| `compute_stage(s, dt, phi_gpu)` | SSP accumulation: `q_rk(:,s) += dt·α(s,ss)·q_rk(:,ss)` for `ss=1…s−1` |
| `compute_stage_ls(s, dt, phi_gpu, dq_gpu, q_gpu)` | Low-storage update: `q = ark·q_n + brk·q + dt·crk·dq` |
| `update_q(dt, phi_gpu, q_gpu)` | Final assembly: `q += dt·β(s)·q_rk(:,s)` for `s=1…nrk` |

### `adam_rk_nvf_kernels` — RK CUDA kernels

All kernels use `!$cuf kernel do(N) <<<*,*>>>` over 5–6D loop nests and mask solid cells via `phi_gpu`.

| Kernel | `do(N)` | Purpose |
|--------|---------|---------|
| `rk_assign_stage_cuf` | `do(5)` | `q_rk(:,s) ← q_gpu` (fluid cells only) |
| `rk_initialize_stages_cuf` | `do(6)` | `q_rk(:,all_s) ← q_gpu` |
| `rk_compute_stage_cuf` | `do(6)` | `q_rk(:,s) += dt·α(s,ss)·q_rk(:,ss)` |
| `rk_compute_stage_ls_cuf` | `do(5)` | `q = ark·q_n + brk·q + dt·crk·dq` (scalar args passed by `value`) |
| `rk_update_q_cuf` | `do(6)` | `q += dt·β(s)·q_rk(:,s)` for `s=1…nrk` |

Scalar time-step and coefficient arguments are passed with `value` intent to avoid global memory reads per thread.

---

## Immersed boundary method

### `adam_ib_nvf_object` — GPU IB wrapper

Holds pointers to the common `ib_object` and to `field_nvf_object` (for grid metrics). Manages the signed-distance field `phi_gpu` on device.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `phi_gpu` | `(nb, 1-ngc:ni+ngc, …, n_solids+1)` | Signed-distance function; last slice = `max` over all solids |
| `q_bcs_vars_gpu` | `(:,:)` | Wall boundary condition reference state per solid |

Sign convention: `phi < 0` inside solid (ghost region), `phi > 0` in fluid.

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(ib, field_gpu)` | Allocate `phi_gpu` (initialised to −1) and `q_bcs_vars_gpu`; associate grid metadata pointers |
| `evolve_eikonal(dq_gpu, q_gpu)` | For each solid: compute gradient-weighted residual `dq`, then apply `q -= dq` inside solid |
| `invert_eikonal(q_gpu)` | Enforce wall BC at surface (φ > 0): `BCS_VISCOUS` mirrors all velocity components; `BCS_EULER` reflects the normal component |

### `adam_ib_nvf_kernels` — IB CUDA kernels

All spatial kernels use `do(4)` over `(k,j,i,b)` with an inner variable loop.

| Kernel | Purpose |
|--------|---------|
| `compute_phi_analytical_sphere_cuf` | `φ = −(‖x−xc‖ − R)` — negative inside sphere |
| `compute_phi_all_solids_cuf` | `φ_all = max(φ₁,…,φ_ns)` — union-of-solids mask |
| `compute_eikonal_dq_phi_cuf` | 1st-order upwind residual: `dq = |n_x|·Δq_x + |n_y|·Δq_y + |n_z|·Δq_z` inside solid; normal direction from centred-difference gradient of `φ`, scaled by 0.9 |
| `evolve_eikonal_q_phi_cuf` | `q -= dq` where `φ > 0` |
| `invert_eikonal_q_phi_cuf` | Viscous: `(u,v,w)→(−u,−v,−w)`; Euler: `u → u − 2(u·n̂)n̂` |
| `move_phi_cuf` | Level-set advection `∂φ/∂t = −v·∇φ` for moving bodies |
| `reduce_cell_order_phi_cuf` | Lower WENO order in cells within `ib_reduction_extent` of the surface, per spatial direction |

---

## Communication maps

### `adam_maps_nvf_object` — GPU maps wrapper

Mirrors all AMR block-to-block and MPI ghost-cell communication index tables to device memory so that buffer packing and unpacking run entirely on the GPU.

**Device arrays:**

| Array | Shape columns | Content |
|-------|--------------|---------|
| `local_map_ghost_cell_gpu` | 9 | `(b_src, b_dst, i_src, j_src, k_src, i_dst, j_dst, k_dst, mode)` |
| `comm_map_send_ghost_cell_gpu` | 7 | `(b_src, i, j, k, v_offset, buf_idx, mode)` |
| `comm_map_recv_ghost_cell_gpu` | 6 | `(buf_idx, b_dst, i, j, k, v_offset)` |
| `send_buffer_ghost_gpu` | — | 1D packed MPI send staging buffer |
| `recv_buffer_ghost_gpu` | — | 1D packed MPI receive staging buffer |
| `local_map_bc_crown_gpu` | — | Boundary condition crown ghost-cell map |

Index arrays use `integer(I8P)` to accommodate large AMR block counts. `mode=1` denotes a one-to-one cell correspondence; `mode=8` denotes an 8-cell average at an AMR refinement interface.

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(maps)` | Store pointer to common maps, call `copy_cpu_gpu(verbose=.true.)` |
| `copy_cpu_gpu(verbose)` | Transfer all map and buffer arrays via `assign_allocatable_gpu` |

---

## Module summary

| Module | Role | Extends / wraps |
|--------|------|----------------|
| `adam_nvf_library` | Aggregate entry point | — |
| `adam_memory_nvf_library` | GPU alloc/copy utilities (`alloc_var_gpu`, `assign_allocatable_gpu`) | — |
| `adam_mpih_nvf_object` | CUDA device binding, kernel dims, error checking | `mpih_object` |
| `adam_field_nvf_object` | GPU field wrapper + host↔device transfer | `field_object` (pointer) |
| `adam_field_nvf_kernels` | Gradient, L2 norm, transpose, ghost-cell pack/unpack | — |
| `adam_weno_nvf_object` | GPU WENO coefficient mirror + ROR tables | `weno_object` (pointer) |
| `adam_weno_nvf_kernels` | `attributes(device)` reconstruction procedures | — |
| `adam_rk_nvf_object` | GPU RK stage storage and update dispatch | `rk_object` (pointer) |
| `adam_rk_nvf_kernels` | Stage assign, accumulate, low-storage, final update | — |
| `adam_ib_nvf_object` | GPU distance field + eikonal BC wrapper | `ib_object` (pointer) |
| `adam_ib_nvf_kernels` | Eikonal evolution, sphere distance, momentum inversion | — |
| `adam_maps_nvf_object` | GPU communication index tables + MPI buffer staging | `maps_object` (pointer) |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
