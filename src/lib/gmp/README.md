# ADAM GMP Library — `src/lib/gmp`

The GMP backend is ADAM's **OpenMP target offloading GPU acceleration layer**, providing vendor-neutral GPU support through OpenMP 5.0+ `!$omp target` directives. It follows the same two-tier pattern as the FNL and NVF backends — a **wrapper object** managing device-resident arrays plus a **kernel module** containing the offloaded computation — but replaces all vendor-specific constructs with standard OpenMP semantics, making it portable to NVIDIA, AMD, and Intel GPU hardware.

> **Status**: experimental. The GMP backend is fully structured and compiles, but has not yet reached the production maturity level of the FNL and NVF backends.

No physics or algorithmic logic is duplicated from `src/lib/common`. All equations, coefficients, and data structures are defined once in the common layer and mirrored to the device by the GMP layer.

The aggregate entry point `adam_gmp_library` re-exports the entire GMP API together with `adam_common_library`; a single `use adam_gmp_library` exposes both layers.

---

## Contents

- [GMP vs FNL vs NVF at a glance](#gmp-vs-fnl-vs-nvf-at-a-glance)
- [OpenMP target directives used](#openmp-target-directives-used)
- [Memory management utilities](#memory-management-utilities)
- [MPI handler and device management](#mpi-handler-and-device-management)
- [Field](#field)
- [WENO reconstruction](#weno-reconstruction)
- [Runge-Kutta integration](#runge-kutta-integration)
- [Immersed boundary method](#immersed-boundary-method)
- [Communication maps](#communication-maps)
- [Module summary](#module-summary)

---

## GMP vs FNL vs NVF at a glance

| Aspect | GMP (OpenMP target) | FNL (OpenACC) | NVF (CUDA Fortran) |
|--------|--------------------|--------------|--------------------|
| Standard | OpenMP 5.0+ | OpenACC 3.x | CUDA Fortran (nvfortran) |
| Vendor portability | NVIDIA, AMD, Intel | NVIDIA, AMD (limited) | NVIDIA only |
| Device arrays | `pointer` + `omp_target_alloc` | `declare device_resident` / OpenACC data | `allocatable, device` |
| Kernel launch | `!$omp target teams distribute parallel do` | `!$acc parallel loop gang vector` | `!$cuf kernel do(N) <<<*,*>>>` |
| Parallelism levels | 3 — teams / distribute / parallel do | 3 — gang / worker / vector | 2 — grid / block |
| Loop collapsing | Explicit `collapse(N)` | Explicit `gang vector` | Implicit in `do(N)` |
| Device memory alloc | `omp_target_alloc` (C interop) | FUNDAL `dev_alloc` | `alloc_var_gpu` (local library) |
| Host↔device copy | `omp_target_memcpy` (C interop) | FUNDAL `dev_assign_to_device` | `assign_allocatable_gpu` (local library) |
| GPU-aware MPI | `!$omp target data use_device_addr` | `!$acc host_data use_device` | Direct (CUDA-aware MPI) |
| Device selection | `omp_set_default_device` | Implicit | `CudaSetDevice` |
| Error checking | Runtime return codes | Implicit | Explicit `check_cuda_error` |
| Preprocessor guards | None — fully runtime | `#ifdef _FNL` | `#ifdef _NVF` |
| Production status | Experimental | Production | Production |

---

## OpenMP target directives used

The GMP backend uses a small, consistent set of OpenMP offloading constructs throughout all kernel modules.

### `!$omp target teams distribute parallel do`

The primary kernel directive. Creates a three-level parallel hierarchy:

```fortran
!$omp target teams distribute parallel do collapse(4) &
!$omp&    has_device_addr(q_gpu, phi_gpu) reduction(max:gradient)
do b = 1, blocks_number
  do k = 1, nk
    do j = 1, nj
      do i = 1, ni
        ...
      end do
    end do
  end do
end do
!$omp end target teams distribute parallel do
```

- `teams` — creates a league of thread teams (equivalent to CUDA grid)
- `distribute` — distributes outer loop iterations across teams
- `parallel do` — parallelises remaining iterations within each team
- `collapse(N)` — fuses N nested loops into a single parallel iteration space

### `has_device_addr(array, …)`

Declares that listed arrays are already resident on the device (allocated via `omp_target_alloc`). Prevents the runtime from generating implicit host↔device copies that would otherwise corrupt device pointers.

### `reduction(op:var)`

Standard OpenMP reduction applied across the parallel region:
- `reduction(max:gradient)` — gradient magnitude for AMR marking
- `reduction(+:norm_gpu)` — L2 norm accumulation

### `!$omp target data use_device_addr(array)`

Opens a data region in which the listed array's **device address** is exposed to host code. Used exclusively for GPU-aware MPI calls:

```fortran
!$omp target data use_device_addr(send_buffer_ghost_gpu)
  call MPI_ISEND(send_buffer_ghost_gpu(start), count, MPI_DOUBLE, ...)
!$omp end target data
```

### `map(from:array)`

Transfers an array from device to host at the end of a target region. Used in `copy_transpose_gpu_cpu_gmp` to return the transposed field after the device-side index reorder.

---

## Memory management utilities

### `adam_gmp_utils` — low-level OpenMP runtime wrappers

Thin Fortran wrappers over the OpenMP C runtime functions, accessed via `iso_c_binding`. All higher-level allocation in the GMP backend routes through these.

**Generic interface `omp_target_alloc_f`** — allocate a device pointer:

| Specialisation | Type | Rank |
|----------------|------|------|
| `omp_target_alloc_R8P_1D` … `_6D` | `real(R8P)` | 1 – 6 |
| `omp_target_alloc_I4P_1D`, `_2D`, `_5D` | `integer(I4P)` | 1, 2, 5 |
| `omp_target_alloc_I8P_1D` … `_3D` | `integer(I8P)` | 1 – 3 |

```fortran
subroutine omp_target_alloc_R8P_1D(fptr_dev, ubounds, omp_dev, ierr, lbounds, init_value)
  real(R8P), pointer,       intent(out) :: fptr_dev(:)
  integer(I4P),             intent(in)  :: ubounds(1)
  integer(I4P), optional,   intent(in)  :: omp_dev       ! defaults to omp_get_default_device()
  integer(I4P), optional,   intent(in)  :: lbounds(1)    ! defaults to 1
  real(R8P),    optional,   intent(in)  :: init_value
  integer(I4P),             intent(out) :: ierr
```

When `init_value` is supplied, the allocated region is initialised via an offloaded loop:
```fortran
!$omp target teams distribute parallel do has_device_addr(fptr_dev)
do i = lbounds(1), ubounds(1)
  fptr_dev(i) = init_value
end do
```

**Generic interface `omp_target_free_f`** — deallocate a device pointer. Matching ranks and types to `omp_target_alloc_f`.

**Generic interface `omp_target_memcpy_f`** — copy between host and device:

| Specialisation | Type |
|----------------|------|
| `omp_target_memcpy_R8P` | `real(R8P)` |
| `omp_target_memcpy_I4P` | `integer(I4P)` |
| `omp_target_memcpy_I8P` | `integer(I8P)` |

```fortran
function omp_target_memcpy_R8P(fptr_dst, fptr_src, dst_off, src_off, omp_dst_dev, omp_src_dev)
  integer(I4P)                        :: omp_target_memcpy_R8P   ! return code
  real(R8P), target, intent(out)      :: fptr_dst(..)
  real(R8P), target, intent(in)       :: fptr_src(..)
  integer(I4P),      intent(in)       :: omp_dst_dev, omp_src_dev, dst_off, src_off
```

### `adam_memory_gmp_library` — high-level allocation and transfer

Wraps `adam_gmp_utils` with bounds checking, optional verbose logging, and a consistent interface matching the NVF `adam_memory_nvf_library`.

**Generic interface `alloc_var_gpu`** — allocate with error checking:

Same type/rank coverage as `omp_target_alloc_f`. Aborts with `error stop` on allocation failure.

```fortran
subroutine alloc_var_gpu_R8P_1D(var, ulb, omp_dev, init_val, msg, verbose)
  real(R8P), pointer,      intent(inout) :: var(:)
  integer(I4P),            intent(in)    :: ulb(2)      ! (lower, upper) bounds
  integer(I4P),            intent(in)    :: omp_dev
  real(R8P),    optional,  intent(in)    :: init_val
  character(*), optional,  intent(in)    :: msg
  logical,      optional,  intent(in)    :: verbose
```

**Generic interface `assign_allocatable_gpu`** — copy host allocatable to device pointer:

| Specialisation | Notes |
|----------------|-------|
| `assign_allocatable_gpu_R8P_1D` … `_4D` | Standard copy |
| `assign_allocatable_gpu_R8P_2D` (transposed) | CPU transpose before copy |
| `assign_allocatable_gpu_I4P_1D`, `_5D` | Integer variants |
| `assign_allocatable_gpu_I8P_2D`, `_3D` | 64-bit index maps |

No-op if the source is unallocated or empty.

---

## MPI handler and device management

### `adam_mpih_gmp_object` — extends `mpih_object` with OpenMP device binding

```fortran
type, extends(mpih_object) :: mpih_gmp_object
  integer(I4P) :: mydev      = 0  ! GPU device index for this MPI rank
  integer(I4P) :: myhos      = 0  ! Host device ID (from omp_get_initial_device)
  integer(I4P) :: local_comm = 0  ! Intra-node MPI communicator
end type
```

**`initialize(do_mpi_init, do_device_init)`** — overrides parent:
1. Calls `mpih_object%initialize`
2. If `do_device_init=.true.`:
   - Creates intra-node communicator via `MPI_COMM_SPLIT_TYPE(MPI_COMM_TYPE_SHARED)`
   - Derives local rank `mydev` from the shared communicator
   - Binds the OpenMP default device: `call omp_set_default_device(self%mydev)`
   - Records the host device ID: `self%myhos = omp_get_initial_device()`

Each MPI rank is assigned a unique GPU by its intra-node rank index, exactly as in the NVF backend.

---

## Field

### `adam_field_gmp_object` — GPU field wrapper

Holds all device-side field and geometry arrays as Fortran pointers to `omp_target_alloc`-allocated device memory.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_gpu` | `(nv, 1-ngc:ni+ngc, …, nb)` | Primary conservative field |
| `q_t_gpu` | same | Scratch for transposed GPU↔CPU copy |
| `x_cell_gpu`, `y_cell_gpu`, `z_cell_gpu` | `(nb, ni/nj/nk)` | Cell centroid coordinates |
| `dxyz_gpu` | `(3, nb)` | Mesh spacing `(dx, dy, dz)` per block |
| `fec_1_6_array_gpu` | `(nb)` | Face enumeration code for IB ghost-cell lookup |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(mpih, field, nv_aux, verbose)` | Allocate all device arrays via `alloc_var_gpu`; allocate CPU transposition buffer `q_t`; call `copy_cpu_gpu` |
| `copy_cpu_gpu(verbose)` | Transfer coordinate and spacing arrays to device via `assign_allocatable_gpu` |
| `copy_transpose_cpu_gpu(nv, q_cpu, q_gpu)` | Transpose and upload: loop-based CPU transposition then `omp_target_memcpy_f` |
| `copy_transpose_gpu_cpu(nv, q_gpu, q_cpu)` | Download and transpose: calls `copy_transpose_gpu_cpu_gmp` kernel (device reorder + `map(from:)`) |
| `update_ghost_local(q_gpu)` | Intra-rank ghost update; calls `update_ghost_local_gmp` |
| `update_ghost_mpi(q_gpu, step)` | Three-step asynchronous MPI ghost exchange; device buffers exposed to MPI via `!$omp target data use_device_addr` |
| `compute_q_gradient(b, ivar, q_gpu, gradient)` | AMR criterion: max `‖∇q‖` over block `b`; calls `compute_q_gradient_gmp` |

### `adam_field_gmp_kernels` — field OpenMP kernels

All kernels use `!$omp target teams distribute parallel do` with `has_device_addr` on every device pointer argument.

| Kernel | `collapse(N)` | Purpose |
|--------|--------------|---------|
| `compute_q_gradient_gmp` | `collapse(3)` + `reduction(max:)` | Centred-difference gradient magnitude for AMR marking |
| `compute_normL2_residuals_gmp` | `collapse(4)` + `reduction(+:)` | Per-variable L2 norm `√(Σ dq²)` |
| `copy_transpose_gpu_cpu_gmp` | `collapse(4)` + `map(from:)` | `q_t_gpu(v,i,j,k,b) ← q_gpu(b,i,j,k,v)` then transfer to host |
| `populate_send_buffer_ghost_gmp` | `collapse(1)` | Pack ghost cells into MPI send buffer; `mode=1` direct copy, `mode=8` 8-cell AMR average |
| `receive_recv_buffer_ghost_gmp` | `collapse(1)` | Unpack MPI receive buffer into ghost cells |
| `update_ghost_local_gmp` | `collapse(1)` | Intra-rank block-to-block ghost update with AMR coarse↔fine averaging |

---

## WENO reconstruction

### `adam_weno_gmp_object` — GPU WENO coefficient wrapper

Holds a pointer to the common `weno_object` and mirrors all coefficient arrays to device memory. No kernel module — WENO reconstruction is called from within equation-solver kernels via device-resident coefficient pointers.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `a_gpu` | `(2, 0:S-1, S)` | Optimal WENO weights per sub-stencil and interface |
| `p_gpu` | `(2, 0:S-1, 0:S-1, S)` | Polynomial reconstruction coefficients |
| `d_gpu` | `(0:S-1, 0:S-1, 0:S-1, S)` | Smoothness indicator coefficients |
| `ror_schemes_gpu` | `(:)` | ROR fallback scheme orders near solid walls |
| `ror_ivar_gpu` | `(:)` | Variable indices checked by ROR |
| `ror_stats_gpu` | `(:,:,:,:,:)` | ROR statistics counters |
| `cell_scheme_gpu` | `(nb, ni, nj, nk, 3)` | Per-cell effective reconstruction order per direction |

`initialize(weno)` copies all arrays via `assign_allocatable_gpu`.

---

## Runge-Kutta integration

### `adam_rk_gmp_object` — GPU RK stage manager

Holds a pointer to the common `rk_object` for scheme metadata and allocates stage storage on the device.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_rk_gpu` | `(nb, ni, nj, nk, nv, nrk_stage)` | Stage values; `nrk_stage=1` for low-storage, `nrk_stage=nrk` for SSP |
| `alph_gpu` | `(nrk, nrk)` | SSP alpha coefficients |
| `beta_gpu` | `(nrk)` | SSP beta coefficients |
| `gamm_gpu` | `(nrk)` | SSP gamma coefficients |

**Scheme allocation strategy:**

| Scheme | `nrk_stage` | Mode |
|--------|------------|------|
| `RK_1`, `RK_2`, `RK_3` | 1 | Low-storage — overwrites stage in place |
| `RK_SSP_22`, `RK_SSP_33` | 2 / 3 | Multi-stage |
| `RK_SSP_54` | 5 | Multi-stage |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(rk, nb, ngc, ni, nj, nk, nv)` | Allocate coefficient and stage arrays sized to scheme |
| `initialize_stages(q_gpu)` | Broadcast `q_gpu` into all stage slots |
| `assign_stage(s, q_gpu, phi_gpu)` | Copy `q_gpu` into stage `s`, skipping solid cells |
| `compute_stage(s, dt, phi_gpu)` | SSP accumulation: `q_rk(:,s) += dt·α(s,ss)·q_rk(:,ss)` |
| `compute_stage_ls(s, dt, phi_gpu, dq_gpu, q_gpu)` | Low-storage update: `q = ark·q_n + brk·q + dt·crk·dq` |
| `update_q(dt, phi_gpu, q_gpu)` | Final assembly: `q += dt·β(s)·q_rk(:,s)` for `s=1…nrk` |

### `adam_rk_gmp_kernels` — RK OpenMP kernels

All kernels use `!$omp target teams distribute parallel do` with `has_device_addr` on device arrays. Solid cells are skipped when `phi_gpu` is present (`phi_gpu(b,i,j,k,all_solids) < 0`).

| Kernel | `collapse(N)` | Purpose |
|--------|--------------|---------|
| `rk_assign_stage_gmp` | `collapse(5)` | `q_rk(:,s) ← q_gpu` (fluid cells only) |
| `rk_initialize_stages_gmp` | `collapse(6)` | `q_rk(:,all_s) ← q_gpu` |
| `rk_compute_stage_gmp` | `collapse(6)` | `q_rk(:,s) += dt·α(s,ss)·q_rk(:,ss)` |
| `rk_compute_stage_ls_gmp` | `collapse(5)` | `q = ark·q_n + brk·q + dt·crk·dq` |
| `rk_update_q_gmp` | `collapse(6)` | `q += dt·β(s)·q_rk(:,s)` for `s=1…nrk` |

---

## Immersed boundary method

### `adam_ib_gmp_object` — GPU IB wrapper

Holds pointers to the common `ib_object` and to `field_gmp_object` for grid metrics. Manages the signed-distance field `phi_gpu` on device.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `phi_gpu` | `(nb, 1-ngc:ni+ngc, …, n_solids+1)` | Signed-distance function; last slice = `max` over all solids |
| `q_bcs_vars_gpu` | `(:,:)` | Wall BC reference state per solid |

Sign convention: `phi < 0` inside solid (ghost region), `phi > 0` in fluid.

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(ib, field_gpu)` | Allocate `phi_gpu` and `q_bcs_vars_gpu`; associate grid metadata pointers via `associate_adam_data` |
| `evolve_eikonal(dq_gpu, q_gpu)` | For each solid: compute gradient-weighted residual `dq`, then apply `q -= dq` inside solid |
| `invert_eikonal(q_gpu)` | Enforce wall BC at surface (φ > 0): `BCS_VISCOUS` mirrors all velocity components; `BCS_EULER` reflects the normal component |

### `adam_ib_gmp_kernels` — IB OpenMP kernels

All spatial kernels use `collapse(4)` over `(b,k,j,i)`.

| Kernel | Purpose |
|--------|---------|
| `compute_phi_analytical_sphere_gmp` | `φ = −(‖x−xc‖ − R)` — negative inside sphere |
| `compute_phi_all_solids_gmp` | `φ_all = max(φ₁,…,φ_ns)` — union-of-solids mask |
| `compute_eikonal_dq_phi_gmp` | 1st-order upwind residual: `dq = |n_x|·Δq_x + |n_y|·Δq_y + |n_z|·Δq_z` inside solid; normal from centred `∇φ` |
| `evolve_eikonal_q_phi_gmp` | `q -= dq` where `φ > 0` |
| `invert_eikonal_q_phi_gmp` | Viscous: `(u,v,w)→(−u,−v,−w)`; Euler: `u → u − 2(u·n̂)n̂` |
| `move_phi_gmp` | Level-set advection `∂φ/∂t = −v·∇φ` for moving bodies (two-kernel: compute `dφ` then update) |
| `reduce_cell_order_phi_gmp` | Lower WENO order in cells within `ib_reduction_extent` of the surface, one kernel per spatial direction |

---

## Communication maps

### `adam_maps_gmp_object` — GPU maps wrapper

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

Index arrays use `integer(I8P)` to accommodate large AMR block counts. `mode=1` — one-to-one cell correspondence; `mode=8` — 8-cell average at AMR refinement interface.

All arrays are transferred via `assign_allocatable_gpu`. The send and receive buffers are exposed to MPI using `!$omp target data use_device_addr` rather than copied back to the host.

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(mpih, maps)` | Store pointers, call `copy_cpu_gpu(verbose=.true.)` |
| `copy_cpu_gpu(verbose)` | Transfer all map and buffer arrays via `assign_allocatable_gpu` |

---

## Module summary

| Module | Role | Extends / wraps |
|--------|------|----------------|
| `adam_gmp_library` | Aggregate entry point | — |
| `adam_gmp_utils` | `omp_target_alloc/free/memcpy` Fortran wrappers | — |
| `adam_memory_gmp_library` | High-level GPU alloc/copy with error checking | `adam_gmp_utils` |
| `adam_mpih_gmp_object` | OpenMP device binding, intra-node GPU assignment | `mpih_object` |
| `adam_field_gmp_object` | GPU field wrapper + host↔device transfer | `field_object` (pointer) |
| `adam_field_gmp_kernels` | Gradient, L2 norm, transpose, ghost-cell pack/unpack | — |
| `adam_weno_gmp_object` | GPU WENO coefficient mirror + ROR tables | `weno_object` (pointer) |
| `adam_rk_gmp_object` | GPU RK stage storage and update dispatch | `rk_object` (pointer) |
| `adam_rk_gmp_kernels` | Stage assign, accumulate, low-storage, final update | — |
| `adam_ib_gmp_object` | GPU distance field + eikonal BC wrapper | `ib_object` (pointer) |
| `adam_ib_gmp_kernels` | Eikonal evolution, sphere distance, momentum inversion | — |
| `adam_maps_gmp_object` | GPU communication index tables + MPI buffer staging | `maps_object` (pointer) |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
