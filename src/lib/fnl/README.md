# ADAM FNL Library — `src/lib/fnl`

The FNL backend is ADAM's **OpenACC GPU acceleration layer**, built on top of the [FUNDAL](https://github.com/szaghi/FUNDAL) GPU memory management library. It follows a consistent two-tier pattern for every subsystem: a **wrapper object** that extends the corresponding common CPU class and manages device-resident arrays, paired with a **kernel module** that contains OpenACC-decorated device subroutines implementing the actual computation.

No physics or algorithmic logic is duplicated — all equations, coefficients, and data structures are defined once in `src/lib/common` and mirrored to the GPU by the FNL layer.

The aggregate entry point `adam_fnl_library` re-exports the entire FNL API together with `adam_common_library`, including all FNL singleton modules; a single `use adam_fnl_library` statement in application code exposes both layers.

---

## Contents

- [Design and memory model](#design-and-memory-model)
- [Program-scope singletons](#program-scope-singletons)
- [Field](#field)
- [WENO reconstruction](#weno-reconstruction)
- [Runge-Kutta integration](#runge-kutta-integration)
- [Immersed boundary method](#immersed-boundary-method)
- [Communication maps](#communication-maps)
- [Finite difference/volume operators](#finite-differencevolume-operators)
- [MPI handler](#mpi-handler)
- [Module summary](#module-summary)

---

## Program-scope singletons

Each FNL object is exposed as a **module-level target variable** in its own singleton module. GPU backend code accesses these singletons via `use`; they are never passed as dummy arguments and never embedded in derived types.

| Module | Variable | Type | Purpose |
|--------|----------|------|---------|
| `adam_fnl_mpih_global` | `mpih_fnl` | `mpih_fnl_object` | GPU-aware MPI handler and device initialization |
| `adam_fnl_field_global` | `field_fnl` | `field_fnl_object` | Device field arrays, coordinates, communication maps |
| `adam_fnl_ib_global` | `ib_fnl` | `ib_fnl_object` | Device distance function `phi_gpu` and wall BC variables |
| `adam_fnl_rk_global` | `rk_fnl` | `rk_fnl_object` | Device RK stage arrays |
| `adam_fnl_weno_global` | `weno_fnl` | `weno_fnl_object` | Device WENO coefficients and ROR tables |

All five are re-exported by `adam_fnl_library`.

**Initialization order** — CPU value singletons must be populated from the solver's owned state before calling FNL `%initialize()`:

```fortran
ib   = self%ib    ! cpu ib_object  → ib  singleton
rk   = self%rk    ! cpu rk_object  → rk  singleton
weno = self%weno  ! cpu weno_object → weno singleton
call mpih_fnl%initialize(do_mpi_init=.true., do_device_init=.true., verbose=.true.)
call field_fnl%initialize(...)
call ib_fnl%initialize()
call rk_fnl%initialize()
call weno_fnl%initialize()
```

Application-level backends may define additional FNL singletons for app-specific GPU objects (e.g. `coil_fnl`, `fwlayer_fnl` in the PRISM application).

---

## Design and memory model

### CPU vs GPU array layout

The common library stores field data in Fortran column-major order with variables as the fastest-varying index:

```
CPU:  q_cpu(nv, ni, nj, nk, nb)   — stride-1 on nv
```

On the GPU the block index is moved to the front so that threads mapped to adjacent blocks access adjacent memory locations, improving warp-level coalescing:

```
GPU:  q_gpu(nb, ni, nj, nk, nv)   — blocks contiguous
```

The transposition between layouts is performed by `copy_transpose_cpu_gpu` / `copy_transpose_gpu_cpu` in `adam_fnl_field_object` and their device kernels in `adam_fnl_field_kernels`. This copy happens only at I/O and AMR update boundaries; during time integration the GPU layout is used exclusively.

### FUNDAL integration

GPU memory allocation and host-device transfers route through FUNDAL utilities:

| Utility | Purpose |
|---------|---------|
| `dev_alloc` | Allocate device array |
| `dev_assign_to_device` | Copy host array to device |
| `dev_memcpy_to_device` | Raw device-to-device copy |
| `DEVICEVAR(array)` macro | Mark array as already device-resident (suppresses implicit copy) |

The `DEVICEVAR` macro is defined in `fundal.H` and appears at the top of every kernel file. It is the primary mechanism by which the compiler is told that a given array lives on the device, preventing spurious host-copy insertion by the OpenACC runtime.

### IB solid masking

Every kernel that modifies field variables carries an implicit guard: cells inside a solid body (`phi_gpu(b,i,j,k,all_solids) < 0`) are skipped. This masking is applied at the innermost loop level so that it incurs no branch divergence in fluid-only regions.

---

## Field

### `adam_fnl_field_object` — GPU field wrapper

Extends `field_object` (common). Holds all device-side arrays needed for field operations and provides the host-device transfer interface.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_gpu` | `(nb, ni, nj, nk, nv)` | Primary field — conservative variables |
| `q_t_gpu` | `(nv, ni, nj, nk, nb)` | Transposed scratch — used during CPU↔GPU copies |
| `x_cell_gpu`, `y_cell_gpu`, `z_cell_gpu` | `(nb, ni, nj, nk)` | Cell centroid coordinates |
| `dxyz_gpu` | `(nb, 3)` | Block mesh spacing `(dx, dy, dz)` |
| `fec_1_6_array_gpu` | `(nb, 26)` | Face enumeration codes for IB ghost-cell lookup |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize` | Allocate all device arrays via `dev_alloc` |
| `copy_cpu_gpu` | Transfer `q`, coordinate arrays, and maps to device |
| `copy_transpose_cpu_gpu(nv, q_cpu, q_gpu)` | Transpose and copy `q_cpu(nv,…,nb)` → `q_gpu(nb,…,nv)` |
| `copy_transpose_gpu_cpu(nv, q_gpu, q_cpu)` | Inverse transpose: `q_gpu(nb,…,nv)` → `q_cpu(nv,…,nb)` |
| `update_ghost_local_gpu` | Apply intra-rank ghost-cell updates entirely on device |
| `update_ghost_mpi_gpu` | Pack send buffer on device, perform MPI exchange, unpack on device |
| `compute_q_gradient(b, ivar, q_gpu, gradient)` | AMR refinement criterion: `max|∇q|/|q|` over block `b` |

### `adam_fnl_field_kernels` — field device kernels

All routines carry `!$acc parallel loop independent` and use `DEVICEVAR` on every device pointer argument.

| Kernel | Purpose |
|--------|---------|
| `compute_q_gradient_dev` | Centred-difference gradient magnitude with `reduction(max:)` |
| `compute_normL2_residuals_dev` | L2 norm `√(Σ dq²)` per variable with `reduction(+:)` |
| `copy_transpose_gpu_cpu_dev` | Transpose `(nb,ni,nj,nk,nv)` → `(nv,ni,nj,nk,nb)` on device |
| `populate_send_buffer_ghost_gpu_dev` | Pack ghost-cell values into MPI send buffer; supports 1-cell and 8-cell AMR averaging |
| `receive_recv_buffer_ghost_gpu_dev` | Unpack MPI receive buffer into ghost cells |
| `update_ghost_local_gpu_dev` | Apply intra-rank block-to-block ghost updates; supports AMR coarse↔fine averaging |

---

## WENO reconstruction

### `adam_fnl_weno_object` — GPU WENO coefficient wrapper

Extends `weno_object` (common). The CPU object computes all WENO coefficients once during initialisation; the FNL object mirrors them to device memory and holds the ROR (Reduced-Order Reconstruction) tables used near solid boundaries.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `a_gpu` | `(2, 0:S-1, S)` | Optimal WENO weights per sub-stencil and face |
| `p_gpu` | `(2, 0:S-1, 0:S-1, S)` | Polynomial reconstruction coefficients |
| `d_gpu` | `(0:S-1, 0:S-1, 0:S-1, S)` | Smoothness indicator coefficients |
| `ror_schemes_gpu` | `(:)` | ROR fallback scheme orders near solid walls |
| `ror_ivar_gpu` | `(:)` | Variable indices checked by ROR |
| `cell_scheme_gpu` | `(nb, ni, nj, nk)` | Per-cell effective reconstruction order |

### `adam_fnl_weno_kernels` — WENO device kernels

| Kernel | Directive | Purpose |
|--------|-----------|---------|
| `weno_reconstruct_upwind_dev(S, a, p, d, zeps, V, VR)` | `!$acc routine seq` | Reconstruct left (`VR(1)`) and right (`VR(2)`) interface values from stencil `V` |

The reconstruction follows the standard three-step algorithm:
1. Compute `S` polynomial reconstructions from overlapping sub-stencils
2. Compute smoothness indicators from second-derivative sums
3. Weight and convolve: `VR(f) = Σ_k w(f,k) · VP(f,k)`

`!$acc routine seq` marks the procedure as callable from within a parallel region without launching a new kernel — one thread per stencil, invoked inside the outer loop over cells.

---

## Runge-Kutta integration

### `adam_fnl_rk_object` — GPU RK stage manager

Extends `rk_object` (common). Manages stage storage on device and drives the per-stage updates.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `q_rk_gpu` | `(nb, ni, nj, nk, nv, nrk)` | Stage values (1 stage for low-storage, `nrk` for SSP) |
| `alph_gpu` | `(nrk, nrk)` | SSP alpha coefficients |
| `beta_gpu` | `(nrk)` | SSP beta coefficients |
| `gamm_gpu` | `(nrk)` | SSP gamma coefficients |

**Supported schemes:**

| Scheme | Storage mode | Stages |
|--------|-------------|--------|
| `RK_1`, `RK_2`, `RK_3` | Low-storage | 1 — overwrites `q_gpu` in place |
| `RK_SSP_22`, `RK_SSP_33` | Multi-stage | 2 / 3 |
| `RK_SSP_54` | Multi-stage | 5 |

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(rk, nb, ngc, ni, nj, nk, nv)` | Allocate `q_rk_gpu` sized to scheme requirements |
| `initialize_stages(q_gpu)` | Broadcast `q_gpu` into all stage slots |
| `assign_stage(s, q_gpu, phi_gpu)` | Copy `q_gpu` into stage `s`, skipping solid cells |
| `compute_stage(s, dt, phi_gpu)` | Accumulate stages 1…s−1 into stage `s` (SSP) |
| `compute_stage_ls(s, dt, phi_gpu, dq_gpu, q_gpu)` | Low-storage update: `q = ark·q_n + brk·q + dt·crk·dq` |
| `update_q(s, dt, phi_gpu)` | Final assembly: `q += dt·beta(s)·q_rk(:,:,:,:,:,s)` |

### `adam_fnl_rk_kernels` — RK device kernels

All kernels carry `!$acc parallel loop independent` and mask solid cells via `phi_gpu`.

| Kernel | Purpose |
|--------|---------|
| `rk_assign_stage_dev` | `q_rk(:,s) ← q_gpu` (fluid cells only) |
| `rk_initialize_stages_dev` | `q_rk(:,all_s) ← q_gpu` |
| `rk_compute_stage_dev` | `q_rk(:,s) += dt·α(s,ss)·q_rk(:,ss)` for `ss = 1…s−1` |
| `rk_compute_stage_ls_dev` | `q = ark·q_n + brk·q + dt·crk·dq` (low-storage) |
| `rk_update_q_dev` | `q += dt·β(s)·q_rk(:,s)` for `s = 1…nrk` |

---

## Immersed boundary method

### `adam_fnl_ib_object` — GPU IB wrapper

Extends `ib_object` (common). Manages the signed-distance field `phi_gpu` on device and drives the eikonal enforcement cycle.

**Device arrays:**

| Array | Shape | Purpose |
|-------|-------|---------|
| `phi_gpu` | `(nb, ni, nj, nk, n_solids+1)` | Signed-distance function; last slice holds `max` over all solids |
| `q_bcs_vars_gpu` | `(:,:)` | Boundary condition state variables per solid |

The sign convention: `phi < 0` inside the solid (ghost region), `phi > 0` in the fluid.

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize` | Allocate `phi_gpu` and `q_bcs_vars_gpu` (reads dimensions from `field_fnl` singleton); copy BCS data to device |
| `evolve_eikonal(dxyz_gpu, dq_gpu, q_gpu)` | Advance eikonal equation inside solid: `q -= ∇φ·(q_bc − q)`; `dxyz_gpu` is passed explicitly by the caller (typically `field_fnl%dxyz_gpu`) |
| `invert_eikonal(q_gpu)` | Enforce wall BC at solid surface (φ > 0): reflect momentum |

Wall BC modes applied by `invert_eikonal`:
- **`BCS_VISCOUS`** (no-slip): `(u, v, w) → (−u, −v, −w)`
- **`BCS_EULER`** (inviscid): `u → u − 2(u·n̂)n̂`

### `adam_fnl_ib_kernels` — IB device kernels

| Kernel | Purpose |
|--------|---------|
| `compute_phi_analytical_sphere_dev` | `φ = −(‖x − xc‖ − R)` — negative inside sphere |
| `compute_phi_all_solids_dev` | `φ_all = max(φ₁, φ₂, …, φ_ns)` — union of all solids |
| `compute_eikonal_dq_phi_dev` | Gradient-weighted residual: `dq = |n_x|·Δq_x + |n_y|·Δq_y + |n_z|·Δq_z` inside solid |
| `evolve_eikonal_q_phi_dev` | `q -= dq` inside solid (φ > 0) |
| `invert_eikonal_q_phi_dev` | Momentum reflection at surface (BCS_VISCOUS or BCS_EULER) |
| `move_phi_dev` | Level-set advection: `∂φ/∂t = −v·∇φ` for moving bodies |
| `reduce_cell_order_phi_dev` | Lower reconstruction order in cells adjacent to solid surface |

---

## Communication maps

### `adam_fnl_maps_object` — GPU maps wrapper

Extends `maps_object` (common). Mirrors all communication index tables to device memory so that ghost-cell packing and unpacking happen entirely on the GPU, eliminating CPU staging for MPI buffers.

**Device arrays:**

| Array | Columns | Content |
|-------|---------|---------|
| `local_map_ghost_cell_gpu` | 9 | `(b_src, b_dst, i_src, j_src, k_src, i_dst, j_dst, k_dst, mode)` |
| `comm_map_send_ghost_cell_gpu` | 7 | `(b_src, i, j, k, v_offset, buf_idx, mode)` |
| `comm_map_recv_ghost_cell_gpu` | 6 | `(buf_idx, b_dst, i, j, k, v_offset)` |
| `send_buffer_ghost_gpu` | — | 1D packed MPI send staging buffer |
| `recv_buffer_ghost_gpu` | — | 1D packed MPI receive staging buffer |
| `local_map_bc_crown_gpu` | — | Boundary condition crown ghost-cell map |

The `mode` column distinguishes two cases:
- `mode = 1` — one-to-one cell correspondence (same refinement level)
- `mode = 8` — eight-cell average (fine block → coarse block at AMR interface)

**Key methods:**

| Method | Purpose |
|--------|---------|
| `initialize(maps)` | Initialise and call `copy_cpu_gpu` |
| `copy_cpu_gpu(verbose)` | Transfer all map arrays from CPU to device via `dev_assign_to_device` |

---

## Finite difference/volume operators

### `adam_fnl_fdv_operators_library` — device-callable spatial operators

Provides the same operators as `adam_fdv_operators_library` (common) in a form callable from within OpenACC parallel regions. All routines carry `!$acc routine seq` — no internal parallelism, one thread per cell.

**Available operators:**

| Operator | FD centred | FV centred |
|----------|-----------|-----------|
| Gradient `∇q` | `compute_gradient_fd_centered_dev` | `compute_gradient_fv_centered_dev` |
| Divergence `∇·q` | `compute_divergence_fd_centered_dev` | `compute_divergence_fv_centered_dev` |
| Curl `∇×q` | `compute_curl_fd_centered_dev` | `compute_curl_fv_centered_dev` |
| Laplacian `∇²q` | `compute_laplacian_fd_centered_dev` | `compute_laplacian_fv_centered_dev` |

Each routine accepts a stencil half-width `s` and the local mesh spacing `dxyz`, allowing the accuracy order to be selected at call time without recompilation.

---

## MPI handler

### `adam_fnl_mpih_object` — FUNDAL MPI alias

```fortran
use :: fundal_mpih_object, only : mpih_fnl_object => mpih_object
```

A direct re-export of FUNDAL's MPI handler under the FNL-namespaced type alias `mpih_fnl_object`. Provides rank/size queries, rank-prefixed console output, and timing utilities. No FNL-specific extensions are needed because FUNDAL's handler already covers GPU-aware MPI requirements.

---

## Module summary

| Module | Role | Extends |
|--------|------|---------|
| `adam_fnl_library` | Aggregate entry point | — |
| `adam_fnl_field_object` | GPU field wrapper + host↔device transfer | `field_object` |
| `adam_fnl_field_kernels` | Gradient, L2 norm, ghost-cell pack/unpack | — |
| `adam_fnl_weno_object` | GPU WENO coefficient mirror + ROR tables | `weno_object` |
| `adam_fnl_weno_kernels` | Upwind WENO reconstruction (`!$acc routine seq`) | — |
| `adam_fnl_rk_object` | GPU RK stage storage and update dispatch | `rk_object` |
| `adam_fnl_rk_kernels` | Stage assign, accumulate, low-storage, final update | — |
| `adam_fnl_ib_object` | GPU distance field + eikonal BC wrapper | `ib_object` |
| `adam_fnl_ib_kernels` | Eikonal evolution, sphere distance, momentum inversion | — |
| `adam_fnl_maps_object` | GPU communication index tables + MPI buffer staging | `maps_object` |
| `adam_fnl_mpih_object` | FUNDAL MPI handler alias | — |
| `adam_fnl_fdv_operators_library` | Device-callable FD/FV spatial operators | — |
| `adam_fnl_mpih_global` | Singleton — `mpih_fnl` | — |
| `adam_fnl_field_global` | Singleton — `field_fnl` | — |
| `adam_fnl_ib_global` | Singleton — `ib_fnl` | — |
| `adam_fnl_rk_global` | Singleton — `rk_fnl` | — |
| `adam_fnl_weno_global` | Singleton — `weno_fnl` | — |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
