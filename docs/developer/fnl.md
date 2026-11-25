# FNL library (GPU)

The FNL backend is ADAM's GPU layer, built on the [FUNDAL](https://github.com/szaghi/FUNDAL) device-memory library. The
same sources compile as OpenACC with nvfortran (`*-fnl-nvf` modes) or as OpenMP target offloading with amdflang
(`*-fnl-omp-amd` modes, `DEV_OMP`): every kernel carries an `!$acc` directive and its `!$omp` twin.

Every subsystem has two parts: a **wrapper object** holding the device arrays, and a **kernel module** with the device
routines. The wrappers are plain types, not extensions of the common classes (the Poisson solvers excepted). Each is a
component of the application's FNL type, initialised from the realm's common objects passed as arguments. Equations,
coefficients and data structures are defined once in `src/lib/common` ([library objects](./library)) and mirrored to
the device.

`use adam_fnl_library` exposes the FNL API together with `adam_common_library` and `adam_fnl_mpih_global` (the
exceptions are listed in the [module summary](#module-summary)).

## Program-scope singletons

The only FNL singleton is `mpih_fnl` (`adam_fnl_mpih_global`, type `mpih_fnl_object`): the GPU-aware MPI handler,
which also owns the device context and is initialised once per process: the forest calls the realm hook
`initialize_process_context` once, before any realm, and the FNL backends override it with `initialize_mpih_fnl()`
(`adam_fnl_mpih_global`: device binding and the ranks sharing the device) followed by the CPU `mpih%initialize`. The
other FNL objects are components of the application's FNL type, one set per realm, initialised after the common
initialisation from the realm's own objects (PRISM, `adam_prism_fnl_object.F90`):

```fortran
call self%prism_common_object%initialize(filename=filename, memory_avail=memory_avail_, verbose=.true.)
call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
call self%ib_fnl%initialize(grid=self%adam%grid, field=self%adam%field, ib=self%ib)
```

The time integrator is not among them: the realm builds it (`ts`, SELAH's stepper on the device fields, see
[Runge-Kutta integration](#runge-kutta-integration)).

Application FNL objects follow the same pattern (PRISM's `coil_fnl`, `pic_fnl`, `pml_fnl`). The WENO reconstruction
needs no helper: SETH's kernels have compile-time coefficients ([#139](https://github.com/szaghi/adam/issues/139)).

## Design and memory model

### CPU and GPU array layout

The common library stores fields with the variable index fastest; on the device the block index moves to the front,
so threads mapped to adjacent blocks read adjacent memory:

```
CPU:  q_cpu(nv, ni, nj, nk, nb)   — stride-1 on nv
GPU:  q_gpu(nb, ni, nj, nk, nv)   — blocks contiguous
```

The device state arrays (`q_gpu`, `dq_gpu`, …) belong to the application's FNL type (e.g. `flume_fnl_object%q_gpu`);
`field_fnl%initialize(..., q_gpu=)` can allocate one in this layout. The application transposes between layouts with
FUNDAL's buffered transposed copies (`dev_memcpy_to_device(..., ij=[1,5], buf=)`) or with the device kernels
`pack_blocks_transposed_dev` / `unpack_blocks_transposed_dev` of `adam_fnl_field_kernels`. The staging buffer of a
buffered copy must have the destination's extent; otherwise nvfortran makes a host copy-in temporary and the copy
fails (use `dev_assign_to_device(..., ij=[1,5])` instead). The transposition happens only at I/O and AMR updates;
time integration uses the device layout alone.

### FUNDAL integration

| Utility | Purpose |
|---------|---------|
| `dev_alloc` | Allocate a device array |
| `dev_assign_to_device` | Copy a host array to the device |
| `dev_memcpy_to_device` | Raw copy to the device |
| `DEVICEVAR(array)` macro | Mark an array as device-resident (no implicit copy) |

`DEVICEVAR` is defined in `fundal.H` and used at the top of every kernel file: it tells the compiler that the array
lives on the device, so the OpenACC runtime inserts no host copy.

### IB solid masking

Kernels that modify field variables skip the cells inside a solid (`phi_gpu(b,i,j,k,all_solids) < 0`), tested in the
innermost loop.

## Field

### `adam_fnl_field_object`

A plain type, initialised from the realm's `grid`, `field` and `maps`. It holds the device geometry and maps (its
`maps` component is a `maps_fnl_object`) and the device ghost exchange.

| Device array | Shape | Purpose |
|--------------|-------|---------|
| `x_cell_gpu`, `y_cell_gpu`, `z_cell_gpu` | `(nb, 1-ngc:ni+ngc)` per axis | Cell centroid coordinates |
| `dxyz_gpu` | `(nb, 3)` | Block spacing |
| `fec_1_6_array_gpu` | `(26)` | Face/edge/corner codes 1–26 → boundary faces 1–6 |

| Method | Purpose |
|--------|---------|
| `initialize(grid, field, maps[, nv_aux][, q_gpu][, verbose])` | Allocate the device arrays (and `q_gpu(nb,…,nv)` when passed) |
| `copy_cpu_gpu` | Copy the coordinates and maps to the device |
| `reduce_ghost_local_gpu`, `reduce_ghost_mpi_gpu` | Fold deposited ghosts back onto their owner cells, locally and across ranks |
| `update_ghost_local_gpu` | Intra-rank ghost update on the device |
| `update_ghost_mpi_gpu` | Pack on the device, MPI exchange, unpack on the device |
| `compute_q_gradient(b, ivar, dx, dy, dz, q_gpu, gradient)` | AMR criterion `max|∇q|/|q|` over block `b` |

`initialize`, `copy_cpu_gpu`, `update_ghost_mpi_gpu` and `compute_q_gradient` take the realm's CPU
`grid`/`field`/`maps` (and the per-rank `comm_map_*_ptr_ghost`) as arguments, so a multi-realm caller always
exchanges with its own decomposition.

### GPU-direct ghost exchange

`update_ghost_mpi_gpu` posts `MPI_Isend`/`MPI_Irecv` directly on the device ghost buffers and waits with
`MPI_Waitall`, with no global barrier before or after (removed 2026-11-25: after the wait the rank's buffers are done,
and MPI's non-overtaking rule matches every message to the receive of its own exchange; a rank with no message posts
nothing, so the exchange is not collective). It posts on the device ghost buffers (`send_buffer_ghost_gpu`,
`recv_buffer_ghost_gpu`), which lets CUDA-aware MPI use GPUDirect RDMA. There is no host-staging fallback. Measured with `ucx_perftest` (33 MB,
[performance](/applications/flume/performance#another-machine-h100-nodes-with-gpudirect-rdma)): GPU to GPU between nodes
52.4 GB/s on ESA SpaceHPC (H100, GPUDirect RDMA through dma-buf), 1.5 GB/s on medooza (A30, no GPUDirect RDMA: UCX
stages the buffers through host memory). How fast the local fill moves its data against the GPU's copy rate, and what a
structured (slab) fill would gain, is measured in the [ghost-fill benchmark](/developer/ghost-fill-benchmark). On the WSL2
development box this path aborts with two or more ranks unless UCX is forced eager (`UCX_RNDV_THRESH=inf`); that
setting is WSL-only and a regression on a cluster. Details: [GPU5](/cookbook/gpu#gpu5-run-on-a-development-box-with-wsl2-not-on-a-cluster).

### `adam_fnl_field_kernels`

All routines carry `!$acc parallel loop independent` and `DEVICEVAR` on every device argument.

| Kernel | Purpose |
|--------|---------|
| `compute_q_gradient_dev` | Centred gradient magnitude, `reduction(max:)` |
| `compute_normL2_residuals_dev` | L2 norm `√(Σ dq²)` per variable, `reduction(+:)` |
| `copy_transpose_gpu_cpu_dev` | Transpose `(nb,ni,nj,nk,nv)` → `(nv,ni,nj,nk,nb)` on the device |
| `pack_blocks_transposed_dev`, `unpack_blocks_transposed_dev` | Layout transposition of the state |
| `populate_send_buffer_ghost_gpu_dev` | Pack ghosts into the MPI send buffer: copies, fine → coarse restrictions, coarse → fine seam interpolations; with the density guard on, a first kernel computes the guard's part once per fine ghost, at the cell's density row ([#157](https://github.com/szaghi/adam/issues/157) B3) |
| `receive_recv_buffer_ghost_gpu_dev` | Unpack the MPI receive buffer into ghosts |
| `update_ghost_local_gpu_dev` | Intra-rank block-to-block ghost update: copies, fine → coarse restrictions (mean, 6-point Lagrange, WENO) and coarse → fine seam interpolations; with the density guard on, a first kernel computes the guard's part once per fine ghost, then the guarded fill; with the guard off, a fill without the guard's code ([#157](https://github.com/szaghi/adam/issues/157)) |
| `seam_guard_density_gpu_dev`, `seam_guard_value_gpu_dev` | The density guard of the seam fill ([#101](https://github.com/szaghi/adam/issues/101)) in two parts, device twins of the CPU's `seam_guard_prepare` and `seam_guard_value`: the interpolated density of a fine ghost clamped to its bound (once per ghost), then its application to one variable |

## WENO reconstruction

The FNL library holds no WENO code: the device kernels of FLUME and RIPPLE call [SETH](https://github.com/szaghi/SETH)'s
WENO kernels, `acc routine seq` and `omp declare target` with compile-time coefficients, so no table is mirrored on the
device and no table argument crosses a kernel boundary ([#139](https://github.com/szaghi/adam/issues/139)). FLUME's face
kernels are one per stencil width (`..._s<N>_dev`, S a named constant), each calling SETH's kernel of that width
(`seth_weno<2S-1>_js`, `_interp_wr`); RIPPLE calls the library's run-time dispatcher `weno_reconstruct`. The former
`adam_fnl_weno_object` (the coefficient mirror and the unused ROR tables), `adam_fnl_weno_kernels` and
`adam_fnl_weno_width_kernels` were deleted in S6; on the device SETH's kernels were measured bitwise equal to them
(nvfortran 26.9) before the swap.

## Runge-Kutta integration

Every FNL backend (FLUME, PRISM, RIPPLE) integrates with [SELAH](https://github.com/szaghi/SELAH)'s device build through
the library adapter `adam_stepper_object` ([#132](https://github.com/szaghi/adam/issues/132), phases 3-5; ADAM's own
device integrator, `adam_fnl_rk_object` and its kernels, was deleted in phase 6): the application sets
`selah_layout = LAYOUT_V_LAST` before the realm initializes, and the realm builds `ts` on the `(b, i, j, k, v)` device
fields (`q_gpu`, `dq_gpu`, `phi_gpu`). The registers are device memory, zeroed at
allocation; the stage state, the residual register and the scratch registers are device pointers the FNL kernels take. The residual of a stage goes into `dq_gpu` (the per-stage reflux corrects it there) and
`ts%set_rhs` copies it into the residual register with SELAH's copy kernel. SELAH's kernels sum the stages in the order
of the deleted ADAM kernels, so FLUME FNL is bitwise the integration it replaced (every regression case and the #132 overlays,
A/B on an A30), except the step residual in the solid cells, now 0 ([#128](https://github.com/szaghi/adam/issues/128)
D4). The end_of_step seam ghosts reach each stage state through `seam_ghosts_to_state_dev` (#132 rule 1).

## Immersed boundary method

### `adam_fnl_ib_object`

A plain type, initialised from the realm's `grid`, `field` and `ib`. It holds the signed distance `phi_gpu` and runs
the eikonal wall enforcement. Sign convention: `phi < 0` inside the solid, `phi > 0` in the fluid.

| Device array | Shape | Purpose |
|--------------|-------|---------|
| `phi_gpu` | `(nb, ni, nj, nk, n_solids+1)` | Signed distance; the last slice is the union of all solids |
| `q_bcs_vars_gpu` | `(:,:)` | Wall state variables per solid |

| Method | Purpose |
|--------|---------|
| `initialize(grid, field, ib)` | Allocate `phi_gpu` and `q_bcs_vars_gpu`; copy the wall data to the device |
| `evolve_eikonal(dxyz_gpu, dq_gpu, q_gpu)` | Advance the eikonal extrapolation inside the solid; `dxyz_gpu` passed by the caller (`field_fnl%dxyz_gpu`) |
| `invert_eikonal(q_gpu)` | Wall condition at the surface: `BCS_VISCOUS` reverses the velocity, `BCS_EULER` reflects its normal component `u → u − 2(u·n̂)n̂` |

### `adam_fnl_ib_kernels`

| Kernel | Purpose |
|--------|---------|
| `compute_phi_analytical_sphere_dev` | `φ = −(‖x − xc‖ − R)` |
| `compute_phi_all_solids_dev` | `φ_all = max(φ₁, …, φ_ns)` |
| `compute_eikonal_dq_phi_dev` | `dq = |n_x|·Δq_x + |n_y|·Δq_y + |n_z|·Δq_z` inside the solid |
| `evolve_eikonal_q_phi_dev` | `q -= dq` inside the solid |
| `invert_eikonal_q_phi_dev` | Momentum reflection at the surface |
| `move_phi_dev`, `reduce_cell_order_phi_dev` | Level-set advection and order reduction near solids; used by no maintained application |

## Communication maps

### `adam_fnl_maps_object`

A plain type, the `maps` component of `field_fnl_object`, initialised from the realm's `maps`. It mirrors the
communication tables to the device, so ghost packing and unpacking need no CPU staging.

| Device array | Columns | Content |
|--------------|---------|---------|
| `local_map_ghost_cell_gpu` | 10 | `(b_src, b_dst, i_src, j_src, k_src, i_dst, j_dst, k_dst, mode, meta)` |
| `local_map_ghost_cell_prol_gpu` | — | The `mode = 4` rows of `local_map_ghost_cell_gpu` (null if none) |
| `seam_rho_gpu` | — | Scratch of the local fill, one value per row: the clamped density of each `mode = 4` row (the guard's part, once per fine ghost) |
| `comm_map_send_ghost_cell_gpu` | 8 | `(b_src, i, j, k, v, buf_idx, mode, meta)`, one row per cell and variable, the variables of a cell contiguous and in order |
| `comm_map_send_density_rows_gpu` | — | The density rows of the `mode = 4` cells of the send map (null if none or the guard is off) |
| `seam_rho_send_gpu` | — | Scratch of the pack, one value per row: the clamped density at those rows |
| `comm_map_recv_ghost_cell_gpu` | 6 | `(buf_idx, b_dst, i, j, k, v_offset)` |
| `send_buffer_ghost_gpu`, `recv_buffer_ghost_gpu` | — | Packed MPI buffers |
| `local_map_bc_crown_gpu` | — | Boundary-condition crown map |

`mode = 1`: one-to-one cells (same level); `mode = 8`: fine → coarse restriction at an AMR interface (the mean of the
2×2×2 fine cells, or the 6-point Lagrange or WENO restriction coded in `meta`); `mode = 4`: coarse → fine interpolation
of the seam ghost (regime of `[amr] seam_ghost_fill`, anchor and octant packed in `meta`). The rows keep the order of
the map build: grouping them by kind was measured and is slower on the H100 with the guard once per ghost (the
prolongations become the tail of the launch; [ghost-fill benchmark](/developer/ghost-fill-benchmark#applied-the-density-guard-once-per-ghost)).

| Method | Purpose |
|--------|---------|
| `initialize(maps)` | Initialise and call `copy_cpu_gpu` |
| `copy_cpu_gpu(verbose)` | Copy every map to the device (`dev_assign_to_device`) |

## Finite difference/volume operators

### `adam_fnl_fdv_operators_library`

The operators of `adam_fdv_operators_library`, callable inside parallel regions (`!$acc routine seq`, one thread per
cell). Each takes the half-width `s` and the spacing `dxyz`, so the order is chosen at call time.

| Operator | FD centred | FV centred |
|----------|-----------|-----------|
| Gradient | `compute_gradient_fd_centered_dev` | `compute_gradient_fv_centered_dev` |
| Divergence | `compute_divergence_fd_centered_dev` | `compute_divergence_fv_centered_dev` |
| Curl | `compute_curl_fd_centered_dev` | `compute_curl_fv_centered_dev` |
| Laplacian | `compute_laplacian_fd_centered_dev` | `compute_laplacian_fv_centered_dev` |

## Turbulence

### `adam_fnl_turbulence_kernels`

Device twins of `adam_turbulence_kernels` ([library objects](./library#adam-turbulence-object-adam-turbulence-library-adam-turbulence-kernels)):
`compute_mu_t_vreman_dev`, `compute_mu_t_wale_dev`, `compute_mu_t_sigma_dev`, on transposed arrays `(b, i, j, k, v)`.
They call the pointwise models of `adam_turbulence_library`; the body is `adam_fnl_turbulence_body.INC`.

## Poisson solvers

The device twins of the [FLAIL](./flail) solvers. Both extend their common parent.

### `adam_fnl_poisson_mg_object`

Extends `poisson_mg_object`. The leaves, the hierarchy and every CSR operator are assembled on the host; `work_place_fnl`
then puts the operators, the level diagonals, colours and leaf maps, the cell volumes and the work vectors on the
device until `destroy`, with halo buffers sized for the largest operator.

| Part | Implementation |
|------|----------------|
| `apply_fnl` | Packs the halo on the device, exchanges it by `MPI_ALLTOALLV` on device buffers (GPU-aware MPI; none on one rank), one row per thread in the host's order |
| `smooth_fnl` | One kernel per colour after its halo exchange. It equals the host's sequential Gauss-Seidel: the colour is the global cell parity and the level operators are 7-point, so no row couples to its own colour (no such entry in 5704 smoother calls over the AMR, isolated, boundary-kind and sub-block suites) |
| Vector hooks | Device kernels; dots and max-norm as device reductions |
| Root | Through the host to the FFT solver, cuFFT by default (`initialize_fnl` passes `poisson_fft_fnl_object` when no `coarse` is given), at most 32^3 cells (`root_cells_max_fnl`) |
| On the host | Right-hand side, James step, ghost-fill placement |

`initialize_fnl` releases and places the device copies itself, since the parent's `initialize` dispatches to the
host's no-op hooks. Results equal the host's to round-off (summation order of the dots, FMA contraction). FLUME's FNL
backend and `bench_poisson_fnl` set `flail%mg_mold` to it; `test-poisson-mg-fnl-nvf` builds `test_poisson_mg` on it.

| Root size on cuFFT, 3.9 M AMR cells | Cost per solve |
|-------------------------------------|----------------|
| 32^3 (default) | 1× |
| 128^3 | 2.4× to 6.8× |

### `adam_fnl_poisson_fft_object`

Extends `poisson_fft_object` with cuFFT. In the distributed slab transform the block planes are redistributed on the
host; the transform runs on each rank's GPU with its arrays resident (`slab_r`, `slab_c`, `pen`, `khat_pen`, the
transpose buffers, created with the plans):

| Hook | Device work |
|------|-------------|
| `slab_update_fnl` | `slab_r` up and down |
| `slab_planes_fnl` | One D2Z/Z2D plan for one plane, executed on every plane |
| `slab_transpose_fnl` | Pack and unpack on the device, `MPI_ALLTOALLV` on device buffers (none on one rank) |
| `slab_lines_fnl` | One Z2Z plan for one row of z-lines, executed on every row |
| `slab_spectral_fnl` | Division by the eigenvalues, or product with the kernel spectrum |

The result does not depend on the rank count. The rank-0 path remains for PRISM's CG solves and the multigrid's small
root: the global source goes to the device, D2Z, a kernel divides by the eigenvalues (host-computed tables), Z2D, and
the potential comes back. cuFFT is bound to the stream of the synchronous OpenACC queue
(`acc_get_cuda_stream(acc_async_sync)`), so transforms and kernels are ordered. Isolated: the padded transforms, with
the kernel spectrum `khat` made once on the device (`enter data`, released by `destroy_plans`).

| Item | Detail |
|------|--------|
| Build | `-D_CUFFT -cudalib=cufft` (mode `test-poisson-fft-fnl-nvf`); `cufft` and `openacc` are intrinsic modules in `fobos.d/templates.fobos`. NVIDIA only: without `_CUFFT` (AMD) making the plans stops the run |
| Memory | Not FUNDAL-agnostic: the arrays live in an OpenACC data region around the cuFFT calls (`host_data use_device`) |
| Verification | SV-2 against the discrete solution at round-off, SV-6 against the exact potentials of two spheres; both against the FFTW solver (`selfgravity/poisson_oracle.py compare`, isolated: `1.6e-15`) |

## MPI handler

### `adam_fnl_mpih_object`

```fortran
use :: fundal_mpih_object, only : mpih_fnl_object => mpih_object
```

FUNDAL's MPI handler under the alias `mpih_fnl_object`: rank and size, rank-prefixed output, timing, GPU-aware MPI.

## Module summary

| Module | Role | Extends |
|--------|------|---------|
| `adam_fnl_library` | Aggregate entry point | — |
| `adam_fnl_field_object` | Device geometry, maps, ghost exchange | — |
| `adam_fnl_field_kernels` | Gradient, L2 norm, ghost pack/unpack, transposition | — |
| `adam_fnl_ib_object` | Signed distance and eikonal wall condition | — |
| `adam_fnl_ib_kernels` | Eikonal evolution, sphere distance, momentum inversion | — |
| `adam_fnl_maps_object` | Communication tables and MPI buffers | — |
| `adam_fnl_turbulence_kernels` | Subgrid-scale eddy viscosity | — |
| `adam_fnl_fdv_operators_library` | Device-callable FD/FV operators | — |
| `adam_fnl_poisson_fft_object` | FFT Poisson solver on cuFFT (not re-exported) | `poisson_fft_object` |
| `adam_fnl_poisson_mg_object` | Multigrid Poisson solver on the device (not re-exported) | `poisson_mg_object` |
| `adam_fnl_mpih_object` | FUNDAL MPI handler alias | — |
| `adam_fnl_mpih_global` | Singleton `mpih_fnl` | — |
