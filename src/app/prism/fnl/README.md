# PRISM — FNL (OpenACC) Backend

> PRISM is part of the ADAM framework, released under the
> [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.en.html) (LGPLv3).
> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.

The FNL backend (`prism-fnl`) is the GPU-accelerated implementation of PRISM using the
**FNL library** (Fortran aNd LAPACK, OpenACC-based GPU framework). It inherits the full PRISM
common infrastructure and extends it with GPU-resident arrays, OpenACC-parallelized FDV operators,
and GPU-aware MPI communication via the FUNDAL library.

## Source files

| File | Module | Description |
|------|--------|-------------|
| `adam_prism_fnl.F90` | entry point | Parses CLI argument, delegates to `prism%simulate` |
| `adam_prism_fnl_object.F90` | `adam_prism_fnl_object` | Main type, all operators and integrators |
| `adam_prism_fnl_library.F90` | `adam_prism_fnl_library` | Barrel re-export of all FNL PRISM modules |
| `adam_prism_fnl_kernels.F90` | `adam_prism_fnl_kernels` | Standalone OpenACC device kernels |
| `adam_prism_fnl_coil_object.F90` | `adam_prism_fnl_coil_object` | GPU coil current source handler |
| `adam_prism_fnl_fWLayer_object.F90` | `adam_prism_fnl_fWLayer_object` | GPU far-wave layer handler |
| `adam_prism_fnl_external_fields_kernels.F90` | `adam_prism_fnl_external_fields_kernels` | GPU external field sub/add kernels |

## Build modes

```bash
FoBiS.py build -mode prism-fnl-nvf-oac        # OpenACC with NVIDIA compiler
FoBiS.py build -mode prism-fnl-nvf-oac-debug  # debug build
```

## GPU memory layout — critical design point

The FNL backend uses the **transposed array layout** `(nb, ni, nj, nk, nv)` for all device
arrays, the inverse of the CPU layout `(nv, ni, nj, nk, nb)`. This is required for coalesced
memory access on NVIDIA GPUs, where adjacent threads execute the innermost loop index:

```fortran
! CPU layout: (nv, ni, nj, nk, nb)  — variable index innermost
q_cpu(v, i, j, k, b)

! GPU layout: (nb, ni, nj, nk, nv)  — variable index outermost
q_gpu(b, i, j, k, v)
```

With `collapse(4)` over `(b, k, j, i)`, threads access `q_gpu(b,i,j,k,v)` with adjacent `i`
values per warp, giving coalesced access in the `ni` dimension.

The `copy_cpu_gpu` / `copy_gpu_cpu` routines perform a full layout transposition on every
host↔device transfer via `field_gpu%copy_transpose_cpu_gpu` / `copy_transpose_gpu_cpu`.

## Type: `prism_fnl_object`

`prism_fnl_object` extends `prism_common_object` and adds:

### GPU-specific ADAM library objects

> **Singleton**: the GPU-aware MPI handler is the program-scope `mpih_fnl` singleton
> (module `adam_fnl_mpih_global`), not an embedded member of `prism_fnl_object`.
> It must be initialized once — via `call mpih_fnl%initialize(do_mpi_init=.true.,
> do_device_init=.true., verbose=.true.)` — at the start of `prism_fnl_object%initialize`,
> before any other FNL object is constructed.

| Member | Type | Description |
|--------|------|-------------|
| `field_gpu` | `field_fnl_object` | Field object with GPU memory and transposed halo exchange |
| `ib_gpu` | `ib_fnl_object` | Immersed Boundary with device-resident `phi_gpu` |
| `rk_gpu` | `rk_fnl_object` | RK stage arrays on device (`q_rk_gpu`) |
| `weno_gpu` | `weno_fnl_object` | WENO coefficients mirrored on device |

### GPU-specific PRISM library objects

| Member | Type | Description |
|--------|------|-------------|
| `coil_gpu` | `prism_fnl_coil_object` | Coil J-vector and flag arrays on device (transposed layout) |
| `fwlayer_gpu` | `prism_fnl_fwlayer_object` | fWLayer damping function `f_gpu` on device |

### Device data pointers

All allocated via FUNDAL `dev_alloc` in `allocate_gpu`:

| Pointer | Shape (GPU transposed) | Description |
|---------|------------------------|-------------|
| `q_gpu` | `(nb, ni, nj, nk, nv)` | Field cell-centred variables (points into `field_gpu%q_gpu`) |
| `dq_gpu` | `(nb, ni, nj, nk, nv)` | Residuals right-hand side |
| `flxyz_c_gpu` | `(nb, 3, 3, ni, nj, nk, nv)` | Cell-centre fluxes with ±-decomposition |
| `flx_f_gpu` | `(nb, 0:ni, nj, nk, nv)` | X-face fluxes |
| `fly_f_gpu` | `(nb, ni, 0:nj, nk, nv)` | Y-face fluxes |
| `flz_f_gpu` | `(nb, ni, nj, 0:nk, nv)` | Z-face fluxes |
| `curl_gpu` | `(nb, ni, nj, nk, nv)` | Curl auxiliary fields |
| `divergence_gpu` | `(nb, ni, nj, nk, nv)` | Divergence auxiliary fields and Poisson buffer |

### Procedure pointers

Same 9 procedure pointer interfaces as the CPU backend; all point to OpenACC-accelerated
implementations operating on GPU arrays.

## Initialization sequence (`initialize`)

```
mpih_fnl%initialize             ← GPU-aware MPI init + CUDA device selection (singleton)
prism_common_object%initialize  ← 34-step common CPU setup; also initializes CPU mpih singleton
field_gpu%initialize(field)     ← mirror global field singleton metadata on device
ib_gpu%initialize               ← copy IB phi to device
rk_gpu%initialize               ← allocate q_rk_gpu stages
weno_gpu%initialize             ← copy WENO coefficients to device
allocate_gpu                    ← allocate all device data pointers
coil_gpu%initialize             ← allocate and transpose coil arrays
fwlayer_gpu%initialize          ← allocate f_gpu
external_fields_initialize_dev  ← copy external field data to device
set procedure pointers          ← based on scheme_space / scheme_time
```

> **Global singletons in scope** (via `adam_fnl_library` → `adam_common_library`):
> `mpih_fnl` (GPU MPI handler), `mpih` (CPU MPI handler), `grid`, `field`, `maps`.
> The global `field` singleton is passed directly to `field_gpu%initialize(field=field, ...)`;
> do not use `self%adam%field` — that pointer is equivalent but less explicit.

## Dispatch table

### Spatial operator pointers

| `scheme_space` | FDV pointer family | `compute_residuals` |
|----------------|-------------------|---------------------|
| `FD_CENTERED` | `*_fd` variants | `compute_residuals_fd_centered` ✓ |
| `FV_CENTERED` | `*_fv` variants | *(not yet implemented)* |
| `WENO` | `*_fv` variants | *(not yet implemented)* |

### Time integration pointer

| `physical_model` | `scheme_time` | RK variant | `integrate` | Status |
|-----------------|--------------|------------|-------------|--------|
| EM | `BLANES_MOAN` | — | `integrate_blanesmoan` | stub |
| EM | `CFM` | — | `integrate_cfm` | stub |
| EM | `LEAPFROG` | — | `integrate_leapfrog` | stub |
| EM | `RK` | `RK_1/2/3` | `integrate_rk_ls` | stub |
| EM | `RK` | `RK_SSP_22/33/54` | `integrate_rk_ssp` | **active** |
| EM | `RK` | `RK_YOSHIDA` | `integrate_rk_yoshida` | stub |
| PIC | `LEAPFROG` | PIC_LEAPFROG/RK | `integrate_leapfrog_pic` | stub |

The only fully implemented time integrator is **SSP RK** (`integrate_rk_ssp`).

## OpenACC parallelization pattern

All FDV operator loops follow the same collapsed parallel pattern:

```fortran
!$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,out_gpu,dxyz_gpu)
do b=1, blocks_number
do k=1, nk
do j=1, nj
do i=1, ni
   ! scalar kernel call
enddo
enddo
enddo
enddo
```

The `DEVICEVAR` clause (FUNDAL macro expanding to OpenACC `present` / `deviceptr` depending on
the SDK version) asserts that the listed arrays are already resident on device. Stencil sub-arrays
are declared `private` when needed to avoid false sharing between threads.

## Methods

### Data management

#### `allocate_gpu`
Allocates all device pointers using `dev_alloc` from FUNDAL. `q_gpu` is a pointer alias into
`field_gpu%q_gpu` (already allocated by `field_gpu%initialize`).

#### `copy_cpu_gpu`
Transfers state from CPU to device after initial conditions or restart load:
- `field_gpu%copy_transpose_cpu_gpu` — transposes and uploads `q`, `curl`, `divergence`, `fWLayer%f`
- `field_gpu%copy_cpu_gpu` — uploads grid metadata (dxyz, cell coordinates, BC maps)
- `coil_gpu%copy_cpu_gpu` — uploads coil amplitude, frequency, phase, J-vector (with transposition)

#### `copy_gpu_cpu`
Downloads state from device to CPU before I/O:
- `field_gpu%copy_transpose_gpu_cpu` — downloads and un-transposes `q`, `curl`, `divergence`, `fWLayer%f`
- `coil_gpu%copy_gpu_cpu` — downloads coil arrays

### Simulation loop (`simulate`)

Mirrors the CPU `simulate` structure with GPU-aware adaptations:

```
initialize
if restart: load_restart_files → copy_cpu_gpu
else:        set_initial_conditions (→ CPU IC + copy_cpu_gpu)

save_simulation_data (→ copy_gpu_cpu before save)
compute_energy (on device)
open residuals file

do
  compute_dt             ← OpenACC reduction on device
  integrate              ← dispatched GPU integrator
  time += dt
  save_simulation_data
  compute_energy
  if converged: exit
enddo

compute_energy_error
save final data
print initial/final energies and RMS errors
MPI + GPU finalize
```

Memory status is dumped for both CPU and GPU when `save_memory_status = .true.`.

### I/O

#### `load_restart_files`
Calls `prism_common_object%load_restart_files` (CPU read), then `copy_cpu_gpu` to upload to device.

#### `save_residuals`
Computes L2 residual norm on device via `compute_normL2_residuals_dev`, reduces with
`MPI_ALLREDUCE(SUM)`, and writes from rank 0. Called only when `is_to_save(residuals_save)`.

#### `save_simulation_data`
Triggered by `it_save`, `restart_save`, or slice monitors. Calls `update_ghost`,
then `copy_gpu_cpu` to bring data to CPU, then `compute_auxiliary_fields`, then saves
XH5F files, restart files, and slices.

### Ghost cell update (`update_ghost`)

Same 3-step async scheme as CPU but on GPU arrays:

| Step | Action |
|------|--------|
| 1 | `field_gpu%update_ghost_local_gpu` — block-local halo fill on device |
| MPI | `field_gpu%update_ghost_mpi_gpu` — GPU-aware MPI halo exchange |
| 3 | `apply_fwl_correction` + `set_boundary_conditions` on device |

#### `update_rk_ghost`
Stub for future implementation of radiative BC ghost updates during RK stages.

### Boundary conditions (`set_boundary_conditions`)

Iterates `local_map_bc_crown_gpu` (device-resident BC map) over `ngc` ghost crowns using a
`!$acc parallel loop independent gang vector` directive. Currently implemented:

| BC type | Status | Implementation |
|---------|--------|----------------|
| `EXTRAPOLATION` | Active | `q_gpu(ghost) = q_gpu(interior)` |
| `NEUMANN` | Active | Mirror reflection formula |
| `DIRICHLET` | Active | Set to zero |
| `Silver_Müller` | Pending | "to be implemented" |
| `PERIODIC` | Pending | "to be implemented" |

### Sources — `compute_coils_current`

Delegates to the `compute_coils_current_dev` kernel in `adam_prism_fnl_kernels` via
`coil_gpu`'s device arrays. The GPU kernel (`adam_prism_fnl_kernels.F90`) applies the same
smooth-start polynomial and branchless step-function arithmetic as the CPU backend, over a
`collapse(4)` gang-vector loop. The optional `gamma` argument shifts the evaluation time for
SSP RK stages.

### Far-wave layer (`apply_fwl_correction`)

Computes per-face index ranges from `fWLayer%C` (layer cell width) and calls
`apply_fwl_correction_dev` for each active face. The device kernel applies the impedance-matching
correction to `(alfa_D, beta_D, alfa_B, beta_B)` component pairs using the fWLayer function `f_gpu`:

$$q_{α_D}' = \frac{1}{2\sqrt{\mu_0}}\!\left[s_2(f-1)\,q_{β_B}\sqrt{\varepsilon_0} + (f+1)\,q_{α_D}\sqrt{\mu_0}\right]$$

and similarly for the other three components (face-dependent sign `s2`).

### Spatial discretization

#### `compute_residuals_fd_centered`
The only fully active space operator. Uses a `!$acc parallel loop collapse(4)` kernel that
explicitly extracts 1D stencil arrays (`qsx_y`, `qsx_z`, etc.) as private temporaries, then calls
the device-side curl kernel `compute_curl_fd_centered_dev`:

$$\frac{\partial \mathbf{D}}{\partial t} = \nabla \times \frac{\mathbf{B}}{\mu_0} - \mathbf{J}, \qquad
\frac{\partial \mathbf{B}}{\partial t} = -\nabla \times \frac{\mathbf{D}}{\varepsilon_0}$$

The stencil arrays are declared as local variables with static bounds at compile time (sized to
`fdv_half_stencils(1)`) to be `private` in the OpenACC region.

**Note**: The NVF SDK 24.x/25.x does not support `firstprivate` with `associate` variables in
`!$acc` regions; the grid scalars (`blocks_number`, `ni`, etc.) are passed via `copyin` instead.

#### FDV operator library

All operators loop `collapse(4)` with `DEVICEVAR` and call the same underlying scalar kernels as
the CPU backend. The `dxyz_gpu(b,:)` array replaces the CPU's `dxyz(:,b)` due to the transposed
layout. Unique to the FNL backend:

- `compute_divergence_fd` implements the finite-difference stencil **inline** using embedded
  coefficient tables `FD1_CC(s_max, s_max)` supporting up to 10th-order schemes, avoiding a
  library call that would require complex array slicing across the transposed layout.
- `compute_divergence_fv` similarly inlines FV reconstruction using
  `compute_reconstruction_r_fv_centered` on 1D stencil temporaries declared `private`.
- Both `divergence` variants take an extra `ovar` parameter selecting the output slot in
  `divergence_gpu(:,:,:,:,ovar)`.
- `compute_laplacian_fv` is fully implemented (unlike the CPU backend stub).

| Operator | FD | FV | Stencil half-width |
|----------|----|----|-------------------|
| `compute_curl` | `compute_curl_fd_centered` | `compute_curl_fv_centered` | `fdv_half_stencils(1)` |
| `compute_derivative1` | `compute_derivative1_fd_centered` | `compute_derivative1_fv_centered` | `fdv_half_stencils(1)` |
| `compute_derivative2` | `compute_derivative2_fd_centered` | `compute_derivative2_fv_centered` | `fdv_half_stencils(1)` |
| `compute_derivative4` | `compute_derivative4_fd_centered` | — | `fdv_half_stencils(1)` |
| `compute_divergence` | inline FD coefficients | inline FV reconstruction | `fdv_half_stencils(1)` |
| `compute_gradient` | `compute_gradient_fd_centered` | `compute_gradient_fv_centered` | `fdv_half_stencils(1)` |
| `compute_laplacian` | `compute_laplacian_fd_centered` | `compute_laplacian_fv_centered` | `fdv_half_stencils(2)` |

### SSP RK time integration (`integrate_rk_ssp`)

The only fully implemented integrator:

```
if external_fields: sub_external_fields_dev
rk_gpu%initialize_stages(q_gpu)
for s = 1..nrk:
  compute_coils_current(gamma=rk%gamm(s))   ← GPU coil kernel
  rk_gpu%compute_stage(s, dt [,phi_gpu])    ← advance q_rk_gpu stage on device
  compute_residuals(q_rk_gpu(:,:,:,:,:,s))  ← OpenACC FD curl kernel
  if s==1: save_residuals
  rk_gpu%assign_stage(s, dq_gpu [,phi_gpu]) ← store stage residuals
rk_gpu%update_q(dt, q_gpu [,phi_gpu])       ← final stage combination
update_rk_ghost(dt)                         ← (stub)
impose_div_free
if external_fields: add_external_fields_dev
```

### `compute_dt`

Same CFL formula as CPU (`dt = CFL × Δx_min / evmax`) but the minimum is computed via an
`!$acc parallel loop reduction(min:dxyz_min)` directly on device, avoiding a device-to-host
transfer. Result is then reduced across MPI ranks via `MPI_ALLREDUCE(MIN)`.

### `compute_energy`

Computes the field energy entirely on device using `!$acc parallel loop reduction(+:energy)`
before reducing across MPI ranks. The energy sums `0.5 × |q|²` over all interior cells for the D
and B fields, exactly as in the CPU backend but without a copy to host.

### `compute_auxiliary_fields`

Currently stubbed — divergence and curl computations for diagnostic fields are commented out
pending stable `compute_divergence_gpu` signatures.

### Divergence cleaning (`impose_div_free` / `impose_ct_correction`)

Same logic as CPU: calls `impose_ct_correction(ivar=1)` for D and/or `ivar=4` for B when Poisson
CT correction is enabled. The inner FLAIL multigrid solver (`compute_smoothing_multigrid`) is
currently commented out pending a GPU port; the gradient correction loop is implemented with an
`!$acc parallel loop collapse(5)` kernel.

## GPU support modules

### `prism_fnl_coil_object`

Manages device-resident coil data. `initialize` allocates device arrays via `dev_alloc`.
`copy_cpu_gpu` manually transposes `coil%J_vec(nv,i,j,k,nb)` → `J_vec_gpu(nb,i,j,k,nv)` and
`coil%coil_flag(i,j,k,nb)` → `coil_flag_gpu(nb,i,j,k)` before uploading.

### `prism_fnl_fwlayer_object`

Stores `f_gpu(nb,i,j,k,3)` (three directional fWLayer function values) allocated on device.
The `apply_fwl_correction_dev` kernel is direction-agnostic: the caller passes the appropriate
index ranges and Barbas-notation field indices for the target face.

### `adam_prism_fnl_kernels`

Standalone module containing public device kernels:

| Subroutine | Description |
|-----------|-------------|
| `compute_coils_current_dev` | Coil current injection with smooth start transient |
| `compute_div_d_b_dev` | Divergence of D and B fields |
| `compute_fluxes_convective_dev` | Convective flux computation (WENO, device) |
| `compute_fluxes_difference_dev` | Flux divergence + source subtraction |
| `set_bc_q_gpu_dev` | Apply boundary conditions on device |
| `set_sir_dev` | Initialize direction versors on device |

### `adam_prism_fnl_external_fields_kernels`

Provides `external_fields_initialize_dev`, `sub_external_fields_dev`, and
`add_external_fields_dev` for GPU-side external field handling in SSP RK integrations.

## Implementation status summary

| Feature | Status |
|---------|--------|
| FD centered space operator | Active |
| FV centered space operator | Not implemented |
| WENO space operator | Not implemented |
| SSP RK time integration | Active |
| Leapfrog / Yoshida / Blanes-Moan / CFM | Stubbed (CPU logic commented out) |
| PIC integration | Stubbed |
| Silver-Müller BC | Pending |
| Periodic BC | Pending |
| FLAIL multigrid (divergence cleaning) | Partially ported |
| Auxiliary fields (curl/divergence diagnostic) | Pending |
