# NASTO — NVF Backend

> GPU-accelerated backend: MPI + CUDA Fortran via the `cudafor` intrinsic module.

The `nvf/` directory provides the CUDA Fortran GPU implementation of the NASTO
solver. It extends `nasto_common_object` with GPU-resident device arrays and
kernels written using native CUDA Fortran language features: `attributes(global)`
global kernels, `attributes(device)` device procedures, `!$cuf kernel` directives,
and the `blockDim` / `blockIdx` / `threadIdx` intrinsics.

## Source layout

| Source file | Module / Program | Role |
|-------------|-----------------|------|
| `adam_nasto_nvf.F90` | program `adam_nasto_nvf` | Entry point — same CLI as CPU backend |
| `adam_nasto_nvf_object.F90` | module `adam_nasto_nvf_object` | Backend object — GPU device arrays, data management, all methods |
| `adam_nasto_nvf_kernels.F90` | module `adam_nasto_nvf_kernels` | Device kernel library — CUDA global kernels and CUF directives |
| `adam_nasto_nvf_cns_kernels.F90` | module `adam_nasto_nvf_cns_kernels` | Device-callable CNS routines (`attributes(device)`) |

---

## CUDA Fortran conventions

The NVF backend uses standard CUDA Fortran (`use cudafor`) rather than an
abstraction library.

**`attributes(global)` kernels** — launched with an explicit grid/block
configuration; threads are identified by the built-in `blockDim`, `blockIdx`,
and `threadIdx` intrinsics:

```fortran
attributes(global) subroutine compute_fluxes_convective_kernel(dir, blocks_number, ni, nj, nk, &
                                                               ngc, nv, weno_s, ..., q_aux_gpu, fluxes_gpu)
integer(I4P), intent(in), value  :: dir, blocks_number, ni, ...  ! scalars by value
real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,...)     ! device arrays

b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
k = blockDim%z * (blockIdx%z - 1) + threadIdx%z
if (b > blocks_number .or. j > nj .or. k > nk) return
```

**`!$cuf kernel do(N) <<<*,*>>>`** — shorthand CUF directive for simpler loops;
the compiler auto-generates the launch configuration:

```fortran
!$cuf kernel do(4) <<<*,*>>>
do b = 1, blocks_number
  do k = 1-ngc, nk+ngc
    do j = 1-ngc, nj+ngc
      do i = 1-ngc, ni+ngc
        dq_gpu(b,i,j,k,iv) = ...
      enddo
    enddo
  enddo
enddo
```

**`attributes(device)` procedures** — called from within device kernels;
scalar dummy arguments use the `value` attribute (CUDA calling convention):

```fortran
attributes(device) subroutine compute_conservatives_device(b, i, j, k, ngc, q_aux_gpu, q)
integer(I4P), intent(in), value  :: b, i, j, k, ngc    ! by value on GPU
real(R8P),    intent(in), device :: q_aux_gpu(1:,...)   ! device array
```

**`alloc_var_gpu` / `assign_allocatable_gpu`** — NVF SDK helper routines from
`adam_nvf_library` for device allocation (replaces FUNDAL's `dev_alloc`):

```fortran
call alloc_var_gpu(var=self%q_aux_gpu, &
                   ulb=reshape([1,nb, 1-ngc,ni+ngc, 1-ngc,nj+ngc, 1-ngc,nk+ngc, 1,nv_aux], [2,5]), &
                   msg=msg_//' q_aux_gpu ')

call assign_allocatable_gpu(lhs=self%q_bc_vars_gpu, rhs=self%bc%q, msg=msg_//' q_bc_vars_gpu ')
```

**`check_cuda_error`** — called after every GPU operation to surface device errors:

```fortran
call self%mpih_gpu%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in compute phi analytical sphere')
```

---

## Array layout transposition

The NVF backend uses the same transposed layout as the FNL backend, optimised
for GPU thread coalescing:

| Backend | Layout | Stride-1 index | Rationale |
|---------|--------|----------------|-----------|
| CPU | `q(nv, ni, nj, nk, nb)` | `nv` (variables) | SIMD vectorisation over variables |
| NVF/GPU | `q_gpu(nb, ni, nj, nk, nv)` | `nb` (blocks) | GPU coalescing: adjacent threads access adjacent blocks |

The innermost kernel loop runs over `b` (blocks), so adjacent GPU threads
access consecutive `q_gpu(b, i, j, k, v)` elements — stride-1 access in the
`nb` dimension. An explicit transpose is performed at every host↔device data
transfer:

```fortran
! CPU → GPU  (transpose q_cpu → q_gpu)
call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%field%q, q_gpu=self%q_gpu)

! GPU → CPU  (transpose q_gpu → q_cpu)
call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
```

---

## `nasto_nvf_object`

Extends `nasto_common_object` with GPU companion objects and device-allocated
working arrays.

### Inheritance

```
nasto_common_object          (common/)
    └── nasto_nvf_object     (nvf/)   ← backend object
```

### GPU companion objects

| CPU object (inherited) | GPU companion | Role |
|------------------------|--------------|------|
| `mpih` (`mpih_object`) | `mpih_gpu` (`mpih_nvf_object`) | Device init + CUDA error checking + GPU-aware MPI barriers |
| `field` (`field_object`) | `field_gpu` (`field_nvf_object`) | Ghost-cell exchange on GPU; `x/y/z_cell_gpu`, `dxyz_gpu` |
| `ib` (`ib_object`) | `ib_gpu` (`ib_nvf_object`) | `phi_gpu` device array; eikonal on device |
| `rk` (`rk_object`) | `rk_gpu` (`rk_nvf_object`) | RK stage arrays on device (`q_rk_gpu`) |
| `weno` (`weno_object`) | `weno_gpu` (`weno_nvf_object`) | WENO coefficient arrays on device (`a_gpu`, `p_gpu`, `d_gpu`) |

### Device data members

All arrays are CUDA Fortran device-resident:

| Member | Attribute | GPU shape | Description |
|--------|-----------|-----------|-------------|
| `q_gpu` | `pointer, device` | `(nb, 1-ngc:ni+ngc, …, nv)` | Conservative variables (alias to `field_gpu%q_gpu`) |
| `q_aux_gpu` | `allocatable, device` | `(nb, 1-ngc:ni+ngc, …, nv_aux)` | Primitive / auxiliary variables |
| `dq_gpu` | `allocatable, device` | `(nb, 1-ngc:ni+ngc, …, nv)` | Residual (RHS) |
| `flx_gpu` | `allocatable, device` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in x |
| `fly_gpu` | `allocatable, device` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in y |
| `flz_gpu` | `allocatable, device` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in z |
| `q_bc_vars_gpu` | `allocatable, device` | `(6, 6)` | Boundary condition primitive states (copied once at init) |

### Methods

**Lifecycle**

| Procedure | Purpose |
|-----------|---------|
| `initialize(filename)` | Calls `mpih_gpu%initialize(do_device_init=.true.)`; `initialize_common`; `mpih_gpu%load_from_file` (reads CUDA settings from INI); initialises GPU companion objects; calls `allocate_gpu` |
| `allocate_gpu(q_gpu)` | Assigns `q_gpu` pointer; calls `alloc_var_gpu` for `q_aux_gpu`, `dq_gpu`, `flx/fly/flz_gpu`; copies BC state to device via `assign_allocatable_gpu` |
| `copy_cpu_gpu()` | Transpose-copies `q` → `q_gpu`; copies grid/field metadata to device |
| `copy_gpu_cpu([compute_copy_q_aux] [,copy_phi])` | Transpose-copies `q_gpu` → `q`; optionally computes and copies `q_aux_gpu`, and copies `phi_gpu` → `phi` |
| `simulate(filename)` | Full simulation driver (identical structure to CPU backend) |

**AMR** — requires a CPU roundtrip because the AMR tree lives on the host:

| Procedure | Purpose |
|-----------|---------|
| `amr_update()` | Ghost sync on GPU → mark blocks → `copy_gpu_cpu` → `adam%amr_update` → `copy_cpu_gpu` → recompute `phi` |
| `mark_by_geo(delta_fine, delta_coarse)` | CPU-side marking using `phi` from the last GPU→CPU copy |
| `mark_by_grad_var(grad_tol, …)` | Ghost sync + `compute_q_auxiliary` on GPU; gradient computed via `field_gpu%compute_q_gradient`; checks `check_cuda_error` |
| `refine_uniform(levels)` | Uniform refinement; delegates to `adam%amr_update` |
| `compute_phi()` | Runs `compute_phi_analytical_sphere_cuf` on device; applies `reduce_cell_order_phi_cuf` near solids; builds aggregate `phi_gpu` via `compute_phi_all_solids_cuf`; checks `check_cuda_error` after each step |
| `move_phi(velocity)` | Calls `move_phi_cuf` — CUDA kernel translating the IB distance field |

**Immersed boundary**

| Procedure | Purpose |
|-----------|---------|
| `integrate_eikonal(q_gpu)` | Runs `n_eikonal` Godunov sweeps on device via `ib_gpu%evolve_eikonal`; applies `invert_eikonal` on device |

**I/O** — save path always copies GPU → CPU before writing:

| Procedure | Purpose |
|-----------|---------|
| `save_hdf5([output_basename])` | Writes `field%q` and `q_aux` (CPU copies) to HDF5; optionally includes `ib%phi` |
| `save_restart_files()` | Saves ADAM binary restart + HDF5 snapshot from CPU data |
| `load_restart_files(t, time)` | Reads restart binary into CPU arrays; rebuilds MPI maps; calls `copy_cpu_gpu` |
| `save_residuals()` | Computes L2 norm of `dq_gpu` via `compute_normL2_residuals_cuf` on GPU; reduces with `MPI_ALLREDUCE` |
| `save_simulation_data()` | Schedules output: syncs ghost cells on GPU → `copy_gpu_cpu(compute_copy_q_aux=.true., copy_phi=.true.)` → saves HDF5/restart/slices |

**Initial and boundary conditions**

| Procedure | Purpose |
|-----------|---------|
| `set_initial_conditions()` | Calls `ic%set_initial_conditions` on CPU; then `copy_cpu_gpu` |
| `set_boundary_conditions(q_gpu)` | Calls `set_bc_q_gpu_cuf` using `field_gpu%maps%local_map_bc_crown_gpu` — all on device |
| `update_ghost(q_gpu [,step])` | GPU-local inter-block copy → GPU-aware MPI exchange → BC application; optional async `step` parameter |

**Numerical methods** — all primary compute on GPU:

| Procedure | Purpose |
|-----------|---------|
| `compute_dt()` | Calls `compute_umax_cuf` on GPU; global minimum via `MPI_ALLREDUCE` |
| `compute_q_auxiliary(q_gpu, q_aux_gpu)` | Calls `compute_q_aux_cuf` on GPU |
| `compute_residuals(q_gpu, dq_gpu)` | Ghost sync → eikonal → aux vars → convective fluxes → diffusive correction → flux divergence; all on GPU |
| `integrate()` | RK time advance on GPU via `rk_gpu`; same scheme dispatch (low-storage / SSP) as CPU backend |

### Residual computation pipeline (GPU)

```
1. update_ghost(q_gpu)                      ! GPU ghost sync + BC
2. integrate_eikonal(q_gpu)                 ! eikonal sweeps on device
3. compute_q_auxiliary(q_gpu, q_aux_gpu)    ! q_gpu → primitive on device
4. compute_fluxes_convective_kernel(x,y,z)  ! WENO + LLF — global CUDA kernel
5. compute_fluxes_diffusive_cuf(…)          ! viscous stress tensor — CUF directive (if μ > 0)
6. compute_fluxes_difference_cuf(…)         ! flux divergence → dq_gpu; IB correction
```

---

## `adam_nasto_nvf_kernels`

Device kernel library. All procedures accept only primitive-type arguments
(no derived types). Uses `cudafor` intrinsics directly.

### Public procedures

**Convective fluxes**

```fortran
attributes(global) subroutine compute_fluxes_convective_kernel(dir, blocks_number, ni, nj, nk, ngc, nv, &
                                                               weno_s, weno_a_gpu, weno_p_gpu, weno_d_gpu, &
                                                               weno_zeps, ror_number, ror_schemes_gpu, &
                                                               ror_threshold, ror_ivar_gpu, ror_stats_gpu, &
                                                               g, q_aux_gpu, fluxes_gpu)
```

A `attributes(global)` kernel launched with a 3D grid over `(b, j_or_i, k)`
dimensions; the remaining direction is looped inside each thread. Calls the
`attributes(device)` routine `compute_fluxes_convective_device` (from
`adam_weno_nvf_kernels`) at each interface, which performs WENO reconstruction
with optional Reduced-Order-Reconstruction (ROR) near shocks.

**Diffusive fluxes**

```fortran
subroutine compute_fluxes_diffusive_cuf(null_x, null_y, null_z, blocks_number, &
                                        ni, nj, nk, ngc, nv, mu, kd,           &
                                        q_aux_gpu, dx_gpu, dy_gpu, dz_gpu,      &
                                        flx_gpu, fly_gpu, flz_gpu)
```

Computes the full viscous stress tensor $\tau_{ij}$ and heat flux at each
face interface, with cross-derivative terms approximated by 4-point stencils.
Adds the result to the convective flux arrays in-place. Implemented via
`!$cuf kernel do(4) <<<*,*>>>` directive. Skipped entirely when $\mu = 0$.

**Flux divergence**

```fortran
subroutine compute_fluxes_difference_cuf(null_x, null_y, null_z, blocks_number, &
                                         ni, nj, nk, ngc, nv, ib_eps,           &
                                         dx_gpu, dy_gpu, dz_gpu,                 &
                                         flx_gpu, fly_gpu, flz_gpu [, phi_gpu],  &
                                         dq_gpu)
```

Computes the conservative flux divergence into `dq_gpu`. When `phi_gpu` is
present (IB run), adjusts the effective cell width at cut cells based on the
signed-distance function.

**Auxiliary variables and CFL**

| Procedure | Purpose |
|-----------|---------|
| `compute_q_aux_cuf(…, q_gpu, q_aux_gpu)` | Full-field conservative → primitive conversion on device |
| `compute_umax_cuf(…, dxyz_gpu, q_aux_gpu, umax)` | CFL + Fourier wave-speed reduction on device |
| `set_bc_q_gpu_cuf(…, l_map_bc_gpu, q_bc_vars_gpu, q_gpu)` | Applies extrapolation and inflow BCs on device using pre-built crown maps |

---

## `adam_nasto_nvf_cns_kernels`

Single-cell device-callable procedures. Each routine has the `attributes(device)`
qualifier so it can be called from within a `attributes(global)` kernel or a
`!$cuf kernel` region. Scalar dummy arguments use the `value` attribute
(CUDA device calling convention).

| Procedure | Purpose |
|-----------|---------|
| `compute_conservatives_device(b,i,j,k,ngc,q_aux_gpu,q)` | Single-cell `q_aux_gpu → q` |
| `compute_conservative_fluxes_device(sir,b,i,j,k,ngc,q_aux_gpu,f)` | Single-cell Euler flux in direction `sir` |
| `compute_max_eigenvalues_device(si,sir,weno_s,…,q_aux_gpu,evmax)` | LLF wave-speed estimate over WENO stencil |
| `compute_eigenvectors_device(…,el,er)` | Flux Jacobian eigenvector matrices (full Roe-average implementation) |
| `compute_q_aux_device(ni,nj,nk,ngc,nb,R,cv,g,q_gpu,q_aux_gpu)` | Full-field `q_gpu → q_aux_gpu` |

---

## Build

```bash
# NVIDIA CUDA Fortran (production)
FoBiS.py build -mode nasto-nvf-cuda

# Debug build
FoBiS.py build -mode nasto-nvf-cuda-debug
```

The executable `adam_nasto_nvf` is placed in `exe/`.

## Running

```bash
# One MPI rank per GPU
mpirun -np 4 exe/adam_nasto_nvf

# With custom INI file
mpirun -np 4 exe/adam_nasto_nvf my_case.ini
```

Set the GPU-to-MPI-rank binding before launch:

```bash
export CUDA_VISIBLE_DEVICES=0,1,2,3   # or use cluster launcher binding
mpirun -np 4 --map-by ppr:1:gpu exe/adam_nasto_nvf
```

---

## Copyrights

NASTO is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
