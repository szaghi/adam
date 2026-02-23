# NASTO — GMP Backend

> GPU-accelerated backend: MPI + OpenMP target offloading via the `adam_gmp_library`.

> **Status**: In active development. API and build modes may change between releases.

The `gmp/` directory provides the OpenMP target GPU implementation of the NASTO
solver. It extends `nasto_common_object` with GPU-resident device arrays and
compute loops written using OpenMP 5.x target directives:
`!$omp target teams distribute parallel do` with `has_device_addr` for
zero-copy access to already-mapped device data, and `!$omp declare target`
for device-callable procedures.

## Source layout

| Source file | Module / Program | Role |
|-------------|-----------------|------|
| `adam_nasto_gmp.F90` | program `adam_nasto_gmp` | Entry point — same CLI as CPU backend |
| `adam_nasto_gmp_object.F90` | module `adam_nasto_gmp_object` | Backend object — device pointers, data management, all methods |
| `adam_nasto_gmp_kernels.F90` | module `adam_nasto_gmp_kernels` | Device kernel library — OpenMP target loops |
| `adam_nasto_gmp_cns_kernels.F90` | module `adam_nasto_gmp_cns_kernels` | Device-callable CNS routines (`!$omp declare target`) |

---

## OpenMP target conventions

The GMP backend uses standard OpenMP 5.x target offloading. No additional
abstraction library is required beyond `adam_gmp_library` for device memory
management.

**`!$omp target teams distribute parallel do collapse(N) has_device_addr(...)`** —
the primary parallel loop construct. `has_device_addr` asserts that the listed
arrays are already resident on the device (allocated via `alloc_var_gpu`),
avoiding implicit `map(tofrom:)` copies at each kernel invocation:

```fortran
!$omp target teams distribute parallel do collapse(4) has_device_addr(weno_a_gpu,weno_p_gpu,  &
!$omp&                                                                q_aux_gpu,fluxes_gpu)
do b = 1, blocks_number
do k = 1-si(3), nk
do j = 1-si(2), nj
do i = 1-si(1), ni
   call compute_fluxes_convective_device(...)
enddo
enddo
enddo
enddo
```

**`!$omp target data map(to:...)`** — used to map small local arrays (e.g., stencil
direction vectors) for a scoped set of target regions before entering a loop.

**`!$omp declare target(...)`** — marks procedures as device-callable so they can
be invoked from within target regions:

```fortran
subroutine compute_conservatives_device(b, i, j, k, ngc, q_aux_gpu, q)
!$omp declare target(compute_conservatives_device)
  ...
```

**`alloc_var_gpu` / `assign_allocatable_gpu`** — GMP SDK helper routines from
`adam_gmp_library` for device-side allocation, taking an explicit `omp_dev`
device ID:

```fortran
call alloc_var_gpu(var=self%q_aux_gpu,                                                         &
                   ulb=reshape([1,nb, 1-ngc,ni+ngc, 1-ngc,nj+ngc, 1-ngc,nk+ngc, 1,nv_aux], [2,5]), &
                   omp_dev=mydev, init_val=0._R8P, msg=ms)

call assign_allocatable_gpu(lhs=self%q_bc_vars_gpu, rhs=self%bc%q, omp_dev=mydev, msg=msg_//' q_bc_vars_gpu ')
```

Unlike the NVF backend, device arrays are plain Fortran `pointer` variables —
there is no `device` attribute. The OpenMP runtime manages the device residency
of the underlying storage.

---

## Array layout transposition

The GMP backend uses the same transposed layout as all GPU backends, optimised
for coalescing:

| Backend | Layout | Stride-1 index | Rationale |
|---------|--------|----------------|-----------|
| CPU | `q(nv, ni, nj, nk, nb)` | `nv` (variables) | SIMD vectorisation over variables |
| GMP/GPU | `q_gpu(nb, ni, nj, nk, nv)` | `nb` (blocks) | GPU coalescing: adjacent threads access adjacent blocks |

An explicit transpose is performed at every host↔device data transfer:

```fortran
! CPU → GPU  (transpose q_cpu → q_gpu)
call self%field_gpu%copy_transpose_cpu_gpu(nv=self%nv, q_cpu=self%field%q, q_gpu=self%q_gpu)

! GPU → CPU  (transpose q_gpu → q_cpu)
call self%field_gpu%copy_transpose_gpu_cpu(nv=self%nv, q_gpu=self%q_gpu, q_cpu=self%field%q)
```

---

## `nasto_gmp_object`

Extends `nasto_common_object` with GPU companion objects and device-allocated
working arrays.

### Inheritance

```
nasto_common_object          (common/)
    └── nasto_gmp_object     (gmp/)   ← backend object
```

### GPU companion objects

| CPU object (inherited) | GPU companion | Role |
|------------------------|--------------|------|
| `mpih` (`mpih_object`) | `mpih_gpu` (`mpih_gmp_object`) | Device initialisation; exposes `mydev` (device ID) |
| `field` (`field_object`) | `field_gpu` (`field_gmp_object`) | Ghost-cell exchange on GPU; `x/y/z_cell_gpu`, `dxyz_gpu` |
| `ib` (`ib_object`) | `ib_gpu` (`ib_gmp_object`) | `phi_gpu` device array; eikonal on device |
| `rk` (`rk_object`) | `rk_gpu` (`rk_gmp_object`) | RK stage arrays on device (`q_rk_gpu`) |
| `weno` (`weno_object`) | `weno_gpu` (`weno_gmp_object`) | WENO coefficient arrays on device (`a_gpu`, `p_gpu`, `d_gpu`) |

All companion objects receive an explicit `mpih=self%mpih_gpu` argument during
initialisation so they can query the device ID and issue device-aware messages.

### Device data members

All arrays are plain Fortran `pointer` variables whose storage is allocated by
`alloc_var_gpu` and mapped to `mpih_gpu%mydev` via the OpenMP runtime:

| Member | GPU shape | Description |
|--------|-----------|-------------|
| `q_gpu` | `(nb, 1-ngc:ni+ngc, …, nv)` | Conservative variables (alias to `field_gpu%q_gpu`) |
| `q_aux_gpu` | `(nb, 1-ngc:ni+ngc, …, nv_aux)` | Primitive / auxiliary variables |
| `dq_gpu` | `(nb, 1-ngc:ni+ngc, …, nv)` | Residual (RHS) |
| `flx_gpu` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in x |
| `fly_gpu` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in y |
| `flz_gpu` | `(nb, 1-ngc:ni+ngc, …, nv)` | Face fluxes in z |
| `q_bc_vars_gpu` | `(6, 6)` | Boundary condition primitive states (copied once at init) |

### Methods

**Lifecycle**

| Procedure | Purpose |
|-----------|---------|
| `initialize(filename)` | Calls `mpih_gpu%initialize(do_device_init=.true.)`; `initialize_common`; initialises all GPU companion objects (passing `mpih=self%mpih_gpu`); calls `allocate_gpu` |
| `allocate_gpu(q_gpu)` | Assigns `q_gpu` pointer; calls `alloc_var_gpu` (with `omp_dev=mydev`) for `q_aux_gpu`, `dq_gpu`, `flx/fly/flz_gpu`; copies BC state via `assign_allocatable_gpu` |
| `copy_cpu_gpu()` | Transpose-copies `q` → `q_gpu`; copies grid/field metadata to device |
| `copy_gpu_cpu([compute_copy_q_aux] [,copy_phi])` | Transpose-copies `q_gpu` → `q`; optionally computes and copies `q_aux_gpu`, and copies `phi_gpu` → `phi` |
| `simulate(filename)` | Full simulation driver (identical structure to CPU backend) |

**AMR** — requires a CPU roundtrip because the AMR tree lives on the host:

| Procedure | Purpose |
|-----------|---------|
| `amr_update()` | Ghost sync on GPU → mark blocks → `copy_gpu_cpu` → `adam%amr_update` → `copy_cpu_gpu` → recompute `phi` |
| `mark_by_geo(delta_fine, delta_coarse)` | CPU-side marking using `phi` from the last GPU→CPU copy |
| `mark_by_grad_var(grad_tol, …)` | Ghost sync + `compute_q_auxiliary` on GPU; gradient computed via `field_gpu%compute_q_gradient` |
| `refine_uniform(levels)` | Uniform refinement; delegates to `adam%amr_update` |
| `compute_phi()` | Runs `compute_phi_analytical_sphere_gmp` on device; applies `reduce_cell_order_phi_gmp` near solids; builds aggregate `phi_gpu` via `compute_phi_all_solids_gmp` |
| `move_phi(velocity)` | Calls `move_phi_gmp` — OpenMP target kernel translating the IB distance field |

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
| `save_residuals()` | Computes L2 norm of `dq_gpu` via `compute_normL2_residuals_gmp` on GPU; reduces with `MPI_ALLREDUCE` |
| `save_simulation_data()` | Schedules output: syncs ghost cells on GPU → `copy_gpu_cpu(compute_copy_q_aux=.true., copy_phi=.true.)` → saves HDF5/restart/slices |

**Initial and boundary conditions**

| Procedure | Purpose |
|-----------|---------|
| `set_initial_conditions()` | Calls `ic%set_initial_conditions` on CPU; then `copy_cpu_gpu` |
| `set_boundary_conditions(q_gpu)` | Calls `set_bc_q_gpu_gmp` using `field_gpu%maps%local_map_bc_crown_gpu` — all on device |
| `update_ghost(q_gpu [,step])` | GPU-local inter-block copy → GPU-aware MPI exchange → BC application; optional async `step` parameter |

**Numerical methods** — all primary compute on GPU:

| Procedure | Purpose |
|-----------|---------|
| `compute_dt()` | Calls `compute_umax_gmp` on GPU; global minimum via `MPI_ALLREDUCE` |
| `compute_q_auxiliary(q_gpu, q_aux_gpu)` | Calls `compute_q_aux_gmp` on GPU |
| `compute_residuals(q_gpu, dq_gpu)` | Ghost sync → eikonal → aux vars → convective fluxes → diffusive correction → flux divergence; all on GPU |
| `integrate()` | RK time advance on GPU via `rk_gpu`; same scheme dispatch (low-storage / SSP) as CPU backend |

### Residual computation pipeline (GPU)

```
1. update_ghost(q_gpu)                      ! GPU ghost sync + BC
2. integrate_eikonal(q_gpu)                 ! eikonal sweeps on device
3. compute_q_auxiliary(q_gpu, q_aux_gpu)    ! q_gpu → primitive on device
4. compute_fluxes_convective_gmp(x,y,z)     ! WENO + LLF — omp target teams collapse(4)
5. compute_fluxes_diffusive_gmp(…)          ! viscous stress tensor — omp target (if μ > 0)
6. compute_fluxes_difference_gmp(…)         ! flux divergence → dq_gpu; IB correction
```

---

## `adam_nasto_gmp_kernels`

Device kernel library. All procedures accept only primitive-type arguments
(no derived types) and annotate parallel regions with OpenMP target directives.

### Public procedures

**Convective fluxes**

```fortran
subroutine compute_fluxes_convective_gmp(dir, blocks_number, ni, nj, nk, ngc, nv, &
                                         weno_s, weno_a_gpu, weno_p_gpu, weno_d_gpu, &
                                         weno_zeps, ror_number, ror_schemes_gpu, &
                                         ror_threshold, ror_ivar_gpu, ror_stats_gpu, &
                                         g, q_aux_gpu, fluxes_gpu)
```

Launches a `!$omp target teams distribute parallel do collapse(4)` region over
the 4D cell space (`b, k, j, i` with `b` innermost for coalescing). Uses
`has_device_addr` to avoid data movement for already-resident arrays. Calls
`compute_fluxes_convective_device` (from `adam_nasto_gmp_cns_kernels`) at each
interface, which performs WENO reconstruction with optional
Reduced-Order-Reconstruction (ROR) near shocks.

**Diffusive fluxes**

```fortran
subroutine compute_fluxes_diffusive_gmp(blocks_number, ni, nj, nk, ngc, nv, mu, kd, &
                                        q_aux_gpu, dx_gpu, dy_gpu, dz_gpu,            &
                                        flx_gpu, fly_gpu, flz_gpu)
```

Computes the full viscous stress tensor $\tau_{ij}$ and heat flux at each
face interface, with cross-derivative terms approximated by 4-point stencils.
Three separate `!$omp target teams distribute parallel do collapse(4)` loops
(one per direction) add the result to the convective flux arrays in-place.
Skipped entirely when $\mu = 0$.

**Flux divergence**

```fortran
subroutine compute_fluxes_difference_gmp(blocks_number, ni, nj, nk, ngc, nv, ib_eps, &
                                         dx_gpu, dy_gpu, dz_gpu,                       &
                                         flx_gpu, fly_gpu, flz_gpu [, phi_gpu],        &
                                         dq_gpu)
```

Computes the conservative flux divergence into `dq_gpu`. When `phi_gpu` is
present (IB run), adjusts the effective cell width at cut cells based on the
signed-distance function, replacing $\Delta x$ with the distance to the
immersed boundary surface.

**Auxiliary variables and CFL**

| Procedure | Purpose |
|-----------|---------|
| `compute_q_aux_gmp(…, q_gpu, q_aux_gpu)` | Full-field conservative → primitive conversion on device |
| `compute_umax_gmp(…, dxyz_gpu, q_aux_gpu, umax)` | CFL + Fourier wave-speed reduction on device |
| `set_bc_q_gpu_gmp(…, l_map_bc_gpu, q_bc_vars_gpu, q_gpu)` | Applies extrapolation and inflow BCs on device using pre-built crown maps |

---

## `adam_nasto_gmp_cns_kernels`

Single-cell device-callable procedures. Each routine is annotated with
`!$omp declare target(...)` so it can be called from within an
`!$omp target teams` region without launching a new kernel.

| Procedure | Purpose |
|-----------|---------|
| `compute_conservatives_device(b,i,j,k,ngc,q_aux_gpu,q)` | Single-cell `q_aux_gpu → q` |
| `compute_conservative_fluxes_device(sir,b,i,j,k,ngc,q_aux_gpu,f)` | Single-cell Euler flux in direction `sir` |
| `compute_max_eigenvalues_device(si,sir,weno_s,…,q_aux_gpu,evmax)` | LLF wave-speed estimate over WENO stencil |
| `compute_eigenvectors_device(…,el,er)` | Flux Jacobian eigenvector matrices (Roe-average implementation) |
| `compute_q_aux_gmp(ni,nj,nk,ngc,nb,R,cv,g,q_gpu,q_aux_gpu)` | Full-field `q_gpu → q_aux_gpu` (also exported by `adam_nasto_gmp_kernels`) |

---

## Build

```bash
# Intel oneAPI (ifx) — OpenMP target offloading
FoBiS.py build -mode nasto-gmp-ifx

# Debug build
FoBiS.py build -mode nasto-gmp-ifx-debug
```

The executable `adam_nasto_gmp` is placed in `exe/`.

## Running

```bash
# One MPI rank per GPU
mpirun -np 4 exe/adam_nasto_gmp

# With custom INI file
mpirun -np 4 exe/adam_nasto_gmp my_case.ini
```

Set the device binding before launch (Intel GPU example):

```bash
export ZE_AFFINITY_MASK=0,1,2,3   # Intel Level Zero device selection
mpirun -np 4 exe/adam_nasto_gmp
```

---

## Copyrights

NASTO is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
