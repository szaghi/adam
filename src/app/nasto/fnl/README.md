# NASTO — FNL Backend

> GPU-accelerated backend: MPI + OpenACC via the FUNDAL device abstraction library.

The `fnl/` directory provides the production GPU implementation of the NASTO
solver. It extends `nasto_common_object` with GPU-resident counterpart objects
and device kernels written through the
[FUNDAL](https://github.com/szaghi/FUNDAL) library, which abstracts OpenACC
and OpenMP target offloading behind a single set of macros and allocation
routines.

## Source layout

| Source file | Module / Program | Role |
|-------------|-----------------|------|
| `adam_nasto_fnl.F90` | program `adam_nasto_fnl` | Entry point — same CLI as CPU backend |
| `adam_nasto_fnl_object.F90` | module `adam_nasto_fnl_object` | Backend object — GPU objects, data management, all methods |
| `adam_nasto_fnl_kernels.F90` | module `adam_nasto_fnl_kernels` | Device kernel library — flux loops with FUNDAL directives |
| `adam_nasto_fnl_cns_kernels.F90` | module `adam_nasto_fnl_cns_kernels` | Device-callable CNS routines (`!$acc routine` / `!$omp declare target`) |
| `adam_nasto_fnl_library.F90` | module `adam_nasto_fnl_library` | Barrel re-export of all FNL and common modules |

---

## FUNDAL and OpenACC conventions

FUNDAL is included via a preprocessor header (`#include "fundal.H"`) that
provides two key abstractions:

**`DEVICEVAR(...)`** — declares a list of dummy arguments as device-resident,
replacing explicit `!$acc copyin`/`copyout` clauses. Used in every parallel
loop:

```fortran
!$acc parallel loop independent DEVICEVAR(q_aux_gpu, fluxes_gpu, weno_a_gpu)
!$omp OMPLOOP              DEVICEVAR(q_aux_gpu, fluxes_gpu, weno_a_gpu)
do b = 1, blocks_number
  ...
```

The same source compiles for OpenACC (nvfortran) and OpenMP target (any
offloading compiler) via macro expansion.

**`dev_alloc` / `dev_assign_to_device`** — FUNDAL API for device memory
allocation and host-to-device copy, returning Fortran pointers to
GPU-resident data:

```fortran
call dev_alloc(fptr_dev=self%q_aux_gpu, &
               ubounds=[nb, ni+ngc, nj+ngc, nk+ngc, nv_aux], &
               lbounds=[1, 1-ngc, 1-ngc, 1-ngc, 1],          &
               init_value=0._R8P, ierr=ierr)

call dev_assign_to_device(dst=self%q_bc_vars_gpu, src=self%bc%q)
```

**Device-callable procedures** are annotated with both directives so they can
be invoked from within device parallel regions:

```fortran
!$acc routine(compute_conservatives_dev)
!$omp declare target(compute_conservatives_dev)
```

---

## Array layout transposition

The FNL backend stores field data in a **transposed layout** relative to the
CPU backend, optimised for GPU thread coalescing:

| Backend | Layout | Stride-1 index | Rationale |
|---------|--------|----------------|-----------|
| CPU | `q(nv, ni, nj, nk, nb)` | `nv` (variables) | SIMD vectorisation over variables |
| FNL/GPU | `q_gpu(nb, ni, nj, nk, nv)` | `nb` (blocks) | GPU coalescing: adjacent threads access adjacent blocks |

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

## `nasto_fnl_object`

Extends `nasto_common_object` with GPU companion objects for every
compute-intensive ADAM SDK component and a set of device pointer working
arrays.

### Inheritance

```
nasto_common_object          (common/)
    └── nasto_fnl_object     (fnl/)   ← backend object
```

### GPU companion objects

Each CPU SDK object has a GPU counterpart holding device-resident data and
device-side methods:

| CPU object (inherited) | GPU companion | Role |
|------------------------|--------------|------|
| `mpih` (`mpih_object`) | `mpih_gpu` (`mpih_fnl_object`) | Device initialisation + GPU-aware MPI barriers |
| `field` (`field_object*`) | `field_gpu` (`field_fnl_object`) | Ghost-cell exchange on GPU; `x/y/z_cell_gpu`, `dxyz_gpu` |
| `ib` (`ib_object`) | `ib_gpu` (`ib_fnl_object`) | `phi_gpu` device array; eikonal on device |
| `rk` (`rk_object`) | `rk_gpu` (`rk_fnl_object`) | RK stage arrays on device (`q_rk_gpu`) |
| `weno` (`weno_object`) | `weno_gpu` (`weno_fnl_object`) | WENO coefficient arrays on device (`a_gpu`, `p_gpu`, `d_gpu`) |

### GPU working arrays (device pointers)

All arrays are `pointer` to FUNDAL-allocated device memory:

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
| `initialize(filename)` | Calls `mpih_gpu%initialize(do_device_init=.true.)`; `initialize_common`; initialises all GPU companion objects; calls `allocate_gpu` |
| `allocate_gpu(q_gpu)` | Assigns `q_gpu` pointer; calls `dev_alloc` for `q_aux_gpu`, `dq_gpu`, `flx/fly/flz_gpu`; copies BC state to device |
| `copy_cpu_gpu()` | Transpose-copies `q` → `q_gpu`; copies grid/field metadata to device |
| `copy_gpu_cpu([compute_copy_q_aux] [,copy_phi])` | Transpose-copies `q_gpu` → `q`; optionally computes and copies `q_aux_gpu`, and copies `phi_gpu` → `phi` |
| `simulate(filename)` | Full simulation driver (identical structure to CPU backend) |

**AMR** — requires a CPU roundtrip because the AMR tree lives on the host:

| Procedure | Purpose |
|-----------|---------|
| `amr_update()` | Ghost sync on GPU → mark blocks (gradient on device) → `copy_gpu_cpu` → `adam%amr_update` → `copy_cpu_gpu` → recompute `phi` |
| `mark_by_geo(delta_fine, delta_coarse)` | CPU-side marking using `phi` from the last GPU→CPU copy |
| `mark_by_grad_var(grad_tol, …)` | Ghost sync + `compute_q_auxiliary` on GPU; gradient computed via `field_gpu%compute_q_gradient` |
| `refine_uniform(levels)` | Uniform refinement; delegates to `adam%amr_update` |
| `compute_phi()` | Runs `compute_phi_analytical_sphere_dev` on device; applies `reduce_cell_order_phi_dev` near solids; builds aggregate `phi_gpu` across all solids |
| `move_phi(velocity)` | Calls `move_phi_dev` — GPU kernel translating the IB distance field |

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
| `save_residuals()` | Computes L2 norm of `dq_gpu` via `compute_normL2_residuals_dev` on GPU; reduces with `MPI_ALLREDUCE` |
| `save_simulation_data()` | Schedules output: syncs ghost cells on GPU → `copy_gpu_cpu(compute_copy_q_aux=.true., copy_phi=.true.)` → saves HDF5/restart/slices |

**Initial and boundary conditions**

| Procedure | Purpose |
|-----------|---------|
| `set_initial_conditions()` | Calls `ic%set_initial_conditions` on CPU; then `copy_cpu_gpu` |
| `set_boundary_conditions(q_gpu)` | Calls `set_bc_q_gpu_dev` using `field_gpu%maps%local_map_bc_crown_gpu` — all on device |
| `update_ghost(q_gpu [,step])` | GPU-local inter-block copy → GPU-aware MPI exchange → BC application; optional async `step` parameter |

**Numerical methods** — all primary compute on GPU:

| Procedure | Purpose |
|-----------|---------|
| `compute_dt()` | Calls `compute_umax_dev` on GPU using `dxyz_gpu`; global minimum via `MPI_ALLREDUCE` |
| `compute_q_auxiliary(q_gpu, q_aux_gpu)` | Calls `compute_q_aux_dev` on GPU |
| `compute_residuals(q_gpu, dq_gpu)` | Ghost sync → eikonal → aux vars → convective fluxes → diffusive correction → flux divergence; all on GPU |
| `integrate()` | RK time advance on GPU via `rk_gpu`; same scheme dispatch (low-storage / SSP) as CPU backend |

### Residual computation pipeline (GPU)

```
1. update_ghost(q_gpu)                      ! GPU ghost sync + BC
2. integrate_eikonal(q_gpu)                 ! eikonal sweeps on device
3. compute_q_auxiliary(q_gpu, q_aux_gpu)    ! q_gpu → primitive on device
4. compute_fluxes_convective_dev(x,y,z)     ! WENO + LLF on device
5. compute_fluxes_diffusive_dev(…)          ! viscous stress tensor on device (if μ > 0)
6. compute_fluxes_difference_dev(…)         ! flux divergence → dq_gpu; IB correction
```

---

## `adam_nasto_fnl_kernels`

Device kernel library. All procedures accept only primitive-type arguments
(no derived types) and annotate parallel regions with FUNDAL macros.

### Public procedures

**Convective fluxes**

```fortran
subroutine compute_fluxes_convective_dev(dir, blocks_number, ni, nj, nk, ngc, nv, &
                                         weno_s, weno_a_gpu, weno_p_gpu, weno_d_gpu, &
                                         weno_zeps, ror_number, ror_schemes_gpu, &
                                         ror_threshold, ror_ivar_gpu, ror_stats_gpu, &
                                         g, q_aux_gpu, fluxes_gpu)
```

Launches a `!$acc parallel loop independent` region over the 4D cell space
(`b, k, j, i` with `b` innermost for coalescing). Calls
`compute_fluxes_convective_ri_dev` (from `adam_fnl_weno_kernels`) at each
interface, which performs WENO reconstruction with optional
Reduced-Order-Reconstruction (ROR) near shocks.

**Diffusive fluxes**

```fortran
subroutine compute_fluxes_diffusive_dev(null_x, null_y, null_z, blocks_number, &
                                        ni, nj, nk, ngc, nv, mu, kd,           &
                                        q_aux_gpu, dx_gpu, dy_gpu, dz_gpu,      &
                                        flx_gpu, fly_gpu, flz_gpu)
```

Computes the full viscous stress tensor $\tau_{ij}$ and heat flux at each
face interface, with cross-derivative terms approximated by 4-point stencils.
Adds the result to the convective flux arrays in-place. Skipped entirely when
$\mu = 0$.

**Flux divergence**

```fortran
subroutine compute_fluxes_difference_dev(null_x, null_y, null_z, blocks_number, &
                                         ni, nj, nk, ngc, nv, ib_eps,           &
                                         dx_gpu, dy_gpu, dz_gpu,                 &
                                         flx_gpu, fly_gpu, flz_gpu [, phi_gpu],  &
                                         dq_gpu)
```

Computes the conservative flux divergence into `dq_gpu`. When `phi_gpu` is
present (IB run), adjusts the effective cell width at cut cells based on the
signed-distance function, replacing $\Delta x$ with the distance to the
immersed boundary surface.

**Auxiliary variables and CFL**

| Procedure | Purpose |
|-----------|---------|
| `compute_q_aux_dev(…, q_gpu, q_aux_gpu)` | Full-field conservative → primitive conversion on device |
| `compute_umax_dev(…, dxyz_gpu, q_aux_gpu, umax)` | CFL + Fourier wave-speed reduction on device |
| `set_bc_q_gpu_dev(…, l_map_bc_gpu, q_bc_vars_gpu, q_gpu)` | Applies extrapolation and inflow BCs on device using pre-built crown maps |

---

## `adam_nasto_fnl_cns_kernels`

Single-cell device-callable procedures. Each routine is annotated with
`!$acc routine(…)` and `!$omp declare target(…)` so it can be called from
within a `!$acc parallel loop` or `!$omp target teams` region without
launching a new kernel.

| Procedure | Purpose |
|-----------|---------|
| `compute_conservatives_dev(b,i,j,k,ngc,q_aux_gpu,q)` | Single-cell `q_aux_gpu → q` |
| `compute_conservative_fluxes_dev(sir,b,i,j,k,ngc,q_aux_gpu,f)` | Single-cell Euler flux in direction `sir` |
| `compute_max_eigenvalues_dev(si,sir,weno_s,…,q_aux_gpu,evmax)` | LLF wave-speed estimate over WENO stencil |
| `compute_eigenvectors_dev(…,el,er)` | Flux Jacobian eigenvector matrices |
| `compute_q_aux_dev(ni,nj,nk,ngc,nb,R,cv,g,q_gpu,q_aux_gpu)` | Full-field `q_gpu → q_aux_gpu` (also exported by `adam_nasto_fnl_kernels`) |

---

## `adam_nasto_fnl_library`

Barrel re-export module. A single `use adam_nasto_fnl_library` statement
brings everything into scope:

```fortran
use adam_nasto_common_library   ! all common/ modules
use adam_nasto_fnl_cns_kernels  ! device-callable CNS routines
use adam_nasto_fnl_kernels      ! device kernel library
use adam_fnl_library            ! FUNDAL SDK: dev_alloc, field_fnl_object, …
```

---

## Build

```bash
# OpenACC — NVIDIA HPC SDK (nvfortran)
FoBiS.py build -mode nasto-fnl-nvf-oac

# Debug build
FoBiS.py build -mode nasto-fnl-nvf-oac-debug
```

The executable `adam_nasto_fnl` is placed in `exe/`.

## Running

```bash
# One MPI rank per GPU
mpirun -np 4 exe/adam_nasto_fnl

# With custom INI file
mpirun -np 4 exe/adam_nasto_fnl my_case.ini
```

Set the GPU-to-MPI-rank binding before launch:

```bash
export CUDA_VISIBLE_DEVICES=0,1,2,3   # or use cluster launcher binding
mpirun -np 4 --map-by ppr:1:gpu exe/adam_nasto_fnl
```

---

## Copyrights

NASTO is part of the ADAM framework, released under the
[GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html)
(LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
