# CHASE — CPU Backend

The CPU backend runs CHASE on multi-core systems using MPI for distributed-memory parallelism and OpenMP for shared-memory parallelism within each rank. It is the reference backend for correctness validation.

## Source Layout

| File | Module / Program | Role |
|------|-----------------|------|
| `adam_chase_cpu.F90` | `program adam_chase_cpu` | Entry point |
| `adam_chase_cpu_object.F90` | `module adam_chase_cpu_object` | `chase_cpu_object` type and implementation |

## Type Hierarchy

```
chase_common_object          (src/app/chase/common/)
└── chase_cpu_object         (src/app/chase/cpu/)
```

`chase_cpu_object` extends `chase_common_object` with CPU-resident scratch arrays and the full simulation pipeline.

## Entry Point

`adam_chase_cpu.F90` follows the standard ADAM CLI pattern:

```fortran
program adam_chase_cpu
  use adam_chase_cpu_object
  type(chase_cpu_object) :: chase
  character(len=:), allocatable :: ini_file

  call get_command_argument(1, ini_file)
  if (.not. allocated(ini_file)) ini_file = 'input.ini'

  call chase%simulate(ini_file)
end program
```

Run as:

```bash
mpirun -np <N> ./adam_chase_cpu [input.ini]
```

If no argument is given, `input.ini` in the working directory is used.

## Data Members

`chase_cpu_object` adds the following scratch arrays to those inherited from `chase_common_object`:

| Member | Shape | Description |
|--------|-------|-------------|
| `dq`  | `(nv, ni, nj, nk, nb)` | Accumulated residuals (RHS) |
| `flx` | `(nv, ni, nj, nk, nb)` | Convective fluxes — x-direction |
| `fly` | `(nv, ni, nj, nk, nb)` | Convective fluxes — y-direction |
| `flz` | `(nv, ni, nj, nk, nb)` | Convective fluxes — z-direction |

All arrays follow the standard 5D layout `(nv, ni, nj, nk, nb)` with Fortran column-major ordering (variable index `nv` is stride-1 in memory).

## Methods

### Lifecycle

| Procedure | Description |
|-----------|-------------|
| `allocate_cpu` | Allocates `dq`, `flx`, `fly`, `flz` after grid dimensions are known |
| `initialize` | Reads INI file, initialises MPI, grid, physics, BC/IC, eikonal, and allocates scratch |

### AMR and Immersed Boundary

| Procedure | Description |
|-----------|-------------|
| `amr_update` | Refines / coarsens the AMR tree, redistributes blocks across ranks |
| `compute_phi` | Computes the signed-distance (level-set) field φ for all IB bodies |
| `mark_by_geo` | Marks AMR cells for refinement based on geometry proximity |
| `mark_by_grad_var` | Marks cells for refinement based on variable-gradient threshold |
| `move_phi` | Updates φ and geometry data after body motion |
| `refine_uniform` | Forces uniform refinement to a specified level |
| `integrate_eikonal` | Advances the eikonal equation to propagate φ |

### I/O

| Procedure | Description |
|-----------|-------------|
| `load_restart_files` | Reads `q`, `phi`, and grid from restart files |
| `save_xh5f` | Writes field data in XH5F format for visualisation |
| `save_residuals` | Appends iteration residuals to the residuals log file |
| `save_restart_files` | Writes checkpoint files (field + grid state) |
| `save_simulation_data` | Orchestrates all periodic I/O calls |

### Initial and Boundary Conditions

| Procedure | Description |
|-----------|-------------|
| `set_initial_conditions` | Applies the selected IC type to `q` |
| `set_boundary_conditions` | Applies physical BCs to domain boundaries each stage |
| `update_ghost` | Exchanges ghost-cell halos via MPI between ranks |

### Numerical Kernel

| Procedure | Description |
|-----------|-------------|
| `compute_dt` | Computes the global CFL time step |
| `compute_q_auxiliary` | Derives auxiliary variables from the conservative state `q` |
| `compute_residuals` | Assembles the full spatial residual into `dq` |
| `integrate` | Advances `q` in time by one step using the selected RK scheme |
| `simulate` | Outer time-loop driver |

## Time Step

CHASE uses an inviscid (CFL-only) time step — there is no Fourier viscous stability limit:

$$\Delta t = \text{CFL} \cdot \min_{b,i,j,k} \frac{\Delta x}{\lvert u \rvert + c}$$

where $c = \sqrt{\gamma p / \rho}$ is the local speed of sound and the minimum is taken over all active cells and directions.

## Residual Pipeline

`compute_residuals` assembles `dq` in five stages:

1. **Ghost exchange** — `update_ghost` synchronises halo cells across MPI ranks.
2. **Eikonal update** — `integrate_eikonal` propagates the signed-distance field φ.
3. **Auxiliary variables** — `compute_q_auxiliary` fills `q_aux` from the current `q`.
4. **Convective fluxes** — `compute_fluxes_convective` (module-private) computes WENO flux vectors for each active spatial direction and stores them in `flx`, `fly`, `flz`.
5. **Flux divergence** — `compute_fluxes_difference` (module-private) accumulates `flx/fly/flz` into `dq`, applying the IB cut-cell width correction via φ.

Directions flagged by `null_xyz` are skipped: their flux arrays are zeroed and the corresponding momentum components are not updated.

## Module-Private Helpers

These non-type-bound procedures implement the innermost numerical kernels:

### `assign_omp_R8P_5D`

OpenMP `collapse(5)` utility to zero-fill (or assign) a 5D real array:

```fortran
!$omp parallel do collapse(5) default(none) ...
do ib=1,nb; do k=1,nk; do j=1,nj; do i=1,ni; do iv=1,nv
  array(iv,i,j,k,ib) = 0.0_R8P
enddo; enddo; enddo; enddo; enddo
```

### `compute_fluxes_convective`

Outer driver for one spatial direction (`dir = 1, 2, 3`). Loops over blocks and cells, calling `compute_fluxes_convective_ri` at each cell interface.

### `compute_fluxes_convective_ri` — WENO Reconstruction Kernel

Per-interface convective flux computation using WENO with optional eigenvector projection:

```fortran
! 1. Compute left/right eigenvector matrices at the interface
call compute_eigenvectors(si, sir, b, i, j, k, ngc, nv, g, q_aux, el, er)

! 2. Project fluxes to characteristic space and apply LLF flux splitting
call decompose_fluxes_convective(weno_s, b, i, j, k, ngc, nv, si, sir, el, q_aux, fmpc)

! 3. WENO reconstruction in characteristic variables
do v = 1, nv
  call weno_reconstruct_upwind(S, weno_a, weno_p, weno_d, weno_zeps, V=fmpc(:,:,v), VR=fpmr(:,v))
enddo

! 4. Back-project to conservative space
do v = 1, nv
  fluxes(v,i,j,k,b) = sum(er(:,v) * (fpmr(1,:) + fpmr(2,:)))
enddo
```

When the WENO basis is `CONSERVATIVE`, `el` and `er` are identity matrices and the projection is bypassed. When `CHARACTERISTICS`, the Euler eigenvectors stored in `adam_chase_parameters` are used.

### `compute_fluxes_difference`

Accumulates the flux divergence into `dq`. For IB cut cells the effective cell width is adjusted using the φ sign-distance field, matching the cut-cell treatment used in NASTO:

```
Δx_eff = Δx * H(φ)    (Heaviside of signed distance)
dq(v,i,j,k,b) += (flx_right - flx_left) / Δx_eff
```

Momentum components in nullified directions (`null_xyz`) are multiplied by zero.

## Runge-Kutta Integration

`integrate` selects the time-advancement scheme at runtime via `self%rk%scheme`:

| Scheme | Identifier | Stages | Type | Storage |
|--------|-----------|--------|------|---------|
| Euler (1st order) | `RK_1` | 1 | Low-storage | In-place `q` update |
| 2nd-order | `RK_2` | 2 | Low-storage | In-place `q` update |
| 3rd-order | `RK_3` | 3 | Low-storage | In-place `q` update |
| SSP-RK(2,2) | `RK_SSP_22` | 2 | SSP | Stages stored in `q_rk` |
| SSP-RK(3,3) | `RK_SSP_33` | 3 | SSP | Stages stored in `q_rk` |
| SSP-RK(5,4) | `RK_SSP_54` | 5 | SSP | Stages stored in `q_rk` |

Residuals (`dq`) are saved at the first stage (`s == 1`) for the convergence log.

## `simulate` — Time-Loop Structure

```
initialize(ini_file)
  ↳ read INI, init MPI, grid, physics, BC/IC, eikonal, allocate scratch
set_initial_conditions  or  load_restart_files
save_simulation_data(t=0)

loop:
  if mod(it, amr%frequency) == 0:  amr_update
  compute_dt
  dt = min(dt, time_max - time)    ! clamp to reach time_max exactly
  integrate(dt)                    ! compute_residuals + RK stages
  time += dt;  it += 1
  save_simulation_data(time, it)
  if time >= time_max  or  it >= it_max:  EXIT

save_simulation_data(final)
close residuals file
mpih%finalize
```

## Building

```bash
# Release build (GNU)
FoBiS.py build -mode chase-gnu

# Debug build
FoBiS.py build -mode chase-gnu-debug
```

The executable is written to `exe/adam_chase_cpu`.

## License

CHASE is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
