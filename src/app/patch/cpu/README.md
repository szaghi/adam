# PATCH — CPU Backend

The CPU backend runs PATCH on multi-core systems using MPI for distributed-memory parallelism and OpenMP (via the SOR-OMP smoother) for shared-memory parallelism within each rank.

## Source Layout

| File | Module / Program | Role |
|------|-----------------|------|
| `adam_patch_cpu.F90` | `program adam_patch_cpu` | Entry point |
| `adam_patch_cpu_object.F90` | `module adam_patch_cpu_object` | `patch_cpu_object` type and implementation |

## Type Hierarchy

```
patch_common_object          (src/app/patch/common/)
└── patch_cpu_object         (src/app/patch/cpu/)
```

`patch_cpu_object` extends `patch_common_object` with the residual scratch array and the full simulation pipeline.

## Entry Point

`adam_patch_cpu.F90` follows the standard ADAM CLI pattern:

```fortran
program adam_patch_cpu
  use adam_patch_cpu_object, only : patch_cpu_object
  type(patch_cpu_object) :: patch

  na = command_argument_count()
  if (na == 0) then
    call patch%simulate(filename='input.ini')
  else
    call get_command_argument(1, input_file_name)
    call patch%simulate(filename=trim(adjustl(input_file_name)))
  endif
end program
```

Run as:

```bash
mpirun -np <N> ./adam_patch_cpu [input.ini]
```

If no argument is given, `input.ini` in the working directory is used.

## Additional Data Member

`patch_cpu_object` adds one scratch array to those inherited from `patch_common_object`:

| Member | Shape | Description |
|--------|-------|-------------|
| `dq` | `(nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb)` | Residual / smoother scratch array |

`dq` serves dual purpose: it is the residual buffer passed to FLAIL smoothing kernels, and in the Multigrid V-cycle it holds the error correction between restriction and prolongation.

## Methods

### Lifecycle

| Procedure | Description |
|-----------|-------------|
| `allocate_cpu` | Allocates `dq` with ghost-cell bounds; zeros it |
| `initialize` | Calls `initialize_common`, allocates CPU data, loads restart or applies IC, saves initial snapshot, opens residuals log |
| `finalize` | Saves final snapshot, closes residuals log, finalises MPI |

### AMR and Immersed Boundary

| Procedure | Description |
|-----------|-------------|
| `amr_update` | Iterative AMR marker loop — runs until the grid stabilises or the AMR iteration budget is exhausted |
| `compute_phi` | Computes signed-distance (level-set) field φ for all IB bodies |
| `mark_by_geo` | Marks blocks for refinement/de-refinement based on IB surface proximity (φ sign change) |
| `mark_by_grad_var` | Marks blocks based on field-gradient threshold (skeleton, not yet active) |
| `move_phi` | Translates φ by a given body velocity vector |
| `refine_uniform` | Forces uniform refinement to a specified number of levels |
| `integrate_eikonal` | Advances the eikonal equation to propagate φ outward from IB surfaces |

### I/O

| Procedure | Description |
|-----------|-------------|
| `load_restart_files` | Loads `q` and grid state from restart files; rebuilds communication maps |
| `save_xh5f` | Writes φ field in XH5F format; output basename is `{basename}-{it:09}` |
| `save_residuals` | Computes MPI-reduced L2 norm of `dq`, normalised by cell count; written by rank 0 |
| `save_restart_files` | Writes checkpoint (grid + field binary + XH5F snapshot) |
| `save_simulation_data` | Orchestrates periodic saves based on `it_save` and `restart_save` frequencies |

### Initial and Boundary Conditions

| Procedure | Description |
|-----------|-------------|
| `set_initial_conditions` | Sets `r` (source distribution) using the selected IC type |
| `set_boundary_conditions` | Applies Dirichlet BCs to ghost cells via the `local_map_bc_crown` map |
| `update_ghost` | Exchanges halos (local + MPI) and applies BCs; supports optional async `step` argument |

### Numerical Kernel

| Procedure | Description |
|-----------|-------------|
| `integrate` | Convergence iteration loop; dispatches to the selected FLAIL smoother |
| `simulate` | Top-level driver: initialise → select smoother → integrate → finalise |

## Key Implementation Details

### `initialize`

```
mpih%initialize(do_mpi_init=.true.)
initialize_common(field=adam%field, filename, memory_avail)
allocate_cpu
if io%restart:
  load_restart_files → time%it, time%time
else:
  set_initial_conditions
  time%time = 0;  time%it = 0
save_simulation_data        ! initial snapshot at t=0
io%open_file_residuals(nv)  ! rank 0 only
```

### `amr_update`

```
outer: do i = 1, amr%iters
  is_grid_changed_all = .false.
  do i_marker = 1, amr%markers_number
    update_ghost(q)
    select case(amr_marker%mode)
      case(AMR_GEO):  mark_by_geo(delta_fine, delta_coarse)
      case(AMR_GRAD): mark_by_grad_var(grad_tol, delta_type, ...)
    endselect
    adam%amr_update(is_marked_by_field=.true., ...)
    if ib%solids_number > 0: compute_phi
    is_grid_changed_all = is_grid_changed_all .or. is_grid_changed
  enddo
  if .not. is_grid_changed_all: EXIT   ! grid stabilised
enddo outer
```

### `update_ghost`

Ghost exchange supports optional asynchronous stepping via the `step` argument:

| `step` value | Action |
|---|---|
| absent | Local update + MPI exchange + apply BCs (fully synchronous) |
| 1 | Local update + start MPI exchange |
| 3 | Apply BCs (called after MPI exchange completes) |

### `integrate` — Convergence Loop

```
integration: do
  time%it += 1

  call compute_smoothing(ni, nj, nk, ngc, nv, blocks_number,
                         dxyz=field%dxyz, f=r, q=q, dq=dq, dq_max=dq_max,
                         iterations_init, iterations_fine, iterations_coarse)

  if dq_max < 1e-6:
    print 'convergence reached at iteration it'
    EXIT integration

  dq = q                           ! store current solution in dq scratch
  time%print_progress(nodes_number=adam%tree%nodes_number)
  save_simulation_data

  if (it_max <= 0 and time >= time_max) or
     (it_max >  0 and it  >= it_max  ): EXIT integration
enddo integration
```

The convergence tolerance (1×10⁻⁶) is hardcoded in `integrate`; the FLAIL `tolerance` field is used only within individual smoothing procedures to tune internal sub-iterations.

### `simulate` — Smoother Dispatch

```fortran
call self%initialize(filename=filename)

select case(self%flail%smoothing)
case(SMOOTHING_MULTIGRID)    ; call self%integrate(compute_smoothing=compute_smoothing_multigrid)
case(SMOOTHING_GAUSS_SEIDEL) ; call self%integrate(compute_smoothing=compute_smoothing_gauss_seidel)
case(SMOOTHING_SOR)          ; call self%integrate(compute_smoothing=compute_smoothing_sor)
case(SMOOTHING_SOR_OMP)      ; call self%integrate(compute_smoothing=compute_smoothing_sor_omp)
endselect

call self%finalize
```

The smoother is passed as a Fortran procedure argument conforming to the `compute_smoothing_interface` abstract interface defined in `adam_flail_object`.

## Building

```bash
# Release build (GNU)
FoBiS.py build -mode patch-gnu

# Debug build
FoBiS.py build -mode patch-gnu-debug
```

The executable is written to `exe/adam_patch_cpu`.

## License

PATCH is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
