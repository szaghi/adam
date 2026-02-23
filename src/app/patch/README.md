# PATCH

**PATCH** — **P**oisson solver with **A**dap**t**ive mes**h** refinement for HPC computing — is a distributed-memory elliptic solver built on the ADAM framework. It solves the Poisson equation for a scalar potential field using iterative smoothing methods on AMR grids.

## Governing Equation

PATCH solves the Poisson boundary-value problem:

$$\nabla^2 \varphi = f$$

where $\varphi$ is the scalar potential (output field) and $f$ is a prescribed source distribution. Convergence is measured by the point-wise residual:

$$r_{i,j,k} = f_{i,j,k} - \nabla^2 \varphi_{i,j,k}$$

The solver iterates until $\max_{i,j,k} |r_{i,j,k}| < \varepsilon$ or a maximum iteration count is reached.

## Source Layout

```
src/app/patch/
├── common/                              # Shared across all backends
│   ├── adam_patch_common_object.F90     # patch_common_object base type
│   ├── adam_patch_bc_object.F90         # Boundary conditions handler
│   ├── adam_patch_ic_object.F90         # Initial conditions handler
│   ├── adam_patch_io_object.F90         # I/O handler
│   ├── adam_patch_time_object.F90       # Time / convergence handler
│   └── adam_patch_common_library.F90    # Barrel re-export of all common modules
└── cpu/                                 # CPU backend
    ├── adam_patch_cpu.F90               # Entry point (program)
    └── adam_patch_cpu_object.F90        # patch_cpu_object implementation
```

## Field Variables

| Array | Shape | Name | Description |
|-------|-------|------|-------------|
| `q`   | `(nv, ni, nj, nk, nb)` | `phi` | Scalar potential (1 variable, `nv=1`) |
| `r`   | `(nv, ni, nj, nk, nb)` | `rho` | Source / forcing distribution |
| `dq`  | `(nv, ni, nj, nk, nb)` | — | Residual scratch array |

All arrays use the standard ADAM 5D layout `(nv, ni, nj, nk, nb)` with Fortran column-major ordering.

## Type Hierarchy

```
patch_common_object          (src/app/patch/common/)
└── patch_cpu_object         (src/app/patch/cpu/)
```

### `patch_common_object` Data Members

| Member | Type | Description |
|--------|------|-------------|
| `mpih` | `mpih_object` | MPI handler |
| `adam` | `adam_object` | Core ADAM object (tree, grid, field, I/O) |
| `field` | `field_object*` | Pointer to `adam%field` |
| `grid` | `grid_object*` | Pointer to `adam%grid` |
| `amr` | `amr_object` | AMR marker handler |
| `ib` | `ib_object` | Immersed boundary handler |
| `flail` | `flail_object` | Linear algebra / smoothing handler (FLAIL) |
| `io` | `patch_io_object` | I/O configuration |
| `ic` | `patch_ic_object` | Initial conditions |
| `bc` | `patch_bc_object` | Boundary conditions |
| `time` | `patch_time_object` | Convergence / time parameters |
| `q` | `real(R8P)(:,:,:,:,:)` | Scalar potential field |
| `r` | `real(R8P)(:,:,:,:,:)` | Source distribution |

### `patch_cpu_object` Additional Member

| Member | Type | Description |
|--------|------|-------------|
| `dq` | `real(R8P)(:,:,:,:,:)` | Residual scratch array |

## Boundary Conditions

The only supported boundary condition is **Dirichlet** — the potential is set to a constant value on each domain face.

Six independent faces are configured in the INI file:

| INI section | Face |
|-------------|------|
| `[bc_x_min]` | x-minimum face |
| `[bc_x_max]` | x-maximum face |
| `[bc_y_min]` | y-minimum face |
| `[bc_y_max]` | y-maximum face |
| `[bc_z_min]` | z-minimum face |
| `[bc_z_max]` | z-maximum face |

```ini
[bc_x_min]
type = dirichlet
phi  = 0.0

[bc_x_max]
type = dirichlet
phi  = 0.0
```

## Initial Conditions

The only supported IC type sets a Gaussian distribution of the source term:

| IC type | `type` key | Description |
|---------|-----------|-------------|
| Gaussian | `gauss` | $r = \exp\!\left(-\tfrac{(x-0.5)^2+(y-0.5)^2+(z-0.5)^2}{0.01}\right)$ centred at $(0.5, 0.5, 0.5)$ |

```ini
[initial_conditions]
type           = gauss
amr_iterations = 2
regions_number = 1

[initial_conditions_region_1]
r      = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
```

## Iterative Smoothing Solvers (FLAIL)

The FLAIL library (`adam_flail_object`) provides the iterative smoothing kernels. The solver is selected in the INI file and dispatched via a Fortran procedure-pointer argument to `integrate`.

| Smoothing method | `[linear-algebra] smoothing` key | Algorithm |
|-----------------|----------------------------------|-----------|
| Gauss-Seidel | `GAUSS-SEIDEL` | Sequential point update: $\varphi^{n+1}_{i,j,k} = \frac{1}{2(1/\Delta x^2 + 1/\Delta y^2 + 1/\Delta z^2)}\left(\text{neighbour sum} - f_{i,j,k}\right)$ |
| SOR | `SOR` | Gauss-Seidel with over-relaxation: $\varphi^{n+1} = \varphi^n + \omega(\varphi^{\text{GS}} - \varphi^n)$, $\omega=1.8$ |
| SOR-OMP | `SOR-OMP` | Parallel SOR with 8-colour red-black ordering (`mod(i+j+k, 8)`) and OpenMP `parallel do` |
| Multigrid | `MULTIGRID` | Two-level V-cycle: pre-smooth (Gauss-Seidel) → restrict residuals (full-weighting) → coarse solve (Gauss-Seidel) → prolongate correction (trilinear) → post-smooth |

FLAIL iteration counts are tunable per smoothing level:

| INI key | Description |
|---------|-------------|
| `iterations` | General iteration count |
| `iterations_init` | Pre-smoothing iterations (Multigrid: initial guess on fine grid) |
| `iterations_fine` | Post-smoothing iterations on fine grid |
| `iterations_coarse` | Smoothing iterations on coarse grid |
| `tolerance` | Maximum residual tolerance for convergence |

```ini
[linear-algebra]
smoothing          = MULTIGRID
iterations         = 1
iterations_init    = 3
iterations_fine    = 3
iterations_coarse  = 10
tolerance          = 1e-6
```

## Methods

### Lifecycle

| Procedure | Description |
|-----------|-------------|
| `allocate_common` | Allocates `r` with ghost-cell bounds |
| `initialize_common` | Reads INI, initialises MPI, grid, AMR tree, field, IB, FLAIL, allocates arrays |
| `allocate_cpu` | Allocates `dq` |
| `initialize` | Calls `initialize_common`, allocates CPU data, loads restart or applies IC, saves initial snapshot |
| `finalize` | Saves final data, closes residuals file, finalises MPI |

### AMR and Immersed Boundary

| Procedure | Description |
|-----------|-------------|
| `amr_update` | Runs AMR marker loop (geometry and gradient criteria), calls `adam%amr_update`, recomputes φ |
| `compute_phi` | Computes signed-distance (level-set) field φ for all IB bodies |
| `mark_by_geo` | Marks blocks for refinement/de-refinement based on IB proximity (φ sign change) |
| `mark_by_grad_var` | Marks blocks based on field-gradient threshold (skeleton, not yet active) |
| `move_phi` | Translates φ by a given body velocity |
| `refine_uniform` | Forces uniform refinement to a given level |
| `integrate_eikonal` | Advances the eikonal equation to propagate φ and inverts it |

### I/O

| Procedure | Description |
|-----------|-------------|
| `load_restart_files` | Loads `q` and grid state from restart files |
| `save_xh5f` | Writes field data in XH5F format |
| `save_residuals` | Computes L2 norm of `dq` via `MPI_ALLREDUCE` and appends to `{basename}-residuals.dat` |
| `save_restart_files` | Writes checkpoint (grid + field + XH5F snapshot) |
| `save_simulation_data` | Orchestrates periodic saves based on `it_save` and `restart_save` frequencies |

### Initial and Boundary Conditions

| Procedure | Description |
|-----------|-------------|
| `set_initial_conditions` | Sets `r` to the selected IC distribution (Gaussian) |
| `set_boundary_conditions` | Applies Dirichlet BCs to ghost cells via the `local_map_bc_crown` map |
| `update_ghost` | Exchanges halos (local + MPI) and applies BCs |

### Numerical Kernel

| Procedure | Description |
|-----------|-------------|
| `integrate` | Convergence iteration loop; dispatches to the selected smoothing procedure |
| `simulate` | Top-level driver: initialise → select smoother → integrate → finalise |

## Convergence Loop

`integrate` runs the selected smoothing procedure until the residual drops below tolerance or the iteration budget is exhausted:

```
loop:
  it += 1
  call compute_smoothing(ni,nj,nk,ngc,nv,blocks_number, dxyz, r, q, dq=dq, dq_max=dq_max,
                         iterations_init, iterations_fine, iterations_coarse)
  if dq_max < tolerance:  EXIT  (convergence)
  dq = q                         ! store current solution for residuals log
  print_progress
  save_simulation_data
  if it >= it_max:  EXIT         (iteration budget)
```

`simulate` selects the smoothing procedure at startup:

```fortran
select case(self%flail%smoothing)
case(SMOOTHING_MULTIGRID)    ; call self%integrate(compute_smoothing=compute_smoothing_multigrid)
case(SMOOTHING_GAUSS_SEIDEL) ; call self%integrate(compute_smoothing=compute_smoothing_gauss_seidel)
case(SMOOTHING_SOR)          ; call self%integrate(compute_smoothing=compute_smoothing_sor)
case(SMOOTHING_SOR_OMP)      ; call self%integrate(compute_smoothing=compute_smoothing_sor_omp)
endselect
```

## Building and Running

```bash
# Release build (GNU)
FoBiS.py build -mode patch-gnu

# Debug build
FoBiS.py build -mode patch-gnu-debug
```

The executable is written to `exe/adam_patch_cpu`.

```bash
mpirun -np <N> exe/adam_patch_cpu [input.ini]
```

If no argument is given, `input.ini` in the working directory is used.

## License

PATCH is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
