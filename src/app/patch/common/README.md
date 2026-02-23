# PATCH — Common Modules

The `common/` directory contains the modules shared by all PATCH backends. They define the problem physics, boundary and initial conditions, I/O, convergence parameters, and the base type that every backend extends.

## Source Layout

| File | Module | Role |
|------|--------|------|
| `adam_patch_common_object.F90` | `adam_patch_common_object` | `patch_common_object` base type and `initialize_common` |
| `adam_patch_bc_object.F90` | `adam_patch_bc_object` | Boundary conditions handler |
| `adam_patch_ic_object.F90` | `adam_patch_ic_object` | Initial conditions handler |
| `adam_patch_io_object.F90` | `adam_patch_io_object` | I/O options and residuals log |
| `adam_patch_time_object.F90` | `adam_patch_time_object` | Convergence and output frequency parameters |
| `adam_patch_common_library.F90` | `adam_patch_common_library` | Barrel re-export of all common modules |

## `patch_common_object`

`patch_common_object` is the base type for all PATCH backends. It aggregates the ADAM framework objects, the FLAIL linear algebra handler, and the PATCH-specific sub-objects.

### Data Members

| Member | Type | Description |
|--------|------|-------------|
| `mpih` | `mpih_object` | MPI handler |
| `adam` | `adam_object` | Core ADAM object (tree, grid, field, I/O) |
| `field` | `field_object*` | Pointer alias to `adam%field` |
| `grid` | `grid_object*` | Pointer alias to `adam%grid` |
| `amr` | `amr_object` | AMR marker handler |
| `ib` | `ib_object` | Immersed boundary handler |
| `flail` | `flail_object` | FLAIL linear algebra / smoothing handler |
| `io` | `patch_io_object` | I/O configuration |
| `ic` | `patch_ic_object` | Initial conditions |
| `bc` | `patch_bc_object` | Boundary conditions |
| `time` | `patch_time_object` | Convergence and output frequency |
| `q` | `real(R8P)(:,:,:,:,:)` | Scalar potential φ — shape `(nv, ni, nj, nk, nb)` |
| `r` | `real(R8P)(:,:,:,:,:)` | Source distribution f — shape `(nv, ni, nj, nk, nb)` |
| `q_name` | `character(3)(1)` | Output label for `q` (fixed to `'phi'`) |
| `ngc, ni, nj, nk, nb` | `integer(I4P)*` | Pointer aliases to grid/field dimension scalars |
| `blocks_number` | `integer(I4P)*` | Pointer alias to `field%blocks_number` |
| `nv` | `integer(I4P)*` | Pointer alias to `field%nv` (always 1 for PATCH) |

### Methods

| Procedure | Description |
|-----------|-------------|
| `allocate_common` | Allocates `r` with ghost-cell bounds `(nv, 1-ngc:ni+ngc, …, nb)` |
| `initialize_common` | Full initialisation sequence (see below) |

### `initialize_common` Sequence

`initialize_common` sets up the entire simulation state:

1. `mpih%initialize` — start MPI, query rank and memory
2. `io%initialize(filename)` — parse INI file
3. `bc%initialize` — load BC configuration
4. `adam%grid%initialize` — set domain bounds, resolution, BC types
5. `adam%compute_blocks_number` — estimate block count from available memory (`nv = 1`, `fields_number = 80`)
6. `adam%initialize` — build AMR tree (11 nodes), field (`nv=1`), maps; allocates `q`
7. `adam%refine_uniform` — apply initial uniform refinement levels from INI
8. `adam%prune` — prune tree according to `ijkl_prune` settings
9. `amr%initialize` — load AMR marker configuration
10. `time%initialize` — load convergence parameters
11. `ic%initialize` — load IC configuration
12. `ib%initialize` — load IB body geometry
13. `flail%initialize` — load linear algebra (smoothing method and iteration counts)
14. `allocate_common` — allocate `r`
15. `adam%io%initialize` — register `r` as `'rho'` field for output

## `patch_bc_object`

Manages Dirichlet boundary conditions on the six domain faces.

### Supported BC Types

| Identifier | INI value | Description |
|------------|-----------|-------------|
| `BC_DIRICHLET` | `dirichlet` | Fix φ to a constant value on the face |

### Data Members

| Member | Description |
|--------|-------------|
| `bc_type(6)` | BC type index for each face (1 = Dirichlet) |
| `q(:,:)` | Prescribed φ value for each face and variable |

### INI Configuration

```ini
[bc_x_min]
type = dirichlet
phi  = 0.0

[bc_x_max]
type = dirichlet
phi  = 0.0

[bc_y_min]
type = dirichlet
phi  = 0.0

[bc_y_max]
type = dirichlet
phi  = 0.0

[bc_z_min]
type = dirichlet
phi  = 0.0

[bc_z_max]
type = dirichlet
phi  = 0.0
```

### Methods

| Procedure | Description |
|-----------|-------------|
| `initialize` | Initialises MPI and calls `load_from_file` |
| `load_from_file` | Reads type and φ value for all six faces |
| `description` | Returns a human-readable summary |

## `patch_ic_object`

Sets the initial distribution of the source term `r` (forcing function f).

### Supported IC Types

| Identifier | `type` key | Description |
|------------|-----------|-------------|
| `IC_TYPE_GAUSS` | `gauss` | Gaussian centred at $(0.5,\,0.5,\,0.5)$ with variance $\sigma^2 = 0.01$ |

The Gaussian IC is defined pointwise as:

$$r_{i,j,k} = \exp\!\left(-\frac{(x_i - 0.5)^2 + (y_j - 0.5)^2 + (z_k - 0.5)^2}{0.01}\right)$$

### Data Members

| Member | Description |
|--------|-------------|
| `ic_type` | Selected IC type string |
| `amr_iterations` | Number of AMR refinement cycles to apply after imposing IC |
| `regions_number` | Number of IC regions |
| `r(1, regions_number)` | Source amplitude per region |
| `emin(3, regions_number)` | Bounding box minimum corner per region |
| `emax(3, regions_number)` | Bounding box maximum corner per region |

### INI Configuration

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

### Methods

| Procedure | Description |
|-----------|-------------|
| `initialize` | Initialises MPI and calls `load_from_file` |
| `load_from_file` | Reads type, region count, and per-region parameters |
| `set_initial_conditions` | Fills `r` array according to the selected IC type |
| `description` | Returns a human-readable summary |

## `patch_io_object`

Manages all I/O options: output frequency, restart, and the residuals log file.

### Data Members and INI Keys

INI section: `[IO]`

| INI key | Member | Default | Description |
|---------|--------|---------|-------------|
| `output_basename` | `output_basename` | — | Basename for output files |
| `it_save` | `it_save` | 100 | Main output iteration frequency |
| `restart` | `restart` | `.false.` | Enable restart from existing files |
| `restart_basename` | `restart_basename` | — | Basename for restart files |
| `restart_save` | `restart_save` | 100 | Restart checkpoint frequency |
| `residuals_save` | `residuals_save` | 10 | Residuals log write frequency |
| `save_memory_status` | `save_memory_status` | `.false.` | Log memory allocation status |

```ini
[IO]
output_basename   = results/patch
it_save           = 50
restart           = .false.
restart_basename  = restart/patch
restart_save      = 100
residuals_save    = 10
save_memory_status = .false.
```

### Residuals Log

The residuals file is written to `{output_basename}-residuals.dat` with columns:

```
VARIABLES="it" "time" "blocks_number" "rq1"
```

One line is written per `residuals_save` iterations. The `rq` value is the MPI-reduced L2 norm of `dq` normalised by the number of cells.

### Methods

| Procedure | Description |
|-----------|-------------|
| `initialize` | Opens and parses the INI file |
| `load_from_file` | Reads all IO options |
| `open_file_residuals` | Opens the residuals log file and writes the header |
| `close_file_residuals` | Closes the residuals log file |
| `save_residuals` | Appends one record to the residuals log |
| `description` | Returns a human-readable summary |

## `patch_time_object`

Holds iteration counters, convergence parameters, and output frequency for the solver loop. PATCH is an iterative (not time-marching) solver; the `time` and `CFL` fields are present for framework compatibility.

### Data Members and INI Keys

INI section: `[time]`

| INI key | Member | Default | Description |
|---------|--------|---------|-------------|
| `it_max` | `it_max` | -1 | Maximum iterations (-1 = no limit) |
| `time_max` | `time_max` | 1.0 | Maximum pseudo-time (used when `it_max ≤ 0`) |
| `CFL` | `CFL` | 0.3 | CFL number (framework compatibility) |

```ini
[time]
it_max   = 500
time_max = 1.0
CFL      = 0.3
```

### Termination Criteria

The solver exits if **either** condition is met:

| Condition | Description |
|-----------|-------------|
| `dq_max < tolerance` | Residual falls below FLAIL tolerance (primary criterion) |
| `it >= it_max` and `it_max > 0` | Iteration budget exhausted |
| `time >= time_max` and `it_max ≤ 0` | Pseudo-time budget exhausted |

### Methods

| Procedure | Description |
|-----------|-------------|
| `initialize` | Initialises MPI and calls `load_from_file` |
| `load_from_file` | Reads `it_max`, `time_max`, `CFL` |
| `is_to_save` | Returns `.true.` when data should be written at current iteration |
| `print_progress` | Prints iteration count, step, pseudo-time, and percentage to stdout |
| `description` | Returns a human-readable summary |

## `adam_patch_common_library`

Barrel re-export of all common modules. Backends use a single `use adam_patch_common_library` to access all common types and constants.

| Re-exported module | Key public symbols |
|--------------------|--------------------|
| `adam_patch_bc_object` | `patch_bc_object`, `BC_DIRICHLET` |
| `adam_patch_common_object` | `patch_common_object` |
| `adam_patch_ic_object` | `patch_ic_object`, `IC_TYPE_GAUSS` |
| `adam_patch_io_object` | `patch_io_object` |
| `adam_patch_time_object` | `patch_time_object` |

## License

PATCH is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
