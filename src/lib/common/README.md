# ADAM Common Library — `src/lib/common`

The common library is the **portable, physics-agnostic core** of the ADAM SDK. Every GPU backend (`fnl`, `nvf`, `gmp`) extends these objects rather than replacing them; every solver application (`nasto`, `prism`, `chase`, …) composes them. No GPU-specific code lives here — the entire directory compiles with any standard Fortran 2003+ compiler.

The aggregate entry point `adam_common_library` re-exports all modules below; a single `use adam_common_library` statement in application code exposes the complete API.

---

## Contents

- [AMR data structures](#amr-data-structures)
- [Grid and field](#grid-and-field)
- [Communication and maps](#communication-and-maps)
- [Numerical methods — spatial](#numerical-methods--spatial)
- [Numerical methods — temporal](#numerical-methods--temporal)
- [Physics support](#physics-support)
- [Infrastructure](#infrastructure)

---

## AMR data structures

### `adam_tree_object` — octree/quadtree hierarchy

The tree is a **CPU-resident hash-map** keyed on 64-bit Morton indices that linearise the four octree coordinates `(level, bx, by, bz)` into a single integer. It provides O(1) block insertion, deletion, and neighbour lookup during AMR update phases. Morton ordering clusters geometrically adjacent blocks in index space, minimising ghost-cell communication volume across MPI ranks.

The refinement ratio is selectable: 4 (quadtree) or 8 (octree). Full 26-neighbourhood support covers all face, edge, and corner adjacencies needed by high-order stencils.

The tree is never accessed by numerical kernels during time integration — it is consulted only during AMR update cycles (marking → refinement/coarsening → load rebalancing).

### `adam_tree_node_object` — individual tree node

Represents a single block in the octree. Key fields:

| Field | Type | Purpose |
|-------|------|---------|
| `code` | `I8P` | Morton key |
| `refinement_needed` | `I4P` | AMR marker (`TO_BE_REFINED`, `TO_BE_DEREFINED`, `TO_NOT_TOUCH`) |
| `myrank` | `I4P` | Owning MPI rank |
| `block_index` | `I8P` | Compact index into the field array |
| `neighbor(26)` | `tree_node_neighbor_object` | Face/edge/corner adjacency data |
| `surface_stl_distance` | `R8P` | Signed distance to nearest IB surface |

### `adam_tree_bucket_object` — hash-map bucket

Double-linked list container underpinning the tree's hash-map. Stores nodes by Morton key and provides O(1) `add_node`, `remove_node`, and `has_code` operations. Supports iterator-based traversal via the `traverse` method and a user-supplied callback conforming to `iterator_interface`.

---

## Grid and field

### `adam_grid_object` — block-structured Cartesian grid

Manages the geometric description of a single structured block and the global domain layout. Provides cell/node coordinate computation, mesh metrics, ghost-cell layout, and boundary condition assignment.

Key configuration (`[grid]` INI section):

| Parameter | Meaning |
|-----------|---------|
| `ni`, `nj`, `nk` | Cells per block per direction |
| `ngc` | Ghost-cell width (typically 3 for WENO5, 4 for WENO7) |
| `emin_*`, `emax_*` | Domain extents |
| `bc_type(6)` | Boundary condition type on each of the 6 domain faces |

Notable methods: `block_emin/emax`, `cell_xyz`, `node_xyz`, `compute_metrics`, `compute_weight_neighbor`, `load_from_ini_file`.

### `adam_field_object` — 5D field array

Dense contiguous array holding all simulation variables for all blocks owned by an MPI rank:

```
q(nv, ni, nj, nk, nb)
```

| Dimension | Meaning |
|-----------|---------|
| `nv` | Physical variables (density, momenta, energy, …) |
| `ni, nj, nk` | Cell indices within a block, including ghost cells |
| `nb` | Compact sequential block index (1 … current block count) |

Fortran column-major order keeps `nv` contiguous in memory, enabling coalesced reads across CUDA threads or OpenACC/OpenMP SIMD lanes. The array is allocated once at simulation start to the maximum block count and reused across all Runge-Kutta stages — no dynamic allocation during time integration.

---

## Communication and maps

### `adam_maps_object` — index translation and MPI communication patterns

Bridges the tree (Morton-key addressed) and the field (compact-index addressed). Rebuilt only when the AMR topology changes; invisible to numerical kernels during time integration.

**Local maps** (updated every AMR cycle):

| Map | Purpose |
|-----|---------|
| `local_map` | Block index ↔ Morton key translation (`b2m` / `m2b`) |
| `local_map_ghost` | Block-to-ghost-block adjacency for exchange |
| `local_map_ghost_cell` | Cell-level ghost stencil offsets |
| `local_map_bc_face/edge/corner/crown` | Boundary condition application stencils |

**MPI communication lists** (per-rank send/receive schedules):

- `comm_map_send/recv` — block-level exchange for AMR rebalancing
- `comm_map_send/recv_ghost` — ghost-cell exchange buffers before stencil sweeps
- `send_buffer_ghost`, `recv_buffer_ghost` — pre-allocated MPI staging buffers

Key methods: `initialize`, `make_comm_local_maps`, `blocks_reorder`.

---

## Numerical methods — spatial

### `adam_weno_object` — WENO reconstruction

Weighted Essentially Non-Oscillatory finite difference reconstruction. Upwind and centered variants are available at orders 1, 3, 5, 7, and 9 (upwind) and 2+ (centered). The stencil half-width `S_MAX` governs the ghost-cell requirement.

Scheme selectors:

| Constant | Scheme |
|----------|--------|
| `WENO_U_1` … `WENO_U_9` | Upwind WENO, orders 1/3/5/7/9 |
| `WENO_C_2` … | Centered WENO |

Reference: Chi-Wang Shu, NASA/CR-97-206253 (1997).

### `adam_fdv_operators_library` — finite difference/volume operators

Spatial derivative operators at configurable order for use in application physics layers. All operators are available in three flavours:

| Flavour | Description |
|---------|-------------|
| Finite difference, centered | Standard FD stencils on cell centres |
| Finite volume, centered | Cell-average-based FV operators |
| Finite volume, upwind | WENO-reconstructed upwind FV fluxes |

Operators: gradient, divergence, curl, Laplacian, and scalar derivatives of orders 1–6. Abstract `compute_*_fdv_interface` procedure pointers allow backend dispatch without branching.

### `adam_riemann_euler_library` — Euler equation Riemann solvers

Convective interface flux computation for the compressible Euler equations. All solvers accept left/right primitive states and return the numerical flux; the optional `lmax` output gives the maximum wave speed for CFL estimation.

| Procedure | Solver |
|-----------|--------|
| `compute_riemann_euler_exact` | Exact solver (Newton iteration on Rankine-Hugoniot) |
| `compute_riemann_euler_hllc` | HLLC (Harten–Lax–van Leer–Contact) |
| `compute_riemann_euler_hllc_lm` | HLLC with low-Mach correction |
| `compute_riemann_euler_hllem` | HLLEM (entropy-fix variant) |
| `compute_riemann_euler_llf` | Local Lax–Friedrichs |
| `compute_riemann_euler_ts` | Trivial (first-order upwind) |

---

## Numerical methods — temporal

### `adam_rk_object` — explicit Runge-Kutta integration

| Constant | Scheme | Stages | Order |
|----------|--------|--------|-------|
| `RK_1` | Forward Euler | 1 | 1 |
| `RK_2` | Heun | 2 | 2 |
| `RK_3` | TVD-RK3 | 3 | 3 |
| `RK_SSP_22` | SSP(2,2) | 2 | 2 |
| `RK_SSP_33` | SSP(3,3) | 3 | 3 |
| `RK_SSP_54` | SSP(5,4) | 5 | 4 |
| `RK_YOSHIDA` | Yoshida symplectic | 3 | 4 |

Stage storage `q_rk(nv, ni, nj, nk, nb, nrk)` is pre-allocated; no allocation occurs during time-stepping. Configuration section: `[runge_kutta]`.

### `adam_leapfrog_object` — leapfrog with RAW filter

Second-order leapfrog scheme with optional Robert–Asselin–Williams noise filter:

```
U^{n+2} = U^{n} + 2·Δt·R(t^{n+1}, U^{n+1})
Δ = ν/2·(U^{n} - 2·U^{n+1} + U^{n+2})    ! RAW correction
```

Filter coefficient `nu` (default 0.01) and Williams parameter `alpha` (default 0.53) are configurable. Configuration section: `[leapfrog]`.

### `adam_blanes_moan_object` — Blanes–Moan splitting integrator

Geometric integrator for separable Hamiltonian systems `H = A + B`. Two schemes available:

| Constant | Order | Stages |
|----------|-------|--------|
| `BM_SCHEME4` | 4 | 4 |
| `BM_SCHEME6` | 6 | 6 |

Coefficients `a(:)` and `b(:)` implement the alternating A-B operator splitting. Configuration section: `[blanes_moan]`.

### `adam_cfm_object` — Commutator-Free Magnus integrators

Exponential integrators for non-autonomous ODEs of the form `dU/dt = A(t)·U`. Computes `U^{n+1} = exp(α·R)·U^n` without evaluating the matrix exponential directly (commutator-free formulation).

| Constant | Order | Stages | Exponentials |
|----------|-------|--------|--------------|
| `CFM_4` | 4 | — | — |
| `CFM_6` | 6 | — | — |
| `CFM_8` | 8 | — | — |

Stage residuals `dq` and field buffer `q` are pre-allocated. Configuration section: `[commutator_free_magnus]`.

---

## Physics support

### `adam_eos_ic_object` — ideal compressible gas EOS

Thermodynamic model for a calorically perfect ideal gas. Stores primary properties and derives all secondary quantities on initialisation.

| Property | Symbol | Description |
|----------|--------|-------------|
| `cp`, `cv` | — | Specific heats at constant pressure/volume |
| `g` | γ | Heat capacity ratio |
| `gm1`, `gp1` | γ−1, γ+1 | Pre-computed combinations |
| `delta`, `eta` | (γ−1)/2, 2γ/(γ−1) | Isentropic flow coefficients |
| `mu` | μ | Dynamic viscosity |
| `kd` | k | Thermal conductivity |

Key procedures: `density`, `pressure`, `temperature`, `internal_energy`, `total_energy`, `speed_of_sound`, `primitive2conservative`, `conservative2primitive`. Configuration section: `[physics_specie_*]`.

### `adam_ib_object` — immersed boundary method

Manages solid bodies immersed in the Cartesian grid. Computes the signed-distance field `phi(nv, ni, nj, nk, nb)` via an eikonal solver and applies wall boundary conditions at solid surfaces.

Supported solid definitions:

| Type | Constant | Description |
|------|----------|-------------|
| Sphere / circle | `IB_ANALYTICAL_SPHERE` | Defined by centre, radius, axis |
| Circle (XY plane) | `IB_ANALYTICAL_CIRCLE` | 2D equivalent |
| Rectangle | `IB_ANALYTICAL_RECTANGLE` | Defined by centre, edge lengths, axis |
| Surface mesh | — | OFF-format triangulated surface |

Wall boundary condition types:

| Constant | Value | Type |
|----------|-------|------|
| `BCS_VISCOUS` | 1 | Viscous (no-slip) wall |
| `BCS_EULER` | 2 | Inviscid (slip) wall |

Configuration section: `[solid]`. Number of eikonal integration steps configurable via `n_eikonal`.

### `adam_amr_object` — AMR refinement markers

Configures the criteria that drive block refinement and coarsening. Multiple markers can be combined; each marker specifies an independent criterion applied at every AMR update cycle.

Marker modes:

| Constant | Value | Criterion |
|----------|-------|-----------|
| `AMR_GEO` | 1 | Geometry-based (distance to IB surface) |
| `AMR_GRAD` | 2 | Field gradient magnitude |

Each marker carries threshold parameters `delta_fine` / `delta_coarse` (cell-spacing ratios), a tolerance `tol`, and a reference to the field variable and direction. AMR is triggered every `frequency` time steps and iterated `iters` times per cycle to enforce 2:1 level balance. Configuration section: `[amr]`.

### `adam_flail_object` — linear algebra interface

Iterative solver and smoother infrastructure for implicit operators (e.g., elliptic sub-problems). Provides a common interface over several smoothing strategies:

| Constant | Method |
|----------|--------|
| `SMOOTHING_MULTIGRID` | Geometric multigrid V-cycle |
| `SMOOTHING_GAUSS_SEIDEL` | Gauss–Seidel iteration |
| `SMOOTHING_SOR` | Successive Over-Relaxation |
| `SMOOTHING_SOR_OMP` | Parallel SOR (OpenMP) |

Convergence is controlled by `tolerance` (default 10⁻⁶) and `iterations`. Configuration section: `[linear-algebra]`.

---

## Infrastructure

### `adam_mpih_object` — MPI wrapper

Thin convenience layer over MPI. Provides rank/size query, tic-toc performance timing, rank-prefixed console output, non-blocking request bookkeeping (`req_send_recv`), and `abort`/`error_stop` helpers. All communication in the maps and field objects routes through this object.

### `adam_io_object` — parallel I/O

Manages all simulation data output and restart files. Holds pointers to the primary field and to up to 30 auxiliary arrays (vectors and scalars at multiple precisions: `R8P`, `R4P`, `I8P`, `I4P`, `I2P`, `I1P`). Supported output formats: HDF5 (parallel, MPI-IO), VTK, XH5F, and MAT (MATLAB).

### `adam_slices_object` — 2D diagnostic slices

Extracts planar 2D slices from the 3D field at configurable spatial bounds and resolution (`emin`, `emax`, `nijk`). Save frequency is controlled by `n_save`. Output format: MAT. Configuration section: `[slices]`.

### `adam_adam_object` — top-level orchestrator

Composer that binds all SDK objects into a single `adam_object` and drives the simulation lifecycle:

- **Initialisation**: `initialize` — sets up grid, tree, maps, field, and I/O from INI file
- **AMR cycle**: `adapt`, `amr_update`, `refine_uniform` — marking, refinement/coarsening, load rebalancing
- **Communication**: `make_comm_local_maps_ghost_bc`, `mpi_redistribute`, `mpi_gather_refinement_needed`, `blocks_reorder`
- **I/O**: `save_hdf5`, `save_restart_files`, `save_vtk`, `save_xh5f`, `save_slice`, `load_restart_files`
- **Geometry**: `interpolate_at_point`, `prune`

### `adam_parameters` — global constants

Named integer constants shared across all modules: boundary condition type flags, AMR marker values (`TO_BE_REFINED`, `TO_BE_DEREFINED`, `TO_NOT_TOUCH`), face enumeration codes (`FEC_1_6_ARRAY`), and coordinate-direction mappings (`FEC_TO_DELTA`, `DELTA_TO_FEC`).

### `adam_common_library` — aggregate entry point

Re-exports every module listed above, **including all seven singleton modules**. Applications need only:

```fortran
use adam_common_library
```

### Program-scope singletons

Each core object is exposed as a **module-level target variable** in its own singleton module. Any Fortran module that `use`-es the singleton module gains direct access to the object without receiving it as a dummy argument or embedding it in a derived type.

| Module | Variable | Type | Re-exported by |
|--------|----------|------|----------------|
| `adam_mpih_global` | `mpih` | `mpih_object` | `adam_common_library` |
| `adam_grid_global` | `grid` | `grid_object` | `adam_common_library` |
| `adam_field_global` | `field` | `field_object` | `adam_common_library` |
| `adam_maps_global` | `maps` | `maps_object` | `adam_common_library` |
| `adam_weno_global` | `weno` | `weno_object` | `adam_common_library` |
| `adam_ib_global` | `ib` | `ib_object` | `adam_common_library` |
| `adam_rk_global` | `rk` | `rk_object` | `adam_common_library` |

**Access pattern:**

```fortran
use :: adam_grid_global,  only: grid
use :: adam_field_global, only: field

associate(ni=>grid%ni, ngc=>grid%ngc, nb=>field%nb)
  ! numerical kernel — no argument passing needed
endassociate
```

Singletons are never passed as dummy arguments and never embedded as members of other derived types.

---

## Module summary

| Module | Category | File |
|--------|----------|------|
| `adam_parameters` | Infrastructure | `adam_parameters.f90` |
| `adam_tree_node_object` | AMR | `adam_tree_node_object.f90` |
| `adam_tree_bucket_object` | AMR | `adam_tree_bucket_object.f90` |
| `adam_tree_object` | AMR | `adam_tree_object.F90` |
| `adam_grid_object` | Grid/Field | `adam_grid_object.F90` |
| `adam_field_object` | Grid/Field | `adam_field_object.F90` |
| `adam_maps_object` | Communication | `adam_maps_object.F90` |
| `adam_weno_object` | Numerics — spatial | `adam_weno_object.F90` |
| `adam_fdv_operators_library` | Numerics — spatial | `adam_fdv_operators_library.F90` |
| `adam_riemann_euler_library` | Numerics — spatial | `adam_riemann_euler_library.F90` |
| `adam_rk_object` | Numerics — temporal | `adam_rk_object.F90` |
| `adam_leapfrog_object` | Numerics — temporal | `adam_leapfrog_object.F90` |
| `adam_blanes_moan_object` | Numerics — temporal | `adam_blanes_moan_object.F90` |
| `adam_cfm_object` | Numerics — temporal | `adam_cfm_object.F90` |
| `adam_eos_ic_object` | Physics | `adam_eos_ic_object.F90` |
| `adam_ib_object` | Physics | `adam_ib_object.F90` |
| `adam_amr_object` | Physics | `adam_amr_object.F90` |
| `adam_flail_object` | Physics | `adam_flail_object.F90` |
| `adam_mpih_object` | Infrastructure | `adam_mpih_object.F90` |
| `adam_io_object` | Infrastructure | `adam_io_object.F90` |
| `adam_slices_object` | Infrastructure | `adam_slices_object.F90` |
| `adam_adam_object` | Infrastructure | `adam_adam_object.F90` |
| `adam_common_library` | Entry point | `adam_common_library.F90` |
| `adam_mpih_global` | Singleton | `adam_mpih_global.F90` |
| `adam_grid_global` | Singleton | `adam_grid_global.F90` |
| `adam_field_global` | Singleton | `adam_field_global.F90` |
| `adam_maps_global` | Singleton | `adam_maps_global.F90` |
| `adam_weno_global` | Singleton | `adam_weno_global.F90` |
| `adam_ib_global` | Singleton | `adam_ib_global.F90` |
| `adam_rk_global` | Singleton | `adam_rk_global.F90` |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
