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
- [Forest orchestration (multi-realm)](#forest-orchestration-multi-realm)

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

Re-exports every module listed above, **including all eight singleton modules**. Applications need only:

```fortran
use adam_common_library
```

For convenience, `adam_globals` bundles only the singleton variables (`mpih`, `grid`, `field`, `maps`, `tree`, `weno`, `ib`, `rk`) without dragging in the type-definition modules.

### Program-scope singletons

Each core object is exposed as a **module-level target variable** in its own singleton module. Any Fortran module that `use`-es the singleton module gains direct access to the object without receiving it as a dummy argument or embedding it in a derived type.

| Module | Variable | Type | Re-exported by |
|--------|----------|------|----------------|
| `adam_mpih_global` | `mpih` | `mpih_object` | `adam_common_library` |
| `adam_grid_global` | `grid` | `grid_object` | `adam_common_library` |
| `adam_field_global` | `field` | `field_object` | `adam_common_library` |
| `adam_maps_global` | `maps` | `maps_object` | `adam_common_library` |
| `adam_tree_global` | `tree` | `tree_object` | `adam_common_library` |
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

## Forest orchestration (multi-realm)

The **forest** machinery coordinates more than one independent simulation domain — a **realm** — through a synchronized time loop. The user-facing conceptual overview, the manifest schema, and the α/β cadence trade-offs are documented in [`/guide/forest.md`](/guide/forest); this section is the library-developer reference for the contract surface.

### `adam_realm_object` — the per-realm contract

`realm_object` is the abstract parent every solver application extends. It owns the conventional per-domain components (`tree`, `grid`, `field`, `maps`, `weno`, `rk`, `io`, `ib`, ...) and declares the **`_forest`-suffixed TBP family** the forest dispatches against. Default implementations error_stop; every app extension that participates in a multi-realm forest must override them.

| TBP | Role |
|---|---|
| `initialize_forest(filename, realms_number)` | Per-realm prologue: load INI, allocate state, place initial conditions, open IO. |
| `finalize_forest()`                          | Per-realm epilogue: close IO, release resources, MPI finalize. |
| `compute_local_dt_forest(dt_local)`          | Report this realm's CFL-limited dt; the forest min-reduces across realms. |
| `advance_one_step_forest(dt)`                | N=1 fast path: realm owns the entire step internally. |
| `open_step_forest(dt)`                       | N>1 path: per-step prologue (external fields, RK stage buffer init, time bookkeeping). |
| `begin_stage_forest(k, K_total, dt)`         | N>1 path: open integrator stage `k`; sets `self%stage_active = k`. |
| `end_stage_forest(k, K_total, dt, ...)`      | N>1 path: residuals + stage assignment for stage `k`. |
| `close_step_forest(dt)`                      | N>1 path: per-step epilogue (state assembly, BC, div-clean, time advance). |
| `fill_seam_from_peer_forest(peer, p_idx)`    | Receive-side roundtrip: copy peer's interior into self's seam ghosts. |
| `apply_reflux_to_stage_forest(stage, dt, flux_register)` | Apply Berger-Colella reflux at the realm's end-of-step (α.r1 gate). |
| `stages_per_step_forest()`                   | Report K (the realm's integrator stride; `rk%nrk` for SSP-RK). |
| `coupling_descriptor_forest(scheme_time, rk_scheme, nv)` | Report `(scheme_time, rk_scheme, nv)` for β admissibility (issue #18). |
| `is_done_forest(done)`                       | Report whether the simulation has reached its termination predicate. |
| `post_step_forest(dt, t, it, realm)`         | Per-step diagnostics / IO block (savers, residual prints). |

Each `_forest` TBP is `pass(self)`-dispatched. The forest receives `class(realm_object), intent(inout) :: realm(:)` as a parameter and never owns realm state; the realm array is allocated by the driver as a concrete app type (e.g. `type(prism_cpu_object), allocatable :: realm(:)`).

The N=1 fast path bypasses every per-stage TBP. Single-realm apps that do not participate in multi-realm forests may leave them at the default error_stop.

### `adam_forest_object` — the orchestrator

`forest_object` is behaviour-only: it owns no derived-type state (only an `n` integer realm count). All forest operations are TBP dispatches against the realm array passed in.

**Public methods**:

| Method | Purpose |
|---|---|
| `initialize(realm, filename)`                 | Shared-INI path: every realm reads the same INI. Used by legacy single-INI drivers (e.g. PRISM's no-manifest form). |
| `initialize_from_manifest(realm, manifest)`   | Manifest path: each realm reads its own INI from the parsed `forest_manifest_t`; topology is wired afterwards. |
| `simulate(realm, filename)`                   | Top-level entry: `initialize` → time loop → `finalize`. |
| `simulate_from_manifest(realm, manifest)`     | Top-level entry: manifest variant. |
| `compute_global_dt(realm, dt)`                | Min-reduce per-realm CFL across realms (intra-process) and ranks (`MPI_ALLREDUCE`). |
| `evolve_one_step(realm)`                      | One global timestep: Phase 0 → K-loop (Phases 1–3) → Phase 4 → Phase 5; per-seam cadence gating inside. |
| `post_step(realm)`                            | Per-step diagnostics block across realms. |
| `is_done(realm, done)`                        | AND-reduce per-realm termination across realms and ranks. |
| `finalize(realm)`                             | Sequence each realm's `finalize_forest`. |

**Private helpers**:

| Helper | Purpose |
|---|---|
| `populate_inter_realm_topology(realm, manifest)` | Translate manifest face-pairs into per-realm `maps%inter_realm_neighbors` and the derived `seam_local_*` arrays. |
| `check_beta_admissibility(realm, manifest)`      | Enforce β admissibility on every `stage_coincident` seam (PRD #18 M2). |
| `apply_reflux_corrections(realm, stage, dt)`     | Dispatch each realm's `apply_reflux_to_stage_forest` (realm-side body gates on end-of-step under α.r1). |

### Manifest data flow

```
forest.ini                          (text manifest)
  │
  ▼  adam_forest_manifest.F90 ─ read_forest_manifest
forest_manifest_t                   (transient struct: realm_ini paths + face_pairs)
  │
  ▼  forest_object%populate_inter_realm_topology
realm(is)%adam%maps%inter_realm_neighbors(:)              (per-realm symmetric topology)
  │
  ▼  build_inter_realm_ghost_cell_map + build_seam_local_map
realm(is)%adam%maps%seam_local_map_ghost_cell(:,:)        (sorted per-cell ghost map)
realm(is)%adam%maps%seam_local_peer_realm(:)              (one entry per distinct peer)
realm(is)%adam%maps%seam_local_peer_row_start(:)
realm(is)%adam%maps%seam_local_peer_row_count(:)
realm(is)%adam%maps%seam_local_cadence(:)                 (CADENCE_END_OF_STEP | CADENCE_STAGE_COINCIDENT)
realm(is)%adam%maps%seam_local_send_buf / recv_buf        (per-peer pack/unpack buffers, allocated once)
```

The orchestrator's per-step Phase 2 (β) and Phase 5 (α) loops consume the `seam_local_*` arrays directly; no manifest reads happen at step time.

### Per-seam coupling cadence

Two integer parameters on `adam_maps_object` tag each peer slot:

- `CADENCE_END_OF_STEP = 0` — α default: seam fill at the end of each global step, after `close_step_forest` (Phase 5).
- `CADENCE_STAGE_COINCIDENT = 1` — β opt-in: seam fill at every RK substage, before `end_stage_forest` (Phase 2 inside the K loop).

The orchestrator gates per-seam: Phase 2 fires only on `STAGE_COINCIDENT` seams; Phase 5 fires only on `END_OF_STEP` seams. A seam is filled exactly once per step under either cadence.

β admissibility — same `scheme_time`, `rk_scheme`, `nv`, and K between both endpoint realms — is verified at `initialize_from_manifest` time via `check_beta_admissibility`. Mismatch → immediate `error_stop`; no silent downgrade.

### Load-bearing invariants

- **Phase 2 → Phase 3 ordering (β)**: Phase 2 must complete on ALL realms before Phase 3 starts on ANY realm; otherwise the read-after-overwrite race between `fill_seam_from_peer_forest` writes and `compute_residuals` reads returns. The serial inner loops within a rank give this for free under the Phase-A replicated-forest layout.
- **`stage_active` discipline**: `begin_stage_forest(k)` MUST set `self%stage_active = k`; `close_step_forest(dt)` MUST clear it back to 0. The `fill_seam_from_peer_forest` buffer-selection logic reads `self%stage_active` and `peer%stage_active` to pick between `q` (committed) and `rk%q_rk(:,...,stage_active)`. Mid-step writes to `q` outside this discipline corrupt peer reads.
- **α.r1 reflux gate**: `apply_reflux_to_stage_forest` and `accumulate_seam_fluxes_fv` MUST gate on `stage == rk%nrk` (the realm's own final substage). `flux_register`'s third axis is collapsed to 1 under α.r1; per-stage RK-weighted reflux (Wang 2018) is deferred. Both α and β share this gate — β does not undo α.r1.
- **N=1 fast path**: when `size(realm) == 1`, `evolve_one_step` routes through `advance_one_step_forest` directly with zero exposure to the per-stage TBPs. App backends MAY leave the per-stage TBPs at the parent's error_stop default for N=1-only use.
- **No singletons embedded on realms**: the program-scope singletons (`mpih`, `grid`, `field`, ...) are aliased per-step from each active realm's components; they are not value members on the realm type. See "Program-scope singletons" above.

### Tracked through

| PRD | Scope | Status |
|---|---|---|
| [#10](https://github.com/szaghi/adam/issues/10)  | Phase D inception — forest orchestrator design          | shipped |
| [#13](https://github.com/szaghi/adam/issues/13)  | Coarse-fine interface machinery (A/B/C)                 | A shipped; B & C deferred |
| [#16](https://github.com/szaghi/adam/issues/16)  | α end-of-step barrier + asymmetric K                    | shipped |
| [#18](https://github.com/szaghi/adam/issues/18)  | β stage-coincident recovery + cross-config oracle       | shipped |

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
| `adam_realm_object` | Forest | `adam_realm_object.F90` |
| `adam_forest_object` | Forest | `adam_forest_object.F90` |
| `adam_forest_manifest` | Forest | `adam_forest_manifest.F90` |
| `adam_flux_register_object` | Forest | `adam_flux_register_object.F90` |
| `adam_common_library` | Entry point | `adam_common_library.F90` |
| `adam_mpih_global` | Singleton | `adam_mpih_global.F90` |
| `adam_grid_global` | Singleton | `adam_grid_global.F90` |
| `adam_field_global` | Singleton | `adam_field_global.F90` |
| `adam_maps_global` | Singleton | `adam_maps_global.F90` |
| `adam_tree_global` | Singleton | `adam_tree_global.F90` |
| `adam_weno_global` | Singleton | `adam_weno_global.F90` |
| `adam_ib_global` | Singleton | `adam_ib_global.F90` |
| `adam_rk_global` | Singleton | `adam_rk_global.F90` |
| `adam_globals` | Singleton aggregator | `adam_globals.F90` |

---

## Copyrights

ADAM is released under the [GNU Lesser General Public License v3.0](/license) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
