# ADAM Architecture Reference

Comprehensive reference for class hierarchy, singletons, key types, library structure, and design patterns.

---

## Inheritance Hierarchy

```
src/lib/common/  (CPU base layer)
├── mpih_object
├── grid_object
├── field_object
├── maps_object
├── tree_object  ←─ tree_node_object, tree_bucket_object
├── weno_object
├── rk_object
├── ib_object    ←─ analytical_sphere_object, analytical_rectangle_object
├── io_object
├── amr_object
├── adam_object  (owns: grid, tree, field, maps as VALUE components — issue #10 Step 1; b40dc451)
├── realm_object (formerly equation_object — issue #10 Step 6 rename, 356cac46; owns: adam, io, amr, ib, rk, weno, blanesmoan, cfm, leapfrog, flail, slices)
└── forest_object (behavior-only orchestrator over realm(:) — issue #10 Step 6 Phase C, 2d74fb6d + a294681a)

src/lib/fnl/  (OpenACC GPU backend)
├── mpih_fnl_object   extends mpih_object
├── field_fnl_object  (composite, wraps field_object data on GPU)
├── ib_fnl_object     (composite, GPU phi and BC vars)
├── weno_fnl_object   (composite, GPU WENO coefficients)
└── rk_fnl_object     (composite, GPU RK stages)

src/app/nasto/  (Navier-Stokes solver)
├── nasto_common_object   (standalone, does NOT extend realm_object)
│   ├── nasto_cpu_object  extends nasto_common_object
│   ├── nasto_nvf_object  extends nasto_common_object  [CUDA Fortran]
│   ├── nasto_gmp_object  extends nasto_common_object  [OpenMP offload]
│   └── nasto_fnl_object  extends nasto_common_object  [OpenACC]

src/app/prism/  (Maxwell equations solver)
├── prism_common_object   extends realm_object
│   ├── prism_cpu_object  extends prism_common_object
│   └── prism_fnl_object  extends prism_common_object  [OpenACC]

src/app/chase/
└── chase_common_object
    └── chase_cpu_object  extends chase_common_object

src/app/patch/
└── patch_common_object
    └── patch_cpu_object  extends patch_common_object
```

---

## Program-Scope Singletons

Singletons are 13-line modules exposing a single `target` module variable. They are accessed via `use`, never passed as dummy arguments, never embedded in derived types.

### CPU Singletons (`src/lib/common/`)

| Module | Variable | Type | Purpose |
|--------|----------|------|---------|
| `adam_mpih_global` | `mpih` | `mpih_object` | MPI handler (rank, comm, timing) |
| `adam_grid_global` | `grid` | `grid_object` | Domain/discretization parameters |
| `adam_field_global` | `field` | `field_object` | Block field data and mesh arrays |
| `adam_maps_global` | `maps` | `maps_object` | Ghost cell communication maps |
| `adam_tree_global` | `tree` | `tree_object` | AMR Morton-ordered octree |
| `adam_weno_global` | `weno` | `weno_object` | WENO scheme coefficients |
| `adam_ib_global` | `ib` | `ib_object` | Immersed boundary / eikonal solver |
| `adam_rk_global` | `rk` | `rk_object` | Runge-Kutta integrator + stage arrays |

All 8 are re-exported via `adam_common_library` (and bundled by the convenience module `adam_globals`).

### FNL GPU Singletons (`src/lib/fnl/`)

| Module | Variable | Type | Purpose |
|--------|----------|------|---------|
| `adam_fnl_mpih_global` | `mpih_fnl` | `mpih_fnl_object` | GPU-aware MPI + device init |
| `adam_fnl_field_global` | `field_fnl` | `field_fnl_object` | GPU field arrays (coords, dxyz) |
| `adam_fnl_ib_global` | `ib_fnl` | `ib_fnl_object` | GPU distance function phi |
| `adam_fnl_rk_global` | `rk_fnl` | `rk_fnl_object` | GPU RK stage arrays |
| `adam_fnl_weno_global` | `weno_fnl` | `weno_fnl_object` | GPU WENO coefficients |

All 5 are re-exported via `adam_fnl_library`.

### PRISM FNL Singletons (`src/app/prism/fnl/`)

| Module | Variable | Type | Purpose |
|--------|----------|------|---------|
| `adam_prism_fnl_coil_global` | `coil_fnl` | `prism_fnl_coil_object` | GPU coil current arrays |
| `adam_prism_fnl_fwlayer_global` | `fwlayer_fnl` | `prism_fnl_fwlayer_object` | GPU Faraday-wall layer arrays |

Both re-exported via `adam_prism_fnl_library`.

### Singleton Template

```fortran
module adam_<name>_global
use :: adam_<name>_object
implicit none
public
type(<name>_object), target :: <name>
endmodule adam_<name>_global
```

### FNL Initialization Order

CPU value singletons (`ib`, `rk`, `weno`) must be copied from the solver state **before** calling FNL `%initialize()`, because FNL objects read them at init time:

```fortran
! In nasto_fnl_object%initialize / prism_fnl_object%initialize_prism:
ib   = self%ib    ! copy CPU ib_object into ib singleton
rk   = self%rk    ! copy CPU rk_object into rk singleton
weno = self%weno  ! copy CPU weno_object into weno singleton
call mpih_fnl%initialize(do_mpi_init=.true., do_device_init=.true.)
call field_fnl%initialize(...)
call ib_fnl%initialize()
call rk_fnl%initialize()
call weno_fnl%initialize()
```

---

## Key Type Members

### `adam_object` (`src/lib/common/adam_adam_object.F90`)

After issue #10 Step 1 (commit `b40dc451`), `adam_object` owns `grid`, `tree`, `field`, `maps` as **value** components (no longer method-only). The program-scope singletons still exist but are now **pointer shims** aliasing into `adam_singleton%grid`, etc. (see "Program-Scope Singletons" below).

```fortran
type :: adam_object
   type(grid_object)  :: grid
   type(tree_object)  :: tree
   type(field_object) :: field
   type(maps_object)  :: maps
end type
```

### `realm_object` (`src/lib/common/adam_realm_object.F90`)

Renamed from `equation_object` in issue #10 Step 6 (commit `356cac46`). A **realm** has both *laws* (PDE, numerics, integrator, BCs) and *borders* (grid, tree, field, maps) — the name reflects the synthesis. Aggregator base class for per-realm solver state:

```fortran
type :: realm_object
   type(io_object)         :: io
   type(adam_object)       :: adam     ! owns grid/tree/field/maps as value components (Step 1)
   type(amr_object)        :: amr
   type(weno_object)       :: weno     ! issue #10 Step 2: moved here from singleton
   type(ib_object)         :: ib       ! issue #10 Step 2
   type(rk_object)         :: rk       ! issue #10 Step 2
   type(slices_object)     :: slices
   type(blanesmoan_object) :: blanesmoan
   type(cfm_object)        :: cfm
   type(leapfrog_object)   :: leapfrog
   type(flail_object)      :: flail
   ! FDV scheme metadata (used by backend dispatch):
   character(:), allocatable :: fdv_scheme
   integer(I4P)              :: fdv_order, fdv_half_stencil, fdv_half_stencils(6)
   ! Forest-position context (set by forest%initialize_from_manifest before the
   ! realm's initialize_forest fires; consumed by realm-side forest TBPs):
   integer(I4P) :: realm_index   = 0_I4P
   integer(I4P) :: stage_active  = 0_I4P  ! 0 between steps; k while RK substage k is open
   ! Scalar replica pointers (point into grid/field components):
   integer(I4P), pointer :: ngc, ni, nj, nk, nb, blocks_number, nv
   ! FDV operator procedure pointers (set at init by backend):
   procedure(compute_curl_interface),        pointer :: compute_curl
   procedure(compute_gradient_interface),    pointer :: compute_gradient
   procedure(compute_divergence_interface),  pointer :: compute_divergence
   procedure(compute_laplacian_interface),   pointer :: compute_laplacian
   procedure(compute_derivative1_interface), pointer :: compute_derivative1
   procedure(compute_derivative2_interface), pointer :: compute_derivative2
   procedure(compute_derivative4_interface), pointer :: compute_derivative4
contains
   ! Forest contract — `_forest`-suffixed TBPs the orchestrator dispatches.
   ! Default bodies are error_stop; every multi-realm participant overrides.
   ! Per-step contract (N=1 fast path uses advance_one_step_forest only):
   procedure, pass(self) :: initialize_forest          ! per-realm init (manifest or shared INI)
   procedure, pass(self) :: finalize_forest            ! per-realm shutdown
   procedure, pass(self) :: compute_local_dt_forest    ! report CFL dt → forest min-reduce
   procedure, pass(self) :: advance_one_step_forest    ! N=1 fast path: realm owns whole step
   procedure, pass(self) :: open_step_forest           ! N>1: per-step prologue
   procedure, pass(self) :: begin_stage_forest         ! N>1: open RK substage k; sets stage_active=k
   procedure, pass(self) :: end_stage_forest           ! N>1: residuals + assign for substage k
   procedure, pass(self) :: close_step_forest          ! N>1: per-step epilogue (clears stage_active)
   procedure, pass(self) :: fill_seam_from_peer_forest ! receive-side: copy peer interior → self ghosts
   procedure, pass(self) :: apply_reflux_to_stage_forest ! Berger-Colella reflux (α.r1 gates on stage==nrk)
   procedure, pass(self) :: stages_per_step_forest     ! report K (rk%nrk for SSP-RK)
   procedure, pass(self) :: coupling_descriptor_forest ! report (scheme_time, rk_scheme, nv) for β admissibility
   procedure, pass(self) :: is_done_forest             ! AND-reduced termination predicate
   procedure, pass(self) :: post_step_forest           ! per-step diagnostics / IO block
end type
```

### `forest_object` (`src/lib/common/adam_forest_object.F90`)

Added in issue #10 Step 6 Phase C (commits `2d74fb6d` + `a294681a`); extended through PRDs [#13](https://github.com/szaghi/adam/issues/13) (interface machinery A), [#16](https://github.com/szaghi/adam/issues/16) (α end-of-step barrier + asymmetric K), and [#18](https://github.com/szaghi/adam/issues/18) (β stage-coincident seam coupling). **Behavior-only orchestrator**: owns no derived-type state, receives the realm array as `class(realm_object), intent(inout) :: realm(:)` and orchestrates inter-realm operations via the `_forest`-suffixed TBPs above.

```fortran
type :: forest_object
   integer(I4P)               :: n = 0_I4P                 ! realm count
   type(flux_register_object) :: flux_register             ! Berger-Colella reflux register (α.r1: third axis = 1)
   ! NO derived-type-pointer components, ever (R1 of issue #10).
contains
   procedure, pass(self) :: initialize                 ! shared-INI: every realm reads the same file
   procedure, pass(self) :: initialize_from_manifest   ! manifest path: each realm reads its own INI
   procedure, pass(self) :: simulate                   ! entry: initialize → loop → finalize
   procedure, pass(self) :: simulate_from_manifest     ! entry: manifest variant
   procedure, pass(self) :: compute_global_dt          ! min-reduce realm(:)%compute_local_dt_forest + MPI_ALLREDUCE
   procedure, pass(self) :: evolve_one_step            ! one global timestep (Phase 0..5; per-seam α/β gating)
   procedure, pass(self) :: post_step                  ! iterate realm(:)%post_step_forest
   procedure, pass(self) :: is_done                    ! AND-reduce realm(:)%is_done_forest + MPI_ALLREDUCE
   procedure, pass(self) :: finalize                   ! sequence realm(:)%finalize_forest
   procedure, pass(self), private :: populate_inter_realm_topology  ! manifest → per-realm seam_local_* arrays
   procedure, pass(self), private :: check_beta_admissibility       ! enforce β contract on stage_coincident seams
   procedure, pass(self), private :: apply_reflux_corrections       ! dispatch realm(:)%apply_reflux_to_stage_forest
end type
```

`evolve_one_step` phase outline (N>1 path):

```
Phase 0 — open_step_forest      (per-realm prologue)
For k = 1..K_max:
  Phase 1 — begin_stage_forest  (per-realm, K-gated by k ≤ K_realm(is))
  Phase 2 — fill_seam_from_peer_forest  (per-seam, fires only for β seams)
  Phase 3 — end_stage_forest    (per-realm, K-gated)
            + flux_register%reduce_fine_sums + apply_reflux_corrections
Phase 4 — close_step_forest     (per-realm epilogue)
Phase 5 — fill_seam_from_peer_forest  (per-seam, fires only for α seams)
```

Each inter-realm seam carries a `coupling_cadence` (cached in `realm(is)%adam%maps%seam_local_cadence(p)`):

- `CADENCE_END_OF_STEP` (α default; AMReX `FillCoarsePatch` convention): mid-step peer ghosts are stale-by-one-step; one synchronized exchange at Phase 5. Admits asymmetric K. Cost: first-order seam coupling in time.
- `CADENCE_STAGE_COINCIDENT` (β; opt-in per seam): per-substage exchange at Phase 2. Admissibility (enforced by `check_beta_admissibility` at init): same `scheme_time`, `rk_scheme`, `physics%nv`, K on both endpoints. Mismatch → `error_stop`, no silent downgrade. Recovers single-realm bit-equivalence when admissible.

Driver pattern (manifest path, e.g. `src/app/prism/cpu/adam_prism_cpu.F90`):

```fortran
type(prism_cpu_object), allocatable :: realm(:)
type(forest_object)                  :: forest
type(forest_manifest_t)              :: manifest

if (is_forest_manifest(input_file_name)) then
   call read_forest_manifest(filename=input_file_name, manifest=manifest)
   allocate(realm(manifest%realms_number))
   call forest%simulate_from_manifest(realm=realm, manifest=manifest)
else
   allocate(realm(1))                                    ! N=1 fast path
   call forest%simulate(realm=realm, filename=input_file_name)
endif
```

**nvfortran 26.1 polymorphic-array workaround** (commit `0062a237`): the polymorphic-array TBP dispatch site `realm(is)%end_stage_forest(...)` segfaults inside libnvf's `pgf90_copy_f90_argl_i8` while marshalling the explicit-bound `q_gpu(1:, 1-self%ngc:, ...)` actual passed deep in the body to `compute_residuals_dev`. The fix is a one-line `associate(r => realm(is))` wrapper at the call site — binds the array element to a polymorphic scalar before the TBP fires, letting nvfortran resolve the element's concrete-type stride correctly. Pre-emptively applied at every polymorphic-array TBP dispatch site in `evolve_one_step` (Phase 1 `begin_stage_forest`, Phase 2 `fill_seam_from_peer_forest`, Phase 3 `end_stage_forest`, Phase 4 `close_step_forest`, Phase 5 `fill_seam_from_peer_forest`). gfortran does not exhibit the issue; CPU builds are unaffected.

Tracked through:

| PRD | Scope | Status |
|---|---|---|
| [#10](https://github.com/szaghi/adam/issues/10) | Phase D inception — forest orchestrator design | shipped |
| [#13](https://github.com/szaghi/adam/issues/13) | Coarse-fine interface machinery (A/B/C) | A shipped; B & C deferred |
| [#16](https://github.com/szaghi/adam/issues/16) | α end-of-step barrier + asymmetric K | shipped (closed 2026-05-28) |
| [#18](https://github.com/szaghi/adam/issues/18) | β stage-coincident recovery + cross-config oracle | shipped (closed 2026-05-28) |

User-facing conceptual overview, manifest schema, and α/β cadence rationale are documented in `docs/guide/forest.md` (the project VitePress site). Library-developer contract reference is in `src/lib/common/README.md` → "Forest orchestration".

### `grid_object` (`src/lib/common/adam_grid_object.F90`)

```fortran
real(R8P)    :: domain_emin(3), domain_emax(3)    ! Domain bounding box
integer(I4P) :: ni, nj, nk                        ! Interior cells per block per direction
integer(I4P) :: ngc                               ! Ghost cell width (typ. 3 or 4)
integer(I4P) :: block_weight, weight_neighbor(26)
integer(I4P) :: bc_type(6)                        ! Face boundary condition types
logical      :: is_ijk_periodic(3)
logical      :: null_xyz(3)
real(R8P), allocatable :: block_dxyz(:,:)         ! [3, 0:MAX_REF_LEVELS]
real(R8P), allocatable :: cell_dxyz(:,:)
integer(I4P), allocatable :: nb_max(:)
real(R8P), allocatable :: lin_space_x(:,:), lin_space_y(:,:), lin_space_z(:,:)
```

### `field_object` (`src/lib/common/adam_field_object.F90`)

```fortran
integer(I4P) :: nv, nv_pic          ! Variable counts
integer(I4P) :: nb, blocks_number   ! Block capacity and actual count
integer(I4P) :: np                  ! Particles per block
! Morton-order block arrays:
integer(I8P), allocatable :: code(:)              ! [nb]
integer(I4P), allocatable :: coordinates(:,:)     ! [4, nb]  IJKL
real(R8P),    allocatable :: emin(:,:), emax(:,:) ! [3, nb]  block bounds
real(R8P),    allocatable :: dxyz(:,:)            ! [3, nb]  cell spacing
real(R8P),    allocatable :: x_cell(:,:), y_cell(:,:), z_cell(:,:)
! MPI:
integer(I4P), allocatable :: blocks_numbers(:), disp_count(:)
integer(I4P), allocatable :: refinements_needed(:), refinements_needed_all(:)
! Working storage:
real(R8P),    allocatable :: q_work(:,:,:,:,:)
real(R8P),    allocatable :: residuals(:)
! Scalar replicas (pointers into grid):
integer(I4P), pointer :: ngc, ni, nj, nk
```

### `weno_object` (`src/lib/common/adam_weno_object.F90`)

```fortran
character(:), allocatable :: scheme            ! e.g. "weno-z", "weno-js"
integer(I4P)              :: S                 ! Stencil half-width; order = 2S-1
logical                   :: is_centered
real(R8P), allocatable :: a(:,:,:)             ! Optimal weights  [1:2, 0:S-1, 1:S]
real(R8P), allocatable :: p(:,:,:,:)           ! Poly coefficients [1:2, 0:S-1, 0:S-1, 1:S]
real(R8P), allocatable :: d(:,:,:,:)           ! Smoothness indicators [0:S-1, 0:S-1, 0:S-1, 1:S]
real(R8P), allocatable :: c(:,:)               ! Centered polynomials [1-S:S, 1:S]
integer(I4P) :: ror_number, wexp
real(R8P)    :: zeps, ror_threshold
```

### `ib_object` (`src/lib/common/adam_ib_object.F90`)

```fortran
integer(I4P)              :: solids_number
character(99), allocatable :: s_name(:), definition(:)
integer(I4P),  allocatable :: bc_type(:)           ! BCS_VISCOUS or BCS_EULER
real(R8P),     allocatable :: q(:,:)               ! BC variables
type(analytical_sphere_object),    allocatable :: sphere(:)
type(analytical_rectangle_object), allocatable :: rectangle(:)
integer(I4P)              :: n_eikonal             ! Eikonal pseudo-time steps
! Distance function (solids+1 for combined):
real(R8P), allocatable :: phi(:,:,:,:,:)           ! [solids+1, ni+2*ngc, nj+2*ngc, nk+2*ngc, nb]
```

### `rk_object` (`src/lib/common/adam_rk_object.F90`)

```fortran
character(:), allocatable :: scheme            ! e.g. "rk-ssp-3", "rk-ls-4"
integer(I4P)              :: nrk               ! Number of RK stages
real(R8P), allocatable :: ark(:), brk(:), crk(:)       ! Low-storage coefficients
real(R8P), allocatable :: alph(:,:), beta(:), gamm(:)  ! Butcher tableau
real(R8P), allocatable :: ssa(:), ssb(:)               ! Symplectic splitting
real(R8P), allocatable :: q_rk(:,:,:,:,:,:)    ! Stage storage [nv, ni+2*ngc, nj+2*ngc, nk+2*ngc, nb, nrk]
integer(I4P), pointer :: ngc, ni, nj, nk
```

### `mpih_object` (`src/lib/common/adam_mpih_object.F90`)

```fortran
integer(I4P)              :: myrank, procs_number
character(:), allocatable :: myrankstr
real(R8P)    :: memory_avail           ! Available host memory (GB)
integer(I4P) :: error
real(R8P)    :: timing(1:2)            ! tic-toc wall time
integer(I4P) :: tictoc
```

### `field_fnl_object` (`src/lib/fnl/adam_fnl_field_object.F90`)

```fortran
type(maps_fnl_object) :: maps           ! GPU communication maps
! Device arrays (OpenACC device_resident):
integer(I4P), pointer :: fec_1_6_array_gpu(:) => null()
real(R8P),    pointer :: x_cell_gpu(:,:), y_cell_gpu(:,:), z_cell_gpu(:,:)
real(R8P),    pointer :: dxyz_gpu(:,:)
! Scalar replicas (point into CPU field singleton):
integer(I4P), pointer :: ngc, ni, nj, nk, nb, blocks_number, nv
```

### `ib_fnl_object` (`src/lib/fnl/adam_fnl_ib_object.F90`)

```fortran
real(R8P), pointer :: q_bcs_vars_gpu(:,:)   ! BC variables on device
real(R8P), pointer :: phi_gpu(:,:,:,:,:)    ! [solids+1, ni+2ngc, nj+2ngc, nk+2ngc, nb]
```

### `nasto_common_object` (`src/app/nasto/common/adam_nasto_common_object.F90`)

Key members (composite, not extending `realm_object`):

```fortran
type(mpih_object)               :: mpih
type(adam_object)               :: adam
type(field_object),   pointer   :: field => null()
type(grid_object),    pointer   :: grid  => null()
type(amr_object)                :: amr
type(ib_object)                 :: ib
type(rk_object)                 :: rk
type(weno_object)               :: weno
type(slices_object)             :: slices
type(nasto_io_object)           :: io
type(nasto_physics_object)      :: physics
type(nasto_ic_object)           :: ic
type(nasto_bc_object)           :: bc
type(nasto_time_object)         :: time
real(R8P), allocatable :: q(:,:,:,:,:)      ! [nv, ni+2ngc, nj+2ngc, nk+2ngc, nb]
real(R8P), allocatable :: q_aux(:,:,:,:,:)
integer(I4P), pointer :: ngc, ni, nj, nk, nb, blocks_number, ns, nv, nv_aux
```

### `nasto_fnl_object` (`src/app/nasto/fnl/adam_nasto_fnl_object.F90`)

Extends `nasto_common_object`. Additional GPU members:

```fortran
real(R8P), pointer :: q_gpu(:,:,:,:,:)
real(R8P), pointer :: q_aux_gpu(:,:,:,:,:)
real(R8P), pointer :: dq_gpu(:,:,:,:,:)      ! eikonal RHS buffer
real(R8P), pointer :: flx_gpu(:,:,:,:,:)     ! x-fluxes
real(R8P), pointer :: fly_gpu(:,:,:,:,:)     ! y-fluxes
real(R8P), pointer :: flz_gpu(:,:,:,:,:)     ! z-fluxes
real(R8P), pointer :: q_bc_vars_gpu(:,:)
```

Key procedures: `initialize(self, filename)`, `simulate`, `integrate`, `compute_residuals`, `compute_dt`, `update_ghost`, `set_boundary_conditions`, `set_initial_conditions`, `load_restart_files`, `save_hdf5`, `save_restart_files`, AMR methods, `integrate_eikonal`.

### `prism_common_object` (`src/app/prism/common/adam_prism_common_object.F90`)

Extends `realm_object`. Additional PRISM-specific members:

```fortran
type(prism_numerics_object)            :: numerics
type(prism_physics_object)             :: physics
type(prism_ic_object)                  :: ic
type(prism_bc_object)                  :: bc
type(prism_rk_bc_object)               :: rk_bc
type(prism_time_object)                :: time
type(prism_fwlayer_object)             :: fwlayer
type(prism_coil_object)                :: coil
type(prism_external_fields_object)     :: external_fields
type(prism_pic_object)                 :: pic
type(prism_leapfrog_pic_object)        :: leapfrog_pic
type(prism_rk_pic_object)              :: rk_pic
real(R8P), allocatable :: q(:,:,:,:,:)       ! [nv, ni+2ngc, nj+2ngc, nk+2ngc, nb]
real(R8P), allocatable :: dq(:,:,:,:,:)
real(R8P), allocatable :: curl(:,:,:,:,:), divergence(:,:,:,:,:)
real(R8P), allocatable :: energy_D(:), energy_B(:), coil_power(:)
integer(I4P), pointer  :: nv_c, nv_s, nv_cl  ! variable-count replicas
```

### `prism_fnl_object` (`src/app/prism/fnl/adam_prism_fnl_object.F90`)

Extends `prism_common_object`. Additional GPU members:

```fortran
real(R8P), pointer :: q_gpu(:,:,:,:,:)
real(R8P), pointer :: dq_gpu(:,:,:,:,:)
real(R8P), pointer :: flxyz_c_gpu(:,:,:,:,:,:,:)  ! cell-center fluxes [nv, ni+2ngc, nj+2ngc, nk+2ngc, nb, 2, 3]
real(R8P), pointer :: flx_f_gpu(:,:,:,:,:)         ! face fluxes x
real(R8P), pointer :: fly_f_gpu(:,:,:,:,:)
real(R8P), pointer :: flz_f_gpu(:,:,:,:,:)
real(R8P), pointer :: curl_gpu(:,:,:,:,:)
real(R8P), pointer :: divergence_gpu(:,:,:,:,:)
! Device FDV procedure pointers:
procedure(compute_curl_interface_dev),      pointer :: compute_curl_dev
procedure(compute_gradient_interface_dev),  pointer :: compute_gradient_dev
procedure(compute_divergence_interface_dev),pointer :: compute_divergence_dev
procedure(compute_laplacian_interface_dev), pointer :: compute_laplacian_dev
procedure(compute_residuals_interface_dev), pointer :: compute_residuals_dev
procedure(integrate_interface_dev),         pointer :: integrate_dev
```

Key procedures: `initialize_prism(self, filename)`, `simulate`, `compute_dt`, `compute_energy`, `compute_energy_error`, `update_ghost`, `update_rk_ghost`, `apply_fwl_correction`, `compute_coils_current`, `impose_ct_correction`, `impose_div_free`, `save_simulation_data`.

---

## Library Aggregation Chain

```
adam_common_library
 ├── adam_adam_object, adam_field_object, adam_grid_object, adam_tree_object
 ├── adam_maps_object, adam_weno_object, adam_rk_object, adam_ib_object, adam_mpih_object
 ├── adam_realm_object, adam_forest_object, adam_io_object, adam_amr_object, adam_refinement_plan_object
 ├── adam_eos_ic_object, adam_cfm_object, adam_blanes_moan_object
 ├── adam_leapfrog_object, adam_flail_object, adam_slices_object
 ├── adam_tree_node_object, adam_tree_bucket_object
 ├── adam_fdv_operators_library, adam_riemann_euler_library
 └── [8 CPU singletons]: mpih_global, grid_global, field_global, maps_global,
                          tree_global, weno_global, ib_global, rk_global
                         (also bundled by `adam_globals` convenience aggregator)

adam_fnl_library
 ├── adam_common_library  (full CPU layer)
 ├── adam_fnl_field_object, adam_fnl_ib_object, adam_fnl_weno_object
 ├── adam_fnl_rk_object, adam_fnl_maps_object, adam_fnl_mpih_object
 ├── adam_fnl_field_kernels, adam_fnl_ib_kernels, adam_fnl_rk_kernels, adam_fnl_weno_kernels
 ├── adam_fnl_fdv_operators_library
 └── [5 FNL singletons]: fnl_mpih_global, fnl_field_global, fnl_weno_global,
                          fnl_ib_global, fnl_rk_global

adam_nasto_common_library
 ├── adam_nasto_bc_object, adam_nasto_common_object, adam_nasto_eos_object
 ├── adam_nasto_ic_object, adam_nasto_io_object, adam_nasto_physics_object, adam_nasto_time_object
 └── adam_nasto_parameters

adam_nasto_fnl_library
 ├── adam_nasto_common_library
 ├── adam_fnl_library
 └── adam_nasto_fnl_cns_kernels, adam_nasto_fnl_kernels

adam_prism_common_library
 ├── adam_prism_common_object, adam_prism_bc_object, adam_prism_coil_object
 ├── adam_prism_fwlayer_object, adam_prism_ic_object, adam_prism_numerics_object
 ├── adam_prism_physics_object, adam_prism_pic_object, adam_prism_particle_injection_object
 ├── adam_prism_external_fields_object, adam_prism_leapfrog_pic_object
 ├── adam_prism_rk_bc_object, adam_prism_rk_pic_object, adam_prism_io_object
 ├── adam_prism_time_object, adam_prism_riemann_library, adam_prism_parameters
 └── (uses adam_common_library via realm_object chain)

adam_prism_fnl_library
 ├── adam_prism_common_library
 ├── adam_fnl_library
 ├── adam_prism_fnl_coil_object,   adam_prism_fnl_coil_global
 ├── adam_prism_fnl_fwlayer_object, adam_prism_fnl_fwlayer_global
 └── adam_prism_fnl_external_fields_kernels
```

---

## Design Patterns

### Singleton Access Pattern

```fortran
! Any module needing the grid dimensions:
use :: adam_grid_global, only: grid
use :: adam_field_global, only: field
...
associate(ni=>grid%ni, nj=>grid%nj, ngc=>grid%ngc, nb=>field%nb)
  ...
endassociate
```

Never: `subroutine foo(self, grid, field)` — singletons are not passed as arguments.

### Scalar Replica Pointers

Several types cache frequently used integers as pointer members (e.g. `integer(I4P), pointer :: ngc`) pointing into the canonical storage in `grid_object` or `field_object`. These are set once at initialization via `ngc => grid%ngc`. Do not read these via `self%ngc` in new code — use the singleton directly.

### 5D Field Array Layout

All physical field arrays use shape `(nv, ni+2*ngc, nj+2*ngc, nk+2*ngc, nb)`:
- `nv` is stride-1 (innermost, contiguous) — variables loop innermost
- Ghost cell bounds: `1-ngc : ni+ngc` in each spatial direction
- `nb` blocks are independent — the outermost loop for block parallelism

### AMR Decoupling via Transfer Object

`refinement_plan_object` is a one-shot DTO produced by `tree%update()` and consumed once by `field%update()`, decoupling the tree and field AMR cycles. It carries block-to-refine / block-to-derefine index arrays.

### Write-Back Composition Exception

`prism_fnl_coil_object` and `prism_fnl_fwlayer_object` retain a CPU-side `pointer` member (`self%coil`, `self%fwlayer`) because their `copy_gpu_cpu` procedures must write results back to the host object. This is the only FNL pattern where a type holds a pointer to another type's data.

### `evolve_eikonal` API

`ib_fnl%evolve_eikonal` takes an explicit `dxyz_gpu(:,:)` parameter (not accessed via `field_fnl` internally) to avoid circular module dependencies. Callers pass `field_fnl%dxyz_gpu`.

---

## Backend Summary

| Backend | Macro | Compiler | GPU API | Applies to |
|---------|-------|----------|---------|------------|
| cpu | (none) | gfortran / ifort | — | NASTO, PRISM, CHASE, PATCH |
| nvf | `_NVF` | nvfortran | CUDA Fortran | NASTO |
| fnl | `_FNL` | nvfortran | OpenACC (FUNDAL) | NASTO, PRISM |
| gmp | `_GMP` | nvfortran / ifort | OpenMP target | NASTO (experimental) |

Macros are mutually exclusive. Only one of `_NVF`, `_FNL`, `_GMP` is defined per build.
