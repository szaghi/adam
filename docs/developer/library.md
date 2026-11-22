# Library objects

`src/lib` holds the ADAM SDK: the physics-agnostic core that every application is built on. This page is the reference
for the common layer (`src/lib/common`) and for the realm/forest contract an application implements. The GPU layer is
on [FNL library](./fnl); the elliptic solvers are on [FLAIL](./flail).

## Directory overview

| Directory | Role | Parallelism | Built by |
|-----------|------|-------------|----------|
| [`common/`](#amr-data-structures) | Portable core objects: grid, tree, field, numerics, I/O, forest | CPU, MPI, OpenMP | every mode |
| [`fnl/`](./fnl) | GPU backend on FUNDAL | OpenACC or OpenMP target, MPI | PRISM and FLUME FNL modes |
| [`nvf/`](https://github.com/szaghi/adam/tree/develop/src/lib/nvf) | CUDA Fortran backend | CUDA Fortran, MPI | no mode |
| [`gmp/`](https://github.com/szaghi/adam/tree/develop/src/lib/gmp) | OpenMP target offloading backend | OpenMP target, MPI | no mode |

The SDK has two axes of variability:

- **Physics** belongs to the application (`src/app`). The SDK knows no PDE.
- **Hardware** is chosen at compile time. `_FNL` selects the GPU backend: OpenACC with nvfortran (`*-fnl-nvf` modes) or
  the same kernels as OpenMP offload with amdflang (`*-fnl-omp-amd` modes, `DEV_OMP`).

The common layer has no GPU code and compiles with any Fortran 2003+ compiler. `use adam_common_library` exposes the
whole API (the exceptions are listed in the [module summary](#module-summary)).

## Program-scope singletons

The only program-scope singletons are the MPI handlers, one per process: `mpih` (`adam_mpih_global`, `mpih_object`)
and, on the FNL backend, `mpih_fnl` (`adam_fnl_mpih_global`), which also owns the device context. Every other object
is a component of a realm (`self%adam%grid`, `self%adam%field`, `self%adam%maps`, `self%rk`, `self%weno`, `self%ib`,
`self%flail`, …) and reaches library routines as an argument:

```fortran
call self%adam%field%update_ghost_local(q=q, grid=self%adam%grid, maps=self%adam%maps)
```

Do not add `*_global` modules. `adam_globals` re-exports `mpih` alone, without the type-definition modules.

## AMR data structures

### `adam_tree_object`

The octree is a CPU-resident hash map keyed on 64-bit Morton indices, which linearise `(level, bx, by, bz)` into one
integer. Insertion, deletion and neighbour lookup are O(1). Morton order clusters adjacent blocks in index space, which
keeps the ghost exchange between ranks small.

The refinement ratio is 2 (binary tree), 4 (quadtree) or 8 (octree). The 26-neighbourhood covers every face, edge and
corner adjacency a high-order stencil needs. Kernels never read the tree during time integration; only the AMR cycle
(marking, refinement and coarsening, rebalancing) does.

### `adam_tree_node_object`

One block of the tree.

| Field | Type | Purpose |
|-------|------|---------|
| `code` | `I8P` | Morton key |
| `refinement_needed` | `I4P` | AMR marker (`TO_BE_REFINED`, `TO_BE_DEREFINED`, `TO_NOT_TOUCH`) |
| `myrank` | `I4P` | Owning MPI rank |
| `block_index` | `I8P` | Compact index into the field array |
| `neighbor(26)` | `tree_node_neighbor_object` | Face, edge and corner adjacency |
| `surface_stl_distance` | `R8P` | Signed distance to the nearest IB surface |

### `adam_tree_bucket_object`

The doubly linked list behind each hash-map bucket: `add_node`, `remove_node` and `has_code` in O(1), and `traverse`
with a callback conforming to `iterator_interface`.

### `adam_refinement_plan_object`

A one-shot transfer object of one AMR cycle. `tree_object%adapt` produces it and `field_object%adapt` consumes it, so
the tree holds no field-block lists and does not depend on the field.

| Component | Content |
|-----------|---------|
| `ratio` | Refinement ratio (2, 4, 8) |
| `block_to_refine(2, n)`, `block_refined(2, ratio n)` | Field blocks to refine, and the children with their Morton codes |
| `block_to_derefine(n)`, `block_derefined(2, n/ratio)` | Field blocks to derefine, and the parents with their Morton codes |
| `has_refinements()`, `has_derefinements()`, `destroy()` | Whether this rank has work of each kind; cleanup |

### `adam_amr_object`

The `[amr]` section and its `[amr_marker_N]` markers. Markers combine; each carries its own thresholds and field
variable. At initialization the markers refine the grid, iterated `iters` times to enforce 2:1 balance. A realm that
regrids during the run (FLUME) re-applies them every `frequency` steps through `regrid_forest`; PRISM refines at
initialization only.

| Constant | Mode | Criterion |
|----------|------|-----------|
| `AMR_GEO` | 1 | Geometry: an IB solid (`AMR_GEO_SOLID`), a primitive box up to `target_level` (`AMR_GEO_PRIMITIVE_BOX`) or an STL surface (`AMR_GEO_STL`, refused by PRISM and FLUME) |
| `AMR_GRAD` | 2 | Field gradient magnitude |
| `AMR_TV` | 3 | Field total variation |
| `AMR_LOHNER` | 4 | Löhner normalised second derivative (`refine_tol`, `derefine_tol`) |
| `AMR_JEANS` | 5 | Jeans length resolution (`jeans_cells`, `jeans_cells_derefine`) |
| `AMR_USER` | 6 | A criterion the program registered (`name`) |

PRISM implements the primitive box and `AMR_TV`; its `AMR_GEO_SOLID` and `AMR_GRAD` branches do nothing. FLUME
implements the box, solid, gradient, Löhner, Jeans and registered criteria, and stops at initialization on any other.

### `adam_amr_criterion_object`

A refinement criterion written as code. It is a type extending the abstract `amr_criterion_object`, registered with
the realm before the simulation starts and selected by `[amr_marker_N] mode = 6, name = <its name>`. Neither the
library nor the app is edited: the program (a copy of the app's driver) registers it.

| Member | Role |
|--------|------|
| `name` | The name an input selects |
| `vote(realm, b)` (deferred) | `TO_BE_REFINED`, `TO_NOT_TOUCH` or `TO_BE_DEREFINED` for block `b` of `realm` |
| `mark(realm, votes, initial)` | The votes of every block of the rank; the default calls `vote` block by block. A criterion that compares neighbours or dilates its marks overrides it |
| `runtime_only` | Skipped by the AMR of the initial condition (default `.true.`) |

Rules in FLUME:

| Rule | Detail |
|------|--------|
| Registration | `register_amr_criterion(criterion)` on every realm of a forest. The realm keeps a pointer: the instance needs the `target` attribute and must outlive the run (a main-program or `save` variable) |
| Unknown name | Stops the run at initialization |
| Combination | Refine if any marker votes refine, derefine if all do, never below `iu_ref_levels` |
| Where it runs | On the host, for the blocks of the rank, at runtime regrids. `realm%adam` gives the grid, the tree and the block metrics; a criterion that reads the app's fields selects the app's type |
| One list | FLUME builds `amr_marker_list` from its `[amr_marker_N]` sections in input order; the built-in markers are criteria of the same kind (`flume_box_marker` and the solid, gradient, Löhner and Jeans markers) |

Example and test: `src/tests/flume/verification/regrid/test_flume_user_criterion.F90` (mode
`test-flume-user-criterion-cpu-gnu`), a criterion that votes as the box marker; the run is bitwise the same on one
realm and on two (`multirealm/check.sh` leg 10).

## Grid and field

### `adam_grid_object`

The geometry of the blocks and of the domain: cell and node coordinates, metrics, ghost layout.

| `[grid]` key | Meaning |
|--------------|---------|
| `ni`, `nj`, `nk` | Cells per block per direction |
| `ngc` | Ghost-cell width (3 for WENO5, 4 for WENO7) |
| `emin_*`, `emax_*` | Domain extents |

The face boundary types `bc_type(6)` are not read from `[grid]`: the application sets them with `set_bc_type`.
Methods: `block_emin`/`block_emax`, `cell_xyz`, `node_xyz`, `compute_metrics`, `compute_weight_neighbor`,
`load_from_ini_file`.

### `adam_field_object`

Per-block geometry (`emin`, `emax`, `dxyz`, node and cell coordinates), the ghost exchange (`update_ghost_local`,
`update_ghost_mpi`, and the reverse folds `reduce_ghost_local`, `reduce_ghost_mpi`), the refinement and derefinement
of block data (`adapt`), and the restart blocks (`save_blocks`, `load_blocks`).

The state arrays belong to the application (e.g. `flume_common_object%q`), dense and contiguous over the blocks of a
rank:

```
q(nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb)
```

| Dimension | Meaning |
|-----------|---------|
| `nv` | Physical variables |
| `ni, nj, nk` | Cell indices in a block, ghosts included |
| `nb` | Compact block index (1 … block count) |

`nv` is the contiguous index. The application allocates `q` once, to the block count that `realm_object%initialize`
derives from the memory budget, and reuses it in every Runge-Kutta stage: nothing is allocated during time
integration.

## Communication and maps

### `adam_maps_object`

The bridge between the tree (Morton keys) and the field (compact indices). Rebuilt only when the AMR topology changes.

| Local map | Purpose |
|-----------|---------|
| `local_map` | Block index changes of the rank's nodes (old → new) |
| `local_map_ghost` | Block-to-ghost-block adjacency |
| `local_map_ghost_cell` | Cell-level ghost stencil offsets |
| `local_map_bc_face/edge/corner/crown` | Boundary-condition stencils |

| MPI list | Purpose |
|----------|---------|
| `comm_map_send/recv` | Block exchange for rebalancing |
| `comm_map_send/recv_ghost` | Ghost-cell exchange |
| `send_buffer_ghost`, `recv_buffer_ghost` | Preallocated MPI staging buffers |

Methods: `initialize`, `make_comm_local_maps`, `blocks_reorder`. The seam maps of the forest (`seam_local_*`,
`seam_mpi_*`, `amr_seam_quadrant`) are described under [forest orchestration](#forest-orchestration-multi-realm).

## Numerical methods: spatial

### `adam_weno_object`

The WENO configuration and the run-time dispatchers over [SETH](https://github.com/szaghi/SETH)'s WENO kernels
([#139](https://github.com/szaghi/adam/issues/139)). SETH is the only source of WENO, TENO and hybrid operators: every
coefficient table, smoothness indicator and nonlinear weight is SETH's, in kernels with compile-time coefficients
(`pure`, `acc routine seq`, `omp declare target`), so there is nothing to mirror on a device. ADAM keeps the policy:
`weno_object` holds `[weno]` (`scheme`, its half-width `S`, `weights` `js` or `si` and their descaler switch `sigma`,
the regulariser `eps`), and the module holds `weno_reconstruct`, `weno_reconstruct_wr` and `weno_reconstruct_interp_wr`,
which select SETH's kernel of width `S` (`seth_weno<2S-1>_js[_wr|_interp_wr]`) for callers whose `S` is a run-time value
per call (RIPPLE). FLUME and PRISM call SETH's kernel of each width from kernels compiled with `S` a named constant.
SETH's stencil convention: cells by offset from the face, `vp(1-S:S)` left-biased, `vm(1-S:S)` right-biased, both
unshifted (ADAM's former packed rows mapped as `vp(m) = v(2,m)`, `vm(m+1) = v(1,m)`).

| Constant | Scheme |
|----------|--------|
| `WENO_U_1` … `WENO_U_9` | Upwind WENO, orders 1/3/5/7/9 (S = 1 … 5) |
| `WENO_C_2` … `WENO_C_8` | Centred names: refused at initialization until an application calls SETH's centred kernels (D6) |

ADAM's own tables and kernels (host, FNL generic, FNL per width, NVF and GMP twins), the centred coefficients, the ROR
and immersed-boundary order reduction and `cell_scheme` were deleted in S6; the substitution was measured before
(`test_weno_seth_parity`, since deleted: bitwise on the host with gfortran for every width, both weights, reconstruction
and interpolation, KEEP-PE and Ducros; on the device SETH's kernels equal ADAM's FNL per-width kernels bitwise), and the
application A/B runs are bitwise. `scripts/check_space_separation.py` (CI leg `space-separation`) keeps a
smoothness-indicator constant, a nonlinear weight, the KEEP two-point mass flux or a name SETH replaced (`weno_a/p/d`,
`weno_reconstruct_upwind*`, `compute_kep_flux`) from being written in `src/app` or `src/lib`; the one waiver left is the
seam WENO-6 transfer, which S7 moves to SETH.

### `adam_fdv_operators_library`

Gradient, divergence, curl, Laplacian and scalar derivatives of orders 1–6, each in three flavours: centred finite
difference, centred finite volume, upwind (WENO) finite volume. The abstract `compute_*_fdv_interface` procedure
pointers dispatch between flavours without branching.

### `adam_riemann_euler_library`

Euler Riemann solvers on left/right primitive states, with the optional maximum wave speed `lmax`. No maintained
application calls them.

| Procedure | Solver |
|-----------|--------|
| `compute_riemann_euler_exact` | Exact (Newton iteration) |
| `compute_riemann_euler_hllc` | HLLC |
| `compute_riemann_euler_hllc_lm` | HLLC with low-Mach correction |
| `compute_riemann_euler_hllem` | HLLEM |
| `compute_riemann_euler_llf` | Local Lax–Friedrichs |
| `compute_riemann_euler_ts` | First-order upwind |

## Numerical methods: temporal

### `adam_stepper_object`

The time integrator adapter on [SELAH](https://github.com/szaghi/SELAH), the library of time integrators written for
ADAM ([#132](https://github.com/szaghi/adam/issues/132), SELAH milestone M7). SELAH owns the coefficients, the
registers and the stage kernels; the adapter is the ADAM side of it:

- **Input.** `[runge_kutta] scheme` (no other key), the names of the table below, defined here.
- **Mapping** (`stepper_map`). `runge-kutta-1/2/3` are the two-register Shu–Osher forms of Euler, SSP-22 and SSP-33
  (ADAM's low-storage coefficients, bitwise); `runge-kutta-ssp-11/22/33/54` are the Butcher forms of Euler, SSP-22,
  SSP-33 and SSP-54 (ADAM's tableaux, bitwise), each stage on its own state register as in ADAM's former
  integrator (`adam_rk_object`, deleted in phase 6 of #132);
  `runge-kutta-yoshida` is SELAH's Yoshida-4 splitting, which needs the variables of its two partitions.
- **Storage** (`stepper_stored_fields`). SELAH's `registers_needed`, the residual register included: 2 for the
  two-register forms, s + 1 for the Butcher forms (one field more than the former integrator's `nrk` stages, whose
  residual went to the application's `dq`). The registers are allocated once, at the pool size, and each block slot is zeroed
  the first time it is active ([#45](https://github.com/szaghi/adam/issues/45)).
- **Field and mask.** Each call takes the field and builds the box of the active blocks (`field%blocks_number` at the
  time of the call: it changes at every regrid). `ib%phi` is the mask: a cell is advanced where the last component (the
  all-solids distance) is negative, as in the former integrator; in the Butcher form the masked cells of each stage state
  start as `q`, where the immersed-boundary extrapolation of the residual starts from.
- **Step.** `begin_step`, then per stage `begin_stage(k)`, the residual of the stage state (`q` itself when
  `state_is_q(k)`, else `stage_state()`) written into `rhs_register()`, `end_stage(k)`; then `finish_step`, which in
  the Butcher form advances `q` and returns the step residual `dq` (0 in the masked cells). `stage_weight(k)` is the
  effective Butcher weight (the flux register's), `stage_time(k)` the abscissa.
- **Scratch.** `scratch(n)` gives register `n` between steps (every Runge–Kutta form has at least two): the D_flux
  diagnostic and the equilibrium residual borrow registers 1 and 2. A step rebuilds what it reads of them.

The same adapter serves the device fields of the FNL backends: `layout=LAYOUT_V_LAST` (`realm_object%selah_layout`), the
registers in device memory (zeroed at allocation), its own masked-cell copy a FUNDAL kernel, and `set_rhs` (a residual
computed into the application's own array copied into the residual register with SELAH's copy kernel). Tests:
`test-stepper-object-gnu` (the mapping, the stored fields, one step per scheme, masked and unmasked, the scratch
registers) and `test-stepper-object-xfail-gnu` (an unknown scheme stops the run).

### Scheme names

The `[runge_kutta] scheme` names (`adam_stepper_object` constants). ADAM's own integrators (`adam_rk_object`, the FNL
`adam_fnl_rk_object` and its kernels, the NVF and GMP twins, PRISM's boundary clone) were deleted in phase 6 of
[#132](https://github.com/szaghi/adam/issues/132): every application integrates with SELAH, and
`scripts/check_time_separation.py` (CI) keeps a Runge–Kutta combination from being written again in `src/app` or
`src/lib`: a time step times a stage coefficient or a stage-indexed element, a literal fraction or multiple of the
time step, a direct read of the stepper's tableau, a literal SSP-RK3/RK4 stage. What remains is waived statement by
statement, each with the issue that removes it ([#134](https://github.com/szaghi/adam/issues/134)): PRISM's PIC particles,
PML auxiliary fields and charge-conserving current ([SELAH#4](https://github.com/szaghi/SELAH/issues/4)), FLUME's sink
kick and positivity limiter ([SELAH#5](https://github.com/szaghi/SELAH/issues/5)), PRISM's t0 centring drift (kept).

| Constant | Scheme | Stages | Order |
|----------|--------|--------|-------|
| `RK_1` | Forward Euler | 1 | 1 |
| `RK_2` | Heun | 2 | 2 |
| `RK_3` | TVD-RK3 | 3 | 3 |
| `RK_SSP_11` | SSP(1,1) | 1 | 1 |
| `RK_SSP_22` | SSP(2,2) | 2 | 2 |
| `RK_SSP_33` | SSP(3,3) | 3 | 3 |
| `RK_SSP_54` | SSP(5,4) | 5 | 4 |
| `RK_YOSHIDA` | Yoshida symplectic | 4 | 4 |

### Removed integrators

The leapfrog with the Robert–Asselin–Williams filter, the Blanes–Moan splitting and the commutator-free Magnus
integrators (`adam_leapfrog_object`, `adam_blanes_moan_object`, `adam_cfm_object`), never initialized by any
application, were deleted ([#134](https://github.com/szaghi/adam/issues/134)) and move to SELAH
([SELAH#3](https://github.com/szaghi/SELAH/issues/3)), which re-implements them against the theory (ADAM's RAW filter,
Blanes–Moan loop and CFM exponentials were defective). PRISM refuses `scheme_time = leapfrog`, `blanes_moan` and
`commutator_free_magnus` at initialization.

## Physics support

### `adam_ib_object`

Solid bodies immersed in the Cartesian grid. The signed distance `phi(nv, ni, nj, nk, nb)` comes from an eikonal
solver; wall conditions are applied at the solid surface. Sections: `[solids]` (`number`, `n_eikonal` eikonal steps)
and one `[solid_N]` per solid.

| Solid | Constant | Definition |
|-------|----------|------------|
| Sphere | `IB_ANALYTICAL_SPHERE` | Centre, radius, axis |
| Circle (xy plane) | `IB_ANALYTICAL_CIRCLE` | 2-D equivalent |
| Rectangle | `IB_ANALYTICAL_RECTANGLE` | Centre, edge lengths, axis |
| Surface mesh | — | `file.off`: not implemented, refused at load |

| Wall constant | Value | Type |
|---------------|-------|------|
| `BCS_VISCOUS` | 1 | No-slip wall |
| `BCS_EULER` | 2 | Slip wall |

### `adam_turbulence_object`, `adam_turbulence_library`, `adam_turbulence_kernels`

Algebraic subgrid-scale models of a realm, `[turbulence]` section. The app owns the object (a realm component), the
eddy-viscosity field and the choice of the model kernel.

| Key | Values | Default |
|-----|--------|---------|
| `sgs_model` | `none`, `vreman`, `wale`, `sigma` (any other value is fatal) | `none` |
| `sgs_constant` | `> 0` | 0.07 Vreman, 0.5 WALE, 1.35 Sigma |
| `prandtl_t` | `> 0` | 0.9 |

| Module | Content |
|--------|---------|
| `adam_turbulence_library` | Pointwise `sgs_vreman`, `sgs_wale`, `sgs_sigma`: the kinematic eddy viscosity from the velocity gradient, the constant, the steps and the filter width. `pure`, `!$acc routine seq` + `!$omp declare target`, shared by CPU and FNL |
| `adam_turbulence_kernels` | CPU `compute_mu_t_vreman`, `_wale`, `_sigma`: `mu_t = rho nu_t` on the interior and on 2 ghost layers per active direction (the 4th-order dissipative fluxes read cells -1 … 2), 2nd-order central gradients (`ngc >= 3`), filter width the geometric mean of the block's active steps (it jumps at a 2:1 face). One kernel per model, selected on the host |

| Model | Eddy viscosity | Vanishes in |
|-------|----------------|-------------|
| Vreman (2004) | `c sqrt(B / (a_ij a_ij))`, `c = 2.5 C_s^2` | rank-1 gradients, pure shear among them |
| WALE (Nicoud & Ducros 1999) | `(c delta)^2 (Sd:Sd)^(3/2) / ((S:S)^(5/2) + (Sd:Sd)^(5/4))` | pure shear |
| Sigma (Nicoud et al. 2011) | `(c delta)^2 s_3 (s_1 - s_2)(s_2 - s_3) / s_1^2` | 2-D flow, pure shear, solid rotation |

The constants were calibrated on other discretizations, not on WENO. The FNL twin is `adam_fnl_turbulence_kernels`.

### `adam_eos_ic_object`

A calorically perfect ideal gas (`cp`, `cv`, `g`, the derived `gm1`, `gp1`, `delta`, `eta`, viscosity `mu`,
conductivity `kd`; conversions between primitive and conservative variables), `[physics_specie_*]` sections.
`realm_object` imports it; no maintained application uses it.

### Elliptic solvers

`adam_flail_object` (the realm component `flail`), `adam_flail_parameters`, `adam_poisson_fft_object`,
`adam_poisson_mg_object` and `adam_lattice_green_library` are described on [FLAIL](./flail).

## Infrastructure

### `adam_mpih_object`

A thin layer over MPI: rank and size, tic-toc timing, rank-prefixed output, barrier, `abort` and `error_stop`. The
ghost exchange itself is in `adam_field_object` (`update_ghost_mpi`).

### `adam_io_object`

The `[IO]` section (output basenames and cadences, auxiliary-field switches, the seam div(B) monitor), the history
files (residuals, energy, divergence, charge conservation) and the XH5F field writer `save_field` (rank 4 and 5; kinds
`R8P`, `R4P`, `I8P`, `I4P`, `I2P`, `I1P`). Field output is one HDF5 file per rank, `<basename>-procNNNNNN.h5`, indexed
by `<basename>.xdmf` (`realm_object%open_file_xh5f`, `save_q_xh5f`).

### `adam_slices_object`

Planar slices of the 3-D field at given bounds and resolution (`emin`, `emax`, `nijk`), every `n_save` steps, section
`[slices]`. Output is raw binary written with MPI-IO (`adam_object%save_slice`, extension `.mat`), converted by ASCOT
or `scripts/slices_converter.py`.

### `adam_adam_object`

`adam_object` binds grid, tree, maps, field and I/O and drives their lifecycle.

| Group | Methods |
|-------|---------|
| Initialization | `initialize` (from the INI file) |
| AMR cycle | `adapt`, `amr_update`, `refine_uniform` |
| Communication | `make_comm_local_maps_ghost_bc`, `mpi_redistribute`, `mpi_gather_refinement_needed`, `blocks_reorder` |
| I/O | `save_restart_files`, `load_restart_files` (stream files: blocks per rank `-procNNNNNN.fbd`, tree nodes `.tnd`, time `.time`), `save_slice`, `save_vtk` |
| Geometry | `interpolate_at_point`, `prune` |

### `adam_parameters`

Named constants shared by every module: boundary type flags, AMR marker values, face enumeration codes
(`FEC_1_6_ARRAY`) and direction mappings (`FEC_TO_DELTA`, `DELTA_TO_FEC`).

## Forest orchestration (multi-realm)

A forest couples several realms in one synchronized time loop. The user view (manifest, α/β cadence) is on
[forest](/concepts/forest); this section is the contract an application implements.

### `adam_realm_object`

`realm_object` is the concrete parent of every application. It owns every library object as a component (`adam` with
`grid`, `tree`, `field`, `maps`; `io`, `amr`, `slices`, `weno`, `ib`, `rk`, `flail`, …) and declares the
`_forest`-suffixed type-bound procedures the forest dispatches on.

| TBP | Role | Default |
|---|---|---|
| `initialize_process_context(verbose)` | Process-global MPI (+ device) init, once, before the first realm (FNL overrides it) | `mpih%initialize` |
| `initialize_forest(filename, realms_number)` (non-overridable) | The initialization template: `initialize_app_forest`, then `restart_state_forest` or the initial state and regrids, then `post_initial_state_forest` (below) | the template |
| `initialize_app_forest(filename, realms_number)` | Template hook: parse, build the library and app objects, set up the backend | `error_stop` |
| `restart_state_forest()` | Template hook: the state from the restart files | refuses `[IO].(restart)` |
| `initial_amr_iterations_forest()`, `set_initial_state_forest()`, `amr_update_initial_forest()` | Template hooks: `n` passes of initial state + regrid, then the final initial state | `error_stop` |
| `rebalance_initial_forest()`, `reset_time_forest()` | Template hooks: `[amr].(balance_static)` on the initial grid; time and step counter to zero | `error_stop` |
| `post_initial_state_forest()` | Template hook: device push, ghosts, t0 diagnostics, history files, AMR lock | `error_stop` |
| `finalize_forest()` | Close I/O, release resources | `error_stop` |
| `finalize_mpi_forest()` | Process-global MPI finalize, once, after every realm's `finalize_forest` (FNL overrides it) | `mpih%finalize` |
| `compute_local_dt_forest(dt_local)` | This realm's CFL-limited dt; the forest min-reduces | `error_stop` |
| `advance_one_step_forest(dt)` | Fast path: the realm runs the whole step | `error_stop` |
| `open_step_forest(dt)` | Staged path: step prologue (external fields, stage buffers, time bookkeeping) | `error_stop` |
| `begin_stage_forest(k, K_total, dt, realm)` | Staged path: open stage `k`; sets `self%stage_active = k` | `error_stop` |
| `end_stage_forest(k, K_total, dt, realm, flux_register)` | Staged path: residuals and stage assignment of stage `k` | `error_stop` |
| `close_step_forest(dt)` | Staged path: step epilogue (state assembly, BC, cleaning, time advance); sets `self%stage_active = 0` | `error_stop` |
| `fill_seam_from_peer_forest(peer, p_idx)` | Copy the peer's interior into own seam ghosts; source buffer by `stage_active` (`> 0`: `rk%q_rk(:,…,stage)`, else committed `q`) | `error_stop` |
| `pack_seam_cells_forest(p_idx, buf)`, `unpack_seam_cells_forest(p_idx, buf)` | Pack own cells of the cross-rank send rows; unpack the receive rows into own ghosts | `error_stop` |
| `after_topology_build_forest()` | Hook after the seam maps are built (FNL refreshes the device maps) | nothing |
| `apply_reflux_to_stage_forest(stage, dt, flux_register)` | Berger–Colella reflux at the realm's end of step | `error_stop` |
| `stages_per_step_forest()` | K, the realm's stage count (`rk%nrk`) | `error_stop` |
| `coupling_descriptor_forest(scheme_time, rk_scheme, nv)` | `(scheme_time, rk_scheme, nv)` for β admissibility | `''`, `''`, `-1` (cannot join β) |
| `is_done_forest(done)` | Termination predicate | `error_stop` |
| `post_step_forest(dt, t, it, do_save_state, do_save_residuals, do_save_restart, do_amr, realm)` | Per-step diagnostics and I/O | `error_stop` |
| `runtime_amr_forest()`, `regrid_forest(regridded)` | Whether the realm regrids during the run; the regrid, called by `forest%regrid` after `post_step_forest` | false; never |

The forest receives `type(realm_pointer), intent(inout) :: realm(:)` and owns no realm state. The driver allocates the
realms as a concrete type (e.g. `type(prism_cpu_object), allocatable, target :: realm(:)`) and points one
`realm_pointer` at each, so realm `is` is `realm(is)%ptr`. Never pass a polymorphic array
`class(realm_object) :: realm(:)`: nvfortran (25.11, 26.9) can address its elements with the declared type's stride
and dispatch on the wrong realm ([#96](https://github.com/szaghi/adam/issues/96)).

The fast path skips every per-stage TBP. Only a single realm without 2:1 faces (`flux_register%nfaces == 0`) takes
it; a single realm with a 2:1 face runs the staged path. An app that refines must implement the staged TBPs; a
single-realm app can leave the seam TBPs at their defaults.

### `adam_forest_object`

`forest_object` holds the realm count `n`, the flux register (the coarse-fine faces of every realm and seam) and the
timing report. Every operation is a TBP dispatch on the realm array passed in.

| Method | Purpose |
|---|---|
| `initialize(realm, filename)` | Every realm reads the same INI (single-INI drivers) |
| `initialize_from_manifest(realm, manifest)` | Each realm reads its own INI from the parsed manifest; the topology is wired afterwards |
| `simulate(realm, filename)`, `simulate_from_manifest(realm, manifest)` | `initialize`, time loop, `finalize` |
| `compute_global_dt(realm, dt)` | Min-reduce the realms' dt within the process and across ranks |
| `evolve_one_step(realm)` | One global step: phase 0 (α seam fill, open), K stage loop (phases 1–3), phase 4; per-seam cadence gating |
| `post_step(realm)` | Diagnostics of every realm |
| `regrid(realm)` | Each realm's `regrid_forest`; rebuild the flux register when a tree changed |
| `is_done(realm, done)` | AND-reduce termination across realms and ranks |
| `finalize(realm)` | Each realm's `finalize_forest` |

| Private helper | Purpose |
|---|---|
| `populate_inter_realm_topology(realm, manifest)` | Manifest face pairs → `maps%inter_realm_neighbors` and the `seam_local_*` arrays |
| `check_beta_admissibility(realm, manifest)` | Same `scheme_time`, `rk_scheme`, `nv` and K on both ends of every `stage_coincident` seam; a mismatch is an `error_stop` |
| `apply_reflux_corrections(realm, stage, dt)` | Each realm's `apply_reflux_to_stage_forest` |
| `register_intra_realm_amr_seams(realm)` | Coarse-fine faces inside each realm into the flux register |

`ADAM_SEAM_TIMING=1` in the environment arms the timing report of the forest loop and of the seam exchange phases.

### Manifest data flow

```
forest.ini                          (text manifest)
  │
  ▼  adam_forest_manifest.F90 ─ read_forest_manifest
forest_manifest_t                   (realm_ini paths + face_pairs)
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
realm(is)%adam%maps%seam_local_send_buf / recv_buf        (per-peer buffers, allocated once)
```

No manifest is read at step time.

### Per-seam coupling cadence

| Constant | Value | Fill |
|----------|-------|------|
| `CADENCE_END_OF_STEP` (α, default) | 0 | Once per step, before `open_step_forest` (phase 0). `begin_stage_forest` copies them into each stage state (`seam_ghosts_to_state`) |
| `CADENCE_STAGE_COINCIDENT` (β) | 1 | Every stage, before `end_stage_forest` (phase 2) |

A 1:1 same-resolution mirror seam in PRISM must use β to stay div-free: PRISM's α stages read seam ghosts that were
not refreshed ([#116](https://github.com/szaghi/adam/issues/116)).

### `adam_seam_exchange`

The inter-realm seam fill across ranks. Realms are partitioned over the ranks independently, so the peer cell holding
a ghost centre can be on another rank. At topology time the forest splits the seam rows of every realm:

| Rows | Maps | Path |
|------|------|------|
| Local | `maps%seam_local_map_ghost_cell` | Copied by the realm's `fill_seam_from_peer_forest` |
| Cross-rank | `maps%seam_mpi_recv_cell` (ghost owner), `maps%seam_mpi_send_cell` (cell owner) | Enumerated in one canonical order on every rank, so send and receive rows match one to one |

Across a 2:1 seam (`coupling = refined`) a send row has a kind: `SEAM_ROW_COPY`, `SEAM_ROW_INTERPOLATE` (the
coarse-to-fine fill `interp_seam_ghost`) or `SEAM_ROW_RESTRICT` (the fine-to-coarse mean). The cell owner computes the
ghost values. Same-rank interpolate and restrict rows travel as messages to self, without MPI.

| Procedure | Role |
|-----------|------|
| `seam_fill(realm, p)` | One (realm, peer slot) pair: local copy, pack (`pack_seam_cells_forest`), point-to-point transfer, unpack |
| `seam_fill_all`, `seam_fill_batch(realm, cadence, wtime)` | Every pair, or every pair of one cadence in one batch |
| `seam_ghosts_to_state` | Copy the seam ghosts into the active stage state (α, #132 rule 1) |
| `pack_seam_rows`, `unpack_seam_rows` | Row packing shared by the realm TBPs |

Every rank calls `seam_fill` for the same pairs in the same order (the peer slots come from the manifest). A rank with
no rows of a pair returns at once; a rank that skips a pair it has rows for leaves its partners waiting.

### Intra-realm AMR seams

A 2:1 jump inside one realm's tree is a second, disjoint seam family. It is not declared in the manifest.

| Part | Implementation |
|------|----------------|
| Registration | `register_intra_realm_amr_seams` registers the coarse-fine faces in the flux register (`SEAM_KIND_INTRA_REALM_AMR`) and fills `maps%amr_seam_quadrant`, allocated for every realm with blocks |
| Fine ghosts | `field_object%update_ghost_local`, map rows of flag 4 → `interp_seam_ghost`, regime `maps%seam_ghost_fill` from `[amr] seam_ghost_fill` |
| Coarse ghosts | Regime `[amr] seam_restriction`: `mean`, `lagrange6`, `weno6` |
| Device path | `update_ghost_local_gpu(..., seam_ghost_fill=...)` |
| Reflux | Berger–Colella at each realm's end of step, the same register as the inter-realm seams |

| `seam_ghost_fill` | Constant | Order |
|-------------------|----------|-------|
| `injection` | `SEAM_FILL_INJECTION` (0) | copy |
| `restriction-compatible` | `SEAM_FILL_COMPATIBLE` (1) | q = 2 |
| `tricubic` (default) | `SEAM_FILL_TRICUBIC` (2) | q = 4 |
| `quintic` | `SEAM_FILL_QUINTIC` (3) | q = 6 |
| `weno6` | `SEAM_FILL_WENO6` (4) | nonlinear |
| `quintic-mp` | `SEAM_FILL_QUINTIC_MP` (5) | q = 6, bounded |

The constants are in `adam_seam_interpolation_library`. `adam_object%initialize` sets `tricubic` when the key is
absent, overriding the `maps` default (injection).

**Seam fill guard.** When an app sets `maps%seam_density_var` to its density index, with the three momenta at the
next indices (FLUME does; 0, the default, turns it off, as in PRISM), `seam_mp_bound` clamps the fill in the manner of
Suresh and Huynh. The density is bounded by the coarse cells bracketing the ghost, widened by the data's curvature, and
the velocity by the same bound on the coarse velocities. A ghost whose clamped density is below `SEAM_DENSITY_GUARD`
= 0.1 times its anchor's takes the anchor's value in every variable. The decision is the same in every row and in any
order: CPU, FNL kernels, MPI pack and inter-realm pack (`pack_seam_rows`). Fields other than the state pass
`seam_guard = .false.` to `update_ghost_local` / `update_ghost_mpi` (FLUME's positivity-limiter factors, whose zero is
legitimate).

| Measured cost of the guard | Value |
|----------------------------|-------|
| `shock-cylinder-ib` (the regression case with the most seam ghosts) | +20 % run time |
| Other regression cases | within noise |

The two families are independent: a single-realm run never calls `fill_seam_from_peer_forest`, a same-resolution
forest never calls `interp_seam_ghost`, and both coexist in a forest of refined realms.

### Seam div(B) guard-rail

A 2:1 seam breaks the identity `div_h(curl_h) = 0` and injects an O(h^p) div(B) source: it converges under refinement
but grows without bound in time at fixed h. No local collocated fix exists ([#29](https://github.com/szaghi/adam/issues/29)).
`save_divergence_history` (`adam_prism_common_object.F90`) carries an opt-in monitor:

| `[IO]` key | Default | Effect |
|------------|---------|--------|
| `seam_divB_tol` | `-1.0` (off) | Arms the monitor when `> 0`; seam `max|div(B)|` above it is reported |
| `seam_divB_error` | `.false.` | `.true.`: `error_stop` instead of a warning |

The keys live on `adam_io_object`. The monitor covers intra-realm 2:1 faces and inter-realm `refined` seams; without
2:1 seams div(B) stays at round-off and the monitor stays silent.

### Load-bearing invariants

| Invariant | Why |
|-----------|-----|
| Phase 2 completes on all realms before phase 3 starts on any (β) | Otherwise `fill_seam_from_peer_forest` writes race with `compute_residuals` reads |
| A div(B) diagnostic refills the inter-realm seam after `update_ghost` | `update_ghost` fills intra-realm ghosts and physical BCs, not the inter-realm seam; omitting the refill reports a spurious seam div(B). The refill is `present(realm)`-guarded |
| `begin_stage_forest(k)` sets `stage_active = k`; `close_step_forest` clears it | The seam fill picks `q` or `rk%q_rk(:,…,stage_active)` from it; a write to `q` outside this discipline corrupts peer reads |
| `apply_reflux_to_stage_forest` and `accumulate_seam_fluxes_fv` act at `stage == rk%nrk` only | The register's stage axis is collapsed to 1; per-stage reflux is not implemented. α and β share this gate |
| Fast path only for one realm with no register faces | A refining app goes through the staged TBPs |
| No library singletons besides `mpih` and `mpih_fnl` | Every object is a realm component passed as an argument |

## Module summary

| Module | Category | File |
|--------|----------|------|
| `adam_parameters` | Infrastructure | `adam_parameters.f90` |
| `adam_tree_node_object` | AMR | `adam_tree_node_object.f90` |
| `adam_tree_bucket_object` | AMR | `adam_tree_bucket_object.f90` |
| `adam_tree_object` | AMR | `adam_tree_object.F90` |
| `adam_refinement_plan_object` | AMR | `adam_refinement_plan_object.F90` |
| `adam_amr_object` | AMR | `adam_amr_object.F90` |
| `adam_amr_criterion_object` | AMR | `adam_amr_criterion_object.F90` |
| `adam_grid_object` | Grid/field | `adam_grid_object.F90` |
| `adam_field_object` | Grid/field | `adam_field_object.F90` |
| `adam_maps_object` | Communication | `adam_maps_object.F90` |
| `adam_weno_object` | Numerics, spatial | `adam_weno_object.F90` |
| `adam_fdv_operators_library` | Numerics, spatial | `adam_fdv_operators_library.F90` |
| `adam_riemann_euler_library` | Numerics, spatial (unused) | `adam_riemann_euler_library.F90` |
| `adam_stepper_object` | Numerics, temporal (SELAH adapter) | `adam_stepper_object.F90` |
| `adam_eos_ic_object` | Physics (unused) | `adam_eos_ic_object.F90` |
| `adam_ib_object` | Physics | `adam_ib_object.F90` |
| `adam_turbulence_object` | Physics | `adam_turbulence_object.F90` |
| `adam_turbulence_library` | Physics | `adam_turbulence_library.F90` (+ `adam_turbulence_body.INC`) |
| `adam_turbulence_kernels` | Physics | `adam_turbulence_kernels.F90` |
| `adam_flail_object` | Elliptic | `adam_flail_object.F90` |
| `adam_flail_parameters` | Elliptic | `adam_flail_parameters.F90` |
| `adam_poisson_fft_object` | Elliptic | `adam_poisson_fft_object.F90` |
| `adam_poisson_mg_object` | Elliptic | `adam_poisson_mg_object.F90` |
| `adam_lattice_green_library` | Elliptic | `adam_lattice_green_library.F90` |
| `adam_mpih_object` | Infrastructure | `adam_mpih_object.F90` |
| `adam_io_object` | Infrastructure | `adam_io_object.F90` |
| `adam_slices_object` | Infrastructure | `adam_slices_object.F90` |
| `adam_adam_object` | Infrastructure | `adam_adam_object.F90` |
| `adam_realm_object` | Forest | `adam_realm_object.F90` |
| `adam_forest_object` | Forest | `adam_forest_object.F90` |
| `adam_forest_manifest` | Forest | `adam_forest_manifest.F90` |
| `adam_flux_register_object` | Forest | `adam_flux_register_object.F90` |
| `adam_seam_exchange` | Forest | `adam_seam_exchange.F90` |
| `adam_seam_interpolation_library` | Forest | `adam_seam_interpolation_library.F90` |
| `adam_common_library` | Entry point | `adam_common_library.F90` |
| `adam_mpih_global` | Singleton | `adam_mpih_global.F90` |
| `adam_globals` | Singleton aggregator | `adam_globals.F90` |

`adam_common_library` does not re-export `adam_amr_criterion_object`, `adam_seam_interpolation_library`,
`adam_flail_parameters` (its `FLAIL_BC_*` kinds come through `adam_flail_object`), `adam_poisson_fft_object`,
`adam_poisson_mg_object` and `adam_lattice_green_library`: use them directly.
