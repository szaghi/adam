# Library input reference

The sections below are read by the ADAM library itself, so every application built on it (FLUME, PRISM) reads them
the same way. Each application adds its own sections:
[FLUME input reference](/applications/flume/input), [PRISM input reference](/applications/prism/input). Where the two
applications use a library key differently, the **Used** column says so.

## Parsing rules

The input is an INI file read by [FiNeR](https://github.com/szaghi/FiNeR):

- **Section and key names are case-sensitive**: `CFL`, not `cfl`; `[IO]`, not `[io]`. When a key or section is
  duplicated, the first occurrence wins.
- **Comments.** A line whose first non-blank character is `;`, `#` or `!` is a comment. After a value, only `;` starts
  an inline comment.
- **Missing key, or a value that does not convert.** FiNeR returns an error code and leaves the variable untouched,
  both when the key is absent and when its value does not convert to the key's type (a non-numeric value for a number,
  `finer_option_t.F90`). Most readers then stop the run (`failed to load [section].(key)`); the readers that ignore the
  error are marked **no check** below, and the variable keeps its declared default, or stays undefined when it has
  none. Integers are accepted for reals (`1` is a valid real).
- **Logical keys.** FiNeR reads the value when its first character is one of `.tTfF` (`.true.`, `T`, `true`,
  `.false.`, ...); any other value (`yes`, `1`) is a conversion error, as above.
- **CRLF files** (Windows line endings) read as LF ones on every compiler: FiNeR drops every CR before splitting the
  file into lines (FiNeR `c55f8d3`). Before that fix, nvfortran kept the CR at the end of every value and every FNL
  build stopped on a CRLF input at its first numeric key; a dependency fetched before 2026-10-26 needs
  `fobis fetch --update`.
- **Removed keys.** FiNeR ignores keys nobody reads, so a removed key is not an error by itself. The lint
  `src/tests/lint/check-removed-input-keys.sh` (run at the top of the PRISM and FLUME regression sweeps) rejects, in
  the committed PRISM, FLUME and AMR test inputs: `[fWLayer] C` (replaced by `width`), `[PIC] elliptic_correction`
  (removed in [#84](https://github.com/szaghi/adam/issues/84) E3d; `.true.` stops the run) and the whole
  `[linear-algebra]` section (removed in #84 E4 with the legacy elliptic smoothers it configured).

## Legend for the tables

- **Req.**: `yes` means a missing key stops the run. `cond.` means it is required only under the stated condition.
  `no` means optional (the default applies). `no check` means the error is ignored; see the notes for what happens.
- **Used**: who depends on the value. `both`, `FLUME`, `PRISM`, or `neither` (the key is read, so it must usually be
  present, but no application uses it).

---

## [IO]

`io_object%load_from_file`. Each application adds keys of its own to this section: FLUME `save_auxiliary_fields`,
PRISM `save_fields`.

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `output_basename` | string | yes | — | any path prefix | both | Prefix of every output file (fields, histories, slices). |
| `it_save` | int | yes | — | `<= 0` disables field saves | both | Field (XH5F) save cadence in steps. |
| `restart` | logical | yes | — | logical | both | `.true.`: skip the initial conditions and load the `restart_basename` restart files. |
| `restart_basename` | string | yes (even when `restart = .false.`) | — | path prefix | both | Basename of the restart files, read and written. |
| `restart_save` | int | yes | — | `<= 0` disables | both | Restart save cadence in steps. |
| `residuals_save` | int | yes | — | `<= 0` disables | both | Residual L2-norm history cadence (`<output_basename>-residuals.dat`). |
| `divergence_history_save` | int | no | 10 | no check | neither | Read, never used: PRISM writes its divergence history at every diagnostic step, FLUME has its own cadence (`[diagnostics]`). |
| `save_memory_status` | logical | yes | — | logical | PRISM | Print the memory status at the diagnostic steps. |
| `save_residual_fields` | logical | yes | — | logical | both | Also write the residuals `dq_*` into the field files. |
| `save_curl_fields` | logical | yes | — | logical | PRISM | Also write the curl fields. |
| `save_divergence_fields` | logical | yes | — | logical | PRISM | Also write the divergence fields. |
| `save_gradient_fields` | logical | yes | — | logical | neither | Read, never used. |
| `save_laplacian_fields` | logical | yes | — | logical | neither | Read, never used. |
| `seam_divB_tol` | real | no | -1 (off) | no check | PRISM | 2:1 seam div(B) monitor: warn (or stop, below) when the seam `max\|div B\|` exceeds it on an AMR run ([#29](https://github.com/szaghi/adam/issues/29)). FLUME uses `[mhd] divb_tol`. |
| `seam_divB_error` | logical | no | `.false.` | no check | PRISM | `.true.`: the seam monitor stops the run instead of warning. |

**Output files.** The library writes the fields as one `<output_basename>-<it:9 digits>-procNNNNNN.h5` per rank plus one
`.xdmf` (XH5F, ghost cells included), the restart files as `<restart_basename>.time`, `.tnd` and
`-procNNNNNN.fbd`, and the slices as raw MPI-IO `.mat` files. There is no output-format key.

---

## [grid]

`grid_object%load_from_ini_file`. Every key is required.

| Key | Type | Req. | Accepted / checks | Used | Meaning |
|-----|------|------|-------------------|------|---------|
| `ni`, `nj`, `nk` | int | yes | no check of `>= 1`. The applications add parity and width rules for AMR (see their pages). | both | Cells per block along i, j, k, ghost cells excluded. |
| `ngc` | int | yes | each application checks it against its stencils: `>= [weno] S`, `>= fdv_order/2`, and more (see their pages) | both | Ghost cells per side. |
| `emin_x`, `emin_y`, `emin_z` | real | yes | no check | both | Domain minimum corner. Divided by `L0` in a PRISM non-dimensional model. |
| `emax_x`, `emax_y`, `emax_z` | real | yes | no check (`emax > emin` is not validated) | both | Domain maximum corner. |
| `null_x`, `null_y`, `null_z` | logical | yes | logical | both | Null (inactive) direction: excluded from the fluxes and the time step. |

---

## [amr]

Read in four places: the grid reads `ratio`, the tree the refinement shape, `adam_object` the seam fills, `amr_object`
the cadence and the markers.

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `ratio` | int | **no check** (effectively required) | declared 8, **not applied** | `2` (binary), `4` (quadtree), `8` (octree). Any other value leaves the block sizes undefined, with no error. | both | Refinement ratio of the tree. |
| `max_level` | int | **no check** | declared 12, **not applied** | Refinements beyond it are cancelled silently. | both | Maximum refinement level. |
| `iu_ref_levels` | int | **no check** | declared -1, **not applied** | `<= 0`: none | both | Uniform refinement levels applied at initialisation, before the initial conditions. The blocks are spread over the ranks by whole sibling families (`ratio` blocks each), so a grid with fewer families than ranks leaves ranks with no block; supported, but the empty ranks do no work. |
| `i_prune`, `j_prune`, `k_prune`, `l_prune` | int | **no check** | declared -1, **not applied** | see the notes below | both | Pruning of a "simple initial forest". No application calls `prune`, but the values enter the neighbour/boundary search: if all four are `>= 0` the tree sizes its search box `2**(l-l_prune)*(ijk_prune+1)`. **Keep at least one negative** (all four `-1` is safest; `0, 0, 0, -1` is equivalent). |
| `frequency` | int | yes | 100 | `>= 0` | FLUME | Runtime regrid cadence ([#74](https://github.com/szaghi/adam/issues/74)): `0` no runtime regridding, `n > 0` regrid after every `n`-th step. PRISM reads it and refines at initialization only. |
| `iters` | int | yes | 5 | any | both | Maximum marker sweeps per AMR update (it stops early when the grid is stable); a sweep refines or coarsens one level. |
| `markers_number` | int | yes | 0 | `>= 0` | both | Number of `[amr_marker_N]` sections. |
| `regrid_prolongation` | string | no | `linear` (FLUME: `conservative`) | `linear`, `conservative`, `conservative-state`; any other value is fatal | FLUME | How a refined block's children are filled from the parent at a runtime regrid; FLUME's [input reference](/applications/flume/input#setting-up-runtime-amr) describes the three. |
| `seam_ghost_fill` | string | no | `tricubic` | `injection`, `restriction-compatible`, `tricubic`, `quintic`, `weno6`, `quintic-mp`; any other value is fatal. The 6-point fills (`quintic`, `weno6`, `quintic-mp`) need at least 6 cells per refined axis in every block (fatal otherwise); on a `coupling = refined` forest seam both realms must declare the same value. | both | Coarse-to-fine ghost fill at 2:1 AMR seams ([#68](https://github.com/szaghi/adam/issues/68) D1, D2). `tricubic`: the tensor 4-point Lagrange interpolation, exact for per-direction degree 3. `quintic`: the 6-point one, exact for degree 5. `weno6`: a WENO-Z combination of the cubics of the 4-point sub-stencils, dropping those that hold a jump. `quintic-mp`: the quintic clamped to a monotonicity-preserving bound. With every regime, an app that sets `maps%seam_density_var` (FLUME) has the density and the velocity of a fine ghost clamped to a monotonicity-preserving bound of the coarse cells around it, and a ghost whose clamped density is below 0.1 times its coarse cell's take that cell's state ([#101](https://github.com/szaghi/adam/issues/101)). |
| `seam_restriction` | string | no | `mean` | `mean`, `lagrange6`, `weno6`; any other value is fatal. `lagrange6` and `weno6` need at least 6 cells per refined axis in every block (fatal otherwise); on a `coupling = refined` forest seam both realms must declare the same value. | both | Fine-to-coarse ghost fill at 2:1 AMR seams ([#68](https://github.com/szaghi/adam/issues/68) C1–C3). `mean`: the mean of the fine cells, second order. `lagrange6`: the tensor 6-point Lagrange interpolation to the coarse centre, exact for degree 5. `weno6`: a WENO-Z combination of the 4-point cubics, falling back to the mean across a jump. |

**Trap:** the tree reader passes one buffer to every key and **ignores the error**. A missing key therefore takes the
value of the key read just before it: a missing `max_level` becomes equal to `ratio`, and a missing `ratio` is garbage.
The declared defaults never take effect. **Always write all seven tree keys.**

**When the markers act.** At initialization, during the `[initial_conditions] amr_iterations` passes (each pass sets the
initial condition and runs one AMR update of up to `iters` sweeps), in both applications; in FLUME also at the runtime
regrids (`frequency > 0`).

### [amr_marker_N] (N = 1..markers_number)

`amr_object%load_from_file` reads the keys a `mode` needs; the applications implement different modes.

| Mode | FLUME | PRISM |
|------|-------|-------|
| `1` geometric (`AMR_GEO`) | `solid` and `primitive-box`; `stl` fatal | `primitive-box`; `solid` does nothing; `stl` fatal |
| `2` gradient (`AMR_GRAD`) | yes | **does nothing** (the branch is empty) |
| `3` total variation (`AMR_TV`) | fatal at initialisation ("not implemented in FLUME") | yes: the total variation of the coils' current density |
| `4` Löhner (`AMR_LOHNER`) | yes, at runtime regrids only (it reads ghost cells the initial AMR does not fill) | **does nothing** (no branch) |
| `5` Jeans length (`AMR_JEANS`) | yes, needs `[gravity] self_gravity = .true.` | **does nothing** (no branch) |
| `6` registered criterion (`AMR_USER`) | yes, at runtime regrids only: the criterion the program registered under `name` ([#110](https://github.com/szaghi/adam/issues/110)); an unregistered name is fatal at the first regrid | **does nothing** (no branch) |

An unknown mode reads no further keys and stops nothing in the library; FLUME stops on it at the first AMR pass.

| Key | Type | Req. | Accepted / invalid | Modes | Meaning |
|-----|------|------|--------------------|-------|---------|
| `mode` | int | yes | see the table above | all | Marker kind. |
| `delta_type` | string | yes (all modes) | `x`, `y`, `z`, `max` (the maximum over the active directions); other values are fatal at marking in FLUME, and in PRISM's mode-3 marker leave the compared spacing undefined (no check) | 1 (solid), 2, 3, 5 | Which block spacing is compared with `delta_fine`/`delta_coarse`. |
| `delta_fine` | real | yes (all modes) | no check | 1 (solid), 2, 3 | Admissible spacing where the criterion fires. |
| `delta_coarse` | real | yes (all modes) | no check | 1 (solid), 2, 3 | Admissible spacing elsewhere. |
| `geo_type` | string | mode 1 | `solid`, `primitive-box`, `stl`; other values are fatal at load | 1 | Geometric marker kind. |
| `solid` | int | mode 1 + `solid` | `1..[solids] number`, otherwise fatal at marking | 1 | Refine the blocks crossed by that solid's surface. |
| `box_xmin`, `box_ymin`, `box_zmin`, `box_xmax`, `box_ymax`, `box_zmax` | real | mode 1 + `primitive-box` | no check | 1 | Box: the blocks whose centroid lies inside it are refined. |
| `target_level` | int | mode 1 + `primitive-box` | no check (capped by `max_level`) | 1 | Refine the in-box blocks until they reach this level. |
| `stl_filename` | string | mode 1 + `stl` | — | 1 | STL surface (both applications refuse `stl`). |
| `field` | int | modes 2, 3, 4 | `1` = conservative state, `2` = auxiliary fields; other values are fatal at marking | 2, 4 | Array holding the marker variable (PRISM's mode 3 reads it, unused). |
| `var` | int | modes 2, 3, 4 | a variable index of that array, otherwise fatal at marking | 2, 4 | Marker variable. |
| `tol` | real | modes 2, 3 | no check | 2, 3 | Threshold: the gradient magnitude (mode 2) or the total variation (mode 3) above which the block gets `delta_fine`. |
| `refine_tol`, `derefine_tol` | real | mode 4 | `0 <= derefine_tol < refine_tol <= 1`, otherwise fatal | 4 | Löhner thresholds: refine above `refine_tol`, coarsen below `derefine_tol`, keep between. |
| `epsilon` | real | no (mode 4) | `>= 0`, default 0.01 | 4 | Noise filter relative to the variable. |
| `floor` | real | no (mode 4) | `>= 0`, default 0 | 4 | Absolute noise filter, for a variable that vanishes in quiet regions. |
| `buffer` | int | no (mode 4) | `>= 0`, default 1 (FLUME: also `<= ngc - 1`) | 4 | Ghost layers the estimator reads. |
| `jeans_cells` | real | no (mode 5) | `> 0`, default 8 | 5 | Refine where some interior cell resolves its Jeans length with fewer cells. |
| `name` | string | mode 6 | the name a criterion was registered under; empty is fatal at load, unregistered at initialisation | 6 | The registered criterion (`adam_amr_criterion_object`). |
| `jeans_cells_derefine` | real | no (mode 5) | `> jeans_cells`, default `2 * jeans_cells`; otherwise fatal | 5 | Coarsen where every interior cell has more cells per Jeans length; keep between. |

FLUME's page lists its variable indices and the [runtime AMR setup](/applications/flume/input#setting-up-runtime-amr).

---

## [field]

| Key | Type | Req. | Used | Meaning |
|-----|------|------|------|---------|
| `nv` | int | **not read** | neither | Both applications pass `nv` from their physics, so `field_object` never reads this key; it can be omitted. |

---

## [runge_kutta]

`rk_object%load_from_file`; each application also reads `scheme` first, to size the per-block memory budget.

| Key | Type | Req. | Accepted / invalid | Used | Meaning |
|-----|------|------|--------------------|------|---------|
| `scheme` | string | yes | see the table below; an unknown name is fatal | both | Time integrator. It also sizes the per-block memory budget (stage fields). |

| Name | Kind | Stored stage fields | FLUME | PRISM |
|------|------|---------------------|-------|-------|
| `runge-kutta-1` | low-storage, forward Euler | 1 | accepted | EM: accepted; PIC, PML, 2:1 faces: refused |
| `runge-kutta-2` | low-storage TVD, 2nd order | 1 | accepted | as above |
| `runge-kutta-3` | low-storage TVD, 3rd order | 1 | accepted | as above |
| `runge-kutta-ssp-11` | SSP, forward Euler | 1 | accepted | EM and PML: accepted; PIC: refused |
| `runge-kutta-ssp-22` | SSP, 2nd order | 2 | accepted | accepted |
| `runge-kutta-ssp-33` | SSP, 3rd order | 3 | accepted | accepted |
| `runge-kutta-ssp-54` | SSP, 5 stages, 4th order | 5 | accepted | accepted |
| `runge-kutta-yoshida` | symplectic, 4th order, 4 stages | 0 | **refused** | EM on CPU: accepted; FNL, PIC, PML, 2:1 faces: refused |

The staged path of the forest (any 2:1 AMR face, any multi-realm forest) needs an SSP scheme in both applications. In
a forest manifest, seams with `coupling_cadence = stage_coincident` require both realms to use the same `scheme`
(fatal otherwise).

---

## [weno]

`weno_object%load_from_file`. FLUME reconstructs with it always; PRISM only with `[numerics] scheme_space = weno`, a CPU
path (the FNL backend refuses it).

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `scheme` | string | yes | — | see the table below; an unknown name is fatal | both | WENO order (stencil half-width `S`). `ngc >= S` is checked by FLUME always, by PRISM with `scheme_space = weno`. |
| `weights` | string | no | `js` | `js` (Jiang–Shu, absolute ε), `si` (scale-invariant); an unknown value is fatal | both | Nonlinear weights (PRISM uses them with `scheme_space = weno`). |
| `ror_number` | int | yes | 0 | `>= 0` | neither | Number of order-reduction (ROR) stages. `S` always comes from `scheme`: `initialize` sets it after reading the ROR keys. |
| `ror_scheme_1` … `ror_scheme_<ror_number>` | int | cond. (`ror_number > 0`) | — | — | neither | ROR stencil sizes. |
| `ror_threshold` | real | yes | 0.9 | — | neither | ROR trigger. |
| `ror_vars_number` | int | yes | 2 | `>= 0` | neither | Number of ROR check variables. |
| `ror_ivar_1` … `ror_ivar_<ror_vars_number>` | int | cond. (`ror_vars_number > 0`) | — | — | neither | ROR check variables. |
| `enable_ror_stats` | logical | yes | `.false.` | logical | neither (`.true.` allocates an `nb x cells x 3` integer array) | ROR statistics. |
| `ib_reduction_extent` | int | yes | 0 | — | neither | Order reduction near immersed solids. |
| `ib_reduced_order` | int | yes | 1 | — | neither | Reduced `S` near solids. |

| Name | S | Order | Kind |
|------|---|-------|------|
| `weno-u-1`, `weno-u-3`, `weno-u-5`, `weno-u-7`, `weno-u-9` | 1–5 | 1, 3, 5, 7, 9 | upwind |
| `weno-c-2`, `weno-c-4`, `weno-c-6`, `weno-c-8` | 1–4 | 2, 4, 6, 8 | centred (FLUME refuses them) |

---

## [fdv]

`realm_object` reads it: the finite-difference/finite-volume operators (gradient, divergence, curl, Laplacian) of
`adam_fdv_operators_library`.

| Key | Type | Req. | Accepted / invalid | Used | Meaning |
|-----|------|------|--------------------|------|---------|
| `fdv_scheme` | string | yes | `FD`/`fd`/`Fd`/`fD` or `FV`/`fv`/`Fv`/`fV`; other values are fatal | PRISM | Operator family. PRISM's curl, divergence and coil stamping follow it; FLUME's div(B) diagnostic is always centred FD. |
| `fdv_order` | int | yes | even, 2 to 10 (half-stencil `fdv_order/2 <= 5`); no range check at load. PRISM checks `fdv_order/2 <= ngc` at initialization; FLUME checks it (and `<= 3` on FNL) at the first div(B) diagnostic of an MHD run. | both | Order of the centred operators. PRISM: its whole spatial operator and the t0 elliptic solves (which accept 2, 4, 6, 8: an order 10 stops at the first solve). FLUME: the div(B) diagnostic of MHD runs only. |

---

## [solids] and [solid_N]

`ib_object`: analytical immersed solids. FLUME applies them (Euler only); PRISM loads them and passes their distance
function to its Runge-Kutta stages, and no committed PRISM input uses them.

### [solids]

| Key | Type | Req. | Default | Meaning |
|-----|------|------|---------|---------|
| `number` | int | yes | 0 | Number of `[solid_N]` sections. `0`: no immersed boundary. |
| `n_eikonal` | int | yes | 2 | Eikonal extrapolation sweeps per residual evaluation, filling the solid cells. |

### [solid_N] (N = 1..number)

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `name` | string | yes | any | Label (printed only). |
| `bc_type` | int | yes | `1` viscous (all momentum negated inside the solid), `2` inviscid (normal momentum reflected); **any other value applies no inversion, silently** | Wall treatment of the solid-cell state. |
| `definition` | string | yes | `analytical_sphere`, `analytical_circle`, `analytical_rectangle`; any other value (`file.off` included, never implemented) is fatal at load | Geometry kind. |
| `sphere_center_x`, `sphere_center_y`, `sphere_center_z`, `sphere_radius` | real | `analytical_sphere` | no check | Sphere. |
| `circle_center_x`, `circle_center_y`, `circle_center_z`, `circle_radius` | real | `analytical_circle` | no check | Infinite cylinder. |
| `circle_axis` | char | `analytical_circle` | `x`, `y`, `z`; other values are fatal at load | Cylinder axis. |
| `rectangle_center_x`, `rectangle_center_y`, `rectangle_center_z`, `rectangle_edge_1`, `rectangle_edge_2` | real | `analytical_rectangle` | no check | Infinite rectangular prism. |
| `rectangle_axis` | char | `analytical_rectangle` | `x`, `y`, `z`; other values are fatal at load | Prism axis. |

---

## [slices] and [slice_N]

`slices_object`: the state sampled on a lattice of points and written as raw MPI-IO `.mat` files,
`<output_basename>-slice_NN-<it:9>.mat`. Both applications write them.

### [slices]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `slices_number` | int | yes | `>= 0`, negative is fatal | Number of `[slice_N]` sections. `0`: no slices. |

### [slice_N] (N = 1..slices_number)

Every key is required (fatal when missing) and checked at initialization.

| Key | Type | Accepted / invalid | Meaning |
|-----|------|--------------------|---------|
| `itype` | string | `trilinear`, `inverse_distance`; other values are fatal | Interpolation from the cells to the slice points. |
| `n_save` | int | `>= 1`; `0` or negative is fatal | Save cadence in steps; the last step is also saved. |
| `ni`, `nj`, `nk` | int | `>= 1`, otherwise fatal | Points of the lattice, cell-centred on the box. |
| `emin_x`, `emin_y`, `emin_z`, `emax_x`, `emax_y`, `emax_z` | real | `emin <= emax`, and every point inside the realm's domain, otherwise fatal | Box sampled by the lattice, in the grid's units (divided by `L0` in a PRISM non-dimensional model). |

---

## Forest manifest (multi-realm)

`adam_forest_manifest` reads it, for both applications: a separate INI whose `[realm.N]` entries point at full realm
inputs (each with the sections above and the application's own). If the file passed to the executable has a `[forest]`
section with `realms_number >= 1`, it is a manifest; otherwise it is itself a single-realm input. Manifest errors use
Fortran `error stop`, not the library's `error_stop`.

| Section | Key | Type | Req. | Accepted / invalid | Meaning |
|---------|-----|------|------|--------------------|---------|
| `[forest]` | `realms_number` | int | yes | `>= 1`, otherwise fatal. Its presence with `>= 1` makes the file a manifest. | Number of realms. |
| `[realm.N]` (N = 1..realms_number) | `ini` | string (≤ 256 chars) | yes | used verbatim: it resolves against the process working directory, **not** the manifest's directory | Realm input. |
| `[forest.topology]` | `inter_realm_faces_number` | int | no | default 0 (no coupling) | Number of glued faces. |
| `[forest.topology.face_N]` | `realm_a`, `realm_b` | int | yes | `1..realms_number`, otherwise fatal | Realms on each side of the seam. |
| | `face_a`, `face_b` | string | yes | `+x -x +y -y +z -z`, also upper case and the suffix forms `x+`, `X-`, …; other values are fatal | Glued faces. |
| | `coupling` | string | no | `mirror` (default), `refined`, `auto`, `periodic`, `interpolate` (lower, Capitalised or UPPER case); other values are fatal. `periodic` and `interpolate` are reserved: the forest stops at setup on them. | `mirror` joins seam cells of the same size (fatal otherwise). `refined` joins a 2:1 jump ([#52](https://github.com/szaghi/adam/issues/52)): the seam cells must differ by exactly 2 along every axis the trees refine and be equal along a null one (a quadtree's $z$), both realms of the same tree kind (`[amr] ratio`), the seam normal to a refined axis ([#71](https://github.com/szaghi/adam/issues/71)). `auto` ([#110](https://github.com/szaghi/adam/issues/110)) is decided at every topology build: either of the two when each side holds one cell size, else face by face (mixed levels along the seam: each face joins equal cells or a 2:1 jump, either realm the fine one), so runtime AMR may change it; it requires, whatever side is fine, the same block cell counts along the seam, the same `seam_ghost_fill` and `seam_restriction`, and even block cell counts along the refined axes in both realms (a pair that turns 2:1 also needs the `refined` conditions). With runtime AMR the forest keeps the seam cells equal (`mirror`), at the declared 2:1 jump (`refined`), one level per side for both, or, for `auto`, face by face equal or in a ratio of 2 between neighbouring leaves (edges and corners included): cell sizes, not levels, since each realm counts levels from its own root. |
| | `coupling_cadence` | string | no | `end_of_step` (default, α) or `stage_coincident` (β), in the same three case forms; other values are fatal. β requires both realms to have equal `[runge_kutta] scheme`, equal `nv` and equal stages per step (fatal otherwise). Conflicting cadences between the same realm pair are fatal. | Seam fill once per step (α) or at every Runge-Kutta stage (β). A 1:1 same-resolution seam needs β to stay divergence-free ([#31](https://github.com/szaghi/adam/issues/31)). |

The [forest concepts](/concepts/forest) and the [forest cookbook](/cookbook/forest) describe the multi-realm machinery.
