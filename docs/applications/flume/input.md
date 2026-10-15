# Input reference

A FLUME run reads one INI file: `mpirun -np N exe/adam_flume_cpu input.ini` (or `adam_flume_fnl`). This page lists
every section and key the application reads, with its type, whether it is required, its default, the accepted values
and what happens on an invalid one. The two complete inputs at the end are verification cases, copied verbatim; the
[verification gallery](./verification) shows their results.

## Program entry and read order

`adam_flume_cpu [file]` / `adam_flume_fnl [file]` (default `input.ini`). If the file has a `[forest]` section with
`realms_number >= 1`, it is a **forest manifest** and every `[realm.N]` names a per-realm INI of the form documented
below; otherwise the file is itself a single-realm input.

A realm input is read in this order: `[IO]` → `[reference]` (if present it converts the file to code units, then
`[IO]` is read again) → `[numerics]` → `[physics]`
(+ `[mhd]`, `[gravity]`) → `[runge_kutta].scheme` (memory budget) → `[grid]` → `[amr]` (tree, seam fill, markers) → `[solids]`
→ `[slices]` → `[runge_kutta]` → `[weno]` → `[linear-algebra]` → `[fdv]` → `[bc_*]` → `[time]` →
`[initial_conditions]` → `[diagnostics]` → `[IO].save_auxiliary_fields`. After reading, FLUME applies
cross-section checks (slices, WENO scheme, Riemann solver per model, ghost cells, even block cells for AMR).

## Units

Any consistent unit system works: the equations carry no dimensionless number
([units and scaling](./models#units-and-scaling)), so SI, cgs and units of order 1 give the same run. Magnetic fields
are always $\mathbf{B}_{SI}/\sqrt{\mu_0}$ (a Gaussian field divided by $\sqrt{4\pi}$). With the classic WENO
weights the data should be of order 1: either write the input in such units, add a
[`[reference]`](#reference-optional-dimensional-input) section that converts a physical input, or set
`[weno] weights = si`.

## Parsing rules

- **Section and key names are case-sensitive**: it is `CFL`,
  not `cfl`; `[IO]`, not `[io]`. When a key or section is duplicated, the first occurrence wins.
- Full-line comments start with `;`, `#` or `!`. A `;` also starts an inline comment after a value
 .
- **Missing key**: FiNeR returns `error = 3` and leaves the variable untouched. Whether that is fatal depends on the
  reader: most readers `error_stop` ("failed to load [sec].(key)"). The readers that ignore the error are marked
  **no check** below. Some of those have no default, so the value is left undefined.
- **Non-numeric value for a numeric key**: StringiFor's `to_number` returns an **uninitialised** result with
  `error = 0`. The parser does not detect this; the run continues with garbage.
  Integers are accepted for reals (`1` is a valid real).
- **Logical keys** are parsed with list-directed `read`: `.true.`, `.false.`, `T`, `F`, `true`, `false` work; any other
  token aborts with a Fortran runtime I/O error.
- **CRLF files**: FLUME's own string keys go through `strip_control` and are CRLF-safe. Library
  string keys (`[weno].scheme`, `[runge_kutta].scheme`, `[fdv].fdv_scheme`, `[amr_marker_N].geo_type`,
  `[linear-algebra].smoothing`, ...) are not; neither FiNeR nor StringiFor strips `\r` (not verified at runtime).
  **Save inputs with LF line endings.**

## Legend for the tables

- **Req.**: `yes` means a missing key causes `error_stop`. `cond.` means it is required only under the stated condition.
  `no` means optional (the default is applied). `no check` means the error is ignored; see the notes for what happens.
- **Used**: `yes` means FLUME's behaviour depends on the value. `no` means the key is read but FLUME never uses it;
  the key usually must still be present, because the library reader requires it.

---

## [IO]

Library `io_object%load_from_file`, plus one FLUME key. `[IO]` is loaded twice (FLUME common, then
the realm), and both loads are identical.

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `output_basename` | string | yes | — | any path prefix | yes | Prefix of every output file (fields, histories, slices). |
| `it_save` | int | yes | (100 unused) | `<= 0` disables field saves | yes | Field (XH5F) save cadence in steps. The first step (`it = 0`, not on restart) and the last step are always saved when `> 0`. |
| `restart` | logical | yes | — | logical | yes | `.true.`: skip the initial conditions and load `<restart_basename>` restart files. |
| `restart_basename` | string | yes (even when `restart = .false.`) | — | path prefix | yes | Basename used to load the restart files and to write them. |
| `restart_save` | int | yes | — | `<= 0` disables | yes | Restart save cadence in steps. The last step is also saved. |
| `residuals_save` | int | yes | — | `<= 0` disables | yes | Residual L2-norm history cadence (`<output_basename>-residuals.dat`). |
| `divergence_history_save` | int | no | 10 | no check | **no** | PRISM divergence history cadence. FLUME never calls `save_divergence_history`; its div(B) history uses `[diagnostics]`. |
| `save_memory_status` | logical | yes | — | logical | **no** | Not used anywhere in the code base. |
| `save_residual_fields` | logical | yes | — | logical | yes | Also write `dq_*` (residuals) into the XH5F field files. |
| `save_curl_fields` | logical | yes | — | logical | **no** | PRISM only. |
| `save_divergence_fields` | logical | yes | — | logical | **no** | PRISM only. |
| `save_gradient_fields` | logical | yes | — | logical | **no** | PRISM only. |
| `save_laplacian_fields` | logical | yes | — | logical | **no** | PRISM only. |
| `seam_divB_tol` | real | no | -1 (off) | no check | **no** | PRISM seam div(B) guard-rail. For FLUME use `[mhd].divb_tol`. |
| `seam_divB_error` | logical | no | `.false.` | no check | **no** | PRISM only. |
| `save_auxiliary_fields` | logical | yes | — | logical | yes | **FLUME key.** Also write the auxiliary fields (`rho,u,v,w,p,T,H,a`, plus `Bx,By,Bz` for MHD). For MHD it also writes the derived fields `pt, beta, bmag, divb`. |

### Output formats (no format key exists)

FLUME has no output-format selector. It always writes these files:

| Output | File name | Controlled by |
|--------|-----------|---------------|
| Fields: XDMF + HDF5 (XH5F), ghost cells included | `<output_basename>-<it:9 digits>.xdmf` + `...-procNNNNNN.h5` | `it_save`, `save_residual_fields`, `save_auxiliary_fields` |
| Restart | `<restart_basename>.time`, `.tnd`, `-procNNNNNN.fbd` + an XH5F snapshot `<restart_basename>.xdmf/.h5` | `restart_save` |
| Residuals history (ASCII, Tecplot header) | `<output_basename>-residuals.dat` | `residuals_save` |
| Conservation history (ASCII) | `<output_basename>-conservation_history.dat` | `[diagnostics].conservation_history_save` |
| div(B) history (MHD only) | `<output_basename>-divb_history.dat` | `[diagnostics].conservation_history_save` |
| Slices (MPI-IO binary) | `<output_basename>-slice_NN-<it:9>.mat` | `[slices]`, `[slice_N]` |

---

## [grid]

Library `grid_object%load_from_ini_file`. Every key is required.

| Key | Type | Req. | Accepted / checks | Used | Meaning |
|-----|------|------|-------------------|------|---------|
| `ni`, `nj`, `nk` | int | yes | See the notes below the table. | yes | Cells per block along i, j, k (ghost cells excluded). |
| `ngc` | int | yes | `>= [weno] S` (fatal). `>= 2` with `weno-riemann` (fatal). MHD: `>= fdv_order/2` (CPU) or `>= fdv_order/2` with `fdv_order/2 <= 3` (FNL), checked at the first div(B) diagnostic. | yes | Ghost cells per side. |
| `emin_x`, `emin_y`, `emin_z` | real | yes | no check | yes | Domain minimum corner. |
| `emax_x`, `emax_y`, `emax_z` | real | yes | no check (`emax > emin` is not validated) | yes | Domain maximum corner. |
| `null_x`, `null_y`, `null_z` | logical | yes | logical | yes | Null (inactive) direction: excluded from the fluxes and the time step. For Euler it also freezes the momentum along that direction. |

Notes on `ni`, `nj`, `nk`:

- There is no check that the value is `>= 1`.
- With `[initial_conditions].amr_iterations > 0`, the value must be **even** along every non-null axis (fatal). See.

---

## [amr]

Read in three places: the grid reads `ratio`, the tree reads the refinement shape,
`adam_object` reads `seam_ghost_fill` and `seam_restriction`, and `amr_object` reads the markers.

| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `ratio` | int | **no check** (effectively required) | declared 8, **not applied** | `2` (binary), `4` (quadtree), `8` (octree). Any other value leaves the block sizes undefined, with no error. `4` with `markers_number > 0` runs on CPU and FNL (issue #46; fatal on the NVF and GMP backends). | yes | Refinement ratio of the tree. |
| `max_level` | int | **no check** | declared 12, **not applied** | Refinements beyond it are silently cancelled. | yes | Maximum refinement level. |
| `iu_ref_levels` | int | **no check** | declared -1, **not applied** | `<= 0`: none | yes | Uniform refinement levels applied at initialisation, before the initial conditions. The blocks are spread over the ranks by whole sibling families (`ratio` blocks each, never split), so a grid with fewer families than ranks leaves ranks with no block: `iu_ref_levels = 1` on a quadtree is one family, all on rank 0. Supported on every backend ([DV-1](./verification#ranks-without-blocks-dv-1), [#82](https://github.com/szaghi/adam/issues/82)), but the empty ranks do no work. |
| `i_prune`, `j_prune`, `k_prune`, `l_prune` | int | **no check** | declared -1, **not applied** | See the notes below the table. | partly | Pruning of a "simple initial forest". FLUME never calls `prune`, but the values still enter the neighbour/boundary detection. |
| `frequency` | int | yes | 100 | `>= 0` | yes | Runtime regrid cadence ([#74](https://github.com/szaghi/adam/issues/74)). `0`: no runtime regridding, the AMR of the initial condition only. `n > 0`: after every `n`-th step the markers are combined and the grid regridded (up to `iters` sweeps); fatal without markers and on a multi-realm forest. The FNL backend regrids by a host round trip (each regrid logs its wall time). Negative: fatal. |
| `iters` | int | yes | 5 | any | yes | Maximum marker sweeps per AMR update, at initialisation and at each runtime regrid (it stops early when the grid is stable); a sweep refines or coarsens one level. |
| `regrid_prolongation` | string | no | `conservative` | `linear`, `conservative` | yes | How a refined block's children are filled from the parent ([#74](https://github.com/szaghi/adam/issues/74)). `conservative` (the FLUME default): limited linear slopes whose children average to the parent, so a regrid keeps the conserved integrals and a positive variable positive. `linear`: the library's tensor linear interpolation (the library default, kept by the other apps), not conservative. Any other value is fatal. Derefinement is always the mean of the children. Used by the runtime regrids; the initial AMR is unaffected (the initial condition is set again on the final grid). |
| `markers_number` | int | yes | 0 | `>= 0` | yes | Number of `[amr_marker_N]` sections. With markers and `max_level > iu_ref_levels` (2:1 refinement possible), every refined non-null axis (x, y; z too on an octree) needs an even block cell count `ni`/`nj`/`nk`: an odd one is fatal at initialisation (issue #39; it used to give NaN silently). The same axes need at least `2 ngc` cells (6 with WENO5): a thinner block made the coarse ghosts read stale fine ghosts, silently; fatal at initialisation (issue #66). |
| `seam_ghost_fill` | string | no | `tricubic` | `injection`, `restriction-compatible`, `tricubic`, `quintic`, `weno6`, `quintic-mp`. Any other value is fatal. The 6-point fills (`quintic`, `weno6`, `quintic-mp`) need at least 6 cells per refined axis in every block (fatal otherwise); on a `coupling = refined` forest seam both realms must declare the same value. CPU and FNL ([#68](https://github.com/szaghi/adam/issues/68) D1, D2). | yes | Coarse-to-fine ghost fill at 2:1 AMR seams. `tricubic`: the tensor 4-point Lagrange interpolation at the fine ghost's centre, exact for per-direction degree 3. `quintic`: the tensor 6-point one, exact for degree 5, its footprint shifted inside the donor block at its edge; with `seam_restriction = lagrange6` and `reflux = .false.` the composite reaches the scheme's order (AO `euler-q6`: orders 6.1, the error below the uniform coarse run's). `weno6`: per direction a WENO-Z combination of the cubics of the 4-point sub-stencils containing the ghost's interval, dropping those that hold a jump, the clamped linear interpolation when none is left (order 6 centred, 5 and 4 towards the block edge). `quintic-mp`: the quintic clamped to a monotonicity-preserving bound (order 6, never outside the interval's values across a jump). Neither limited fill improves a shock crossing over `tricubic` (D2 in the verification); they are options for smooth flows that need a fill which cannot overshoot. |
| `seam_restriction` | string | no | `mean` | `mean`, `lagrange6`, `weno6`. Any other value is fatal. `lagrange6` and `weno6` need at least 6 cells per refined axis in every block, the null axes an octree refines included (fatal otherwise); on a `coupling = refined` forest seam both realms must declare the same value. CPU and FNL ([#68](https://github.com/szaghi/adam/issues/68) C1, C2, C3). | yes | Fine-to-coarse ghost fill at 2:1 AMR seams. `mean`: the mean of the 2 (4, 8) fine cells, second order. `lagrange6`: the tensor 6-point Lagrange interpolation to the coarse centre, its footprint shifted inside the donor block at its edge, exact for per-direction degree 5 ([#68](https://github.com/szaghi/adam/issues/68) C1). `weno6`: per direction a WENO-Z combination of the cubic interpolants of the 4-point sub-stencils containing the coarse centre, dropping those that hold a jump and falling back to the mean of the fine pair when none is left (a jump inside the pair); exact for cubics, design order 6 centred, 5 one cell from the donor block's edge and 4 at it (no sub-stencil extrapolates) ([#68](https://github.com/szaghi/adam/issues/68) C3). With the reflux on the composite order stays 2 whatever the restriction; with `reflux = .false.` it rises to the prolongation's (4 with `tricubic`, AO `euler-l6` in the verification). |

Notes on the four `*_prune` keys:

- If all four are `>= 0`, the tree uses them to size the neighbour/BC search box (`2**(l-l_prune)*(ijk_prune+1)`).
- The verification inputs use `0, 0, 0, -1`, which is equivalent to "off".
- **Keep at least one of them negative** (the safest choice is all four `-1`).

**Trap:** the tree reader passes the same uninitialised buffer to every key and **ignores the
error**. A missing key therefore takes the value of the key read just before it. For example, a missing `max_level`
becomes equal to `ratio`, and a missing `ratio` is garbage. The declared defaults never take effect. **Always write
all seven tree keys.**

**Two conditions for the markers to take effect:**

- The markers act only while the initial conditions are imposed, and only if `[initial_conditions].amr_iterations > 0`.
  Each of those passes runs one AMR update of up to `iters` sweeps.
- Runtime AMR does not exist in FLUME.

### [amr_marker_N] (N = 1..markers_number)

Keys that are read but unused by FLUME's handling of a mode must still be present when marked `yes`.

**Runtime regrids** (`[amr] frequency > 0`, [#74](https://github.com/szaghi/adam/issues/74)): the markers are combined into one set of flags. Each marker votes per block, refine, keep or coarsen; a block is refined when any marker asks for it, coarsened only when every marker agrees, and kept otherwise. No block is coarsened below the base level `iu_ref_levels`. As a vote, the box marker asks to refine inside the box below `target_level`, to keep inside at it, and to coarsen elsewhere. The Löhner marker (mode 4) reads ghost cells, which the initial AMR does not fill: it acts at the runtime regrids only, and the initial AMR, with the markers applied one after the other as before, skips it.

| Key | Type | Req. | Accepted / invalid | Used | Meaning |
|-----|------|------|--------------------|------|---------|
| `mode` | int | yes | `1` = geometric (`AMR_GEO`), `2` = gradient (`AMR_GRAD`), `3` = total variation (`AMR_TV`), `4` = Löhner (`AMR_LOHNER`). FLUME supports 1, 2 and 4 (4 at runtime regrids only, see below). `3` or any other value → fatal "mode ... is not supported by FLUME" at the first AMR pass. An unknown mode reads no further keys. | yes | Marker kind. |
| `delta_type` | string | yes (all modes) | `x`, `y`, `z`, `max` (max over active directions). Other → fatal at marking. | modes 1 (solid) and 2 | Which block spacing is compared with `delta_fine`/`delta_coarse`. |
| `delta_fine` | real | yes (all modes) | no check | modes 1 (solid) and 2 | Admissible spacing where the criterion fires. |
| `delta_coarse` | real | yes (all modes) | no check | modes 1 (solid) and 2 | Admissible spacing elsewhere. A block is derefined when twice its spacing is still `<=` this. |
| `geo_type` | string | mode 1 | `solid`, `primitive-box`, `stl`. Other → fatal. `stl` → fatal in FLUME ("not supported"). | yes | Geometric marker kind. |
| `solid` | int | mode 1 + `solid` | `1..[solids].number`, otherwise fatal at marking | yes | Refine blocks crossed by that solid's surface. |
| `box_xmin`, `box_ymin`, `box_zmin`, `box_xmax`, `box_ymax`, `box_zmax` | real | mode 1 + `primitive-box` | no check | yes | Box: blocks whose centroid lies inside it are refined. |
| `target_level` | int | mode 1 + `primitive-box` | no check (capped by `max_level`) | yes | Refine the in-box blocks until they reach this level. The marker is additive and never derefines. |
| `stl_filename` | string | mode 1 + `stl` | — | no (`stl` is fatal in FLUME) | STL surface. |
| `field` | int | modes 2, 3, 4 | `1` = conservative `q`, `2` = auxiliary `q_aux`. Other → fatal at marking. | modes 2, 4 | Array holding the marker variable. |
| `var` | int | modes 2, 3, 4 | 1..nv (field 1) or 1..nv_aux (field 2), otherwise fatal | modes 2, 4 | Variable index (see the table below). |
| `tol` | real | modes 2, 3 | no check | mode 2 | Gradient-magnitude threshold: `max\|grad var\| > tol` gives `delta_fine`. The gradient is centred, across the faces between blocks too (a jump lying on a block face is seen, [#76](https://github.com/szaghi/adam/issues/76)), one-sided at the domain boundary. |
| `refine_tol`, `derefine_tol` | real | mode 4 | `0 <= derefine_tol < refine_tol <= 1`, otherwise fatal | yes | Löhner estimator thresholds: a block whose largest estimator exceeds `refine_tol` is refined, one below `derefine_tol` coarsened, one between kept (hysteresis). |
| `epsilon` | real | no (mode 4) | `>= 0`, default 0.01 | yes | Noise filter of the estimator: ripples smaller than about `epsilon` times the variable do not mark. |
| `floor` | real | no (mode 4) | `>= 0`, default 0 | yes | Absolute noise filter: variations of the variable well below `floor` do not mark. Needed on a variable that vanishes in quiet regions (the field outside a magnetic loop), where `epsilon`, relative to the variable, filters nothing and round-off noise reads as a jump; a fraction of the feature's jump is a fair value. |
| `buffer` | int | no (mode 4) | `0..ngc-1`, default 1 | yes | Ghost layers the estimator reads, so a feature at a neighbour's edge marks the block too. |

Variable indices for the gradient marker:

| Index | Conservative (`field = 1`) | Auxiliary (`field = 2`) |
|-------|----------------------------|-------------------------|
| 1 | `r` | `rho` |
| 2 | `ru` | `u` |
| 3 | `rv` | `v` |
| 4 | `rw` | `w` |
| 5 | `rE` | `p` |
| 6 | `bx` (MHD) | `T` |
| 7 | `by` (MHD) | `H` |
| 8 | `bz` (MHD) | `a` |
| 9 | `psi` (MHD + GLM or EGLM) | `Bx` (MHD) |
| 10 | — | `By` (MHD) |
| 11 | — | `Bz` (MHD) |

### Setting up runtime AMR

Runtime AMR ([#74](https://github.com/szaghi/adam/issues/74)) regrids during the run, on both backends. A minimal setup
tracks a feature from a base level to a finest level (the spacings below are for a unit box with 16 cells per block):

```ini
[amr]
max_level      = 3          ; finest level
iu_ref_levels  = 1          ; base level: no block is ever coarsened below it
frequency      = 5          ; regrid after every 5th step
iters          = 2          ; up to 2 sweeps per regrid (a sweep moves a block by one level)
markers_number = 2
; regrid_prolongation = conservative   (the FLUME default; the refined children average to the parent)

[amr_marker_1]               ; refines the initial grid too (Loehner does not act at initialisation)
mode         = 2             ; gradient
delta_type   = max
delta_fine   = 0.0079        ; ~1 % above the finest spacing: the blocks where the gradient fires reach it
delta_coarse = 0.0316        ; ~1 % above the base spacing
field        = 1
var          = 1             ; density
tol          = 0.05

[amr_marker_2]
mode         = 4             ; Loehner: scale-free, ~0 where smooth, ~1 at a jump
delta_type   = max           ; read but unused by mode 4
delta_fine   = 0.0
delta_coarse = 0.0
field        = 1
var          = 1
refine_tol   = 0.5           ; refine above
derefine_tol = 0.2           ; coarsen below (keep between: the hysteresis band)
buffer       = 2             ; ghost layers read: a feature at a neighbour's edge marks the block
; floor      = 1e-4          ; needed when the variable vanishes in quiet regions (e.g. B outside a loop)
```

How the pieces combine, and what to expect:

- **Markers vote.** A block is refined when any marker asks for it and coarsened only when all agree; the box marker
  votes to coarsen outside its box. Use a gradient (or solid, or box) marker to shape the initial grid and Loehner to
  follow the solution: without an initial marker the first `frequency` steps run on the base grid, and their error
  persists.
- **Thresholds.** Loehner's `refine_tol` around 0.5 follows shocks and contacts; smooth structures (a vortex, a
  rarefaction) score low and may stay coarse, so mark them with a gradient marker. On a variable that is zero in quiet
  regions set `floor` to a fraction of the feature's jump, or round-off noise marks the whole grid.
- **Cost.** The time step is global (no subcycling), so a run regridding to level L advances at the step of level L
  everywhere; the savings come from the cells, and on the verification cases the wall time dropped by 3–19 %
  ([verification](./verification#av-accuracy-of-runtime-amr)). On the FNL backend each regrid is a host round trip
  whose copies scale with the blocks in use, about 0.02 s per regrid on a 50-block grid on the development box
  ([#75](https://github.com/szaghi/adam/issues/75)).
- **Limits.** Single-realm runs only (a multi-realm forest with `frequency > 0` is refused); with EGLM the 2:1 seams
  raise the divergence error B_z well above its uniform-grid level ([#78](https://github.com/szaghi/adam/issues/78)).
- **Restart.** A run restarted from any step regrids on the same steps and continues bitwise (the regrid happens before
  the step's output and restart files).
- **Log.** Each regrid prints `flume: regrid at step N: blocks a -> b, R refined, C coarsened, S sweeps` (plus its wall
  time on FNL); a state with non-positive density or pressure after a regrid stops the run with its location.

---

## [field]

| Key | Type | Req. | Used | Meaning |
|-----|------|------|------|---------|
| `nv` | int | **not read** | **no** | FLUME passes `nv` from `[physics]` (5, 8 or 9), so `field_object` never reads this key. It can be omitted, and a wrong value is harmless. |

---

## [runge_kutta]

| Key | Type | Req. | Accepted / invalid | Used | Meaning |
|-----|------|------|--------------------|------|---------|
| `scheme` | string | yes | See the table below. Unknown → fatal. `runge-kutta-yoshida` is accepted by the library but fatal in FLUME ("no CPU/FNL time integrator"). | yes | Time integrator. It also sizes the per-block memory budget (stage fields). |

Scheme names accepted by the library:

| Name | Kind | Stored stage fields | FLUME |
|------|------|---------------------|-------|
| `runge-kutta-1` | low-storage, 1 stage (forward Euler) | 1 | accepted |
| `runge-kutta-2` | low-storage TVD, 2nd order | 1 | accepted |
| `runge-kutta-3` | low-storage TVD, 3rd order | 1 | accepted |
| `runge-kutta-ssp-11` | SSP, 1 stage (forward Euler) | 1 | accepted |
| `runge-kutta-ssp-22` | SSP, 2nd order | 2 | accepted |
| `runge-kutta-ssp-33` | SSP, 3rd order | 3 | accepted |
| `runge-kutta-ssp-54` | SSP, 5 stages, 4th order | 5 | accepted |
| `runge-kutta-yoshida` | symplectic, 4th order | 0 | **refused** (fatal) |

In a forest manifest, seams with `coupling_cadence = stage_coincident` require both realms to use the same `scheme`
(fatal otherwise;).

---

## [weno]



| Key | Type | Req. | Default | Accepted / invalid | Used | Meaning |
|-----|------|------|---------|--------------------|------|---------|
| `scheme` | string | yes | — | See the table below. Unknown → fatal. `weno-c-*` is accepted by the library but **fatal in FLUME**. | yes | WENO order (stencil half-width `S`). |
| `weights` | string | no | `js` | `js` (Jiang–Shu, absolute $\varepsilon$), `si` (scale-invariant). Unknown → fatal. | yes | Nonlinear weights ([numerics](./numerics#weno-reconstruction)). `si` makes the scheme independent of the units of the data; the default `js` keeps every existing result bitwise. `si` with `weno-riemann` and `reconstruction_variables = primitive` is fatal. |
| `ror_number` | int | yes | 0 | `>= 0` | **no** | Number of order-reduction (ROR) stages. FLUME does not implement ROR, and the scheme overwrites `S`. |
| `ror_scheme_1` … `ror_scheme_<ror_number>` | int | cond. (`ror_number > 0`) | — | — | **no** | ROR stencil sizes. |
| `ror_threshold` | real | yes | 0.9 | — | **no** | ROR trigger. |
| `ror_vars_number` | int | yes | 2 | `>= 0` | **no** | Number of ROR check variables. |
| `ror_ivar_1` … `ror_ivar_<ror_vars_number>` | int | cond. (`ror_vars_number > 0`) | — | — | **no** | ROR check variables. |
| `enable_ror_stats` | logical | yes | `.false.` | logical | **no** (but `.true.` allocates an `nb x cells x 3` integer array) | ROR statistics. |
| `ib_reduction_extent` | int | yes | 0 | — | **no** | Order reduction near immersed solids (FLUME does not read `cell_scheme`). |
| `ib_reduced_order` | int | yes | 1 | — | **no** | Reduced `S` near solids. |

Scheme names:

| Name | S | Order | Kind | FLUME |
|------|---|-------|------|-------|
| `weno-u-1` | 1 | 1 | upwind | accepted (`weno-riemann` still needs `ngc >= 2`) |
| `weno-u-3` | 2 | 3 | upwind | accepted |
| `weno-u-5` | 3 | 5 | upwind | accepted |
| `weno-u-7` | 4 | 7 | upwind | accepted |
| `weno-u-9` | 5 | 9 | upwind | accepted |
| `weno-c-2` | 1 | 2 | centred | **refused** (fatal) |
| `weno-c-4` | 2 | 4 | centred | **refused** (fatal) |
| `weno-c-6` | 3 | 6 | centred | **refused** (fatal) |
| `weno-c-8` | 4 | 8 | centred | **refused** (fatal) |

Two further constraints:

- The ghost-cell requirement is `ngc >= S`.
- With `scheme_space = weno-riemann`, the same `scheme` sets the order of the WENO *interpolation* tables
 .

---

## [linear-algebra]

FLAIL elliptic-solver settings. **FLUME never runs an elliptic solve, so every key here is unused.**
The keys are still required: every key except `smoothing` has an `error_stop` on absence. `smoothing` also errors when
missing.

| Key | Type | Req. | Default | Accepted / invalid | Used |
|-----|------|------|---------|--------------------|------|
| `smoothing` | string | yes | (unallocated) | See the notes below the table. | no |
| `iterations` | int | yes | 1 | — | no |
| `iterations_init` | int | yes | 1 | — | no |
| `iterations_fine` | int | yes | 1 | — | no |
| `iterations_coarse` | int | yes | 1 | — | no |
| `tolerance` | real | yes | 1e-6 | — | no |

Accepted spellings for `smoothing`:

| Method | Accepted spellings |
|--------|--------------------|
| Multigrid | `MULTIGRID`, `multigrid`, `Multigrid` |
| Gauss-Seidel | `GAUSS-SEIDEL`, `gauss-seidel`, `Gauss-Seidel` |
| SOR | `SOR`, `sor`, `Sor` |
| SOR with OpenMP | `SOR-OMP`, `sor-omp`, `Sor-omp` |

**Trap:** an unrecognised `smoothing` value is **not** rejected. The select has no `case default`, so the string stays
unallocated and is then printed by `description()`. That is undefined behaviour and in practice
usually crashes at initialisation.

---

## [fdv]

Read by the realm.

| Key | Type | Req. | Accepted / invalid | Used | Meaning |
|-----|------|------|--------------------|------|---------|
| `fdv_scheme` | string | yes | `FD`/`fd`/`Fd`/`fD` or `FV`/`fv`/`Fv`/`fV`. Other → fatal. | **no** | FLUME always uses centred FD for its div(B) diagnostic, whatever this value. |
| `fdv_order` | int | yes | Even values: 2, 4, 6. FNL: `fdv_order/2 <= min(ngc, 3)`. CPU: `fdv_order/2 <= ngc`. Violations are fatal at the first div(B) diagnostic (MHD only). An odd value is rounded down by the `/2`. | MHD only | Order of the centred div(B) stencil used for `divb_history.dat`, the `divb` output field and the `[mhd].divb_tol` monitor. Unused for Euler. |

---

## [solids] and [solid_N]

Library immersed boundary. Solids are **Euler only**: `number > 0` with `mhd-ideal` is fatal.

### [solids]

| Key | Type | Req. | Default | Used | Meaning |
|-----|------|------|---------|------|---------|
| `number` | int | yes | 0 | yes | Number of `[solid_N]` sections. `0` means no immersed boundary. |
| `n_eikonal` | int | yes | 2 | yes (with solids) | Eikonal extrapolation sweeps per residual evaluation, filling the solid cells. |

### [solid_N] (N = 1..number)

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `name` | string | yes | any | Label (printed only). |
| `bc_type` | int | yes | See the notes below the table. | Wall treatment of the solid-cell state. |
| `definition` | string | yes | `analytical_sphere`, `analytical_circle`, `analytical_rectangle`. See the notes below the table. | Geometry kind. |
| `sphere_center_x`, `sphere_center_y`, `sphere_center_z`, `sphere_radius` | real | `definition = analytical_sphere` | no check | Sphere. |
| `circle_center_x`, `circle_center_y`, `circle_center_z`, `circle_radius` | real | `analytical_circle` | no check | Infinite cylinder. |
| `circle_axis` | char | `analytical_circle` | `x`, `y`, `z`. Any other value computes no distance function, silently. | Cylinder axis. |
| `rectangle_center_x`, `rectangle_center_y`, `rectangle_center_z`, `rectangle_edge_1`, `rectangle_edge_2` | real | `analytical_rectangle` | no check | Infinite rectangular prism. |
| `rectangle_axis` | char | `analytical_rectangle` | `x`, `y`, `z`. Any other value is silently ignored. | Prism axis. |

Notes on `bc_type`:

- `1` = viscous (all momentum negated inside the solid).
- `2` = inviscid (normal momentum reflected).
- Any other value applies no inversion, silently.

Notes on `definition`:

- `file.off` appears in the definitions list but is not implemented.
- An unknown value reads no geometry keys and **silently produces no solid** (the distance function stays -1).

---

## [slices] and [slice_N]

### [slices]

| Key | Type | Req. | Used | Meaning |
|-----|------|------|------|---------|
| `slices_number` | int | yes | yes | Number of `[slice_N]` sections. `0` means no slices. |

### [slice_N] (N = 1..slices_number)

**None of these keys is error-checked, and none has a default.** A missing key leaves the value undefined.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `itype` | string | **no check** | `trilinear`, `inverse_distance`. Other → fatal (FLUME check). | Interpolation from the cells to the slice points. |
| `n_save` | int | **no check** | Must be `> 0`. `0` makes `mod(it, 0)` crash at runtime. | Save cadence in steps. The last step is also saved. |
| `ni`, `nj`, `nk` | int | **no check** | `>= 1` | Points of the slice lattice (cell-centred on the box). |
| `emin_x`, `emin_y`, `emin_z`, `emax_x`, `emax_y`, `emax_z` | real | **no check** | — | Box sampled by the lattice. |

Only the conservative variables are sliced.

---

## [reference] (optional: dimensional input)

FLUME's equations carry no dimensionless number (`B` is `B_SI/sqrt(mu0)`), so a dimensional input runs once every
value is divided by the reference of its dimension (issue #49). With this section present, FLUME does that once, on
the loaded file, before any other section is parsed; without it nothing changes. References: density `rho0`, length
`L0`, velocity `u0`; derived: time `L0/u0`, pressure and energy density `rho0 u0^2`, field `u0 sqrt(rho0)`, the GLM
`psi` `u0^2 sqrt(rho0)` (the EGLM one as the field). `cp`, `cv` are replaced by `gamma = cp/cv` (the gas constant
becomes 1, the temperature unit `u0^2/R`). Restart files and logs stay in code units; the other outputs follow
`output_units`. Every restart save also writes `<restart_basename>.reference` (the three references, 1 without this section); a restart under different
references is refused, and restart files without that record (written before it existed, so in code units) restart only
when the references are all 1.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `density` | real | no (1) | `> 0`, otherwise fatal | `rho0`. |
| `length` | real | no (1) | `> 0`, otherwise fatal | `L0`. |
| `velocity` | real or string | no (1) | `> 0`; `acoustic` (`sqrt(gamma pressure / density)`, needs `pressure`); `alfvenic` (`field / sqrt(density)`, needs `field` and `mhd-ideal`). Other → fatal. | `u0`. |
| `pressure` | real | cond. | `> 0` | Reference pressure of the `acoustic` preset. |
| `field` | real | cond. | `> 0` | Reference field of the `alfvenic` preset. |
| `output_units` | string | no (`code`) | `code`, `dimensional`; other → fatal | Units of the fields, grid, time, slices and histories. |

With `output_units = dimensional` every written value is multiplied, at write time only, by the reference of its
dimension: the fields (conservative, residuals, auxiliary, MHD derived: `beta` is dimensionless, `divb` a field per
length), the block origins and spacings, the time, the slice points and values, the conservation integrals (times
`L0^3`, the cell volume includes the null directions), the div(B) and residual histories. The temperature is written
as `p / (rho R)` with the gas constant of the input (`cp - cv`, or 1 with `gamma`). `<output_basename>.units` (written
whenever the section is present) records the output units, the references, the derived ones and the factor of every
written variable. The state, the restart files and the solver are unaffected.

Any other key is fatal. With the section active **every** section and key of the file must be known to the layer:
an unknown one is fatal, so no dimensional value can pass unconverted. Keys whose dimension depends on the context are
resolved from it: `[amr_marker_N].tol` of a gradient marker takes the dimension of the marked variable per length,
`floor` of a Löhner marker that of the marked variable (a marker on the temperature is refused; the Löhner
`refine_tol`, `derefine_tol` and `epsilon` are dimensionless), `[initial_conditions].wave_amplitude` that of the
eigenvector it multiplies.
`type = orszag-tang` (a state hard-coded in code units) and multi-realm runs are refused. The log lists every
converted key with its old and new value.

The dissipative coefficients and the wall keys (issue #49 N4, issue #65) convert as follows; the numbers are
dimensionless and pass unchanged. $R$ is the gas constant of the input ($c_p - c_v$, or 1 with `gamma`), since the layer
replaces `cp`, `cv` by `gamma` and the code temperature becomes $p/\rho$:

| Key | Divided by |
|---|---|
| `viscosity` | $\rho_0 u_0 L_0$ |
| `conductivity` | $\rho_0 u_0 L_0 R$ |
| `resistivity` | $u_0 L_0$ |
| `wall_u`, `wall_v`, `wall_w` | $u_0$ |
| `wall_temperature`, `reference_temperature` | $u_0^2 / R$ |
| `reynolds`, `prandtl`, `magnetic_reynolds`, `viscosity_exponent` | — (dimensionless) |
| `lundquist` | — (refused unless `velocity = alfvenic`) |
| `[gravity]` `g_x`, `g_y`, `g_z` | $u_0^2 / L_0$ |
| `[gravity]` `G` | $u_0^2 / (\rho_0 L_0^2)$ |
| `[gravity]` `froude`, `G_star`; `[initial_conditions]` `polytropic_gamma` | — (dimensionless) |

The conversion divides by the power-of-two reference first and applies $R$ last, a single rounding, so the
verification DC can check every converted value exactly.

---

## [physics]

`physical_model` is required, with either `gamma` or both `cp` and `cv`.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `physical_model` | string | yes | `euler` (nv = 5) or `mhd-ideal` (nv = 8, or 9 with GLM or EGLM). Other → fatal. | Equation set. `mhd-ideal` also loads `[mhd]`. |
| `gamma` | real | cond. | `> 1`; with `cp` or `cv` → fatal | Specific heats ratio; the gas constant is 1 (code units). |
| `cp` | real | cond. | `cp > cv > 0`, otherwise fatal; required without `gamma` | Specific heat at constant pressure. |
| `cv` | real | cond. | as above | Specific heat at constant volume. `gamma = cp/cv`, `R = cp - cv`. |

### Dissipative terms (issue #65, M4)

Each term is given either as its coefficient or as the dimensionless number it stands for. In code units the
references are 1, so a number is the reciprocal coefficient. Giving two keys of a term is fatal; giving none leaves the
term off, so an input without these keys stays ideal and bitwise unchanged. The viscosity and the conductivity drive
the Navier–Stokes fluxes of every model, the resistivity the Ohmic ones of `mhd-ideal`
([numerics](./numerics#dissipative-fluxes-navier-stokes)).

| Key | Type | Accepted / invalid | Meaning |
|-----|------|--------------------|---------|
| `viscosity` | real | `>= 0`; with `reynolds` → fatal | Dynamic viscosity $\mu$ (code units, or dimensional with [`[reference]`](#reference-optional-dimensional-input)). |
| `reynolds` | real | `> 0` | $Re$: $\mu = 1/Re$. |
| `conductivity` | real | `>= 0`; with `prandtl` → fatal | Thermal conductivity $k$; the heat flux is $-k\nabla T$ with $T = p/(\rho R)$. |
| `prandtl` | real | `> 0`; needs `viscosity` or `reynolds` | $Pr$: $k = \mu c_p / Pr$ (follows $\mu(T)$ under the power law). |
| `resistivity` | real | `>= 0`; MHD only; with a number → fatal | Magnetic diffusivity $\eta$ (Ohmic resistivity). |
| `magnetic_reynolds` | real | `> 0`; MHD only | $Rm$: $\eta = 1/Rm$. |
| `lundquist` | real | `> 0`; MHD only; with `[reference]` it needs `velocity = alfvenic` | $S$: $\eta = 1/S$, the magnetic Reynolds number at the Alfvén speed. |
| `viscosity_law` | string | `constant` (default), `power-law`; other → fatal; `power-law` needs a viscosity | Temperature law of the viscosity. |
| `viscosity_exponent` | real | `>= 0`; required by `power-law`, fatal without it | $\omega$ in $\mu(T) = \mu\,(T/T_\mathrm{ref})^\omega$ (non-negative: the diffusive time step bound evaluates the laws at the largest temperature a face can see). |
| `reference_temperature` | real | `> 0`; required by `power-law`, fatal without it | $T_\mathrm{ref}$ (code temperature $p/(\rho R)$). |

`[numerics] positivity_limiter = cell` is refused with any non-zero coefficient: the limiter's first-order backbone
with a central dissipative flux is not admissible without the extension of Zhang (2017, *J. Comput. Phys.* 328), which
is not implemented (D-M4-5). A coefficient is also fatal with immersed solids (their walls are inviscid), and
`dissipative_order = 4` with `[grid] ngc < 3`.

---

## [mhd] (only with `physical_model = mhd-ideal`)

The section is not read for `euler`.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `divergence_control` | string | yes | `glm` (nv = 9, psi added), `eglm` (nv = 9, psi also in the energy; the `glm_*` keys apply) or `none` (nv = 8). Other → fatal. | div(B) control. Under `eglm` the AMR reflux leaves B alone (ρ, ρu, E and psi are refluxed): with the Derigs sources B is not a conserved integral, and refluxing it broke the sources' cancellation of the flux-generated $B_z$ at 2:1 seams ([#68](https://github.com/szaghi/adam/issues/68) G3b, #78). |
| `glm_ch` | real | `glm`, `eglm` | `> 0`, otherwise fatal | Constant cleaning speed `c_h`. It also bounds the time step (`c_h * sum 1/dx`). Keep it just above the fastest signal speed: a larger one speeds up the growth of small asymmetries in sheared, magnetised flows ([#83](https://github.com/szaghi/adam/issues/83)). |
| `glm_alpha` | real | `glm`, `eglm` | `>= 0`, otherwise fatal | Damping parameter: `c_h^2/c_p^2 = glm_alpha * c_h / L`. |
| `glm_damping_length` | real or `min-cell` | `glm`, `eglm` | A positive real, or the literal `min-cell` (the minimum cell spacing over the active directions, MPI-reduced; with runtime regridding, `[amr] frequency > 0`, the spacing of the finest level `max_level` allows, fixed for the run). A non-number or a value `<= 0` is fatal. | Damping length `L`. |
| `glm_ch_check` | string | `glm`, `eglm` | `warning` or `error`. Other → fatal. Without GLM it is forced to `warning` and not read. | What happens when `max(\|u\|+c_f) > glm_ch`: warn (logged each time a new maximum appears) or stop. |
| `divb_tol` | real | yes | `>= 0`, otherwise fatal. `0` disables the monitor; `> 0` needs `[diagnostics] divb_flux = .true.` (hence GLM or EGLM), fatal otherwise. | Monitor on the normalised D_flux ([#68](https://github.com/szaghi/adam/issues/68) G1): `max h\|D_flux\| / max\|B\| > divb_tol` warns or stops. Dimensionless (no `[reference]` conversion): the divergence of the interface field the cleaning acts on, scaled by the spacing and the field, bounded on a resolved run with an O(H) seam floor (A2). |
| `divb_error` | logical | yes | logical | `.true.`: exceeding `divb_tol` is fatal. |
| `rho_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Density positivity floor. |
| `p_floor` | real | yes | `>= 0`, otherwise fatal. `0` disables it. | Pressure positivity floor. With **both** floors at 0, any non-positive density or pressure stops the run. |

---

## [gravity] (optional: external gravity and self-gravity)

A constant acceleration $\mathbf g$ ([#73](https://github.com/szaghi/adam/issues/73), milestone M7) and the
self-gravity of the gas ([#81](https://github.com/szaghi/adam/issues/81), milestone M8). Without the section there is
no gravity and nothing changes.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `g_x`, `g_y`, `g_z` | real | one of the two forms | any; each optional (default 0) | The acceleration. Dimensional (length / time²) with [`[reference]`](#reference-optional-dimensional-input). |
| `froude` | real | the other form, with `axis` | `> 0`, otherwise fatal | The number form: $g = 1/\mathrm{Fr}^2$ along `axis` ($g L_0/u_0^2$ in the reference units; not converted by `[reference]`). |
| `axis` | string | with `froude` | `+x`, `-x`, `+y`, `-y`, `+z`, `-z`; other → fatal | The signed direction of $\mathbf g$. |
| `equilibrium` | string | no (default `none`) | `none` or `atmosphere`; other → fatal. `atmosphere` needs an `isothermal-atmosphere`, `polytropic-atmosphere` or `rayleigh-taylor` initial condition, and is fatal with `[numerics] reflux_policy = troubled` and on multi-realm runs | The well-balanced mode (M7-G2b): the initial atmosphere at rest (zero velocity, region 1's field, no pulse or seed) is kept to round-off. See below. |
| `potential` | string | no (default `uniform`) | `uniform` or `linear`; other → fatal. `linear` is fatal with `equilibrium = atmosphere` and with `self_gravity` | How g reaches the sources. `uniform`: the constant acceleration in every kernel, with the fourth-order energy source. `linear` ([#81](https://github.com/szaghi/adam/issues/81), M8-P1): the same field written as the cell potential Φ = −g·x and applied through the potential path of self-gravity (face gravity −ΔΦ/Δx, second order); a test of that path ([SV-1](./verification#sv-1-p1-gravity-from-a-cell-potential)). |
| `self_gravity` | logical | no (default `.false.`) | with it: `equilibrium = atmosphere`, `potential = linear`, multi-realm runs, `[amr] frequency > 0` (runtime AMR: #81 P5c), immersed solids and a build without an FFT library (the `*-fft` modes) are fatal, and so are the boundaries `poisson_bc` cannot hold (below) | Self-gravity (M8-P3): $L_h\Phi = 4\pi G(\rho - \bar\rho)$ (periodic) or $L_h\Phi = 4\pi G\rho$, $\Phi \to 0$ at infinity (isolated, M8-P4), solved every Runge–Kutta stage, applied through the potential path. See below. |
| `G` | real | one of `G`, `G_star` with `self_gravity` | `> 0`, otherwise fatal; fatal without `self_gravity`, and together with `G_star` | The gravitational constant. Dimensional with [`[reference]`](#reference-optional-dimensional-input). |
| `G_star` | real | the other form | as `G` | The number $G\rho_0L_0^2/u_0^2$; not converted by `[reference]`. |
| `poisson_solver` | string | no (default `auto`) | `auto`, `fft` or `multigrid`; other → fatal. `fft` is fatal on an AMR grid | `auto`: the FFT solver on a uniform grid, the multigrid solver on an AMR grid (leaves at more than one level, M8-P5b). `fft`: the library FFT solver only (FFTW on the CPU, cuFFT on FNL). `multigrid`: the multigrid solver on any grid (periodic only), a check against `fft` on a uniform one. |
| `poisson_bc` | string | no (default `periodic`) | `periodic` or `isolated`; other → fatal. `periodic` is fatal with a non-periodic active direction; `isolated` is fatal in 2-D, with a periodic direction, with non-cubic cells and with the multigrid solver (on AMR grids: [#87](https://github.com/szaghi/adam/issues/87)) | `periodic`: every active direction periodic, the mean density removed. `isolated` (M8-P4): the free-space potential of the mass on the grid (Hockney–Eastwood zero padding with the lattice Green's function of $L_h$); the fluid boundaries are the `[bc_*]` ones, e.g. `wall-inviscid` for a closed box. |

Giving `g_*` together with `froude` is fatal, as is a section with neither (unless `self_gravity`, where $\mathbf g$
defaults to 0).

$\mathbf g$ sets the `isothermal-atmosphere` and `polytropic-atmosphere` initial conditions and the `hydrostatic`
boundary condition (M7-G0), and drives the sources (M7-G1):

- momentum $\rho\mathbf g$; energy from the face mass fluxes the continuity equation uses,
  $\tfrac1{24}\sum_d g_d\,(-F^\rho_{d,i-3/2} + 13F^\rho_{d,i-1/2} + 13F^\rho_{d,i+1/2} - F^\rho_{d,i+3/2})$, so that with
  this static potential the total energy $E + \rho\Phi$ ($\Phi = -\mathbf g\cdot\mathbf x$) is conserved exactly
  (Mikami et al. 2008, Mullen et al. 2021: each face's weights sum to one over the cells) with periodic boundaries and
  across block interfaces; at a `wall-inviscid` wall along $\mathbf g$ the source reads the mirrored face beyond the
  wall, and the exactly conserved total is $\sum(E + \rho\Phi)V \pm \sum g_n\Delta x_n\rho V/12$ over the cells
  beside the wall (+ at the lower wall, − at the upper), which differs from the plain sum by $O(\Delta x^2)$. The fluxes of a
  finite-difference scheme are point values of the primitive function of $\rho\mathbf u$, not of $\rho\mathbf u$
  itself, so the two-face average $\tfrac12(F^\rho_{i-1/2} + F^\rho_{i+1/2})$ alone is a second-order cell value and
  takes the whole scheme to second order; the four-face form is fourth order (M7-G2a,
  [verification](./verification#gv-g2a-the-order-of-the-energy-source)). It reads the face beyond each block end, so
  **gravity needs `[grid] ngc >= S + 1`** ($S$ the WENO half-stencil, 3 for WENO-5; with `weno-riemann`
  $\max(S, 2) + 1$), fatal otherwise;
- with `positivity_limiter = cell`, the limiter's first-order backbone takes the time-centred energy
  $\rho\mathbf g\cdot(\mathbf u + \mathbf g\,\Delta t/2)$, which leaves its internal energy unchanged, and limits the rest
  of the face form with the cell factor $\Lambda$: gravity never makes the backbone inadmissible (near vacuum, or with
  $g\Delta t$ above the sound speed), and $E + \rho\Phi$ stays exact wherever the limiter does not act;
- with `equilibrium = atmosphere` (M7-G2b, the deviation method of Edelmann et al. 2021, §4.3), the residual of the
  atmosphere at rest is computed once per grid (at the first step, after a restart, after every regrid) with the
  run's own space operator and subtracted from every residual: the atmosphere is then an exact discrete steady state
  (round-off), on `weno` and `weno-riemann`, Euler and MHD, across AMR 2:1 seams (the equilibrium's reflux is folded
  into its residual), and a perturbation of it evolves without the stationary truncation flow of the unbalanced
  scheme ([verification](./verification#gv-g2b-the-well-balanced-mode)). It costs one more field (the stored
  residual) and one residual evaluation per grid. The subtracted residual is outside the positivity limiter's factor:
  it is the truncation error of the equilibrium, small against the state wherever the atmosphere is resolved;
- no time-step limit: the first step logs $|\mathbf g|\,\Delta t^2/(\mathrm{CFL}\,h_{\min})$, at least the velocity
  gravity adds in a step over the fastest wave speed (near 1 or above, gravity is stiff);
- refused (fatal) along a null direction and with immersed solids.

With `self_gravity` ([models](./models#self-gravity)):

- every Runge–Kutta stage solves $\Phi$ from that stage's density after the ghost update, with the library FFT solver:
  rank 0 gathers the grid and scatters $\Phi$ back with every ghost layer filled. Periodic: it divides each Fourier
  mode by the eigenvalue of the 7-point Laplacian $L_h$, so the discrete residual is round-off. Isolated: the density
  sits in a zero-padded array of $2(N + n_{gc})$ cells per direction and is convolved with the lattice Green's function
  of $L_h$, whose residual is about $10^{-13}$ per offset (`scripts/lattice_green.py`); the residual of $\Phi$ is that
  summed over the source, $2\cdot10^{-11}$ of $4\pi G\rho$ at $64^3$, and the ghosts beyond the domain hold the
  free-space potential. The solve is serial on rank 0 and, on FNL, the density and $\Phi$ cross the host every stage:
  correct, not scalable;
- on an AMR grid (M8-P5b) the multigrid solver (`adam_poisson_mg_object`) solves the composite Laplacian, whose 2:1
  faces carry the mimetic mean flux (symmetric and consistent; [SV-9](./verification#sv-9-p5a-the-poisson-solver-on-amr-grids)):
  flexible CG preconditioned by a V-cycle, its coarsest level solved by the FFT solver of the backend; on FNL it runs
  on the host. The face gravity at a 2:1 face is the mean-flux one on both sides (first order there, the potential
  second order), and the gravity work across each seam rides on the reflux (below), so $E + W$ stays exact in the
  semi-discrete form;
- across a 2:1 face the reflux, which corrects the coarse cells once at the end of the step, also hands the coarse cell
  the seam's gravity work, so that summed over the cells the energy sources equal $-\dot W$ exactly. Its stages,
  though, advance the coarse cells with their own seam fluxes: on AMR the drift of $E + W$ is first order in
  $\Delta t$, not of the Runge–Kutta order (SV-10). `reflux = .false.` or `reflux_policy = troubled` leave the work
  uncorrected where they leave the flux;
- the conservation history gains the columns `W` $= \tfrac12\sum\rho\Phi\,dV$ and `E_plus_W`, the energy the scheme
  conserves (each row solves $\Phi$ from the committed state);
- the first step logs $G$, the free-fall time $t_{ff} = \sqrt{3\pi/(32G\rho_{\max})}$ of the densest cell, $\Delta t/t_{ff}$
  (a warning above 0.1) and the wall time of one solve;
- a constant $\mathbf g$ may be given too; it adds through its own kernels;
- with `positivity_limiter = cell`, the limiter's gravity energy correction (above) acts on the potential path too:
  $E + W$ stays exact only where the limiter does not act. A star in a light atmosphere ($10^{-4}$ of its central
  density) needs the limiter at its surface and drifts by $2\cdot10^{-3}$ of $E + W$ over a free-fall time at $32^3$,
  against $10^{-9}$ without it (the SV-7 and SV-8 cases use a $10^{-2}$ atmosphere, which runs without it);
- **build:** `flume-cpu-gnu-fft` (FFTW 3, executable `adam_flume_cpu_fft`) or `flume-fnl-nvf-fft` (cuFFT,
  `adam_flume_fnl_fft`); the default modes refuse `self_gravity` at the first solve.

## [numerics]

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `scheme_space` | string | yes | `weno` (WENO flux splitting) or `weno-riemann` (WENO interpolation + Riemann flux + high-order correction). Other → fatal. | Spatial operator. MHD at very low β (≲ 1e-3, e.g. the Balsara–Spicer blast): `weno` needs the positivity limiter at every step and breaks point symmetry at the grid scale; use `weno-riemann` with `hlld` where symmetry or backend agreement matter ([#77](https://github.com/szaghi/adam/issues/77), [verification](./verification#the-mhd-blast-loses-its-point-symmetry-on-the-flux-splitting-path-77)). |
| `reconstruction_variables` | string | yes | With `weno`: `characteristic` or `conservative`. With `weno-riemann`: `characteristic` or `primitive`. Other → fatal. | Variables reconstructed (or interpolated) at the faces. |
| `riemann_solver` | string | `weno-riemann` | `llf`, `hll`, `hllc`, `hlld`. Other → fatal. By model: `euler` allows `llf`, `hll`, `hllc`; `mhd-ideal` allows `llf`, `hll`, `hlld` (fatal otherwise). Not read with `weno`. | Face Riemann solver. |
| `flux_correction` | string | `weno-riemann` | `6th`, `4th`, `none` (2nd order). Other → fatal. | Face-flux correction `F^ = c1 F + c2 (f_i+f_i+1) + c3 (f_i-1+f_i+2)`. |
| `flux_correction_sensor` | string | `weno-riemann` | `weno` or `none`. Other → fatal. | `weno`: the correction is switched off at faces where `min_k w_k/d_k < 0.2` (fixed threshold `FLUX_CORRECTION_SENSOR_TAU`). `none`: always on. |
| `reflux` | logical | yes | logical | Berger–Colella reflux at AMR coarse-fine faces (and inter-realm seams). `.false.` is a diagnostic only: the run is then not conservative across 2:1 faces. |
| `reflux_policy` | string | no | `always`, `troubled` (default `always`; any other value is fatal) | With `reflux = .true.`, which seam face cells the reflux corrects ([#68](https://github.com/szaghi/adam/issues/68) E1). `always`: every one, conservative to round-off. `troubled`: only those flagged in some stage of the step, by the positivity limiter (a blended seam face, issue #50) or by the WENO-5 sensor (the minimum weight ratio $\min_k w_k/d_k$ of density, normal velocity or pressure across the face below 0.2, the weno-riemann flux-correction sensor); elsewhere the coarse flux stays, so smooth regions keep the accuracy of the high-order seam operators at the cost of conservation there (Sod c2f: mass drift $1.5\cdot10^{-7}$, shock unchanged). The log reports the flagged face cells per step. CPU and FNL (E1b: the FNL backend gathers the seam cells' auxiliary variables to the host and flags them there). The default stays `always` (owner decision after E2: exact conservation); `troubled` is the switch that lets `seam_restriction = lagrange6` with `seam_ghost_fill = quintic` keep sixth order across 2:1 seams with the reflux on (AO `euler-tq`), at the price of a bounded conservation drift (about $10^{-5}$ per period on the smooth vortex at N = 64, the reflux-off rate). |
| `positivity_limiter` | string | no (default `none`) | `none` or `cell`. Other → fatal. `cell` is fatal with `[mhd] divergence_control = glm`, a non-SSP `[runge_kutta] scheme`, immersed solids, multi-realm runs and any dissipative coefficient (issue #65, D-M4-5). | `cell`: the cell-based positivity limiter ([numerics](./numerics#positivity-limiter)): every face flux blended with the first-order Lax–Friedrichs backbone so that each stage keeps the density and the pressure positive; the limited faces are logged per stage. |

| `dissipative_order` | int | no (default 4) | `2` or `4`. Other → fatal; `4` needs `[grid] ngc >= 3` (fatal otherwise). | Order of the dissipative face fluxes (issue #65, D-M4-2): `4` is the conservative 4th-order flux with the Shu–Osher correction, `2` the compact central one ([numerics](./numerics#dissipative-fluxes-navier-stokes)). An ideal run ignores it. |

---

## [bc_x_min], [bc_x_max], [bc_y_min], [bc_y_max], [bc_z_min], [bc_z_max]

All six sections are required, even along null directions.

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `type` | string | yes (every face) | `extrapolation`, `inflow`, `wall-inviscid`, `wall-noslip`, `wall-isothermal`, `periodic`, `hydrostatic`. Other → fatal. The two faces of an axis must be **both** periodic or **neither** (fatal). | See the list below. |
| `r`, `u`, `v`, `w`, `p` | real | `type = inflow` | no range check | Prescribed primitive inflow state. |
| `bx`, `by`, `bz` | real | `type = inflow` **and** `mhd-ideal` | no range check | Prescribed inflow magnetic field. `psi` is set to 0. |
| `wall_u`, `wall_v`, `wall_w` | real | no (default 0) with `wall-noslip`, `wall-isothermal` | the component normal to the face must be 0 (fatal) | Velocity of the wall (tangential): a moving lid, a Couette wall. |
| `wall_temperature` | real | `type = wall-isothermal` | `> 0`, otherwise fatal | Wall temperature, $T = p/(\rho R)$ in code units. |

What each `type` does:

- `extrapolation`: zeroth-order extrapolation.
- `inflow`: prescribed state.
- `wall-inviscid`: slip wall (mirror state with the normal momentum negated; for MHD the normal B is negated too).
- `wall-noslip`: adiabatic no-slip wall (the velocity reflected about the wall velocity, the temperature mirrored).
- `wall-isothermal`: isothermal no-slip wall (the velocity as above, the ghost temperature $T_w^2/T$, the geometric mirror). See
  [Boundary conditions](./boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls).
- `periodic`: library periodicity, so true periodic neighbours across blocks and ranks.
- `hydrostatic` ([#73](https://github.com/szaghi/adam/issues/73)): closes a stratified atmosphere under `[gravity]`.
  The ghost density and pressure continue the isothermal hydrostatic equilibrium from the first interior cell along
  the normal, $\rho_g = \rho_1 \exp(g_d s / T_1)$, $p_g = \rho_g T_1$ ($T_1 = p_1/\rho_1$, $s$ the signed distance
  from that cell to the ghost), so an isothermal atmosphere has exact ghosts. The velocity is mirrored with its normal
  component negated (no mass flux). The field and psi are mirrored **even**: field lines pass through the boundary, so
  a uniform field along or across gravity stays uniform (a reflecting conductor is `wall-inviscid`). Without gravity
  along the normal it is a slip wall that copies the thermodynamic state of the first interior cell.

In a forest, a face glued to another realm by the manifest still needs a `type` here. The seam overrides the ghost
fill.

---

## [initial_conditions]

`amr_iterations` and `type` are always required. The other keys depend on `type`, and **every key the
selected type uses is required** (fatal "failed to load" otherwise).

| Key | Type | Req. | Accepted / invalid | Meaning |
|-----|------|------|--------------------|---------|
| `amr_iterations` | int | yes | Negative values are clamped to 0. | Number of initial-condition + AMR passes (each pass: set IC, run the markers). `> 0` requires even `ni/nj/nk`. |
| `type` | string | yes | See the table below. Other → fatal. | Initial condition kind. |
| `regions_number` | int | `riemann-problem` | `>= 1`, otherwise fatal. Ignored by the other types (they fix the count). | Number of `[initial_conditions_region_N]` sections. |

### Accepted `type` values and their extra keys

| `type` | Model | Regions read | Extra `[initial_conditions]` keys (all real unless stated) and checks |
|--------|-------|--------------|-----------------------------------------------------------------------|
| `uniform` | any | 1 | `s`: signed relative amplitude of a deterministic hashed perturbation of `rho` and `p`. |
| `isentropic-vortex` | **euler** (fatal otherwise) | 1 (free stream) | `x0`, `y0`, `radius` (`> 0`), `strength`. |
| `riemann-problem` | any | `regions_number` | `regions_number`. Each region also needs extents (see the next section). Every cell must be covered by a region (fatal otherwise, at set time), and the **first** matching region wins. |
| `glm-pulse` | **mhd-ideal** | 1 | `pulse_axis` (string `x`/`y`/`z`, fatal otherwise), `pulse_center`, `pulse_width` (`> 0`), `pulse_amplitude`. |
| `divb-peak` | **mhd-ideal** | 1 | `peak_x0`, `peak_y0`, `peak_radius` (`> 0`), `peak_amplitude`. |
| `mhd-linear-wave` | **mhd-ideal** | 1 | `wave` (string `fast`/`alfven`/`slow`/`entropy`, fatal otherwise), `wave_angle` (degrees, x-y plane), `wave_amplitude`, `wavelength` (`> 0`). |
| `mhd-cpaw` | **mhd-ideal** | 1 (u, v, w, bx, by, bz must be 0, fatal otherwise) | `polarisation` (string `right`/`left`, fatal otherwise), `wave_angle`, `wave_amplitude`, `wavelength` (`> 0`), `b_par` (`/= 0`). |
| `mhd-vortex` | **mhd-ideal** | 1 (bx, by must be 0, fatal otherwise) | `x0`, `y0`, `radius` (`> 0`), `kappa`, `mu`. |
| `orszag-tang` | **mhd-ideal** | **none** | No extra keys. The unit-period Orszag–Tang state is hard-coded (`gamma` comes from `[physics]`). |
| `mhd-rotor` | **mhd-ideal** | 1 (ambient) | `x0`, `y0`, `r0`, `r1` (`0 < r0 < r1`), `rho_in` (`> 0`), `v0`. |
| `field-loop` | **mhd-ideal** | 1 | `x0`, `y0`, `loop_radius` (`> 0`), `loop_amplitude`. |
| `rotated-riemann` | any | 2 (states given in the **normal frame**: `u`, `bx` normal; `v`, `by` tangential) | `normal_x`, `normal_y` (not both 0), `interface_1`, `interface_2`, `period` (`> 0`), `interface_2_width` (`>= 0`). Requires `interface_1 < interface_2` and `interface_2 + interface_2_width < interface_1 + period`, fatal otherwise. |
| `shu-osher` | **euler** | 2 | `axis` (string `x`/`y`/`z`), `interface`, `rho_amplitude`, `rho_wavenumber`. Requires `rho_amplitude < region_2 r`. Region 1 applies where `s <= interface`; region 2 applies elsewhere, with `r + rho_amplitude * sin(rho_wavenumber * s)`. |
| `linear` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z` (inverse length). Every conservative variable is `q_1 (1 + g . x)`, `q_1` the state of region 1: the verification field of the ghost probe (issue #65 P0), not a flow. |
| `isothermal-atmosphere` | any | 1 (the state at the origin; its velocity and field uniform) | Needs `[gravity]` (fatal otherwise). Optional, both atmospheres: a relative pressure pulse $p(1 + A e^{-|\mathbf x - \mathbf x_0|^2/w^2})$, `pulse_amplitude` $A$ (default 0), `pulse_x`, `pulse_y`, `pulse_z` ($\mathbf x_0$), `pulse_width` $w$ (`> 0` when $A \ne 0$, fatal otherwise; M7-G2b). $\rho = \rho_1 \exp(\mathbf g\cdot\mathbf x / T_1)$, $p = \rho T_1$, $T_1 = p_1/\rho_1$: the hydrostatic equilibrium of the constant gravity ([#73](https://github.com/szaghi/adam/issues/73)). |
| `polytropic-atmosphere` | any | 1 (as above) | `polytropic_gamma` $G > 1$ (fatal otherwise); needs `[gravity]`. $p = p_1(\rho/\rho_1)^G$, $\rho = \rho_1 s^{1/(G-1)}$, $s = 1 + \frac{G-1}{G}\frac{\rho_1}{p_1}\mathbf g\cdot\mathbf x$; fatal where $s \le 0$ (the top of the atmosphere inside the domain). |
| `rayleigh-taylor` | any | 1 (the light layer: $T_l = p_1/\rho_1$, $p_1$ the pressure at the interface; its field uniform) | Needs `[gravity]` along one axis (fatal otherwise). `rt_height` $h_0$ (along $-\mathbf g$, $h = -\mathbf g\cdot\mathbf x/|\mathbf g|$), `rt_width` $\delta$, `rt_density_ratio` $R$ (heavy over light at equal pressure; $R < 1$ is the stable control), `rt_wavelength` $\lambda$, `rt_amplitude` $A$, `rt_seed_width` $w$; every key required, $\delta, R, \lambda, w > 0$ (fatal otherwise). Two isothermal layers joined by $T = a + b\tanh((h - h_0)/\delta)$, $a = (T_l + T_h)/2$, $b = (T_h - T_l)/2$, $T_h = T_l/R$, hydrostatic in closed form, at rest; seeded with the solenoidal single mode $v_h = A\cos(k x_t) f(h)$, $u_t = -\frac{A}{k}\sin(k x_t) f'(h)$, $f = e^{-((h - h_0)/w)^2}$, $k = 2\pi/\lambda$, $x_t$ the first axis other than $\mathbf g$'s (M7-G3; its linear growth rate: `verification/gravity/rt_oracle.py`). |
| `polynomial` | any | 1 | `gradient_x`, `gradient_y`, `gradient_z`, `degree` (an integer $n \ge 1$). Every conservative variable is `q_1 (1 + sum_d g_d x_d^n)`: the ghost probe's field for the high-order seam operators ([#68](https://github.com/szaghi/adam/issues/68) D1), not a flow. Fatal under a non-dimensional reference (`g_d` has dimension $L^{-n}$). |

A model mismatch (for example `glm-pulse` with `euler`) is fatal. The complete formulas are in the module header.

### [initial_conditions_region_N]

N runs over the regions the selected type reads (see the table above).

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `r`, `u`, `v`, `w`, `p` | real | yes | Primitive state (density, velocity, pressure). No positivity check. |
| `bx`, `by`, `bz` | real | `mhd-ideal` | Magnetic field. `psi` starts at 0. |
| `emin_x`, `emin_y`, `emin_z`, `emax_x`, `emax_y`, `emax_z` | real | `type = riemann-problem` only | Region box. A cell belongs to it iff `emin < centre <= emax` on every axis. |

---

## [time]

Every key is required, and none has a range check.

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `it_max` | int | yes | `> 0`: stop after `it_max` steps (iteration-driven). `<= 0`: time-driven run. |
| `time_max` | real | yes | End time (time-driven runs). The last step is clipped to hit it. |
| `CFL` | real | yes | CFL number (note the upper-case key). |

---

## [diagnostics]

| Key | Type | Req. | Meaning |
|-----|------|------|---------|
| `conservation_history_save` | int | yes | Cadence (steps) of `<output_basename>-conservation_history.dat` (volume integrals of every conservative variable). For MHD it is also the cadence of `-divb_history.dat` and of the `divb_tol` monitor. `<= 0` disables them. |
| `ghost_poison` | logical | no (`.false.`) | Verification instrument (issue #65): before each field write every ghost cell is set to NaN and refilled in the order of a stage (seams, exchange, boundary conditions), so a ghost the fill misses, or reads before its donor is written, is written as NaN instead of a stale value. Used by the ghost probe GP; it costs one pass over the ghosts per write. |
| `divb_flux` | logical | no (`.false.`) | MHD with GLM or EGLM ([#68](https://github.com/szaghi/adam/issues/68) A2): adds four columns to `-divb_history.dat` (`max_dflux`, `l1_dflux`, `seam_max_dflux`, `max_h_dflux_over_bmax`). D_flux is the divergence the cleaning acts on: that of the interface normal field $B_n = F_\psi/s$ read from the scheme's own $\psi$ face fluxes ($s = c_h^2$ GLM, $c_h$ EGLM), the coarse 2:1 faces taking the fine fluxes (the reflux), so its integral over a composite periodic domain is round-off. The last column is the scale-free $\max(h\,|D_\text{flux}|)/\max|\mathbf B|$. Each history row costs one residual evaluation of the committed state, on the idle RK stage buffers: fatal without GLM/EGLM or with a one-stage scheme. An inter-realm seam keeps each side's own flux; the step-0 row, written before the forest builds its register, has no reflux substitution and no seam column (as the centred one). |

---

## Forest manifest (multi-realm)

The manifest is a separate INI; every `[realm.N]` points at a full realm input as documented
above. Manifest errors use Fortran `error stop`, not `mpih%error_stop`.

| Section | Key | Type | Req. | Accepted / invalid | Meaning |
|---------|-----|------|------|--------------------|---------|
| `[forest]` | `realms_number` | int | yes | `>= 1`, otherwise fatal. Its presence with `>= 1` is what makes the file a manifest. | Number of realms. |
| `[realm.N]` (N = 1..realms_number) | `ini` | string (≤ 256 chars) | yes | See the note below the table. | Realm input. |
| `[forest.topology]` | `inter_realm_faces_number` | int | no | Default 0 (no coupling). | Number of glued faces. |
| `[forest.topology.face_N]` | `realm_a`, `realm_b` | int | yes | `1..realms_number`, otherwise fatal. | Realms on each side of the seam. |
| | `face_a`, `face_b` | string | yes | `+x -x +y -y +z -z`, also `X`/upper case and suffix forms `x+`, `X-`, ... Other → fatal. | Glued faces. |
| | `coupling` | string | no | `mirror` (default), `periodic`, `interpolate` (each in lower, Capitalised or UPPER case). Other → fatal. | Declared coupling kind. **Stored, but no code dispatches on it**: `periodic` and `interpolate` are "reserved, not implemented" and behave as the pass-through seam. |
| | `coupling_cadence` | string | no | `end_of_step` (default, α) or `stage_coincident` (β), in the same three case forms. Other → fatal. β requires both realms to have equal `[runge_kutta].scheme`, equal `nv` (same physics/divergence control) and equal stages per step (fatal otherwise). Conflicting cadences between the same realm pair are fatal. | Seam fill once per step (α) or at every RK stage (β). For a 1:1 same-resolution seam use β. |

Note on `[realm.N] ini`: the path is **used verbatim**, so it resolves against the process working directory, **not**
the manifest's directory (see the inconsistencies section).

---

## Complete example 1: Euler, Sod shock tube (verbatim)

The file is `src/tests/flume/verification/sod/sod-x.ini`, copied without changes.

```ini
; FLUME verification V1 (issue #35, section 11): Sod shock tube along x, 1-D (the other directions are null).
; 64 blocks, four along x (200 cells); the 8 sibling groups split 32/32 over two ranks (the library never splits a
; sibling group, so one refinement level would leave a rank empty). WENO-5 characteristic, SSP-33, CFL 0.5, t = 0.2.
; check.sh compares density with the exact Riemann solution (L1 bound) and the three directions with each other (bitwise).
[IO]
output_basename        = sod-x
it_save                = 100000
restart                = .false.
restart_basename       = sod-x-restart
restart_save           = 0
residuals_save         = 1
save_memory_status     = .false.
save_residual_fields   = .false.
save_auxiliary_fields  = .false.
save_curl_fields       = .false.
save_divergence_fields = .false.
save_gradient_fields   = .false.
save_laplacian_fields  = .false.

[grid]
ni     = 50
nj     = 4
nk     = 4
ngc    = 3
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
null_x = .false.
null_y = .true.
null_z = .true.

[amr]
max_level      = 2
ratio          = 8
iu_ref_levels  = 2
i_prune        = 0
j_prune        = 0
k_prune        = 0
l_prune        = -1
frequency      = 0
iters          = 1
markers_number = 0

[field]
nv = 5

[runge_kutta]
scheme = runge-kutta-ssp-33

[weno]
scheme              = weno-u-5
ror_number          = 0
ror_threshold       = 0.9
ror_vars_number     = 0
enable_ror_stats    = .false.
ib_reduction_extent = 0
ib_reduced_order    = 2

[linear-algebra]
smoothing         = gauss-seidel
iterations_init   = 3
iterations_coarse = 10
iterations_fine   = 3
iterations        = 10
tolerance         = 1.e-20

[fdv]
fdv_scheme = fd
fdv_order  = 2

[solids]
number    = 0
n_eikonal = 0

[slices]
slices_number = 0

[physics]
physical_model = euler
cp             = 1040.004
cv             = 742.86

[numerics]
scheme_space             = weno
reconstruction_variables = characteristic
reflux                   = .true.

[bc_x_min]
type = extrapolation
[bc_x_max]
type = extrapolation
[bc_y_min]
type = extrapolation
[bc_y_max]
type = extrapolation
[bc_z_min]
type = extrapolation
[bc_z_max]
type = extrapolation

[initial_conditions]
type           = riemann-problem
regions_number = 2
amr_iterations = 0

[initial_conditions_region_1]
r      = 1.0
u      = 0.0
v      = 0.0
w      = 0.0
p      = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 0.5
emax_y = 1.0
emax_z = 1.0

[initial_conditions_region_2]
r      = 0.125
u      = 0.0
v      = 0.0
w      = 0.0
p      = 0.1
emin_x = 0.5
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0

[time]
it_max   = -1
time_max = 0.2
CFL      = 0.5

[diagnostics]
conservation_history_save = 1
```

### Annotations for the Sod input

- **Grid and block count.** With `ratio = 8` and `iu_ref_levels = 2` the domain is split into `8^2 = 64` blocks of
  `50 x 4 x 4` cells. The run is 1-D: `null_y` and `null_z` make y and z inactive, so there are 4 blocks along x and
  200 cells in total along x.
- **`it_save = 100000`.** This is effectively "first and last step only", because the last step is always saved.
- **`restart_save = 0` and `residuals_save = 1`.** Restarts are disabled, and the residual history is written every
  step.
- **Keys FLUME ignores.** `[field] nv`, `[linear-algebra]`, `fdv_scheme` and the `ror_*`/`ib_*` keys of `[weno]`
  are ignored by FLUME. All of them except `[field]` must still be present.
- **Riemann-problem regions.** The regions tile x at 0.5: region 1 owns `0 < x <= 0.5` and region 2 owns
  `0.5 < x <= 1`. The domain minimum corner is covered because `emin` is the domain minimum and cell centres are
  strictly larger.

---

## Complete example 2: ideal MHD with GLM, Orszag–Tang vortex with a 2:1 AMR box (verbatim)

There is **no committed INI** under `src/tests/flume/verification/mhd/orszag-tang/`. `check.sh` generates the input
with `make_orszag_tang.py` from `verification/vortex/vortex-n064.ini`. The file below is copied without changes from
the generated (git-ignored) `src/tests/flume/verification/mhd/orszag-tang/work-adam_flume_cpu-np2-amr/orszag-tang.ini`
(written 2026-09-29, `--cells 32 --refine-box 0.25 0.25 0.75 0.75`). It is written by Python `configparser`, which
explains the section order and the spacing.

```ini
[IO]
output_basename = orszag-tang
it_save = 1000000
restart = .false.
restart_basename = vortex-n064-restart
restart_save = 0
residuals_save = 1
save_memory_status = .false.
save_residual_fields = .false.
save_auxiliary_fields = .false.
save_curl_fields = .false.
save_divergence_fields = .false.
save_gradient_fields = .false.
save_laplacian_fields = .false.

[grid]
ni = 8
nj = 8
nk = 4
ngc = 3
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
null_x = .false.
null_y = .false.
null_z = .true.

[amr]
max_level = 3
ratio = 8
iu_ref_levels = 2
i_prune = 0
j_prune = 0
k_prune = 0
l_prune = -1
frequency = 0
iters = 1
markers_number = 1

[field]
nv = 9

[runge_kutta]
scheme = runge-kutta-ssp-54

[weno]
scheme = weno-u-5
ror_number = 0
ror_threshold = 0.9
ror_vars_number = 0
enable_ror_stats = .false.
ib_reduction_extent = 0
ib_reduced_order = 2

[linear-algebra]
smoothing = gauss-seidel
iterations_init = 3
iterations_coarse = 10
iterations_fine = 3
iterations = 10
tolerance = 1.e-20

[fdv]
fdv_scheme = fd
fdv_order = 2

[solids]
number = 0
n_eikonal = 0

[slices]
slices_number = 0

[physics]
physical_model = mhd-ideal
cp = 2.5
cv = 1.5

[numerics]
scheme_space = weno
reconstruction_variables = characteristic
reflux = .true.

[bc_x_min]
type = periodic

[bc_x_max]
type = periodic

[bc_y_min]
type = periodic

[bc_y_max]
type = periodic

[bc_z_min]
type = extrapolation

[bc_z_max]
type = extrapolation

[initial_conditions]
type = orszag-tang
amr_iterations = 1

[time]
it_max = -1
time_max = 0.5
CFL = 0.4

[diagnostics]
conservation_history_save = 1

[mhd]
divergence_control = glm
divb_tol = 0.0
divb_error = .false.
rho_floor = 0.0
p_floor = 0.0
glm_ch = 4.0
glm_alpha = 0.18
glm_damping_length = 1.0
glm_ch_check = error

[amr_marker_1]
mode = 1
geo_type = primitive-box
delta_type = max
delta_fine = 0.0
delta_coarse = 0.0
box_xmin = 0.25
box_ymin = 0.25
box_zmin = -1e+30
box_xmax = 0.75
box_ymax = 0.75
box_zmax = 1e+30
target_level = 3
```

### Annotations for the Orszag–Tang input

- **Physics.** `cp/cv = 5/3`. `divergence_control = glm` gives `nv = 9` (`psi`); `[field] nv = 9` is ignored.
- **Time step and cleaning speed.** `glm_ch = 4` enters the time step. `glm_ch_check = error` stops the run if a fast
  wave ever outruns `c_h`.
- **Monitors and floors.** `divb_tol = 0` disables the div(B) monitor (on the normalised D_flux, which needs
  `[diagnostics] divb_flux`), but `-divb_history.dat` is still written every step. Both floors are 0, so a non-positive density or pressure is fatal.
- **Grid.** The base grid is `4 x 4` blocks (`iu_ref_levels = 2` on a quasi-2-D octree, with z null) of `8 x 8` cells,
  i.e. 32² cells. `ni` and `nj` are even, as `amr_iterations = 1` requires.
- **Refinement.** `amr_iterations = 1` runs one pass of the single marker. The marker is `mode = 1` +
  `geo_type = primitive-box`: it refines the blocks whose centroid lies in `[0.25, 0.75]^2` to level 3, one level
  above `iu_ref_levels`, capped by `max_level = 3`. `delta_type`, `delta_fine` and `delta_coarse` are required but
  unused for a box.
- **Boundaries and IC.** x and y are periodic (both faces of each axis). The IC has no region sections.
- **Restart basename.** `restart_basename = vortex-n064-restart` is inherited from the base input. It is harmless here
  because `restart_save = 0`.

## Derived output fields

Both groups are written to the XH5F output only when `[IO] save_auxiliary_fields = .true.`. The key is **required**. The auxiliary variables are recomputed on the host from the saved `q`,
including ghost cells. The derived MHD group is written only for the MHD models.

`save_simulation_data` calls `save_xh5f(with_ghost=.true.)`, so **field files include the ghost cells**. In a forest
the inter-realm seam ghosts are refilled before the write, then the intra-realm ghosts and the boundary conditions
(the order of a Runge-Kutta stage), so the written ghosts are those the stencils read; the step-0 file, written while
each realm initialises and before the forest connects them, holds unfilled seam ghosts.

### Auxiliary variables

| name | index | Euler | MHD |
|---|---|---|---|
| `rho` | `IA_R` | $\rho$ | $\rho$ |
| `u`, `v`, `w` | `IA_U..IA_W` | $(\rho u)/\rho$, … | same |
| `p` | `IA_P` | $(\gamma-1)\big(E-\tfrac12\rho\vert \mathbf u\vert ^2\big)$ | $(\gamma-1)\big(E-\tfrac12\rho\vert \mathbf u\vert ^2-\tfrac12\vert \mathbf B\vert ^2\big)$ |
| `T` | `IA_T` | $p/(\rho R)$, with $R=c_p-c_v$ | same |
| `H` | `IA_H` | $(E+p)/\rho$ | $(E+p+\tfrac12\vert \mathbf B\vert ^2)/\rho$, which includes the magnetic pressure |
| `a` | `IA_A` | $\sqrt{\gamma p/\rho}$ | $\sqrt{\gamma p/\rho}$ (the **sound** speed, not the fast speed) |
| `Bx`, `By`, `Bz` | `IA_BX..IA_BZ` | — | copies of $B_x, B_y, B_z$ (capitalised to avoid clashing with the conservative `bx, by, bz`) |

The conservative variables are named `r, ru, rv, rw, rE [, bx, by, bz [, psi]]`, and the residuals `dq_<name>`.

### MHD derived fields

| name | definition | where |
|---|---|---|
| `pt` | $p + \tfrac12\vert \mathbf B\vert ^2$ (total pressure, rationalised units) | every cell, ghosts included |
| `beta` | $2p/\vert \mathbf B\vert ^2$, or `huge(1._R8P)` $\approx 1.8\times10^{308}$ where $\vert \mathbf B\vert ^2 = 0$ | every cell |
| `bmag` | $\vert \mathbf B\vert  = \sqrt{B_x^2+B_y^2+B_z^2}$ | every cell |
| `divb` | $\sum_d w_d\,\delta_d B_d$ (signed; see below) | interior only, **zero in the ghost cells** |

Here $|\mathbf B|^2$ is summed in plain $x,y,z$ order, not through `mhd_sum3`.

**Discretisation of `divb`.** Each $\delta_d$ is the library centred first derivative,
$\delta B = \frac{1}{\Delta}\sum_{m=1}^{s}c_{m,s}\,(B_{i+m}-B_{i-m})$.
The coefficients are:

| $s$ | order | $c_{m,s}$ |
|---|---|---|
| 1 | 2 | $\tfrac12$ |
| 2 | 4 | $\tfrac{8}{12}, -\tfrac{1}{12}$ |
| 3 | 6 | $\tfrac{45}{60}, -\tfrac{9}{60}, \tfrac1{60}$ |
. The half stencil is $s$ = `fdv_half_stencils(1)` $=$ `[fdv] fdv_order / 2`,
so the formal order equals `fdv_order`. Every FLUME test input sets `fdv_order = 2`, which gives
$(B_{i+1}-B_{i-1})/(2\Delta x)$.

Two further rules apply:
- Null directions are weighted $w_d = 0$.
- The scheme is **always FD-centred**, even when `[fdv] fdv_scheme = fv`.

The operator is the same one used by the div(B) history kernels, `compute_divb_norms_dev`. Those kernels report $\max|\nabla\!\cdot\!\mathbf B|$,
$\sum|\nabla\!\cdot\!\mathbf B|\,\Delta V$ and a seam-band maximum, and cap $s$ at `HS_MAX = 3`. The CPU history also
requires $s\le n_{gc}$. The output path has no such cap or check.
