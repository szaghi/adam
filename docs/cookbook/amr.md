# Adaptive mesh refinement

Recipes for refining the grid: a static region at initialization, markers that follow the solution, and FLUME's
runtime regridding. The keys are in the [library input reference](/reference/library-input#amr); the tree and the 2:1
faces are described in [Concepts](/concepts/architecture).

## Which marker in which application

| `[amr_marker_N] mode` | FLUME | PRISM |
|---|---|---|
| `1` geometric: `primitive-box` | yes | yes |
| `1` geometric: `solid` (an immersed solid's surface) | yes | does nothing |
| `2` gradient of a variable | yes | does nothing |
| `3` total variation | refused | yes (the coils' current density) |
| `4` Löhner estimator | yes, at runtime regrids only | does nothing |
| `5` Jeans length | yes, with `[gravity] self_gravity = .true.` | does nothing |

Markers act during the `[initial_conditions] amr_iterations` passes (set the initial condition, run one AMR update of
up to `[amr] iters` sweeps, repeat) in both applications, and in FLUME also at every runtime regrid.

## AMR1. Refine a fixed region at initialization

**Problem.** Put finer blocks where you know the solution needs them, for the whole run.

**Input** (`src/tests/flume/regression/sod-amr/input.ini`):

```ini
[amr]
max_level      = 3
ratio          = 8
iu_ref_levels  = 2
i_prune        = 0
j_prune        = 0
k_prune        = 0
l_prune        = -1
frequency      = 0
iters          = 1
markers_number = 1

[amr_marker_1]
mode         = 1             ; AMR_GEO (geometric marker)
geo_type     = primitive-box
delta_type   = max           ; unused by the box marker; required by the common marker parse
delta_fine   = 0.0
delta_coarse = 0.0
box_xmin     = 0.75          ; the blocks whose centroid is in x > 0.75 are refined once more
box_ymin     = 0.0
box_zmin     = 0.0
box_xmax     = 1.0
box_ymax     = 1.0
box_zmax     = 1.0
target_level = 3

[initial_conditions]
amr_iterations = 1
```

**What you get.** Two uniform levels everywhere (`iu_ref_levels = 2`), and the blocks whose centroid lies in the box
refined to level 3. The faces between level 2 and level 3 are **2:1 faces**: the solver fills the fine ghosts by
interpolation (`[amr] seam_ghost_fill`) and the coarse ghosts by restriction (`seam_restriction`). FLUME (`[numerics] reflux`, on by default), and PRISM
with the FV operators, then correct the coarse cells along the face after every step with the fine fluxes (the
**reflux**, Berger–Colella), so the conserved totals do not drift; PRISM's FD operators do not reflux. A run with a 2:1 face always takes the forest's staged path.

PRISM refines the same way (`src/tests/prism/regression/rmf-amr/input.ini`: the `x < 0` half refined one level).

**Reference.** [`[amr]`, `[amr_marker_N]`](/reference/library-input#amr); the
[seam-conservation verification](/applications/flume/verification/euler/amr-conservation).

## AMR2. Follow the solution: FLUME runtime regridding

**Problem.** Refine where a feature is and coarsen behind it, during the run.

**Input** (`src/tests/flume/regression/vortex-regrid/input.ini`):

```ini
[amr]
max_level      = 3
ratio          = 4
iu_ref_levels  = 1
i_prune        = 0
j_prune        = 0
k_prune        = 0
l_prune        = -1
frequency      = 5
regrid_prolongation = conservative
iters          = 2
markers_number = 1

[amr_marker_1]
mode         = 4             ; Loehner estimator (issue #74)
delta_type   = max           ; unused by the Loehner marker; required by the common marker parse
delta_fine   = 0.0
delta_coarse = 0.0
field        = 1             ; 1 conservative, 2 auxiliary variables
var          = 1
refine_tol   = 0.3
derefine_tol = 0.1
floor        = 0.0
buffer       = 2
```

**What you get.** After every fifth step the markers vote, blocks move by up to `iters` levels between
`iu_ref_levels` (never coarsened below) and `max_level`, and refined children are filled from their parent by
`regrid_prolongation` (`conservative`: the children average to the parent). The Löhner marker does not act at
initialization; add a gradient marker (mode 2) to refine the initial grid too, as in the
[FLUME runtime-AMR guide](/applications/flume/input#setting-up-runtime-amr). Runtime regridding works on both FLUME
backends; PRISM reads `frequency` and refines at initialization only. On a multi-realm forest each realm regrids by its
own `[amr]` and the forest balances the levels across the seams ([forest cookbook R12](/cookbook/forest#r12-runtime-amr-on-a-forest)).

**Reference.** [Setting up runtime AMR](/applications/flume/input#setting-up-runtime-amr);
[runtime AMR verification](/applications/flume/verification/runtime-amr/).

## AMR3. Refine on the Jeans length (self-gravity)

**Problem.** Keep a collapsing cloud resolved by a fixed number of cells per Jeans length.

**Input** (the runtime-AMR collapse written by `src/tests/flume/verification/sinks/make_si.py`):

```ini
[amr]
frequency            = 4
regrid_prolongation  = conservative-state
markers_number       = 1

[amr_marker_1]
mode         = 5
delta_type   = max
delta_fine   = 0.0
delta_coarse = 0.0
jeans_cells  = 8.0
```

**What you get.** A block refines where some interior cell resolves its Jeans length with fewer than `jeans_cells`
cells, and coarsens where every cell has more than `jeans_cells_derefine` (default `2 jeans_cells`); the sound speed is
$c_s^2 = dP/d\rho$ (isothermal and barotropic gases included). It needs `[gravity] self_gravity = .true.`. A collapse
on runtime AMR needs `regrid_prolongation = conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)).

**Reference.** [`[amr_marker_N]`](/reference/library-input#amr-marker-n-n-1-markers-number);
[self-gravity under runtime AMR](/applications/flume/verification/self-gravity/runtime-amr).

## Pitfalls

- **Write all seven tree keys** (`max_level`, `ratio`, `iu_ref_levels`, `i_prune`, `j_prune`, `k_prune`, `l_prune`).
  The tree reader ignores a missing key's error and passes the same buffer to every key: a missing `max_level`
  silently takes the value of `ratio`, and a missing `ratio` is garbage. Keep at least one `*_prune` negative.
- **Even cells per refined axis.** With markers and `max_level > iu_ref_levels`, every refined non-null axis needs an
  even `ni`/`nj`/`nk` of at least `2 ngc` cells; the library checks it at initialization and stops, in both apps
  (`check_amr_block_cells`; [#39](https://github.com/szaghi/adam/issues/39), [#66](https://github.com/szaghi/adam/issues/66)).
- **`ratio`**: `8` for octrees, `4` for quadtrees (one null axis), `2` for 1-D; any other value leaves the block sizes
  undefined with no error.
- **Marker keys that are unused are still required**: every marker reads `delta_type`, `delta_fine` and
  `delta_coarse`.
- **The 6-point seam fills** (`quintic`, `weno6`, `quintic-mp`, and the `lagrange6`/`weno6` restrictions) need at least 6
  cells per refined axis in every block.

**Reference.** [`[amr]`](/reference/library-input#amr) (the trap note under the table).
