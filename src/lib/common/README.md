# ADAM Common Library — `src/lib/common`

The common library is the **portable, physics-agnostic core** of the ADAM SDK. Every GPU backend (`fnl`, `nvf`, `gmp`) extends these objects rather than replacing them; every solver application (`nasto`, `prism`, `chase`, …) composes them. No GPU-specific code lives here — the entire directory compiles with any standard Fortran 2003+ compiler.

The aggregate entry point `adam_common_library` re-exports all modules below; a single `use adam_common_library` statement in application code exposes the complete API.

---

## Contents

- [AMR data structures](#amr-data-structures)
- [Grid and field](#grid-and-field)
- [Communication and maps](#communication-and-maps)
- [Numerical methods — spatial](#numerical-methods-spatial)
- [Numerical methods — temporal](#numerical-methods-temporal)
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

## Numerical methods: spatial

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

## Numerical methods: temporal

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

The next three integrators are library objects that `realm_object` declares but does **not** initialize (the calls
are commented out); no application can use them today, and PRISM refuses them at init (issue #55).

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

### `adam_flail_object` — FLAIL, the elliptic solver library (issue #84)

One entry point for the elliptic solves of every app, the realm component `flail` (`self%flail`). It is built on the two
solvers below, verified in #81/#87: the FFT solver on uniform grids and the multigrid-preconditioned CG solver on
AMR grids. The legacy smoothers (next section) retire once PRISM moves here (#84, E3-E4).

| Procedure | Role |
|-----------|------|
| `setup(grid, field, isolated, backend[, option_name][, bc][, dg_order])` | Choose the solver for the current grid (`FLAIL_BACKEND_AUTO`: the FFT solver on a uniform grid, the multigrid otherwise; `_FFT`, `_MULTIGRID`); call it on the initial grid and after every regrid. `bc(2,3)` (E1c): the kind of every face, forwarded (`FLAIL_BC_PERIODIC`, `_DIRICHLET`, `_NEUMANN`, `_DIRICHLET_GHOST`, re-exported here); without it, periodic or isolated as in E0; isolated with `bc` is fatal. `dg_order` (E2a): the wide D(G) of order 2, 4, 6, 8, on the FFT path only (fatal on AMR: not defined across 2:1 faces yet, #84 E7) |
| `solve(field, q, iv, factor, phi, energy[, ghost])`, `solve_packed(field, buf, factor, phi, energy[, ghost])` | Forward to the chosen solver; `energy = sum f phi dV`; `ghost` the boundary value of each boundary face (`mg%nbnd`, `mg%bnd_kind`, `mg%bnd_dir`, `mg%bnd_x`), multigrid path only (fatal on the FFT path, which solves the homogeneous problem) |
| `iterations()`, `is_multigrid()`, `destroy()` | The last FCG iteration count (0 for the FFT solver), the solver in use, cleanup |

A vector whose components take different kinds (PRISM's vector potential at a PEC face: Neumann for the normal
component, Dirichlet for the others) takes one FLAIL per set of kinds: the operators differ (the FFT extensions, the
singularity), and the instances coexist (`test_poisson_mg --flail`, bitwise). The mapping from an app's boundaries to
the kinds belongs to the app (PRISM's lands with its migration, #84 E3).

A backend chooses its FFT engine by allocating the `fft` component before `setup` (FLUME's FNL backend allocates the
cuFFT one). The boundary kinds per face, `FLAIL_BC_PERIODIC | _DIRICHLET | _NEUMANN | _DIRICHLET_GHOST`, live in
`adam_flail_parameters` (a module of their own: the solvers use them, and the facade uses the solvers). Phases and
status: [issue #84](https://github.com/szaghi/adam/issues/84).

### `adam_flail_legacy_object` — legacy FLAIL smoothers (retired in issue #84, E4)

PRISM's elliptic smoothers until they move to the FLAIL library (E3); the realm component is `flail_legacy`, and it
owns the `[linear-algebra]` section. No new code should use it.

Iterative solver and smoother infrastructure for implicit operators (e.g., elliptic sub-problems). Provides a common interface over several smoothing strategies:

| Constant | Method |
|----------|--------|
| `SMOOTHING_MULTIGRID` | Geometric multigrid V-cycle |
| `SMOOTHING_GAUSS_SEIDEL` | Gauss–Seidel iteration |
| `SMOOTHING_SOR` | Successive Over-Relaxation |
| `SMOOTHING_SOR_OMP` | Parallel SOR (OpenMP) |

Convergence is controlled by `tolerance` (default 10⁻⁶) and `iterations`. Configuration section: `[linear-algebra]`.

### `adam_poisson_fft_object` — FFT Poisson solver (uniform grids: periodic, Dirichlet/Neumann faces, or isolated)

Solves `L_h phi = factor (f - mean f)` (periodic) or `L_h phi = factor f` with `phi -> 0` at infinity (isolated) for
one variable `f` of a conservative array `q(nv,i,j,k,b)`, with `L_h` the second-order 7-point Laplacian (5-point in
2-D; a null direction has no derivative) ([#81](https://github.com/szaghi/adam/issues/81)).

- **Periodic (P2):** each Fourier mode is divided by the eigenvalue of `L_h` itself,
  `lambda_k = sum_d (2 cos(2 pi k_d/N_d) - 2)/dx_d^2`, so the discrete residual is at round-off.
- **Isolated (P4), Hockney–Eastwood:** the source sits in a zero-padded array of `m = 2 (n + ngc)` cells per direction
  and is convolved (one r2c, a product, one c2r) with the **lattice** Green's function of `L_h`, `h^2 G(n)`, whose
  spectrum is computed once. The circular convolution equals the free-space one on every cell the scatter reads,
  ghosts included. `G` comes from `adam_lattice_green_library` (below): its residual is about `1e-13` per offset, so
  the residual of `phi` is that summed over the source (`2e-11` of `factor f` at `64^3`), not round-off.

- **Dirichlet/Neumann faces (issue #84, E1a):** `bc(2,3)` gives each face a homogeneous kind; an axis is periodic, or
  each of its faces is Dirichlet (zero at the face) or Neumann (zero normal derivative), possibly different at the two
  ends. The source is extended odd about a Dirichlet face and even about a Neumann one, over `2n` cells when the faces
  agree and `4n` when they differ, and the periodic transform of the extension solves the 7-point operator with those
  conditions exactly (round-off residual, face placement of the wall: second order). Without a Dirichlet face the source
  mean is dropped, as on a periodic grid. The energy sums the interior only.
- **The wide operator D(G) (issue #84, E2a):** `dg_order = p` (2, 4, 6, 8; 0 the compact Laplacian, the default) solves
  PRISM's elliptic operator, the order-p centred first derivative composed with itself per axis (a stencil of
  `2p + 1` cells per axis), eigenvalue `-sum_d s_p(theta_d)^2 / dx_d^2`, `s_p(theta) = sum_j 2 a_j sin(j theta)`.
  Every centred derivative annihilates the checkerboard, `s_p(pi) = 0`, set exactly (`sin(pi)` is `1e-16`): the modes
  with each `theta_d` in {0, pi} are in the kernel and dropped, Dirichlet faces or not (G of them is zero). The
  extensions are reflections, which the wide stencil commutes with, so the kinds per face hold, the scatter filling
  `p` ghost layers; isolated boundaries take the compact operator only. Why a direct solve: on the prototype
  (`src/tests/poisson/wide_dg_proto.py`) a compact preconditioner, even inverted exactly, needs O(N) PCG iterations
  for D(G) (86 → 677 at p = 2, N = 32 → 256: the checkerboard branch it does not see), while the wide p = 2 operator
  (2^d compact Laplacians on the sub-lattices of spacing 2h) is spectrally equivalent (17–30 iterations, independent of
  N). PRISM's zero-ghost Dirichlet face is not a reflection: `wide_dg_zero_ghost.py` measures it (#84 E2b).
- **Zero-ghost faces (issue #84, E2b):** `FLAIL_BC_ZERO_GHOSTS` (both faces of an axis) is PRISM's elliptic Dirichlet
  face, every ghost layer zero. With the compact operator it is the ghost-Dirichlet kind (one layer read; the multigrid
  maps it the same way). With the wide D(G) the operator is the stencil truncated, symmetric and nonsingular, and no
  transform diagonalises it: rank 0 solves it on the gathered array by conjugate gradients (`tolerance` 1e-12 on the
  max-norm residual, `iterations`), preconditioned by the wide p = 2 operator with the same faces inverted exactly.
  Along a zero-ghost axis that operator splits into the two parity sub-lattices, each the compact Laplacian at spacing
  2h with zero at the ghost centres; the `2^z` sub-problems are FFT solves by sub-solvers of the solver's own type
  (cuFFT on FNL). Spectrally equivalent: 1, 17, 24, 30 iterations at p = 2, 4, 6, 8, independent of the grid and of
  the companion axes' kinds (`wide_dg_zero_ghost.py`, `wide_dg_e2b.py`). The ghost layers returned carry the faces'
  rules. A zero-ghost face against another kind on one axis is fatal (a reflection would couple the parity classes).

All are symmetric linear maps of the source, which is what exact energy and momentum conservation of self-gravity
needs (FLUME [models](/applications/flume/models#self-gravity)).

| Procedure | Role |
|-----------|------|
| `initialize(grid, field, isolated, bc, dg_order)` | Check the scope, size the solver, make the plans (and, isolated, the kernel's spectrum) on rank 0; `bc(2,3)` optional, consistent with the grid's periodicity; `dg_order` the wide D(G) (E2a) |
| `solve(field, q, iv, factor, phi, energy)` | Gather `q(iv)` to rank 0, transform, scatter `phi(1-ngc:n+ngc, ..., b)` with every ghost layer filled from the global solution; optionally `energy = sum f phi` |
| `solve_packed(field, buf, factor, phi, energy)` | The same from the interior source already packed (block-major, `i` fastest: what a device kernel produces) |
| `initialize_blocks(nblk, nijk, ngc, dxyz, level, null_xyz, isolated, dirichlet, bc)` | Size the solver for an explicit uniform block lattice, no grid object and no scope check (the caller owns it): the multigrid's coarsest level; `dirichlet` (#87) zero at the first ghost centres, by the odd extension to `2 (n + 1)` cells per active direction (the type-I sine transform: `FLAIL_BC_DIRICHLET_GHOST` on every axis); `bc` the kinds per face (E1a) |
| `solve_blocks(nbl, coords, buf, factor, gw, out, energy)` | Solve on an explicit list of this rank's blocks (integer coordinates, packed source), returning each block's potential with `gw` ghost layers; `initialize`/`solve_packed` are wrappers over these two |
| `destroy()` | Free the arrays and the plans |

- **Data path:** `MPI_Gatherv` of the block interiors (with their integer block coordinates) to rank 0, which holds the
  whole grid, solves and `MPI_Scatterv`s every block with its ghosts, so no ghost exchange follows. Correct for any rank
  count, ranks without blocks included; **not scalable**.
- **Scope (fatal otherwise):** every block at one level and the blocks covering the domain. Without `bc`, every active
  direction periodic; one periodic face against a non-periodic one, or a ghost-Dirichlet face against a face kind, on
  one axis are fatal. Isolated: 3-D, every direction non-periodic, cubic cells (the kernel is the unit lattice's). AMR
  grids: `adam_poisson_mg_object` (below).
- **Build:** `-D_FFTW` and `-lfftw3` (FFTW 3, `FFTW_ESTIMATE` plans: deterministic). Without an FFT library,
  `initialize` stops with "built without an FFT library". The FNL twin is `adam_fnl_poisson_fft_object` (cuFFT).
  Every mode builds it since the FLAIL facade (#84, E0) is a realm component; only the `*-fft` modes link an FFT library.
- **Verification:** SV-2 (periodic) and SV-6 (isolated), `src/tests/poisson/test_poisson_fft.F90` (`--isolated`)
  driven by `src/tests/flume/verification/selfgravity/poisson.sh`; the face kinds by
  `src/tests/poisson/test_flail_bc.F90` (E1a: per-axis discrete eigenfunctions of every kind mix, the ghost layer
  included, at round-off on 1-3 ranks; the D/D order against the continuum 1.99; refusals).

### `adam_poisson_mg_object` — Poisson solver on AMR grids (multigrid-preconditioned CG)

Solves `A phi = factor V f - B g` on the leaf cells of a 2:1-balanced octree (quadtree with `null_z`), periodic
([#81](https://github.com/szaghi/adam/issues/81), P5a) or with a boundary kind per face (#87, #84 E1b), the mean of
the right-hand side removed when `A` is singular. `A` is the composite second-order Laplacian, V-weighted
(finite-volume form): the 7-point stencil between leaves of one level and, across a 2:1 face, the **mimetic mean
flux** — every fine cell on a coarse cell's face gets the flux `(2/3) w (mean of the face's fine cells - phi_C)`, `w`
its face area over its spacing, and the coarse cell their sum. That coupling is the hybrid mimetic finite-difference
flux of the coarse cell with its face unknowns eliminated: exact for linear `phi` face by face (**consistent**: bounded
truncation, `phi` second order, the 2:1 face gradient first order) and **symmetric** (`phi = K f` with `V K`
symmetric: the exact energy exchange of self-gravity). Tomida & Stone's flux matching is consistent but not symmetric;
its symmetrisation is symmetric but inconsistent (O(1/h) truncation, a non-converging seam gravity): both measured in
`selfgravity/seam_operators.py`.

| Procedure | Role |
|-----------|------|
| `initialize(grid, field, tolerance, max_iterations, coarse, isolated, bc)` | Gather the leaf list (every rank), derive the level hierarchy, assemble the operator, the ghost fill and the V-cycle operators from the geometry; size the coarsest-level FFT solver, of the dynamic type of `coarse` when given (FLUME's FNL backend passes its cuFFT solver), FFTW otherwise; `isolated` also sizes the ghost-to-ghost Green's function of the James step (an isolated FFT solver of the same type on the coarsest box widened by the ghost layer); `bc` the kind of every face (E1b; default ghost-Dirichlet on the bounded axes) |
| `solve(field, q, iv, factor, phi, energy, ghost)` | FCG on `A`, preconditioned by one V-cycle, from the previous solution when the grid has not changed; `phi` returned with every ghost filled; optionally `energy = sum V f phi`; `ghost` the boundary value of each boundary face, of its kind (below); isolated: the James step (no `ghost`) |
| `solve_packed(field, buf, factor, phi, energy, ghost)` | The same from the packed interior source |
| `symmetry_defect(x, y)`, `apply_a(x, y)` | Verification: `|x.Ay - y.Ax| / (|x| |Ay|)` and `y = A x` on the leaves |
| `apply_b(g, y)`, `apply_bt(x, s)` | `y = B g` (the boundary values onto the boundary rows) and `s = B^T x`; the boundary faces are `nbnd`, `bnd_row`, `bnd_w` (the weight in `B`), `bnd_kind`, `bnd_dir` (the outward axis, signed), `bnd_x` (where the value lives: the ghost or the face centre) |
| `destroy()` | Free everything |

- **Solver:** flexible CG (Notay's FCG(1)), stopping on the relative max-norm residual: `tolerance` (default `1e-10`)
  or, when the source's mean dominates it, its round-off floor `100 eps max|factor f| / max|factor (f - mean f)|` (a
  `1e-6` perturbation of a uniform density is known to `1e-10` of itself; without the floor every solve ran to the cap,
  #81 P5b); `tolerance_used` keeps the one applied, and reaching `max_iterations` (200) above it prints a warning.
  The V-cycle runs over levels holding the leaves and a **shadow** block per refined node (ADAM keeps leaves only);
  red-black Gauss-Seidel (2 sweeps down, 2 up), restriction by sum, piecewise-constant prolongation, the coarsest level
  (the minimum leaf level: the whole domain) solved by `adam_poisson_fft_object%solve_blocks`. Diagnostics after a
  solve: `iterations`, `residual`, `factor_mean`.
- **Ghosts:** the first layer per face from the solution (across a 2:1 face, the value whose one-sided difference is
  the operator's own face gradient), deeper layers copying the first.
- **Distributed, assembled:** CSR rows over global cell keys, each operator with a halo-exchange plan
  (`MPI_ALLTOALLV`); every rank owns its leaves' rows and the shadows whose first leaf it owns. Correct for any rank
  count (results equal to round-off across rank counts); not tuned for speed.
- **Boundary faces (#87 P1; per face, #84 E1b):** an active axis is periodic or bounded; each face of a bounded axis
  has a kind (`adam_flail_parameters`), its row homogeneous (symmetric) and its value `g` on the right-hand side
  (`A phi = factor V f - B g`), `w` the face area over the spacing:

  | Kind | `g` | Row | `B` | Ghost |
  |------|-----|-----|-----|-------|
  | `FLAIL_BC_DIRICHLET_GHOST` (#87, the default) | value at the ghost centre | `-w phi_C` | `w` | `g` |
  | `FLAIL_BC_DIRICHLET` | value at the face | `-2 w phi_C` | `2 w` | `2 g - phi_C` |
  | `FLAIL_BC_NEUMANN` | outward derivative at the face | none | `w h` | `phi_C + h g` |

  Without a Dirichlet face `A` is singular as on a periodic grid: the mean of `factor V f - B g` and of the solution
  are removed. The V-cycle's levels take the same kinds, homogeneous, and the coarsest level the FFT solver with them
  (E1a): on a uniform grid one V-cycle is the exact solve.
- **Isolated boundaries (#87, P2), the James step:** `phi0 = A^-1 factor V f` (zero ghosts), the screening charge
  `s = B^T phi0` on the ghost cells, their free-space values `g = -G s` (`G` the lattice Green's function over the
  coarsest spacing, from ghost cell to ghost cell, by the padded FFT on the coarsest box widened by the ghost layer,
  gathered on rank 0), then `phi = A^-1 (factor V f - B g)`: the map `(A^-1 + A^-1 B G B^T A^-1) factor V` is
  symmetric in the V-weighted form, the exact discrete free-space solution on a uniform grid, second order on AMR.
  Two solves per call, each warm-started from its own last solution (`iterations0`, `iterations`).
- **Scope (fatal otherwise):** without `bc` every active direction periodic or none; with it, the grid's periodicity
  and `bc` agree, the two faces of an axis both ghost-Dirichlet or
  neither (the FFT solver's scope); isolated: 3-D, no periodic direction, no `bc`, cubic cells, every leaf on the
  domain boundary at the coarsest level; 2:1 balance across
  faces; fewer than `2^19` cells per direction at the finest level. Host code: FLUME's FNL backend runs it on the host
  with a cuFFT coarsest level (#81 P5b).
- **Build:** as `adam_poisson_fft_object` (`-D_FFTW`, not re-exported by `adam_common_library`; mode
  `test-poisson-mg-gnu`).
- **Verification:** SV-9, SV-12 (Dirichlet), SV-13 (isolated) and the boundary kinds (E1b, `--bc`: the FFT solver's
  result on uniform grids, a manufactured potential on AMR), `src/tests/poisson/test_poisson_mg.F90`, driven by
  `src/tests/flume/verification/selfgravity/poisson.sh` (legs `sv9`, `sv13`, `bcmg`).

### `adam_lattice_green_library` — lattice Green's function of the 7-point Laplacian (generated)

`lattice_green(n)` returns `G(n)` on the unit cubic lattice, `L_1 G = delta_(n,0)`, `G -> -1/(4 pi |n|)`: the kernel of
the isolated FFT solver. **Generated, do not edit:** `scripts/lattice_green.py` writes it and validates it
(`--check` compares the committed file byte for byte). For `max |n_d| <= 32` it is a table (the wedge
`n1 >= n2 >= n3 >= 0`, 6545 values): the Dirichlet problem on the cube `[-32, 32]^3`, solved exactly by a type-I sine
transform, with the asymptotic expansion on the faces. Beyond, the expansion itself:
`G0 + G1 + G2 + G3`, `G0 = -1/(4 pi r)`, `G1 = -5 (n1^4 + n2^4 + n3^4 - 3 r^4/5)/(32 pi r^7)`, `G2 = P8/r^13`,
`G3 = P12/r^19`, each order of the Taylor series of `L_1` solved symbolically; the two coefficients the series
leaves free (one cubic harmonic of degree 4 in `G2`, one of degree 6 in `G3`) are fitted to the table. Checked by the
generator: the table's residual `2.2e-16`; the residual of the composite kernel on `[-64, 64]^3`, across the faces
where table and expansion meet, `1.0e-13`; `G(0)` against Watson's closed form for the simple cubic lattice,
`2.0e-13` relative.

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
| `fill_seam_from_peer_forest(peer, p_idx)`    | Receive-side roundtrip: copy peer's interior into self's seam ghosts. Source/target buffer chosen by `self%stage_active`/`peer%stage_active` (`>0` ⇒ `rk%q_rk(:,...,stage)`, else committed `q`). |
| `after_topology_build_forest()`              | Backend hook fired once after the seam maps are built — FNL uses it to refresh the device-resident seam maps (`field_fnl%maps%copy_cpu_gpu`). |
| `apply_reflux_to_stage_forest(stage, dt, flux_register)` | Apply Berger-Colella reflux at the realm's end-of-step (α.r1 gate). |
| `stages_per_step_forest()`                   | Report K (the realm's integrator stride; `rk%nrk` for SSP-RK). |
| `coupling_descriptor_forest(scheme_time, rk_scheme, nv)` | Report `(scheme_time, rk_scheme, nv)` for β admissibility (issue #18). Default returns `''/''/-1` — a realm can't join β unless it overrides. |
| `is_done_forest(done)`                       | Report whether the simulation has reached its termination predicate. |
| `post_step_forest(dt, t, it, realm)`         | Per-step diagnostics / IO block (savers, residual prints). |

Each `_forest` TBP is `pass(self)`-dispatched. The forest receives `type(realm_pointer), intent(inout) :: realm(:)` as a parameter and never owns realm state: the driver allocates the realms as a concrete app type (e.g. `type(prism_cpu_object), allocatable, target :: realm(:)`) and points one `realm_pointer` at each, so element `is` is `realm(is)%ptr`. The forest never takes a polymorphic array `class(realm_object) :: realm(:)`: nvfortran (25.11, 26.9) can address its elements with the declared type's stride and dispatch on the wrong realm ([#96](https://github.com/szaghi/adam/issues/96)).

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

### Intra-realm AMR seams (a second, independent seam family)

The inter-realm machinery above couples two *realms*. A separate family — the **intra-realm AMR seam** — handles a 2:1 resolution jump *inside a single realm's* Morton octree (a level-`ℓ` block abutting a level-`ℓ+1` block, produced by AMR markers). It is not manifest-declared and uses different code:

- **Registration**: `forest_object%register_intra_realm_amr_seams(realm)` (private) does a two-pass registration of the coarse-fine faces into the flux register with `SEAM_KIND_INTRA_REALM_AMR` (`adam_flux_register_object.F90`), and allocates/fills `maps%amr_seam_quadrant`. That array's allocation is also the **gate** for the seam `div(B)` guard-rail (below).
- **Ghost fill**: the fine block's coarse-overlapping ghosts are filled through `field_object%update_ghost_local`'s **flag-4** map-row path → `interp_seam_ghost`. The interpolation *regime* is `maps%seam_ghost_fill`, one of `SEAM_FILL_INJECTION` (0, copy), `SEAM_FILL_COMPATIBLE` (1, q=2 restriction-compatible), `SEAM_FILL_TRICUBIC` (2, q=4) — declared in `adam_seam_interpolation_library.F90`, selected by `[amr] seam_ghost_fill`. **Effective default is `tricubic`**: `adam_object%initialize` overwrites the `maps` struct default (`SEAM_FILL_INJECTION`) with `SEAM_FILL_TRICUBIC` when the key is absent. The FNL device path carries the regime via `update_ghost_local_gpu(..., seam_ghost_fill=..., ...)`.
- **Reflux**: Berger-Colella reflux at the 2:1 face, dispatched by `apply_reflux_corrections` at each realm's own end-of-step (α.r1 gate), identical machinery to the inter-realm register.

The two families are structurally disjoint — a single-realm run never touches `fill_seam_from_peer_forest`; a same-resolution forest never touches `interp_seam_ghost`. They can coexist (a forest of internally-AMR-refined realms).

### Seam div(B) guard-rail (`[IO]`, issue #29)

The 2:1 AMR seam breaks the mimetic `div_h(curl_h) ≡ 0` identity (coarse Δx vs fine Δx/2 stencils no longer commute) and injects an O(h^p) `div(B)` source that is refinement-convergent but unbounded in `t` at fixed `h` — a constraint violation, not an energy instability. No local collocated-FD fix exists (five classes ruled out); the resolution is **accept-truncation**. To make the otherwise-silent growth visible, `save_divergence_history` (`adam_prism_common_object.F90`) carries an opt-in monitor: `[IO] seam_divB_tol` (default `-1.0`, off) and `seam_divB_error` (default `.false.`, warn-only). It arms only when `seam_divB_tol > 0` **and** `allocated(maps%amr_seam_quadrant)` (an AMR seam is present); on exceedance it `error_stop`s (`seam_divB_error = .true.`) or warns. `seam_divB_tol`/`seam_divB_error` live on the base `adam_io_object`.

### Load-bearing invariants

- **Phase 2 → Phase 3 ordering (β)**: Phase 2 must complete on ALL realms before Phase 3 starts on ANY realm; otherwise the read-after-overwrite race between `fill_seam_from_peer_forest` writes and `compute_residuals` reads returns. The serial inner loops within a rank give this for free under the Phase-A replicated-forest layout.
- **Inter-realm seam diagnostic refill (issue #31)**: any div(B) diagnostic that runs `self%update_ghost` before reading the divergence stencil must ALSO refill the inter-realm seam (`fill_seam_from_peer_forest` per peer), because `update_ghost` fills intra-realm ghosts + physical BCs but NOT the inter-realm seam (the forest owns it). Omitting it reports a spurious seam-skin div(B) on a div-free field. `present(realm)`-guarded so single-realm runs are untouched. A 1:1 same-resolution inter-realm mirror seam must use `stage_coincident` (β) to stay div-free — α leaves the seam ghosts unfilled during RK substages.
- **`stage_active` discipline**: `begin_stage_forest(k)` MUST set `self%stage_active = k`; `close_step_forest(dt)` MUST clear it back to 0. The `fill_seam_from_peer_forest` buffer-selection logic reads `self%stage_active` and `peer%stage_active` to pick between `q` (committed) and `rk%q_rk(:,...,stage_active)`. Mid-step writes to `q` outside this discipline corrupt peer reads.
- **α.r1 reflux gate**: `apply_reflux_to_stage_forest` and `accumulate_seam_fluxes_fv` MUST gate on `stage == rk%nrk` (the realm's own final substage). `flux_register`'s third axis is collapsed to 1 under α.r1; per-stage RK-weighted reflux (Wang 2018) is deferred. Both α and β share this gate — β does not undo α.r1.
- **N=1 fast path**: when `size(realm) == 1`, `evolve_one_step` routes through `advance_one_step_forest` directly with zero exposure to the per-stage TBPs. App backends MAY leave the per-stage TBPs at the parent's error_stop default for N=1-only use.
- **No singletons embedded on realms**: the program-scope singletons (`mpih`, `grid`, `field`, ...) are aliased per-step from each active realm's components; they are not value members on the realm type. See "Program-scope singletons" above.

### Tracked through

| PRD | Scope | Status |
|---|---|---|
| [#10](https://github.com/szaghi/adam/issues/10)  | Phase D inception — forest orchestrator design          | shipped |
| [#13](https://github.com/szaghi/adam/issues/13)  | Coarse-fine interface machinery (A/B/C)                 | A + C shipped; B characterized (→ #29) |
| [#16](https://github.com/szaghi/adam/issues/16)  | α end-of-step barrier + asymmetric K                    | shipped |
| [#18](https://github.com/szaghi/adam/issues/18)  | β stage-coincident recovery + cross-config oracle       | shipped |
| [#21](https://github.com/szaghi/adam/issues/21)  | tricubic coarse→fine seam ghost fill (`seam_ghost_fill`) | shipped |
| [#22](https://github.com/szaghi/adam/issues/22)/[#23](https://github.com/szaghi/adam/issues/23) | FNL seam fill + FV residual/reflux device twins | shipped |
| [#28](https://github.com/szaghi/adam/issues/28)  | np>1 reflux register ownership-gating                   | shipped |
| [#29](https://github.com/szaghi/adam/issues/29)  | 2:1 seam div(B) — accept-truncation + guard-rail        | shipped (resolution) |
| [#31](https://github.com/szaghi/adam/issues/31)  | inter-realm 1:1 seam div-free (β + diagnostic refill); FNL fWLayer HtoD | shipped |

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
