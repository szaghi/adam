# FLAIL elliptic solver

FLAIL is the one entry point for the elliptic solves of every application: the realm component `flail`
(`self%flail`, type `flail_object`). It wraps two solvers: a direct FFT solver on uniform grids and a
multigrid-preconditioned conjugate-gradient solver on AMR grids. FLUME uses it for self-gravity, PRISM for its
time-zero field solves.

| Module | Content |
|--------|---------|
| `adam_flail_object` | The facade: solver choice, forwarding, boundary kinds |
| `adam_flail_parameters` | The boundary kinds `FLAIL_BC_*` (a module of their own: the solvers use them, the facade uses the solvers) |
| `adam_poisson_fft_object` | [FFT solver](#fft-solver) on uniform grids |
| `adam_poisson_mg_object` | [Multigrid solver](#multigrid-solver) on 2:1-balanced AMR grids |
| `adam_lattice_green_library` | [Lattice Green's function](#lattice-green-s-function) of the 7-point Laplacian (generated) |
| `adam_fnl_poisson_fft_object`, `adam_fnl_poisson_mg_object` | Device twins ([FNL library](./fnl#poisson-solvers)) |

## The facade

| Procedure | Role |
|-----------|------|
| `setup(grid, field, isolated, backend[, option_name][, bc][, dg_order])` | Choose the solver for the current grid; call it on the initial grid and after every regrid |
| `solve(field, q, iv, factor, phi, energy[, ghost])` | Solve for variable `iv` of `q`; `energy = sum f phi dV` |
| `solve_packed(field, buf, factor, phi, energy[, ghost])` | The same from the interior source already packed |
| `iterations()` | Last FCG iteration count (0 for the FFT solver) |
| `is_multigrid()`, `destroy()` | The solver in use; cleanup |

| `setup` argument | Values |
|------------------|--------|
| `backend` | `FLAIL_BACKEND_AUTO` (FFT on a uniform grid, multigrid otherwise), `FLAIL_BACKEND_FFT`, `FLAIL_BACKEND_MULTIGRID` |
| `isolated` | `phi -> 0` at infinity; fatal together with `bc` |
| `bc(2,3)` | The kind of every face ([boundary kinds](#boundary-kinds)); without it, periodic or isolated |
| `dg_order` | 0 (compact Laplacian, default) or 2, 4, 6, 8: PRISM's wide operator D(G), FFT path only (fatal on AMR) |
| `ghost` (in `solve`) | The boundary value of each boundary face, multigrid only (fatal on the FFT path, which solves the homogeneous problem) |

A backend picks its FFT engine by allocating the `fft` component before `setup` (FLUME's FNL backend allocates the
cuFFT one), and its multigrid through `flail%mg_mold`. A vector whose components take different kinds (PRISM's vector
potential at a PEC face: Neumann for the normal component, Dirichlet for the others) uses one FLAIL per set of kinds;
the instances coexist. The mapping from an application's boundaries to the kinds belongs to the application.

### Input keys

FLAIL has no input section. The applications expose it as follows.

| Application | Key | Values |
|-------------|-----|--------|
| FLUME | `[gravity] poisson_solver` | `auto` (default), `fft` (an AMR grid is fatal), `multigrid` |
| FLUME | `[gravity] poisson_bc` | `periodic` (default), `isolated` |
| PRISM | none | the t0 solves map `[bc_*] type` to the kinds and set `dg_order = fdv_order` |

## Boundary kinds

| Kind | Meaning | FFT | Multigrid |
|------|---------|-----|-----------|
| `FLAIL_BC_PERIODIC` | Periodic axis | yes | yes |
| `FLAIL_BC_DIRICHLET` | Zero (or `g`) at the face | yes (odd extension) | yes |
| `FLAIL_BC_NEUMANN` | Zero (or `g`) normal derivative at the face | yes (even extension) | yes |
| `FLAIL_BC_DIRICHLET_GHOST` | Zero (or `g`) at the first ghost centre | yes (type-I sine transform) | yes (default on bounded axes) |
| `FLAIL_BC_ZERO_GHOSTS` | Every ghost layer zero (PRISM's elliptic Dirichlet) | yes (CG with D(G), both faces of an axis) | as ghost-Dirichlet |
| `FLAIL_BC_FREE_GHOSTS` | Ghosts hold the free-space potential of the source (PRISM's open face) | wide D(G) only | refused |

Rules on one axis: both faces periodic or neither; a ghost-Dirichlet or zero-ghost face against another kind is fatal.

## FFT solver

`adam_poisson_fft_object` solves `L_h phi = factor (f - mean f)` (periodic) or `L_h phi = factor f` for one variable
of `q(nv,i,j,k,b)`. `L_h` is the second-order 7-point Laplacian (5-point in 2-D; a null direction has no derivative)
or, with `dg_order`, the wide D(G).

| Problem | Method | Residual |
|---------|--------|----------|
| Periodic | Each Fourier mode divided by the eigenvalue of `L_h`, `lambda_k = sum_d (2 cos(2 pi k_d/N_d) - 2)/dx_d^2` | round-off |
| Dirichlet/Neumann faces | Odd extension about a Dirichlet face, even about a Neumann one, over `2n` cells (faces agree) or `4n` (they differ); periodic transform of the extension. Without a Dirichlet face the source mean is dropped | round-off; wall placement second order |
| Isolated (Hockney–Eastwood) | Source zero-padded to `m >= 2 (n + ngc)` cells per direction (the smallest even size with no prime factor above 13), convolved with the lattice Green's function, whose spectrum is computed once | `~1e-13` per offset, summed over the source (`2e-11` of `factor f` at `64^3`) |
| Wide D(G), order p | Eigenvalue `-sum_d s_p(theta_d)^2 / dx_d^2`, `s_p(theta) = sum_j 2 a_j sin(j theta)`; the checkerboard modes (`theta_d` in {0, pi}) are in the kernel and dropped; reflections commute with the stencil, so Dirichlet/Neumann kinds hold and `p` ghost layers are filled | round-off |
| Zero-ghost faces with D(G) | The truncated stencil, solved on rank 0 by CG (tolerance `1e-12`, max-norm), preconditioned by the wide p = 2 operator inverted exactly (`2^z` FFT sub-solves on the parity sub-lattices) | `1e-12` |
| Free-space ghosts with D(G) | Ghosts hold `g(x) = -sum_y factor f(y) dV / (4 pi |x - y|)`, summed on rank 0 over the gathered source; the zero-ghost solve with `D(G) u_g` on the right-hand side. Edges and corners of the ghost region stay zero | `1e-12` |

Preconditioner iteration counts for D(G) with zero-ghost faces: 1, 17, 24, 30 at p = 2, 4, 6, 8, independent of the
grid. A compact preconditioner needs O(N) iterations instead (86 → 677 at p = 2, N = 32 → 256), because it does not
see the checkerboard branch (`src/tests/poisson/wide_dg_proto.py`, `wide_dg_zero_ghost.py`).

Every solve is a symmetric linear map of the source, which exact energy and momentum conservation of self-gravity
needs (FLUME [models](/applications/flume/models#self-gravity)).

| Procedure | Role |
|-----------|------|
| `initialize(grid, field, isolated, bc, dg_order)` | Check the scope, size the solver, make the plans (and the kernel spectrum when isolated) |
| `solve(field, q, iv, factor, phi, energy)` | Transform and return `phi(1-ngc:n+ngc, …, b)` with every ghost layer from the global solution; optional `energy = sum f phi` |
| `solve_packed(field, buf, factor, phi, energy)` | The same from the packed interior source (block-major, `i` fastest) |
| `initialize_blocks(nblk, nijk, ngc, dxyz, level, null_xyz, isolated, dirichlet, bc, redundant)` | Size the solver for an explicit uniform block lattice, without a grid object (the multigrid's root). `redundant`: every rank holds the arrays and solves, no scatter. `dirichlet`: ghost-Dirichlet on every axis |
| `solve_blocks(nbl, coords, buf, factor, gw, out, energy)` | Solve on an explicit list of this rank's blocks, returning `gw` ghost layers |
| `destroy()` | Free arrays and plans |

**Data path.** The transform is distributed over z-slabs. Each rank owns z-planes of the transform array; block planes
go to their owners, planes are transformed in 2-D, transposed to z-lines (`MPI_ALLTOALLV`) and transformed in 1-D,
divided by the eigenvalues (or multiplied by the kernel spectrum), and sent back. Each block fetches its cells with
their ghosts from the owners, so no ghost exchange follows. The stages are backend hooks (`slab_update`,
`slab_planes`, `slab_transpose`, `slab_lines`, `slab_spectral`, `slab_keep_kernel`). The potential is bitwise the same
on every rank count (`poisson.sh` leg `slab`: CPU np 1–3, FNL np 1, 2, 4). At most `m(3)` ranks hold planes. The CG
solves, the multigrid's redundant root and the CG sub-problems gather to rank 0 instead. `poisson_fft_object%wall(3)`
accumulates the gather, solve and scatter times of `solve_blocks`.

**Scope (fatal otherwise).** Every block at one level, covering the domain. Without `bc`, every active direction
periodic. Isolated: 3-D, no periodic direction, cubic cells.

**Build.** `-D_FFTW -lfftw3` (`FFTW_ESTIMATE` plans: deterministic). Without an FFT library `initialize` stops with
"built without an FFT library". Every PRISM mode and FLUME's `*-fft` modes link one (FFTW on the host, cuFFT on NVIDIA
FNL).

**Verification.** SV-2 (periodic) and SV-6 (isolated): `src/tests/poisson/test_poisson_fft.F90`, driven by
`src/tests/flume/verification/selfgravity/poisson.sh`. Face kinds: `src/tests/poisson/test_flail_bc.F90` (per-axis
discrete eigenfunctions of every kind mix at round-off on 1–3 ranks; Dirichlet/Dirichlet order against the continuum
1.99; refusals).

## Multigrid solver

`adam_poisson_mg_object` solves `A phi = factor V f - B g` on the leaves of a 2:1-balanced octree (quadtree with
`null_z`), periodic or with a kind per face; the mean of the right-hand side is removed when `A` is singular.

`A` is the composite second-order Laplacian in finite-volume form: the 7-point stencil between leaves of one level
and, across a 2:1 face, the mimetic mean flux. Every fine cell on a coarse cell's face gets the flux
`(2/3) w (mean of the face's fine cells - phi_C)`, `w` the face area over the spacing; the coarse cell gets their sum.

| Property | Value |
|----------|-------|
| Consistency | Exact for linear `phi` face by face: `phi` second order, the 2:1 face gradient first order |
| Symmetry | `V K` symmetric for `phi = K f`: the exact energy exchange of self-gravity |
| Alternatives rejected | Flux matching (Tomida & Stone): consistent, not symmetric. Its symmetrisation: symmetric, O(1/h) truncation (`selfgravity/seam_operators.py`) |

| Procedure | Role |
|-----------|------|
| `initialize(grid, field, tolerance, max_iterations, coarse, isolated, bc, root_cells)` | Gather the leaf list, derive the level hierarchy, assemble the operators, size the root FFT solver (of the dynamic type of `coarse`, FFTW otherwise) and, isolated, the ghost-to-ghost Green's function |
| `solve(field, q, iv, factor, phi, energy, ghost)` | FCG on `A` preconditioned by one V-cycle, warm-started when the grid has not changed; `phi` with every ghost filled |
| `solve_packed(field, buf, factor, phi, energy, ghost)` | The same from the packed source |
| `symmetry_defect(x, y)`, `apply_a(x, y)` | Verification: `|x.Ay - y.Ax| / (|x| |Ay|)` and `y = A x` |
| `apply_b(g, y)`, `apply_bt(x, s)` | `y = B g` and `s = B^T x` over the boundary faces (`nbnd`, `bnd_row`, `bnd_w`, `bnd_kind`, `bnd_dir`, `bnd_x`) |
| `destroy()` | Free everything |

| Component | Choice |
|-----------|--------|
| Krylov solver | Flexible CG (Notay's FCG(1)), relative max-norm residual: `tolerance` (default `1e-10`) or, when the source mean dominates, the round-off floor `100 eps max|factor f| / max|factor (f - mean f)|`; `tolerance_used` keeps the one applied; reaching `max_iterations` (200) prints a warning |
| V-cycle | Leaves plus one shadow block per refined node; red-black Gauss-Seidel, 2 sweeps down and 2 up; restriction by sum; piecewise-constant prolongation |
| Below the coarsest tree level | Levels coarsen inside the blocks (same blocks, same ranks, half the cells per level) while the cell counts are even and above `root_cells` |
| Root | `solve_blocks` of the FFT solver; redundant on every rank (`MPI_ALLGATHERV`) up to 32^3 cells, distributed above. Default `root_cells`: a host transform of 128^3 cells (periodic root 128^3, ghost-Dirichlet root 64^3); 32^3 on cuFFT |
| Ghosts | First layer from the solution (across a 2:1 face, the value matching the operator's face gradient); deeper layers copy the first |
| Distribution | CSR rows over global cell keys, one halo plan per operator (`MPI_ALLTOALLV`); results equal to round-off across rank counts |
| Backend hooks | `apply`, `smooth`, `vec_fill`, `vec_copy`, `vec_lincomb`, `vec_axpy`, `vec_dot`, `vec_relmax`, `vec_permute`, `vec_update`, `work_place`, on vectors allocated once by `initialize`. The right-hand side, the root, the James step and the ghost fill stay on the host |

The root design follows Athena++ (Tomida & Stone 2023, §3). Coarsening inside the blocks costs FCG iterations
(`test_poisson_mg --sub`: 1 → 14 on a uniform box of 1:2:3 cells, from point smoothing on anisotropic cells; +1 to +13
on AMR) and removes the rank-0 gather of the base grid in every V-cycle.

**Boundary faces.** An active axis is periodic or bounded. Each bounded face has a kind; its row is homogeneous
(symmetric) and its value `g` goes to the right-hand side (`w` the face area over the spacing, `h` the spacing):

| Kind | `g` | Row | `B` | Ghost |
|------|-----|-----|-----|-------|
| `FLAIL_BC_DIRICHLET_GHOST` (default) | value at the ghost centre | `-w phi_C` | `w` | `g` |
| `FLAIL_BC_DIRICHLET` | value at the face | `-2 w phi_C` | `2 w` | `2 g - phi_C` |
| `FLAIL_BC_NEUMANN` | outward derivative at the face | none | `w h` | `phi_C + h g` |

Without a Dirichlet face `A` is singular: the means of `factor V f - B g` and of the solution are removed. The V-cycle
levels take the same kinds, homogeneous; on a uniform grid that is its own root, one V-cycle is the exact solve.

**Isolated boundaries (James step).** `phi0 = A^-1 factor V f` (zero ghosts), the screening charge `s = B^T phi0` on
the ghosts, their free-space values `g = -G s` (`G` the lattice Green's function on the coarsest spacing, by the
padded FFT on the coarsest box widened by the ghost layer, on rank 0), then `phi = A^-1 (factor V f - B g)`. The map is
symmetric in the V-weighted form: the exact discrete free-space solution on a uniform grid, second order on AMR. Two
solves per call, each warm-started (`iterations0`, `iterations`).

**Scope (fatal otherwise).** Without `bc`, every active direction periodic or none. With it, `bc` agrees with the
grid's periodicity. Isolated: 3-D, no periodic direction, no `bc`, cubic cells, every leaf on the domain boundary at
the coarsest level. 2:1 balance across faces; fewer than `2^19` cells per direction at the finest level.

**Build and verification.** As the FFT solver (`-D_FFTW`; mode `test-poisson-mg-gnu`). SV-9, SV-12 (Dirichlet), SV-13
(isolated) and the boundary kinds (`--bc`: the FFT solution on uniform grids, a manufactured potential on AMR):
`src/tests/poisson/test_poisson_mg.F90`, driven by `selfgravity/poisson.sh` (legs `sv9`, `sv13`, `bcmg`).

## Lattice Green's function

`lattice_green(n)` returns `G(n)` on the unit cubic lattice, `L_1 G = delta_(n,0)`, `G -> -1/(4 pi |n|)`: the kernel
of the isolated FFT solver. **Generated, do not edit:** `scripts/lattice_green.py` writes it and `--check` compares the
committed file byte for byte.

| Range | Content |
|-------|---------|
| `max |n_d| <= 32` | Table of the wedge `n1 >= n2 >= n3 >= 0` (6545 values): the Dirichlet problem on `[-32, 32]^3`, solved by a type-I sine transform with the asymptotic expansion on the faces |
| Beyond | `G0 + G1 + G2 + G3`, `G0 = -1/(4 pi r)`, `G1 = -5 (n1^4 + n2^4 + n3^4 - 3 r^4/5)/(32 pi r^7)`, `G2 = P8/r^13`, `G3 = P12/r^19`; the two free coefficients fitted to the table |

| Check (by the generator) | Value |
|--------------------------|-------|
| Table residual | `2.2e-16` |
| Composite kernel residual on `[-64, 64]^3` | `1.0e-13` |
| `G(0)` against Watson's closed form | `2.0e-13` relative |

## PRISM time-zero solves

`impose_pic_fields_time_zero` (`adam_prism_common_object.F90`) solves `D(G) phi = -rho/eps` (PIC; then
`D <- D - eps G phi`) and `D(G) A_d = -mu J_d` per component (then `B <- B + curl A`), one `setup` per component,
`dg_order = fdv_order`, on the FFT path.

| PRISM face | FLAIL kind |
|------------|------------|
| `periodic` | `FLAIL_BC_PERIODIC` |
| `Dirichlet` | `FLAIL_BC_ZERO_GHOSTS` |
| `PEC` | `FLAIL_BC_DIRICHLET` (potential, tangential A), `FLAIL_BC_NEUMANN` (normal A) |
| `Neumann`, `Silver_Muller` (exact/open) | `FLAIL_BC_FREE_GHOSTS` |
| null axis | periodic |

| Limit | Detail |
|-------|--------|
| Refined grids | Refused (D(G) is not defined across 2:1 faces) |
| Memory of the CG solves | Zero-ghost and free-space solves gather the grid on rank 0 (and its GPU on FNL), about 24 B per transform cell; PEC solves are distributed (8 N cells with reflected faces) |
| Checkerboard | Centred D(G) annihilates `(-1)^(i+j+k)`; with PEC faces that component of `rho` stays in `div D - rho` |
| Free-space ghost sum | Direct, on rank 0, O(ghost-slab cells × cells) per component |

Regression cases `pic-pec-t0`, `rmf-pec-t0`, `pic-open-t0` and `rmf-open-t0` run these solves.

## Performance

`exe/bench_poisson` (modes `bench-poisson-gnu`, `bench-poisson-fnl-nvf`) times FLAIL on the cube `[-1,1]^3` with a
Gaussian source, uniform or with the central box refined. `src/tests/poisson/bench.sh [cpu|fnl] [--weak] [--strong]
[--np-max N]` runs the ladders and prints one row per run. Solve time: minimum over the repetitions, maximum over the
ranks. Machine: medooza (CPU nodes; A30 GPUs for FNL).

| Solver | Grid | Backend | Problem | 1 | 2 | 4 | 16 ranks / GPUs |
|--------|------|---------|---------|---|---|---|-----------------|
| FFT | 2.1 M uniform cells | CPU | periodic | 0.115 s | | | 0.026 s |
| FFT | 2.1 M uniform cells | CPU | isolated | 0.729 s | | | 0.108 s |
| FFT | 2.1 M uniform cells | FNL | periodic | 0.067 s | | 0.030 s | |
| FFT | 2.1 M uniform cells | FNL | isolated | 0.138 s | | 0.112 s | |
| Multigrid | 3.9 M AMR cells | CPU, 128^3 root | periodic | 8.6 s | | 2.7 s | 1.6 s |
| Multigrid | 3.9 M AMR cells | CPU, 64^3 root | isolated | 26.8 s | | 8.5 s | 5.2 s |
| Multigrid | 3.9 M AMR cells | FNL, 32^3 cuFFT root | periodic | 0.33 s | 0.19 s | 0.12 s | |
| Multigrid | 3.9 M AMR cells | FNL, 32^3 cuFFT root | isolated | 0.85 s | 0.52 s | 0.43 s | |

Weak scaling of the FFT solver, CPU, 8 ranks, 16.7 M uniform cells: periodic 0.24 s, isolated 1.5 s per solve.
