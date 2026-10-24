# The multigrid Poisson solver on AMR grids

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-9, SV-12, SV-13 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5a; [#87](https://github.com/szaghi/adam/issues/87), P1, P2 | CPU only | yes | — |

## Case and scope

On an AMR grid the exact $E + W$ exchange needs the composite Laplacian $A$ (finite-volume form, rows weighted by the
cell volume $V$) to be **symmetric**: $\Phi = 4\pi G\,K\rho$ with $VK$ symmetric, the AMR form of $K = K^T$ now that
the volumes differ. Accuracy needs it **consistent**: a bounded truncation error at the 2:1 faces and a converging
face gradient there, the face gravity the momentum and energy sources use. The coupling across a 2:1 face decides
both, and the usual choices give one or the other.

The case has two parts:
- **the operator study** (`seam_operators.py`): three couplings of a coarse cell C and the fine cells $f_k$ on its
  face, measured in a 2-D model, to choose the one the solver uses;
- **the solver leg** (`test_poisson_mg`): the multigrid solver on that operator, on uniform grids against the FFT
  solver (SV-9a) and on box-refined grids against the exact potential (SV-9b);
- **the Dirichlet legs** (`test_poisson_mg`, [#87](https://github.com/szaghi/adam/issues/87) P1): the same solver with
  Dirichlet domain faces, the ground the isolated solve on AMR is built on: on a uniform grid the coarsest-level solve
  (the FFT solver's Dirichlet mode) is the whole solve (SV-12a); on box-refined grids against the exact potential
  (SV-12b);
- **the isolated legs** (`test_poisson_mg --isolated`, [#87](https://github.com/szaghi/adam/issues/87) P2): free-space
  boundaries on AMR by a symmetric James step, SV-6's stars: on uniform grids against the FFT solver's isolated mode
  (SV-13a), on grids with the star refined against the exact potentials (SV-13b), the symmetry of the solution map
  measured in situ;
- **the boundary-kind legs** (`test_poisson_mg --bc`, [#84](https://github.com/szaghi/adam/issues/84) E1b): a kind per
  face, periodic, Dirichlet or Neumann at the face, or ghost-Dirichlet, mixed across and within axes: on uniform grids
  against the FFT solver with the same kinds, on grids with a refined corner (leaves of two levels on the boundary)
  against a manufactured potential whose boundary values of each kind are the ghost values, the domain-face ghosts
  included.

What the case catches:
- a non-symmetric assembled operator (the symmetry defect of the assembled rows, which the assembly does not
  symmetrise);
- an inconsistent seam coupling: a truncation error that grows as $1/h$ under refinement;
- a multigrid that is not the FFT solve on a uniform grid, a ghost the fill leaves empty, a stalled FCG;
- a boundary kind whose row, right-hand side or ghost is wrong (the order against the manufactured potential, the
  truncation, the agreement with the FFT solver);
- a mixed grid accepted without `bc`, boundary kinds that disagree with the grid's periodicity, or that the
  coarsest-level FFT solver cannot take, accepted;
- an isolated solve that is not the free-space solution on a uniform grid, or not a symmetric map ($VK \ne (VK)^T$:
  $E + W$ would not be conserved);
- an isolated solve accepted with a refined leaf on the domain boundary.

## Models and numerics

**The three couplings** (the operator study):
- **B**, Tomida and Stone's flux matching: the fine ghost interpolated from C and its tangential neighbours, the coarse
  flux the sum of the fine ones. Consistent, not symmetric.
- **S1** $= (B + B^T)/2$, the symmetrisation the P5-0 design proposed. The transposed tangential terms land in the
  coarse rows with no matching Taylor terms.
- **M**, the **mimetic mean flux** (`adam_poisson_mg_object`): every fine face carries
  $\tfrac23 w\,(\bar\Phi_f - \Phi_C)$, $\bar\Phi_f$ the mean of the face's fine cells and $w$ the fine face area over
  the fine spacing; C gets their sum. It is the hybrid mimetic finite-difference flux of the coarse cell, seen as a
  polygon with hanging nodes, with its face unknowns eliminated: consistency makes the coarse cell's sub-face block
  rank one (it sees only the mean of its face unknowns), and the elimination gives $\bar\Phi_f$. Exact for linear
  $\Phi$ face by face, and symmetric by construction.

| Item | Choice | Code |
|---|---|---|
| Operator | the composite Laplacian of the leaves, $V$-weighted, with M across 2:1 faces; assembled from the geometry as distributed sparse rows over global cell keys, each with its halo-exchange plan | `initialize` (`assemble_a`), `a_row`, `face_cells`, `face_weight`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| 2:1 ghosts | the first ghost layer across a 2:1 face is the value whose one-sided difference equals M's face gradient | `initialize` (`assemble_ghost_fill`), [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| Iteration | flexible conjugate gradients (Notay's FCG(1)) on $A$, stopping at a relative max-norm residual of $10^{-10}$ | `solve`, `solve_packed`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| Preconditioner | one V-cycle. Levels: each holds its leaves plus a **shadow** block for every refined node (ADAM keeps leaves only). Smoothing: red-black Gauss–Seidel, 2 sweeps down and 2 up. Transfers: restriction by sum, piecewise-constant prolongation | `vcycle`, `smooth`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| Coarsest level and root ([#84](https://github.com/szaghi/adam/issues/84) E5c) | the minimum leaf level, which covers the domain; below it the levels coarsen inside the blocks (the same blocks and ranks, half the cells per block per level along the active directions) while those are even and the level has more than `root_cells` cells (default the root solver's `root_cells_max`: $32^3$ on the host, $128^3$ on cuFFT, which solves faster than the host smooths); the last, the root, is solved exactly by the FFT solver, on every rank when it has at most $32^3$ cells (`redundant`: the source gathered by `MPI_ALLGATHERV`, no scatter), on rank 0 above that (every rank transforming a $128^3$ root contends for the node's memory). A coarsest level of at most that many cells is its own root | `initialize` (`lbot`, `assemble_levels`), `vcycle`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90); `solve_everywhere`, [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |
| Boundary faces ([#87](https://github.com/szaghi/adam/issues/87) P1; per face, [#84](https://github.com/szaghi/adam/issues/84) E1b) | an active axis is periodic or bounded, each bounded face of a kind: ghost-Dirichlet (the #87 default; $g$ at the ghost centre: $-w\Phi_C$ in the row, $wg$ to the right-hand side), Dirichlet at the face ($-2w\Phi_C$, $2wg$) or Neumann at the face ($g$ the outward derivative: no row term, $whg$); without a Dirichlet face the operator is singular and the means are removed; the V-cycle's levels take the same kinds and the coarsest level the FFT solver with them (odd and even extensions, E1a) | `a_row`, `assemble_ghost_fill`, `solve(..., ghost)`, `apply_b`, `apply_bt`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| Isolated, the James step ([#87](https://github.com/szaghi/adam/issues/87) P2) | `initialize(..., isolated=.true.)`: the zero-ghost solve $\Phi_0 = A^{-1}4\pi G V\rho$; the screening charge $s = B^T\Phi_0$ on the ghost cells; their values $g = -Gs$, $G$ the inverse of the $V$-weighted coarsest operator on the infinite lattice (the lattice Green's function over $H$), by the zero-padded FFT on the coarsest box widened by the ghost layer; then $\Phi = A^{-1}(4\pi G V\rho - Bg)$. The map $K = A^{-1} + A^{-1}BGB^TA^{-1}$ is symmetric since $A$ and $G$ are. Each solve warm-starts from its own last solution | `solve_packed`, `ghost_green`, `solve_fcg`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90); `solve_blocks` (isolated), [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |
| Scope | without `bc` every active direction periodic or none (a mixed grid would get a silent zero on its bounded faces); with it, the grid's periodicity and `bc` agree, an axis's two faces both ghost-Dirichlet or neither; isolated: no periodic direction, no `bc`, 3-D, cubic cells, and every leaf on the domain boundary at the coarsest level (D-87-2; the ghost cells then form the coarse lattice $G$ is exact on); 2:1-balanced faces; anything else is fatal | `initialize` (boundaries), `assemble_ghost_fill` (the boundary layer), `classify` (2:1 balance), [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
| Verification hooks | the symmetry defect $\lvert x\cdot Ay - y\cdot Ax\rvert/(\lvert x\rvert\lvert Ay\rvert)$, $y = Ax$ on the leaves | `symmetry_defect`, `apply_a`, [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |

The operator on AMR grids is described in [models](../../models#self-gravity).

## Initial and boundary conditions

| Leg | Domain and source | Grid | Ranks |
|---|---|---|---|
| Operator study (`seam_operators.py`) | the doubly periodic unit square; $\Phi = \sin 2\pi x\cos 2\pi y + \tfrac12\cos(2\pi x + 4\pi y)$ | coarse ($H = 1/N$) on $x < 1/2$, fine ($h = H/2$) on $x > 1/2$: two seams ($x = 1/2$ and $x = 0 \equiv 1$); $N$ = 8, 16, 32, 48; dense solves | — (NumPy) |
| SV-9a | SV-2's: $[0,1]\times[0,2]\times[0,3]$, periodic, the two oblique modes, $G = 1$, the exact potential in closed form | uniform, octree and quadtree at level 2 | 1, 2, 3 |
| SV-9b | the same | the central box $[L/4, 3L/4]$ refined once over base levels 2 and 3: octree (120 and 960 leaves) and quadtree (28 and 112); 8 cells per block side, 3 ghosts | 1, 2, 3 |
| SV-12a | the same source and box without periodic boundaries, the exact potential as the ghost values (Dirichlet at the centres of the first ghost layer): $A\Phi = 4\pi G V(\rho - \rho_0) - Bg$ | uniform octree at levels 1 and 2 | 1, 2, 3 |
| SV-12b | the same | SV-9b's box-refined grids, octree and quadtree | 1, 2, 3 |
| SV-13a | SV-6's: the cube $[-1,1]^3$, isolated, the $n = 1$ polytrope $\rho = \sin(kr)/(kr)$, $k = \pi/a$, radius $a = 1/2$, $G = 1$, the exact potential | uniform octree at levels 1 and 2 ($16^3$, $32^3$) | 1, 2, 3 |
| SV-13b | the same polytrope, and the uniform sphere ($\rho = 1$, cell averages of $8^3$ samples) | base levels 2 and 3 with the star's box $[-1/2, 1/2]^3$ refined once (120 and 960 leaves; finest $64^3$ and $128^3$ equivalent) | 1, 2, 3 |
| E1b, uniform ([#84](https://github.com/szaghi/adam/issues/84)) | SV-2's box and source, homogeneous, the face kinds per axis (x, y, z): P/P D/N N/D; D/D N/N G/G; P/P N/N N/N (singular); 2-D N/D P/P | uniform octree at level 1, quadtree at level 2 | 1, 2, 3 |
| E1b, AMR | the same kinds; the manufactured $\Phi_m = \prod_d\cos(2\pi x_d/L_d + \theta_d)$, $\theta = (0.3, 0.7, 1.1)$, its Laplacian the source and its boundary values of each kind the ghost values | base levels 2 and 3 with the corner octant $[0, L/2]$ refined once (leaves of both levels on the low faces), octree and quadtree | 1, 2, 3 |
| E5c, levels inside the blocks ([#84](https://github.com/szaghi/adam/issues/84)) | SV-9a's and SV-9b's box and source (periodic, D/D N/N G/G, P/P D/N N/D) and SV-13b's polytrope (isolated) | base level 2: uniform octree, the box refined once (octree and quadtree); every level forced down to one cell per block (`root_cells = 1`: 3 levels for $8^3$ blocks) against none (`root_cells = huge`) | 1, 2, 3 |
| Refusals | `--mixed`: a grid periodic along $x$ only and no `bc` (without it, every active direction periodic or none); `--refuse-bc-mismatch`: the same grid, `bc` periodic along $y$ as well; `--refuse-bc-mixed-ghost`: the same grid, a ghost-Dirichlet face and a Dirichlet one along $y$; `--refuse-boundary-refined`: the isolated cube at level 1 with the octant $x, y, z > 0$ refined | — | 2 |

The FNL path comes with FLUME on AMR (P5b, P5c: [SV-10](./amr-jeans), [SV-11](./runtime-amr)).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/seam_operators.py` | the operator study (`[N ...]`, default 8 16 32 48) |
| `src/tests/poisson/test_poisson_mg.F90` | the unit test (mode `test-poisson-mg-gnu`), SV-9a, SV-9b, SV-12a, SV-12b; `--isolated` SV-13a, SV-13b; `--bc` the boundary kinds (E1b); `--mixed`, `--refuse-bc-mismatch`, `--refuse-bc-mixed-ghost` and `--refuse-boundary-refined` the refusals |
| `src/tests/flume/verification/selfgravity/poisson.sh` | legs `sv9`, `sv13`, `bcmg` and `sub`: the test on every rank count of `--np-list`, and the refusals; `POISSON_MG_EXE` selects the build |
| `src/tests/flume/verification/make_doc_runs.sh`, group `amr` | `sv9-poisson` and `sv13-poisson`: the test (and `--isolated`) with its log |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_amr`, `fig_selfgravity_amr_isolated` | draws the figures (`--only selfgravity-amr selfgravity-amr-isolated`), the first running `seam_operators.py` |

```bash
python3 src/tests/flume/verification/selfgravity/seam_operators.py      # the three couplings
cd src/tests/flume/verification/selfgravity && ./poisson.sh --legs "sv9 sv13 bcmg"   # the solver, np 1, 2, 3
```

## Results

![SV-9: three 2:1 couplings against N, and the multigrid solver's errors on AMR](/flume/selfgravity-amr.png)

*Left: the truncation error of S1 grows as $1/h$; M's stays at about 40, an $O(1)$ seam truncation that is still
consistent (the classical supra-convergence of cell-centred schemes); B's decays in this model. Middle: the potentials
of the three couplings coincide at second order; the seam face gradient converges at second order for B, at first
order for M, and stalls for S1. Right: the multigrid solver on M in 3-D and 2-D: the potential second order, the 2:1
face gravity first order, the truncation bounded.*

**The operator study** (measured, no bounds):

| Coupling, $N$ = 8 → 48 | $\lVert A - A^T\rVert/\lVert A\rVert$ | $\max\lvert\tau\rvert$ | seam face gradient error, order (32 → 48) | $\Phi$ order (32 → 48) |
|---|---|---|---|---|
| B | $6.8\cdot10^{-3}$ | 18.1 → 3.2 | $4.1\cdot10^{-2}$ at $N = 48$, 1.93 | 1.99 |
| S1 | 0 | 30.7 → 107, growing as $1/h$ | 0.98 at $N = 48$, **0.11: does not converge** | 1.99 |
| M | 0 | 37.4 → 41.0, bounded | 0.42 at $N = 48$, 1.04 | 2.00 |

**The solver** (CPU only):

| Asserted | Bound | CPU |
|---|---|---|
| SV-9a: multigrid against the FFT solver, $\max\lvert\Phi_{mg} - \Phi_{fft}\rvert / \max\lvert\Phi\rvert$; FCG iterations | $10^{-8}$; $\le 1$ | $\le 3.0\cdot10^{-15}$; 1 |
| SV-9b: FCG residual, relative max norm | $10^{-10}$ | below it |
| symmetry defect of the assembled rows (not symmetrised by the assembly) | $10^{-12}$ | $\le 7.9\cdot10^{-18}$ |
| ghost cells left unfilled | 0 | 0 |
| order of $\Phi$ against the exact potential, base level 2 → 3 | $\ge 1.7$ | 1.92 (3-D), 1.92 (2-D) |
| growth of $\max\lvert\tau\rvert$ under refinement | $\le 1.5$ | 1.006 (3-D), 1.005 (2-D) |
| SV-12a (Dirichlet, uniform): FCG iterations; order of $\Phi$, level 1 → 2 | $\le 1$; $\ge 1.7$ | 1; 1.95 |
| SV-12a: growth of $\max\lvert A\Phi_{exact} + Bg - f\rvert$ | $\le 1.5$ | 0.27 (it falls at second order, 0.64 → 0.17) |
| SV-12b (Dirichlet, AMR): FCG residual; symmetry defect of the Dirichlet rows; ghosts unfilled | $10^{-10}$; $10^{-12}$; 0 | below it; $\le 6.7\cdot10^{-18}$; 0 |
| SV-12b: order of $\Phi$; growth of $\max\lvert\tau\rvert$ | $\ge 1.7$; $\le 1.5$ | 1.91 (3-D), 1.93 (2-D); 1.006, 1.005 |
| SV-13a (isolated, uniform): $\max\lvert\Phi_{mg} - \Phi_{fft}\rvert / \max\lvert\Phi\rvert$ against the FFT solver's isolated mode; FCG iterations of each solve | $10^{-9}$; $\le 1$ | $\le 2.6\cdot10^{-13}$; 1 + 1 |
| SV-13a: order of $\Phi$ against the exact potential, $16^3 \to 32^3$ | $\ge 1.7$ | 2.14 |
| SV-13b (isolated, AMR): FCG residual (tolerance $10^{-12}$); ghosts unfilled | below it; 0 | below it; 0 |
| SV-13b: order of $\Phi$, base level 2 → 3 | $\ge 1.7$ (polytrope), $\ge 1.2$ (uniform sphere, as SV-6) | 2.01, 1.86 |
| SV-13: symmetry of the solution map in situ, $\lvert x\cdot VKy - y\cdot VKx\rvert / \lvert x\cdot VKy\rvert$ for two pseudo-random sources, first grid of each leg | $10^{-9}$ | $\le 1.7\cdot10^{-11}$ |
| E1b (boundary kinds, uniform): $\max\lvert\Phi_{mg} - \Phi_{fft}\rvert / \max\lvert\Phi\rvert$ against the FFT solver with the same kinds; FCG iterations | $10^{-8}$; $\le 1$ | $\le 1.2\cdot10^{-14}$; 1 |
| E1b (AMR): FCG residual; symmetry defect; ghosts unfilled | below the tolerance; $10^{-12}$; 0 | below it; $\le 2.5\cdot10^{-17}$; 0 |
| E1b: order of $\Phi$ against $\Phi_m$, base level 2 → 3; order of the domain-face ghosts against $\Phi_m$ at the ghost centres; growth of $\max\lvert\tau\rvert$ | $\ge 1.7$; $\ge 1.7$; $\le 1.5$ | $\Phi$ 1.95 (P/P D/N N/D), 2.00 (D/D N/N G/G), 1.95 (P/P N/N N/N, singular), 2.00 (2-D N/D P/P); ghosts 1.95, 2.03, 1.95, 2.00; 1.00–1.02 |
| mixed periodic and bounded directions without `bc`; `bc` against the grid's periodicity; a ghost-Dirichlet face against a face kind on one axis | refused, exit 1 | refused |
| E5c (`--sub`): both solves below their tolerance; $\max\lvert\Phi_{sub} - \Phi_{exact}\rvert / \max\lvert\Phi\rvert$; FCG iterations against the exact coarsest level | below it; $10^{-7}$; $\le 1.5\times$ (20 where the exact needs 1) | below it; $1.6\cdot10^{-14}$ - $1.1\cdot10^{-11}$; uniform 14 (1), AMR periodic 47 (36), quadtree 25 (24), D/D N/N G/G 41 (39), P/P D/N N/D 53 (40), isolated 52 (47) |
| isolated, a leaf above the coarsest level on the boundary | refused, exit 1 | refused |

Reported, not bounded:
- **2:1 face gravity order:** 1.03 (3-D) and 1.11 (2-D).
- **FCG iterations:** 36–39 in 3-D (mean factor 0.52–0.54) and 24–25 in 2-D (0.37–0.40); with Dirichlet faces
  35–36 (3-D) and 22–23 (2-D); isolated (SV-13b, tolerance $10^{-12}$) 24–25 for the zero-ghost solve and 23 for the
  second; with the E1b boundary kinds 31–40 (3-D) and 21 (2-D).
- **E1b truncation:** $\max\lvert\tau\rvert = 12.0$ (3-D) and 10.3 → 10.5 (2-D), the same in every configuration,
  with or without Dirichlet faces: the maximum sits on the 2:1 seams (the $O(1)$ seam truncation of SV-9b for this
  potential), not on the boundary rows, so the growth bound does not test them; the boundary rows and weights are
  pinned by the orders of $\Phi$ and of the ghosts. Mutations, each alone, turn the suite red: the Neumann weight in
  $B$ negated (every AMR configuration, orders 0.02–0.14), the Neumann ghost's $h g$ negated (ghost orders 0.99–1.01),
  the Dirichlet ghost $2g + \Phi_C$ instead of $2g - \Phi_C$ (the three configurations with a Dirichlet face, ghost
  orders 0.02–0.03). Not covered: the preconditioner's boundary diagonal, which can change only the iteration count.
- **The James ghost values** against the exact potential at the ghost centres: $9.1\cdot10^{-4} \to 6.2\cdot10^{-5}$
  (SV-13a), $1.0\cdot10^{-4} \to 2.5\cdot10^{-5}$ (SV-13b polytrope), $2.8\cdot10^{-4} \to 8.9\cdot10^{-5}$
  (sphere).
- **AMR against the uniform grid at the same base level** ($32^3$, polytrope): $1.18\cdot10^{-3}$ against
  $2.31\cdot10^{-3}$.

Every printed number is the same on 1, 2 and 3 ranks, up to round-off. The bound on the potential's order is 1.7, not
1.9, because two grids give a single pair.

![SV-13: the isolated solve on AMR, errors against the finest cells per side, and its exactness measures](/flume/selfgravity-amr-isolated.png)

*Left: the potential (solid) and the James ghost values (dashed) against the exact potential, uniform and with the
star refined: second order, the ghost values converging with the potential. Right: the multigrid against the FFT
solver's isolated mode on uniform grids (the same discrete free-space solution, to $2.6\cdot10^{-13}$, the lattice
Green's function's own accuracy), and the symmetry defect of the solution map, within 20 times the FCG tolerance, against
the bound. CPU, 2 ranks.*

Related cases: the FFT solver the coarsest level uses is [SV-2](./fft-poisson), its isolated mode
[SV-6](./isolated-solver); the solver coupled to FLUME on a static AMR grid is [SV-10](./amr-jeans), under runtime AMR
[SV-11](./runtime-amr).

## Findings

- **S1 is refuted.** The P5-0 design comment claimed the symmetrisation costs nothing. Its potential is second order
  like the others, but its truncation grows as $1/h$ and the gravity at the seam does not converge. The potential's
  max-norm hides the defect; the gravity, the quantity the scheme uses, shows it.
- **M pays for symmetry with one order of seam gravity.** Its seam face gradient is first order (B's is second order
  in this model), and 10 times B's at $N = 48$. The error sits on the one layer of seam faces; the potential stays
  second order.
- **M's truncation is $O(1)$ at the seam but consistent:** bounded under refinement, the classical supra-convergence
  of cell-centred schemes.
- The multigrid runs on the host only (CPU); the FNL backend calls it with its coarsest level on cuFFT (SV-10).
- **The levels inside the blocks trade iterations for the rank-0 gather** (E5c). Forced to one cell per block they
  add 1 to 13 FCG iterations on AMR and 13 on a uniform grid, where the exact coarsest level needs 1: the point
  smoother on SV-2's $1:2:3$ cells (the prototype `src/tests/poisson/mg_subblock_e5c.py` measures the same, 26 for
  white noise, against 8 on cubic cells; trilinear prolongation gains nothing over piecewise constant there). The
  default root of $32^3$ cells on the host keeps the exact coarsest level wherever it is cheap; above it the gather of
  the whole base grid per V-cycle, 70% of the isolated solve at $128^3$, goes: on 3.9 M AMR cells (medooza CPU) the
  isolated solve falls from 39.9 s to 6.5 s at 16 ranks and from 50.7 s to 31.2 s at one, the periodic one from 2.6 s
  to 2.3 s at 16 ranks but rises from 9.2 s to 12.7 s at one: there one FFT of the $128^3$ base level beats smoothing
  it with the generic CSR sweep (a stencil smoother is #84 E6). On FNL the multigrid smooths on the host and cuFFT
  solves the coarsest level fast, so the cuFFT root keeps the exact coarsest level up to $128^3$ cells. A mutation with no off-diagonal
  in the levels below the coarsest fails 5 of the 6 configurations on the iteration bound.
- **Dirichlet boundaries cost nothing in symmetry or order** (SV-12): a boundary face adds $-w\Phi_C$ to its row and
  moves $wg$ to the right-hand side, so the rows stay symmetric ($6.7\cdot10^{-18}$) and consistent (the uniform
  truncation falls at second order). The coarsest level's odd extension is exact: one iteration on a uniform grid.
  Mixed boundaries are refused; the periodic legs are unchanged to every printed digit.
- **The James step is symmetric and second order on AMR** (SV-13), the P0 prototype's result in the library: on a
  uniform grid it is the discrete free-space solution (the FFT solver's, to $2.6\cdot10^{-13}$), on AMR the potential
  is second order (2.01 for the polytrope) and the solution map is symmetric to the solve tolerance, which the exact
  $E + W$ exchange needs. The cost is two composite solves and one padded FFT of the coarsest box per Poisson solve.
  The price is D-87-2: every leaf on the domain boundary at the coarsest level, fatal otherwise; the FLUME coupling
  (P3) lifts the run-time refusals.
