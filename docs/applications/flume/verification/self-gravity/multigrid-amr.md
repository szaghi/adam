# The multigrid Poisson solver on AMR grids

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-9, SV-12 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P5a; [#87](https://github.com/szaghi/adam/issues/87), P1 | CPU only | yes | — |

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
  (SV-12b).

What the case catches:
- a non-symmetric assembled operator (the symmetry defect of the assembled rows, which the assembly does not
  symmetrise);
- an inconsistent seam coupling: a truncation error that grows as $1/h$ under refinement;
- a multigrid that is not the FFT solve on a uniform grid, a ghost the fill leaves empty, a stalled FCG;
- a solver that accepts a non-periodic grid.

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
| Coarsest level | the minimum leaf level, which covers the domain, solved exactly by the FFT solver | `solve_blocks`, [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |
| Scope | periodic only (P5a deferred isolated boundaries on AMR to P5b; they remain refused, [#87](https://github.com/szaghi/adam/issues/87)); 2:1-balanced faces; anything else is fatal | `initialize` (periodicity), `classify` (2:1 balance), [`adam_poisson_mg_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_mg_object.F90) |
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
| Refusal | `--mixed`: one uniform grid periodic along $x$ only (every active direction periodic, or none) | — | 2 |

The FNL path comes with FLUME on AMR (P5b, P5c: [SV-10](./amr-jeans), [SV-11](./runtime-amr)).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/selfgravity/seam_operators.py` | the operator study (`[N ...]`, default 8 16 32 48) |
| `src/tests/poisson/test_poisson_mg.F90` | the unit test (mode `test-poisson-mg-gnu`), SV-9a, SV-9b, SV-12a, SV-12b; `--mixed` the refusal |
| `src/tests/flume/verification/selfgravity/poisson.sh` | leg `sv9`: the test on every rank count of `--np-list`, and the refusal; `POISSON_MG_EXE` selects the build |
| `src/tests/flume/verification/make_doc_runs.sh`, group `amr` | `sv9-poisson`: the test with its log |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_amr` | draws the figure (`--only selfgravity-amr`), running `seam_operators.py` |

```bash
python3 src/tests/flume/verification/selfgravity/seam_operators.py      # the three couplings
cd src/tests/flume/verification/selfgravity && ./poisson.sh --legs sv9   # the solver, np 1, 2, 3
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
| mixed periodic and non-periodic directions | refused, exit 1 | refused |

Reported, not bounded:
- **2:1 face gravity order:** 1.03 (3-D) and 1.11 (2-D).
- **FCG iterations:** 36–39 in 3-D (mean factor 0.52–0.54) and 24–25 in 2-D (0.37–0.40); with Dirichlet faces
  35–36 (3-D) and 22–23 (2-D).

Every printed number is the same on 1, 2 and 3 ranks, up to round-off. The bound on the potential's order is 1.7, not
1.9, because two grids give a single pair.

Related cases: the FFT solver the coarsest level uses is [SV-2](./fft-poisson); the solver coupled to FLUME on a static
AMR grid is [SV-10](./amr-jeans), under runtime AMR [SV-11](./runtime-amr).

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
- **Dirichlet boundaries cost nothing in symmetry or order** (SV-12): a boundary face adds $-w\Phi_C$ to its row and
  moves $wg$ to the right-hand side, so the rows stay symmetric ($6.7\cdot10^{-18}$) and consistent (the uniform
  truncation falls at second order). The coarsest level's odd extension is exact: one iteration on a uniform grid.
  Mixed boundaries are refused; the periodic legs are unchanged to every printed digit.
