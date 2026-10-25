# The FFT Poisson solver

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SV-2 | Self-gravity | M8 ([#81](https://github.com/szaghi/adam/issues/81)), P2 | CPU, FNL | yes | — |

## Case and scope

Self-gravity needs a Poisson solver that is exact for the *discrete* equation $L_h\Phi = 4\pi G(\rho - \bar\rho)$,
with $L_h$ the 7-point Laplacian the gravity stencils telescope against, so that the scheme is the one P0 analysed.
P2 wrote this as the condition of exact conservation; P4 found the condition weaker: any linear, symmetric solve
conserves $E + \tfrac12\rho\Phi$ and momentum exactly in the semi-discrete form (see [SV-6](./isolated-solver)). The
round-off residual remains a property of this solver, verified here.

What the case catches (the unit test `test_poisson_fft`):
- a wrong eigenvalue of $L_h$ (a missing factor, a swapped axis, a wrong spacing): the error against the discrete
  solution leaves round-off. The domain is anisotropic, $[0,1]\times[0,2]\times[0,3]$, so an axis mix-up cannot cancel;
- a misplaced block in the gather, or a wrong ghost wrap in the scatter: the error over every cell and ghost, and the
  residual from the scattered ghosts;
- a ghost layer the scatter leaves unfilled (a sentinel value);
- a solver that accepts a grid it cannot solve: a non-periodic grid must stop the run.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equation | $L_h\Phi = 4\pi G(\rho - \bar\rho)$, $L_h$ the 7-point Laplacian (5-point in 2-D) | [models](../../models#self-gravity) |
| Solver, CPU | FFTW, the distributed slab transform ([#84](https://github.com/szaghi/adam/issues/84) E5d): each rank owns z-planes of the transform array, transforms them in 2-D, the planes are transposed to z-lines by `MPI_ALLTOALLV` and transformed in 1-D, and back; one plan per plane and per row of lines, so the potential is bitwise the same on every rank count (the `slab` leg of `poisson.sh`); each Fourier mode is divided by the eigenvalue of $L_h$ itself, $\lambda_k = \sum_d (2\cos(2\pi k_d/N_d) - 2)/\Delta x_d^2$; the mean ($\lambda_k = 0$) is set to zero; the inverse transform follows, then every block fetches its cells with all their ghost layers from the owners of its planes, by the periodic wrap | `initialize`, `solve`, `transform`, [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |
| Solver, FNL | the same slab transform on each rank's GPU, its arrays resident there: cuFFT for the planes and the lines (one plan per plane and per row: bitwise equal across rank counts too), the transposes packed on the device and exchanged by GPU-aware `MPI_ALLTOALLV`, the division on the device | `slab_transpose_fnl`, `slab_spectral_fnl`, `transform_fnl`, `spectral_solve_dev`, [`adam_fnl_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_poisson_fft_object.F90) |
| Scope | uniform, periodic grids; anything else stops the run. Until #84 E5d the solve gathered the grid to rank 0 (serial, not scalable); since then it is distributed over the z-planes (at most $m_z$ ranks take part). Isolated boundaries come in P4 ([SV-6](./isolated-solver)), AMR in P5 ([SV-9](./multigrid-amr)) | `initialize`, [`adam_poisson_fft_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_poisson_fft_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0,1]\times[0,2]\times[0,3]$, periodic on every face |
| Source | $\rho = 1.5 + \sin(k_1x)\cos(k_1y)\cos(k_1z) + 0.5\cos(\mathbf k_2\cdot\mathbf x)$, $\mathbf k_1 = 2\pi(1, 1, \tfrac13)$, $\mathbf k_2 = 2\pi(2, \tfrac12, \tfrac13)$; in 2-D the z wave numbers are 0 |
| Factor | $4\pi G$ with $G = 1$ |
| Grids | 3-D octree at levels 1–3 ($16^3$ to $64^3$ cells); 2-D quadtree (`null_z`) at levels 1–3 ($16^2$ to $64^2$); a single $8^3$ block (level 0: the wrap inside one block; on 2+ ranks, [ranks without blocks](../infrastructure/ranks-without-blocks)). 8 cells per block side, 3 ghosts |
| Ranks | CPU: 1, 2, 3 (the `poisson.sh` default list), and 4; FNL: 1, 2 |
| Exact solutions | each mode is an eigenfunction of $L_h$: the **discrete** solution divides it by its $\lambda_k$, the **continuum** one by $-\lvert\mathbf k\rvert^2$ |
| Refusal | `--no-periodic`: one 3-D configuration without periodic boundaries; the run must stop with exit 1, naming the direction (`direction 1 is not periodic`) |

## Input files

| File | Role |
|---|---|
| `src/tests/poisson/test_poisson_fft.F90` | the unit test (modes `test-poisson-fft-gnu`, `test-poisson-fft-fnl-nvf`); `--dump DIR` writes rank 0's potential of every level |
| `src/tests/flume/verification/selfgravity/poisson.sh` | runs the test on every rank count of `--np-list` (leg `sv2`), and the refusal; `POISSON_EXE` selects the build |
| `src/tests/flume/verification/selfgravity/poisson_oracle.py` | `compare DIR_CPU DIR_FNL`: the CPU = FNL check on the dumped potentials, and the figure's reader |
| `src/tests/flume/verification/make_doc_runs.sh`, group `selfgravity` | `sv2-poisson`: the test with its log and dumps |
| `src/tests/flume/verification/make_doc_figures.py`, `fig_selfgravity_poisson` | draws the figure (`--only selfgravity-poisson`) |

```bash
cd src/tests/flume/verification/selfgravity
./poisson.sh --legs sv2 --dump <dir-cpu>                                    # CPU (FFTW)
POISSON_EXE=$PWD/../../../../../exe/test_poisson_fft_fnl ./poisson.sh --legs sv2 --np-list "1 2" --dump <dir-fnl>
python3 poisson_oracle.py compare <dir-cpu> <dir-fnl>                       # CPU = FNL
```

## Results

![SV-2: the errors against N, the 2-D potential and its error against the continuum solution](/flume/selfgravity-poisson.png)

*Left: against N, the error against the continuum solution falls as $\Delta x^2$, while the error against the discrete
solution and the residual stay at round-off. Middle and right: the 2-D potential at $N = 64$ and its error against
the continuum solution, $4.3\cdot10^{-4}$, largest where $\lvert\Phi\rvert$ is.*

| Asserted | Bound | CPU (FFTW) | FNL (cuFFT, A30 and the WSL GPU) |
|---|---|---|---|
| $\max\lvert\Phi - \Phi_h\rvert / \max\lvert\Phi_h\rvert$, every cell and ghost | $10^{-12}$ | $\le 2.8\cdot10^{-15}$ | $\le 2.8\cdot10^{-15}$ |
| residual of $L_h$ from the scattered ghosts, relative | $10^{-10}$ | $\le 3.7\cdot10^{-14}$ | $\le 4.3\cdot10^{-14}$ |
| order against the continuum solution, $N = 32 \to 64$ | $\ge 1.9$ | 1.98 (3-D), 1.98 (2-D) | 1.98, 1.98 |
| ghost cells the scatter left unfilled | 0 | 0 | 0 |
| a non-periodic grid | refused, exit 1 | refused | refused |
| CPU against FNL, $\max\lvert\Phi_{\text{CPU}} - \Phi_{\text{FNL}}\rvert / \max\lvert\Phi\rvert$ (`poisson_oracle.py compare`) | $10^{-13}$ | $\le 9.8\cdot10^{-15}$ | $\le 9.8\cdot10^{-15}$ |

Related cases: the isolated variant of the same solver is [SV-6](./isolated-solver); the multigrid solver on AMR grids,
checked against this one on uniform grids, is [SV-9](./multigrid-amr); the coupling to the scheme is
[SV-3 to SV-5](./jeans).

## Findings

- The printed numbers do not depend on the rank count (np 1 and np 4 agree to every digit), and the A30 potential is
  bitwise that of the WSL GPU.
- The solver is exact for $L_h$, and the scheme is second order through $L_h$.
- The residual grows from $10^{-15}$ to $4\cdot10^{-14}$ between $N = 16$ and 64: $L_h$ amplifies the round-off of
  $\Phi$ by $(k\Delta x)^{-2}$. That is the floor the conservation of P3 inherits.
- The solve is distributed over the z-planes since #84 E5d (it was serial on rank 0) and its result does not depend on
  the rank count, bitwise.
