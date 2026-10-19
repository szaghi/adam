# Self-gravity (M8)

Self-gravity ([#81](https://github.com/szaghi/adam/issues/81)) solves $\nabla^2\Phi = 4\pi G(\rho - \bar\rho)$ every
Runge–Kutta stage and applies gravity through the cell potential. The design (the P0 comment on #81) uses the spatial
discretisation of Mullen et al. (2021), with the 7-point Laplacian $L_h$:
- **face gravity:** $g_{d,i+1/2} = -(\Phi_{i+1} - \Phi_i)/\Delta x_d$;
- **cell gravity:** the mean of the cell's two faces;
- **momentum:** $\rho$ times the cell gravity;
- **energy:** $\tfrac12\sum_d (g_{d,i+1/2} F^\rho_{d,i+1/2} + g_{d,i-1/2} F^\rho_{d,i-1/2})$, from the mass fluxes the
  continuity equation uses.

Summed over the cells, the energy source is $\sum_i \Phi_i (F_{i+1/2} - F_{i-1/2})/\Delta x$. So $E + \rho\Phi$ is
conserved exactly with a static potential, and $E + \tfrac12\rho\Phi$ with a self-consistent one: in the semi-discrete
form, given a Poisson solve that is linear and symmetric, $\Phi = 4\pi G\,K\rho$ with $K = K^T$. P4 corrected the first
claim, which asked for a round-off residual of the Poisson equation (see [models](../../models#self-gravity)). The
scheme is second order.

The phases of #81 and the cases that verify them:

| Phase | What it adds | Cases |
|---|---|---|
| P1 | the potential path (`[gravity] potential = linear`) | SV-1 |
| P2 | the FFT Poisson solver, periodic | SV-2 |
| P3 | the coupling, one solve per Runge–Kutta stage (`[gravity] self_gravity`) | SV-3, SV-4, SV-5 |
| P4 | isolated boundaries (`poisson_bc = isolated`) | SV-6, SV-7, SV-8 |
| P5a | the multigrid Poisson solver on AMR grids | SV-9 |
| P5b | self-gravity on static AMR grids, the seam gravity work in the reflux | SV-10 |
| P5c | the per-stage reflux, self-gravity under runtime AMR | SV-10 (`reflux_cadence = stage`), SV-11 |
| P6 | the regression cases `selfgravity-jeans`, `selfgravity-amr-stage`, `selfgravity-star` | — |
| [#87](https://github.com/szaghi/adam/issues/87) P1 | Dirichlet faces in the multigrid solver, a DST-I coarsest level | SV-12 |
| [#87](https://github.com/szaghi/adam/issues/87) P2 | isolated boundaries on AMR, the symmetric James step | SV-13 |
| [#87](https://github.com/szaghi/adam/issues/87) P3 | isolated self-gravity on AMR in FLUME: the refusals lifted, a regrid that refines a boundary block fatal | SV-14 |

The figures are regenerated with `make_doc_runs.sh <dir> <group>` and `make_doc_figures.py --only <figure>`
(`src/tests/flume/verification/`):

| `make_doc_runs.sh` group | Executables it needs | `make_doc_figures.py --only` |
|---|---|---|
| `selfgravity` | `exe/adam_flume_cpu`, `exe/test_poisson_fft` (mode `test-poisson-fft-gnu`) | `selfgravity-potential selfgravity-poisson` |
| `jeans` | `exe/adam_flume_cpu_fft` (mode `flume-cpu-gnu-fft`) | `selfgravity-jeans` |
| `isolated` | `exe/test_poisson_fft`, `exe/adam_flume_cpu_fft` | `selfgravity-isolated selfgravity-star` |
| `amr` | `exe/test_poisson_mg` (mode `test-poisson-mg-gnu`), `exe/adam_flume_cpu_fft` | `selfgravity-amr selfgravity-amr-jeans selfgravity-amr-regrid` |
| `sinks` | `exe/adam_flume_cpu_fft` | `sink-particles` |

| Case | ID | What it verifies |
|---|---|---|
| [Gravity from a cell potential](./potential) | SV-1 | the potential path reproduces the constant-g oracles of M7 |
| [The FFT Poisson solver](./fft-poisson) | SV-2 | the periodic solver is exact for $L_h$, second order against the continuum |
| [Jeans modes](./jeans) | SV-3 to SV-5 | the coupling: dispersion relation, growth rate, symmetry, exact $E + W$ exchange, scaling covariance |
| [The isolated solver](./isolated-solver) | SV-6 | zero padding with the lattice Green's function: free-space potential of the mass on the grid |
| [The n = 1 polytrope](./polytrope) | SV-7, SV-8, SV-14 | a star in equilibrium and a collapsing one under isolated self-gravity, uniform and on AMR |
| [The multigrid solver on AMR grids](./multigrid-amr) | SV-9, SV-12, SV-13 | the mimetic 2:1 coupling (symmetric and consistent), the FCG solver, Dirichlet faces, isolated boundaries by the symmetric James step |
| [Jeans modes on a static AMR grid](./amr-jeans) | SV-10 | the multigrid coupling, the seam gravity work, the per-stage reflux |
| [Self-gravity under runtime AMR](./runtime-amr) | SV-11 | the solver and $\Phi$ follow every regrid; $E + W$ exact between regrids |
| [Sink particles](./sink-particles) | SI | point masses through the Poisson source: rest, Kepler binary, creation, accretion at Shu's rate, mergers, runtime AMR ([#94](https://github.com/szaghi/adam/issues/94)) |
