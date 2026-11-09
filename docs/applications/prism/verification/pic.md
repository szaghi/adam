# PIC: periodic faces and the Fresnel slab

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `pic-periodic`, `pic-slab-pulse` | [#111](https://github.com/szaghi/adam/issues/111), [#112](https://github.com/szaghi/adam/issues/112) ([#105](https://github.com/szaghi/adam/issues/105) open) | P | CPU, FNL | no | `golden-check/{cpu,fnl}` (digest + residuals, compared by `pic_check.sh`) |

## Case and scope

Two cases run the PIC path of the production runs: TSC particle weighting, the charge-conserving current
`esirkepov_Q10` with the coherent gather, the t0 FLAIL solves of $\mathbf D$ and $\mathbf B$, FD6, SSP (5,4), two
macroparticles per cell per species, a fixed seed and no pairing.

- `pic-periodic`: a uniform plasma periodic on every face, with a plane wave of one box length. It exercises the
  periodic current fold and the particle wrap on three axes ([particle boundaries](../pic#particle-boundaries)).
- `pic-slab-pulse`: the Fresnel pilot at $\Delta x = c/(8\omega_{pe})$. A plasma slab $3\,c/\omega_{pe}$ thick, PEC
  along $x$ and periodic in $y$, $z$; a single-cycle `gaussian_derivative_pulse` crosses the slab.

The charge-conserving CPU current is single-block and single-rank, so each case runs as legs: the CPU on one block for
the first 4 steps; FNL at 1 rank on one block; FNL at 2 ranks on the same grid as 64 blocks. On the FNL virtual block
the multi-block result equals the single-block one by construction, up to the order of the atomic adds, so the 2-rank
leg checks the halos, the migration and the decomposition against the 1-rank leg.

What the cases catch:
- a deposit that does not conserve charge: the charge-conservation history must stay at round-off;
- a t0 field that violates Gauss's law;
- a change of the particle path beyond the digest tolerance (the push, the deposit, the gather, the wrap);
- a decomposition-dependent result on FNL.

`pic-slab-pulse` also carries the one PRISM measurement against an exact physical solution: run to
$t\,\omega_{pe} = 32$ with three probe planes, the same input matches the exact cold slab.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `physical_model = PIC`, electrons and ions, Maxwellian ($T_e = 30000$ K, $T_i = 300$ K), density $5.663678\cdot10^{31}$ m⁻³ | [particle-in-cell](../pic#particles) |
| Weighting | TSC charge; conserving current `esirkepov_Q10`; field gather `esirkepov_Q10` (coherent) | `deposit_matched_current`, [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90); [`adam_prism_fnl_pic_conserving_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/fnl/adam_prism_fnl_pic_conserving_object.F90); [charge-conserving currents](../pic#charge-conserving-currents) |
| Particle faces | periodic → fold and end-of-step wrap; other kinds → legacy | [`adam_prism_pic_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_pic_bc_object.F90) |
| t0 fields | FLAIL solves of $\phi$ and $\mathbf A$; closure `grounded` (`pic-periodic`), `isolated` (`pic-slab-pulse`) | [the fields at t0](../pic#the-fields-at-t0) |
| Time | SSP Runge–Kutta (5,4); CFL 0.5 (`pic-periodic`), 0.25 (`pic-slab-pulse`: at CFL 0.5 $\omega_{\max}\Delta t \approx 0.07$ hits the plasma-frequency limit) | [`adam_prism_plasma_frequencies_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_plasma_frequencies_object.F90) |
| Diagnostic | `charge_conservation_save = 1`: `Q_deposition_error`, `dQ_num_net` per step | [charge conservation diagnostic](../pic#charge-conservation-diagnostic) |

## Initial and boundary conditions

| Item | `pic-periodic` | `pic-slab-pulse` |
|---|---|---|
| Domain (m) | $x \in [0, 2\cdot10^{-9}]$, $y, z \in [-7.5, 7.5]\cdot10^{-10}$ | $x \in [0, 4.8\cdot10^{-8}]$, $y, z \in [-1.5, 1.5]\cdot10^{-9}$ |
| Grid | $32\times24\times24$, one block ($\Delta x = c/(16\omega_{pe})$); 2-rank leg: 64 blocks of $8\times6\times6$ | $384\times24\times24$, one block ($\Delta x = c/(8\omega_{pe})$); 2-rank leg: 64 blocks |
| Plasma | uniform | slab normal to $x$, thickness $3\cdot10^{-9}$ m, centre $x = 2.45\cdot10^{-8}$ m |
| Field initial condition | plane wave along $x$, $\lambda = 2\cdot10^{-9}$ m (one box length) | [Gaussian-derivative pulse](../initial-conditions#gaussian-derivative-pulse) along $+x$, `r0 = 1.5e-8`, `sigma = 1e-9` |
| Faces | `periodic` on all six | `PEC` in $x$, `periodic` in $y$, $z$ |
| Steps | 40 (CPU leg: 4) | 1152, $t\,\omega_{pe} = 18$ (CPU leg: 4) |
| Macroparticles | 2 per cell per species, weight 3.457, seed 12345 | 2 per cell per species, weight 27.65, seed 12345 |

Run times: CPU 1 rank about 2.5 / 2 min; FNL 1 rank about 35 s / 10 min (RTX 4070 laptop); FNL 2 ranks about
1.5 / 10 min (two A30).

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/pic-{periodic,slab-pulse}/input.ini` | one block (CPU and FNL at 1 rank) |
| `src/tests/prism/regression/pic-{periodic,slab-pulse}/input-np2.ini` | the same grid as 64 blocks (FNL at 2 ranks) |
| `src/tests/prism/regression/pic-{periodic,slab-pulse}/check.sh <cpu\|fnl> [--capture]` | call the shared driver |
| `src/tests/prism/regression/pic_check.sh` | the legs: charge, Gauss at t0, digest, residuals |
| `src/tests/prism/regression/pic-*/golden-check/{cpu,fnl}/` | the goldens (not `golden/`: `run.sh` runs a goldened case at 2 ranks, which the CPU leg refuses) |
| `src/tests/prism/fnl/fd/PIC/Fresnel_test/pilot/pulse/input.ini` | the Fresnel pilot the slab case derives from |

```bash
cd src/tests/prism/regression/pic-periodic && ./check.sh cpu
FNL_NP=1 ./check.sh fnl                     # FNL, 1 rank (WSL); FNL_NP="1 2" on a cluster node
PIC_CHECK_LAUNCH="srun --mpi=pmix -u -n" ./check.sh fnl   # under Slurm
```

## Results

Each leg asserts, as coded:

| Check | Bound | `pic-periodic` CPU | `pic-periodic` FNL | `pic-slab-pulse` CPU | `pic-slab-pulse` FNL |
|---|---|---|---|---|---|
| run exits 0 | — | pass, np 1 | pass, np 1, 2 | pass, np 1 | pass, np 1, 2 |
| charge: max $\lvert$`Q_deposition_error`$\rvert$, $\lvert$`dQ_num_net`$\rvert$ over the run, relative to one species' charge | $\le 10^{-12}$ | $1.0\cdot10^{-16}$ | $7.3\cdot10^{-18}$ | $9.3\cdot10^{-17}$ | $1.3\cdot10^{-17}$ |
| Gauss at t0: max $\lvert\nabla\cdot\mathbf D - \rho + c\chi\rvert/\max\lvert\rho\rvert$ | $\le 10^{-12}$ | $8.9\cdot10^{-15}$ | $5.0\cdot10^{-15}$ | $6.3\cdot10^{-15}$ | $5.2\cdot10^{-15}$ |
| digest against `golden-check/<backend>` | `rtol = 1e-6`, `atol = 1e-3`; sums on the scale $\sqrt{n\sum q^2}$ | 48 rows, 0 outside | 48 + 42 rows (np 1, 2), 0 outside | 48 rows, 0 outside | 96 + 84 rows (np 1, 2), 0 outside |
| residuals against the golden (1-rank legs) | `rtol = 1e-6`, `atol = 1e-3` | 1 row, match | 1 row, match | 1 row, match | 3 rows, match |

At 2 ranks the digest skips the per-block metadata rows (`origin`, `dxdydz`, `time_iteration`) and the residuals log,
which depend on the block layout, not on the fields.

The Fresnel slab against the exact cold slab, plane-mean $L_2$ error (`pic-slab-pulse` input run to
$t\,\omega_{pe} = 32$ with three probe planes; not a leg):

| Probe plane | `pic-slab-pulse` input | Pilot plane study, same $\Delta x$ |
|---|---|---|
| reflected | 3.3 % | 2.7 % |
| inside the slab | 3.3 % | 4.2 % |
| transmitted | 6.2 % | 7.0 % |

The error is second order in $\Delta x$.

Related pages: [t0 elliptic solves](./t0-elliptic) for the t0 fields on their own; the comparison with the exact
spectral solution is open as [#105](https://github.com/szaghi/adam/issues/105).

## Findings

- The 2-rank FNL leg is compared with the 1-rank golden: the virtual block makes the multi-block result equal to the
  single-block one up to the order of the atomic adds, so any halo, migration or decomposition defect shows.
- The digest compares the sum of a zero-mean field ($\mathbf J$, $\rho$, the residuals) on its own scale,
  $\sqrt{n\sum q^2}$: its last digits follow the order of the atomic adds.
- The slab pulse runs at CFL 0.25: at CFL 0.5 the plasma-frequency time-step limit would bind. The case doubled its
  step count for it (576 → 1152).
- The Fresnel slab matches the exact cold slab to 3.3 / 3.3 / 6.2 %, second order in $\Delta x$.
