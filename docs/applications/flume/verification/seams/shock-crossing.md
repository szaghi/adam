# Shock crossing a static 2:1 face

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| AS (A3, C3, C3b, C3c) | Seam accuracy and divergence | M6 ([#68](https://github.com/szaghi/adam/issues/68)) | CPU, FNL (Sod); FNL only (RJ2a, C3, C3c) | yes | — |

## Case and scope

The Sod problem ([V1](../euler/sod)) and RJ2a ([MV-4](../mhd/ryu-jones-2a)) run on the V1 grid with one 2:1 face at
$x = 0.75$. The right shock (RJ2a: the fast shock) crosses the face at $t = 0.143$ (Sod) and $t = 0.110$ (RJ2a) and is
past it at $t = 0.2$, so the face then sits in the constant state behind the shock, the plateau, where every deviation
from the exact solution is error with a known reference.

The committed shock-across-a-seam checks (the [multirealm](../euler/multi-realm) legs, the `sod-amr` goldens) compare
FLUME runs with each other, so an error common to all of them passes. AS compares with the exact solution. A seam can
fail a shock in three ways a run-to-run comparison does not see.

What the case catches:
- **a conservation defect:** the mass integral against the exact one, and the shock landing in the wrong place;
- **a reflection:** acoustic waves left in the plateau behind the shock (velocity and pressure);
- **an entropy glitch:** a density and entropy defect at constant pressure, advected away from the face by the flow;
- **a transverse dependence** of the 1-D solution (the spread among the transverse copies must stay at round-off).

AS is the shock baseline of the seam remedies of #68 C, D and E (D2's "no worse than tricubic"). C3 reruns it under the
6-point restriction, C3b profiles the FNL cost, C3c adds the WENO restriction.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations, Sod | compressible Euler, ideal gas, $\gamma = 1.4$ | [models](../../models#compressible-euler) |
| Equations, RJ2a | ideal MHD, $\gamma = 5/3$, no divergence control (`divergence_control = none`; in 1-D $B_n$ is uniform) | [models](../../models#ideal-mhd) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `compute_face_split_fluxes`, [`adam_flume_euler_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_library.F90); `mhd_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.5 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Fine-to-coarse ghosts | `mean` (A3); `lagrange6` and `weno6` (C3, C3c) | `update_ghost_local`, `restrict_seam_high_order`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); `seam_restriction_weno6`, `seam_weno6_core`, `seam_restriction_meta_pack`, [`adam_seam_interpolation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_interpolation_library.F90) |
| Coarse-to-fine ghosts | `tricubic` | `interp_seam_ghost`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90) |
| Reflux | Berger–Colella, on; off only as a measurement (`--no-reflux`) | `accumulate_seam_fluxes`, `apply_reflux_to_stage_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Boundaries | `extrapolation` on all six faces | [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant) |
| FNL ghost kernel | the local ghost update, its map loops collapsed (C3b) | `update_ghost_local_gpu_dev`, [`adam_fnl_field_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/fnl/adam_fnl_field_kernels.F90) |

**C3c, the WENO restriction** (`[amr] seam_restriction = weno6`, #68 C3). Per direction, the WENO-Z combination
(Borges et al. 2008, $p = 2$) of the cubic interpolants of the 4-point sub-stencils of the 6 nodes that contain the
coarse centre, with Jiang–Shu indicators and $\tau$ the highest undivided difference of the candidates' union. An
interval whose difference exceeds twice its neighbours' is a jump; a sub-stencil holding one is dropped, and with none
left (a jump inside the fine pair) the restriction is the mean of the pair, the troubled-cell fallback. No sub-stencil
extrapolates (across a jump it would import the far side's state), so the design order is 6 centred, 5 one cell from
the donor block's edge and 4 at it. It is applied direction by direction ($x$, $y$, $z$) on CPU, FNL and the
inter-realm rows; the regime travels in the row's meta (1..216 Lagrange, 217..432 WENO), so the inter-realm packers need
no new argument. The unit test `test_seam_interpolation`
([`src/tests/amr/test_seam_interpolation.F90`](https://github.com/szaghi/adam/blob/develop/src/tests/amr/test_seam_interpolation.F90))
checks it: exact on cubics at any scale and offset, design orders on smooth data, the fine pair's own state across a
jump elsewhere and the mean of a jump inside the pair, mirror-symmetric on random data.

## Initial and boundary conditions

| Item | Sod | RJ2a |
|---|---|---|
| Domain | $[0, 1]$ along $x$; $y$ and $z$ null | the same |
| Base grid | octree (`ratio = 8`), the V1 grid: four blocks along $x$ of $N/4 \times 4 \times 4$ cells (64 base blocks), `ngc = 3`; $N = 200$ | the same, $N = 256$ |
| Refined side, `c2f` | $x > 0.75$ (one block column): the shock goes from coarse to fine | the same |
| Refined side, `f2c` | $x < 0.75$ (three block columns): the shock goes from fine to coarse | the same |
| Refinement | one 2:1 level at initialisation (`primitive-box` marker, `max_level = 3`, no regridding); the octree refines the null $y$ and $z$ too | the same |
| Uniform controls | $N$ and $2N$ (the control of `c2f` is $2N$, of `f2c` $N$: the spacing of the side that holds the shock at $t = 0.2$) | the same |
| `--long` | adds the refined runs at $2N$ against the controls at $2N$ and $4N$ | the same |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions), $(\rho, u, p) = (1, 0, 1)$ for $x \le 0.5$, $(0.125, 0, 0.1)$ for $x > 0.5$ | `riemann-problem`, the RJ2a states of [MV-4](../mhd/ryu-jones-2a), discontinuity at $x = 0.5$ |
| Boundaries | `extrapolation` on all faces | the same |
| Final time | $t = 0.2$ | $t = 0.2$ |
| Ranks | 2 | 2 |
| C3 variant | `--seam-restriction R`: the null $y$, $z$ blocks widened from 4 to 8 cells (the 6-point footprint needs 6), for every value of R so that a comparison runs on one grid | the same |

Cost, 2 ranks: Sod about 7 min on the CPU and 8 min on FNL; RJ2a about 25 min on FNL and hours on the CPU (the octree
refinement holds 88 and 200 blocks per rank against 32), so the RJ2a leg runs on FNL; the refined RJ2a runs alone take 5
to 14 minutes there.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/amr-shock/make_amr_shock.py` | writes the inputs from `sod/sod-x.ini` (Sod) or `mhd/rj2a/make_rj2a.py` (RJ2a); `--side c2f` or `f2c`, `--uniform`, `--no-reflux`, `--seam-restriction`, `--reflux-policy`, `--it-save` |
| `src/tests/flume/verification/mhd/amr_box.py` | `refine_box`: the refined side, one 2:1 level at initialisation |
| `src/tests/flume/verification/amr-shock/amr_shock_oracle.py` | mass, shock position, plateau, band and spread against the exact solution |
| `src/tests/flume/verification/make_doc_runs.sh <runs> seams`, `make_doc_figures.py --only shock-crossing` | re-runs the Sod crossings `as-c2f`, `as-f2c` and the uniform controls with checkpoints kept (CPU); `rj2a-shock` (FNL) the RJ2a crossings `as-rj2a-{c2f,f2c,uniform}` |
| `src/tests/flume/verification/amr-shock/check.sh` | runs the legs (`--legs sod,rj2a`, `--long`, `--np N`, `--keep`) |
| `src/tests/flume/verification/mhd/positivity/check.sh` | the PV-5 blast across 2:1 seams (C3; `--legs pv5`, not a default leg) |
| `src/tests/amr/test_seam_interpolation.F90` | the unit test of the seam interpolants and restrictions (C3c) |

```bash
cd src/tests/flume/verification/amr-shock && ./check.sh --legs sod   # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh          # FNL, both legs
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --long   # adds the refined runs at 2N
```

The C3 and C3c runs reuse `make_amr_shock.py --seam-restriction R`; the C3 PV-5 runs are the
[positivity](../mhd/positivity-limiter) PV-5 blast with `nk` widened to 8.

**Oracle** (`amr_shock_oracle.py`):
- **mass:** the relative defect of the mass integral against the initial mass plus $t$ times the net inflow at the
  ends (their states never change by $t = 0.2$);
- **shock:** the mass-equivalent shock position (the step holding the mass of the cells within $W = 0.03$ of the
  exact shock, so a smeared profile counts as its integral) minus the exact one, in local cells;
- **plateau:** between the left wave (Sod: the contact; RJ2a: the rotational discontinuity) and the shock, $W$ away
  from both: $\max|q - q_2| / |q_1 - q_2|$ for $\rho$, $u_n$ and $p$ ($q_1$ ahead of the shock, $q_2$ behind it), and
  $\max|s/s_2 - 1|$ of $s = p/\rho^\gamma$ and its mean;
- **band:** the same maxima over the cells within `ngc` of their own spacing from the face;
- **spread:** the largest relative difference among the transverse copies of a column.

$W$ is a length, the same for every run (6 cells of the Sod $N = 200$ grid): a margin in cells of the run would let the
contact, wider in cells on the fine controls, into their plateau.

| Bound | Value |
|---|---|
| mass defect | $\le 10^{-8}$ (on RJ2a at $N = 256$ the uniform run itself leaves $3\cdot10^{-10}$: numerical precursors of the left fast rarefaction reach the extrapolated end $x = 0$; $9\cdot10^{-15}$ at $N = 512$) |
| $\lvert$shock$\rvert$ | $\le 0.5$ local cells |
| spread | $\le 10^{-12}$ |
| every plateau and band metric | at most the baseline ($N$ = 200 / 256) plus 10 % |

## Results

![AS: the Sod shock crossing a 2:1 face](/flume/shock-crossing.png)

*Sod at $t = 0.2$, $N = 200$, today's operators, CPU: the density against the exact solution, coarse to fine and fine to coarse, and the plateau behind the shock with the uniform controls at the coarse ($N = 200$) and the fine ($N = 400$) spacing: the c2f crossing leaves a bump behind the face at $x \approx 0.81$, the f2c one a dip on both sides; neither control has it.*

![AS: the RJ2a fast shock crossing a 2:1 face](/flume/shock-crossing-rj2a.png)

*RJ2a at $t = 0.2$, $N = 256$, today's operators, FNL (`make_doc_runs.sh <runs> rj2a-shock`): the density against the
exact solution and the plateau behind the fast shock with the uniform control; the same signature as on Sod, a bump
behind the face after c2f and a dip after f2c, which the uniform run does not have.*

**A3** (2 ranks; Sod: CPU and FNL agree to every printed digit; RJ2a: FNL). Plateau maxima of the refined run and, in
brackets, of its uniform control at the spacing of the side that holds the shock:

| Case | Shock (cells) | $\rho$ | $u_n$ | $p$ | $s$ |
|---|---|---|---|---|---|
| Sod c2f, N = 200 | 0.16 (0.17) | $2.5\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $4.5\cdot10^{-3}$ ($1.7\cdot10^{-5}$) | $6.9\cdot10^{-3}$ ($2.4\cdot10^{-5}$) | $1.8\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| Sod c2f, N = 400 | 0.16 (0.17) | $2.4\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $4.7\cdot10^{-3}$ ($1.2\cdot10^{-6}$) | $7.2\cdot10^{-3}$ ($1.5\cdot10^{-6}$) | $1.7\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| Sod f2c, N = 200 | 0.18 (0.17) | $1.6\cdot10^{-2}$ ($9.5\cdot10^{-4}$) | $6.7\cdot10^{-3}$ ($1.1\cdot10^{-4}$) | $1.0\cdot10^{-2}$ ($1.2\cdot10^{-4}$) | $1.2\cdot10^{-2}$ ($7.1\cdot10^{-4}$) |
| Sod f2c, N = 400 | 0.18 (0.17) | $1.5\cdot10^{-2}$ ($1.3\cdot10^{-4}$) | $6.6\cdot10^{-3}$ ($1.7\cdot10^{-5}$) | $1.0\cdot10^{-2}$ ($2.4\cdot10^{-5}$) | $1.1\cdot10^{-2}$ ($9.4\cdot10^{-5}$) |
| RJ2a c2f, N = 256 | 0.11 (0.11) | $1.4\cdot10^{-2}$ ($1.9\cdot10^{-4}$) | $8.1\cdot10^{-4}$ ($8.7\cdot10^{-5}$) | $6.2\cdot10^{-3}$ ($7.4\cdot10^{-5}$) | $5.9\cdot10^{-3}$ ($5.2\cdot10^{-5}$) |
| RJ2a f2c, N = 256 | 0.11 (0.11) | $1.9\cdot10^{-2}$ ($2.1\cdot10^{-4}$) | $9.1\cdot10^{-4}$ ($2.1\cdot10^{-4}$) | $6.5\cdot10^{-3}$ ($2.0\cdot10^{-4}$) | $7.6\cdot10^{-3}$ ($5.8\cdot10^{-5}$) |

With the reflux turned off (`make_amr_shock.py --no-reflux`, a measurement only): mass defect $1.2$ to
$2.2\cdot10^{-4}$ (Sod) and $0.9$ to $1.5\cdot10^{-4}$ (RJ2a), against at most $3\cdot10^{-10}$ with it; the shock off
by 0.2 to 0.35 cells; the plateau maxima at most 30 % larger than with the reflux.

### C3: the 6-point restriction at shocks

`lagrange6` is linear and unlimited, so it could overshoot where a shock lies in its footprint
([#68](https://github.com/szaghi/adam/issues/68) C3). The AS cases rerun under either restriction on the widened grid
(the `mean` run on it reproduces the AS baseline to every printed digit). FNL, 2 ranks:

| Case | Plateau $\rho$ | Plateau $p$ | Plateau $s$ | Band $\rho$ | Band $p$ | Shock (cells) |
|---|---|---|---|---|---|---|
| Sod c2f, mean / lagrange6 | $2.463$ / $2.457\cdot10^{-2}$ | $6.936$ / $6.911\cdot10^{-3}$ | $1.792$ / $1.787\cdot10^{-2}$ | $6.352$ / $6.330\cdot10^{-3}$ | $6.936$ / $6.911\cdot10^{-3}$ | 0.158 / 0.158 |
| Sod f2c, mean / lagrange6 | $1.576$ / $1.575\cdot10^{-2}$ | $1.018$ / $1.015\cdot10^{-2}$ | $1.180$ / $1.179\cdot10^{-2}$ | $3.481$ / $3.430\cdot10^{-3}$ | $3.718$ / $3.663\cdot10^{-3}$ | 0.177 / 0.177 |
| RJ2a c2f, mean / lagrange6 | $1.413$ / $1.410\cdot10^{-2}$ | $6.220$ / $6.218\cdot10^{-3}$ | $5.906$ / $5.895\cdot10^{-3}$ | $1.005$ / $1.025\cdot10^{-4}$ | $1.402$ / $1.300\cdot10^{-4}$ | 0.109 / 0.109 |

The PV-5 blast (the limited Balsara–Spicer blast across 2:1 seams, `nk` widened to 8) passes under both restrictions
and both schemes: $t = 0.01$ reached, no inadmissible stage, no non-finite flux, mass and seam flux mismatch unchanged;
the largest number of faces limited in a stage is 46336 against 48512 with the mean (split) and 17920 for both (HLLD).
The PV-5 runs take 2.4 to 2.5 times the wall time of the mean on FNL (4900 s against 1955 s, split; 4395 s against
1841 s, HLLD, with the same 433 steps; WSL, indicative only), against 16 % on the quadtree
[AO](./order-across-seams) case.

### C3b: the cost was the ghost kernel

Profiled with Nsight Systems (40 PV-5 steps, one rank), the local ghost update was 98 % of the FNL GPU time even with
the mean: 0.70 s per call. nvfortran parallelised its loop over the variables only (9 threads) and ran the map rows in
sequence (a "loop carried dependence" through the indirect writes): the schedule
[#27](https://github.com/szaghi/adam/issues/27) specifies the fix for, and that
[#53](https://github.com/szaghi/adam/issues/53) measured as the dominant FNL cost. With the two loops collapsed, as the
other map kernels are:

| Run | Before | After |
|---|---|---|
| local ghost kernel over 444 calls, `mean` | 310 s | 0.84 s |
| local ghost kernel over 444 calls, `lagrange6` | 815 s | 1.39 s |
| full PV-5 split run, `mean` | 1955 s | 276 s |
| full PV-5 split run, `lagrange6` | 4900 s | 288 s (premium 4 %, with 5 % more steps) |

The results are bitwise identical. WSL timings: #27 asks for the benchmark on a real node.

### C3c: the WENO restriction

| Check (FNL, 2 ranks) | mean | lagrange6 | weno6 |
|---|---|---|---|
| AO euler, reflux off, N = 256: L1 all / band, Linf | — | $8.55\cdot10^{-8}$ / $2.61\cdot10^{-7}$, $6.36\cdot10^{-6}$ | $8.43\cdot10^{-8}$ / $2.26\cdot10^{-7}$, $6.39\cdot10^{-6}$ |
| AO orders L1 / Linf (128 → 256) | — | 5.33 / 5.53 | 5.34 / 5.53 |
| Sod c2f plateau $\rho$, band $\rho$ | $2.463\cdot10^{-2}$, $6.352\cdot10^{-3}$ | $2.457\cdot10^{-2}$, $6.330\cdot10^{-3}$ | $2.464\cdot10^{-2}$, $6.365\cdot10^{-3}$ |
| Sod f2c plateau $\rho$, band $\rho$ | $1.576\cdot10^{-2}$, $3.481\cdot10^{-3}$ | $1.575\cdot10^{-2}$, $3.430\cdot10^{-3}$ | $1.576\cdot10^{-2}$, $3.516\cdot10^{-3}$ |
| RJ2a c2f plateau $\rho$, band $p$ | $1.413\cdot10^{-2}$, $1.402\cdot10^{-4}$ | $1.410\cdot10^{-2}$, $1.300\cdot10^{-4}$ | $1.414\cdot10^{-2}$, $1.360\cdot10^{-4}$ |
| PV-5 split / HLLD: most faces limited in a stage | 48512 / 17920 | 46336 / 17920 | 48064 / 17920 |
| PV-5 split wall (WSL, after C3b) | 276 s | 288 s | 278 s |

The AO rows are the legs `euler-l6` and `euler-w6` of [AO](./order-across-seams), which CPU and FNL reproduce to every
printed digit.

Related cases: [AO](./order-across-seams) is the smooth baseline of the same operators; the seam studies
[D1](../seams-studies/quintic-fill), [D2](../seams-studies/limited-fills), [E0](../seams-studies/capturing-error) and
[E2](../seams-studies/reflux-policy) rerun these crossings; the uniform references are [Sod](../euler/sod) and
[Ryu–Jones 2a](../mhd/ryu-jones-2a); the blast is the [positivity limiter](../mhd/positivity-limiter) case PV-5.

## Findings

- **The reflux does what it is for.** Mass is exact (to the floor the uniform RJ2a run shares), and the shock lands
  where the uniform run puts it. Without the reflux the mass drifts by $10^{-4}$ and the shock moves by a fraction of a
  cell.
- **A shock crossing a 2:1 face leaves an $O(1)$ error behind it.** The plateau carries a density and entropy defect of
  0.6 to 2.5 % of the shock jump and acoustic waves of 0.08 to 1 % in velocity and pressure: 4 to 4800 times the
  uniform level. This error **does not converge**: from $N = 200$ to 400 the maxima stay the same while the uniform
  level falls up to 14 times; only the width of the defect shrinks (the mean entropy error halves).
- **It is not the reflux.** Without the reflux the plateau defect changes by at most 30 %. A3 placed its source in the
  coarse-fine transfer of a discontinuity, the mean restriction and the prolongation of the ghosts, which Phases C and D
  replace; [E0](../seams-studies/capturing-error) later showed that no seam operator controls it: it is a capturing
  error.
- **The entropy defect leaves the face.** It moves with the flow, $u_2 t$ downstream of the face by $t = 0.2$, so the
  band entropy maxima are 100 to 1900 times below the plateau ones; the acoustic maxima of Sod c2f still sit at the
  face.
- **C3: the 6-point restriction does not worsen a shock crossing.** Every plateau and band defect is within 2 % of the
  mean's (the RJ2a band density 2 % above, the others equal or below), and the shock lands in the same place with the
  mass exact. The defect behind a crossing shock does not come from the restriction's order.
- **C3: it does not worsen positivity.** The limited blast passes as with the mean, with no more faces limited.
- **C3: on FNL it was slow on an octree,** 2.4 to 2.5 times the mean's wall time: a restricted ghost reads 216 fine
  cells, one after the other, in the loop that copies the other ghosts. The nonlinear `weno6` followed as an option, not
  as a fix.
- **C3b: the cost was the ghost kernel, for every regime.** Collapsing its two loops cut the full PV-5 run from 1955 s
  to 276 s (mean) and from 4900 s to 288 s (`lagrange6`), bitwise identical.
- **C3c: on the smooth vortex `weno6` equals `lagrange6`** (the band L1 13 % lower): no jump is detected, and the
  fourth-order rows at the donor block's edge do not show at these resolutions. The AO leg `euler-w6` asserts it.
- **C3c: at shocks it is as neutral as `lagrange6`:** shock positions identical, mass exact, every defect within 1.5 %
  of the mean's; PV-5 passes under both schemes with no inadmissible stage. Neither 6-point restriction makes a shock
  crossing worse, so the default stays `mean` until Phase D/E evidence; `weno6` is the option for data with jumps
  inside the restriction footprint that these cases do not produce.
- [Multirealm](../euler/multi-realm) leg 6 runs both regimes (bitwise unions), and [GP](../infrastructure/ghost-probe)
  probes `weno6` on `box3d`, `channel2d` and `refined3d`.
