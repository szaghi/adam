# MHD rotor

| ID | Group | Milestone | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|---|
| MV-13 | Ideal MHD | M2 ([#41](https://github.com/szaghi/adam/issues/41)) · P6b | P | CPU, FNL | yes | `rotor` |

## Case and scope

A dense disc spins in a magnetised ambient gas (Balsara & Spicer 1999; the first rotor of Tóth 2000) and launches
torsional Alfvén waves that brake it, compresses the field and drives strong shocks: a stringent 2-D test of positivity
and of the multi-dimensional field evolution. There is no exact solution; the oracle checks:

- **symmetry:** the 180° rotation composed with $\mathbf{B} \to -\mathbf{B}$ (the ambient field is uniform) maps the
  solution to itself. A symmetry defect of the scheme is at truncation level ($10^{-3}$) from the first steps; a correct
  scheme leaves a round-off seed amplified by the flow;
- **positivity:** the floors are disabled, so a non-positive density or pressure stops the run.

There is no conservation check: the outflow boundaries do not conserve.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | ideal MHD, $\gamma = 1.4$ ($c_p = 3.5$, $c_v = 2.5$) | [models](../../models#ideal-mhd) |
| Space | WENO5 flux splitting on characteristic variables; `weno-riemann` under `--numerics` | `compute_residuals_weno`, `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Divergence control | GLM, $c_h = 6$ (`glm_ch_check = error`), `glm_alpha = 0.18`, `glm_damping_length = 1.0` | `mhd_glm_face_split_fluxes`, [`adam_flume_mhd_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_library.F90); [models](../../models#divergence-control-mixed-glm) |
| Floors | disabled | [models](../../models#positivity-floors) |
| Time | SSP Runge–Kutta, 5 stages (`runge-kutta-ssp-54`), CFL 0.4 | [`adam_stepper_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_stepper_object.F90) (SELAH) |
| Boundaries | zeroth-order extrapolation (outflow) | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [boundary conditions](../../boundary-conditions#extrapolation-zeroth-order-constant) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | the vortex grid of [V2](../euler/isentropic-vortex) (quadtree, $4 \times 4$ blocks), $256^2$ cells ($64^2$ per block); a power of two, so the cell centres are exact binary fractions and the initial state is bitwise symmetric (Tóth uses $200^2$) |
| Ranks | 2 |
| Initial condition | [`mhd-rotor`](../../initial-conditions#mhd-rotor): ambient $\rho = 1$ at rest, $p = 1$, $\mathbf{B} = (5/\sqrt{4\pi}, 0, 0)$; disc $\rho = 10$ for $r < 0.1$ around $(0.5, 0.5)$, rigid rotation at speed 2 at $r = 0.1$, linear taper to $r = 0.115$ |
| Boundaries | `extrapolation` on $x$ and $y$ |
| Final time | $t = 0.15$ (1152 steps) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/mhd/rotor/make_rotor.py` | writes the input from `vortex/vortex-n064.ini` (`--cells`, `--cfl`, `--time-max`, `--glm-ch`, `--it-save`) |
| `src/tests/flume/verification/mhd/rotor/check.sh` | runs the case (`--np N`, `--numerics SPEC`, `--legs mv13`) |
| `src/tests/flume/verification/mhd/orszag-tang/orszag_tang_oracle.py` | the symmetry check, `--field-parity even` (the $\mathbf{B} \to -\mathbf{B}$ composition) |
| `src/tests/flume/verification/make_doc_figures.py` `fig_rotor` | draws `rotor.png` from `<runs>/rotor` (`--only rotor`) |

```bash
cd src/tests/flume/verification/mhd
rotor/check.sh                                                     # CPU (about 10 min)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl rotor/check.sh    # FNL
rotor/check.sh --numerics hlld:primitive                           # weno-riemann, HLLD, primitive interpolation
```

`make_doc_runs.sh <runs> mhd` writes the figure run (`<runs>/rotor`: $256^2$); `make_doc_figures.py --runs <runs>` draws it.

## Results

![MHD rotor](/flume/rotor.png)

*Density and magnetic pressure $\lvert B\rvert^2/2$ of the rotor.*

| Check | Bound | CPU | FNL |
|---|---|---|---|
| symmetry (180° rotation with $\mathbf{B} \to -\mathbf{B}$), splitting, $t = 0.15$ | $10^{-8}$ | $6.5\cdot10^{-10}$ (2026-11-22; $1.8\cdot10^{-9}$ at M3-P3c; $10^{-11}$ at step 100) | $6.4\cdot10^{-10}$ |
| positivity (floors disabled) | run completes | pass | pass |
| `--numerics`: symmetry at step 100 | $10^{-8}$ | HLL $5.0\cdot10^{-13}$, HLLD primitive $3.6\cdot10^{-12}$ (without the sensor $3.0\cdot10^{-12}$), HLLD characteristic $3.8\cdot10^{-12}$ | $3.9\cdot10^{-13}$, $4.2\cdot10^{-12}$ ($4.3\cdot10^{-12}$), $3.2\cdot10^{-12}$ |
| `--numerics`: symmetry at the end, HLL | $5\cdot10^{-10}$ | $4.4\cdot10^{-11}$ | $1.2\cdot10^{-11}$ |
| `--numerics`: symmetry at the end, HLLD primitive | $7\cdot10^{-4}$ (one sensor flip, see below) | $6.08\cdot10^{-5}$ | $6.08\cdot10^{-5}$ |
| `--numerics`: symmetry at the end, HLLD primitive without the sensor | $4\cdot10^{-6}$ | $3.4\cdot10^{-7}$ | $1.8\cdot10^{-8}$ |
| `--numerics`: symmetry at the end, HLLD characteristic | $8\cdot10^{-9}$ | $7.4\cdot10^{-10}$ | $3.6\cdot10^{-10}$ |
| `--numerics`: symmetry at the end, any other spec | $10^{-4}$ (gross failures only) | not measured | not measured |

*Measured: ADAM `0bcdf870`, CPU gfortran 16.2.0 + Open MPI 5.0.11 (`iac_gnu`) on medooza, FNL nvfortran 26.9, CUDA
12.9 (`iac_nvf`) on NVIDIA A30, 2026-11-22 (`rotor/check.sh --numerics <spec>`, job 1369317): the `--numerics` rows at
the end. Repeated on ADAM `0bcdf870` + the night-2026-11-22 tree (the paper commit, tag pending), medooza (gnode22) and
A30, 2026-11-23 (`rotor/check.sh`, `rotor/check.sh --numerics <spec>` for the four specs, jobs 1369376 CPU, 1369378
FNL): every end value the same to the printed digits, plus the splitting row and the step-100 row. Repeated again on
ADAM `95ffd728` + the night-2026-11-22 tree (part 3), medooza and A30, 2026-11-23 (`rotor/check.sh --numerics <spec>`
for the four specs, job 1369637): the same to the printed digits. The positivity row is from M3-P3c (2026-10-01).*

$w$ and $B_z$ stay exactly 0 on the CPU.

Related cases: [Orszag–Tang](./orszag-tang) shares the symmetry oracle; [MHD on `weno-riemann`](./weno-riemann)
explains the two-point symmetry check; [multi-realm runs](../euler/multi-realm) split this rotor on a quadtree.

## Findings

- The final symmetry defect grows from round-off and orders by the dissipation of the solver and by the backend's
  round-off: the least dissipative interpolations amplify the seed most, and FNL contracts differently into FMAs. The
  exception is a flip of a switch, below.
- A single end-of-run tolerance cannot separate a round-off seed from a defect across solvers, so under `--numerics`
  the symmetry is asserted twice: at step 100, before the flow amplifies the seed, within $10^{-8}$, and at the end
  within a bound per spec (`rotor/check.sh`): 10 times the larger of the CPU and FNL defects measured on 2026-11-22,
  rounded up to one digit; a spec never measured keeps $10^{-4}$, which catches only a truncation-level defect
  ($10^{-3}$).
- **HLLD primitive with the sensor flips a switch.** Its end defect is $6.08\cdot10^{-5}$ on both backends, against
  $2.2\cdot10^{-9}$ before `a2258f41` ("FLUME CPU integrates with SELAH", [#132](https://github.com/szaghi/adam/issues/132)).
  The commit changed only the order of the Runge–Kutta sums: SELAH adds $\Delta t\sum_k a_{sk}R_k$ to $q$ once, the
  order of the former FNL kernels, where `rk_object` added each term to $q$ in turn. With the former order restored in
  SELAH's two kernels (a scratch build), the run is bitwise that of `a8f3fbd2`, $2.22\cdot10^{-9}$ at the end. The
  defect is a jump: the time series reads $6\cdot10^{-10}$ at step 910 and $5\cdot10^{-7}$ at step 920. At step
  912, stage 5, the WENO sensor of the flux correction ($\min_k w_k/d_k \ge \tau = 0.2$) reads
  $\tau - 7.2\cdot10^{-10}$ at the $x$ face $(0.371, 0.693)$ and $\tau + 3.4\cdot10^{-10}$ at its image: the round-off
  seed, $10^{-9}$ by then, is larger than the distance to $\tau$, so the correction is dropped on one side only. Both
  summation orders are a lottery on the seed: with the CFL moved by $10^{-8}$ (relative) the end defect is
  $3.4\cdot10^{-7}$ or $2.5\cdot10^{-4}$ with SELAH's order, $3.1\cdot10^{-7}$ or $1.9\cdot10^{-4}$ with the former one.
  CPU and FNL give the same value because they now sum in the same order and meet the same face. The end bound of
  this spec is therefore the rule's $7\cdot10^{-4}$, which admits a flip; the step-100 check, before the seed reaches
  $\tau$, is the one that catches a defect of the scheme. The other specs with the sensor flipped no face in the runs
  measured, and could on another seed. A sensor without a jump (a ramp of the correction weight around $\tau$) would
  remove the lottery; it changes every case with the sensor and is not done here (medooza CPU, gfortran 16.2,
  2026-11-23, jobs 1369623, 1369628).
- The end defects moved since M3-P3c (the SSP(5,4) tableau fix #104 and SELAH #132 changed the round-off seed): HLLD
  characteristic $1.1\cdot10^{-7} \to 7.4\cdot10^{-10}$, splitting $1.8\cdot10^{-9} \to 6.5\cdot10^{-10}$, HLLD primitive
  without the sensor $7.4\cdot10^{-8} \to 3.4\cdot10^{-7}$ (CPU).
