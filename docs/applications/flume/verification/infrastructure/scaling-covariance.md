# Scaling covariance and the reference layer

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| NV-0, NV-4, NV-5, NV-6, NV-8, NV-9 | Infrastructure | — ([#49](https://github.com/szaghi/adam/issues/49), non-dimensional form, N0–N2) | CPU, FNL | no | — |

## Case and scope

Ideal Euler and MHD in FLUME's units carry no dimensionless number, so an input rescaled by powers of two (lengths
$\times 2^j$, velocities $\times 2^k$, density $\times 4^m$) is the same problem in other units, and a scheme without
absolute constants reproduces it **bit for bit** once the outputs are scaled back. The scales are exact powers of two,
so any difference is a constant of the scheme that does not scale.

What the legs catch:
- **absolute constants in the scheme** (NV-0, NV-4): the absolute WENO regulariser $\varepsilon$ (`zeps = 1e-6`), an
  absolute positivity floor, a guarded immersed-boundary cut spacing. Each rescaling moves one quantity at a time
  (length only, velocity only, density only, all together), so a mismatch is attributed to the quantity whose
  magnitude moved;
- **a wrong reference layer** (NV-5): a dimensional input with a [`[reference]`](../../input#reference-optional-dimensional-input)
  section must run as its hand-normalised twin; an option the layer misclassifies or leaves unconverted shows as a
  mismatch;
- **restart files that lose the units** (NV-8) and **dimensional output that is not the base output times its scale**
  (NV-9);
- **a scheme that stops limiting at physical magnitudes** (NV-6): Sod in SI interstellar units, where the classic
  weights collapse to the linear ones.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Models | Euler and ideal MHD | [models](../../models#units-and-scaling) |
| WENO weights | `js` (classic, absolute $\varepsilon$, the default) and `si` (scale-invariant, `[weno] weights = si`, N1) | `weno_descaling`, `weno_weights_exponent`, [`adam_weno_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_weno_object.F90); [numerics](../../numerics#weno-reconstruction) |
| Positivity floor | relative only (N1b) | `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Reference layer | classifies every option, converts the dimensional ones, logs the conversions, records the references with the restart files, writes dimensional output on request | `classify`, `convert_options`, `check_restart`, `save_restart`, `variable_output`, [`adam_flume_reference_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_reference_object.F90) |
| Independent classification | the oracle's own table of every option (it refuses an unknown one) | `src/tests/flume/verification/scaling/scaling.py` |

## Initial and boundary conditions

Defaults of every leg: $j = 2$, $k = -1$, $m = -2$ (`--j`, `--k`, `--m`), 2 ranks. The inputs are those of existing
cases, unchanged except for the units.

| Leg | Script | Cases | Rescalings or references |
|---|---|---|---|
| NV-0 / NV-4 | `check.sh` | `sod-x` (Euler, fast path), the isentropic vortex (`vortex-periodic`, periodic Euler), RJ2a (`rj2a-x`, MHD without cleaning), the limited Balsara–Spicer blast (`blast-limiter`, EGLM, splitting scheme); `orszag-tang` reads no initial-condition parameter and is listed as not rescalable | length only $(j, 0, 0)$, velocity only $(0, k, 0)$, density only $(0, 0, m)$, all $(j, k, m)$; NV-0 prints the baseline and exits 0, `--expect-bitwise` (NV-4) fails on any non-bitwise comparison |
| NV-5 | `reference.sh` | every single-realm regression case (the forest manifests and `orszag-tang` are refused by the layer); `--inputs` for the generated inputs | $L_0 = 2^j$, $u_0 = 2^k$, $\rho_0 = 4^m$ |
| NV-8 | `reference.sh --restart` | `sod-x` (fast path), `shock-cylinder-ib` (immersed boundary, AMR), `amr-periodic-reflux` (staged path, reflux), `blast-amr-limiter` (MHD EGLM, AMR, positivity limiter) | as NV-5; 20 steps against 10 steps plus a restart to 20 |
| NV-9 | `reference.sh --output` | `sod-x`, `shock-cylinder-ib`, `field-loop` (GLM), `blast-amr-limiter` (EGLM, AMR), `uniform-amr-mhd` (3-D, AMR), each with residual and auxiliary fields and a 64-point slice along $x$ through the domain centre | as NV-5, `output_units = dimensional` |
| NV-6 | `interstellar.sh` | the V1 Sod (`sod-x`) on `weno` and on `weno-riemann` HLLC (6th-order correction, `weno` sensor) | $L_0 = 1$ pc, $u_0 = 10$ km/s, $\rho_0 = 10^{-21}$ kg/m³ (not powers of two) |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/scaling/scaling.py` | `rescale`, `dimensionalize`, `physical`, `compare` (`--conservative`, `--output`), `check-log` |
| `src/tests/flume/verification/scaling/check.sh` | NV-0 / NV-4 (`--weights js\|si`, `--expect-bitwise`, `--legs LIST`: `sod-x`, `vortex`, `rj2a`, `blast`, `orszag-tang`, default all) |
| `src/tests/flume/verification/scaling/reference.sh` | NV-5, NV-8 (`--restart`), NV-9 (`--output`); `--legs LIST` names regression cases, default every case for NV-5 and the cases of the table for NV-8 and NV-9; `--inputs` takes no leg selection |
| `src/tests/flume/verification/scaling/interstellar.sh`, `interstellar_oracle.py` | NV-6 (`--legs LIST`: `nv6-weno`, `nv6-weno-riemann`, default both) |
| `src/tests/flume/verification/scaling/weights-exe.sh` | runs any verification script with `[weno] weights = si` (`FLUME_EXE` pointed at it, `WEIGHTS_EXE` at the executable) |
| `src/tests/flume/regression/*/input.ini` | the inputs of NV-0, NV-5, NV-8, NV-9 |

```bash
cd src/tests/flume/verification/scaling
./check.sh                                   # NV-0 baseline, js weights
./check.sh --weights si --expect-bitwise     # NV-4
./reference.sh                               # NV-5
./reference.sh --restart ; ./reference.sh --output   # NV-8, NV-9
./interstellar.sh                            # NV-6
FLUME_EXE=$PWD/weights-exe.sh ../vortex/check.sh     # any oracle with the si weights
```

## Results

No figure: every leg is a bitwise identity or a contract (#85 classes NV-5 / NV-8 as setup table only).

| Weights | Bound | Result (CPU and FNL) |
|---|---|---|
| `js` (default) | baseline (NV-0) | length-only bitwise; every velocity or density rescaling differs by $10^{-4}$ to $10^{-1}$ (the absolute $\varepsilon$) |
| `si` (`--weights si --expect-bitwise`) | bitwise (NV-4) | 16/16 bitwise; and the blast at density $\times 4^{-24}$, where the limiter floor acts, bitwise |

The `si` row needs the relative-only positivity floor and the unguarded immersed-boundary cut spacing (#49, N1b). Run
through `weights-exe.sh`, the orders and the shock errors of the `si` weights are judged by the same oracles and bounds
as the default ones. One tolerance differs: RJ2a's GLM-against-no-cleaning and EGLM-against-GLM legs are bitwise with
`js` but within $10^{-11}$ with `si` (measured $\le 2.5\cdot10^{-13}$), the round-off of a uniform field reconstructed
with face-varying weights ([numerics](../../numerics#weno-reconstruction)).

**Reference layer (NV-5).** The references are powers of two, so the solver sees the base numbers exactly and the run
must equal the base run bit for bit: every conservative field of the last checkpoint, the history and slice files byte
for byte. `scaling.py check-log` compares every conversion the layer logs with the base value (exact equality, every
dimensional option converted), which also covers the options whose conversion does not show in the solution.

| Leg | Cases | Bound | Result (CPU and FNL) |
|---|---|---|---|
| regression cases | the 21 single-realm ones at the time of N2 (`orszag-tang` and the 5 multi-realm ones refused by the layer) | bitwise, conversions exact and complete | 21/21 bitwise, conversions exact and complete |
| generated inputs | the MHD linear waves (fast GLM, Alfvén EGLM, slow), CPAW, the magnetised vortex, the div(B) peak, the EGLM pulse, Shu–Osher, the isentropic vortex, Sod with a slice, Sod with two gradient AMR markers | bitwise, conversions exact and complete | 11/11 bitwise, conversions exact and complete |

**Restart (NV-8).** On the four cases a run of 20 steps equals a run of 10 steps plus a restart to 20 bit for bit
(conservative fields, residual and conservation histories). Restarting with another density reference, or without the
`.reference` record, is refused.

**Dimensional output (NV-9).** The twin computes the base numbers and writes them times the references: every field,
the block spacing, the XDMF times, every history, the slices and the `.units` record equal the base output times its
exact power of two, bit for bit (`scaling.py compare --output`).

| Exception | Tolerance | Measured |
|---|---|---|
| temperature: the base computes $p/(\rho R)$ with the regression inputs' $c_p$, $c_v$, the twin $(p/\rho)/R$ with $R = 1$ in the solver | $2\cdot10^{-15}$ relative | $2.3\cdot10^{-16}$ on `sod-x` |
| underflow: the FNL build flushes subnormal results to zero (`-fast`), so a written 0 matches an expected subnormal | — | one `rv` of $3\cdot10^{-308}$ on the `shock-cylinder-ib` slice, scaled by $2^{-5}$ |

**Interstellar units (NV-6).** With the classic weights the absolute `zeps = 1e-6` dominates smoothness indicators of
order $10^{-42}$: the weights collapse to the linear ones and the scheme stops limiting, as #48 (A3) predicted from the
formula.

| Run | Bound | `weno` | `weno-riemann` HLLC |
|---|---|---|---|
| code units, `js`: total variation of $\rho/\rho_0$ | — | 0.877 | 0.884 |
| raw SI, `js` (negative control) | TV $> 1.1\times$ code units, or a stop on a non-finite state | TV 1.162 | TV 2.455 ($L_1\times20$, an 8% undershoot); survives on the CPU, stops at step 90 on FNL |
| raw SI, `si`, against code units `si` | $10^{-12}$ in $\rho/\rho_0$ | $\le 1.1\cdot10^{-13}$ | $\le 1.1\cdot10^{-13}$ |
| SI with `[reference]`, `js`, against code units `js` | $10^{-12}$ in $\rho/\rho_0$ | $\le 1.1\cdot10^{-13}$ | $\le 1.1\cdot10^{-13}$ |

Related cases: [Sod](../euler/sod), [Isentropic vortex](../euler/isentropic-vortex) and the MHD cases provide the
inputs; [Ranks without blocks](./ranks-without-blocks) is the other identity leg of this group.

## Findings

- The classic weights break covariance through the absolute $\varepsilon$ only: length-only rescalings are bitwise
  with `js`, and the `si` weights with the relative floor make every rescaling bitwise.
- `weno-riemann` reads the same weights in its flux-correction sensor, which goes blind with them: in raw SI it leaves
  the unlimited correction on and undershoots towards vacuum; round-off decides whether the run survives (CPU) or stops
  (FNL, step 90).
- Both cures restore the code-units solution within round-off: the scale-invariant weights on the raw SI input, and
  the reference layer with the classic weights.
- The first NV-5 run found two options the scaling tool had misclassified since N0: `rho_amplitude` taken as
  dimensionless (it is a density) and `loop_amplitude` as a vector potential (it is the field of the loop). The two
  independent tables disagreed, and the field-loop case was not bitwise.
- `reference.sh` runs every single-realm case of the suite by default; the suite has grown since N2 (39 cases, 34 of
  them single-realm), so the 21/21 count is the measurement of that time, not of the current list.
- A self-gravity case runs on the FFT twin of the executable under test (`<exe>_fft`, or `FLUME_EXE_FFT`), as in the
  regression harness, so one run covers the whole suite. The polytrope's `star_radius`, `star_x0`, `star_y0`,
  `star_z0` (lengths), `star_density` (a density) and `star_pressure_factor` (a pure number) are classified in both
  tables since 2026-10: before, the two star cases stopped on an unclassified option, loudly, and NV-5 did not cover
  them.
