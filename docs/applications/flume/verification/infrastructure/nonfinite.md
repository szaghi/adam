# Non-finite state guard

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| — (`nonfinite`) | Infrastructure | — ([#45](https://github.com/szaghi/adam/issues/45)) | CPU, FNL | no | — |

## Case and scope

A positive control of the non-finite guard. Before [#45](https://github.com/szaghi/adam/issues/45) a run whose state
went NaN finished normally and exited 0: the CPU `sod-amr` event of #45 was caught only by the regression digest. Every
step now counts the non-finite (NaN or infinite) values of the committed state over the interior cells, sums the count
over the ranks and, on a non-zero total, stops the run, every rank naming its first non-finite value (block Morton
code, cell, cell centre, variable).

The case plants a NaN in the right state of the Sod tube of [V1](../euler/sod) (`[initial_conditions_region_2]`
pressure `p = nan`) and requires the run to fail **on the guard's message**.

What the case catches:
- **a guard that never fires:** it would pass every other test unnoticed, since those tests run finite states only;
- **a run that stops for another reason:** a non-zero exit is not enough, the log must hold the guard's message;
- **a test folded away by the floating-point mode:** the count tests the bit pattern (the 11 exponent bits of a
  binary64 all set), not `x /= x`, so no compiler mode (`-fast` included) can remove it.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Base case | the V1 Sod tube: compressible Euler, WENO5 characteristic flux splitting, SSP-33, CFL 0.5 | [Sod shock tube](../euler/sod) |
| Count, CPU | per-model kernel over the interior cells, exponent field all ones, OpenMP reduction | `check_nonfinite`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `count_nonfinite`, [`adam_flume_cpu_aux_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_aux_kernels_agnostic.INC) |
| Count, FNL | the same count on the device; the state is copied to the host only to locate a failure | `check_nonfinite`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); `count_nonfinite_dev`, [`adam_flume_fnl_aux_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_aux_kernels_agnostic.INC) |
| Reduction and stop | sum over the ranks; every rank reports its first non-finite value and stops through `error_stop` | `nonfinite_total`, `stop_nonfinite`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Call site | first action of every step's post-processing, both backends | `post_step_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90), [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Input | `sod/sod-x.ini` with the line `p = 0.1` of `[initial_conditions_region_2]` replaced by `p = nan` (written by `check.sh` to `sod-nan.ini`) |
| Domain, grid | as [V1](../euler/sod): $[0, 1]$ along $x$, 64 blocks of $50 \times 4 \times 4$ cells, `ngc = 3` |
| Ranks | 2 |
| Left state, $x \le 0.5$ | $(\rho, u, p) = (1, 0, 1)$ |
| Right state, $x > 0.5$ | $(\rho, u, p) = (0.125, 0, \text{NaN})$ |
| Initial condition | [`riemann-problem`](../../initial-conditions#riemann-problem-piecewise-constant-axis-aligned-regions) |
| Boundaries | `extrapolation` on all six faces |
| Expected | non-zero exit and the log line `non-finite (NaN or infinite) values in the state at step` |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/sod/sod-x.ini` | the base input |
| `src/tests/flume/verification/nonfinite/check.sh` | plants the NaN, runs the case, requires the failure on the guard (`--np N`, `--legs nonfinite`, its one leg; `FLUME_EXE` selects the backend) |

```bash
cd src/tests/flume/verification/nonfinite && ./check.sh              # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh          # FNL
```

## Results

A contract leg: no figure. Values read from the kept logs `nonfinite/work-adam_flume_{cpu,fnl}-np2/log.txt`:

| Check | Expected | CPU | FNL |
|---|---|---|---|
| exit status of the run | non-zero | non-zero (`MPI_ABORT`) | non-zero (`MPI_ABORT`) |
| guard message in the log | present | present | present |
| step of the stop | — | 1 | 1 |
| non-finite values counted (all ranks) | > 0 | 81408 | 81408 |
| first value, rank 0 | — | block code 9, cell (45, 1, 1), centre $(0.4725, 0.03125, 0.03125)$, variable `r` | the same |
| first value, rank 1 | — | block code 41, cell (45, 1, 1), centre $(0.4725, 0.03125, 0.53125)$, variable `r` | the same |

Related cases: [Sod shock tube](../euler/sod) is the base input; [ranks without blocks](./ranks-without-blocks) is the
other infrastructure leg that runs a contract rather than a flow.

## Findings

- The guard fires at the first step on both backends with the same count and the same first site: one step of WENO5
  carries the NaN of the right state to $x = 0.4725$, left of the diaphragm, into the density.
- Before the guard (commit `4f7824b0`) the same input ran to the end and exited 0, per that commit's record.
- The guard costs one pass over the interior per step and leaves the numerics untouched; on FNL the state crosses to
  the host only when a failure has to be located.
- The check deletes the checkpoints after the run; it keeps only the log and the histories.
