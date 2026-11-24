# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-24 18:14 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/ftz-windows` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | cpu | 64 | 8.917 | 9.1% (!) | 8.907 | 0.1781 | 0 | 56,623,104 | 884,736 | 1.01e+04 |
| 2 | 2a-sphere-uniform | U | cpu | 512 | 1.224 | 0.5% | 1.222 | 0.04437 | 0 | 56,623,104 | 110,592 | 1.11e+04 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | 2a-sphere-uniform | 64 | 56,623,104 | 884,736 | 8.917 | 1.01e+04 | 1 (64) | - | - |
| cpu | - | 2a-sphere-uniform | 512 | 56,623,104 | 110,592 | 1.224 | 1.11e+04 | 0.91 (64) | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
| 2a-sphere-uniform | cpu np 64 | cpu | 512 | 7.28 |

## Warnings

- cpu np 64 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 9.1% above 2%
- cpu np 512 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
