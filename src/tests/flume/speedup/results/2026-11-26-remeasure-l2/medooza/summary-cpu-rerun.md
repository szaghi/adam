# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 11:18 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/l2-cpu2` (commits: unknown). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | cpu | 64 | 7.816 | 5.4% (!) | 7.794 | 0.1961 | 0 | 56,623,104 | 884,736 | 8.83e+03 |
| 2 | 2a-sphere-uniform | U | cpu | 512 | 1.169 | 12.6% (!) | 1.114 | 0.249 | 0 | 56,623,104 | 110,592 | 1.06e+04 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | 2a-sphere-uniform | 64 | 56,623,104 | 884,736 | 7.816 | 8.83e+03 | 1 (64) | - | - |
| cpu | - | 2a-sphere-uniform | 512 | 56,623,104 | 110,592 | 1.169 | 1.06e+04 | 0.836 (64) | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
| 2a-sphere-uniform | cpu np 64 | cpu | 512 | 6.68 |

## Warnings

- cpu np 64 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 5.4% above 2%
- cpu np 512 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 12.6% above 2%
