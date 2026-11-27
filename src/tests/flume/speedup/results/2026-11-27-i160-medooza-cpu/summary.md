# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-27 08:58 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i160/work` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | cpu | 64 | 1.223 | 4.6% (!) | 1.219 | 0.004113 | 0.00073 | 4,128,768 | 64,512 | 1.9e+04 |
| 3 | 3a-channel-uniform | U | cpu | 64 | 3.59 | 6.7% (!) | 3.536 | 0.1288 | 0 | 6,291,456 | 98,304 | 3.65e+04 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | 1a-step-uniform | 64 | 4,128,768 | 64,512 | 1.223 | 1.9e+04 | - | - | - |
| cpu | - | 3a-channel-uniform | 64 | 6,291,456 | 98,304 | 3.59 | 3.65e+04 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- cpu np 64 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 4.6% above 2%
- cpu np 64 3a-channel-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 6.7% above 2%
