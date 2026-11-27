# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-27 10:43 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/cpu` (commits: 07a414ba). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | cpu | 96 | 0.5779 | 4.4% (!) | 0.5768 | 0.001033 | 0.000311 | 4,128,768 | 43,008 | 1.34e+04 |
| 2 | 2a-sphere-uniform | U | cpu | 96 | 3.902 | 4.5% (!) | 3.895 | 0.4148 | 0 | 56,623,104 | 589,824 | 6.61e+03 |
| 3 | 3a-channel-uniform | U | cpu | 96 | 1.564 | 6.0% (!) | 1.563 | 0.0534 | 0 | 6,291,456 | 65,536 | 2.39e+04 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | 1a-step-uniform | 96 | 4,128,768 | 43,008 | 0.5779 | 1.34e+04 | - | - | - |
| cpu | - | 2a-sphere-uniform | 96 | 56,623,104 | 589,824 | 3.902 | 6.61e+03 | - | - | - |
| cpu | - | 3a-channel-uniform | 96 | 6,291,456 | 65,536 | 1.564 | 2.39e+04 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- cpu np 96 1a-step-uniform: spread 4.4% above 2%
- cpu np 96 2a-sphere-uniform: spread 4.5% above 2%
- cpu np 96 3a-channel-uniform: spread 6.0% above 2%
