# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 18:17 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/i157n/A` (commits: 80582712). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.1953 | 0.1% | 0.1948 | 0.000199 | 0 | 56,623,104 | 14,155,776 | 13.8 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.18 | 0.3% | 0.1201 | 0.05997 | 0 | 14,068,408 | 3,517,102 | 51.2 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.07985 | 3.8% (!) | 0.07966 | 0.000205 | 0 | 14,045,184 | 3,511,296 | 22.7 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 2 | fnl-nvf | 4 | 1.09 | 4.02 | 0.607 | 125.4% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.1953 | 13.8 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.18 | 51.2 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.07985 | 22.7 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2c-sphere-fixed: spread 3.8% above 2%
