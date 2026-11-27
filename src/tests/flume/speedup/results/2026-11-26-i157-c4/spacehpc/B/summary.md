# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 18:17 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/i157n/B` (commits: i157n-s0). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.1874 | 0.2% | 0.187 | 0.000324 | 0 | 56,623,104 | 14,155,776 | 13.2 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.1718 | 1.0% | 0.111 | 0.06127 | 0 | 14,068,408 | 3,517,102 | 48.8 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.08333 | 10.5% (!) | 0.08314 | 0.000299 | 0 | 14,045,184 | 3,511,296 | 23.7 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 2 | fnl-nvf | 4 | 1.09 | 4.02 | 0.558 | 106.1% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.1874 | 13.2 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.1718 | 48.8 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.08333 | 23.7 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2c-sphere-fixed: spread 10.5% above 2%
