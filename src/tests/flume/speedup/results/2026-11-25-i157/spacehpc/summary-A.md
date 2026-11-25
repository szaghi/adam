# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-25 16:09 UTC from `/tmp/claude-1000/-home-stefano-fortran-adam/427178af-df51-4839-bbb7-e1d56d892c3d/scratchpad/h100/A` (commits: b8641a59). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.2536 | 0.3% | 0.2529 | 0.000366 | 0 | 56,623,104 | 14,155,776 | 17.9 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.2177 | 0.6% | 0.1587 | 0.05895 | 0 | 14,068,408 | 3,517,102 | 61.9 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.1364 | 2.0% (!) | 0.1361 | 0.000125 | 0 | 14,045,184 | 3,511,296 | 38.8 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 2 | fnl-nvf | 4 | 1.17 | 4.02 | 0.461 | 59.6% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.2536 | 17.9 | - | sphere 14.2 M/dev | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.2177 | 61.9 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.1364 | 38.8 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2c-sphere-fixed: spread 2.0% above 2%
