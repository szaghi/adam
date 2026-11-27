# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 22:11 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157c/B` (commits: i157n-s0). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2156 | 0.2% | 0.2154 | 0.000622 | 0.000974 | 4,128,768 | 1,032,192 | 209 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1552 | 0.1% | 0.1548 | 0.000585 | 0.001674 | 2,603,520 | 650,880 | 239 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2078 | 1.0% | 0.1598 | 0.04798 | 0.001603 | 2,602,783 | 650,696 | 319 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.5544 | 0.4% | 0.5535 | 0.001997 | 0 | 56,623,104 | 14,155,776 | 39.2 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.406 | 0.0% | 0.2639 | 0.142 | 0 | 14,068,408 | 3,517,102 | 115 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.1926 | 0.1% | 0.1923 | 0.000154 | 0 | 14,045,184 | 3,511,296 | 54.9 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.2219 | 0.8% | 0.2213 | 0.000478 | 0 | 6,291,456 | 1,572,864 | 141 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.1786 | 0.4% | 0.1372 | 0.04115 | 0 | 3,538,944 | 884,736 | 202 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.1376 | 0.1% | 0.1373 | 4.5e-05 | 0 | 3,538,944 | 884,736 | 156 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.876 | 33.8% |
| 2 | fnl-nvf | 4 | 1.37 | 4.02 | 0.714 | 110.8% |
| 3 | fnl-nvf | 4 | 1.24 | 1.78 | 0.907 | 29.8% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.5544 | 39.2 | - | sphere 14.2 M/dev | - |
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2156 | 209 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1552 | 239 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2078 | 319 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.406 | 115 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.1926 | 54.9 | - | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.2219 | 141 | - | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.1786 | 202 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.1376 | 156 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series
- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
