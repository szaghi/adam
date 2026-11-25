# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-25 15:48 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/B` (commits: 10c69ae2). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2164 | 0.5% | 0.2163 | 0.000619 | 0.000755 | 4,128,768 | 1,032,192 | 210 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1553 | 0.3% | 0.1549 | 0.000579 | 0.001469 | 2,603,520 | 650,880 | 239 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2083 | 0.4% | 0.16 | 0.04835 | 0.001332 | 2,602,783 | 650,696 | 320 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.7777 | 0.4% | 0.7768 | 0.001227 | 0 | 56,623,104 | 14,155,776 | 54.9 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.4503 | 0.3% | 0.3277 | 0.1224 | 0 | 14,068,408 | 3,517,102 | 128 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.3278 | 0.1% | 0.3275 | 0.000103 | 0 | 14,045,184 | 3,511,296 | 93.4 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.2279 | 0.4% | 0.2274 | 0.000388 | 0 | 6,291,456 | 1,572,864 | 145 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.1853 | 0.7% | 0.1406 | 0.04446 | 0 | 3,538,944 | 884,736 | 209 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.141 | 0.2% | 0.1407 | 8e-05 | 0 | 3,538,944 | 884,736 | 159 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.878 | 34.1% |
| 2 | fnl-nvf | 4 | 1.73 | 4.02 | 0.589 | 37.4% |
| 3 | fnl-nvf | 4 | 1.23 | 1.78 | 0.909 | 31.4% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.7777 | 54.9 | - | sphere 14.2 M/dev | - |
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2164 | 210 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1553 | 239 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2083 | 320 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.4503 | 128 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.3278 | 93.4 | - | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.2279 | 145 | - | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.1853 | 209 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.141 | 159 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series
- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
