# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 22:11 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157c/A` (commits: 80582712). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2161 | 0.5% | 0.2159 | 0.00062 | 0.000751 | 4,128,768 | 1,032,192 | 209 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1552 | 0.5% | 0.1548 | 0.000571 | 0.001448 | 2,603,520 | 650,880 | 238 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2072 | 0.9% | 0.1596 | 0.04772 | 0.001353 | 2,602,783 | 650,696 | 318 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.5743 | 0.6% | 0.5734 | 0.002336 | 0 | 56,623,104 | 14,155,776 | 40.6 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.4219 | 0.7% | 0.2808 | 0.1408 | 0 | 14,068,408 | 3,517,102 | 120 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.1986 | 0.1% | 0.1983 | 0.000143 | 0 | 14,045,184 | 3,511,296 | 56.5 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.2286 | 1.1% | 0.2281 | 0.000447 | 0 | 6,291,456 | 1,572,864 | 145 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.183 | 0.2% | 0.1403 | 0.04243 | 0 | 3,538,944 | 884,736 | 207 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.141 | 0.1% | 0.1406 | 5.2e-05 | 0 | 3,538,944 | 884,736 | 159 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.878 | 33.5% |
| 2 | fnl-nvf | 4 | 1.36 | 4.02 | 0.717 | 112.5% |
| 3 | fnl-nvf | 4 | 1.25 | 1.78 | 0.912 | 29.8% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.5743 | 40.6 | - | sphere 14.2 M/dev | - |
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2161 | 209 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1552 | 238 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2072 | 318 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.4219 | 120 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.1986 | 56.5 | - | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.2286 | 145 | - | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.183 | 207 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.141 | 159 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series
- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
