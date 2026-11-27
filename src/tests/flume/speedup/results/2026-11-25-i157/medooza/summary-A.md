# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-25 15:48 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157/A` (commits: b8641a59). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2167 | 0.7% | 0.2163 | 0.000626 | 0.000747 | 4,128,768 | 1,032,192 | 210 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1547 | 0.1% | 0.1543 | 0.000582 | 0.001394 | 2,603,520 | 650,880 | 238 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2077 | 0.1% | 0.1594 | 0.04833 | 0.001317 | 2,602,783 | 650,696 | 319 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.7797 | 0.9% | 0.7789 | 0.000985 | 0 | 56,623,104 | 14,155,776 | 55.1 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.5054 | 0.5% | 0.3816 | 0.1237 | 0 | 14,068,408 | 3,517,102 | 144 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.382 | 0.1% | 0.3817 | 0.000105 | 0 | 14,045,184 | 3,511,296 | 109 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.228 | 0.5% | 0.2275 | 0.0004 | 0 | 6,291,456 | 1,572,864 | 145 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.1883 | 0.2% | 0.143 | 0.04484 | 0 | 3,538,944 | 884,736 | 213 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.1432 | 0.1% | 0.1429 | 5e-05 | 0 | 3,538,944 | 884,736 | 162 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.883 | 34.2% |
| 2 | fnl-nvf | 4 | 1.54 | 4.02 | 0.506 | 32.3% |
| 3 | fnl-nvf | 4 | 1.21 | 1.78 | 0.896 | 31.5% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.7797 | 55.1 | - | sphere 14.2 M/dev | - |
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2167 | 210 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1547 | 238 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2077 | 319 | - | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.5054 | 144 | - | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.382 | 109 | - | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.228 | 145 | - | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.1883 | 213 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.1432 | 162 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series
- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
