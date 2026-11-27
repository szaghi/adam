# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-27 06:05 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/i157h/B` (commits: i157n-s0). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.08629 | 1.7% | 0.08616 | 0.000643 | 0.000817 | 4,128,768 | 1,032,192 | 83.6 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.06701 | 0.4% | 0.06692 | 0.000475 | 0.0014 | 2,603,520 | 650,880 | 103 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.08845 | 2.3% (!) | 0.0685 | 0.02062 | 0.001326 | 2,602,783 | 650,696 | 136 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.06087 | 0.7% | 0.06056 | 8.1e-05 | 0 | 6,291,456 | 1,572,864 | 38.7 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.07126 | 0.5% | 0.04301 | 0.02796 | 0 | 3,538,944 | 884,736 | 80.5 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.04367 | 0.3% | 0.04344 | 4.8e-05 | 0 | 3,538,944 | 884,736 | 49.4 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 0.976 | 1.59 | 0.812 | 32.0% |
| 3 | fnl-nvf | 4 | 0.854 | 1.78 | 0.784 | 63.2% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.08629 | 83.6 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.06701 | 103 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.08845 | 136 | - | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.06087 | 38.7 | - | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.07126 | 80.5 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.04367 | 49.4 | - | - | - |

Series:

- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 1c-step-runtime: spread 2.3% above 2%
