# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 23:31 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157d/C` (commits: i157n-s0). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2165 | 0.4% | 0.2161 | 0.000624 | 0.000898 | 4,128,768 | 1,032,192 | 210 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1552 | 0.5% | 0.155 | 0.000578 | 0.001666 | 2,603,520 | 650,880 | 238 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2073 | 0.2% | 0.1603 | 0.04742 | 0.001578 | 2,602,783 | 650,696 | 319 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.88 | 33.5% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2165 | 210 | - | step 1.03 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1552 | 238 | - | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2073 | 319 | - | - | - |

Series:

- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
