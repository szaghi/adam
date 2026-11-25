# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-25 13:22 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/matrix` (commits: unknown). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.2541 | 0.1% | 0.2537 | 0.000452 | 0 | 56,623,104 | 14,155,776 | 18 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 8 | 0.1371 | 0.7% | 0.1367 | 0.00025 | 0 | 56,623,104 | 7,077,888 | 19.4 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 16 | 0.07803 | 0.9% | 0.07783 | 0.000446 | 0 | 56,623,104 | 3,538,944 | 22 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 32 | 0.04624 | 0.2% | 0.0461 | 0.00044 | 0 | 56,623,104 | 1,769,472 | 26.1 |
| 2 | w2-sphere-l5 | U | fnl-nvf | 32 | 0.2776 | 0.1% | 0.2771 | 0.001294 | 0 | 452,984,832 | 14,155,776 | 19.6 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.2541 | 18 | 1 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 8 | 56,623,104 | 7,077,888 | 0.1371 | 19.4 | 0.927 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 16 | 56,623,104 | 3,538,944 | 0.07803 | 22 | 0.814 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 32 | 56,623,104 | 1,769,472 | 0.04624 | 26.1 | 0.687 (4) | - | - |
| fnl-nvf | - | w2-sphere-l5 | 32 | 452,984,832 | 14,155,776 | 0.2776 | 19.6 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 w2-sphere-l5: built from sources that differ from the recorded commit (src_modified=yes)
