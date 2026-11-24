# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-24 11:55 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/gpu-weak-noib` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.2606 | 0.1% | 0.2601 | 0.000482 | 0 | 28,094,464 | 7,023,616 | 37.1 |
| 2 | wg-l4b38 | U | fnl-nvf | 32 | 0.2901 | 0.1% | 0.2895 | 0.001578 | 0 | 224,755,712 | 7,023,616 | 41.3 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.2606 | 37.1 | - | - | - |
| fnl-nvf | - | wg-l4b38 | 32 | 224,755,712 | 7,023,616 | 0.2901 | 41.3 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
