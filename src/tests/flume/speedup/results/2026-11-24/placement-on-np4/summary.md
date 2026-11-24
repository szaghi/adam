# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-24 07:43 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/bind/on-np4` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.7808 | 0.1% | 0.7799 | 0.002146 | 0 | 56,623,104 | 14,155,776 | 55.2 |
| 2 | wg-b24 | U | fnl-nvf | 4 | 0.1137 | 0.9% | 0.1135 | 5.6e-05 | 0 | 7,077,888 | 1,769,472 | 64.3 |
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.3552 | 0.1% | 0.3546 | 0.000349 | 0 | 28,094,464 | 7,023,616 | 50.6 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.7808 | 55.2 | - | - | - |
| fnl-nvf | - | wg-b24 | 4 | 7,077,888 | 1,769,472 | 0.1137 | 64.3 | - | - | - |
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.3552 | 50.6 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
