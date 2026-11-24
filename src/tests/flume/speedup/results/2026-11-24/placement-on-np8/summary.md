# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-24 07:57 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/bind/on-np8` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 8 | 0.4908 | 0.0% | 0.4902 | 0.00058 | 0 | 56,623,104 | 7,077,888 | 69.3 |
| 2 | wg-b24 | U | fnl-nvf | 8 | 0.08738 | 0.5% | 0.08725 | 6.6e-05 | 0 | 7,077,888 | 884,736 | 98.8 |
| 2 | wg-b38 | U | fnl-nvf | 8 | 0.2356 | 0.1% | 0.2353 | 0.000101 | 0 | 28,094,464 | 3,511,808 | 67.1 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 8 | 56,623,104 | 7,077,888 | 0.4908 | 69.3 | - | - | - |
| fnl-nvf | - | wg-b24 | 8 | 7,077,888 | 884,736 | 0.08738 | 98.8 | - | - | - |
| fnl-nvf | - | wg-b38 | 8 | 28,094,464 | 3,511,808 | 0.2356 | 67.1 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 8 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
