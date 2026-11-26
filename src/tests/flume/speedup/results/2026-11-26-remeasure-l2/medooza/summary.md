# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 10:00 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/l2` (commits: unknown). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | cpu | 64 | 8.004 | 7.0% (!) | 7.994 | 0.3951 | 0 | 56,623,104 | 884,736 | 9.05e+03 |
| 2 | 2a-sphere-uniform | U | cpu | 512 | 1.181 | 6.6% (!) | 1.124 | 0.2528 | 0 | 56,623,104 | 110,592 | 1.07e+04 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.5745 | 1.0% | 0.5736 | 0.001243 | 0 | 56,623,104 | 14,155,776 | 40.6 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 8 | 0.3111 | 0.2% | 0.3105 | 0.000649 | 0 | 56,623,104 | 7,077,888 | 43.9 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 16 | 0.1604 | 3.2% (!) | 0.1552 | 0.005282 | 0 | 56,623,104 | 3,538,944 | 45.3 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 32 | 0.08524 | 1.1% | 0.08442 | 0.001627 | 0 | 56,623,104 | 1,769,472 | 48.2 |
| 2 | w2-sphere-l5 | U | fnl-nvf | 32 | 0.6044 | 3.0% (!) | 0.6013 | 0.006608 | 0 | 452,984,832 | 14,155,776 | 42.7 |
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.2698 | 0.1% | 0.2692 | 0.000295 | 0 | 28,094,464 | 7,023,616 | 38.4 |
| 2 | wg-l4b38 | U | fnl-nvf | 32 | 0.2909 | 0.1% | 0.2899 | 0.0047 | 0 | 224,755,712 | 7,023,616 | 41.4 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | 2a-sphere-uniform | 64 | 56,623,104 | 884,736 | 8.004 | 9.05e+03 | 1 (64) | - | - |
| cpu | - | 2a-sphere-uniform | 512 | 56,623,104 | 110,592 | 1.181 | 1.07e+04 | 0.847 (64) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.5745 | 40.6 | 1 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 8 | 56,623,104 | 7,077,888 | 0.3111 | 43.9 | 0.923 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 16 | 56,623,104 | 3,538,944 | 0.1604 | 45.3 | 0.896 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 32 | 56,623,104 | 1,769,472 | 0.08524 | 48.2 | 0.842 (4) | - | - |
| fnl-nvf | - | w2-sphere-l5 | 32 | 452,984,832 | 14,155,776 | 0.6044 | 42.7 | - | - | - |
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.2698 | 38.4 | - | - | - |
| fnl-nvf | - | wg-l4b38 | 32 | 224,755,712 | 7,023,616 | 0.2909 | 41.4 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
| 2a-sphere-uniform | cpu np 64 | cpu | 512 | 6.78 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 4 | 13.9 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 8 | 25.7 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 16 | 49.9 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 32 | 93.9 |

## Warnings

- cpu np 64 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 7.0% above 2%
- cpu np 512 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 6.6% above 2%
- fnl-nvf np 4 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 3.2% above 2%
- fnl-nvf np 32 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 w2-sphere-l5: built from sources that differ from the recorded commit (src_modified=yes); spread 3.0% above 2%
- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
