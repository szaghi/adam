# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-26 17:46 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/l2` (commits: unknown). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.1952 | 0.2% | 0.1947 | 0.000207 | 0 | 56,623,104 | 14,155,776 | 13.8 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 8 | 0.09889 | 0.2% | 0.09853 | 0.000252 | 0 | 56,623,104 | 7,077,888 | 14 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 16 | 0.05163 | 0.9% | 0.05109 | 0.000603 | 0 | 56,623,104 | 3,538,944 | 14.6 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 32 | 0.02662 | 8.1% (!) | 0.02641 | 0.000389 | 0 | 56,623,104 | 1,769,472 | 15 |
| 2 | w2-sphere-l5 | U | fnl-nvf | 32 | 0.2002 | 6.5% (!) | 0.1992 | 0.00191 | 0 | 452,984,832 | 14,155,776 | 14.1 |
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.09282 | 0.1% | 0.09252 | 0.000221 | 0 | 28,094,464 | 7,023,616 | 13.2 |
| 2 | wg-l4b38 | U | fnl-nvf | 4 | 0.727 | 0.1% | 0.7258 | 0.000481 | 0 | 224,755,712 | 56,188,928 | 12.9 |
| 2 | wg-l4b38 | U | fnl-nvf | 8 | 0.3683 | 0.3% | 0.3676 | 0.000842 | 0 | 224,755,712 | 28,094,464 | 13.1 |
| 2 | wg-l4b38 | U | fnl-nvf | 16 | 0.1857 | 1.1% | 0.1842 | 0.00136 | 0 | 224,755,712 | 14,047,232 | 13.2 |
| 2 | wg-l4b38 | U | fnl-nvf | 32 | 0.09586 | 1.5% | 0.09502 | 0.000737 | 0 | 224,755,712 | 7,023,616 | 13.6 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.1952 | 13.8 | 1 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 8 | 56,623,104 | 7,077,888 | 0.09889 | 14 | 0.987 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 16 | 56,623,104 | 3,538,944 | 0.05163 | 14.6 | 0.945 (4) | - | - |
| fnl-nvf | - | 2a-sphere-uniform | 32 | 56,623,104 | 1,769,472 | 0.02662 | 15 | 0.916 (4) | - | - |
| fnl-nvf | - | w2-sphere-l5 | 32 | 452,984,832 | 14,155,776 | 0.2002 | 14.1 | - | - | - |
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.09282 | 13.2 | - | - | - |
| fnl-nvf | - | wg-l4b38 | 4 | 224,755,712 | 56,188,928 | 0.727 | 12.9 | 1 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 8 | 224,755,712 | 28,094,464 | 0.3683 | 13.1 | 0.987 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 16 | 224,755,712 | 14,047,232 | 0.1857 | 13.2 | 0.979 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 32 | 224,755,712 | 7,023,616 | 0.09586 | 13.6 | 0.948 (4) | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 4 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 8.1% above 2%
- fnl-nvf np 32 w2-sphere-l5: built from sources that differ from the recorded commit (src_modified=yes); spread 6.5% above 2%
- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
