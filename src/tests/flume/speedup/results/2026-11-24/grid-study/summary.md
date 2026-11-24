# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-24 08:25 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/grid-study` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | wg-b16 | U | fnl-nvf | 1 | 0.1324 | 1.0% | 0.1322 | 2e-06 | 0 | 2,097,152 | 2,097,152 | 63.1 |
| 2 | wg-b16 | U | fnl-nvf | 2 | 0.07212 | 0.4% | 0.07198 | 5e-05 | 0 | 2,097,152 | 1,048,576 | 68.8 |
| 2 | wg-b16 | U | fnl-nvf | 4 | 0.0469 | 0.1% | 0.04675 | 8.1e-05 | 0 | 2,097,152 | 524,288 | 89.5 |
| 2 | wg-b20 | U | fnl-nvf | 1 | 0.223 | 0.5% | 0.2227 | 2e-06 | 0 | 4,096,000 | 4,096,000 | 54.4 |
| 2 | wg-b20 | U | fnl-nvf | 2 | 0.117 | 0.7% | 0.1168 | 2.6e-05 | 0 | 4,096,000 | 2,048,000 | 57.1 |
| 2 | wg-b20 | U | fnl-nvf | 4 | 0.07295 | 0.4% | 0.07282 | 6.7e-05 | 0 | 4,096,000 | 1,024,000 | 71.2 |
| 2 | wg-b24 | U | cpu | 64 | 1.366 | 5.3% (!) | 1.346 | 0.06357 | 0 | 7,077,888 | 110,592 | 1.24e+04 |
| 2 | wg-b24 | U | fnl-nvf | 1 | 0.3611 | 0.3% | 0.3606 | 2e-06 | 0 | 7,077,888 | 7,077,888 | 51 |
| 2 | wg-b24 | U | fnl-nvf | 2 | 0.1876 | 0.5% | 0.1873 | 0.000317 | 0 | 7,077,888 | 3,538,944 | 53 |
| 2 | wg-b24 | U | fnl-nvf | 4 | 0.1147 | 0.2% | 0.1145 | 4.4e-05 | 0 | 7,077,888 | 1,769,472 | 64.8 |
| 2 | wg-b28 | U | fnl-nvf | 1 | 0.5361 | 0.2% | 0.5353 | 2e-06 | 0 | 11,239,424 | 11,239,424 | 47.7 |
| 2 | wg-b28 | U | fnl-nvf | 2 | 0.2764 | 0.2% | 0.2759 | 8.1e-05 | 0 | 11,239,424 | 5,619,712 | 49.2 |
| 2 | wg-b28 | U | fnl-nvf | 4 | 0.1636 | 0.0% | 0.1634 | 0.000147 | 0 | 11,239,424 | 2,809,856 | 58.2 |
| 2 | wg-b32 | U | fnl-nvf | 1 | 0.7751 | 0.1% | 0.7741 | 2e-06 | 0 | 16,777,216 | 16,777,216 | 46.2 |
| 2 | wg-b32 | U | fnl-nvf | 2 | 0.3939 | 0.2% | 0.3933 | 0.000257 | 0 | 16,777,216 | 8,388,608 | 47 |
| 2 | wg-b32 | U | fnl-nvf | 4 | 0.2335 | 0.1% | 0.2332 | 0.000122 | 0 | 16,777,216 | 4,194,304 | 55.7 |
| 2 | wg-b36 | U | fnl-nvf | 1 | 1.061 | 0.6% | 1.06 | 3e-06 | 0 | 23,887,872 | 23,887,872 | 44.4 |
| 2 | wg-b36 | U | fnl-nvf | 2 | 0.5417 | 0.2% | 0.541 | 0.000488 | 0 | 23,887,872 | 11,943,936 | 45.4 |
| 2 | wg-b36 | U | fnl-nvf | 4 | 0.3128 | 0.2% | 0.3123 | 0.000169 | 0 | 23,887,872 | 5,971,968 | 52.4 |
| 2 | wg-b38 | U | cpu | 64 | 4.458 | 0.7% | 4.447 | 0.08341 | 0 | 28,094,464 | 438,976 | 1.02e+04 |
| 2 | wg-b38 | U | fnl-nvf | 1 | 1.234 | 0.1% | 1.233 | 4e-06 | 0 | 28,094,464 | 28,094,464 | 43.9 |
| 2 | wg-b38 | U | fnl-nvf | 2 | 0.6281 | 0.3% | 0.6272 | 0.00043 | 0 | 28,094,464 | 14,047,232 | 44.7 |
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.355 | 0.1% | 0.3544 | 0.000186 | 0 | 28,094,464 | 7,023,616 | 50.5 |
| 2 | wg-l2b48 | U | fnl-nvf | 1 | 0.3111 | 0.2% | 0.3106 | 1e-06 | 0 | 7,077,888 | 7,077,888 | 44 |
| 2 | wg-l4b12 | U | fnl-nvf | 1 | 0.4673 | 0.1% | 0.4668 | 2e-06 | 0 | 7,077,888 | 7,077,888 | 66 |
| 2 | wg-l4b38 | U | fnl-nvf | 32 | 0.4703 | 0.3% | 0.4697 | 0.001623 | 0 | 224,755,712 | 7,023,616 | 67 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | - | wg-b24 | 64 | 7,077,888 | 110,592 | 1.366 | 1.24e+04 | - | - | - |
| cpu | - | wg-b38 | 64 | 28,094,464 | 438,976 | 4.458 | 1.02e+04 | - | - | - |
| fnl-nvf | - | wg-b16 | 1 | 2,097,152 | 2,097,152 | 0.1324 | 63.1 | 1 (1) | - | - |
| fnl-nvf | - | wg-b16 | 2 | 2,097,152 | 1,048,576 | 0.07212 | 68.8 | 0.918 (1) | - | - |
| fnl-nvf | - | wg-b16 | 4 | 2,097,152 | 524,288 | 0.0469 | 89.5 | 0.706 (1) | - | - |
| fnl-nvf | - | wg-b20 | 1 | 4,096,000 | 4,096,000 | 0.223 | 54.4 | 1 (1) | - | - |
| fnl-nvf | - | wg-b20 | 2 | 4,096,000 | 2,048,000 | 0.117 | 57.1 | 0.953 (1) | - | - |
| fnl-nvf | - | wg-b20 | 4 | 4,096,000 | 1,024,000 | 0.07295 | 71.2 | 0.764 (1) | - | - |
| fnl-nvf | - | wg-b24 | 1 | 7,077,888 | 7,077,888 | 0.3611 | 51 | 1 (1) | - | - |
| fnl-nvf | - | wg-b24 | 2 | 7,077,888 | 3,538,944 | 0.1876 | 53 | 0.962 (1) | - | - |
| fnl-nvf | - | wg-b24 | 4 | 7,077,888 | 1,769,472 | 0.1147 | 64.8 | 0.787 (1) | - | - |
| fnl-nvf | - | wg-b28 | 1 | 11,239,424 | 11,239,424 | 0.5361 | 47.7 | 1 (1) | - | - |
| fnl-nvf | - | wg-b28 | 2 | 11,239,424 | 5,619,712 | 0.2764 | 49.2 | 0.97 (1) | - | - |
| fnl-nvf | - | wg-b28 | 4 | 11,239,424 | 2,809,856 | 0.1636 | 58.2 | 0.819 (1) | - | - |
| fnl-nvf | - | wg-b32 | 1 | 16,777,216 | 16,777,216 | 0.7751 | 46.2 | 1 (1) | - | - |
| fnl-nvf | - | wg-b32 | 2 | 16,777,216 | 8,388,608 | 0.3939 | 47 | 0.984 (1) | - | - |
| fnl-nvf | - | wg-b32 | 4 | 16,777,216 | 4,194,304 | 0.2335 | 55.7 | 0.83 (1) | - | - |
| fnl-nvf | - | wg-b36 | 1 | 23,887,872 | 23,887,872 | 1.061 | 44.4 | 1 (1) | - | - |
| fnl-nvf | - | wg-b36 | 2 | 23,887,872 | 11,943,936 | 0.5417 | 45.4 | 0.979 (1) | - | - |
| fnl-nvf | - | wg-b36 | 4 | 23,887,872 | 5,971,968 | 0.3128 | 52.4 | 0.848 (1) | - | - |
| fnl-nvf | - | wg-b38 | 1 | 28,094,464 | 28,094,464 | 1.234 | 43.9 | 1 (1) | - | - |
| fnl-nvf | - | wg-b38 | 2 | 28,094,464 | 14,047,232 | 0.6281 | 44.7 | 0.983 (1) | - | - |
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.355 | 50.5 | 0.869 (1) | - | - |
| fnl-nvf | - | wg-l2b48 | 1 | 7,077,888 | 7,077,888 | 0.3111 | 44 | - | - | - |
| fnl-nvf | - | wg-l4b12 | 1 | 7,077,888 | 7,077,888 | 0.4673 | 66 | - | - | - |
| fnl-nvf | - | wg-l4b38 | 32 | 224,755,712 | 7,023,616 | 0.4703 | 67 | - | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
| wg-b24 | cpu np 64 | fnl-nvf | 1 | 3.78 |
| wg-b24 | cpu np 64 | fnl-nvf | 2 | 7.28 |
| wg-b24 | cpu np 64 | fnl-nvf | 4 | 11.9 |
| wg-b38 | cpu np 64 | fnl-nvf | 1 | 3.61 |
| wg-b38 | cpu np 64 | fnl-nvf | 2 | 7.1 |
| wg-b38 | cpu np 64 | fnl-nvf | 4 | 12.6 |

## Warnings

- fnl-nvf np 1 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b20: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b20: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b20: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 wg-b24: built from sources that differ from the recorded commit (src_modified=yes); spread 5.3% above 2%
- fnl-nvf np 1 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b28: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b28: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b28: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b32: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b32: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b32: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b36: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b36: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b36: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-l2b48: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-l4b12: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
