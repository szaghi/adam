# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-25 13:09 UTC from `/lustre/projects/1026/CALCOLI/adam/speedup-133/grid-study` (commits: unknown). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | wg-b16 | U | fnl-nvf | 1 | 0.05236 | 0.1% | 0.05225 | 4e-06 | 0 | 2,097,152 | 2,097,152 | 25 |
| 2 | wg-b16 | U | fnl-nvf | 2 | 0.03249 | 0.2% | 0.03239 | 5.1e-05 | 0 | 2,097,152 | 1,048,576 | 31 |
| 2 | wg-b16 | U | fnl-nvf | 4 | 0.02193 | 0.9% | 0.02187 | 8.6e-05 | 0 | 2,097,152 | 524,288 | 41.8 |
| 2 | wg-b20 | U | fnl-nvf | 1 | 0.08441 | 0.0% | 0.08423 | 6e-06 | 0 | 4,096,000 | 4,096,000 | 20.6 |
| 2 | wg-b24 | U | fnl-nvf | 1 | 0.131 | 0.2% | 0.1307 | 6e-06 | 0 | 7,077,888 | 7,077,888 | 18.5 |
| 2 | wg-b24 | U | fnl-nvf | 2 | 0.07204 | 0.2% | 0.07186 | 3.3e-05 | 0 | 7,077,888 | 3,538,944 | 20.4 |
| 2 | wg-b24 | U | fnl-nvf | 4 | 0.04221 | 0.1% | 0.0421 | 8.6e-05 | 0 | 7,077,888 | 1,769,472 | 23.9 |
| 2 | wg-b28 | U | fnl-nvf | 1 | 0.1896 | 0.2% | 0.1892 | 7e-06 | 0 | 11,239,424 | 11,239,424 | 16.9 |
| 2 | wg-b32 | U | fnl-nvf | 1 | 0.2572 | 0.2% | 0.2568 | 7e-06 | 0 | 16,777,216 | 16,777,216 | 15.3 |
| 2 | wg-b36 | U | fnl-nvf | 1 | 0.3784 | 0.2% | 0.3778 | 9e-06 | 0 | 23,887,872 | 23,887,872 | 15.8 |
| 2 | wg-b38 | U | fnl-nvf | 1 | 0.4395 | 0.2% | 0.4389 | 8e-06 | 0 | 28,094,464 | 28,094,464 | 15.6 |
| 2 | wg-b38 | U | fnl-nvf | 2 | 0.2246 | 0.1% | 0.2241 | 0.000305 | 0 | 28,094,464 | 14,047,232 | 16 |
| 2 | wg-b38 | U | fnl-nvf | 4 | 0.1189 | 0.5% | 0.1185 | 0.000149 | 0 | 28,094,464 | 7,023,616 | 16.9 |
| 2 | wg-b48 | U | fnl-nvf | 1 | 0.8462 | 0.0% | 0.845 | 9e-06 | 0 | 56,623,104 | 56,623,104 | 14.9 |
| 2 | wg-b48 | U | fnl-nvf | 2 | 0.4238 | 0.2% | 0.4231 | 0.000436 | 0 | 56,623,104 | 28,311,552 | 15 |
| 2 | wg-b48 | U | fnl-nvf | 4 | 0.2184 | 0.1% | 0.218 | 0.000258 | 0 | 56,623,104 | 14,155,776 | 15.4 |
| 2 | wg-b56 | U | fnl-nvf | 1 | 1.305 | 0.1% | 1.303 | 1e-05 | 0 | 89,915,392 | 89,915,392 | 14.5 |
| 2 | wg-b64 | U | fnl-nvf | 1 | 1.908 | 0.0% | 1.906 | 1e-05 | 0 | 134,217,728 | 134,217,728 | 14.2 |
| 2 | wg-b64 | U | fnl-nvf | 2 | 0.9651 | 0.0% | 0.9637 | 0.000548 | 0 | 134,217,728 | 67,108,864 | 14.4 |
| 2 | wg-b64 | U | fnl-nvf | 4 | 0.4869 | 0.1% | 0.486 | 0.001079 | 0 | 134,217,728 | 33,554,432 | 14.5 |
| 2 | wg-l2b48 | U | fnl-nvf | 1 | 0.1162 | 0.1% | 0.116 | 6e-06 | 0 | 7,077,888 | 7,077,888 | 16.4 |
| 2 | wg-l4b12 | U | fnl-nvf | 1 | 0.1741 | 0.2% | 0.1738 | 7e-06 | 0 | 7,077,888 | 7,077,888 | 24.6 |
| 2 | wg-l4b38 | U | fnl-nvf | 4 | 0.8613 | 0.0% | 0.8601 | 0.000268 | 0 | 224,755,712 | 56,188,928 | 15.3 |
| 2 | wg-l4b38 | U | fnl-nvf | 8 | 0.4508 | 0.1% | 0.45 | 0.000544 | 0 | 224,755,712 | 28,094,464 | 16 |
| 2 | wg-l4b38 | U | fnl-nvf | 16 | 0.2387 | 0.2% | 0.2383 | 0.000805 | 0 | 224,755,712 | 14,047,232 | 17 |
| 2 | wg-l4b38 | U | fnl-nvf | 32 | 0.1303 | 0.1% | 0.13 | 0.000855 | 0 | 224,755,712 | 7,023,616 | 18.6 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| fnl-nvf | - | wg-b16 | 1 | 2,097,152 | 2,097,152 | 0.05236 | 25 | 1 (1) | - | - |
| fnl-nvf | - | wg-b16 | 2 | 2,097,152 | 1,048,576 | 0.03249 | 31 | 0.806 (1) | - | - |
| fnl-nvf | - | wg-b16 | 4 | 2,097,152 | 524,288 | 0.02193 | 41.8 | 0.597 (1) | - | - |
| fnl-nvf | - | wg-b20 | 1 | 4,096,000 | 4,096,000 | 0.08441 | 20.6 | - | - | - |
| fnl-nvf | - | wg-b24 | 1 | 7,077,888 | 7,077,888 | 0.131 | 18.5 | 1 (1) | - | - |
| fnl-nvf | - | wg-b24 | 2 | 7,077,888 | 3,538,944 | 0.07204 | 20.4 | 0.909 (1) | - | - |
| fnl-nvf | - | wg-b24 | 4 | 7,077,888 | 1,769,472 | 0.04221 | 23.9 | 0.776 (1) | - | - |
| fnl-nvf | - | wg-b28 | 1 | 11,239,424 | 11,239,424 | 0.1896 | 16.9 | - | - | - |
| fnl-nvf | - | wg-b32 | 1 | 16,777,216 | 16,777,216 | 0.2572 | 15.3 | - | - | - |
| fnl-nvf | - | wg-b36 | 1 | 23,887,872 | 23,887,872 | 0.3784 | 15.8 | - | - | - |
| fnl-nvf | - | wg-b38 | 1 | 28,094,464 | 28,094,464 | 0.4395 | 15.6 | 1 (1) | - | - |
| fnl-nvf | - | wg-b38 | 2 | 28,094,464 | 14,047,232 | 0.2246 | 16 | 0.979 (1) | - | - |
| fnl-nvf | - | wg-b38 | 4 | 28,094,464 | 7,023,616 | 0.1189 | 16.9 | 0.925 (1) | - | - |
| fnl-nvf | - | wg-b48 | 1 | 56,623,104 | 56,623,104 | 0.8462 | 14.9 | 1 (1) | - | - |
| fnl-nvf | - | wg-b48 | 2 | 56,623,104 | 28,311,552 | 0.4238 | 15 | 0.998 (1) | - | - |
| fnl-nvf | - | wg-b48 | 4 | 56,623,104 | 14,155,776 | 0.2184 | 15.4 | 0.969 (1) | - | - |
| fnl-nvf | - | wg-b56 | 1 | 89,915,392 | 89,915,392 | 1.305 | 14.5 | - | - | - |
| fnl-nvf | - | wg-b64 | 1 | 134,217,728 | 134,217,728 | 1.908 | 14.2 | 1 (1) | - | - |
| fnl-nvf | - | wg-b64 | 2 | 134,217,728 | 67,108,864 | 0.9651 | 14.4 | 0.989 (1) | - | - |
| fnl-nvf | - | wg-b64 | 4 | 134,217,728 | 33,554,432 | 0.4869 | 14.5 | 0.98 (1) | - | - |
| fnl-nvf | - | wg-l2b48 | 1 | 7,077,888 | 7,077,888 | 0.1162 | 16.4 | - | - | - |
| fnl-nvf | - | wg-l4b12 | 1 | 7,077,888 | 7,077,888 | 0.1741 | 24.6 | - | - | - |
| fnl-nvf | - | wg-l4b38 | 4 | 224,755,712 | 56,188,928 | 0.8613 | 15.3 | 1 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 8 | 224,755,712 | 28,094,464 | 0.4508 | 16 | 0.955 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 16 | 224,755,712 | 14,047,232 | 0.2387 | 17 | 0.902 (4) | - | - |
| fnl-nvf | - | wg-l4b38 | 32 | 224,755,712 | 7,023,616 | 0.1303 | 18.6 | 0.826 (4) | - | - |

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|

## Warnings

- fnl-nvf np 1 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b16: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b20: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b24: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b28: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b32: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b36: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b48: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b48: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b48: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b56: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-b64: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 wg-b64: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-b64: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-l2b48: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 wg-l4b12: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 wg-l4b38: built from sources that differ from the recorded commit (src_modified=yes)
