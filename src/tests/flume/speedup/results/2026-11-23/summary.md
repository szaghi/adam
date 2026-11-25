# FLUME speedup matrix: summary

Written by `analyze.py` on 2026-11-23 19:12 UTC from `/ifs/hpc/home/szaghi/adam-tests/speedup-133/work` (commits: HEAD). Protocol and ratios: `README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).

## Step costs

| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | Seam fills (s/step) | N cells | Cells/device | ns/cell-step |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1a-step-uniform | U | cpu | 64 | 1.26 | 2.9% (!) | 1.255 | 0.004649 | 0.000733 | 4,128,768 | 64,512 | 1.95e+04 |
| 1 | 1a-step-uniform | U | cpu | 128 | 0.6415 | 2.3% (!) | 0.6371 | 0.005093 | 0.000117 | 4,128,768 | 32,256 | 1.99e+04 |
| 1 | 1a-step-uniform | U | cpu | 256 | 0.3275 | 2.8% (!) | 0.3227 | 0.004971 | 0.000178 | 4,128,768 | 16,128 | 2.03e+04 |
| 1 | 1a-step-uniform | U | fnl-nvf | 4 | 0.2162 | 0.7% | 0.216 | 0.000626 | 0.000758 | 4,128,768 | 1,032,192 | 209 |
| 1 | 1a-step-uniform | U | fnl-nvf | 8 | 0.1188 | 0.3% | 0.1186 | 0.000474 | 0.00071 | 4,128,768 | 516,096 | 230 |
| 1 | 1a-step-uniform | U | fnl-nvf | 16 | 0.07043 | 0.5% | 0.07011 | 0.000351 | 0.0006535 | 4,128,768 | 258,048 | 273 |
| 1 | 1b-step-fixed | F | cpu | 64 | 0.8865 | 2.8% (!) | 0.8818 | 0.005701 | 0.006582 | 2,603,520 | 40,680 | 2.18e+04 |
| 1 | 1b-step-fixed | F | cpu | 128 | 0.5443 | 0.7% | 0.5421 | 0.005115 | 0.006434 | 2,603,520 | 20,340 | 2.68e+04 |
| 1 | 1b-step-fixed | F | cpu | 256 | 0.3339 | 1.6% | 0.3287 | 0.004571 | 0.006029 | 2,603,520 | 10,170 | 3.28e+04 |
| 1 | 1b-step-fixed | F | fnl-nvf | 4 | 0.1541 | 0.0% | 0.1539 | 0.000509 | 0.001386 | 2,603,520 | 650,880 | 237 |
| 1 | 1b-step-fixed | F | fnl-nvf | 8 | 0.08898 | 0.1% | 0.08878 | 0.000517 | 0.001313 | 2,603,520 | 325,440 | 273 |
| 1 | 1b-step-fixed | F | fnl-nvf | 16 | 0.05725 | 1.5% | 0.05701 | 0.0004965 | 0.00107 | 2,603,520 | 162,720 | 352 |
| 1 | 1c-step-runtime | R | cpu | 64 | 1.187 | 3.2% (!) | 1.175 | 0.01344 | 0.009504 | 2,602,783 | 40,668 | 2.92e+04 |
| 1 | 1c-step-runtime | R | cpu | 128 | 1.182 | 3.8% (!) | 1.17 | 0.01362 | 0.009689 | 2,602,783 | 20,334 | 5.81e+04 |
| 1 | 1c-step-runtime | R | cpu | 256 | 1.175 | 1.3% | 1.163 | 0.0153 | 0.009412 | 2,602,783 | 10,167 | 1.16e+05 |
| 1 | 1c-step-runtime | R | fnl-nvf | 4 | 0.2073 | 0.8% | 0.1593 | 0.04825 | 0.001342 | 2,602,783 | 650,696 | 319 |
| 1 | 1c-step-runtime | R | fnl-nvf | 8 | 0.1194 | 0.5% | 0.0926 | 0.02685 | 0.00141 | 2,602,783 | 325,348 | 367 |
| 1 | 1c-step-runtime | R | fnl-nvf | 16 | 0.07565 | 1.9% | 0.06057 | 0.01542 | 0.001212 | 2,602,783 | 162,674 | 465 |
| 1 | w1-step-m22 | U | fnl-nvf | 8 | 0.2065 | 0.6% | 0.206 | 0.000502 | 0.000765 | 7,805,952 | 975,744 | 212 |
| 1 | w1-step-m32 | U | fnl-nvf | 16 | 0.2211 | 0.1% | 0.2207 | 0.000344 | 0.000787 | 16,515,072 | 1,032,192 | 214 |
| 1 | w1-step-m8 | U | fnl-nvf | 1 | 0.209 | 0.1% | 0.2088 | 3e-06 | 0.000219 | 1,032,192 | 1,032,192 | 202 |
| 2 | 2a-sphere-uniform | U | cpu | 64 | 12.68 | 7.5% (!) | 12.67 | 1.438 | 0 | 56,623,104 | 884,736 | 1.43e+04 |
| 2 | 2a-sphere-uniform | U | cpu | 128 | 6.481 | 0.5% | 6.471 | 0.7354 | 0 | 56,623,104 | 442,368 | 1.47e+04 |
| 2 | 2a-sphere-uniform | U | cpu | 256 | 3.244 | 0.2% | 3.242 | 0.3643 | 0 | 56,623,104 | 221,184 | 1.47e+04 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 4 | 0.7799 | 0.1% | 0.779 | 0.001484 | 0 | 56,623,104 | 14,155,776 | 55.1 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 8 | 0.4911 | 0.1% | 0.4905 | 0.000603 | 0 | 56,623,104 | 7,077,888 | 69.4 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 16 | 0.2599 | 0.1% | 0.2595 | 0.000821 | 0 | 56,623,104 | 3,538,944 | 73.4 |
| 2 | 2a-sphere-uniform | U | fnl-nvf | 32 | 0.1674 | 0.4% | 0.1671 | 0.000751 | 0 | 56,623,104 | 1,769,472 | 94.6 |
| 2 | 2b-sphere-runtime | R | cpu | 64 | 5.759 | 6.3% (!) | 5.696 | 0.07472 | 0 | 14,068,408 | 219,819 | 2.62e+04 |
| 2 | 2b-sphere-runtime | R | cpu | 128 | 3.945 | 0.5% | 3.892 | 0.05786 | 0 | 14,068,408 | 109,909 | 3.59e+04 |
| 2 | 2b-sphere-runtime | R | cpu | 256 | 3.064 | 1.2% | 3.021 | 0.04808 | 0 | 14,068,408 | 54,955 | 5.58e+04 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 4 | 0.526 | 0.8% | 0.3819 | 0.1438 | 0 | 14,068,408 | 3,517,102 | 150 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 8 | 0.3499 | 0.8% | 0.2792 | 0.07107 | 0 | 14,068,408 | 1,758,551 | 199 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 16 | 0.3118 | 2.0% | 0.2648 | 0.04682 | 0 | 14,068,408 | 879,276 | 355 |
| 2 | 2b-sphere-runtime | R | fnl-nvf | 32 | 0.1875 | 0.6% | 0.1575 | 0.03112 | 0 | 14,068,408 | 439,638 | 427 |
| 2 | 2c-sphere-fixed | F | cpu | 64 | 4.638 | 1.1% | 4.633 | 0.003907 | 0 | 14,045,184 | 219,456 | 2.11e+04 |
| 2 | 2c-sphere-fixed | F | cpu | 128 | 3.13 | 1.0% | 3.126 | 0.004708 | 0 | 14,045,184 | 109,728 | 2.85e+04 |
| 2 | 2c-sphere-fixed | F | cpu | 256 | 1.968 | 0.2% | 1.965 | 0.004421 | 0 | 14,045,184 | 54,864 | 3.59e+04 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 4 | 0.3814 | 0.1% | 0.3811 | 0.000147 | 0 | 14,045,184 | 3,511,296 | 109 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 8 | 0.2841 | 0.2% | 0.2839 | 0.000596 | 0 | 14,045,184 | 1,755,648 | 162 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 16 | 0.2567 | 0.2% | 0.2565 | 0.001203 | 0 | 14,045,184 | 877,824 | 292 |
| 2 | 2c-sphere-fixed | F | fnl-nvf | 32 | 0.1615 | 0.3% | 0.1612 | 0.000904 | 0 | 14,045,184 | 438,912 | 368 |
| 2 | w2-sphere-l3 | U | fnl-nvf | 1 | 0.3611 | 0.4% | 0.3606 | 2e-06 | 0 | 7,077,888 | 7,077,888 | 51 |
| 2 | w2-sphere-l3 | U | fnl-nvf | 2 | 0.1865 | 0.2% | 0.1861 | 3.3e-05 | 0 | 7,077,888 | 3,538,944 | 52.7 |
| 2 | w2-sphere-l3 | U | fnl-nvf | 4 | 0.1146 | 0.3% | 0.1144 | 2.8e-05 | 0 | 7,077,888 | 1,769,472 | 64.8 |
| 2 | w2-sphere-l5 | U | fnl-nvf | 32 | 0.9814 | 0.2% | 0.9805 | 0.00262 | 0 | 452,984,832 | 14,155,776 | 69.3 |
| 3 | 3a-channel-uniform | U | cpu | 64 | 3.171 | 6.3% (!) | 3.167 | 0.007702 | 0 | 6,291,456 | 98,304 | 3.23e+04 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 1 | 0.8 | 2.8% (!) | 0.7988 | 3e-06 | 0 | 6,291,456 | 6,291,456 | 127 |
| 3 | 3a-channel-uniform | U | fnl-nvf | 4 | 0.2296 | 0.2% | 0.229 | 0.000213 | 0 | 6,291,456 | 1,572,864 | 146 |
| 3 | 3b-channel-runtime | R | cpu | 64 | 3.082 | 4.6% (!) | 3.069 | 0.01402 | 0 | 3,538,944 | 55,296 | 5.57e+04 |
| 3 | 3b-channel-runtime | R | fnl-nvf | 4 | 0.1842 | 0.2% | 0.1428 | 0.04108 | 0 | 3,538,944 | 884,736 | 208 |
| 3 | 3c-channel-fixed | F | cpu | 64 | 2.217 | 0.9% | 2.212 | 0.003698 | 0 | 3,538,944 | 55,296 | 4.01e+04 |
| 3 | 3c-channel-fixed | F | fnl-nvf | 4 | 0.1432 | 0.5% | 0.1428 | 6.2e-05 | 0 | 3,538,944 | 884,736 | 162 |

## AMR ratios

U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency [T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.

| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |
|---|---|---|---|---|---|---|
| 1 | cpu | 64 | 1.06 | 1.59 | 0.896 | 33.9% |
| 1 | cpu | 128 | 0.543 | 1.59 | 0.743 | 117.3% |
| 1 | cpu | 256 | 0.279 | 1.59 | 0.618 | 251.9% |
| 1 | fnl-nvf | 4 | 1.04 | 1.59 | 0.885 | 34.5% |
| 1 | fnl-nvf | 8 | 0.995 | 1.59 | 0.842 | 34.2% |
| 1 | fnl-nvf | 16 | 0.931 | 1.59 | 0.776 | 32.1% |
| 2 | cpu | 64 | 2.2 | 4.02 | 0.678 | 24.2% |
| 2 | cpu | 128 | 1.64 | 4.02 | 0.514 | 26.0% |
| 2 | cpu | 256 | 1.06 | 4.02 | 0.409 | 55.7% |
| 2 | fnl-nvf | 4 | 1.48 | 4.02 | 0.507 | 37.9% |
| 2 | fnl-nvf | 8 | 1.4 | 4.02 | 0.429 | 23.2% |
| 2 | fnl-nvf | 16 | 0.834 | 4.02 | 0.251 | 21.5% |
| 2 | fnl-nvf | 32 | 0.892 | 4.02 | 0.257 | 16.2% |
| 3 | cpu | 64 | 1.03 | 1.78 | 0.805 | 39.0% |
| 3 | fnl-nvf | 4 | 1.25 | 1.78 | 0.902 | 28.6% |

## Scaling and grind time

g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the group's smallest device count. Rows sorted by backend, series and load: read g against cells/device as the grind curve.

| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | Strong eff. (ref np) | Load group | Weak eff. (ref np) |
|---|---|---|---|---|---|---|---|---|---|---|
| cpu | sphere | 2a-sphere-uniform | 64 | 56,623,104 | 884,736 | 12.68 | 1.43e+04 | 1 (64) | sphere 0.885 M/dev | - |
| cpu | sphere | 2a-sphere-uniform | 128 | 56,623,104 | 442,368 | 6.481 | 1.47e+04 | 0.978 (64) | sphere 0.442 M/dev | - |
| cpu | sphere | 2a-sphere-uniform | 256 | 56,623,104 | 221,184 | 3.244 | 1.47e+04 | 0.977 (64) | sphere 0.221 M/dev | - |
| cpu | step | 1a-step-uniform | 64 | 4,128,768 | 64,512 | 1.26 | 1.95e+04 | 1 (64) | step 0.0645 M/dev | - |
| cpu | step | 1a-step-uniform | 128 | 4,128,768 | 32,256 | 0.6415 | 1.99e+04 | 0.982 (64) | step 0.0323 M/dev | - |
| cpu | step | 1a-step-uniform | 256 | 4,128,768 | 16,128 | 0.3275 | 2.03e+04 | 0.962 (64) | step 0.0161 M/dev | - |
| cpu | - | 1b-step-fixed | 64 | 2,603,520 | 40,680 | 0.8865 | 2.18e+04 | 1 (64) | - | - |
| cpu | - | 1b-step-fixed | 128 | 2,603,520 | 20,340 | 0.5443 | 2.68e+04 | 0.814 (64) | - | - |
| cpu | - | 1b-step-fixed | 256 | 2,603,520 | 10,170 | 0.3339 | 3.28e+04 | 0.664 (64) | - | - |
| cpu | - | 1c-step-runtime | 64 | 2,602,783 | 40,668 | 1.187 | 2.92e+04 | 1 (64) | - | - |
| cpu | - | 1c-step-runtime | 128 | 2,602,783 | 20,334 | 1.182 | 5.81e+04 | 0.502 (64) | - | - |
| cpu | - | 1c-step-runtime | 256 | 2,602,783 | 10,167 | 1.175 | 1.16e+05 | 0.252 (64) | - | - |
| cpu | - | 2b-sphere-runtime | 64 | 14,068,408 | 219,819 | 5.759 | 2.62e+04 | 1 (64) | - | - |
| cpu | - | 2b-sphere-runtime | 128 | 14,068,408 | 109,909 | 3.945 | 3.59e+04 | 0.73 (64) | - | - |
| cpu | - | 2b-sphere-runtime | 256 | 14,068,408 | 54,955 | 3.064 | 5.58e+04 | 0.47 (64) | - | - |
| cpu | - | 2c-sphere-fixed | 64 | 14,045,184 | 219,456 | 4.638 | 2.11e+04 | 1 (64) | - | - |
| cpu | - | 2c-sphere-fixed | 128 | 14,045,184 | 109,728 | 3.13 | 2.85e+04 | 0.741 (64) | - | - |
| cpu | - | 2c-sphere-fixed | 256 | 14,045,184 | 54,864 | 1.968 | 3.59e+04 | 0.589 (64) | - | - |
| cpu | - | 3a-channel-uniform | 64 | 6,291,456 | 98,304 | 3.171 | 3.23e+04 | - | - | - |
| cpu | - | 3b-channel-runtime | 64 | 3,538,944 | 55,296 | 3.082 | 5.57e+04 | - | - | - |
| cpu | - | 3c-channel-fixed | 64 | 3,538,944 | 55,296 | 2.217 | 4.01e+04 | - | - | - |
| fnl-nvf | sphere | 2a-sphere-uniform | 4 | 56,623,104 | 14,155,776 | 0.7799 | 55.1 | 1 (4) | sphere 14.2 M/dev | 1 (4) |
| fnl-nvf | sphere | w2-sphere-l5 | 32 | 452,984,832 | 14,155,776 | 0.9814 | 69.3 | - | sphere 14.2 M/dev | 0.795 (4) |
| fnl-nvf | sphere | w2-sphere-l3 | 1 | 7,077,888 | 7,077,888 | 0.3611 | 51 | 1 (1) | sphere 7.08 M/dev | 1 (1) |
| fnl-nvf | sphere | 2a-sphere-uniform | 8 | 56,623,104 | 7,077,888 | 0.4911 | 69.4 | 0.794 (4) | sphere 7.08 M/dev | 0.735 (1) |
| fnl-nvf | sphere | w2-sphere-l3 | 2 | 7,077,888 | 3,538,944 | 0.1865 | 52.7 | 0.968 (1) | sphere 3.54 M/dev | 1 (2) |
| fnl-nvf | sphere | 2a-sphere-uniform | 16 | 56,623,104 | 3,538,944 | 0.2599 | 73.4 | 0.75 (4) | sphere 3.54 M/dev | 0.717 (2) |
| fnl-nvf | sphere | w2-sphere-l3 | 4 | 7,077,888 | 1,769,472 | 0.1146 | 64.8 | 0.787 (1) | sphere 1.77 M/dev | 1 (4) |
| fnl-nvf | sphere | 2a-sphere-uniform | 32 | 56,623,104 | 1,769,472 | 0.1674 | 94.6 | 0.582 (4) | sphere 1.77 M/dev | 0.685 (4) |
| fnl-nvf | step | w1-step-m8 | 1 | 1,032,192 | 1,032,192 | 0.209 | 202 | - | step 1.03 M/dev | 1 (1) |
| fnl-nvf | step | 1a-step-uniform | 4 | 4,128,768 | 1,032,192 | 0.2162 | 209 | 1 (4) | step 1.03 M/dev | 0.967 (1) |
| fnl-nvf | step | w1-step-m32 | 16 | 16,515,072 | 1,032,192 | 0.2211 | 214 | - | step 1.03 M/dev | 0.945 (1) |
| fnl-nvf | step | w1-step-m22 | 8 | 7,805,952 | 975,744 | 0.2065 | 212 | - | step 1.03 M/dev | 0.957 (1) |
| fnl-nvf | step | 1a-step-uniform | 8 | 4,128,768 | 516,096 | 0.1188 | 230 | 0.91 (4) | step 0.516 M/dev | - |
| fnl-nvf | step | 1a-step-uniform | 16 | 4,128,768 | 258,048 | 0.07043 | 273 | 0.767 (4) | step 0.258 M/dev | - |
| fnl-nvf | - | 1b-step-fixed | 4 | 2,603,520 | 650,880 | 0.1541 | 237 | 1 (4) | - | - |
| fnl-nvf | - | 1b-step-fixed | 8 | 2,603,520 | 325,440 | 0.08898 | 273 | 0.866 (4) | - | - |
| fnl-nvf | - | 1b-step-fixed | 16 | 2,603,520 | 162,720 | 0.05725 | 352 | 0.673 (4) | - | - |
| fnl-nvf | - | 1c-step-runtime | 4 | 2,602,783 | 650,696 | 0.2073 | 319 | 1 (4) | - | - |
| fnl-nvf | - | 1c-step-runtime | 8 | 2,602,783 | 325,348 | 0.1194 | 367 | 0.868 (4) | - | - |
| fnl-nvf | - | 1c-step-runtime | 16 | 2,602,783 | 162,674 | 0.07565 | 465 | 0.685 (4) | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 4 | 14,068,408 | 3,517,102 | 0.526 | 150 | 1 (4) | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 8 | 14,068,408 | 1,758,551 | 0.3499 | 199 | 0.752 (4) | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 16 | 14,068,408 | 879,276 | 0.3118 | 355 | 0.422 (4) | - | - |
| fnl-nvf | - | 2b-sphere-runtime | 32 | 14,068,408 | 439,638 | 0.1875 | 427 | 0.351 (4) | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 4 | 14,045,184 | 3,511,296 | 0.3814 | 109 | 1 (4) | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 8 | 14,045,184 | 1,755,648 | 0.2841 | 162 | 0.671 (4) | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 16 | 14,045,184 | 877,824 | 0.2567 | 292 | 0.371 (4) | - | - |
| fnl-nvf | - | 2c-sphere-fixed | 32 | 14,045,184 | 438,912 | 0.1615 | 368 | 0.295 (4) | - | - |
| fnl-nvf | - | 3a-channel-uniform | 1 | 6,291,456 | 6,291,456 | 0.8 | 127 | 1 (1) | - | - |
| fnl-nvf | - | 3a-channel-uniform | 4 | 6,291,456 | 1,572,864 | 0.2296 | 146 | 0.871 (1) | - | - |
| fnl-nvf | - | 3b-channel-runtime | 4 | 3,538,944 | 884,736 | 0.1842 | 208 | - | - | - |
| fnl-nvf | - | 3c-channel-fixed | 4 | 3,538,944 | 884,736 | 0.1432 | 162 | - | - | - |

Series:

- `sphere`: shock over a sphere, 3-D octree, WENO-5, 24^3 blocks: the weak-scaling and grind-time series
- `step`: Mach 3 step, 2-D forest, WENO-7: memory-capped 2-D, launch-bound regime (H1), ~1 M cells per A30; not a grind measurement (the block shape changes with m)

## Backend speedups

| Case | Reference | Backend | Ranks | Speedup |
|---|---|---|---|---|
| 1a-step-uniform | cpu np 64 | cpu | 128 | 1.96 |
| 1a-step-uniform | cpu np 64 | cpu | 256 | 3.85 |
| 1a-step-uniform | cpu np 64 | fnl-nvf | 4 | 5.83 |
| 1a-step-uniform | cpu np 64 | fnl-nvf | 8 | 10.6 |
| 1a-step-uniform | cpu np 64 | fnl-nvf | 16 | 17.9 |
| 1b-step-fixed | cpu np 64 | cpu | 128 | 1.63 |
| 1b-step-fixed | cpu np 64 | cpu | 256 | 2.65 |
| 1b-step-fixed | cpu np 64 | fnl-nvf | 4 | 5.75 |
| 1b-step-fixed | cpu np 64 | fnl-nvf | 8 | 9.96 |
| 1b-step-fixed | cpu np 64 | fnl-nvf | 16 | 15.5 |
| 1c-step-runtime | cpu np 64 | cpu | 128 | 1 |
| 1c-step-runtime | cpu np 64 | cpu | 256 | 1.01 |
| 1c-step-runtime | cpu np 64 | fnl-nvf | 4 | 5.72 |
| 1c-step-runtime | cpu np 64 | fnl-nvf | 8 | 9.94 |
| 1c-step-runtime | cpu np 64 | fnl-nvf | 16 | 15.7 |
| 2a-sphere-uniform | cpu np 64 | cpu | 128 | 1.96 |
| 2a-sphere-uniform | cpu np 64 | cpu | 256 | 3.91 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 4 | 16.3 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 8 | 25.8 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 16 | 48.8 |
| 2a-sphere-uniform | cpu np 64 | fnl-nvf | 32 | 75.7 |
| 2b-sphere-runtime | cpu np 64 | cpu | 128 | 1.46 |
| 2b-sphere-runtime | cpu np 64 | cpu | 256 | 1.88 |
| 2b-sphere-runtime | cpu np 64 | fnl-nvf | 4 | 10.9 |
| 2b-sphere-runtime | cpu np 64 | fnl-nvf | 8 | 16.5 |
| 2b-sphere-runtime | cpu np 64 | fnl-nvf | 16 | 18.5 |
| 2b-sphere-runtime | cpu np 64 | fnl-nvf | 32 | 30.7 |
| 2c-sphere-fixed | cpu np 64 | cpu | 128 | 1.48 |
| 2c-sphere-fixed | cpu np 64 | cpu | 256 | 2.36 |
| 2c-sphere-fixed | cpu np 64 | fnl-nvf | 4 | 12.2 |
| 2c-sphere-fixed | cpu np 64 | fnl-nvf | 8 | 16.3 |
| 2c-sphere-fixed | cpu np 64 | fnl-nvf | 16 | 18.1 |
| 2c-sphere-fixed | cpu np 64 | fnl-nvf | 32 | 28.7 |
| 3a-channel-uniform | cpu np 64 | fnl-nvf | 1 | 3.96 |
| 3a-channel-uniform | cpu np 64 | fnl-nvf | 4 | 13.8 |
| 3b-channel-runtime | cpu np 64 | fnl-nvf | 4 | 16.7 |
| 3c-channel-fixed | cpu np 64 | fnl-nvf | 4 | 15.5 |

## Warnings

- cpu np 64 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 2.9% above 2%
- cpu np 128 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 2.3% above 2%
- cpu np 256 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 2.8% above 2%
- fnl-nvf np 4 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 1a-step-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes); spread 2.8% above 2%
- cpu np 128 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 256 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 1b-step-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes); spread 3.2% above 2%
- cpu np 128 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes); spread 3.8% above 2%
- cpu np 256 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 1c-step-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 w1-step-m22: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 w1-step-m32: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 w1-step-m8: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 7.5% above 2%
- cpu np 128 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 256 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 2a-sphere-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes); spread 6.3% above 2%
- cpu np 128 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 256 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 2b-sphere-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 128 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 256 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 8 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 16 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 2c-sphere-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 1 w2-sphere-l3: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 2 w2-sphere-l3: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 w2-sphere-l3: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 32 w2-sphere-l5: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 3a-channel-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 6.3% above 2%
- fnl-nvf np 1 3a-channel-uniform: built from sources that differ from the recorded commit (src_modified=yes); spread 2.8% above 2%
- fnl-nvf np 4 3a-channel-uniform: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 3b-channel-runtime: built from sources that differ from the recorded commit (src_modified=yes); spread 4.6% above 2%
- fnl-nvf np 4 3b-channel-runtime: built from sources that differ from the recorded commit (src_modified=yes)
- cpu np 64 3c-channel-fixed: built from sources that differ from the recorded commit (src_modified=yes)
- fnl-nvf np 4 3c-channel-fixed: built from sources that differ from the recorded commit (src_modified=yes)
