#!/usr/bin/env bash
# pml-restart: restart leg of the classic PML (issue #129 R3), see ../pic_check.sh (restart-only mode).
# A Gaussian pulse (EM, no particles) runs into a classic PML on x+ (PEC behind it): 200x6x6 cells as 8 blocks of
# 100x3x3 (one uniform refinement), so at np 2 each rank holds 2 of the 4 blocks on the layer. The run is stopped at
# step 1400, with the pulse inside the layer (peak at x = 0.0865, layer 0.08..0.1), restarted, and compared at 1700
# with the uninterrupted run. A restart without the PML auxiliary state (.pml) gives Bz off by 2.3 on 2.9 at 1700.
# Legs: CPU np 1 and np 2, FNL np 1 and np 2 (FNL_NP). No golden: the uninterrupted run is the reference.
# Times of one run (2026-11-24, the restart leg about one more): CPU np 1 / np 2 ~16 / ~20 s, FNL np 1 ~3.5 min (WSL,
# RTX 4070 laptop: 8 small blocks, kernel-launch bound), FNL np 2 ~1 min (two A30, medooza). All bit-identical.
# Usage: ./check.sh <cpu|fnl>
PIC_CHECK_RESTART_ONLY=1 PIC_CHECK_RESTART_AT=1400 \
   exec "$(cd "$(dirname "$0")" && pwd)/../pic_check.sh" "$(cd "$(dirname "$0")" && pwd)" pml_restart 1700 "$@"
