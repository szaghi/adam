#!/usr/bin/env bash
# Build ghost_fill_bench for one GPU architecture and run the size matrix (README "Runs"); one result file.
#   ./run.sh <cc> <out>      e.g. ./run.sh cc90,cuda12.9 results/h100-<job>.txt
# The caller loads the toolchain (nvfortran on PATH). One GPU is used (device 0 of what is visible).
set -euo pipefail
cd "$(dirname "$0")"
CC=${1:?run.sh: the -gpu target, e.g. cc90,cuda12.9}
OUT=${2:?run.sh: the result file}
mkdir -p "$(dirname "$OUT")" build
nvfortran -acc -gpu="$CC" -fast -Minfo=accel -o build/ghost_fill_bench ghost_fill_bench.F90 > build/minfo.txt 2>&1
{
   echo "# host $(hostname) $(date -Is); $(nvfortran --version | grep -m1 nvfortran); -gpu=$CC"
   nvidia-smi --query-gpu=index,name,memory.total,clocks.max.sm,clocks.sm,clocks_event_reasons.active,power.draw \
      --format=csv,noheader
   # nbx nby nbz b ngc nv nrep: ~7 M cells per GPU at five block edges, ~56 M at three, WENO-7's ngc 4 once
   for cfg in "12 12 12 16 3" "8 8 8 24 3" "6 6 6 32 3" "4 4 8 38 3" "4 4 4 48 3" \
              "16 16 16 24 3" "8 8 16 38 3" "8 8 8 48 3" "4 4 8 38 4"; do
      ./build/ghost_fill_bench $cfg 5 50
      echo
   done
   nvidia-smi --query-gpu=index,clocks.sm,clocks_event_reasons.active --format=csv,noheader
} > "$OUT" 2>&1
echo "run.sh: $OUT"
