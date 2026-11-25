#!/usr/bin/env bash
# Build test_amr_ghost_bench (FoBiS mode test-amr-ghost-bench-fnl-nvf) and run its matrix on one GPU; one result file.
#   ./run.sh <varset> <out>      e.g. ./run.sh iac_nvf results/a30-<job>.txt (from any directory)
# The caller loads the toolchain (nvfortran and its MPI on PATH). Runs happen in a scratch directory next to <out>.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../../.." && pwd)"
VARSET=${1:?run.sh: the varset, e.g. iac_nvf}
OUT=${2:?run.sh: the result file}
mkdir -p "$(dirname "$OUT")"           # before readlink -f, which fails on a missing parent
OUT=$(readlink -f "$OUT")
(cd "$REPO" && fobis build --mode test-amr-ghost-bench-fnl-nvf --varset "$VARSET" -j 16 > "$OUT.build.log" 2>&1) \
   || { tail -30 "$OUT.build.log"; exit 3; }
EXE="$REPO/exe/test_amr_ghost_bench_fnl"
WORK="$(dirname "$OUT")/work-$$" ; mkdir -p "$WORK" ; cd "$WORK"
run() { echo "## $*" ; mpirun -np 1 "$EXE" "$@" 2>&1 | grep -E '^#|^V0|^D[0-9]|differ|rror|STOP' ; echo ; }
{
   echo "# host $(hostname) $(date -Is); $(nvfortran --version | grep -m1 nvfortran); varset $VARSET"
   nvidia-smi --query-gpu=index,name,memory.total,clocks.max.sm,clocks_event_reasons.active --format=csv,noheader
   # b base box fill restriction nrep mem_gb
   if [[ ${RUN_ONLY:-} == large ]]; then
      for pair in "tricubic mean" "weno6 weno6"; do run 24 3 0.5 $pair 20 16.0 40; done
   else
   # 1) the production-like grid (2c: 24^3 blocks, ~254 per GPU): every regime pair
   for fill in tricubic quintic weno6 quintic-mp; do
      for restr in mean lagrange6 weno6; do run 24 2 0.75 $fill $restr 30 16.0; done
   done
   # 2) the refined fraction, the default pair and the WENO pair
   for box in 0.25 0.5; do
      for pair in "tricubic mean" "weno6 weno6"; do run 24 2 $box $pair 30 16.0; done
   done
   # 3) the block edge, at the refined half-box
   for b in 16 32; do
      for pair in "tricubic mean" "quintic lagrange6"; do run $b 2 0.5 $pair 30 16.0; done
   done
   # 4) a larger load (base level 3: 960 blocks of 24^3, 13 M cells)
   for pair in "tricubic mean" "weno6 weno6"; do run 24 3 0.5 $pair 20 16.0 40; done
   fi
   nvidia-smi --query-gpu=index,clocks.sm,clocks_event_reasons.active --format=csv,noheader
} > "$OUT" 2>&1
cd .. && rm -rf "$WORK"
echo "run.sh: $OUT"
