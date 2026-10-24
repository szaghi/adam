#!/usr/bin/env bash
# Scaling ladders of FLAIL's solves (issue #84, E5), on exe/bench_poisson (CPU, FFTW) or exe/bench_poisson_fnl (FNL,
# cuFFT on rank 0's GPU).
#
# Why: both FLAIL solvers gather to rank 0 (the FFT solver the whole uniform grid, the multigrid its coarsest level in
# every V-cycle and the James step's ghost charges in every isolated solve). The ladders measure what that costs as the
# ranks grow, phase by phase (gather, rank-0 solve, scatter), before and after E5 replaces the gathers.
#
# Ladders (octree, the cube [-1,1]^3, 8^3 or 16^3 cells per block; see bench_poisson.F90):
#   - weak: the cells per rank fixed, (np, level) = (1, 3), (8, 4), (64, 5), on four paths: FFT periodic, FFT
#     isolated, multigrid on AMR (the box [-1/2, 1/2]^3 refined once) periodic and isolated, with ni = 8 (262 144
#     uniform cells per rank); the two FFT paths also with ni = 16 (2 097 152; the multigrid there would take hours:
#     its coarsest level, the whole base grid, is solved on rank 0 in every V-cycle);
#   - strong: level 4, ni = 8 (2 097 152 cells uniform), np = 1, 2, 4, 8, 16, 32, 64, the same four paths.
# One BENCH line per run goes to $BENCH_DIR/bench.txt; the table printed at the end has the solve time (min over the
# repetitions, max over the ranks) and rank 0's phases per solve.
#
# Usage: src/tests/poisson/bench.sh [cpu|fnl] [--weak] [--strong] [--np-max N]
#   BENCH_DIR (default ./bench-poisson-work), BENCH_EXE, BENCH_REPS (default 5); MPIRUN (default mpirun).
#   FNL: one GPU per rank, so np-max is the GPUs available (default 4); the weak ladder stops at np-max.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
LEG=cpu ; WEAK=0 ; STRONG=0 ; NPMAX=64
while [[ $# -gt 0 ]]; do
   case "$1" in
      cpu|fnl)  LEG="$1" ; shift ;;
      --weak)   WEAK=1 ; shift ;;
      --strong) STRONG=1 ; shift ;;
      --np-max) NPMAX="$2" ; shift 2 ;;
      *) echo "bench.sh: unknown argument '$1' (accepted: cpu|fnl --weak --strong --np-max N)" >&2 ; exit 2 ;;
   esac
done
[[ $WEAK -eq 0 && $STRONG -eq 0 ]] && { WEAK=1 ; STRONG=1 ; }
if [[ $LEG == fnl ]]; then
   EXE="${BENCH_EXE:-$REPO_ROOT/exe/bench_poisson_fnl}" ; [[ $NPMAX -eq 64 ]] && NPMAX=4
else
   EXE="${BENCH_EXE:-$REPO_ROOT/exe/bench_poisson}"
fi
DIR="${BENCH_DIR:-$PWD/bench-poisson-work}"
REPS="${BENCH_REPS:-5}"
MPIRUN="${MPIRUN:-mpirun}"
[[ -x "$EXE" ]] || { echo "bench.sh: $EXE not found (build bench-poisson-gnu or bench-poisson-fnl-nvf)" >&2 ; exit 2 ; }
mkdir -p "$DIR" || exit 2
OUT="$DIR/bench.txt"
PATHS=("" "--isolated" "--amr" "--amr --isolated")
status=0

run() { # run <np> <level> <ni> <path flags>: one benchmark, its BENCH line appended to $OUT
   local np=$1 level=$2 ni=$3 flags=$4 log
   log="$DIR/log-$LEG-np$np-l$level-n$ni${flags// /}.txt"
   if (cd "$DIR" && $MPIRUN -np "$np" "$EXE" --level "$level" --ni "$ni" --reps "$REPS" $flags < /dev/null > "$log" 2>&1) \
      && grep -aq '^BENCH ' "$log"; then
      grep -a '^BENCH ' "$log" | sed "s/^BENCH /BENCH leg=$LEG /" >> "$OUT"
      echo "  ok   np $np level $level ni $ni $flags"
   else
      echo "  FAIL np $np level $level ni $ni $flags (see $log)"
      status=1
   fi
}

if [[ $WEAK -eq 1 ]]; then
   echo "weak ladder ($LEG)"
   for ni in 8 16; do
      for f in "${PATHS[@]}"; do
         [[ $ni -eq 16 && $f == *--amr* ]] && continue
         for pair in "1 3" "8 4" "64 5"; do
            set -- $pair
            [[ $1 -le $NPMAX ]] && run "$1" "$2" "$ni" "$f"
         done
      done
   done
fi
if [[ $STRONG -eq 1 ]]; then
   echo "strong ladder ($LEG)"
   for f in "${PATHS[@]}"; do
      for np in 1 2 4 8 16 32 64; do
         [[ $np -le $NPMAX ]] && run "$np" 4 8 "$f"
      done
   done
fi

echo
printf '%-4s %-4s %-3s %-3s %-4s %-4s %-11s %-5s %-10s %-10s %-10s %-10s %-10s\n' leg path iso amr np lev cells its \
   solve_min gather rank0 scatter james_r0
awk '{
   for (i = 2; i <= NF; i++) { split($i, kv, "=") ; v[kv[1]] = kv[2] }
   g = v["fft_gather"] + v["coarse_gather"] ; r = v["fft_rank0"] + v["coarse_rank0"] ; s = v["fft_scatter"] + v["coarse_scatter"]
   printf "%-4s %-4s %-3s %-3s %-4s %-4s %-11s %-5s %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n", v["leg"], v["path"], \
      v["isolated"], v["amr"], v["np"], v["level"], v["cells"], v["iterations"], v["solve_min"], g, r, s, v["james_rank0"]
}' "$OUT"
exit $status
