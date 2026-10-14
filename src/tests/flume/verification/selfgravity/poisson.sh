#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8). SV-2 (P2): the FFT Poisson solver on uniform periodic grids.
#
# Why: self-gravity conserves E + rho Phi and momentum exactly only if Phi solves the discrete Poisson equation
# L_h Phi = 4 pi G (rho - mean rho) to round-off (#81 P0). The library solver (adam_poisson_fft_object, FFTW; its FNL
# twin, cuFFT on rank 0's GPU) divides each Fourier mode by the eigenvalue of L_h, gathers the blocks to rank 0 and
# scatters Phi back with every ghost layer filled by the periodic wrap. The unit test test_poisson_fft holds it to:
#   - the discrete solution (two oblique modes on [0,1] x [0,2] x [0,3], each an eigenfunction of L_h), every cell and
#     ghost, at round-off (1e-12 relative);
#   - the discrete residual from the scattered ghosts (1e-10 relative);
#   - the continuum solution, second order (>= 1.9 on the finest pair, N = 32 -> 64);
# on a 3-D octree, a 2-D quadtree and a single block (ranks without blocks when np > 1). Each rank count of the list
# must pass; a non-periodic grid must be refused (exit 1, the message naming the direction).
# Measured (CPU and FNL): discrete <= 2.3e-15, residual <= 4.3e-14, orders 1.98 (3-D) and 1.98 (2-D); the printed
# numbers do not depend on the rank count.
#
# CPU = FNL: run this script once per build with --dump DIR, then
#   poisson_oracle.py compare DIR_CPU DIR_FNL      (max |phi_CPU - phi_FNL| / max |phi|, 1e-13)
#
# Usage: ./poisson.sh [--np-list "1 2 3"] [--dump DIR] [--keep]
#
# POISSON_EXE overrides the executable under test, e.g. POISSON_EXE=$REPO/exe/test_poisson_fft_fnl ./poisson.sh
# (built by `fobis build --mode test-poisson-fft-gnu` or `test-poisson-fft-fnl-nvf`). The caller owns the matching
# environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${POISSON_EXE:-$REPO_ROOT/exe/test_poisson_fft}"
NP_LIST="1 2 3"
DUMP=""
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np-list) NP_LIST="$2" ; shift 2 ;;
      --dump)    DUMP="$2" ; shift 2 ;;
      --keep)    KEEP=1 ; shift ;;
      *)         echo "poisson.sh: unknown argument '$1' (accepted: --np-list LIST, --dump DIR, --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "poisson.sh: executable '$EXE' not found (build it or set POISSON_EXE)" >&2 ; exit 2
fi
EXE="$(realpath "$EXE")"

TAG="$(basename "$EXE")"
work="$CASE_DIR/work-$TAG-sv2"
rm -rf "$work" ; mkdir -p "$work"
FAILED=0
echo ">> SV-2: the FFT Poisson solver ($TAG)"
for np in $NP_LIST; do
   log="$work/log-np$np.txt"
   if ! (cd "$work" && mpirun -np "$np" "$EXE" > "$log" 2>&1); then
      echo "   np $np: the run failed, see $log  FAIL"
      FAILED=1
      continue
   fi
   grep '^   ' "$log" | sed "s/^/   np $np:/"
   if ! grep -q 'test_poisson_fft: PASSED' "$log"; then FAILED=1 ; fi
done
log="$work/log-refusal.txt"
if (cd "$work" && mpirun -np 2 "$EXE" --no-periodic > "$log" 2>&1); then
   echo "   refusal: a non-periodic grid was accepted  FAIL"
   FAILED=1
elif grep -q 'direction 1 is not periodic' "$log"; then
   echo "   refusal: a non-periodic grid stops the run (direction 1 is not periodic)  PASS"
else
   echo "   refusal: the run failed without the expected message, see $log  FAIL"
   FAILED=1
fi
if [[ -n "$DUMP" ]]; then
   mkdir -p "$DUMP"
   DUMP="$(cd "$DUMP" && pwd)"
   if ! (cd "$work" && mpirun -np 2 "$EXE" --dump "$DUMP" > "$work/log-dump.txt" 2>&1); then
      echo "   dump: the run failed, see $work/log-dump.txt  FAIL"
      FAILED=1
   else
      echo "   dump: $(ls "$DUMP"/phi-*.bin | wc -l) potentials in $DUMP"
   fi
fi
if [[ $KEEP -eq 0 && $FAILED -eq 0 ]]; then rm -rf "$work" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "SV-2 FAILED ($TAG)" ; exit 1
fi
echo "SV-2 PASSED ($TAG)"
