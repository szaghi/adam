#!/usr/bin/env bash
# FLUME verification of self-gravity (issue #81, M8). SV-2 (P2): the FFT Poisson solver on uniform periodic grids;
# SV-6 (P4): the same solver with isolated boundaries.
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
# SV-6 (test_poisson_fft --isolated): Hockney-Eastwood zero padding with the lattice Green's function of L_h
# (adam_lattice_green_library, scripts/lattice_green.py), phi -> 0 at infinity, on the cube [-1,1]^3 (16^3 to 64^3):
#   - an n = 1 polytrope (smooth) and a uniform sphere (cell averages), radius 1/2, against their exact potentials,
#     every cell and ghost: least-squares order >= 1.9 (polytrope) and >= 1.2 (sphere: its surface cuts cells);
#   - the residual of L_h, 1e-10 relative (the kernel's own residual summed over the source, not round-off);
#   - refused, each with its message: a periodic grid, non-cubic cells, a 2-D grid.
# Measured (CPU): orders 2.03 (polytrope) and 1.35 (sphere), residual <= 2.1e-11 at 64^3.
#
# CPU = FNL: run this script once per build with --dump DIR, then
#   poisson_oracle.py compare DIR_CPU DIR_FNL      (max |phi_CPU - phi_FNL| / max |phi|, 1e-13)
#
# Usage: ./poisson.sh [--np-list "1 2 3"] [--legs "sv2 sv6"] [--dump DIR] [--keep]
#
# POISSON_EXE overrides the executable under test, e.g. POISSON_EXE=$REPO/exe/test_poisson_fft_fnl ./poisson.sh
# (built by `fobis build --mode test-poisson-fft-gnu` or `test-poisson-fft-fnl-nvf`). The caller owns the matching
# environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${POISSON_EXE:-$REPO_ROOT/exe/test_poisson_fft}"
NP_LIST="1 2 3"
LEGS="sv2 sv6"
DUMP=""
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np-list) NP_LIST="$2" ; shift 2 ;;
      --legs)    LEGS="$2" ; shift 2 ;;
      --dump)    DUMP="$2" ; shift 2 ;;
      --keep)    KEEP=1 ; shift ;;
      *)         echo "poisson.sh: unknown argument '$1' (accepted: --np-list LIST, --legs LIST, --dump DIR, --keep)" >&2
                 exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "poisson.sh: executable '$EXE' not found (build it or set POISSON_EXE)" >&2 ; exit 2
fi
EXE="$(realpath "$EXE")"

TAG="$(basename "$EXE")"
work="$CASE_DIR/work-$TAG-poisson"
rm -rf "$work" ; mkdir -p "$work"
FAILED=0
suite() { # suite <label> <flags>: the test on every rank count of the list, its lines and verdict
   local label="$1" flags="$2" np log
   for np in $NP_LIST; do
      log="$work/log-$label-np$np.txt"
      # shellcheck disable=SC2086
      if ! (cd "$work" && mpirun -np "$np" "$EXE" $flags > "$log" 2>&1); then
         echo "   np $np: the run failed, see $log  FAIL"
         FAILED=1
         continue
      fi
      grep '^   ' "$log" | sed "s/^/   np $np:/"
      if ! grep -q 'test_poisson_fft: PASSED' "$log"; then FAILED=1 ; fi
   done
}
refusal() { # refusal <label> <flag> <expected message>: the run must stop, naming the cause
   local label="$1" flag="$2" msg="$3" log="$work/log-refusal-$1.txt"
   if (cd "$work" && mpirun -np 2 "$EXE" "$flag" > "$log" 2>&1); then
      echo "   refusal, $label: accepted  FAIL"
      FAILED=1
   elif grep -q "$msg" "$log"; then
      echo "   refusal, $label: stops the run ($msg)  PASS"
   else
      echo "   refusal, $label: the run failed without the expected message, see $log  FAIL"
      FAILED=1
   fi
}
for leg in $LEGS; do
   case "$leg" in
      sv2)
         echo ">> SV-2: the FFT Poisson solver, periodic ($TAG)"
         suite periodic ""
         refusal "a non-periodic grid" --no-periodic 'direction 1 is not periodic' ;;
      sv6)
         echo ">> SV-6: the FFT Poisson solver, isolated ($TAG)"
         suite isolated "--isolated"
         refusal "isolated on a periodic grid" --refuse-isolated-periodic 'direction 1 is periodic'
         refusal "isolated with non-cubic cells" --refuse-isolated-anisotropic 'isolated boundaries need cubic cells'
         refusal "isolated in 2-D" --refuse-isolated-2d 'isolated boundaries are 3-D only' ;;
      *) echo "poisson.sh: unknown leg '$leg' (sv2, sv6)" >&2 ; exit 2 ;;
   esac
done
if [[ -n "$DUMP" ]]; then
   mkdir -p "$DUMP"
   DUMP="$(cd "$DUMP" && pwd)"
   for flags in "" "--isolated"; do
      # shellcheck disable=SC2086
      if ! (cd "$work" && mpirun -np 2 "$EXE" $flags --dump "$DUMP" > "$work/log-dump.txt" 2>&1); then
         echo "   dump ${flags:-periodic}: the run failed, see $work/log-dump.txt  FAIL"
         FAILED=1
      fi
   done
   echo "   dump: $(ls "$DUMP"/phi-*.bin | wc -l) potentials in $DUMP"
fi
if [[ $KEEP -eq 0 && $FAILED -eq 0 ]]; then rm -rf "$work" ; fi
if [[ $FAILED -ne 0 ]]; then
   echo "Poisson legs ($LEGS) FAILED ($TAG)" ; exit 1
fi
echo "Poisson legs ($LEGS) PASSED ($TAG)"
