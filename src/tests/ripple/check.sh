#!/usr/bin/env bash
# RIPPLE oracles (issue #59): the toy app of the "build your own solver" tutorial against the exact solution.
#
# RIPPLE advects u0 = sin(2 pi x) sin(2 pi y) at a = (1, 0.5) across the periodic unit square; at t = 1 the exact
# solution is u0(x - a t). oracle.py reads the checkpoints (XH5F, point values at the cell centres) and the history of
# the total of u dV. The grids are 4x4 blocks (a quadtree, z null) of N/4 x N/4 cells. Legs:
#   r1  periodic.ini, the uniform grid (the forest's fast path), N = 32, 64, 128 on 2 ranks: the L1 order of the
#       finest pair at least 4.8 (WENO-5 with SSP-54 at CFL 0.4; measured 5.00, 5.00) and the finest L1 error at most
#       the measured 1.81e-7 plus 2 %;
#   r2  amr.ini, the four central blocks refined one level (the staged path: flux register and end-of-step reflux), the
#       same ladder: the L1 order of the finest pair at least 1.9 and the finest L1 at most the measured 3.12e-5 plus
#       2 %; the total of u dV constant to 1e-14 of max|u0| times the volume on every grid (the reflux conserves). The
#       order is 2, not 5, by design of the library's 2:1 transfers: the mean restriction of the coarse ghosts and the
#       reflux of point-value fluxes are second order (as FLUME's AO leg measures, issue #68);
#   r3  the AMR case at N = 64 on 1 and 2 ranks: bitwise identical;
#   r5  CRLF input: the periodic case at N = 64 on 1 rank with LF and with CR+LF line endings, bitwise identical, on
#       RIPPLE_R5_EXE (default the CPU executable; the leg that matters runs the FNL one: nvfortran keeps the CR at the
#       end of a formatted record, gfortran drops it; FiNeR drops every CR since c55f8d3);
#   r4  CPU and FNL agree: the periodic and the AMR cases at N = 64 on 2 ranks with the FNL executable, against the
#       CPU runs of r1 and r2 (run those first, in the CPU environment: the two executables need their own MPI), within
#       1e-12; the FNL total of u dV constant as in r2.
#
# Usage: ./check.sh [--legs LIST]   (legs r1 r2 r3 r4 r5, comma-separated; default r1,r2,r3,r5)
#
# RIPPLE_CPU_EXE, RIPPLE_FNL_EXE override the executables (default exe/adam_ripple_cpu, exe/adam_ripple_fnl), MPIRUN
# the launcher (default "mpirun"), PYTHON the interpreter (default the regression venv, which has h5py); the caller owns
# the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../.." && pwd)"
CPU_EXE="${RIPPLE_CPU_EXE:-$REPO_ROOT/exe/adam_ripple_cpu}"
FNL_EXE="${RIPPLE_FNL_EXE:-$REPO_ROOT/exe/adam_ripple_fnl}"
MPIRUN="${MPIRUN:-mpirun}"
PYTHON="${PYTHON:-$REPO_ROOT/exe/.regression-venv/bin/python}"
LEGS="r1,r2,r3,r5"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) LEGS="$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
for leg in ${LEGS//,/ }; do
   case "$leg" in
      r1|r2|r3|r4|r5) ;;
      *) echo "check.sh: unknown leg '$leg' (accepted: r1 r2 r3 r4 r5)" >&2 ; exit 2 ;;
   esac
done
ORACLE=("$PYTHON" "$CASE_DIR/oracle.py")
FAILED=0

# run EXE INPUT NI NP DIR: write INPUT with ni = nj = NI into DIR, run it on NP ranks
run() {
   local exe="$1" input="$2" ni="$3" np="$4" dir="$5"
   rm -rf "$dir" && mkdir -p "$dir"
   sed -e "s/^ni .*/ni     = $ni/" -e "s/^nj .*/nj     = $ni/" "$CASE_DIR/$input" > "$dir/input.ini"
   if ! (cd "$dir" && $MPIRUN -np "$np" "$exe" input.ini > log.txt 2>&1); then
      echo "FAIL: $exe $input ni=$ni np=$np exited non-zero (see $dir/log.txt)" ; return 1
   fi
}

# verdict LEG COMMAND...: run an oracle check, report PASS or FAIL
verdict() {
   local leg="$1" ; shift
   if "$@"; then echo "PASS [$leg]" ; else echo "FAIL [$leg]" ; FAILED=1 ; fi
}

has_leg() { [[ ",$LEGS," == *",$1,"* ]]; }

WORK="$CASE_DIR/work"
mkdir -p "$WORK"
if has_leg r1; then
   for ni in 8 16 32; do run "$CPU_EXE" periodic.ini $ni 2 "$WORK/r1-$ni" || FAILED=1 ; done
   verdict r1 "${ORACLE[@]}" convergence "$WORK"/r1-{8,16,32} --min-order 4.8 --max-l1 1.85e-7
fi
if has_leg r2; then
   for ni in 8 16 32; do run "$CPU_EXE" amr.ini $ni 2 "$WORK/r2-$ni" || FAILED=1 ; done
   verdict r2 "${ORACLE[@]}" convergence "$WORK"/r2-{8,16,32} --min-order 1.9 --max-l1 3.19e-5
   for ni in 8 16 32; do verdict "r2 conservation N=$((4*ni))" "${ORACLE[@]}" conservation "$WORK/r2-$ni" --tol 1e-14 ; done
fi
if has_leg r3; then
   run "$CPU_EXE" amr.ini 16 1 "$WORK/r3-np1" || FAILED=1
   run "$CPU_EXE" amr.ini 16 2 "$WORK/r3-np2" || FAILED=1
   verdict r3 "${ORACLE[@]}" compare "$WORK/r3-np1" "$WORK/r3-np2" --tol 0
fi
if has_leg r4; then
   for pair in periodic:r1 amr:r2; do
      input="${pair%%:*}" ; cpu="$WORK/${pair##*:}-16"
      if [[ ! -d "$cpu" ]]; then echo "FAIL [r4 $input]: no CPU run $cpu (run leg ${pair##*:} first)" ; FAILED=1 ; continue ; fi
      run "$FNL_EXE" $input.ini 16 2 "$WORK/r4-$input-fnl" || FAILED=1
      verdict "r4 $input" "${ORACLE[@]}" compare "$cpu" "$WORK/r4-$input-fnl" --tol 1e-12
   done
   verdict "r4 conservation" "${ORACLE[@]}" conservation "$WORK/r4-amr-fnl" --tol 1e-14
fi
if has_leg r5; then
   R5_EXE="${RIPPLE_R5_EXE:-$CPU_EXE}"
   run "$R5_EXE" periodic.ini 16 1 "$WORK/r5-lf" || FAILED=1
   rm -rf "$WORK/r5-crlf" && mkdir -p "$WORK/r5-crlf" && sed 's/$/\r/' "$WORK/r5-lf/input.ini" > "$WORK/r5-crlf/input.ini"
   if ! (cd "$WORK/r5-crlf" && $MPIRUN -np 1 "$R5_EXE" input.ini > log.txt 2>&1); then
      echo "FAIL: $R5_EXE on the CRLF input exited non-zero (see $WORK/r5-crlf/log.txt)" ; FAILED=1
   fi
   verdict r5 "${ORACLE[@]}" compare "$WORK/r5-lf" "$WORK/r5-crlf" --tol 0
fi
if [[ $FAILED -eq 0 ]]; then echo "RIPPLE oracles ($LEGS) PASSED" ; else echo "RIPPLE oracles ($LEGS) FAILED" ; fi
exit $FAILED
