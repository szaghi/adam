#!/usr/bin/env bash
# FLUME verification VV (issue #65, P2): the Navier-Stokes dissipative fluxes.
#
# Why: a dissipative flux can be consistent and still wrong in its order (a wrong correction term), its coefficient
# (a missing 4/3, a k without R) or its direction (a cross term read with the wrong stencil). Each leg has an exact
# solution and an asserted bound:
#   vv1  shear wave (waves.py, linearised NS mode): measured order of the finest pair >= ORDER4_MIN at
#        dissipative_order = 4 along x, along y and along the diagonal (cross derivatives), >= ORDER2_MIN at order 2;
#   vv2  viscous-thermal acoustic wave (mu and k together): the same order bound at order 4, on the odd part of the
#        +A / -A twin runs (an acoustic wave steepens at O(A^2), which would otherwise floor the error);
#   vv3  compressible Couette flow between an isothermal wall and a moving adiabatic one (couette.py);
#   vv4  Becker's viscous shock (becker.py).
#
# Usage: ./check.sh [--np N] [--leg vv1|vv2|vv3|vv4 ...]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEGS=()
ORDER4_MIN="3.8"
ORDER2_MIN="1.9"
BECKER_ORDER_MIN="3.0"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)  NP="$2" ; shift 2 ;;
      --leg) LEGS+=("$2") ; shift 2 ;;
      *)     echo "check.sh: unknown argument '$1' (accepted: --np N, --leg L)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(vv1 vv2 vv3 vv4)
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

run() { (cd "$1" && mpirun -np "$NP" "$EXE" input.ini < /dev/null > log.txt 2>&1) ; } # mpirun reads stdin

ladder() { # ladder <tool> <name> <order min> "<resolutions>" <make options...>: run, judge the finest pair
   # LADDER_TWIN=1 also runs each resolution with the amplitude -A into <work>-neg (waves.py: the odd part is judged)
   local tool="$CASE_DIR/$1" name="$2" pmin="$3" ns=($4) works=() w v n
   shift 4
   echo "-- $name"
   for n in "${ns[@]}"; do
      w="$CASE_DIR/work-$TAG-$name-n$n"
      for v in "$w" ${LADDER_TWIN:+"$w-neg"}; do
         rm -rf "$v" ; mkdir -p "$v"
         if [[ $v == *-neg ]]; then
            "$VENV_PY" "$tool" make "$v/input.ini" --n "$n" "$@" --sign -1
         else
            "$VENV_PY" "$tool" make "$v/input.ini" --n "$n" "$@"
         fi
         if ! run "$v"; then echo "   run failed, see $v/log.txt" ; FAILED=1 ; return 0 ; fi
         rm -f "$v"/*-000000000-*.h5
      done
      works+=("$w")
   done
   # pipefail: a failing oracle fails the pipeline, which marks the leg failed instead of aborting the check (set -e)
   if ! "$VENV_PY" "$tool" oracle "${works[@]}" --order-min "$pmin" | sed 's/^/   /'; then FAILED=1 ; fi
}

echo ">> VV: Navier-Stokes dissipative fluxes ($TAG)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      vv1)
         ladder waves.py shear-x-o4 "$ORDER4_MIN" "32 64 128" --mode shear --angle 0 --order 4
         ladder waves.py shear-x-o2 "$ORDER2_MIN" "32 64 128" --mode shear --angle 0 --order 2
         ladder waves.py shear-y-o4 "$ORDER4_MIN" "32 64 128" --mode shear --angle 90 --order 4
         ladder waves.py shear-xy-o4 "$ORDER4_MIN" "24 48 96" --mode shear --angle 45 --order 4 --time 0.1 ;;
      vv2)
         LADDER_TWIN=1 ladder waves.py acoustic-x-o4 "$ORDER4_MIN" "32 64 128" --mode acoustic --angle 0 --order 4 \
                --mu 0.01 --kappa 0.02 ;;
      vv3)
         ladder profiles.py couette-o4 "$ORDER2_MIN" "32 64 128" --case couette --order 4 ;;
      vv4)
         ladder profiles.py becker-m2-o4 "$BECKER_ORDER_MIN" "64 128 256" --case becker --mach 2 --order 4 --time 0.05
         ladder profiles.py becker-m3-o4 "$BECKER_ORDER_MIN" "64 128 256" --case becker --mach 3 --order 4 --time 0.05 \
                --mu 0.01 ;;
      *) echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "VV PASSED ($TAG)"
else
   echo "VV FAILED ($TAG)"
   exit 1
fi
