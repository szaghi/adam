#!/usr/bin/env bash
# FLUME verification PS: passive scalars, `[physics] passive_scalars = N` (issue #92).
#
# Why: the scalar models carry NS_MAX = 4 more slots, rho s_k, advected by the same flux machinery as the flow; a scalar
# must not act on the flow, must be conserved, and the slots beyond N must stay empty. The legs:
#   ps1  Sod (sod-x) with two scalars, 1 left and 0 right, on euler and on euler-de, each against the same run without
#        scalars (make_ps.py): the flow fields (and the entropy of euler-de) bitwise equal, each scalar integral
#        constant to round-off (1e-13), the two unused slots identically zero; the range of s_k reported;
#   ps2  the contract: passive_scalars = 5 (above NS_MAX) and passive_scalars with mhd-ideal refused.
# On FNL only ps2's refusal of the scalar models is checked until #92 P2.
#
# Usage: ./check.sh [--legs LIST]   (legs ps1 ps2, all; default all)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "ps1 ps2"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ ! -x "$EXE" ]]; then echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; fi
EXE="$(realpath "$EXE")"
TAG="$(basename "$EXE")"
FNL=0 ; [[ $TAG == *fnl* ]] && FNL=1
FAILED=0

run() { # run <work dir>: run the input of the work dir; 0 on success
   (cd "$1" && mpirun -np 2 "$EXE" input.ini < /dev/null > log.txt 2>&1)
}
make_case() { # make_case <name> [make_ps.py options]: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_ps.py" "$W/input.ini" "$@"
}
expect_refused() { # expect_refused <work dir> <message regex>
   if run "$1"; then
      echo "   accepted, expected a refusal: see $1/log.txt  FAIL" ; FAILED=1
   elif grep -aqE "$2" "$1/log.txt"; then
      echo "   refused as expected"
   else
      echo "   failed without the expected message ($2): see $1/log.txt  FAIL" ; FAILED=1
   fi
}

echo ">> PS: passive scalars ($TAG)"
for leg in $LEGS; do
   case $leg in
      ps1)
         if [[ $FNL -eq 1 ]]; then echo "-- ps1: CPU only until #92 P2 (skipped)" ; continue ; fi
         for model in euler euler-de; do
            echo "-- ps1: Sod, $model, two scalars against none"
            make_case "ps1-$model-base" --model "$model" --scalars 0 ; base="$W"
            make_case "ps1-$model-ps" --model "$model" --scalars 2 ; ps="$W"
            if run "$base" && run "$ps"; then
               "$VENV_PY" "$CASE_DIR/ps_oracle.py" "$base" "$ps" --scalars 2 || FAILED=1
            else
               echo "   a run failed: see $base/log.txt, $ps/log.txt  FAIL" ; FAILED=1
            fi
         done ;;
      ps2)
         if [[ $FNL -eq 1 ]]; then
            echo "-- ps2: passive scalars on FNL"
            make_case ps2-fnl --scalars 1
            expect_refused "$W" "passive_scalars\) runs on the CPU backend only"
         else
            echo "-- ps2: passive_scalars = 5"
            make_case ps2-five --scalars 5
            expect_refused "$W" "passive_scalars\) must lie in 0\.\.4"
            echo "-- ps2: passive_scalars with mhd-ideal"
            make_case ps2-mhd --model mhd-ideal --scalars 1
            expect_refused "$W" "passive_scalars\) needs a hydrodynamic physical_model"
         fi ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "PS FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "PS PASSED ($TAG, legs: $LEGS)"
