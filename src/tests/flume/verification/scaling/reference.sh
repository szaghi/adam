#!/usr/bin/env bash
# FLUME reference layer verification NV-5 (issue #49, N2).
#
# Why: a dimensional input with a [reference] section must run as its hand-normalised twin. scaling.py dimensionalize
# writes a case in the units L0 = 2^j, u0 = 2^k, rho0 = 4^m and adds the [reference] section that converts it back;
# the references are powers of two, so the reference layer (adam_flume_reference_object) hands the solver the base
# numbers exactly and the run must equal the base run BIT FOR BIT: every conservative field of the last checkpoint
# (scaling.py compare --conservative; the temperature changes unit, R* = 1) and every history file (.dat, byte for
# byte). Any option the layer misclassifies, or leaves unconverted, shows up as a mismatch; scaling.py classifies the
# options independently (its table mirrors the Fortran one) and refuses an unknown one. The conversions the layer logs
# are checked too (scaling.py check-log): every converted value equals the base value exactly and every dimensional
# option was converted, which covers the options whose conversion does not show in the solution.
#
# Cases: every single-realm regression case (src/tests/flume/regression), which between them exercise region states
# and boxes, inflow boundaries, the immersed solid, geometric AMR markers, the GLM/EGLM options and floors, the
# initial-condition parameters; the refused ones (orszag-tang, forest manifests) are listed as skipped. --inputs runs
# any input files instead (each as input.ini in its own work directory), e.g. the generated verification inputs that
# exercise the initial conditions, gradient AMR markers and slices no regression case uses.
#
# Usage: ./reference.sh [--np N] [--j J --k K --m M] [--cases "sod-x rotor ..." | --inputs "a.ini b.ini ..."]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./reference.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
REG_DIR="$REPO_ROOT/src/tests/flume/regression"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TOOL="$CASE_DIR/scaling.py"
NP=2 ; J=2 ; K=-1 ; M=-2 ; INPUTS=""
CASES="$(cd "$REG_DIR" && for d in */; do [[ -f "$d/input.ini" ]] && echo "${d%/}"; done)"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --j)     J="$2" ; shift 2 ;;
      --k)     K="$2" ; shift 2 ;;
      --m)     M="$2" ; shift 2 ;;
      --cases) CASES="$2" ; shift 2 ;;
      --inputs) INPUTS="$2" ; shift 2 ;;
      *) echo "reference.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
[[ -x "$EXE" ]] || { echo "reference.sh: executable '$EXE' not found" >&2 ; exit 2 ; }

TAG="$(basename "$EXE")-np$NP"
fails=0 ; passed=0 ; skipped=0

run() { # run <work-dir>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "reference.sh: run failed, see $1/log.txt" >&2
      exit 1
   fi
}

copy_case() { # copy_case <case or input file> <work-dir>: the case inputs, without goldens and outputs
   rm -rf "$2" ; mkdir -p "$2"
   if [[ -f "$1" ]]; then
      cp "$1" "$2/input.ini"
   else
      find "$REG_DIR/$1" -maxdepth 1 -type f ! -name '*.h5' ! -name '*.xdmf' ! -name '*.dat' ! -name '*.log' \
           -exec cp {} "$2" \;
   fi
}

echo ">> reference layer NV-5 ($(basename "$EXE"), np $NP): L0 = 2^$J, u0 = 2^$K, rho0 = 4^$M"
[[ -n $INPUTS ]] && CASES="$INPUTS"
for item in $CASES; do
   c="$(basename "$item" .ini)"
   base="$CASE_DIR/work-ref-$TAG-$c-base"
   ref="$CASE_DIR/work-ref-$TAG-$c-reference"
   copy_case "$item" "$base"
   copy_case "$item" "$ref"
   rc=0
   "$VENV_PY" "$TOOL" dimensionalize "$base/input.ini" "$ref/input.ini" --j "$J" --k "$K" --m "$M" \
      > "$ref/dimensionalize.txt" 2>&1 || rc=$?
   if [[ $rc -eq 3 ]]; then
      printf '   %-22s skipped (refused by the reference layer)\n' "$c"
      skipped=$(( skipped + 1 ))
      rm -rf "$base" "$ref"
      continue
   elif [[ $rc -ne 0 ]]; then
      printf '   %-22s FAIL: %s\n' "$c" "$(tail -1 "$ref/dimensionalize.txt")"
      fails=$(( fails + 1 ))
      continue
   fi
   run "$base"
   run "$ref"
   ngc="$(awk -F= '/^\[/{s=$0} s=="[grid]" && $1~/^ *ngc *$/{print $2+0}' "$base/input.ini")"
   status=0
   result="$("$VENV_PY" "$TOOL" compare "$base" "$ref" --ngc "${ngc:-3}" --conservative --psi glm)" || status=1
   for h in "$base"/*.dat "$base"/*.mat; do
      [[ -f "$h" ]] || continue
      cmp -s "$h" "$ref/$(basename "$h")" || { status=1 ; result="$result; $(basename "$h") differs" ; }
   done
   if ! log_check="$("$VENV_PY" "$TOOL" check-log "$base/input.ini" "$ref/log.txt")"; then
      status=1 ; result="$result; conversions: $(echo "$log_check" | grep -v 'conversions:' | tr -s ' ' | tr '\n' ';')"
   else
      result="$result; $(echo "$log_check" | tail -1 | sed 's/^ *//')"
   fi
   printf '   %-22s %s\n' "$c" "$(echo "$result" | sed 's/^ *//')"
   if [[ $status -eq 0 ]]; then passed=$(( passed + 1 )) ; else fails=$(( fails + 1 )) ; fi
   find "$base" "$ref" -name '*.h5' -delete
done
if [[ $fails -gt 0 ]]; then
   echo "reference layer NV-5 FAILED: $fails cases not bitwise, $passed bitwise, $skipped skipped ($TAG)"
   exit 1
fi
echo "reference layer NV-5: $passed cases bitwise, $skipped skipped ($TAG)"
