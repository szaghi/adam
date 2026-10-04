#!/usr/bin/env bash
# FLUME scaling-covariance verification NV-0 / NV-4 (issue #49, N0).
#
# Why: ideal Euler and MHD in FLUME's units carry no dimensionless number, so an input rescaled by powers of two
# (lengths 2^j, velocities 2^k, density 4^m) is the same problem in other units, and a code without absolute constants
# reproduces the base run BIT FOR BIT once the outputs are scaled back (scaling.py: the scales are exact powers of two).
# Each case runs at its base scale and at four rescalings, each compared with the base:
#   length only (j, 0, 0), velocity only (0, k, 0), density only (0, 0, m), all together (j, k, m),
# so a mismatch is attributed to the quantity whose magnitude moved. Cases: Sod (sod-x, Euler, fast path), the
# isentropic vortex (periodic, Euler), RJ2a (MHD, no cleaning), the Balsara-Spicer blast with the positivity limiter
# (EGLM, splitting scheme). Orszag-Tang reads no initial-condition parameter and cannot be rescaled from its input: it
# is listed as such.
#
# Default: the BASELINE (NV-0) is printed and the script exits 0 whatever it finds; #49 expects the absolute WENO
# regulariser (zeps = 1e-6) and the absolute positivity floor to break covariance until N1. --expect-bitwise turns
# it into the gate NV-4: any non-bitwise comparison fails.
#
# Usage: ./check.sh [--np N] [--j J --k K --m M] [--expect-bitwise] [--cases "sod-x vortex rj2a blast orszag-tang"]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TOOL="$CASE_DIR/scaling.py"
NP=2 ; J=2 ; K=-1 ; M=-2 ; GATE=0
CASES="sod-x vortex rj2a blast orszag-tang"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)             NP="$2" ; shift 2 ;;
      --j)              J="$2" ; shift 2 ;;
      --k)              K="$2" ; shift 2 ;;
      --m)              M="$2" ; shift 2 ;;
      --cases)          CASES="$2" ; shift 2 ;;
      --expect-bitwise) GATE=1 ; shift ;;
      *) echo "check.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
[[ -x "$EXE" ]] || { echo "check.sh: executable '$EXE' not found" >&2 ; exit 2 ; }

declare -A INPUT=([sod-x]="$REPO_ROOT/src/tests/flume/verification/sod/sod-x.ini"
                  [vortex]="$REPO_ROOT/src/tests/flume/regression/vortex-periodic/input.ini"
                  [rj2a]="$REPO_ROOT/src/tests/flume/regression/rj2a-x/input.ini"
                  [blast]="$REPO_ROOT/src/tests/flume/regression/blast-limiter/input.ini"
                  [orszag-tang]="$REPO_ROOT/src/tests/flume/regression/orszag-tang/input.ini")
declare -A PSI=([blast]=eglm)
TAG="$(basename "$EXE")-np$NP"
fails=0

run() { # run <work-dir> <ini-name>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: run failed, see $1/log.txt" >&2
      exit 1
   fi
}

echo ">> scaling covariance ($(basename "$EXE"), np $NP): j = $J, k = $K, m = $M"
for c in $CASES; do
   ini="${INPUT[$c]}" ; name="$(basename "$ini")"
   base="$CASE_DIR/work-$TAG-$c-base"
   rm -rf "$base" ; mkdir -p "$base"
   if ! "$VENV_PY" "$TOOL" rescale "$ini" "$base/$name" 2> "$base/rescale.txt"; then
      echo "== $c: not rescalable: $(tail -1 "$base/rescale.txt")"
      continue
   fi
   cp "$ini" "$base/$name"
   echo "== $c"
   run "$base" "$name"
   for v in "length $J 0 0" "velocity 0 $K 0" "density 0 0 $M" "all $J $K $M"; do
      read -r label j k m <<< "$v"
      w="$CASE_DIR/work-$TAG-$c-$label"
      rm -rf "$w" ; mkdir -p "$w"
      "$VENV_PY" "$TOOL" rescale "$ini" "$w/$name" --j "$j" --k "$k" --m "$m"
      run "$w" "$name"
      printf '   %-8s (j %2s, k %2s, m %2s)' "$label" "$j" "$k" "$m"
      if ! "$VENV_PY" "$TOOL" compare "$base" "$w" --ngc 3 --j "$j" --k "$k" --m "$m" --psi "${PSI[$c]:-glm}" \
           | sed 's/^ */ /'; then
         fails=$(( fails + 1 ))
      fi
   done
done
for w in "$CASE_DIR"/work-"$TAG"-*; do find "$w" -name '*.h5' -delete; done
if [[ $GATE -eq 1 && $fails -gt 0 ]]; then
   echo "scaling covariance FAILED: $fails comparisons not bitwise ($TAG)"
   exit 1
fi
echo "scaling covariance: $fails comparisons not bitwise ($TAG)$( [[ $GATE -eq 0 ]] && echo ', baseline mode')"
