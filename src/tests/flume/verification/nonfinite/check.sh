#!/usr/bin/env bash
# FLUME verification: the non-finite state guard (issue #45), a positive control.
#
# Why: a run whose state goes NaN used to finish normally and exit 0 (the sod-amr event of issue #45 was caught only by
# the regression digest). Every step now counts the non-finite values of the committed state and stops the run,
# locating the first one on every rank. This check plants a NaN in the right state of the Sod tube of V1 (region 2
# pressure = nan) and requires the run to fail with the guard's message; a guard that never fires would pass every
# other test unnoticed.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
TAG="$(basename "$EXE")-np$NP"

echo ">> non-finite guard: Sod with a NaN pressure in region 2 ($(basename "$EXE"), np $NP)"
w="$CASE_DIR/work-$TAG"
rm -rf "$w" ; mkdir -p "$w"
awk '/^\[initial_conditions_region_2\]/{r2=1} /^\[/{if ($0 !~ /region_2/) r2=0} r2 && /^p[[:space:]]*=/{print "p      = nan"; next} {print}' \
   "$VERIF_DIR/sod/sod-x.ini" > "$w/sod-nan.ini"
if ! grep -q "^p      = nan" "$w/sod-nan.ini"; then
   echo "check.sh: could not plant the NaN in $w/sod-nan.ini" >&2
   exit 2
fi
if (cd "$w" && timeout 600 mpirun -np "$NP" "$EXE" sod-nan.ini > log.txt 2>&1); then
   echo "non-finite guard FAILED ($TAG): the run with a NaN state exited 0, see $w/log.txt"
   exit 1
fi
if ! grep -a -q "non-finite (NaN or infinite) values in the state at step" "$w/log.txt"; then
   echo "non-finite guard FAILED ($TAG): the run stopped, but not on the guard, see $w/log.txt"
   exit 1
fi
grep -a -m2 "non-finite (NaN or infinite)" "$w/log.txt"
find "$w" -name '*.h5' -delete
echo "non-finite guard PASSED ($TAG)"
