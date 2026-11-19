#!/usr/bin/env bash
# A/B check of a refactoring (SELAH M7, ADAM #132): one regression case run by two executables, compared bit for bit.
#
# Why this script exists: run.sh compares a run with its goldens within rtol = 1e-6, so it cannot tell an identical
# result from one that moved at round-off. The A/B runs the same case with executable A (e.g. the develop build) and
# executable B (the refactored build) and requires every history log byte for byte and every checkpoint dataset bit for
# bit (regression_ab.py); on a difference it prints the first divergent row and the largest deviation.
#
# --scheme overrides `[runge_kutta] scheme` in the run directory's copies of the case inputs (every INI with a
# [runge_kutta] section; a realm INI that extends another inherits it): the old integrators no committed case runs,
# without committing an input. Nothing is written in the case directory except work-ab-*/.
#
# Usage:
#   regression-ab.sh --case DIR --exe-a EXE --exe-b EXE [--scheme NAME] [--np N] [--work DIR] [--timeout SECONDS]
# Environment: PYTHON (default exe/.regression-venv/bin/python, with h5py and numpy), ADAM_MEMORY_FRACTION and
# ADAM_UNUSED_KEYS (defaults as run.sh), MPIRUN_ARGS (extra mpirun arguments).
# Exit status: 0 bitwise, 1 different, 2 usage error, 3 a run did not complete.

set -uo pipefail

CASE="" ; EXE_A="" ; EXE_B="" ; SCHEME="" ; NP=2 ; WORK="" ; TIMEOUT=900
while [[ $# -gt 0 ]]; do
   case "$1" in
      --case)    CASE="$2" ; shift ;;
      --exe-a)   EXE_A="$2" ; shift ;;
      --exe-b)   EXE_B="$2" ; shift ;;
      --scheme)  SCHEME="$2" ; shift ;;
      --np)      NP="$2" ; shift ;;
      --work)    WORK="$2" ; shift ;;
      --timeout) TIMEOUT="$2" ; shift ;;
      *) echo "ERROR: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
   shift
done
[[ -n "$CASE" && -n "$EXE_A" && -n "$EXE_B" ]] || { sed -n 2,17p "$0" >&2 ; exit 2 ; }
[[ -f "$CASE/input.ini" ]] || { echo "ERROR: $CASE/input.ini not found" >&2 ; exit 2 ; }
for e in "$EXE_A" "$EXE_B"; do [[ -x "$e" ]] || { echo "ERROR: executable $e not found" >&2 ; exit 2 ; } ; done
CASE="$(cd "$CASE" && pwd)"
EXE_A="$(cd "$(dirname "$EXE_A")" && pwd)/$(basename "$EXE_A")"
EXE_B="$(cd "$(dirname "$EXE_B")" && pwd)/$(basename "$EXE_B")"
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PYTHON="${PYTHON:-$REPO_ROOT/exe/.regression-venv/bin/python}"
[[ -x "$PYTHON" ]] || { echo "ERROR: no Python with h5py at $PYTHON (set PYTHON)" >&2 ; exit 2 ; }
label="${SCHEME:-as-committed}"
WORK="${WORK:-$CASE/work-ab-$label}"
export ADAM_UNUSED_KEYS="${ADAM_UNUSED_KEYS:-error}"
if [[ -z "${ADAM_MEMORY_FRACTION:-}" ]]; then # the cap of run.sh: 5 %, at least 2.5 GB of the node
   ADAM_MEMORY_FRACTION=$(awk '/^MemTotal:/ {f = 2.5e6 / $2; if (f < 0.05) f = 0.05; if (f > 0.85) f = 0.85;
                                            printf "%.4f", f; found = 1} END {if (!found) printf "0.05"}' /proc/meminfo 2>/dev/null)
fi
export ADAM_MEMORY_FRACTION

prepare() { # $1 = run directory: the case inputs, the scheme overlay
   rm -rf "$1" && mkdir -p "$1" && cp "$CASE"/*.ini "$1/" || return 1
   [[ -z "$SCHEME" ]] && return 0
   local f n=0
   for f in "$1"/*.ini; do
      grep -qiE '^[[:space:]]*\[runge_kutta\]' "$f" || continue
      awk -v s="$SCHEME" '
         /^[[:space:]]*\[/ { sec = tolower($0); gsub(/[[:space:]]/, "", sec) }
         sec == "[runge_kutta]" && /^[[:space:]]*scheme[[:space:]]*=/ { print "scheme = " s " ; A/B overlay (regression-ab.sh)"; next }
         { print }' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
      n=$((n + 1))
   done
   [[ $n -gt 0 ]] || { echo "ERROR: no [runge_kutta] section in the inputs of $CASE: --scheme has nothing to override" >&2 ; return 1 ; }
}

run() { # $1 = run directory, $2 = executable
   local t0=$SECONDS rc
   ( cd "$1" && timeout --signal=TERM --kill-after=30s "$TIMEOUT" mpirun -np "$NP" ${MPIRUN_ARGS:-} "$2" input.ini > run.log 2>&1 )
   rc=$?
   if ! grep -q 'ADAM run complete' "$1/run.log"; then
      echo "AB $(basename "$CASE") $label: RUN FAILED in $1 (mpirun rc $rc, no 'ADAM run complete'); tail of run.log:"
      tail -n 8 "$1/run.log" | sed 's/^/   /'
      return 3
   fi
   echo "   $(basename "$1"): $(( SECONDS - t0 )) s, mpirun rc $rc"
}

prepare "$WORK/a" || exit 2
prepare "$WORK/b" || exit 2
echo "AB $(basename "$CASE") $label: A = $EXE_A, B = $EXE_B, np $NP"
run "$WORK/a" "$EXE_A" || exit 3
run "$WORK/b" "$EXE_B" || exit 3
"$PYTHON" "$REPO_ROOT/src/tests/regression_ab.py" "$WORK/a" "$WORK/b" > "$WORK/compare.txt" 2>&1
rc=$?
sed 's/^/   /' "$WORK/compare.txt"
case $rc in
   0) echo "AB $(basename "$CASE") $label: BITWISE" ;;
   1) echo "AB $(basename "$CASE") $label: DIFFERENT" ;;
   *) echo "AB $(basename "$CASE") $label: COMPARE ERROR (rc $rc)" ; exit 2 ;;
esac
exit $rc
