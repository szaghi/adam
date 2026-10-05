#!/usr/bin/env bash
# Run a FLUME executable with the scale-invariant WENO weights (issue #49, N1: NV-2, NV-3), whatever the input says.
#
# Why: the verification scripts (vortex, sod, riemann-flux, mhd/*) build their inputs themselves and honour FLUME_EXE.
# Pointed at this wrapper, every one of them runs with `[weno] weights = si` and judges the result with its own oracle
# and bounds, unchanged, so the orders and the shock errors of the new weights are measured by the same checks as the
# default ones. The wrapper writes `<name>.si.ini`, a copy of the input with the key set (each MPI rank writes the same
# file through an atomic rename, so the ranks do not race), and execs the executable on it.
#
# Usage: FLUME_EXE=$REPO/src/tests/flume/verification/scaling/weights-exe.sh [WEIGHTS_EXE=<exe>] [WEIGHTS=si] ./check.sh
#   WEIGHTS_EXE  the FLUME executable to run (default exe/adam_flume_cpu of the repository)
#   WEIGHTS      the value of [weno] weights (default si)
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../../../../.." && pwd)"
exe="${WEIGHTS_EXE:-$repo/exe/adam_flume_cpu}"
weights="${WEIGHTS:-si}"
args=()
for a in "$@"; do
   if [[ $a == *.ini && -f $a ]]; then
      out="${a%.ini}.si.ini"
      awk -v w="$weights" '
         /^\[/ { if (inweno && !done) { print "weights = " w; done = 1 } inweno = ($0 ~ /^\[weno\]/) }
         inweno && /^[ \t]*weights[ \t]*=/ { next }
         { print }
         END { if (inweno && !done) print "weights = " w }' "$a" > "$out.$$"
      grep -q '^\[weno\]' "$out.$$" || { echo "weights-exe.sh: $a has no [weno] section" >&2; exit 2; }
      mv -f "$out.$$" "$out"
      args+=("$out")
   else
      args+=("$a")
   fi
done
exec "$exe" "${args[@]}"
