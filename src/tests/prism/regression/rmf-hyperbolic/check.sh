#!/usr/bin/env bash
# rmf-hyperbolic: the divergence-correction refusals of issue #61 (manual legs; run.sh runs only the golden digest).
#
# The case itself (constrained_transport = DB, divergence_correction = hyperbolic) is the golden-digest anchor of the
# GLM layout: nv = 11, phi and psi in slots 7 and 8, J in 9..11. These legs check the two inputs PRISM must refuse at
# initialization, each with its own message, instead of running:
#   1. hyperbolic without constrained transport: no phi/psi slot exists, so var_J* and the eigensystem were left
#      undefined (q(var_Jx,...) indexed out of bounds);
#   2. poisson: the CT projection is not implemented on either backend (the FNL kernel is a stub).
# A crash also exits non-zero, so each leg requires the message in the log as well.
#
# Usage: PRISM_EXE=/path/to/adam_prism_cpu ./check.sh     (default: <repo>/exe/adam_prism_cpu)
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"

fail=0

refusal_leg() {
   # $1 leg name, $2 constrained_transport, $3 divergence_correction, $4 expected message
   local work="$CASE_DIR/work-refusal-${EXE##*_}-$1" rc=0 # per backend (cpu, fnl): the two may run at once
   rm -rf "$work" && mkdir -p "$work"
   sed -e "s/^constrained_transport *=.*/constrained_transport = $2/" \
       -e "s/^divergence_correction *=.*/divergence_correction = $3/" "$CASE_DIR/input.ini" > "$work/input.ini"
   ( cd "$work" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 ) || rc=$?
   find "$work" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' \) -delete
   if [[ $rc -eq 0 ]]; then
      echo "FAIL [rmf-hyperbolic] $1: constrained_transport = $2 with divergence_correction = $3 RAN instead of refusing"
      fail=1
   elif ! grep -q "$4" "$work/run.log"; then
      echo "FAIL [rmf-hyperbolic] $1: died (exit $rc) WITHOUT the refusal message '$4': a crash, not a refusal"
      fail=1
   else
      echo ">> [rmf-hyperbolic] $1: clean refusal (exit $rc)"
   fi
}

refusal_leg hyperbolic-no-ct no hyperbolic "hyperbolic needs"
refusal_leg poisson          B  poisson    "poisson is not implemented"

if [[ $fail -ne 0 ]]; then
   echo "FAIL [rmf-hyperbolic] refusal legs"
   exit 1
fi
echo "PASS [rmf-hyperbolic] refusal legs"
