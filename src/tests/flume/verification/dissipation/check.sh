#!/usr/bin/env bash
# FLUME verification DC (issue #65, P1, P2): the input contract of the dissipative terms and the no-slip walls.
#
# Why: every coefficient of M4 is given either as a coefficient or as its dimensionless number (issue #49, N4), the
# positivity limiter is refused with any of them (D-M4-5), the no-slip walls take a tangential wall velocity and the
# isothermal one a temperature (D-M4-7). A wrong combination must stop the run with a message naming it, never run
# with a coefficient silently dropped. make_contract.py lists the cases; for each:
#   refuse          the run stops with the expected message;
#   log             the run logs the coefficients its keys imply (mu = 1/Re, k = mu cp/Pr, eta = 1/Rm), exactly, then
#                   completes and logs its diffusive dt limit (Euler, P2) or stops at the P3 guard (MHD);
#   ideal           zero coefficients and dissipative_order = 2 leave the run bit for bit equal to the base run;
#   convert         the case dimensionalised by ../scaling/scaling.py ([reference] with powers of two) logs every
#                   dimensional key converted back exactly (scaling.py check-log: coefficients, temperatures, wall
#                   velocity), then completes (Euler) or stops at the guard (MHD);
#   convert-refuse  the dimensionalised case is refused (lundquist without the Alfvenic preset).
# Base inputs: the regression sod-x (Euler) and uniform-amr-mhd (MHD) inputs.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
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
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TOOL="$CASE_DIR/make_contract.py"
SCALING="$CASE_DIR/../scaling/scaling.py"
REGRESSION="$REPO_ROOT/src/tests/flume/regression"
TAG="$(basename "$EXE")-np$NP"
FAILED=0
GUARD_MHD="are not yet computed with [physics].(physical_model)=mhd-ideal"
DT_LOG="dissipative dt limit (issue #65)"

base_of() { [[ $1 == euler ]] && echo "$REGRESSION/sod-x/input.ini" || echo "$REGRESSION/uniform-amr-mhd/input.ini" ; }
run() { (cd "$1" && mpirun -np "$NP" "$EXE" input.ini < /dev/null > log.txt 2>&1) ; } # mpirun reads stdin
runs_through() { # runs_through <work> <model>: Euler completes and logs its dt limit, MHD stops at the P3 guard
   if [[ $2 == euler ]]; then
      run "$1" && grep -aqF -- "$DT_LOG" "$1/log.txt"
   else
      ! run "$1" && grep -aqF -- "$GUARD_MHD" "$1/log.txt"
   fi
}
verdict() { # verdict <name> <ok: 0 pass> <note>
   if [[ $2 -eq 0 ]]; then
      printf '   %-26s PASS\n' "$1"
   else
      printf '   %-26s FAIL%s\n' "$1" "${3:+: $3}"
      FAILED=1
   fi
}
same_fields() { # same_fields <work a> <work b>: every dataset of the last checkpoint bit for bit equal
   "$VENV_PY" - "$1" "$2" <<'EOF'
import sys
from pathlib import Path
import h5py
import numpy as np
a, b = Path(sys.argv[1]), Path(sys.argv[2])
fa = sorted(p.name for p in a.glob("*-proc*.h5") if "-restart-" not in p.name)
last = max(int(n.split("-")[-2]) for n in fa)
names = [n for n in fa if int(n.split("-")[-2]) == last]
bad = 0
for n in names:
    with h5py.File(a / n) as ha, h5py.File(b / n) as hb:
        bad += sum(int(not np.array_equal(ha[k][()], hb[k][()])) for k in ha)
sys.exit(int(bad > 0 or not names))
EOF
}

echo ">> DC: input contract of the dissipative terms ($TAG)"
base_done=""
while IFS=$'\t' read -r name kind model expected; do
   w="$CASE_DIR/work-$TAG-$name"
   rm -rf "$w" ; mkdir -p "$w"
   base="$(base_of "$model")"
   "$VENV_PY" "$TOOL" "$base" "$w/input.ini" "$name"
   ok=0
   case "$kind" in
      refuse)
         run "$w" && ok=1
         grep -aqF -- "$expected" "$w/log.txt" || ok=1
         verdict "$name" "$ok" "expected the refusal '$expected', see $w/log.txt" ;;
      log)
         runs_through "$w" "$model" || ok=1
         "$VENV_PY" "$TOOL" --check-log "$w/log.txt" "$base" "$name" || ok=1
         verdict "$name" "$ok" "see $w/log.txt" ;;
      ideal)
         b="$CASE_DIR/work-$TAG-base-$model"
         if [[ " $base_done " != *" $model "* ]]; then
            rm -rf "$b" ; mkdir -p "$b" ; cp "$base" "$b/input.ini"
            if ! run "$b"; then verdict "base-$model" 1 "the base run failed, see $b/log.txt" ; continue ; fi
            base_done="$base_done $model"
         fi
         run "$w" || ok=1
         if [[ $ok -eq 0 ]]; then same_fields "$b" "$w" || ok=1 ; fi
         verdict "$name" "$ok" "not bit for bit equal to the base run, see $w" ;;
      convert|convert-refuse)
         mv "$w/input.ini" "$w/base.ini"
         "$VENV_PY" "$SCALING" dimensionalize "$w/base.ini" "$w/input.ini" --j 1 --k -1 --m 1 > /dev/null
         if [[ $kind == convert ]]; then
            runs_through "$w" "$model" || ok=1
            "$VENV_PY" "$SCALING" check-log "$w/base.ini" "$w/log.txt" | sed 's/^/   /' || ok=1
         else
            run "$w" && ok=1
            grep -aqF -- "$expected" "$w/log.txt" || ok=1
         fi
         verdict "$name" "$ok" "see $w/log.txt" ;;
   esac
   find "$w" -name '*.h5' -delete
done < <("$VENV_PY" "$TOOL" --list)

if [[ $FAILED -eq 0 ]]; then
   echo "DC PASSED ($TAG)"
else
   echo "DC FAILED ($TAG)"
   exit 1
fi
