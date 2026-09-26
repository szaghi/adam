#!/usr/bin/env bash
# FLUME MHD verification, M2-P3 exit (issue #41): Brio-Wu and RJ4d run clean, and the positivity floors work.
#
# Why: Brio-Wu (compound wave, strong field) and RJ4d (switch-on/off waves) stress the MHD eigensystem near its
# degeneracies; a scheme that needs the floors on them is not clean. The floors themselves (issue #41, section 3.8)
# are exercised on a state with a negative pressure, which no clean solver produces. Legs (x, N = 512, np 2):
#   1. Brio-Wu (no cleaning and GLM) and RJ4d, floors armed at 1e-10: the run reaches t_end, no stage floors a cell
#      (no "MHD floors:" line in the log), and the final density and pressure are positive and finite;
#   2. the negative-pressure state with p_floor = 1e-3: the run completes and the log reports floored cells;
#   3. the same state with the floors disabled (0): the run stops on the non-positive state with its message.
# One run at a time; checkpoints deleted after use.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
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
TAG="$(basename "$EXE")-np$NP"
FAILED=0

prepare() { # prepare <label> <make_riemann args...>: fresh work dir with input.ini, prints the work dir
   local w="$CASE_DIR/work-$TAG-$1" ; shift
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_riemann.py" "$VERIF_DIR/sod/sod-x.ini" "$w/input.ini" "$@"
   echo "$w"
}

positive() { # positive <work>: final density and pressure positive and finite (auxiliary fields, interior cells)
   "$VENV_PY" - "$1" <<'PYEOF'
import sys
from pathlib import Path
import h5py
import numpy as np
work = Path(sys.argv[1])
files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
last = max(int(p.name.split("-")[-2]) for p in files)
lo = {"rho": np.inf, "p": np.inf}
finite = True
for path in (p for p in files if int(p.name.split("-")[-2]) == last):
    with h5py.File(path, "r") as h5:
        for key in h5:
            name = key.rsplit("-", 1)[1]
            if name in lo:
                a = h5[key][()][3:-3, 3:-3, 3:-3]
                finite &= bool(np.all(np.isfinite(a)))
                lo[name] = min(lo[name], float(a.min()))
ok = finite and lo["rho"] > 0.0 and lo["p"] > 0.0
print(f"   {work.name}: final min rho {lo['rho']:.4e}, min p {lo['p']:.4e}, finite {finite}  {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
}

for spec in "bw none" "bw glm" "rj4d none"; do
   set -- $spec
   echo ">> $1 ($2), N = 512, floors armed at 1e-10 ($(basename "$EXE"), np $NP)"
   w="$(prepare "$1-$2" --case "$1" --cells 512 --divergence-control "$2" --rho-floor 1e-10 --p-floor 1e-10)"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "   $1 ($2): run failed, see $w/log.txt  FAIL" ; FAILED=1 ; continue
   fi
   n=$(grep -c "MHD floors:" "$w/log.txt" || true)
   if [[ "$n" -eq 0 ]]; then echo "   $1 ($2): no floored cell in any stage  PASS"
   else echo "   $1 ($2): $n stages floored cells  FAIL" ; FAILED=1 ; fi
   positive "$w" || FAILED=1
   find "$w" -name '*.h5' -delete
done

echo ">> negative pressure, p_floor = 1e-3: the floors fix the state"
w="$(prepare negative-floored --case negative --cells 512 --divergence-control none --p-floor 1e-3)"
if (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) && grep -q "MHD floors:" "$w/log.txt"; then
   echo "   floors engaged: $(grep -m1 "MHD floors:" "$w/log.txt" | sed 's/^.*MHD floors: //')  PASS"
else
   echo "   floors not engaged or run failed, see $w/log.txt  FAIL" ; FAILED=1
fi
find "$w" -name '*.h5' -delete

echo ">> negative pressure, floors disabled: the run stops on the non-positive state"
w="$(prepare negative-stop --case negative --cells 512 --divergence-control none)"
(cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || true
if grep -q "cells with a non-positive density or pressure" "$w/log.txt"; then
   echo "   stops with: $(grep -m1 "non-positive density" "$w/log.txt" | sed 's/^.*: \([0-9]* cells\)/\1/' | cut -c1-110)  PASS"
else
   echo "   the expected stop message was not found in $w/log.txt  FAIL" ; FAILED=1
fi
find "$w" -name '*.h5' -delete

if [[ $FAILED -eq 0 ]]; then
   echo "MHD Riemann/floors verification PASSED ($TAG)"
else
   echo "MHD Riemann/floors verification FAILED ($TAG)"
   exit 1
fi
