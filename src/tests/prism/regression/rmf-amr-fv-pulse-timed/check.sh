#!/usr/bin/env bash
# rmf-amr-fv-pulse-timed oracle (issue #38): the FV reflux keeps a time-driven AMR run conservative on its capped step.
#
# THE CLAIM UNDER TEST: on a time-driven run (it_max <= 0) PRISM caps its last step to land on time_max, and the
# Berger-Colella reflux correction must be scaled by that capped step, the one the update used. Before 4c1dd1a0 it was
# scaled by the forest's uncapped dt, and the coarse-fine interface was not conservative on that step.
#
# The case (input.ini) is source-free: a gaussian EM pulse crossing two 2:1 seam planes on the FV path, 4.5 steps of
# time, so the fifth step is half a step. With J = 0 the FV update with reflux keeps the domain integrals of D and B
# constant up to the fluxes through the domain boundary, which this geometry keeps small (see input.ini).
#
# ACCEPTANCE:
#   1. the run completes and its last time step is shorter than the first (the capped step happened);
#   2. the reflux diagnostic is printed (the seams are active);
#   3. at every saved step, the drift of the integral of each D component relative to the largest D integral, and of
#      each B component relative to the largest B integral, is at most DRIFT_TOL.
#
# Measured (CPU, np 2): int(Dz) drifts 1.8e-15, 1.4e-13, 1.8e-12, 1.0e-11 over the four full steps (the seam's
# disturbance reaching the x boundaries) and 2.2e-11 after the capped one; int(Bx) stays at 1e-16. NEGATIVE CONTROL
# (the fix reverted, the correction scaled by the forest dt): the same drift up to step 4, then 1.4e-7 on the capped
# step, 6500 times larger. DRIFT_TOL = 1e-9 sits 45x above the fixed run and 140x below the reverted one. Note that
# int(Bx) cannot see the defect (Bx has no x-flux, the reflux at the x-normal seams never touches it): int(Dz) does.
#
# Usage: ./check.sh [--np N]
#
# PRISM_EXE overrides the executable under test, e.g. PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH plus the WSL UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
DRIFT_TOL="1.0e-9"
NP=2

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set PRISM_EXE)" >&2
   exit 2
fi
W="$CASE_DIR/work-$(basename "$EXE")-np$NP"
rm -rf "$W" ; mkdir -p "$W"
cp "$CASE_DIR/input.ini" "$W/"

echo ">> rmf-amr-fv-pulse-timed: FV pulse across 2:1 seams, last step capped ($(basename "$EXE"), np $NP)"
STATUS=0
if ! (cd "$W" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
   echo "   run failed, see $W/log.txt  FAIL"
   exit 1
fi
mapfile -t DTS < <(grep -a '^\[mpi-00000\]time step:' "$W/log.txt" | awk '{print $NF}')
if "$VENV_PY" -c "import sys; sys.exit(0 if float('${DTS[-1]}') < 0.99 * float('${DTS[0]}') else 1)"; then
   echo "   ${#DTS[@]} steps, last step ${DTS[-1]} against ${DTS[0]} (capped)  PASS"
else
   echo "   ${#DTS[@]} steps, last step ${DTS[-1]} against ${DTS[0]}: the last step is not capped  FAIL"
   STATUS=1
fi
if grep -aq 'forest: reflux max|F_coarse-F_fine_sum|' "$W/log.txt"; then
   echo "   reflux diagnostic printed (seams active)  PASS"
else
   echo "   no reflux diagnostic: the run has no 2:1 seam  FAIL"
   STATUS=1
fi
"$VENV_PY" - "$W" "$DRIFT_TOL" <<'EOF' || STATUS=1
import glob
import sys
from collections import defaultdict

import h5py
import numpy as np

work, tol, ngc = sys.argv[1], float(sys.argv[2]), 3
groups = {"D": ("Dx", "Dy", "Dz"), "B": ("Bx", "By", "Bz")}
names = groups["D"] + groups["B"]
files = defaultdict(list)
for f in glob.glob(f"{work}/*-[0-9]*-proc*.h5"):
    if "-restart-" not in f:
        files[int(f.split("-")[-2])].append(f)
integrals = {}
for step, paths in sorted(files.items()):
    total = np.zeros(len(names))
    for path in paths:
        with h5py.File(path, "r") as h5:
            for blk in {k.rsplit("-", 1)[0] for k in h5}:
                dv = float(np.prod(h5[f"{blk}-dxdydz"][()]))
                for v, name in enumerate(names):
                    total[v] += float(np.sum(h5[f"{blk}-{name}"][()][ngc:-ngc, ngc:-ngc, ngc:-ngc])) * dv
    integrals[step] = total
steps = sorted(integrals)
if len(steps) < 2:
    sys.exit(f"   fewer than 2 checkpoints in {work}  FAIL")
first = integrals[steps[0]]
scale = np.zeros(len(names))
for group in groups.values():
    idx = [names.index(n) for n in group]
    s = max(float(np.max(np.abs(integrals[t][idx]))) for t in steps)
    scale[idx] = s if s > 0.0 else 1.0
worst = 0.0
for t in steps[1:]:
    drift = np.abs(integrals[t] - first) / scale
    worst = max(worst, float(np.max(drift)))
    print(f"   step {t}: drift " + " ".join(f"{n} {d:.1e}" for n, d in zip(names, drift, strict=True)))
ok = worst <= tol
print(f"   max drift {worst:.2e} (tolerance {tol:.0e})  {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
find "$W" -name '*.h5' -delete

if [[ $STATUS -eq 0 ]]; then
   echo "rmf-amr-fv-pulse-timed PASSED"
else
   echo "rmf-amr-fv-pulse-timed FAILED"
   exit 1
fi
