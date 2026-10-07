#!/usr/bin/env bash
# FLUME verification EV-step (issue #46): the Woodward-Colella Mach 3 forward-facing step, a three-realm forest with
# quadtree AMR.
#
# Why: the step tunnel is L-shaped, so it is three rectangular realms glued by 1:1 mirror seams with the step faces as
# physical walls (make_step.py), and its 2-D flow runs on quadtrees (ratio 4, nk = 1) since issue #46. The case puts
# together inter-realm mirror seams (two, on different axes), intra-realm quadtree 2:1 seams (refined boxes off the
# inter-realm seams) and the walls, inflow and outflow of a supersonic tunnel. Legs:
#   trees  N = 40 (2x2 blocks per realm, one refined level on two boxes), t = 0.5: the forest on quadtrees (nk = 1)
#          against the same forest on octrees with a null z axis (nk = 4), quadtree_oracle.py (../mhd/quadtree): every
#          run z-invariant, the quadtree within TOL of the octree at the same (x, y), momentum scaled as one vector.
#          Both runs must also keep density and pressure positive and finite;
#   full   N = 80, the coarse grid of Woodward and Colella, quadtree, t = 4 (their reference time), with the refined
#          boxes and, as the control of the figure's zoom, on the uniform 1/80 grid (no box): both positive and finite;
#          the checkpoints are kept (--keep) for the documentation figures (make_doc_figures.py --only step step-trees).
#          Measured (CPU, np 2): 7392 steps, min density 0.356, min pressure 0.393 (refined). About 30 + 11 min on the
#          CPU. Not a default leg.
#
# Measured (trees, np 2, CPU): quadtree against octree 6.5e-13, octree z spread 7.6e-13, quadtree 0; 4 min.
#
# Usage: ./check.sh [--np N] [--legs trees,full] [--keep]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEGS="trees"
KEEP=0
TOL="1.0e-10"    # quadtree against octree, relative (measured 6.5e-13)
Z_TOL="1.0e-10"  # z invariance, relative (measured: octree 7.6e-13, quadtree 0)

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) LEGS="$2" ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs trees,full, --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

has_leg() { [[ ",$LEGS," == *",$1,"* ]]; }
run() { # run <work> <make_step.py arguments>: write the forest, run it
   local w=$1
   shift
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_step.py" "$w" "$@"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" step.ini > log.txt 2>&1); then
      echo "check.sh: $(basename "$w") failed, see $w/log.txt" >&2
      FAILED=1
      return 1
   fi
}
positive() { # positive <work>: density and pressure of the last checkpoint positive and finite in every realm
   "$VENV_PY" - "$1" <<'EOF'
import sys
from pathlib import Path

import h5py
import numpy as np

work, gamma = Path(sys.argv[1]), 1.4
files = [p for p in work.glob("*-proc*.h5") if "-restart-" not in p.name]
last = max(int(p.name.split("-")[-2]) for p in files)
rho_min = p_min = np.inf
bad = 0
for path in (p for p in files if int(p.name.split("-")[-2]) == last):
    with h5py.File(path, "r") as h5:
        for blk in {k.rsplit("-", 1)[0] for k in h5}:
            q = {v: h5[f"{blk}-{v}"][()][3:-3, 3:-3, 3:-3] for v in ("r", "ru", "rv", "rw", "rE")}
            p = (gamma - 1.0) * (q["rE"] - 0.5 * (q["ru"] ** 2 + q["rv"] ** 2 + q["rw"] ** 2) / q["r"])
            bad += int(np.count_nonzero(~np.isfinite(p)) + np.count_nonzero(~np.isfinite(q["r"])))
            rho_min, p_min = min(rho_min, float(q["r"].min())), min(p_min, float(p.min()))
ok = bad == 0 and rho_min > 0.0 and p_min > 0.0
print(f"   {work.name}: step {last}, min density {rho_min:.3e}, min pressure {p_min:.3e}, {bad} non-finite"
      f"  {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
}

if has_leg trees; then
   echo ">> EV-step trees: 3-realm forest, N = 40 + refined boxes, t = 0.5, quadtree vs octree ($TAG)"
   works=()
   for tree in oct quad; do
      w="$CASE_DIR/work-$TAG-trees-$tree"
      run "$w" --tree "$tree" --cells 40 --refine --time-max 0.5 && works+=("$w")
   done
   if [[ ${#works[@]} -eq 2 ]]; then
      for w in "${works[@]}"; do positive "$w" || FAILED=1; done
      "$VENV_PY" "$VERIF_DIR/mhd/quadtree/quadtree_oracle.py" "${works[@]}" --ngc 3 --tol "$TOL" --z-tol "$Z_TOL" \
         --groups ru,rv,rw | sed 's/^/   /' || FAILED=1
   fi
   [[ $KEEP -eq 1 ]] || for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done
fi

if has_leg full; then
   echo ">> EV-step full: 3-realm forest, N = 80, quadtree, t = 4, refined boxes and uniform control"
   for variant in refined uniform; do
      w="$CASE_DIR/work-$TAG-full" ; opts=(--refine)
      [[ $variant == uniform ]] && { w="$w-uniform" ; opts=() ; }
      if run "$w" --tree quad --cells 80 --time-max 4.0 "${opts[@]}"; then
         positive "$w" || FAILED=1
      fi
      [[ $KEEP -eq 1 ]] || find "$w" -name '*.h5' -delete
   done
fi

if [[ $FAILED -eq 0 ]]; then
   echo "EV-step PASSED ($TAG)"
else
   echo "EV-step FAILED ($TAG)"
   exit 1
fi
