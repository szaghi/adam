#!/usr/bin/env bash
# PRISM verification of the [slices] output (slices_restauration.md, F1/F5).
#
# Why: the slices are outputs no run reads back, so a defect in them is silent. slices.ini is the source-free gaussian
# pulse of rmf-amr-fd-pulse (2:1 AMR seam, 36 blocks) with curl and divergence saved and three slices (a trilinear
# line through the seam, an inverse-distance plane, a trilinear box around the block corner at the origin). For every
# process count:
#   1. slices_oracle.py recomputes every slice value of q, curl and div from the XH5F checkpoint of the slice step
#      (saved with ghost cells) and checks the intra-realm ghost cells of every field set against the same-level
#      neighbour interior (the auxiliary fields are computed on the interior cells only; save_slices refreshes them);
#   2. the slice files are compared, bitwise, with those of the first process count: each point is written by the one
#      rank owning it, so the output must not depend on the decomposition.
#
# Usage: ./check.sh [--build] [--np "1 2 4"]
#
# PRISM_EXE overrides the executable under test, e.g. PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
NPS="1 2 4"
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NPS="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np \"N ...\")" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "check.sh: executable $EXE not found" >&2 ; exit 2 ; }

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c "import h5py, numpy" 2>/dev/null; then
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")"
first=""
status=0
for np in $NPS; do
   work="$CASE_DIR/work-$TAG-np$np"
   rm -rf "$work" ; mkdir -p "$work"
   cp "$CASE_DIR/slices.ini" "$work/"
   echo ">> $TAG, np $np"
   if ! (cd "$work" && mpirun -np "$np" "$EXE" slices.ini > log.txt 2>&1); then
      echo "   run failed, see $work/log.txt" ; status=1 ; continue
   fi
   "$VENV_PY" -I "$CASE_DIR/slices_oracle.py" "$work" || status=1
   if [[ -z "$first" ]]; then
      first="$work"
   else
      for mat in "$first"/*.mat; do
         if cmp -s "$mat" "$work/$(basename "$mat")"; then
            echo "   $(basename "$mat"): identical to np $(basename "$first" | sed 's/.*-np//')  PASS (bitwise)"
         else
            echo "   $(basename "$mat"): differs from np $(basename "$first" | sed 's/.*-np//')  FAIL" ; status=1
         fi
      done
   fi
done
if [[ $status -eq 0 ]]; then echo "slices verification PASSED ($TAG)"; else echo "slices verification FAILED ($TAG)"; fi
exit $status
