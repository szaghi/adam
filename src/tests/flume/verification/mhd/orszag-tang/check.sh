#!/usr/bin/env bash
# FLUME MHD verification MV-12 (issue #41, M2-P6): the Orszag-Tang vortex, symmetry, conservation and positivity.
#
# Why: the standard 2-D nonlinear MHD test (Stone et al. 2008, section 8.4), shocks and a current sheet from a smooth
# start. 128^2, t = 0.5, GLM, np 2. Through orszag_tang_oracle.py:
#   - the 180 degrees rotational symmetry of the last checkpoint, relative defect at most SYM_TOL. The initial state is
#     bitwise symmetric (the IC evaluates odd functions of 2x - 1, 2y - 1); the scheme is symmetric in exact arithmetic
#     but its left- and right-biased reconstructions round differently, so the defect is a round-off seed amplified by
#     the flow: measured (CPU, M2-P6a) 1.7e-13 at step 100, 3e-11 at step 300 when the shocks form, 1.3e-9 at t = 0.5
#     (1281 steps), against a truncation-level defect (~1e-3) for a wrong fill or a sign error; w and bz stay exactly 0
#     on the CPU and carry round-off noise on FNL (~1e-12 of the in-plane field), scaled by the vector magnitude;
#   - the integrals of rho, rho u, E and B constant to CONS_TOL relative (measured 1.2e-14);
#   - positivity with zero floored cells: the floors are disabled, so a non-positive density or pressure stops the run.
# About 4.5 min on the CPU.
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
SYM_TOL="1.0e-8"
CONS_TOL="1.0e-13"

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

echo ">> MV-12 Orszag-Tang 128^2 ($(basename "$EXE"), np $NP)"
w="$CASE_DIR/work-$TAG"
rm -rf "$w" ; mkdir -p "$w"
"$VENV_PY" "$CASE_DIR/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/orszag-tang.ini" --cells 128
if ! (cd "$w" && mpirun -np "$NP" "$EXE" orszag-tang.ini > log.txt 2>&1); then
   echo "check.sh: Orszag-Tang run failed (a floor stop is a positivity failure), see $w/log.txt" >&2
   exit 1
fi
if "$VENV_PY" "$CASE_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$SYM_TOL" --cons-tol "$CONS_TOL"; then
   STATUS=0
else
   STATUS=1
fi
find "$w" -name '*.h5' -delete

if [[ $STATUS -eq 0 ]]; then
   echo "MV-12 PASSED ($TAG)"
else
   echo "MV-12 FAILED ($TAG)"
   exit 1
fi
