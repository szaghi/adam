#!/usr/bin/env bash
# FLUME MHD verification MV-13 (issue #41, M2-P6): the MHD rotor, symmetry and positivity.
#
# Why: a dense disc spun up in a magnetised ambient (Balsara and Spicer 1999, Toth 2000 first rotor): torsional Alfven
# waves, field compression and strong shocks, a stringent 2-D positivity test. 256^2 (a power of two, so the initial
# state is bitwise symmetric; Toth uses 200^2), t = 0.15, GLM, np 2. Through orszag_tang_oracle.py --field-parity even:
#   - the symmetry of the 180 degrees rotation composed with B -> -B (the ambient field is uniform), relative defect at
#     most SYM_TOL: as for MV-12 a round-off seed amplified by the flow, measured (CPU, M2-P6b) 1e-11 at step 100,
#     1.8e-9 at t = 0.15 (1152 steps); w and bz stay exactly 0 on the CPU;
#   - positivity with zero floored cells: the floors are disabled, so a non-positive density or pressure stops the run.
# No conservation check: the outflow boundaries do not conserve. About 10 min on the CPU.
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
OT_DIR="$(cd "$CASE_DIR/../orszag-tang" && pwd)"

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

echo ">> MV-13 rotor 256^2 ($(basename "$EXE"), np $NP)"
w="$CASE_DIR/work-$TAG"
rm -rf "$w" ; mkdir -p "$w"
"$VENV_PY" "$CASE_DIR/make_rotor.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/rotor.ini" --cells 256
if ! (cd "$w" && mpirun -np "$NP" "$EXE" rotor.ini > log.txt 2>&1); then
   echo "check.sh: rotor run failed (a floor stop is a positivity failure), see $w/log.txt" >&2
   exit 1
fi
if "$VENV_PY" "$OT_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$SYM_TOL" --field-parity even; then
   STATUS=0
else
   STATUS=1
fi
find "$w" -name '*.h5' -delete

if [[ $STATUS -eq 0 ]]; then
   echo "MV-13 PASSED ($TAG)"
else
   echo "MV-13 FAILED ($TAG)"
   exit 1
fi
