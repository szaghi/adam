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
# Usage: ./check.sh [--np N] [--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]]
#
# --numerics runs the legs on `scheme_space = weno-riemann` (mhd/numerics.sh, issue #47 M3-P3c). The final defect
# grows smoothly from round-off (time series, M3-P3c: no jump, so no switch flip) and orders by the dissipation of the
# solver and by the backend's round-off: final 3.4e-11 HLL, 1.8e-9 split, 5.4e-9 HLLD primitive, 7.4e-8 without the
# sensor, 1.1e-7 HLLD characteristic on the CPU, 2.3e-7 and 6.1e-6 for the two HLLD on FNL (FMA contraction). So the
# symmetry is asserted twice: at step 100, before the flow amplifies the seed, within EARLY_SYM_TOL 1e-8 (measured CPU /
# FNL 3.9e-12 / 2.1e-12 HLLD primitive, 4.9e-10 / 7.9e-10 HLLD characteristic; a symmetry defect of the scheme is at
# truncation level, 1e-3, from the first steps), and at the end within SYM_TOL 1e-4 (gross failures only).
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
source "$CASE_DIR/../numerics.sh"
SYM_TOL="1.0e-8"
EARLY_SYM_TOL="1.0e-8"
OT_DIR="$(cd "$CASE_DIR/../orszag-tang" && pwd)"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      --numerics) NUMERICS="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N, --numerics SPEC)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
numerics_check
[[ -n $NUMERICS ]] && SYM_TOL="1.0e-4"
MAKE_OPTS=()
[[ -n $NUMERICS ]] && MAKE_OPTS=(--it-save 100)
TAG="$(basename "$EXE")-np$NP$(numerics_tag)"

echo ">> MV-13 rotor 256^2 ($(basename "$EXE"), np $NP)"
w="$CASE_DIR/work-$TAG"
rm -rf "$w" ; mkdir -p "$w"
"$VENV_PY" "$CASE_DIR/make_rotor.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/rotor.ini" --cells 256 "${MAKE_OPTS[@]}"
numerics_apply "$w/rotor.ini"
if ! (cd "$w" && mpirun -np "$NP" "$EXE" rotor.ini > log.txt 2>&1); then
   echo "check.sh: rotor run failed (a floor stop is a positivity failure), see $w/log.txt" >&2
   exit 1
fi
if "$VENV_PY" "$OT_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$SYM_TOL" --field-parity even; then
   STATUS=0
else
   STATUS=1
fi
if [[ -n $NUMERICS ]]; then
   "$VENV_PY" "$OT_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$EARLY_SYM_TOL" --field-parity even --step 100 || STATUS=1
fi
find "$w" -name '*.h5' -delete

if [[ $STATUS -eq 0 ]]; then
   echo "MV-13 PASSED ($TAG)"
else
   echo "MV-13 FAILED ($TAG)"
   exit 1
fi
