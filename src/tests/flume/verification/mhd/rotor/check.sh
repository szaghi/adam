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
# Usage: ./check.sh [--np N] [--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]] [--legs LIST]   LIST: mv13 or all
#
# --numerics runs the legs on `scheme_space = weno-riemann` (mhd/numerics.sh, issue #47 M3-P3c). The final defect
# grows smoothly from round-off (time series, M3-P3c: no jump, so no switch flip) and orders by the dissipation of the
# solver and by the backend's round-off (at M3-P3c, October: final 3.4e-11 HLL, 1.8e-9 split, 5.4e-9 HLLD primitive,
# 7.4e-8 without the sensor, 1.1e-7 HLLD characteristic on the CPU, 2.3e-7 and 6.1e-6 for the two HLLD on FNL, FMA
# contraction). So the
# symmetry is asserted twice: at step 100, before the flow amplifies the seed, within EARLY_SYM_TOL 1e-8 (measured CPU /
# FNL 3.9e-12 / 2.1e-12 HLLD primitive, 4.9e-10 / 7.9e-10 HLLD characteristic; a symmetry defect of the scheme is at
# truncation level, 1e-3, from the first steps), and at the end within a bound per spec (issue #142 A9): 10 times the
# larger of the CPU and FNL defects measured on medooza (develop 0bcdf870, CPU gfortran 16.2 / FNL nvfortran 26.9 on
# an A30, job 1369317, 2026-11-22), rounded up to one digit:
#   hll:*:6th:weno                  5e-10 (CPU 4.4e-11, FNL 1.2e-11)
#   hlld:primitive:6th:weno         1e-4  (CPU 6.08e-5, FNL 6.08e-5: the SAME value on both backends, so not round-off
#                                          amplified by the flow but a deterministic asymmetry, 4 orders above the
#                                          5.4e-9 of M3-P3c; under investigation, the bound stays the former 1e-4)
#   hlld:primitive:6th:none         4e-6  (CPU 3.4e-7, FNL 1.8e-8)
#   hlld:characteristic:6th:weno    8e-9  (CPU 7.4e-10, FNL 3.6e-10)
#   any other spec                  1e-4  (not measured: gross failures only)
# Without --numerics (the splitting scheme) the end bound stays SYM_TOL 1e-8 (CPU 6.5e-10, FNL 6.4e-10).
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
source "$VERIF_DIR/legs.sh"
legs_declare "mv13"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      --numerics) NUMERICS="$2" ; shift 2 ;;
      --legs) legs_select "$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N, --numerics SPEC, --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
numerics_check
if [[ -n $NUMERICS ]]; then
   case "$NUMERICS" in # the end-of-run bound per spec (see the header)
      hll:*:6th:weno)               SYM_TOL="5.0e-10" ;;
      hlld:primitive:6th:weno)      SYM_TOL="1.0e-4" ;;
      hlld:primitive:6th:none)      SYM_TOL="4.0e-6" ;;
      hlld:characteristic:6th:weno) SYM_TOL="8.0e-9" ;;
      *)                            SYM_TOL="1.0e-4" ;;
   esac
fi
MAKE_OPTS=()
[[ -n $NUMERICS ]] && MAKE_OPTS=(--it-save 100)
TAG="$(basename "$EXE")-np$NP$(numerics_tag)"

STATUS=0
if has_leg mv13; then
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
      "$VENV_PY" "$OT_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$EARLY_SYM_TOL" --field-parity even --step 100 \
         || STATUS=1
   fi
   find "$w" -name '*.h5' -delete

fi
if [[ $STATUS -eq 0 ]]; then
   echo "MV-13 PASSED (legs: $LEGS) ($TAG)"
else
   echo "MV-13 FAILED (legs: $LEGS) ($TAG)"
   exit 1
fi
