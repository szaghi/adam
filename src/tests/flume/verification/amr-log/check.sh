#!/usr/bin/env bash
# FLUME verification AL (issue #68, M6-A4): integral drift without the reflux, and the WENO sensor near 2:1 seams.
#
# Why: the reflux policy of #68 E turns the Berger-Colella reflux off at the seam faces a WENO sensor judges smooth.
# Whether that is safe, and where its threshold sits, rests on two logs no check kept (amr_log.py):
#   drift   the relative drift of the domain integrals with the reflux on and off: its maximum, and its growth (mean
#           |drift| over the last quarter over the second quarter; about 1 bounded, 2 or more accumulating);
#   sensor  the WENO-5 sensor r = min_k w_k/d_k on (rho, u_n, p) of every face, per checkpoint, in the seam band and
#           the interior: the minimum and the fraction below tau = 0.2 (the weno-riemann threshold).
# Legs:
#   vortex  the AO Euler vortex (verification/amr-order), quadtree, refined quadrant [0.5, 1]^2, N = 64, run to t = 1
#           (one period of the free stream (1, 1): the vortex crosses every seam), periodic: every integral is
#           conserved;
#   sod     the AS Sod (verification/amr-shock), c2f, N = 200, t = 0.2: the shock crosses the 2:1 face at t = 0.143;
#           mass and energy are conserved (no flux through the ends before t = 0.2).
# Asserted: with the reflux, every logged integral drifts at most DRIFT_MAX (round-off). Everything else is a log.
# Cost, 2 ranks: about 11 min on the CPU, 5 min on FNL (the two backends log the same values).
# --reflux-policy P (#68 E2, the A4 rows under `troubled`) passes [numerics] reflux_policy = P to the reflux-on runs,
# whose work directories then end in -reflux-P; under `troubled` the reflux is partial, so their drift is logged, not
# asserted.
#
# Usage: ./check.sh [--np N] [--legs vortex,sod|all] [--reflux-policy always|troubled] [--keep]   (default: vortex,sod;
#        --leg L is an older spelling of --legs, repeatable)
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF="$(cd "$CASE_DIR/.." && pwd)"
VERIF_DIR="$VERIF"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "vortex sod" "vortex sod"
NP=2
KEEP=0
POLICY=""
DRIFT_MAX="1.0e-13"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) legs_select "$2" ; shift 2 ;;
      --leg)  legs_select "$2" ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      --reflux-policy) POLICY="$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs L, --leg L, --reflux-policy P," \
                   "--keep)" >&2
              exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
TAG="$(basename "$EXE")-np$NP"
FAILED=0
POLICY_OPTS=() ; ON="reflux"
if [[ -n $POLICY ]]; then
   [[ $POLICY =~ ^(always|troubled)$ ]] || { echo "check.sh: --reflux-policy always|troubled" >&2 ; exit 2 ; }
   POLICY_OPTS=(--reflux-policy "$POLICY") ; ON="reflux-$POLICY"
fi

run() { # run <work> <case> <generator> [generator options]: write the input, run it, report the wall time
   local work=$1 case=$2 gen=$3 ; shift 3
   rm -rf "${work:?}" ; mkdir -p "$work"
   "$VENV_PY" "$gen" "$case" "$work/$case.ini" "$@"
   local t0=$SECONDS
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$case.ini" < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $work/log.txt" ; FAILED=1 ; return 1
   fi
   echo "   $(basename "$work"): $(( SECONDS - t0 )) s"
}

echo ">> AL: integral drift and WENO sensor near 2:1 seams ($(basename "$EXE"), np $NP)"
for leg in $LEGS; do
   pre="$CASE_DIR/work-$TAG-$leg"
   case "$leg" in
      vortex)
         echo "-- vortex: AO euler, N = 64, t = 1, reflux on and off"
         run "$pre-$ON" euler "$VERIF/amr-order/make_amr_order.py" --cells 64 --time-max 1.0 --it-save 10 \
             "${POLICY_OPTS[@]}" || continue
         run "$pre-noreflux" euler "$VERIF/amr-order/make_amr_order.py" --cells 64 --time-max 1.0 --it-save 10 \
             --no-reflux || continue
         vars=(r ru rv rE) ; where=(--box 0.5 0.5 1.0 1.0) ;;
      sod)
         echo "-- sod: AS sod c2f, N = 200, t = 0.2, reflux on and off"
         run "$pre-$ON" sod "$VERIF/amr-shock/make_amr_shock.py" --cells 200 --side c2f --it-save 10 \
             "${POLICY_OPTS[@]}" || continue
         run "$pre-noreflux" sod "$VERIF/amr-shock/make_amr_shock.py" --cells 200 --side c2f --it-save 10 \
             --no-reflux || continue
         vars=(r rE) ; where=(--seam-x 0.75) ;;
   esac
   if [[ $POLICY == troubled ]]; then
      "$VENV_PY" "$CASE_DIR/amr_log.py" drift "$pre-$ON" --vars "${vars[@]}"
   else
      "$VENV_PY" "$CASE_DIR/amr_log.py" drift "$pre-$ON" --vars "${vars[@]}" --max "$DRIFT_MAX" || FAILED=1
   fi
   "$VENV_PY" "$CASE_DIR/amr_log.py" drift "$pre-noreflux" --vars "${vars[@]}"
   "$VENV_PY" "$CASE_DIR/amr_log.py" sensor "$pre-$ON" "${where[@]}"
done
if [[ $KEEP -eq 0 ]]; then
   find "$CASE_DIR" -path "$CASE_DIR/work-$TAG-*" -name '*.h5' -delete
fi

if [[ $FAILED -eq 0 ]]; then
   echo "AL PASSED ($TAG, legs: $LEGS)"
else
   echo "AL FAILED ($TAG, legs: $LEGS)"
   exit 1
fi
