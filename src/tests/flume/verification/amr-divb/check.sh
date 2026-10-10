#!/usr/bin/env bash
# FLUME verification AD (issue #68, M6-A2): D_flux, the divergence the cleaning acts on, across static 2:1 seams.
#
# Why: the div(B) column of every MHD check is the centred divergence of the cell values, not the divergence GLM and
# EGLM act on (#68 F4). `[diagnostics] divb_flux` adds D_flux: the divergence of the interface normal field read from
# the scheme's own psi face fluxes, the coarse 2:1 faces taking the fine fluxes (the reflux), so that it sums exactly
# over a composite grid (Gauss: on a periodic domain its signed integral is round-off; without the fine fluxes it is
# not). make_amr_divb.py writes the cases, amr_divb_oracle.py measures them.
#
# Legs:
#   cpaw  the source-free CPAW (MV-6, GLM, inclined, two x seams) for 10 periods, refined N = 64/128 (--long adds 256)
#         and the uniform controls: the seam D_flux floor must stay bounded in time (growth, last quarter over second
#         quarter, at most GROWTH_MAX) and converge (finest-pair order at least FLOOR_ORDER_MIN). Measured on FNL
#         (M6-A2): floors 5.39e-4, 5.43e-5, 2.12e-5 (orders +3.31, +1.36: first order asymptotically, the O(H) divergence
#         of an O(H^2) flux inconsistency at the 2:1 face, #68 F1), growth 0.95, 0.72, 0.80; uniform floors 9.97e-6,
#         4.57e-7, 6.74e-9, so the seam floor is 54x, 119x, 3146x the uniform one: the seams are the divergence error.
#         Cost: about 12 min on FNL (N = 256: 30 more), about 2 h on the CPU (N = 128 dominates): run it on FNL.
#   loop  the field loop of MV-9 under GLM and EGLM (#78): refined N = 64 (box [0, 1] x [-0.5, 0.5]) against the
#         uniform N = 128, t = 1: the D_flux columns and <|B_z|>/A0 reported, no assertion (#78 is fixed or explained
#         in #68 G3).
#
# Usage: ./check.sh [--np N] [--leg cpaw|loop ...] [--long] [--keep]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NP=2
KEEP=0
LONG=0
LEGS=()
GROWTH_MAX="1.5"
FLOOR_ORDER_MIN="1.0"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --leg)  LEGS+=("$2") ; shift 2 ;;
      --long) LONG=1 ; shift ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --leg L, --long, --keep)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(cpaw loop)
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
TAG="$(basename "$EXE")-np$NP"
FAILED=0

run() { # run <work> <case> [make_amr_divb options]: write the input, run it, report the wall time
   local work=$1 case=$2 ; shift 2
   rm -rf "${work:?}" ; mkdir -p "$work"
   "$VENV_PY" "$CASE_DIR/make_amr_divb.py" "$case" "$work/$case.ini" "$@"
   local t0=$SECONDS
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$case.ini" < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $work/log.txt" ; FAILED=1 ; return 1
   fi
   echo "   $(basename "$work"): $(( SECONDS - t0 )) s"
}

echo ">> AD: D_flux across static 2:1 seams ($(basename "$EXE"), np $NP)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      cpaw)
         echo "-- cpaw: source-free CPAW, GLM, two x seams, 10 periods"
         ladder=(64 128) ; [[ $LONG -eq 1 ]] && ladder+=(256)
         refined=() ; uniform=()
         for n in "${ladder[@]}"; do
            w="$CASE_DIR/work-$TAG-cpaw-n$n" ; run "$w" cpaw --cells "$n" && refined+=("$w")
            w="$CASE_DIR/work-$TAG-cpaw-uniform-n$n" ; run "$w" cpaw --cells "$n" --uniform && uniform+=("$w")
         done
         "$VENV_PY" "$CASE_DIR/amr_divb_oracle.py" growth --refined "${refined[@]}" --uniform "${uniform[@]}" \
            --ratio-max "$GROWTH_MAX" --floor-order-min "$FLOOR_ORDER_MIN" || FAILED=1 ;;
      loop)
         echo "-- loop: field loop, GLM and EGLM, refined N = 64 against uniform N = 128, t = 1 (#78)"
         works=()
         for c in glm eglm; do
            w="$CASE_DIR/work-$TAG-loop-$c-n64" ; run "$w" loop --cells 64 --cleaning "$c" --history 1 && works+=("$w")
            w="$CASE_DIR/work-$TAG-loop-$c-uniform-n128"
            run "$w" loop --cells 128 --cleaning "$c" --uniform --history 1 && works+=("$w")
         done
         "$VENV_PY" "$CASE_DIR/amr_divb_oracle.py" loop "${works[@]}" || FAILED=1 ;;
      *) echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done
if [[ $KEEP -eq 0 ]]; then
   find "$CASE_DIR" -path "$CASE_DIR/work-$TAG-*" -name '*.h5' -delete
fi

if [[ $FAILED -eq 0 ]]; then
   echo "AD PASSED ($TAG)"
else
   echo "AD FAILED ($TAG)"
   exit 1
fi
