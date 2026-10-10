#!/usr/bin/env bash
# FLUME verification AO (issue #68, M6-A1): the order of accuracy of a smooth flow across static 2:1 seams.
#
# Why: every smooth order oracle (V2, MV-5..7) runs on a uniform grid. Here a vortex starts centred on the corner of a
# refined quadrant (make_amr_order.py: an x seam, a y seam and their corner through its core, convected into the fine
# quadrant), and amr_order_oracle.py measures the error per region against the exact solution: the seam band (cells
# within ngc of their own spacing from the refined box boundary) and the interior. The #68 audit (F1, Appendix A)
# predicts a second-order band, from the mean restriction of the coarse ghosts and the Berger-Colella reflux of
# point-value fluxes, and design order in the interior; the uniform control (the same case without the box) is the
# design-order reference.
#
# Measured (CPU baseline, M6-A1, finest pair N = 128 -> 256): the whole composite solution is second order, not only
# the band: euler L1 order all 2.01, band 1.71, interior 1.91 (Linf 1.46); mhd all 2.07, band 1.85, interior 1.92 (Linf
# 1.89). The error made at the seams is carried into the interior by the flow (the vortex crosses into the fine
# quadrant), so a seam caps every region a feature reaches after crossing it. At N = 256 the refined composite error is
# 48x (euler) and 70x (mhd) the error of the uniform run at the coarse spacing, whose L1 orders are 5.79 and 5.34: on
# this smooth flow the refinement makes the solution worse. Bounds: the finest refined errors at most the CPU baseline
# plus 2 % (an improvement of the transfers or of the reflux passes; a degradation fails); the uniform control at
# least ORDER_UNIFORM_MIN (the instrument itself at design order).
#
# Legs (each: the refined ladder and the uniform control, base N = 64/128/256):
#   euler         the V2 isentropic vortex, quadtree, t = 0.1, CFL 0.4; density errors
#   mhd           the MV-7 magnetised vortex, quadtree, GLM, t = 0.5, CFL 0.4; the 8 conservative variables
#   euler-octree  the euler leg on an octree (null z, nk = 4), N = 64/128 only (not in the default set: an octree
#                 refines along the null z too, about 16x the cost of the quadtree for the same 2-D data)
#
# Usage: ./check.sh [--np N] [--leg euler|mhd|euler-octree ...] [--keep]
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
LEGS=()
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --leg)  LEGS+=("$2") ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --leg L, --keep)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(euler mhd)
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
TAG="$(basename "$EXE")-np$NP"
FAILED=0
ORDER_UNIFORM_MIN="4.7"
# CPU baseline (M6-A1) plus 2 %, the finest refined run (N = 256): region:norm:max
MAX_EULER=(all:l1:4.129e-06 band:l1:1.607e-05 all:linf:4.772e-04)
MAX_MHD=(all:l1:1.871e-06 band:l1:8.102e-06 all:linf:1.295e-04)
NO_BOUNDS=()                                                         # euler-octree: a measurement

run() { # run <work> <case> <cells> [make_amr_order options]: write the input, run it, report the wall time
   local work=$1 case=$2 cells=$3 ; shift 3
   rm -rf "${work:?}" ; mkdir -p "$work"
   "$VENV_PY" "$CASE_DIR/make_amr_order.py" "$case" "$work/$case.ini" --cells "$cells" "$@"
   local t0=$SECONDS
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$case.ini" < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $work/log.txt" ; FAILED=1 ; return 1
   fi
   echo "   $(basename "$work"): $(( SECONDS - t0 )) s"
}

ladder() { # ladder <leg> <case> <box> <ratio> <bounds-array-name> <cells...>: the refined and the uniform runs, then
   # the oracle on each
   local leg=$1 case=$2 box=$3 ratio=$4 ; shift 4
   local -n bounds=$1 ; shift
   local refined=() uniform=()
   for n in "$@"; do
      run "$CASE_DIR/work-$TAG-$leg-n$n" "$case" "$n" --ratio "$ratio" && refined+=("$CASE_DIR/work-$TAG-$leg-n$n")
      run "$CASE_DIR/work-$TAG-$leg-uniform-n$n" "$case" "$n" --ratio "$ratio" --uniform && \
         uniform+=("$CASE_DIR/work-$TAG-$leg-uniform-n$n")
   done
   echo "   refined ($leg):"
   # shellcheck disable=SC2086
   local maxima=() ; [[ ${#bounds[@]} -gt 0 ]] && maxima=(--max "${bounds[@]}")
   # shellcheck disable=SC2086
   "$VENV_PY" "$CASE_DIR/amr_order_oracle.py" "$case" "${refined[@]}" --box $box "${maxima[@]}" || FAILED=1
   echo "   uniform control ($leg):"
   "$VENV_PY" "$CASE_DIR/amr_order_oracle.py" "$case" "${uniform[@]}" --assert "interior:l1:$ORDER_UNIFORM_MIN" || FAILED=1
   if [[ $KEEP -eq 0 ]]; then
      for w in "${refined[@]}" "${uniform[@]}"; do find "$w" -name '*.h5' -delete; done
   fi
}

echo ">> AO: order across static 2:1 seams ($(basename "$EXE"), np $NP)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      euler)        echo "-- euler: V2 vortex on a quadtree, quadrant [0.5, 1]^2 refined"
                    ladder euler euler "0.5 0.5 1.0 1.0" 4 MAX_EULER 64 128 256 ;;
      mhd)          echo "-- mhd: MV-7 vortex on a quadtree, quadrant [0, 7]^2 refined"
                    ladder mhd mhd "0.0 0.0 7.0 7.0" 4 MAX_MHD 64 128 256 ;;
      euler-octree) echo "-- euler-octree: V2 vortex on an octree, quadrant [0.5, 1]^2 refined"
                    ladder euler-octree euler "0.5 0.5 1.0 1.0" 8 NO_BOUNDS 64 128 ;;
      *)            echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "AO PASSED ($TAG)"
else
   echo "AO FAILED ($TAG)"
   exit 1
fi
