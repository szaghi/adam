#!/usr/bin/env bash
# FLUME MHD verification MV-8 (issue #41, M2-P6d): the rotated shock tube (Toth 2000, J. Comput. Phys. 161, 6.3.2).
#
# Why: Ryu-Jones 1a rotated in the x-y plane makes the dimension-by-dimension scheme solve a 1-D problem on 2-D data;
# the exact normal field stays uniform, B_par = 5/sqrt(4 pi), so its deviation measures the divergence error of the
# base scheme with GLM, and the whole solution must match the problem solved in 1-D. Periodic square (FLUME has no
# shifted-periodic boundary): the jump at s = 1/4 and a resolved ramp back to the left state (make_rotated_shock_tube.py
# explains why a second jump cannot be used). Through rotated_shock_tube_oracle.py, gamma = 5/3, GLM (c_h 24,
# glm_ch_check = error), WENO-5 characteristic, SSP-54, CFL 0.4, floors off, np 2, for each normal:
#   * tan^-1 2 (Toth's 63.4 degrees, normal (1, 2)) and 45 degrees (normal (1, 1); not degenerate for this scheme, unlike
#     Toth's CT and CD schemes: B_par is conserved to truncation, not to round-off);
#   * the 1-D reference: the same problem along x, 1024 cells per period (4 times the normal resolution of N = 256);
#   * N = 128, 256: dB_par (Toth eq. 45 against the analytic value) and the mean relative L1 error of
#     rho, v_par, v_perp, p, B_par, B_perp against the reference below the per-N bounds, both decreasing with N.
# For reference, Toth's second-order base scheme at 63.4 degrees, N = 256 on his strip (twice the cells across the fan
# of N = 256 here): dB_par 0.0037, mean error 0.0238. One run at a time; checkpoints deleted after use; ~30 min.
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
LADDER=(128 256)
REFERENCE=256                      # the 1-D reference has 4 x 256 cells per period
DBPAR_MAX_63=(6.146e-03 2.842e-03) # CPU baseline (M2-P6d) plus 2 %, per N of the ladder
DMEAN_MAX_63=(8.989e-02 4.095e-02)
DBPAR_MAX_45=(6.424e-03 4.252e-03)
DMEAN_MAX_45=(5.440e-02 2.600e-02)

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
STATUS=0

case_run() { # case_run <angle> <cells> <2d|ref>: run one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$1-$3-$2" base="$VERIF_DIR/vortex/vortex-n064.ini" extra=()
   if [[ "$3" == ref ]]; then base="$VERIF_DIR/sod/sod-x.ini" ; extra=(--reference) ; fi
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_rotated_shock_tube.py" "$base" "$w/rst.ini" --angle "$1" --cells "$2" "${extra[@]}"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" rst.ini > log.txt 2>&1); then
      echo "check.sh: rotated shock tube $1 degrees $3 N=$2 run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

for angle in 63 45; do
   echo ">> MV-8 rotated shock tube, normal angle $angle, N = ${LADDER[*]} ($(basename "$EXE"), np $NP)"
   ref="$(case_run "$angle" "$REFERENCE" ref)" || exit 1
   works=()
   for n in "${LADDER[@]}"; do
      w="$(case_run "$angle" "$n" 2d)" || exit 1
      works+=("$w")
   done
   dbpar="DBPAR_MAX_$angle[@]" ; dmean="DMEAN_MAX_$angle[@]"
   if ! "$VENV_PY" "$CASE_DIR/rotated_shock_tube_oracle.py" "${works[@]}" --reference "$ref" \
         --dbpar-max "${!dbpar}" --dmean-max "${!dmean}"; then
      STATUS=1
   fi
   for w in "$ref" "${works[@]}"; do find "$w" -name '*.h5' -delete; done
done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-8 PASSED ($TAG)"
else
   echo "MV-8 FAILED ($TAG)"
   exit 1
fi
