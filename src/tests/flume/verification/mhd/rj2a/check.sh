#!/usr/bin/env bash
# FLUME MHD verification MV-4 (issue #41, M2-P3): the RJ2a Riemann problem against its exact solution.
#
# Why: RJ2a (Ryu & Jones 1995 Fig. 2a; exact 7-state solution of Dai & Woodward 1994, Tables Ia/Ib, as in Athena++
# shock_tube.cpp) develops all seven MHD waves from one discontinuity: it checks the whole MHD path (Roe-Balsara
# eigensystem at the face average, characteristic splitting, fluxes, rotation of the transverse components). Through
# rj2a_oracle.py, with gamma = 5/3, t = 0.2, WENO-5 characteristic, SSP-33, CFL 0.5, np 2 (the Sod grid of V1):
#   1. N = 256 along x, y, z (no cleaning): L1 of each conservative variable against the exact solution, their sum
#      below L1_MAX_256, and the y and z runs equal the x run in the rotated frame BITWISE (DIR_TOL = 0: the 3-term
#      |u|^2, |B|^2, u.B sums are order-independent, mhd_sum3, and the projections sum in the frame order, M2-P3d);
#   2. N = 512 along x: the L1 sum below L1_MAX_512 and below the N = 256 one (the solution converges; at first order,
#      as discontinuities dominate);
#   3. GLM (N = 256, x): equal to the run without cleaning BITWISE on the 8 shared variables, psi exactly zero (in 1-D
#      B_n is uniform, so the (B_n, psi) block is inert and block-diagonal). Under --numerics (weno-riemann) within
#      PAIR_TOL = 1e-11: the face values of the interpolated uniform B_n differ by round-off, which seeds psi (measured
#      9.2e-16, HLLD primitive) and, through the shocks, the shared variables (1.3e-12), the scale of the x/y/z
#      round-off before M2-P3d;
#   4. (--numerics with HLLD) RV-6 of issue #47: the HLLD L1 sum at N = 256 at most the HLL one;
#   5. EV-3 of issue #47 (M3-P4c): EGLM (N = 256, x) equal to GLM BITWISE on the 8 shared variables, psi exactly zero
#      (psi = 0 makes the EGLM energy coupling and sources vanish); under --numerics within PAIR_TOL (measured, M3-P4b,
#      HLLD characteristic: 6.3e-13 CPU, 0 FNL).
# One run at a time; checkpoints deleted after use.
#
# Usage: ./check.sh [--np N] [--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]]
#
# --numerics runs the legs on `scheme_space = weno-riemann` (mhd/numerics.sh, issue #47 M3-P3c). The specs
# hlld:primitive and hlld:characteristic (6th, weno) assert the bounds measured on them (CPU, M3-P3c);
# any other spec asserts only the scheme-independent checks (a measurement).
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
L1_MAX_256="3.919e-02" # CPU and FNL baseline 3.842043e-02 (M2-P3b), plus 2%
L1_MAX_512="2.159e-02" # CPU and FNL baseline 2.116856e-02 (M2-P3b), plus 2%
DIR_TOL="0.0"          # x/y/z in the rotated frame, bitwise (1.4e-12 before M2-P3d)
PAIR_TOL="0.0"         # GLM vs no cleaning, bitwise; weno-riemann (--numerics): 1e-11, see the GLM leg

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
ORACLE="$CASE_DIR/rj2a_oracle.py"
numerics_check
TAG="$(basename "$EXE")-np$NP$(numerics_tag)"
[[ -n $NUMERICS ]] && PAIR_TOL="1.0e-11"
# [weno] weights = si (FLUME_EXE = scaling/weights-exe.sh, WEIGHTS_EXE set; issue #49): the smoothness indicator of a
# uniform field is round-off, not zero (a quadratic form), and the si descaler varies face to face across the shocks,
# so the uniform B_n is reconstructed with face values that differ by round-off, as under --numerics: PAIR_TOL 1e-11
# (measured: GLM vs no cleaning 1.1e-13, psi 1.3e-15; EGLM vs GLM 2.5e-13). The default js weights stay bitwise.
[[ -n ${WEIGHTS_EXE:-} && ${WEIGHTS:-si} == si ]] && PAIR_TOL="1.0e-11"
# CPU weno-riemann HLLD (M3-P3c), N = 256, 512: primitive 5.132654e-02, 2.926153e-02; characteristic 3.764077e-02,
# 2.122612e-02; plus 2 %
case "$NUMERICS" in
   "")                           ;;
   hlld:primitive:6th:weno)      L1_MAX_256="5.235e-02" ; L1_MAX_512="2.985e-02" ;;
   hlld:characteristic:6th:weno) L1_MAX_256="3.840e-02" ; L1_MAX_512="2.166e-02" ;;
   *)                            L1_MAX_256="" ; L1_MAX_512="" ;;
esac
B256=(--l1-max "$L1_MAX_256") ; B512=(--l1-max "$L1_MAX_512")
[[ -z $L1_MAX_256 ]] && B256=() && B512=()
FAILED=0

case_run() { # case_run <axis> <cells> <variant>: run one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$3-$2-$1"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_rj2a.py" "$VERIF_DIR/sod/sod-$1.ini" "$w/rj2a-$1.ini" --axis "$1" --cells "$2" \
      --divergence-control "$3"
   numerics_apply "$w/rj2a-$1.ini"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" "rj2a-$1.ini" > log.txt 2>&1); then
      echo "check.sh: RJ2a $1 N=$2 $3 run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

echo ">> MV-4 RJ2a, N = 256, x/y/z, no cleaning ($(basename "$EXE"), np $NP)"
W256=()
for axis in x y z; do w="$(case_run $axis 256 none)" || exit 1 ; W256+=("$w"); done
"$VENV_PY" "$ORACLE" "${W256[@]}" "${B256[@]}" --dir-tol "$DIR_TOL" || FAILED=1
echo ">> MV-4 RJ2a, N = 512, x, no cleaning"
W512="$(case_run x 512 none)" || exit 1
"$VENV_PY" "$ORACLE" "$W512" "${B512[@]}" || FAILED=1
sum256=$("$VENV_PY" "$ORACLE" "${W256[0]}" | sed -n 's/.*sum \([0-9.e+-]*\).*/\1/p')
sum512=$("$VENV_PY" "$ORACLE" "$W512" | sed -n 's/.*sum \([0-9.e+-]*\).*/\1/p')
if "$VENV_PY" -c "import sys; sys.exit(0 if $sum512 < $sum256 else 1)"; then
   echo "   convergence: L1 sum N=512 $sum512 < N=256 $sum256  PASS"
else
   echo "   convergence: L1 sum N=512 $sum512 not below N=256 $sum256  FAIL" ; FAILED=1
fi
find "$W512" -name '*.h5' -delete
if [[ $NUMERICS_SOLVER == hlld ]]; then
   echo ">> RV-6 RJ2a, N = 256, x, no cleaning, HLL: the HLLD L1 sum at most the HLL one"
   NUMERICS_SOLVER=hll ; TAG_HLLD="$TAG" ; TAG="${TAG/-hlld-/-hll-}"
   WHLL="$(case_run x 256 none)" || exit 1
   NUMERICS_SOLVER=hlld ; TAG="$TAG_HLLD"
   sumhll=$("$VENV_PY" "$ORACLE" "$WHLL" | sed -n 's/.*sum \([0-9.e+-]*\).*/\1/p')
   if "$VENV_PY" -c "import sys; sys.exit(0 if $sum256 <= $sumhll else 1)"; then
      echo "   HLLD L1 sum $sum256 <= HLL $sumhll  PASS"
   else
      echo "   HLLD L1 sum $sum256 > HLL $sumhll  FAIL" ; FAILED=1
   fi
   find "$WHLL" -name '*.h5' -delete
fi
echo ">> MV-4 RJ2a, N = 256, x, GLM against no cleaning"
WGLM="$(case_run x 256 glm)" || exit 1
"$VENV_PY" "$ORACLE" --pair "${W256[0]}" "$WGLM" --pair-tol "$PAIR_TOL" || FAILED=1
echo ">> EV-3 RJ2a, N = 256, x, EGLM against GLM"
WEGLM="$(case_run x 256 eglm)" || exit 1
"$VENV_PY" "$ORACLE" --pair "$WGLM" "$WEGLM" --pair-tol "$PAIR_TOL" || FAILED=1
for w in "${W256[@]}" "$WGLM" "$WEGLM"; do find "$w" -name '*.h5' -delete; done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-4 PASSED ($TAG)"
else
   echo "MV-4 FAILED ($TAG)"
   exit 1
fi
