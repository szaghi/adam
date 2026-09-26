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
#      B_n is uniform, so the (B_n, psi) block is inert and block-diagonal).
# One run at a time; checkpoints deleted after use.
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
L1_MAX_256="3.919e-02" # CPU and FNL baseline 3.842043e-02 (M2-P3b), plus 2%
L1_MAX_512="2.159e-02" # CPU and FNL baseline 2.116856e-02 (M2-P3b), plus 2%
DIR_TOL="0.0"          # x/y/z in the rotated frame, bitwise (1.4e-12 before M2-P3d)

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
ORACLE="$CASE_DIR/rj2a_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

case_run() { # case_run <axis> <cells> <variant>: run one case, print its work directory
   local w="$CASE_DIR/work-$TAG-$3-$2-$1"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_rj2a.py" "$VERIF_DIR/sod/sod-$1.ini" "$w/rj2a-$1.ini" --axis "$1" --cells "$2" \
      --divergence-control "$3"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" "rj2a-$1.ini" > log.txt 2>&1); then
      echo "check.sh: RJ2a $1 N=$2 $3 run failed, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

echo ">> MV-4 RJ2a, N = 256, x/y/z, no cleaning ($(basename "$EXE"), np $NP)"
W256=()
for axis in x y z; do w="$(case_run $axis 256 none)" || exit 1 ; W256+=("$w"); done
"$VENV_PY" "$ORACLE" "${W256[@]}" --l1-max "$L1_MAX_256" --dir-tol "$DIR_TOL" || FAILED=1
echo ">> MV-4 RJ2a, N = 512, x, no cleaning"
W512="$(case_run x 512 none)" || exit 1
"$VENV_PY" "$ORACLE" "$W512" --l1-max "$L1_MAX_512" || FAILED=1
sum256=$("$VENV_PY" "$ORACLE" "${W256[0]}" | sed -n 's/.*sum \([0-9.e+-]*\).*/\1/p')
sum512=$("$VENV_PY" "$ORACLE" "$W512" | sed -n 's/.*sum \([0-9.e+-]*\).*/\1/p')
if "$VENV_PY" -c "import sys; sys.exit(0 if $sum512 < $sum256 else 1)"; then
   echo "   convergence: L1 sum N=512 $sum512 < N=256 $sum256  PASS"
else
   echo "   convergence: L1 sum N=512 $sum512 not below N=256 $sum256  FAIL" ; FAILED=1
fi
find "$W512" -name '*.h5' -delete
echo ">> MV-4 RJ2a, N = 256, x, GLM against no cleaning"
WGLM="$(case_run x 256 glm)" || exit 1
"$VENV_PY" "$ORACLE" --pair "${W256[0]}" "$WGLM" || FAILED=1
for w in "${W256[@]}" "$WGLM"; do find "$w" -name '*.h5' -delete; done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-4 PASSED ($TAG)"
else
   echo "MV-4 FAILED ($TAG)"
   exit 1
fi
