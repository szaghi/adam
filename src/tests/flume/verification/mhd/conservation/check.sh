#!/usr/bin/env bash
# FLUME MHD verification MV-11 (issue #41, M2-P7a): conservation across AMR coarse-fine faces with reflux, MHD.
#
# Why: the V3 case (verification/conservation: triply periodic box, one refined octant, 6 coarse-fine faces, uniform
# flow with a 1% seeded density/pressure perturbation) through the MHD path with the uniform field B = (0.8, 0.5, -0.3):
# the perturbation drives the induction, so the flux register carries the B and psi rows at every seam. Legs:
#   1. GLM without damping (glm_alpha = 0), reflux on : all 9 volume integrals (rho, rho u, E, each of B, psi) constant
#      within MAX_DRIFT; psi starts at zero, so its drift is relative to the largest initial integral;
#   2. the same, reflux off: the negative control must drift by at least MIN_DRIFT (the seams are exercised);
#   3. GLM with damping (glm_alpha = 0.18), reflux on: the 8 physical integrals within MAX_DRIFT; psi is reported, not
#      bounded. The damping source -k psi couples the Runge-Kutta stages: the seam-flux mismatch of stage s reaches
#      psi^(n+1) through the later stages' damping, so it enters with weight b_s (1 - O(k dt)) while the register
#      corrects b_s; the remainder is O(k dt) of the uncorrected leak (measured, M2-P7a: int psi 2.7e-11, 2.7e-10,
#      1.3e-9 for glm_alpha 0.018, 0.18, 0.9, linear; without the seams it stays at 1e-21). With damping int psi is not
#      a conserved quantity in the first place;
#   4. divergence control none (nv = 8), reflux on: the 8 integrals within MAX_DRIFT.
# Every run must register the 6 coarse-fine faces. Measured (M2-P7a, CPU): drift <= 4.7e-14 with reflux (B_z the
# largest), 7e-6 to 1.7e-5 without. Inputs are derived from the committed V3 one (mhd/zero-field/make_input.py).
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
MAX_DRIFT="1.0e-13"
MIN_DRIFT="1.0e-10"

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
ORACLE="$VERIF_DIR/conservation/conservation_oracle.py"
TAG="$(basename "$EXE")-np$NP"

case_run() { # case_run <name> <none|glm> <glm_alpha> <reflux true|false>: run one leg, print its work directory
   local w="$CASE_DIR/work-$TAG-$1" extra=()
   [[ "$2" == glm ]] && extra=(--set "mhd.glm_alpha=$3")
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$VERIF_DIR/mhd/zero-field/make_input.py" "$VERIF_DIR/conservation/amr-periodic.ini" "$w/input.ini" \
      --divergence-control "$2" --b 0.8 0.5 -0.3 "${extra[@]}" --set "numerics.reflux=.$4."
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      echo "check.sh: MV-11 leg $1 run failed, see $w/log.txt" >&2
      return 1
   fi
   if ! grep -aq "registered intra-realm AMR seam faces: +6" "$w/log.txt"; then
      echo "check.sh: MV-11 leg $1 must register 6 coarse-fine faces, see $w/log.txt" >&2
      return 1
   fi
   echo "$w"
}

echo ">> MV-11 MHD conservation across 2:1 faces ($(basename "$EXE"), np $NP)"
glm=$(case_run glm-reflux glm 0.0 true) || exit 1
leaky=$(case_run glm-noreflux glm 0.0 false) || exit 1
damped=$(case_run glm-damped-reflux glm 0.18 true) || exit 1
none=$(case_run none-reflux none 0 true) || exit 1
STATUS=0
"$VENV_PY" "$ORACLE" --conserved "$glm" --max-drift "$MAX_DRIFT" --leaky "$leaky" --min-drift "$MIN_DRIFT" || STATUS=1
"$VENV_PY" "$ORACLE" --conserved "$damped" --max-drift "$MAX_DRIFT" --exclude psi || STATUS=1
"$VENV_PY" "$ORACLE" --conserved "$none" --max-drift "$MAX_DRIFT" || STATUS=1
for w in "$glm" "$leaky" "$damped" "$none"; do find "$w" -name '*.h5' -delete; done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-11 PASSED ($TAG)"
else
   echo "MV-11 FAILED ($TAG)"
   exit 1
fi
