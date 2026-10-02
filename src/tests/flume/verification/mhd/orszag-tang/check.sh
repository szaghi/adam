#!/usr/bin/env bash
# FLUME MHD verification MV-12 (issue #41, M2-P6): the Orszag-Tang vortex, symmetry, conservation and positivity.
#
# Why: the standard 2-D nonlinear MHD test (Stone et al. 2008, section 8.4), shocks and a current sheet from a smooth
# start. 128^2, t = 0.5, GLM, np 2. Through orszag_tang_oracle.py:
#   - the 180 degrees rotational symmetry of the last checkpoint, relative defect at most SYM_TOL. The initial state is
#     bitwise symmetric (the IC evaluates odd functions of 2x - 1, 2y - 1); the scheme is symmetric in exact arithmetic
#     but its left- and right-biased reconstructions round differently, so the defect is a round-off seed amplified by
#     the flow: measured (CPU, M2-P6a) 1.7e-13 at step 100, 3e-11 at step 300 when the shocks form, 1.3e-9 at t = 0.5
#     (1281 steps), against a truncation-level defect (~1e-3) for a wrong fill or a sign error; w and bz stay exactly 0
#     on the CPU and carry round-off noise on FNL (~1e-12 of the in-plane field), scaled by the vector magnitude;
#   - the integrals of rho, rho u, E and B constant to CONS_TOL relative (measured 1.2e-14);
#   - positivity with zero floored cells: the floors are disabled, so a non-positive density or pressure stops the run.
# About 4.5 min on the CPU.
#
# --amr runs the AMR leg instead (M2-P7b): 32^2 base with [0.25, 0.75]^2 refined once (2:1 seams the shocks cross,
# octree with null z, nk = 4, see ../amr_box.py), t = 0.5, np 2, plus the uniform 32^2 run as reference:
#   - symmetry and conservation as above on the AMR run (measured 9.1e-13 and 1.1e-16: the refined box is symmetric
#     about the centre, and the reflux keeps the integrals);
#   - the seam div(B) (seam_divb_oracle.py): the peak of max|div B| beside the seam faces at most SEAM_REF_RATIO times
#     the peak of max|div B| of the uniform reference (measured 0.887 / 1.170 = 0.76; the shocks set both);
# about 6 min on the CPU.
#
# --divergence-control eglm runs the uniform leg with EGLM (issue #47, M3-P4c, EV-4), checkpoints every 25 steps:
#   - symmetry asserted twice, as for the rotor: at step 100 within EARLY_SYM_TOL (measured, CPU: 5.4e-13 splitting,
#     1.2e-12 weno-riemann HLLD) and at the end within 1e-4 (the EGLM sources amplify the round-off seed more than GLM
#     on the splitting scheme: 4.3e-7 against 5.9e-10 at t = 0.5; weno-riemann HLLD 1.1e-11);
#   - EV-4 (../eglm_conservation_oracle.py) instead of CONS_TOL: rho to round-off, the momentum, B and energy drifts
#     within twice the time integral of their source magnitudes (div B from the div(B) history, psi u.grad psi from the
#     checkpoints); a necessary bound, not a tight one (OT: energy drift 2.8e-3 of it; rho u and B ~1e-11, the 180
#     degrees symmetry cancels their sources).
# --amr is not combined with eglm (the EGLM sources at 2:1 seams are out of scope).
#
# Usage: ./check.sh [--np N] [--amr] [--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]] [--divergence-control glm|eglm]
#
# --numerics runs the legs on `scheme_space = weno-riemann` (mhd/numerics.sh, issue #47 M3-P3c): the bounds
# measured on the flux-splitting scheme are then not asserted (measurement), the scheme-independent checks are.
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
CONS_TOL="1.0e-13"
SEAM_REF_RATIO="1.0"
EARLY_SYM_TOL="1.0e-8"
AMR=0
DC="glm"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)  NP="$2" ; shift 2 ;;
      --amr) AMR=1 ; shift ;;
      --numerics) NUMERICS="$2" ; shift 2 ;;
      --divergence-control) DC="$2" ; shift 2 ;;
      *)     echo "check.sh: unknown argument '$1' (accepted: --np N, --amr, --numerics SPEC," \
                  "--divergence-control glm|eglm)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
numerics_check
if [[ ! $DC =~ ^(glm|eglm)$ || ( $DC == eglm && $AMR -eq 1 ) ]]; then
   echo "check.sh: --divergence-control '$DC' is not glm or eglm, or eglm with --amr" >&2
   exit 2
fi
TAG="$(basename "$EXE")-np$NP$(numerics_tag)"
[[ $DC == eglm ]] && TAG="$TAG-eglm"

# run_ot <work> <make_orszag_tang.py options>: write the input and run it.
run_ot() {
   local w="$1" ; shift
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$CASE_DIR/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$w/orszag-tang.ini" "$@"
   numerics_apply "$w/orszag-tang.ini"
   if ! (cd "$w" && mpirun -np "$NP" "$EXE" orszag-tang.ini > log.txt 2>&1); then
      echo "check.sh: Orszag-Tang run failed (a floor stop is a positivity failure), see $w/log.txt" >&2
      exit 1
   fi
}

STATUS=0
if [[ $AMR -eq 0 && $DC == eglm ]]; then
   echo ">> MV-12 / EV-4 Orszag-Tang 128^2, EGLM ($(basename "$EXE"), np $NP)"
   w="$CASE_DIR/work-$TAG"
   run_ot "$w" --cells 128 --divergence-control eglm --it-save 25
   "$VENV_PY" "$CASE_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$EARLY_SYM_TOL" --step 100 || STATUS=1
   "$VENV_PY" "$CASE_DIR/orszag_tang_oracle.py" "$w" --sym-tol "1.0e-4" || STATUS=1
   "$VENV_PY" "$CASE_DIR/../eglm_conservation_oracle.py" "$w" || STATUS=1
   works=("$w")
elif [[ $AMR -eq 0 ]]; then
   echo ">> MV-12 Orszag-Tang 128^2 ($(basename "$EXE"), np $NP)"
   w="$CASE_DIR/work-$TAG"
   run_ot "$w" --cells 128
   "$VENV_PY" "$CASE_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$SYM_TOL" --cons-tol "$CONS_TOL" || STATUS=1
   works=("$w")
else
   TAG="$TAG-amr"
   echo ">> MV-12 Orszag-Tang 32^2 + [0.25, 0.75]^2 refined 2:1 ($(basename "$EXE"), np $NP)"
   ref="$CASE_DIR/work-$TAG-ref" ; w="$CASE_DIR/work-$TAG"
   run_ot "$ref" --cells 32
   run_ot "$w" --cells 32 --refine-box 0.25 0.25 0.75 0.75
   "$VENV_PY" "$CASE_DIR/orszag_tang_oracle.py" "$w" --sym-tol "$SYM_TOL" --cons-tol "$CONS_TOL" || STATUS=1
   "$VENV_PY" "$CASE_DIR/../seam_divb_oracle.py" "$w" --ref "$ref" --ref-ratio-max "$SEAM_REF_RATIO" || STATUS=1
   works=("$ref" "$w")
fi
for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done

if [[ $STATUS -eq 0 ]]; then
   echo "MV-12 PASSED ($TAG)"
else
   echo "MV-12 FAILED ($TAG)"
   exit 1
fi
