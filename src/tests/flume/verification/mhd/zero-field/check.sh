#!/usr/bin/env bash
# FLUME MHD verification MV-2 (issue #41, M2-P3): the MHD path with B = 0, and a uniform MHD state across a 2:1 patch.
#
# Why: the MHD face kernels (7x7 Roe-Balsara core at the face average, decoupled B_n or the GLM pair) replace the
# Euler ones whenever physical_model = mhd-ideal. With B = 0 the MHD system is the Euler one plus inert fields, which
# gives exact expectations before any MHD-specific test (MV-4):
#   1. Sod along x, y, z through the MHD path (both variants): L1(rho) against the exact Riemann solution within the
#      M1 bound (not bitwise to the Euler path: at B = 0 the MHD eigenspaces are degenerate and WENO acts in a
#      different basis), the y and z runs equal the x run BITWISE after the axes permutation (sod_oracle.py), and B,
#      psi stay exactly zero;
#   2. a uniform MHD state (rho, u, p and B non-zero) on the periodic box of the conservation V3 case:
#      a. single level: the last checkpoint equals the step-0 state BITWISE (the MHD fluxes of equal states cancel);
#      b. with the 2:1 refined octant: equal to round-off (STEADY_TOL relative). Not bitwise: the coarse-fine ghost
#         fill reproduces a constant only to round-off for general values, and the Euler path shows the same (a
#         uniform Euler state with u, v, w = 0.8, 0.5, -0.3 drifts by 1.3e-15 on this grid, while the committed
#         all-ones state happens to be exact), so the criterion is the interface's, not the MHD kernels'.
# Inputs are derived from the committed Euler ones (make_input.py). One run at a time; checkpoints deleted after use.
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
L1_MAX="3.31e-03"     # the M1 V1 bound (verification/sod/check.sh)
STEADY_TOL="1.0e-14"  # uniform state across the 2:1 patch: a few ulp

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
MHD_ORACLE="$VERIF_DIR/mhd/plumbing/mhd_plumbing_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

run() { # run <work-dir> <ini>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: run failed, see $1/log.txt" >&2
      exit 1
   fi
}

for variant in none glm; do
   names=(bx by bz)
   [[ $variant == glm ]] && names+=(psi)
   echo ">> MV-2 Sod with B = 0 through the MHD path, divergence_control = $variant ($(basename "$EXE"), np $NP)"
   works=()
   for axis in x y z; do
      w="$CASE_DIR/work-$TAG-$variant-sod-$axis"
      rm -rf "$w" ; mkdir -p "$w"
      "$VENV_PY" "$CASE_DIR/make_input.py" "$VERIF_DIR/sod/sod-$axis.ini" "$w/sod-$axis.ini" --divergence-control $variant
      run "$w" "sod-$axis.ini"
      "$VENV_PY" "$MHD_ORACLE" --zero "$w" --names "${names[@]}" || FAILED=1
      works+=("$w")
   done
   "$VENV_PY" "$VERIF_DIR/sod/sod_oracle.py" "${works[0]}/sod-x.ini" "${works[@]}" --l1-max "$L1_MAX" || FAILED=1
   for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done
   for grid in single amr; do
      echo ">> MV-2 uniform MHD state, $grid grid, divergence_control = $variant"
      w="$CASE_DIR/work-$TAG-$variant-uniform-$grid"
      rm -rf "$w" ; mkdir -p "$w"
      extra=() ; tol=0.0
      if [[ $grid == single ]]; then extra=(--set amr.markers_number=0) ; else tol=$STEADY_TOL ; fi
      "$VENV_PY" "$CASE_DIR/make_input.py" "$VERIF_DIR/conservation/amr-periodic.ini" "$w/input.ini" \
         --divergence-control $variant --b 0.8 0.5 -0.3 --set initial_conditions.s=0.0 --set time.it_max=10 \
         --set IO.it_save=10 --set initial_conditions_region_1.u=0.3 --set initial_conditions_region_1.v=-0.2 \
         --set initial_conditions_region_1.w=0.1 "${extra[@]}"
      run "$w" input.ini
      "$VENV_PY" "$MHD_ORACLE" --steady "$w" --steady-tol $tol || FAILED=1
      find "$w" -name '*.h5' -delete
   done
done

if [[ $FAILED -eq 0 ]]; then
   echo "MV-2 PASSED ($TAG)"
else
   echo "MV-2 FAILED ($TAG)"
   exit 1
fi
