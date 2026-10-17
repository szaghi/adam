#!/usr/bin/env bash
# FLUME verification AS (issue #68, M6-A3): a shock crossing a static 2:1 face, against the exact solution.
#
# Why: the committed shock-across-a-seam checks (multirealm legs, sod-amr goldens) compare FLUME runs with each other,
# so an error common to all of them passes. Here the Sod (V1) and RJ2a (MV-4) problems run on the V1 grid with one 2:1
# face at x = 0.75, which the right shock (RJ2a: the fast shock) crosses well before t = 0.2; make_amr_shock.py writes
# the cases, amr_shock_oracle.py compares them with the exact solution: the mass defect, the mass-equivalent shock
# position, and the error of the constant state behind the shock (the plateau, where the face sits at t = 0.2) and of
# the cells beside the face (the band).
#
# Legs (each: c2f, the shock goes from the coarse to the fine side; f2c, from the fine to the coarse side; and the
# uniform controls at the coarse spacing N and at the fine spacing 2N):
#   sod   N = 200 (V1's grid);
#   rj2a  N = 256 (MV-4's grid), no divergence control (1-D: B_n is uniform).
# Asserted on every refined run: the mass defect at most MASS_MAX (the reflux conserves; without it the defect is
# 1e-4; the floor is not round-off on RJ2a N = 256, 2.9e-10 on the uniform control as well: numerical precursors of
# the left fast rarefaction reach the extrapolated end x = 0 at that level, 9e-15 at N = 512), the shock position
# within SHOCK_MAX cells of the local spacing, the transverse spread at round-off (SPREAD_MAX), and every plateau and
# band metric at most the baseline plus 10 % (the bounds below; the baselines of the seam remedies of #68 C, D, E: an
# improvement passes, a degradation fails).
# --long adds the refined runs at 2N against the controls at 2N and 4N: a report of how the seam error scales.
# Cost, 2 ranks: sod about 7 min on the CPU, 8 min on FNL; rj2a about 25 min on FNL and hours on the CPU (the octree
# refines the null y and z too: 88 and 200 blocks against 32), so run the rj2a leg on FNL.
#
# Usage: ./check.sh [--np N] [--legs sod,rj2a|all] [--long] [--keep]   (default: sod,rj2a; --leg L is an older spelling
#        of --legs, repeatable)
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
source "$VERIF_DIR/legs.sh"
legs_declare "sod rj2a" "sod rj2a"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NP=2
KEEP=0
LONG=0
MASS_MAX="1.0e-8"
SHOCK_MAX="0.5"
SPREAD_MAX="1.0e-12"
declare -A CELLS=([sod]=200 [rj2a]=256)
# the baselines plus 10 % (M6-A3; Sod: CPU and FNL agree to the printed digits; RJ2a: FNL), per leg and side
declare -A BOUNDS=(
   [sod-c2f]="plateau_r:2.8e-2 plateau_u:5.0e-3 plateau_p:7.7e-3 plateau_s:2.0e-2 plateau_s_mean:2.6e-3 band_r:7.0e-3
              band_u:5.0e-3 band_p:7.7e-3 band_s:1.6e-4"
   [sod-f2c]="plateau_r:1.8e-2 plateau_u:7.4e-3 plateau_p:1.2e-2 plateau_s:1.3e-2 plateau_s_mean:2.5e-3 band_r:3.9e-3
              band_u:2.7e-3 band_p:4.1e-3 band_s:1.2e-4"
   [rj2a-c2f]="plateau_r:1.6e-2 plateau_u:8.9e-4 plateau_p:6.9e-3 plateau_s:6.5e-3 plateau_s_mean:9.6e-4 band_r:1.2e-4
               band_u:6.3e-5 band_p:1.6e-4 band_s:2.1e-5"
   [rj2a-f2c]="plateau_r:2.1e-2 plateau_u:1.0e-3 plateau_p:7.2e-3 plateau_s:8.5e-3 plateau_s_mean:1.1e-3 band_r:4.8e-4
               band_u:6.1e-5 band_p:3.8e-4 band_s:4.9e-5"
)
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) legs_select "$2" ; shift 2 ;;
      --leg)  legs_select "$2" ; shift 2 ;;
      --long) LONG=1 ; shift ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs L, --leg L, --long, --keep)" >&2
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

run() { # run <work> <case> [make_amr_shock options]: write the input, run it, report the wall time
   local work=$1 case=$2 ; shift 2
   rm -rf "${work:?}" ; mkdir -p "$work"
   "$VENV_PY" "$CASE_DIR/make_amr_shock.py" "$case" "$work/$case.ini" "$@"
   local t0=$SECONDS
   if ! (cd "$work" && mpirun -np "$NP" "$EXE" "$case.ini" < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $work/log.txt" ; FAILED=1 ; return 1
   fi
   echo "   $(basename "$work"): $(( SECONDS - t0 )) s"
}

oracle() { # oracle <control> <refined> [bounds]: compare a refined run with its control, assert the bounds
   local control=$1 refined=$2 bounds=${3:-}
   local args=(--max "mass:$MASS_MAX" --max "shock:$SHOCK_MAX" --max "spread:$SPREAD_MAX")
   for b in $bounds ; do args+=(--max "$b") ; done
   "$VENV_PY" "$CASE_DIR/amr_shock_oracle.py" "$refined" --control "$control" "${args[@]}" || FAILED=1
}

echo ">> AS: a shock crossing a static 2:1 face ($(basename "$EXE"), np $NP)"
for leg in $LEGS; do
   ladder=("${CELLS[$leg]}") ; [[ $LONG -eq 1 ]] && ladder+=($(( 2 * CELLS[$leg] )))
   for n in "${ladder[@]}"; do
      echo "-- $leg, N = $n"
      pre="$CASE_DIR/work-$TAG-$leg"
      run "$pre-uniform-n$n" "$leg" --cells "$n" --uniform || continue
      run "$pre-uniform-n$(( 2 * n ))" "$leg" --cells "$(( 2 * n ))" --uniform || continue
      for side in c2f f2c; do
         run "$pre-$side-n$n" "$leg" --cells "$n" --side "$side" || continue
         # the shock is on the refined side at t = 0.2 for c2f, on the coarse side for f2c
         control="$pre-uniform-n$n" ; [[ $side == c2f ]] && control="$pre-uniform-n$(( 2 * n ))"
         bounds="" ; [[ $n -eq ${CELLS[$leg]} ]] && bounds="${BOUNDS[$leg-$side]}"
         oracle "$control" "$pre-$side-n$n" "$bounds"
      done
   done
done
if [[ $KEEP -eq 0 ]]; then
   find "$CASE_DIR" -path "$CASE_DIR/work-$TAG-*" -name '*.h5' -delete
fi

if [[ $FAILED -eq 0 ]]; then
   echo "AS PASSED ($TAG, legs: $LEGS)"
else
   echo "AS FAILED ($TAG, legs: $LEGS)"
   exit 1
fi
