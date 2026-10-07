#!/usr/bin/env bash
# FLUME MHD verification MV-15 (issue #46): 2:1 AMR seams on a quadtree.
#
# Why: a quadtree (ratio 4) refines x and y only, so its 2:1 seams are 2:1 in x and y and 1:1 in z. Until #46 the seam
# machinery (ghost restriction and interpolation, reflux, the limiter's seam flux synchronisation) assumed an octree:
# the coarse cells beside a seam picked up a spurious z dependence (M2-P7b: 1e-2 at the first step, negative pressure in
# Orszag-Tang at step 26), and quadtree AMR was refused at initialisation. The same 2-D problem now runs on three trees
# (../amr_box.py):
#   - octree with null z, nk = 4: the reference, the layout of every AMR verification before #46;
#   - quadtree, nk = 1: true 2-D, 4 times fewer cells per block;
#   - quadtree, nk = 4 with an active z axis (Orszag-Tang only): the z fluxes are computed, the solution must stay
#     z-invariant.
# quadtree_oracle.py asserts, per variable, the z invariance of every run and the agreement of the nk = 1 quadtree run
# with the reference at the same (x, y): to round-off, not bit for bit (an octree restriction averages 8 fine cells, two
# identical z layers, a quadtree 4; the summation orders differ, and even the reference's identical z layers differ at
# round-off, its tricubic weights depending on the z sub-position). The nk = 4 quadtree has its own time steps (its z
# cells enter the CFL condition): only its z invariance is asserted, bit for bit. Legs (CPU, np 2, measured):
#   ot     Orszag-Tang (../orszag-tang/make_orszag_tang.py, GLM), 32^2 base with [0.25, 0.75]^2 refined, t = 0.2:
#          quadtree nk = 1 against the octree 4.1e-12 (psi; 3.0e-13 or less elsewhere), octree z spread 5.7e-12,
#          quadtree nk = 4 z spread 0;
#   blast  the Balsara-Spicer blast (../positivity-probe/make_blast.py, EGLM), 32^2 base with [0.25, 0.75]^2 refined,
#          positivity limiter on, t = 0.006 (the outer shock crosses the seams): the limiter synchronises the seam flux
#          at every stage (2x1 fine faces per coarse face cell on a quadtree). Its switches amplify round-off: the
#          octree's identical z layers drift apart by up to 2.8e-5 (psi; r 9.3e-6), and the quadtree differs from the
#          reference by no more than that, per variable (ratios 0.40 to 0.80 on the CPU, 0.87 to 1.02 on FNL; same
#          limiter activity, 579 stages). Before #46 a quadtree seam gave 1e-2 at the first step.
# About 6 min on the CPU, 9 on FNL (WSL).
#
# Usage: ./check.sh [--np N] [--legs ot,blast]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEGS="ot,blast"
TOL="1.0e-10"               # agreement with the octree reference, relative (measured 4.1e-12, psi)
Z_TOL="1.0e-10"             # z invariance of the octree reference (measured 5.7e-12, psi)
BLAST_Z_TOL="1.0e-4"        # z invariance of the blast reference (measured 2.8e-5, psi)
BLAST_SPREAD_FACTOR="2.0"   # blast agreement over the reference z spread (measured 0.40-0.80 CPU, 0.87-1.02 FNL)
KEEP=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) LEGS="$2" ; shift 2 ;;
      --keep) KEEP=1 ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs ot,blast, --keep)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

has_leg() { [[ ",$LEGS," == *",$1,"* ]]; }
run() { # run <work> <ini>: run one case, fail the check if the run fails
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: $(basename "$1") failed, see $1/log.txt" >&2
      FAILED=1
      return 1
   fi
}

if has_leg ot; then
   echo ">> MV-15 Orszag-Tang 32^2 + [0.25, 0.75]^2 refined 2:1, octree vs quadtree ($(basename "$EXE"), np $NP)"
   works=()
   for tree in octree quad-nk1 quad-nk4; do
      w="$CASE_DIR/work-$TAG-ot-$tree"
      rm -rf "$w" ; mkdir -p "$w"
      case $tree in
         octree)   opts=(--ratio 8) ;;
         quad-nk1) opts=(--ratio 4) ;;
         quad-nk4) opts=(--ratio 4 --nk 4) ;;
      esac
      "$VENV_PY" "$VERIF_DIR/mhd/orszag-tang/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" \
         "$w/orszag-tang.ini" --cells 32 --time-max 0.2 --refine-box 0.25 0.25 0.75 0.75 "${opts[@]}"
      run "$w" orszag-tang.ini && works+=("$w")
   done
   if [[ ${#works[@]} -eq 3 ]]; then
      "$VENV_PY" "$CASE_DIR/quadtree_oracle.py" "${works[0]}" "${works[1]}" --ngc 3 --tol "$TOL" --z-tol "$Z_TOL" \
         | sed 's/^/   /' || FAILED=1
      "$VENV_PY" "$CASE_DIR/quadtree_oracle.py" "${works[2]}" --ngc 3 --z-tol 0 | sed 's/^/   /' || FAILED=1
   fi
   [[ $KEEP -eq 1 ]] || for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done
fi

if has_leg blast; then
   echo ">> MV-15 Balsara-Spicer blast 32^2 + [0.25, 0.75]^2 refined 2:1, limiter on, octree vs quadtree"
   works=()
   for tree in octree quad-nk1; do
      w="$CASE_DIR/work-$TAG-blast-$tree"
      rm -rf "$w" ; mkdir -p "$w"
      "$VENV_PY" "$VERIF_DIR/mhd/positivity-probe/make_blast.py" "$w/blast.ini" --eglm --cells 32 --time-max 0.006 \
         --refine-box 0.25 0.25 0.75 0.75 --ratio "$([[ $tree == octree ]] && echo 8 || echo 4)" > /dev/null
      sed -i -E 's/^(reflux\s*=.*)$/\1\npositivity_limiter       = cell/' "$w/blast.ini"
      grep -q '^positivity_limiter *= *cell' "$w/blast.ini"
      run "$w" blast.ini && works+=("$w")
   done
   if [[ ${#works[@]} -eq 2 ]]; then
      "$VENV_PY" "$CASE_DIR/quadtree_oracle.py" "${works[@]}" --ngc 3 --tol "$TOL" --z-tol "$BLAST_Z_TOL" \
         --spread-factor "$BLAST_SPREAD_FACTOR" | sed 's/^/   /' || FAILED=1
      for w in "${works[@]}"; do
         n=$(grep -ac 'positivity limiter:' "$w/log.txt" || true)
         echo "   $(basename "$w"): $n limiter log lines (stages with limited faces)"
      done
   fi
   [[ $KEEP -eq 1 ]] || for w in "${works[@]}"; do find "$w" -name '*.h5' -delete; done
fi

if [[ $FAILED -eq 0 ]]; then
   echo "MV-15 PASSED ($TAG)"
else
   echo "MV-15 FAILED ($TAG)"
   exit 1
fi
