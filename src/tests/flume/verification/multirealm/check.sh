#!/usr/bin/env bash
# FLUME verification of the multi-realm coupling (issue #37): 2-realm Sod against the single-realm Sod.
#
# Why: a forest manifest glues realms through inter-realm seams. sod-2realm.ini splits the sod-x domain at x = 0.5
# (the diaphragm, so the seam carries the Riemann problem from the first step) into two realms with the sod-x cell
# size; with a mirror seam filled at every Runge-Kutta stage (beta cadence) the seam is a block interface like any
# other, so the union of the two realms must reproduce the single-realm sod-x bit for bit. multirealm_oracle.py
# compares the interior cells by their centres (tolerance 0) and the sum of the realms' conservation histories with
# the single-realm one (round-off: the realms sum their cells separately).
#
# Leg 2 (seam + AMR): the same split with ni = 48 / 24 (the 2:1 refinement needs even block cells) and realm 2
# refined on x > 0.75, against the single-realm sod-amr with the same refinement. Realm 2 then holds both kinds of
# flux register face: the intra-realm 2:1 face at x = 0.75, crossed by the shock at t ~ 0.14, and the fine side of the
# inter-realm seam. Measured: CPU and FNL bitwise on all 135168 cells at t = 0.2 (216 steps), once the forest composed
# the two registrations (before, the manifest path registered the inter-realm faces only: the 2:1 face got no reflux).
#
# Measured (issue #37, np 2, t = 0.2, 174 steps): CPU and FNL bitwise on all 51200 cells, the summed conservation
# histories within 2e-13. The first runs found three defects: the CPU BC routine stopped on the forest's BC_SEAM crown
# rows; the fine side of the inter-realm reflux register was 2:1-restricted like an AMR seam, which wrote F_coarse - 0
# into three quarters of the realm-1 seam skin at every step; the FNL backend never copied the forest-built seam and
# BC maps to the device (illegal address in the seam fill kernel).
#
# Leg 3 (issue #40): sod-2realm-z.ini is the leg-1 split along z. On N >= 2 ranks each realm's blocks are split by z,
# so the two sides of the seam land on different ranks and every seam ghost crosses ranks (np 2: 3468 rows each way,
# no local row): the seam ghosts travel through the cross-rank seam exchange and the inter-realm reflux register is
# completed across ranks. The leg requires the union to reproduce sod-z bit for bit on one rank and on N ranks. Before
# issue #40 the N-rank run ended normally with a wrong solution (205 steps instead of 174, realm 1 mass frozen at 0.5,
# reflux mismatch 69).
#
# Leg 4 (issue #52): sod-amr-refined.ini splits sod-amr at its 2:1 face into a coarse realm (x < 0.75) and a fine realm
# (x > 0.75) glued by `coupling = refined`: the 2:1 seam ghosts (coarse->fine interpolation, fine->coarse 2x2x2 mean,
# the intra-realm formulas) and the 2:1 inter-realm reflux register. The union must reproduce the single-realm sod-amr
# of leg 2 bit for bit, and the same split turned to z (sod-amr-refined-z.ini, whose seam the partition splits across
# ranks) sod-amr-z; the same pair declared `mirror` must be refused at initialization for its different cell sizes,
# and so must the pair with the positivity limiter (refused on every multi-realm run: the inter-realm seam faces carry
# no limiting factor; the limiter's 2:1 seam flux synchronisation of issue #50 is intra-realm).
#
# Leg 5 (issue #51): the seam blocks of the two realms do not line up (same cells, different block sizes along the
# seam). sod-2realm-blocks.ini is the leg-1 split with realm 2 on 2 blocks per axis (each of its seam blocks faces 2x2
# blocks of realm 1, so its skin is scattered into four register faces); sod-2realm-blocks-z.ini the leg-3 split with
# realm 1, the register side, on the larger blocks (each register face holds 2x2 realm-2 faces; on N ranks the seam
# crosses ranks). Both must reproduce sod-x / sod-z bit for bit. Sod is 1-D, so these legs see the coverage of the
# register faces, not where the overlaps land along the face: the MHD rotor of mhd/multirealm/check.sh (leg 5) does.
#
# Leg 6 (issue #68 C2, C3, D1, D2): leg 4 (along x and z) with `[amr] seam_restriction = lagrange6`, then `weno6`,
# then `[amr] seam_ghost_fill = quintic`, `weno6`, `quintic-mp`, in every input: the coarse ghosts of the refined seam take the
# fine owner's 6-point restriction, the fine ghosts the coarse owner's 6-point fill, as the intra-realm 2:1 face does,
# so the unions must still reproduce the single-realm runs bit for bit (the shock crosses the seam: another operator
# on either side differs).
#
# Leg 7 (issue #110, P1): runtime AMR on a multi-realm forest. Leg 2's inputs with [amr] frequency = 4 in realm 2 (realm
# 1 has no marker, and FLUME refuses runtime regridding without one): its box marker re-marks the same blocks, so every
# regrid leaves the grid as it is, while the forest fills the seams
# before marking and, under ADAM_FOREST_REBUILD_TOPOLOGY=1, rebuilds the inter-realm topology after every step. The
# union must still reproduce the single-realm sod-amr of leg 2 bit for bit (the refusal of issue #74, D-M5-4, is gone).
#
# The rebuild itself is tested by running legs 1-6 again with ADAM_FOREST_REBUILD_TOPOLOGY=1 exported (a single-node
# mpirun or srun passes it to the ranks; across nodes Open MPI needs `-x ADAM_FOREST_REBUILD_TOPOLOGY`; the forest
# OR-reduces it, so one rank seeing it is enough): the forest then rebuilds the seam rows, the register faces, the BC_SEAM crown rows and their device
# copies after every step and refills every seam, and every leg must still reproduce its single-realm run bit for bit.
#
# Leg 8 (issue #110, P2): runtime AMR moving refined regions across an inter-realm seam. The single-realm Sod with a
# Loehner marker on the density (refine 0.3, derefine 0.1, a regrid every 4 steps, levels 2 to 3), against the same
# Sod split at the diaphragm x = 0.5 into two realms glued by `coupling = auto`, both regridding by the same marker.
# The realms use blocks of 48 x 8 x 8 cells one level lower (levels 1 to 2): the same x-slabs and the same cells as the
# single run (48 x 4 x 4), so the block-wise marks, the sibling families and every cell agree, and Sod is uniform in y
# and z, so every leaf along the seam votes alike. The forest balances the marks across the seam (one level per side,
# at most one level apart), which the seam switches between mirror and refined as the shock and the rarefaction cross
# it. With `[amr] seam_ghost_fill = injection` in every input the union reproduces the single-realm run bit for bit
# (measured: 221184 cells at step 330, regrids at steps 4, 216 and 240 in both). With the default tricubic fill the two
# differ by round-off (6e-14 to 2e-13 relative next to the 2:1 face at x = 0.25): one level lower, the realms' coarse
# blocks are twice as large in y and z, so the 4-point footprints of the 2:1 face sit elsewhere in their block and the
# same constant-in-y data is summed with other weights; that variant is checked to 1e-12 with the same regrid steps.
#
# Leg 9 (issue #110, P2): a mirror seam between realms whose equal cells sit at different levels (leg 5's forest, with
# 24 x 4 x 4 cells per block at levels 2 to 3 in realm 1 and 48 x 8 x 8 at levels 1 to 2 in realm 2: dx = 0.5/96 and
# dy = 1/16 on both sides), both regridding by the Loehner marker of leg 8. The balance must keep the seam cells equal
# (it compares cell sizes, not levels, which count from each realm's own root): the run completes with regrids, the
# seam stays a mirror (never logged as a 2:1 seam), and the mass of the two realms is conserved to round-off while
# the waves stay inside the box (t = 0.2).
#
# Leg 10 (issue #110, P4): a refinement criterion the program registers (`register_amr_criterion`), selected by
# `[amr_marker_1] mode = 6, name = box-criterion`. The test driver exe/test_flume_user_criterion (fobis mode
# test-flume-user-criterion-cpu-gnu; CRIT_EXE overrides it) is the FLUME CPU driver plus a criterion that votes as the
# built-in primitive-box marker does. The initial grid is not refined (amr_iterations = 0), so the first runtime regrid
# refines the box: the registered criterion must give the same run as the built-in marker, bit for bit, on one realm
# (sod-amr) and on the 2-realm forest of leg 2 (the box in realm 2). The same input with the plain FLUME driver, which
# registers nothing, must stop at initialisation. Not covered: a derefine vote of a criterion (with these inputs every
# derefine vote falls on the base level) and its combination with another marker; leg 8 covers that path with Loehner.
#
# Usage: ./check.sh [--build] [--np N] [--legs LIST]   LIST: leg1 ... leg10 or all (default all); leg4 compares against
#        the sod-amr run of leg2, leg5 against the sod-x and sod-z runs of leg1 and leg3, so they require those legs
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
BUILD=0
source "$VERIF_DIR/legs.sh"
legs_declare "leg1 leg2 leg3 leg4 leg5 leg6 leg7 leg8 leg9 leg10"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      --legs)  legs_select "$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N, --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
if has_leg leg4 && ! has_leg leg2; then
   echo "check.sh: leg4 compares against the sod-amr run of leg2: select leg2 too" >&2 ; exit 2
fi
if has_leg leg5 && ! { has_leg leg1 && has_leg leg3; }; then
   echo "check.sh: leg5 compares against the sod-x and sod-z runs of leg1 and leg3: select them too" >&2 ; exit 2
fi

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode flume-cpu-gnu)
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")-np$NP"

run() { # run <work-dir> <input> <files...>; RUN_NP overrides the rank count, RUN_MUST_FAIL=1 expects a failure
   local work="$1" input="$2" np="${RUN_NP:-$NP}" rc=0
   shift 2
   rm -rf "$work" ; mkdir -p "$work" ; cp "$@" "$work/"
   echo ">> $(basename "$work"): mpirun -np $np $(basename "$EXE") $input"
   (cd "$work" && mpirun -np "$np" "$EXE" "$input" > log.txt 2>&1) || rc=$?
   if [[ "${RUN_MUST_FAIL:-0}" == 1 ]]; then
      if [[ $rc -eq 0 ]]; then
         echo "check.sh: run succeeded but must fail, see $work/log.txt" >&2
         exit 1
      fi
   elif [[ $rc -ne 0 ]]; then
      echo "check.sh: run failed, see $work/log.txt" >&2
      exit 1
   fi
}

# l8_input <in> <out> <ni> <nj> <nk> <base level> <max level>: runtime AMR by a Loehner marker on the density in
# place of the input's markers, no initial AMR
l8_input() {
   python3 - "$@" <<'PYEOF'
import re
import sys

src, out, ni, nj, nk, base, top = sys.argv[1:8]
s = open(src).read()
s = re.sub(r"\n\[amr_marker_\d+\][^\[]*", "\n", s)  # drop the input's markers
lohner = ("[amr_marker_1]\nmode         = 4\ndelta_type   = max\ndelta_fine   = 0.0\ndelta_coarse = 0.0\n"
       "field        = 1\nvar          = 1\nrefine_tol   = 0.3\nderefine_tol = 0.1\nfloor        = 0.0\nbuffer       = 2\n\n")
s = s.replace("\n[field]", "\n" + lohner + "[field]", 1)
for key, val in (("ni", ni), ("nj", nj), ("nk", nk), ("iu_ref_levels", base), ("max_level", top),
              ("frequency", "4"), ("markers_number", "1"), ("amr_iterations", "0")):
 s, n = re.subn(rf"^{key}\s*=.*$", f"{key} = {val}", s, count=1, flags=re.M)
 assert n == 1, (src, key)
open(out, "w").write(s)
PYEOF
}

if has_leg leg1; then
   echo "== leg 1: 2-realm Sod vs sod-x"
   single="$CASE_DIR/work-$TAG-single"
   multi="$CASE_DIR/work-$TAG-2realm"
   run "$single" sod-x.ini "$VERIF_DIR/sod/sod-x.ini"
   run "$multi" sod-2realm.ini "$CASE_DIR/sod-2realm.ini" "$CASE_DIR/sod-2realm-r1.ini" "$CASE_DIR/sod-2realm-r2.ini"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
fi

if has_leg leg2; then
   echo "== leg 2: seam + AMR, 2-realm Sod with realm 2 refined on x > 0.75 vs sod-amr"
   single="$CASE_DIR/work-$TAG-amr-single"
   multi="$CASE_DIR/work-$TAG-amr-2realm"
   run "$single" sod-amr.ini "$CASE_DIR/sod-amr.ini"
   run "$multi" sod-amr-2realm.ini "$CASE_DIR/sod-amr-2realm.ini" "$CASE_DIR/sod-amr-2realm-r1.ini" \
       "$CASE_DIR/sod-amr-2realm-r2.ini"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
fi
if has_leg leg3; then
   echo "== leg 3: z-split 2-realm Sod vs sod-z, on one rank and on $NP (the seam crosses ranks, issue #40)"
   zfiles=("$CASE_DIR/sod-2realm-z.ini" "$CASE_DIR/sod-2realm-z-r1.ini" "$CASE_DIR/sod-2realm-z-r2.ini")
   single="$CASE_DIR/work-$TAG-z-single"
   run "$single" sod-z.ini "$VERIF_DIR/sod/sod-z.ini"
   multi="$CASE_DIR/work-$TAG-z-2realm"
   run "$multi" sod-2realm-z.ini "${zfiles[@]}"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   if [[ $NP -ge 2 ]]; then
      single="$CASE_DIR/work-$TAG-z-single-np1"
      RUN_NP=1 run "$single" sod-z.ini "$VERIF_DIR/sod/sod-z.ini"
      multi="$CASE_DIR/work-$TAG-z-2realm-np1"
      RUN_NP=1 run "$multi" sod-2realm-z.ini "${zfiles[@]}"
      "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   fi
fi
if has_leg leg4; then
   echo "== leg 4: refined (2:1) 2-realm Sod vs sod-amr (issue #52)"
   rfiles=("$CASE_DIR/sod-amr-refined.ini" "$CASE_DIR/sod-amr-refined-r1.ini" "$CASE_DIR/sod-amr-refined-r2.ini")
   multi="$CASE_DIR/work-$TAG-refined"
   run "$multi" sod-amr-refined.ini "${rfiles[@]}"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$CASE_DIR/work-$TAG-amr-single" --ngc 3 --tol 0
   single="$CASE_DIR/work-$TAG-amr-z-single"
   run "$single" sod-amr-z.ini "$CASE_DIR/sod-amr-z.ini"
   multi="$CASE_DIR/work-$TAG-refined-z"
   run "$multi" sod-amr-refined-z.ini "$CASE_DIR/sod-amr-refined-z.ini" "$CASE_DIR/sod-amr-refined-z-r1.ini" \
       "$CASE_DIR/sod-amr-refined-z-r2.ini"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   work="$CASE_DIR/work-$TAG-refined-as-mirror"
   mirror_dir="$CASE_DIR/work-$TAG-refined-as-mirror-src" ; rm -rf "$mirror_dir" ; mkdir -p "$mirror_dir"
   sed 's/^coupling *= *refined/coupling         = mirror/' "${rfiles[0]}" > "$mirror_dir/sod-amr-refined.ini"
   RUN_MUST_FAIL=1 run "$work" sod-amr-refined.ini "$mirror_dir/sod-amr-refined.ini" "${rfiles[1]}" "${rfiles[2]}"
   if ! grep -aq \
        'error stop forest_object%populate_inter_realm_topology: .*coupling = mirror joins cells of the same size' \
        "$work/log.txt"; then
      echo "check.sh: the refined pair declared mirror was not refused, see $work/log.txt" >&2
      exit 1
   fi
   echo "   $(basename "$work"): refused as expected"
   work="$CASE_DIR/work-$TAG-refined-limiter"
   limiter_dir="$CASE_DIR/work-$TAG-refined-limiter-src" ; rm -rf "$limiter_dir" ; mkdir -p "$limiter_dir"
   for f in "${rfiles[@]:1}"; do
      sed 's/^\[numerics\]/&\npositivity_limiter = cell/' "$f" > "$limiter_dir/$(basename "$f")"
   done
   RUN_MUST_FAIL=1 run "$work" sod-amr-refined.ini "${rfiles[0]}" "$limiter_dir"/sod-amr-refined-r?.ini
   if ! grep -aq 'positivity_limiter)=cell is not supported on multi-realm runs' "$work/log.txt"; then
      echo "check.sh: the positivity limiter on the refined pair was not refused, see $work/log.txt" >&2
      exit 1
   fi
   echo "   $(basename "$work"): refused as expected"
fi
if has_leg leg5; then
   echo "== leg 5: seam blocks that do not line up (issue #51)"
   bfiles=("$CASE_DIR/sod-2realm-blocks.ini" "$CASE_DIR/sod-2realm-blocks-r1.ini" "$CASE_DIR/sod-2realm-blocks-r2.ini")
   multi="$CASE_DIR/work-$TAG-blocks"
   run "$multi" sod-2realm-blocks.ini "${bfiles[@]}"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$CASE_DIR/work-$TAG-single" --ngc 3 --tol 0
   bfiles=("$CASE_DIR/sod-2realm-blocks-z.ini" "$CASE_DIR/sod-2realm-blocks-z-r1.ini"
           "$CASE_DIR/sod-2realm-blocks-z-r2.ini")
   multi="$CASE_DIR/work-$TAG-blocks-z"
   run "$multi" sod-2realm-blocks-z.ini "${bfiles[@]}"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$CASE_DIR/work-$TAG-z-single" --ngc 3 --tol 0
fi
if has_leg leg6; then
   for cfg in "seam_restriction = lagrange6|lagrange6 *restriction" "seam_restriction = weno6|weno6 *restriction" \
              "seam_ghost_fill = quintic|seam ghost fill regime: *+3" \
              "seam_ghost_fill = weno6|seam ghost fill regime: *+4" \
              "seam_ghost_fill = quintic-mp|seam ghost fill regime: *+5"; do
   key="${cfg%%|*}" ; pat="${cfg#*|}" ; rs="${key%% *}-${key##* }" # e.g. seam_ghost_fill-weno6
   echo "== leg 6: leg 4 with [amr] $key (issue #68 C2, C3, D1)"
   # the refined seam's coarse ghosts take the 6-point restriction of the fine owner (the codes and the regime packed in
   # the send rows' meta), its fine ghosts the quintic fill of the coarse owner (the 6-point anchor shifts in the meta),
   # the intra-realm 2:1 face the same: the unions must still reproduce the single-realm runs bit for bit, which they
   # cannot if either side uses another operator (the shock crosses the seam); the single run must log the regime
   l6_dir="$CASE_DIR/work-$TAG-$rs-src" ; rm -rf "$l6_dir" ; mkdir -p "$l6_dir"
   for f in sod-amr.ini sod-amr-z.ini sod-amr-refined-r1.ini sod-amr-refined-r2.ini sod-amr-refined-z-r1.ini \
            sod-amr-refined-z-r2.ini; do
      # the octree refines the null axes too, so their 4-cell blocks would not fit the 6-point footprint (fatal): 8
      sed -e "s/^\[amr\]/&\n$key/" -e 's/^\(n[ijk]\) *= *4 *$/\1 = 8/' \
          "$CASE_DIR/$f" > "$l6_dir/$f"
   done
   for axis in x z; do
      sfx="" ; [[ $axis == z ]] && sfx="-z"
      single="$CASE_DIR/work-$TAG-$rs-amr$sfx-single"
      run "$single" "sod-amr$sfx.ini" "$l6_dir/sod-amr$sfx.ini"
      if ! grep -aq "$pat" "$single/log.txt"; then
         echo "check.sh: [amr] $key is not active in $single/log.txt" >&2
         exit 1
      fi
      multi="$CASE_DIR/work-$TAG-$rs-refined$sfx"
      run "$multi" "sod-amr-refined$sfx.ini" "$CASE_DIR/sod-amr-refined$sfx.ini" "$l6_dir/sod-amr-refined$sfx-r1.ini" \
          "$l6_dir/sod-amr-refined$sfx-r2.ini"
      "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   done
   done
fi
if has_leg leg7 && ! has_leg leg2; then
   echo "check.sh: leg7 compares against the sod-amr run of leg2: select leg2 too" >&2 ; exit 2
fi
if has_leg leg7; then
   echo "== leg 7: leg 2 with runtime AMR in realm 2 ([amr] frequency = 4, issue #110): unchanged grid, bitwise"
   rt_dir="$CASE_DIR/work-$TAG-runtime-amr-src" ; rm -rf "$rt_dir" ; mkdir -p "$rt_dir"
   cp "$CASE_DIR/sod-amr-2realm-r1.ini" "$rt_dir/"
   sed 's/^frequency *= *0/frequency      = 4/' "$CASE_DIR/sod-amr-2realm-r2.ini" > "$rt_dir/sod-amr-2realm-r2.ini"
   grep -q '^frequency *= *4' "$rt_dir/sod-amr-2realm-r2.ini" || { echo "check.sh: no [amr] frequency" >&2 ; exit 1 ; }
   multi="$CASE_DIR/work-$TAG-runtime-amr"
   run "$multi" sod-amr-2realm.ini "$CASE_DIR/sod-amr-2realm.ini" "$rt_dir/sod-amr-2realm-r1.ini" \
       "$rt_dir/sod-amr-2realm-r2.ini"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$CASE_DIR/work-$TAG-amr-single" --ngc 3 --tol 0
fi
if has_leg leg8; then
   echo "== leg 8: runtime AMR across a coupling = auto seam (Loehner marker in both realms) vs the single realm (#110)"
   l8_dir="$CASE_DIR/work-$TAG-l8-src" ; rm -rf "$l8_dir" ; mkdir -p "$l8_dir"
   for fill in injection tricubic; do
   echo "   [amr] seam_ghost_fill = $fill"
   src="$l8_dir/$fill" ; mkdir -p "$src"
   l8_input "$CASE_DIR/sod-amr.ini" "$src/sod-amr.ini" 48 4 4 2 3
   l8_input "$CASE_DIR/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r1.ini" 48 8 8 1 2
   l8_input "$CASE_DIR/sod-amr-2realm-r2.ini" "$src/sod-amr-2realm-r2.ini" 48 8 8 1 2
   sed -i "/^\[amr\]/a seam_ghost_fill = $fill" "$src/sod-amr.ini" "$src/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r2.ini"
   sed 's/^coupling *= *mirror/coupling         = auto/' "$CASE_DIR/sod-amr-2realm.ini" > "$src/sod-amr-2realm.ini"
   grep -q '^coupling *= *auto' "$src/sod-amr-2realm.ini" || { echo "check.sh: no coupling in the manifest" >&2 ; exit 1 ; }
   single="$CASE_DIR/work-$TAG-l8-$fill-single"
   run "$single" sod-amr.ini "$src/sod-amr.ini"
   multi="$CASE_DIR/work-$TAG-l8-$fill-2realm"
   run "$multi" sod-amr-2realm.ini "$src/sod-amr-2realm.ini" "$src/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r2.ini"
   # the same regrid steps (each realm logs its own regrids; the single run, every one)
   steps_single=$(grep -ao 'regrid at step [0-9]*' "$single/log.txt" | sort -u | tr '\n' ' ')
   steps_multi=$(grep -ao 'regrid at step [0-9]*' "$multi/log.txt" | sort -u | tr '\n' ' ')
   if [[ -z "$steps_single" || "$steps_single" != "$steps_multi" ]]; then
      echo "check.sh: regrid steps differ: single [$steps_single], 2 realms [$steps_multi]" >&2 ; exit 1
   fi
   echo "   the same regrid steps: ${steps_single//regrid at step /}"
   tol=0 ; [[ $fill == tricubic ]] && tol=1e-12
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol $tol
   done
fi
if has_leg leg9; then
   echo "== leg 9: runtime AMR across a mirror seam between realms at different levels, equal cells (#110)"
   l9_dir="$CASE_DIR/work-$TAG-l9-src" ; rm -rf "$l9_dir" ; mkdir -p "$l9_dir"
   l8_input "$CASE_DIR/sod-2realm-blocks-r1.ini" "$l9_dir/sod-2realm-blocks-r1.ini" 24 4 4 2 3
   l8_input "$CASE_DIR/sod-2realm-blocks-r2.ini" "$l9_dir/sod-2realm-blocks-r2.ini" 48 8 8 1 2
   cp "$CASE_DIR/sod-2realm-blocks.ini" "$l9_dir/"
   multi="$CASE_DIR/work-$TAG-l9-2realm"
   run "$multi" sod-2realm-blocks.ini "$l9_dir/sod-2realm-blocks.ini" "$l9_dir/sod-2realm-blocks-r1.ini" \
       "$l9_dir/sod-2realm-blocks-r2.ini"
   n_regrid=$(grep -ac 'regrid at step' "$multi/log.txt" || true)
   if [[ $n_regrid -eq 0 ]]; then echo "check.sh: leg9 never regridded, see $multi/log.txt" >&2 ; exit 1 ; fi
   if grep -aq 'is a 2:1 seam' "$multi/log.txt"; then echo "check.sh: leg9's mirror seam turned 2:1" >&2 ; exit 1 ; fi
   python3 - "$multi" <<'PYEOF'
import glob
import sys

rows = []
for path in sorted(glob.glob(sys.argv[1] + "/*conservation_history*")):
    data = [line.split() for line in open(path).read().splitlines()[1:] if line.strip()]
    rows.append({int(r[0]): float(r[2]) for r in data})
steps = sorted(set.intersection(*(set(r) for r in rows)))
mass = [sum(r[s] for r in rows) for s in steps]
drift = max(abs(m - mass[0]) for m in mass) / abs(mass[0])
ok = len(rows) == 2 and drift < 1e-12
print(f"   mass of the 2 realms over {len(steps)} steps: largest relative drift {drift:.2e} (bound 1e-12) "
      f"{'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
   echo "   $(basename "$multi"): $n_regrid regrid log lines, the seam stayed a mirror"
fi
if has_leg leg10; then
   echo "== leg 10: a registered refinement criterion against the built-in box marker it reproduces (#110)"
   crit_exe="${CRIT_EXE:-$REPO_ROOT/exe/test_flume_user_criterion}"
   if [[ ! -x "$crit_exe" ]]; then
      echo "check.sh: $crit_exe not found (fobis build --mode test-flume-user-criterion-cpu-gnu, or CRIT_EXE)" >&2 ; exit 2
   fi
   l10_dir="$CASE_DIR/work-$TAG-l10-src" ; rm -rf "$l10_dir" ; mkdir -p "$l10_dir/builtin" "$l10_dir/criterion"
   for f in sod-amr.ini sod-amr-2realm-r2.ini; do
      # runtime AMR every 4 steps from an unrefined grid: the box marker acts at the first regrid
      sed -e 's/^frequency *= *0/frequency      = 4/' -e 's/^amr_iterations *= *[0-9]*/amr_iterations = 0/' \
          "$CASE_DIR/$f" > "$l10_dir/builtin/$f"
      python3 - "$l10_dir/builtin/$f" "$l10_dir/criterion/$f" <<'PYEOF'
import re
import sys

s = open(sys.argv[1]).read()
crit = ("[amr_marker_1]\nmode         = 6\nname         = box-criterion\ndelta_type   = max\ndelta_fine   = 0.0\n"
        "delta_coarse = 0.0\n\n")
s, n = re.subn(r"\[amr_marker_1\][^\[]*", crit, s, count=1)
assert n == 1
open(sys.argv[2], "w").write(s)
PYEOF
   done
   for d in builtin criterion; do
      cp "$CASE_DIR/sod-amr-2realm.ini" "$CASE_DIR/sod-amr-2realm-r1.ini" "$l10_dir/$d/"
   done
   exe_save="$EXE"
   for case in single 2realm; do
      if [[ $case == single ]]; then inp=sod-amr.ini ; files=("sod-amr.ini") ; else
         inp=sod-amr-2realm.ini ; files=("sod-amr-2realm.ini" "sod-amr-2realm-r1.ini" "sod-amr-2realm-r2.ini") ; fi
      EXE="$exe_save"
      run "$CASE_DIR/work-$TAG-l10-$case-builtin" "$inp" "${files[@]/#/$l10_dir/builtin/}"
      EXE="$crit_exe"
      run "$CASE_DIR/work-$TAG-l10-$case-criterion" "$inp" "${files[@]/#/$l10_dir/criterion/}"
      EXE="$exe_save"
      if ! grep -aq "regrid at step" "$CASE_DIR/work-$TAG-l10-$case-criterion/log.txt"; then
         echo "check.sh: leg10 ($case) never regridded" >&2 ; exit 1
      fi
      "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$CASE_DIR/work-$TAG-l10-$case-criterion" \
         "$CASE_DIR/work-$TAG-l10-$case-builtin" --ngc 3 --tol 0
   done
   # the plain FLUME driver registers nothing: the same input must stop at initialisation, naming the criterion
   work="$CASE_DIR/work-$TAG-l10-unregistered"
   RUN_MUST_FAIL=1 run "$work" sod-amr.ini "$l10_dir/criterion/sod-amr.ini"
   if ! grep -aq 'name = "box-criterion" is not a registered criterion' "$work/log.txt"; then
      echo "check.sh: an unregistered criterion was not refused, see $work/log.txt" >&2 ; exit 1
   fi
   echo "   $(basename "$work"): refused as expected"
fi
echo "multi-realm verification PASSED (legs: $LEGS) ($TAG)"
