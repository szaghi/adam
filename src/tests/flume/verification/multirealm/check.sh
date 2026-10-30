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
#   leg11  (issue #110 P3) mixed levels along a seam: leg 8's grids, the markers two boxes next to the seam, one in each
#          realm over half of it (realm 2: x in [0.5, 0.75], y < 0.5; realm 1: x in [0.25, 0.5], y > 0.5), so the
#          seam holds 2:1 faces in both directions (realm 2 fine below y = 0.5, realm 1 fine above), and runtime
#          regridding on (frequency 4: the box votes keep the grid, the per-face balance and the mixed-seam topology
#          run every 4 steps). Against the single realm with the same boxes: bitwise with seam_ghost_fill = injection,
#          within 1e-12 with tricubic (as leg 8); and with the topology rebuilt every step (the hook), bitwise.
#   leg12  (issue #110 P3) a feature moving along part of a seam, the balance face by face: the isentropic vortex of
#          the vortex-periodic regression case (radius 0.07) started at (0.47, 0.3), across the seam, in a free stream
#          (0, 1) along it to t = 0.3, the box [0, 1]^2 x [0, 0.5] split at x = 0.5 into two coupling = auto realms (an
#          octree, 8x16x8 cells per block, cubic cells of 1/32 at base level 1, max_level 3, y periodic, x inviscid
#          walls, z extrapolated with w = 0: a closed box), a Loehner marker on the density in both, a regrid every 5 steps. The
#          refined patch covers part of the seam and moves along it, so the seam turns mixed. No single realm
#          reproduces these blocks (the realm roots are not square), so the leg asserts invariants: the seam mixed at
#          some regrid, the grid refined and coarsened, the mass of the two realms constant to 1e-12, the run on 1 and
#          on 2 ranks bitwise, and the run with the topology rebuilt every step (the hook) bitwise.
#   leg13  (issue #110 P5) restart across multi-realm regrids: leg 12's case, 195 steps against 155 + restart + 40. The
#          restart is saved at step 155, a regrid that coarsens and leaves 2:1 faces on the seam, so the restarted run
#          builds a mixed seam from the restored trees, and it regrids again into a new grid at step 190. The last
#          fields bitwise and every history file identical (the rg4 pattern on two realms); a restart on fewer ranks
#          than saved stops (issue #109). The seams are stage_coincident: an end_of_step seam carries the seam ghosts of
#          the previous step's stage buffers, which the restart files do not hold, so it does not restart bitwise.
#   leg14  (issue #113) the end_of_step (alpha) cadence: leg 12's forest without regridding (the base level only, the
#          seam a 1:1 face everywhere), the vortex moving along the seam, 120 steps at CFL 0.4 (about t = 0.3) and
#          240 at CFL 0.2. Alpha fills the seam ghosts at the top of every step and all the stages of the step read
#          that state, so it lags the stage_coincident (beta) run by a first-order-in-time error: the largest relative
#          difference between the alpha and the beta run is bounded (L14_BOUND at CFL 0.4) and halves with the time
#          step (the CFL 0.4 / CFL 0.2 ratio in L14_RATIO). Alpha's restart is bitwise: 40 steps against 25 + restart
#          + 15, the last fields and every history file. Before #113 the stages read seam ghosts no fill reached.
#   leg15  (issue #71) inter-realm 2:1 seams on quadtrees: leg 4's inputs with [amr] ratio = 4 (z not refined), split
#          along x, then turned to y (the Sod along y, the split at y = 0.75). Along z the 2:1 seam ghosts copy the
#          coincident cell, the coarse ghosts are the mean of 2x2 fine cells, and the register face holds two fine
#          blocks: each union must reproduce its single quadtree realm bit for bit. The z split (a seam normal to the
#          axis the quadtree does not refine) must be refused at initialization as a 2:1 seam, and accepted as a 1:1
#          one: leg 3's z split on quadtrees with coupling = auto resolves to a mirror (its 25 cells along z are odd,
#          which an octree refuses) and must reproduce the single quadtree sod-z bit for bit.
#
# Usage: ./check.sh [--build] [--np N] [--legs LIST]   LIST: leg1 ... leg15 or all (default all); leg4 compares against
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
legs_declare "leg1 leg2 leg3 leg4 leg5 leg6 leg7 leg8 leg9 leg10 leg11 leg12 leg13 leg14 leg15"

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

# l11_input <in> <out> <ni> <nj> <nk> <base level> <max level> <target level>: two box markers next to the seam
# (leg 11), runtime regridding every 4 steps.
l11_input() {
   python3 - "$@" <<'PYEOF'
import re
import sys

src, out, ni, nj, nk, base, top, target = sys.argv[1:9]
s = open(src).read()
s = re.sub(r"\n\[amr_marker_\d+\][^\[]*", "\n", s)  # drop the input's markers
box = ("[amr_marker_{m}]\nmode         = 1\ngeo_type     = primitive-box\ndelta_type   = max\ndelta_fine   = 0.0\n"
       "delta_coarse = 0.0\nbox_xmin     = {x0}\nbox_ymin     = {y0}\nbox_zmin     = 0.0\nbox_xmax     = {x1}\n"
       "box_ymax     = {y1}\nbox_zmax     = 1.0\ntarget_level = {t}\n\n")
boxes = (box.format(m=1, x0=0.5, x1=0.75, y0=0.0, y1=0.5, t=target) +
         box.format(m=2, x0=0.25, x1=0.5, y0=0.5, y1=1.0, t=target))
s = s.replace("\n[field]", "\n" + boxes + "[field]", 1)
for key, val in (("ni", ni), ("nj", nj), ("nk", nk), ("iu_ref_levels", base), ("max_level", top),
                 ("frequency", "4"), ("markers_number", "2"), ("amr_iterations", "1")):
    s, n = re.subn(rf"^{key}\s*=.*$", f"{key} = {val}", s, count=1, flags=re.M)
    assert n == 1, (src, key)
open(out, "w").write(s)
PYEOF
}

# l12_inputs <dir>: the manifest and the two realm inputs of leg 12 (the vortex moving along part of the seam).
l12_inputs() {
   python3 "$CASE_DIR/make_vortex_forest.py" "$1"
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
if has_leg leg11; then
   echo "== leg 11: mixed levels along a coupling = auto seam (a box next to it in each realm) vs the single realm (#110 P3)"
   l11_dir="$CASE_DIR/work-$TAG-l11-src" ; rm -rf "$l11_dir" ; mkdir -p "$l11_dir"
   for fill in injection tricubic hook; do
   echo "   [amr] seam_ghost_fill = ${fill/hook/injection, topology rebuilt every step}"
   src="$l11_dir/$fill" ; mkdir -p "$src"
   l11_input "$CASE_DIR/sod-amr.ini" "$src/sod-amr.ini" 48 4 4 2 3 3
   l11_input "$CASE_DIR/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r1.ini" 48 8 8 1 2 2
   l11_input "$CASE_DIR/sod-amr-2realm-r2.ini" "$src/sod-amr-2realm-r2.ini" 48 8 8 1 2 2
   sed -i "/^\[amr\]/a seam_ghost_fill = ${fill/hook/injection}" "$src/sod-amr.ini" "$src/sod-amr-2realm-r1.ini" \
       "$src/sod-amr-2realm-r2.ini"
   sed 's/^coupling *= *mirror/coupling         = auto/' "$CASE_DIR/sod-amr-2realm.ini" > "$src/sod-amr-2realm.ini"
   grep -q '^coupling *= *auto' "$src/sod-amr-2realm.ini" || { echo "check.sh: no coupling in the manifest" >&2 ; exit 1 ; }
   single="$CASE_DIR/work-$TAG-l11-$fill-single"
   run "$single" sod-amr.ini "$src/sod-amr.ini"
   multi="$CASE_DIR/work-$TAG-l11-$fill-2realm"
   if [[ $fill == hook ]]; then
      ADAM_FOREST_REBUILD_TOPOLOGY=1 run "$multi" sod-amr-2realm.ini "$src/sod-amr-2realm.ini" \
         "$src/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r2.ini"
   else
      run "$multi" sod-amr-2realm.ini "$src/sod-amr-2realm.ini" "$src/sod-amr-2realm-r1.ini" "$src/sod-amr-2realm-r2.ini"
   fi
   if ! grep -aq 'mixed levels along the seam' "$multi/log.txt"; then
      echo "check.sh: the seam is not mixed: see $multi/log.txt" >&2 ; exit 1
   fi
   grep -a -m1 'mixed levels along the seam' "$multi/log.txt" | sed 's/^[^f]*forest: /   /'
   tol=0 ; [[ $fill == tricubic ]] && tol=1e-12
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol $tol
   done
fi
if has_leg leg12; then
   echo "== leg 12: a vortex moving along part of a coupling = auto seam, both realms regridding (#110 P3)"
   l12_dir="$CASE_DIR/work-$TAG-l12-src" ; rm -rf "$l12_dir" ; mkdir -p "$l12_dir"
   l12_inputs "$l12_dir"
   l12_files=("$l12_dir/vortex-2realm.ini" "$l12_dir/vortex-r1.ini" "$l12_dir/vortex-r2.ini")
   for variant in np2 np1 hook; do
      multi="$CASE_DIR/work-$TAG-l12-$variant"
      case $variant in
         np1)  RUN_NP=1 run "$multi" vortex-2realm.ini "${l12_files[@]}" ;;
         hook) ADAM_FOREST_REBUILD_TOPOLOGY=1 run "$multi" vortex-2realm.ini "${l12_files[@]}" ;;
         *)    run "$multi" vortex-2realm.ini "${l12_files[@]}" ;;
      esac
   done
   multi="$CASE_DIR/work-$TAG-l12-np2"
   # every rank prints the topology line: count rank 0's
   n_mixed=$(grep -ac '^\[mpi-00000\].*mixed levels along the seam' "$multi/log.txt" || true)
   n_regrid=$(grep -ac 'regrid at step' "$multi/log.txt" || true)
   n_refined=$(grep -aEc 'regrid at step .* [1-9][0-9]* refined' "$multi/log.txt" || true)
   n_coarsened=$(grep -aEc 'regrid at step .* [1-9][0-9]* coarsened' "$multi/log.txt" || true)
   if [[ $n_mixed -eq 0 ]]; then echo "check.sh: leg12's seam never turned mixed, see $multi/log.txt" >&2 ; exit 1 ; fi
   if [[ $n_regrid -eq 0 ]]; then echo "check.sh: leg12 never regridded, see $multi/log.txt" >&2 ; exit 1 ; fi
   if [[ $n_refined -eq 0 ]]; then echo "check.sh: leg12 never refined, see $multi/log.txt" >&2 ; exit 1 ; fi
   if [[ $n_coarsened -eq 0 ]]; then echo "check.sh: leg12 never coarsened, see $multi/log.txt" >&2 ; exit 1 ; fi
   echo "   $n_regrid realm regrids ($n_refined refining, $n_coarsened coarsening), the seam mixed at $n_mixed topology builds"
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
   for variant in np1 hook; do
      "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$CASE_DIR/work-$TAG-l12-$variant" "$multi" --ngc 3 --tol 0
   done
fi

if has_leg leg13; then
   echo "== leg 13: restart across multi-realm regrids, 195 steps vs 155 + restart + 40 (#110 P5)"
   l13_src="$CASE_DIR/work-$TAG-l13-src" ; rm -rf "$l13_src" ; mkdir -p "$l13_src"
   l12_inputs "$l13_src"
   # it_max > 0: the step count, not time_max, ends the runs
   sed -i 's/^it_max *=.*/it_max   = 195/' "$l13_src"/vortex-r[12].ini
   l13_files=("$l13_src/vortex-2realm.ini" "$l13_src/vortex-r1.ini" "$l13_src/vortex-r2.ini")
   a="$CASE_DIR/work-$TAG-l13-A" ; b="$CASE_DIR/work-$TAG-l13-B"
   run "$a" vortex-2realm.ini "${l13_files[@]}"
   sed -i 's/^it_max *=.*/it_max   = 155/; s/^restart_save *=.*/restart_save           = 155/' "$l13_src"/vortex-r[12].ini
   run "$b" vortex-2realm.ini "${l13_files[@]}"
   # the topology line printed last before the step-155 regrid line is the one its rebuild built
   seam155=$(awk '/mixed levels along the seam|is a 2:1 seam/ {t = $0} /flume: regrid at step 155:/ {print t ; exit}' \
             "$b/log.txt")
   if grep -aq 'flume: regrid at step 155: .* [1-9][0-9]* coarsened' "$b/log.txt" &&
      [[ $seam155 =~ "and "[1-9][0-9]*" 2:1 faces" ]]; then
      echo "   step 155 coarsens and leaves 2:1 faces on the seam"
   else
      echo "check.sh: leg13's step 155 does not coarsen or leaves no 2:1 face on the seam, see $b/log.txt" >&2
      exit 1
   fi
   mv "$b/log.txt" "$b/log-1.txt"
   sed -i -e 's/^restart *=.*/restart                = .true./' -e 's/^it_max *=.*/it_max   = 195/' \
          -e 's/^restart_save *=.*/restart_save           = 0/' "$b"/vortex-r[12].ini
   # a restart on fewer ranks than saved stops (issue #109); run on a copy, the restart files stay untouched
   if [[ $NP -gt 1 ]]; then
      fewer="$CASE_DIR/work-$TAG-l13-fewer" ; rm -rf "$fewer" ; cp -r "$b" "$fewer"
      echo ">> $(basename "$fewer"): mpirun -np $((NP - 1)) $(basename "$EXE") vortex-2realm.ini (must stop)"
      if (cd "$fewer" && mpirun -np $((NP - 1)) "$EXE" vortex-2realm.ini > log.txt 2>&1); then
         echo "check.sh: leg13's restart on $((NP - 1)) ranks did not stop, see $fewer/log.txt" >&2 ; exit 1
      fi
      grep -aq 'saved on more ranks' "$fewer/log.txt" ||
         { echo "check.sh: leg13's restart on fewer ranks stopped without the rank message, see $fewer/log.txt" >&2 ; exit 1 ; }
      echo "   a restart on $((NP - 1)) ranks stops (saved on $NP)"
   fi
   echo ">> $(basename "$b"): mpirun -np $NP $(basename "$EXE") vortex-2realm.ini (restart)"
   (cd "$b" && mpirun -np "$NP" "$EXE" vortex-2realm.ini > log.txt 2>&1) ||
      { echo "check.sh: leg13's restarted run failed, see $b/log.txt" >&2 ; exit 1 ; }
   # the first topology line of the restarted run is the build from the restored trees
   seam0=$(grep -a -m1 'mixed levels along the seam\|is a 2:1 seam' "$b/log.txt" || true)
   [[ $seam0 =~ "and "[1-9][0-9]*" 2:1 faces" ]] ||
      { echo "check.sh: leg13's restarted run did not start on a mixed seam with 2:1 faces, see $b/log.txt" >&2 ; exit 1 ; }
   grep -aq 'flume: regrid at step 1[6-9][05]: .* \([1-9][0-9]* refined\|[1-9][0-9]* coarsened\)' "$b/log.txt" ||
      { echo "check.sh: leg13's restarted run never changed the grid, see $b/log.txt" >&2 ; exit 1 ; }
   echo "   the restarted run starts on the mixed seam and changes the grid again"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$b" "$a" --ngc 3 --tol 0
   for f in $(cd "$a" && ls ./*.dat); do
      if cmp -s "$a/$f" "$b/$f"; then echo "   ${f#./}: identical"
      else echo "check.sh: leg13's ${f#./} differs between the uninterrupted and the restarted run" >&2 ; exit 1 ; fi
   done
fi
if has_leg leg14; then
   echo "== leg 14: the end_of_step (alpha) cadence against stage_coincident (beta), and alpha's restart (#113)"
   # measured on medooza (#113): 4.5e-2 at CFL 0.4, 1.5e-2 at CFL 0.2, ratio 3.05 (CPU); before #113 3.2 to 6.9
   L14_BOUND=1e-1      # the largest relative alpha-beta difference at CFL 0.4
   L14_RATIO="1.6 4.5" # the range of the CFL 0.4 / CFL 0.2 ratio: the gap must shrink with the step (first order: 2)
   declare -A l14_diff
   for cfl in 0.4 0.2; do
      for cadence in stage_coincident end_of_step; do
         src="$CASE_DIR/work-$TAG-l14-src-$cadence-$cfl" ; rm -rf "$src" ; mkdir -p "$src"
         python3 "$CASE_DIR/make_vortex_forest.py" "$src" --static --cadence $cadence --cfl $cfl
         # a step count, not time_max, ends the runs (about t = 0.3): alpha and beta keep a common last step even if
         # their time steps differ (before #113 alpha took 140 steps to t = 0.3, beta 132)
         sed -i "s/^it_max *=.*/it_max   = $([[ $cfl == 0.4 ]] && echo 120 || echo 240)/" "$src"/vortex-r[12].ini
         run "$CASE_DIR/work-$TAG-l14-$cadence-$cfl" vortex-2realm.ini "$src"/vortex-{2realm,r1,r2}.ini
      done
      out=$("$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$CASE_DIR/work-$TAG-l14-end_of_step-$cfl" \
            "$CASE_DIR/work-$TAG-l14-stage_coincident-$cfl" --ngc 3 --tol 1 --fields-only) ||
         { echo "$out" ; echo "check.sh: leg14's alpha and beta runs at CFL $cfl differ by more than their largest" \
                "values (or have no common step)" >&2 ; exit 1 ; }
      echo "$out" | sed 's/  PASS$//; s/^/   alpha vs beta, CFL '"$cfl"': /'
      l14_diff[$cfl]=$(echo "$out" | sed -n 's/.*max relative difference \([0-9.e+-]*\).*/\1/p' | sort -g | tail -1)
   done
   python3 - "${l14_diff[0.4]}" "${l14_diff[0.2]}" "$L14_BOUND" $L14_RATIO <<'PYEOF'
import sys

d4, d2, bound, lo, hi = (float(x) for x in sys.argv[1:6])
ratio = d4 / d2 if d2 > 0.0 else float("inf")
ok = d4 <= bound and lo <= ratio <= hi
print(f"   largest alpha-beta difference {d4:.3e} at CFL 0.4 (bound {bound:.1e}), {d2:.3e} at CFL 0.2: ratio "
      f"{ratio:.2f} (first order: in [{lo}, {hi}])  {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
PYEOF
   src="$CASE_DIR/work-$TAG-l14-src-restart" ; rm -rf "$src" ; mkdir -p "$src"
   python3 "$CASE_DIR/make_vortex_forest.py" "$src" --static --cadence end_of_step
   sed -i 's/^it_max *=.*/it_max   = 40/' "$src"/vortex-r[12].ini
   a="$CASE_DIR/work-$TAG-l14-restart-A" ; b="$CASE_DIR/work-$TAG-l14-restart-B"
   run "$a" vortex-2realm.ini "$src"/vortex-{2realm,r1,r2}.ini
   sed -i 's/^it_max *=.*/it_max   = 25/; s/^restart_save *=.*/restart_save           = 25/' "$src"/vortex-r[12].ini
   run "$b" vortex-2realm.ini "$src"/vortex-{2realm,r1,r2}.ini
   mv "$b/log.txt" "$b/log-1.txt"
   sed -i -e 's/^restart *=.*/restart                = .true./' -e 's/^it_max *=.*/it_max   = 40/' \
          -e 's/^restart_save *=.*/restart_save           = 0/' "$b"/vortex-r[12].ini
   echo ">> $(basename "$b"): mpirun -np $NP $(basename "$EXE") vortex-2realm.ini (restart)"
   (cd "$b" && mpirun -np "$NP" "$EXE" vortex-2realm.ini > log.txt 2>&1) ||
      { echo "check.sh: leg14's restarted run failed, see $b/log.txt" >&2 ; exit 1 ; }
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$b" "$a" --ngc 3 --tol 0
   for f in $(cd "$a" && ls ./*.dat); do
      if cmp -s "$a/$f" "$b/$f"; then echo "   ${f#./}: identical"
      else echo "check.sh: leg14's ${f#./} differs between the uninterrupted and the restarted run" >&2 ; exit 1 ; fi
   done
fi
if has_leg leg15; then
   echo "== leg 15: refined (2:1) 2-realm Sod on quadtrees vs the single quadtree realm (issue #71)"
   qdir="$CASE_DIR/work-$TAG-l15-src" ; rm -rf "$qdir" ; mkdir -p "$qdir"/{x,y,z}
   qfiles=(sod-amr.ini sod-amr-refined.ini sod-amr-refined-r1.ini sod-amr-refined-r2.ini)
   for f in "${qfiles[@]}"; do
      sed 's/^ratio *= *8/ratio          = 4/' "$CASE_DIR/$f" > "$qdir/x/$f"
      # turned to y: x <-> y in every key, cell count, velocity and seam face
      python3 - "$qdir/x/$f" "$qdir/y/$f" <<'PYEOF'
import re
import sys

s = open(sys.argv[1]).read()
swap = {"_x": "_y", "_y": "_x", "xmin": "ymin", "ymin": "xmin", "xmax": "ymax", "ymax": "xmax"}
s = re.sub(r"_x|_y|xmin|ymin|xmax|ymax", lambda m: swap[m.group(0)], s)
s = re.sub(r"^(ni|nj)(\s*=)", lambda m: ("nj" if m.group(1) == "ni" else "ni") + m.group(2), s, flags=re.M)
s = re.sub(r"^(u|v)(\s*=)", lambda m: ("v" if m.group(1) == "u" else "u") + m.group(2), s, flags=re.M)
s = re.sub(r"^(face_[ab]\s*=\s*[+-])x", r"\1y", s, flags=re.M)
open(sys.argv[2], "w").write(s)
PYEOF
   done
   for f in sod-amr-refined-z{,-r1,-r2}.ini; do
      sed 's/^ratio *= *8/ratio          = 4/' "$CASE_DIR/$f" > "$qdir/z/$f"
   done
   for ax in x y; do
      single="$CASE_DIR/work-$TAG-l15-$ax-single"
      run "$single" sod-amr.ini "$qdir/$ax/sod-amr.ini"
      multi="$CASE_DIR/work-$TAG-l15-$ax-refined"
      run "$multi" sod-amr-refined.ini "$qdir/$ax"/sod-amr-refined{,-r1,-r2}.ini
      "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   done
   work="$CASE_DIR/work-$TAG-l15-z-refused"
   RUN_MUST_FAIL=1 run "$work" sod-amr-refined-z.ini "$qdir/z"/sod-amr-refined-z{,-r1,-r2}.ini
   if ! grep -aq 'a 2:1 seam must be normal to an axis the trees refine' "$work/log.txt"; then
      echo "check.sh: the z split of the quadtrees was not refused, see $work/log.txt" >&2
      exit 1
   fi
   echo "   $(basename "$work"): refused as expected"
   for f in sod-2realm-z{,-r1,-r2}.ini; do
      sed 's/^ratio *= *8/ratio          = 4/; s/^coupling *= *mirror/coupling         = auto/' "$CASE_DIR/$f" > "$qdir/z/$f"
   done
   sed 's/^ratio *= *8/ratio          = 4/' "$VERIF_DIR/sod/sod-z.ini" > "$qdir/z/sod-z.ini"
   single="$CASE_DIR/work-$TAG-l15-z-single"
   run "$single" sod-z.ini "$qdir/z/sod-z.ini"
   multi="$CASE_DIR/work-$TAG-l15-z-auto"
   run "$multi" sod-2realm-z.ini "$qdir/z"/sod-2realm-z{,-r1,-r2}.ini
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
fi
echo "multi-realm verification PASSED (legs: $LEGS) ($TAG)"
