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
# Usage: ./check.sh [--build] [--np N] [--legs LIST]   LIST: leg1 ... leg7 or all (default all); leg4 compares against
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
legs_declare "leg1 leg2 leg3 leg4 leg5 leg6 leg7"

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
echo "multi-realm verification PASSED (legs: $LEGS) ($TAG)"
