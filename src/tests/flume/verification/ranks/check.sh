#!/usr/bin/env bash
# FLUME verification DV-1: a rank that owns no block (issue #82).
#
# Why: the tree never splits a family of sibling blocks across ranks (tree%mpi_redistribute keeps the ratio siblings
# together), so a grid with fewer families than ranks leaves ranks empty. On FNL every MHD run died on such a rank at
# initialization: the div(B) history copied a blocks_number-sized array to the device, zero bytes there, which the
# OpenACC runtime returns as a null pointer and FUNDAL reports as a failed allocation. Small multi-rank tests hide this
# class of defect, because they rarely leave a rank empty; this leg makes the empty rank the test.
#
# The grid: 4 blocks of 16^2 cells at level 1 on [0,1]^2, one sibling family, so at np N >= 2 rank 0 owns all of them.
# Nine cases (gravity/make_atmosphere.py inputs, the strat2d family with a 10% pressure pulse, 20 steps at CFL 0.4):
#   repro      the issue's input: a uniform EGLM state, x periodic, y walls, no gravity
#   euler      the Euler atmosphere
#   mhd-none, mhd-glm, mhd-eglm   the MHD atmosphere (field across gravity) under each divergence control
#   limiter    EGLM with the positivity limiter
#   balanced   EGLM with [gravity] equilibrium = atmosphere
#   hlld       EGLM on weno-riemann with HLLD
#   amr        EGLM with runtime regrids every 5 steps to level 2 (a Loehner marker on the density, which acts at the
#              runtime regrids only): rank 1 starts empty and the first regrid repopulates it
#   split      the mhd-eglm input with the static balancing on (issue #123, [amr] balance_static, the default): the
#              family is cut by block count, so at np N <= 4 every rank owns blocks and the children of one parent sit on
#              different ranks through the whole run (the other cases set balance_static = .false. to keep their layout)
#   partition  the library unit test tests/amr/test_tree_balance (a host program, np 1, any EXE): the contiguous split of the
#              Morton-ordered units (complete sibling sets kept whole) that tree%mpi_redistribute uses at runtime regrids, on
#              S4c's tree, a uniform tree and 2000 random ones against the exact optimum (issue #123)
#   repart-uniform, repart-refined, repart-forest   the restart repartition tool speedup/repartition.py (issue #133
#              phase 1), not selected by default: a run on --dump-np ranks (default 2) saves a restart at step 10; the
#              restart is continued to step 30 on its own rank count (the reference) and, repartitioned, on every rank
#              count of --to-np (default "1 3 4"), with each split rule of the leg. Cases: uniform = the V2 vortex
#              (verification/vortex/vortex-n064.ini, 16 blocks, static), rules units and siblings; refined = the
#              runtime-regrid vortex (regression/vortex-regrid, levels 1..3, a regrid every 5 steps), restart at step
#              40, continued to 60 through the regrid of step 50 that coarsens a sibling set on the new layout, rule
#              units; forest = the 2-realm forest with realm 2 refined (verification/multirealm/sod-amr-2realm.ini,
#              static, one repartition per realm), rule siblings.
#              Asserted: the dump's layout is the one the leg's first rule writes (repartition.py --info), and the tool
#              repartitioning it to its own rank count rewrites its .fbd and .tnd byte for byte; every repartitioned
#              continuation equals the reference bitwise at its last step, every realm, every field, every interior
#              cell.
#              FNL on WSL (two GPUs): --legs repart-uniform --dump-np 1 --to-np 2.
# Each case runs on 1 and on N ranks. Asserted (ranks_oracle.py):
#   - ownership: in the np N run every rank but 0 owns no block at step 0 (else the leg would exercise nothing), and
#     in the amr case every rank owns blocks at the end; in the split case every rank owns blocks at step 0;
#   - identity: the np N run equals the np 1 run bitwise at the last checkpoint, every field, every cell.
# Measured (P0/P2, FNL on 2 A30 GPUs, nvhpc 25.5): before the fix every MHD case aborted at np 2
# ("dev_assign_to_device: device memory not allocated"), Euler ran; after it every case runs and is bitwise.
#
# Usage: ./check.sh [--np N] [--legs LIST] [--keep] [--it-save N] [--dump-np N] [--to-np "LIST"]
# (--legs: the cases above, default all but the repart-* legs; --cases "repro euler ..." is an older spelling of --legs)
# (--dump-np, --to-np: the rank counts of the repart-* legs, see above)
# (--it-save: a checkpoint every N steps, for the docs figure; the default saves the first and the last only)
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
GEN="$CASE_DIR/../gravity/make_atmosphere.py"
NP=2
source "$VERIF_DIR/legs.sh"
DV1_LEGS="repro euler mhd-none mhd-glm mhd-eglm limiter balanced hlld amr split partition"
legs_declare "$DV1_LEGS repart-uniform repart-refined repart-forest" "$DV1_LEGS"
KEEP=0
IT_SAVE=""
DUMP_NP=2
TO_NP="1 3 4"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --legs)  legs_select "$2" ; shift 2 ;;
      --cases) legs_select "$2" ; shift 2 ;;
      --keep)  KEEP=1 ; shift ;;
      --it-save) IT_SAVE="$2" ; shift 2 ;;
      --dump-np) DUMP_NP="$2" ; shift 2 ;;
      --to-np) TO_NP="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST," \
                    "--cases LIST, --keep, --it-save N, --dump-np N, --to-np LIST)" >&2
               exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ $NP -lt 2 ]]; then
   echo "check.sh: --np must be at least 2 (an empty rank needs a second rank)" >&2 ; exit 2
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

case_ini() { # case_ini <name> <ini>: write the input of a case
   local name=$1 ini=$2 gen
   gen=(python3 -I "$GEN" "$ini")
   case $name in
      repro)    "${gen[@]}" iso2d-mhd --steps 20 --cfl 0.4 ;;
      euler)    "${gen[@]}" strat2d --steps 20 --cfl 0.4 --pulse 0.1 ;;
      mhd-none) "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --divergence-control none ;;
      mhd-glm)  "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --divergence-control glm ;;
      mhd-eglm|limiter|amr|split) "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 ;;
      balanced) "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --balanced ;;
      hlld)     "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --riemann hlld ;;
      *)        echo "check.sh: unknown case '$name'" >&2 ; exit 2 ;;
   esac
   # one sibling family: 4 blocks of 16^2 at level 1
   sed -i -e 's/^ni = .*/ni = 16/' -e 's/^nj = .*/nj = 16/' -e 's/^iu_ref_levels = .*/iu_ref_levels = 1/' \
          -e 's/^max_level = .*/max_level = 1/' -e 's/^ghost_poison = .*/ghost_poison = .false./' "$ini"
   if [[ -n $IT_SAVE ]]; then sed -i "s/^it_save = .*/it_save = $IT_SAVE/" "$ini" ; fi
   # the static grid keeps its sibling-preserving layout (rank 0 owns the family) unless the case is split (issue #123)
   if [[ $name != split ]]; then sed -i 's/^frequency = \(.*\)$/frequency = \1\nbalance_static = .false./' "$ini" ; fi
   case $name in
      repro)
         sed -i -e 's/^type = hydrostatic/type = wall-inviscid/' -e 's/^type = isothermal-atmosphere/type = uniform/' \
                -e '/^\[gravity\]/,/^g_y/d' -e 's/^amr_iterations = 1$/amr_iterations = 1\ns = 0.0/' "$ini" ;;
      limiter)
         sed -i 's/^reflux = .true./reflux = .true.\npositivity_limiter = cell/' "$ini" ;;
      amr)
         sed -i -e 's/^max_level = .*/max_level = 2/' -e 's/^frequency = .*/frequency = 5/' \
                -e 's/^markers_number = .*/markers_number = 1/' "$ini"
         printf '%s\n' "[amr_marker_1]" "mode = 4" "delta_type = max" "delta_fine = 0.0" "delta_coarse = 0.0" \
                       "field = 1" "var = 1" "refine_tol = 0.1" "derefine_tol = 0.02" "buffer = 2" >> "$ini" ;;
   esac
}

TAG="$(basename "$EXE")"

ini_put() { # ini_put <ini> <section> <key> <value>: set a key of a section, added when absent (in place)
   awk -v sec="$2" -v key="$3" -v val="$4" '
      function emit() { if (insec && !done) { print key " = " val ; done = 1 } }
      /^[[:space:]]*\[/ { emit() ; name = $0 ; gsub(/[][[:space:]]/, "", name) ; insec = (name == sec) ; if (insec) seen = 1 }
      insec && $0 ~ "^[[:space:]]*" key "[[:space:]]*=" { if (!done) print key " = " val ; done = 1 ; next }
      { print }
      END { emit() ; if (!seen) { print "[" sec "]" ; print key " = " val } }' "$1" > "$1.tmp" && mv "$1.tmp" "$1"
}

ini_value() { # ini_value <ini> <key>: the first value of a key (comments stripped)
   sed -n "s/^[[:space:]]*$2[[:space:]]*=[[:space:]]*\([^ ;]*\).*/\1/p" "$1" | head -1
}

REPART="$REPO_ROOT/src/tests/flume/speedup/repartition.py"

repart_inputs() { # repart_inputs <case> <dir> <it_max> <restart_save> <restart>: the inputs of one run of a repart leg
   local case=$1 dir=$2 ini
   mkdir -p "$dir"
   case $case in
      uniform) cp "$VERIF_DIR/vortex/vortex-n064.ini" "$dir/input.ini" ;;
      refined) cp "$VERIF_DIR/../regression/vortex-regrid/input.ini" "$dir/input.ini" ;;
      forest)  cp "$VERIF_DIR/multirealm/sod-amr-2realm.ini" "$dir/input.ini"
               cp "$VERIF_DIR/multirealm/sod-amr-2realm-r1.ini" "$VERIF_DIR/multirealm/sod-amr-2realm-r2.ini" "$dir/" ;;
   esac
   for ini in $(repart_realms "$case"); do
      ini_put "$dir/$ini" time it_max "$3"
      ini_put "$dir/$ini" time time_max 100.0
      ini_put "$dir/$ini" IO it_save 100000
      ini_put "$dir/$ini" IO restart_save "$4"
      ini_put "$dir/$ini" IO restart "$5"
   done
}

repart_realms() { # repart_realms <case>: the realm inputs of a case
   if [[ $1 == forest ]]; then echo "sod-amr-2realm-r1.ini sod-amr-2realm-r2.ini" ; else echo "input.ini" ; fi
}

repart_run() { # repart_run <dir> <np>: run the case of a directory
   echo ">> $(basename "$1"): mpirun -np $2 $TAG"
   if ! (cd "$1" && mpirun -np "$2" "$EXE" input.ini > log.txt 2>&1); then
      echo "   $(basename "$1"): the run failed, see $1/log.txt  FAIL" ; return 1
   fi
}

repart_leg() { # repart_leg <case>: the restart repartition leg of a case (issue #133 phase 1)
   local case=$1 base="$CASE_DIR/work-$TAG-repart-$1" rules ini rb basenames=() ok=0 m rule dir f n0=10 n1=30
   case $case in uniform) rules="units siblings" ;; refined) rules="units" ;; forest) rules="siblings" ;; esac
   # refined: the window holds the regrid of step 50, the first to coarsen (a sibling set merged on the new layout)
   [[ $case == refined ]] && { n0=40 ; n1=60 ; }
   rm -rf "$base" ; mkdir -p "$base"
   # the dump: n0 steps on DUMP_NP ranks, a restart saved at step n0
   repart_inputs "$case" "$base/dump" $n0 $n0 .false.
   repart_run "$base/dump" "$DUMP_NP" || return 1
   for ini in $(repart_realms "$case"); do
      basenames+=("$base/dump/$(ini_value "$base/dump/$ini" restart_basename)")
   done
   # the layout: the leg's first rule writes the dump's layout, and rewrites its files byte for byte
   if ! python3 -I "$REPART" "${basenames[@]}" --info | tee "$base/info.txt" | sed 's/^/   /'; then
      echo "   repart-$case: repartition.py --info failed  FAIL" ; return 1
   fi
   rule="${rules%% *}"
   if [[ $(grep -c "layout reproduced by rule: .*$rule" "$base/info.txt") -ne ${#basenames[@]} ]]; then
      echo "   repart-$case: the dump's layout is not the one rule $rule writes  FAIL" ; ok=1
   fi
   python3 -I "$REPART" "${basenames[@]}" --to "$DUMP_NP" --out "$base/same-np" --rule "$rule" > /dev/null || ok=1
   for rb in "${basenames[@]}"; do
      for f in "$rb".tnd "$rb"-proc*.fbd "$rb"-proc*.wmf; do
         [[ -f $f ]] || continue
         if ! cmp -s "$f" "$base/same-np/$(basename "$f")"; then
            echo "   repart-$case: $(basename "$f") rewritten for np $DUMP_NP by rule $rule differs  FAIL" ; ok=1
         fi
      done
   done
   [[ $ok -eq 0 ]] && echo "   repart-$case: np $DUMP_NP -> np $DUMP_NP by rule $rule rewrites the dump byte for byte  PASS"
   # the reference: the restart continued to step n1 on its own rank count
   repart_inputs "$case" "$base/ref" $n1 0 .true.
   for rb in "${basenames[@]}"; do
      cp "$rb".* "$rb"-proc*.fbd "$base/ref/"
      cp "$rb"-proc*.wmf "$base/ref/" 2> /dev/null || true
   done
   repart_run "$base/ref" "$DUMP_NP" || return 1
   # the repartitioned continuations, each compared with the reference realm by realm
   for rule in $rules; do
      for m in $TO_NP; do
         dir="$base/np$m-$rule"
         repart_inputs "$case" "$dir" $n1 0 .true.
         if ! python3 -I "$REPART" "${basenames[@]}" --to "$m" --out "$dir" --rule "$rule" | sed 's/^/   /'; then
            echo "   repart-$case: repartition.py to np $m ($rule) failed  FAIL" ; ok=1 ; continue
         fi
         repart_run "$dir" "$m" || { ok=1 ; continue ; }
         for ini in $(repart_realms "$case"); do
            "$VENV_PY" "$CASE_DIR/ranks_oracle.py" same "$base/ref" "$dir" \
               --basename "$(ini_value "$dir/$ini" output_basename)" --ini "$ini" || ok=1
         done
      done
   done
   if [[ $KEEP -eq 0 && $ok -eq 0 ]]; then rm -rf "$base" ; fi
   return $ok
}

FAILED=0
for name in $LEGS; do
   if [[ $name == repart-* ]]; then
      repart_leg "${name#repart-}" || FAILED=1
      continue
   fi
   if [[ $name == partition ]]; then
      unit="$REPO_ROOT/exe/test_tree_balance"
      echo ">> partition: mpirun -np 1 test_tree_balance"
      if [[ -x $unit ]] && mpirun -np 1 "$unit" > "$CASE_DIR/partition.log" 2>&1 && grep -q "test_tree_balance: PASSED" "$CASE_DIR/partition.log"; then
         echo "   partition: PASS" ; rm -f "$CASE_DIR/partition.log"
      else
         echo "   partition: FAIL (build it: fobis build --mode test-tree-balance-gnu; see $CASE_DIR/partition.log)" ; FAILED=1
      fi
      continue
   fi
   for np in 1 "$NP"; do
      work="$CASE_DIR/work-$TAG-$name-np$np"
      rm -rf "$work" ; mkdir -p "$work"
      case_ini "$name" "$work/input.ini"
      echo ">> $name: mpirun -np $np $TAG"
      if ! (cd "$work" && mpirun -np "$np" "$EXE" input.ini > log.txt 2>&1); then
         echo "   $name np $np: the run failed, see $work/log.txt  FAIL"
         FAILED=1
         continue 2
      fi
   done
   one="$CASE_DIR/work-$TAG-$name-np1" ; many="$CASE_DIR/work-$TAG-$name-np$NP"
   own=(--empty-at-start) ; [[ $name == amr ]] && own+=(--filled-at-end)
   [[ $name == split ]] && own=(--split-at-start)
   "$VENV_PY" "$CASE_DIR/ranks_oracle.py" ownership "$many" "${own[@]}" || FAILED=1
   "$VENV_PY" "$CASE_DIR/ranks_oracle.py" same "$one" "$many" || FAILED=1
   if [[ $KEEP -eq 0 ]]; then
      rm -rf "$one" "$many"
   fi
done
if [[ $FAILED -ne 0 ]]; then
   echo "DV-1 FAILED ($TAG, np $NP)" ; exit 1
fi
echo "DV-1 PASSED ($TAG, np $NP, legs: $LEGS)"
