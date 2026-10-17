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
# Each case runs on 1 and on N ranks. Asserted (ranks_oracle.py):
#   - ownership: in the np N run every rank but 0 owns no block at step 0 (else the leg would exercise nothing), and
#     in the amr case every rank owns blocks at the end;
#   - identity: the np N run equals the np 1 run bitwise at the last checkpoint, every field, every cell.
# Measured (P0/P2, FNL on 2 A30 GPUs, nvhpc 25.5): before the fix every MHD case aborted at np 2
# ("dev_assign_to_device: device memory not allocated"), Euler ran; after it every case runs and is bitwise.
#
# Usage: ./check.sh [--np N] [--legs LIST] [--keep] [--it-save N]
# (--legs: the cases above, default all; --cases "repro euler ..." is an older spelling of --legs)
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
legs_declare "repro euler mhd-none mhd-glm mhd-eglm limiter balanced hlld amr"
KEEP=0
IT_SAVE=""

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --legs)  legs_select "$2" ; shift 2 ;;
      --cases) legs_select "$2" ; shift 2 ;;
      --keep)  KEEP=1 ; shift ;;
      --it-save) IT_SAVE="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST," \
                    "--cases LIST, --keep, --it-save N)" >&2
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
      mhd-eglm|limiter|amr) "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 ;;
      balanced) "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --balanced ;;
      hlld)     "${gen[@]}" strat2d-mhdx --steps 20 --cfl 0.4 --pulse 0.1 --riemann hlld ;;
      *)        echo "check.sh: unknown case '$name'" >&2 ; exit 2 ;;
   esac
   # one sibling family: 4 blocks of 16^2 at level 1
   sed -i -e 's/^ni = .*/ni = 16/' -e 's/^nj = .*/nj = 16/' -e 's/^iu_ref_levels = .*/iu_ref_levels = 1/' \
          -e 's/^max_level = .*/max_level = 1/' -e 's/^ghost_poison = .*/ghost_poison = .false./' "$ini"
   if [[ -n $IT_SAVE ]]; then sed -i "s/^it_save = .*/it_save = $IT_SAVE/" "$ini" ; fi
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
FAILED=0
for name in $LEGS; do
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
