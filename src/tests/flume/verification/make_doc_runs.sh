#!/usr/bin/env bash
# Re-run, with their checkpoints kept, the verification cases whose figures make_doc_figures.py draws from `--runs`
# (issue #47, M3-P6): the check.sh scripts delete the checkpoints after their oracles, so the fields have to be
# produced again. One sub-directory per case in <runs-dir>, np 2, one run at a time:
#   blast   blast-{split,hlld}, blast-wushu-{split,hlld}: PV-1 and PV-2, EGLM, positivity limiter (4 runs, 8 min)
#   vacuum  {leblanc,double-rarefaction,sedov}-{split,hllc}: PV-3, positivity limiter (6 runs, 55 min: Sedov takes 45)
#   rj2a    rj2a, rj2a-hlld: Ryu-Jones 2a, 256 cells, no cleaning (2 runs, 1 min)
#   order   linear-wave-{split,hlld}-{fast,alfven,slow,entropy}-{16,32,64}: MV-5, 1-D (24 runs, 5 min)
#   gravity the M7 gravity cases (#73), as gravity/check.sh writes them: iso2d, poly2d (G0); freefall2d, box2d,
#           box2d-mhd (G1); wave-{weno,riemann}-{8,16,32} (G2a, GV-9); {strat2d,polystrat2d,strat2d-amr}-{unb,bal}
#           saved every 5 steps, lin-{unb,bal}-{6,4} (G2b) (21 runs, 12 min)
# Then: make_doc_figures.py --out docs/public/flume --runs <runs-dir> --only blast near-vacuum mhd-riemann order
#       (gravity: --only gravity-atmospheres gravity-force gravity-order gravity-balance)
#
# Usage: ./make_doc_runs.sh <runs-dir> [blast] [vacuum] [rj2a] [order] [gravity]     (default: all five)
#
# FLUME_EXE overrides the executable (default exe/adam_flume_cpu); the caller owns the matching MPI environment.
set -euo pipefail

VERIF_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$VERIF_DIR/../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NP=2
source "$VERIF_DIR/mhd/numerics.sh"

if [[ $# -lt 1 ]]; then
   echo "make_doc_runs.sh: usage: make_doc_runs.sh <runs-dir> [blast] [vacuum] [rj2a] [order] [gravity]" >&2 ; exit 2
fi
mkdir -p "$1"
RUNS="$(cd "$1" && pwd)" ; shift
GROUPS_="${*:-blast vacuum rj2a order gravity}"
for g in $GROUPS_; do
   if [[ ! $g =~ ^(blast|vacuum|rj2a|order|gravity)$ ]]; then
      echo "make_doc_runs.sh: unknown group '$g' (accepted: blast, vacuum, rj2a, order, gravity)" >&2 ; exit 2
   fi
done
if [[ ! -x "$EXE" ]]; then
   echo "make_doc_runs.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2
fi
FAILED=0

fresh() { # fresh <case>: an empty case directory, its path in W
   W="$RUNS/$1"
   rm -rf "$W" ; mkdir -p "$W"
}
solver() { # solver <ini> <split|llf|hll|hllc|hlld>: select weno-riemann with that solver (characteristic, 6th, weno)
   NUMERICS="" ; NUMERICS_SOLVER=""
   [[ $2 == split ]] && return 0
   NUMERICS=$2 NUMERICS_SOLVER=$2 NUMERICS_RECON=characteristic NUMERICS_CORRECTION=6th NUMERICS_SENSOR=weno
   numerics_apply "$1"
}
limiter_on() { # limiter_on <ini>: add positivity_limiter = cell to the [numerics] block (after reflux)
   sed -i -E 's/^(reflux\s*=.*)$/\1\npositivity_limiter       = cell/' "$1"
   grep -q '^positivity_limiter *= *cell' "$1"
}
run() { # run <ini>: run the case of W
   if (cd "$W" && mpirun -np "$NP" "$EXE" "$1" > log.txt 2>&1); then
      echo "   $(basename "$W"): done"
   else
      echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
   fi
}

for g in $GROUPS_; do
   echo ">> $g"
   case $g in
      blast)
         for scheme in split hlld; do
            fresh "blast-$scheme"
            "$VENV_PY" "$VERIF_DIR/mhd/positivity-probe/make_blast.py" "$W/blast.ini" --eglm --cells 128 > /dev/null
            solver "$W/blast.ini" $scheme ; limiter_on "$W/blast.ini" ; run blast.ini
            fresh "blast-wushu-$scheme"
            "$VENV_PY" "$VERIF_DIR/mhd/positivity-probe/make_blast.py" "$W/blast.ini" --eglm --cells 128 --b-axis x \
               --p-in 1.0e4 --b0 282.0947917738782 --time-max 0.001 --glm-ch 400.0 > /dev/null
            solver "$W/blast.ini" $scheme ; limiter_on "$W/blast.ini" ; run blast.ini
         done ;;
      vacuum)
         for prob in leblanc double-rarefaction sedov; do
            for scheme in split hllc; do
               fresh "$prob-$scheme"
               "$VENV_PY" "$VERIF_DIR/mhd/positivity/make_vacuum.py" "$VERIF_DIR/sod/sod-x.ini" "$W/$prob.ini" \
                  --problem "$prob"
               solver "$W/$prob.ini" $scheme ; limiter_on "$W/$prob.ini" ; run "$prob.ini"
            done
         done ;;
      rj2a)
         for scheme in split hlld; do
            if [[ $scheme == split ]]; then fresh rj2a ; else fresh rj2a-hlld ; fi
            "$VENV_PY" "$VERIF_DIR/mhd/rj2a/make_rj2a.py" "$VERIF_DIR/sod/sod-x.ini" "$W/rj2a-x.ini" --axis x --cells 256 \
               --divergence-control none
            solver "$W/rj2a-x.ini" $scheme ; run rj2a-x.ini
         done ;;
      order)
         for scheme in split hlld; do
            for wave in fast alfven slow entropy; do
               for n in 16 32 64; do
                  fresh "linear-wave-$scheme-$wave-$n"
                  "$VENV_PY" "$VERIF_DIR/mhd/linear-wave/make_linear_wave.py" "$VERIF_DIR/sod/sod-x.ini" \
                     "$W/linear-wave.ini" --wave $wave --geometry 1d --cells $n --amplitude 1.0e-7 \
                     --divergence-control none
                  solver "$W/linear-wave.ini" $scheme ; run linear-wave.ini
               done
            done
         done ;;
      gravity)
         atm() { # atm <case dir> <make_atmosphere arguments>: write input.ini of a gravity case into a fresh W
            fresh "$1" ; shift
            python3 -I "$VERIF_DIR/gravity/make_atmosphere.py" "$W/input.ini" "$@"
         }
         for c in iso2d poly2d freefall2d box2d box2d-mhd; do
            atm "$c" "$c" ; run input.ini
         done
         for scheme in weno riemann; do
            extra=() ; [[ $scheme == riemann ]] && extra=(--riemann hllc)
            for cells in 8 16 32; do
               atm "wave-$scheme-$cells" wave2d --cells "$cells" --time-max 0.2 \
                  --cfl "$(python3 -c "print(0.4 * (8 / $cells) ** (2 / 3))")" "${extra[@]}"
               sed -i 's/^it_save = .*/it_save = 100000/' "$W/input.ini" ; run input.ini
            done
         done
         for c in strat2d polystrat2d strat2d-amr; do
            for mode in unb bal; do
               opt=() ; [[ $mode == bal ]] && opt=(--balanced)
               atm "$c-$mode" "$c" "${opt[@]}"
               sed -i 's/^it_save = .*/it_save = 5/' "$W/input.ini" ; run input.ini
            done
         done
         for mode in unb bal; do
            opt=() ; [[ $mode == bal ]] && opt=(--balanced)
            for a in 6 4; do
               atm "lin-$mode-$a" strat2d --pulse "1e-$a" --time-max 0.1 "${opt[@]}" ; run input.ini
            done
         done ;;
   esac
done
exit $FAILED
