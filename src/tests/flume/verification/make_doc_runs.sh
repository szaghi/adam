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
#           saved every 5 steps, lin-{unb,bal}-{6,4} (G2b); rt-{8,16} (G3, GV-4) (23 runs, 27 min)
#   rt      rt2d-nl, rt2d-mhd: the nonlinear single-mode Rayleigh-Taylor, hydro and EGLM, N = 128 to t = 5 (2 runs)
#   rt-sym  rt2d-mhd under other cleanings (#83), rt-sym-{none,glm,ch1.5,ch6} and the EGLM damping rt-sym-{a0,a1,a5,
#           mincell}; the reference (EGLM, glm_ch 3, glm_alpha 0.18, length 1) is the rt group's rt2d-mhd (8 runs,
#           about 2 min each with the FNL build on 2 A30 GPUs)
#   ranks   ranks-amr-np{1,2}: the DV-1 runtime-AMR case of ranks/check.sh (#82), a checkpoint every step: rank 1 owns
#           no block until the first regrid (2 runs, under a minute)
#   selfgravity sv1-box2d-{uniform,linear}: the closed box (M7 box2d, 16 cells per block side) to t = 0.2, about 20
#           checkpoints; sv1-wave-{uniform,linear}-{8,16,32}: the falling wave (GV-9) to t = 0.2; each through the
#           constant-g path and through the potential path (#81, SV-1) (8 runs, 2 min); sv2-poisson: the SV-2 unit
#           test of the FFT Poisson solver (POISSON_EXE, default exe/test_poisson_fft, mode test-poisson-fft-gnu) with
#           its log and the potential of every level dumped (1 run, seconds)
#   jeans   the Jeans cases of selfgravity/jeans.sh (#81, SV-3, SV-4) with the FFT executable (JEANS_EXE, default
#           exe/adam_flume_cpu_fft, mode flume-cpu-gnu-fft): sv3-2d-{8,16,32}, sv3-3d-{8,16}, sv4-linear (a checkpoint
#           every 5 steps), sv4-collapse-{0.4,0.2,0.1} (8 runs, about 15 min)
#   isolated the isolated self-gravity cases (#81, P4): sv6-poisson, the SV-6 unit test (`--isolated`, POISSON_EXE) with
#           its log and dumps; sv7-l{2,3}, the equilibrium star of selfgravity/star.sh at 32^3 and 64^3 (a light
#           atmosphere, the positivity limiter) to t = 2; sv8-collapse-{0.4,0.2,0.1}, the collapsing star, a checkpoint
#           every 10 steps (JEANS_EXE, default exe/adam_flume_cpu_fft) (6 runs, about 60 min: sv7-l3 takes 45)
# Then: make_doc_figures.py --out docs/public/flume --runs <runs-dir> --only blast near-vacuum mhd-riemann order
#       (gravity: --only gravity-atmospheres gravity-force gravity-order gravity-balance gravity-rt-linear;
#        rt: --only gravity-rt; rt and rt-sym: --only gravity-rt-symmetry; ranks: --only ranks;
#        selfgravity: --only selfgravity-potential selfgravity-poisson; jeans: --only selfgravity-jeans;
#        isolated: --only selfgravity-isolated selfgravity-star)
#
# Usage: ./make_doc_runs.sh <runs-dir> [blast] [vacuum] [rj2a] [order] [gravity] [rt] [rt-sym] [ranks] [selfgravity]
#        [jeans] [isolated] (default: all eleven)
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
   echo "make_doc_runs.sh: usage: make_doc_runs.sh <runs-dir> [blast] [vacuum] [rj2a] [order] [gravity] [rt]" \
        "[rt-sym] [ranks] [selfgravity] [jeans] [isolated]" >&2
   exit 2
fi
mkdir -p "$1"
RUNS="$(cd "$1" && pwd)" ; shift
GROUPS_="${*:-blast vacuum rj2a order gravity rt rt-sym ranks selfgravity jeans isolated}"
for g in $GROUPS_; do
   if [[ ! $g =~ ^(blast|vacuum|rj2a|order|gravity|rt|rt-sym|ranks|selfgravity|jeans|isolated)$ ]]; then
      echo "make_doc_runs.sh: unknown group '$g' (accepted: blast, vacuum, rj2a, order, gravity, rt, rt-sym, ranks," \
           "selfgravity, jeans, isolated)" >&2
      exit 2
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
         done
         for cells in 8 16; do
            atm "rt-$cells" rt2d --cells "$cells" --it-save $((5 * cells)) ; run input.ini
         done ;;
      rt)
         for c in rt2d-nl rt2d-mhd; do
            fresh "$c"
            python3 -I "$VERIF_DIR/gravity/make_atmosphere.py" "$W/input.ini" "$c" --cells 32 --it-save 100
            run input.ini
         done ;;
      rt-sym)
         for spec in "none:--divergence-control none" "glm:--divergence-control glm" "ch1.5:--glm-ch 1.5" \
                     "ch6:--glm-ch 6.0" "a0:--glm-alpha 0.0" "a1:--glm-alpha 1.0" "a5:--glm-alpha 5.0" \
                     "mincell:--glm-damping-length min-cell"; do
            fresh "rt-sym-${spec%%:*}"
            read -ra opt <<< "${spec#*:}"
            python3 -I "$VERIF_DIR/gravity/make_atmosphere.py" "$W/input.ini" rt2d-mhd --cells 32 --it-save 100 \
               "${opt[@]}"
            run input.ini
         done ;;
      ranks)
         if FLUME_EXE="$EXE" "$VERIF_DIR/ranks/check.sh" --cases amr --np "$NP" --keep --it-save 1; then
            for np in 1 "$NP"; do
               fresh "ranks-amr-np$np" ; rmdir "$W"
               mv "$VERIF_DIR/ranks/work-$(basename "$EXE")-amr-np$np" "$W"
            done
         else
            echo "make_doc_runs.sh: ranks/check.sh (amr) failed" >&2 ; FAILED=1
         fi ;;
      selfgravity)
         for pot in uniform linear; do
            fresh "sv1-box2d-$pot"
            python3 -I "$VERIF_DIR/gravity/make_atmosphere.py" "$W/input.ini" box2d --cells 16 --time-max 0.2 \
               --potential "$pot"
            sed -i 's/^it_save = .*/it_save = 2/' "$W/input.ini"
            run input.ini
            for cells in 8 16 32; do
               fresh "sv1-wave-$pot-$cells"
               python3 -I "$VERIF_DIR/gravity/make_atmosphere.py" "$W/input.ini" wave2d --cells "$cells" \
                  --time-max 0.2 --cfl "$(python3 -c "print(0.4 * (8 / $cells) ** (2 / 3))")" --potential "$pot"
               sed -i 's/^it_save = .*/it_save = 100000/' "$W/input.ini"
               run input.ini
            done
         done
         fresh sv2-poisson
         POISSON_EXE="${POISSON_EXE:-$REPO_ROOT/exe/test_poisson_fft}"
         if (cd "$W" && mpirun -np "$NP" "$POISSON_EXE" --dump "$W" > log.txt 2>&1); then
            echo "   $(basename "$W"): done"
         else
            echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
         fi ;;
      jeans)
         JEANS_EXE="${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
         jeans() { # jeans <case dir> <make_jeans arguments>: write and run one Jeans case with the FFT executable
            fresh "$1" ; shift
            python3 -I "$VERIF_DIR/selfgravity/make_jeans.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$NP" "$JEANS_EXE" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         for cells in 8 16 32; do jeans "sv3-2d-$cells" jeans2d --cells "$cells" ; done
         for cells in 8 16; do jeans "sv3-3d-$cells" jeans3d --cells "$cells" ; done
         jeans sv4-linear jeansu2d --cells 16 --it-save 5
         for cfl in 0.4 0.2 0.1; do jeans "sv4-collapse-$cfl" jeansc2d --cells 16 --cfl "$cfl" ; done ;;
      isolated)
         fresh sv6-poisson
         POISSON_EXE="${POISSON_EXE:-$REPO_ROOT/exe/test_poisson_fft}"
         if (cd "$W" && mpirun -np "$NP" "$POISSON_EXE" --isolated --dump "$W" > log.txt 2>&1); then
            echo "   $(basename "$W"): done"
         else
            echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
         fi
         STAR_EXE="${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
         star() { # star <case dir> <make_polytrope arguments>: write and run one polytrope case with the FFT executable
            fresh "$1" ; shift
            python3 -I "$VERIF_DIR/selfgravity/make_polytrope.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$NP" "$STAR_EXE" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         for lev in 2 3; do
            star "sv7-l$lev" star --levels "$lev" --atmosphere 1e-4 --limiter --time-max 2.0 --it-save 10
         done
         for cfl in 0.4 0.2 0.1; do star "sv8-collapse-$cfl" starc --cfl "$cfl" --it-save 10 ; done ;;
   esac
done
exit $FAILED
