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
#           every 10 steps (JEANS_EXE, default exe/adam_flume_cpu_fft); sv14-{uniform,box},
#           sv14-collapse-{0.4,0.2,0.1}, the smaller star (radius 0.3) on AMR (#87, P3: isolated self-gravity by the
#           multigrid's James step) as selfgravity/star.sh runs them (11 runs)
#   amr     sv9-poisson: the SV-9 unit test of the multigrid Poisson solver on AMR grids (#81, P5a; POISSON_MG_EXE,
#           default exe/test_poisson_mg, mode test-poisson-mg-gnu) with its log; sv13-poisson, the same with
#           `--isolated` (#87, P2: the James step, SV-13); the SV-10 Jeans cases on the
#           box-refined grid (#81, P5b; JEANS_EXE): sv10-2d-{8,16,32}, the stable wave, sv10-uniform-{8,16,32}, the
#           same on the uniform grid, sv10-collapse-{0.4,0.2,0.1}, the collapse, sv10-stage-{0.4,0.2,0.1}, the same
#           with reflux_cadence = stage (P5c), and sv11-regrid, the collapse (base 32^2) under runtime AMR with a
#           checkpoint every 10 steps (14 runs, about 35 min)
#   mhd     the MHD figures of the gallery (issue #85, R2a): brio-wu (512 cells, GLM), rj4d (512, no cleaning), glm-pulse
#           (128, GLM, amplitude 0.1), orszag-tang (128^2, t = 0.5), orszag-tang-amr (32^2 + [0.25, 0.75]^2 refined),
#           rotor (256^2), field-loop (128 x 64), rj2a-hlld-prim (RJ2a on HLLD with primitive interpolation) (8 runs,
#           about 40 min)
#   rotated rst-{63,45}-256 and their 1-D references rst-{63,45}-ref: the rotated shock tube MV-8 (4 runs, about 40 min)
#   quadtree ot-octree, ot-quad-nk1: MV-15, Orszag-Tang 32^2 + [0.25, 0.75]^2 refined, t = 0.2 (2 runs, about 3 min)
#   regrid  regrid-vortex-{conservative,linear}: the rg3 vortex regridding every 5 steps, a checkpoint every 5 steps
#           (2 runs, about 2 min)
#   collapse-regrid sv15-{conservative-state,conservative}: the runtime-AMR collapse of selfgravity/star.sh SV-15
#           (#88); sv16-{jeans,control}: SV-16, the smaller star with and without the Jeans marker (#91); a
#           checkpoint every 2 steps, with the FFT executable (JEANS_EXE); the conservative and control runs go
#           non-finite by design and are kept; sv17-de: SV-17, the SV-15 star with euler-de (#91 P2) (5 runs, 6 min)
#   dualenergy noh-{euler,euler-de}: the cold planar Noh shock of dualenergy/check.sh de2 (#91 P2) (2 runs, 1 min)
#   scalars ps-sod (Sod with two scalars, scalars/check.sh ps1) and ps-vortex (the regridded vortex with a uniform
#           scalar, ps4, a checkpoint every 5 steps) (#92 P3) (2 runs, 2 min)
#   cooling co-box (the uniform box of cooling/check.sh co1 from T = 5), co-limited (co1 with dt_factor 0.01) and
#           co-shock (the radiative shock of co7) (#93 P3) (3 runs, 2 min)
#   eos     eo-tube (the isothermal shock tube of eos/check.sh eo1), eo-tube-hllc (eo6, weno-riemann hllc),
#           eo-collide (the Mach-5 collision of eo2), eo-wave-{12,24,48} (the acoustic wave of eo3) and eo-steep (the
#           wave of amplitude 0.3 of eo5 on 200 cells, a checkpoint every 20 steps) (#100 P1) (7 runs, 2 min); and
#           the barotropic table POLY (P3): bo-tube (eb1, 1600 cells), bo-collide-hll and bo-collide-hll-slow (eb7 hll in
#           the shifted frame and in the slow-shock frame), bo-wave-{12,24,48} (eb3), bo-mf2021-{200,400,800,1600} (eb6)
#           (10 more runs, 4 min)
#   sinks   si-binary (the circular binary of sinks/check.sh si2, one orbit), si-shu-{32,64} (the overdense singular
#           isothermal sphere of si11 at 32^3 and 64^3) and si-collapse-amr (the si10 collapse on runtime AMR, two levels
#           above the base) (#94 P4), with the FFT executable (JEANS_EXE) (4 runs, 6 min)
#   eosgrav the barotropic gases under gravity (#100 P4, eos/gravity.sh) with the FFT executable (JEANS_EXE):
#           eg-jeans-{isothermal,barotropic,box}-{8,16,32} (eg2, eg4), eg-jeansu-{isothermal,barotropic} (eg3, a
#           checkpoint every 5 steps), eg-be-{3,6,9} (eg7, every 10 steps), eg-shu (eg6, 32^3) and eg-bb, the
#           Boss-Bodenheimer cloud on a 32^3 base with two Jeans-refined levels and sinks to 1.25 t_ff (8 ranks, a
#           checkpoint every 20 steps)
#   multirealm sod-x, sod-2realm, sod-amr, sod-amr-refined: the Sod unions of multirealm/check.sh legs 1 and 4 and
#           their single-realm twins (4 runs, about 4 min)
#   seams   the Sod shock crossing a static 2:1 face at x = 0.75 (amr-shock/make_amr_shock.py, #68; issue #85 R2b):
#           as-{c2f,f2c} (AS, today's operators), as-uniform-{200,400}; e0-{fill}-{restriction}-{on,off} for the fills
#           tricubic and restriction-compatible, the restrictions mean and lagrange6, the reflux on and off, and
#           e0-{tricubic,compatible}-400 (E0); d-{quintic,weno6,quintic-mp}-{c2f,f2c} on lagrange6 (D1, D2);
#           e2-troubled-{c2f,f2c} (lagrange6 + quintic under reflux_policy = troubled, E2); N = 200 unless named
# The groups below run on FNL (issue #85, R3; FLUME_EXE=exe/adam_flume_fnl, the caller's nvhpc environment):
#   blast77 blast77-{split,hlld}: the av3b/av3c MHD blast (EGLM, positivity limiter), uniform N = 128, to t = 0.01,
#           a checkpoint every 20 steps, splitting and weno-riemann HLLD (#77) (2 runs)
#   divb    ad-cpaw-n64: the AD CPAW refined at N = 64, its last checkpoint kept; g2-cpaw-{coarse,mincell}-{64,128}: the
#           same with glm_damping_length = the coarse spacing sqrt(5)/N and min-cell (G2); loop-{glm,eglm}-{n64,
#           n64-noreflux,uniform-n128}: the A2 field loop, the divergence history every step and a checkpoint every
#           20 steps (G3a, G3b) (11 runs)
#   rj2a-shock as-rj2a-{c2f,f2c,uniform}: the AS RJ2a crossing at N = 256, t = 0.2 (3 runs)
# Then: make_doc_figures.py --out docs/public/flume --runs <runs-dir> --only blast near-vacuum mhd-riemann order
#       (gravity: --only gravity-atmospheres gravity-force gravity-order gravity-balance gravity-rt-linear;
#        rt: --only gravity-rt; rt and rt-sym: --only gravity-rt-symmetry; ranks: --only ranks;
#        selfgravity: --only selfgravity-potential selfgravity-poisson; jeans: --only selfgravity-jeans;
#        isolated: --only selfgravity-isolated selfgravity-star selfgravity-star-amr;
#        collapse-regrid: --only selfgravity-star-regrid; dualenergy: --only dual-energy; scalars: --only passive-scalars;
#        cooling: --only radiative-cooling; sinks: --only sink-particles; eos: --only isothermal-eos barotropic-eos;
#        amr: --only selfgravity-amr selfgravity-amr-jeans
#        selfgravity-amr-regrid selfgravity-amr-isolated)
#
# Usage: ./make_doc_runs.sh <runs-dir> [blast] [vacuum] [rj2a] [order] [gravity] [rt] [rt-sym] [ranks] [selfgravity]
#        [jeans] [isolated] [amr] (default: all twelve)
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
        "[rt-sym] [ranks] [selfgravity] [jeans] [isolated] [amr] [mhd] [rotated] [quadtree] [regrid] [multirealm] [seams] [blast77] [divb] [rj2a-shock] [collapse-regrid] [dualenergy] [scalars] [cooling] [sinks] [eos]" >&2
   exit 2
fi
mkdir -p "$1"
RUNS="$(cd "$1" && pwd)" ; shift
GROUPS_="${*:-blast vacuum rj2a order gravity rt rt-sym ranks selfgravity jeans isolated amr mhd rotated quadtree regrid multirealm seams blast77 divb rj2a-shock collapse-regrid dualenergy scalars cooling sinks eos eosgrav}"
for g in $GROUPS_; do
   if [[ ! $g =~ ^(blast|vacuum|rj2a|order|gravity|rt|rt-sym|ranks|selfgravity|jeans|isolated|amr|mhd|rotated|quadtree|regrid|multirealm|seams|blast77|divb|rj2a-shock|collapse-regrid|dualenergy|scalars|cooling|sinks|eos|eosgrav)$ ]]; then
      echo "make_doc_runs.sh: unknown group '$g' (accepted: blast, vacuum, rj2a, order, gravity, rt, rt-sym, ranks," \
           "selfgravity, jeans, isolated, amr, mhd, rotated, quadtree, regrid, multirealm, seams, blast77, divb," \
           "rj2a-shock, collapse-regrid, dualenergy, scalars, cooling, sinks, eos, eosgrav)" >&2
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
         for cfl in 0.4 0.2 0.1; do star "sv8-collapse-$cfl" starc --cfl "$cfl" --it-save 10 ; done
         # SV-14 (issue #87, P3), the smaller star on AMR as selfgravity/star.sh runs it
         for grid in uniform box; do
            star "sv14-$grid" star --levels 2 --radius 0.3 $([[ $grid == box ]] && echo --box) --atmosphere 1e-4 \
               --limiter --time-max 2.0 --it-save 10
         done
         for cfl in 0.4 0.2 0.1; do
            star "sv14-collapse-$cfl" starc --radius 0.3 --box --reflux-cadence stage --cfl "$cfl"
         done ;;
      collapse-regrid)
         # SV-15 (issue #88): the runtime-AMR collapse as selfgravity/star.sh runs it; conservative stops by design
         for prol in conservative-state conservative; do
            fresh "sv15-$prol"
            python3 -I "$VERIF_DIR/selfgravity/make_polytrope.py" "$W/input.ini" starc --radius 0.3 --reflux-cadence stage \
               --regrid 5 --it-save 2 --prolongation "$prol"
            if (cd "$W" && mpirun -np "$NP" "${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            elif [[ $prol == conservative ]]; then
               echo "   $(basename "$W"): stopped, as expected (#88)"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         done
         # SV-16 (issue #91): the radius-0.22 star with and without the Jeans marker; the control stops by design
         for v in jeans control; do
            fresh "sv16-$v"
            python3 -I "$VERIF_DIR/selfgravity/make_polytrope.py" "$W/input.ini" starc --radius 0.22 --reflux-cadence stage \
               --regrid 5 --it-save 2 $([[ $v == jeans ]] && echo --jeans 8)
            if (cd "$W" && mpirun -np "$NP" "${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            elif [[ $v == control ]]; then
               echo "   $(basename "$W"): stopped, as expected (#91)"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         done
         # SV-17 (issue #91, P2): the SV-15 star with dual energy
         fresh sv17-de
         python3 -I "$VERIF_DIR/selfgravity/make_polytrope.py" "$W/input.ini" starc --radius 0.3 --reflux-cadence stage \
            --regrid 5 --it-save 2 --physical-model euler-de
         if (cd "$W" && mpirun -np "$NP" "${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}" input.ini > log.txt 2>&1); then
            echo "   $(basename "$W"): done"
         else
            echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
         fi ;;
      scalars)
         fresh ps-sod
         python3 -I "$VERIF_DIR/scalars/make_ps.py" "$W/input.ini" --scalars 2
         run input.ini
         fresh ps-vortex
         python3 -I "$VERIF_DIR/scalars/make_ps.py" "$W/input.ini" --case vortex --scalars 1
         sed -i -E 's/^(\s*it_save\s*=).*/\1 5/' "$W/input.ini"
         run input.ini ;;
      cooling)
         fresh co-box
         python3 -I "$VERIF_DIR/cooling/make_co.py" "$W/input.ini" --p 5.0 --time-max 6.5
         run input.ini
         fresh co-limited
         python3 -I "$VERIF_DIR/cooling/make_co.py" "$W/input.ini" --time-max 2.0 --dt-factor 0.01
         run input.ini
         fresh co-shock
         python3 -I "$VERIF_DIR/cooling/make_co.py" "$W/input.ini" --case shock
         run input.ini ;;
      eos)
         iso() { # iso <case dir> <ranks> <make_eos arguments>: write and run one isothermal case (one block: 1 rank)
            local np=$2
            fresh "$1" ; shift 2
            python3 -I "$VERIF_DIR/eos/make_eos.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$np" "$EXE" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         iso eo-tube 1 --case tube
         iso eo-tube-hllc 1 --case tube --scheme weno-riemann --solver hllc
         iso eo-collide 1 --case collide --mach 5.0 --time-max 1.0
         for n in 12 24 48; do iso "eo-wave-$n" 2 --case wave --ni "$n" --cfl 0.1 --rk runge-kutta-ssp-54 ; done
         iso eo-steep 2 --case wave --ni 200 --set "initial_conditions:wave_amplitude=0.3" --set "IO:it_save=20"
         poly=(--eos barotropic --densities "2.0 8.0" --exponents "1.0 1.4 1.6666666666666667")
         iso bo-tube 1 --case tube "${poly[@]}" --left "10.0 0.0 0.0" --right "0.5 0.0 0.7" --time-max 0.15 --ni 1600
         iso bo-collide-hll 1 --case collide "${poly[@]}" --left "1.0 3.4" --right "1.5 -2.6" --time-max 0.3 \
            --scheme weno-riemann --solver hll
         iso bo-collide-hll-slow 1 --case collide "${poly[@]}" --left "1.0 4.0" --right "1.5 -2.0" --time-max 0.3 \
            --scheme weno-riemann --solver hll
         period=$("$VENV_PY" -I -c 'print(repr(1.0 / (1.4 * 2.0**-0.4 * 5.0**0.4) ** 0.5))')
         for n in 12 24 48; do
            iso "bo-wave-$n" 2 --case wave "${poly[@]}" --rho0 5.0 --ni "$n" --cfl 0.1 --rk runge-kutta-ssp-54 \
               --time-max "$period"
         done
         for n in 200 400 800 1600; do
            iso "bo-mf2021-$n" 2 --case tube --eos barotropic --c0 2.0e4 --left "1.0e-7 0.0 0.0" --right "1.0e-8 0.0 0.0" \
               --time-max 3.0e-7 --ni "$n" --level 1 --null-y --set bc_x_min:type=periodic \
               --set bc_x_max:type=periodic --set weno:weights=si
         done ;;
      eosgrav)
         grav() { # grav <case dir> <ranks> <make_eos_gravity arguments>: write and run one case with the FFT executable
            local np=$2
            fresh "$1" ; shift 2
            "$VENV_PY" -I "$VERIF_DIR/eos/make_eos_gravity.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$np" "${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}" input.ini > log.txt 2>&1)
            then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         for n in 8 16 32; do
            for eos in isothermal barotropic; do grav "eg-jeans-$eos-$n" 4 jeans --eos "$eos" --cells "$n" ; done
            grav "eg-jeans-box-$n" 4 jeans --eos isothermal --box --cells "$n"
         done
         for eos in isothermal barotropic; do grav "eg-jeansu-$eos" 4 jeansu --eos "$eos" --cells 16 --it-save 5 ; done
         for xi in 3 6 9; do grav "eg-be-$xi" 4 be --xi "$xi" --time-max 4.0 --it-save 10 ; done
         grav eg-shu 4 shu
         grav eg-bb 8 bb --cells 32 --amr 2 --it-save 20 ;;
      sinks)
         sink() { # sink <case dir> <make_si arguments>: write and run one sink case with the FFT executable
            fresh "$1" ; shift
            python3 -I "$VERIF_DIR/sinks/make_si.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$NP" "${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}" input.ini > log.txt 2>&1)
            then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         sink si-binary --case binary --time-max 1.5708
         for cells in 32 64; do sink "si-shu-$cells" --case shu --cells "$cells" ; done
         sink si-collapse-amr --case collapse --cells 16 --amr 2 ;;
      dualenergy)
         for model in euler euler-de; do
            fresh "noh-$model"
            python3 -I "$VERIF_DIR/dualenergy/make_de.py" "$W/input.ini" noh --model "$model"
            run input.ini
         done ;;
      amr)
         fresh sv9-poisson
         POISSON_MG_EXE="${POISSON_MG_EXE:-$REPO_ROOT/exe/test_poisson_mg}"
         if (cd "$W" && mpirun -np "$NP" "$POISSON_MG_EXE" > log.txt 2>&1); then
            echo "   $(basename "$W"): done"
         else
            echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
         fi
         fresh sv13-poisson
         if (cd "$W" && mpirun -np "$NP" "$POISSON_MG_EXE" --isolated > log.txt 2>&1); then
            echo "   $(basename "$W"): done"
         else
            echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
         fi
         JEANS_EXE="${JEANS_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
         amr_jeans() { # amr_jeans <case dir> <make_jeans arguments>: one Jeans case (with --box: the box-refined grid)
            fresh "$1" ; shift
            python3 -I "$VERIF_DIR/selfgravity/make_jeans.py" "$W/input.ini" "$@"
            if (cd "$W" && mpirun -np "$NP" "$JEANS_EXE" input.ini > log.txt 2>&1); then
               echo "   $(basename "$W"): done"
            else
               echo "   $(basename "$W"): FAILED, see $W/log.txt" ; FAILED=1
            fi
         }
         for cells in 8 16 32; do amr_jeans "sv10-2d-$cells" jeans2d --cells "$cells" --box ; done
         for cells in 8 16 32; do amr_jeans "sv10-uniform-$cells" jeans2d --cells "$cells" ; done
         for cfl in 0.4 0.2 0.1; do amr_jeans "sv10-collapse-$cfl" jeansc2d --cells 16 --cfl "$cfl" --box ; done
         for cfl in 0.4 0.2 0.1; do
            amr_jeans "sv10-stage-$cfl" jeansc2d --cells 16 --cfl "$cfl" --box --reflux-cadence stage
         done
         amr_jeans sv11-regrid jeansc2d --cells 8 --regrid 5 --reflux-cadence stage --it-save 10 ;;
      mhd)
         M="$VERIF_DIR/mhd"
         fresh brio-wu
         "$VENV_PY" "$M/riemann/make_riemann.py" "$VERIF_DIR/sod/sod-x.ini" "$W/input.ini" --case bw --cells 512 \
            --divergence-control glm --rho-floor 1e-10 --p-floor 1e-10 ; run input.ini
         fresh rj4d
         "$VENV_PY" "$M/riemann/make_riemann.py" "$VERIF_DIR/sod/sod-x.ini" "$W/input.ini" --case rj4d --cells 512 \
            --divergence-control none --rho-floor 1e-10 --p-floor 1e-10 ; run input.ini
         fresh glm-pulse
         "$VENV_PY" "$M/glm-pulse/make_glm_pulse.py" "$VERIF_DIR/sod/sod-x.ini" "$W/glm-pulse.ini" --amplitude 0.1 \
            --cells 128 --divergence-control glm ; run glm-pulse.ini
         fresh orszag-tang
         "$VENV_PY" "$M/orszag-tang/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$W/orszag-tang.ini" \
            --cells 128 ; run orszag-tang.ini
         fresh orszag-tang-amr
         "$VENV_PY" "$M/orszag-tang/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$W/orszag-tang.ini" \
            --cells 32 --refine-box 0.25 0.25 0.75 0.75 ; run orszag-tang.ini
         fresh rotor
         "$VENV_PY" "$M/rotor/make_rotor.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$W/rotor.ini" --cells 256
         run rotor.ini
         fresh field-loop
         "$VENV_PY" "$M/field-loop/make_field_loop.py" "$VERIF_DIR/vortex/vortex-n064.ini" "$W/field-loop.ini" \
            --cells 128 ; run field-loop.ini
         fresh rj2a-hlld-prim
         "$VENV_PY" "$M/rj2a/make_rj2a.py" "$VERIF_DIR/sod/sod-x.ini" "$W/rj2a-x.ini" --axis x --cells 256 \
            --divergence-control none
         NUMERICS=hlld NUMERICS_SOLVER=hlld NUMERICS_RECON=primitive NUMERICS_CORRECTION=6th NUMERICS_SENSOR=weno \
            numerics_apply "$W/rj2a-x.ini"
         run rj2a-x.ini ;;
      rotated)
         for angle in 63 45; do
            fresh "rst-$angle-256"
            "$VENV_PY" "$VERIF_DIR/mhd/rotated-shock-tube/make_rotated_shock_tube.py" "$VERIF_DIR/vortex/vortex-n064.ini" \
               "$W/rst.ini" --angle "$angle" --cells 256 ; run rst.ini
            fresh "rst-$angle-ref"
            "$VENV_PY" "$VERIF_DIR/mhd/rotated-shock-tube/make_rotated_shock_tube.py" "$VERIF_DIR/sod/sod-x.ini" \
               "$W/rst.ini" --angle "$angle" --cells 256 --reference ; run rst.ini
         done ;;
      quadtree)
         for tree in octree quad-nk1; do
            fresh "ot-$tree"
            ratio=8 ; [[ $tree == quad-nk1 ]] && ratio=4
            "$VENV_PY" "$VERIF_DIR/mhd/orszag-tang/make_orszag_tang.py" "$VERIF_DIR/vortex/vortex-n064.ini" \
               "$W/orszag-tang.ini" --cells 32 --time-max 0.2 --refine-box 0.25 0.25 0.75 0.75 --ratio "$ratio"
            run orszag-tang.ini
         done ;;
      regrid)
         for prol in conservative linear; do
            fresh "regrid-vortex-$prol"
            python3 -I "$VERIF_DIR/regrid/make_regrid.py" vortex "$W/input.ini" --prolongation "$prol"
            sed -i -E 's/^(\s*it_save\s*=).*/\1 5/' "$W/input.ini"
            grep -qE '^\s*it_save\s*=\s*5\b' "$W/input.ini" || { echo "   it_save not set" ; FAILED=1 ; }
            run input.ini
         done ;;
      multirealm)
         mr() { # mr <case dir> <manifest or input> [<realm inputs>]: copy the inputs of multirealm/ and run the first
            fresh "$1" ; shift
            for f in "$@"; do cp "$VERIF_DIR/multirealm/$f" "$W/" 2>/dev/null || cp "$VERIF_DIR/sod/$f" "$W/" ; done
            run "$1"
         }
         mr sod-x sod-x.ini
         mr sod-2realm sod-2realm.ini sod-2realm-r1.ini sod-2realm-r2.ini
         mr sod-amr sod-amr.ini
         mr sod-amr-refined sod-amr-refined.ini sod-amr-refined-r1.ini sod-amr-refined-r2.ini ;;
      seams)
         as_run() { # as_run <case dir> <make_amr_shock arguments>: one Sod crossing
            fresh "$1" ; shift
            "$VENV_PY" "$VERIF_DIR/amr-shock/make_amr_shock.py" sod "$W/sod.ini" "$@" ; run sod.ini
         }
         as_run as-c2f --cells 200 --side c2f
         as_run as-f2c --cells 200 --side f2c
         as_run as-uniform-200 --cells 200 --uniform
         as_run as-uniform-400 --cells 400 --uniform
         for fill in tricubic restriction-compatible; do
            for restr in mean lagrange6; do
               as_run "e0-$fill-$restr-on" --cells 200 --side c2f --seam-ghost-fill "$fill" --seam-restriction "$restr"
               as_run "e0-$fill-$restr-off" --cells 200 --side c2f --seam-ghost-fill "$fill" --seam-restriction "$restr" \
                  --no-reflux
            done
            as_run "e0-$fill-400" --cells 400 --side c2f --seam-ghost-fill "$fill"
         done
         for fill in quintic weno6 quintic-mp; do
            for side in c2f f2c; do
               as_run "d-$fill-$side" --cells 200 --side "$side" --seam-ghost-fill "$fill" --seam-restriction lagrange6
            done
         done
         for side in c2f f2c; do
            as_run "e2-troubled-$side" --cells 200 --side "$side" --seam-ghost-fill quintic --seam-restriction lagrange6 \
               --reflux-policy troubled
         done ;;
      blast77)
         for scheme in split hlld; do
            fresh "blast77-$scheme"
            python3 -I "$VERIF_DIR/regrid/make_regrid.py" blast-av "$W/blast-av.ini" --base-level 3 --max-level 3 \
               --frequency 0 --time-max 0.01
            sed -i -E 's/^(\s*it_save\s*=).*/\1 20/' "$W/blast-av.ini"
            if [[ $scheme == hlld ]]; then # as regrid/accuracy.sh riemann(): weno-riemann, HLLD, 6th-order correction
               sed -i -e 's/^scheme_space *=.*/scheme_space = weno-riemann/' \
                      -e 's/^reconstruction_variables *=.*/reconstruction_variables = characteristic\nriemann_solver = hlld\nflux_correction = 6th\nflux_correction_sensor = weno/' \
                      "$W/blast-av.ini"
            fi
            run blast-av.ini
         done ;;
      divb)
         DIVB="$VERIF_DIR/amr-divb/make_amr_divb.py"
         fresh ad-cpaw-n64
         "$VENV_PY" "$DIVB" cpaw "$W/cpaw.ini" --cells 64 ; run cpaw.ini
         for n in 64 128; do
            coarse=$(python3 -c "import math; print(repr(math.sqrt(5.0) / $n))")
            fresh "g2-cpaw-coarse-$n"
            "$VENV_PY" "$DIVB" cpaw "$W/cpaw.ini" --cells "$n" --glm-damping-length "$coarse" ; run cpaw.ini
            fresh "g2-cpaw-mincell-$n"
            "$VENV_PY" "$DIVB" cpaw "$W/cpaw.ini" --cells "$n" --glm-damping-length min-cell ; run cpaw.ini
         done
         for cl in glm eglm; do
            fresh "loop-$cl-n64"
            "$VENV_PY" "$DIVB" loop "$W/loop.ini" --cells 64 --cleaning "$cl" --history 1 --it-save 20 ; run loop.ini
            fresh "loop-$cl-n64-noreflux"
            "$VENV_PY" "$DIVB" loop "$W/loop.ini" --cells 64 --cleaning "$cl" --history 1 --it-save 20 --no-reflux
            run loop.ini
            fresh "loop-$cl-uniform-n128"
            "$VENV_PY" "$DIVB" loop "$W/loop.ini" --cells 128 --cleaning "$cl" --uniform --history 1 --it-save 20
            run loop.ini
         done ;;
      rj2a-shock)
         for side in c2f f2c; do
            fresh "as-rj2a-$side"
            "$VENV_PY" "$VERIF_DIR/amr-shock/make_amr_shock.py" rj2a "$W/rj2a.ini" --cells 256 --side "$side"
            run rj2a.ini
         done
         fresh as-rj2a-uniform
         "$VENV_PY" "$VERIF_DIR/amr-shock/make_amr_shock.py" rj2a "$W/rj2a.ini" --cells 256 --uniform ; run rj2a.ini ;;
   esac
done
exit $FAILED
