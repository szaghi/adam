#!/usr/bin/env bash
# FLUME verification of the turbulence diagnostics and of the Taylor-Green initial condition (issue #119, T0.1, T0.2).
#
# Why: the ILES baseline of #119 (T0.3-T0.6) rests on the integrals of `[diagnostics] turbulence_history_save` and on the
# `taylor-green` initial condition, so both are checked first against what is known exactly.
#
#   tt1  the first history row at M 0.1 and M 1.25, n = 16, 32, 64, against the exact integrals (tgv_oracle.py t0): K within
#        1e-12 (the midpoint rule is exact on its trigonometric integrand), Z and eps at the 4th order of the gradients, eps_d
#        and pd zero (the field is solenoidal, and the same difference of x and y cancels), eps = eps_s + eps_d.
#   tt2  the same 32^3 run, 10 steps, on one rank and on N: the histories agree to 1e-11 of each column's scale (the
#        sums are in another order, so not bitwise).
#   tt3  the energy budget, Re 1600, M 0.1, to t = 2 (laminar), n = 32, 64 (--budget-ladder; 128 on the GPU: about 14 h on two
#        CPU ranks): the K drop against the integral of eps - pd. What is left is the scheme's numerical dissipation: it
#        must fall with the resolution.
#   tt4  spectra.py on the initial field: Parseval to 1e-12, the energy of the single-mode field in shell 2.
#   tt6  the pressure work: at M 0.8, n = 64, pd is large against what the scheme dissipates, and the budget to t = 2 must close
#        better with it than without it or with its sign flipped (at M 0.1 neither tt1, zero by symmetry, nor tt3 sees it).
#   ht1  the isotropic turbulence initial field (T0.5), M_t 0.6, k_0 4, n = 32, 64: every shell against E(k) within 1e-10 of the
#        peak, the mean K exact, the dilatation falling with n (hit_oracle.py).
#   ht2  10 steps of it, n = 32, on one rank and on N: unlike tt2 the blocks differ, so this does compare sums across ranks.
#   sg1  the subgrid-scale formulas (T1), test_sgs (mode test-sgs-gnu, SGS_TEST_EXE) against numpy, and the zeros of each model.
#   sg2  the budget with the model, Vreman, M 0.4, Re 1600, t = 2, n = 32, 64, 128 (--sgs-ladder): what eps_sgs leaves falls with n,
#        and the budget closes better with eps_sgs than without it or flipped. WALE and Sigma on 64^3 for the sign check.
#   sg3  no viscosity with each model, 32^3 to t = 1: eps zero, eps_sgs positive, mass and energy conserved.
#   sg4  one rank against N with a model on the random field (the same history, within 1e-11).
#   hy1  the hybrid central/WENO scheme (T2), inviscid Taylor-Green, M 0.1, to t = 2.5, n = 32 (--hy-n; 64 takes 100 minutes on two CPU ranks):
#        K less the pressure work against plain WENO, which must dissipate at least 10 times more; mass and energy conserved.
#   hy2  one rank against N with the hybrid scheme on the random field, where the sensor fires.
#   hy3  the central flux itself (test_kep, mode test-kep-gnu, KEP_TEST_EXE): order, kinetic energy identity, pressure equilibrium. hy1 cannot
#        tell a KEP flux from a non-KEP one: at 32^3 the plain central flux loses less K than this one.
#   hy4  a contact of density ratio 8 and 100 advected on a periodic square, t = 0.3: finite, with u and p within 1e-3 of 1. The Ducros ratio cannot
#        see it, the density jump in the sensor does.
#   hy5  T2.4, 2:1 box on the Taylor-Green vortex, 32^3 base, M 0.1, to t = 1, reflux on: conserved to round-off.
#   hy6  T2.4, physical faces and a shock: the regression sod-x input on the hybrid scheme against plain WENO: density within --hy-sod-tol L1, same shock.
#   hy7  T2.4, the regression inputs of runtime regridding (vortex-regrid, closed: 1e-11) and of a forest (sod-2realm) on the hybrid scheme: they
#        run and conserve. sod-2realm has open ends and the central flux sends a ripple ahead of the waves that reaches them (|u| 3e-8 there, 0 on WENO), so
#        real mass leaves: the bound is 1e-8, and the sum of the realms must drift as the single realm of hy6 does, within 1 %.
#   hy8  T2.5, a density wave of relative amplitude 0.2 advected once around, n = 32, 64, 128: pressure and velocity stay uniform; WENO reported.
#   hy9  T2.6, the model only where the scheme is central: on the random field, every cell flagged gives eps_sgs exactly 0, the default threshold
#        does not; one rank against N with a model on a 2:1 box.
#   hy10 T2.4 measure, no bound: the skin relaxed on 2:1 seams (hy5 case, central-seams). Prints what conservation does and the last history rows
#        of both, or that the run failed; never fails the script. On physical faces: the difference hy11 prints.
#   hy11 T2.4, positive test of the skin (the conservation legs pass with no flag at all): on a 2:1 box and with waves at open ends, relaxing the
#        skin must change the density.
#   hy12 T2.4, restart on a refined grid: from step 10 to 20 against the straight run, every cell's density within 1e-13.
#   hy13 T2.7 robustness, hybrid against WENO, density L1: a vortex after Jiang & Shu through a shock and sod-amr within 5e-3, Shu-Osher within 2e-4 (set
#        after the 1.7e-5 result). In 1-D there is no vorticity, so the Ducros ratio flags any varying flow: only the 2-D case exercises the sensor. A vortex
#        through a shock with a 2:1 box, skin on and off, runs and reports the difference.
#   wm1  T3.1 laminar channel, isothermal walls, Re 100, M 0.1, from the laminar profile: the driving force against 3 mu U / h^2 with order. Mass within 1e-6 only:
#        the no-slip isothermal wall of #65 lets mass through its face at O(h^3) (1.6e-7, 2.1e-8, 2.0e-9 at 16, 32, 64), the modelled wall sets it to 0;
#        one rank against N on the perturbed flow.
#   wm2  T3.3, T3.4: test_wall_model (laminar limit, log law, Crocco); the laminar channel on modelled walls: the force with order above 1.5 (the mirror cell and the linear profile are second order), mass conserved. The time 2 is before the diffusive steady state, t about 40, so the same transient on both sides.
#   wm3  T3.6: the same with 2:1 faces along the walls, half the wall refined: mass conserved, the force strictly between those of the uniform base and of the
#        uniform grid with the fine wall spacing (the flow does not depend on x, z). Runtime regridding is not run.
#   wm4  T3.3: the energy flux: a wall hotter than the gas, modelled against resolved isothermal walls: the same change of the total energy within 15 %; one rank against N on modelled walls.
#   tt5  (--agree TAG) this executable's tt2 history against the executable TAG (tt2 run): backends within 1e-11.
#
# Usage: ./check.sh [--np N] [--legs LIST] [--budget-ladder "32 64"] [--agree EXE-NAME]   FLUME_EXE overrides the executable.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
AGREE=""
BUDGET_LADDER="32 64"
SGS_LADDER="32 64"
HY_N=32
source "$VERIF_DIR/legs.sh"
legs_declare "tt1 tt2 tt3 tt4 tt6 ht1 ht2 sg1 sg2 sg3 sg4 hy1 hy2 hy3 hy4 hy5 hy6 hy7 hy8 hy9 hy10 hy11 hy12 hy13 wm1 wm2 wm3 wm4"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)    NP="$2" ; shift 2 ;;
      --legs)  legs_select "$2" ; shift 2 ;;
      --agree) AGREE="$2" ; shift 2 ;;
      --budget-ladder) BUDGET_LADDER="$2" ; shift 2 ;;
      --sgs-ladder) SGS_LADDER="$2" ; shift 2 ;;
      --hy-n) HY_N="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST, --agree EXE-NAME)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
[[ -x "$EXE" ]] || { echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; }
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")"

run() { # run <name> <np> <make_tgv options>: echo the work directory
   local name="$1" np="$2" work="$CASE_DIR/work-$TAG-$1" ; shift 2
   rm -rf "${work:?}" ; mkdir -p "$work"
   python3 -I "$CASE_DIR/make_tgv.py" "$work/input.ini" "$@"
   (cd "$work" && mpirun -np "$np" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: $name failed, see $work/log.txt" >&2 ; exit 1 ; }
   echo "$work"
}

loaded() { # loaded <work>: stop unless the last rank holds cells (one root block is not split on an unrefined one-level grid)
   local last ; last=$(ls "$1"/tgv-000000000-proc*.h5 | tail -1)
   (( $(stat -c %s "$last") > 1000 )) || { echo "check.sh: $last holds no block: the comparison across ranks would be vacuous" >&2 ; exit 1 ; }
}

if has_leg tt1; then
   for mach in 0.1 1.25; do
      echo "== tt1: first history row against the exact integrals, M $mach"
      ladder=()
      for n in 16 32 64; do
         w=$(run "tt1-m$mach-n$n" "$NP" --n $n --levels 1 --mach $mach --it-max 1 --cfl 0.1)
         ladder+=("$n:$w")
      done
      "$VENV_PY" "$CASE_DIR/tgv_oracle.py" t0 "${ladder[@]}" --mach $mach
   done
fi
if has_leg tt2; then
   echo "== tt2: one rank against $NP"
   a=$(run tt2-np1 1 --n 32 --levels 2 --it-max 10)
   b=$(run tt2-npN "$NP" --n 32 --levels 2 --it-max 10)
   loaded "$b"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
fi
if has_leg tt3; then
   echo "== tt3: energy budget to t = 2, Re 1600, M 0.1"
   ladder=()
   for n in $BUDGET_LADDER; do
      w=$(run "tt3-n$n" "$NP" --n $n --levels 2 --time-max 2.0)
      ladder+=("$n:$w")
   done
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" budget "${ladder[@]}" --time 2.0
fi
if has_leg tt4; then
   echo "== tt4: spectrum of the initial field"
   w=$(run tt4 "$NP" --n 32 --levels 1 --it-max 1 --it-save 1)
   "$VENV_PY" "$CASE_DIR/spectra.py" "$w" --step 0 --peak-share 0.99 --peak-shell 2
fi
if has_leg tt6; then
   echo "== tt6: pressure work, M 0.8"
   w=$(run tt6 "$NP" --n 64 --levels 2 --mach 0.8 --time-max 2.0)
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" pd "$w" --time 2.0
fi
HIT=(--mach 0.6 --hit-k0 4 --reynolds 100)
if has_leg ht1; then
   echo "== ht1: isotropic turbulence, the initial field"
   ladder=()
   for n in 32 64; do
      w=$(run "ht1-n$n" "$NP" --n $n --levels 1 "${HIT[@]}" --it-max 1 --it-save 1 --cfl 0.1)
      ladder+=("$n:$w")
   done
   "$VENV_PY" "$CASE_DIR/hit_oracle.py" "${ladder[@]}" --k0 4
fi
if has_leg ht2; then
   echo "== ht2: isotropic turbulence, one rank against $NP"
   a=$(run ht2-np1 1 --n 32 --levels 2 "${HIT[@]}" --it-max 10)
   b=$(run ht2-npN "$NP" --n 32 --levels 2 "${HIT[@]}" --it-max 10)
   loaded "$b"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
fi
if has_leg sg1; then
   echo "== sg1: subgrid-scale formulas"
   "$VENV_PY" "$CASE_DIR/sgs_oracle.py" formula "${SGS_TEST_EXE:-$REPO_ROOT/exe/test_sgs}"
fi
if has_leg sg2; then
   echo "== sg2: budget with the model"
   ladder=()
   for n in $SGS_LADDER; do
      w=$(run "sg2-vreman-n$n" "$NP" --n $n --levels 2 --mach 0.4 --sgs-model vreman --time-max 2.0)
      ladder+=("$n:$w")
   done
   "$VENV_PY" "$CASE_DIR/sgs_oracle.py" budget "${ladder[@]}"
   for m in wale sigma; do
      w=$(run "sg2-$m" "$NP" --n 64 --levels 2 --mach 0.4 --sgs-model $m --time-max 2.0)
      "$VENV_PY" "$CASE_DIR/sgs_oracle.py" budget "64:$w"
   done
fi
if has_leg sg3; then
   for m in vreman wale sigma; do
      echo "== sg3: no viscosity, $m"
      w=$(run "sg3-$m" "$NP" --n 32 --levels 2 --mach 0.4 --reynolds 0 --sgs-model $m --time-max 1.0)
      "$VENV_PY" "$CASE_DIR/sgs_oracle.py" inviscid "$w"
   done
fi
if has_leg sg4; then
   echo "== sg4: one rank against $NP with a model"
   a=$(run sg4-np1 1 --n 32 --levels 2 "${HIT[@]}" --sgs-model vreman --it-max 10)
   b=$(run sg4-npN "$NP" --n 32 --levels 2 "${HIT[@]}" --sgs-model vreman --it-max 10)
   loaded "$b"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
fi
if has_leg hy1; then
   echo "== hy1: inviscid kinetic energy, hybrid against WENO"
   a=$(run hy1-hybrid "$NP" --n $HY_N --levels 2 --reynolds 0 --scheme-space hybrid-kep-weno --time-max 2.5)
   b=$(run hy1-weno "$NP" --n $HY_N --levels 2 --reynolds 0 --time-max 2.5)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" "$a" "$b"
fi
if has_leg hy2; then
   echo "== hy2: one rank against $NP with the hybrid scheme"
   a=$(run hy2-np1 1 --n 32 --levels 2 "${HIT[@]}" --scheme-space hybrid-kep-weno --it-max 10)
   b=$(run hy2-npN "$NP" --n 32 --levels 2 "${HIT[@]}" --scheme-space hybrid-kep-weno --it-max 10)
   loaded "$b"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
fi
if has_leg hy3; then
   echo "== hy3: the central flux"
   "${KEP_TEST_EXE:-$REPO_ROOT/exe/test_kep}"
fi
if has_leg hy4; then
   for r in 8 100; do
      echo "== hy4: contact, density ratio $r"
      w=$(run "hy4-$r" "$NP" --n 64 --levels 2 --reynolds 0 --scheme-space hybrid-kep-weno --contact $r --cfl 0.3 --time-max 0.3 --it-save 100000)
      "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --contact "$w"
   done
fi
if has_leg hy5; then
   echo "== hy5: a 2:1 box with the hybrid scheme"
   w=$(run hy5 "$NP" --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --time-max 1.0)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --conserve "$w/tgv-conservation_history.dat"
fi
regression_run() { # regression_run <case> <work>: the case's inputs on the hybrid scheme, in <work>
   local src="$CASE_DIR/../../regression/$1" w="$2" ; rm -rf "${w:?}" ; mkdir -p "$w"
   cp "$src"/*.ini "$w"/
   sed -i -E 's/^(scheme_space\s*=\s*)weno\s*$/\1hybrid-kep-weno/' "$w"/*.ini
   grep -q "hybrid-kep-weno" "$w"/*.ini || { echo "check.sh: no weno scheme_space in $1" >&2 ; exit 1 ; }
   (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: $1 failed on the hybrid scheme, see $w/log.txt" >&2 ; exit 1 ; }
}
if has_leg hy6; then
   echo "== hy6: Sod on the hybrid scheme against WENO"
   w="$CASE_DIR/work-$TAG-hy6" ; regression_run sod-x "$w"
   u="$CASE_DIR/work-$TAG-hy6-weno" ; rm -rf "${u:?}" ; mkdir -p "$u" ; cp "$CASE_DIR/../../regression/sod-x/input.ini" "$u"/
   (cd "$u" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: sod-x failed, see $u/log.txt" >&2 ; exit 1 ; }
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --compare "$w" "$u"
fi
if has_leg hy7; then
   for c in vortex-regrid sod-2realm; do
      echo "== hy7: $c on the hybrid scheme"
      w="$CASE_DIR/work-$TAG-hy7-$c" ; regression_run "$c" "$w"
      if [[ $c == sod-2realm ]]; then tol=1e-8 ; else tol=1e-11 ; fi
      "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --conserve "$w"/*-conservation_history.dat --conserve-tol $tol
      if [[ $c == sod-2realm ]]; then
         [[ -f "$CASE_DIR/work-$TAG-hy6/sod-x-conservation_history.dat" ]] || { echo "check.sh: hy7 sod-2realm needs hy6 in the same legs" >&2 ; exit 1 ; }
         "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --same-drift "$CASE_DIR/work-$TAG-hy6/sod-x-conservation_history.dat" "$w"
      fi
   done
fi
if has_leg hy8; then
   for n in 32 64 128; do
      echo "== hy8: density wave, $n cells"
      a=$(run "hy8-$n" "$NP" --n $n --levels 2 --reynolds 0 --density-wave 0.2 --scheme-space hybrid-kep-weno --cfl 0.3 --time-max 1.0 --it-save 100000)
      b=$(run "hy8-$n-weno" "$NP" --n $n --levels 2 --reynolds 0 --density-wave 0.2 --cfl 0.3 --time-max 1.0 --it-save 100000)
      "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --density-wave "$a" "$b"
   done
fi
if has_leg hy9; then
   echo "== hy9: the model only where the scheme is central"
   a=$(run hy9-flagged "$NP" --n 32 --levels 2 "${HIT[@]}" --scheme-space hybrid-kep-weno --sgs-model vreman --hybrid-threshold 1e-30 --it-max 10)
   b=$(run hy9-free "$NP" --n 32 --levels 2 "${HIT[@]}" --scheme-space hybrid-kep-weno --sgs-model vreman --it-max 10)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --masked "$a" "$b"
   a=$(run hy9-np1 1 --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --sgs-model vreman --it-max 10)
   b=$(run hy9-npN "$NP" --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --sgs-model vreman --it-max 10)
   loaded "$b"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$a" "$b"
   # a mask missing the ghost cells gives the two blocks of an interface two SGS fluxes: flagged and free cells mixed, conserved
   w=$(run hy9-mix "$NP" --n 32 --levels 2 --amr-box 0.5 "${HIT[@]}" --scheme-space hybrid-kep-weno --sgs-model vreman --it-max 10)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --conserve "$w/tgv-conservation_history.dat"
fi
if has_leg hy10; then
   echo "== hy10: the skin relaxed on 2:1 seams"
   w="$CASE_DIR/work-$TAG-hy10-seams" ; rm -rf "${w:?}" ; mkdir -p "$w"
   python3 -I "$CASE_DIR/make_tgv.py" "$w/input.ini" --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --hybrid-skin central-seams --time-max 1.0
   if (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1); then
      "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --conserve "$w/tgv-conservation_history.dat" || true
      tail -1 "$w/tgv-turbulence_history.dat"; [[ -f "$CASE_DIR/work-$TAG-hy5/tgv-turbulence_history.dat" ]] && tail -1 "$CASE_DIR/work-$TAG-hy5/tgv-turbulence_history.dat"
   else echo "hy10: central-seams run FAILED, see $w/log.txt"; fi
fi
if has_leg hy11; then
   echo "== hy11: the skin is really WENO: relaxing it must change the result"
   a=$(run hy11-seams-on "$NP" --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --it-max 20)
   b=$(run hy11-seams-off "$NP" --n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno --hybrid-skin central-seams --it-max 20)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --differ "$a" "$b"
   for s in weno central-bcs; do
      w="$CASE_DIR/work-$TAG-hy11-bcs-$s" ; rm -rf "${w:?}" ; mkdir -p "$w"
      cp "$CASE_DIR/../../regression/sod-x/input.ini" "$w"/
      # past t = 0.2 the waves reach the open ends, where the skin acts
      sed -i -E "s/^(scheme_space\s*=\s*)weno\s*$/\1hybrid-kep-weno\nhybrid_skin = $s/; s/^(time_max\s*=\s*).*/\11.0/" "$w/input.ini"
      (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: hy11 $s failed, see $w/log.txt" >&2 ; exit 1 ; }
   done
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --differ "$CASE_DIR/work-$TAG-hy11-bcs-weno" "$CASE_DIR/work-$TAG-hy11-bcs-central-bcs"
fi
if has_leg hy12; then
   echo "== hy12: restart of the hybrid scheme on a refined grid"
   # one restart file, overwritten: step 10 of a first run, continued to 20, against 20 straight
   A=(--n 32 --levels 2 --amr-box 0.5 --scheme-space hybrid-kep-weno)
   a=$(run hy12-straight "$NP" "${A[@]}" --it-max 20)
   b=$(run hy12-restart "$NP" "${A[@]}" --it-max 10 --restart-save 10)
   python3 -I "$CASE_DIR/make_tgv.py" "$b/input.ini" "${A[@]}" --it-max 20 --restart-save 10 --restart
   (cd "$b" && mpirun -np "$NP" "$EXE" input.ini > log2.txt 2>&1) || { echo "check.sh: hy12 restart failed, see $b/log2.txt" >&2 ; exit 1 ; }
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --compare --compare-max --compare-tol 1e-13 "$b" "$a"
fi
if has_leg hy13; then
   SV=(--n 256 --levels 3 --reynolds 0 --shock-vortex --cfl 0.3 --time-max 0.35 --it-save 100000)
   echo "== hy13: shock-vortex, hybrid against WENO"
   a=$(run hy13-hybrid "$NP" "${SV[@]}" --scheme-space hybrid-kep-weno)
   b=$(run hy13-weno "$NP" "${SV[@]}")
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --compare "$a" "$b"
   for scheme in hybrid-kep-weno weno; do
      w="$CASE_DIR/work-$TAG-hy13-so-$scheme" ; rm -rf "${w:?}" ; mkdir -p "$w"
      "$VENV_PY" "$CASE_DIR/make_shu_osher.py" "$CASE_DIR/../../regression/sod-x/input.ini" "$w/input.ini" "$scheme"
      (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: Shu-Osher $scheme failed, see $w/log.txt" >&2 ; exit 1 ; }
   done
   echo "== hy13: Shu-Osher, hybrid against WENO"
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --compare --compare-tol 2e-4 "$CASE_DIR/work-$TAG-hy13-so-hybrid-kep-weno" "$CASE_DIR/work-$TAG-hy13-so-weno"
   echo "== hy13: a shock through the 2:1 seams of sod-amr, hybrid against WENO"
   w="$CASE_DIR/work-$TAG-hy13-amr" ; regression_run sod-amr "$w"
   u="$CASE_DIR/work-$TAG-hy13-amr-weno" ; rm -rf "${u:?}" ; mkdir -p "$u" ; cp "$CASE_DIR/../../regression/sod-amr/input.ini" "$u"/
   (cd "$u" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: sod-amr failed, see $u/log.txt" >&2 ; exit 1 ; }
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --compare "$w" "$u"
   echo "== hy13: shock-vortex through a 2:1 box, skin on and off"
   a=$(run hy13-box-on "$NP" "${SV[@]}" --scheme-space hybrid-kep-weno --amr-box 0.5)
   b=$(run hy13-box-off "$NP" "${SV[@]}" --scheme-space hybrid-kep-weno --amr-box 0.5 --hybrid-skin central-seams)
   "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --differ "$a" "$b"
fi
channel_run() { # channel_run <name> <np> <make_channel options>: echo the work directory
   local name="$1" np="$2" work="$CASE_DIR/work-$TAG-$1" ; shift 2
   rm -rf "${work:?}" ; mkdir -p "$work"
   python3 -I "$CASE_DIR/make_channel.py" "$work/input.ini" "$@"
   (cd "$work" && mpirun -np "$np" "$EXE" input.ini > log.txt 2>&1) || { echo "check.sh: $name failed, see $work/log.txt" >&2 ; exit 1 ; }
   echo "$work"
}
if has_leg wm1; then
   echo "== wm1: laminar channel"
   a=$(channel_run wm1-16 "$NP" --ny 16 --time-max 2.0)
   b=$(channel_run wm1-32 "$NP" --ny 32 --time-max 2.0)
   c=$(channel_run wm1-64 "$NP" --ny 64 --time-max 2.0)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" laminar --mass-tol 1e-6 "$a" "$b" "$c"
   a=$(channel_run wm1-np1 1 --ny 32 --perturbation 0.1 --it-max 10)
   b=$(channel_run wm1-npN "$NP" --ny 32 --perturbation 0.1 --it-max 10)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" agree "$a" "$b"
fi
if has_leg wm2; then
   echo "== wm2: the wall model"
   [[ "${WM_TEST_EXE:-}" == skip ]] || "${WM_TEST_EXE:-$REPO_ROOT/exe/test_wall_model}"  # a host program: WM_TEST_EXE=skip on the GPU side
   a=$(channel_run wm2-16 "$NP" --ny 16 --wall wall-modelled --time-max 2.0)
   b=$(channel_run wm2-32 "$NP" --ny 32 --wall wall-modelled --time-max 2.0)
   c=$(channel_run wm2-64 "$NP" --ny 64 --wall wall-modelled --time-max 2.0)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" laminar --min-order 1.5 --tol 0.05 "$a" "$b" "$c"
fi
if has_leg wm3; then
   echo "== wm3: 2:1 faces along the wall"
   a=$(channel_run wm3 "$NP" --nx 32 --ny 32 --nz 32 --wall wall-modelled --amr-box --time-max 0.5)
   b=$(channel_run wm3-coarse "$NP" --nx 32 --ny 32 --nz 32 --wall wall-modelled --time-max 0.5)
   c=$(channel_run wm3-fine "$NP" --nx 32 --ny 64 --nz 32 --wall wall-modelled --time-max 0.5)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" between "$a" "$b" "$c"
fi
if has_leg wm4; then
   echo "== wm4: wall heat flux and ranks on modelled walls"
   a=$(channel_run wm4-model "$NP" --ny 64 --wall wall-modelled --wall-heat 1.5 --time-max 1.0)
   b=$(channel_run wm4-resolved "$NP" --ny 64 --wall-heat 1.5 --time-max 1.0)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" energy --tol 0.15 "$a" "$b"
   a=$(channel_run wm4-np1 1 --ny 32 --wall wall-modelled --perturbation 0.1 --it-max 10)
   b=$(channel_run wm4-npN "$NP" --ny 32 --wall wall-modelled --perturbation 0.1 --it-max 10)
   "$VENV_PY" "$CASE_DIR/channel_oracle.py" agree "$a" "$b"
fi
if [[ -n "$AGREE" ]]; then
   echo "== tt5: $TAG against $AGREE"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-tt2-npN" "$CASE_DIR/work-$AGREE-tt2-npN"
   [[ -d "$CASE_DIR/work-$AGREE-hy2-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-hy2-npN" "$CASE_DIR/work-$AGREE-hy2-npN"
   [[ -d "$CASE_DIR/work-$AGREE-sg4-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-sg4-npN" "$CASE_DIR/work-$AGREE-sg4-npN"
   [[ -d "$CASE_DIR/work-$AGREE-ht2-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-ht2-npN" "$CASE_DIR/work-$AGREE-ht2-npN"
fi
echo "turbulence verification PASSED (legs: $LEGS) ($TAG)"
