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
legs_declare "tt1 tt2 tt3 tt4 tt6 ht1 ht2 sg1 sg2 sg3 sg4 hy1 hy2 hy3 hy4 hy5 hy6 hy7"
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
      if [[ $c == sod-2realm && -f "$CASE_DIR/work-$TAG-hy6/sod-x-conservation_history.dat" ]]; then
         "$VENV_PY" "$CASE_DIR/hybrid_oracle.py" --same-drift "$CASE_DIR/work-$TAG-hy6/sod-x-conservation_history.dat" "$w"
      fi
   done
fi
if [[ -n "$AGREE" ]]; then
   echo "== tt5: $TAG against $AGREE"
   "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-tt2-npN" "$CASE_DIR/work-$AGREE-tt2-npN"
   [[ -d "$CASE_DIR/work-$AGREE-hy2-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-hy2-npN" "$CASE_DIR/work-$AGREE-hy2-npN"
   [[ -d "$CASE_DIR/work-$AGREE-sg4-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-sg4-npN" "$CASE_DIR/work-$AGREE-sg4-npN"
   [[ -d "$CASE_DIR/work-$AGREE-ht2-npN" ]] && "$VENV_PY" "$CASE_DIR/tgv_oracle.py" agree "$CASE_DIR/work-$TAG-ht2-npN" "$CASE_DIR/work-$AGREE-ht2-npN"
fi
echo "turbulence verification PASSED (legs: $LEGS) ($TAG)"
