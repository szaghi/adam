#!/usr/bin/env bash
# FLUME verification EO, P4: the isothermal and barotropic gases under gravity (issue #100, M13 P4).
#
# Why: the barotropic gases have no energy equation (D-M13-1), so gravity must move their momentum only, and the
# sinks and the Jeans marker must take c^2 = dP/drho and the thermal energy from the equation of state, not from the
# energy slot (D-M13-4). make_eos_gravity.py writes the cases (the gamma-law Jeans, sinks and box inputs with the gas
# of [eos] in place of gamma); the barotropic cases use POLY of eos_oracle.py on a background where dP/drho and P / rho
# differ. The legs:
#   eg1  free fall: a uniform periodic box under the constant g_x = 0.5 (isothermal; barotropic rho0 = 5): rho uniform,
#        rho u = rho g t (1e-12, relative), rE exactly 0 at every checkpoint;
#   eg2  the stable Jeans wave (lambda = lambda_J / 2, w^2 = c^2 k^2 - 4 pi G rho0 with c^2 = dP/drho), one period,
#        32^2, 64^2, 128^2: the L1 order of the finest pair >= 1.8 (the potential path is second order), isothermal and
#        barotropic (rho0 = 5: c^2 = P'(5) = 2.0198, where P / rho = 1.4427 would miss the phase by 1.3 rad);
#   eg3  the unstable Jeans mode (lambda = 1.5 lambda_J), 64^2: the growth rate within 1e-2 of sigma = sqrt(5/4) c k,
#        and the mirror symmetry within 1e-10, both gases;
#   eg4  static AMR (the central box refined, 2:1 seams, the multigrid solver, the seam reflux): the isothermal stable
#        wave at 32^2 .. 128^2 base, order >= 1.8; the barotropic collapse (amplitude 1e-2, 64^2 base): the momentum
#        within 1e-14, the mass within 1e-13, the energy integral exactly 0 at every step (no seam gravity work);
#   eg5  the Jeans marker (5.4 cells per Jeans length, a regrid at step 4, 32^3 periodic box at rest): isothermal at
#        5 cells per lambda_J, refined; barotropic at 6 by dP/drho (4.65 by P / rho), not refined;
#   eg6  sinks, isothermal, with the density-only positivity limiter: a Jeans-stable box at 1.2 rho_thr forms none
#        (the virial thermal term E_th = 3/2 int P dV; with the energy slot, E_th = 0, every cell would be a candidate);
#        accretion from a cold background (mass to 1e-13, momentum within 5 % of the exchange); the overdense singular
#        isothermal sphere (Federrath et al. 2010, A = 29.3, now truly isothermal) accreting within 6 % of Shu's
#        m_0(A) a^3 / G; a supercritical Bonnor-Ebert sphere (xi = 10, density doubled) forming one sink at the centre;
#   eg7  Bonnor-Ebert stability (closed cube [-1, 1]^3, 32^3, isolated self-gravity, to t = 4, about six free-fall
#        times of the centre): the sub-critical spheres xi = 3 and 6 (at the walls) stay within 0.25 of their central
#        density and below Mach 0.25; the super-critical xi = 9 departs from it by at least 0.5. (In a cube the
#        Lane-Emden profile is no exact equilibrium: the corner gas, 46 % of the sphere's mass at xi = 3, pulls through
#        its cubic multipoles, 3 % of the gravity at r = 0.6 and 10 % at 0.9: a stability test, not a static one.)
#   eg8  the contract: bonnor-ebert with the barotropic gas, bonnor-ebert without self-gravity, boss-bodenheimer with
#        the gamma-law gas, the singular isothermal sphere with the barotropic gas are refused;
#   eg9  CPU and FNL agree (run eg1, eg2, eg4 and eg6 with each executable first): the free fall (eos_oracle agree), the
#        64^2 Jeans waves and the AMR collapse (density within 1e-6 of the mode amplitude, jeans_oracle compare), the
#        singular isothermal sphere and the Bonnor-Ebert collapse (si_oracle agree, the centred sinks by their masses
#        and distance from the centre). Not the cold accretion: the limiter clips cells there every step, and the two
#        backends, equal to round-off for five steps, part at 1e-6 by step 15 (the limited runs are not bitwise across
#        backends, as the M3-P6 limited regression cases, which have per-backend goldens).
#   eg10 the Boss-Bodenheimer cloud (make_doc_runs.sh's eg-bb: 32 cells per block, the Jeans marker to 2 levels more,
#        sinks), 300 steps (it_max-driven, so a stalled run ends too): still advancing at the end (the last 50 steps at
#        least 0.1 of the mean step) and a sink formed. Its envelope, expanding into an ambient 100 times thinner, meets
#        the 2:1 seams at the corners of the refined cube; a limiter that extrapolates the velocity m / rho of a cell or
#        ghost it takes to the density floor drives dt to 0 there: the seam ghost blend stalled the run at t = 1.3721,
#        the cell limiter without its kinetic chord at 1.46 (0 and 0.007 of the mean step; 0.37 now). And, from the
#        checkpoints every 20 steps, the corners stay like the gas around (issue #101): in the shell 1.55 < r < 1.85
#        the cells near an edge of the refined cube at most twice as fast and at least half as dense as the others
#        (the tricubic seam ghosts of non-positive density, made zero-gradient copies by the blend, drove them to
#        23 c0 and 1/25 of the density while the run kept advancing).
#
# Usage: ./gravity.sh [--legs LIST]   (legs eg1..eg10, all; default all but eg9)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu_fft: self-gravity needs the FFT build); the
# caller owns the matching environment. EG_NP sets the ranks of the 2-D and 3-D cases (default 4; the free fall runs 2).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu_fft}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
EOS="$CASE_DIR/eos_oracle.py"
JEANS="$VERIF_DIR/selfgravity/jeans_oracle.py"
STAR="$VERIF_DIR/selfgravity/star_oracle.py"
SINKS="$VERIF_DIR/sinks/si_oracle.py"
source "$VERIF_DIR/legs.sh"
legs_declare "eg1 eg2 eg3 eg4 eg5 eg6 eg7 eg8 eg9 eg10" "eg1 eg2 eg3 eg4 eg5 eg6 eg7 eg8 eg10"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "gravity.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve gravity.sh
if [[ ! -x "$EXE" ]]; then echo "gravity.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; fi
EXE="$(realpath "$EXE")"
TAG="$(basename "$EXE")"
FAILED=0
NP="${EG_NP:-4}"

run() { # run <work dir> [ranks]: run the input of the work dir; 0 on success
   (cd "$1" && mpirun -np "${2:-$NP}" "$EXE" input.ini < /dev/null > log.txt 2>&1)
}
make_case() { # make_case <name> <make_eos_gravity.py case and options>: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   "$VENV_PY" -I "$CASE_DIR/make_eos_gravity.py" "$W/input.ini" "$@"
}
check() { # check <work dir> <ranks> <oracle> <arguments>: run the case, then the oracle
   local w=$1 np=$2 oracle=$3 ; shift 3
   if run "$w" "$np"; then
      "$VENV_PY" "$oracle" "$@" || FAILED=1
   else
      echo "   the run failed: see $w/log.txt  FAIL" ; FAILED=1
   fi
}
expect_refused() { # expect_refused <name> <message regex> <make_eos_gravity.py case and options>: one contract case
   local name=$1 msg=$2 ; shift 2
   echo "-- eg8: $name"
   make_case "eg8-$name" "$@"
   if run "$W" 1; then
      echo "   accepted, expected a refusal: see $W/log.txt  FAIL" ; FAILED=1
   elif grep -aqE "$msg" "$W/log.txt"; then
      echo "   refused as expected"
   else
      echo "   failed without the expected message ($msg): see $W/log.txt  FAIL" ; FAILED=1
   fi
}
ladder() { # ladder <label> <make_eos_gravity.py case and options>: the 8, 16, 32 cells-per-block runs, then the order
   local label=$1 n ; shift
   local works=()
   for n in 8 16 32; do
      echo "-- $label, $((4 * n))^2"
      make_case "$label-$n" "$@" --cells "$n"
      if run "$W"; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
   done
   if [[ ${#works[@]} -eq 3 ]]; then "$VENV_PY" "$JEANS" wave "${works[@]}" --min-order 1.8 || FAILED=1 ; fi
}
agree_runs() { # agree_runs <case> <oracle> <arguments before the two runs> [-- <arguments after them>]
   local case_=$1 oracle=$2 ; shift 2
   local cpu="$CASE_DIR/work-adam_flume_cpu_fft-$case_" fnl="$CASE_DIR/work-adam_flume_fnl_fft-$case_"
   local pre=() post=()
   while [[ $# -gt 0 && $1 != "--" ]]; do pre+=("$1") ; shift ; done
   [[ $# -gt 0 ]] && shift
   post=("$@")
   echo "-- eg9: CPU and FNL agree, $case_"
   if [[ -f "$cpu/log.txt" && -f "$fnl/log.txt" ]]; then
      "$VENV_PY" "$oracle" "${pre[@]}" "$cpu" "$fnl" "${post[@]}" || FAILED=1
   else
      echo "   missing a run: run its leg with both executables first" ; FAILED=1
   fi
}
echo ">> EO P4: the barotropic gases under gravity ($TAG)"
for leg in $LEGS; do
   case $leg in
      eg1)
         for eos in isothermal barotropic; do
            echo "-- eg1: free fall, $eos"
            make_case "eg1-freefall-$eos" freefall --eos "$eos"
            check "$W" 2 "$EOS" freefall "$W"
         done ;;
      eg2)
         ladder eg2-jeans-isothermal jeans --eos isothermal
         ladder eg2-jeans-barotropic jeans --eos barotropic ;;
      eg3)
         for eos in isothermal barotropic; do
            echo "-- eg3: the unstable Jeans mode, $eos"
            make_case "eg3-jeansu-$eos" jeansu --eos "$eos" --cells 16 --it-save 5
            if run "$W"; then
               "$VENV_PY" "$JEANS" growth "$W" --tol 1e-2 || FAILED=1
               "$VENV_PY" "$JEANS" symmetry "$W" --tol 1e-10 || FAILED=1
            else
               echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1
            fi
         done ;;
      eg4)
         ladder eg4-jeans-box jeans --eos isothermal --box
         echo "-- eg4: the barotropic collapse on the static AMR box"
         make_case eg4-collapse-box jeansu --eos barotropic --cells 16 --box --amplitude 1e-2
         check "$W" "$NP" "$JEANS" conservation "$W" --tol-momentum 1e-14 --tol-mass 1e-13
         "$VENV_PY" "$EOS" inert-history "$W" || FAILED=1 ;;
      eg5)
         echo "-- eg5: the Jeans marker, isothermal at 5 cells per Jeans length"
         make_case eg5-marker-isothermal marker --eos isothermal
         check "$W" "$NP" "$EOS" refined "$W" --expect yes
         echo "-- eg5: the Jeans marker, barotropic at 6 cells by dP/drho (4.65 by P / rho)"
         make_case eg5-marker-barotropic marker --eos barotropic
         check "$W" "$NP" "$EOS" refined "$W" --expect no ;;
      eg6)
         echo "-- eg6: a Jeans-stable isothermal box above the threshold forms no sink"
         make_case eg6-stable stable --it-max 20
         check "$W" "$NP" "$SINKS" none "$W"
         echo "-- eg6: accretion from a cold isothermal background"
         make_case eg6-accrete accrete
         check "$W" "$NP" "$SINKS" accrete "$W"
         echo "-- eg6: the overdense singular isothermal sphere (Shu's A-family)"
         make_case eg6-shu shu
         check "$W" "$NP" "$SINKS" shu "$W"
         echo "-- eg6: a supercritical Bonnor-Ebert sphere forms one central sink"
         make_case eg6-collapse collapse
         check "$W" "$NP" "$SINKS" collapse "$W" ;;
      eg7)
         for xi in 3 6; do
            echo "-- eg7: the sub-critical Bonnor-Ebert sphere, xi = $xi at the walls"
            make_case "eg7-be-$xi" be --xi "$xi" --time-max 4.0 --it-save 10
            check "$W" "$NP" "$STAR" equilibrium "$W" --tol-centre 0.25 --tol-mach 0.25 --tol-energy inf
         done
         echo "-- eg7: the super-critical Bonnor-Ebert sphere, xi = 9 at the walls"
         make_case eg7-be-9 be --xi 9 --time-max 4.0 --it-save 10
         check "$W" "$NP" "$STAR" unstable "$W" --min-departure 0.5 ;;
      eg8)
         expect_refused be-barotropic "bonnor-ebert needs \[eos\]\.\(type\) = isothermal" be --cells 16 \
            --set "eos:type=barotropic"
         expect_refused be-no-gravity "bonnor-ebert needs \[gravity\]\.\(self_gravity\) = \.true\." be --cells 16 \
            --set "gravity:self_gravity=.false." --set "gravity:g_z=-1.0" --unset "gravity:G" \
            --unset "gravity:poisson_bc"
         expect_refused bb-ideal "boss-bodenheimer needs \[eos\]\.\(type\) = isothermal or barotropic" bb --cells 16 \
            --set "eos:type=ideal" --set "physics:gamma=1.4"
         expect_refused sis-barotropic "singular-isothermal-sphere is not defined for" shu --cells 16 \
            --set "eos:type=barotropic" ;;
      eg9)
         for eos in isothermal barotropic; do
            agree_runs "eg1-freefall-$eos" "$EOS" agree
            agree_runs "eg2-jeans-$eos-16" "$JEANS" compare -- --tol 1e-6
         done
         agree_runs eg4-collapse-box "$JEANS" compare -- --tol 1e-6
         agree_runs eg6-shu "$SINKS" agree -- --mass-only
         agree_runs eg6-collapse "$SINKS" agree -- --mass-only ;;
      eg10)
         echo "-- eg10: the Boss-Bodenheimer cloud still advances after 300 steps"
         make_case eg10-bb bb --cells 32 --amr 2 --it-max 300 --it-save 20
         check "$W" "$NP" "$EOS" advancing "$W" --last 50 --min-ratio 0.1
         "$VENV_PY" "$EOS" corners "$W" --max-speed-ratio 2 --min-density-ratio 0.5 || FAILED=1 ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "EO P4 FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "EO P4 PASSED ($TAG, legs: $LEGS)"
