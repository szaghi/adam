#!/usr/bin/env bash
# FLUME verification EO: the isothermal and barotropic equations of state, `[eos] type = isothermal | barotropic`
# (issue #100, M13 P1, P2, P3; gravity, self-gravity, the Jeans marker and sinks, P4: gravity.sh).
#
# Why: the isothermal gas P = c0^2 rho drops the energy equation (D-M13-1: the state keeps its layout, the energy slot
# rE is inert, zero flux and zero source, held at 0), so its Riemann problem has an exact solution with no pressure
# iteration on the energy: every run is checked against eos_oracle.py, which builds that solution from the jump
# conditions independently of the code. make_eos.py writes the cases from the sod-x input. The legs:
#   eo1  the shock tube (rho 1 | 0.125, v 0.3 | -0.2, 400 cells, WENO-5 characteristic, SSP-33, CFL 0.5, t = 0.2): the
#        star plateau mean within 1e-3 of rho*, the right shock (mid-density crossing) within 1.5 cells, the contact
#        (mid transverse-velocity crossing) within 2 cells, rE exactly 0; the L1 errors reported;
#   eo2  two Mach-5 streams colliding (rho = 1, u = -+5, 400 cells, to t = 1): the dense layer rho* = M_1^2 = 26.96
#        within 1e-3 (its mean: the strong shocks leave a decaying, zero-mean oscillation train), both shocks within
#        1.5 cells;
#   eo3  the order: a right-running acoustic wave of amplitude 1e-6 on the periodic box, one period, N = 12, 24, 48
#        (two blocks, two ranks), SSP-54 at CFL 0.1 (the time error below the WENO-5 space error): the observed L1
#        order of the finest pair above 4.5;
#   eo4  the contract: what the isothermal gas does not define yet is refused with its own message ([cooling],
#        viscosity, [reference], immersed solids, the no-slip and hydrostatic walls, an IC whose density needs gamma),
#        and so are a malformed [eos] (euler-de or passive scalars with it, the gas keys cp/cv, a missing or
#        non-positive sound_speed, an unknown type); gravity, sinks and the positivity limiter are defined since P4
#        (gravity.sh), mhd-ideal since P5 (mhd.sh);
#   eo5  the inert slot and the mass: a wave of amplitude 0.3 on the periodic box, which breaks into a shock at
#        t ~ 1 / (c0 A k) = 0.53 (N = 48, two ranks, to t = 1): rE exactly 0 in every checkpoint, the total mass
#        constant to 1e-13;
#   eo6  the Riemann solvers: eo1 with weno-riemann and llf, hll, hllc (characteristic reconstruction, the 6th-order
#        flux correction), each against the exact solution as in eo1;
#   eo7  (P2) CPU and FNL agree (run eo1, eo2, eo3 and eo6 with each executable first): on the last checkpoint of the
#        tube, the collision, the N = 48 wave and the hllc tube, every field within 1e-10 of its scale (the floor of
#        AV av8): the density its largest magnitude, a momentum at least max(rho) c0 (its flux carries the pressure
#        c0^2 rho, so a 1e-6 wave's momentum has the round-off of an O(1) flux); rE exactly 0 on both.
# The barotropic legs (P3) on the convex test table of eos_oracle.py, POLY: c0 = 1, densities 2 8, exponents
# 1 1.4 5/3 (P continuous, c^2 = dP/drho jumping up at both breaks):
#   eb1  the tube rho 10 | 0.5, v 0 | 0.7, 1600 cells, t = 0.15: rho* = 2.602, the left fan crossing the break 8 splits
#        into two fans with a plateau at rho = 8; as eo1, plus that plateau (its interior cells, mean) within 1e-3;
#   eb2  the collision (1, 3.4) | (1.5, -2.6), 400 cells, t = 0.3: rho* = 8.056, both shocks from segment 1 or 2 into
#        segment 3 (shock speeds -0.466 and 0.612); the layer mean within 1e-3, both shocks within 1.5 cells. The frame
#        is (1, 4) | (1.5, -2) shifted by -0.6: there the left shock moves at 0.134 relative to the grid, and that slow
#        shock sheds the persistent post-shock noise of slowly moving shocks (Arora & Roe 1997), 4.6 % under HLL and
#        HLLC; in this frame every shock moves at least 0.36, and the noise is 0.47 % (a solver defect would not depend
#        on the frame);
#   eb3  the order on the Gamma = 1.4 segment: the eo3 wave on the background rho0 = 5 (c^2 = P'(5) = 2.0198), one
#        period, N = 12, 24, 48: the finest-pair order above 4.5;
#   eb4  the contract of a table: decreasing densities, exponents of another count, an exponent below 1, more than 8
#        segments, a non-positive reference_density, the gas keys, [reference] and no sound_speed are refused;
#   eb5  the isothermal limit: the one-segment table exponents = 1 with c0 = 1.3 against the isothermal gas with
#        c0 = 1.3 on the eo1 tube, t = 0.15: the two runs agree within 1e-10 (the eo7 rule; the barotropic path takes
#        P, c and the Roe speed from the table, the isothermal one in closed form);
#   eb6  the default table (Mathew & Federrath 2021, cgs, c0 = 0.2 km/s), non-convex at 3.84e-8 g cm^-3 (no exact
#        solution: composite waves): the periodic tube rho 1e-7 | 1e-8 g cm^-3 at rest on 1 cm, t = 3e-7 s, its left
#        fan running down through the non-convex break, N = 200, 400, 800 against N = 1600 (`si` WENO weights: the
#        data are far from order 1); the L1 distance to the reference shrinking at an order >= 0.5, rho > 0, rE = 0,
#        the total mass constant to 1e-13;
#   eb7  the Riemann solvers: eb2 with weno-riemann and llf, hll, hllc (the HLL speeds use the Roe secant speed), each
#        against the exact solution as in eb2;
#   eb8  CPU and FNL agree (run eb1, eb2, eb3 and eb7 with each executable first), as eo7.
#
# Usage: ./check.sh [--legs LIST]   (legs eo1..eo7, eb1..eb8, all; default all but eo7 and eb8)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "eo1 eo2 eo3 eo4 eo5 eo6 eo7 eb1 eb2 eb3 eb4 eb5 eb6 eb7 eb8" \
             "eo1 eo2 eo3 eo4 eo5 eo6 eb1 eb2 eb3 eb4 eb5 eb6 eb7"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh
if [[ ! -x "$EXE" ]]; then echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; fi
EXE="$(realpath "$EXE")"
TAG="$(basename "$EXE")"
FAILED=0

run() { # run <work dir> [ranks]: run the input of the work dir; 0 on success
   (cd "$1" && mpirun -np "${2:-1}" "$EXE" input.ini < /dev/null > log.txt 2>&1)
}
make_case() { # make_case <name> [make_eos.py options]: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   python3 -I "$CASE_DIR/make_eos.py" "$W/input.ini" "$@"
}
run_oracle() { # run_oracle <work dir> <ranks> <eos_oracle.py arguments>: run the case, then the oracle
   local w=$1 np=$2 ; shift 2
   if run "$w" "$np"; then
      "$VENV_PY" "$CASE_DIR/eos_oracle.py" "$@" || FAILED=1
   else
      echo "   the run failed: see $w/log.txt  FAIL" ; FAILED=1
   fi
}
expect_refused() { # expect_refused <work dir> <message regex>
   if run "$1"; then
      echo "   accepted, expected a refusal: see $1/log.txt  FAIL" ; FAILED=1
   elif grep -aqE "$2" "$1/log.txt"; then
      echo "   refused as expected"
   else
      echo "   failed without the expected message ($2): see $1/log.txt  FAIL" ; FAILED=1
   fi
}
refuse() { # refuse <name> <message regex> [make_eos.py options]: one contract case
   local name=$1 msg=$2 ; shift 2
   echo "-- eo4: $name"
   make_case "eo4-$name" "$@"
   expect_refused "$W" "$msg"
}

ISO="with \[eos\]\.\(type\)=isothermal"
POLY=(--eos barotropic --densities "2.0 8.0" --exponents "1.0 1.4 1.6666666666666667") # eos_oracle.py POLY
agree() { # agree <label> <case>...: the CPU and FNL work directories of each case
   local case_ label=$1 ; shift
   for case_ in "$@"; do
      echo "-- $label: CPU and FNL agree, $case_"
      cpu="$CASE_DIR/work-adam_flume_cpu-$case_" ; fnl="$CASE_DIR/work-adam_flume_fnl-$case_"
      if [[ -f "$cpu/log.txt" && -f "$fnl/log.txt" ]]; then
         "$VENV_PY" "$CASE_DIR/eos_oracle.py" agree "$cpu" "$fnl" || FAILED=1
      else
         echo "   missing a run: run the legs of $case_ with both executables first" ; FAILED=1
      fi
   done
}
echo ">> EO: the isothermal equation of state ($TAG)"
for leg in $LEGS; do
   case $leg in
      eo1)
         echo "-- eo1: the shock tube, WENO-5 characteristic"
         make_case eo1-tube --case tube
         run_oracle "$W" 1 riemann "$W" ;;
      eo2)
         echo "-- eo2: the Mach-5 collision"
         make_case eo2-collide --case collide --mach 5.0 --time-max 1.0
         run_oracle "$W" 1 riemann "$W" ;;
      eo3)
         works=()
         for n in 12 24 48; do
            echo "-- eo3: the acoustic wave, N = $n"
            make_case "eo3-wave-$n" --case wave --ni "$n" --cfl 0.1 --rk runge-kutta-ssp-54
            if run "$W" 2; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
         done
         if [[ ${#works[@]} -eq 3 ]]; then "$VENV_PY" "$CASE_DIR/eos_oracle.py" wave "${works[@]}" || FAILED=1 ; fi ;;
      eo4)
         refuse cooling "\[cooling\] is refused $ISO" --set "cooling:particle_mass=1.0"
         refuse viscosity "dissipative terms .* are refused with" --set "physics:viscosity=0.01"
         refuse reference "\[reference\] is refused $ISO" --set "reference:length=1.0"
         refuse solids "immersed solids are refused $ISO" --set "solids:number=1" --set "solid_1:name=cylinder" \
            --set "solid_1:bc_type=2" --set "solid_1:definition=analytical_circle" \
            --set "solid_1:circle_center_x=0.5" --set "solid_1:circle_center_y=0.5" \
            --set "solid_1:circle_center_z=0.5" --set "solid_1:circle_radius=0.1" --set "solid_1:circle_axis=z"
         refuse noslip "\(type\)=wall-noslip is refused $ISO" --set "bc_x_min:type=wall-noslip"
         refuse hydrostatic "\(type\)=hydrostatic is refused $ISO" --set "bc_x_min:type=hydrostatic"
         refuse vortex "\(type\)=isentropic-vortex is not defined for" \
            --set "initial_conditions:type=isentropic-vortex" \
            --set "initial_conditions:x0=0.5" --set "initial_conditions:y0=0.5" \
            --set "initial_conditions:radius=0.1" --set "initial_conditions:strength=1.0"
         refuse euler-de "isothermal needs \[physics\]\.\(physical_model\)=euler or mhd-ideal" \
            --set "physics:physical_model=euler-de"
         refuse scalars "isothermal needs \[physics\]\.\(physical_model\)=euler or mhd-ideal without passive_scalars" \
            --set "physics:passive_scalars=1"
         refuse gas "gamma, cp, cv are not used with \[eos\]\.\(type\) = isothermal" --keep-gas
         refuse no-speed "failed to load \[eos\]\.\(sound_speed\)" --unset "eos:sound_speed"
         refuse negative-speed "\[eos\]\.\(sound_speed\) must be positive" --set "eos:sound_speed=-1.0"
         refuse unknown "unknown \[eos\]\.\(type\) \"polytropic\"" --set "eos:type=polytropic" ;;
      eo5)
         echo "-- eo5: the steepened wave, the inert slot and the mass"
         make_case eo5-steep --case wave --ni 48 --set "initial_conditions:wave_amplitude=0.3"
         run_oracle "$W" 2 inert "$W" ;;
      eo6)
         for solver in llf hll hllc; do
            echo "-- eo6: the shock tube, weno-riemann $solver"
            make_case "eo6-$solver" --case tube --scheme weno-riemann --solver "$solver"
            run_oracle "$W" 1 riemann "$W"
         done ;;
      eo7)
         agree eo7 eo1-tube eo2-collide eo3-wave-48 eo6-hllc ;;
      eb1)
         echo "-- eb1: the barotropic tube across the breaks"
         make_case eb1-tube --case tube "${POLY[@]}" --left "10.0 0.0 0.0" --right "0.5 0.0 0.7" --time-max 0.15 \
            --ni 1600
         run_oracle "$W" 1 riemann "$W" ;;
      eb2)
         echo "-- eb2: the barotropic collision"
         make_case eb2-collide --case collide "${POLY[@]}" --left "1.0 3.4" --right "1.5 -2.6" --time-max 0.3
         run_oracle "$W" 1 riemann "$W" ;;
      eb3)
         period=$("$VENV_PY" -I -c 'print(repr(1.0 / (1.4 * 2.0**-0.4 * 5.0**0.4) ** 0.5))') # 1 / c(5) of POLY
         works=()
         for n in 12 24 48; do
            echo "-- eb3: the acoustic wave on rho0 = 5, N = $n"
            make_case "eb3-wave-$n" --case wave "${POLY[@]}" --rho0 5.0 --ni "$n" --cfl 0.1 \
               --rk runge-kutta-ssp-54 --time-max "$period"
            if run "$W" 2; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
         done
         if [[ ${#works[@]} -eq 3 ]]; then "$VENV_PY" "$CASE_DIR/eos_oracle.py" wave "${works[@]}" || FAILED=1 ; fi ;;
      eb4)
         eb4() { local name=$1 msg=$2 ; shift 2 ; refuse "$name" "$msg" --eos barotropic "$@" ; }
         eb4 decreasing "\[eos\]\.\(densities\) must be strictly increasing" --densities "8.0 2.0" \
            --exponents "1.0 1.4 1.6"
         eb4 count "\[eos\]\.\(exponents\) has 2 values, one per segment is needed: 3" --densities "2.0 8.0" \
            --exponents "1.0 1.4"
         eb4 below-one "\[eos\]\.\(exponents\) must all be >= 1" --densities "2.0 8.0" --exponents "0.9 1.4 1.6"
         eb4 segments "\[eos\] has 9 segments, at most 8 are accepted" --densities "1 2 3 4 5 6 7 8" \
            --exponents "1 1 1 1 1 1 1 1 1"
         eb4 reference-density "\[eos\]\.\(reference_density\) must be positive" --rho-ref -1.0
         eb4 gas "gamma, cp, cv are not used with \[eos\]\.\(type\) = barotropic" --keep-gas
         eb4 reference "\[reference\] is refused with \[eos\]\.\(type\)=barotropic" --set "reference:length=1.0"
         eb4 no-speed "failed to load \[eos\]\.\(sound_speed\) \(needed by type = barotropic\)" \
            --unset "eos:sound_speed" ;;
      eb5)
         echo "-- eb5: the one-segment table against the isothermal gas, c0 = 1.3"
         make_case eb5-iso --case tube --c0 1.3 --time-max 0.15 ; iso="$W"
         make_case eb5-baro --case tube --c0 1.3 --time-max 0.15 --eos barotropic --exponents "1.0" ; baro="$W"
         if run "$iso" && run "$baro"; then
            "$VENV_PY" "$CASE_DIR/eos_oracle.py" agree "$iso" "$baro" || FAILED=1
         else
            echo "   a run failed: see $iso/log.txt, $baro/log.txt  FAIL" ; FAILED=1
         fi ;;
      eb6)
         works=()
         for n in 200 400 800 1600; do
            echo "-- eb6: the default table (Mathew & Federrath 2021), N = $n"
            make_case "eb6-mf2021-$n" --case tube --eos barotropic --c0 2.0e4 --left "1.0e-7 0.0 0.0" \
               --right "1.0e-8 0.0 0.0" --time-max 3.0e-7 --ni "$n" --level 1 --null-y \
               --set bc_x_min:type=periodic --set bc_x_max:type=periodic --set weno:weights=si
            if run "$W" 2; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
         done
         if [[ ${#works[@]} -eq 4 ]]; then
            "$VENV_PY" "$CASE_DIR/eos_oracle.py" selfconv "${works[@]}" || FAILED=1
         fi ;;
      eb7)
         for solver in llf hll hllc; do
            echo "-- eb7: the barotropic collision, weno-riemann $solver"
            make_case "eb7-$solver" --case collide "${POLY[@]}" --left "1.0 3.4" --right "1.5 -2.6" \
               --time-max 0.3 --scheme weno-riemann --solver "$solver"
            run_oracle "$W" 1 riemann "$W"
         done ;;
      eb8)
         agree eb8 eb1-tube eb2-collide eb3-wave-48 eb7-hllc ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "EO FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "EO PASSED ($TAG, legs: $LEGS)"
