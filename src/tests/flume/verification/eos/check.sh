#!/usr/bin/env bash
# FLUME verification EO: the isothermal equation of state, `[eos] type = isothermal` (issue #100, M13 P1, P2).
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
#        [sinks], [gravity], viscosity, positivity_limiter = cell, [reference], immersed solids, the no-slip and
#        hydrostatic walls, an IC whose density needs gamma), and so are a malformed [eos] (mhd-ideal or passive
#        scalars with it, the gas keys cp/cv, a missing or non-positive sound_speed, an unknown type);
#   eo5  the inert slot and the mass: a wave of amplitude 0.3 on the periodic box, which breaks into a shock at
#        t ~ 1 / (c0 A k) = 0.53 (N = 48, two ranks, to t = 1): rE exactly 0 in every checkpoint, the total mass
#        constant to 1e-13;
#   eo6  the Riemann solvers: eo1 with weno-riemann and llf, hll, hllc (characteristic reconstruction, the 6th-order
#        flux correction), each against the exact solution as in eo1;
#   eo7  (P2) CPU and FNL agree (run eo1, eo2, eo3 and eo6 with each executable first): on the last checkpoint of the
#        tube, the collision, the N = 48 wave and the hllc tube, every field within 1e-10 of its scale (the floor of
#        AV av8): the density its largest magnitude, a momentum at least max(rho) c0 (its flux carries the pressure
#        c0^2 rho, so a 1e-6 wave's momentum has the round-off of an O(1) flux); rE exactly 0 on both.
#
# Usage: ./check.sh [--legs LIST]   (legs eo1 eo2 eo3 eo4 eo5 eo6 eo7, all; default all but eo7)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the matching environment.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
source "$VERIF_DIR/legs.sh"
legs_declare "eo1 eo2 eo3 eo4 eo5 eo6 eo7" "eo1 eo2 eo3 eo4 eo5 eo6"
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
echo ">> EO: the isothermal equation of state ($TAG)"
for leg in $LEGS; do
   case $leg in
      eo1)
         echo "-- eo1: the shock tube, WENO-5 characteristic"
         make_case eo1-tube --case tube
         run_oracle "$W" 1 riemann "$W" --case tube ;;
      eo2)
         echo "-- eo2: the Mach-5 collision"
         make_case eo2-collide --case collide --mach 5.0 --time-max 1.0
         run_oracle "$W" 1 riemann "$W" --case collide --mach 5.0 ;;
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
         refuse sinks "\[sinks\] is refused $ISO" --set "sinks:mass=1.0"
         refuse gravity "\[gravity\] is refused $ISO" --set "gravity:g_x=0.0"
         refuse viscosity "dissipative terms .* are refused with" --set "physics:viscosity=0.01"
         refuse limiter "\(positivity_limiter\)=cell is refused $ISO" --set "numerics:positivity_limiter=cell"
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
         refuse mhd "isothermal needs \[physics\]\.\(physical_model\)=euler" --set "physics:physical_model=mhd-ideal"
         refuse scalars "isothermal needs \[physics\]\.\(physical_model\)=euler without passive_scalars" \
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
            run_oracle "$W" 1 riemann "$W" --case tube
         done ;;
      eo7)
         for case_ in eo1-tube eo2-collide eo3-wave-48 eo6-hllc; do
            echo "-- eo7: CPU and FNL agree, $case_"
            cpu="$CASE_DIR/work-adam_flume_cpu-$case_" ; fnl="$CASE_DIR/work-adam_flume_fnl-$case_"
            if [[ -f "$cpu/log.txt" && -f "$fnl/log.txt" ]]; then
               "$VENV_PY" "$CASE_DIR/eos_oracle.py" agree "$cpu" "$fnl" || FAILED=1
            else
               echo "   missing a run: run eo1, eo2, eo3 and eo6 with both executables first" ; FAILED=1
            fi
         done ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "EO FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "EO PASSED ($TAG, legs: $LEGS)"
