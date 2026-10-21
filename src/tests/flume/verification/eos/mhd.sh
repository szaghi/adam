#!/usr/bin/env bash
# FLUME verification EO, P5: the isothermal and barotropic MHD models (issue #100, M13 P5a, CPU).
#
# Why: the barotropic gases of P1-P4 gain a magnetic field: ideal MHD without an energy equation (D-M13-1), the six
# waves of Kim et al. (1999) for the characteristic projection and the isothermal HLLD of Mignone (2007), with
# a^2 = dP/drho (the isothermal eigensystem is the barotropic one). mhd_oracle.py checks the two papers first (selftest:
# L R = I and R Lambda L = dF/dU, HLLD consistency and exactness on a rotational discontinuity, its hydro limit, the
# exact CPAW); make_eos_mhd.py writes the cases. The legs:
#   em1  linear waves (one wavelength, amplitude 1e-7, the background of Stone et al. 2008, B = (1, sqrt(2), 1/2)):
#        fast, Alfven, slow, isothermal (c0 = 1: c_f = 2, c_a = 1, c_s = 1/2) along x (16, 32 cells) and inclined
#        (tan a = 2, 32, 64), barotropic along x on rho0 = 5 (the Gamma = 1.4 segment, c^2 = dP/drho = 2.0198), and
#        the fast wave with GLM: the order of the finest pair >= 4.7 (WENO-5; a wrong eigenvector leaves an O(amplitude)
#        error that does not converge);
#   em2  the em1 isothermal waves along x on weno-riemann (HLLD, characteristic, 6th-order correction): order >= 4.7;
#   em3  the circularly polarised Alfven wave (Toth 2000, B_perp = 0.1, an exact nonlinear solution), isothermal,
#        inclined, right and left, 32 and 64 cells: order >= 4.7, the two polarisations bitwise equal;
#   em4  Mignone's tubes 1-3 (Table 1, a = 1, weno-riemann HLLD): 200, 400, 800 cells against 3200, the L1 distance of
#        rho shrinking at an order >= 0.5 over the ladder, rho > 0, rE exactly 0; tube 3 (B_x = 0, two Mach-5 streams)
#        the degenerate HLLD (S*_L = S*_R = u*), also with the positivity limiter (positivity_limiter = cell, the
#        kinetic chord): rho > 0 and rE 0 asserted, its distance to the reference reported (the limiter acts at the
#        collision and leaves a start-up error in the zero-speed tangential mode: the doc page);
#   em5  rotation invariance: tube 2 along x, y and z (the states permuted cyclically), weno and weno-riemann: the
#        last checkpoints bitwise equal in the x frame (issue #41, MV-4: mhd_sum3 and the cyclic tangents);
#   em6  HLLD resolves what HLL smears: tube 2 (rotational discontinuities) at 400 cells, the L1 distance of B_y to the
#        3200-cell reference smaller with HLLD than with HLL;
#   em7  the isothermal limit: the one-segment table exponents = 1 with c0 = 1.3 against the isothermal gas with
#        c0 = 1.3 on tube 2: agreement within 1e-10 (the barotropic path takes P and c from the table);
#   em8  the contract: EGLM, the [mhd] floors, the Orszag-Tang vortex (a gamma-law density) and the entropy wave are
#        refused with the barotropic gases;
#   em9  (FNL executable only) the inert MHD models refused on the FNL backend until P5b;
#   em10 gravity on the magnetised gases: the uniform box of em1 (amplitude 0, B = (1, sqrt(2), 1/2)) under g_x = 0.5,
#        isothermal and barotropic (rho0 = 5), with GLM (c_h = 4: the gas accelerates to 0.2), SSP-33: rho uniform,
#        rho u = rho g t (eos_oracle.py freefall), rE exactly 0 (gravity moves the momentum only). SSP-33 as in eg1:
#        the library's SSP-54 weights sum to 1 - 8.8e-11 (14 printed digits), which a constant force integrates as a
#        relative drift of that size (issue #104).
#
# Usage: ./mhd.sh [--legs LIST]   (legs em1..em10; default em1..em8, em10)
#
# FLUME_EXE overrides the executable under test (default exe/adam_flume_cpu); the caller owns the environment. EM_NP sets
# the ranks (default 2).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
ORACLE="$CASE_DIR/mhd_oracle.py"
LW="$VERIF_DIR/mhd/linear-wave/linear_wave_oracle.py"
CPAW="$VERIF_DIR/mhd/cpaw/cpaw_oracle.py"
source "$VERIF_DIR/legs.sh"
legs_declare "em1 em2 em3 em4 em5 em6 em7 em8 em9 em10" "em1 em2 em3 em4 em5 em6 em7 em8 em10"
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      *)      echo "mhd.sh: unknown argument '$1' (accepted: --legs LIST)" >&2 ; exit 2 ;;
   esac
done
legs_resolve mhd.sh
if [[ ! -x "$EXE" ]]; then echo "mhd.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2 ; exit 2 ; fi
EXE="$(realpath "$EXE")"
TAG="$(basename "$EXE")"
FAILED=0
NP="${EM_NP:-2}"

run() { # run <work dir>: run the input of the work dir; 0 on success (an error stop exits 0 under mpirun: the log)
   (cd "$1" && mpirun -np "$NP" "$EXE" input.ini < /dev/null > log.txt 2>&1) && ! grep -aq "error stop" "$1/log.txt"
}
make_case() { # make_case <name> <make_eos_mhd.py case and options>: write one case, its dir in W
   local name=$1 ; shift
   W="$CASE_DIR/work-$TAG-$name"
   rm -rf "$W" ; mkdir -p "$W"
   "$VENV_PY" -I "$CASE_DIR/make_eos_mhd.py" "$W/input.ini" "$@"
}
ladder() { # ladder <label> <cells...> -- <make_eos_mhd.py case and options>: the runs, then the linear-wave order
   local label=$1 n ; shift
   local cells=() works=()
   while [[ $1 != "--" ]]; do cells+=("$1") ; shift ; done
   shift
   for n in "${cells[@]}"; do
      echo "-- $label, $n cells"
      make_case "$label-$n" "$@" --cells "$n"
      if run "$W"; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
   done
   if [[ ${#works[@]} -eq ${#cells[@]} ]]; then "$VENV_PY" "$LW" "${works[@]}" --order-min 4.7 || FAILED=1 ; fi
}
expect_refused() { # expect_refused <name> <message regex> <make_eos_mhd.py case and options>: one contract case
   local name=$1 msg=$2 ; shift 2
   echo "-- $name"
   make_case "$name" "$@"
   (cd "$W" && mpirun -np 1 "$EXE" input.ini < /dev/null > log.txt 2>&1) || true
   if grep -aqE "$msg" "$W/log.txt"; then
      echo "   refused as expected"
   else
      echo "   no refusal matching ($msg): see $W/log.txt  FAIL" ; FAILED=1
   fi
}
echo ">> EO P5: the barotropic MHD models ($TAG)"
for leg in $LEGS; do
   case $leg in
      em1)
         for wv in fast alfven slow; do
            ladder "em1-iso-1d-$wv" 16 32 -- wave --wave "$wv"
            ladder "em1-iso-2d-$wv" 32 64 -- wave --wave "$wv" --geometry 2d
            ladder "em1-baro-1d-$wv" 16 32 -- wave --wave "$wv" --eos barotropic --rho0 5
         done
         ladder em1-iso-glm-1d-fast 16 32 -- wave --wave fast --dc glm ;;
      em2)
         for wv in fast alfven slow; do
            ladder "em2-iso-hlld-1d-$wv" 16 32 -- wave --wave "$wv" --scheme weno-riemann --solver hlld
         done ;;
      em3)
         rights=() ; lefts=()
         for n in 32 64; do
            for pol in right left; do
               echo "-- em3: CPAW, $pol, $n cells"
               make_case "em3-cpaw-$pol-$n" cpaw --geometry 2d --polarisation "$pol" --cells "$n"
               if run "$W"; then
                  if [[ $pol == right ]]; then rights+=("$W") ; else lefts+=("$W") ; fi
               else
                  echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1
               fi
            done
         done
         if [[ ${#rights[@]} -eq 2 && ${#lefts[@]} -eq 2 ]]; then
            "$VENV_PY" "$CPAW" --right "${rights[@]}" --left "${lefts[@]}" --order-min 4.7 --lr-tol 0 || FAILED=1
         fi ;;
      em4)
         for t in 1 2 3; do
            works=()
            for n in 200 400 800 3200; do
               echo "-- em4: Mignone tube $t, HLLD, $n cells"
               make_case "em4-tube$t-$n" tube --tube "$t" --scheme weno-riemann --solver hlld --cells "$n"
               if run "$W"; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
            done
            if [[ ${#works[@]} -eq 4 ]]; then
               "$VENV_PY" "$ORACLE" selfconv "${works[@]:0:3}" --ref "${works[3]}" || FAILED=1
            fi
         done
         echo "-- em4: Mignone tube 3 with the positivity limiter"
         make_case em4-tube3-limited tube --tube 3 --scheme weno-riemann --solver hlld --cells 400 \
            --set numerics:positivity_limiter=cell
         if run "$W"; then
            "$VENV_PY" "$ORACLE" selfconv "$W" --ref "$CASE_DIR/work-$TAG-em4-tube3-3200" --order-min -99 || FAILED=1
         else
            echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1
         fi ;;
      em5)
         for sch in weno weno-riemann; do
            works=()
            for ax in x y z; do
               echo "-- em5: tube 2 along $ax, $sch"
               make_case "em5-$sch-$ax" tube --tube 2 --axis "$ax" --scheme "$sch" --cells 200
               if run "$W"; then works+=("$W") ; else echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1 ; fi
            done
            if [[ ${#works[@]} -eq 3 ]]; then "$VENV_PY" "$ORACLE" rotation "${works[@]}" || FAILED=1 ; fi
         done ;;
      em6)
         echo "-- em6: tube 2, HLL against HLLD"
         make_case em6-hlld tube --tube 2 --scheme weno-riemann --solver hlld --cells 400 ; hlld="$W"
         make_case em6-hll tube --tube 2 --scheme weno-riemann --solver hll --cells 400 ; hll="$W"
         make_case em6-ref tube --tube 2 --scheme weno-riemann --solver hlld --cells 3200 ; ref="$W"
         if run "$hlld" && run "$hll" && run "$ref"; then
            "$VENV_PY" "$ORACLE" hlld-gain "$hlld" "$hll" --ref "$ref" || FAILED=1
         else
            echo "   a run failed  FAIL" ; FAILED=1
         fi ;;
      em7)
         echo "-- em7: the one-segment table against the isothermal gas, c0 = 1.3"
         make_case em7-iso tube --tube 2 --c0 1.3 ; iso="$W"
         make_case em7-baro tube --tube 2 --c0 1.3 --eos barotropic --one-segment ; baro="$W"
         if run "$iso" && run "$baro"; then
            "$VENV_PY" "$ORACLE" agree "$iso" "$baro" --tol 1e-10 || FAILED=1
         else
            echo "   a run failed  FAIL" ; FAILED=1
         fi ;;
      em8)
         expect_refused em8-eglm "divergence_control\)=eglm is refused with \[eos\]" tube --cells 40 --dc glm \
            --set mhd:divergence_control=eglm
         expect_refused em8-floors "rho_floor, p_floor are refused with \[eos\]" tube --cells 40 --set mhd:p_floor=1e-6
         expect_refused em8-orszag-tang "orszag-tang is not defined for \[eos\]" tube --cells 40 \
            --set initial_conditions:type=orszag-tang
         expect_refused em8-entropy "no entropy mode" wave --cells 16 --set initial_conditions:wave=entropy ;;
      em9)
         expect_refused em9-fnl "runs on the CPU backend only" tube --cells 40 ;;
      em10)
         for gas in isothermal barotropic; do
            echo "-- em10: free fall, $gas, GLM"
            make_case "em10-$gas" wave --cells 48 --amplitude 0.0 --eos "$gas" --dc glm --time-max 0.4 --it-save 10 \
               --rho0 "$([[ $gas == barotropic ]] && echo 5 || echo 1)" --set gravity:g_x=0.5 --set mhd:glm_ch=4.0 \
               --set runge_kutta:scheme=runge-kutta-ssp-33
            if run "$W"; then
               "$VENV_PY" "$CASE_DIR/eos_oracle.py" freefall "$W" || FAILED=1
            else
               echo "   the run failed: see $W/log.txt  FAIL" ; FAILED=1
            fi
         done ;;
   esac
done
if [[ $FAILED -ne 0 ]]; then echo "EO P5 FAILED ($TAG, legs: $LEGS)" ; exit 1 ; fi
echo "EO P5 PASSED ($TAG, legs: $LEGS)"
