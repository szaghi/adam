#!/usr/bin/env bash
# FLUME verification VV (issue #65, P2): the Navier-Stokes dissipative fluxes.
#
# Why: a dissipative flux can be consistent and still wrong in its order (a wrong correction term), its coefficient
# (a missing 4/3, a k without R) or its direction (a cross term read with the wrong stencil). Each leg has an exact
# solution and an asserted bound:
#   vv1  shear wave (waves.py, linearised NS mode): measured order of the finest pair >= ORDER4_MIN at
#        dissipative_order = 4 along x, along y and along the diagonal (cross derivatives), >= ORDER2_MIN at order 2;
#   vv2  viscous-thermal acoustic wave (mu and k together): the same order bound at order 4, on the odd part of the
#        +A / -A twin runs (an acoustic wave steepens at O(A^2), which would otherwise floor the error);
#   vv3  compressible Couette flow between an isothermal wall and a moving adiabatic one (couette.py);
#   vv4  Becker's viscous shock (profiles.py);
#   vv5  MHD (issue #65 P3): Ohmic decay of a transverse field (eta k^2), along x at orders 4 and 2 and along the
#        diagonal (the curl form's cross terms), and the visco-resistive Alfven wave (B_n = 1, mu and eta together),
#        against the exact linearised (v_t, b_t) system, and vv2's acoustic wave on mhd-none with B = 0 (the
#        viscous and heat fluxes on MHD, read through its own auxiliary temperature): the same order bounds as vv1.
#        The diagonal runs without divergence control: GLM's upwinding (c_h = 3) damps the static diagonal field at
#        fifth order (8.5e-12 at N = 48, measured with eta = 0), opposite in sign to the Ohmic error and as large, so
#        the sum cancels near N = 48 and a 24/48/96 ladder is not asymptotic (+1.72; +3.32 at 96/192). That is the
#        ideal scheme's error; vv6 covers GLM and EGLM with eta;
#   vv6  MHD with GLM and EGLM: the Ohmic energy budget (the field's energy loss reappears as heat, total energy
#        conserved) and div(B) not raised above the ideal twin's (waves.py budget).
#   vv7  (P4) conservation across 2:1 seams: a diagonal wave of amplitude 0.01 in the periodic box with its centre
#        refined 2:1 (quadtree), Euler (mu, k) and MHD-GLM (mu, eta, B_n = 1), t = 0.05: with reflux every integral
#        holds within MAX_DRIFT of int|q| (waves.py conserve); without reflux it drifts by at least MIN_DRIFT, and by
#        at least MIN_DRIFT differently from the ideal twin's, so the dissipative fluxes do cross the seams unmatched.
#        The forest: the same wave with outflow x faces split at x = 0.5 into two realms, mirror (1:1, quadtree) and
#        refined (2:1, octree; the refined coupling needs ratio 2 along every axis), against the single realm:
#        fields within FOREST_TOL (the sine is evaluated from each realm's block origins, one ulp apart at step 0).
#   vv8  (P4) accuracy across 2:1 seams: the diagonal shear (mu) and Ohmic (eta, mhd-none) waves on the refined
#        quadtree, order >= SEAM_ORDER_MIN: a 2:1 seam of point values caps the order at 2 (issue #21), so the
#        composite error is second order and its ratio to the uniform runs grows with N. Without reflux the Ohmic
#        ladder falls to first order (<= SEAM_NOREFLUX_MAX): the reflux carries the dissipative flux.
#
#   agree (P5, runs nothing) every run of the legs above on the CPU against the same run on FNL, cell by cell on the
#        last checkpoint, each difference over the physical scale of its variable (waves.py agree): <= AGREE_TOL.
#        Needs the work directories of both executables (run the legs with each first); not in the default set.
#
# Usage: ./check.sh [--np N] [--leg vv1|vv2|vv3|vv4|vv5|vv6|vv7|vv8|agree ...]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEGS=()
ORDER4_MIN="3.8"
ORDER2_MIN="1.9"
BECKER_ORDER_MIN="3.0"
MAX_DRIFT="1.0e-13"
MIN_DRIFT="1.0e-10"
FOREST_TOL="1.0e-10"
SEAM_ORDER_MIN="1.8"
SEAM_NOREFLUX_MAX="1.3"
AGREE_TOL="1.0e-10"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)  NP="$2" ; shift 2 ;;
      --leg) LEGS+=("$2") ; shift 2 ;;
      *)     echo "check.sh: unknown argument '$1' (accepted: --np N, --leg L)" >&2 ; exit 2 ;;
   esac
done
[[ ${#LEGS[@]} -eq 0 ]] && LEGS=(vv1 vv2 vv3 vv4 vv5 vv6 vv7 vv8)
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

run() { (cd "$1" && mpirun -np "$NP" "$EXE" "${2:-input.ini}" < /dev/null > log.txt 2>&1) ; } # mpirun reads stdin

ladder() { # ladder <tool> <name> <order min> "<resolutions>" <make options...>: run, judge the finest pair
   # LADDER_TWIN=1 also runs each resolution with the amplitude -A into <work>-neg (waves.py: the odd part is judged)
   local tool="$CASE_DIR/$1" name="$2" pmin="$3" ns=($4) works=() w v n
   shift 4
   echo "-- $name"
   for n in "${ns[@]}"; do
      w="$CASE_DIR/work-$TAG-$name-n$n"
      for v in "$w" ${LADDER_TWIN:+"$w-neg"}; do
         rm -rf "$v" ; mkdir -p "$v"
         if [[ $v == *-neg ]]; then
            "$VENV_PY" "$tool" make "$v/input.ini" --n "$n" "$@" --sign -1
         else
            "$VENV_PY" "$tool" make "$v/input.ini" --n "$n" "$@"
         fi
         if ! run "$v"; then echo "   run failed, see $v/log.txt" ; FAILED=1 ; return 0 ; fi
         rm -f "$v"/*-000000000-*.h5
      done
      works+=("$w")
   done
   # pipefail: a failing oracle fails the pipeline, which marks the leg failed instead of aborting the check (set -e)
   # LADDER_ORDER_MAX bounds the order from above instead (a negative control), the minimum then 0
   local bound=(--order-min "$pmin")
   [[ -n "${LADDER_ORDER_MAX:-}" ]] && bound=(--order-max "$LADDER_ORDER_MAX")
   if ! "$VENV_PY" "$tool" oracle "${works[@]}" "${bound[@]}" | sed 's/^/   /'; then FAILED=1 ; fi
}

conserve() { # conserve <name> <waves.py make options...>: reflux on, off and the ideal twin off; waves.py conserve
   local name="$1" w="$CASE_DIR/work-$TAG-$1" v
   shift
   echo "-- $name"
   for v in "$w" "$w-off" "$w-ideal"; do
      rm -rf "$v" ; mkdir -p "$v"
      case "$v" in
         *-off)   "$VENV_PY" "$CASE_DIR/waves.py" make "$v/input.ini" "$@" --reflux false ;;
         *-ideal) "$VENV_PY" "$CASE_DIR/waves.py" make "$v/input.ini" "$@" --reflux false --mu 0.0 --kappa 0.0 \
                     --eta 0.0 ;;
         *)       "$VENV_PY" "$CASE_DIR/waves.py" make "$v/input.ini" "$@" ;;
      esac
      if ! run "$v"; then echo "   run failed, see $v/log.txt" ; FAILED=1 ; return 0 ; fi
   done
   if ! "$VENV_PY" "$CASE_DIR/waves.py" conserve "$w" --leaky "$w-off" --ideal "$w-ideal" --max-drift "$MAX_DRIFT" \
        --min-drift "$MIN_DRIFT" | sed 's/^/   /'; then FAILED=1 ; fi
}

forest() { # forest <name> mirror|refined <waves.py make options...>: single realm vs its 2-realm split
   local name="$1" kind="$2" w="$CASE_DIR/work-$TAG-$1" flag=()
   shift 2
   [[ $kind == refined ]] && flag=(--refined)
   echo "-- $name"
   rm -rf "$w-single" "$w-forest" ; mkdir -p "$w-single" "$w-forest"
   "$VENV_PY" "$CASE_DIR/waves.py" make "$w-single/input.ini" --bc-x extrapolation "$@"
   "$VENV_PY" "$CASE_DIR/waves.py" split "$w-single/input.ini" "$w-forest" wave "${flag[@]}"
   if ! run "$w-single"; then echo "   run failed, see $w-single/log.txt" ; FAILED=1 ; return 0 ; fi
   if ! run "$w-forest" wave.ini; then echo "   run failed, see $w-forest/log.txt" ; FAILED=1 ; return 0 ; fi
   # --fields-only: the momentum integrals of a sine wave are zero, so their relative difference measures nothing
   if ! "$VENV_PY" "$CASE_DIR/../multirealm/multirealm_oracle.py" "$w-forest" "$w-single" --ngc 3 --tol "$FOREST_TOL" \
        --fields-only | sed 's/^/   /'; then FAILED=1 ; fi
}

budget() { # budget <name> <waves.py make options...>: the resistive run and its ideal twin, then the VV-6 oracle
   local name="$1" w="$CASE_DIR/work-$TAG-$1" v
   shift
   echo "-- $name"
   for v in "$w" "$w-ideal"; do
      rm -rf "$v" ; mkdir -p "$v"
      if [[ $v == *-ideal ]]; then
         "$VENV_PY" "$CASE_DIR/waves.py" make "$v/input.ini" "$@" --eta 0.0
      else
         "$VENV_PY" "$CASE_DIR/waves.py" make "$v/input.ini" "$@"
      fi
      if ! run "$v"; then echo "   run failed, see $v/log.txt" ; FAILED=1 ; return 0 ; fi
   done
   if ! "$VENV_PY" "$CASE_DIR/waves.py" budget "$w" "$w-ideal" | sed 's/^/   /'; then FAILED=1 ; fi
}

echo ">> VV: Navier-Stokes and Ohmic dissipative fluxes ($TAG)"
for leg in "${LEGS[@]}"; do
   case "$leg" in
      vv1)
         ladder waves.py shear-x-o4 "$ORDER4_MIN" "32 64 128" --mode shear --angle 0 --order 4
         ladder waves.py shear-x-o2 "$ORDER2_MIN" "32 64 128" --mode shear --angle 0 --order 2
         ladder waves.py shear-y-o4 "$ORDER4_MIN" "32 64 128" --mode shear --angle 90 --order 4
         ladder waves.py shear-xy-o4 "$ORDER4_MIN" "24 48 96" --mode shear --angle 45 --order 4 --time 0.1 ;;
      vv2)
         LADDER_TWIN=1 ladder waves.py acoustic-x-o4 "$ORDER4_MIN" "32 64 128" --mode acoustic --angle 0 --order 4 \
                --mu 0.01 --kappa 0.02 ;;
      vv3)
         ladder profiles.py couette-o4 "$ORDER2_MIN" "32 64 128" --case couette --order 4 ;;
      vv4)
         ladder profiles.py becker-m2-o4 "$BECKER_ORDER_MIN" "64 128 256" --case becker --mach 2 --order 4 --time 0.05
         ladder profiles.py becker-m3-o4 "$BECKER_ORDER_MIN" "64 128 256" --case becker --mach 3 --order 4 --time 0.05 \
                --mu 0.01 ;;
      vv5)
         ladder waves.py ohmic-x-o4 "$ORDER4_MIN" "32 64 128" --model mhd-none --mode magnetic --angle 0 --order 4 \
                --mu 0.0 --eta 0.01
         ladder waves.py ohmic-x-o2 "$ORDER2_MIN" "32 64 128" --model mhd-none --mode magnetic --angle 0 --order 2 \
                --mu 0.0 --eta 0.01
         ladder waves.py ohmic-xy-o4 "$ORDER4_MIN" "24 48 96" --model mhd-none --mode magnetic --angle 45 --order 4 \
                --mu 0.0 --eta 0.01 --time 0.1
         ladder waves.py alfven-x-o4 "$ORDER4_MIN" "32 64 128" --model mhd-glm --mode magnetic --angle 0 --order 4 \
                --mu 0.01 --eta 0.005 --b0 1.0
         LADDER_TWIN=1 ladder waves.py acoustic-mhd-x-o4 "$ORDER4_MIN" "32 64 128" --model mhd-none --mode acoustic \
                --angle 0 --order 4 --mu 0.01 --kappa 0.02 ;;
      vv6)
         # A = 1e-4: the budget compares the field energy (A^2) with internal-energy differences, so it needs the signal
         budget budget-glm --n 48 --model mhd-glm --mode magnetic --angle 45 --order 4 --mu 0.0 --eta 0.01 --time 0.2 \
                --amplitude 1.0e-4
         budget budget-eglm --n 48 --model mhd-eglm --mode magnetic --angle 45 --order 4 --mu 0.0 --eta 0.01 \
                --time 0.2 --amplitude 1.0e-4 ;;
      vv7)
         box=(--angle 45 --order 4 --time 0.05 --amplitude 0.01 --refine 0.25 0.25 0.75 0.75 --ratio 4)
         euler=(--mode acoustic --mu 0.01 --kappa 0.02)
         mhd=(--model mhd-glm --mode magnetic --b0 1.0 --mu 0.01 --eta 0.01)
         # quadtree only: with z null an octree solves the same problem (measured: identical drifts to the printed
         # digits, at 4 times the cells); the octree's 2:1 seams are exercised by the refined forest below
         conserve conserve-euler-quadtree --n 48 "${euler[@]}" "${box[@]}"
         conserve conserve-mhd-quadtree --n 48 "${mhd[@]}" "${box[@]}"
         wave=(--n 48 --angle 45 --order 4 --time 0.05 --amplitude 0.01)
         forest forest-euler-mirror mirror "${wave[@]}" "${euler[@]}"
         forest forest-mhd-mirror mirror "${wave[@]}" "${mhd[@]}"
         forest forest-euler-refined refined "${wave[@]}" "${euler[@]}" --refine 0.5 0.0 1.0 1.0 --ratio 8
         forest forest-mhd-refined refined "${wave[@]}" "${mhd[@]}" --refine 0.5 0.0 1.0 1.0 --ratio 8 ;;
      vv8)
         seam=(--angle 45 --order 4 --time 0.1 --refine 0.25 0.25 0.75 0.75 --ratio 4)
         ladder waves.py seam-shear-o4 "$SEAM_ORDER_MIN" "24 48 96" --mode shear --mu 0.01 "${seam[@]}"
         ladder waves.py seam-ohmic-o4 "$SEAM_ORDER_MIN" "24 48 96" --model mhd-none --mode magnetic --mu 0.0 \
                --eta 0.01 "${seam[@]}"
         LADDER_ORDER_MAX="$SEAM_NOREFLUX_MAX" ladder waves.py seam-ohmic-o4-noreflux 0 "24 48 96" --model mhd-none \
                --mode magnetic --mu 0.0 --eta 0.01 "${seam[@]}" --reflux false ;;
      agree)
         echo "-- agree: every work-adam_flume_cpu-np$NP-* run against its work-adam_flume_fnl-np$NP-* twin"
         pairs=0
         for c in "$CASE_DIR"/work-adam_flume_cpu-np"$NP"-*; do
            f="${c/adam_flume_cpu/adam_flume_fnl}"
            [[ -d $f ]] || continue
            pairs=$((pairs + 1))
            if ! "$VENV_PY" "$CASE_DIR/waves.py" agree "$c" "$f" --tol "$AGREE_TOL" | sed 's/^/   /'; then FAILED=1 ; fi
         done
         if [[ $pairs -eq 0 ]]; then echo "   no CPU/FNL pair of work directories: run the legs on both first" ; FAILED=1
         else echo "   $pairs pairs compared"; fi ;;
      *) echo "check.sh: unknown leg '$leg'" >&2 ; exit 2 ;;
   esac
done

if [[ $FAILED -eq 0 ]]; then
   echo "VV PASSED ($TAG)"
else
   echo "VV FAILED ($TAG)"
   exit 1
fi
