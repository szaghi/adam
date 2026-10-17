#!/usr/bin/env bash
# FLUME verification GV (issue #73, M7): stratified atmospheres under constant gravity.
#
# M7-G0 (the inputs, before the gravity source of G1): every case of make_atmosphere.py runs one negligible step
# (CFL = 1e-30) with the ghosts poisoned, and atmosphere_oracle.py profile checks the checkpoint:
#   - the interior against the analytic atmosphere (isothermal or polytropic, region 1 at the origin);
#   - the ghosts beyond the `hydrostatic` faces against the boundary rule (the isothermal continuation from the first
#     interior cell along the normal; for the isothermal atmosphere the analytic profile itself, edges and corners too),
#     the velocity zero, the field and psi even;
#   - the periodic ghosts against the analytic profile at the wrapped position.
# Cases: iso2d, poly2d, tilted2d (every face hydrostatic, a refined block), tilted3d, iso2d-mhd (EGLM, a field oblique
# to gravity), froude2d (the number form of [gravity]), on each of --np.
# Negative control: iso2d run with extrapolation on the y faces and judged as if they were hydrostatic must FAIL (the
# oracle tells the boundary rule from a copy).
# NV-5 (--reference): ../scaling/reference.sh --inputs on iso2d, iso2d-mhd and froude2d: the [gravity] components
# convert as an acceleration (u0^2 / L0), froude does not, and the dimensional run equals the base one bit for bit.
#
# M7-G1 (the force; np 2, skipped with --no-force):
#   GV-1 freefall2d: rho u = rho g t and E = E_0 + rho g^2 t^2 / 2 to round-off (atmosphere_oracle.py freefall);
#   GV-3 box2d, box2d-mhd: a closed box, mass and the total energy W conserved to round-off (atmosphere_oracle.py
#        energy; W = sum (E + rho Phi) V plus the O(dx^2) wall-row term of the fourth-order energy source, M7-G2a;
#        the plain sum E + rho Phi is reported);
#   limiter box2d-lim (p = 0.05, the limiter acting) and box2d-stiff (p = 1e-3, g dt ~ 4.5 c): the runs complete,
#        the mass is conserved, E + rho Phi is reported (exact where the limiter does not act, option B);
#   GV-2 iso2d run with the force to t = 0.5 at N = 32, 64, 128: the spurious max |u| / c reported, with its order
#        (the unbalanced baseline of M7-G2);
#   GV-6 (--scaling): box2d with the scale-invariant weights against its input rescaled by powers of two
#        (../scaling/scaling.py rescale --j 2 --k -1 --m -2), bitwise (scaling.py compare);
#   GV-7 (--agree TAG): box2d and box2d-mhd against the runs of the executable TAG (run that backend first, --keep).
#
# M7-G2a (the order of the energy source; with the force):
#   GV-9 wave2d, a density wave falling freely in a periodic box, at N = 32, 64, 128 with the CFL scaled as N^(-2/3)
#        (the RK3 error below the spatial one), on `weno` and on `weno-riemann` (HLLC, sixth-order correction): the
#        orders of rho and p of the finest pair at least 3.5 (atmosphere_oracle.py wave). The face form of M7-G1 alone
#        gives 2 (the energy source is then second order); its fourth-order correction restores the scheme's order.
#
# M7-G2b (the well-balanced mode, [gravity] equilibrium = atmosphere; with the force, np 2, N = 32):
#   GV-2b strat2d (isothermal, ten scale heights) and polystrat2d (polytropic, G = 1.2) at rest to t = 0.5: balanced
#         max |u| / c below 1e-12 (asserted); the unbalanced strat2d reported;
#   GV-5  strat2d-mhdy, strat2d-mhdx (EGLM, a uniform field along and across g) to t = 0.2, balanced, on weno and on
#         weno-riemann (HLLD): below 1e-12;
#   GV-8  strat2d-amr (a refined box, 2:1 seams across the stratification) to t = 0.5: balanced below 1e-12 (the
#         equilibrium's seam mismatch is subtracted from the reflux), the unbalanced one reported;
#   GV-10 strat2d with a pressure pulse of 1e-6 and of 1e-4 to t = 0.1: the balanced responses per unit amplitude
#         agree to 1e-2 (asserted; the unbalanced ones are reported: their stationary truncation flow does not scale).
#   GV-7 (--agree) also compares the balanced 1e-4 pulse and the balanced strat2d-amr runs (momentum scaled by at
#         least max rho c: a balanced momentum is the perturbation, round-off at rest).
#
# M7-G3 (Rayleigh-Taylor, the `rayleigh-taylor` initial condition: rt2d, two isothermal layers 2 over 1 under
# g = (0, -1), interface width 0.05, on [0, 0.5] x [0.25, 0.75]; with the force, np 2):
#   GV-4a rt2d unseeded, one negligible step: the layers and their hydrostatic ghosts against the closed form of
#         rt_oracle.py (1e-13);
#   GV-4b rt2d unseeded to t = 1: unbalanced (reported), balanced (max |u| / c below 1e-12);
#   GV-4c rt2d seeded (1e-5, solenoidal) at N = 32 and 64 to t = 4: the growth rate of the mode, the slope of
#         ln a over t in [2, 4], against the compressible linear rate of rt_oracle.py (sigma = 1.5233) to 1e-3 (the
#         fit window's share; the two resolutions agree to 1.2e-4);
#   GV-7 (--agree) also compares the seeded rt2d at N = 32.
#
# Usage: ./check.sh [--np "1 2"] [--legs LIST] [--cases "iso2d ..."] [--reference] [--no-force] [--scaling]
#        [--agree TAG] [--keep]
# Legs, in run order: the G0 cases iso2d poly2d tilted2d tilted3d iso2d-mhd froude2d (on each --np), negative (the
# negative control), the force legs gv1 gv3 limiter gv2 gv9 gv2b gv5 gv8 gv10 gv4a gv4b gv4c, then gv6 and nv5.
# Default: every leg but gv6 and nv5; --no-force drops the force legs from it (selecting one with it is an error),
# --scaling adds gv6 and --reference nv5. --cases "..." selects G0 cases and keeps the other default legs, as
# before --legs. GV-7 is not a leg: --agree TAG runs it and needs gv3, gv8, gv10 and gv4c selected.
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
source "$VERIF_DIR/legs.sh"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NPS="1 2"
CASES="iso2d poly2d tilted2d tilted3d iso2d-mhd froude2d"
FORCE_LEGS="gv1 gv3 limiter gv2 gv9 gv2b gv5 gv8 gv10 gv4a gv4b gv4c"
legs_declare "$CASES negative $FORCE_LEGS gv6 nv5" "$CASES negative $FORCE_LEGS"
CASES_SEL=""
REFERENCE=0
KEEP=0
FORCE=1
SCALING=0
AGREE=""
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)        NPS="$2" ; shift 2 ;;
      --legs)      legs_select "$2" ; shift 2 ;;
      --cases)     CASES_SEL+=" ${2//,/ }" ; shift 2 ;;
      --reference) REFERENCE=1 ; shift ;;
      --no-force)  FORCE=0 ; shift ;;
      --scaling)   SCALING=1 ; shift ;;
      --agree)     AGREE="$2" ; shift 2 ;;
      --keep)      KEEP=1 ; shift ;;
      *) echo "check.sh: unknown argument '$1' (accepted: --np, --legs, --cases, --reference, --no-force, --scaling, --agree, --keep)" \
                >&2 ; exit 2 ;;
   esac
done
if [[ $FORCE -eq 0 ]]; then
   legs_declare "$CASES negative $FORCE_LEGS gv6 nv5" "$CASES negative"
fi
if [[ -n $CASES_SEL ]]; then # --cases: the G0 cases named, with the other legs of the default
   for case in $CASES_SEL; do
      [[ " $CASES " == *" $case "* ]] || { echo "check.sh: unknown case '$case' (accepted: $CASES)" >&2 ; exit 2 ; }
   done
   legs_select "$CASES_SEL negative"
   if [[ $FORCE -eq 1 ]]; then legs_select "$FORCE_LEGS" ; fi
fi
if [[ -z $LEGS_SELECTED && ( $SCALING -eq 1 || $REFERENCE -eq 1 ) ]]; then
   legs_select "$LEGS_DEFAULT"
fi
if [[ $SCALING -eq 1 ]]; then legs_select gv6 ; fi
if [[ $REFERENCE -eq 1 ]]; then legs_select nv5 ; fi
legs_resolve check.sh
if [[ $FORCE -eq 0 ]]; then
   for leg in $FORCE_LEGS; do
      if has_leg "$leg"; then echo "check.sh: --no-force excludes the force leg '$leg'" >&2 ; exit 2 ; fi
   done
fi
if [[ -n "$AGREE" ]]; then # GV-7 compares the runs of these legs
   for leg in gv3 gv8 gv10 gv4c; do
      if ! has_leg "$leg"; then echo "check.sh: --agree (GV-7) needs the leg '$leg'" >&2 ; exit 2 ; fi
   done
fi
TAG="$(basename "$EXE")"
FAILED=0

run_case() { # run_case <work> <np>: run the input of <work>, report a failure
   if ! (cd "$1" && mpirun -np "$2" "$EXE" input.ini < /dev/null > log.txt 2>&1); then
      echo "   run failed: see $1/log.txt" ; FAILED=1 ; return 1
   fi
}

echo ">> GV: stratified atmospheres ($TAG)"
for np in $NPS; do
   for case in $CASES; do
      has_leg "$case" || continue
      work="$CASE_DIR/work-$TAG-np$np-$case"
      rm -rf "${work:?}" ; mkdir -p "$work"
      python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" "$case"
      echo "-- $case, np $np"
      run_case "$work" "$np" || continue
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" profile "$work" || FAILED=1
      [[ $KEEP -eq 1 ]] || find "$work" -name '*.h5' -delete
   done
done

if has_leg negative; then
   echo "-- negative control: iso2d with extrapolation on the y faces, judged as hydrostatic (must FAIL)"
   work="$CASE_DIR/work-$TAG-negative"
   rm -rf "${work:?}" ; mkdir -p "$work"
   python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" iso2d
   sed -i '/^\[bc_y_m[ai][nx]\]/{n;s/^type = hydrostatic/type = extrapolation/}' "$work/input.ini"
   if run_case "$work" 1; then
      sed -i '/^\[bc_y_m[ai][nx]\]/{n;s/^type = extrapolation/type = hydrostatic/}' "$work/input.ini"
      if "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" profile "$work" > "$work/oracle.txt"; then
         echo "   the oracle passed copied ghosts: FAIL" ; FAILED=1
      else
         echo "   the oracle rejects copied ghosts" \
              "($(grep -a 'hydrostatic ghosts rho' "$work/oracle.txt" | sed 's/^ *//')): PASS"
      fi
   fi
   [[ $KEEP -eq 1 ]] || find "$work" -name '*.h5' -delete
fi

if [[ $FORCE -eq 1 ]]; then
   for case in freefall2d box2d box2d-mhd box2d-lim box2d-stiff; do
      case $case in freefall2d) leg=gv1 ;; box2d|box2d-mhd) leg=gv3 ;; *) leg=limiter ;; esac
      has_leg "$leg" || continue
      work="$CASE_DIR/work-$TAG-np2-$case"
      rm -rf "${work:?}" ; mkdir -p "$work"
      python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" "$case"
      echo "-- $case, np 2"
      run_case "$work" 2 || continue
      echo "   limiter: $(grep -a '^\[mpi-00000\]positivity limiter:' "$work/log.txt" | awk '{f+=$3; b+=$6} END{print f+0" faces limited, "b+0" inadmissible backbones"}')"
      case $case in
         freefall2d)          "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" freefall "$work" || FAILED=1 ;;
         box2d|box2d-mhd)     "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" energy "$work" || FAILED=1 ;;
         *)                   "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" energy "$work" --report-total || FAILED=1 ;;
      esac
      [[ $KEEP -eq 1 || -n "$AGREE" ]] || find "$work" -name '*.h5' -delete
   done
   if has_leg gv2; then
      echo "-- GV-2: iso2d at rest under the force to t = 0.5, N = 32, 64, 128 (reported)"
      rest=()
      for cells in 8 16 32; do
         work="$CASE_DIR/work-$TAG-np2-rest-$cells"
         rm -rf "${work:?}" ; mkdir -p "$work"
         python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" iso2d --cells "$cells" --cfl 0.4 --time-max 0.5
         sed -i 's/^it_save = .*/it_save = 100000/' "$work/input.ini"
         run_case "$work" 2 && rest+=("$work")
      done
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "${rest[@]}"
      [[ $KEEP -eq 1 ]] || for w in "${rest[@]}"; do find "$w" -name '*.h5' -delete; done
   fi
   if has_leg gv9; then
      for scheme in weno riemann; do
         echo "-- GV-9: wave2d on $scheme, N = 32, 64, 128 (orders >= 3.5)"
         wave=()
         for cells in 8 16 32; do
            work="$CASE_DIR/work-$TAG-np2-wave-$scheme-$cells"
            rm -rf "${work:?}" ; mkdir -p "$work"
            extra=() ; [[ $scheme == riemann ]] && extra=(--riemann hllc)
            python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" wave2d --cells "$cells" --time-max 0.2 \
               --cfl "$(python3 -c "print(0.4 * (8 / $cells) ** (2 / 3))")" "${extra[@]}"
            sed -i 's/^it_save = .*/it_save = 100000/' "$work/input.ini"
            run_case "$work" 2 && wave+=("$work")
         done
         "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" wave "${wave[@]}" || FAILED=1
         [[ $KEEP -eq 1 ]] || for w in "${wave[@]}"; do find "$w" -name '*.h5' -delete; done
      done
   fi
   g2b() { # g2b <name> <case> [make_atmosphere options]: run one G2b case at np 2, echo its work directory
      local work="$CASE_DIR/work-$TAG-np2-$1" ; shift
      rm -rf "${work:?}" ; mkdir -p "$work"
      python3 -I "$CASE_DIR/make_atmosphere.py" "$work/input.ini" "$@"
      run_case "$work" 2 > /dev/null && echo "$work"
   }
   if has_leg gv2b; then
      echo "-- GV-2b: strat2d unbalanced (reported), strat2d and polystrat2d balanced (max |u| / c < 1e-12)"
      unb=$(g2b strat2d-unb strat2d) || FAILED=1
      [[ -n "$unb" ]] && "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "$unb"
      bal=()
      for case in strat2d polystrat2d; do
         w=$(g2b "$case-bal" "$case" --balanced) && bal+=("$w") || FAILED=1
      done
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "${bal[@]}" --max 1e-12 || FAILED=1
   fi
   if has_leg gv5; then
      echo "-- GV-5: strat2d-mhdy, strat2d-mhdx balanced, weno and weno-riemann (HLLD) (max |u| / c < 1e-12)"
      mhd=()
      for case in strat2d-mhdy strat2d-mhdx; do
         w=$(g2b "$case-bal" "$case" --balanced --time-max 0.2) && mhd+=("$w") || FAILED=1
         w=$(g2b "$case-bal-hlld" "$case" --balanced --time-max 0.2 --riemann hlld) && mhd+=("$w") || FAILED=1
      done
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "${mhd[@]}" --max 1e-12 || FAILED=1
   fi
   if has_leg gv8; then
      echo "-- GV-8: strat2d-amr unbalanced (reported), balanced (max |u| / c < 1e-12)"
      w=$(g2b amr-unb strat2d-amr) && "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "$w"
      w=$(g2b amr-bal strat2d-amr --balanced) && "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "$w" --max 1e-12 \
         || FAILED=1
   fi
   if has_leg gv10; then
      echo "-- GV-10: a pressure pulse of 1e-6 and 1e-4 on strat2d, unbalanced (reported) and balanced (linear to 1e-2)"
      for mode in unb bal; do
         opt=() ; [[ $mode == bal ]] && opt=(--balanced)
         a=$(g2b "lin-$mode-6" strat2d --pulse 1e-6 --time-max 0.1 "${opt[@]}")
         b=$(g2b "lin-$mode-4" strat2d --pulse 1e-4 --time-max 0.1 "${opt[@]}")
         if [[ $mode == bal ]]; then
            "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" linear "$a" "$b" || FAILED=1
         else
            "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" linear "$a" "$b" | \
               sed 's/: FAIL (tol/: reported, not asserted (tol/' || true
         fi
      done
   fi
   if has_leg gv2b || has_leg gv5 || has_leg gv8 || has_leg gv10; then
      [[ $KEEP -eq 1 || -n "$AGREE" ]] || find "$CASE_DIR"/work-"$TAG"-np2-*-bal* "$CASE_DIR"/work-"$TAG"-np2-*-unb* \
         "$CASE_DIR"/work-"$TAG"-np2-lin-* -name '*.h5' -delete 2> /dev/null || true # a glob of a leg not run
   fi
   if has_leg gv4a; then
      echo "-- GV-4a: rt2d unseeded, one negligible step: the layers and their hydrostatic ghosts (1e-13)"
      w=$(g2b rt-profile rt2d --seed 0 --steps 1 --cfl 1e-30) && \
         "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" profile "$w" \
         || FAILED=1
   fi
   if has_leg gv4b; then
      echo "-- GV-4b: rt2d unseeded to t = 1, unbalanced (reported) and balanced (max |u| / c < 1e-12)"
      w=$(g2b rt-unb rt2d --seed 0 --time-max 1.0) && "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "$w"
      w=$(g2b rt-bal rt2d --seed 0 --time-max 1.0 --balanced) && \
         "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" rest "$w" --max 1e-12 || FAILED=1
   fi
   if has_leg gv4c; then
      echo "-- GV-4c: rt2d growth rate at N = 32, 64 against the compressible linear theory (1e-3)"
      rt=()
      for cells in 8 16; do
         w=$(g2b "rt-$cells" rt2d --cells "$cells" --it-save $((5 * cells))) && rt+=("$w") || FAILED=1
      done
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" growth "${rt[@]}" --t0 2 --t1 4 --tol 1e-3 || FAILED=1
   fi
   if has_leg gv4a || has_leg gv4b || has_leg gv4c; then
      [[ $KEEP -eq 1 || -n "$AGREE" ]] || find "$CASE_DIR"/work-"$TAG"-np2-rt-* -name '*.h5' -delete
   fi
fi

if has_leg gv6; then
   echo "-- GV-6: box2d (scale-invariant weights) against its input rescaled by 2^2, 2^-1, 4^-2"
   base="$CASE_DIR/work-$TAG-np2-scaling-base" ; scaled="$CASE_DIR/work-$TAG-np2-scaling-scaled"
   rm -rf "${base:?}" "${scaled:?}" ; mkdir -p "$base" "$scaled"
   python3 -I "$CASE_DIR/make_atmosphere.py" "$base/input.ini" box2d --weights si --steps 100
   "$VENV_PY" "$CASE_DIR/../scaling/scaling.py" rescale "$base/input.ini" "$scaled/input.ini" --j 2 --k -1 --m -2
   if run_case "$base" 2 && run_case "$scaled" 2; then
      "$VENV_PY" "$CASE_DIR/../scaling/scaling.py" compare "$base" "$scaled" --j 2 --k -1 --m -2 --ngc 4 || FAILED=1
   fi
   [[ $KEEP -eq 1 ]] || find "$base" "$scaled" -name '*.h5' -delete
fi

if [[ -n "$AGREE" ]]; then
   echo "-- GV-7: box2d, box2d-mhd, lin-bal-4, amr-bal, rt-8 against $AGREE"
   for case in box2d box2d-mhd lin-bal-4 amr-bal rt-8; do
      "$VENV_PY" "$CASE_DIR/atmosphere_oracle.py" agree "$CASE_DIR/work-$TAG-np2-$case" "$CASE_DIR/work-$AGREE-np2-$case" \
         || FAILED=1
   done
fi

if has_leg nv5; then
   echo "-- NV-5: [reference] on iso2d, iso2d-mhd, froude2d"
   inputs=()
   mkdir -p "$CASE_DIR/work-$TAG-ref"
   for case in iso2d iso2d-mhd froude2d; do
      python3 -I "$CASE_DIR/make_atmosphere.py" "$CASE_DIR/work-$TAG-ref/$case.ini" "$case"
      inputs+=("$CASE_DIR/work-$TAG-ref/$case.ini")
   done
   (cd "$CASE_DIR/../scaling" && FLUME_EXE="$EXE" ./reference.sh --inputs "${inputs[*]}") || FAILED=1
fi

if [[ $FAILED -eq 0 ]]; then
   echo "GV PASSED ($TAG, legs: $LEGS)"
else
   echo "GV FAILED ($TAG, legs: $LEGS)"
fi
exit $FAILED
