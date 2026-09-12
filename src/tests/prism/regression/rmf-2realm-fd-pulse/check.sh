#!/usr/bin/env bash
# rmf-2realm-fd-pulse divergence oracle (issue #31, inter-realm 1:1 mirror seam).
#
# THE CLAIM UNDER TEST: a 1:1 same-resolution INTER-realm mirror seam must be
# div-free to machine precision. It is numerically identical to a 1:1
# intra-block interface -- fill_seam_from_peer_forest is a pure peer-interior
# COPY, so the mimetic div_h(curl_h) = 0 identity holds across it. With J = 0
# (source-free gaussian pulse, zero coils) BOTH discrete divergences are
# conserved, so any drift above round-off is a defect in the seam machinery.
#
# This is NOT the #29 2:1-jump case. There the coarse and fine operators differ
# and an O(h^p) div(B) source is expected and accepted; here there is no
# resolution jump, so nothing may leak.
#
# ACCEPTANCE
#   1. beta cadence (as committed): max|div(D)| AND max|div(B)| <= DIV_TOL, on
#      BOTH realms. This is the anchor -- it is what #31 fixed.
#   2. both realms reach 100% and report no error/abort/NaN.
#   3. alpha cadence (negative control): flipping coupling_cadence to
#      end_of_step must push div(B) ABOVE DIV_TOL. Without this leg the oracle
#      could pass vacuously -- e.g. if the seam silently stopped being
#      exercised at all -- so the leg proves the measurement has teeth. It
#      asserts only THAT alpha leaks, deliberately not a pinned magnitude:
#      alpha's drift is not a quantity this suite has any reason to hold fixed.
#
#      Why alpha leaks (#31): fill_seam_from_peer_forest runs once per step
#      under alpha, so RK substages 2..N read UNFILLED seam ghosts out of the
#      stage buffer q_rk(:,...,k) and the curl evolves the seam skin from
#      garbage. Beta fills every substage. Alpha's end-of-step lag is a
#      Berger-Oliger AMR-subcycling convention, not a div-preserving one.
#
# Usage: ./check.sh [--build]
#
# PRISM_EXE: override the executable under test, e.g.
#   PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH plus the
# WSL UCX knobs of issue #12); --build always builds the CPU default.
#
# Runs at -np 1: the seam is then a same-rank peer copy, which is the mechanism
# under test. (-np 2 exercises the MPI seam path instead; the committed digest
# goldens of the sibling rmf-2realm* cases cover that.)
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"

DIV_TOL="1.0E-13"   # round-off ceiling: a 1:1 seam conserves both invariants

do_build=0
for arg in "$@"; do
   case "$arg" in
      --build) do_build=1 ;;
      *) echo "ERROR: unknown flag $arg (use --build)" >&2; exit 2 ;;
   esac
done
if [[ $do_build -eq 1 ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# divergence-history columns:
#   it  blocks_number  time  D_divergence  B_divergence  J_divergence
max_div_d() { grep '^+' "$1" | awk '{v=$4; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }
max_div_b() { grep '^+' "$1" | awk '{v=$5; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }

R1_HIST="rmf_2realm_fd_pulse_r1-divergence_history.dat"
R2_HIST="rmf_2realm_fd_pulse_r2-divergence_history.dat"

run_in() { # workdir, sed-expr applied to the forest manifest
   local wd="$1" sed_expr="$2"
   rm -rf "$wd" && mkdir -p "$wd"
   sed "$sed_expr" "$CASE_DIR/input.ini" > "$wd/input.ini"
   cp "$CASE_DIR/realm_1.ini" "$CASE_DIR/realm_2.ini" "$wd/"
   ( cd "$wd" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 )
   find "$wd" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' -o -name '*.tnd' \) -delete
}

fail=0

# --- legs 1+2: beta cadence, the committed configuration ----------------------
BETA="$CASE_DIR/work-cpu"
echo ">> [rmf-2realm-fd-pulse] running beta cadence (stage_coincident, as committed)"
run_in "$BETA" 's/^$/&/'   # identity: use the manifest as-is

if grep -qiE 'error|abort| nan |segfault' "$BETA/run.log"; then
   echo "FAIL [rmf-2realm-fd-pulse] beta run reported an error/abort/NaN"; fail=1
fi
if ! grep -qE 'progress:[[:space:]]*100%' "$BETA/run.log"; then
   echo "FAIL [rmf-2realm-fd-pulse] beta time loop did not reach 100%"; fail=1
fi
for spec in "r1:$R1_HIST" "r2:$R2_HIST"; do
   realm="${spec%%:*}" ; hist="$BETA/${spec#*:}"
   if [[ ! -f "$hist" ]]; then
      echo "FAIL [rmf-2realm-fd-pulse] beta $realm divergence history missing: $(basename "$hist")"
      fail=1 ; continue
   fi
   divd="$(max_div_d "$hist")" ; divb="$(max_div_b "$hist")"
   echo ">> [rmf-2realm-fd-pulse] beta $realm max|div(D)| = $divd  max|div(B)| = $divb  (both <= $DIV_TOL)"
   if ! awk "BEGIN{exit !($divd <= $DIV_TOL)}"; then
      echo "FAIL [rmf-2realm-fd-pulse] beta $realm div(D) above round-off — source-free d/dt divD = -divJ = 0 broken"
      fail=1
   fi
   if ! awk "BEGIN{exit !($divb <= $DIV_TOL)}"; then
      echo "FAIL [rmf-2realm-fd-pulse] beta $realm div(B) above round-off — a 1:1 inter-realm mirror seam is"
      echo "                           leaking. It is a pure peer-interior copy with no resolution jump, so"
      echo "                           the mimetic identity must hold exactly across it (issue #31)."
      fail=1
   fi
done

# --- leg 3: alpha cadence, the negative control -------------------------------
ALPHA="$CASE_DIR/work-cpu-alpha"
echo ">> [rmf-2realm-fd-pulse] running alpha cadence (end_of_step) — must LEAK, proving the test has teeth"
run_in "$ALPHA" 's/^coupling_cadence = stage_coincident.*$/coupling_cadence = end_of_step/'

if ! grep -qE 'coupling_cadence[[:space:]]*=[[:space:]]*end_of_step' "$ALPHA/input.ini"; then
   echo "FAIL [rmf-2realm-fd-pulse] alpha leg did not rewrite coupling_cadence — the sed anchor drifted"
   fail=1
elif grep -qiE 'error|abort| nan |segfault' "$ALPHA/run.log"; then
   echo "FAIL [rmf-2realm-fd-pulse] alpha run reported an error/abort/NaN (expected a clean run that leaks)"
   fail=1
else
   a_divb_r1="$(max_div_b "$ALPHA/$R1_HIST")"
   a_divb_r2="$(max_div_b "$ALPHA/$R2_HIST")"
   echo ">> [rmf-2realm-fd-pulse] alpha max|div(B)| = $a_divb_r1 (r1), $a_divb_r2 (r2) — expect > $DIV_TOL"
   if awk "BEGIN{exit !($a_divb_r1 <= $DIV_TOL && $a_divb_r2 <= $DIV_TOL)}"; then
      echo "FAIL [rmf-2realm-fd-pulse] alpha cadence did NOT leak div(B). Either alpha became"
      echo "                           div-preserving (good news, but the beta leg then proves nothing and"
      echo "                           this oracle must be re-derived), or the seam is no longer being"
      echo "                           exercised at all and the beta leg is passing vacuously."
      fail=1
   fi
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-2realm-fd-pulse] 1:1 inter-realm seam div-free under beta; alpha leaks as expected"
   exit 0
else
   echo "FAIL [rmf-2realm-fd-pulse] inter-realm 1:1 seam divergence check failed"
   exit 1
fi
