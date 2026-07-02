#!/usr/bin/env bash
# rmf-amr-fd Phase-B seam div(B) check (issue #13 Phase B).
#
# This is the fail-before / pass-after oracle for Phase B. It runs the
# fd_centered single-realm case with a static intra-realm 2:1 AMR jump and
# measures max|div_h(B)| over the domain per step:
#
#   - In the FD-centered scheme the B-update is dq(B) = -curl_h(D) and the div(B)
#     diagnostic is div_h(B); both use the SAME centered stencil
#     (compute_derivative1_fd_centered / FD1_CC), so div_h(curl_h D) = 0 to
#     round-off on a uniform grid and div(B) is conserved exactly by the interior
#     scheme. The marker-OFF control confirms this: div(B) stays at round-off.
#
#   - With the 2:1 seam present (marker ON), the coarse and fine curl_h(D) read
#     each other's ghosts via 0th-order injection -> the two one-sided stencils
#     are not the same operator -> a dt-rate div(B) source is injected at seam
#     cells and div(B) grows. Phase A reflux does NOT fix this (different
#     operator; verified empirically at M4).
#
# Phase B shares a canonical D across the seam so curl_h(D) is single-valued
# there and the cancellation is restored. ACCEPTANCE: seam-case max|div(B)| drops
# to round-off (<= DIVB_TOL) across the run.
#
# Status of THIS check:
#   - BEFORE Phase B: FAILS (seam div(B) ~ 1e-6, far above tol) — this is the
#     intended fail-before state proving the leak is real and reflux-immune.
#   - AFTER  Phase B: PASSES (seam div(B) <= tol).
# A control run (marker OFF) must ALWAYS pass — it proves the leak is the seam,
# not the fd_centered scheme itself.
#
# Usage: ./check.sh            (expects exe/adam_prism_cpu already built)
#        ./check.sh --build    (build prism-cpu-gnu first)
#
# mpirun and the GNU MPI toolchain must be on PATH (see rmf-amr/run.sh header).

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="$REPO_ROOT/exe/adam_prism_cpu"

# Phase-B acceptance tolerance on seam max|div(B)|. The interior scheme holds
# div(B) at round-off (~1e-18 in the control); 1e-13 is a generous round-off
# ceiling that the seam leak (~1e-6) violates by 7 orders of magnitude today.
DIVB_TOL="1.0E-13"

if [[ "${1:-}" == "--build" ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# Largest B_divergence (column 5) over all rows of a divergence_history file.
max_div_b() { grep '^+' "$1" | awk '{v=$5; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }

run_in() { # workdir, ini-transform-sed
   local wd="$1" sed_expr="$2"
   rm -rf "$wd" && mkdir -p "$wd"
   sed "$sed_expr" "$CASE_DIR/input.ini" > "$wd/input.ini"
   ( cd "$wd" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 )
}

fail=0
HIST="rmf_amr_fd_regression-divergence_history.dat"

# --- Control run: marker OFF, no seam. div(B) must stay at round-off. ---
CTRL="$CASE_DIR/work-cpu-ctrl"
echo ">> [rmf-amr-fd] running marker-disabled control (no seam)"
run_in "$CTRL" 's/markers_number = 1/markers_number = 0/'
if grep -qiE 'error|abort| nan |segfault' "$CTRL/run.log"; then
   echo "FAIL [rmf-amr-fd] control reported an error/abort/NaN"; fail=1
fi
ctrl_divb="$(max_div_b "$CTRL/$HIST")"
echo ">> [rmf-amr-fd] control max|div(B)| = $ctrl_divb (expect <= $DIVB_TOL)"
if ! awk "BEGIN{exit !($ctrl_divb <= $DIVB_TOL)}"; then
   echo "FAIL [rmf-amr-fd] control div(B) above round-off — fd matched-stencil identity broken without a seam"; fail=1
fi

# --- Main run: marker ON, 2:1 seam present. ---
WORK="$CASE_DIR/work-cpu"
echo ">> [rmf-amr-fd] running marker-active case (2:1 seam)"
run_in "$WORK" 's/^$/&/'   # identity (copy as-is)
if grep -qiE 'error|abort| nan |segfault' "$WORK/run.log"; then
   echo "FAIL [rmf-amr-fd] run reported an error/abort/NaN"; fail=1
fi
if ! grep -qE 'progress:[[:space:]]*100%' "$WORK/run.log"; then
   echo "FAIL [rmf-amr-fd] fd_centered time loop did not reach 100%"; fail=1
fi
seam_divb="$(max_div_b "$WORK/$HIST")"
echo ">> [rmf-amr-fd] seam   max|div(B)| = $seam_divb (Phase-B acceptance: <= $DIVB_TOL)"

# Phase-B acceptance: seam div(B) at round-off. FAILS before Phase B (the leak),
# PASSES after the canonical-D seam sharing lands.
if ! awk "BEGIN{exit !($seam_divb <= $DIVB_TOL)}"; then
   echo "FAIL [rmf-amr-fd] PHASE-B NOT MET: seam max|div(B)| = $seam_divb > $DIVB_TOL"
   echo "                  (this is the intended fail-before state until the canonical-D"
   echo "                   seam sharing restores div_h(curl_h D) = 0 across the 2:1 jump)"
   fail=1
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr-fd] Phase B: seam max|div(B)| at round-off across the 2:1 AMR jump"
   exit 0
else
   echo "FAIL [rmf-amr-fd] Phase B seam div(B) check failed"
   exit 1
fi
