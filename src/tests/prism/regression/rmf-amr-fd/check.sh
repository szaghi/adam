#!/usr/bin/env bash
# rmf-amr-fd seam div(B) regression check (issue #21, N4 acceptance).
#
# The fd_centered single-realm case with a static intra-realm 2:1 AMR jump.
# On a uniform grid the matched-stencil identity div_h(curl_h D) = 0 holds to
# round-off (control run). At the seam, coarse<->fine ghosts are filled by the
# tricubic seam interpolation ([amr] seam_ghost_fill, default tricubic; #21
# N2) over the full ghost slab (#21 N3.5), so the residual seam div(B) is a
# TRUNCATION-ORDER quantity, not a defect.
#
# ACCEPTANCE (redefined per #21 — the historical "seam div(B) <= 1e-13"
# pointwise criterion is RETIRED, see the regression README and #21 §1):
#   1. control (marker OFF) max|div(B)| <= DIVB_TOL  — the interior identity;
#   2. seam (marker ON) max|div(B)| within RTOL of the pinned baseline
#      (golden-style: catches regressions AND silent improvements).
#
# CONVERGENCE of the seam leak under refinement is asserted by
# rmf-amr-fd-pulse/check.sh --convergence, NOT here: the coil filament makes
# div(J) scale ~ h^-2, so matched-time seam metrics on THIS case are
# contaminated by the source and are not a valid order measurement (#21 N3).
#
# Baseline provenance: N3.5 binary (stale-ghost-layer fix), default tricubic
# fill, ni=16, it_max=5, -np 1 — max over the divergence history (= step 5).
#
# Usage: ./check.sh            (expects exe/adam_prism_cpu already built)
#        ./check.sh --build    (build prism-cpu-gnu first)
#
# mpirun and the GNU MPI toolchain must be on PATH (see rmf-amr/run.sh header).

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="$REPO_ROOT/exe/adam_prism_cpu"

DIVB_TOL="1.0E-13"           # control: round-off ceiling for the interior identity
SEAM_BASELINE="1.555802E-07" # pinned seam max|div(B)| (see provenance above)
SEAM_RTOL="0.05"             # 5% relative band around the baseline

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
   # keep the small .dat/.log evidence, drop the heavy field dumps
   find "$wd" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' \) -delete
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
echo ">> [rmf-amr-fd] seam   max|div(B)| = $seam_divb (baseline $SEAM_BASELINE, rtol $SEAM_RTOL)"
if ! awk "BEGIN{d=($seam_divb-$SEAM_BASELINE)/$SEAM_BASELINE; if(d<0)d=-d; exit !(d<=$SEAM_RTOL)}"; then
   echo "FAIL [rmf-amr-fd] seam max|div(B)| moved off the pinned baseline"
   echo "                  (a DROP means an improvement landed — rebaseline deliberately;"
   echo "                   a RISE means the seam fill or the exchange regressed)"
   fail=1
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr-fd] control at round-off, seam div(B) on the pinned truncation baseline"
   exit 0
else
   echo "FAIL [rmf-amr-fd] seam div(B) regression check failed"
   exit 1
fi
