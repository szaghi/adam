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
# div(J) truthfulness leg (issue #26 G1.a/G2, DEFAULT-ON since G2): asserts the
# seam max|div(J)| within a FACTOR-3 band of the CPU-pinned baseline.
# Why factor 3 and not the ±5% used for div(B): div(J) of the coil filament is
# a TEN-DECADE cancellation — the stencil cancels |J|/h ~ 1e+6-scale gradients
# down to ~3e-5 — so any relative J noise n reappears in div(J) as ~n*1e+6.
# The irreducible cross-compiler noise in the coil-geometry initialization
# (gfortran vs nvfortran host code, ~1e-11 relative on J after the amplitude)
# puts a ~1e-5 ABSOLUTE noise floor under this diagnostic at this resolution:
# cross-backend agreement tighter than a small factor is mathematically
# impossible, while the defect class this leg guards against (#26 G2: unstamped
# J ghosts + one-dt-lagged stamp) sat SIX ORDERS off. Band history: G1.a red at
# 2.35E+01; after the G2 stamp fixes, FNL reads 2.62E-05 vs CPU 3.22E-05.
# The legacy --divj flag is accepted as a no-op (the leg always runs).
#
# mpirun and the GNU MPI toolchain must be on PATH (see rmf-amr/run.sh header).
#
# PRISM_EXE (issue #22, GA4): override the executable under test, e.g.
#   PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH + the
# WSL UCX knobs of issue #12) and the build (--build always builds the CPU
# default, never the override). Baselines are CPU-pinned; the ±5% band is the
# cross-backend/compiler-noise allowance.

set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"

DIVB_TOL="1.0E-13"           # control: round-off ceiling for the interior identity
SEAM_BASELINE="1.555802E-07" # pinned seam max|div(B)| (see provenance above)
SEAM_RTOL="0.05"             # 5% relative band around the baseline
SEAM_J_BASELINE="3.217958E-05" # pinned seam max|div(J)| (CPU N5-era binary,
                               # captured 2026-07-04 at #26 G1.a: ni=16, it_max=5, -np 1)
SEAM_J_FACTOR="3.0"            # factor band around the div(J) baseline (see header)

do_build=0
for arg in "$@"; do
   case "$arg" in
      --build) do_build=1 ;;
      --divj)  ;; # legacy opt-in flag (issue #26 G1.a), now a no-op: the leg is default-on
      *) echo "ERROR: unknown flag $arg (use --build)" >&2; exit 2 ;;
   esac
done
if [[ $do_build -eq 1 ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# Largest B_divergence (column 5) over all rows of a divergence_history file.
max_div_b() { grep '^+' "$1" | awk '{v=$5; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }
# Largest J_divergence (column 6) over all rows of a divergence_history file.
max_div_j() { grep '^+' "$1" | awk '{v=$6; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }

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

seam_divj="$(max_div_j "$WORK/$HIST")"
echo ">> [rmf-amr-fd] seam   max|div(J)| = $seam_divj (baseline $SEAM_J_BASELINE, factor band x/$SEAM_J_FACTOR)"
if ! awk "BEGIN{r=$seam_divj/$SEAM_J_BASELINE; exit !(r>=1.0/$SEAM_J_FACTOR && r<=$SEAM_J_FACTOR)}"; then
   echo "FAIL [rmf-amr-fd] seam max|div(J)| off the CPU-pinned baseline beyond the noise-floor"
   echo "                  band — the div(J) diagnostic is not truthful on this backend"
   echo "                  (J stamping extent/order or the divergence path; see issue #26 G2)"
   fail=1
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr-fd] control at round-off, seam div(B) on the pinned truncation baseline"
   exit 0
else
   echo "FAIL [rmf-amr-fd] seam div(B) regression check failed"
   exit 1
fi
