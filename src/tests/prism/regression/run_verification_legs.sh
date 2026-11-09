#!/usr/bin/env bash
# Run the PRISM verification legs documented in docs/applications/prism/verification/, one case after the other, and
# log each case to <outdir>/<case>.log. The logs carry the measured values the pages quote (the "<run>" cells).
#
# Usage:
#   run_verification_legs.sh [--fnl] [--only "case ..."] [--refusals] <repo-root> <exe> <outdir>
#
#   <repo-root>  the ADAM tree whose src/tests/prism/regression/ scripts are run; the check.sh scripts write their
#                work-* directories inside it, so pass a copy of the tree on a machine where the source tree must stay
#                clean
#   <exe>        the PRISM executable under test (exported as PRISM_EXE): exe/adam_prism_cpu, or exe/adam_prism_fnl
#                with --fnl
#   <outdir>     where <case>.log (and the work directory of periodic-pulse-blocks, which has no check.sh) go
#   --fnl        FNL backend: the PIC cases run their fnl legs (FNL_NP, default "1 2"), the slices run on
#                SLICES_NP (default "1 2"); the caller owns the NVHPC environment (and the WSL UCX settings)
#   --only LIST  run only the named cases (blank-separated), in page order
#   --refusals   also run refusals.sh (contract legs, no measured value) into <outdir>/refusals.log
#
# Cases, in the order of the verification pages:
#   1  periodic-pulse-blocks           run at 2 ranks; energy drift and divergence summary
#   2  rmf-amr                         check.sh
#   3  rmf-amr-fd, rmf-amr-fd-pulse    check.sh; rmf-amr-fd-pulse with --convergence (its log also holds page 4, the
#                                      quadtree leg)
#   5  rmf-2realm-fd-pulse             check.sh
#   6  rmf-2realm-fv-pulse-blocks      check.sh
#   7  rmf-2realm-fd-pulse-refined     check.sh
#   8  rmf-amr-fv-pulse-timed          check.sh
#   11 rmf-pec-t0 rmf-open-t0 pic-pec-t0 pic-open-t0   check.sh
#   12 pic-periodic pic-slab-pulse     check.sh <cpu|fnl>
#   13 slices                          check.sh --convergence
# Pages 9 and 10 are digest and refusal contracts (run.sh, refusals.sh): no measured value.
#
# Needs mpirun on PATH and python3 (a venv with h5py and numpy is created at <repo-root>/exe/.regression-venv when
# missing; PRISM_PY overrides the python the t0 and PIC scripts use). Exit status: 0 when every leg passes, 1 when one
# fails, 2 on a usage error.
set -uo pipefail

usage() { sed -n '2,/^set -uo/{/^set -uo/d;s/^# \{0,1\}//;p}' "$0"; }

BACKEND=cpu
ONLY=""
REFUSALS=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --fnl)      BACKEND=fnl ; shift ;;
      --only)     [[ $# -ge 2 ]] || { echo "ERROR: --only needs a list" >&2; exit 2; }; ONLY="$2" ; shift 2 ;;
      --refusals) REFUSALS=1 ; shift ;;
      -h|--help)  usage ; exit 0 ;;
      --*)        echo "ERROR: unknown option $1" >&2 ; usage >&2 ; exit 2 ;;
      *)          break ;;
   esac
done
[[ $# -eq 3 ]] || { usage >&2 ; exit 2 ; }
REPO_ROOT="$(cd "$1" 2>/dev/null && pwd)" || { echo "ERROR: repo root '$1' not found" >&2; exit 2; }
EXE="$(cd "$(dirname "$2")" 2>/dev/null && pwd)/$(basename "$2")"
OUTDIR="$3"
REG="$REPO_ROOT/src/tests/prism/regression"
[[ -d "$REG" ]] || { echo "ERROR: $REG not found" >&2; exit 2; }
[[ -x "$EXE" ]] || { echo "ERROR: executable $EXE not found" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }
[[ "$(basename "$EXE")" == *_"$BACKEND" ]] || echo "WARNING: $(basename "$EXE") does not end in _$BACKEND (--fnl selects fnl)"
mkdir -p "$OUTDIR" && OUTDIR="$(cd "$OUTDIR" && pwd)"
export PRISM_EXE="$EXE"

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c "import h5py, numpy" 2>/dev/null; then
   echo ">> creating the python venv $REPO_ROOT/exe/.regression-venv (h5py, numpy)"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv" && "$VENV_PY" -m pip install --quiet h5py numpy \
      || { echo "ERROR: cannot create the python venv" >&2; exit 2; }
fi

CASES=(periodic-pulse-blocks rmf-amr rmf-amr-fd rmf-amr-fd-pulse rmf-2realm-fd-pulse rmf-2realm-fv-pulse-blocks
       rmf-2realm-fd-pulse-refined rmf-amr-fv-pulse-timed rmf-pec-t0 rmf-open-t0 pic-pec-t0 pic-open-t0
       pic-periodic pic-slab-pulse slices)
for name in $ONLY; do
   [[ " ${CASES[*]} " == *" $name "* ]] || { echo "ERROR: --only names '$name', not a verification case" >&2; exit 2; }
done

# periodic-pulse-blocks has no check.sh: run it at 2 ranks (as run.sh does) and summarise its energy and divergence
# histories (columns: it blocks time D_energy B_energy ...; it blocks time D_divergence B_divergence J_divergence)
uniform_leg() {
   local wd="$OUTDIR/work-periodic-pulse-blocks-$BACKEND" base=periodic_pulse_blocks
   rm -rf "$wd" && mkdir -p "$wd" && cp "$REG/periodic-pulse-blocks/input.ini" "$wd/"
   ( cd "$wd" && mpirun -np 2 "$EXE" input.ini > run.log 2>&1 ) || { echo "FAIL [periodic-pulse-blocks] run failed, see $wd/run.log"; return 1; }
   [[ -f "$wd/$base-energy_history.dat" && -f "$wd/$base-divergence_history.dat" ]] \
      || { echo "FAIL [periodic-pulse-blocks] energy or divergence history missing in $wd"; return 1; }
   echo ">> energy history, first and last rows:"
   grep '^+' "$wd/$base-energy_history.dat" | sed -n '1p;$p'
   grep '^+' "$wd/$base-energy_history.dat" | awk '
      NR == 1 { w0 = $4 + $5 }
      { w = $4 + $5; d = $4 - $5; if (d < 0) d = -d; r = d / w; if (r > m) m = r; last = w; it = $1 }
      END { dr = (last - w0) / w0; if (dr < 0) dr = -dr
            printf ">> total energy relative drift, first row to it %d: %.3e\n", it, dr
            printf ">> max |W_D - W_B| / (W_D + W_B): %.3e\n", m }'
   grep '^+' "$wd/$base-divergence_history.dat" | awk '
      { for (c = 4; c <= 5; c++) { v = $c; if (v < 0) v = -v; if (v > m[c]) m[c] = v; l[c] = v } }
      END { printf ">> max |div D| %.3e, max |div B| %.3e; last row |div D| %.3e, |div B| %.3e\n", m[4], m[5], l[4], l[5] }'
   find "$wd" -type f \( -name '*.h5' -o -name '*.xdmf' \) -delete
   echo "PASS [periodic-pulse-blocks] run completed (no bound: the values are reported)"
}

case_cmd() { # the command of one case, as an array in CMD
   local c="$1"
   case "$c" in
      periodic-pulse-blocks)          CMD=(uniform_leg) ;;
      rmf-amr-fd-pulse)               CMD=("$REG/$c/check.sh" --convergence) ;;
      pic-periodic|pic-slab-pulse)    CMD=("$REG/$c/check.sh" "$BACKEND") ;;
      slices)
         local nps="${SLICES_NP:-}"
         [[ -n "$nps" ]] || { [[ $BACKEND == fnl ]] && nps="1 2" || nps="1 2 4"; }
         CMD=("$REG/$c/check.sh" --np "$nps" --convergence) ;;
      *)                              CMD=("$REG/$c/check.sh") ;;
   esac
}

declare -a SUMMARY=()
failed=0
run_case() { # case name, log name
   local c="$1" log="$OUTDIR/$2.log" rc start
   case_cmd "$c"
   echo "============================================================"
   echo "== [$c/$BACKEND] ${CMD[*]}  -> $log"
   start=$(date +%s)
   {
      echo "# case $c, backend $BACKEND, exe $EXE, repo $REPO_ROOT, $(date -Is)"
      echo "# command: ${CMD[*]}"
      "${CMD[@]}"
   } > "$log" 2>&1
   rc=$?
   echo "# exit status $rc, $(( $(date +%s) - start )) s" >> "$log"
   if [[ $rc -eq 0 ]]; then SUMMARY+=("PASS  $c ($(( $(date +%s) - start )) s)")
   else SUMMARY+=("FAIL  $c (exit $rc, see $log)"); failed=1; fi
   grep -aE '^(PASS|FAIL)|(PASSED|FAILED)' "$log" | tail -n 1 | sed 's/^/   /'
}

for c in "${CASES[@]}"; do
   [[ -z "$ONLY" || " $ONLY " == *" $c "* ]] || continue
   run_case "$c" "$c"
done
if [[ $REFUSALS -eq 1 ]]; then
   CMD=("$REG/refusals.sh")
   { echo "# refusals.sh, backend $BACKEND, exe $EXE"; "${CMD[@]}"; } > "$OUTDIR/refusals.log" 2>&1
   rc=$?
   [[ $rc -eq 0 ]] && SUMMARY+=("PASS  refusals") || { SUMMARY+=("FAIL  refusals (exit $rc)"); failed=1; }
fi

echo "============================================================"
echo "== PRISM verification legs, backend $BACKEND, logs in $OUTDIR"
printf '   %s\n' "${SUMMARY[@]}"
exit $failed
