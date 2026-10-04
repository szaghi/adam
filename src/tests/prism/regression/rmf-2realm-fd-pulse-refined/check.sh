#!/usr/bin/env bash
# rmf-2realm-fd-pulse-refined oracle (issue #52, inter-realm 2:1 seam).
#
# THE CLAIM UNDER TEST: rmf-amr-fd-pulse split at its 2:1 face x = 0 into a fine realm
# (x < 0) and a coarse realm (x > 0) glued by `coupling = refined` reproduces the
# single-realm run bit for bit. The inter-realm 2:1 seam reuses the intra-realm 2:1
# formulas (coarse->fine ghost interpolation, fine->coarse 2x2x2 mean, the 2:1 flux
# register), so nothing about the split may show in the solution.
#
# This includes the #29 seam div(B) source: a 2:1 jump injects an O(h^p) div(B) at the
# seam (accept-truncation, see the PRISM gotchas in CLAUDE.md), and the split must carry
# the SAME source, not a smaller or a larger one.
#
# ACCEPTANCE
#   1. both runs reach 100% and report no error/abort/NaN;
#   2. every checkpoint field (B, D, J, div_*, res_*) of the union of the realms equals
#      the single-realm field on every cell (multirealm_oracle.py --tol 0);
#   3. the div(D) and div(B) histories: at every step the single-realm maximum equals
#      the larger of the two realms' maxima, to the printed digits;
#   4. the seam div(B) monitor ([IO].seam_divB_tol, warn-only, set to 1 by this script)
#      fires on the split: the guard-rail of #29 covers inter-realm 2:1 seams too.
#
# Measured (issue #52, CPU and FNL, np 2, 5 steps): 147456 cells bitwise on all 27
# fields; max|div(B)| 6.86 -> 12.65 in both runs, the coarse realm holding the maximum
# at steps 1-2 and the fine realm at steps 3-5. On 2 ranks every seam row crosses ranks.
#
# Usage: ./check.sh [--build] [--np N]
#
# PRISM_EXE: override the executable under test, e.g.
#   PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH plus the
# WSL UCX knobs of issue #12); --build always builds the CPU default.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
ORACLE="$REPO_ROOT/src/tests/flume/verification/multirealm/multirealm_oracle.py"
SINGLE_INI="$CASE_DIR/../rmf-amr-fd-pulse/input.ini"
NP=2

do_build=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) do_build=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *) echo "ERROR: unknown flag $1 (use --build, --np N)" >&2; exit 2 ;;
   esac
done
if [[ $do_build -eq 1 ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

TAG="$(basename "$EXE")-np$NP"
SINGLE="$CASE_DIR/work-$TAG-single"
SPLIT="$CASE_DIR/work-$TAG-refined"
fail=0

run_in() { # workdir
   ( cd "$1" && timeout 600 mpirun -np "$NP" "$EXE" > run.log 2>&1 ) || true
   if grep -qiE 'error|abort| nan |segfault' "$1/run.log"; then
      echo "FAIL [rmf-2realm-fd-pulse-refined] $(basename "$1") reported an error/abort/NaN"; fail=1
   fi
   if ! grep -qE 'progress:[[:space:]]*100%' "$1/run.log"; then
      echo "FAIL [rmf-2realm-fd-pulse-refined] $(basename "$1") did not reach 100%"; fail=1
   fi
}

rm -rf "$SINGLE" "$SPLIT" ; mkdir -p "$SINGLE" "$SPLIT"
cp "$SINGLE_INI" "$SINGLE/input.ini"
cp "$CASE_DIR/input.ini" "$SPLIT/"
for r in 1 2; do # leg 4: arm the seam div(B) monitor, warn-only
   sed 's/^divergence_history_save = 1 .*/&\nseam_divB_tol           = 1.0/' "$CASE_DIR/realm_$r.ini" > "$SPLIT/realm_$r.ini"
done
echo ">> [rmf-2realm-fd-pulse-refined] single-realm rmf-amr-fd-pulse ($TAG)"
run_in "$SINGLE"
echo ">> [rmf-2realm-fd-pulse-refined] refined 2-realm split ($TAG)"
run_in "$SPLIT"

echo ">> [rmf-2realm-fd-pulse-refined] fields: union of the realms vs single realm, bitwise"
"$VENV_PY" "$ORACLE" "$SPLIT" "$SINGLE" --ngc 3 --tol 0 --fields-only | sed 's/^/   /' || fail=1

echo ">> [rmf-2realm-fd-pulse-refined] div histories: single max == max over the realms, every step"
# columns: it  blocks_number  time  D_divergence  B_divergence  J_divergence
if ! awk 'FNR == 1 { f++ } /^\+/ { for (c = 4; c <= 5; c++) { v = $c + 0; if (f == 1) s[$1, c] = v;
             else if (!(($1, c) in m) || v > m[$1, c]) m[$1, c] = v } }
          END { n = 0; bad = 0
                for (k in s) { n++; if (s[k] != m[k]) { bad++; split(k, p, SUBSEP)
                   printf "   it %d col %d: single %.17g, realms %.17g\n", p[1], p[2], s[k], m[k] } }
                printf "   %d (step, column) pairs, %d differ\n", n, bad; exit (bad > 0 || n == 0) }' \
        "$SINGLE/rmf_amr_fd_pulse_regression-divergence_history.dat" \
        "$SPLIT/rmf_2realm_fd_pulse_refined_r1-divergence_history.dat" \
        "$SPLIT/rmf_2realm_fd_pulse_refined_r2-divergence_history.dat"; then
   echo "FAIL [rmf-2realm-fd-pulse-refined] the split does not carry the single-realm div(D)/div(B)"; fail=1
fi

echo ">> [rmf-2realm-fd-pulse-refined] seam div(B) monitor (tol 1.0) on the split"
if grep -aq 'WARNING: seam div(B) .* exceeds \[IO\].seam_divB_tol' "$SPLIT/run.log"; then
   echo "   fired: $(grep -ac 'WARNING: seam div(B)' "$SPLIT/run.log") warnings"
else
   echo "FAIL [rmf-2realm-fd-pulse-refined] the seam div(B) monitor did not fire on the inter-realm 2:1 seam"; fail=1
fi

for w in "$SINGLE" "$SPLIT"; do
   find "$w" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' -o -name '*.tnd' \) -delete
done
if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-2realm-fd-pulse-refined] refined 2-realm split == single-realm rmf-amr-fd-pulse ($TAG)"
   exit 0
fi
echo "FAIL [rmf-2realm-fd-pulse-refined] ($TAG)"
exit 1
