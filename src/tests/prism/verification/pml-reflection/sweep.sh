#!/usr/bin/env bash
# PRISM absorbing boundaries at normal incidence: the reflected amplitude of a 1-D Gaussian pulse (issue #145, item 3).
#
# Why: the paper's 20-cell classical PML is one point; this sweep measures the layer against its design reflection
# R = exp(-2 gamma_max L / ((n + 1) c)) (continuum, PEC-backed, frequency-independent at normal incidence):
#   - thickness N = 10, 20, 40, 80 cells at gamma_max = 7.8e10 1/s, n = 2 (the paper's layer is N = 20);
#   - gamma_max = 3.9e10 .. 6.24e11 1/s at N = 20;
#   - N = 10, 20, 40, 80 at the same design R ~ 9.4e-7 (gamma_max L fixed): what is left is the discrete reflection of
#     the grading;
#   - the CFS layer at N = 20 (kappa_max = 1, alpha_max = 0: the classical layer through the CFS path; and the
#     regression's kappa_max = 5, alpha_max = 2.5e8 1/s of rmf-pml-amr);
#   - the boundaries without a layer: first-order Silver-Mueller, Neumann, PEC;
#   - `far`: no boundary within reach, the noise floor of the measurement (the dispersive trail of the pulse).
# make_pml_reflection.py writes each input, pml_reflection.py measures A_r/A_0 and (E_r/E_0)^(1/2) in the window
# [-0.04, 0.04] m; results.txt collects one line per case, results.json the measurements.
#
# Usage: ./sweep.sh <out> [--jobs J] [--cases "name ..."] [--keep]
#   PRISM_EXE overrides the executable (default <repo>/exe/adam_prism_cpu); the runs are serial-in-space (np 1).
#   --keep keeps the checkpoints (the figures read the profiles.npz each case writes anyway).
# About 30 s per case on one core (the high-gamma cases take the PML step cap 0.1/gamma_max, up to 5000 steps).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
PY="$REPO_ROOT/exe/.regression-venv/bin/python"
OUT="${1:?usage: sweep.sh <out> [--jobs J] [--cases LIST] [--keep]}"; shift
JOBS=4 ; KEEP=0 ; ONLY=""
while [[ $# -gt 0 ]]; do
   case "$1" in
      --jobs)  JOBS="$2" ; shift 2 ;;
      --cases) ONLY="$2" ; shift 2 ;;
      --keep)  KEEP=1 ; shift ;;
      *) echo "sweep.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
[[ -x "$EXE" ]] || { echo "sweep.sh: executable '$EXE' not found (set PRISM_EXE)" >&2 ; exit 2 ; }

# name: make_pml_reflection.py options
CASES=(
   "pml-n10:--boundary pml --cells 10"
   "pml-n20:--boundary pml --cells 20"
   "pml-n40:--boundary pml --cells 40"
   "pml-n80:--boundary pml --cells 80"
   "pml-n20-g3.9e10:--boundary pml --cells 20 --gamma-max 3.9E10"
   "pml-n20-g1.56e11:--boundary pml --cells 20 --gamma-max 1.56E11"
   "pml-n20-g3.12e11:--boundary pml --cells 20 --gamma-max 3.12E11"
   "pml-n20-g6.24e11:--boundary pml --cells 20 --gamma-max 6.24E11"
   "pml-n10-g6.24e11:--boundary pml --cells 10 --gamma-max 6.24E11"
   "pml-n40-g1.56e11:--boundary pml --cells 40 --gamma-max 1.56E11"
   "cfs-n20-k1:--boundary cfs --cells 20 --k-max 1.0 --alfa-max 0.0"
   "cfs-n20-k5:--boundary cfs --cells 20 --k-max 5.0 --alfa-max 2.5E8"
   "silver-muller:--boundary silver-muller"
   "neumann:--boundary neumann"
   "pec:--boundary pec"
   "far:--boundary far"
)

mkdir -p "$OUT"
run_case() { # name options...
   local name="$1" ; shift
   local w="$OUT/$name"
   rm -rf "$w" && mkdir -p "$w"
   "$PY" "$CASE_DIR/make_pml_reflection.py" "$w/input.ini" "$@"
   local t0 ; t0=$(date +%s.%N)
   if ! (cd "$w" && mpirun -np 1 "$EXE" input.ini > run.log 2>&1); then
      echo "$name: RUN FAILED (see $w/run.log)" > "$w/result.txt" ; return 0
   fi
   local wall ; wall=$(echo "$(date +%s.%N) - $t0" | bc)
   "$PY" "$CASE_DIR/pml_reflection.py" "$w" --json "$w/result.json" --profiles "$w/profiles.npz" > "$w/result.txt"
   echo "   wall ${wall} s" >> "$w/result.txt"
   [[ $KEEP -eq 1 ]] || find "$w" -name '*.h5' -o -name '*.xdmf' | xargs -r rm -f
}
export -f run_case
export OUT PY CASE_DIR EXE KEEP

for c in "${CASES[@]}"; do
   name="${c%%:*}" ; opts="${c#*:}"
   if [[ -n $ONLY && " $ONLY " != *" $name "* ]]; then continue; fi
   echo "$name $opts"
done | xargs -P "$JOBS" -L 1 bash -c 'run_case "$@"' _

: > "$OUT/results.txt"
for c in "${CASES[@]}"; do
   name="${c%%:*}"
   [[ -f "$OUT/$name/result.txt" ]] && cat "$OUT/$name/result.txt" >> "$OUT/results.txt"
done
"$PY" - "$OUT" <<'EOF'
import json, sys
from pathlib import Path
out = Path(sys.argv[1])
res = {p.parent.name: json.loads(p.read_text()) for p in sorted(out.glob("*/result.json"))}
(out / "results.json").write_text(json.dumps(res, indent=1) + "\n")
EOF
cat "$OUT/results.txt"
