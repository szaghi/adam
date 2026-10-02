#!/usr/bin/env bash
# FLUME verification PV-1 and PV-4 (issue #47, M3-P5a, M3-P5b): the positivity limiter, `[numerics] positivity_limiter = cell`.
#
# Why: the limiter (D-9) blends every face flux with the first-order Lax-Friedrichs backbone so that each stage keeps
# the density and the pressure positive; it must make the hard problems run without floors, and leave the smooth ones
# untouched. np 2, the floors disabled (a non-positive density or pressure stops the run):
#   1. PV-1, the Balsara-Spicer blast (beta 2.5e-4, 128^2, t = 0.01, EGLM, ../positivity-probe/make_blast.py) with the
#      limiter, on the splitting scheme (which fails at step 3 without it, M3-P4c) and on weno-riemann HLLD
#      characteristic: the run reaches t = 0.01, and the limiter never meets an inadmissible backbone (the log line
#      `positivity limiter: N faces limited, M inadmissible backbones` has M = 0 at every stage);
#   2. PV-4, inactive on smooth problems: with the limiter on, the isentropic vortex (Euler, N = 64, 128) and the 1-D
#      fast wave (EGLM, N = 16, 32) are BITWISE equal to the runs without it (interior cells) (a face whose two cells keep Lambda = 1 is
#      not touched, and the EGLM sources with Lambda = 1 take the unlimited arithmetic), so their orders are unchanged.
# About 8 min on the CPU. One run at a time; checkpoints deleted after use.
#
# Usage: ./check.sh [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh (M3-P5b)
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
source "$CASE_DIR/../numerics.sh"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np) NP="$2" ; shift 2 ;;
      *)    echo "check.sh: unknown argument '$1' (accepted: --np N)" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

limiter_on() { # limiter_on <ini>: add positivity_limiter = cell to the [numerics] block (after reflux)
   sed -i -E 's/^(reflux\s*=.*)$/\1\npositivity_limiter       = cell/' "$1"
   grep -q '^positivity_limiter *= *cell' "$1"
}
run() { # run <work> <ini>: run one case, fail the check if the run fails
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: $(basename "$1") failed (a non-positive state stops the run), see $1/log.txt" >&2
      FAILED=1
      return 1
   fi
}
same() { # same <work-a> <work-b>: the interiors of the last checkpoints of the two runs bitwise equal (FNL ghost
         # cells are not reproducible run to run, M3-P3b)
   "$VENV_PY" - "$1" "$2" <<'EOF'
import glob
import sys

import h5py
import numpy as np

a, b = sys.argv[1:]
fa = sorted(f for f in glob.glob(a + "/*-proc*.h5") if "restart" not in f)
last = max(int(f.split("-")[-2]) for f in fa)
diff = 0.0
for f in (f for f in fa if int(f.split("-")[-2]) == last):
    with h5py.File(f) as x, h5py.File(f.replace(a, b)) as y:
        for k in x:
            if x[k].ndim == 3:
                diff = max(diff, float(np.max(np.abs(x[k][()][3:-3, 3:-3, 3:-3] - y[k][()][3:-3, 3:-3, 3:-3]))))
ok = diff == 0.0
print(f"   {b.split('/')[-1]} vs {a.split('/')[-1]}: step {last}, max |difference| {diff:.3e}  "
      f"{'PASS (bitwise)' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
}

echo ">> PV-1 Balsara-Spicer blast, 128^2, EGLM, limiter on ($(basename "$EXE"), np $NP)"
for spec in split hlld:characteristic; do
   w="$CASE_DIR/work-$TAG-blast-${spec//:/-}"
   rm -rf "$w" ; mkdir -p "$w"
   "$VENV_PY" "$VERIF_DIR/mhd/positivity-probe/make_blast.py" "$w/blast.ini" --cells 128 --eglm > /dev/null
   NUMERICS="" ; NUMERICS_SOLVER=""
   if [[ $spec != split ]]; then NUMERICS=$spec; numerics_check; numerics_apply "$w/blast.ini"; fi
   limiter_on "$w/blast.ini"
   if run "$w" blast.ini; then
      last=$(tail -1 "$w/blast-conservation_history.dat" | awk '{print $2}')
      limited=$(grep -a 'positivity limiter:' "$w/log.txt" | sed -E 's/.*: \+([0-9]+) faces.*/\1/' | sort -n | tail -1 \
                || true) # no line when the limiter never acts
      badrows=$(grep -a 'positivity limiter:' "$w/log.txt" | grep -vc '+0 inadmissible' || true)
      if "$VENV_PY" -c "import sys; sys.exit(0 if abs($last - 0.01) < 1e-12 else 1)" && [[ $badrows -eq 0 ]]; then
         echo "   $spec: t = $last, max ${limited:-0} faces limited per stage, 0 inadmissible backbones  PASS"
      else
         echo "   $spec: t = $last, $badrows stages with inadmissible backbones  FAIL" ; FAILED=1
      fi
   fi
   find "$w" -name '*.h5' -delete
done

echo ">> PV-4 the limiter is inactive on smooth problems: bitwise equal runs"
for n in 064 128; do
   for v in off on; do
      w="$CASE_DIR/work-$TAG-vortex-n$n-$v"
      rm -rf "$w" ; mkdir -p "$w"
      cp "$VERIF_DIR/vortex/vortex-n$n.ini" "$w/"
      [[ $v == on ]] && limiter_on "$w/vortex-n$n.ini"
      run "$w" "vortex-n$n.ini" || true
   done
   same "$CASE_DIR/work-$TAG-vortex-n$n-off" "$CASE_DIR/work-$TAG-vortex-n$n-on" || FAILED=1
done
for n in 16 32; do
   for v in off on; do
      w="$CASE_DIR/work-$TAG-wave-n$n-$v"
      rm -rf "$w" ; mkdir -p "$w"
      "$VENV_PY" "$VERIF_DIR/mhd/linear-wave/make_linear_wave.py" "$VERIF_DIR/sod/sod-x.ini" "$w/linear-wave.ini" \
         --wave fast --geometry 1d --cells "$n" --amplitude 1.0e-7 --divergence-control eglm
      [[ $v == on ]] && limiter_on "$w/linear-wave.ini"
      run "$w" linear-wave.ini || true
   done
   same "$CASE_DIR/work-$TAG-wave-n$n-off" "$CASE_DIR/work-$TAG-wave-n$n-on" || FAILED=1
done
find "$CASE_DIR" -path "$CASE_DIR/work-$TAG-*" -name '*.h5' -delete

if [[ $FAILED -eq 0 ]]; then
   echo "PV-1 PV-4 PASSED ($TAG)"
else
   echo "PV-1 PV-4 FAILED ($TAG)"
   exit 1
fi
