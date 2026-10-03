#!/usr/bin/env bash
# FLUME verification of the multi-realm coupling (issue #37): 2-realm Sod against the single-realm Sod.
#
# Why: a forest manifest glues realms through inter-realm seams. sod-2realm.ini splits the sod-x domain at x = 0.5
# (the diaphragm, so the seam carries the Riemann problem from the first step) into two realms with the sod-x cell
# size; with a mirror seam filled at every Runge-Kutta stage (beta cadence) the seam is a block interface like any
# other, so the union of the two realms must reproduce the single-realm sod-x bit for bit. multirealm_oracle.py
# compares the interior cells by their centres (tolerance 0) and the sum of the realms' conservation histories with
# the single-realm one (round-off: the realms sum their cells separately).
#
# Leg 2 (seam + AMR): the same split with ni = 48 / 24 (the 2:1 refinement needs even block cells) and realm 2
# refined on x > 0.75, against the single-realm sod-amr with the same refinement. Realm 2 then holds both kinds of
# flux register face: the intra-realm 2:1 face at x = 0.75, crossed by the shock at t ~ 0.14, and the fine side of the
# inter-realm seam. Measured: CPU and FNL bitwise on all 135168 cells at t = 0.2 (216 steps), once the forest composed
# the two registrations (before, the manifest path registered the inter-realm faces only: the 2:1 face got no reflux).
#
# Measured (issue #37, np 2, t = 0.2, 174 steps): CPU and FNL bitwise on all 51200 cells, the summed conservation
# histories within 2e-13. The first runs found three defects: the CPU BC routine stopped on the forest's BC_SEAM crown
# rows; the fine side of the inter-realm reflux register was 2:1-restricted like an AMR seam, which wrote F_coarse - 0
# into three quarters of the realm-1 seam skin at every step; the FNL backend never copied the forest-built seam and
# BC maps to the device (illegal address in the seam fill kernel).
#
# Leg 3 (issue #40, needs N >= 2): the inter-realm seam is rank-local, its fill and reflux never cross ranks.
# sod-2realm-z.ini is the leg-1 split along z: valid on one rank (the union reproduces sod-z bit for bit), but on two
# ranks each realm's blocks are split by z, so the two sides of the seam land on different ranks. Before the guard it
# ran to the end with exit code 0 and a wrong solution (205 steps instead of 174, realm 1 mass frozen at 0.5, reflux
# mismatch 69); the leg requires the one-rank run to match sod-z and the N-rank run to stop with the issue #40 message.
#
# Usage: ./check.sh [--build] [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N)" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode flume-cpu-gnu)
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")-np$NP"

run() { # run <work-dir> <input> <files...>; RUN_NP overrides the rank count, RUN_MUST_FAIL=1 expects a failure
   local work="$1" input="$2" np="${RUN_NP:-$NP}" rc=0
   shift 2
   rm -rf "$work" ; mkdir -p "$work" ; cp "$@" "$work/"
   echo ">> $(basename "$work"): mpirun -np $np $(basename "$EXE") $input"
   (cd "$work" && mpirun -np "$np" "$EXE" "$input" > log.txt 2>&1) || rc=$?
   if [[ "${RUN_MUST_FAIL:-0}" == 1 ]]; then
      if [[ $rc -eq 0 ]]; then
         echo "check.sh: run succeeded but must fail, see $work/log.txt" >&2
         exit 1
      fi
   elif [[ $rc -ne 0 ]]; then
      echo "check.sh: run failed, see $work/log.txt" >&2
      exit 1
   fi
}

echo "== leg 1: 2-realm Sod vs sod-x"
single="$CASE_DIR/work-$TAG-single"
multi="$CASE_DIR/work-$TAG-2realm"
run "$single" sod-x.ini "$VERIF_DIR/sod/sod-x.ini"
run "$multi" sod-2realm.ini "$CASE_DIR/sod-2realm.ini" "$CASE_DIR/sod-2realm-r1.ini" "$CASE_DIR/sod-2realm-r2.ini"
"$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0

echo "== leg 2: seam + AMR, 2-realm Sod with realm 2 refined on x > 0.75 vs sod-amr"
single="$CASE_DIR/work-$TAG-amr-single"
multi="$CASE_DIR/work-$TAG-amr-2realm"
run "$single" sod-amr.ini "$CASE_DIR/sod-amr.ini"
run "$multi" sod-amr-2realm.ini "$CASE_DIR/sod-amr-2realm.ini" "$CASE_DIR/sod-amr-2realm-r1.ini" \
    "$CASE_DIR/sod-amr-2realm-r2.ini"
"$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
if [[ $NP -ge 2 ]]; then
   echo "== leg 3: cross-rank seam, z-split 2-realm Sod: matches sod-z on one rank, refused on $NP (issue #40)"
   zfiles=("$CASE_DIR/sod-2realm-z.ini" "$CASE_DIR/sod-2realm-z-r1.ini" "$CASE_DIR/sod-2realm-z-r2.ini")
   single="$CASE_DIR/work-$TAG-z-single"
   multi="$CASE_DIR/work-$TAG-z-2realm-np1"
   RUN_NP=1 run "$single" sod-z.ini "$VERIF_DIR/sod/sod-z.ini"
   RUN_NP=1 run "$multi" sod-2realm-z.ini "${zfiles[@]}"
   "$VENV_PY" "$CASE_DIR/multirealm_oracle.py" "$multi" "$single" --ngc 3 --tol 0
   refused="$CASE_DIR/work-$TAG-z-2realm"
   RUN_MUST_FAIL=1 run "$refused" sod-2realm-z.ini "${zfiles[@]}"
   if ! grep -aq 'error stop forest_object%populate_inter_realm_topology: .*(issue #40)' "$refused/log.txt"; then
      echo "check.sh: the cross-rank seam run failed without the issue #40 message, see $refused/log.txt" >&2
      exit 1
   fi
   echo "   refused: $(grep -a 'issue #40' "$refused/log.txt" | head -1 | cut -c1-120)..."
else
   echo "== leg 3 skipped: the cross-rank seam needs --np 2 or more"
fi
echo "multi-realm verification PASSED ($TAG)"
