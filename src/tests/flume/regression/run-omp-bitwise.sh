#!/usr/bin/env bash
# FLUME OpenMP-CPU determinism check (issue #35): the OpenMP build must reproduce the serial build BITWISE.
#
# Why: every OpenMP region of a FLUME CPU run writes only its own cell or face, and its one reduction is a max, so the
# result must not depend on the thread count. One event on record says otherwise (issue #35, P6 note): a single-ulp
# change of one residual at iteration 195 of a 4-thread shock-cylinder run, grown to 2.3e-12 by step 268, never
# reproduced since (0 of 40 runs, 2026-09-24) and not explained by a review of the 21 regions. run.sh cannot see such
# an event (it compares within rtol 1e-6), so this check runs the immersed-boundary case, the only one that ever showed
# it, with the serial and the OpenMP executables and requires byte-identical residuals and conservation histories and an
# identical field digest. Run on every push, it turns a rare race into a reproducer with its first divergent iteration.
#
# Usage:
#   ./run-omp-bitwise.sh                     # build flume-cpu-gnu and flume-cpu-gnu-omp, np 2, 2 threads per rank
#   ./run-omp-bitwise.sh --no-build          # use the existing exe/adam_flume_cpu and exe/adam_flume_cpu_omp
#   ./run-omp-bitwise.sh --threads 4 --np 2  # threads per rank, MPI ranks
#
# Exits 0 when the two runs are bitwise identical, 1 on any difference (printing the first divergent iteration),
# 2 on a usage or setup error.
set -euo pipefail

NO_BUILD=0
THREADS=2
NP=2
while [[ $# -gt 0 ]]; do
   case "$1" in
      --no-build) NO_BUILD=1 ; shift ;;
      --threads)  THREADS="$2" ; shift 2 ;;
      --np)       NP="$2" ; shift 2 ;;
      *) echo "ERROR: unknown argument '$1' (accepted: --no-build, --threads N, --np N)" >&2 ; exit 2 ;;
   esac
done

REGRESSION_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$REGRESSION_DIR/../../../.." && pwd)"
CASE_DIR="$REGRESSION_DIR/shock-cylinder-ib"
DIGEST_PY="$REPO_ROOT/src/tests/prism/regression/digest.py"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"

if ! command -v mpirun >/dev/null 2>&1; then
   echo "ERROR: mpirun not on PATH" >&2
   exit 2
fi
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating digest venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

cd "$REPO_ROOT"
if [[ $NO_BUILD -eq 0 ]]; then
   fobis build --mode flume-cpu-gnu
   fobis build --mode flume-cpu-gnu-omp
fi
for exe in exe/adam_flume_cpu exe/adam_flume_cpu_omp; do
   [[ -x "$exe" ]] || { echo "ERROR: executable $exe not found" >&2 ; exit 2 ; }
done

ob="$(sed -n -E 's/^[[:space:]]*output_basename[[:space:]]*=[[:space:]]*([^[:space:];]+).*/\1/p' "$CASE_DIR/input.ini" | head -n 1)"

run() { # run <work-dir> <exe> <threads>
   rm -rf "$1" ; mkdir -p "$1" ; cp "$CASE_DIR/input.ini" "$1/"
   echo ">> $(basename "$2"): mpirun -np $NP, OMP_NUM_THREADS=$3"
   if ! (cd "$1" && OMP_NUM_THREADS="$3" mpirun -np "$NP" "$REPO_ROOT/$2" input.ini > run.log 2>&1); then
      echo "ERROR: run failed, see $1/run.log" >&2
      exit 2
   fi
   shopt -s nullglob
   local h5=("$1/$ob"-*-proc*.h5)
   shopt -u nullglob
   "$VENV_PY" "$DIGEST_PY" write "$1/digest.txt" "${h5[@]}" --case-dir "$CASE_DIR" > /dev/null
}

serial="$CASE_DIR/work-omp-serial"
threaded="$CASE_DIR/work-omp-threads"
run "$serial" exe/adam_flume_cpu 1
run "$threaded" exe/adam_flume_cpu_omp "$THREADS"

failed=0
for f in "$ob-residuals.dat" "$ob-conservation_history.dat" digest.txt; do
   if cmp -s "$serial/$f" "$threaded/$f"; then
      echo "   $f: identical"
   else
      failed=1
      line="$(cmp "$serial/$f" "$threaded/$f" | sed -n -E 's/.*line ([0-9]+).*/\1/p')"
      echo "   $f: DIFFERS from line $line"
      diff <(sed -n "${line}p" "$serial/$f") <(sed -n "${line}p" "$threaded/$f") | sed 's/^/      /' || true
   fi
done
if [[ $failed -ne 0 ]]; then
   echo "FAIL [omp-bitwise] OpenMP ($NP x $THREADS threads) differs from serial: $serial and $threaded are the reproducer"
   exit 1
fi
rm -f "$serial"/*.h5 "$threaded"/*.h5
echo "PASS [omp-bitwise] OpenMP ($NP x $THREADS threads) bitwise identical to serial"
