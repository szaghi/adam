#!/usr/bin/env bash
# FLUME verification of the I/O channels (issue #35, P7): V7 restart round trip, slices, auxiliary fields.
#
# Why: restart files, slices and auxiliary fields are outputs no run reads back in normal operation, so a defect in
# them is silent. Two legs:
#   1. V7, restart round trip, on three cases covering the three integration paths: sod-x (fast path),
#      amr-periodic (staged path, AMR, reflux) and shock-cylinder (immersed boundary, AMR). Run A goes N steps;
#      run B goes N/2 steps saving a restart, then restarts and completes N. The final fields must be bitwise
#      identical on the interior cells (edge/corner ghost cells that no map fills keep stale values the directional
#      stencils never read), and the residuals and conservation histories byte-identical (the restarted run appends
#      to them and does not re-save its starting step).
#   2. sod-x with [IO] save_auxiliary_fields and one trilinear slice on the x cell centres: io_oracle.py checks the
#      auxiliary fields against the conservative ones and the slice against the cell values.
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

run() { # run <work-dir> <ini> <log>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > "$3" 2>&1); then
      echo "check.sh: run failed, see $1/$3" >&2
      exit 1
   fi
}

restart_round_trip() { # restart_round_trip <case-ini> <N>
   local ini="$1" n="$2" h=$(( $2 / 2 )) name base
   name="$(basename "$ini")" ; base="${name%.ini}"
   local a="$CASE_DIR/work-$TAG-$base-A" b="$CASE_DIR/work-$TAG-$base-B"
   local common=(-e "s/^it_max   = .*/it_max   = $n/" -e "s/^it_save                = .*/it_save                = 100000/")
   rm -rf "$a" "$b" ; mkdir -p "$a" "$b"
   echo ">> V7 $base: A $n steps, B $h + restart + $h ($(basename "$EXE"), np $NP)"
   sed "${common[@]}" -e "s/^restart_save           = .*/restart_save           = 0/" "$ini" > "$a/$name"
   run "$a" "$name" log.txt
   sed "${common[@]}" -e "s/^it_max   = .*/it_max   = $h/" -e "s/^restart_save           = .*/restart_save           = $h/" \
       "$ini" > "$b/$name"
   run "$b" "$name" log-1.txt
   sed "${common[@]}" -e "s/^restart_save           = .*/restart_save           = 0/" \
       -e "s/^restart                = .*/restart                = .true./" "$ini" > "$b/$name"
   run "$b" "$name" log-2.txt
   "$VENV_PY" "$VERIF_DIR/conservation/conservation_oracle.py" --compare "$a" "$b" --tol 0 --ngc 3
   for hist in residuals conservation_history; do
      if cmp -s "$a/$base-$hist.dat" "$b/$base-$hist.dat"; then
         echo "   $hist history: identical  PASS"
      else
         echo "   $hist history: differs  FAIL" >&2
         exit 1
      fi
   done
}

restart_round_trip "$VERIF_DIR/sod/sod-x.ini" 40
restart_round_trip "$VERIF_DIR/conservation/amr-periodic.ini" 20
restart_round_trip "$VERIF_DIR/shock-cylinder/shock-cylinder.ini" 40

work="$CASE_DIR/work-$TAG-sod-x-io"
rm -rf "$work" ; mkdir -p "$work"
echo ">> slices and auxiliary fields: sod-x, 40 steps"
sed -e "s/^it_max   = .*/it_max   = 40/" -e "s/^save_auxiliary_fields  = .*/save_auxiliary_fields  = .true./" \
    -e "s/^slices_number = 0/slices_number = 1\n\n[slice_1]\nitype  = trilinear\nn_save = 40\nni     = 200\nnj     = 1\nnk     = 1\nemin_x = 0.0\nemin_y = 0.46875\nemin_z = 0.46875\nemax_x = 1.0\nemax_y = 0.53125\nemax_z = 0.53125/" \
    "$VERIF_DIR/sod/sod-x.ini" > "$work/sod-x.ini"
run "$work" sod-x.ini log.txt
"$VENV_PY" "$CASE_DIR/io_oracle.py" "$work"
echo "I/O verification PASSED ($TAG)"
