#!/usr/bin/env bash
# PRISM verification of the [slices] output (slices_restauration.md, F1/F5).
#
# Why: the slices are outputs no run reads back, so a defect in them is silent. slices.ini is the source-free gaussian
# pulse of rmf-amr-fd-pulse (2:1 AMR seam, 36 blocks) with curl and divergence saved and three slices (a trilinear
# line through the seam, an inverse-distance plane, a trilinear box around the block corner at the origin). For every
# process count:
#   1. slices_oracle.py recomputes every slice value of q, curl and div from the XH5F checkpoint of the slice step
#      (saved with ghost cells) and checks the intra-realm ghost cells of every field set against the same-level
#      neighbour interior (the auxiliary fields are computed on the interior cells only; save_slices refreshes them);
#   2. the slice files are compared, bitwise, with those of the first process count: each point is written by the one
#      rank owning it, so the output must not depend on the decomposition.
#
# Usage: ./check.sh [--build] [--np "1 2 4"]
#
# PRISM_EXE overrides the executable under test, e.g. PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
NPS="1 2 4"
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NPS="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np \"N ...\")" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "check.sh: executable $EXE not found" >&2 ; exit 2 ; }

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c "import h5py, numpy" 2>/dev/null; then
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet h5py numpy
fi

TAG="$(basename "$EXE")"
first=""
status=0
for np in $NPS; do
   work="$CASE_DIR/work-$TAG-np$np"
   rm -rf "$work" ; mkdir -p "$work"
   cp "$CASE_DIR/slices.ini" "$work/"
   echo ">> $TAG, np $np"
   if ! (cd "$work" && mpirun -np "$np" "$EXE" slices.ini > log.txt 2>&1); then
      echo "   run failed, see $work/log.txt" ; status=1 ; continue
   fi
   "$VENV_PY" -I "$CASE_DIR/slices_oracle.py" "$work" || status=1
   if [[ -z "$first" ]]; then
      first="$work"
   else
      for mat in "$first"/*.mat; do
         if cmp -s "$mat" "$work/$(basename "$mat")"; then
            echo "   $(basename "$mat"): identical to np $(basename "$first" | sed 's/.*-np//')  PASS (bitwise)"
         else
            echo "   $(basename "$mat"): differs from np $(basename "$first" | sed 's/.*-np//')  FAIL" ; status=1
         fi
      done
   fi
done

# adimensional model (F3): the [grid] and [slice_N] extents are dimensional, the grid divides them by L0 and so must
# the slices (before, the points stayed dimensional: in the wrong place, or outside the domain and an abort at the first
# save); the oracle expects the extents divided by L0. Run on the last process count.
last_np="${NPS##* }"
if [[ -n "$last_np" ]]; then
   work="$CASE_DIR/work-$TAG-adim-np$last_np"
   rm -rf "$work" ; mkdir -p "$work"
   # a plane wave along y: unlike the gaussian pulse, its initial conditions are adimensionalized (B, D and lambda)
   sed -e 's/^physical_model = EM$/physical_model = AdimEM\nL0             = 0.01\nB0             = 10.0/' \
       -e 's/^type            = gaussian_pulse$/type            = plane_wave\nkx              = 0.0\nky              = 1.0\nkz              = 0.0\nlambda          = 0.04/' \
       "$CASE_DIR/slices.ini" > "$work/slices.ini"
   echo ">> $TAG, adimensional model (AdimEM, L0 = 0.01, plane wave), np $last_np"
   if (cd "$work" && mpirun -np "$last_np" "$EXE" slices.ini > log.txt 2>&1); then
      "$VENV_PY" -I "$CASE_DIR/slices_oracle.py" "$work" || status=1
   else
      echo "   run failed, see $work/log.txt" ; status=1
   fi
fi

# input validation (F2): every [slice_N] option is required and checked at init, an invalid one stops the run with a
# message naming it (before, a missing option left the slice undefined and n_save = 0 crashed on mod(it, 0)); an itype
# with a trailing CR (INI authored on Windows) is accepted
reject() { # reject <label> <expected message> <sed expression on slices.ini>
   local work="$CASE_DIR/work-$TAG-input-$1"
   rm -rf "$work" ; mkdir -p "$work"
   sed -e "$3" "$CASE_DIR/slices.ini" > "$work/slices.ini"
   (cd "$work" && mpirun -np 1 "$EXE" slices.ini > log.txt 2>&1) || true
   if grep -qF "$2" "$work/log.txt" && ! ls "$work"/*.mat >/dev/null 2>&1; then
      echo "   input $1: stops with \"$2\"  PASS"
   else
      echo "   input $1: expected a stop with \"$2\", see $work/log.txt  FAIL" ; status=1
   fi
}
echo ">> $TAG, input validation, np 1"
reject missing-n_save "failed to load [slice_1].(n_save)" '/^\[slice_1\]/,/^\[slice_2\]/{/^n_save/d}'
reject n_save-zero    "[slice_1].(n_save) = 0 must be >= 1" '/^\[slice_1\]/,/^\[slice_2\]/s/^n_save = .*/n_save = 0/'
reject itype-unknown  'unknown [slice_1].(itype) "bilinear"' '/^\[slice_1\]/,/^\[slice_2\]/s/^itype  = .*/itype  = bilinear/'
reject nj-zero        "[slice_1].(ni, nj, nk)" '/^\[slice_1\]/,/^\[slice_2\]/s/^nj     = .*/nj     = 0/'
reject emin-above-max "[slice_1] emin" '/^\[slice_1\]/,/^\[slice_2\]/s/^emin_y = .*/emin_y = 0.002/'
reject number-negative "[slices].(slices_number) = -1 must be >= 0" 's/^slices_number = .*/slices_number = -1/'
reject outside-domain "outside the realm domain" '/^\[slice_1\]/,/^\[slice_2\]/s/^emax_x = .*/emax_x = 0.07/'
work="$CASE_DIR/work-$TAG-input-itype-crlf"
rm -rf "$work" ; mkdir -p "$work"
sed -e '/^\[slice_1\]/,/^\[slice_2\]/s/^itype  = trilinear$/itype  = trilinear\r/' "$CASE_DIR/slices.ini" > "$work/slices.ini"
if (cd "$work" && mpirun -np 1 "$EXE" slices.ini > log.txt 2>&1) && [[ -f "$work/slices-slice_01-000000005.mat" ]]; then
   echo "   input itype-crlf: trailing CR accepted, slice saved  PASS"
else
   echo "   input itype-crlf: run failed or no slice, see $work/log.txt  FAIL" ; status=1
fi

if [[ $status -eq 0 ]]; then echo "slices verification PASSED ($TAG)"; else echo "slices verification FAILED ($TAG)"; fi
exit $status
