#!/usr/bin/env bash
# FLUME regression harness (issue #35, P8): runs every regression case under one backend (cpu or fnl) and compares the
# produced outputs against committed golden references. A copy of the PRISM harness (src/tests/prism/regression/run.sh)
# that reuses its digest.py unchanged.
#
# Three references per case (raw HDF5 checkpoints are NOT committed):
#   golden/<backend>/digest.txt                     per-variable field digest of every checkpoint (digest.py)
#   golden/<backend>/*-residuals.dat                per-iteration residuals log
#   golden/<backend>/*-conservation_history.dat     per-iteration volume integrals of the conservative variables
# The two logs are compared with digest.py compare-residuals (header verbatim, integer columns exact, float columns
# within the digest tolerances).
#
# Usage:
#   ./run.sh cpu                       # build flume-cpu-gnu (default varset), run, diff
#   ./run.sh fnl --varset local_nvf    # build flume-fnl-nvf with that varset, run, diff (see run-fnl-local.sh)
#   ./run.sh cpu --no-build            # skip the build step (use existing exe/)
#
# Exits 0 on full pass, non-zero on any case failure.
#
# Cases (the input.ini of each is a frozen copy of a verification input; the verification oracles live in
# src/tests/flume/verification and assert the physics, this suite only detects changes):
#   sod-x, sod-y, sod-z      V1 Sod shock tube along each direction, WENO-5 characteristic, SSP-33, fast path
#   vortex-periodic          V2 isentropic vortex, periodic quadtree, SSP-54, 100 steps
#   amr-periodic-reflux      V3 init-time AMR, 2:1 faces on every side of the refined octant, reflux, staged path
#   shock-cylinder-ib        V6 Mach 2 shock over a cylinder, immersed boundary, solid AMR marker, 120 steps
#   sod-2realm               sod-x split in two realms at the diaphragm, mirror seam, beta cadence (issue #37); a forest
#                            manifest (input.ini) plus one INI per realm
#
# A private Python venv (exe/.regression-venv/, gitignored) is created on first run to provide h5py for digest.py.

set -euo pipefail

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
BACKEND="${1:-}"
NO_BUILD=0
VARSET=""
shift || true
while [[ $# -gt 0 ]]; do
   case "$1" in
      --no-build) NO_BUILD=1 ;;
      --varset)
         shift || { echo "ERROR: --varset needs a value" >&2; exit 2; }
         VARSET="$1"
         ;;
      *) echo "ERROR: unknown argument '$1' (accepted: --no-build, --varset <name>)" >&2; exit 2 ;;
   esac
   shift
done

case "$BACKEND" in
   cpu)
      MODE="flume-cpu-gnu"
      EXE="exe/adam_flume_cpu"
      ;;
   fnl)
      MODE="flume-fnl-nvf"
      EXE="exe/adam_flume_fnl"
      ;;
   *)
      echo "Usage: $0 {cpu|fnl} [--no-build] [--varset <name>]" >&2
      exit 2
      ;;
esac

# ---------------------------------------------------------------------------
# Paths and environment checks
# ---------------------------------------------------------------------------
REGRESSION_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$REGRESSION_DIR/../../../.." && pwd)" || {
   echo "ERROR: cannot resolve repo root from $REGRESSION_DIR" >&2
   exit 2
}

# GPU race-shape lint gate (issue #26 G1.d): the check scans all of src/, FLUME device code included.
"$REPO_ROOT/src/tests/lint/check-gpu-race-shapes.sh"

if ! command -v mpirun >/dev/null 2>&1; then
   echo "ERROR: mpirun not on PATH" >&2
   exit 2
fi

# ---------------------------------------------------------------------------
# Digest tooling: the PRISM digest.py, reused unchanged (FLUME checkpoints share the block_*-proc*-<var> dataset layout
# and the <basename>-<step>-proc<rank>.h5 file names). The venv lives under exe/, not src/, so that fobis and formal
# never scan the Fortran fixtures shipped inside numpy/h5py (see the PRISM harness).
# ---------------------------------------------------------------------------
DIGEST_PY="$REPO_ROOT/src/tests/prism/regression/digest.py"
if [[ ! -f "$DIGEST_PY" ]]; then
   echo "ERROR: digest.py not found at $DIGEST_PY" >&2
   exit 2
fi
VENV_DIR="$REPO_ROOT/exe/.regression-venv"
VENV_PY="$VENV_DIR/bin/python"
if [[ ! -x "$VENV_PY" ]]; then
   echo ">> creating digest venv at $VENV_DIR"
   python3 -m venv "$VENV_DIR"
   "$VENV_DIR/bin/pip" install --quiet --upgrade pip
   "$VENV_DIR/bin/pip" install --quiet h5py numpy
fi
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> repairing digest venv (h5py/numpy import failed)"
   "$VENV_DIR/bin/pip" install --quiet h5py numpy
fi

# ---------------------------------------------------------------------------
# Build (skippable)
# ---------------------------------------------------------------------------
cd "$REPO_ROOT"

if [[ $NO_BUILD -eq 0 ]]; then
   build_args=(build --mode "$MODE")
   if [[ -n "$VARSET" ]]; then
      build_args+=(--varset "$VARSET")
   fi
   echo ">> building: fobis ${build_args[*]}"
   fobis "${build_args[@]}"
fi

if [[ ! -x "$EXE" ]]; then
   echo "ERROR: executable $EXE not found after build" >&2
   exit 2
fi
EXE_ABS="$REPO_ROOT/$EXE"

# ---------------------------------------------------------------------------
# Run every regression case
# ---------------------------------------------------------------------------
fail_count=0
pass_count=0
skip_count=0
declare -a failed_cases=()

for case_dir in "$REGRESSION_DIR"/*/; do
   case_name="$(basename "$case_dir")"
   [[ -f "$case_dir/input.ini" ]] || { skip_count=$((skip_count + 1)); continue; }

   # A case with no golden/<backend>/ is not a regression anchor yet: skip it, unless REGRESSION_RUN_GOLDENLESS=1 (the
   # initial golden-capture workflow, which produces work-<backend>/ for promotion into golden/).
   golden_dir="${case_dir%/}/golden/$BACKEND"
   if [[ ! -d "$golden_dir" ]]; then
      if [[ "${REGRESSION_RUN_GOLDENLESS:-0}" == "1" ]]; then
         echo "!! [$case_name/$BACKEND] no golden at $golden_dir — running anyway (REGRESSION_RUN_GOLDENLESS=1)"
      else
         echo ">> [$case_name/$BACKEND] no golden at $golden_dir — skipping case"
         echo "   (set REGRESSION_RUN_GOLDENLESS=1 to run it for initial golden capture)"
         skip_count=$((skip_count + 1))
         continue
      fi
   fi

   # A multi-realm case is a forest manifest (input.ini) plus one INI per realm, each with its own output_basename:
   # collect them all, copy every INI (issue #37).
   output_basenames=()
   while IFS= read -r ob; do
      [[ -n "$ob" ]] && output_basenames+=("$ob")
   done < <(sed -n -E 's/^[[:space:]]*output_basename[[:space:]]*=[[:space:]]*([^[:space:];]+).*/\1/p' "$case_dir"/*.ini)
   if [[ ${#output_basenames[@]} -eq 0 ]]; then
      echo "FAIL [$case_name/$BACKEND] could not parse output_basename from any .ini in $case_dir"
      fail_count=$((fail_count + 1))
      failed_cases+=("$case_name")
      continue
   fi

   workdir="$case_dir/work-$BACKEND"
   rm -rf "$workdir"
   mkdir -p "$workdir"
   cp "$case_dir"/*.ini "$workdir/"

   echo
   echo "============================================================"
   echo "== [$case_name/$BACKEND] running"
   echo "============================================================"
   # Run the case, tolerant of a hung MPI/GPU finalize (the same two bounds as the PRISM harness): the end-of-run marker
   # detector (REGRESSION_TEARDOWN_DETECT, default 1; grace REGRESSION_TEARDOWN_GRACE, default 15 s) terminates a
   # process still alive after the driver printed "ADAM run complete", and REGRESSION_MPIRUN_TIMEOUT (default 300 s on
   # cpu, 900 s on fnl, whose shock-cylinder-ib runs ~300 s on the WSL box; 0 disables) is the hard backstop. The
   # digest/residual comparison below is the pass/fail gate.
   default_timeout=300
   [[ "$BACKEND" == "fnl" ]] && default_timeout=900
   mpirun_timeout="${REGRESSION_MPIRUN_TIMEOUT:-$default_timeout}"
   teardown_detect="${REGRESSION_TEARDOWN_DETECT:-1}"
   teardown_grace="${REGRESSION_TEARDOWN_GRACE:-15}"
   mpirun_cmd=(mpirun -np 2 "$EXE_ABS" input.ini)
   if [[ "$mpirun_timeout" != "0" ]]; then
      mpirun_cmd=(timeout --signal=TERM --kill-after=30s "$mpirun_timeout" "${mpirun_cmd[@]}")
   fi
   pushd "$workdir" >/dev/null
   t0=$(date +%s)
   : >run.log
   set +e
   "${mpirun_cmd[@]}" >run.log 2>&1 &
   run_pid=$!
   tail -n +1 -f --pid="$run_pid" run.log 2>/dev/null &
   tail_pid=$!
   detector_pid=""
   if [[ "$teardown_detect" != "0" ]]; then
      (
         while kill -0 "$run_pid" 2>/dev/null; do
            sleep 3
            if grep -q 'ADAM run complete' run.log 2>/dev/null; then
               sleep "$teardown_grace"
               if kill -0 "$run_pid" 2>/dev/null; then
                  echo "[detector] end-of-run marker seen; still alive after ${teardown_grace}s — terminating hung MPI_Finalize" >>run.log
                  kill -TERM "$run_pid" 2>/dev/null
                  ( sleep 30 ; kill -KILL "$run_pid" 2>/dev/null ) &
               fi
               exit 0
            fi
         done
      ) &
      detector_pid=$!
   fi
   wait "$run_pid"
   mpirun_rc=$?
   [[ -n "$detector_pid" ]] && { kill "$detector_pid" 2>/dev/null ; wait "$detector_pid" 2>/dev/null ; }
   kill "$tail_pid" 2>/dev/null ; wait "$tail_pid" 2>/dev/null
   set -e
   t1=$(date +%s)
   detector_fired=0
   grep -q '^\[detector\]' run.log 2>/dev/null && detector_fired=1
   popd >/dev/null
   if [[ "$detector_fired" -eq 1 ]]; then
      echo ">> [$case_name/$BACKEND] run complete; hung MPI_Finalize terminated after ${teardown_grace}s — continuing to digest"
   elif [[ "$mpirun_rc" -eq 124 ]]; then
      echo ">> [$case_name/$BACKEND] mpirun hit the ${mpirun_timeout}s hard timeout — continuing to digest produced output"
   elif [[ "$mpirun_rc" -ne 0 ]]; then
      echo ">> [$case_name/$BACKEND] mpirun exited non-zero ($mpirun_rc) — continuing to digest (checkpoints decide pass/fail)"
   fi
   echo ">> [$case_name/$BACKEND] runtime: $((t1 - t0))s"

   # ----- Compare outputs against golden -----
   case_failed=0

   # Digest only <output_basename>-<step>-proc<rank>.h5: the filter is structural, so a restart dump never enters. The
   # checkpoints of every realm at the same step aggregate into one digest row (digest.py keys rows on the step).
   shopt -s nullglob
   produced_h5=()
   for ob in "${output_basenames[@]}"; do
      for h5 in "$workdir/$ob"-*.h5; do
         [[ "$(basename "$h5")" =~ -[0-9]+-proc[0-9]+\.h5$ ]] || continue
         produced_h5+=("$h5")
      done
   done
   shopt -u nullglob
   if [[ ${#produced_h5[@]} -eq 0 ]]; then
      echo "FAIL [$case_name/$BACKEND] no '<output_basename>-<step>-proc<rank>.h5' checkpoints for: ${output_basenames[*]}"
      case_failed=1
   elif ! "$VENV_PY" "$DIGEST_PY" write "$workdir/digest.txt" "${produced_h5[@]}" --case-dir "$case_dir"; then
      echo "FAIL [$case_name/$BACKEND] digest computation failed"
      case_failed=1
   fi

   if [[ -d "$golden_dir" ]]; then
      golden_digest="$golden_dir/digest.txt"
      if [[ ! -f "$golden_digest" ]]; then
         echo "FAIL [$case_name/$BACKEND] golden digest missing: $golden_digest"
         case_failed=1
      elif [[ -f "$workdir/digest.txt" ]]; then
         if ! "$VENV_PY" "$DIGEST_PY" compare "$workdir/digest.txt" "$golden_digest"; then
            echo "FAIL [$case_name/$BACKEND] field digest mismatch"
            case_failed=1
         fi
      fi
      shopt -s nullglob
      for log_golden in "$golden_dir"/*-residuals.dat "$golden_dir"/*-conservation_history.dat; do
         fname="$(basename "$log_golden")"
         produced="$workdir/$fname"
         if [[ ! -f "$produced" ]]; then
            echo "FAIL [$case_name/$BACKEND] missing log: $fname"
            case_failed=1
            continue
         fi
         if ! "$VENV_PY" "$DIGEST_PY" compare-residuals "$produced" "$log_golden"; then
            echo "FAIL [$case_name/$BACKEND] log differs: $fname"
            case_failed=1
         fi
      done
      shopt -u nullglob
   fi

   # Cross-configuration oracle (issue #37): a case whose cells are the union of another case's cells (e.g. a
   # multi-realm split of a single-realm case) names that case in `equivalent_to`; its digest must match the other
   # case's golden too (the per-block metadata rows, which count blocks, are skipped by digest.py).
   if [[ -f "$case_dir/equivalent_to" && -f "$workdir/digest.txt" ]]; then
      equivalent="$(tr -d '[:space:]' < "$case_dir/equivalent_to")"
      equivalent_golden="$REGRESSION_DIR/$equivalent/golden/$BACKEND/digest.txt"
      if [[ ! -f "$equivalent_golden" ]]; then
         echo "FAIL [$case_name/$BACKEND] cross-configuration oracle: missing reference $equivalent_golden"
         case_failed=1
      elif ! "$VENV_PY" "$DIGEST_PY" compare "$workdir/digest.txt" "$equivalent_golden"; then
         echo "FAIL [$case_name/$BACKEND] cross-configuration oracle: digest differs from the $equivalent golden"
         case_failed=1
      else
         echo ">> [$case_name/$BACKEND] cross-configuration oracle: matches the $equivalent golden"
      fi
   fi

   if [[ $case_failed -eq 0 ]]; then
      echo "PASS [$case_name/$BACKEND]"
      pass_count=$((pass_count + 1))
   else
      fail_count=$((fail_count + 1))
      failed_cases+=("$case_name")
   fi
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo
echo "============================================================"
echo "== Summary [$BACKEND]"
echo "============================================================"
echo "PASS: $pass_count"
echo "FAIL: $fail_count"
echo "SKIP: $skip_count"
if [[ $fail_count -gt 0 ]]; then
   echo "Failed cases: ${failed_cases[*]}"
   exit 1
fi
exit 0
