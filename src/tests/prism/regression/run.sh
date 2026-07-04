#!/usr/bin/env bash
# PRISM regression harness — runs every regression case under one backend
# (cpu or fnl) and compares the produced outputs against committed golden
# references.
#
# Two references per case (raw HDF5 checkpoints are ~100 MB/rank and are NOT
# committed — see README "Reference data"):
#   golden/<backend>/digest.txt        — per-variable field digest (digest.py)
#   golden/<backend>/*-residuals.dat   — per-iteration residuals log (byte-exact)
#
# Usage:
#   ./run.sh cpu                       # build prism-cpu-gnu (default varset), run, diff
#   ./run.sh fnl --varset local_nvf    # build prism-fnl-nvf with that varset, run, diff
#   ./run.sh cpu --no-build            # skip the build step (use existing exe/)
#
# Exits 0 on full pass, non-zero on any case failure.
#
# HDF5 paths and NVF_CC are fobos *variables* resolved from the active varset
# (see the repo `fobos` [varset:*] sections) — they are NOT shell environment
# variables and this harness does not set them. The CPU backend uses fobos's
# default varset; the FNL backend needs `--varset local_nvf`, which the
# run-fnl-local.sh wrapper supplies.
#
# A private Python venv (exe/.regression-venv/, gitignored) is created on
# first run to provide h5py for digest.py.

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
      *) echo "ERROR: unknown argument '$1'" >&2; exit 2 ;;
   esac
   shift
done

case "$BACKEND" in
   cpu)
      MODE="prism-cpu-gnu"
      EXE="exe/adam_prism_cpu"
      ;;
   fnl)
      MODE="prism-fnl-nvf"
      EXE="exe/adam_prism_fnl"
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

# GPU race-shape lint gate (issue #26 G1.d): fail the sweep BEFORE any run if a
# forbidden device-code pattern reappears — the gate fires at the first firewall
# run, where FNL development actually happens (no CI GPU runners).
"$REGRESSION_DIR/../../lint/check-gpu-race-shapes.sh"

if ! command -v mpirun >/dev/null 2>&1; then
   echo "ERROR: mpirun not on PATH" >&2
   exit 2
fi

# ---------------------------------------------------------------------------
# Digest tooling
# ---------------------------------------------------------------------------
# Regression goldens are field *digests* (digest.py), not raw HDF5 — a full
# checkpoint is ~100 MB/rank and cannot be committed. digest.py needs h5py;
# WSL2 / CI system Python is externally managed, so use a private venv that
# is created on first run and reused thereafter.
DIGEST_PY="$REGRESSION_DIR/digest.py"
if [[ ! -f "$DIGEST_PY" ]]; then
   echo "ERROR: digest.py not found at $DIGEST_PY" >&2
   exit 2
fi

# The venv lives under exe/ — NOT under src/ — on purpose: every tool that
# scans src/ (fobis for the build, formal for the API docs) would otherwise
# pick up the Fortran test fixtures shipped inside numpy/h5py and either
# break the link step or generate junk doc pages. exe/ is already gitignored
# and ignored by those tools, so one location keeps the venv invisible to
# all of them.
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
   # Skip non-case directories (golden refs aren't cases)
   [[ "$case_name" == "golden" ]] && continue
   [[ -f "$case_dir/input.ini" ]] || { skip_count=$((skip_count + 1)); continue; }

   # A case with no golden/<backend>/ is not a regression anchor yet: it has
   # nothing to diff against, and running it can hang/OOM the whole sweep (e.g.
   # an in-progress N>1 case), taking down committed cases that sort after it.
   # So SKIP goldenless cases by default. To run one anyway — the initial
   # golden-capture workflow — set REGRESSION_RUN_GOLDENLESS=1, which runs the
   # case and produces work-<backend>/digest.txt for promotion into golden/.
   golden_dir="${case_dir%/}/golden/$BACKEND"
   if [[ ! -d "$golden_dir" ]]; then
      if [[ "${REGRESSION_RUN_GOLDENLESS:-0}" == "1" ]]; then
         echo "!! [$case_name/$BACKEND] no golden at $golden_dir — running anyway"
         echo "   (REGRESSION_RUN_GOLDENLESS=1; initial golden-capture mode)"
      else
         echo ">> [$case_name/$BACKEND] no golden at $golden_dir — skipping case"
         echo "   (set REGRESSION_RUN_GOLDENLESS=1 to run it for initial golden capture)"
         skip_count=$((skip_count + 1))
         continue
      fi
   fi

   # The checkpoint file prefix(es). PRISM also writes a restart dump
   # (restart_basename) which is NOT a regression checkpoint and must be
   # excluded from the digest — filter strictly on output_basename.
   #
   # Multi-realm (Phase D) regression cases use a forest manifest at
   # input.ini that points at per-realm INIs (e.g. realm_1.ini, realm_2.ini),
   # each with its own output_basename. Collect every output_basename from
   # the case directory's .ini files so the produced-h5 glob covers all
   # realms.
   output_basenames=()
   while IFS= read -r ob; do
      [[ -n "$ob" ]] && output_basenames+=("$ob")
   done < <(sed -n -E 's/^[[:space:]]*output_basename[[:space:]]*=[[:space:]]*([^[:space:];]+).*/\1/p' \
                      "$case_dir"/*.ini)
   if [[ ${#output_basenames[@]} -eq 0 ]]; then
      echo "FAIL [$case_name/$BACKEND] could not parse output_basename from any .ini in $case_dir"
      fail_count=$((fail_count + 1))
      failed_cases+=("$case_name")
      continue
   fi

   workdir="$case_dir/work-$BACKEND"
   rm -rf "$workdir"
   mkdir -p "$workdir"
   # Copy every .ini from the case directory (input.ini plus any per-realm
   # INIs the manifest references). The PRISM driver auto-detects whether
   # input.ini is a manifest or a single-realm config.
   cp "$case_dir"/*.ini "$workdir/"

   echo
   echo "============================================================"
   echo "== [$case_name/$BACKEND] running"
   echo "============================================================"
   pushd "$workdir" >/dev/null
   t0=$(date +%s)
   mpirun -np 2 "$EXE_ABS" 2>&1 | tail -20
   t1=$(date +%s)
   popd >/dev/null
   echo ">> [$case_name/$BACKEND] runtime: $((t1 - t0))s"

   # ----- Compare outputs against golden -----
   # Two committed references per case (see README "Reference data"):
   #   golden/<backend>/digest.txt          — field digest of every HDF5 checkpoint
   #   golden/<backend>/*-residuals.dat     — per-iteration residuals log (byte-exact)
   case_failed=0

   # Always (re)compute the produced digest, even with no golden — during
   # initial setup this is what gets copied into golden/. Only the
   # output_basename checkpoints are digested; the restart dump is excluded.
   # For multi-realm cases we glob across every realm's output_basename.
   shopt -s nullglob
   produced_h5=()
   for ob in "${output_basenames[@]}"; do
      for h5 in "$workdir/$ob"-*.h5; do
         produced_h5+=("$h5")
      done
   done
   shopt -u nullglob
   if [[ ${#produced_h5[@]} -eq 0 ]]; then
      echo "FAIL [$case_name/$BACKEND] no '*-*.h5' checkpoints produced for any of: ${output_basenames[*]}"
      case_failed=1
   else
      if ! "$VENV_PY" "$DIGEST_PY" write "$workdir/digest.txt" "${produced_h5[@]}" --case-dir "$case_dir"; then
         echo "FAIL [$case_name/$BACKEND] digest computation failed"
         case_failed=1
      fi
   fi

   if [[ -d "$golden_dir" ]]; then
      # Field digest comparison (tolerance-aware; see digest.py).
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
      # Residuals log (line-by-line text diff).
      shopt -s nullglob
      for residuals_golden in "$golden_dir"/*-residuals.dat; do
         fname="$(basename "$residuals_golden")"
         produced="$workdir/$fname"
         if [[ ! -f "$produced" ]]; then
            echo "FAIL [$case_name/$BACKEND] missing residuals: $fname"
            case_failed=1
            continue
         fi
         # Tolerance-based residuals comparison: bit-exact diff is too strict
         # for cross-runner CI (different x86 microarchitectures, MPI rank
         # allocations, libm versions produce last-digit roundoff). Uses the
         # same (rtol, atol) as digest.py compare — calibrated for cross-
         # compiler runs. See digest.py header for tolerance rationale.
         if ! "$VENV_PY" "$DIGEST_PY" compare-residuals "$produced" "$residuals_golden"; then
            echo "FAIL [$case_name/$BACKEND] residuals differ: $fname"
            case_failed=1
         fi
      done
      shopt -u nullglob
   fi

   # --- β cross-config invariant (issue #18) ---
   # The rmf-2realm-stagesync case exercises β (stage-coincident seam fill)
   # on a geometry whose union is bit-aligned with the single-realm rmf
   # case. β's design oracle is: under same-K + same-physics + same-ODE +
   # different/identical spatial operators, the multi-realm digest must
   # match the monolithic single-realm digest within the usual tolerance
   # (per-block metadata `dxdydz`, `origin`, `time_iteration` differ by
   # block-count scaling and are auto-downgraded to SKIP_METADATA by
   # digest.py). Enforced here as a continuously-checked invariant
   # (hardcoded for this single case; if a second case ever needs cross-
   # config equivalence, generalise to a per-case `equiv_to` annotation).
   if [[ "$case_name" == "rmf-2realm-stagesync" && -f "$workdir/digest.txt" ]]; then
      rmf_golden="${REGRESSION_DIR}/rmf/golden/${BACKEND}/digest.txt"
      if [[ ! -f "$rmf_golden" ]]; then
         echo "FAIL [$case_name/$BACKEND] β cross-config oracle: missing single-realm reference $rmf_golden"
         case_failed=1
      elif ! "$VENV_PY" "$DIGEST_PY" compare "$workdir/digest.txt" "$rmf_golden"; then
         echo "FAIL [$case_name/$BACKEND] β cross-config oracle: digest does not match single-realm rmf golden (issue #18)"
         case_failed=1
      else
         echo ">> [$case_name/$BACKEND] β cross-config oracle: matches single-realm rmf golden"
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
