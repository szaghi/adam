#!/usr/bin/env bash
# GPU race-shape lint gate (issue #26 G1.d).
#
# Static grep gate for device-code patterns that nvfortran turns into
# UNPRIVATIZED compiler temporaries — threads race on them; benign on uniform
# data, destructive at 2:1 AMR level mixes, invisible to compute-sanitizer
# (in-bounds) and to uniform-grid validation (value-coincident). Full rules:
# CLAUDE-gpu "Device Code Pitfalls" (#22 F1-bis/F2, #26 G1).
#
# Pattern 1: strided sections of the shared device dxyz array passed from
# inside kernels — `dxyz_gpu(b,1:3)`. Exact by construction: `dxyz_gpu` exists
# only in device-kernel scope; the sweep (#26 G1.b/c) hoists every use to a
# private `dxyz_b(3)`. Comment lines (leading `!`) are exempt: the rationale
# comments legitimately name the forbidden pattern.
#
# NOTE (#26 G1.b lesson): this gate catches the dxyz pattern ONLY. The race
# FAMILY is wider — any strided section of a shared device array as a device
# call argument (e.g. the OUT-section `curl_gpu(b,i,j,k,ivar:)` fixed in
# G1.b) — and is not grep-able in general: reviews and -Minfo audits
# ("Copy in/out of X in call to Y", "pgf90_alloc") remain mandatory.
#
# Invoked from the top of the regression sweeps (run.sh); standalone use:
#   ./src/tests/lint/check-gpu-race-shapes.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
fail=0

hits="$(grep -rn 'dxyz_gpu(b,1:3)' "$REPO_ROOT/src" --include='*.F90' | grep -vE ':[0-9]+:[[:space:]]*!' || true)"
if [[ -n "$hits" ]]; then
   echo "FAIL [lint/gpu-race-shapes] forbidden strided device section 'dxyz_gpu(b,1:3)' (hoist to a private dxyz_b(3), see CLAUDE-gpu):"
   echo "$hits"
   fail=1
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [lint/gpu-race-shapes] no forbidden GPU race shapes"
   exit 0
fi
exit 1
