#!/usr/bin/env bash
# Removed-input-key lint gate (issue #58).
#
# Static check of the committed PRISM inputs for keys the parser no longer reads. A removed key is not an error by
# itself (FiNeR ignores unknown keys), but its replacement is required, so an input that still carries only the old
# key error_stops at load, far from the commit that removed it. Rules:
#
#   [fWLayer] C  -> removed by 8e05d363; replaced by the physical [fWLayer] width (required even when no layer is
#                   used). 48 research inputs kept C until issue #58.
#   [PIC] elliptic_correction -> removed by issue #84 E3d (the t0 solve on FLAIL leaves div D - rho at round-off up to
#                   the checkerboard component no D can carry); .true. stops the run, so no replacement. 29 inputs
#                   carried it (all .false.) until E3d.
#
# Section-aware (a C key in another section is not the fWLayer one) and CRLF-tolerant. Git-ignored files (the
# harnesses' work-* copies) are skipped.
#
# Invoked from the top of the PRISM regression sweep (run.sh); standalone use:
#   ./src/tests/lint/check-removed-input-keys.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

hits="$(cd "$REPO_ROOT" && git ls-files -z 'src/tests/prism/*.ini' | xargs -0 awk '
   FNR == 1 { section = "" }
   { line = $0; sub(/\r$/, "", line) }
   line ~ /^[[:space:]]*\[/ { section = line; gsub(/[[:space:]]/, "", section); next }
   section == "[fWLayer]" && line ~ /^[[:space:]]*C[[:space:]]*=/ { print FILENAME ":" FNR ": [fWLayer] C" }
   section == "[PIC]" && line ~ /^[[:space:]]*elliptic_correction[[:space:]]*=/ {
      print FILENAME ":" FNR ": [PIC] elliptic_correction" }
')"
if [[ -n "$hits" ]]; then
   echo "FAIL [lint/removed-input-keys] removed keys in committed PRISM inputs (replace [fWLayer] C with width; delete"
   echo "     [PIC] elliptic_correction):"
   echo "$hits"
   exit 1
fi
echo "PASS [lint/removed-input-keys] no removed keys in the committed PRISM inputs"
