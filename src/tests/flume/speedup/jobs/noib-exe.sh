#!/bin/bash
# A run.sh --exe wrapper that runs the sphere without its immersed body (#133): every rank writes the case input without
# the [solids] and [solid_N] sections next to it (an atomic rename, the same bytes from every rank) and execs the real
# executable on it. The grid, the dump and the flow are the leg's; only the body, its eikonal extension and its 9 extra
# ghost exchanges per step go.   ./run.sh windows ... --exe "$PWD/jobs/noib-exe.sh"   with NOIB_EXE=<the executable>
set -euo pipefail
: "${NOIB_EXE:?noib-exe.sh: set NOIB_EXE to the FLUME executable}"
in=${1:?noib-exe.sh: the input file}
tmp="$in.noib.$$.$RANDOM"
awk '/^\[/ {skip = ($0 ~ /^\[solids\]/ || $0 ~ /^\[solid_[0-9]+\]/)} !skip' "$in" > "$tmp"
mv -f "$tmp" "$in.noib"
exec "$NOIB_EXE" "$in.noib" "${@:2}"
