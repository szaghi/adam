#!/usr/bin/env bash
# FLUME #68 Phase B: the 2-D NumPy prototype of the 2:1 transfer operators and reflux (proto2d.py). No ADAM runs.
#
# B1: linear advection across a 2:1 box, 5th-order linear upwind flux, arms restriction (mean | lag6) x prolongation
# (cubic | quintic) x reflux (on | off), N = 32 / 64 / 128 (finest pair asserted; --long runs 64 / 128 / 256).
# Asserted, lower bounds only (a transfer that recovers order passes): the uniform control at order 4.8,
# lag6-quintic-off (order 5 recovered with the tangential direction and the corners present) at 4.7 in L1 and Linf,
# lag6-cubic-off (the cubic prolongation's order 4) at 3.7. The caps of the mean restriction and of the reflux (Linf
# 1.9) are reported. Cost: about 1 min (--long: about 4 min).
#
# Usage: ./check.sh [--long]
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
CELLS=(32 64 128)
while [[ $# -gt 0 ]]; do
   case "$1" in
      --long) CELLS=(64 128 256) ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --long)" >&2 ; exit 2 ;;
   esac
done

if "$VENV_PY" "$CASE_DIR/proto2d.py" b1 --cells "${CELLS[@]}" \
      --assert uniform:linf:4.8 --assert lag6-quintic-off:l1:4.7 --assert lag6-quintic-off:linf:4.7 \
      --assert lag6-cubic-off:linf:3.7; then
   echo "PROTO PASSED"
else
   echo "PROTO FAILED"
   exit 1
fi
