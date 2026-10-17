#!/usr/bin/env bash
# FLUME #68 Phase B: the 2-D NumPy prototype of the 2:1 transfer operators and reflux (proto2d.py). No ADAM runs.
#
# B1: linear advection across a 2:1 box, 5th-order linear upwind flux, arms restriction (mean | lag6) x prolongation
# (cubic | quintic) x reflux (on | off), N = 32 / 64 / 128 (finest pair asserted; --long runs 64 / 128 / 256).
# Asserted, lower bounds only (a transfer that recovers order passes): the uniform control at order 4.8,
# lag6-quintic-off (order 5 recovered with the tangential direction and the corners present) at 4.7 in L1 and Linf,
# lag6-cubic-off (the cubic prolongation's order 4) at 3.7. The caps of the mean restriction and of the reflux (Linf
# 1.9) are reported. Cost: about 1 min (--long: about 4 min).
# B2: WENO5-JS with the reflux policies always | troubled | never (proto2d.py b2 --check): the smooth order of troubled
# at least 4.7; on the Burgers shocks crossing coarse to fine and fine to coarse, troubled within 0.05 cells of always
# and its mass defect at most 1e-3 of the pulse, never's at least 1e-3 (the negative control). Shocks at N = 64
# (--long: 64 and 128). Cost: about 3 min (--long: about 5 min).
# B3: the 6-point restriction footprints against the fine block size (proto2d.py b3, reflux off): shift-inward at 6
# fine cells per block and centred two-pass, each with the cubic (order 3.7) and quintic (order 4.7) prolongation;
# shift-inward at 4 fine cells is reported infeasible. N = 24 / 48 / 96 (--long: 48 / 96 / 192). About 20 s.
#
# Usage: ./check.sh [--legs b1,b2,b3|all] [--long]   (default: b1,b2,b3; a leg runs only when the ones before it pass)
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
source "$VERIF_DIR/legs.sh"
legs_declare "b1 b2 b3" "b1 b2 b3"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
CELLS=(32 64 128)
SHOCK_CELLS=(64)
B3_CELLS=(24 48 96)
while [[ $# -gt 0 ]]; do
   case "$1" in
      --legs) legs_select "$2" ; shift 2 ;;
      --long) CELLS=(64 128 256) ; SHOCK_CELLS=(64 128) ; B3_CELLS=(48 96 192) ; shift ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --legs L, --long)" >&2 ; exit 2 ;;
   esac
done
legs_resolve check.sh

FAILED=0
if has_leg b1; then
   "$VENV_PY" "$CASE_DIR/proto2d.py" b1 --cells "${CELLS[@]}" \
      --assert uniform:linf:4.8 --assert lag6-quintic-off:l1:4.7 --assert lag6-quintic-off:linf:4.7 \
      --assert lag6-cubic-off:linf:3.7 || FAILED=1
fi
if [[ $FAILED -eq 0 ]] && has_leg b2; then
   "$VENV_PY" "$CASE_DIR/proto2d.py" b2 --cells 32 64 128 --shock-cells "${SHOCK_CELLS[@]}" --check || FAILED=1
fi
if [[ $FAILED -eq 0 ]] && has_leg b3; then
   "$VENV_PY" "$CASE_DIR/proto2d.py" b3 --cells "${B3_CELLS[@]}" --blocks 4 6 8 0 \
      --assert lag6-b6-quintic:linf:4.7 --assert lag6c-quintic:linf:4.7 --assert lag6-b6-cubic:linf:3.7 \
      --assert lag6c-cubic:linf:3.7 || FAILED=1
fi
if [[ $FAILED -eq 0 ]]; then
   echo "PROTO PASSED (legs: $LEGS)"
else
   echo "PROTO FAILED (legs: $LEGS)"
   exit 1
fi
