#!/usr/bin/env bash
# ADAM, check of adam_ini_library%ini_load (issue #80, P4): the merge checks of test_ini_extends, then the two inputs that
# must stop the run (a missing base file, a cycle), each with its message. Usage: ./check.sh (after
# `fobis build --mode test-ini-extends-gnu`); INI_EXE overrides the executable.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO_ROOT="$(cd "$HERE/../../.." && pwd)"
EXE="${INI_EXE:-$REPO_ROOT/exe/test_ini_extends}"
WORK="$HERE/work"; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK" || exit 2
rc=0
"$EXE" > merge.log 2>&1 || rc=1
cat merge.log
for c in missing:"does not exist" cycle:"chain longer than"; do
   leg=${c%%:*}; msg=${c#*:}
   "$EXE" $leg > $leg.log 2>&1; st=$?
   if [[ $st -ne 0 ]] && grep -q "$msg" $leg.log; then echo "PASS $leg: stopped (exit $st): $(grep -m1 "$msg" $leg.log)"
   else echo "FAIL $leg: exit $st, message \"$msg\" $(grep -q "$msg" $leg.log && echo found || echo missing)"; rc=1; fi
done
[[ $rc -eq 0 ]] && echo "ini check PASSED" || echo "ini check FAILED"
exit $rc
