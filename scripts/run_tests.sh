#!/usr/bin/env bash
# Run the test suite.
# Each *test* binary exits 0 on success, 1 on failure.
#
# Usage: run_tests.sh [--np N]
#   -n, --np N   Run tests via mpirun with N ranks (default: serial)

NP=0

while [[ $# -gt 0 ]]; do
  case $1 in
    --np | -n ) NP="$2"; shift 2 ;;
    * ) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

if [[ -t 1 ]]; then
  RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; BOLD=$'\033[1m'; RESET=$'\033[0m'
else
  RED=''; GREEN=''; BOLD=''; RESET=''
fi

if [[ $NP -gt 0 ]]; then
  runner=(mpirun -np "$NP")
else
  runner=()
fi

pass=0; fail=0
tmpout=$(mktemp)
trap 'rm -f "$tmpout"' EXIT

for exe in exe/*test*; do
  [[ -x "$exe" ]] || continue
  name=$(basename "$exe")
  if "${runner[@]}" "$exe" > "$tmpout" 2>&1; then
    printf "  ${GREEN}PASS${RESET}  %s\n" "$name"
    pass=$((pass + 1))
  else
    printf "  ${RED}FAIL${RESET}  %s\n" "$name"
    sed 's/^/       /' "$tmpout"
    fail=$((fail + 1))
  fi
done

total=$((pass + fail))
printf "\n${BOLD}%d/%d passed${RESET}\n" "$pass" "$total"
exit $fail
