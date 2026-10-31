#!/usr/bin/env bash
# The Taylor-Green ILES campaign of issue #119 (T0.3, T0.4): run the case list, cheapest first, with restart chains because
# the GPU partition stops a job after an hour. Each case lives in <root>/<name>; `done` marks a finished one, a killed run
# restarts from its last restart save (tgv_report.py keeps the last copy of the rows run twice). Exit 0: all cases done;
# exit 3: the wall-clock budget ran out, submit again.
#
# Usage: tgv_campaign.sh <root> <executable> <np> <budget seconds> [name filter, a grep pattern] [case file]
# A case line: name n mach scheme_space scheme [more make_tgv.py options, a later one overrides].
set -uo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$1" ; EXE="$2" ; NP="$3" ; BUDGET="$4" ; FILTER="${5:-.}" ; CASE_FILE="${6:-}"
T_END=$(( $(date +%s) + BUDGET ))

# name n mach scheme_space scheme: 64 and 128 first, the supersonic ones (6 times fewer steps) before 256^3
CASES="
m0.1-weno-5-n64 64 0.1 weno weno-u-5
m0.1-riemann-5-n64 64 0.1 weno-riemann weno-u-5
m1.25-weno-5-n64 64 1.25 weno weno-u-5
m1.25-riemann-5-n64 64 1.25 weno-riemann weno-u-5
m0.1-weno-5-n128 128 0.1 weno weno-u-5
m0.1-riemann-5-n128 128 0.1 weno-riemann weno-u-5
m1.25-weno-5-n128 128 1.25 weno weno-u-5
m1.25-riemann-5-n128 128 1.25 weno-riemann weno-u-5
m0.1-weno-7-n128 128 0.1 weno weno-u-7
m0.1-weno-9-n128 128 0.1 weno weno-u-9
m1.25-weno-5-n256 256 1.25 weno weno-u-5
m1.25-riemann-5-n256 256 1.25 weno-riemann weno-u-5
m0.1-weno-5-n256 256 0.1 weno weno-u-5
"

[[ -n "$CASE_FILE" ]] && CASES="$(cat "$CASE_FILE")"
while read -r name n mach space scheme extra; do
   [[ -n "$name" && "$name" != \#* ]] || continue
   grep -q "$FILTER" <<< "$name" || continue
   work="$ROOT/$name"
   [[ -f "$work/done" ]] && continue
   left=$(( T_END - $(date +%s) ))
   (( left > 120 )) || exit 3
   mkdir -p "$work"
   # about one checkpoint per unit of time, from dt = CFL h / (3 (V + c)), and a restart save every 3 % of the run
   read -r it_save restart_save < <(python3 -I -c "import math; d = 0.5 * 2 * math.pi / $n / (3 * (1 + 1 / $mach)); print(round(1 / d), max(100, round(0.6 / d)))")
   restart=() ; ls "$work"/tgv-restart.time > /dev/null 2>&1 && restart=(--restart)
   python3 -I "$CASE_DIR/make_tgv.py" "$work/input.ini" --n "$n" --levels 2 --mach "$mach" --scheme "$scheme" --scheme-space "$space" \
           --time-max 20.0 --it-save "$it_save" --restart-save "$restart_save" --history 10 "${restart[@]}" $extra
   echo "$(date +%T) $name: ${restart[*]:-start}, ${left} s left"
   (cd "$work" && timeout "$left" mpirun -np "$NP" "$EXE" input.ini >> log.txt 2>&1 < /dev/null)
   rc=$?
   if (( rc == 0 )); then
      touch "$work/done" ; echo "$(date +%T) $name: done"
   elif (( rc == 124 )); then
      exit 3
   else
      echo "$(date +%T) $name: FAILED, exit $rc, see $work/log.txt" ; exit 1
   fi
done <<< "$CASES"
echo "$(date +%T) campaign done"
