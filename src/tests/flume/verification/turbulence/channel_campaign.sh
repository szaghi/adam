#!/usr/bin/env bash
# The channel campaign of issue #119 T3 (T3.2, T3.7): run the case list with restart chains, as tgv_campaign.sh. A case line:
# name nx ny nz mach reynolds [more make_channel.py options]. Exit 0: all done; 3: the budget ran out, submit again.
#
# Usage: channel_campaign.sh <root> <executable> <np> <budget seconds> <case file> [name filter]
set -uo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$1" ; EXE="$2" ; NP="$3" ; BUDGET="$4" ; CASE_FILE="$5" ; FILTER="${6:-.}"
T_END=$(( $(date +%s) + BUDGET ))

while read -r name nx ny nz mach re extra; do
   [[ -n "$name" && "$name" != \#* ]] || continue
   grep -q "$FILTER" <<< "$name" || continue
   work="$ROOT/$name"
   [[ -f "$work/done" ]] && continue
   left=$(( T_END - $(date +%s) ))
   (( left > 120 )) || exit 3
   mkdir -p "$work"
   restart=() ; ls "$work"/channel-restart.time > /dev/null 2>&1 && restart=(--restart)
   # a checkpoint every 10 time units, a restart save every 2, histories every 10 steps
   python3 -I "$CASE_DIR/make_channel.py" "$work/input.ini" --nx "$nx" --ny "$ny" --nz "$nz" --mach "$mach" --reynolds "$re" --power-law 0.7 \
           --it-save 100000000 --restart-save 500 --history 10 "${restart[@]}" $extra
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
done < "$CASE_FILE"
echo "$(date +%T) campaign done"
