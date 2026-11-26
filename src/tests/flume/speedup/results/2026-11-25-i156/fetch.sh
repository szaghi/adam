#!/bin/bash
# #156 lever 2 (2026-11-25): mirror the reduced SpaceHPC results into this directory: the analyze.py summaries of the
# sphere windows per executable (A = d68d5f97, B = A + no ghost update after the last eikonal sweep), the phase-timer
# reports (the PHASE TIMING lines of each 110-step run of wg-l4b38 and 2a), the job logs.
# Run from WSL: ./fetch.sh   (needs `ssh spacehpc`).
set -euo pipefail
cd "$(dirname "$0")"
W=/lustre/projects/1026/CALCOLI/adam/speedup-133
mkdir -p spacehpc/phase spacehpc/logs
for t in A B; do for e in csv md; do rsync -a "spacehpc:$W/i156/$t/results/summary.$e" "spacehpc/summary-$t.$e"; done; done
# shellcheck disable=SC2016
ssh spacehpc "cd $W/phase && ls -d i156*/windows/*/*/r?-n110 2> /dev/null" | while read -r d; do
   t=${d%%/*}; np=${d#*/windows/fnl-nvf-np}; np=${np%%/*}; c=$(basename "$(dirname "$d")"); r=$(basename "$d")
   ssh -n spacehpc "grep -h 'PHASE TIMING' $W/phase/$d/log.txt" > "spacehpc/phase/$t-$c-np$np-$r.txt" || true
done
rsync -a "spacehpc:$W/logs/"{184502,184501}.pbs.OU "spacehpc:$W/i156/i156-amr.184502.pbs.out" spacehpc/logs/ 2> /dev/null || true
ls -R spacehpc | head -40
