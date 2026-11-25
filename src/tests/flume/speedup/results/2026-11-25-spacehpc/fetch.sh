#!/bin/bash
# Mirror the reduced SpaceHPC results of #133 (2026-11-25) into this directory: analyze.py summaries, the phase-timer
# reports (the PHASE TIMING lines of each 110-step run), the job logs and their GPU-health summaries, the bandwidth tests.
# Run from WSL: ./fetch.sh   (needs `ssh spacehpc`).
set -euo pipefail
cd "$(dirname "$0")"
R=/lustre/projects/1026/CALCOLI/adam/speedup-133
for t in matrix grid-study; do mkdir -p "$t" ; rsync -a "spacehpc:$R/$t/results/summary."{md,csv} "$t/" ; done
mkdir -p phase/reports logs
# shellcheck disable=SC2016
ssh spacehpc "cd $R/phase && for d in */windows/*/*/r1-n110; do echo \"\$d\"; done" | while read -r d; do
   t=${d%%/*}; np=${d#*/windows/fnl-nvf-np}; np=${np%%/*}; c=$(basename "$(dirname "$d")")
   ssh -n spacehpc "grep -h 'PHASE TIMING' $R/phase/$d/log.txt" > "phase/reports/$t-$c-np$np.txt" || true
done
rsync -a "spacehpc:$R/logs/"{*.pbs.OU,netbw-*.txt} logs/ 2> /dev/null || true
ls -R | head -60
