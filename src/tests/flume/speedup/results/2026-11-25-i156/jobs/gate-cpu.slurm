#!/bin/bash
#SBATCH -A iacroma -p prod-gn -N 1 -n 4 -c 8 --mem=60G -t 02:00:00 -J i156-gcpu
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/i156/gate-cpu.%j.out
# Night 2026-11-25-i156-l2, step 01 gate, CPU: shock-cylinder-ib against the goldens. medooza's gnu.sh has no FFTW, so
# the -fft executable (self-gravity cases only) cannot link: the plain executable stands in for it (the case does not
# use self-gravity), built fresh here, then run.sh --no-build.
set -uo pipefail
G=/ifs/hpc/home/szaghi/adam-tests/i156/adam
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1
cd $G && fobis build --mode flume-cpu-gnu --varset iac_gnu > exe/build-cpu-gate.log 2>&1; echo "BUILD-RC $?"
ls -la exe/adam_flume_cpu && cp -p exe/adam_flume_cpu exe/adam_flume_cpu_fft
cd $G/src/tests/flume/regression
REGRESSION_CASES="shock-cylinder-ib" ./run.sh cpu --no-build; echo "CPU-RC $?"
echo "GATE-DONE $(date -Is)"
