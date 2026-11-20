#!/usr/bin/env bash
# FLUME CPU regression suite (issue #35, P8), then the OpenMP build against the serial one, bit for bit, on the
# immersed-boundary case (np 2 x 2 threads), run whatever the suite's outcome. The FNL goldens are checked on medooza:
# the hosted runners have no GPU.
source "$(dirname "$0")/common.sh"
REGRESSION_JOBS=2 ./src/tests/flume/regression/run.sh cpu --varset ci_gnu
rc=$?
./src/tests/flume/regression/run-omp-bitwise.sh --threads 2 --np 2 --varset ci_gnu || rc=1
collect src/tests/flume/regression/**/work-cpu/digest.txt \
        src/tests/flume/regression/**/work-cpu/*-residuals.dat \
        src/tests/flume/regression/**/work-cpu/*-conservation_history.dat \
        src/tests/flume/regression/**/work-omp-*/digest.txt \
        src/tests/flume/regression/**/work-omp-*/*.dat \
        src/tests/flume/regression/**/work-omp-*/run.log
exit $rc
