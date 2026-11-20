#!/usr/bin/env bash
# PRISM CPU regression suite: two cases at a time (np 2 each) fill the 4-core runner (#80 P5 dispatcher). The produced
# digests and residuals are uploaded (a golden is captured from a known-good CI run; a failing run is diagnosed).
source "$(dirname "$0")/common.sh"
REGRESSION_JOBS=2 ./src/tests/prism/regression/run.sh cpu --varset ci_gnu
rc=$?
collect src/tests/prism/regression/**/work-cpu/digest.txt src/tests/prism/regression/**/work-cpu/*-residuals.dat
exit $rc
