#!/usr/bin/env bash
# RIPPLE, the toy app of the "build your own solver" tutorial (issue #59): the CPU build and its oracles against the
# exact translated profile (uniform order 5, AMR order 2 with conservation, rank invariance, a registered refinement
# criterion, issue #110). The FNL leg r4 runs on medooza.
source "$(dirname "$0")/common.sh"
set -e
fobis build --mode ripple-cpu-gnu --varset ci_gnu
fobis build --mode test-ripple-crest-criterion-cpu-gnu --varset ci_gnu
./src/tests/ripple/check.sh --legs r1,r2,r3,r5,r6
