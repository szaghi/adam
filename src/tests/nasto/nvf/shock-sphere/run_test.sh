#!/bin/bash

if [ $# -eq 0 ] ; then
   mpirun -n 1 adam_nasto_nvf adam-nasto-shock-sphere.ini | tee run_test.log
else
   mpirun -n 1 compute-sanitizer adam_nasto_nvf adam-nasto-shock-sphere.ini | tee run_test.log
fi
