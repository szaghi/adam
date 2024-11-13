#!/bin/bash

backend=${PWD##*/}
backend=${backend:-/}
dbg=''
clean='no'
while [ $# -gt 0 ]; do
  case "$1" in
    "-c")
      shift; clean='yes'
      ;;
    "-d")
      shift; dbg=''
      ;;
    *)
      echo; echo "Unknown switch $1"; exit 1
      ;;
  esac
  shift
done
if [ "$clean" == "yes" ] ; then
   echo "clean directory"
   ./clean.sh
fi
echo "mpirun -n 1 $dbg adam_nasto_$backend adam-nasto*.ini | tee run_test.log"
mpirun -n 1 $dbg adam_nasto_$backend adam-nasto*.ini | tee run_test.log
exit 0
