#!/bin/bash
#SBATCH -A iacroma -p prod-gn -N 1 --gres=gpu:4 -n 4 -c 8 --mem=120G -t 04:00:00 -J l2b-gD
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/i156/gate-D.%j.out
# Night 2026-11-26-i156-l2b (#156 L2 part 2): tree D = the part-1 gate tree (adam/, 4a029436's code, its
# shock-cylinder-ib outputs of jobs 1369857/1369858 are A) with the night's changed files, built clean and run like A
# (shock-cylinder-ib, CPU and FNL), then D's outputs compared byte for byte with A's.
set -uo pipefail
G=/ifs/hpc/home/szaghi/adam-tests/i156
rm -rf $G/D && mkdir -p $G/D && rsync -a --exclude /exe --exclude 'work-*' $G/adam/ $G/D/ && rsync -a $G/D-src/ $G/D/
# the compute nodes' clock runs behind the login node's: stamp every source with this node's clock, so the fresh
# executables are newer than their sources (run.sh's stale-executable guard)
find $G/D/src -type f -exec touch {} +
grep -c 'ib_rows=.true.' $G/D/src/app/flume/*/adam_flume_*_object.F90
cd $G/D
( source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1
  fobis build --mode flume-cpu-gnu --varset iac_gnu > build-cpu.log 2>&1; echo "D-CPU-BUILD-RC $?"
  cp -p exe/adam_flume_cpu exe/adam_flume_cpu_fft
  cd src/tests/flume/regression && REGRESSION_CASES="shock-cylinder-ib" ./run.sh cpu --no-build; echo "D-CPU-RC $?" )
( source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1; export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
  cd src/tests/flume/regression && REGRESSION_CASES="shock-cylinder-ib" ./run.sh fnl --varset iac_nvf > $G/D/fnl-run.log 2>&1; echo "D-FNL-RC $?"
  grep -aE "PASS|FAIL" $G/D/fnl-run.log | tail -3 )
C=src/tests/flume/regression/shock-cylinder-ib
for b in cpu fnl; do for f in $G/adam/$C/work-$b/*.dat $G/adam/$C/work-$b/digest.txt; do
   n=$(basename $f); cmp -s $f $G/D/$C/work-$b/$n && echo "AD-$b same $n" || echo "AD-$b DIFF $n"; done; done
echo "GATE-D-DONE $(date -Is)"
