#!/bin/bash
#SBATCH -A iacroma -p prod-gn -N 1 --gres=gpu:4 -n 4 -c 8 --mem=120G -t 04:00:00 -J i156-gA
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/i156/gate-A.%j.out
# Night 2026-11-25-i156-l2, step 01 gate, the control: tree A = the gate tree with the two FLUME objects of d68d5f97
# (the base), built and run like B (shock-cylinder-ib, CPU and FNL), so B's outputs can be compared byte for byte with
# A's on the same machine and toolchain (the goldens differ from a medooza build in the last digit already at step 0).
set -uo pipefail
G=/ifs/hpc/home/szaghi/adam-tests/i156
mkdir -p $G/A && rsync -a --exclude /exe --exclude 'work-*' $G/adam/ $G/A/
cp $G/adam_flume_cpu_object.F90 $G/A/src/app/flume/cpu/ && cp $G/adam_flume_fnl_object.F90 $G/A/src/app/flume/fnl/
grep -c 'e < self%ib%n_eikonal' $G/A/src/app/flume/*/adam_flume_*_object.F90
cd $G/A
( source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1
  fobis build --mode flume-cpu-gnu --varset iac_gnu > build-cpu.log 2>&1; echo "A-CPU-BUILD-RC $?"
  cp -p exe/adam_flume_cpu exe/adam_flume_cpu_fft
  cd src/tests/flume/regression && REGRESSION_CASES="shock-cylinder-ib" ./run.sh cpu --no-build; echo "A-CPU-RC $?" )
( source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1; export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
  cd src/tests/flume/regression && REGRESSION_CASES="shock-cylinder-ib" ./run.sh fnl --varset iac_nvf > $G/A/fnl-run.log 2>&1; echo "A-FNL-RC $?" )
C=src/tests/flume/regression/shock-cylinder-ib
for b in cpu fnl; do for f in $G/A/$C/work-$b/*.dat $G/A/$C/work-$b/digest.txt; do
   n=$(basename $f); cmp -s $f $G/adam/$C/work-$b/$n && echo "AB-$b same $n" || echo "AB-$b DIFF $n"; done; done
echo "GATE-A-DONE $(date -Is)"
