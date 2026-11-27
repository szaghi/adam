#!/bin/bash
# Night 2026-11-27-i160-spacehpc-cpu: why the medooza CPU 3a is 13 % above its 2026-11-23 value. A/B on one node of
# the 3a windows (the dumps of i160/work, repartitioned for 64 ranks): A = the 2026-11-23 executable (sha 5b9bd687,
# no -mdaz-ftz), B = HEAD 5b200ffa with -mdaz-ftz (i160/tree), alternated run by run, 3 x (10, 110) steps.
#SBATCH -J i160-ab3a
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 04:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i160/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i160; R=$D/ab3a; SRC=$D/work/windows/cpu-np64/3a-channel-uniform
declare -A EXE=([A]=/ifs/hpc/home/szaghi/adam-tests/speedup-133/adam/src/tests/flume/speedup/exe/cpu/adam_flume_cpu [B]=$D/tree/src/tests/flume/speedup/exe/cpu/adam_flume_cpu)
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
export OMPI_MCA_coll_hcoll_enable=0 ADAM_SEAM_TIMING=1
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID"
for rep in 1 2 3; do for n in 10 110; do for t in A B; do
   d=$R/$t/r$rep-n$n; rm -rf $d; mkdir -p $d; cp $SRC/r$rep-n$n/*.ini $d/; for f in $SRC/dump/*-restart*; do ln -s $f $d/; done
   t0=$(date +%s.%N); (cd $d && srun --mpi=pmix -n 64 ${EXE[$t]} input.ini > log.txt 2>&1); rc=$?; t1=$(date +%s.%N)
   echo "$(date +%T) $t r$rep-n$n rc $rc wall $(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}') $(grep -a 'forest timing' $d/log.txt | head -1 | awk '{print $7, $8, $9, $10}') $(grep -a 'regrid timing' $d/log.txt | head -1 | awk '{print $8, $9}')"
done; done; done
echo "I160-AB3A-DONE $(date -Is)"
