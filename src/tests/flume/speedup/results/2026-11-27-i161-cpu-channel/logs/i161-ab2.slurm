#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161: A/B on one node of the FLUME CPU windows (the dumps of i160/work, 64 ranks):
# the executables of LABELS (i161/trees/<label>, or A = the #133 matrix's CPU build sha256 5b9bd687) alternated run by
# run, the order rotated each repetition, PAIRS x (10, 110) steps. The step cost is (loop(110) - loop(10)) / 100.
#   sbatch --export=ALL,[SRCROOT=<windows dir>,]LABELS="A b516 82c7",PAIRS=3,CASE=3a-channel-uniform,TAGNAME=t11 i161-ab.slurm
#SBATCH -J i161-ab
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 04:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; R=$D/ab/${TAGNAME:?}; SRC=${SRCROOT:-$D/../i160/work/windows/cpu-np64}/${CASE:?}
read -ra L <<< "${LABELS:?}"; PAIRS=${PAIRS:-3}
exe() { [[ $1 == A ]] && echo /ifs/hpc/home/szaghi/adam-tests/speedup-133/adam/src/tests/flume/speedup/exe/cpu/adam_flume_cpu || echo $D/trees/$1/src/tests/flume/speedup/exe/cpu/adam_flume_cpu; }
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
export OMPI_MCA_coll_hcoll_enable=0 ADAM_SEAM_TIMING=1
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID, case $CASE, labels ${L[*]}, pairs $PAIRS"
for t in "${L[@]}"; do echo "exe $t $(exe $t) sha256 $(sha256sum $(exe $t) | cut -c1-16)"; done
for rep in $(seq 1 $PAIRS); do
   k=$(( (rep - 1) % ${#L[@]} )); order=("${L[@]:k}" "${L[@]:0:k}")
   for n in 10 110; do for t in "${order[@]}"; do
      d=$R/$t/r$rep-n$n; rm -rf "${d:?}"; mkdir -p $d; cp $SRC/r$rep-n$n/*.ini $d/; for f in $SRC/dump/*-restart*; do ln -s $f $d/; done
      t0=$(date +%s.%N); (cd $d && srun --mpi=pmix -n 64 $(exe $t) input.ini > log.txt 2>&1); rc=$?; t1=$(date +%s.%N)
      echo "$(date +%T) $t r$rep-n$n rc $rc wall $(awk -v a=$t0 -v b=$t1 'BEGIN{printf "%.1f", b-a}') $(grep -a 'forest timing' $d/log.txt | head -1 | awk '{print $7, $8, $9, $10}') $(grep -a 'regrid timing' $d/log.txt | head -1 | awk '{print $8, $9}')"
   done; done
done
echo "I161-AB-DONE $(date -Is)"
