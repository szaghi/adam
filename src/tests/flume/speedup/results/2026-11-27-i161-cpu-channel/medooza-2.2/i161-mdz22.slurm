#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 2.2 (medooza): the 3a CPU windows with the fix (trees/fixA, built by the
# i161-fix job), run.sh windows --cases 3a on 64 ranks from the matrix's warm-up dumps (work/dumps/leg3, linked read
# only) into its own work tree, then analyze.py. As i160-mdz-cpu.slurm (#160 1.3), without the build: the kit of the i160
# tree (HEAD's), the executable trees/fixA's (--exe).
#SBATCH -J i161-mdz22
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 03:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; T=$D/trees/fixA; W=$D/work22
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
export OMPI_MCA_coll_hcoll_enable=0
[[ -f $T/READY ]] || { echo "fixA not built"; exit 4; }
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID; exe sha256 $(sha256sum $T/src/tests/flume/speedup/exe/cpu/adam_flume_cpu | cut -c1-16)"
rm -rf "${W:?}" && mkdir -p $W/dumps && ln -sfn /ifs/hpc/home/szaghi/adam-tests/speedup-133/work/dumps/leg3 $W/dumps/leg3
cd $D/../i160/tree/src/tests/flume/speedup
t0=$(date +%s); ./run.sh windows --backend cpu --np 64 --work $W --launcher "srun --mpi=pmix -n" --cases 3a --exe $T/src/tests/flume/speedup/exe/cpu/adam_flume_cpu; r=$?
echo "WINDOWS rc $r after $(( $(date +%s) - t0 )) s"
python3 -I analyze.py --work $W --out $W/results > $D/analyze22.log 2>&1; echo "ANALYZE rc $?"
grep -E "^\| 3 \|" $W/results/summary.md | head -1
echo "I161-MDZ22-DONE $(date -Is)"
