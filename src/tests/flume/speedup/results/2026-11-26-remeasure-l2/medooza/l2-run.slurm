#!/bin/bash
# Night 2026-11-26-remeasure-l2: windows of the sphere at HEAD 0984b5f5 (#156 lever 2) in the work tree l2 (the matrix
# and grid-study dumps linked; outputs only there). BK = fnl-nvf | cpu; the arguments are cases, `case@grid` for the
# grid-study cases (their conf). Nodes, GPUs and ranks come from sbatch's options.
#SBATCH -A iacroma -p prod-gn
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/%x.%j.out
set -uo pipefail
W=/ifs/hpc/home/szaghi/adam-tests/speedup-133; WK=$W/l2; S=$W/adam-l2/src/tests/flume/speedup; BK=${BK:?BK}
mkdir -p $WK/dumps
for d in leg2 w2-sphere-l5; do ln -sfn $W/work/dumps/$d $WK/dumps/$d; done
for d in wg-b38 wg-l4b38; do ln -sfn $W/grid-study/dumps/$d $WK/dumps/$d; done
if [[ $BK == cpu ]]; then source /ifs/hpc/home/szaghi/gnu/gnu.sh; else
   source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5; export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
   nvidia-smi --query-gpu=index,memory.free --format=csv,noheader; fi
echo "job $SLURM_JOB_ID $SLURM_JOB_NODELIST np $SLURM_NTASKS $BK cases $* $(date -Is)"; cd "$S" || exit 3
rc=0
for a in "$@"; do
   c=${a%@*}; conf=(); [[ $a == *@grid ]] && conf=(--conf $W/grid-study.conf)
   ./run.sh windows --backend $BK --np $SLURM_NTASKS --work $WK --cases "$c" "${conf[@]}" --launcher "srun --mpi=pmix -n"
   r=$?; (( r == 0 )) || rc=1
   echo "L2 $BK np $SLURM_NTASKS $c rc $r"
done
python3 -I analyze.py --work $WK --out $WK/results-$SLURM_JOB_ID > /dev/null 2>&1
echo "L2-RUN-DONE rc $rc $(date -Is)"
