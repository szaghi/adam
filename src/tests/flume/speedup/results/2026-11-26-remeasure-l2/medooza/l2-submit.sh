#!/bin/bash
# submit the l2 runs after the build job (at most 4 jobs of the night queued at once; the rest by hand later)
W=/ifs/hpc/home/szaghi/adam-tests/speedup-133; B=${1:?build job id}
cd $W
F="--gres=gpu:4 --ntasks-per-node=4 --cpus-per-task=8 --mem=96G"
sbatch --parsable -J l2-f1 --dependency=afterok:$B -N 1 $F -t 03:00:00 --export=ALL,BK=fnl-nvf l2-run.slurm 2a wg-b38@grid
sbatch --parsable -J l2-c1 --dependency=afterok:$B -N 1 --ntasks-per-node=64 --mem=256G -t 04:00:00 --export=ALL,BK=cpu l2-run.slurm 2a
sbatch --parsable -J l2-f8 --dependency=afterok:$B -N 8 $F -t 04:00:00 --export=ALL,BK=fnl-nvf l2-run.slurm 2a wg-l4b38@grid w2-sphere-l5
