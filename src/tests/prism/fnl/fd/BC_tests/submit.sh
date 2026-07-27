#!/bin/bash -l

#SBATCH -p prod-gn
#SBATCH --time 23:59:59
#SBATCH --job-name "PRISM"
#SBATCH --exclusive
#SBATCH --nodes=2                    # Number of Nodes
#SBATCH --ntasks=8                   # Total task
#SBATCH --ntasks-per-node=4          # Task MPI per nodo
#SBATCH --cpus-per-task=1            # CPU per task
#SBATCH --gres=gpu:4                 # GPU per nodo (A30)
#SBATCH --gpus-per-node=4            # GPU per nodo (A30)
#SBATCH --mem=0
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

#SBATCH --mail-type=FAIL ## When send mail

## Modules (same used to compile)
#module unload intel/nvidia/cuda-12.3.2 ; module load intel/nvidia_hpc_sdk/nvhpc/25.5
module purge 
module load intel/slurm intel/nvidia_hpc_sdk/nvhpc/25.5 
#module load intel/slurm intel/nvidia_hpc_sdk/nvhpc-hpcx/25.5
#module load intel/slurm intel/nvidia_hpc_sdk/nvhpc-hpcx-cuda12/25.5
## Task execution
cd $SLURM_SUBMIT_DIR

echo "Job ID: $SLURM_JOB_ID"
echo "Nodi allocati: $SLURM_JOB_NODELIST"
echo "Numero di nodi: $SLURM_JOB_NUM_NODES"
echo "Tasks totali: $SLURM_NTASKS"

mpirun -np 8 ./adam_prism_fnl input.ini

