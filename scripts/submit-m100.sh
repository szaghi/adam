#!/bin/bash
#SBATCH -A IscrB_ICY-ADAM
#SBATCH --output=output-%j.log
#SBATCH --error=error-%j.log

#SBATCH --partition=m100_usr_prod

#SBATCH --qos=m100_qos_bprod
#SBATCH --time=24:00:00

#SBATCH --cpus-per-task 32
#SBATCH --exclusive

#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4

#SBATCH --ntasks-per-socket 2

#SBATCH -N 64

module purge
module load profile/global
module load git python/3.8.2 gnu/8.4.0 hpc-sdk/2021--binary spectrum_mpi/10.4.0--binary

mpirun -n 256 -gpu adam_nasto_gpu

## ogni nodo ha 128 core/ 4 GPU
## QOS DBG   ha massimo 2 ore, 2 nodi
## QOS BPROD ha minimo 17 nodi, massimo 256 nodi e massimo 24 ore
## QOS NORMAL (senza specificare) ha massimo 16 nodi e massimo 24 ore
