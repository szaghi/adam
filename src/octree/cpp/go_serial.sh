#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00

make cleandata 
make 
./miniamr >& miniamr.log
tar cvf output.tar *.vt[mr]
