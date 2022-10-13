#!/bin/bash
#SBATCH -A IscrB_ICY-ADAM
#SBATCH --output=output-pvpython-%j.log
#SBATCH --error=error-pvpython-%j.log

#SBATCH --partition=m100_usr_prod
##SBATCH --qos=m100_qos_dbg

#SBATCH --time=08:00:00

#SBATCH --cpus-per-task 32
#SBATCH --exclusive

#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1

#SBATCH --ntasks-per-socket 1

#SBATCH -N 1

module purge
module load profile/global
module load cuda/11.0 spectrum_mpi/10.4.0--binary  gnu/8.4.0 paraview_egl/5.8.1

# $1 = paraview (pv)python macro file
# $2 = basename of XDMF file(s)

template=`echo $1 | sed -e 's/\.template//'`
for xdmf in $2*.xdmf ; do
   it=`echo $xdmf | awk -F '-' '{print $3}' | sed -e 's/^0*//' | sed -e 's/\.xdmf//'`
   if [ -f "z-slice-$it.vtm" ]; then
      echo "$it already processed"
   else
      echo $xdmf
      cat $template.template | sed -e "s/AAAIT/-$it/" | sed -e "s/BBFILENAMEXDMF/$xdmf/" > $template-$it.py
      mpirun -n 1 pvpython --force-offscreen-rendering $template-$it.py
   fi
done

## ogni nodo ha 128 core/ 4 GPU
## QOS DBG   ha massimo 2 ore, 2 nodi
## QOS BPROD ha minimo 17 nodi, massimo 256 nodi e massimo 24 ore
## QOS NORMAL (senza specificare) ha massimo 16 nodi e massimo 24 ore
