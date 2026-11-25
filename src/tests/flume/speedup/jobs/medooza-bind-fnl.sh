#!/bin/bash
# FLUME FNL with each rank pinned to its A30's NUMA node on medooza prod-gn (README, "Rank placement in a GPU node").
# nvidia-smi topo -m (job 1369679): GPU0 NUMA 0 (cores 0-15), GPU1 NUMA 1 (16-31), GPU2 and GPU3 NUMA 3 (48-63).
# The rank of local id l sees only GPU l, runs on 8 cores of that GPU's NUMA node and allocates host memory there.
# Pass it as run.sh --exe and launch with a step that owns every core (the default step's cpuset holds 0-15, 32-47):
#   ./run.sh windows ... --exe jobs/medooza-bind-fnl.sh --launcher "srun --mpi=pmix -c 16 --cpu-bind=none -n"
# in an --exclusive job. Measured: no speedup (within 0.8 % on 1 and 2 nodes); GPU 0 keeps the 0.9 GB the other
# ranks' contexts take there when unpinned. BIND_FNL_EXE overrides the executable.
l=${SLURM_LOCALID:?not under srun}
cores=(0-7 16-23 48-55 56-63)
numa=(0 1 3 3)
X=${BIND_FNL_EXE:-$(dirname "$(readlink -f "$0")")/../exe/fnl-nvf/adam_flume_fnl}
export CUDA_VISIBLE_DEVICES=$l
if command -v numactl > /dev/null; then
   exec numactl --physcpubind="${cores[$l]}" --membind="${numa[$l]}" "$X" "$@"
fi
exec taskset -c "${cores[$l]}" "$X" "$@"
