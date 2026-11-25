#!/bin/bash
# GPU to GPU and host to host bandwidth between the first two nodes of a PBS job (ucx_perftest tag_bw, 33 MB and 8 MB,
# 60 iterations: the medooza protocol, results/2026-11-24/phase/network.txt), with the site's UCX_NET_DEVICES and with
# the GPU's own HCA. Run from a spacehpc-fnl.pbs plan (`sh` line): needs $PBS_NODEFILE and the nvhpc-hpcx-cuda12 module.
set -uo pipefail
mapfile -t N < <(sort -u "$PBS_NODEFILE")
(( ${#N[@]} >= 2 )) || { echo "netbw.sh: needs 2 nodes" ; exit 0 ; }
P=$(dirname "$(command -v ucx_info)")/ucx_perftest
echo "ucx_perftest $P; nodes ${N[0]} -> ${N[1]}; site UCX_NET_DEVICES=${UCX_NET_DEVICES:-unset}"
echo "columns of the final line: iterations, latency (typical, avg, overall usec), bandwidth (average, overall MB/s), message rate"
one() { # label size mem gpu env...
   local lab=$1 s=$2 m=$3 g=$4 port=$((15000 + RANDOM % 2000)); shift 4
   mpirun -np 1 --host "${N[1]}" -x CUDA_VISIBLE_DEVICES="$g" $(for e in "$@"; do echo "-x $e"; done) env "$@" "$P" -t tag_bw -m "$m" -s "$s" -n 60 -p "$port" > /dev/null 2>&1 & local sp=$!
   sleep 6
   echo "$lab: $(env CUDA_VISIBLE_DEVICES="$g" "$@" "$P" "${N[1]}" -t tag_bw -m "$m" -s "$s" -n 60 -p "$port" 2>&1 | grep Final)"
   wait $sp
}
one "inter-node host->host 33 MB [site]" 33554432 host 0
one "inter-node cuda->cuda 33 MB GPU0 [site]" 33554432 cuda 0
one "inter-node cuda->cuda 33 MB GPU1 [site]" 33554432 cuda 1
one "inter-node cuda->cuda 8 MB GPU0 [site]" 8388608 cuda 0
one "inter-node cuda->cuda 33 MB GPU0 [own HCA mlx5_1]" 33554432 cuda 0 UCX_NET_DEVICES=mlx5_1:1
one "inter-node cuda->cuda 33 MB GPU1 [own HCA mlx5_0]" 33554432 cuda 1 UCX_NET_DEVICES=mlx5_0:1
one "inter-node cuda->cuda 33 MB GPU0 [no GPUDirect RDMA]" 33554432 cuda 0 UCX_IB_GPU_DIRECT_RDMA=no
echo NETBW-DONE
