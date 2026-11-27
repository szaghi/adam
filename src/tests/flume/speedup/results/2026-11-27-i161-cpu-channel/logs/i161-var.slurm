#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 1.3: two instrumented variants of 82c70653 (= 5b200ffa on the CPU path),
# i0 (as is) and i1 (b8641a59's two MPI_Barrier restored), with per-rank timers of the ghost exchange (pack, post,
# waitall, barriers, unpack; instrument.py), built here, then the 3a A/B of b516ftz / i0 / i1 (i161-ab.slurm, 3 pairs).
#SBATCH -J i161-var
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 04:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; T=$D/trees
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID"
for l in i0 i1; do rm -rf "${T:?}/$l"; done
tar -xzf $D/vtrees.tgz -C $T || exit 3
build() { local l=$1; tar -xzf $D/../i160/deps.tgz -C $T/$l && ln -s /home/szaghi/adam/lib $T/$l/lib; cd $T/$l/src/tests/flume/speedup || return 9
   local t0=$(date +%s); ./build.sh cpu --varset iac_gnu --jobs 32 > $D/build-$l.log 2>&1; local r=$?
   local e=exe/cpu/adam_flume_cpu
   echo "$(date +%T) BUILD $l rc $r $(( $(date +%s) - t0 )) s; $(tail -1 $T/$l/COMMIT); lflags: $(grep -a -m1 'Linking flags' builds/cpu/build.txt | sed 's/.*Linking flags: //'); sha256 $( [[ -x $e ]] && sha256sum $e | cut -c1-16)"; }
build i0 & build i1 & wait
for l in i0 i1; do [[ -x $T/$l/src/tests/flume/speedup/exe/cpu/adam_flume_cpu ]] || { echo "BUILD $l FAILED"; grep -a -m5 -i error $D/build-$l.log; exit 4; }; done
LABELS="b516ftz i0 i1" PAIRS=3 CASE=3a-channel-uniform TAGNAME=var bash $D/i161-ab.slurm
for d in $D/ab/var/*/r*-n110; do echo "$d: $(grep -a 'ghost mpi timing' $d/log.txt | head -1)"; done
echo "I161-VAR-DONE $(date -Is)"
