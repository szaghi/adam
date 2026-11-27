#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 2.1 timing: build fixA (HEAD + one MPI_Barrier before the posts) and fixB
# (HEAD + b8641a59's two barriers restored) with build.sh cpu --varset iac_gnu, then the 3a A/B of
# b516ftz / fixA / fixB / 82c7 (= HEAD on the CPU path), 3 pairs (i161-ab2.slurm). Writes trees/fix*/READY when built.
#SBATCH -J i161-fix
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 05:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; T=$D/trees
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID"
for l in fixA fixB; do rm -rf "${T:?}/$l"; done
tar -xzf $D/ftrees.tgz -C $T || exit 3
build() { local l=$1; tar -xzf $D/../i160/deps.tgz -C $T/$l && ln -s /home/szaghi/adam/lib $T/$l/lib; cd $T/$l/src/tests/flume/speedup || return 9
   local t0=$(date +%s); ./build.sh cpu --varset iac_gnu --jobs 32 > $D/build-$l.log 2>&1; local r=$?
   local e=exe/cpu/adam_flume_cpu
   echo "$(date +%T) BUILD $l rc $r $(( $(date +%s) - t0 )) s; $(tail -1 $T/$l/COMMIT); lflags: $(grep -a -m1 'Linking flags' builds/cpu/build.txt | sed 's/.*Linking flags: //'); unresolved $( [[ -x $e ]] && ldd $e | grep -c 'not found'); sha256 $( [[ -x $e ]] && sha256sum $e | cut -c1-16)"
   [[ $r -eq 0 && -x $e ]] && touch $T/$l/READY; }
build fixA & build fixB & wait
for l in fixA fixB; do [[ -f $T/$l/READY ]] || { echo "BUILD $l FAILED"; grep -a -m5 -i error $D/build-$l.log; exit 4; }; done
LABELS="b516ftz fixA fixB 82c7" PAIRS=3 CASE=3a-channel-uniform TAGNAME=fix bash $D/i161-ab2.slurm
echo "I161-FIX-DONE $(date -Is)"
