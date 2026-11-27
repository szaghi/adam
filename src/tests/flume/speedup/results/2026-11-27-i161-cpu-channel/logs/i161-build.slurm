#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 1.1/1.2: build flume-cpu-gnu (build.sh cpu --varset iac_gnu) of each probe tree
# (trees.tgz: git archive of the commit; *ftz trees have template-gnu lflags patched to -O2 -mdaz-ftz as 2a67b16c), with
# the deps of deps.tgz (unchanged since 2026-11-22), four builds at a time on 16 cores each.
#SBATCH -J i161-build
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=128G
#SBATCH -t 01:30:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; T=$D/trees
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID; $(gfortran --version | head -1); $(fobis --version 2>&1 | head -1)"
rm -rf "${T:?}" && mkdir -p $T && tar -xzf $D/trees.tgz -C $T || exit 3
LABELS=$(cd $T && ls)
for l in $LABELS; do tar -xzf /ifs/hpc/home/szaghi/adam-tests/speedup-133/i160/deps.tgz -C $T/$l && ln -s /home/szaghi/adam/lib $T/$l/lib; done
build() { local l=$1; cd $T/$l/src/tests/flume/speedup || return 9; local t0=$(date +%s)
   ./build.sh cpu --varset iac_gnu --jobs 16 > $D/build-$l.log 2>&1; local r=$?
   local e=exe/cpu/adam_flume_cpu; local u=$( [[ -x $e ]] && ldd $e | grep -c 'not found')
   echo "$(date +%T) BUILD $l rc $r $(( $(date +%s) - t0 )) s; commit $(head -1 $T/$l/COMMIT); lflags: $(grep -a -m1 'Linking flags' builds/cpu/build.txt | sed 's/.*Linking flags: //'); unresolved ${u:-noexe}; sha256 $( [[ -x $e ]] && sha256sum $e | cut -c1-16)"; }
n=0; for l in $LABELS; do build $l & n=$((n+1)); (( n % 4 == 0 )) && wait; done; wait
echo "I161-BUILD-DONE $(date -Is)"
