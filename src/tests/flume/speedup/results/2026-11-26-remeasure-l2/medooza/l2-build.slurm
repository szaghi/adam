#!/bin/bash
# Night 2026-11-26-remeasure-l2: the timing tree adam-l2 = HEAD 0984b5f5 (code of 44cfdfa8, #156 lever 2), copied from
# the final gate tree i156/E (its sources are HEAD's: md5 list checked first), then build.sh cpu and fnl-nvf.
#SBATCH -A iacroma -p prod-gn -N 1 -n 1 -c 32 --mem=64G -t 03:00:00 -J l2-build
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/l2-build.%j.out
set -uo pipefail
W=/ifs/hpc/home/szaghi/adam-tests/speedup-133; A=$W/adam-l2; E=/ifs/hpc/home/szaghi/adam-tests/i156/E
rm -rf $A && mkdir -p $A && rsync -a --exclude /exe --exclude 'work-*' --exclude /.git $E/ $A/
ln -sfn /home/szaghi/adam/lib $A/lib 2> /dev/null || true
cd $A && md5sum -c --quiet $W/l2-head-md5.txt > $W/l2-md5-check.txt 2>&1; echo "MD5 mismatches: $(grep -c FAILED $W/l2-md5-check.txt)"
grep FAILED $W/l2-md5-check.txt | head
echo "HEAD 0984b5f5 (code 44cfdfa8): copied from i156/E, md5-checked against git ls-files of HEAD" > $A/PROVENANCE.txt
# the compute node's clock runs behind the login node's: stamp the sources here so the executables are newer
find $A/src -type f -exec touch {} +
S=$A/src/tests/flume/speedup; t0=$(date +%s)
( source /ifs/hpc/home/szaghi/gnu/gnu.sh; cd "$S" && ./build.sh cpu --varset iac_gnu --jobs 16 > $W/l2-build-cpu.log 2>&1; echo $? > $W/l2-build-cpu.rc ) &
( source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh; cd "$S" && ./build.sh fnl-nvf --varset iac_nvf --jobs 16 > $W/l2-build-fnl.log 2>&1; echo $? > $W/l2-build-fnl.rc ) &
wait
echo "$(date +%T) builds cpu rc $(cat $W/l2-build-cpu.rc), fnl rc $(cat $W/l2-build-fnl.rc) after $(( $(date +%s) - t0 )) s"
ls -la $S/exe/*/
echo "L2-BUILD-DONE"
