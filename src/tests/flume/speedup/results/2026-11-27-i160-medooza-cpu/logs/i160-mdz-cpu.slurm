#!/bin/bash
# Night 2026-11-27-i160-spacehpc-cpu, #160 step 1.3: the medooza CPU windows of 3a (channel) and 1a (step) with
# -mdaz-ftz. Builds flume-cpu-gnu of HEAD 5b200ffa (tree.tgz + deps.tgz) with build.sh cpu --varset iac_gnu, then
# run.sh windows --cases 3a,1a on 64 ranks from the matrix's warm-up dumps (work/dumps/leg{1,3}, linked read only) into
# its own work tree (never the matrix's: its results stay), then analyze.py.
#   sbatch i160-mdz-cpu.slurm
#SBATCH -J i160-mdz-cpu
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=256G
#SBATCH -t 12:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i160/%x.%j.out
set -uo pipefail
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i160; T=$D/tree; W=$D/work
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
export OMPI_MCA_coll_hcoll_enable=0
echo "start $(date -Is) on $(hostname), job $SLURM_JOB_ID; $(gfortran --version | head -1)"
lscpu | grep -E '^(Model name|Socket|Core|Thread)' || true
rm -rf $T && mkdir -p $T && tar -xzf $D/tree.tgz -C $T && tar -xzf $D/deps.tgz -C $T || exit 3
ln -s /home/szaghi/adam/lib $T/lib
cd $T && git init -q && git add -A src fobos fobos.d > /dev/null 2>&1 && git -c user.name=night -c user.email=night@localhost commit -qm "tree of ADAM 5b200ffa" && echo "tree commit $(git rev-parse --short HEAD) = ADAM 5b200ffa"
cd $T/src/tests/flume/speedup
t0=$(date +%s); ./build.sh cpu --varset iac_gnu --jobs 32 > $D/build.log 2>&1; r=$?
echo "BUILD rc $r after $(( $(date +%s) - t0 )) s; mdaz-ftz in build.txt: $(grep -c -- -mdaz-ftz builds/cpu/build.txt)"
ldd exe/cpu/adam_flume_cpu | grep 'not found' && { echo "BUILD unresolved libraries"; exit 3; }
[[ $r -eq 0 ]] || exit $r
mkdir -p $W/dumps
for g in leg1 leg3; do ln -sfn /ifs/hpc/home/szaghi/adam-tests/speedup-133/work/dumps/$g $W/dumps/$g; done
t0=$(date +%s)
./run.sh windows --backend cpu --np 64 --work $W --launcher "srun --mpi=pmix -n" --cases 3a,1a; r=$?
echo "WINDOWS rc $r after $(( $(date +%s) - t0 )) s"
python3 -I analyze.py --work $W --out $W/results > $D/analyze.log 2>&1; echo "ANALYZE rc $?"
grep -E "^\| [13] \|" $W/results/summary.md
echo "I160-MDZ-DONE $(date -Is)"
