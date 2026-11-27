#!/bin/bash
# Night 2026-11-27-i161-cpu-channel, #161 2.1 gate (CPU): clean GNU builds of base (HEAD e8ac01f5 = 5b200ffa on the
# CPU path) and fix (base + one MPI_Barrier before the posts of update_ghost_mpi); the FLUME and PRISM CPU suites and
# RIPPLE on fix; DV-1 (ranks, empty ranks: the exchange is collective again) at np 2 and 4 on fix; the A/B of every FLUME
# and PRISM case base against fix, bit for bit (ab-all.sh of the #157 night). Also the Open MPI sm transport settings.
#SBATCH -J i161-gate
#SBATCH -A iacroma
#SBATCH -p prod-gn
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH -t 04:00:00
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/speedup-133/i161/%x.%j.out
D=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i161; G=$D/gate; F=/ifs/hpc/home/szaghi/adam-tests/opt/fftw-3.3.10
source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1 || exit 1
export LIBRARY_PATH=$F/lib:${LIBRARY_PATH:-} LD_LIBRARY_PATH=$F/lib:${LD_LIBRARY_PATH:-} CPATH=$F/include:${CPATH:-}
export OMPI_MCA_coll_hcoll_enable=0 PRTE_MCA_hwloc_default_binding_policy=none
echo "start $(date +%T) on $(hostname), job $SLURM_JOB_ID; $(gfortran --version | head -1)"
ompi_info --param btl sm --level 9 2>/dev/null | grep -E 'single_copy|eager_limit|max_send_size|rndv' | head -8
ompi_info --param smsc all --level 9 2>/dev/null | grep -E 'MCA smsc|priority' | head -8
rc=0
for t in base fix; do
  A=$G/$t; rm -rf "${A:?}" && mkdir -p $A && tar -xzf $D/$t-tree.tgz -C $A && tar -xzf $D/../i160/deps.tgz -C $A || exit 3
  ln -s /home/szaghi/adam/lib $A/lib; mkdir -p $A/exe; ln -s /ifs/hpc/home/szaghi/adam-tests/oracle-venv $A/exe/.regression-venv
  cd $A && git init -q && git add -A src fobos fobos.d > /dev/null 2>&1
  for m in flume-cpu-gnu flume-cpu-gnu-fft prism-cpu-gnu ripple-cpu-gnu test-ripple-crest-criterion-cpu-gnu; do
    t0=$(date +%s); fobis build --mode $m --varset iac_gnu -j 32 > $G/$t-build-$m.log 2>&1; r=$?
    echo "$(date +%T) BUILD $t $m exit $r after $(( $(date +%s) - t0 )) s"; rc=$((rc|r))
  done
done
[[ $rc -eq 0 ]] || { echo "GATE-CPU-DONE rc $rc (build)"; exit $rc; }
A=$G/fix; cd $A
find $A/src $A/fobos $A/fobos.d -type f -newermt '2026-11-18 12:00' -exec touch -h -d '2026-11-18 12:00' {} +
for app in flume prism; do
  t0=$(date +%s); REGRESSION_JOBS=8 src/tests/$app/regression/run.sh cpu --no-build --varset iac_gnu > $G/fix-cpu-$app-suite.log 2>&1; r=$?
  echo "SUITE_${app} rc $r after $(( $(date +%s) - t0 )) s: $(tail -2 $G/fix-cpu-$app-suite.log | tr '\n' ' ')"
  grep -h '^FAIL \[' $G/fix-cpu-$app-suite.log | sort -u
done
t0=$(date +%s); (cd src/tests/ripple && ./check.sh) > $G/fix-cpu-ripple.log 2>&1; echo "RIPPLE rc $? after $(( $(date +%s) - t0 )) s: $(grep -E 'PASS|FAIL' $G/fix-cpu-ripple.log | tail -3 | tr '\n' ' ')"
for np in 2 4; do
  t0=$(date +%s); (cd src/tests/flume/verification/ranks && FLUME_EXE=$A/exe/adam_flume_cpu ./check.sh --np $np --legs "repro euler mhd-none mhd-glm mhd-eglm limiter balanced hlld amr split") > $G/fix-ranks-np$np.log 2>&1
  echo "RANKS np $np rc $? after $(( $(date +%s) - t0 )) s: $(grep -c -E '^ *PASS' $G/fix-ranks-np$np.log) PASS, $(grep -c -E '^ *FAIL' $G/fix-ranks-np$np.log) FAIL; $(tail -1 $G/fix-ranks-np$np.log)"
done
t0=$(date +%s); bash $D/ab-all.sh cpu $G/base/exe $A/exe $G/ab-cpu 6 > $G/fix-ab-cpu.log 2>&1
echo "AB after $(( $(date +%s) - t0 )) s: $(grep -A50 AB-ALL $G/fix-ab-cpu.log | tr '\n' ' ')"
echo "GATE-CPU-DONE $(date +%T)"
