#!/bin/bash
#SBATCH -A iacroma -p prod-gn -N 1 --gres=gpu:4 -n 4 -c 8 --mem=120G -t 04:00:00 -J i156-gate
#SBATCH -o /ifs/hpc/home/szaghi/adam-tests/i156/gate-ib.%j.out
# Night 2026-11-25-i156-l2, step 01 gate: the FLUME regression case with an immersed body (shock-cylinder-ib: IB +
# solid AMR marker) on both backends against the committed goldens; every other case skips the changed branch.
# Fresh tree (the night working tree's tracked files + third_party, no exe/), so the builds are clean.
set -uo pipefail
G=/ifs/hpc/home/szaghi/adam-tests/i156/adam
# the removed-keys lint lists the committed inputs with git ls-files: a throwaway index of the scratch tree's inputs
[[ -d $G/.git ]] || ( cd $G && git init -q && find src/tests -name '*.ini' -print0 | xargs -0 git add -f -- )
cd $G/src/tests/flume/regression
export REGRESSION_CASES="shock-cylinder-ib"
echo "=== $(date -Is) CPU on $(hostname)"
( source /ifs/hpc/home/szaghi/gnu/gnu.sh > /dev/null 2>&1; ./run.sh cpu --varset iac_gnu ); echo "CPU-RC $?"
echo "=== $(date -Is) FNL"
( source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1; export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc
  ./run.sh fnl --varset iac_nvf ); echo "FNL-RC $?"
echo "GATE-DONE $(date -Is)"
