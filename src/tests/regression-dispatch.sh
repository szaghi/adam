#!/usr/bin/env bash
# Run the cases of a regression suite concurrently (issue #80 P5), sourced by src/tests/{flume,prism}/regression/run.sh
# after its lints, digest venv and build, when REGRESSION_JOBS > 1:
#
#   regression_dispatch <run.sh> <backend> <case>...
#
# REGRESSION_JOBS workers take the cases from one queue; each case is a child `run.sh <backend> --no-build` restricted to
# it (REGRESSION_CASES=<case>, REGRESSION_CHILD=1: no lints, no build). The cases already run in their own work
# directories, so the children share nothing but the executables. On fnl each worker owns its own GPUs
# (REGRESSION_GPUS_PER_JOB, default 2: a case runs 2 ranks; CUDA_VISIBLE_DEVICES set per worker, from
# REGRESSION_GPUS, default every GPU nvidia-smi lists), and there are never more workers than GPU groups. MPI core
# binding is off in the children, else every concurrent `mpirun -np 2` binds its ranks to cores 0 and 1 and the cases
# pile onto two cores (measured on medooza: 8 workers ran each case ~8x slower, 18 hit the timeout):
# PRTE_MCA_hwloc_default_binding_policy=none for Open MPI 5, OMPI_MCA_hwloc_base_binding_policy=none for Open MPI 4
# (the HPC-X of the NVIDIA HPC SDK). Under Slurm every mpirun is a job step: allocate 2 tasks per worker. The memory
# budget is capped by the harness (ADAM_MEMORY_FRACTION), so concurrent CPU cases fit a node.
# At the end every case's log is printed in the case order, then one summary; the exit status is 1 on any failure.

regression_dispatch() {
   local runsh="$1" backend="$2"; shift 2
   local -a cases=("$@")
   local jobs="${REGRESSION_JOBS:-1}" per="${REGRESSION_GPUS_PER_JOB:-2}"
   local -a groups=()
   if [[ "$backend" == "fnl" ]]; then
      local -a gpus=()
      if [[ -n "${REGRESSION_GPUS:-}" ]]; then
         read -r -a gpus <<< "${REGRESSION_GPUS//,/ }"
      else
         mapfile -t gpus < <(nvidia-smi --query-gpu=index --format=csv,noheader 2>/dev/null)
      fi
      local i
      for ((i = 0; i + per <= ${#gpus[@]}; i += per)); do
         groups+=("$(IFS=,; echo "${gpus[*]:i:per}")")
      done
      [[ ${#groups[@]} -eq 0 ]] && { echo "ERROR: REGRESSION_JOBS=$jobs needs groups of $per GPU(s); found ${#gpus[@]}" >&2; return 2; }
      (( jobs > ${#groups[@]} )) && jobs=${#groups[@]}
   fi
   local logdir; logdir="$(mktemp -d "${TMPDIR:-/tmp}/regression-XXXXXX")"
   printf '%s\n' "${cases[@]}" > "$logdir/queue"
   echo ">> running ${#cases[@]} case(s) on $jobs worker(s)${groups:+ (GPU groups: ${groups[*]:0:jobs})}, logs in $logdir"
   local t0; t0=$(date +%s)
   local w
   for ((w = 0; w < jobs; w++)); do
      (
         while :; do
            local c
            c="$(flock "$logdir/queue.lock" bash -c 'q="$1"; head -n 1 "$q"; sed -i 1d "$q"' _ "$logdir/queue")"
            [[ -z "$c" ]] && break
            local -a env_case=(REGRESSION_CHILD=1 REGRESSION_JOBS=1 REGRESSION_CASES="$c" PRTE_MCA_hwloc_default_binding_policy=none
                                OMPI_MCA_hwloc_base_binding_policy=none)
            [[ ${#groups[@]} -gt 0 ]] && env_case+=(CUDA_VISIBLE_DEVICES="${groups[w]}")
            env "${env_case[@]}" "$runsh" "$backend" --no-build > "$logdir/$c.log" 2>&1
            echo "$?" > "$logdir/$c.rc"
            echo ">> [$c/$backend] done (worker $w, rc $(cat "$logdir/$c.rc"), $(( $(date +%s) - t0 )) s since start)"
         done
      ) &
   done
   wait
   local pass=0 fail=0 skip=0 c n
   local -a failed=()
   # one verdict per case (a child's own SKIP count also counts the directories that are no case)
   for c in "${cases[@]}"; do
      cat "$logdir/$c.log"
      n=$(sed -n 's/^PASS: //p' "$logdir/$c.log" | tail -n 1)
      if [[ "$(cat "$logdir/$c.rc" 2>/dev/null)" != "0" ]]; then
         fail=$((fail + 1)); failed+=("$c")
      elif [[ "${n:-0}" -gt 0 ]]; then
         pass=$((pass + 1))
      else
         skip=$((skip + 1))
      fi
   done
   echo
   echo "============================================================"
   echo "== Summary [$backend] ($jobs workers, $(( $(date +%s) - t0 )) s)"
   echo "============================================================"
   echo "PASS: $pass"
   echo "FAIL: $fail"
   echo "SKIP: $skip"
   if [[ $fail -gt 0 ]]; then
      echo "Failed cases: ${failed[*]}"
      return 1
   fi
   return 0
}
