#!/bin/bash
# Night 2026-11-26-i157-complete (#157 B4 fix, step 07): the #133 step windows (leg 1, the only leg with seam
# packs) re-run with HEAD 80582712 (A), the s03 build (B: B4 allocating its guard scratch per pack) and the s07 build (C:
# the persistent scratch), alternated run by run on the same node. Reads the matrix windows read only; writes under $R.
W=/ifs/hpc/home/szaghi/adam-tests/speedup-133/work
WIN=$W/windows
R=/ifs/hpc/home/szaghi/adam-tests/speedup-133/i157d
declare -A EXE=([A]=/ifs/hpc/home/szaghi/adam-tests/i157n/ab-base/adam_flume_fnl [B]=/ifs/hpc/home/szaghi/adam-tests/i157n/s03-fnl/exe/adam_flume_fnl [C]=/ifs/hpc/home/szaghi/adam-tests/i157n/s07-fnl/exe/adam_flume_fnl)
declare -A COMMIT=([A]=80582712 [B]=i157n-s03 [C]=i157n-s07)
log() { echo "$(date +%T) $*"; }
fnl_env() { source /ifs/hpc/home/szaghi/nvidia/nvhpc.sh ompi5 > /dev/null 2>&1; export OMPI_MCA_coll_ucc_enable=0 OMPI_MCA_coll=^ucc; }
# win <A|B> <np> <case> <rep> <n> [VAR=value ...]: one window run into $R/<tag>/windows/fnl-nvf-np<np>/<case>/r<rep>-n<n>
win() {
   local tag=$1 np=$2 case=$3 rep=$4 n=$5 src d inp=input.ini t0 t1 rc f exe
   shift 5
   src=$WIN/fnl-nvf-np$np/$case; d=$R/$tag/windows/fnl-nvf-np$np/$case/r$rep-n$n; exe=${EXE[$tag]}
   rm -rf "$d"; mkdir -p "$d"; cp "$src/r$rep-n$n"/*.ini "$d"/
   for f in "$src"/dump/*-restart*; do ln -s "$f" "$d"/; done
   [[ -f $d/step.ini ]] && inp=step.ini
   t0=$(date +%s.%N)
   ( fnl_env; cd "$d" && env ADAM_SEAM_TIMING=1 "$@" srun --mpi=pmix -n "$np" "$exe" "$inp" > log.txt 2>&1 ); rc=$?
   t1=$(date +%s.%N)
   awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f\n", b - a}' > "$d/wall.txt"
   grep -vE '^(date|host|slurm_job|command|exe|exe_sha256|commit|src_modified)=' "$src/r$rep-n$n/run.info" > "$d/run.info"
   printf 'date=%s\nhost=%s\nslurm_job=%s\ncommand=srun --mpi=pmix -n %s %s %s\nexe=%s\nexe_sha256=%s\ncommit=%s\nsrc_modified=no\nenv=%s\nrc=%s\n' \
      "$(date -Is)" "$(hostname)" "$SLURM_JOB_ID" "$np" "$exe" "$inp" "$exe" "$(sha256sum "$exe" | cut -c1-16)" \
      "${COMMIT[$tag]}" "$*" "$rc" >> "$d/run.info"
   log "$tag $case np$np r$rep-n$n: rc $rc wall $(cat "$d/wall.txt") $(grep -a 'forest timing' "$d/log.txt" | head -1 | awk '{print $7, $8, $9, $10}')"
   return $rc
}
# windows <np> <case...>: 3 repetitions x (10, 110) steps, A and B alternated
windows() {
   local np=$1 c rep n tag
   shift
   for c in "$@"; do
      for rep in 1 2 3; do for n in 10 110; do for tag in A B C; do win $tag "$np" "$c" $rep $n; done; done; done
   done
}
prepare() {
   mkdir -p "$R/A" "$R/B" "$R/C"
   ln -sfn "$W/cases" "$R/A/cases"; ln -sfn "$W/cases" "$R/B/cases"; ln -sfn "$W/cases" "$R/C/cases"
   log "start on $(hostname), job $SLURM_JOB_ID"
   nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.free --format=csv,noheader
   while IFS= read -r free; do
      [[ ${free// /} -lt 20000 ]] && { log "a GPU has only $free MiB free: resubmit with -x $(hostname)"; exit 4; }
   done < <(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits)
   n=0; until [[ -f /ifs/hpc/home/szaghi/adam-tests/i157n/s07-fnl/BUILT && -f /ifs/hpc/home/szaghi/adam-tests/i157n/ab-base/READY-fnl ]] || (( n > 480 )); do sleep 60; n=$((n+1)); done
   for t in A B C; do [[ -x ${EXE[$t]} ]] || { log "no executable ${EXE[$t]}"; exit 3; }; done
}
