#!/usr/bin/env bash
# FLUME speedup matrix: run the chain of a leg (warm-up, settle, windows) on one backend and rank count (issue #133).
#
# Why: every backend and rank count must time the same grid and state (D3): each leg warms up once (its uniform case to
# t_w, a restart dump), settles once (its runtime case restarted from the warm-up dump, the first regrids derefine what
# the marker does not need, a second dump), and every window starts from one of those two dumps, repartitioned for the
# rank count that reads it (repartition.py). A window is REPETITIONS interleaved pairs of runs of WINDOW_SHORT and
# WINDOW_LONG steps, every output off, ADAM_SEAM_TIMING=1 (the forest timing lines at the end of each log); analyze.py
# turns the pairs into step costs and ratios.
#
# Stages:
#   cases     write the cases from matrix.conf into <work>/cases (make_cases.py; every other stage does it first too)
#   warmup    per leg of the selected cases: its uniform case from the cold start to t_w (--reduced: WARMUP_STEPS),
#             restart dump in <work>/dumps/leg<N>/warmup/; --use-dump DIR links an existing dump instead of running
#             (leg 3: a surviving T3 resolved-hybrid dump, the same basename channel-restart)
#   settle    per leg of the selected runtime and frozen cases: the warm-up dump repartitioned for --np (rule units),
#             the runtime case restarted from it for SETTLE_STEPS steps, dump in <work>/dumps/leg<N>/settle/
#   windows   per selected case: the source dump (uniform: warm-up; runtime and frozen: settle) repartitioned for --np
#             into <work>/windows/<backend>-np<N>/<case>/dump/ (rule siblings for the uniform and frozen grids, the
#             layout of a static cold start; units for the runtime grid, the layout of a regrid), then for each
#             repetition, for each case, the short and the long run in r<rep>-n<steps>/
#   chain     warmup, settle (when a runtime or frozen case is selected) and windows, in this order
# The warm-up and the settle run once per leg, on the fastest backend (D3); run them with the backend and rank count
# of your choice, then the windows of every backend start from their dumps. A stage that finds its directory already
# there moves it to <work>/attic/ first.
#
# Usage: ./run.sh <stage> [--backend cpu|fnl-nvf|fnl-amd] [--np N] [--cases LIST] [--work DIR] [--launcher "CMD"]
#                 [--exe PATH] [--reduced] [--use-dump DIR] [--conf FILE] [--dry-run]
#   --backend   cpu (default), fnl-nvf or fnl-amd: the executable ../exe/<backend>/ of build.sh, and the window
#               directory name
#   --np        MPI ranks (default $SLURM_NTASKS, or 2); FNL: one rank per GPU
#   --cases     comma-separated case names or prefixes: 1a, 2b, 3 (every case of leg 3), all (default)
#   --work      the run tree (default runs/full, or runs/reduced with --reduced, next to this script); on medooza a
#               directory under /ifs/hpc/home/szaghi/adam-tests, never the repository
#   --launcher  the launcher followed by the rank count (default "mpirun -x ADAM_SEAM_TIMING -np"; Open MPI forwards
#               only the variables -x names; srun forwards the environment: "srun --mpi=pmix -n")
#   --exe       the executable (default ../exe/<backend>/adam_flume_{cpu,fnl})
#   --reduced   the reduced sizes of matrix.conf (make_cases.py --reduced): the WSL plumbing test, seconds per run
#   --use-dump  warmup only: link the restart files found in DIR instead of running the warm-up
#   --conf      the matrix (default matrix.conf)
#   --dry-run   print the commands, run nothing (the cases are still written)
#
# Output: each run directory holds the inputs of the run, log.txt, wall.txt (seconds) and run.info (host, job, stage,
# case, role, repetition, steps, start step, blocks per realm at the start, cells per block, command, executable and
# its digest, commit, whether src/ differs from it). WSL: one run at a time, plumbing only (never timings), FNL np <= 2
# and the UCX knobs of issue #12 (CLAUDE.md, WSL2 GPU+MPI caveats).
set -euo pipefail

SPEEDUP_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SPEEDUP_DIR/../../../.." && pwd)"
STAGE="${1:-}"
[[ $# -gt 0 ]] && shift
BACKEND=cpu
NP="${SLURM_NTASKS:-2}"
CASES=all
WORK=""
LAUNCHER="mpirun -x ADAM_SEAM_TIMING -np"
EXE=""
REDUCED=0
USE_DUMP=""
CONF="$SPEEDUP_DIR/matrix.conf"
DRY=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --backend)  BACKEND="$2" ; shift 2 ;;
      --np)       NP="$2" ; shift 2 ;;
      --cases)    CASES="$2" ; shift 2 ;;
      --work)     WORK="$2" ; shift 2 ;;
      --launcher) LAUNCHER="$2" ; shift 2 ;;
      --exe)      EXE="$2" ; shift 2 ;;
      --reduced)  REDUCED=1 ; shift ;;
      --use-dump) USE_DUMP="$2" ; shift 2 ;;
      --conf)     CONF="$2" ; shift 2 ;;
      --dry-run)  DRY=1 ; shift ;;
      *) echo "run.sh: unknown argument '$1' (see the usage at the top of the script)" >&2 ; exit 2 ;;
   esac
done

die() { echo "run.sh: $*" >&2 ; exit 2 ; }

case "$STAGE" in cases|warmup|settle|windows|chain) ;; *) die "the stage must be cases, warmup, settle, windows or chain" ;; esac
case "$BACKEND" in
   cpu)             EXE="${EXE:-$SPEEDUP_DIR/exe/cpu/adam_flume_cpu}" ;;
   fnl-nvf|fnl-amd) EXE="${EXE:-$SPEEDUP_DIR/exe/$BACKEND/adam_flume_fnl}" ;;
   *) die "--backend must be cpu, fnl-nvf or fnl-amd, got '$BACKEND'" ;;
esac
[[ $NP =~ ^[1-9][0-9]*$ ]] || die "--np must be a positive integer, got '$NP'"
[[ -n $USE_DUMP && $STAGE != warmup ]] && die "--use-dump belongs to the warmup stage"
if [[ -z $WORK ]]; then
   if [[ $REDUCED -eq 1 ]]; then WORK="$SPEEDUP_DIR/runs/reduced" ; else WORK="$SPEEDUP_DIR/runs/full" ; fi
fi
mkdir -p "$WORK"
WORK="$(cd "$WORK" && pwd)"
if [[ $DRY -eq 0 && $STAGE != cases && ! -x $EXE ]]; then
   die "executable '$EXE' not found (./build.sh $BACKEND, or --exe)"
fi
export ADAM_SEAM_TIMING=1
PY=(python3 -I)

# --- the cases ---------------------------------------------------------------------------------------------------------
MAKE_ARGS=("$WORK/cases" --conf "$CONF")
[[ $REDUCED -eq 1 ]] && MAKE_ARGS+=(--reduced)
"${PY[@]}" "$SPEEDUP_DIR/make_cases.py" "${MAKE_ARGS[@]}" > /dev/null
cp "$CONF" "$WORK/cases/matrix.conf"
mapfile -t ALL_CASES < <(find "$WORK/cases" -mindepth 2 -maxdepth 2 -name case.conf -printf '%h\n' | xargs -n1 basename | sort)
SELECTED=()
IFS=',' read -r -a TOKENS <<< "$CASES"
for c in "${ALL_CASES[@]}"; do
   for t in "${TOKENS[@]}"; do
      if [[ $t == all || $c == "$t"* ]]; then SELECTED+=("$c") ; break ; fi
   done
done
[[ ${#SELECTED[@]} -gt 0 ]] || die "--cases '$CASES' selects none of: ${ALL_CASES[*]}"
if [[ $STAGE == cases ]]; then
   printf 'run.sh: cases in %s:\n' "$WORK/cases"
   printf '   %s\n' "${ALL_CASES[@]}"
   exit 0
fi

declare -A C # the keys of one case.conf
load_case() {
   local file="$WORK/cases/$1/case.conf" key value
   [[ -f $file ]] || die "no case '$1' ($file)"
   C=()
   while IFS='=' read -r key value; do
      [[ $key =~ ^[A-Z_0-9]+$ ]] || continue
      eval "C[\$key]=$value" # make_cases.py writes shell-quoted values
   done < "$file"
}

# the case of a leg with a role (uniform, runtime, frozen): <leg> <role>
case_of() {
   local c
   for c in "${ALL_CASES[@]}"; do
      load_case "$c"
      if [[ ${C[LEG]} == "$1" && ${C[ROLE]} == "$2" ]]; then echo "$c" ; return 0 ; fi
   done
   die "no $2 case in leg $1"
}

# the step of a restart (int32 step, float64 time in <basename>.time), and the block count of its tree (int64 in .tnd)
dump_step() { "${PY[@]}" -c 'import struct, sys; print(struct.unpack("=id", open(sys.argv[1], "rb").read(12))[0])' "$1.time" ; }
dump_blocks() { "${PY[@]}" -c 'import struct, sys; print(struct.unpack("=q", open(sys.argv[1], "rb").read(8))[0])' "$1.tnd" ; }

# move an existing directory aside, into <work>/attic/
fresh_dir() {
   local dir="$1" rel
   if [[ -e $dir ]]; then
      rel="${dir#"$WORK"/}"
      mkdir -p "$WORK/attic/$(dirname "$rel")"
      mv "$dir" "$WORK/attic/$rel-$(date +%Y%m%d-%H%M%S)"
   fi
   mkdir -p "$dir"
}

# copy the inputs of the loaded case into a run directory
copy_inputs() {
   local dir="$1" ini
   cp "$WORK/cases/${C[CASE]}/${C[INPUT]}" "$dir/"
   for ini in ${C[INIS]}; do cp "$WORK/cases/${C[CASE]}/$ini" "$dir/" ; done
}

# rewrite options (section.key=value) of every realm input of the loaded case in a run directory
set_options() {
   local dir="$1" ini
   shift
   for ini in ${C[INIS]}; do "${PY[@]}" "$SPEEDUP_DIR/ini_set.py" "$dir/$ini" "$@" ; done
}

# repartition the dump of every realm of the loaded case: <src dir> <dst dir> <rule>
repartition() {
   local src="$1" dst="$2" rule="$3" b paths=()
   for b in ${C[BASENAMES]}; do paths+=("$src/$b") ; done
   "${PY[@]}" "$SPEEDUP_DIR/repartition.py" "${paths[@]}" --to "$NP" --out "$dst" --rule "$rule" --ratio "${C[RATIO]}" \
      --force | sed 's/^/   /'
}

# the step of the dump of the loaded case in a directory, checked equal across the realms
case_dump_step() {
   local dir="$1" b s step=""
   for b in ${C[BASENAMES]}; do
      [[ -f $dir/$b.time ]] || die "no restart dump $dir/$b.time"
      s="$(dump_step "$dir/$b")"
      [[ -z $step || $s == "$step" ]] || die "$dir: the realms' dumps are at different steps ($step, $s)"
      step="$s"
   done
   echo "$step"
}

# run the loaded case in a directory: <dir> <stage> [run.info lines...]
run_case() {
   local dir="$1" stage="$2" cmd status t0 t1
   shift 2
   cmd="$LAUNCHER $NP $EXE ${C[INPUT]}"
   if [[ $DRY -eq 1 ]]; then
      echo "   cd $dir && ADAM_SEAM_TIMING=1 $cmd"
      return 0
   fi
   {
      echo "date=$(date -Is)"
      echo "host=$(hostname)"
      echo "slurm_job=${SLURM_JOB_ID:-}"
      echo "stage=$stage"
      echo "backend=$BACKEND"
      echo "np=$NP"
      echo "case=${C[CASE]}"
      echo "leg=${C[LEG]}"
      echo "role=${C[ROLE]}"
      echo "reduced=${C[REDUCED]}"
      echo "block_cells=${C[BLOCK_CELLS]// /,}"
      echo "regrid_period=${C[REGRID_PERIOD]}"
      printf '%s\n' "$@"
      echo "command=$cmd"
      echo "exe=$EXE"
      echo "exe_sha256=$(sha256sum "$EXE" | cut -c1-16)"
      echo "commit=$(git -C "$REPO_ROOT" rev-parse HEAD 2> /dev/null || echo unknown)"
      # the sources differ from that commit (the scripts under src/tests do not enter the executable)
      if git -C "$REPO_ROOT" diff --quiet HEAD -- src ':!src/tests' 2> /dev/null; then
         echo "src_modified=no"
      else
         echo "src_modified=yes"
      fi
   } > "$dir/run.info"
   t0=$(date +%s.%N)
   status=0
   (cd "$dir" && eval "$cmd" > log.txt 2>&1) || status=$?
   t1=$(date +%s.%N)
   awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f\n", b - a}' > "$dir/wall.txt"
   echo "   ${dir#"$WORK"/}: wall $(cat "$dir/wall.txt") s, exit $status"
   grep -a "forest timing" "$dir/log.txt" | sed 's/^\[mpi-[0-9]*\]/      /' || true
   [[ $status -eq 0 ]] || die "the run exited $status, see $dir/log.txt"
}

# the legs of the selected cases, and whether a selected case of a leg needs the settle dump
declare -A LEGS=() SETTLE=()
for c in "${SELECTED[@]}"; do
   load_case "$c"
   LEGS[${C[LEG]}]=1
   [[ ${C[SOURCE]} == settle ]] && SETTLE[${C[LEG]}]=1
done

stage_warmup() {
   local leg="$1" dir b f n
   load_case "$(case_of "$leg" uniform)"
   dir="$WORK/dumps/leg$leg/warmup"
   echo ">> warm-up, leg $leg: ${C[CASE]} -> ${dir#"$WORK"/}"
   if [[ $DRY -eq 0 ]]; then fresh_dir "$dir" ; fi
   if [[ -n $USE_DUMP ]]; then
      for b in ${C[BASENAMES]}; do
         [[ -f $USE_DUMP/$b.time ]] || die "--use-dump: no $USE_DUMP/$b.time"
         n=0
         for f in "$USE_DUMP/$b".* "$USE_DUMP/$b"-proc*; do
            [[ -e $f ]] || continue
            if [[ $DRY -eq 1 ]]; then echo "   ln -s $f $dir/" ; else ln -s "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" "$dir/" ; fi
            n=$((n + 1))
         done
         echo "   $b: $n files linked from $USE_DUMP"
      done
      return 0
   fi
   local opts=("IO.restart=.false." "IO.restart_save=${C[HUGE]}" "IO.it_save=${C[HUGE]}")
   if [[ ${C[WARMUP_STEPS]} -gt 0 ]]; then
      opts+=("time.it_max=${C[WARMUP_STEPS]}")
   else
      opts+=("time.it_max=-1" "time.time_max=${C[T_W]}")
   fi
   if [[ $DRY -eq 0 ]]; then copy_inputs "$dir" ; set_options "$dir" "${opts[@]}" ; fi
   run_case "$dir" warmup "steps=${C[WARMUP_STEPS]}" "t_w=${C[T_W]}"
   if [[ $DRY -eq 0 ]]; then echo "   dump at step $(case_dump_step "$dir")" ; fi
}

stage_settle() {
   local leg="$1" warm dir s_w it_max
   load_case "$(case_of "$leg" runtime)"
   warm="$WORK/dumps/leg$leg/warmup"
   dir="$WORK/dumps/leg$leg/settle"
   echo ">> settle, leg $leg: ${C[CASE]}, $NP ranks -> ${dir#"$WORK"/}"
   if [[ $DRY -eq 1 ]]; then
      echo "   repartition.py <warm-up dump> --to $NP --out $dir --rule units"
      echo "   it_max = <warm-up step> + ${C[SETTLE_STEPS]}"
      run_case "$dir" settle
      return 0
   fi
   s_w="$(case_dump_step "$warm")"
   it_max=$((s_w + C[SETTLE_STEPS]))
   fresh_dir "$dir"
   copy_inputs "$dir"
   repartition "$warm" "$dir" units
   set_options "$dir" "IO.restart=.true." "IO.restart_save=${C[HUGE]}" "IO.it_save=${C[HUGE]}" "time.it_max=$it_max"
   run_case "$dir" settle "start_step=$s_w" "steps=${C[SETTLE_STEPS]}"
   [[ $(case_dump_step "$dir") == "$it_max" ]] || die "$dir: the settle dump is not at step $it_max"
   echo "   dump at step $it_max"
}

stage_windows() {
   local root="$WORK/windows/$BACKEND-np$NP" c wdir s0 b blocks rep n rdir f
   echo ">> windows, $BACKEND, $NP ranks: ${SELECTED[*]} -> ${root#"$WORK"/}"
   for c in "${SELECTED[@]}"; do # the dump of each window, repartitioned for NP
      load_case "$c"
      wdir="$root/$c"
      if [[ $DRY -eq 1 ]]; then
         echo "   $c: repartition.py <dumps/leg${C[LEG]}/${C[SOURCE]}> --to $NP --out $wdir/dump --rule ${C[RULE]}"
         continue
      fi
      fresh_dir "$wdir"
      repartition "$WORK/dumps/leg${C[LEG]}/${C[SOURCE]}" "$wdir/dump" "${C[RULE]}"
      s0="$(case_dump_step "$wdir/dump")"
      blocks=""
      for b in ${C[BASENAMES]}; do blocks+="${blocks:+,}$(dump_blocks "$wdir/dump/$b")" ; done
      printf 'S0=%s\nREALM_BLOCKS=%s\n' "$s0" "$blocks" > "$wdir/window.conf"
   done
   for ((rep = 1; rep <= C[REPETITIONS]; rep++)); do
      for c in "${SELECTED[@]}"; do
         load_case "$c"
         wdir="$root/$c"
         s0="<dump step>" ; blocks="?"
         if [[ $DRY -eq 0 ]]; then
            s0="$(sed -n 's/^S0=//p' "$wdir/window.conf")"
            blocks="$(sed -n 's/^REALM_BLOCKS=//p' "$wdir/window.conf")"
         fi
         for n in "${C[WINDOW_SHORT]}" "${C[WINDOW_LONG]}"; do
            rdir="$wdir/r$rep-n$n"
            if [[ $DRY -eq 0 ]]; then
               mkdir -p "$rdir"
               copy_inputs "$rdir"
               for b in ${C[BASENAMES]}; do
                  for f in "$wdir/dump/$b".* "$wdir/dump/$b"-proc*; do [[ -e $f ]] && ln -s "../dump/$(basename "$f")" "$rdir/" ; done
               done
               # shellcheck disable=SC2086 # OFF_KEYS is a list of section.key=0 words
               set_options "$rdir" "IO.restart=.true." "time.it_max=$((s0 + n))" ${C[OFF_KEYS]}
            fi
            run_case "$rdir" window "rep=$rep" "steps=$n" "start_step=$s0" "realm_blocks=$blocks" \
               "window_short=${C[WINDOW_SHORT]}" "window_long=${C[WINDOW_LONG]}" "spread_warning=${C[SPREAD_WARNING]}"
         done
      done
   done
}

if [[ $STAGE == warmup || $STAGE == chain ]]; then
   for leg in $(printf '%s\n' "${!LEGS[@]}" | sort -n); do stage_warmup "$leg" ; done
fi
if [[ $STAGE == settle || $STAGE == chain ]]; then
   [[ ${#SETTLE[@]} -gt 0 || $STAGE == chain ]] || die "settle: no runtime or frozen case selected (--cases)"
   for leg in $(printf '%s\n' "${!SETTLE[@]}" | sort -n); do stage_settle "$leg" ; done
fi
if [[ $STAGE == windows || $STAGE == chain ]]; then
   stage_windows
fi
echo "run.sh: done ($STAGE, $BACKEND, $NP ranks, work $WORK)"
