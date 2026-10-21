#!/usr/bin/env bash
# FLUME speedup matrix: run the measurements (data/) and the full simulations (full/) of one leg.
#
# Why this protocol: FLUME has no step timer, and a run's wall time also holds its start-up (MPI, the device, the
# initial AMR, the distance function of the immersed body) and its output. The measurement removes both:
#   1. warmup: the case runs to T_WARM with every output off but one restart dump at its end. All legs restart from
#      the same physical state (each from its own dump: the restart files are per rank);
#   2. window: REPEATS times, interleaved, a run restarted from the dump for N1 steps ("short") and one for N1 + NWIN
#      steps ("long"), every output off (it_save = restart_save = residuals_save = conservation_history_save = 0, so
#      not even the last step is written). Both pay the same start-up, restart read and first-step transients, so
#      (T_long - T_short) / NWIN is the cost of a step on the developed flow. N1 and NWIN are multiples of the regrid
#      cadence and the cadence keys on the step count (restored by the restart), so the two runs differ by exactly
#      NWIN / AMR_FREQUENCY regrids, which belong to the cost of a dynamic-AMR step.
# The window is NWIN steps of a flow at T_WARM: a few percent of the run, so the grid and the solution barely change
# across it and every leg measures the same work.
#
# Every run must exit 0 (medooza CPU: export OMPI_MCA_coll_hcoll_enable=0, the gnu.sh toolchain does); a failed run is
# recorded with its status and analyze.py ignores its repetition.
#
# Usage: ./run.sh <leg> [--np N] [--cases LIST] [--what data|full|all] [--launcher "CMD"] [--exe PATH] [--reuse-warmup]
#   leg             cpu, fnl-nvf or fnl-amd (the executable is exe/<leg>/, see build.sh)
#   --np            MPI ranks (default $SLURM_NTASKS); FNL: one rank per GPU
#   --cases         case names or substrings, commas or blanks (default: every case); e.g. --cases cylinder-amr,sphere
#   --what          data (default), full (only the FULL_CASES of matrix.conf), or all
#   --launcher      the launcher, followed by the rank count (default "mpirun -np"; medooza: "srun --mpi=pmix -n")
#   --reuse-warmup  keep an existing warmup dump (same leg and np) instead of running it again
set -euo pipefail

SPEEDUP_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SPEEDUP_DIR/../../../.." && pwd)"
TOOL="$SPEEDUP_DIR/make_inputs.py"
# shellcheck source=matrix.conf
source "$SPEEDUP_DIR/matrix.conf"

LEG="${1:-}"
[[ $# -gt 0 ]] && shift
NP="${SLURM_NTASKS:-}"
SELECT=""
WHAT=data
LAUNCHER="mpirun -np"
EXE=""
REUSE=0
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)           NP="$2" ; shift 2 ;;
      --cases)        SELECT+=" ${2//,/ }" ; shift 2 ;;
      --what)         WHAT="$2" ; shift 2 ;;
      --launcher)     LAUNCHER="$2" ; shift 2 ;;
      --exe)          EXE="$2" ; shift 2 ;;
      --reuse-warmup) REUSE=1 ; shift ;;
      *) echo "run.sh: unknown argument '$1' (see the usage at the top of the script)" >&2 ; exit 2 ;;
   esac
done
case "$LEG" in
   cpu)             EXE="${EXE:-$SPEEDUP_DIR/exe/cpu/adam_flume_cpu}" ;;
   fnl-nvf|fnl-amd) EXE="${EXE:-$SPEEDUP_DIR/exe/$LEG/adam_flume_fnl}" ;;
   *) echo "run.sh: the leg must be cpu, fnl-nvf or fnl-amd, got '$LEG'" >&2 ; exit 2 ;;
esac
case "$WHAT" in data|full|all) ;; *) echo "run.sh: --what must be data, full or all" >&2 ; exit 2 ;; esac
[[ -n $NP ]] || { echo "run.sh: give --np (or run inside a Slurm allocation)" >&2 ; exit 2 ; }
[[ -x $EXE ]] || { echo "run.sh: executable '$EXE' not found (./build.sh $LEG, or --exe)" >&2 ; exit 2 ; }
if (( N1 % AMR_FREQUENCY || NWIN % AMR_FREQUENCY )); then
   echo "run.sh: N1 = $N1 and NWIN = $NWIN must be multiples of AMR_FREQUENCY = $AMR_FREQUENCY" >&2
   exit 2
fi

CASES=()
while read -r c; do
   if [[ -z ${SELECT// } ]]; then CASES+=("$c"); continue; fi
   for s in $SELECT; do
      if [[ $c == *"$s"* ]]; then CASES+=("$c"); break; fi
   done
done < <(python3 "$TOOL" cases)
[[ ${#CASES[@]} -gt 0 ]] || { echo "run.sh: --cases '$SELECT' selects no case" >&2 ; exit 2 ; }
LABEL="$LEG-np$NP"
FAILED=0

ini_value() { # ini_value <ini> <section> <key>: the value of an option (comments stripped)
   awk -F= -v s="[$2]" -v k="$3" '/^\[/{sec=$0} sec==s && $1~"^ *"k" *$" {gsub(/[ \t]|;.*/, "", $2); print $2; exit}' "$1"
}
launch() { # launch <dir> <csv> <rep> <kind> <steps>: run the input of <dir>, time it, append a row to <csv>
   local dir=$1 csv=$2 t0 t1 status=0
   t0=$(date +%s.%N)
   (cd "$dir" && $LAUNCHER "$NP" "$EXE" input.ini > log.txt 2>&1) || status=$?
   t1=$(date +%s.%N)
   echo "$3,$4,$5,$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b - a}'),$status" >> "$csv"
   if [[ $status -ne 0 ]]; then
      echo "run.sh: $dir exited $status, see $dir/log.txt" >&2
      FAILED=1
   fi
}
run_info() { # run_info <file> <key=value>...: what the timings depend on
   {
      echo "date=$(date -Is)"
      echo "host=$(hostname)"
      echo "slurm_job=${SLURM_JOB_ID:-}"
      echo "leg=$LEG"
      echo "np=$NP"
      echo "launcher=$LAUNCHER"
      echo "exe=$EXE"
      echo "exe_sha256=$(sha256sum "$EXE" | cut -c1-16)"
      echo "commit=$(git -C "$REPO_ROOT" rev-parse HEAD 2> /dev/null || echo unknown)"
      if command -v nvidia-smi > /dev/null; then
         echo "gpu=$(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader | sort | uniq -c | xargs)"
      fi
      if command -v rocm-smi > /dev/null; then
         echo "gpu=$(rocm-smi --showproductname 2> /dev/null | grep -m1 -i 'card series' | xargs)"
      fi
      printf '%s\n' "$@"
   } > "$1"
}

for case in "${CASES[@]}"; do
   cdir="$SPEEDUP_DIR/$case"
   if grep -qw TBD "$cdir/input.ini"; then
      echo "run.sh: $case/input.ini still holds TBD: set the grid in matrix.conf, then ./make_inputs.py inputs" >&2
      exit 2
   fi
   rb="$(ini_value "$cdir/input.ini" IO restart_basename)"

   if [[ $WHAT != full ]]; then
      d="$cdir/data/$LABEL"
      csv="$d/timings.csv"
      echo ">> $case, $LABEL: warmup to t = $T_WARM, then $REPEATS x ($N1, $N1 + $NWIN) steps"
      mkdir -p "$d"
      if [[ $REUSE -eq 0 || ! -f "$d/warmup/$rb.time" ]]; then
         rm -rf "$d"/warmup "$d"/rep*-short "$d"/rep*-long
         echo "rep,kind,steps,seconds,status" > "$csv"
         mkdir -p "$d/warmup"
         cp "$cdir/input.ini" "$d/warmup/input.ini"
         python3 "$TOOL" set "$d/warmup/input.ini" IO.it_save=0 IO.restart_save=1000000000 IO.residuals_save=0 \
            diagnostics.conservation_history_save=0 time.it_max=-1 "time.time_max=$T_WARM"
         launch "$d/warmup" "$csv" 0 warmup 0
         # the warmup also dumps the restart at step 0: a warmup that died later leaves that one behind
         if [[ "$(tail -n1 "$csv" | cut -d, -f5)" != 0 || ! -f "$d/warmup/$rb.time" ]]; then
            echo "run.sh: the warmup of $case failed or wrote no restart, skipping its window" >&2
            FAILED=1
            continue
         fi
      else
         rm -rf "$d"/rep*-short "$d"/rep*-long
         grep -E '^(rep|0,warmup)' "$csv" > "$csv.tmp" && mv "$csv.tmp" "$csv"
      fi
      read -r it_w t_w < <(python3 "$TOOL" restart-step "$d/warmup/$rb.time")
      run_info "$d/run.info" "it_warm=$it_w" "t_warm=$t_w" "n1=$N1" "nwin=$NWIN" "repeats=$REPEATS" \
         "amr_frequency=$AMR_FREQUENCY"
      for (( r = 1; r <= REPEATS; r++ )); do
         for kind in short long; do
            steps=$N1
            [[ $kind == long ]] && steps=$(( N1 + NWIN ))
            w="$d/rep$r-$kind"
            mkdir -p "$w"
            cp "$d/warmup/input.ini" "$w/input.ini"
            python3 "$TOOL" set "$w/input.ini" IO.restart=.true. IO.restart_save=0 "time.it_max=$(( it_w + steps ))" \
               time.time_max=1.0e30
            for f in "$d/warmup/$rb".* "$d/warmup/$rb"-proc*; do
               if [[ -e $f ]]; then ln -sf "$f" "$w/"; fi
            done
            launch "$w" "$csv" "$r" "$kind" "$steps"
            echo "   rep $r $kind ($steps steps): $(tail -n1 "$csv" | cut -d, -f4) s"
         done
      done
   fi

   if [[ $WHAT != data && " $FULL_CASES " == *" $case "* ]]; then
      w="$cdir/full/$LABEL"
      echo ">> $case, $LABEL: full simulation to t = $T_END"
      rm -rf "$w"
      mkdir -p "$w"
      cp "$cdir/input.ini" "$w/input.ini"
      echo "rep,kind,steps,seconds,status" > "$w/timings.csv"
      run_info "$w/run.info" "t_end=$T_END"
      launch "$w" "$w/timings.csv" 0 full 0
      echo "   $(tail -n1 "$w/timings.csv" | cut -d, -f4) s"
   fi
done

if [[ $FAILED -eq 0 ]]; then
   echo "run.sh: $LABEL done (${#CASES[@]} cases, --what $WHAT); next: ./analyze.sh"
else
   echo "run.sh: $LABEL had failed runs (see the messages above)" >&2
   exit 1
fi
