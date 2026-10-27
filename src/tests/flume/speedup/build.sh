#!/usr/bin/env bash
# FLUME speedup matrix: build the executable of one leg from a clean tree, keep it, and record the toolchain.
#
# Why: a speedup is only as good as the record of what was built. Each leg builds from scratch (fobis clean first: the
# apps share obj/mod per compiler, and an incremental build after a library change has left stale .mod files before),
# copies the executable to exe/<leg>/ so that a later build of another mode cannot replace it under a running matrix
# (both FNL legs build exe/adam_flume_fnl), and writes builds/<leg>/build.txt (the full fobis log) and
# builds/<leg>/toolchain.md (versions, varset, flags, commit, executable digest), which every case README embeds.
#
# Legs:
#   cpu      flume-cpu-gnu       gfortran + Open MPI, MPI only          e.g. --varset iac_gnu (medooza), local_gnu (WSL)
#   fnl-nvf  flume-fnl-nvf       nvfortran OpenACC, NVIDIA GPUs         e.g. --varset iac_nvf (medooza), local_nvf (WSL)
#   fnl-amd  flume-fnl-omp-amd   amdflang OpenMP offload, AMD GPUs      e.g. --varset "local_amd gfx90a"
#
# The caller loads the toolchain first (medooza: source /ifs/hpc/home/szaghi/gnu/gnu.sh or .../nvidia/nvhpc.sh); the
# script records whatever is on PATH.
#
# Variants (profiling experiments, never the matrix legs): --tag NAME builds into exe/<leg>-<NAME>/ and records into
# builds/<leg>-<NAME>/, and leaves the case READMEs alone, so a variant never replaces the executable or the toolchain
# record of a leg. --gpu OPTS (fnl-nvf only, needs --tag) appends -gpu=OPTS to the compile and link lines of the mode
# (fobis appends --cflags/--lflags to the mode's flags; nvfortran merges the -gpu= options), e.g. --gpu lto or
# --gpu maxregcount:128: the fobos files and the default build stay as they are. --gpu refuses --no-clean (objects
# keep the flags they were compiled with). The options of this script are checked before the clean.
#
# Usage: ./build.sh <leg> [--varset "V ..."] [--no-clean] [--jobs N] [--tag NAME [--gpu OPTS]]
set -euo pipefail

SPEEDUP_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SPEEDUP_DIR/../../../.." && pwd)"
LEG="${1:-}"
[[ $# -gt 0 ]] && shift
VARSET=""
CLEAN=1
JOBS=""
TAG=""
GPU_OPTS=""
while [[ $# -gt 0 ]]; do
   case "$1" in
      --varset)   VARSET="$2" ; shift 2 ;;
      --no-clean) CLEAN=0 ; shift ;;
      --jobs)     JOBS="$2" ; shift 2 ;;
      --tag)      TAG="$2" ; shift 2 ;;
      --gpu)      GPU_OPTS="$2" ; shift 2 ;;
      *)          echo "build.sh: unknown argument '$1' (see the usage at the top of the script)" >&2 ; exit 2 ;;
   esac
done
case "$LEG" in
   cpu)     MODE=flume-cpu-gnu     ; EXE=adam_flume_cpu ; COMPILERS="gfortran" ;;
   fnl-nvf) MODE=flume-fnl-nvf     ; EXE=adam_flume_fnl ; COMPILERS="nvfortran" ;;
   fnl-amd) MODE=flume-fnl-omp-amd ; EXE=adam_flume_fnl ; COMPILERS="amdflang flang" ;;
   *)       echo "build.sh: the leg must be cpu, fnl-nvf or fnl-amd, got '$LEG'" >&2 ; exit 2 ;;
esac
if [[ $LEG != cpu && -z $VARSET ]]; then
   echo "build.sh: leg $LEG needs --varset (it sets \$HDF5_PREFIX and the GPU target)" >&2
   exit 2
fi

if [[ -n $GPU_OPTS && $LEG != fnl-nvf ]]; then
   echo "build.sh: --gpu sets nvfortran -gpu= options, for the fnl-nvf leg only" >&2
   exit 2
fi
if [[ -n $GPU_OPTS && -z $TAG ]]; then
   echo "build.sh: --gpu needs --tag: a variant must not replace the executable of the $LEG leg" >&2
   exit 2
fi
if [[ -n $GPU_OPTS && $CLEAN -eq 0 ]]; then
   echo "build.sh: --gpu needs a clean build (the objects keep the flags they were compiled with)" >&2
   exit 2
fi
if [[ -n $JOBS && ! $JOBS =~ ^[1-9][0-9]*$ ]]; then
   echo "build.sh: --jobs takes a positive integer, got '$JOBS'" >&2
   exit 2
fi
if [[ -n $TAG && ! $TAG =~ ^[A-Za-z0-9_.-]+$ ]]; then
   echo "build.sh: --tag takes letters, digits, '.', '_' and '-' only, got '$TAG'" >&2
   exit 2
fi
NAME="$LEG${TAG:+-$TAG}"
OUT="$SPEEDUP_DIR/builds/$NAME"
mkdir -p "$OUT" "$SPEEDUP_DIR/exe/$NAME"
LOG="$OUT/build.txt"
fobis_opt=()
[[ -n $VARSET ]] && fobis_opt+=(--varset "$VARSET")
[[ -n $JOBS ]] && fobis_opt+=(-j "$JOBS")
[[ -n $GPU_OPTS ]] && fobis_opt+=(--cflags "-gpu=$GPU_OPTS" --lflags "-gpu=$GPU_OPTS")

cd "$REPO_ROOT"
start=$(date +%s)
status=0
{
   echo ">> $(date -Is) $(hostname): fobis build --mode $MODE ${fobis_opt[*]} (clean: $CLEAN)"
   if [[ $CLEAN -eq 1 ]]; then fobis clean --mode "$MODE"; fi
   fobis build --mode "$MODE" "${fobis_opt[@]}"
} > "$LOG" 2>&1 || status=$?
seconds=$(( $(date +%s) - start ))
if [[ $status -ne 0 || ! -x "exe/$EXE" ]]; then
   echo "build.sh: $MODE failed (status $status), see $LOG" >&2
   exit 1
fi
cp "exe/$EXE" "$SPEEDUP_DIR/exe/$NAME/$EXE"

first_line() { "$@" 2>&1 | grep -m1 -v '^\s*$' || true; }
{
   echo "| Item | Value |"
   echo "|---|---|"
   echo "| Built | $(date -Is) on \`$(hostname)\`, ${seconds} s |"
   cmd="fobis build --mode $MODE${VARSET:+ --varset $VARSET}"
   [[ -z $GPU_OPTS ]] || cmd+=" --cflags -gpu=$GPU_OPTS --lflags -gpu=$GPU_OPTS"
   echo "| Mode | \`$cmd\` (clean: $CLEAN) |"
   if [[ -n $TAG ]]; then echo "| Variant | \`$TAG\` (not a matrix leg) |"; fi
   # the sources of the executable: the scripts and inputs under src/tests do not enter it
   dirty="" ; git diff --quiet HEAD -- src ':!src/tests' || dirty=", **src/ modified**"
   echo "| Commit | \`$(git rev-parse HEAD)\`$dirty |"
   echo "| FoBiS | $(first_line fobis --version) |"
   for c in $COMPILERS mpif90; do
      if command -v "$c" > /dev/null; then echo "| \`$c\` | $(first_line "$c" --version) |"; fi
   done
   if command -v mpirun > /dev/null; then echo "| \`mpirun\` | $(first_line mpirun --version) |"; fi
   cflags="$(grep -a -m1 'Compiling flags:' "$LOG" | sed -E 's/.*Compiling flags: "(.*)"/\1/' || true)"
   lflags="$(grep -a -m1 'Linking flags:' "$LOG" | sed -E 's/.*Linking flags: "(.*)"/\1/' || true)"
   echo "| Compile flags | \`${cflags:-n/a}\` |"
   echo "| Link flags | \`${lflags:-n/a}\` |"
   echo "| Executable | \`exe/$NAME/$EXE\`, sha256 \`$(sha256sum "$SPEEDUP_DIR/exe/$NAME/$EXE" | cut -c1-16)\` |"
} > "$OUT/toolchain.md"

# a variant is not a matrix leg: the case READMEs embed the toolchain of the legs only
if [[ -z $TAG ]]; then python3 "$SPEEDUP_DIR/make_inputs.py" readmes; fi
echo "build.sh: $MODE built in ${seconds} s -> $SPEEDUP_DIR/exe/$NAME/$EXE; toolchain in $OUT/toolchain.md"
