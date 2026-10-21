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
# Usage: ./build.sh <leg> [--varset "V ..."] [--no-clean]
set -euo pipefail

SPEEDUP_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SPEEDUP_DIR/../../../.." && pwd)"
LEG="${1:-}"
[[ $# -gt 0 ]] && shift
VARSET=""
CLEAN=1
while [[ $# -gt 0 ]]; do
   case "$1" in
      --varset)   VARSET="$2" ; shift 2 ;;
      --no-clean) CLEAN=0 ; shift ;;
      *)          echo "build.sh: unknown argument '$1' (accepted: --varset V, --no-clean)" >&2 ; exit 2 ;;
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

OUT="$SPEEDUP_DIR/builds/$LEG"
mkdir -p "$OUT" "$SPEEDUP_DIR/exe/$LEG"
LOG="$OUT/build.txt"
varset_opt=()
[[ -n $VARSET ]] && varset_opt=(--varset "$VARSET")

cd "$REPO_ROOT"
start=$(date +%s)
status=0
{
   echo ">> $(date -Is) $(hostname): fobis build --mode $MODE ${VARSET:+--varset \"$VARSET\"} (clean: $CLEAN)"
   if [[ $CLEAN -eq 1 ]]; then fobis clean --mode "$MODE"; fi
   fobis build --mode "$MODE" "${varset_opt[@]}"
} > "$LOG" 2>&1 || status=$?
seconds=$(( $(date +%s) - start ))
if [[ $status -ne 0 || ! -x "exe/$EXE" ]]; then
   echo "build.sh: $MODE failed (status $status), see $LOG" >&2
   exit 1
fi
cp "exe/$EXE" "$SPEEDUP_DIR/exe/$LEG/$EXE"

first_line() { "$@" 2>&1 | grep -m1 -v '^\s*$' || true; }
{
   echo "| Item | Value |"
   echo "|---|---|"
   echo "| Built | $(date -Is) on \`$(hostname)\`, ${seconds} s |"
   echo "| Mode | \`fobis build --mode $MODE${VARSET:+ --varset $VARSET}\` (clean: $CLEAN) |"
   dirty="" ; git diff --quiet HEAD -- src || dirty=", **src/ modified**"
   echo "| Commit | \`$(git rev-parse HEAD)\`$dirty |"
   echo "| FoBiS | $(first_line fobis --version) |"
   for c in $COMPILERS mpif90; do
      if command -v "$c" > /dev/null; then echo "| \`$c\` | $(first_line "$c" --version) |"; fi
   done
   if command -v mpirun > /dev/null; then echo "| \`mpirun\` | $(first_line mpirun --version) |"; fi
   flags="$(grep -m1 -E ' -c .*\.F90' "$LOG" | sed -E 's/ [^ ]+\.F90.*//' || true)"
   echo "| First compile line | \`${flags:-n/a}\` |"
   echo "| Executable | \`exe/$LEG/$EXE\`, sha256 \`$(sha256sum "$SPEEDUP_DIR/exe/$LEG/$EXE" | cut -c1-16)\` |"
} > "$OUT/toolchain.md"

python3 "$SPEEDUP_DIR/make_inputs.py" readmes
echo "build.sh: $MODE built in ${seconds} s -> $SPEEDUP_DIR/exe/$LEG/$EXE; toolchain in $OUT/toolchain.md"
