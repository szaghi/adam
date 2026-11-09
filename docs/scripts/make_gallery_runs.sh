#!/usr/bin/env bash
# Run the four cases of the landing-page gallery and render them (docs/scripts/make_gallery.py).
#
# Every case is a committed verification input or generator, with only the checkpoint period changed so that each run
# leaves about 100 frames:
#   orszag-tang  FLUME MV-12 (src/tests/flume/verification/mhd/orszag-tang), 128^2, GLM, t = 0.5, a checkpoint
#                every 13 steps (1281 steps, about 1.5 min on 4 ranks);
#   rt-hydro     FLUME G3 rt2d-nl (src/tests/flume/verification/gravity/make_atmosphere.py), 128^2, t = 5;
#   rt-mhd       FLUME G3 rt2d-mhd, the same with EGLM and B = (0.1, 0, 0), half the critical field;
#   step         FLUME EV-step full (src/tests/flume/verification/step/make_step.py), the Woodward-Colella Mach 3
#                step on the three-realm forest, quadtrees, 1/80 base with the static refined boxes, t = 4;
#   pulse-pec, pulse-sm, pulse-pml
#                PRISM tutorials 2 and 3 (docs/tutorials/prism/faces.md, pml.md): the Gaussian pulse of
#                src/tests/prism/regression/periodic-pulse-blocks started at the centre against PEC, Silver-Muller and
#                Silver-Muller + classic PML x faces, all three run to the tutorials' final time 6.67e-10 s at CFL
#                0.18, below the PML cap of the time step, so the three runs share one time step and one set of
#                checkpoint times.
# One run at a time, NP ranks each.
#
# Usage: make_gallery_runs.sh <runs-dir> [--np N] [--cases LIST] [--no-render]
#   LIST: any of orszag-tang rt step pulse (default all); the renderer then writes docs/public/gallery.
#
# FLUME_EXE and PRISM_EXE override the executables (default exe/adam_flume_cpu, exe/adam_prism_cpu); the caller owns
# the matching MPI environment (mpirun on PATH).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
VERIF="$REPO_ROOT/src/tests/flume/verification"
FLUME_EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
PRISM_EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
NP=4
CASES="orszag-tang rt step pulse"
RENDER=1

[[ $# -ge 1 ]] || { echo "usage: make_gallery_runs.sh <runs-dir> [--np N] [--cases LIST] [--no-render]" >&2 ; exit 2 ; }
RUNS="$(mkdir -p "$1" && cd "$1" && pwd)"
shift
while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)        NP="$2" ; shift 2 ;;
      --cases)     CASES="$2" ; shift 2 ;;
      --no-render) RENDER=0 ; shift ;;
      *) echo "make_gallery_runs.sh: unknown argument '$1' (accepted: --np N, --cases LIST, --no-render)" >&2 ; exit 2 ;;
   esac
done

fresh() { # fresh <name>: an empty work directory, echoed
   local w="$RUNS/$1"
   rm -rf "${w:?}" ; mkdir -p "$w"
   echo "$w"
}

run() { # run <work> <exe> <input>: run one case, stop on failure
   echo "-- $(basename "$1"): $(basename "$2") $3, np $NP"
   local t0=$SECONDS
   if ! (cd "$1" && mpirun -np "$NP" "$2" "$3" < /dev/null > log.txt 2>&1); then
      echo "make_gallery_runs.sh: $(basename "$1") failed, see $1/log.txt" >&2
      exit 1
   fi
   echo "   done in $((SECONDS - t0)) s"
}

for case in $CASES; do
   case $case in
      orszag-tang)
         w=$(fresh orszag-tang)
         "$VENV_PY" "$VERIF/mhd/orszag-tang/make_orszag_tang.py" "$VERIF/vortex/vortex-n064.ini" \
            "$w/orszag-tang.ini" --cells 128 --it-save 13
         run "$w" "$FLUME_EXE" orszag-tang.ini ;;
      rt)
         w=$(fresh rt-hydro)
         python3 -I "$VERIF/gravity/make_atmosphere.py" "$w/input.ini" rt2d-nl --cells 32 --it-save 20
         run "$w" "$FLUME_EXE" input.ini
         w=$(fresh rt-mhd)
         python3 -I "$VERIF/gravity/make_atmosphere.py" "$w/input.ini" rt2d-mhd --cells 32 --it-save 40
         run "$w" "$FLUME_EXE" input.ini ;;
      step)
         w=$(fresh step)
         "$VENV_PY" "$VERIF/step/make_step.py" "$w" --tree quad --cells 80 --time-max 4 --refine --it-save 75
         run "$w" "$FLUME_EXE" step.ini ;;
      pulse)
         base="$REPO_ROOT/src/tests/prism/regression/periodic-pulse-blocks/input.ini"
         for face in pec sm pml; do
            w=$(fresh "pulse-$face")
            kind=$([[ $face == pec ]] && echo PEC || echo Silver_Muller)
            sed -e 's/^r0 *=.*/r0 = 0.0/' \
                -e "/^\[bc_x_m[ai][nx]\]/{n;s/^type *=.*/type = $kind/}" \
                -e 's/^output_basename *=.*/output_basename = pulse/' \
                -e 's/^it_save *=.*/it_save = 6/' \
                -e 's/^time_max *=.*/time_max = 6.6712814328639916E-10/' \
                -e 's/^it_max *=.*/it_max = -1/' \
                -e 's/^CFL *=.*/CFL = 0.18/' "$base" > "$w/input.ini"
            if [[ $face == pml ]]; then
               cat >> "$w/input.ini" <<'EOF'

[PML]
PML_type       = classic
width          = 0.04
x_minus_layer  = .true.
x_plus_layer   = .true.
y_minus_layer  = .false.
y_plus_layer   = .false.
z_minus_layer  = .false.
z_plus_layer   = .false.
gamma_max      = 7.8E10
gamma_exponent = 2.0
EOF
            fi
            run "$w" "$PRISM_EXE" input.ini
         done ;;
      *) echo "make_gallery_runs.sh: unknown case '$case' (accepted: orszag-tang rt step pulse)" >&2 ; exit 2 ;;
   esac
done

if [[ $RENDER -eq 1 ]]; then
   "$VENV_PY" "$REPO_ROOT/docs/scripts/make_gallery.py" --runs "$RUNS" --out "$REPO_ROOT/docs/public/gallery"
fi
