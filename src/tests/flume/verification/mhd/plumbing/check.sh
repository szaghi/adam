#!/usr/bin/env bash
# FLUME MHD plumbing verification (issue #41, M2-P1): the MHD model end to end with a zero residual.
#
# Why: M2-P1 wired the MHD model (input predicate, state width and names, MHD conversions and fast speed, the per-model
# kernel instances, IC/BC keys and wall rule, outputs, restart). For both variants, divergence_control = none (nv = 8)
# and glm (nv = 9, psi):
#   1. dt: one step of the unperturbed uniform state; the step must be CFL / sum_d (|u_d| + c_{f,d}) / dx_d with the
#      fast speed computed independently by the oracle (exercises the device aux + fast-speed kernels on FNL), with GLM
#      bounded by CFL / (c_h sum_d 1 / dx_d) (here the c_h bound is the active one, M2-P4);
#   2. auxiliary fields: the saved MHD auxiliaries equal the values recomputed from the conservative fields;
#   3. restart: 5 steps + restart + 5 steps of a seeded (s = 0.05) state equals the continuous 10-step run bitwise,
#      histories byte-identical, and every model variable (bx, by, bz [, psi]) is saved.
# Since M2-P3 the MHD face fluxes are live, so the P1 zero-residual legs (constant integrals, 10 steps = 1 step) are
# retired: the exact steady state is checked by MV-2 (verification/mhd/zero-field, uniform state, bitwise).
# Plus two refused configurations (variant none only): an unknown divergence_control and isentropic-vortex + MHD
# must stop with their error message.
#
# Usage: ./check.sh [--build] [--np N]
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12);
# --build always builds the CPU default, never the override.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
BUILD=0

while [[ $# -gt 0 ]]; do
   case "$1" in
      --build) BUILD=1 ; shift ;;
      --np)    NP="$2" ; shift 2 ;;
      *)       echo "check.sh: unknown argument '$1' (accepted: --build, --np N)" >&2 ; exit 2 ;;
   esac
done

if [[ $BUILD -eq 1 ]]; then
   (cd "$REPO_ROOT" && fobis build --mode flume-cpu-gnu)
fi
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi

VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
if ! "$VENV_PY" -c 'import h5py, numpy' 2>/dev/null; then
   echo ">> creating the oracle venv at $REPO_ROOT/exe/.regression-venv"
   python3 -m venv "$REPO_ROOT/exe/.regression-venv"
   "$VENV_PY" -m pip install --quiet --upgrade pip
   "$VENV_PY" -m pip install --quiet h5py numpy
fi
ORACLE="$CASE_DIR/mhd_plumbing_oracle.py"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

run() { # run <work-dir> <ini> <log>
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > "$3" 2>&1); then
      echo "check.sh: run failed, see $1/$3" >&2
      exit 1
   fi
}

oracle() { "$VENV_PY" "$ORACLE" "$@" || FAILED=1 ; }

setup() { # setup <work-dir> <source-ini> [sed expressions...]: fresh work dir with the edited input as input.ini
   local w="$1" src="$2" ; shift 2
   rm -rf "$w" ; mkdir -p "$w"
   sed -e "" "$@" "$src" > "$w/input.ini"
}

refused() { # refused <label> <expected message> <source-ini> [sed expressions...]
   local label="$1" msg="$2" src="$3" w="$CASE_DIR/work-$TAG-refused-$1" ; shift 3
   setup "$w" "$src" "$@"
   (cd "$w" && mpirun -np "$NP" "$EXE" input.ini > log.txt 2>&1) || true
   if grep -qF "$msg" "$w/log.txt"; then
      echo "   refused $label: stops with \"$msg\"  PASS"
   else
      echo "   refused $label: expected message \"$msg\" not found in $w/log.txt  FAIL"
      FAILED=1
   fi
}

for variant in none glm; do
   ini="$CASE_DIR/mhd-uniform-$variant.ini"
   names=(r ru rv rw rE bx by bz)
   [[ $variant == glm ]] && names+=(psi)
   base="$CASE_DIR/work-$TAG-$variant"
   echo ">> MHD plumbing, divergence_control = $variant ($(basename "$EXE"), np $NP)"
   # 1. dt of the unperturbed uniform state
   setup "$base-dt" "$ini" -e "s/^s              = .*/s              = 0.0/" -e "s/^it_max   = .*/it_max   = 1/"
   run "$base-dt" input.ini log.txt
   oracle --dt "$base-dt" "$base-dt/input.ini"
   # 2. auxiliary fields of the seeded 10-step run
   setup "$base-A" "$ini"
   run "$base-A" input.ini log.txt
   oracle --aux "$base-A" "$base-A/input.ini"
   # 3. restart round trip
   setup "$base-B" "$ini" -e "s/^it_max   = .*/it_max   = 5/" -e "s/^restart_save           = .*/restart_save           = 5/"
   run "$base-B" input.ini log-1.txt
   sed -i -e "s/^it_max   = .*/it_max   = 10/" -e "s/^restart_save           = .*/restart_save           = 0/" \
          -e "s/^restart                = .*/restart                = .true./" "$base-B/input.ini"
   run "$base-B" input.ini log-2.txt
   oracle --compare "$base-A" "$base-B" --names "${names[@]}"
   for hist in residuals conservation_history; do
      if cmp -s "$base-A/mhd-uniform-$hist.dat" "$base-B/mhd-uniform-$hist.dat"; then
         echo "   restart $hist history: identical  PASS"
      else
         echo "   restart $hist history: differs  FAIL" ; FAILED=1
      fi
   done
   find "$base"-* -name '*.h5' -delete
done

echo ">> refused configurations"
refused divergence-control 'unknown [mhd].(divergence_control) "bogus"' "$CASE_DIR/mhd-uniform-none.ini" \
        -e "s/^divergence_control = .*/divergence_control = bogus/"
refused vortex-mhd 'isentropic-vortex requires [physics].(physical_model) = euler' "$CASE_DIR/mhd-uniform-none.ini" \
        -e "s/^type           = uniform/type           = isentropic-vortex/"

if [[ $FAILED -eq 0 ]]; then
   echo "MHD plumbing verification PASSED ($TAG)"
else
   echo "MHD plumbing verification FAILED ($TAG)"
   exit 1
fi
