#!/usr/bin/env bash
# FLUME verification PV-1 to PV-4 (issue #47, M3-P5): the positivity limiter, `[numerics] positivity_limiter = cell`.
#
# Why: the limiter (D-9) blends every face flux with the first-order Lax-Friedrichs backbone so that each stage keeps
# the density and the pressure positive; it must make the hard problems run without floors, and leave the smooth ones
# untouched. np 2, the floors disabled (a non-positive MHD density or pressure stops the run). A leg with the limiter
# passes when the run reaches its final time and the limiter never meets an inadmissible backbone (the log line
# `positivity limiter: N faces limited, M inadmissible backbones` has M = 0 at every stage); the Euler legs, which have
# no floors to stop the run, also assert a positive density and pressure in the last checkpoint. Each problem runs on
# the splitting scheme and on weno-riemann (HLLD or HLLC, characteristic). Legs:
#   pv1     the Balsara-Spicer blast (beta 2.5e-4, 128^2, t = 0.01, EGLM, ../positivity-probe/make_blast.py); the
#           splitting scheme fails at step 3 without the limiter (M3-P4c);
#   pv1-3d  the same blast as a sphere on 64^3 (not in the default legs: about 1 h per scheme on the CPU);
#   pv2     the second blast of Wu & Shu (2018, doi:10.1137/18M1168042, Example 4.4): p 1e4 in the disc, B =
#           1000/sqrt(4 pi) along x, beta 2.51e-6, t = 0.001, 128^2, EGLM with c_h = 400; without the limiter the
#           splitting scheme fails at step 1 and HLLD at step 41;
#   pv3     Euler near vacuum (make_vacuum.py): the LeBlanc shock tube, the double rarefaction and the planar Sedov
#           blast. The splitting scheme needs no limiting on any of them; HLLC needs none on LeBlanc and fails without
#           the limiter on the double rarefaction (step 79) and on Sedov (step 69);
#   pv4     inactive on smooth problems: with the limiter on, the isentropic vortex (Euler, N = 64, 128) and the 1-D
#           fast wave (EGLM, N = 16, 32) are BITWISE equal to the runs without it (interior cells; a face whose two
#           cells keep Lambda = 1 is not touched, and the EGLM sources with Lambda = 1 take the unlimited arithmetic).
# Sedov on HLLC is the case of the limiter's relative floor (POSITIVITY_LIMITER_KAPPA): with the absolute floor alone
# one stage drains the centre cell to rho = 1e-13 with its energy kept, the sound speed grows by 1e6 and the next stage
# meets inadmissible backbones (NaN at step 10, at any CFL).
# One run at a time; checkpoints deleted after use. On the CPU: pv1 4 min, pv2 4 min, pv3 55 min, pv4 9 min.
#
# Usage: ./check.sh [--np N] [--legs pv1,pv1-3d,pv2,pv3,pv4]     (default legs: pv1,pv2,pv3,pv4)
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh (M3-P5b)
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VERIF_DIR="$(cd "$CASE_DIR/../.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEGS="pv1,pv2,pv3,pv4"
source "$CASE_DIR/../numerics.sh"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)   NP="$2" ; shift 2 ;;
      --legs) LEGS="$2" ; shift 2 ;;
      *)      echo "check.sh: unknown argument '$1' (accepted: --np N, --legs LIST)" >&2 ; exit 2 ;;
   esac
done
for leg in ${LEGS//,/ }; do
   if [[ ! $leg =~ ^(pv1|pv1-3d|pv2|pv3|pv4)$ ]]; then
      echo "check.sh: unknown leg '$leg' (accepted: pv1, pv1-3d, pv2, pv3, pv4)" >&2 ; exit 2
   fi
done
has_leg() { [[ ",$LEGS," == *",$1,"* ]]; }
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP"
FAILED=0

limiter_on() { # limiter_on <ini>: add positivity_limiter = cell to the [numerics] block (after reflux)
   sed -i -E 's/^(reflux\s*=.*)$/\1\npositivity_limiter       = cell/' "$1"
   grep -q '^positivity_limiter *= *cell' "$1"
}
run() { # run <work> <ini>: run one case, fail the check if the run fails
   if ! (cd "$1" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: $(basename "$1") failed (a non-positive state stops the run), see $1/log.txt" >&2
      FAILED=1
      return 1
   fi
}
same() { # same <work-a> <work-b>: the interiors of the last checkpoints of the two runs bitwise equal (FNL ghost
         # cells are not reproducible run to run, M3-P3b)
   "$VENV_PY" - "$1" "$2" <<'EOF'
import glob
import sys

import h5py
import numpy as np

a, b = sys.argv[1:]
fa = sorted(f for f in glob.glob(a + "/*-proc*.h5") if "restart" not in f)
last = max(int(f.split("-")[-2]) for f in fa)
diff = 0.0
for f in (f for f in fa if int(f.split("-")[-2]) == last):
    with h5py.File(f) as x, h5py.File(f.replace(a, b)) as y:
        for k in x:
            if x[k].ndim == 3:
                diff = max(diff, float(np.max(np.abs(x[k][()][3:-3, 3:-3, 3:-3] - y[k][()][3:-3, 3:-3, 3:-3]))))
ok = diff == 0.0
print(f"   {b.split('/')[-1]} vs {a.split('/')[-1]}: step {last}, max |difference| {diff:.3e}  "
      f"{'PASS (bitwise)' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
}

verdict() { # verdict <work> <basename> <t_final> <label>: final time reached, no inadmissible backbone
   local last limited badrows
   last=$(tail -1 "$1/$2-conservation_history.dat" | awk '{print $2}')
   limited=$(grep -a 'positivity limiter:' "$1/log.txt" | sed -E 's/.*: \+([0-9]+) faces.*/\1/' | sort -n | tail -1 \
             || true) # no line when the limiter never acts
   badrows=$(grep -a 'positivity limiter:' "$1/log.txt" | grep -vc '+0 inadmissible' || true)
   if "$VENV_PY" -c "import sys; sys.exit(0 if abs($last - $3) < 1e-12 else 1)" && [[ $badrows -eq 0 ]]; then
      echo "   $4: t = $last, max ${limited:-0} faces limited per stage, 0 inadmissible backbones  PASS"
   else
      echo "   $4: t = $last, $badrows stages with inadmissible backbones  FAIL" ; FAILED=1
   fi
}
blast() { # blast <tag> <t_final> <make_blast.py arguments>: the limited blast on the two schemes
   local tag=$1 tend=$2 spec w
   shift 2
   for spec in split hlld:characteristic; do
      w="$CASE_DIR/work-$TAG-$tag-${spec//:/-}"
      rm -rf "$w" ; mkdir -p "$w"
      "$VENV_PY" "$VERIF_DIR/mhd/positivity-probe/make_blast.py" "$w/blast.ini" --eglm "$@" > /dev/null
      NUMERICS="" ; NUMERICS_SOLVER=""
      if [[ $spec != split ]]; then NUMERICS=$spec; numerics_check; numerics_apply "$w/blast.ini"; fi
      limiter_on "$w/blast.ini"
      if run "$w" blast.ini; then verdict "$w" blast "$tend" "$spec"; fi
      find "$w" -name '*.h5' -delete
   done
}
positive() { # positive <work> <gamma>: density and pressure of the last Euler checkpoint are positive (interior cells)
   "$VENV_PY" - "$1" "$2" <<'EOF'
import glob
import sys

import h5py
import numpy as np

w, gamma = sys.argv[1], float(sys.argv[2])
fs = sorted(f for f in glob.glob(w + "/*-proc*.h5") if "restart" not in f)
last = max(int(f.split("-")[-2]) for f in fs)
rho_min = p_min = np.inf
for f in (f for f in fs if int(f.split("-")[-2]) == last):
    with h5py.File(f) as x:
        blocks = {}
        for k in x:
            if x[k].ndim == 3:
                name, v = k.rsplit("-", 1)
                blocks.setdefault(name, {})[v] = x[k][()][3:-3, 3:-3, 3:-3]
        for q in blocks.values():
            p = (gamma - 1.0) * (q["rE"] - 0.5 * (q["ru"] ** 2 + q["rv"] ** 2 + q["rw"] ** 2) / q["r"])
            rho_min, p_min = min(rho_min, float(q["r"].min())), min(p_min, float(p.min()))
ok = rho_min > 0.0 and p_min > 0.0
print(f"      step {last}: min density {rho_min:.3e}, min pressure {p_min:.3e}  {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
EOF
}

if has_leg pv1; then
   echo ">> PV-1 Balsara-Spicer blast, 128^2, EGLM, limiter on ($(basename "$EXE"), np $NP)"
   blast blast 0.01 --cells 128
fi
if has_leg pv1-3d; then
   echo ">> PV-1 Balsara-Spicer blast as a sphere, 64^3, EGLM, limiter on"
   blast blast3d 0.01 --cells 64 --3d
fi
if has_leg pv2; then
   echo ">> PV-2 Wu-Shu blast, beta 2.51e-6, 128^2, EGLM, limiter on"
   blast wushu 0.001 --cells 128 --b-axis x --p-in 1.0e4 --b0 282.0947917738782 --time-max 0.001 --glm-ch 400.0
fi
if has_leg pv3; then
   echo ">> PV-3 Euler near vacuum, limiter on"
   for item in leblanc:6.0:1.6666666666666667 double-rarefaction:0.6:1.4 sedov:0.001:1.4; do
      IFS=: read -r prob tend gamma <<< "$item"
      for solver in split hllc; do
         w="$CASE_DIR/work-$TAG-$prob-$solver"
         rm -rf "$w" ; mkdir -p "$w"
         "$VENV_PY" "$CASE_DIR/make_vacuum.py" "$VERIF_DIR/sod/sod-x.ini" "$w/$prob.ini" --problem "$prob"
         NUMERICS="" ; NUMERICS_SOLVER=""
         if [[ $solver == hllc ]]; then # the Euler solver, which numerics_check (MHD) does not know
            NUMERICS=hllc NUMERICS_SOLVER=hllc NUMERICS_RECON=characteristic NUMERICS_CORRECTION=6th NUMERICS_SENSOR=weno
            numerics_apply "$w/$prob.ini"
         fi
         limiter_on "$w/$prob.ini"
         if run "$w" "$prob.ini"; then
            verdict "$w" "$prob" "$tend" "$prob $solver"
            positive "$w" "$gamma" || FAILED=1
         fi
         find "$w" -name '*.h5' -delete
      done
   done
fi

if has_leg pv4; then
   echo ">> PV-4 the limiter is inactive on smooth problems: bitwise equal runs"
   for n in 064 128; do
      for v in off on; do
         w="$CASE_DIR/work-$TAG-vortex-n$n-$v"
         rm -rf "$w" ; mkdir -p "$w"
         cp "$VERIF_DIR/vortex/vortex-n$n.ini" "$w/"
         [[ $v == on ]] && limiter_on "$w/vortex-n$n.ini"
         run "$w" "vortex-n$n.ini" || true
      done
      same "$CASE_DIR/work-$TAG-vortex-n$n-off" "$CASE_DIR/work-$TAG-vortex-n$n-on" || FAILED=1
   done
   for n in 16 32; do
      for v in off on; do
         w="$CASE_DIR/work-$TAG-wave-n$n-$v"
         rm -rf "$w" ; mkdir -p "$w"
         "$VENV_PY" "$VERIF_DIR/mhd/linear-wave/make_linear_wave.py" "$VERIF_DIR/sod/sod-x.ini" "$w/linear-wave.ini" \
            --wave fast --geometry 1d --cells "$n" --amplitude 1.0e-7 --divergence-control eglm
         [[ $v == on ]] && limiter_on "$w/linear-wave.ini"
         run "$w" linear-wave.ini || true
      done
      same "$CASE_DIR/work-$TAG-wave-n$n-off" "$CASE_DIR/work-$TAG-wave-n$n-on" || FAILED=1
   done
fi
find "$CASE_DIR" -path "$CASE_DIR/work-$TAG-*" -name '*.h5' -delete

if [[ $FAILED -eq 0 ]]; then
   echo "legs $LEGS PASSED ($TAG)"
else
   echo "legs $LEGS FAILED ($TAG)"
   exit 1
fi
