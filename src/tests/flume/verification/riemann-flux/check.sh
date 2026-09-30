#!/usr/bin/env bash
# FLUME verification RV-2, RV-3: the weno-riemann scheme on the V2 vortex and the V1 Sod problems (issue #47).
#
# Why: `scheme_space = weno-riemann` replaces the WENO flux splitting by the WENO interpolation of the face states, a
# Riemann flux and a high-order correction (Chen, Toth & Gombosi 2016). Its oracles are those of the splitting scheme,
# on the same inputs with the [numerics] block rewritten:
#   RV-2 vortex  the isentropic vortex ladder 64/128/256 (vortex/): L1 order of the finest pair >= ORDER_MIN;
#   RV-3 sod     Sod along x, y, z (sod/): L1(rho) against the exact solution <= L1_MAX, y and z BITWISE equal to x
#                after the permutation, every transverse copy bitwise identical;
#   RV-3 lax     the Lax problem on the Sod inputs (left 0.445, 0.698, 3.528, right 0.5, 0, 0.571, t = 0.13 on [0, 1]),
#                same checks, L1 bound LAX_L1_MAX;
#   RV-3 shu-osher  the Shu-Osher problem on the Sod inputs (400 cells on [-5, 5], t = 1.8): L1(rho) against a 16x
#                finer split-scheme reference (shu_osher_oracle.py) <= SO_L1_MAX, positivity, x/y/z bitwise;
#   RV-4 conservation  the V3 periodic AMR box (conservation/): the five integrals constant within 1e-13 with reflux,
#                drifting by at least 1e-10 without (negative control);
#   V6 cylinder  the Mach 2 shock over the immersed cylinder with solid AMR (shock-cylinder/): refined surface blocks,
#                mirror symmetry within 1e-10, positive density and pressure.
#
# Usage: ./check.sh [--np N] [--leg vortex|sod|lax|shu-osher|conservation|cylinder|all] [--solver S] [--correction C]
#                   [--sensor S] [--recon R]
#                   [--order-min X] [--l1-max X] [--lax-l1-max X] [--so-l1-max X]
# Defaults: all legs, hllc, 6th, weno, characteristic; ORDER_MIN 4.8 (#47 RV-2); the RV-3 bounds are 1.05 x the split
# scheme's L1(rho) (#47 RV-3, owner decision 2026-09-30): L1_MAX (Sod) 3.41e-03, LAX_L1_MAX 1.12e-02, SO_L1_MAX 0.274;
# LLF exceeds the Lax and Shu-Osher bounds (see README) and needs explicit ones.
# `--solver split` keeps the inputs' flux-splitting scheme (`scheme_space = weno`): the baseline of the same legs.
#
# FLUME_EXE overrides the executable under test, e.g. FLUME_EXE=$REPO/exe/adam_flume_fnl ./check.sh
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH and, on WSL, the UCX knobs of issue #12).
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
VER_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${FLUME_EXE:-$REPO_ROOT/exe/adam_flume_cpu}"
NP=2
LEG=all
SOLVER=hllc
CORRECTION=6th
SENSOR=weno
RECON=characteristic
ORDER_MIN="4.8"
L1_MAX="3.41e-03"
LAX_L1_MAX="1.12e-02"
SO_L1_MAX="0.274"

while [[ $# -gt 0 ]]; do
   case "$1" in
      --np)         NP="$2" ; shift 2 ;;
      --leg)        LEG="$2" ; shift 2 ;;
      --solver)     SOLVER="$2" ; shift 2 ;;
      --correction) CORRECTION="$2" ; shift 2 ;;
      --sensor)     SENSOR="$2" ; shift 2 ;;
      --recon)      RECON="$2" ; shift 2 ;;
      --order-min)  ORDER_MIN="$2" ; shift 2 ;;
      --l1-max)     L1_MAX="$2" ; shift 2 ;;
      --lax-l1-max) LAX_L1_MAX="$2" ; shift 2 ;;
      --so-l1-max)  SO_L1_MAX="$2" ; shift 2 ;;
      *) echo "check.sh: unknown argument '$1'" >&2 ; exit 2 ;;
   esac
done
if [[ ! -x "$EXE" ]]; then
   echo "check.sh: executable '$EXE' not found (build it or set FLUME_EXE)" >&2
   exit 2
fi
VENV_PY="$REPO_ROOT/exe/.regression-venv/bin/python"
TAG="$(basename "$EXE")-np$NP-$SOLVER-$CORRECTION-$SENSOR-$RECON"

set_numerics() { # rewrite the [numerics] block of an input for weno-riemann (split: keep it)
   [[ $SOLVER == split ]] && return 0
   "$VENV_PY" - "$1" "$SOLVER" "$CORRECTION" "$SENSOR" "$RECON" <<'EOF'
import re, sys
path, solver, corr, sensor, recon = sys.argv[1:]
text = open(path).read()
new = (f"scheme_space             = weno-riemann\nreconstruction_variables = {recon}\nriemann_solver           = {solver}\n"
       f"flux_correction          = {corr}\nflux_correction_sensor   = {sensor}\n")
text, n = re.subn(r"^scheme_space\s*=.*\n^reconstruction_variables\s*=.*\n", new, text, flags=re.M)
if n != 1:
    sys.exit(f"{path}: [numerics] block not found")
open(path, "w").write(text)
EOF
}

set_lax() { # turn a Sod input into the Lax problem: states of the two regions, normal velocity by the input's axis
   "$VENV_PY" - "$1" <<'EOF'
import re, sys
path = sys.argv[1]
axis = re.search(r"sod-([xyz])\.ini$", path).group(1)
un = {"x": "u", "y": "v", "z": "w"}[axis]
text = open(path).read()
for region, state in (("1", {"r": "0.445", un: "0.698", "p": "3.528"}), ("2", {"r": "0.5", un: "0.0", "p": "0.571"})):
    head = f"[initial_conditions_region_{region}]"
    start = text.index(head)
    end = text.find("\n[", start + len(head))
    block = text[start:end]
    for key, value in state.items():
        block, n = re.subn(rf"^{key}(\s*)=.*$", rf"{key}\g<1>= {value}", block, flags=re.M)
        if n != 1:
            sys.exit(f"{path}: key {key} not found in {head}")
    text = text[:start] + block + text[end:]
text, n = re.subn(r"^time_max(\s*)=.*$", r"time_max\g<1>= 0.13", text, flags=re.M)
if n != 1:
    sys.exit(f"{path}: time_max not found")
open(path, "w").write(text)
EOF
}

set_shu_osher() { # turn a Sod input into the Shu-Osher problem along the input's axis
   "$VENV_PY" - "$1" <<'EOF'
import re, sys
path = sys.argv[1]
axis = re.search(r"sod-([xyz])\.ini$", path).group(1)
un = {"x": "u", "y": "v", "z": "w"}[axis]
nc = {"x": "ni", "y": "nj", "z": "nk"}[axis]
text = open(path).read()


def section(text, head, subs):
    start = text.index(head)
    end = text.find("\n[", start + len(head))
    end = len(text) if end < 0 else end
    block = text[start:end]
    for key, value in subs.items():
        block, n = re.subn(rf"^{key}(\s*)=.*$", rf"{key}\g<1>= {value}", block, flags=re.M)
        if n != 1:
            sys.exit(f"{path}: key {key} not found in {head}")
    return text[:start] + block + text[end:]


text = section(text, "[grid]", {nc: "100", f"emin_{axis}": "-5.0", f"emax_{axis}": "5.0"})
text = section(text, "[initial_conditions]", {"type": "shu-osher"})
text = text.replace("[initial_conditions]\n", "[initial_conditions]\n"
                    f"axis           = {axis}\ninterface      = -4.0\nrho_amplitude  = 0.2\nrho_wavenumber = 5.0\n", 1)
text = section(text, "[initial_conditions_region_1]", {"r": "3.857143", un: "2.629369", "p": "10.33333"})
text = section(text, "[initial_conditions_region_2]", {"r": "1.0", un: "0.0", "p": "1.0"})
text = section(text, "[time]", {"time_max": "1.8"})
open(path, "w").write(text)
EOF
}

set_reflux_off() { # the negative control of the conservation leg
   sed -i "s/^reflux                   = .true./reflux                   = .false./" "$1"
   grep -q "^reflux                   = .false." "$1"
}

run() { # case-dir input work [input-hook]
   rm -rf "$3"
   mkdir -p "$3"
   cp "$1/$2" "$3/"
   set_numerics "$3/$2"
   if [[ -n ${4:-} ]]; then "$4" "$3/$2"; fi
   echo ">> $2: mpirun -np $NP $(basename "$EXE") ($SOLVER, $CORRECTION, $SENSOR, $RECON)"
   local start
   start=$(date +%s)
   if ! (cd "$3" && mpirun -np "$NP" "$EXE" "$2" > log.txt 2>&1); then
      echo "check.sh: $2 run failed, see $3/log.txt" >&2
      exit 1
   fi
   echo "   done in $(( $(date +%s) - start )) s"
}

if [[ $LEG == vortex || $LEG == all ]]; then
   WORK=()
   for n in 064 128 256; do
      run "$VER_DIR/vortex" "vortex-n$n.ini" "$CASE_DIR/work-$TAG-vortex-n$n"
      WORK+=("$CASE_DIR/work-$TAG-vortex-n$n")
   done
   "$VENV_PY" "$VER_DIR/vortex/vortex_oracle.py" "${WORK[@]}" --order-min "$ORDER_MIN"
   echo "RV-2 PASSED ($TAG)"
fi
if [[ $LEG == sod || $LEG == all ]]; then
   WORK=()
   for d in x y z; do
      run "$VER_DIR/sod" "sod-$d.ini" "$CASE_DIR/work-$TAG-sod-$d"
      WORK+=("$CASE_DIR/work-$TAG-sod-$d")
   done
   "$VENV_PY" "$VER_DIR/sod/sod_oracle.py" "$VER_DIR/sod/sod-x.ini" "${WORK[@]}" --l1-max "$L1_MAX"
   echo "RV-3 sod PASSED ($TAG)"
fi
if [[ $LEG == lax || $LEG == all ]]; then
   WORK=()
   for d in x y z; do
      run "$VER_DIR/sod" "sod-$d.ini" "$CASE_DIR/work-$TAG-lax-$d" set_lax
      WORK+=("$CASE_DIR/work-$TAG-lax-$d")
   done
   "$VENV_PY" "$VER_DIR/sod/sod_oracle.py" "${WORK[0]}/sod-x.ini" "${WORK[@]}" ${LAX_L1_MAX:+--l1-max "$LAX_L1_MAX"}
   echo "RV-3 lax PASSED ($TAG)"
fi
if [[ $LEG == shu-osher || $LEG == all ]]; then
   WORK=()
   for d in x y z; do
      run "$VER_DIR/sod" "sod-$d.ini" "$CASE_DIR/work-$TAG-shu-osher-$d" set_shu_osher
      WORK+=("$CASE_DIR/work-$TAG-shu-osher-$d")
   done
   "$VENV_PY" "$CASE_DIR/shu_osher_oracle.py" "${WORK[0]}/sod-x.ini" "${WORK[@]}" ${SO_L1_MAX:+--l1-max "$SO_L1_MAX"}
   echo "RV-3 shu-osher PASSED ($TAG)"
fi
if [[ $LEG == conservation || $LEG == all ]]; then
   run "$VER_DIR/conservation" amr-periodic.ini "$CASE_DIR/work-$TAG-conservation-reflux-true"
   run "$VER_DIR/conservation" amr-periodic.ini "$CASE_DIR/work-$TAG-conservation-reflux-false" set_reflux_off
   for reflux in true false; do
      log="$CASE_DIR/work-$TAG-conservation-reflux-$reflux/log.txt"
      if ! grep -q "registered intra-realm AMR seam faces: +6" "$log"; then
         echo "check.sh: the case must register 6 coarse-fine faces (reflux .$reflux.)" >&2
         exit 1
      fi
   done
   "$VENV_PY" "$VER_DIR/conservation/conservation_oracle.py"                             \
      --conserved "$CASE_DIR/work-$TAG-conservation-reflux-true" --max-drift 1.0e-13    \
      --leaky "$CASE_DIR/work-$TAG-conservation-reflux-false" --min-drift 1.0e-10
   echo "RV-4 PASSED ($TAG)"
fi
if [[ $LEG == cylinder || $LEG == all ]]; then
   run "$VER_DIR/shock-cylinder" shock-cylinder.ini "$CASE_DIR/work-$TAG-cylinder"
   "$VENV_PY" "$VER_DIR/shock-cylinder/shock_cylinder_oracle.py" "$CASE_DIR/work-$TAG-cylinder" --mirror-tol 1.0e-10
   echo "V6 cylinder PASSED ($TAG)"
fi
