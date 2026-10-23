#!/usr/bin/env bash
# PRISM init refusals: inputs that must stop at initialization with their own message instead of running (manual;
# run.sh does not call it). Each leg rewrites the rmf input with sed and requires a non-zero exit AND its message in the
# log: a crash also exits non-zero, and only the message tells a refusal from a crash.
#
# Usage: PRISM_EXE=/path/to/adam_prism_<cpu|fnl> ./refusals.sh     (default: <repo>/exe/adam_prism_cpu)
# The backend is the executable's suffix; a leg that names other backends is skipped.
set -euo pipefail

REGRESSION_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$REGRESSION_DIR/../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"
BACKEND="${EXE##*_}"
BASE="$REGRESSION_DIR/rmf/input.ini"

# name | backends | sed script applied to rmf/input.ini | expected message (grep -F)
LEGS=(
   # issue #61: hyperbolic cleaning has no phi/psi slot without constrained transport (var_J* were left undefined)
   "hyperbolic-no-ct|cpu fnl|s/^divergence_correction *=.*/divergence_correction = hyperbolic/|hyperbolic needs"
   # issue #61 / 68457200: the Poisson projection is not implemented on either backend
   "poisson|cpu fnl|s/^constrained_transport *=.*/constrained_transport = B/;s/^divergence_correction *=.*/divergence_correction = poisson/|poisson is not implemented"
   # issue #55: field integrators whose objects realm_object never initializes (CPU) or whose device routines are empty
   "blanes-moan|cpu fnl|s/^scheme_time *=.*/scheme_time = blanes_moan/|Blanes-Moan time integration"
   "cfm|cpu fnl|s/^scheme_time *=.*/scheme_time = commutator_free_magnus/|Commutator-Free Magnus time integration"
   "leapfrog|cpu fnl|s/^scheme_time *=.*/scheme_time = leapfrog/|leapfrog time integration"
   "yoshida|fnl|s/^scheme = runge-kutta-ssp-54/scheme = runge-kutta-yoshida/|Yoshida Runge-Kutta time integration"
   # issue #56: initial conditions that set_initial_conditions never applied (the run started from the zero field), and
   # an RMF axis that left the component indices undefined
   "ic-riemann-problem|cpu fnl|s/^type *= *vacuum/type = riemann-problem/|riemann-problem is not implemented"
   "ic-magnetic-nozzle|cpu fnl|s/^type *= *vacuum/type = magnetic_nozzle/|is not implemented as an initial condition"
   "ic-unknown|cpu fnl|s/^type *= *vacuum/type = no_such_ic/|unknown [initial_conditions].(type)"
   "ic-rmf-axis|cpu fnl|s/^type *= *vacuum/type = rmf_field/;s/^external_fields_applied *=.*/&\nRMF_frequency = 1.0E6\nRMF_B_amplitude = 1.0E-3\nRMF_rotation_axis = W/|unknown [external_fields].(RMF_rotation_axis)"
   # issue #57: a BC no backend implements (FNL applied none), and values that used to fall back silently or with a
   # warning; a typo now stops, naming the option and the value
   "bc-radiative|cpu fnl|s/^type = Neumann/type = radiative/|radiative is not implemented"
   "scheme-time-typo|cpu fnl|s/^scheme_time *=.*/scheme_time = runge-kutta/|(scheme_time) \"runge-kutta\""
   "scheme-space-typo|cpu fnl|s/^scheme_space *=.*/scheme_space = fd_centred/|(scheme_space) \"fd_centred\""
   "reconstruction-typo|cpu fnl|s/^reconstruction_variables *=.*/reconstruction_variables = characteristichs/|(reconstruction_variables) \"characteristichs\""
   "ct-typo|cpu fnl|s/^constrained_transport *=.*/constrained_transport = DD/|(constrained_transport) \"DD\""
   "dc-typo|cpu fnl|s/^divergence_correction *=.*/divergence_correction = hyperbolc/|(divergence_correction) \"hyperbolc\""
   "ef-typo|cpu fnl|s/^external_fields_applied *=.*/external_fields_applied = rmff/|(external_fields_applied) \"rmff\""
   "fwl-flag-typo|cpu fnl|s/^x_minus_layer *=.*/x_minus_layer = yse/|(x_minus_layer) \"yse\""
   "ib-file-off|cpu fnl|s/^number    = 0/number    = 1/;s/^n_eikonal = 2/&\n[solid_1]\nname = s\nbc_type = 1\ndefinition = file.off/|(definition) \"file.off\""
   "ib-axis|cpu fnl|s/^number    = 0/number    = 1/;s/^n_eikonal = 2/&\n[solid_1]\nname = s\nbc_type = 1\ndefinition = analytical_circle\ncircle_center_x = 0.0\ncircle_center_y = 0.0\ncircle_center_z = 0.0\ncircle_radius = 0.01\ncircle_axis = w/|(circle_axis) \"w\""
)

fail=0
ran=0
for leg in "${LEGS[@]}"; do
   IFS='|' read -r name backends script message <<< "$leg"
   if [[ " $backends " != *" $BACKEND "* ]]; then
      echo ">> [refusals/$BACKEND] $name: skipped (applies to $backends)"
      continue
   fi
   # inside the rmf case, like every case's work-*: a directory with an input.ini right under the suite would be
   # picked up by run.sh as a case
   work="$REGRESSION_DIR/rmf/work-refusal-$BACKEND-$name"
   rm -rf "$work" && mkdir -p "$work"
   sed -e "$script" "$BASE" > "$work/input.ini"
   if cmp -s "$BASE" "$work/input.ini"; then
      echo "FAIL [refusals/$BACKEND] $name: the sed script changed nothing (rmf/input.ini drifted?)"
      fail=1; continue
   fi
   rc=0
   ( cd "$work" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 ) || rc=$?
   find "$work" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' \) -delete
   ran=$((ran + 1))
   if [[ $rc -eq 0 ]]; then
      echo "FAIL [refusals/$BACKEND] $name: RAN instead of refusing"
      fail=1
   elif ! grep -aqF "$message" "$work/run.log"; then
      echo "FAIL [refusals/$BACKEND] $name: died (exit $rc) WITHOUT the message '$message': a crash, not a refusal"
      fail=1
   else
      echo ">> [refusals/$BACKEND] $name: clean refusal (exit $rc)"
   fi
done

if [[ $fail -ne 0 ]]; then
   echo "FAIL [refusals/$BACKEND] init refusals"
   exit 1
fi
echo "PASS [refusals/$BACKEND] init refusals ($ran legs)"
