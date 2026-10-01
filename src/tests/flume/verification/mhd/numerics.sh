# Shared `--numerics` option of the MHD check.sh scripts (issue #47, M3-P3c): run a verification on the `weno-riemann`
# scheme instead of the inputs' WENO flux splitting, with the oracles of the splitting scheme (RV-5, RV-6, RV-7).
#
# Sourced by the check.sh of linear-wave, cpaw, riemann, rj2a, orszag-tang, rotor and field-loop, after VENV_PY is set:
#   NUMERICS                 the spec `solver[:recon[:correction[:sensor]]]` (e.g. `hlld`, `hll:primitive`); empty
#                            (the default) keeps the inputs untouched, scheme `weno`;
#   numerics_check           validate NUMERICS (exit 2 on a malformed spec), set the defaults of the omitted fields:
#                            recon `characteristic` (the MHD default, #47 D-5), correction `6th`, sensor `weno`;
#   numerics_tag             the work-directory suffix of the spec (`-hlld-characteristic-6th-weno`; empty for `weno`);
#   numerics_apply <ini>     rewrite the [numerics] block of a generated input (no-op without NUMERICS).
# shellcheck shell=bash

NUMERICS=""
NUMERICS_SOLVER="" NUMERICS_RECON="" NUMERICS_CORRECTION="" NUMERICS_SENSOR=""

numerics_check() {
   [[ -z $NUMERICS ]] && return 0
   local IFS=:
   # shellcheck disable=SC2086
   set -- $NUMERICS
   NUMERICS_SOLVER="$1" NUMERICS_RECON="${2:-characteristic}" NUMERICS_CORRECTION="${3:-6th}" NUMERICS_SENSOR="${4:-weno}"
   if [[ $# -gt 4 || ! $NUMERICS_SOLVER =~ ^(llf|hll|hlld)$ || ! $NUMERICS_RECON =~ ^(primitive|characteristic)$ \
         || ! $NUMERICS_CORRECTION =~ ^(none|4th|6th)$ || ! $NUMERICS_SENSOR =~ ^(weno|none)$ ]]; then
      echo "check.sh: --numerics '$NUMERICS' is not solver[:recon[:correction[:sensor]]] with solver llf|hll|hlld," \
           "recon primitive|characteristic, correction none|4th|6th, sensor weno|none" >&2
      exit 2
   fi
   NUMERICS="$NUMERICS_SOLVER:$NUMERICS_RECON:$NUMERICS_CORRECTION:$NUMERICS_SENSOR"
}

numerics_tag() {
   [[ -z $NUMERICS ]] && return 0
   echo "-$NUMERICS_SOLVER-$NUMERICS_RECON-$NUMERICS_CORRECTION-$NUMERICS_SENSOR"
}

numerics_apply() {
   [[ -z $NUMERICS ]] && return 0
   "$VENV_PY" - "$1" "$NUMERICS_SOLVER" "$NUMERICS_RECON" "$NUMERICS_CORRECTION" "$NUMERICS_SENSOR" <<'EOF'
import re, sys
path, solver, recon, corr, sensor = sys.argv[1:]
text = open(path).read()
new = (f"scheme_space             = weno-riemann\nreconstruction_variables = {recon}\n"
       f"riemann_solver           = {solver}\nflux_correction          = {corr}\n"
       f"flux_correction_sensor   = {sensor}\n")
text, n = re.subn(r"^scheme_space\s*=.*\n^reconstruction_variables\s*=.*\n", new, text, flags=re.M)
if n != 1:
    sys.exit(f"{path}: [numerics] block not found")
open(path, "w").write(text)
EOF
}
