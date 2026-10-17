# Leg selection shared by the FLUME verification scripts: the `--legs` option.
#
# Why: the verification layer is run with a focused scope (the legs a change touches), not as a whole; the regression
# suite is the layer that runs on every change. Every verification script therefore selects its legs the same way, and
# a selection that names no real leg must fail: a misspelt leg that silently ran nothing would report a PASS.
#
# A script declares its legs and its default selection, collects the selection while parsing its arguments, then
# resolves it before it runs anything:
#
#   source "$VERIF_DIR/legs.sh"
#   legs_declare "sv3-2d sv3-3d sv4 sv5" "sv3-2d sv4"   # every leg (in run order), then the default selection
#   while [[ $# -gt 0 ]]; do
#      case "$1" in
#         --legs) legs_select "$2" ; shift 2 ;;          # older spellings (--leg, --cases) call legs_select too
#         ...
#      esac
#   done
#   legs_resolve check.sh                               # default, `all`, validation; exits 2 on a bad selection
#   if has_leg sv4; then ... fi                         # or: for leg in $LEGS; do ... done
#
# A selection is a list separated by commas or blanks; repeated options add up; `all` selects every leg. The legs run
# in declaration order whatever the order of the selection. On success LEGS holds the resolved legs, blank-separated,
# and the selection is printed (">> legs: ...") so every log records what it covered.

LEGS_ALL=""      # every leg of the script, in run order
LEGS_DEFAULT=""  # the legs run when no selection is given
LEGS_SELECTED="" # the selection collected from the command line
LEGS=""          # the resolved selection (legs_resolve)

legs_declare() { # legs_declare "<every leg>" ["<default legs>"]: the default is every leg when omitted
   LEGS_ALL="$(echo ${1//,/ })"
   LEGS_DEFAULT="$(echo ${2:-$1})"
   LEGS_DEFAULT="$(echo ${LEGS_DEFAULT//,/ })"
}

legs_select() { # legs_select "<legs>": add legs (commas or blanks) to the selection
   LEGS_SELECTED+=" ${1//,/ }"
}

legs_resolve() { # legs_resolve <script name>: resolve the selection into LEGS, or exit 2
   local script=$1 leg picked=""
   local accepted="${LEGS_ALL// /, }, all"
   if [[ -z $LEGS_ALL ]]; then
      echo "$script: no legs declared (legs_declare must come before legs_resolve)" >&2 ; exit 2
   fi
   for leg in ${LEGS_SELECTED:-$LEGS_DEFAULT}; do
      if [[ $leg == all ]]; then
         picked+=" $LEGS_ALL"
      elif [[ " $LEGS_ALL " == *" $leg "* ]]; then
         picked+=" $leg"
      else
         echo "$script: unknown leg '$leg' (accepted: $accepted)" >&2 ; exit 2
      fi
   done
   LEGS=""
   for leg in $LEGS_ALL; do
      [[ " $picked " == *" $leg "* ]] && LEGS+="${LEGS:+ }$leg"
   done
   if [[ -z $LEGS ]]; then
      echo "$script: no leg selected (accepted: $accepted)" >&2 ; exit 2
   fi
   echo ">> legs: $LEGS"
}

has_leg() { # has_leg <leg>: true when <leg> is selected
   [[ " $LEGS " == *" $1 "* ]]
}
