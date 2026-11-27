#!/usr/bin/env bash
# Device array-constructor lint gate (issue #157, Giacomo Rossi's MI300X finding).
#
# flang (AMD AFAR drop 23.2.1, OpenMP offload) lowers an array constructor with variable elements, assigned to a
# fixed-size local (`base = [i, j, k]`), through a heap temporary: in device code a malloc/free pair of the OpenMP
# device runtime, whose heap is global and serialised by atomics. One such line in the seam fill's MP bound made the
# MI300X 2:1 ghost fill 7.8x slower than the H100; one in `compute_eigenvectors` doubled the face-flux kernel.
# nvfortran allocates nothing for the same line, so NVIDIA timings never show it. Write the elements one by one
# (`base(1) = i ; base(2) = j ; base(3) = k`): the same values, bitwise. Constructors of literals and named constants
# (`[IA_U, IA_V, IA_W]`, `[1_I4P, 0_I4P]`) are folded and allocate nothing: they pass.
#
# Device code, for this gate: the procedures that carry `!$acc routine` or `!$omp declare target`, and the loop nests
# opened by a kernel directive (`!$acc parallel|kernels|serial`, `!$omp OMPLOOP|target`) in any procedure. A
# constructor fails when its elements contain a lowercase name (a variable, an array element, a function result).
# Comment lines and quoted strings are ignored. Scope: the maintained library and applications (src/lib,
# src/app/{flume,prism,ripple}), the CUDA Fortran and OpenMP-offload twins of NASTO excluded (nvf/, gmp/).
#
# Invoked from the top of the regression sweeps (run.sh); standalone use:
#   ./src/tests/lint/check-device-array-constructors.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

mapfile -t files < <(find "$REPO_ROOT/src/lib" "$REPO_ROOT/src/app/flume" "$REPO_ROOT/src/app/prism" \
                          "$REPO_ROOT/src/app/ripple" \( -name '*.F90' -o -name '*.INC' \) \
                          -not -path '*/nvf/*' -not -path '*/gmp/*' 2>/dev/null | sort)

hits="$(awk '
function strip(s) {                          # drop quoted strings and the trailing comment of a code line
   gsub(/"[^"]*"/, "", s) ; gsub(/\047[^\047]*\047/, "", s)
   sub(/![^$].*$/, "", s) ; sub(/!$/, "", s)
   return s
}
function flush_proc(   n) {                  # report the routine-scope hits of a device procedure
   if (dev) for (n = 1; n <= nh; n++) print hit[n]
   nh = 0 ; dev = 0
}
FNR == 1 { flush_proc() ; depth = 0 ; armed = 0 ; region = 0 }
{
   raw = $0 ; low = tolower(raw)
   if (low ~ /^[[:space:]]*!\$(acc|omp)/) {
      if (low ~ /^[[:space:]]*!\$acc[[:space:]]+routine/ || low ~ /^[[:space:]]*!\$omp[[:space:]]+declare[[:space:]]+target/) dev = 1
      if (low ~ /^[[:space:]]*!\$acc[[:space:]]+(parallel|kernels|serial)/ || low ~ /^[[:space:]]*!\$omp[[:space:]]+(omploop|target)/) {
         if (low !~ /end/ && depth == 0) armed = 1
      }
      next
   }
   if (low ~ /^[[:space:]]*!/) next
   code = strip(raw) ; lcode = tolower(code)
   if (lcode ~ /^[[:space:]]*((pure|elemental|impure|recursive|module|simple)[[:space:]]+)*(subroutine|function)[[:space:]]/ || \
       lcode ~ /^[[:space:]]*((pure|elemental|impure|recursive|module)[[:space:]]+)*(real|integer|logical)(\([^)]*\))?[[:space:]]+function[[:space:]]/) {
      flush_proc() ; depth = 0 ; armed = 0
   }
   # loop depth of an armed kernel nest
   ndo = gsub(/(^|;)[[:space:]]*(([a-z_][a-z0-9_]*)[[:space:]]*:[[:space:]]*)?do([[:space:]]|$)/, "&", lcode)
   nend = gsub(/end[[:space:]]*do/, "&", lcode)
   inkernel = (armed || depth > 0)
   s = code
   while (match(s, /\[[^][]*\]/)) {
      body = substr(s, RSTART + 1, RLENGTH - 2) ; s = substr(s, RSTART + RLENGTH)
      sub(/^.*::/, "", body)
      gsub(/[0-9.][eEdD][+-]?[0-9]+/, "0", body) ; gsub(/\.[a-zA-Z]+\./, "", body) # exponents, .true., operators
      if (body ~ /[a-z]/) {
         line = FILENAME ":" FNR ": " raw
         if (inkernel) print line ; else hit[++nh] = line
      }
   }
   if (armed && ndo > 0) { depth = ndo - nend ; armed = 0 }           # the kernel nest opens
   else if (depth > 0) { depth += ndo - nend ; if (depth < 0) depth = 0 } # inside it; host loops are not counted
}
/^[[:space:]]*end[[:space:]]*(subroutine|function)/ { flush_proc() }
END { flush_proc() }
' "${files[@]}" | sed "s|$REPO_ROOT/||" | sort -u)"

if [[ -n "$hits" ]]; then
   echo "FAIL [lint/device-array-constructors] array constructor with variable elements in device code (flang allocates"
   echo "     it on the device heap: assign the elements one by one, issue #157):"
   echo "$hits"
   exit 1
fi
echo "PASS [lint/device-array-constructors] no variable array constructors in device code"
