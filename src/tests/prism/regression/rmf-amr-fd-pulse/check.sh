#!/usr/bin/env bash
# rmf-amr-fd-pulse structural + divergence check (issue #21, N4 acceptance).
#
# The SOURCE-FREE two-invariant oracle: gaussian EM pulse IC, zero coils, same
# 2:1 AMR seam as rmf-amr-fd. With J = 0 the marker-off control must conserve
# BOTH discrete divergences to round-off. With the seam present:
#   - div(D) stays at round-off STRUCTURALLY: the +y transverse pulse has
#     div(D) = dDz/dz with no z-variation, and any partition-of-unity ghost
#     fill (injection, compatible, tricubic) preserves that exactly;
#   - div(B) is a TRUNCATION-ORDER quantity of the tricubic seam ghost fill
#     ([amr] seam_ghost_fill, #21 N2) over the full ghost slab (#21 N3.5),
#     pinned to a baseline and, on demand, asserted to CONVERGE under
#     refinement — this case is the valid convergence instrument (no source;
#     rmf-amr-fd's coil contaminates matched-time metrics, see #21 N3).
#
# ACCEPTANCE (redefined per #21 — the historical "seam div(B) <= 1e-13"
# pointwise criterion is RETIRED, see the regression README and #21 §1):
#   1. control: max|div(D)| and max|div(B)| <= DIV_TOL (no seam at all);
#   1b. [#29] JUMP-ISOLATION: uniform-refined (whole domain -> 1:1 same-level
#      boundaries, NO 2:1 jump) must also keep max|div(D)| AND max|div(B)| <=
#      DIV_TOL. This is the CONTROLLED COMPARISON: identical pulse physics +
#      identical level-2 resolution + the same many-block ghost-exchange as the
#      AMR case, differing ONLY in whether a resolution jump exists. 1:1 must be
#      machine-zero; the 2:1 seam (leg 3) is ~1.26E+01 — isolating the JUMP as
#      the sole div(B) source, not the seam machinery. (issue #29 accept-truncation.)
#   2. seam: max|div(D)| <= DIV_TOL (structural invariant);
#   3. seam: max|div(B)| within RTOL of the pinned baseline;
#   4. [--convergence] seam div(B) at matched time over a 2-rung ladder
#      (ni=16 it_max=5 vs ni=32 it_max=10; dt scales with CFL so final times
#      coincide): p_obs = log2(E16/E32) >= P_MIN. Measured +1.15 at N3.5
#      (second ratio +1.62 at ni=64, trending to the predicted +2).
#
# Baseline provenance: N3.5 binary (stale-ghost-layer fix), default tricubic
# fill, fdv_order 6, ni=16, it_max=5, -np 1 — max over the history (= step 5).
#
# Usage: ./check.sh [--build] [--convergence]
#
# PRISM_EXE (issue #22, GA4): override the executable under test, e.g.
#   PRISM_EXE=$REPO/exe/adam_prism_fnl ./check.sh --convergence
# The caller owns the matching environment (FNL: nvhpc mpirun on PATH + the
# WSL UCX knobs of issue #12) and the build (--build always builds the CPU
# default, never the override). Baselines are CPU-pinned; the ±5% band is the
# cross-backend/compiler-noise allowance; p_obs is compiler-independent.
#
# RK-scheme contract legs (issue #25, DEFAULT-ON):
#   1. LS refusal: a low-storage RK scheme (runge-kutta-1) on this SEAM case
#      must REFUSE CLEANLY — the staged forest path is SSP-only by contract
#      (stages_per_step_forest). The leg asserts non-zero exit AND the
#      contract message in the log: a segfault also exits non-zero, and the
#      leg must distinguish clean refusal from the pre-#25 crash/UB.
#   2. SSP-11 instrument: runge-kutta-ssp-11 (forward Euler in SSP form, the
#      1-stage debugging discriminator) must run the seam case to 100% with
#      div(D) structural zero and div(B) on its own CPU-pinned baseline.
set -euo pipefail

CASE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$CASE_DIR/../../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_cpu}"

DIV_TOL="1.0E-13"            # round-off ceiling for the conserved invariants
SEAM_B_BASELINE="1.265004E+01" # pinned seam max|div(B)| (see provenance above)
SEAM_RTOL="0.05"             # 5% relative band around the baseline
P_MIN="0.8"                  # convergence-order floor (measured ~ +1.15)
SEAM_B_SSP11_BASELINE="1.181251E+01" # pinned seam max|div(B)| under runge-kutta-ssp-11 (#25 F2; CPU
                                     # F3-fix binary, 2026-07-04: ni=16, it_max=5, -np 1, tricubic fill)
LS_REFUSAL_MSG="not stage-splittable" # the stages_per_step_forest contract message (#25 F1)

do_build=0 ; do_convergence=0
for arg in "$@"; do
   case "$arg" in
      --build)       do_build=1 ;;
      --convergence) do_convergence=1 ;;
      *) echo "ERROR: unknown flag $arg (use --build / --convergence)" >&2; exit 2 ;;
   esac
done
if [[ $do_build -eq 1 ]]; then
   echo ">> building prism-cpu-gnu"
   (cd "$REPO_ROOT" && fobis build --mode prism-cpu-gnu)
fi
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found — run with --build" >&2; exit 2; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not on PATH" >&2; exit 2; }

# divergence-history columns: it nblk time divD divB divJ
max_div_d()  { grep '^+' "$1" | awk '{v=$4; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }
max_div_b()  { grep '^+' "$1" | awk '{v=$5; if(v<0)v=-v; if(v>m)m=v} END{printf "%.6E", m+0}'; }
last_div_b() { grep '^+' "$1" | awk 'END{v=$5; if(v<0)v=-v; printf "%.6E", v+0}'; }

run_in() { # workdir, ini-transform-sed
   local wd="$1" sed_expr="$2"
   rm -rf "$wd" && mkdir -p "$wd"
   sed "$sed_expr" "$CASE_DIR/input.ini" > "$wd/input.ini"
   ( cd "$wd" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 )
   find "$wd" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' \) -delete
}

fail=0
HIST="rmf_amr_fd_pulse_regression-divergence_history.dat"

CTRL="$CASE_DIR/work-cpu-ctrl"
echo ">> [rmf-amr-fd-pulse] running marker-disabled control (no seam, no sources)"
run_in "$CTRL" 's/markers_number = 1/markers_number = 0/'
if grep -qiE 'error|abort| nan |segfault' "$CTRL/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] control reported an error/abort/NaN"; fail=1
fi
ctrl_divd="$(max_div_d "$CTRL/$HIST")"
ctrl_divb="$(max_div_b "$CTRL/$HIST")"
echo ">> [rmf-amr-fd-pulse] control max|div(D)| = $ctrl_divd (expect <= $DIV_TOL, J=0)"
echo ">> [rmf-amr-fd-pulse] control max|div(B)| = $ctrl_divb (expect <= $DIV_TOL)"
if ! awk "BEGIN{exit !($ctrl_divd <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] control div(D) above round-off — source-free d/dt divD = -divJ = 0 broken"; fail=1
fi
if ! awk "BEGIN{exit !($ctrl_divb <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] control div(B) above round-off — matched-stencil identity broken without a seam"; fail=1
fi

# --- JUMP-ISOLATION leg (issue #29): 1:1 same-resolution boundaries do NOT leak; only the 2:1 JUMP does.
# Refine the WHOLE domain (box_xmax 0.0 -> +0.06) so every block reaches level 2 UNIFORMLY: many
# same-level (1:1) block boundaries, ghost-exchanged every step, but NO 2:1 resolution jump. The
# mimetic div_h(curl_h)=0 identity survives 1:1 adjacency (both sides same Delta-x -> operators
# commute), so div(B) must stay at round-off. Paired with the 2:1 seam leg below (~1.26E+01), this
# ISOLATES the resolution jump as the sole div(B) source — the controlled comparison that proves the
# leak is the JUMP, not the seam machinery/ghost-exchange (measured 64 blocks, div(B)=0.0 to 100 steps).
UNIF="$CASE_DIR/work-cpu-uniform"
echo ">> [rmf-amr-fd-pulse] running uniform-refined leg (whole domain -> 1:1 boundaries, NO 2:1 jump)"
run_in "$UNIF" 's/^box_xmax     =  0.0/box_xmax     =  0.06/'
if grep -qiE 'error|abort| nan |segfault' "$UNIF/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] uniform-refined leg reported an error/abort/NaN"; fail=1
fi
unif_divd="$(max_div_d "$UNIF/$HIST")"
unif_divb="$(max_div_b "$UNIF/$HIST")"
echo ">> [rmf-amr-fd-pulse] 1:1 uniform max|div(D)| = $unif_divd (expect <= $DIV_TOL)"
echo ">> [rmf-amr-fd-pulse] 1:1 uniform max|div(B)| = $unif_divb (expect <= $DIV_TOL — 1:1 boundaries must NOT leak)"
if ! awk "BEGIN{exit !($unif_divd <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] uniform-refined div(D) above round-off — source-free invariant broke"; fail=1
fi
if ! awk "BEGIN{exit !($unif_divb <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] uniform-refined div(B) above round-off — a 1:1 SAME-LEVEL boundary is"
   echo "                        leaking div(B). The mimetic identity must hold at 1:1 adjacency; a"
   echo "                        leak here means the seam machinery/ghost-exchange itself regressed"
   echo "                        (NOT the 2:1-jump truncation source, which is separate and expected)."
   fail=1
fi

WORK="$CASE_DIR/work-cpu"
echo ">> [rmf-amr-fd-pulse] running marker-active case (2:1 seam)"
run_in "$WORK" 's/^$/&/'   # identity (copy as-is)
if grep -qiE 'error|abort| nan |segfault' "$WORK/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] run reported an error/abort/NaN"; fail=1
fi
if ! grep -qE 'progress:[[:space:]]*100%' "$WORK/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] fd_centered time loop did not reach 100%"; fail=1
fi
seam_divd="$(max_div_d "$WORK/$HIST")"
seam_divb="$(max_div_b "$WORK/$HIST")"
echo ">> [rmf-amr-fd-pulse] seam   max|div(D)| = $seam_divd (structural invariant: <= $DIV_TOL)"
echo ">> [rmf-amr-fd-pulse] seam   max|div(B)| = $seam_divb (baseline $SEAM_B_BASELINE, rtol $SEAM_RTOL)"
if ! awk "BEGIN{exit !($seam_divd <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] seam div(D) above round-off — the partition-of-unity structural"
   echo "                        invariant broke (ghost fill no longer reproduces constants?)"
   fail=1
fi
if ! awk "BEGIN{d=($seam_divb-$SEAM_B_BASELINE)/$SEAM_B_BASELINE; if(d<0)d=-d; exit !(d<=$SEAM_RTOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] seam max|div(B)| moved off the pinned baseline"
   echo "                        (a DROP means an improvement landed — rebaseline deliberately;"
   echo "                         a RISE means the seam fill or the exchange regressed)"
   fail=1
fi

# --- RK-scheme contract legs (issue #25, default-on) ---
RK1="$CASE_DIR/work-rk1-refusal"
echo ">> [rmf-amr-fd-pulse] running LS-refusal leg (runge-kutta-1 on the seam: must refuse cleanly)"
rk1_rc=0
rm -rf "$RK1" && mkdir -p "$RK1"
sed 's/^scheme = runge-kutta-ssp-54/scheme = runge-kutta-1/' "$CASE_DIR/input.ini" > "$RK1/input.ini"
( cd "$RK1" && timeout 300 mpirun -np 1 "$EXE" > run.log 2>&1 ) || rk1_rc=$?
find "$RK1" -type f \( -name '*.h5' -o -name '*.fbd' -o -name '*.xdmf' \) -delete
if [[ $rk1_rc -eq 0 ]]; then
   echo "FAIL [rmf-amr-fd-pulse] runge-kutta-1 RAN on the staged (seam) path — the SSP-only"
   echo "                        contract of stages_per_step_forest is not enforced"
   fail=1
elif ! grep -q "$LS_REFUSAL_MSG" "$RK1/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] runge-kutta-1 died WITHOUT the contract message ('$LS_REFUSAL_MSG')"
   echo "                        — crash/UB instead of the clean stages_per_step_forest refusal"
   fail=1
else
   echo ">> [rmf-amr-fd-pulse] LS refusal: clean error-stop with contract message (exit $rk1_rc)"
fi

SSP11="$CASE_DIR/work-ssp11"
echo ">> [rmf-amr-fd-pulse] running SSP-11 instrument leg (forward Euler in SSP form, seam case)"
run_in "$SSP11" 's/^scheme = runge-kutta-ssp-54/scheme = runge-kutta-ssp-11/'
if grep -qiE 'error|abort| nan |segfault' "$SSP11/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] SSP-11 run reported an error/abort/NaN"; fail=1
fi
if ! grep -qE 'progress:[[:space:]]*100%' "$SSP11/run.log"; then
   echo "FAIL [rmf-amr-fd-pulse] SSP-11 time loop did not reach 100%"; fail=1
fi
ssp11_divd="$(max_div_d "$SSP11/$HIST")"
ssp11_divb="$(max_div_b "$SSP11/$HIST")"
echo ">> [rmf-amr-fd-pulse] SSP-11 seam max|div(D)| = $ssp11_divd (structural invariant: <= $DIV_TOL)"
echo ">> [rmf-amr-fd-pulse] SSP-11 seam max|div(B)| = $ssp11_divb (baseline $SEAM_B_SSP11_BASELINE, rtol $SEAM_RTOL)"
if ! awk "BEGIN{exit !($ssp11_divd <= $DIV_TOL)}"; then
   echo "FAIL [rmf-amr-fd-pulse] SSP-11 seam div(D) above round-off"; fail=1
fi
if ! awk "BEGIN{d=($ssp11_divb-$SEAM_B_SSP11_BASELINE)/$SEAM_B_SSP11_BASELINE; if(d<0)d=-d; exit !(d<=$SEAM_RTOL)}" 2>/dev/null; then
   echo "FAIL [rmf-amr-fd-pulse] SSP-11 seam max|div(B)| off its pinned baseline"; fail=1
fi

if [[ $do_convergence -eq 1 ]]; then
   echo ">> [rmf-amr-fd-pulse] --convergence: 2-rung matched-time ladder (ni=16 vs ni=32)"
   LAD16="$CASE_DIR/work-conv-16"
   LAD32="$CASE_DIR/work-conv-32"
   trim='s/^it_save                = 5/it_save                = 999999/;s/^restart_save           = 500/restart_save           = 999999/;s/^save_residual_fields   = .true./save_residual_fields   = .false./;s/^save_divergence_fields = .true./save_divergence_fields = .false./'
   run_in "$LAD16" "$trim"
   run_in "$LAD32" "$trim;s/^ni     = 16/ni     = 32/;s/^nj     = 16/nj     = 32/;s/^nk     = 16/nk     = 32/;s/^it_max   = 5/it_max   = 10/"
   e16="$(last_div_b "$LAD16/$HIST")"
   e32="$(last_div_b "$LAD32/$HIST")"
   p_obs="$(awk "BEGIN{printf \"%.3f\", log($e16/$e32)/log(2)}")"
   echo ">> [rmf-amr-fd-pulse] matched-t seam div(B): E16 = $e16, E32 = $e32, p_obs = $p_obs (floor $P_MIN)"
   if ! awk "BEGIN{exit !($p_obs >= $P_MIN)}"; then
      echo "FAIL [rmf-amr-fd-pulse] seam div(B) no longer converges under refinement (p_obs = $p_obs < $P_MIN)"
      fail=1
   fi
fi

if [[ $fail -eq 0 ]]; then
   echo "PASS [rmf-amr-fd-pulse] two-invariant oracle: control at round-off, seam div(D) structural zero, div(B) on baseline"
   exit 0
else
   echo "FAIL [rmf-amr-fd-pulse] source-free seam divergence check failed"
   exit 1
fi
