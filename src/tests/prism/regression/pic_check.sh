#!/usr/bin/env bash
# Shared driver of the PIC regression cases (pic-periodic, pic-slab-pulse), called by each case's check.sh.
#
# The cases run the PIC path of the production runs: TSC, conserving current esirkepov_Q10 with the coherent gather, the
# t0 FLAIL solves of D and B, FD6, RK SSP54, macroparticles (particle_weight, 2 per cell per species), fixed seed, no
# pairing.
#   pic-periodic   uniform plasma, periodic on every face (closure grounded), plane wave of one box length,
#                  32x24x24 cells (dx = c/(16 w_pe)), 40 steps: periodic current fold and particle wrap on three axes.
#   pic-slab-pulse the Fresnel pilot (src/tests/prism/fnl/fd/PIC/Fresnel_test/pilot/pulse) at dx = c/(8 w_pe): slab
#                  3 c/w_pe thick, PEC along x, periodic y,z (closure isolated), gaussian_derivative_pulse through the
#                  slab, 384x24x24 cells, 1152 steps at CFL 0.25 (t w_pe 18; CFL 0.5 would give omega_max dt ~0.07 and
#                  the plasma_frequencies dt limit, see adam_prism_plasma_frequencies_object). Run to t w_pe 32 with
#                  three probe planes, the same
#                  input matches the exact cold slab like the pilot plane study at this dx (plane-mean L2 3.3 / 3.3 /
#                  6.2 % reflected / inside / transmitted; study 2.7 / 4.2 / 7.0 %; second order in dx).
# Times of one run (pic-periodic / pic-slab-pulse, 2026-11-23): CPU np 1 ~15 / ~20 s, FNL np 1 ~20 s / ~3 min (RTX 4070
# laptop), FNL np 2 ~1.2 / ~3 min (two A30, medooza); the restart legs about double them.
#
# The charge-conserving CPU current is single-block, single-rank (check_conserving_single_block), so a case runs as legs:
#   cpu  np 1, input.ini (one block), the first CPU_STEPS steps only  -> golden-check/cpu
#   fnl  np 1, input.ini (one block)                                  -> golden-check/fnl
#   fnl  np 2, input-np2.ini (the same grid as 64 blocks)             -> golden-check/fnl as well: on the FNL virtual block
#        the multi-block result equals the single-block one by construction (adam_prism_fnl_pic_conserving_object), up
#        to the order of the atomic adds, so the np 2 leg checks halos, migration and decomposition against np 1.
# Every leg also asserts Gauss's law at t0 (div D - rho on the step-0 checkpoint) and the charge conservation of the
# run (charge_conservation_history: Q_deposition_error and dQ_num_net relative to the injected charge, at round-off).
# Restart (issue #129, PIC_CHECK_RESTART=0 skips it): every leg runs again as it_max/2 steps, a restart, and the rest,
# compared with the uninterrupted run just made: the first dt after the restart (its omega_max dt, also where the CFL
# sets dt), the final .pbd of every rank (particles and plasma-frequency maxima) and the fields of the last checkpoint
# (div_* left out, #129 R5). CPU bit-identical; FNL within PIC_CHECK_RESTART_RTOL for the fields (default 1e-9) and
# PIC_CHECK_RESTART_PRTOL for the dt and the particles, variable by variable (default 1e-11): the atomic deposit is not
# reproducible, measured restart vs uninterrupted 4e-12 (fields, residuals) and 1.7e-14 (particles, pic-slab-pulse).
# Restart-only mode (PIC_CHECK_RESTART_ONLY=1, set by the case's check.sh; pml-restart): a case without particles runs
# only the uninterrupted run and its restart leg, cut at step PIC_CHECK_RESTART_AT instead of it_max/2; no charge, Gauss
# or golden check, no .pbd compared; the CPU legs are np 1 and np 2 (no single-block conserving current), every leg on
# input.ini when the case has no input-np2.ini.
# The goldens live in golden-check/, not golden/: run.sh runs every case with a golden/ at -np 2, which the CPU leg
# refuses. Their residuals are force-added (the repo-root .gitignore ignores src/tests/**/*.dat).
#
# Usage: check.sh <cpu|fnl> [--capture]
#   cpu: the gnu toolchain on PATH (mpirun of the CPU build); fnl: the nvhpc one. FNL_NP selects the FNL legs (default
#   "1 2"; on WSL use FNL_NP=1: np 2 shares one GPU and is 10-40x slower, run it on the cluster). --capture writes the
#   np 1 outputs into golden-check/<backend>/ instead of comparing (np 2 is never captured).
#   PRISM_EXE overrides exe/adam_prism_<backend>, PRISM_PY the python with h5py (default exe/.regression-venv),
#   PIC_CHECK_LAUNCH the launcher taking the rank count last (default "mpirun -np"; on Slurm "srun --mpi=pmix -u -n"),
#   PIC_CHECK_RESTART=0 skips the restart legs (they double the run time).
set -uo pipefail

CASE_DIR="$1"; BASE="$2"; CPU_STEPS="$3"; shift 3
BACKEND="${1:-}"; CAPTURE="${2:-}"
[[ "$BACKEND" == cpu || "$BACKEND" == fnl ]] || { echo "usage: check.sh <cpu|fnl> [--capture]" >&2; exit 2; }
CASE="$(basename "$CASE_DIR")"
REG_DIR="$(cd "$CASE_DIR/.." && pwd)"
REPO_ROOT="$(cd "$REG_DIR/../../../.." && pwd)"
EXE="${PRISM_EXE:-$REPO_ROOT/exe/adam_prism_$BACKEND}"
PY="${PRISM_PY:-$REPO_ROOT/exe/.regression-venv/bin/python}"
DIGEST="$REG_DIR/digest.py"
TIMEOUT="${PIC_CHECK_TIMEOUT:-1800}"
[[ -x "$EXE" ]] || { echo "ERROR: $EXE not found" >&2; exit 2; }
[[ -x "$PY" ]] || { echo "ERROR: $PY not found (run run.sh once, or set PRISM_PY)" >&2; exit 2; }
LAUNCH="${PIC_CHECK_LAUNCH:-mpirun -np}"
RESTART="${PIC_CHECK_RESTART:-1}"
RESTART_ONLY="${PIC_CHECK_RESTART_ONLY:-0}"
RESTART_AT="${PIC_CHECK_RESTART_AT:-}"
[[ "$RESTART_ONLY" == 1 && "$RESTART" == 0 ]] && { echo "SKIP [$(basename "$CASE_DIR")/$BACKEND] restart-only case, PIC_CHECK_RESTART=0"; exit 0; }
RESTART_RTOL=0; RESTART_PRTOL=0
[[ "$BACKEND" == fnl ]] && { RESTART_RTOL="${PIC_CHECK_RESTART_RTOL:-1e-9}"; RESTART_PRTOL="${PIC_CHECK_RESTART_PRTOL:-1e-11}"; }
command -v ${LAUNCH%% *} >/dev/null 2>&1 || { echo "ERROR: ${LAUNCH%% *} not on PATH" >&2; exit 2; }

if [[ "$BACKEND" == cpu ]]; then [[ "$RESTART_ONLY" == 1 ]] && nps="1 2" || nps=1; else nps="${FNL_NP:-1 2}"; fi
legs=(); for np in $nps; do
   [[ $np == 1 || ! -f "$CASE_DIR/input-np2.ini" ]] && legs+=("$np input.ini") || legs+=("$np input-np2.ini"); done

fail=0
for leg in "${legs[@]}"; do
   read -r np ini <<<"$leg"
   tag="$CASE/$BACKEND/np$np"
   wd="$CASE_DIR/work-check-$BACKEND-np$np"
   rm -rf "$wd" && mkdir -p "$wd"
   cp "$CASE_DIR/$ini" "$wd/input.ini"
   if [[ "$BACKEND" == cpu ]]; then   # the CPU leg: the first CPU_STEPS steps, one checkpoint at the end
      sed -i -E "s/^(it_max\s*=\s*)[0-9]+/\1$CPU_STEPS/; s/^((it_save|restart_save|residuals_save)\s*=\s*)[0-9]+/\1$CPU_STEPS/" \
         "$wd/input.ini"
   fi
   start=$(date +%s)
   ( cd "$wd" && timeout "$TIMEOUT" $LAUNCH "$np" "$EXE" > run.log 2>&1 )
   rc=$?
   echo ">> [$tag] run rc=$rc, $(( $(date +%s) - start )) s"
   if [[ $rc -ne 0 ]]; then echo "FAIL [$tag] the run exited non-zero, see $wd/run.log"; fail=1; continue; fi

   if [[ "$RESTART_ONLY" != 1 ]]; then   # the PIC checks and the golden (not in restart-only mode)
   # charge conservation of the run: |Q_deposition_error| and |dQ_num_net| <= 1e-12 of the injected charge
   "$PY" - "$wd/$BASE-charge_conservation_history.dat" "$wd/run.log" <<'EOF' || { echo "FAIL [$tag] charge conservation"; fail=1; }
import re, sys
import numpy as np
log = open(sys.argv[2]).read()
ne = float(re.search(r'Electrons number:\s*([-+0-9.Ee]+)', log).group(1))
w = float(re.search(r'Particle weight \(real particles per macroparticle\):\s*([-+0-9.Ee]+)', log).group(1))
q = 1.6e-19*ne*w                                   # charge of one species
names = [l[1:].strip() for l in open(sys.argv[1]) if l.startswith('#')]
h = np.loadtxt(sys.argv[1], comments='#', ndmin=2)
err = max(np.abs(h[:, names.index(c)]).max() for c in ('Q_deposition_error', 'dQ_num_net'))
print(f'>> charge: max |Q_deposition_error|, |dQ_num_net| = {err:.3e} C = {err/q:.2e} of one species')
sys.exit(0 if err <= 1e-12*q else 1)
EOF

   # Gauss law at t0 after the FLAIL solve, on the step-0 checkpoint: div D = rho to round-off, up to the checkerboard
   # chi = (-1)^(i+j+k) (out of reach of PRISM's centred div of any D), which the matched deposit leaves out of rho
   "$PY" - "$wd" "$BASE" <<'EOF' || { echo "FAIL [$tag] Gauss law at t0"; fail=1; }
import glob, sys
import h5py
import numpy as np
NGC, TOL = 6, 1e-12
blocks = []
for fn in sorted(glob.glob(f'{sys.argv[1]}/{sys.argv[2]}-000000000-proc*.h5')):
    f = h5py.File(fn, 'r')
    for b in sorted({k.rsplit('-', 1)[0] for k in f.keys()}):
        cut = lambda a: a[tuple(slice(NGC, n - NGC) for n in a.shape)]          # arrays (z, y, x): ghosts per axis
        h = f[f'{b}-dxdydz'][()][::-1]
        blocks.append((f[f'{b}-origin'][()][::-1] + NGC*h, h, cut(f[f'{b}-rho'][()]), cut(f[f'{b}-div_D'][()])))
x0 = np.min([b[0] for b in blocks], axis=0)
num = den = 0.0
chis = []
for start, h, rho, _ in blocks:
    off = np.rint((start - x0)/h).astype(int).sum()
    k, j, i = np.indices(rho.shape)
    chis.append((-1.0)**(i + j + k + off))
    num += (rho*chis[-1]).sum(); den += (chis[-1]**2).sum()
c = num/den
rmax = max(np.abs(b[2]).max() for b in blocks)
res = max(np.abs(b[3] - b[2] + c*chi).max() for b, chi in zip(blocks, chis))/rmax if rmax > 0 else np.inf
print(f'>> gauss t0: max|rho| {rmax:.4e}, checkerboard c/max|rho| {c/rmax:.2e}, max|div D - rho + c chi|/max|rho| {res:.2e}')
sys.exit(0 if res <= TOL else 1)
EOF

   shopt -s nullglob; h5=("$wd/$BASE"-*-proc*.h5); shopt -u nullglob
   h5=($(printf '%s\n' "${h5[@]}" | grep -E -- "-[0-9]+-proc[0-9]+\.h5$"))
   "$PY" "$DIGEST" write "$wd/digest.txt" "${h5[@]}" --case-dir "$CASE_DIR" >/dev/null || { echo "FAIL [$tag] digest"; fail=1; continue; }
   gdir="$CASE_DIR/golden-check/$BACKEND"
   if [[ "$CAPTURE" == --capture ]]; then
      [[ $np == 1 ]] || continue
      mkdir -p "$gdir" && cp "$wd/digest.txt" "$wd/$BASE-residuals.dat" "$gdir/"
      echo ">> [$tag] golden captured in $gdir"
      continue
   fi
   [[ -f "$gdir/digest.txt" ]] || { echo "FAIL [$tag] no golden at $gdir (run with --capture)"; fail=1; continue; }
   # digest.py compare, but the sum of a zero-mean field (J, rho, residuals) is a cancellation whose last digits follow
   # the order of the atomic adds: compared on its own scale, sqrt(count sum_sq) (Cauchy-Schwarz bound of |sum|)
   # np 2 runs 64 blocks against the one-block golden: the per-block metadata rows (origin, dxdydz, time_iteration) and
   # the residuals log (summed per block, with the block count) depend on the layout, not on the fields: skipped there
   "$PY" - "$wd/digest.txt" "$gdir/digest.txt" "$np" <<'EOF' || { echo "FAIL [$tag] field digest"; fail=1; }
import sys
import numpy as np
RTOL, ATOL = 1e-6, 1e-3                            # digest.py defaults
def rows(fn):
    skip = ('origin', 'dxdydz', 'time_iteration') if sys.argv[3] != '1' else ()
    return {tuple(l.split()[:2]): [float(x) for x in l.split()[2:]] for l in open(fn)
            if not l.startswith('#') and l.strip() and l.split()[1] not in skip}
p, g = rows(sys.argv[1]), rows(sys.argv[2])
bad = sorted(set(p) ^ set(g))
for k in sorted(set(p) & set(g)):
    (n, mn, mx, s, ss), (gn, gmn, gmx, gs, gss) = p[k], g[k]
    ok = n == gn and all(np.isclose(a, b, rtol=RTOL, atol=ATOL) for a, b in ((mn, gmn), (mx, gmx), (ss, gss)))
    ok = ok and abs(s - gs) <= RTOL*np.sqrt(gn*gss) + ATOL
    if not ok:
        bad.append(k); print('  MISMATCH', k, p[k], g[k])
print(f'>> digest: {len(g)} rows, {len(bad)} outside rtol {RTOL:g} (sum on the scale sqrt(count sum_sq))')
sys.exit(1 if bad else 0)
EOF
   [[ $np == 1 ]] && { "$PY" "$DIGEST" compare-residuals "$wd/$BASE-residuals.dat" "$gdir/$BASE-residuals.dat" \
      || { echo "FAIL [$tag] residuals"; fail=1; }; }
   fi

   # restart leg: half the steps, a restart, the rest; against the uninterrupted run of this leg
   [[ "$RESTART" == 0 || "$CAPTURE" == --capture ]] && continue
   rtag="$tag/restart"
   rd="$wd-restart"
   rm -rf "$rd" && mkdir -p "$rd"
   it_max=$(sed -nE 's/^it_max\s*=\s*([0-9]+).*/\1/p' "$wd/input.ini")
   half=${RESTART_AT:-$(( it_max/2 ))}
   cp "$wd/input.ini" "$rd/input.ini"
   sed -i -E "s/^(it_max\s*=\s*)[0-9]+/\1$half/; s/^(restart_save\s*=\s*)[0-9]+/\1$half/; /^restart\s*=/d" "$rd/input.ini"
   start=$(date +%s)
   ( cd "$rd" && timeout "$TIMEOUT" $LAUNCH "$np" "$EXE" > run1.log 2>&1 )
   rc=$?
   if [[ $rc -eq 0 ]]; then
      sed -i -E "s/^(it_max\s*=\s*)[0-9]+/\1$it_max/; s/^\[IO\]\s*$/[IO]\nrestart = .true./" "$rd/input.ini"
      ( cd "$rd" && timeout "$TIMEOUT" $LAUNCH "$np" "$EXE" > run2.log 2>&1 )
      rc=$?
   fi
   echo ">> [$rtag] $half + $(( it_max - half )) steps, rc=$rc, $(( $(date +%s) - start )) s"
   if [[ $rc -ne 0 ]]; then echo "FAIL [$rtag] the run exited non-zero, see $rd/run1.log, run2.log"; fail=1; continue; fi
   particles=(); [[ "$RESTART_ONLY" == 1 ]] && particles=(--no-particles)
   "$PY" "$REG_DIR/pic_restart_compare.py" "$wd" "$rd" "$BASE" "$half" "$it_max" "$RESTART_RTOL" "$RESTART_PRTOL" \
      "${particles[@]}" || { echo "FAIL [$rtag] the restarted run differs from the uninterrupted one"; fail=1; }
done
[[ $fail -eq 0 ]] && echo "PASS [$CASE/$BACKEND]" || echo "FAIL [$CASE/$BACKEND]"
exit $fail
