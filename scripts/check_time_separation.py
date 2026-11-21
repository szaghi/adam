#!/usr/bin/env python3
"""Check that no Runge-Kutta combination is written in ADAM's sources: the time integrators are SELAH's (#132, 6.2).

Why this script exists: M7 (#132) moved every time integrator ADAM ran onto SELAH and deleted ADAM's own (rk_object,
rk_fnl_object and its kernels, the NVF and GMP twins, PRISM's boundary clone). A Runge-Kutta stage or update written
again in src/app or src/lib would be a second source of truth beside SELAH's tableaux, drifting from them as the old
ones did (#128). This check stops one from coming back unnoticed, as SELAH's scripts/check_separation.py does for the
space operators in GIDEON and SETH. A finding is an assignment that

- multiplies the time step by a stage coefficient: dt * alph(...), beta(s) * dt, ... (alph, beta, gamm, the Butcher
  arrays; ark, brk, crk, the low-storage ones; ssa, ssb, the splitting ones; any of them as a device copy, *_gpu),
- multiplies the time step by any array element indexed first by a stage (dt * c(r), b(s) * dt, dt*s_coeffs(s,ss):
  a tableau row under another name, #134 D3),
- applies a literal fraction or multiple to the time step (0.5 * dt, dt / 2, 2 * dt: a half kick, a Strang
  sub-step, a leapfrog step),
- reads the stepper's tableau directly (%st%bt) outside the adapter (adam_stepper_object), or
- writes an SSP-RK3 stage or the classic RK4 sum with literal coefficients (0.75 u + 0.25 ..., u/3 + 2/3 ...,
  k1 + 2 k2 + 2 k3 + k4).

The patterns are regular expressions, not a Fortran parser: a statement is the joined continuation lines, comments
cut; the sources are the .F90 and .INC files of src/app and src/lib (the archived apps excluded). A remaining
statement is waived in WAIVERS by file name and a pattern it matches, with the reason and the issue that removes it.
`--self-test` plants each kind in a scratch file and checks that the rule catches it, and that a plain dt product (a
CFL, a flux) does not.

Usage: python3 scripts/check_time_separation.py [--root DIR] [--self-test]     (exit status 1 on a finding)
"""

from __future__ import annotations

import argparse
import re
import sys
import tempfile
from pathlib import Path

COEFF = r'(?:\w+%)?(?:alph|beta|gamm|ark|brk|crk|ssa|ssb)(?:_gpu)?\s*\('
STAGE = re.compile(r'\bdt\w*\s*\*\s*' + COEFF + r'|' + COEFF + r'[^()]*\)\s*\*\s*dt\b', re.IGNORECASE)
DT = r'(?:\w+%)*dt\b'
STAGE_INDEXED = r'[\w%]+\s*\(\s*(?:s|ss|r|stage|istage)\s*[,)]'
STAGE_ROW = re.compile(DT + r'\s*\*\s*' + STAGE_INDEXED + r'|' + STAGE_INDEXED + r'[^=]*?\*\s*' + DT, re.IGNORECASE)
SPLIT = re.compile(r'\b(?:0\.5|2\.?|2\.0)\w*\s*\*\s*' + DT + r'|' + DT + r'\s*\*\s*0\.5\w*\b|' + DT + r'\s*/\s*2(?:\.\w*)?\b', re.IGNORECASE)
TABLEAU = re.compile(r'%st%bt\b', re.IGNORECASE)
ADAPTER = 'adam_stepper_object.F90'  # the one file that may read the stepper's tableau
RK3 = re.compile(r'0\.75_R\w*\s*\*.*0\.25_R\w*|/\s*3\._R\w*\s*\+\s*2\._R\w*\s*/\s*3\._R\w*', re.IGNORECASE)
RK4 = re.compile(r'\bk1\b\s*\+\s*2\w*\s*\*\s*k2\b', re.IGNORECASE)

# file name -> [(pattern a waived statement matches, reason and the issue that removes it)]
_PIC = "PRISM's particles (PIC) advance in step with the field stepper on its SELAH tableau (szaghi/SELAH#4)"
_CC = 'the charge-conserving PIC current combines the tableau rows (szaghi/SELAH#4)'
_LIM = ("the positivity limiter's gravity work rho g . (u + g h/2) of one forward-Euler step h, the stage's "
        "(SELAH's stage_euler_fraction, #134 D4): the exact work of a constant-g kick over h, not a split (kept)")
_T0 = 'the PIC t0 centring drift x -+ dt/2 v builds J(t0), not a step (kept)'
WAIVERS: dict[str, list[tuple[str, str]]] = {
    'adam_prism_rk_pic_object.F90': [('.', _PIC)],
    'adam_prism_fnl_rk_pic_object.F90': [('.', _PIC + ', device twin')],
    'adam_prism_cpu_object.F90': [(r'%st%bt%[ab]\b', _CC), (r'q_pic_next.*coeff\(r\)', _CC),
                                  (r'q_pic_(minus|plus)|0\.5_R8P\s*\*\s*self%time%dt\s*\*\s*self%q_pic', _T0)],
    'adam_prism_fnl_object.F90': [(r'%st%bt%[ab]\b', _CC)],
    'adam_prism_fnl_pic_conserving_object.F90': [(r'coeff|\bc\(r\)', _CC), (r'half_dt|0\.5_R8P\s*\*\s*dt', _T0)],
    'adam_flume_cpu_positivity_kernels_agnostic.INC': [(r'0\.5_R8P\s*\*\s*dt\s*\*\s*q', _LIM)],
    'adam_flume_fnl_positivity_kernels_agnostic.INC': [(r'0\.5_R8P\s*\*\s*dt\s*\*\s*q', _LIM)],
}


def statements(path: Path):
    """(line number, statement) of a Fortran source, continuation lines joined, comments cut."""
    buf, start = '', 0
    for n, raw in enumerate(path.read_text(errors='replace').splitlines(), 1):
        code = raw.split('!')[0].rstrip()
        if not buf:
            start = n
        if code.endswith('&'):
            buf += code[:-1] + ' '
            continue
        stmt = (buf + code.lstrip('& ')).strip()
        buf = ''
        for part in stmt.split(';'):
            if part.strip():
                yield start, part.strip()


def waived(path: Path, stmt: str) -> bool:
    return any(re.search(p, stmt) for p, _ in WAIVERS.get(path.name, []))


def combinations(path: Path) -> list[str]:
    out = []
    for n, s in statements(path):
        tableau = TABLEAU.search(s) and path.name != ADAPTER
        if not tableau and ('=' not in s or re.match(r'(?i)(do|if|call|where)\b', s)):
            continue
        rhs = s.split('=', 1)[1] if '=' in s else s
        kind = ('tableau read' if tableau else
                'Runge-Kutta combination' if STAGE.search(rhs) or RK3.search(rhs) or RK4.search(rhs) else
                'stage-row combination' if STAGE_ROW.search(rhs) else
                'splitting fraction' if SPLIT.search(rhs) else '')
        if kind and not waived(path, s):
            out.append(f'{path}:{n}: {kind}: {s[:110]}')
    return out


ARCHIVED = ('nasto', 'chase', 'patch', 'ascot')  # archived apps, not maintained


def sources(root: Path) -> list[Path]:
    return sorted(p for d in ('src/app', 'src/lib') for p in (root / d).rglob('*')
                  if p.suffix in ('.F90', '.f90', '.INC') and not set(ARCHIVED) & set(p.parts))


def self_test() -> int:
    """Plant each kind of combination in scratch files: each must be caught, a plain dt product must not."""
    planted = {
        'stage.F90': ('      q_rk(v,i,j,k,b,s) = q_rk(v,i,j,k,b,s) + dt * alph(s,ss) * q_rk(v,i,j,k,b,ss)\n', 1),
        'update.F90': ('      q(v,i,j,k,b) = q(v,i,j,k,b) + &\n                     beta_gpu(s) * dt * dq(v,i,j,k,b)\n', 1),
        'lowstorage.F90': ('      q = self%ark(s) * q0 + self%brk(s) * q + self%crk(s) * dt * dq\n', 1),
        'rk3.F90': ('      u2 = 0.75_R8P*u + 0.25_R8P*(u1 + dt*r)\n', 1),
        'rk4.F90': ('      u = u + dt*(k1 + 2*k2 + 2*k3 + k4)/6\n', 1),
        'row.F90': ('      q_next(v,p) = q_next(v,p) + self%time%dt * coeff(r) * k_stages(v,p,r)\n', 1),
        'kick.F90': ('      v = v + 0.5_R8P * dt * a\n', 1),
        'leapfrog.F90': ('      q_new(v,i) = q_old(v,i) + 2._R8P * dt * dq(v,i)\n', 1),
        'tableau.F90': ('      alph = ts%st%bt%a\n', 1),
        'cfl.F90': ('      dt = cfl * dx / (abs(u) + c)\n      flux = dt * f(v,i) / dx\n      x(i) = x0 + dt * v(i)\n', 0),
    }
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for name, (src, expect) in planted.items():
            p = Path(d) / name
            p.write_text('program p\n' + src + 'endprogram p\n')
            got = len(combinations(p))
            ok = got == expect
            bad += not ok
            print(f'self-test {name:15s} findings {got} (expected {expect}): {"ok" if ok else "FAIL"}')
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--root', default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument('--self-test', action='store_true')
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    root = Path(args.root)
    files = sources(root)
    findings = [f for p in files for f in combinations(p)]
    for f in findings:
        print(f)
    print(f'check_time_separation: {len(files)} sources, {len(findings)} finding(s), '
          f'{sum(len(w) for w in WAIVERS.values())} waiver(s) in {len(WAIVERS)} file(s)')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
