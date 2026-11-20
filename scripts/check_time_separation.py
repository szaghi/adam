#!/usr/bin/env python3
"""Check that no Runge-Kutta combination is written in ADAM's sources: the time integrators are SELAH's (#132, 6.2).

Why this script exists: M7 (#132) moved every time integrator ADAM ran onto SELAH and deleted ADAM's own (rk_object,
rk_fnl_object and its kernels, the NVF and GMP twins, PRISM's boundary clone). A Runge-Kutta stage or update written
again in src/app or src/lib would be a second source of truth beside SELAH's tableaux, drifting from them as the old
ones did (#128). This check stops one from coming back unnoticed, as SELAH's scripts/check_separation.py does for the
space operators in GIDEON and SETH. A finding is an assignment that

- multiplies the time step by a stage coefficient: dt * alph(...), beta(s) * dt, ... (alph, beta, gamm, the Butcher
  arrays; ark, brk, crk, the low-storage ones; ssa, ssb, the splitting ones; any of them as a device copy, *_gpu), or
- writes an SSP-RK3 stage or the classic RK4 sum with literal coefficients (0.75 u + 0.25 ..., u/3 + 2/3 ...,
  k1 + 2 k2 + 2 k3 + k4).

The patterns are regular expressions, not a Fortran parser: a statement is the joined continuation lines, comments
cut. A remaining combination is waived by file name in WAIVERS, with the reason. `--self-test` plants each kind in a
scratch file and checks that the rule catches it, and that a plain dt product (a CFL, a flux) does not.

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
RK3 = re.compile(r'0\.75_R\w*\s*\*.*0\.25_R\w*|/\s*3\._R\w*\s*\+\s*2\._R\w*\s*/\s*3\._R\w*', re.IGNORECASE)
RK4 = re.compile(r'\bk1\b\s*\+\s*2\w*\s*\*\s*k2\b', re.IGNORECASE)

# file name -> reason: every combination of the file is waived
WAIVERS: dict[str, str] = {
    'adam_prism_rk_pic_object.F90': "PRISM's particles (PIC) advance in step with the field stepper, with its SELAH "
                                    'tableau: SELAH has no stepper for the particle arrays (#132 phase 6)',
    'adam_prism_fnl_rk_pic_object.F90': 'the device twin of adam_prism_rk_pic_object.F90 (its coefficients, copied)',
    'adam_prism_rk_pml_object.F90': "PRISM's PML auxiliary face fields advance in step with the field stepper, with its "
                                    'SELAH tableau: SELAH has no stepper for the face arrays (#132 phase 4)',
    'adam_prism_fnl_rk_pml_object.F90': 'the device twin of adam_prism_rk_pml_object.F90 (its coefficients, copied)',
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


def combinations(path: Path) -> list[str]:
    if path.name in WAIVERS:
        return []
    out = []
    for n, s in statements(path):
        if '=' not in s or re.match(r'(?i)(do|if|call|where)\b', s):
            continue
        rhs = s.split('=', 1)[1]
        if STAGE.search(rhs) or RK3.search(rhs) or RK4.search(rhs):
            out.append(f'{path}:{n}: Runge-Kutta combination: {s[:110]}')
    return out


def sources(root: Path) -> list[Path]:
    return sorted(p for d in ('src/app', 'src/lib') for p in (root / d).rglob('*.[fF]90'))


def self_test() -> int:
    """Plant each kind of combination in scratch files: each must be caught, a plain dt product must not."""
    planted = {
        'stage.F90': ('      q_rk(v,i,j,k,b,s) = q_rk(v,i,j,k,b,s) + dt * alph(s,ss) * q_rk(v,i,j,k,b,ss)\n', 1),
        'update.F90': ('      q(v,i,j,k,b) = q(v,i,j,k,b) + &\n                     beta_gpu(s) * dt * dq(v,i,j,k,b)\n', 1),
        'lowstorage.F90': ('      q = self%ark(s) * q0 + self%brk(s) * q + self%crk(s) * dt * dq\n', 1),
        'rk3.F90': ('      u2 = 0.75_R8P*u + 0.25_R8P*(u1 + dt*r)\n', 1),
        'rk4.F90': ('      u = u + dt*(k1 + 2*k2 + 2*k3 + k4)/6\n', 1),
        'cfl.F90': ('      dt = cfl * dx / (abs(u) + c)\n      flux = dt * f(v,i) / dx\n', 0),
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
    print(f'check_time_separation: {len(files)} sources, {len(findings)} finding(s), {len(WAIVERS)} file(s) waived')
    return 1 if findings else 0


if __name__ == '__main__':
    sys.exit(main())
