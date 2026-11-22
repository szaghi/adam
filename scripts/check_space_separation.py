#!/usr/bin/env python3
"""Check that no WENO, TENO or hybrid operator is written in ADAM's sources: the space operators are SETH's (#139).

Why this script exists: #139 makes SETH the only source of WENO, TENO and hybrid operators in ADAM. ADAM held one WENO
algorithm in four copies plus two more hand-written ones (finding F1 of #139), and the copies drifted (the host and the
device ran different exponents until #95, F6). This check stops a coefficient table, a smoothness indicator, a
nonlinear weight, a central split flux or a deleted kernel name from coming back unnoticed, as check_time_separation.py
does for the time integrators (#132). A finding is a statement that

- writes a smoothness-indicator constant: 13/12, 781/720, or two or more of the Jiang-Shu WENO-5 rows' thirds
  (10/3, 31/3, 25/3, 19/3, 11/3, 13/3) in one statement,
- forms a nonlinear weight: a division by (eps + IS)**p, a regulariser named *eps* or written as a literal 1e-N
  (Jiang-Shu), or a WENO-Z ratio (tau / (IS + eps))**p,
- writes the KEEP two-point mass flux, a quarter of a product of two sums (0.25 (a_l + a_r)(b_l + b_r)), or
- names a kernel or table that SETH replaced: weno_a, weno_p, weno_d, weno_reconstruct_upwind*, compute_kep_flux.

The patterns are regular expressions, not a Fortran parser: a statement is the joined continuation lines, comments cut;
the sources are the .F90 and .INC files of src/app and src/lib (the archived apps excluded). A remaining file is waived
in WAIVERS by name, with the step of #139 that migrates it: the waiver list is the progress meter of the migration, and
nothing new can be written while it runs (D9). `--self-test` plants each kind in a scratch file and checks that the rule
catches it, and that a plain stencil (a finite difference, a SETH call) does not.

Usage: python3 scripts/check_space_separation.py [--root DIR] [--self-test]     (exit status 1 on a finding)
"""

from __future__ import annotations

import argparse
import fnmatch
import re
import sys
import tempfile
from pathlib import Path

R = r'(?:\._R\w*|\.\d*_R\w*|\.d0|\.0?)?'  # the kind suffix of a literal: 13._R8P, 13.0, 13
IS_CONST = re.compile(r'\b13' + R + r'\s*/\s*12' + R + r'(?!\d)|\b781' + R + r'\s*/\s*720' + R + r'(?!\d)', re.IGNORECASE)
THIRDS = [re.compile(r'\b' + n + R + r'\s*/\s*3' + R + r'(?![\d.])', re.IGNORECASE) for n in ('10', '31', '25', '19', '11', '13')]
EPS = r'(?:[\w%]*eps\w*|1\.?\d*e-\d+\w*)'
JS_WEIGHT = re.compile(r'/\s*\(\s*' + EPS + r'\s*\+[^()]*(?:\([^()]*\)[^()]*)*\)\s*\*\*', re.IGNORECASE)
Z_WEIGHT = re.compile(r'\(\s*\w+\s*/\s*\([^=]*' + EPS + r'\s*\)\s*\)\s*\*\*', re.IGNORECASE)
SUM = r'\((?:[^()]|\([^()]*\))+\+(?:[^()]|\([^()]*\))+\)'  # a parenthesised sum, one level of inner parentheses
KEEP = re.compile(r'0\.25_R\w*\s*\*\s*' + SUM + r'\s*\*\s*' + SUM, re.IGNORECASE)
NAMES = re.compile(r'\b(?:weno_[apd](?:_gpu)?\b|weno_reconstruct_upwind\w*|compute_kep_flux)\b', re.IGNORECASE)

# file name -> (the step of #139 that migrates it, what it still holds)
_S6 = "S6: ADAM's WENO deleted"
_S7 = 'S7: the seam WENO-6 transfer onto SETH (E2)'
WAIVERS: dict[str, str] = {  # a file name, or a glob of names
    'adam_weno_object.F90': _S6 + ' (tables, host kernels)',
    'adam_fnl_weno_kernels.F90': _S6 + ' (the device copy)',
    'adam_fnl_weno_width_kernels.F90': _S6 + ' (the per-width device copy)',
    'adam_fnl_weno_width_body.INC': _S6 + ' (the per-width device body)',
    'adam_weno_nvf_kernels.F90': _S6 + ' (the CUDA Fortran twin, NASTO only)',
    'adam_flume_hybrid_library.F90': _S6 + ' (compute_kep_flux, compute_ducros: test_kep and the parity test still call them)',
    'adam_seam_interpolation_library.F90': _S7,
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


def waived(path: Path) -> bool:
    return any(fnmatch.fnmatchcase(path.name, w) for w in WAIVERS)


def operators(path: Path) -> list[str]:
    out = []
    for n, s in statements(path):
        if re.match(r'(?i)#', s):  # preprocessor lines
            continue
        kind = ('deleted name' if NAMES.search(s) else
                'smoothness-indicator constant' if IS_CONST.search(s) or sum(bool(t.search(s)) for t in THIRDS) >= 2 else
                'nonlinear weight' if JS_WEIGHT.search(s) or Z_WEIGHT.search(s) else
                'KEEP split flux' if KEEP.search(s) else '')
        if kind:
            out.append(f'{path}:{n}: {kind}: {s[:110]}')
    return out


ARCHIVED = ('nasto', 'chase', 'patch', 'ascot')  # archived apps, not maintained


def sources(root: Path) -> list[Path]:
    return sorted(p for d in ('src/app', 'src/lib') for p in (root / d).rglob('*')
                  if p.suffix in ('.F90', '.f90', '.INC') and not set(ARCHIVED) & set(p.parts))


def self_test() -> int:
    """Plant each kind in scratch files: each must be caught; a plain stencil and a SETH call must not."""
    planted = {
        'is13.F90': ('      is0 = 13._R8P/12._R8P*(v(1) - 2._R8P*v(2) + v(3))**2 + 0.25_R8P*(v(1) - v(3))**2\n', 1),
        'is781.F90': ('      b = c*c + 781._R8P / 720._R8P * t * t\n', 1),
        'thirds.F90': ('      is0 = 10._R8P/3._R8P * v(c)**2 - 31._R8P/3._R8P * v(c-1) * v(c) + &\n'
                       '            25._R8P/3._R8P * v(c-1)**2\n', 1),
        'jsweight.F90': ('      a0 = d0 / (1.e-6_R8P + is0)**3\n', 1),
        'jsweight2.F90': ('      a(f) = weno_aa(f)*(1._R8P/(zeps+IS(f)*rmu*rmu)**wexp)\n', 1),
        'zweight.F90': ('      alpha = dlin(s) * (1._R8P + (tau / (c * c + SEAM_EPS))**2)\n', 1),
        'keep.F90': ('      mf = 0.25_R8P * (w(1,c) + w(1,r)) * (w(1+d,c) + w(1+d,r))\n', 1),
        'name.F90': ('      call weno_reconstruct_upwind_s3_dev(a, p, d, zeps, rmu, v, vr)\n', 1),
        'table.F90': ('      x = weno_a(1,0,3) + weno_p(1,0,0,3)\n', 1),
        'plain.F90': ('      du = (u(i+1) - u(i-1)) / (2._R8P*dx)\n'
                      '      d2 = (4._R8P/3._R8P)*(u(i+1) - 2*u(i) + u(i-1)) - (1._R8P/12._R8P)*u(i+2)\n'
                      '      call seth_weno5_js(vp=vp, vm=vm, eps=eps, rmu=rmu, rp=rp, rm=rm)\n'
                      '      f = 0.5_R8P * (fl + fr) - 0.5_R8P * a * (qr - ql)\n'
                      '      x = 1._R8P / (eps + y)\n', 0),
    }
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for name, (src, expect) in planted.items():
            p = Path(d) / name
            p.write_text('program p\n' + src + 'endprogram p\n')
            got = len(operators(p))
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
    findings = [f for p in files if not waived(p) for f in operators(p)]
    stale = sorted(w for w in WAIVERS if not any(fnmatch.fnmatchcase(p.name, w) and operators(p) for p in files))
    for f in findings:
        print(f)
    for w in stale:
        print(f'stale waiver: {w} holds no finding: delete it from WAIVERS ({WAIVERS[w]})')
    print(f'check_space_separation: {len(files)} sources, {len(findings)} finding(s), {len(WAIVERS)} waived file(s), '
          f'{len(stale)} stale waiver(s)')
    return 1 if findings or stale else 0


if __name__ == '__main__':
    sys.exit(main())
