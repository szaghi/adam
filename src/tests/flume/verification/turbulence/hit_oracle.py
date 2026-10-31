#!/usr/bin/env python3
"""Check the isotropic turbulence initial field (issue #119, T0.5) of each run (argument n:work), coarse to fine.

Every shell of the spectrum against E(s) ~ (s/k_0)^4 exp(-2 (s/k_0)^2), s <= 3 k_0, normalised to a mean K of 3/2 (u_rms = 1,
rho = 1); the K of the first history row; the dilatational dissipation over the solenoidal one, which vanishes only as the
gradients converge, so it must fall with n: at these n the top shells are not in the asymptotic range, hence a factor, not
an order.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

VOLUME = (2.0 * np.pi) ** 3


def main() -> int:
    """Check the ladder."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--k0", type=float, default=4.0)
    parser.add_argument("--shell-tol", type=float, default=1e-10, help="relative to the largest shell")
    parser.add_argument("--kin-tol", type=float, default=1e-12)
    parser.add_argument("--fall-min", type=float, default=4.0)
    args = parser.parse_args()
    spectra = Path(__file__).resolve().parent / "spectra.py"
    ok = True
    ratios = []
    for arg in args.runs:
        n, work = arg.split(":")
        if subprocess.run([sys.executable, str(spectra), work, "--step", "0"], check=False).returncode:
            ok = False
        got = np.loadtxt(Path(work) / "spectrum-000000000.dat")[:, 1]
        shells = np.arange(len(got))
        x = shells / args.k0
        want = np.where(shells <= int(3 * args.k0), x**4 * np.exp(-2.0 * x**2), 0.0)
        want *= 1.5 / want.sum()
        err = float(np.max(np.abs(got - want)) / want.max())
        row = np.loadtxt(Path(work) / "tgv-turbulence_history.dat", skiprows=1, ndmin=2)[0]
        kin = abs(row[2] / VOLUME - 1.5) / 1.5
        ratios.append(row[5] / row[4])
        good = err < args.shell_tol and kin < args.kin_tol
        ok = ok and good
        print(f"n {int(n):4d}: largest shell error {err:.2e} (bound {args.shell_tol}), K rel {kin:.2e}, eps_d/eps_s {ratios[-1]:.3e}  "
              f"{'PASS' if good else 'FAIL'}")
    falls = [ratios[i] / ratios[i + 1] for i in range(len(ratios) - 1)]
    good = bool(falls) and min(falls) >= args.fall_min
    ok = ok and good
    print(f"eps_d/eps_s falls by {' '.join(f'{f:.1f}' for f in falls)} per doubling (min {args.fall_min})  {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
