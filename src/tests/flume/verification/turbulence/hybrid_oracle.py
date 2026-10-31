#!/usr/bin/env python3
"""Check the inviscid Taylor-Green run of the hybrid scheme against the same run on WENO (issue #119, T2): what it dissipates, K less the
integral of the pressure work, at most a tenth of WENO's (--factor), and mass and total energy conserved. Arguments: the two work directories.
With --contact, one directory: an advected contact.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "multirealm"))
from multirealm_oracle import basenames, last_step, load_cells  # noqa: E402


def lost(work: str) -> float:
    """Return what the scheme dissipated: the K change less the pressure work, over the initial K. The bare K change is no measure: at low
    Mach the pressure work is of the size of the hybrid's dissipation, of either sign."""
    h = np.loadtxt(Path(work) / "tgv-turbulence_history.dat", skiprows=1, ndmin=2)
    return float(abs(h[-1, 2] - h[0, 2] - np.trapezoid(h[:, 7], h[:, 1])) / h[0, 2])


def drift(work: str, tol: float) -> tuple[bool, str]:
    """Return whether the mass and the energy are conserved within the bound, with the text."""
    path = Path(work) / "tgv-conservation_history.dat"
    with open(path) as fh:
        names = [w.strip('"') for w in fh.readline().split("=")[1].split()]
    c = np.loadtxt(path, skiprows=1, ndmin=2)
    vals = [float(np.max(np.abs(c[:, names.index(v)] - c[0, names.index(v)])) / abs(c[0, names.index(v)])) for v in ("int_r", "int_rE")]
    return max(vals) < tol, f"mass {vals[0]:.2e}, energy {vals[1]:.2e}"


def contact(work: str, tol: float) -> int:
    """Check the velocity and the pressure of the last checkpoint of an advected contact."""
    path = Path(work)
    names = basenames(path)
    variables, cells = load_cells(path, names, last_step(path, names), 3)
    q = np.array(list(cells.values()))
    r, u, v, e = (q[:, variables.index(n)] for n in ("r", "ru", "rv", "rE"))
    u, v = u / r, v / r
    p = 0.4 * (e - 0.5 * r * (u**2 + v**2))
    dev = max(float(np.max(np.abs(u - 1.0))), float(np.max(np.abs(p - 1.0))))
    ok = dev < tol and float(r.min()) > 0.0
    print(f"{len(r)} cells, rho in [{r.min():.4f}, {r.max():.4f}], largest deviation of u and p {dev:.2e} (bound {tol})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Check."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hybrid")
    parser.add_argument("weno", nargs="?")
    parser.add_argument("--contact", action="store_true", help="the one directory is an advected contact")
    parser.add_argument("--contact-tol", type=float, default=1e-3)
    parser.add_argument("--factor", type=float, default=10.0)
    parser.add_argument("--tol", type=float, default=1e-12)
    args = parser.parse_args()
    if args.contact:
        return contact(args.hybrid, args.contact_tol)
    lh, lw = lost(args.hybrid), lost(args.weno)
    good, text = drift(args.hybrid, args.tol)
    ok = good and lh * args.factor <= lw
    print(f"K dissipated, hybrid {lh:.3e}, WENO {lw:.3e}, ratio {lw / lh:.1f} (min {args.factor}); drift {text}  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
