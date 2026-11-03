#!/usr/bin/env python3
"""The Taylor-Green campaign of issue #119 against the 512^3 spectral DNS (High-Order CFD Workshop, case C3.5, spectral_Re1600_512.gdiag): one row per case
directory. Peak -dK/dt over a window of 0.5, as the page; L2 errors over 0.3 < t < 19.7, relative to the DNS: of -dK/dt (both over the window), of Z and of K.

Usage: tgv_dns.py <gdiag> <directory holding the cases> <case>...   (the DNS file: turbulence/reference/spectral_Re1600_512.gdiag)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

VOLUME = (2.0 * np.pi) ** 3


def main() -> int:
    """Print the table."""
    tdns, kdns, ddns, zdns = np.loadtxt(sys.argv[1], comments="#").T
    print(f"DNS: peak -dK/dt {ddns.max():.5f} at t {tdns[ddns.argmax()]:.2f}, Z {zdns.max():.3f} at t {tdns[zdns.argmax()]:.2f}")
    print(f"{'case':30s} {'-dK/dt':>8s} {'at t':>6s} {'Z':>7s} {'at t':>6s} {'L2 dK':>7s} {'L2 Z':>7s} {'L2 K':>7s}")
    keep = (tdns > 0.3) & (tdns < 19.7)
    t = tdns[keep]
    window = np.ones(50) / 50.0
    dref = np.convolve(ddns, window, "same")[keep]

    def err(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.sqrt(np.trapezoid((a - b) ** 2, t) / np.trapezoid(b**2, t)))

    for case in sys.argv[3:]:
        h = np.loadtxt(Path(sys.argv[2]) / case / "tgv-turbulence_history.dat", skiprows=1, ndmin=2)
        _, last = np.unique(h[::-1, 0], return_index=True)  # the last copy of the rows run twice
        h = h[::-1][last]
        k, z = h[:, 2] / VOLUME, h[:, 3] / VOLUME
        dk = -(np.interp(t + 0.25, h[:, 1], k) - np.interp(t - 0.25, h[:, 1], k)) / 0.5
        print(f"{case:30s} {dk.max():8.5f} {t[dk.argmax()]:6.2f} {z.max():7.3f} {h[z.argmax(), 1]:6.2f} {err(dk, dref):7.3f} "
              f"{err(np.interp(t, h[:, 1], z), zdns[keep]):7.3f} {err(np.interp(t, h[:, 1], k), kdns[keep]):7.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
