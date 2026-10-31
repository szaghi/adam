#!/usr/bin/env python3
"""The table of the Taylor-Green campaign (issue #119, T0.3, T0.4): one row per case directory of tgv_campaign.sh.

All columns are per unit volume, as the literature gives them. `-dK/dt` is the total dissipation, `eps - pd` the resolved
part; the difference is the scheme's numerical dissipation, reported as its share of the K lost until the peak of -dK/dt.
Rows run twice by a restart are kept once, the last copy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

VOLUME = (2.0 * np.pi) ** 3


def history(work: Path) -> np.ndarray:
    """Return the history rows, unique in the step column, the last copy of each, per unit volume."""
    rows = np.loadtxt(work / "tgv-turbulence_history.dat", comments="VARIABLES", ndmin=2)
    _, last = np.unique(rows[::-1, 0], return_index=True)
    rows = rows[::-1][last]
    rows[:, 2:] /= VOLUME
    return rows


def main() -> None:
    """Print the table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    print("| case | t reached | peak -dK/dt | at t | peak eps - pd | at t | peak Z | at t | numerical share to the peak | K(t end) |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for work in sorted(p for p in args.root.iterdir() if (p / "tgv-turbulence_history.dat").exists()):
        h = history(work)
        t, kin, zeta, eps, pd = h[:, 1], h[:, 2], h[:, 3], h[:, 6], h[:, 7]
        rate = -np.gradient(kin, t)
        i, j, m = np.argmax(rate), np.argmax(eps - pd), np.argmax(zeta)
        phys = np.trapezoid(eps[: i + 1] - pd[: i + 1], t[: i + 1])
        lost = kin[0] - kin[i]
        print(f"| {work.name} | {t[-1]:.2f} | {rate[i]:.5f} | {t[i]:.2f} | {(eps - pd)[j]:.5f} | {t[j]:.2f} | {zeta[m]:.3f} | "
              f"{t[m]:.2f} | {(lost - phys) / lost:.3f} | {kin[-1]:.5f} |")


if __name__ == "__main__":
    main()
