#!/usr/bin/env python3
"""Write the inputs of the FLUME dual-energy verification (issue #91, P2) from the sod-x regression input.

Cases:
    sod    the sod-x shock tube as it is (1-D along x, WENO5 characteristic, SSP-33, CFL 0.5, t = 0.2)
    noh    the planar Noh problem in a cold gas: two streams rho = 1, u = +-1, p = 1e-6 (Mach 845, gamma = 1.4)
           colliding at x = 0.5, to t = 0.6: the exact solution has a shock moving out at (gamma - 1) / 2 = 0.2,
           rho = (gamma + 1) / (gamma - 1) = 6 and p = rho0 u0 (u0 + 0.2) = 1.2 behind it, the cold streams ahead

Usage:
    make_de.py <out.ini> <case> [--model euler|euler-de] [--eta ETA]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "sod-x" / "input.ini"


def case_ini(case: str, model: str, eta: float | None) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    s = re.sub(r"(?m)^physical_model\s*=.*$", f"physical_model = {model}", s)
    if eta is not None:
        s = re.sub(r"(?m)^(physical_model\s*=.*)$", rf"\1\ndual_energy_eta = {eta!r}", s)
    if case == "noh":
        s = re.sub(r"(?m)^cp\s*=.*\n^cv\s*=.*\n", "gamma = 1.4\n", s)
        regions = s.split("[initial_conditions_region_")
        for k, sign in ((1, 1.0), (2, -1.0)):
            r = regions[k]
            r = re.sub(r"(?m)^r\s*=.*$", "r      = 1.0", r)
            r = re.sub(r"(?m)^u\s*=.*$", f"u      = {sign!r}", r)
            r = re.sub(r"(?m)^p\s*=.*$", "p      = 1.0e-6", r)
            regions[k] = r
        s = "[initial_conditions_region_".join(regions)
        s = re.sub(r"(?m)^(\s*time_max\s*=).*$", r"\1 0.6", s)
    s = re.sub(r"(?m)^(\s*it_save\s*=).*$", r"\1 100000", s)
    return s


def main() -> int:
    """Write the requested case."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path)
    parser.add_argument("case", choices=("sod", "noh"))
    parser.add_argument("--model", choices=("euler", "euler-de"), default="euler-de")
    parser.add_argument("--eta", type=float, default=None, help="[physics].(dual_energy_eta) (default: absent)")
    args = parser.parse_args()
    args.out.write_text(case_ini(args.case, args.model, args.eta))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
