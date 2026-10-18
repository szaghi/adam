#!/usr/bin/env python3
"""Write the inputs of the FLUME passive-scalar verification (issue #92) from the sod-x regression input.

    make_ps.py <out.ini> [--model euler|euler-de|mhd-ideal] [--scalars N] [--left S] [--right S]

The sod-x shock tube (1-D along x, WENO5 characteristic, SSP-33, CFL 0.5, t = 0.2) with `passive_scalars = N` (0: the
key absent) and every scalar set to `--left` in region 1 (x < 0.5) and `--right` in region 2.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "sod-x" / "input.ini"


def case_ini(model: str, scalars: int, left: float, right: float) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    s = re.sub(r"(?m)^physical_model\s*=.*$", f"physical_model = {model}", s)
    if scalars > 0:
        s = re.sub(r"(?m)^(physical_model\s*=.*)$", rf"\1\npassive_scalars = {scalars}", s)
        for region, value in ((1, left), (2, right)):
            keys = "".join(f"s{k} = {value!r}\n" for k in range(1, scalars + 1))
            s = s.replace(f"[initial_conditions_region_{region}]\n", f"[initial_conditions_region_{region}]\n{keys}")
    return re.sub(r"(?m)^(\s*it_save\s*=).*$", r"\1 100000", s)


def main() -> int:
    """Write the requested case."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--model", choices=("euler", "euler-de", "mhd-ideal"), default="euler")
    ap.add_argument("--scalars", type=int, default=1)
    ap.add_argument("--left", type=float, default=1.0)
    ap.add_argument("--right", type=float, default=0.0)
    a = ap.parse_args()
    a.out.write_text(case_ini(a.model, a.scalars, a.left, a.right))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
