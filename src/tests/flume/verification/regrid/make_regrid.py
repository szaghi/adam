#!/usr/bin/env python3
"""Write the inputs of the FLUME runtime regridding checks (issue #74, RG legs rg3-rg5).

Why this generator exists: the runtime regrid legs need cases where the grid really changes during the run, built
from inputs the suite already verifies, so that only the regrid is new:

- ``vortex``: the isentropic vortex in a periodic box (``regression/vortex-periodic``) on a quadtree, base level 1,
  ``max_level`` 3, 8x8 cells per block, a Loehner marker on the density, a regrid every 5 steps, reflux on. The vortex
  travels across blocks, so regrids refine ahead of it and derefine behind it; nothing crosses the boundary, so every
  volume integral must stay constant (rg3) and a restart taken at a regrid step must continue bitwise (rg4).
- ``cylinder``: the shock over the cylinder (``verification/shock-cylinder``) with the solid marker and a Loehner
  marker on the density, a regrid every 5 steps (rg5).

Usage:
    make_regrid.py vortex <out.ini> [--prolongation linear|conservative] [--it-max N] [--restart-save N] [--restart]
    make_regrid.py cylinder <out.ini> [--it-max N]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
VORTEX = REPO / "src/tests/flume/regression/vortex-periodic/input.ini"
CYLINDER = REPO / "src/tests/flume/verification/shock-cylinder/shock-cylinder.ini"

LOHNER = """[amr_marker_{n}]
mode         = 4             ; Loehner estimator (issue #74)
delta_type   = max           ; unused by the Loehner marker; required by the common marker parse
delta_fine   = 0.0
delta_coarse = 0.0
field        = 1             ; conservative variables
var          = 1             ; density
refine_tol   = {refine}
derefine_tol = {derefine}
buffer       = 2

"""


def setkey(text: str, key: str, value: str) -> str:
    """Return ``text`` with every ``key = ...`` line set to ``value`` (the key must exist).

    Parameters
    ----------
    text : str
        INI text.
    key : str
        Option name.
    value : str
        New value.

    Returns
    -------
    str
        The edited text.

    Raises
    ------
    KeyError
        If the key is not in the text.
    """
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=).*$")
    if not pattern.search(text):
        raise KeyError(key)
    return pattern.sub(lambda m: f"{m.group(1)} {value}", text)


def vortex(args: argparse.Namespace) -> str:
    """Return the input of the travelling vortex case."""
    text = VORTEX.read_text()
    for key, value in (
        ("ni", "8"),
        ("nj", "8"),
        ("iu_ref_levels", "1"),
        ("max_level", "3"),
        ("iters", "2"),
        ("frequency", f"5\nregrid_prolongation = {args.prolongation}"),
        ("markers_number", "1"),
        ("output_basename", "vortex-rg"),
        ("restart_basename", "vortex-rg-restart"),
        ("restart_save", str(args.restart_save)),
        ("restart", ".true." if args.restart else ".false."),
        ("it_max", str(args.it_max)),
        ("time_max", "10.0"),
        ("CFL", "0.4"),
    ):
        text = setkey(text, key, value)
    return text.replace("[field]", LOHNER.format(n=1, refine="0.3", derefine="0.1") + "[field]", 1)


def cylinder(args: argparse.Namespace) -> str:
    """Return the input of the shock over the cylinder with runtime regridding."""
    text = CYLINDER.read_text()
    for key, value in (
        ("frequency", "5"),
        ("markers_number", "2"),
        ("it_max", str(args.it_max)),
    ):
        text = setkey(text, key, value)
    return text.replace("[field]", LOHNER.format(n=2, refine="0.6", derefine="0.2") + "[field]", 1)


def main() -> None:
    """Parse the command line and write the input."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", choices=("vortex", "cylinder"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--prolongation", choices=("linear", "conservative"), default="conservative")
    parser.add_argument("--it-max", type=int, default=200)
    parser.add_argument("--restart-save", type=int, default=0)
    parser.add_argument("--restart", action="store_true")
    args = parser.parse_args()
    text = vortex(args) if args.case == "vortex" else cylinder(args)
    args.out.write_text(text)


if __name__ == "__main__":
    main()
