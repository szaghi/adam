#!/usr/bin/env python3
"""Write an Euler near-vacuum input for the positivity verification PV-3 (issue #47, M3-P5c).

Why: the three 1-D Euler problems that the positivity literature uses to drive density or pressure to (near) zero; a
scheme without positivity control produces a negative pressure in them. Grid and numerics of verification/sod/sod-x.ini
(4 blocks along x, y and z null, WENO-5 characteristic, SSP-33, CFL 0.5), extrapolation boundaries:
  * leblanc: gamma = 5/3, [0, 9], (rho, u, e) = (1, 0, 0.1) for x < 3, (1e-3, 0, 1e-7) beyond, t = 6, 400 cells
    (pressure jump 1e9, the internal energy one 1e6);
  * double-rarefaction: gamma = 1.4, [-1, 1], (rho, u, p) = (7, -1, 0.2) for x < 0, (7, 1, 0.2) beyond, t = 0.6,
    h = 0.005 (400 cells): the two rarefactions leave a near-vacuum at the centre;
  * sedov: gamma = 1.4, the planar Sedov blast on [0, 4] with h = 0.005 (800 cells; the domain is shifted by h/2 so
    that a cell is centred at x = 2), (rho, u, p) = (1, 0, 4e-13) except p = 2.56e8 in the centre cell, t = 0.001.
Setups as Zhang & Shu (2010, doi:10.1016/j.jcp.2010.08.016) and Guo et al. (arXiv:1402.5618, Examples 8-10).
Usage: make_vacuum.py <sod-x.ini> <out.ini> --problem leblanc|double-rarefaction|sedov
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

H = 0.005  # cell size of double-rarefaction and sedov
PROBLEMS = {  # cells per block, domain, final time, (cp, cv), regions as (emin_x, emax_x, r, u, p), first match wins
    "leblanc": (100, (0.0, 9.0), "6.0", ("2.5", "1.5"),
                [(-1e30, 3.0, 1.0, 0.0, 2.0 / 3.0 * 0.1), (3.0, 1e30, 1.0e-3, 0.0, 2.0 / 3.0 * 1.0e-3 * 1.0e-7)]),
    "double-rarefaction": (100, (-1.0, 1.0), "0.6", ("1.4", "1.0"),
                           [(-1e30, 0.0, 7.0, -1.0, 0.2), (0.0, 1e30, 7.0, 1.0, 0.2)]),
    "sedov": (200, (-0.5 * H, 4.0 - 0.5 * H), "0.001", ("1.4", "1.0"),
              [(2.0 - 0.5 * H, 2.0 + 0.5 * H, 1.0, 0.0, 2.56e8), (-1e30, 1e30, 1.0, 0.0, 4.0e-13)]),
}


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("base", type=Path, help="verification/sod/sod-x.ini")
    parser.add_argument("out", type=Path)
    parser.add_argument("--problem", choices=tuple(PROBLEMS), required=True)
    args = parser.parse_args()
    ni, (x0, x1), t_max, (cp, cv), regions = PROBLEMS[args.problem]
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(args.base)
    ini["grid"].update({"ni": str(ni), "emin_x": repr(x0), "emax_x": repr(x1)})
    ini["physics"].update({"cp": cp, "cv": cv})
    for sec in [k for k in ini.sections() if k.startswith("initial_conditions_region_")]:
        ini.remove_section(sec)
    ini["initial_conditions"]["regions_number"] = str(len(regions))
    for i, (e0, e1, r, u, p) in enumerate(regions, start=1):
        ini[f"initial_conditions_region_{i}"] = {
            "r": repr(r), "u": repr(u), "v": "0.0", "w": "0.0", "p": repr(p), "emin_x": repr(e0), "emax_x": repr(e1),
            "emin_y": "-1e30", "emax_y": "1e30", "emin_z": "-1e30", "emax_z": "1e30"}
    ini["time"].update({"it_max": "-1", "time_max": t_max})
    ini["IO"].update({"output_basename": args.problem, "restart_basename": f"{args.problem}-restart"})
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
