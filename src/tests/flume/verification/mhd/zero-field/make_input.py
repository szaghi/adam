#!/usr/bin/env python3
"""Derive an MHD input from a committed Euler input (issue #41, M2-P3, MV-2).

Why: MV-2 runs committed Euler cases through the MHD path, so the MHD inputs are generated from the Euler ones
(one source of truth): physical_model becomes mhd-ideal, an [mhd] section is added (divergence control none or glm),
every region and inflow state gains bx, by, bz. Options override single keys (section.key=value).

Usage:
    make_input.py <euler.ini> <out.ini> --divergence-control none|glm [--b BX BY BZ] [--set section.key=value ...]
"""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path


def main() -> None:
    """Write the MHD input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("euler", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--divergence-control", choices=("none", "glm"), required=True)
    parser.add_argument("--b", type=float, nargs=3, default=(0.0, 0.0, 0.0), metavar=("BX", "BY", "BZ"))
    parser.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    args = parser.parse_args()
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str  # FLUME keys are case-sensitive (e.g. [time] CFL)
    ini.read(args.euler)
    ini["physics"]["physical_model"] = "mhd-ideal"
    ini["mhd"] = {"divergence_control": args.divergence_control, "divb_tol": "0.0", "divb_error": ".false.",
                  "rho_floor": "0.0", "p_floor": "0.0"}
    if args.divergence_control == "glm":
        ini["mhd"].update({"glm_ch": "3.0", "glm_alpha": "0.18", "glm_damping_length": "1.0",
                           "glm_ch_check": "warning"})
    ini["field"]["nv"] = "9" if args.divergence_control == "glm" else "8"
    for name in ini.sections():
        if name.startswith("initial_conditions_region_") or (name.startswith("bc_") and ini[name]["type"] == "inflow"):
            for key, val in zip(("bx", "by", "bz"), args.b, strict=True):
                ini[name][key] = repr(val)
    for item in args.set:
        key, val = item.split("=", 1)
        section, option = key.split(".", 1)
        ini[section][option] = val
    with args.out.open("w") as out:
        ini.write(out)


if __name__ == "__main__":
    main()
