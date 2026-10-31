#!/usr/bin/env python3
"""Write the input of the compressible Taylor-Green vortex (issue #119, T0).

The box is [0, 2 pi]^3, periodic, with rho_0 = V_0 = L = 1, gamma = 1.4 and R = 1, Re = 1 / mu, on 2^(3 levels) blocks of
(n / 2^levels)^3 cells. The other sections come from the vortex-periodic regression input. Standard library only.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parent.parent.parent / "regression" / "vortex-periodic" / "input.ini"


def section(text: str, name: str, entries: dict[str, str]) -> str:
    """Return `text` with the keys of a section set, a key that is absent added at its top."""
    match = re.search(rf"^\[{re.escape(name)}\][^\[]*", text, flags=re.M)
    if match is None:
        raise SystemExit(f"make_tgv: no section [{name}] in {TEMPLATE}")
    body = match.group(0)
    for key, val in entries.items():
        line = f"{key} = {val}"
        body, count = re.subn(rf"^{re.escape(key)}\s*=.*$", line, body, flags=re.M)
        if not count:
            body = body.replace(f"[{name}]\n", f"[{name}]\n{line}\n", 1)
    return text[: match.start()] + body + text[match.end() :]


def main() -> None:
    """Parse the arguments and write the input."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--n", type=int, default=64, help="cells per direction")
    parser.add_argument("--levels", type=int, default=2, help="refinement levels: 2^levels blocks per direction")
    parser.add_argument("--mach", type=float, default=0.1)
    parser.add_argument("--reynolds", type=float, default=1600.0, help="0: inviscid")
    parser.add_argument("--scheme", default="weno-u-5")
    parser.add_argument("--scheme-space", default="weno")
    parser.add_argument("--riemann", default="hllc")
    parser.add_argument("--cfl", type=float, default=0.5)
    parser.add_argument("--time-max", type=float, default=20.0)
    parser.add_argument("--it-max", type=int, default=-1)
    parser.add_argument("--it-save", type=int, default=100000)
    parser.add_argument("--history", type=int, default=1, help="cadence of the conservation and turbulence histories")
    args = parser.parse_args()

    blocks = 2**args.levels
    if args.n % blocks:
        raise SystemExit(f"make_tgv: {args.n} cells do not split into {blocks} blocks")
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)  # the template's own header
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    box = 2.0 * math.pi
    text = section(text, "IO", {"output_basename": "tgv", "it_save": str(args.it_save)})
    grid = {"ni": str(args.n // blocks), "nj": str(args.n // blocks), "nk": str(args.n // blocks),
            "ngc": str(max(3, (int(args.scheme.rsplit("-", 1)[1]) + 1) // 2))}
    grid |= {f"emin_{a}": "0.0" for a in "xyz"} | {f"emax_{a}": f"{box:.16e}" for a in "xyz"}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": str(args.levels), "ratio": "8", "iu_ref_levels": str(args.levels),
                                 "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": args.scheme, "weights": "si"})
    physics = {"cp": "3.5", "cv": "2.5"}
    if args.reynolds > 0.0:
        physics |= {"reynolds": str(args.reynolds), "prandtl": "0.71"}
    text = section(text, "physics", physics)
    numerics = {"scheme_space": args.scheme_space}
    if args.scheme_space == "weno-riemann":
        numerics["riemann_solver"] = args.riemann
    text = section(text, "numerics", numerics)
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    ic = {"type": "taylor-green", "amr_iterations": "0", "tg_density": "1.0", "tg_velocity": "1.0",
          "tg_mach": str(args.mach)}
    text = section(text, "initial_conditions", ic)
    text = section(text, "time", {"it_max": str(args.it_max), "time_max": str(args.time_max), "CFL": str(args.cfl)})
    text = section(text, "diagnostics", {"conservation_history_save": str(args.history),
                                         "turbulence_history_save": str(args.history)})
    args.out.write_text(text)


if __name__ == "__main__":
    main()
