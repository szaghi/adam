#!/usr/bin/env python3
"""Write an input of the shock-generator verification (issue #121, S3), from make_tbl.py's template.

A box [0, 2] x [0, 1] x [0, 0.125] of 128 x 64 x 8 cells (8 blocks, periodic z), an inviscid free stream rho = 1, p = 1 / gamma (c = 1)
at Mach `--mach` along x, entering through x_min (`inflow`) and leaving through x_max and y_min (`extrapolation`). The top face y_max is
the `shock-generator`: the oblique shock of deflection `--deflection` (or angle `--angle`) enters at `--shock-x`; `--top extrapolation`
replaces it for the bit-for-bit leg. `--bottom wall-inviscid` makes y_min a slip wall, where the shock reflects (S4a).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "turbulence"))
from make_tbl import TEMPLATE, section  # noqa: E402

GAMMA = 1.4


def main() -> int:
    """Write the input."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--mach", type=float, default=2.0)
    ap.add_argument("--deflection", type=float, help="shock_deflection, degrees (default 10 when --angle is absent)")
    ap.add_argument("--angle", type=float, help="shock_angle, degrees")
    ap.add_argument("--shock-x", type=float, default=0.5)
    ap.add_argument("--shock-x-end", type=float)
    ap.add_argument("--top", default="shock-generator", choices=["shock-generator", "extrapolation"])
    ap.add_argument("--bottom", default="extrapolation", choices=["extrapolation", "wall-inviscid"])
    ap.add_argument("--time-max", type=float, default=4.0)
    ap.add_argument("--nx", type=int, default=128)
    args = ap.parse_args()
    p_inf = 1.0 / GAMMA
    nx, ny, nz = args.nx, args.nx // 2, 8
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    text = section(text, "IO", {"output_basename": "sg", "it_save": "100000", "restart_basename": "sg-restart", "restart_save": "0",
                                "restart": ".false."})
    grid = {"ni": str(nx // 2), "nj": str(ny // 2), "nk": str(nz // 2), "ngc": "3", "emin_x": "0.0", "emax_x": "2.0", "emin_y": "0.0",
            "emax_y": "1.0", "emin_z": "0.0", "emax_z": repr(nz * 2.0 / nx)}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": "1", "ratio": "8", "iu_ref_levels": "1", "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": "weno-u-5", "weights": "si"})
    text = re.sub(r"^\[physics\]\n[^\[]*", "", text, flags=re.M)
    text = section(text, "physics", {"physical_model": "euler", "gamma": repr(GAMMA)})
    text = section(text, "numerics", {"scheme_space": "weno"})
    free = {"r": "1.0", "u": repr(args.mach), "v": "0.0", "w": "0.0", "p": repr(p_inf)}
    text = section(text, "bc_x_min", {"type": "inflow"} | free)
    text = section(text, "bc_x_max", {"type": "extrapolation"})
    text = section(text, "bc_y_min", {"type": args.bottom})
    top = {"type": args.top}
    if args.top == "shock-generator":
        top |= free | {"shock_x": repr(args.shock_x)}
        top |= {"shock_angle": repr(args.angle)} if args.angle is not None else {"shock_deflection": repr(args.deflection or 10.0)}
        if args.shock_x_end is not None:
            top |= {"shock_x_end": repr(args.shock_x_end)}
    text = section(text, "bc_y_max", top)
    for face in ("z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    text = re.sub(r"^\[initial_conditions\]\n[^\[]*", "", text, flags=re.M)
    text = section(text, "initial_conditions", {"type": "uniform", "amr_iterations": "0", "s": "0.0"})
    text = section(text, "initial_conditions_region_1", free)
    text = section(text, "time", {"it_max": "-1", "time_max": repr(args.time_max), "CFL": "0.5"})
    text = section(text, "diagnostics", {"conservation_history_save": "10", "turbulence_history_save": "0"})
    args.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
