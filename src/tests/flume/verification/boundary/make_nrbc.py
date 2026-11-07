#!/usr/bin/env python3
"""Write an input of the non-reflecting boundary verification (issue #121, S2), from make_tbl.py's template.

A box [0, 2] x [0, 0.0625]^2 of 256 x 8 x 8 cells (8 blocks, periodic y and z), a free stream rho = 1, p = 1 / gamma (c = 1) at speed u
along x, and either a right-running acoustic wave packet (the `sine-wave` IC with a Gaussian envelope) or a uniform overpressure. The x
faces take `--xmin` and `--xmax`: `non-reflecting` (target the free stream), `extrapolation` or `inflow` (the free stream).
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
    ap.add_argument("--u", type=float, default=0.0, help="free-stream speed along x (c = 1)")
    ap.add_argument("--xmin", default="non-reflecting", choices=["non-reflecting", "extrapolation", "inflow"])
    ap.add_argument("--xmax", default="non-reflecting", choices=["non-reflecting", "extrapolation", "inflow"])
    ap.add_argument("--ic", default="packet", choices=["packet", "overpressure"])
    ap.add_argument("--amplitude", type=float, default=1.0e-3, help="packet: relative amplitude; overpressure: p / p_target - 1")
    ap.add_argument("--sigma", type=float, default=0.28)
    ap.add_argument("--left", action="store_true", help="packet: left-running (toward x_min), else right-running")
    ap.add_argument("--time-max", type=float, default=1.5)
    ap.add_argument("--it-save", type=int, default=100000)
    ap.add_argument("--nx", type=int, default=256)
    args = ap.parse_args()
    p_inf = 1.0 / GAMMA
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    text = section(text, "IO", {"output_basename": "nrbc", "it_save": str(args.it_save), "restart_basename": "nrbc-restart",
                                "restart_save": "0", "restart": ".false."})
    grid = {"ni": str(args.nx // 2), "nj": "4", "nk": "4", "ngc": "3", "emin_x": "0.0", "emax_x": "2.0", "emin_y": "0.0",
            "emax_y": "0.0625", "emin_z": "0.0", "emax_z": "0.0625"}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": "1", "ratio": "8", "iu_ref_levels": "1", "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": "weno-u-5", "weights": "si"})
    text = re.sub(r"^\[physics\]\n[^\[]*", "", text, flags=re.M)
    text = section(text, "physics", {"physical_model": "euler", "gamma": repr(GAMMA)})
    text = section(text, "numerics", {"scheme_space": "weno"})
    free = {"r": "1.0", "u": repr(args.u), "v": "0.0", "w": "0.0", "p": repr(p_inf)}
    for face, kind in (("x_min", args.xmin), ("x_max", args.xmax)):
        keys = {"type": kind}
        if kind == "non-reflecting":
            keys |= free | {"nr_sigma": repr(args.sigma), "nr_length": "2.0"}
        elif kind == "inflow":
            keys |= free
        text = section(text, f"bc_{face}", keys)
    for face in ("y_min", "y_max", "z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    if args.ic == "packet":
        ic = {"type": "sine-wave", "amr_iterations": "0", "wave_mode": "acoustic", "wave_angle": "180.0" if args.left else "0.0",
              "wave_amplitude": repr(args.amplitude), "wavelength": "0.2", "wave_envelope_width": "0.1",
              "wave_envelope_center": "-1.0" if args.left else "1.0"}  # along n
        p0 = p_inf
    else:
        ic = {"type": "uniform", "amr_iterations": "0", "s": "0.0"}
        p0 = p_inf * (1.0 + args.amplitude)
    text = re.sub(r"^\[initial_conditions\]\n[^\[]*", "", text, flags=re.M)
    text = section(text, "initial_conditions", ic)
    text = section(text, "initial_conditions_region_1", {"r": "1.0", "u": repr(args.u), "v": "0.0", "w": "0.0", "p": repr(p0)})
    text = section(text, "time", {"it_max": "-1", "time_max": repr(args.time_max), "CFL": "0.5"})
    text = section(text, "diagnostics", {"conservation_history_save": "10", "turbulence_history_save": "0"})
    args.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
