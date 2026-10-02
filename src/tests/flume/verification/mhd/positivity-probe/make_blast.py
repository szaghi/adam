#!/usr/bin/env python3
"""Write the Balsara & Spicer (1999) strong MHD blast as a FLUME input (issue #41, M2 stretch S0-S2).

rho = 1, p = 0.1 outside, p = 1000 for r < 0.1 around (0.5, 0.5), B = 100/sqrt(4 pi) at 45 degrees (ambient
beta = 2.5e-4), gamma = 1.4, t = 0.01. The disk is the cell-centre sampling of r < 0.1: one riemann-problem strip
region per cell row (first match wins), the ambient region last. The grid is the Orszag-Tang one (2-D periodic
[0, 1]^2, make_orszag_tang.py). --limiter writes [mhd] positivity_limiter, read only by a build with
one-sided-limiter-cpu.patch applied (the key is ignored otherwise). --eglm selects EGLM cleaning (issue #47, M3-P4c),
--none no cleaning (GLM by default).
Issue #47, M3-P5c: --3d makes the blast a sphere on the 3-D periodic [0, 1]^3 (octree, 4x4x4 blocks, one strip region
per cell row of the (y, z) plane); --b-axis x lays the field along x (Wu & Shu 2018, Example 4.4: with --p-in 1e4,
--b0 1000/sqrt(4 pi), --time-max 0.001 their second, beta 2.51e-6 blast); --time-max sets the final time.
Usage: make_blast.py <out.ini> --cells N [--p-floor F] [--rho-floor F] [--glm-ch C] [--p-in P] [--b0 B] [--cfl C]
                     [--limiter] [--none | --eglm] [--3d] [--b-axis diagonal|x] [--time-max T]
"""

from __future__ import annotations

import argparse
import configparser
import math
import subprocess
import sys
from pathlib import Path

V = Path(__file__).resolve().parents[2]  # src/tests/flume/verification

parser = argparse.ArgumentParser()
parser.add_argument("out", type=Path)
parser.add_argument("--cells", type=int, required=True)
parser.add_argument("--p-floor", default="0.0")
parser.add_argument("--rho-floor", default="0.0")
parser.add_argument("--glm-ch", default="60.0")
parser.add_argument("--p-in", default="1000.0")
parser.add_argument("--cfl", default=None)
parser.add_argument("--limiter", action="store_true")
parser.add_argument("--none", action="store_true")
parser.add_argument("--eglm", action="store_true")
parser.add_argument("--b0", type=float, default=100.0 / (4.0 * 3.141592653589793) ** 0.5)
parser.add_argument("--3d", dest="three_d", action="store_true")
parser.add_argument("--b-axis", choices=("diagonal", "x"), default="diagonal")
parser.add_argument("--time-max", default="0.01")
args = parser.parse_args()
subprocess.run([sys.executable, str(V / "mhd/orszag-tang/make_orszag_tang.py"), str(V / "vortex/vortex-n064.ini"),
                str(args.out), "--cells", str(args.cells)], check=True)
ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
ini.optionxform = str
ini.read(args.out)
for s in [s for s in ini.sections() if s.startswith("initial_conditions")]:
    ini.remove_section(s)
n, h, b0 = args.cells, 1.0 / args.cells, args.b0
bx, by = (b0, 0.0) if args.b_axis == "x" else (b0 / math.sqrt(2.0), b0 / math.sqrt(2.0))
if args.three_d:
    ini["grid"].update({"nk": str(n // 4), "null_z": ".false.", "emin_z": "0.0", "emax_z": "1.0"})
    ini["amr"]["ratio"] = "8"
    ini["bc_z_min"]["type"] = ini["bc_z_max"]["type"] = "periodic"
strips = []
for k in range(n) if args.three_d else [None]:
    dz = (k + 0.5) * h - 0.5 if args.three_d else 0.0
    z0, z1 = (k * h, (k + 1) * h) if args.three_d else (-1e30, 1e30)
    for j in range(n):
        yc = (j + 0.5) * h
        dy = yc - 0.5
        if dy**2 + dz**2 >= 0.1**2:
            continue
        half = math.sqrt(0.1**2 - dy**2 - dz**2)
        i_lo = [i for i in range(n) if abs((i + 0.5) * h - 0.5) < half]
        if i_lo:
            strips.append((i_lo[0] * h, (i_lo[-1] + 1) * h, j * h, (j + 1) * h, z0, z1))
ini["initial_conditions"] = {"type": "riemann-problem", "amr_iterations": "0", "regions_number": str(len(strips) + 1)}
state = {"u": "0.0", "v": "0.0", "w": "0.0", "bx": repr(bx), "by": repr(by), "bz": "0.0", "r": "1.0"}
for r, (x0, x1, y0, y1, z0, z1) in enumerate(strips, start=1):
    extent = {"emin_x": repr(x0), "emax_x": repr(x1), "emin_y": repr(y0), "emax_y": repr(y1),
              "emin_z": "-1e30" if z0 == -1e30 else repr(z0), "emax_z": "1e30" if z1 == 1e30 else repr(z1)}
    ini[f"initial_conditions_region_{r}"] = {**state, "p": args.p_in, **extent}
ini[f"initial_conditions_region_{len(strips) + 1}"] = {**state, "p": "0.1", "emin_x": "-1e30", "emax_x": "1e30",
                                                       "emin_y": "-1e30", "emax_y": "1e30", "emin_z": "-1e30",
                                                       "emax_z": "1e30"}
ini["physics"].update({"cp": "1.4", "cv": "1.0"})
ini["mhd"].update({"glm_ch": args.glm_ch, "glm_ch_check": "warning", "p_floor": args.p_floor,
                   "rho_floor": args.rho_floor})
ini["time"].update({"it_max": "-1", "time_max": args.time_max})
if args.cfl:
    ini["time"]["CFL"] = args.cfl
if args.none:
    ini["mhd"]["divergence_control"] = "none"
if args.eglm:
    ini["mhd"]["divergence_control"] = "eglm"
if args.limiter:
    ini["mhd"]["positivity_limiter"] = ".true."
ini["IO"].update({"output_basename": "blast", "it_save": "100000"})
with open(args.out, "w") as f:
    ini.write(f)
print(f"{len(strips)} disk strips, ambient beta {2 * 0.1 / b0**2:.2e}")
