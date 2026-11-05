#!/usr/bin/env python3
"""Write a plane channel input of issue #119 T3, from make_tgv.py's template: walls at y = +-1, periodic x and z, mass-flow forcing, rho = U = 1, so the bulk
Reynolds number is the `reynolds` key (Coleman, Kim & Moser define it on the half height)."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

TEMPLATE = Path(__file__).resolve().parent.parent.parent / "regression" / "vortex-periodic" / "input.ini"


def section(text: str, name: str, keys: dict[str, str]) -> str:
    """Set keys in a section, adding those absent."""
    m = re.search(rf"^\[{re.escape(name)}\]\n(.*?)(?=^\[|\Z)", text, flags=re.M | re.S)
    if not m:
        raise SystemExit(f"make_channel: no section [{name}]")
    body = m.group(1)
    for key, value in keys.items():
        body, n = re.subn(rf"^{key}\s*=.*$", f"{key} = {value}", body, flags=re.M)
        if n == 0:
            body = f"{key} = {value}\n" + body
    return text[: m.start(1)] + body + text[m.end(1) :]


def main() -> int:
    """Write the input."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", type=Path)
    ap.add_argument("--nx", type=int, default=16)
    ap.add_argument("--ny", type=int, default=32)
    ap.add_argument("--nz", type=int, default=16)
    ap.add_argument("--levels", type=int, default=2)
    ap.add_argument("--lx", type=float, default=4.0 * math.pi)
    ap.add_argument("--lz", type=float, default=4.0 * math.pi / 3.0)
    ap.add_argument("--mach", type=float, default=0.1)
    ap.add_argument("--reynolds", type=float, default=100.0)
    ap.add_argument("--power-law", type=float, default=0.0, help="> 0: viscosity (T / T_wall)^this, as Coleman (0.7)")
    ap.add_argument("--perturbation", type=float, default=0.0)
    ap.add_argument("--scheme-space", default="weno")
    ap.add_argument("--wall", default="wall-isothermal")
    ap.add_argument("--wall-heat", type=float, default=1.0, help="the wall temperature over the gas temperature of the start")
    ap.add_argument("--wm-height", type=float, default=0.0, help="[bc] wall_model_height, 0: the 3rd cell")
    ap.add_argument("--wm-damping", default="", help="[bc] wall_model_damping")
    ap.add_argument("--wm-filter", type=float, default=0.0, help="[bc] wall_model_filter, the time scale of the input filter, 0: none")
    ap.add_argument("--wm-model", default="", help="[bc] wall_model: equilibrium or inverse-transformation")
    ap.add_argument("--wm-eddy-viscosity", action="store_true", help="[bc] wall_model_eddy_viscosity: the model's eddy viscosity below the matching height")
    ap.add_argument("--wm-edge-velocity", type=float, default=1.5, help="[bc] wall_model_edge_velocity (inverse model), the bulk speed being 1")
    ap.add_argument("--wm-edge-heat", type=float, default=1.0, help="[bc] wall_model_edge_temperature (inverse model) over the gas temperature of the start")
    ap.add_argument("--amr-box", action="store_true", help="one more level on the middle half along x, over the whole height: 2:1 faces along the walls")
    ap.add_argument("--sgs-model", default="none")
    ap.add_argument("--gas", default="cp-cv", choices=["cp-cv", "gamma"], help="the gas as cp, cv (R = 1) or as gamma = 1.4 (the reference-layer legs: a cp, cv input is not bitwise across it)")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--time-max", type=float, default=1.0)
    ap.add_argument("--it-max", type=int, default=-1)
    ap.add_argument("--it-save", type=int, default=100000)
    ap.add_argument("--restart-save", type=int, default=0)
    ap.add_argument("--restart", action="store_true")
    ap.add_argument("--history", type=int, default=1)
    ap.add_argument("--wall-history", type=int, default=0, help="[diagnostics] wall_history_save: the wall model's tau_w and q_w in nx bins along x")
    args = ap.parse_args()
    blocks = 2**args.levels
    for n in (args.nx, args.ny, args.nz):
        if n % blocks:
            raise SystemExit(f"make_channel: {n} cells do not split into {blocks} blocks")
    gamma = 1.4
    t_wall = 1.0 / (gamma * args.mach**2)  # R = 1, rho = U = 1
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    text = section(text, "IO", {"output_basename": "channel", "it_save": str(args.it_save), "restart_basename": "channel-restart",
                                "restart_save": str(args.restart_save), "restart": ".true." if args.restart else ".false."})
    grid = {"ni": str(args.nx // blocks), "nj": str(args.ny // blocks), "nk": str(args.nz // blocks), "ngc": "3",
            "emin_x": "0.0", "emax_x": repr(args.lx), "emin_y": "-1.0", "emax_y": "1.0", "emin_z": "0.0", "emax_z": repr(args.lz)}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": str(args.levels), "ratio": "8", "iu_ref_levels": str(args.levels), "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": "weno-u-5", "weights": "si"})
    physics = {"cp": "3.5", "cv": "2.5", "reynolds": str(args.reynolds), "prandtl": "0.7"} if args.gas == "cp-cv" else {"gamma": "1.4", "reynolds": str(args.reynolds), "prandtl": "0.7"}
    if args.power_law > 0.0:
        physics |= {"viscosity_law": "power-law", "viscosity_exponent": str(args.power_law), "reference_temperature": repr(t_wall)}
    if args.gas == "gamma":  # the template's cp, cv go: the physics refuses gamma with them
        text = re.sub(r"^(cp|cv)\s*=.*\n", "", text, flags=re.M)
    text = section(text, "physics", physics)
    text = section(text, "numerics", {"scheme_space": args.scheme_space})
    for face in ("x_min", "x_max", "z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    for face in ("y_min", "y_max"):
        keys = {"type": args.wall}
        if args.wall == "wall-isothermal" or args.wall == "wall-modelled":
            keys["wall_temperature"] = repr(t_wall * args.wall_heat)
        if args.wall == "wall-modelled":
            keys |= {"wall_model_height": str(args.wm_height)} if args.wm_height > 0.0 else {}
            keys |= {"wall_model_damping": args.wm_damping} if args.wm_damping else {}
            keys |= {"wall_model_filter": repr(args.wm_filter)} if args.wm_filter > 0.0 else {}
            keys |= {"wall_model_eddy_viscosity": ".true."} if args.wm_eddy_viscosity else {}
            if args.wm_model:
                keys["wall_model"] = args.wm_model
            if args.wm_model == "inverse-transformation":
                keys |= {"wall_model_edge_velocity": repr(args.wm_edge_velocity), "wall_model_edge_temperature": repr(t_wall * args.wm_edge_heat)}
        text = section(text, f"bc_{face}", keys)
    text = section(text, "initial_conditions", {"type": "channel", "amr_iterations": "0", "ch_rho": "1.0", "ch_ub": "1.0", "ch_mach": str(args.mach), "ch_h": "1.0",
                                                 "ch_lx": repr(args.lx), "ch_lz": repr(args.lz), "ch_pert": str(args.perturbation)})
    text = section(text, "time", {"it_max": str(args.it_max), "time_max": str(args.time_max), "CFL": str(args.cfl)})
    text = section(text, "diagnostics", {"conservation_history_save": str(args.history), "turbulence_history_save": str(args.history)})
    if args.wall_history > 0:
        text = section(text, "diagnostics", {"wall_history_save": str(args.wall_history), "wall_history_bins": str(args.nx), "wall_history_axis": "x"})
    if args.amr_box:
        text = section(text, "amr", {"max_level": str(args.levels + 1), "markers_number": "1"})
        text = section(text, "initial_conditions", {"amr_iterations": "1"})
        text = section(text, "numerics", {"reflux": ".true."})
        text += f"""
[amr_marker_1]
mode = 1
geo_type = primitive-box
delta_type = max
delta_fine = 0.0
delta_coarse = 0.0
box_xmin = {0.25 * args.lx}
box_ymin = -2.0
box_zmin = -1.0
box_xmax = {0.75 * args.lx}
box_ymax = 2.0
box_zmax = {args.lz + 1.0}
target_level = {args.levels + 1}
"""
    text += "\n[forcing]\ntype = mass-flow\naxis = x\nbulk_momentum = 1.0\n"
    if args.sgs_model != "none":
        text += f"\n[turbulence]\nsgs_model = {args.sgs_model}\n"
    args.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
