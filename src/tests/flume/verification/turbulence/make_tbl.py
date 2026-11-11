#!/usr/bin/env python3
"""Write a compressible turbulent boundary layer input of issue #119 T3.7, from make_channel.py's template: a wall at y = 0, inflow at x = 0 through the
digital-filter face, extrapolated outflow and top, periodic z, the inflow thickness delta_i = 1 and rho = U = 1 in the free stream, so `reynolds` is
U_inf delta_i / nu_inf (the viscosity a power law of the temperature, or Sutherland's law, with the free-stream reference) and T_inf = 1 / (gamma M^2) (R = 1)."""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent.parent / "regression" / "vortex-periodic" / "input.ini"


def section(text: str, name: str, keys: dict[str, str]) -> str:
    """Set keys in a section, adding the section or the keys when absent."""
    m = re.search(rf"^\[{re.escape(name)}\]\n(.*?)(?=^\[|\Z)", text, flags=re.M | re.S)
    if not m:
        return text.rstrip("\n") + f"\n\n[{name}]\n" + "".join(f"{k} = {v}\n" for k, v in keys.items())
    body = m.group(1)
    for key, value in keys.items():
        body, n = re.subn(rf"^{key}\s*=.*$", f"{key} = {value}", body, flags=re.M)
        if n == 0:
            body = f"{key} = {value}\n" + body
    return text[: m.start(1)] + body + text[m.end(1) :]


def main() -> int:
    """Write the input and copy the profile table next to it."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out", type=Path)
    ap.add_argument("--profile", type=Path, required=True, help="profile table (tbl_oracle.py profile)")
    ap.add_argument("--nx", type=int, default=64)
    ap.add_argument("--ny", type=int, default=16)
    ap.add_argument("--nz", type=int, default=16)
    ap.add_argument("--levels", type=int, default=2)
    ap.add_argument("--lx", type=float, default=8.0)
    ap.add_argument("--ly", type=float, default=2.0)
    ap.add_argument("--lz", type=float, default=2.0)
    ap.add_argument("--mach", type=float, default=2.5)
    ap.add_argument("--reynolds", type=float, default=19370.0, help="U_inf delta_i / nu_inf (M2p5: 4.84e6 /m times delta_i 4 mm)")
    ap.add_argument("--power-law", type=float, default=0.76, help="viscosity (T / T_inf)^this")
    ap.add_argument("--sutherland", type=float, default=0.0,
                    help="Sutherland's law instead of the power law, S / T_inf (air: 110.4 K over T_inf in K; 0: the power law)")
    ap.add_argument("--wall-heat", type=float, default=2.1037, help="T_w / T_inf (M2p5: 568 / 270)")
    ap.add_argument("--wall", default="wall-modelled")
    ap.add_argument("--wm-height", type=float, default=0.0)
    ap.add_argument("--wm-damping", default="semi-local")
    ap.add_argument("--wm-eddy-viscosity", action="store_true")
    ap.add_argument("--wm-pressure-gradient", action="store_true", help="[bc] wall_model_pressure_gradient: dp/dx in the equilibrium model (#121, S4d)")
    ap.add_argument("--wm-pg-filter", type=float, default=0.0, help="[bc] wall_model_pressure_gradient_filter, the time scale of the gradient's filter")
    ap.add_argument("--wm-pg-ypg", type=float, default=0.0, help="[bc] wall_model_pressure_gradient_ypg: the Bernoulli-consistent y_pg (Hickel et al. 2012: 4; #121 P1)")
    ap.add_argument("--df-lengths", default="0.8,0.3,0.3", help="integral lengths streamwise, normal, spanwise in units of delta_i")
    ap.add_argument("--df-seed", type=int, default=1)
    ap.add_argument("--scheme-space", default="hybrid-kep-weno")
    ap.add_argument("--sgs-model", default="none")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--time-max", type=float, default=1.0)
    ap.add_argument("--it-max", type=int, default=-1)
    ap.add_argument("--it-save", type=int, default=100000)
    ap.add_argument("--restart-save", type=int, default=0)
    ap.add_argument("--restart", action="store_true")
    ap.add_argument("--history", type=int, default=1)
    ap.add_argument("--wall-history", type=int, default=0, help="[diagnostics] wall_history_save: the wall model's tau_w and q_w in nx bins along x")
    ap.add_argument("--top-shock", type=float, default=0.0,
                    help="the top face a shock-generator of this deflection (degrees, issue #121 S4), 0: extrapolation")
    ap.add_argument("--shock-x", type=float, default=0.0, help="with --top-shock: where the shock enters through the top")
    ap.add_argument("--slice-z", default="", help="x-y slices at these comma-separated z, cell-centre points every cell along x and y (trilinear)")
    ap.add_argument("--slice-every", type=int, default=0, help="steps between slice saves")
    ap.add_argument("--refine-box", default="", help="XMIN,XMAX,YMAX: one more level (2:1) on the blocks meeting x in [XMIN, XMAX], y < YMAX, every z")
    args = ap.parse_args()
    blocks = 2**args.levels
    for n in (args.nx, args.ny, args.nz):
        if n % blocks:
            raise SystemExit(f"make_tbl: {n} cells do not split into {blocks} blocks")
    gamma = 1.4
    t_inf = 1.0 / (gamma * args.mach**2)
    p_inf = t_inf  # rho = 1, R = 1
    lengths = args.df_lengths.split(",")
    shutil.copyfile(args.profile, args.out.parent / "profile.txt")
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    text = section(text, "IO", {"output_basename": "tbl", "it_save": str(args.it_save), "restart_basename": "tbl-restart",
                                "restart_save": str(args.restart_save), "restart": ".true." if args.restart else ".false."})
    grid = {"ni": str(args.nx // blocks), "nj": str(args.ny // blocks), "nk": str(args.nz // blocks), "ngc": "3",
            "emin_x": "0.0", "emax_x": repr(args.lx), "emin_y": "0.0", "emax_y": repr(args.ly), "emin_z": "0.0", "emax_z": repr(args.lz)}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": str(args.levels), "ratio": "8", "iu_ref_levels": str(args.levels), "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": "weno-u-5", "weights": "si"})
    law = {"viscosity_law": "power-law", "viscosity_exponent": str(args.power_law)}
    if args.sutherland > 0.0:
        law = {"viscosity_law": "sutherland", "sutherland_temperature": repr(args.sutherland * t_inf)}
    text = section(text, "physics", {"cp": "3.5", "cv": "2.5", "reynolds": str(args.reynolds), "prandtl": "0.71", **law,
                                     "reference_temperature": repr(t_inf)})
    text = section(text, "numerics", {"scheme_space": args.scheme_space})
    dx = args.lx / args.nx
    dy = args.ly / args.ny
    text = section(text, "bc_x_min", {"type": "inflow-digital-filter", "r": "1.0", "u": "1.0", "p": repr(p_inf), "df_profile": "profile.txt",
                                      "df_delta": "1.0", "df_wall_axis": "y", "df_spacing": repr(dy), "df_length_streamwise": lengths[0],
                                      "df_length_normal": lengths[1], "df_length_spanwise": lengths[2], "df_seed": str(args.df_seed)})
    text = section(text, "bc_x_max", {"type": "extrapolation"})
    if args.top_shock > 0.0:  # the incident shock of an interaction (issue #121, S4): the free stream upstream of it
        text = section(text, "bc_y_max", {"type": "shock-generator", "r": "1.0", "u": "1.0", "v": "0.0", "w": "0.0", "p": repr(p_inf),
                                          "shock_deflection": repr(args.top_shock), "shock_x": repr(args.shock_x)})
    else:
        text = section(text, "bc_y_max", {"type": "extrapolation"})
    for face in ("z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    keys = {"type": args.wall}
    if args.wall in ("wall-isothermal", "wall-modelled"):
        keys["wall_temperature"] = repr(t_inf * args.wall_heat)
    if args.wall == "wall-modelled":
        keys |= {"wall_model_height": str(args.wm_height)} if args.wm_height > 0.0 else {}
        keys |= {"wall_model_damping": args.wm_damping} if args.wm_damping else {}
        keys |= {"wall_model_eddy_viscosity": ".true."} if args.wm_eddy_viscosity else {}
        keys |= {"wall_model_pressure_gradient": ".true."} if args.wm_pressure_gradient else {}
        keys |= {"wall_model_pressure_gradient_filter": repr(args.wm_pg_filter)} if args.wm_pg_filter > 0.0 else {}
        keys |= {"wall_model_pressure_gradient_ypg": repr(args.wm_pg_ypg)} if args.wm_pg_ypg > 0.0 else {}
    text = section(text, "bc_y_min", keys)
    text = section(text, "initial_conditions", {"type": "boundary-layer", "amr_iterations": "0", "bl_profile": "profile.txt", "bl_delta": "1.0",
                                                 "bl_rho": "1.0", "bl_u": "1.0", "bl_p": repr(p_inf)})
    text = section(text, "time", {"it_max": str(args.it_max), "time_max": str(args.time_max), "CFL": str(args.cfl)})
    diag = {"conservation_history_save": str(args.history), "turbulence_history_save": "0"}
    if args.wall_history > 0:
        diag |= {"wall_history_save": str(args.wall_history), "wall_history_bins": str(args.nx), "wall_history_axis": "x"}
    text = section(text, "diagnostics", diag)
    if args.sgs_model != "none":
        text += f"\n[turbulence]\nsgs_model = {args.sgs_model}\n"
    if args.refine_box:  # a static 2:1 box at the start, as the channel's wm3 (make_channel.py --amr-box); no runtime regridding
        xa, xb, yb = (float(v) for v in args.refine_box.split(","))
        text = section(text, "amr", {"max_level": str(args.levels + 1), "markers_number": "1"})
        text = section(text, "initial_conditions", {"amr_iterations": "1"})
        text = section(text, "numerics", {"reflux": ".true."})
        text += (f"\n[amr_marker_1]\nmode = 1\ngeo_type = primitive-box\ndelta_type = max\ndelta_fine = 0.0\ndelta_coarse = 0.0\n"
                 f"box_xmin = {xa!r}\nbox_ymin = -1.0\nbox_zmin = -1.0\nbox_xmax = {xb!r}\nbox_ymax = {yb!r}\nbox_zmax = {args.lz + 1.0!r}\n"
                 f"target_level = {args.levels + 1}\n")
    zs = [float(z) for z in args.slice_z.split(",")] if args.slice_z else []
    text = re.sub(r"^\[slices\]\n[^\[]*", "", text, flags=re.M)
    text += f"\n[slices]\nslices_number = {len(zs)}\n"
    dz = args.lz / args.nz
    for s, z in enumerate(zs, start=1):
        text += (f"\n[slice_{s}]\nitype = trilinear\nn_save = {args.slice_every}\nni = {args.nx}\nnj = {args.ny}\nnk = 1\n"
                 f"emin_x = 0.0\nemax_x = {args.lx!r}\nemin_y = 0.0\nemax_y = {args.ly!r}\nemin_z = {z - 0.5 * dz!r}\nemax_z = {z + 0.5 * dz!r}\n")
    args.out.write_text(text)
    print(f"make_tbl: {args.out}, dx {dx:.4g}, dy {dy:.4g}, T_inf {t_inf:.6g}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
