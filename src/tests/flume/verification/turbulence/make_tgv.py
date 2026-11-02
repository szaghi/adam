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
    parser.add_argument("--hit-k0", type=float, default=0.0, help="> 0: isotropic turbulence of peak k0, Mach = M_t, Reynolds = Re_lambda")
    parser.add_argument("--hit-seed", type=int, default=1)
    parser.add_argument("--sgs-model", default="none")
    parser.add_argument("--sgs-constant", type=float, default=0.0, help="0: the model's default")
    parser.add_argument("--amr-box", type=float, default=0.0, help="> 0: one more level on the centred box of this side over 2 pi")
    parser.add_argument("--seam-fill", default="")
    parser.add_argument("--seam-restriction", default="")
    parser.add_argument("--contact", type=float, default=0.0, help="> 0: instead, a contact of this density ratio advected on a periodic square")
    parser.add_argument("--density-wave", type=float, default=0.0, help="> 0: instead, a density wave of this relative amplitude advected at u = p = 1")
    parser.add_argument("--hybrid-threshold", type=float, default=0.0, help="> 0: [numerics] hybrid_threshold")
    parser.add_argument("--hybrid-skin", default="", help="[numerics] hybrid_skin")
    parser.add_argument("--shock-vortex", action="store_true", help="instead, a vortex after Jiang & Shu (1996), the code's profile, through a stationary Mach 1.1 shock, on the unit square")
    parser.add_argument("--restart-save", type=int, default=0, help="cadence of the restart files")
    parser.add_argument("--restart", action="store_true", help="restart from them")
    parser.add_argument("--history", type=int, default=1, help="cadence of the conservation and turbulence histories")
    args = parser.parse_args()

    blocks = 2**args.levels
    if args.n % blocks:
        raise SystemExit(f"make_tgv: {args.n} cells do not split into {blocks} blocks")
    text = TEMPLATE.read_text()
    text = re.sub(r"\A(?:;.*\n)+", "", text)  # the template's own header
    text = re.sub(r"^\[initial_conditions_region_\d+\][^\[]*", "", text, flags=re.M)
    box = 2.0 * math.pi
    text = section(text, "IO", {"output_basename": "tgv", "it_save": str(args.it_save), "restart_basename": "tgv-restart",
                                 "restart_save": str(args.restart_save), "restart": ".true." if args.restart else ".false."})
    grid = {"ni": str(args.n // blocks), "nj": str(args.n // blocks), "nk": str(args.n // blocks),
            "ngc": str(max(3, (int(args.scheme.rsplit("-", 1)[1]) + 1) // 2))}
    grid |= {f"emin_{a}": "0.0" for a in "xyz"} | {f"emax_{a}": f"{box:.16e}" for a in "xyz"}
    grid |= {f"null_{a}": ".false." for a in "xyz"}
    text = section(text, "grid", grid)
    text = section(text, "amr", {"max_level": str(args.levels), "ratio": "8", "iu_ref_levels": str(args.levels),
                                 "frequency": "0", "markers_number": "0"})
    text = section(text, "weno", {"scheme": args.scheme, "weights": "si"})
    physics = {"cp": "3.5", "cv": "2.5"}
    if args.hit_k0 > 0.0:
        # Johnsen et al. (2010): u_rms = rho_0 = 1, lambda_0 = 2 / k_0, mu = (T / T_0)^0.76, Pr 0.7; code T = p / rho (R = 1)
        physics |= {"reynolds": str(args.reynolds * args.hit_k0 / 2.0), "prandtl": "0.7", "viscosity_law": "power-law",
                    "viscosity_exponent": "0.76", "reference_temperature": str((3.0**0.5 / args.mach) ** 2 / 1.4)}
    elif args.reynolds > 0.0:
        physics |= {"reynolds": str(args.reynolds), "prandtl": "0.71"}
    text = section(text, "physics", physics)
    numerics = {"scheme_space": args.scheme_space}
    numerics |= {"hybrid_threshold": str(args.hybrid_threshold)} if args.hybrid_threshold > 0.0 else {}
    numerics |= {"hybrid_skin": args.hybrid_skin} if args.hybrid_skin else {}
    if args.scheme_space == "weno-riemann":
        numerics["riemann_solver"] = args.riemann
        numerics |= {"flux_correction": "6th", "flux_correction_sensor": "weno"}
    text = section(text, "numerics", numerics)
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        text = section(text, f"bc_{face}", {"type": "periodic"})
    ic = {"type": "taylor-green", "amr_iterations": "0", "tg_density": "1.0", "tg_velocity": "1.0",
          "tg_mach": str(args.mach)}
    if args.hit_k0 > 0.0:
        ic = {"type": "isotropic-turbulence", "amr_iterations": "0", "hit_density": "1.0", "hit_urms": "1.0",
              "hit_mach": str(args.mach), "hit_k0": str(args.hit_k0), "hit_seed": str(args.hit_seed)}
    text = section(text, "initial_conditions", ic)
    text = section(text, "time", {"it_max": str(args.it_max), "time_max": str(args.time_max), "CFL": str(args.cfl)})
    text = section(text, "diagnostics", {"conservation_history_save": str(args.history),
                                         "turbulence_history_save": str(args.history)})
    if args.contact > 0.0:
        text = section(text, "grid", {"nk": "1", "emax_x": "1.0", "emax_y": "1.0", "emax_z": "1.0", "null_z": ".true."})
        for face in ("z_min", "z_max"):
            text = section(text, f"bc_{face}", {"type": "extrapolation"})
        for key in ("tg_density", "tg_velocity", "tg_mach"):
            text = re.sub(rf"^{key} = .*\n", "", text, flags=re.M)
        text = section(text, "initial_conditions", {"type": "riemann-problem", "regions_number": "3"})
        for i, (rho, lo, hi) in enumerate(((1.0, 0.0, 0.25), (args.contact, 0.25, 0.75), (1.0, 0.75, 1.0)), start=1):
            text += f"""
[initial_conditions_region_{i}]
r = {rho}
u = 1.0
v = 0.0
w = 0.0
p = 1.0
emin_x = {lo}
emin_y = 0.0
emin_z = 0.0
emax_x = {hi}
emax_y = 1.0
emax_z = 1.0
"""
    if args.density_wave > 0.0:
        text = section(text, "grid", {"nk": "1", "emax_x": "1.0", "emax_y": "1.0", "emax_z": "1.0", "null_z": ".true."})
        for face in ("z_min", "z_max"):
            text = section(text, f"bc_{face}", {"type": "extrapolation"})
        for key in ("tg_density", "tg_velocity", "tg_mach"):
            text = re.sub(rf"^{key} = .*\n", "", text, flags=re.M)
        text = section(text, "initial_conditions", {"type": "sine-wave", "wave_mode": "entropy", "wave_angle": "0.0",
                                                     "wave_amplitude": str(args.density_wave), "wavelength": "1.0"})
        text += """
[initial_conditions_region_1]
r = 1.0
u = 1.0
v = 0.0
w = 0.0
p = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
"""
    if args.shock_vortex:
        box = 1.0  # for --amr-box too
        text = section(text, "amr", {"ratio": "4"})  # a quadtree: seams cannot refine the one cell along z
        text = section(text, "grid", {"nk": "1", "emax_x": "1.0", "emax_y": "1.0", "emax_z": "1.0", "null_z": ".true."})
        for face, kind in (("x_min", "extrapolation"), ("x_max", "extrapolation"), ("z_min", "extrapolation"), ("z_max", "extrapolation")):
            text = section(text, f"bc_{face}", {"type": kind})
        for key in ("tg_density", "tg_velocity", "tg_mach"):
            text = re.sub(rf"^{key} = .*\n", "", text, flags=re.M)
        text = section(text, "initial_conditions", {"type": "isentropic-vortex", "x0": "0.25", "y0": "0.5", "radius": "0.05",
                                                     "strength": "1.885", "shock_x": "0.5"})
        text += f"""
[initial_conditions_region_1]
r = 1.0
u = {1.1 * 1.4**0.5:.16e}
v = 0.0
w = 0.0
p = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
"""
    if args.amr_box > 0.0:
        lo, hi = 0.5 * box * (1.0 - args.amr_box), 0.5 * box * (1.0 + args.amr_box)
        amr = {"max_level": str(args.levels + 1), "markers_number": "1"}
        amr |= {"seam_ghost_fill": args.seam_fill} if args.seam_fill else {}
        amr |= {"seam_restriction": args.seam_restriction} if args.seam_restriction else {}
        text = section(text, "amr", amr)
        text = section(text, "numerics", {"reflux": ".true."})
        text = section(text, "initial_conditions", {"amr_iterations": "1"})  # 0 would leave the box unrefined
        text += f"""
[amr_marker_1]
mode = 1
geo_type = primitive-box
delta_type = max
delta_fine = 0.0
delta_coarse = 0.0
box_xmin = {lo:.16e}
box_ymin = {lo:.16e}
box_zmin = {lo:.16e}
box_xmax = {hi:.16e}
box_ymax = {hi:.16e}
box_zmax = {hi:.16e}
target_level = {args.levels + 1}
"""
    if args.sgs_model != "none":
        text += f"\n[turbulence]\nsgs_model = {args.sgs_model}\n"
        if args.sgs_constant > 0.0:
            text += f"sgs_constant = {args.sgs_constant}\n"
    args.out.write_text(text)


if __name__ == "__main__":
    main()
