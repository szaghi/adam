#!/usr/bin/env python3
"""Write the FLUME ghost-probe cases (issue #65, P0): realms holding a linear field, run for one negligible step.

Every case uses `[initial_conditions] type = linear` (`q = q_1 (1 + g . x)`, with a non-zero velocity so that every
wall sign flip shows) and `CFL = 1e-30`, `time_max = 0`: one step that moves the field by ~1e-30 relative, after which
the post-step save writes the fields with their ghosts filled in the order of a Runge-Kutta stage. `ghost_probe.py`
then checks every face and edge ghost against the value it must hold.

Cases (each exercises the exchange paths a cross-derivative stencil reads):

    box3d      one realm [0,1]^3, octree, 2x2x2 blocks of 8^3, the block at the origin refined: 2:1 seams that meet
               three faces (x_min inflow, y_min wall, z_min wall); x_max and z_max extrapolation, y_max wall;
               --model mhd for the MHD wall rule (the wall-normal field negated too);
    channel2d  one realm [0,2] x [0,1], quadtree (nk = 1, null z), x periodic, y_min wall, y_max extrapolation, the
               block at the origin refined: 2:1 seams across the periodic boundary at the wall;
    mirror3d   two realms [0,1]^3 and [1,2] x [0,1]^2 glued by a 1:1 mirror seam, octree, walls and extrapolation
               on the other faces: seam edges in every plane;
    refined3d  as mirror3d with the second realm one level finer and a `coupling = refined` (2:1) seam.

The Woodward-Colella step forest (`../step/make_step.py`) is probed by `check.sh` through `--linearize`.

Usage:
    make_probe.py <out-dir> <case> [--model euler|mhd]
    make_probe.py <out-dir> --linearize    # rewrite the realm INIs found in <out-dir> into probe runs
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

GRADIENT = (0.05, 0.11, 0.07)
STATE = {"r": 1.0, "u": 0.3, "v": -0.2, "w": 0.1, "p": 1.0, "bx": 0.8, "by": 0.5, "bz": -0.3}
INFLOW = {"r": 1.2, "u": 0.4, "v": 0.1, "w": -0.1, "p": 1.1, "bx": 0.7, "by": 0.4, "bz": -0.2}


def realm_ini(name: str, extent: tuple[float, ...], cells: int, levels: int, ratio: int, null_z: bool,
              faces: dict[str, str], model: str, box: tuple[float, ...] | None = None) -> str:
    """Return a realm INI: `cells` per block along the active axes, `levels` uniform levels, one more level on `box`."""
    keys = ("r", "u", "v", "w", "p") + (("bx", "by", "bz") if model == "mhd" else ())
    lines = [
        f"; ghost probe realm {name} (make_probe.py, issue #65 P0)",
        "[IO]", f"output_basename = {name}", "it_save = 1000000", "restart = .false.",
        f"restart_basename = {name}-restart", "restart_save = 0", "residuals_save = 1", "save_memory_status = .false.",
        "save_residual_fields = .false.", "save_auxiliary_fields = .false.", "save_curl_fields = .false.",
        "save_divergence_fields = .false.", "save_gradient_fields = .false.", "save_laplacian_fields = .false.",
        "[grid]", f"ni = {cells}", f"nj = {cells}", f"nk = {1 if null_z else cells}", "ngc = 3",
        f"emin_x = {extent[0]!r}", f"emin_y = {extent[2]!r}", f"emin_z = {extent[4]!r}",
        f"emax_x = {extent[1]!r}", f"emax_y = {extent[3]!r}", f"emax_z = {extent[5]!r}",
        "null_x = .false.", "null_y = .false.", f"null_z = {'.true.' if null_z else '.false.'}",
        "[amr]", f"max_level = {levels + (1 if box else 0)}", f"ratio = {ratio}", f"iu_ref_levels = {levels}",
        "i_prune = 0", "j_prune = 0", "k_prune = 0", "l_prune = -1", "frequency = 999999", "iters = 1",
        f"markers_number = {1 if box else 0}",
    ]
    if box:
        lines += ["[amr_marker_1]", "mode = 1", "geo_type = primitive-box", "delta_type = max", "delta_fine = 0.0",
                  "delta_coarse = 0.0", f"box_xmin = {box[0]!r}", f"box_ymin = {box[2]!r}", f"box_zmin = {box[4]!r}",
                  f"box_xmax = {box[1]!r}", f"box_ymax = {box[3]!r}", f"box_zmax = {box[5]!r}",
                  f"target_level = {levels + 1}"]
    lines += [
        "[field]", f"nv = {9 if model == 'mhd' else 5}", "[runge_kutta]", "scheme = runge-kutta-ssp-33",
        "[weno]", "scheme = weno-u-5", "ror_number = 0", "ror_threshold = 0.9", "ror_vars_number = 0",
        "enable_ror_stats = .false.", "ib_reduction_extent = 0", "ib_reduced_order = 2",
        "[linear-algebra]", "smoothing = gauss-seidel", "iterations_init = 3", "iterations_coarse = 10",
        "iterations_fine = 3", "iterations = 10", "tolerance = 1.e-20",
        "[fdv]", "fdv_scheme = fd", "fdv_order = 2", "[solids]", "number = 0", "n_eikonal = 0",
        "[slices]", "slices_number = 0",
        "[physics]", f"physical_model = {'mhd-ideal' if model == 'mhd' else 'euler'}", "gamma = 1.4",
        "[numerics]", "scheme_space = weno", "reconstruction_variables = characteristic", "reflux = .true.",
    ]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        bc = faces.get(face, "extrapolation")
        lines += [f"[bc_{face}]", f"type = {'extrapolation' if bc == 'seam' else bc}"]
        if bc == "inflow":
            lines += [f"{k} = {INFLOW[k]!r}" for k in keys]
    gz = 0.0 if null_z else GRADIENT[2]
    lines += ["[initial_conditions]", "type = linear", f"gradient_x = {GRADIENT[0]!r}", f"gradient_y = {GRADIENT[1]!r}",
              f"gradient_z = {gz!r}", f"amr_iterations = {1 if box else 0}", "[initial_conditions_region_1]"]
    lines += [f"{k} = {STATE[k]!r}" for k in keys]
    lines += ["[time]", "it_max = 0", "time_max = 0.0", "CFL = 1.0e-30", "[diagnostics]",
              "conservation_history_save = 1"]
    if model == "mhd":
        lines += ["[mhd]", "divergence_control = glm", "divb_tol = 0.0", "divb_error = .false.", "rho_floor = 0.0",
                  "p_floor = 0.0", "glm_ch = 3.0", "glm_alpha = 0.18", "glm_damping_length = 1.0",
                  "glm_ch_check = warning"]
    return "\n".join(lines) + "\n"


def manifest(names: list[str], coupling: str) -> str:
    """Return a two-realm manifest glued along x (+x of the first, -x of the second)."""
    lines = ["; ghost probe forest (make_probe.py, issue #65 P0)", "[forest]", f"realms_number = {len(names)}"]
    for i, name in enumerate(names, start=1):
        lines += [f"[realm.{i}]", f"ini = {name}.ini"]
    lines += ["[forest.topology]", "inter_realm_faces_number = 1", "[forest.topology.face_1]", "realm_a = 1",
              "face_a = +x", "realm_b = 2", "face_b = -x", f"coupling = {coupling}",
              "coupling_cadence = stage_coincident"]
    return "\n".join(lines) + "\n"


def write_case(out: Path, case: str, model: str) -> None:
    """Write the INIs of `case` into `out`; the run input is always `probe.ini` (a realm INI or a manifest)."""
    out.mkdir(parents=True, exist_ok=True)
    walls = {"y_min": "wall-inviscid", "y_max": "wall-inviscid", "z_min": "wall-inviscid"}
    if case == "box3d":
        ini = realm_ini("box3d", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), 8, 1, 8, False, {"x_min": "inflow", **walls}, model,
                        box=(0.0, 0.5, 0.0, 0.5, 0.0, 0.5))
        (out / "probe.ini").write_text(ini)
    elif case == "channel2d":
        faces = {"x_min": "periodic", "x_max": "periodic", "y_min": "wall-inviscid"}
        ini = realm_ini("channel2d", (0.0, 2.0, 0.0, 1.0, 0.0, 1.0), 8, 2, 4, True, faces, model,
                        box=(0.0, 0.5, 0.0, 0.25, -1.0e30, 1.0e30))
        (out / "probe.ini").write_text(ini)
    elif case in ("mirror3d", "refined3d"):
        fine = 1 if case == "refined3d" else 0
        r1 = realm_ini("probe-r1", (0.0, 1.0, 0.0, 1.0, 0.0, 1.0), 8, 1, 8, False,
                       {"x_min": "inflow", "x_max": "seam", **walls}, model)
        r2 = realm_ini("probe-r2", (1.0, 2.0, 0.0, 1.0, 0.0, 1.0), 8, 1 + fine, 8, False,
                       {"x_min": "seam", **walls}, model)
        (out / "probe-r1.ini").write_text(r1)
        (out / "probe-r2.ini").write_text(r2)
        (out / "probe.ini").write_text(manifest(["probe-r1", "probe-r2"], "refined" if fine else "mirror"))
    else:
        raise SystemExit(f"make_probe: unknown case '{case}'")


def linearize(out: Path) -> None:
    """Rewrite the realm INIs of `out` into probe runs: linear field, non-zero velocity, one negligible step."""
    gradient = "".join(f"gradient_{a} = {g!r}\n" for a, g in zip("xyz", (*GRADIENT[:2], 0.0), strict=True))
    for path in out.glob("*.ini"):
        text = path.read_text()
        if "[initial_conditions]" not in text:
            continue
        text = re.sub(r"\[initial_conditions\]\ntype\s*=\s*uniform\ns\s*=.*\n", "[initial_conditions]\ntype = linear\n"
                      + gradient, text)
        text = re.sub(r"(\[initial_conditions_region_1\]\n(?:.*\n)*?)v\s*=\s*0\.0\nw\s*=\s*0\.0\n",
                      r"\1v      = 0.5\nw      = 0.3\n", text)
        text = re.sub(r"it_max\s*=.*\n", "it_max   = 0\n", text)
        text = re.sub(r"time_max\s*=.*\n", "time_max = 0.0\n", text)
        text = re.sub(r"CFL\s*=.*\n", "CFL      = 1.0e-30\n", text)
        if "type = linear" not in text:
            raise SystemExit(f"make_probe: {path} has no uniform initial condition to linearise")
        path.write_text(text)


def main() -> None:
    """Write a case, or linearise the INIs already in the directory."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path, help="output directory")
    parser.add_argument("case", nargs="?", help="box3d, channel2d, mirror3d or refined3d")
    parser.add_argument("--model", choices=("euler", "mhd"), default="euler")
    parser.add_argument("--linearize", action="store_true", help="rewrite the realm INIs of <out> into probe runs")
    args = parser.parse_args()
    if args.linearize:
        linearize(args.out)
    elif args.case:
        write_case(args.out, args.case, args.model)
    else:
        parser.error("give a case or --linearize")


if __name__ == "__main__":
    main()
