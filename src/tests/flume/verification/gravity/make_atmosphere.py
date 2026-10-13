#!/usr/bin/env python3
"""Write the FLUME stratified-atmosphere cases of issue #73 (M7): a hydrostatic atmosphere under `[gravity]`.

M7-G0 (no gravity source yet) runs each case for one negligible step (`CFL = 1e-30`) with the ghosts poisoned
(`[diagnostics] ghost_poison`), so the post-step checkpoint holds the initial atmosphere and its ghosts filled in the
order of a Runge-Kutta stage; `atmosphere_oracle.py profile` checks both. Cases:

    iso2d        isothermal, g = (0, -1), [0,1]^2 quadtree 4x4 blocks of 8^2, x periodic, y hydrostatic (Euler);
    poly2d       polytropic, G = 1.2, otherwise as iso2d;
    tilted2d     isothermal, g = (-0.6, -0.8), every face hydrostatic (the realm corners too), one refined block;
    tilted3d     isothermal, g = (-0.48, -0.64, -0.6), [0,1]^3 octree 4x4x4 blocks of 8^3, every face hydrostatic;
    iso2d-mhd    iso2d with EGLM and a field oblique to gravity, B = (0.5, 0.7, 0.2) (even at hydrostatic faces);
    froude2d     iso2d written with the number form, froude = 1, axis = -y.

Two uniform levels, so the blocks spread over the ranks (one level leaves every block on rank 0).

Usage:
    make_atmosphere.py <out.ini> <case> [--steps N --cfl C]   # N steps at CFL C instead of the negligible step
"""

from __future__ import annotations

import argparse
from pathlib import Path

CASES = ("iso2d", "poly2d", "tilted2d", "tilted3d", "iso2d-mhd", "froude2d")
GRAVITY = {"iso2d": (0.0, -1.0, 0.0), "poly2d": (0.0, -1.0, 0.0), "tilted2d": (-0.6, -0.8, 0.0),
           "tilted3d": (-0.48, -0.64, -0.6), "iso2d-mhd": (0.0, -1.0, 0.0)}


def case_ini(case: str, steps: int, cfl: float) -> str:
    """Return the INI of a case."""
    is3d = case == "tilted3d"
    mhd = case.endswith("-mhd")
    cells = 8
    tilted = case.startswith("tilted")
    box = (0.0, 0.5, 0.0, 0.5, 0.0, 1.0) if case == "tilted2d" else None
    lines = [
        f"; FLUME stratified atmosphere {case} (verification/gravity/make_atmosphere.py, issue #73 M7-G0)",
        "[IO]", f"output_basename = {case}", "it_save = 1000000", "restart = .false.",
        f"restart_basename = {case}-restart", "restart_save = 0", "residuals_save = 1", "save_memory_status = .false.",
        "save_residual_fields = .false.", "save_auxiliary_fields = .false.", "save_curl_fields = .false.",
        "save_divergence_fields = .false.", "save_gradient_fields = .false.", "save_laplacian_fields = .false.",
        "[grid]", f"ni = {cells}", f"nj = {cells}", f"nk = {cells if is3d else 1}", "ngc = 3",
        "emin_x = 0.0", "emin_y = 0.0", "emin_z = 0.0", "emax_x = 1.0", "emax_y = 1.0", "emax_z = 1.0",
        "null_x = .false.", "null_y = .false.", f"null_z = {'.false.' if is3d else '.true.'}",
        "[amr]", f"max_level = {3 if box else 2}", f"ratio = {8 if is3d else 4}", "iu_ref_levels = 2",
        "i_prune = 0", "j_prune = 0", "k_prune = 0", "l_prune = -1", "frequency = 0", "iters = 1",
        f"markers_number = {1 if box else 0}",
    ]
    if box:
        lines += ["[amr_marker_1]", "mode = 1", "geo_type = primitive-box", "delta_type = max", "delta_fine = 0.0",
                  "delta_coarse = 0.0", f"box_xmin = {box[0]!r}", f"box_ymin = {box[2]!r}", f"box_zmin = {box[4]!r}",
                  f"box_xmax = {box[1]!r}", f"box_ymax = {box[3]!r}", f"box_zmax = {box[5]!r}", "target_level = 3"]
    lines += [
        "[field]", f"nv = {9 if mhd else 5}", "[runge_kutta]", "scheme = runge-kutta-ssp-33",
        "[weno]", "scheme = weno-u-5", "ror_number = 0", "ror_threshold = 0.9", "ror_vars_number = 0",
        "enable_ror_stats = .false.", "ib_reduction_extent = 0", "ib_reduced_order = 2",
        "[linear-algebra]", "smoothing = gauss-seidel", "iterations_init = 3", "iterations_coarse = 10",
        "iterations_fine = 3", "iterations = 10", "tolerance = 1.e-20",
        "[fdv]", "fdv_scheme = fd", "fdv_order = 2", "[solids]", "number = 0", "n_eikonal = 0",
        "[slices]", "slices_number = 0",
        "[physics]", f"physical_model = {'mhd-ideal' if mhd else 'euler'}", "gamma = 1.4",
        "[numerics]", "scheme_space = weno", "reconstruction_variables = characteristic", "reflux = .true.",
    ]
    if mhd:
        lines += ["[mhd]", "divergence_control = eglm", "divb_tol = 0.0", "divb_error = .false.", "rho_floor = 0.0",
                  "p_floor = 0.0", "glm_ch = 3.0", "glm_alpha = 0.18", "glm_damping_length = 1.0",
                  "glm_ch_check = warning"]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        axis = face[0]
        if axis == "z" and not is3d:
            bc = "extrapolation"
        elif axis == "x" and not tilted:
            bc = "periodic"
        else:
            bc = "hydrostatic"
        lines += [f"[bc_{face}]", f"type = {bc}"]
    lines += ["[time]", f"it_max = {steps}", "time_max = 1.0e30", f"CFL = {cfl!r}",
              "[diagnostics]", "ghost_poison = .true.", "conservation_history_save = 1"]
    if case == "froude2d":
        lines += ["[gravity]", "froude = 1.0", "axis = -y"]
    else:
        lines += ["[gravity]"] + [f"g_{x} = {v!r}" for x, v in zip("xyz", GRAVITY[case], strict=True) if v != 0.0]
    ic_type = "polytropic-atmosphere" if case == "poly2d" else "isothermal-atmosphere"
    lines += ["[initial_conditions]", f"type = {ic_type}", "amr_iterations = 1"]
    if case == "poly2d":
        lines += ["polytropic_gamma = 1.2"]
    state = {"r": 1.0, "u": 0.0, "v": 0.0, "w": 0.0, "p": 1.0}
    if mhd:
        state |= {"bx": 0.5, "by": 0.7, "bz": 0.2}
    lines += ["[initial_conditions_region_1]"] + [f"{k} = {v!r}" for k, v in state.items()]
    return "\n".join(lines) + "\n"


def main() -> int:
    """Parse the command line and write the case."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--cfl", type=float, default=1.0e-30)
    args = parser.parse_args()
    args.out.write_text(case_ini(args.case, args.steps, args.cfl))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
