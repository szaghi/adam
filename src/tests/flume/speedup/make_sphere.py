#!/usr/bin/env python3
"""Write the 3-D shock diffraction over a sphere, the second leg of the FLUME speedup matrix (issue #133, stdlib only).

The flow is the 3-D twin of the V6 verification (verification/shock-cylinder): a planar shock of Mach --mach (2 by
default) travelling along +x enters the unit cube from x = --shock-x (0.2) and diffracts over a sphere of radius
--radius (0.13) at the centre of the cube. Gas gamma = 1.4 (cp, cv of the V6 input), pre-shock rho = 1.4, p = 1 (sound
speed 1), the post-shock state from the Rankine-Hugoniot relations (Mach 2: rho 3.7333, u 1.25, p 4.5), which is also
the inflow at x = 0; outflow (extrapolation) at x = 1, inviscid walls on the four lateral faces. The sphere is an
immersed boundary of the `analytical_sphere` kind, an Euler wall (`bc_type = 2`); the rectangle kind is not used (its
distance function is wrong, issue #143). WENO-5 characteristic flux splitting (--weno), SSP-33, CFL 0.5, reflux on.

The grid is an octree of cubic blocks of --block cells (24 by default). The uniform grid is level --base-level (the
uniform case of the matrix sets it to the finest level); with --runtime the realm regrids every --frequency steps by
the density-gradient marker (`[amr_marker_1] mode = 2`) from --base-level up to --max-level: a block where
max |grad rho| > tol takes the finest spacing h(max), every other block coarsens to the base spacing h(base), one
level per regrid; tol = --tol-factor * drho / h(base), drho the jump of the incident shock (0.1 drho / h(base) in the
matrix, the tolerance of the removed matrix, #133). The spacings passed to the marker carry a relative margin of 1e-9,
so that a block on the admissible level is neither refined nor coarsened by round-off.

The initial condition is the two-state Riemann problem (shock at x = --shock-x) on the initial uniform grid; the matrix
never cold-starts the runtime case (it restarts from the warm-up dump), so no initial AMR iteration is asked.

Usage:
    make_sphere.py <out.ini> [--block N] [--base-level LB] [--max-level LF] [--runtime] [--frequency F]
                   [--tol-factor C] [--mach M] [--shock-x X] [--radius R] [--weno 5|7|9] [--cfl C] [--time-max T]
                   [--it-save N] [--restart-save N] [--residuals-save N] [--history N]
"""

from __future__ import annotations

import argparse
from pathlib import Path

GAMMA = 1.4
CP, CV = 1040.004, 742.86  # the V6 input's gas, gamma = 1.4
PRE = {"r": 1.4, "u": 0.0, "v": 0.0, "w": 0.0, "p": 1.0}  # sound speed 1
MARGIN = 1.0 + 1.0e-9  # relative margin on the marker spacings


def post_shock(mach: float) -> dict[str, float]:
    """Return the state behind a shock of Mach `mach` moving along +x into `PRE` (Rankine-Hugoniot relations).

    Parameters
    ----------
    mach : float
        Shock Mach number relative to the gas at rest ahead of it (> 1).

    Returns
    -------
    dict[str, float]
        The post-shock primitive state: r, u, v, w, p.
    """
    if mach <= 1.0:
        raise SystemExit(f"make_sphere: --mach {mach} is not a shock (needs > 1)")
    a = (GAMMA * PRE["p"] / PRE["r"]) ** 0.5
    m2 = mach * mach
    ratio = (GAMMA + 1.0) * m2 / ((GAMMA - 1.0) * m2 + 2.0)
    return {"r": PRE["r"] * ratio, "u": mach * a * (1.0 - 1.0 / ratio), "v": 0.0, "w": 0.0,
            "p": PRE["p"] * (1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (m2 - 1.0))}


def sphere_ini(args: argparse.Namespace) -> str:
    """Return the input of the sphere case.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed command line.

    Returns
    -------
    str
        The INI text.
    """
    s = {5: 3, 7: 4, 9: 5}[args.weno]
    ngc = max(3, s)
    if args.block % 2 or args.block < 2 * ngc:
        raise SystemExit(f"make_sphere: --block {args.block}: a refined axis needs an even count >= {2 * ngc}")
    lb = args.base_level
    lf = args.max_level if args.runtime else lb
    if lf < lb:
        raise SystemExit("make_sphere: --max-level below --base-level")
    post = post_shock(args.mach)
    h_base, h_fine = 1.0 / (args.block * 2**lb), 1.0 / (args.block * 2**lf)
    drho = post["r"] - PRE["r"]
    tol = args.tol_factor * drho / h_base
    lines = [
        f"; FLUME speedup matrix, leg 2 (make_sphere.py): Mach {args.mach} shock over a sphere of radius {args.radius},",
        f"; immersed boundary, octree of {args.block}^3-cell blocks, levels {lb}..{lf}: h {h_base:.6g} .. {h_fine:.6g}"
        + (f", gradient marker tol {tol:.6g} every {args.frequency} steps" if args.runtime else ", uniform"),
        "[IO]", "output_basename        = sphere", f"it_save                = {args.it_save}",
        "restart                = .false.", "restart_basename       = sphere-restart",
        f"restart_save           = {args.restart_save}", f"residuals_save         = {args.residuals_save}",
        "save_residual_fields   = .false.", "save_auxiliary_fields  = .false.",
        "[grid]", f"ni     = {args.block}", f"nj     = {args.block}", f"nk     = {args.block}", f"ngc    = {ngc}",
        "emin_x = 0.0", "emin_y = 0.0", "emin_z = 0.0", "emax_x = 1.0", "emax_y = 1.0", "emax_z = 1.0",
        "null_x = .false.", "null_y = .false.", "null_z = .false.",
        "[amr]", f"max_level      = {lf}", "ratio          = 8", f"iu_ref_levels  = {lb}", "l_prune        = -1",
        f"frequency      = {args.frequency if args.runtime else 0}", "iters          = 1",
        f"markers_number = {1 if args.runtime else 0}",
    ]
    if args.runtime:
        lines += ["regrid_prolongation = conservative",
                  "[amr_marker_1]", "mode         = 2             ; density gradient",
                  "delta_type   = max", f"delta_fine   = {h_fine * MARGIN!r}   ; finest spacing where |grad rho| > tol",
                  f"delta_coarse = {h_base * MARGIN!r}   ; base spacing elsewhere", "field        = 1",
                  "var          = 1", f"tol          = {tol!r}   ; {args.tol_factor} * drho / h(base), drho = {drho:.6g}"]
    lines += [
        "[runge_kutta]", "scheme = runge-kutta-ssp-33",
        "[weno]", f"scheme              = weno-u-{args.weno}",
        "[solids]", "number    = 1", "n_eikonal = 2",
        "[solid_1]", "name            = sphere", "bc_type         = 2                 ; Euler wall",
        "definition      = analytical_sphere", "sphere_center_x = 0.5", "sphere_center_y = 0.5",
        "sphere_center_z = 0.5", f"sphere_radius   = {args.radius!r}",
        "[slices]", "slices_number = 0",
        "[physics]", "physical_model = euler", f"cp             = {CP!r}", f"cv             = {CV!r}",
        "[numerics]", "scheme_space             = weno", "reconstruction_variables = characteristic",
        "reflux                   = .true.",
        "[bc_x_min]", "type = inflow", *(f"{k} = {v!r}" for k, v in post.items()),
        "[bc_x_max]", "type = extrapolation",
    ]
    for face in ("y_min", "y_max", "z_min", "z_max"):
        lines += [f"[bc_{face}]", "type = wall-inviscid"]
    lines += ["[initial_conditions]", "type           = riemann-problem", "regions_number = 2",
              "amr_iterations = 0"]
    for region, (state, x0, x1) in enumerate(((post, 0.0, args.shock_x), (PRE, args.shock_x, 1.0)), start=1):
        lines += [f"[initial_conditions_region_{region}]", *(f"{k:6s} = {v!r}" for k, v in state.items()),
                  f"emin_x = {x0!r}", f"emax_x = {x1!r}", "emin_y = 0.0", "emax_y = 1.0", "emin_z = 0.0",
                  "emax_z = 1.0"]
    lines += ["[time]", "it_max   = -1", f"time_max = {args.time_max}", f"CFL      = {args.cfl}",
              "[diagnostics]", f"conservation_history_save = {args.history}"]
    return "\n".join(lines) + "\n"


def main() -> None:
    """Parse the command line and write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path, help="the input file to write")
    parser.add_argument("--block", type=int, default=24, help="cells per block along each axis")
    parser.add_argument("--base-level", type=int, default=2, help="level of the initial uniform grid")
    parser.add_argument("--max-level", type=int, default=4, help="finest level (with --runtime)")
    parser.add_argument("--runtime", action="store_true", help="runtime AMR, density-gradient marker")
    parser.add_argument("--frequency", type=int, default=10, help="regrid period in steps")
    parser.add_argument("--tol-factor", type=float, default=0.1, help="marker tol = this * drho / h(base)")
    parser.add_argument("--mach", type=float, default=2.0, help="shock Mach number")
    parser.add_argument("--shock-x", type=float, default=0.2, help="initial shock position")
    parser.add_argument("--radius", type=float, default=0.13, help="sphere radius (centre of the unit cube)")
    parser.add_argument("--weno", type=int, choices=(5, 7, 9), default=5, help="WENO order")
    parser.add_argument("--cfl", default="0.5")
    parser.add_argument("--time-max", default="0.25")
    parser.add_argument("--it-save", default="1000000", help="checkpoint period (the first and last are saved)")
    parser.add_argument("--restart-save", default="0", help="restart dump period (0: none)")
    parser.add_argument("--residuals-save", default="1", help="residuals period (0: none)")
    parser.add_argument("--history", default="1", help="conservation history period (0: none)")
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(sphere_ini(args))


if __name__ == "__main__":
    main()
