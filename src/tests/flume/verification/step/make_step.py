#!/usr/bin/env python3
"""Write the Woodward-Colella Mach 3 forward-facing step as a three-realm FLUME forest (issue #46).

Woodward & Colella (1984, J. Comput. Phys. 54, 115-173, section IV b): a Mach 3 wind tunnel [0, 3] x [0, 1] with a step
of height 0.2 at x = 0.6; the gas (gamma = 1.4) enters from the left with rho = 1.4, u = 3, p = 1 (sound speed 1), and
leaves on the right; every other boundary, the step faces included, is a reflecting wall. The reference results are
taken at t = 4.

The tunnel is L-shaped, three rectangles glued by 1:1 mirror seams, so the step is body-fitted: its faces are physical
walls of the realms, with no immersed boundary.

    y = 1   +--------+-----------------------------+   wall
            |   B    |              C              |
     inflow |        | seam                        | outflow
    y = 0.2 +--------+   +------ wall (step top) --+
            |   A    | wall (step front)
    y = 0   +--------+
          x = 0    x = 0.6                       x = 3

A seam A-B at y = 0.2 (x in [0, 0.6]), a seam B-C at x = 0.6 (y in [0.2, 1]). All three realms share the base cell
1/N (--cells N, 80 by default, the coarse grid of Woodward and Colella) and N/20 blocks per axis (iu_ref_levels =
log2(N/20): 4x4 blocks at N = 80, 2x2 at N = 40; the step realm A keeps 4 cells per block along y). --refine adds one
level on the blocks whose centroid lies in two boxes kept off the seams (a refined block may not touch an inter-realm
seam): C on x in [1.3, 3], y in [0.2, 0.55] (2:1 seams at x = 1.2 and y = 0.6 at N = 80, x = 1.8 and y = 0.6 at
N = 40, crossed by the corner's expansion and the reflected shocks) and B on x in [0, 0.4], y in [0.55, 1].

--tree quad builds every realm as a quadtree with nk = 1 (true 2-D, issue #46); --tree oct as an octree with a null z
axis and nk = 4, the layout of the AMR verifications before #46. The z extent is 1 in both (null axis).

The scheme is WENO-5 characteristic flux splitting with SSP-54 (no positivity limiter: it is refused on multi-realm
runs). The corner is left untreated (Woodward and Colella fix the entropy near it): its numerical boundary layer is
part of the reference pictures of every scheme without the fix.

Usage:
    make_step.py <out-dir> [--tree quad|oct] [--cells N] [--time-max T] [--refine] [--it-save N] [--cfl C]
"""

from __future__ import annotations

import argparse
from pathlib import Path

GAMMA_CP, GAMMA_CV = 3.5, 2.5  # gamma = 1.4 (R = 1: only the ratio matters)
INFLOW = {"r": 1.4, "u": 3.0, "v": 0.0, "w": 0.0, "p": 1.0}
# realm: (xmin, xmax, ymin, ymax, faces {face: bc}, refine box (xmin, xmax, ymin, ymax) or None)
REALMS = {
    "A": (0.0, 0.6, 0.0, 0.2, {"x_min": "inflow", "x_max": "wall-inviscid", "y_min": "wall-inviscid",
                               "y_max": "seam"}, None),
    "B": (0.0, 0.6, 0.2, 1.0, {"x_min": "inflow", "x_max": "seam", "y_min": "seam", "y_max": "wall-inviscid"},
          (0.0, 0.4, 0.55, 1.0)),
    "C": (0.6, 3.0, 0.2, 1.0, {"x_min": "seam", "x_max": "extrapolation", "y_min": "wall-inviscid",
                               "y_max": "wall-inviscid"}, (1.3, 3.0, 0.2, 0.55)),
}
SEAMS = (("A", "+y", "B", "-y"), ("B", "+x", "C", "-x"))


def realm_ini(name: str, tree: str, cells: int, refine: bool, time_max: str, it_save: str, cfl: str) -> str:
    """Return the INI of one realm."""
    xmin, xmax, ymin, ymax, faces, box = REALMS[name]
    if cells % 20 or (cells // 20) & (cells // 20 - 1):
        raise SystemExit(f"make_step: --cells {cells} is not 20 times a power of 2")
    levels = (cells // 20).bit_length() - 1  # N / 20 blocks per axis
    nx, ny = round((xmax - xmin) * cells), round((ymax - ymin) * cells)
    ni, nj = nx // 2**levels, ny // 2**levels
    if ni * 2**levels != nx or nj * 2**levels != ny or ni % 2 or nj % 2 or min(ni, nj) < 4:
        raise SystemExit(f"make_step: --cells {cells} gives realm {name} blocks of {ni} x {nj} cells (need even, >= 4)")
    refined = refine and box is not None
    lines = [
        f"; Woodward-Colella step, realm {name} (make_step.py, --tree {tree}, --cells {cells})",
        "[IO]", f"output_basename        = step-{name}", f"it_save                = {it_save}",
        "restart                = .false.", f"restart_basename       = step-{name}-restart",
        "restart_save           = 0",
        "residuals_save         = 1", "save_memory_status     = .false.", "save_residual_fields   = .false.",
        "save_auxiliary_fields  = .false.", "save_curl_fields       = .false.", "save_divergence_fields = .false.",
        "save_gradient_fields   = .false.", "save_laplacian_fields  = .false.",
        "[grid]", f"ni     = {ni}", f"nj     = {nj}", f"nk     = {1 if tree == 'quad' else 4}", "ngc    = 3",
        f"emin_x = {xmin!r}", f"emin_y = {ymin!r}", "emin_z = 0.0", f"emax_x = {xmax!r}", f"emax_y = {ymax!r}",
        "emax_z = 1.0", "null_x = .false.", "null_y = .false.", "null_z = .true.",
        "[amr]", f"max_level      = {levels + (1 if refined else 0)}", f"ratio          = {4 if tree == 'quad' else 8}",
        f"iu_ref_levels  = {levels}", "i_prune        = 0", "j_prune        = 0", "k_prune        = 0",
        "l_prune        = -1", "frequency      = 999999", "iters          = 1",
        f"markers_number = {1 if refined else 0}",
    ]
    if refined:
        lines += ["[amr_marker_1]", "mode         = 1", "geo_type     = primitive-box", "delta_type   = max",
                  "delta_fine   = 0.0", "delta_coarse = 0.0", f"box_xmin     = {box[0]!r}",
                  f"box_ymin     = {box[2]!r}",
                  "box_zmin     = -1.0e30", f"box_xmax     = {box[1]!r}", f"box_ymax     = {box[3]!r}",
                  "box_zmax     = 1.0e30", f"target_level = {levels + 1}"]
    lines += [
        "[field]", "nv = 5", "[runge_kutta]", "scheme = runge-kutta-ssp-54",
        "[weno]", "scheme              = weno-u-5", "ror_number          = 0", "ror_threshold       = 0.9",
        "ror_vars_number     = 0", "enable_ror_stats    = .false.", "ib_reduction_extent = 0",
        "ib_reduced_order    = 2",
        "[linear-algebra]", "smoothing         = gauss-seidel", "iterations_init   = 3", "iterations_coarse = 10",
        "iterations_fine   = 3", "iterations        = 10", "tolerance         = 1.e-20",
        "[fdv]", "fdv_scheme = fd", "fdv_order  = 2", "[solids]", "number    = 0", "n_eikonal = 0",
        "[slices]", "slices_number = 0",
        "[physics]", "physical_model = euler", f"cp             = {GAMMA_CP!r}", f"cv             = {GAMMA_CV!r}",
        "[numerics]", "scheme_space             = weno", "reconstruction_variables = characteristic",
        "reflux                   = .true.",
    ]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        bc = faces.get(face, "extrapolation")
        lines += [f"[bc_{face}]", f"type = {'extrapolation' if bc == 'seam' else bc}"]
        if bc == "inflow":
            lines += [f"{k} = {v!r}" for k, v in INFLOW.items()]
    lines += ["[initial_conditions]", "type           = uniform", "s              = 0.0",
              f"amr_iterations = {1 if refined else 0}", "[initial_conditions_region_1]"]
    lines += [f"{k:6s} = {v!r}" for k, v in INFLOW.items()]
    lines += ["[time]", "it_max   = -1", f"time_max = {time_max}", f"CFL      = {cfl}",
              "[diagnostics]", "conservation_history_save = 1"]
    return "\n".join(lines) + "\n"


def manifest() -> str:
    """Return the forest manifest."""
    lines = ["; Woodward-Colella step: three realms glued by mirror seams (make_step.py)", "[forest]",
             f"realms_number = {len(REALMS)}"]
    for i, name in enumerate(REALMS, start=1):
        lines += [f"[realm.{i}]", f"ini = step-{name}.ini"]
    lines += ["[forest.topology]", f"inter_realm_faces_number = {len(SEAMS)}"]
    index = {name: i for i, name in enumerate(REALMS, start=1)}
    for f, (a, fa, b, fb) in enumerate(SEAMS, start=1):
        lines += [f"[forest.topology.face_{f}]", f"realm_a          = {index[a]}", f"face_a           = {fa}",
                  f"realm_b          = {index[b]}", f"face_b           = {fb}", "coupling         = mirror",
                  "coupling_cadence = stage_coincident"]
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write the manifest and the realm inputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path, help="output directory")
    parser.add_argument("--tree", choices=("quad", "oct"), default="quad")
    parser.add_argument("--cells", type=int, default=80, help="base cells per unit length")
    parser.add_argument("--time-max", default="4.0")
    parser.add_argument("--refine", action="store_true", help="refine one level on two boxes off the seams")
    parser.add_argument("--it-save", default="1000000", help="checkpoint period (the first and last are always saved)")
    parser.add_argument("--cfl", default="0.5")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "step.ini").write_text(manifest())
    for name in REALMS:
        (args.out / f"step-{name}.ini").write_text(
            realm_ini(name, args.tree, args.cells, args.refine, args.time_max, args.it_save, args.cfl))


if __name__ == "__main__":
    main()
