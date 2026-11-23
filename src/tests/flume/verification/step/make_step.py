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
1/N (--cells N, 80 by default, the coarse grid of Woodward and Colella); B and C have N/10 blocks per axis (8x8 at
N = 80), narrow enough to refine close to the seams, A N/40 (2x2 at N = 80, a single block at N = 40), so that its
blocks keep 8 cells along y: a refined block needs at least 2 ngc = 6 cells along every active axis (realm_object
check_amr_block_cells). The seams A-B and B-C are misaligned 1:1 mirror seams (issue #51).

--refine adds one level on the blocks whose centroid lies in a box of each realm. The refinement is STATIC and
GEOMETRIC: the boxes are placed by hand on the features of the reference solution at t = 4 and the grid is frozen
after initialisation (this legacy layout predates runtime AMR, M5, and its mirror seams keep the blocks next to them
on the base level: the strongest features, the Mach stem and the triple point at x ~ 0.6 on the B-C seam, stay on
the base grid). The boxes, at N = 80:
  A  x in [0, 0.6], y in [0, 0.1]: the bow shock standing on the lower wall at x ~ 0.3;
  B  x in [0.3, 0.525], y in [0.3, 0.7]: the bow shock bending up towards the Mach stem (the column x > 0.525
     touches the seam);
  C  x in [0.9, 3], y in [0.2, 1]: everything but the column on the seam, i.e. the reflected shock reaching the step
     top at x ~ 1.2, the slip line along y ~ 0.8 and the re-reflection reaching the top wall at x ~ 2.4.
At N = 40 the same boxes select nothing in A (its single block touches the seam); B x in [0.3, 0.45], y in
[0.4, 0.8]; C x in [1.2, 3].

--tree quad builds every realm as a quadtree with nk = 1 (true 2-D, issue #46); --tree oct as an octree with a null z
axis, the layout of the AMR verifications before #46: nk = 4 in B and C and 16 in A (a block level halves the
octree's z blocks too, the z cells must match across the mirror seams, and the 2:1 seam fill needs nk >= 4 along a
refined axis). The z extent is 1 in both (null axis).

The scheme is WENO-5 characteristic flux splitting with SSP-54 (no positivity limiter: it is refused on multi-realm
runs). The corner is left untreated (Woodward and Colella fix the entropy near it): its numerical boundary layer is
part of the reference pictures of every scheme without the fix.

--layout aligned (issue #133 leg 1, #142's 3-D step) lines the blocks up across the seams so that they can be
`coupling = auto` seams, along which the runtime grid refines and coarsens freely (one level per seam side, equal
across, forest_object balance_seam_marks): `auto` needs the same block cell counts along a seam
(forest_object populate_inter_realm_topology). With one root block per realm and square cells the block shapes are
fixed by the realms' aspect ratios: A (3m, m), B (3m, 4m), C (12m, 4m) cells (--m, 16 by default), so that at level L
every realm has the spacing h = 0.2/(m 2^L). --dims 3 extrudes the tunnel in z on octrees: nk cells per block (--nk,
10 by default), z extent Lz = 0.2 nk/m (cubic cells), reflecting walls on the z faces (the flow stays z-invariant: the
3-D case costs, it adds no physics). --weno 5|7|9 picks the scheme (ngc = max(3, S)); a refined axis needs an even
count of at least 2 ngc cells (m and nk >= 10 at WENO-9). --base-level Lb is the uniform level of the initial grid,
--max-level Lf the finest; --runtime regrids every --frequency steps by a Loehner marker on the density (refine 0.3,
derefine 0.1, buffer ngc - 1), from Lb up to Lf; without --runtime the grid is uniform at Lb (Lf ignored).
The legacy layout (default) keeps the EV-step legs and their figures unchanged.

Usage:
    make_step.py <out-dir> [--tree quad|oct] [--cells N] [--time-max T] [--refine] [--it-save N] [--cfl C]
    make_step.py <out-dir> --layout aligned [--dims 2|3] [--m M] [--nk NK] [--weno 5|7|9] [--base-level LB]
                 [--max-level LF] [--runtime] [--frequency F] [--time-max T] [--it-save N] [--cfl C]
                 [--coupling auto|mirror] [--restart-save N] [--residuals-save N]
"""

from __future__ import annotations

import argparse
from pathlib import Path

GAMMA_CP, GAMMA_CV = 3.5, 2.5  # gamma = 1.4 (R = 1: only the ratio matters)
INFLOW = {"r": 1.4, "u": 3.0, "v": 0.0, "w": 0.0, "p": 1.0}
# realm: (xmin, xmax, ymin, ymax, faces {face: bc}, refine box on block centroids (xmin, xmax, ymin, ymax), extra
# block levels beyond log2(N/20))
REALMS = {
    "A": (0.0, 0.6, 0.0, 0.2, {"x_min": "inflow", "x_max": "wall-inviscid", "y_min": "wall-inviscid",
                               "y_max": "seam"}, (0.1, 0.5, 0.0, 0.09), -1),
    "B": (0.0, 0.6, 0.2, 1.0, {"x_min": "inflow", "x_max": "seam", "y_min": "seam", "y_max": "wall-inviscid"},
          (0.3, 0.52, 0.31, 0.7), 1),
    "C": (0.6, 3.0, 0.2, 1.0, {"x_min": "seam", "x_max": "extrapolation", "y_min": "wall-inviscid",
                               "y_max": "wall-inviscid"}, (0.95, 3.0, 0.2, 1.0), 1),
}
SEAMS = (("A", "+y", "B", "-y"), ("B", "+x", "C", "-x"))
EXTRA_MAX = max(realm[-1] for realm in REALMS.values())


def realm_ini(name: str, tree: str, cells: int, refine: bool, time_max: str, it_save: str, cfl: str) -> str:
    """Return the INI of one realm."""
    xmin, xmax, ymin, ymax, faces, box, extra = REALMS[name]
    if cells % 20 or (cells // 20) & (cells // 20 - 1):
        raise SystemExit(f"make_step: --cells {cells} is not 20 times a power of 2")
    levels = max(0, (cells // 20).bit_length() - 1 + extra)  # N / 20 blocks per axis, times 2**extra
    nx, ny = round((xmax - xmin) * cells), round((ymax - ymin) * cells)
    ni, nj = nx // 2**levels, ny // 2**levels
    refined = refine and box is not None
    least = 6 if refined else 4  # 2 ngc along a refined axis
    if ni * 2**levels != nx or nj * 2**levels != ny or ni % 2 or nj % 2 or min(ni, nj) < least:
        raise SystemExit(f"make_step: --cells {cells} gives realm {name} blocks of {ni} x {nj} cells (need even, "
                         f">= {least})")
    lines = [
        f"; Woodward-Colella step, realm {name} (make_step.py, --tree {tree}, --cells {cells})",
        "[IO]", f"output_basename        = step-{name}", f"it_save                = {it_save}",
        "restart                = .false.", f"restart_basename       = step-{name}-restart",
        "restart_save           = 0",
        "residuals_save         = 1","save_residual_fields   = .false.",
        "save_auxiliary_fields  = .false.",
        "[grid]", f"ni     = {ni}", f"nj     = {nj}",
        f"nk     = {1 if tree == 'quad' else 4 * 2 ** (EXTRA_MAX - extra)}", "ngc    = 3",
        f"emin_x = {xmin!r}", f"emin_y = {ymin!r}", "emin_z = 0.0", f"emax_x = {xmax!r}", f"emax_y = {ymax!r}",
        "emax_z = 1.0", "null_x = .false.", "null_y = .false.", "null_z = .true.",
        "[amr]", f"max_level      = {levels + (1 if refined else 0)}", f"ratio          = {4 if tree == 'quad' else 8}",
        f"iu_ref_levels  = {levels}", "l_prune        = -1", "frequency      = 0", "iters          = 1",
        f"markers_number = {1 if refined else 0}",
    ]
    if refined:
        lines += ["[amr_marker_1]", "mode         = 1", "geo_type     = primitive-box", "delta_type   = max",
                  "delta_fine   = 0.0", "delta_coarse = 0.0", f"box_xmin     = {box[0]!r}",
                  f"box_ymin     = {box[2]!r}",
                  "box_zmin     = -1.0e30", f"box_xmax     = {box[1]!r}", f"box_ymax     = {box[3]!r}",
                  "box_zmax     = 1.0e30", f"target_level = {levels + 1}"]
    lines += [
        "[runge_kutta]", "scheme = runge-kutta-ssp-54",
        "[weno]", "scheme              = weno-u-5", 
        
        
        "[solids]", "number    = 0",
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


# aligned layout: realm -> (block cells along x, along y) in units of m
ALIGNED = {"A": (3, 1), "B": (3, 4), "C": (12, 4)}


def aligned_realm_ini(name: str, args: argparse.Namespace) -> str:
    """Return the INI of one realm of the aligned layout (--layout aligned)."""
    xmin, xmax, ymin, ymax, faces, _, _ = REALMS[name]
    s = {5: 3, 7: 4, 9: 5}[args.weno]
    ngc = max(3, s)
    least = 2 * ngc
    if args.m % 2 or args.m < least or (args.dims == 3 and (args.nk % 2 or args.nk < least)):
        raise SystemExit(f"make_step: --m {args.m} / --nk {args.nk}: a refined axis needs an even count >= {least}")
    lb = args.base_level
    lf = args.max_level if args.runtime else lb
    if lf < lb:
        raise SystemExit("make_step: --max-level below --base-level")
    ni, nj = (k * args.m for k in ALIGNED[name])
    nk = args.nk if args.dims == 3 else 1
    lz = 0.2 * nk / args.m if args.dims == 3 else 1.0
    h = 0.2 / (args.m * 2**lb)
    lines = [
        f"; Woodward-Colella step, realm {name} (make_step.py --layout aligned, --dims {args.dims}, --m {args.m}, "
        f"--weno {args.weno}, levels {lb}..{lf}: h {h:.6g} .. {0.2 / (args.m * 2**lf):.6g})",
        "[IO]", f"output_basename        = step-{name}", f"it_save                = {args.it_save}",
        "restart                = .false.", f"restart_basename       = step-{name}-restart",
        f"restart_save           = {args.restart_save}",
        f"residuals_save         = {args.residuals_save}", "save_residual_fields   = .false.",
        "save_auxiliary_fields  = .false.",
        "[grid]", f"ni     = {ni}", f"nj     = {nj}", f"nk     = {nk}", f"ngc    = {ngc}",
        f"emin_x = {xmin!r}", f"emin_y = {ymin!r}", "emin_z = 0.0", f"emax_x = {xmax!r}", f"emax_y = {ymax!r}",
        f"emax_z = {lz!r}", "null_x = .false.", "null_y = .false.",
        f"null_z = {'.true.' if args.dims == 2 else '.false.'}",
        "[amr]", f"max_level      = {lf}", f"ratio          = {4 if args.dims == 2 else 8}",
        f"iu_ref_levels  = {lb}", "l_prune        = -1", f"frequency      = {args.frequency if args.runtime else 0}",
        "iters          = 1", f"markers_number = {1 if args.runtime else 0}",
    ]
    if args.runtime:
        lines += ["regrid_prolongation = conservative",
                  "[amr_marker_1]", "mode         = 4", "delta_type   = max", "delta_fine   = 0.0",
                  "delta_coarse = 0.0", "field        = 1", "var          = 1", "refine_tol   = 0.3",
                  "derefine_tol = 0.1", "floor        = 0.0", f"buffer       = {ngc - 1}"]
    lines += [
        "[runge_kutta]", "scheme = runge-kutta-ssp-54",
        "[weno]", f"scheme              = weno-u-{args.weno}",
        "[solids]", "number    = 0",
        "[slices]", "slices_number = 0",
        "[physics]", "physical_model = euler", f"cp             = {GAMMA_CP!r}", f"cv             = {GAMMA_CV!r}",
        "[numerics]", "scheme_space             = weno", "reconstruction_variables = characteristic",
        "reflux                   = .true.",
    ]
    for face in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max"):
        bc = faces.get(face, "wall-inviscid" if args.dims == 3 else "extrapolation")
        lines += [f"[bc_{face}]", f"type = {'extrapolation' if bc == 'seam' else bc}"]
        if bc == "inflow":
            lines += [f"{k} = {v!r}" for k, v in INFLOW.items()]
    lines += ["[initial_conditions]", "type           = uniform", "s              = 0.0",
              "amr_iterations = 0", "[initial_conditions_region_1]"]
    lines += [f"{k:6s} = {v!r}" for k, v in INFLOW.items()]
    lines += ["[time]", "it_max   = -1", f"time_max = {args.time_max}", f"CFL      = {args.cfl}",
              "[diagnostics]", "conservation_history_save = 1"]
    return "\n".join(lines) + "\n"


def manifest(coupling: str = "mirror") -> str:
    """Return the forest manifest."""
    lines = [f"; Woodward-Colella step: three realms glued by {coupling} seams (make_step.py)", "[forest]",
             f"realms_number = {len(REALMS)}"]
    for i, name in enumerate(REALMS, start=1):
        lines += [f"[realm.{i}]", f"ini = step-{name}.ini"]
    lines += ["[forest.topology]", f"inter_realm_faces_number = {len(SEAMS)}"]
    index = {name: i for i, name in enumerate(REALMS, start=1)}
    for f, (a, fa, b, fb) in enumerate(SEAMS, start=1):
        lines += [f"[forest.topology.face_{f}]", f"realm_a          = {index[a]}", f"face_a           = {fa}",
                  f"realm_b          = {index[b]}", f"face_b           = {fb}", f"coupling         = {coupling}",
                  "coupling_cadence = stage_coincident"]
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write the manifest and the realm inputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path, help="output directory")
    parser.add_argument("--tree", choices=("quad", "oct"), default="quad")
    parser.add_argument("--cells", type=int, default=80, help="base cells per unit length")
    parser.add_argument("--time-max", default="4.0")
    parser.add_argument("--refine", action="store_true", help="refine one level on a static box of each realm")
    parser.add_argument("--it-save", default="1000000", help="checkpoint period (the first and last are always saved)")
    parser.add_argument("--cfl", default="0.5")
    parser.add_argument("--layout", choices=("legacy", "aligned"), default="legacy")
    parser.add_argument("--dims", type=int, choices=(2, 3), default=2, help="aligned: 2-D quadtrees or 3-D octrees")
    parser.add_argument("--m", type=int, default=16, help="aligned: block cells along A's y")
    parser.add_argument("--nk", type=int, default=10, help="aligned, 3-D: block cells along z")
    parser.add_argument("--weno", type=int, choices=(5, 7, 9), default=5, help="aligned: WENO order")
    parser.add_argument("--base-level", type=int, default=1, help="aligned: level of the initial uniform grid")
    parser.add_argument("--max-level", type=int, default=1, help="aligned: finest level (with --runtime)")
    parser.add_argument("--runtime", action="store_true", help="aligned: runtime AMR, Loehner marker on density")
    parser.add_argument("--frequency", type=int, default=10, help="aligned: regrid period in steps")
    parser.add_argument("--coupling", choices=("auto", "mirror"), default=None,
                        help="seam coupling (aligned: auto by default; legacy: mirror)")
    parser.add_argument("--restart-save", default="0", help="aligned: restart dump period (0: none)")
    parser.add_argument("--residuals-save", default="1", help="aligned: residuals period")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.layout == "aligned":
        (args.out / "step.ini").write_text(manifest(args.coupling or "auto"))
        for name in REALMS:
            (args.out / f"step-{name}.ini").write_text(aligned_realm_ini(name, args))
        return
    if args.coupling not in (None, "mirror"):
        raise SystemExit("make_step: the legacy layout's blocks do not line up: its seams are mirror seams")
    (args.out / "step.ini").write_text(manifest())
    for name in REALMS:
        (args.out / f"step-{name}.ini").write_text(
            realm_ini(name, args.tree, args.cells, args.refine, args.time_max, args.it_save, args.cfl))


if __name__ == "__main__":
    main()
