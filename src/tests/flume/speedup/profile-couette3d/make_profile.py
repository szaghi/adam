#!/usr/bin/env python3
"""Write the input of the FLUME profiling case: a 3-D compressible Couette channel with refined walls (stdlib only).

Why this case: a profile is only useful if the run is cheap, reproducible and goes through the code that production
runs spend their time in. The exact compressible Couette flow (regression/couette, issue #65 VV-3) extended to 3-D is
steady, so every profile of it sees the same work, and it exercises in one short run:

- the convective residual: WENO-5 flux splitting on characteristic variables, along x, y and z (octree, nothing null);
- the dissipative residual: viscous stress and heat flux, 4th order, of the same magnitude as the convective time
  step bound at these sizes (both enter dt);
- the time integration: SSP Runge-Kutta 5 stages (stage buffers, stage ghost exchanges);
- the boundary conditions: an isothermal no-slip wall (y = 0), a moving adiabatic no-slip wall (y = 1), periodic x
  and z (edges and corners composed);
- the ghost exchange: intra-rank and MPI, every stage;
- AMR: the base level plus one level on the block layer along each wall (two primitive-box markers), so every stage
  fills 2:1 seam ghosts (restriction and prolongation) and the reflux corrects the coarse faces; a runtime regrid
  every 10 steps (marking, ghost refresh, the library regrid, on FNL the host round trip), which leaves the grid
  unchanged (the box markers vote "keep");
- the reductions: the global time step and the conservation integrals, every step.

Not exercised: the immersed boundary, the prolongation of a regrid that changes the grid, inter-realm seams (but see
--realms 2), MHD, gravity, the positivity limiter.

Sizes (blocks = 8^L + 14 4^L: the base level L, its two wall layers refined once; memory ~65 block-sized arrays):
    small   NI 16, base L2:  288 blocks,  1.2 M cells, ~1.6 GB   WSL (CPU np 2, one GPU), one A30
    medium  NI 16, base L3: 1408 blocks,  5.8 M cells, ~8 GB     medooza CPU node (np 64), one A30, 4 A30
    large   NI 32, base L3: 1408 blocks, 46 M cells,   ~40 GB    one AMD 192 GB GPU, 4 A30 (~10 GB each)

--uniform writes the simplest variant instead: one realm on a uniform grid, no AMR marker, no regrid, no 2:1 seam
(convective and dissipative residuals, the stages, the boundaries, the ghost exchange and the reductions only):
    small   NI 16, level 3:  512 blocks,  2.1 M cells, ~2.8 GB   WSL
    medium  NI 32, level 3:  512 blocks, 16.8 M cells, ~15 GB    one A30, the medooza CPU node
    large   NI 32, level 4: 4096 blocks,  134 M cells, ~117 GB   one AMD 192 GB GPU

--realms 2 writes the refined-walls variant as a two-realm forest (issue #110 P6): the channel split at x = 0.5, each
realm NI/2 cells along x (cubic cells, the same cells as one realm), x extrapolated outside (the Couette flow does not
vary along x, so it stays steady; a seam cannot wrap the periodic x axis), one coupling = auto seam at x = 0.5 that
crosses both refined wall layers: a seam with mixed levels (equal cells face by face, finer along the walls). <out.ini>
is then the manifest, next to <stem>-r1.ini and <stem>-r2.ini. The box markers keep the grid, so the topology is never
rebuilt unless the run sets ADAM_FOREST_REBUILD_TOPOLOGY=1 (profile.sh --rebuild): a rebuild after every step.

Usage:
    make_profile.py <small|medium|large> <out.ini> [--uniform | --realms 2] [--steps N] [--frequency N]
                    [--scheme weno-u-<2S-1>]
"""

from __future__ import annotations

import argparse
from pathlib import Path

SIZES = {"small": (16, 2), "medium": (16, 3), "large": (32, 3)}  # NI, base level
SCHEMES = ("weno-u-1", "weno-u-3", "weno-u-5", "weno-u-7", "weno-u-9")  # FLUME's upwind WENO, S = 1 ... 5
UNIFORM_SIZES = {"small": (16, 3), "medium": (32, 3), "large": (32, 4)}  # NI, level


def profile_ini(size: str, steps: int, frequency: int, uniform: bool = False, scheme: str = "weno-u-5",
                realm: tuple[float, float, str] | None = None) -> str:
    """Return the input of one size: the refined-walls variant, or the uniform one, with the WENO scheme `scheme`; with
    `realm` = (xmin, xmax, basename), the realm of a forest split along x, extrapolated at its x faces."""
    half_width = (int(scheme.rsplit("-", 1)[1]) + 1) // 2  # weno-u-<2S-1>
    ngc = max(3, half_width)  # the stencil reaches S cells beyond a face; 3 (WENO-5) is the committed inputs' ngc
    ni, base = (UNIFORM_SIZES if uniform else SIZES)[size]
    layer = 1.0 / 2**base  # the height of a base block: the wall layer is the first block row along y
    xmin, xmax, basename = realm if realm else (0.0, 1.0, "couette3d")
    ni_x = ni if realm is None else int(round(ni * (xmax - xmin)))  # cubic cells in every realm
    if ni_x < 2 * ngc:  # realm_object check_amr_block_cells would stop the run at initialisation
        raise SystemExit(f"make_profile.py: --realms 2 with {scheme} on size {size} gives blocks of {ni_x} cells along "
                         f"x, fewer than 2 ngc = {2 * ngc}: use a smaller scheme or a larger size")
    x_bc = "periodic" if realm is None else "extrapolation"
    if uniform:
        what = f"NI {ni}, uniform at level {base}: one realm, no AMR marker, no regrid; {steps} steps"
    else:
        what = (f"NI {ni}, base level {base}, the wall block layers refined to level {base + 1}; {steps} steps, a "
                f"regrid every {frequency}")
    lines = [
        f"; FLUME profiling case: 3-D Couette channel, size {size}{' uniform' if uniform else ''} (make_profile.py: "
        "rerun it, do not edit)",
        f"; {what}",
        "[IO]", f"output_basename        = {basename}", "it_save                = 0",
        "restart                = .false.",
        f"restart_basename       = {basename}-restart", "restart_save           = 0", "residuals_save         = 0",
        "save_memory_status     = .false.", "save_residual_fields   = .false.", "save_auxiliary_fields  = .false.",
        "save_curl_fields       = .false.", "save_divergence_fields = .false.", "save_gradient_fields   = .false.",
        "save_laplacian_fields  = .false.",
        "", "[grid]", f"ni     = {ni_x}", f"nj     = {ni}", f"nk     = {ni}", f"ngc    = {ngc}",
        f"emin_x = {xmin!r}", "emin_y = 0.0", "emin_z = 0.0", f"emax_x = {xmax!r}", "emax_y = 1.0", "emax_z = 1.0",
        "null_x = .false.", "null_y = .false.", "null_z = .false.",
        "", "[amr]", f"max_level      = {base if uniform else base + 1}", "ratio          = 8",
        f"iu_ref_levels  = {base}", "i_prune        = 0", "j_prune        = 0", "k_prune        = 0",
        "l_prune        = -1", f"frequency      = {0 if uniform else frequency}", "iters          = 1",
        f"markers_number = {0 if uniform else 2}",
    ]
    for n, (lo, hi) in enumerate(() if uniform else ((-1.0e30, layer), (1.0 - layer, 1.0e30)), start=1):
        lines += [
            "", f"[amr_marker_{n}]", "mode         = 1           ; geometric", "geo_type     = primitive-box",
            "delta_type   = max", "delta_fine   = 0.0", "delta_coarse = 0.0",
            "box_xmin     = -1.0e30", f"box_ymin     = {lo!r}", "box_zmin     = -1.0e30",
            "box_xmax     = 1.0e30", f"box_ymax     = {hi!r}", "box_zmax     = 1.0e30", f"target_level = {base + 1}",
        ]
    lines += [
        "", "[field]", "nv = 5", "", "[runge_kutta]", "scheme = runge-kutta-ssp-54",
        "", "[weno]", f"scheme              = {scheme}", "ror_number          = 0", "ror_threshold       = 0.9",
        "ror_vars_number     = 0", "enable_ror_stats    = .false.", "ib_reduction_extent = 0", "ib_reduced_order    = 2",
        "", "[fdv]", "fdv_scheme = fd", "fdv_order  = 2", "", "[solids]", "number    = 0", "n_eikonal = 0",
        "", "[slices]", "slices_number = 0",
        "", "[physics]", "physical_model = euler", "gamma          = 1.4", "viscosity      = 0.01", "prandtl        = 0.75",
        "", "[numerics]", "scheme_space             = weno", "reconstruction_variables = characteristic",
        "reflux                   = .true.", "dissipative_order        = 4",
        "", "[time]", f"it_max   = {steps}", "time_max = 1.0e30", "CFL      = 0.2",
        "", "[diagnostics]", "conservation_history_save = 0",
        "", "[bc_x_min]", f"type = {x_bc}", "[bc_x_max]", f"type = {x_bc}",
        "[bc_y_min]", "type             = wall-isothermal", "wall_temperature = 1.0",
        "[bc_y_max]", "type   = wall-noslip", "wall_u = 1.0",
        "[bc_z_min]", "type = periodic", "[bc_z_max]", "type = periodic",
        "", "[initial_conditions]", "type             = couette", f"amr_iterations   = {0 if uniform else 1}",
        "wall_y0          = 0.0",
        "height           = 1.0", "wall_velocity    = 1.0", "wall_temperature = 1.0",
        "", "[initial_conditions_region_1]", "r = 1.0", "u = 0.0", "v = 0.0", "w = 0.0", "p = 1.0",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    """Write the input."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("size", choices=tuple(SIZES))
    parser.add_argument("out", type=Path)
    parser.add_argument("--uniform", action="store_true", help="one realm, uniform grid, no AMR marker, no regrid")
    parser.add_argument("--realms", type=int, default=1, choices=(1, 2),
                        help="2: the refined-walls variant as a two-realm forest split at x = 0.5 (out = the manifest)")
    parser.add_argument("--steps", type=int, default=20, help="time steps of the run (default 20)")
    parser.add_argument("--frequency", type=int, default=10, help="runtime regrid cadence in steps (default 10)")
    parser.add_argument("--scheme", default="weno-u-5", choices=SCHEMES,
                        help="the WENO scheme (default weno-u-5); ngc follows its stencil half-width")
    args = parser.parse_args()
    if args.realms == 1:
        args.out.write_text(profile_ini(args.size, args.steps, args.frequency, args.uniform, args.scheme))
        return
    if args.uniform:
        parser.error("--realms 2 is a variant of the refined walls, not of --uniform")
    stem = args.out.stem
    for n, (xmin, xmax) in enumerate(((0.0, 0.5), (0.5, 1.0)), start=1):
        realm = (xmin, xmax, f"couette3d-r{n}")
        realm_ini = profile_ini(args.size, args.steps, args.frequency, False, args.scheme, realm)
        (args.out.parent / f"{stem}-r{n}.ini").write_text(realm_ini)
    args.out.write_text(
        f"; FLUME profiling case: 3-D Couette channel, size {args.size}, two realms split at x = 0.5 (make_profile.py "
        "--realms 2: rerun it, do not edit)\n"
        f"[forest]\nrealms_number = 2\n\n[realm.1]\nini = {stem}-r1.ini\n\n[realm.2]\nini = {stem}-r2.ini\n\n"
        "[forest.topology]\ninter_realm_faces_number = 1\n\n[forest.topology.face_1]\nrealm_a          = 1\n"
        "face_a           = +x\nrealm_b          = 2\nface_b           = -x\ncoupling         = auto\n"
        "coupling_cadence = stage_coincident\n")


if __name__ == "__main__":
    main()
