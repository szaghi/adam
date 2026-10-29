#!/usr/bin/env python3
"""Write the vortex forest of multirealm/check.sh legs 12 and 13 (issue #110 P3, P5): a feature along an auto seam.

The isentropic vortex of the vortex-periodic regression case (radius 0.07), started at (0.47, 0.3) across the seam in
a free stream (0, 1) along it, in the box [0, 1]^2 x [0, 0.5] split at x = 0.5 into two `coupling = auto` realms.
Each realm is an octree of 8x16x8 cells per block (inter-realm 2:1 faces need even block cell counts on every axis, so
a quadtree is refused), cubic cells of 1/32 at base level 1, `max_level` 3. y is periodic, the x faces are inviscid
walls (mirrored ghosts, no mass flux) and z is extrapolated with w = 0, so the box is closed and the mass of the two
realms is conserved. Both realms regrid every 5 steps by a Loehner marker on the density (refine 0.3, derefine 0.1).
The refined patch covers part of the seam and moves along it, so the seam turns mixed.

Usage:
    make_vortex_forest.py <out-dir> [--it-save N]

Writes <out-dir>/vortex-2realm.ini (the manifest), vortex-r1.ini and vortex-r2.ini. --it-save N keeps a checkpoint
every N steps (the documentation figures); by default the inputs of check.sh, which save the first and last steps.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "regrid"))
import make_regrid as mr  # noqa: E402


def realm_text(name: str, x0: str, x1: str, bmin: str, bmax: str, it_save: int | None) -> str:
    """Return the input of one realm of the forest: x in [x0, x1], walls on the x faces bmin/bmax of the box."""
    text = mr.VORTEX.read_text()
    for key, value in (
        ("ni", "8"),
        ("nj", "16"),
        ("nk", "8"),
        ("ratio", "8"),
        ("emax_z", "0.5"),
        ("null_z", ".false."),
        ("emin_x", x0),
        ("emax_x", x1),
        ("iu_ref_levels", "1"),
        ("max_level", "3"),
        ("iters", "2"),
        ("frequency", "5"),
        ("markers_number", "1"),
        ("x0", "0.47"),
        ("y0", "0.3"),
        ("u", "0.0"),
        ("v", "1.0"),
        ("output_basename", f"vortex-{name}"),
        ("restart_basename", f"vortex-{name}-restart"),
        ("it_max", "-1"),
        ("time_max", "0.3"),
        ("CFL", "0.4"),
    ):
        text = mr.setkey(text, key, value)
    if it_save is not None:
        text = mr.setkey(text, "it_save", str(it_save))
    text = text.replace("[bc_x_min]\ntype = periodic", f"[bc_x_min]\ntype = {bmin}", 1)
    text = text.replace("[bc_x_max]\ntype = periodic", f"[bc_x_max]\ntype = {bmax}", 1)
    if text.count("type = periodic") != 2:  # y only
        sys.exit(f"make_vortex_forest.py: realm {name}: the x faces of the vortex input are not periodic")
    return text.replace(
        "[field]",
        mr.LOHNER.format(n=1, field=1, var=1, refine="0.3", derefine="0.1", floor="0.0") + "[field]",
        1,
    )


def main() -> None:
    """Write the manifest and the two realm inputs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", type=Path)
    parser.add_argument("--it-save", type=int, default=None, help="a checkpoint every N steps")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, x0, x1, bmin, bmax in (
        ("r1", "0.0", "0.5", "wall-inviscid", "extrapolation"),
        ("r2", "0.5", "1.0", "extrapolation", "wall-inviscid"),
    ):
        (args.out / f"vortex-{name}.ini").write_text(realm_text(name, x0, x1, bmin, bmax, args.it_save))
    (args.out / "vortex-2realm.ini").write_text(
        "; legs 12 and 13 of multirealm/check.sh (issue #110 P3, P5): the vortex moving along part of an auto seam\n"
        "[forest]\nrealms_number = 2\n\n[realm.1]\nini = vortex-r1.ini\n\n[realm.2]\nini = vortex-r2.ini\n\n"
        "[forest.topology]\ninter_realm_faces_number = 1\n\n[forest.topology.face_1]\nrealm_a          = 1\n"
        "face_a           = +x\nrealm_b          = 2\nface_b           = -x\ncoupling         = auto\n"
        "coupling_cadence = stage_coincident\n"
    )


if __name__ == "__main__":
    main()
