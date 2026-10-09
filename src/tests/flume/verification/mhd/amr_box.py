"""Refine a box of an MHD verification input by one 2:1 level (issue #41, M2-P7b).

Why: the AMR variants of the 2-D MHD verifications (MV-9 field loop, MV-12 Orszag-Tang) reuse their uniform inputs and
refine the blocks whose centroid lies in a box, at initialisation only (the V3 pattern: geometric marker
primitive-box, no regridding afterwards): max_level grows by one, and the marker targets that level. The 2:1 seam
machinery (restriction, ghost interpolation, reflux) assumed an octree until issue #46: on a quadtree (ratio 4) the
coarse cells beside a seam picked up a spurious z dependence (M2-P7b: 1e-2 at the first step, negative pressure in
Orszag-Tang at step 26, the same on the Euler path), so the variants run on an octree (ratio 8) with null z, where the
coarse-fine ghost interpolation needs at least 4 cells per block on every axis, null ones included (M2-P4): nk grows to
4. Both are costs of the AMR variant, not physical changes. Since #46 a quadtree seam is 2:1 in x and y and 1:1 in z,
so ratio = 4 keeps the base nk (1 for a true 2-D run, 4 times fewer cells than the octree variant).
"""

from __future__ import annotations

import configparser


def refine_box(ini: configparser.ConfigParser, box: list[float], ratio: int = 8) -> None:
    """Add one 2:1 level over the blocks whose centroid lies in box = [xmin, ymin, xmax, ymax] (all z).

    ratio is the tree ratio: 8 (octree, nk raised to 4) or 4 (quadtree, nk kept).
    """
    if ratio not in (4, 8):
        raise ValueError(f"refine_box: ratio must be 4 or 8, not {ratio}")
    level = int(ini["amr"]["max_level"]) + 1
    ini["amr"].update({"max_level": str(level), "ratio": str(ratio), "markers_number": "1", "frequency": "0"})
    ini["amr_marker_1"] = {"mode": "1", "geo_type": "primitive-box", "delta_type": "max", "delta_fine": "0.0",
                           "delta_coarse": "0.0", "box_xmin": repr(box[0]), "box_ymin": repr(box[1]),
                           "box_zmin": repr(-1.0e30), "box_xmax": repr(box[2]), "box_ymax": repr(box[3]),
                           "box_zmax": repr(1.0e30), "target_level": str(level)}
    ini["initial_conditions"]["amr_iterations"] = "1"
    if ratio == 8:
        ini["grid"]["nk"] = str(max(4, int(ini["grid"]["nk"])))
