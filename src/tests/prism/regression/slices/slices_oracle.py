#!/usr/bin/env python3
"""PRISM verification oracle of the [slices] output (slices_restauration.md, F1/F5).

Why this oracle exists: the slices (library slices_object, MPI-IO .mat files) are derived outputs nothing reads back,
so a wrong value goes unnoticed. Two checks on the XH5F checkpoint saved at the slice step (with ghost cells):

* slices: every point of every slice of every field set (q, curl, div) is recomputed from the checkpoint blocks: the
  leaf block containing the written point, the eight cells around it (ghost cells included) and the trilinear or inverse
  distance weights, as adam_object%interpolate_at_point defines them. The point coordinates are checked against the
  [slice_N] box, which also catches a point no rank wrote.
* auxiliary ghosts: curl and divergence are computed on the interior cells only and save_slices refreshes their
  intra-realm ghost cells; on every face shared by two blocks of the same level the ghost slab of one block must equal
  the interior slab of the other, bitwise, for every field set. Faces across the 2:1 seam are interpolated, not
  copied, and are counted but not checked.

.mat layout (adam_object%save_slice): 4 int32 (ni, nj, nk, nvar = 3 + nq), then per point 3 float64 coordinates and
nq float64 variables, points in (i, j, k) order, i fastest.

Usage:
    slices_oracle.py <work-dir> [--ini NAME] [--tol T]
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

# field sets in the order of the state vector (io_initialize, EM model without phi/psi/rho)
SETS = {
    "": ("Dx", "Dy", "Dz", "Bx", "By", "Bz", "Jx", "Jy", "Jz"),
    "curl": ("curl_Dx", "curl_Dy", "curl_Dz", "curl_Bx", "curl_By", "curl_Bz", "curl_Jx", "curl_Jy", "curl_Jz"),
    "div": ("div_D", "div_B", "div_J", "totvar", "div_05", "div_06", "div_09", "div_10", "div_11"),
}


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read a PRISM INI file (`;` comments, no interpolation)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(path)
    return ini


def load_blocks(work: Path, basename: str, it: int, ngc: int) -> list[dict]:
    """Load the blocks of the checkpoint at iteration `it`: interior box, spacing, field sets with ghosts [v,i,j,k]."""
    blocks = []
    for path in sorted(work.glob(f"{basename}-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                dx = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + ngc * dx
                fields = {s: np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in names])
                          for s, names in SETS.items() if f"{blk}-{names[0]}" in h5}
                n = np.array(next(iter(fields.values())).shape[1:]) - 2 * ngc
                blocks.append({"name": f"{path.name}:{blk}", "lo": lo, "dx": dx, "n": n, "hi": lo + n * dx,
                               "q": fields})
    return blocks


def owner(blocks: list[dict], point: np.ndarray, domain: tuple) -> dict:
    """Return the leaf block owning `point` by the library rule (grid%get_closest_block, tree%get_closest_block): the
    block index on the finest level is int((point - domain_emin) / finest_block_size), then the living ancestor. A point
    within a few ulps of a block face can so belong to the block on its + side and use that block's ghost cells."""
    emin, size_max = domain
    idx = np.floor((point - emin) / size_max).astype(int)
    found = []
    for b in blocks:
        ratio = int(round(float(b["n"][0] * b["dx"][0] / size_max[0])))
        bidx = np.rint((b["lo"] - emin) / (b["n"] * b["dx"])).astype(int)
        if np.all(idx // ratio == bidx):
            found.append(b)
    if len(found) != 1:
        raise RuntimeError(f"point {point}: {len(found)} leaf blocks own it")
    return found[0]


def interpolate(blocks: list[dict], point: np.ndarray, itype: str, ngc: int, domain: tuple) -> dict[str, tuple]:
    """Interpolate every field set at `point` as adam_object%interpolate_at_point does; return, per field set, the
    values and the largest magnitude of the eight stencil values (the round-off scale of the convex combination)."""
    b = owner(blocks, point, domain)
    dx = b["dx"]
    ijk1 = np.floor((point - 0.5 * dx - b["lo"]) / dx).astype(int) + 1  # 1-based, may be a ghost (0 or n + 1)
    corners = [np.array([i, j, k]) for k in (0, 1) for j in (0, 1) for i in (0, 1)]
    cells = [ijk1 + c for c in corners]
    xyz = [b["lo"] + (c - 0.5) * dx for c in cells]
    out = {}
    for s, q in b["q"].items():
        val = [q[:, c[0] - 1 + ngc, c[1] - 1 + ngc, c[2] - 1 + ngc] for c in cells]
        if itype == "trilinear":
            t = (point - xyz[0]) / (xyz[7] - xyz[0])
            p = 1.0 - t
            qx = [p[0] * val[2 * m] + t[0] * val[2 * m + 1] for m in range(4)]
            qxy = [p[1] * qx[0] + t[1] * qx[1], p[1] * qx[2] + t[1] * qx[3]]
            out[s] = (p[2] * qxy[0] + t[2] * qxy[1], np.max(np.abs(val), axis=0))
        elif itype == "inverse_distance":
            w = np.array([1.0 / float(np.sum((x - point) ** 2)) ** 3 for x in xyz])
            out[s] = (sum(wi * vi for wi, vi in zip(w, val)) / w.sum(), np.max(np.abs(val), axis=0))
        else:
            raise RuntimeError(f"unknown itype {itype}")
    return out


def check_slices(work: Path, ini: configparser.ConfigParser, blocks: list[dict], it: int, ngc: int, domain: tuple,
                 tol: float) -> int:
    """Check every slice file against the recomputed values; return the failures count."""
    basename = ini["IO"]["output_basename"]
    failures = 0
    for s in range(1, int(ini["slices"]["slices_number"]) + 1):
        sec = ini[f"slice_{s}"]
        nijk = np.array([int(sec["ni"]), int(sec["nj"]), int(sec["nk"])])
        emin = np.array([float(sec[f"emin_{a}"]) for a in "xyz"])
        emax = np.array([float(sec[f"emax_{a}"]) for a in "xyz"])
        expected = [emin + (np.array([i, j, k]) + 0.5) * (emax - emin) / nijk
                    for k in range(nijk[2]) for j in range(nijk[1]) for i in range(nijk[0])]
        for fset in blocks[0]["q"]:
            tag = f"{basename}-{fset}-" if fset else f"{basename}-"
            mat = work / f"{tag}slice_{s:02d}-{it:09d}.mat"
            if not mat.exists():
                print(f"slice {mat.name}: missing  FAIL")
                failures += 1
                continue
            raw = mat.read_bytes()
            head = np.frombuffer(raw[:16], dtype=np.int32)
            nq = len(SETS[fset])
            data = np.frombuffer(raw[16:], dtype=np.float64)
            ok_head = list(head) == [*nijk, 3 + nq] and data.size == nijk.prod() * (3 + nq)
            if not ok_head:
                print(f"slice {mat.name}: header {list(head)}, {data.size} values  FAIL")
                failures += 1
                continue
            data = data.reshape(nijk.prod(), 3 + nq)
            xyz_err = float(np.max(np.abs(data[:, :3] - np.array(expected))))
            # recompute at the written coordinates: a point on a block face (x = 0 here, on the seam) is owned by one
            # side or the other by its last bits (gfortran writes 0, nvfortran -2.7e-18)
            ref = [interpolate(blocks, p, sec["itype"].strip(), ngc, domain) for p in data[:, :3]]
            val = np.array([r[fset][0] for r in ref])
            # each value is a convex combination of its eight stencil values: its round-off scales with their largest
            # magnitude (a variable can be round-off noise along a slice, where any wider scale is meaningless)
            scale = np.array([r[fset][1] for r in ref])
            scale = np.where(scale > 0.0, scale, 1.0)
            err = float(np.max(np.abs(data[:, 3:] - val) / scale))
            ok = xyz_err <= 1.0e-15 and err <= tol and np.all(np.isfinite(data))
            failures += 0 if ok else 1
            print(f"slice {mat.name}: {'x'.join(map(str, nijk))} {sec['itype'].strip()}, coordinates error "
                  f"{xyz_err:.1e}, max relative difference {err:.3e}  {'PASS' if ok else 'FAIL'} (tol {tol:.1e})")
    return failures


def check_ghosts(blocks: list[dict], ngc: int) -> int:
    """Check the face ghost slabs against the same-level neighbour interior slabs; return the failures count."""
    checked, skipped, worst = 0, 0, {s: 0.0 for s in blocks[0]["q"]}
    for a in blocks:
        for b in blocks:
            if a is b:
                continue
            for ax in range(3):
                if not np.isclose(a["hi"][ax], b["lo"][ax], rtol=0.0, atol=1.0e-12 * a["dx"][ax]):
                    continue
                other = [o for o in range(3) if o != ax]
                overlap = all(min(a["hi"][o], b["hi"][o]) - max(a["lo"][o], b["lo"][o]) > 0.5 * a["dx"][o]
                              for o in other)
                if not overlap:
                    continue
                if not np.allclose(a["dx"], b["dx"], rtol=1.0e-12, atol=0.0):
                    skipped += 1
                    continue
                checked += 1
                n = a["n"]
                for s in a["q"]:
                    # a's ghost slab beyond its + face along ax == b's first interior slab
                    ga = [slice(ngc, ngc + n[o]) for o in range(3)]
                    gb = [slice(ngc, ngc + n[o]) for o in range(3)]
                    ga[ax] = slice(ngc + n[ax], 2 * ngc + n[ax])
                    gb[ax] = slice(ngc, 2 * ngc)
                    da = a["q"][s][(slice(None), *ga)]
                    db = b["q"][s][(slice(None), *gb)]
                    # and b's ghost slab before its - face == a's last interior slab
                    ha = [slice(ngc, ngc + n[o]) for o in range(3)]
                    hb = [slice(ngc, ngc + n[o]) for o in range(3)]
                    ha[ax] = slice(n[ax], ngc + n[ax])
                    hb[ax] = slice(0, ngc)
                    worst[s] = max(worst[s], float(np.max(np.abs(da - db))),
                                   float(np.max(np.abs(a["q"][s][(slice(None), *ha)] - b["q"][s][(slice(None), *hb)]))))
    ok = all(v == 0.0 for v in worst.values()) and checked > 0
    detail = ", ".join(f"{s or 'q'} {v:.1e}" for s, v in worst.items())
    print(f"face ghosts vs same-level neighbour interior: {checked} faces checked, {skipped} across levels skipped; "
          f"max difference {detail}  {'PASS' if ok else 'FAIL'} (bitwise)")
    return 0 if ok else 1


def main() -> int:
    """Run the checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="work directory of the run")
    parser.add_argument("--ini", default="slices.ini", help="input file name inside the work directory")
    parser.add_argument("--tol", type=float, default=1.0e-12, help="relative tolerance of the slice values")
    args = parser.parse_args()

    ini = read_ini(args.work / args.ini)
    ngc = int(ini["grid"]["ngc"])
    it = int(ini["time"]["it_max"])
    blocks = load_blocks(args.work, ini["IO"]["output_basename"], it, ngc)
    if not blocks:
        sys.exit(f"slices_oracle: no checkpoint at iteration {it} in {args.work}")
    # the finest block size, as grid%block_dxyz(:, max_level)
    emin = np.array([float(ini["grid"][f"emin_{a}"]) for a in "xyz"])
    emax = np.array([float(ini["grid"][f"emax_{a}"]) for a in "xyz"])
    refined = {2: 1, 4: 2, 8: 3}[int(ini["amr"].get("ratio", "8"))]  # binary, quad, octree: refined axes x[, y[, z]]
    size_max = (emax - emin) / np.array([2 ** int(ini["amr"]["max_level"]) if a < refined else 1 for a in range(3)])
    domain = (emin, size_max)
    failures = check_slices(args.work, ini, blocks, it, ngc, domain, args.tol)
    failures += check_ghosts(blocks, ngc)
    print(f"slices verification {'PASSED' if failures == 0 else 'FAILED'} ({failures} failures)")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
