#!/usr/bin/env python3
"""Draw the figures of the PRISM verification page (docs/applications/prism/verification.md) from regression runs.

Why this script exists: the documentation shows real PRISM results. Every figure is drawn from the checkpoints of a
regression check run with --keep, so it can be regenerated after any change of the solver.

Figures:
    quadtree-pulse  rmf-amr-fd-pulse/check.sh --keep (quadtree leg, issue #46): the source-free pulse across the 2:1
                    seam on an octree and on a quadtree, Bx and By planes with blocks outlined, the column-by-column
                    difference of the two trees, and their divergence histories.

Usage:
    make_doc_figures.py --out DIR [--only NAME ...]

Needs numpy, h5py and matplotlib.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
DPI = 110
NGC = 3
FIELDS = ("Dx", "Dy", "Dz", "Bx", "By", "Bz")


def last_blocks(work: Path) -> tuple[int, list[dict]]:
    """Return the last saved step and its blocks: interior fields (x, y, z order), spacing and interior corner."""
    files = [p for p in work.glob("*-proc*.h5") if "restart" not in p.name and p.name.split("-")[-2].isdigit()]
    if not files:
        raise FileNotFoundError(f"no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    out = []
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                dxyz = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + NGC * dxyz
                f = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[NGC:-NGC, NGC:-NGC, NGC:-NGC] for v in FIELDS}
                out.append({"lo": lo, "d": dxyz, "f": f})
    return last, out


def slab(bl: list[dict], axis: int, at: float) -> list[tuple[dict, int]]:
    """Return the blocks crossed by the plane `coordinate[axis] = at`, with the cell index of the plane."""
    sel = []
    for b in bl:
        n = b["f"]["Bx"].shape[axis]
        k = int(np.floor((at - b["lo"][axis]) / b["d"][axis]))
        if 0 <= k < n:
            sel.append((b, k))
    return sel


def plane_map(ax, bl: list[dict], var: str, normal: int, at: float, cmap: str,  # noqa: ANN001, ANN201
              lim: tuple | None = None):
    """Draw `var` on the plane normal to axis `normal` at `at`, blocks outlined; return the image."""
    axes = [a for a in range(3) if a != normal]
    sel = slab(bl, normal, at)
    vals = [np.take(b["f"][var], k, axis=normal) for b, k in sel]
    lo, hi = lim if lim else (min(float(v.min()) for v in vals), max(float(v.max()) for v in vals))
    im = None
    for (b, _), v in zip(sel, vals, strict=True):
        e0 = b["lo"][axes[0]] + np.arange(v.shape[0] + 1) * b["d"][axes[0]]
        e1 = b["lo"][axes[1]] + np.arange(v.shape[1] + 1) * b["d"][axes[1]]
        im = ax.pcolormesh(e0, e1, v.T, cmap=cmap, vmin=lo, vmax=hi, shading="flat")
        ax.add_patch(Rectangle((e0[0], e1[0]), e0[-1] - e0[0], e1[-1] - e1[0], fill=False, lw=0.4, ec="k"))
    ax.set_aspect("equal")
    ax.set_xlabel("xyz"[axes[0]])
    ax.set_ylabel("xyz"[axes[1]])
    return im


def column_key(b: dict, i: int, j: int) -> tuple[float, float]:
    """Return the rounded (x, y) centre of cell (i, j) of a block."""
    return (round(float(b["lo"][0] + (i + 0.5) * b["d"][0]), 12), round(float(b["lo"][1] + (j + 0.5) * b["d"][1]), 12))


def columns(bl: list[dict]) -> dict[tuple[float, float], np.ndarray]:
    """Return the first z cell of every (x, y) column: its fields, keyed by the rounded cell centre."""
    cols = {}
    for b in bl:
        nx, ny, _ = b["f"]["Bx"].shape
        for i in range(nx):
            for j in range(ny):
                cols.setdefault(column_key(b, i, j), np.array([b["f"][v][i, j, 0] for v in FIELDS]))
    return cols


def history(work: Path) -> np.ndarray:
    """Return the divergence history rows (it, nblk, time, div(D), div(B), div(J))."""
    path = next(work.glob("*-divergence_history.dat"))
    rows = [line.split() for line in path.read_text().splitlines() if line.startswith("+")]
    return np.array([[float(x) for x in r[:6]] for r in rows])


def fig_quadtree_pulse(out: Path) -> None:
    """rmf-amr-fd-pulse quadtree leg: octree and quadtree, fields, difference, divergence histories."""
    case = HERE / "rmf-amr-fd-pulse"
    oct_w, quad_w = case / "work-cpu-quad-octree", case / "work-cpu-quadtree"
    step, oct_b = last_blocks(oct_w)
    _, quad_b = last_blocks(quad_w)
    z0 = y0 = 1.0e-4  # just off the block faces at 0, through the By peak
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.6), constrained_layout=True)
    bmax = max(float(np.abs(b["f"]["Bx"]).max()) for b in oct_b)
    lim = (-bmax, bmax)
    im = plane_map(axs[0, 0], quad_b, "Bx", 2, z0, "RdBu_r", lim)
    fig.colorbar(im, ax=axs[0, 0], label=r"$B_x$", shrink=0.8)
    axs[0, 0].set_title(f"quadtree, $B_x$ on z = 0 (step {step}): 2:1 seam at x = 0")
    bymax = max(float(np.abs(b["f"]["By"]).max()) for b in oct_b)
    for ax, bl in ((axs[0, 1], oct_b), (axs[0, 2], quad_b)):
        im = plane_map(ax, bl, "By", 1, y0, "PuOr", (-bymax, bymax))
        fig.colorbar(im, ax=ax, label=r"$B_y$", shrink=0.8)
    axs[0, 1].set_title("octree, $B_y$ on y = 0: x < 0 blocks split in z too")
    axs[0, 2].set_title("quadtree, $B_y$ on y = 0: one z layer of blocks")
    oc, qc = columns(oct_b), columns(quad_b)
    scale = np.array([max(abs(c[i]) for c in oc.values()) for i in range(len(FIELDS))])
    for grp in ((0, 1, 2), (3, 4, 5)):
        scale[list(grp)] = scale[list(grp)].max()
    diff = {k: float(np.max(np.abs(qc[k] - oc[k]) / scale)) for k in qc}
    floor = 1.0e-18
    for b, _ in slab(quad_b, 2, z0):
        nx, ny = b["f"]["Bx"].shape[:2]
        v = np.array([[max(diff[column_key(b, i, j)], floor) for j in range(ny)] for i in range(nx)])
        xe = b["lo"][0] + np.arange(nx + 1) * b["d"][0]
        ye = b["lo"][1] + np.arange(ny + 1) * b["d"][1]
        im = axs[1, 0].pcolormesh(xe, ye, v.T, cmap="viridis", norm=LogNorm(vmin=floor, vmax=1e-15), shading="flat")
    axs[1, 0].set_aspect("equal")
    fig.colorbar(im, ax=axs[1, 0], label="max relative difference over D, B (0 shown as 1e-18)", shrink=0.8)
    axs[1, 0].set_title(f"quadtree vs octree, per (x, y) column: max {max(diff.values()):.1e}")
    axs[1, 0].set_xlabel("x")
    axs[1, 0].set_ylabel("y")
    im = plane_map(axs[1, 1], quad_b, "By", 2, z0, "PuOr")
    fig.colorbar(im, ax=axs[1, 1], label=r"$B_y$", shrink=0.8)
    axs[1, 1].set_title(r"quadtree, $B_y$ on z = 0: generated at the seam")
    ho, hq = history(oct_w), history(quad_w)
    axs[1, 2].semilogy(ho[:, 0], np.abs(ho[:, 4]), "s-", mfc="none", label="octree max|div(B)|")
    axs[1, 2].semilogy(hq[:, 0], np.abs(hq[:, 4]), "x--", label="quadtree max|div(B)|")
    axs[1, 2].semilogy(ho[:, 0], np.maximum(np.abs(ho[:, 3]), 1e-17), "o:", mfc="none",
                       label="octree max|div(D)| (0 shown as 1e-17)")
    axs[1, 2].semilogy(hq[:, 0], np.maximum(np.abs(hq[:, 3]), 1e-17), "+:", label="quadtree max|div(D)|")
    axs[1, 2].set_xlabel("step")
    axs[1, 2].grid(alpha=0.3, which="both")
    axs[1, 2].legend(fontsize=8)
    axs[1, 2].set_title("divergence histories: the same seam source on both trees")
    fig.suptitle("PRISM source-free pulse across a 2:1 seam (rmf-amr-fd-pulse): octree against quadtree")
    fig.savefig(out / "quadtree-pulse.png", dpi=DPI)
    plt.close(fig)


FIGURES = {"quadtree-pulse": fig_quadtree_pulse}


def main() -> int:
    """Draw the requested figures; report the missing data."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True, help="output directory of the PNG files")
    ap.add_argument("--only", nargs="*", default=None, help="figures to draw (default: all)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    status = 0
    for name, fn in FIGURES.items():
        if a.only and name not in a.only:
            continue
        try:
            fn(a.out)
            print(f"{name}: done")
        except (FileNotFoundError, StopIteration, KeyError) as err:
            print(f"{name}: skipped ({err!r})")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
