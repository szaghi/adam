#!/usr/bin/env python3
"""Planes of a FLUME checkpoint for the turbulence figures (#119, #121): a few 2-D cuts of a 3-D run, small enough to
keep next to the run and plot anywhere.

The 3-D checkpoints of the boundary-layer, channel and interaction runs are hundreds of MB each; the figures need three
cuts and a spanwise average. This reads one checkpoint block by block on the machine that holds it and writes an npz:

- ``xy_*``  the x-y plane nearest the spanwise position ``--z`` (default mid-span),
- ``xz_*``  the wall-parallel x-z plane nearest ``--y`` (default the first cell above the wall),
- ``yz_*``  the cross-flow y-z plane nearest ``--x`` (default mid-length),
- ``avg_*`` the x-y plane averaged over z (all cells),

each with the density ``r``, the velocity ``u, v, w`` and the pressure ``p`` (gamma from the run's input, default 1.4),
plus the cell-centre coordinates ``x, y, z`` of the uniform grid. One refinement level only (the runs it serves are
uniform); a mixed grid is refused.

Usage: field_slices.py <work> [--it N | --time T] [--x X] [--y Y] [--z Z] [--out FILE] (default <work>/slices-<it>.npz)
"""

from __future__ import annotations

import argparse
import configparser
import re
import sys
from pathlib import Path

import h5py
import numpy as np

NAMES = ("r", "ru", "rv", "rw", "rE")


def run_setup(work: Path) -> tuple[int, float]:
    """Return the ghost width and gamma of a run (its input.ini)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=("!", "#", ";"), strict=False)
    ini.read(work / "input.ini")
    ngc = ini.getint("grid", "ngc", fallback=3)
    gamma = 1.4
    for key in ("gamma", "gamma_gas"):
        if ini.has_option("physics", key):
            gamma = float(ini.get("physics", key).replace("d", "e").replace("D", "e").split()[0])
    return ngc, gamma


def checkpoint(work: Path, it: int | None) -> tuple[int, list[Path]]:
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        raise FileNotFoundError(f"no checkpoint in {work}")
    its = sorted({int(p.name.split("-")[-2]) for p in files})
    it = its[-1] if it is None else it
    return it, [p for p in files if int(p.name.split("-")[-2]) == it]


def block_names(h5: h5py.File) -> list[str]:
    return sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")})


def primitives(q: dict[str, np.ndarray], gamma: float) -> dict[str, np.ndarray]:
    r = q["r"]
    u, v, w = q["ru"] / r, q["rv"] / r, q["rw"] / r
    return {"r": r, "u": u, "v": v, "w": w, "p": (gamma - 1.0) * (q["rE"] - 0.5 * r * (u * u + v * v + w * w))}


def extract(work: Path, it: int | None, x: float | None, y: float | None, z: float | None) -> dict[str, np.ndarray]:
    ngc, gamma = run_setup(work)
    it, files = checkpoint(work, it)
    geo = []  # pass 1: the grid
    for path in files:
        with h5py.File(path, "r") as h5:
            for blk in block_names(h5):
                d = h5[f"{blk}-dxdydz"][()][::-1]
                lo = h5[f"{blk}-origin"][()][::-1] + ngc * d
                n = np.array(h5[f"{blk}-r"].shape[::-1]) - 2 * ngc
                geo.append((lo, d, n))
    ds = np.array([g[1] for g in geo])
    if np.ptp(ds, axis=0).max() > 1e-9 * ds.max():
        raise SystemExit(f"{work}: more than one cell size, field_slices handles one level only")
    d = ds[0]
    lo = np.min([g[0] for g in geo], axis=0)
    hi = np.max([g[0] + g[2] * g[1] for g in geo], axis=0)
    n = np.rint((hi - lo) / d).astype(int)
    c = [lo[a] + (np.arange(n[a]) + 0.5) * d[a] for a in range(3)]
    pos = [x if x is not None else 0.5 * (lo[0] + hi[0]), y if y is not None else c[1][0], z if z is not None else 0.5 * (lo[2] + hi[2])]
    cut = [int(np.clip(np.floor((pos[a] - lo[a]) / d[a]), 0, n[a] - 1)) for a in range(3)]
    keys = ("r", "u", "v", "w", "p")
    out = {f"{p}_{k}": np.full(shape, np.nan) for p, shape in (("xy", (n[0], n[1])), ("xz", (n[0], n[2])), ("yz", (n[1], n[2])))
           for k in keys}
    out |= {f"avg_{k}": np.zeros((n[0], n[1])) for k in keys}
    count = np.zeros((n[0], n[1]))
    for path in files:  # pass 2: the planes
        with h5py.File(path, "r") as h5:
            for blk in block_names(h5):
                b0 = np.rint((h5[f"{blk}-origin"][()][::-1] + ngc * d - lo) / d).astype(int)
                q = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc] for v in NAMES}
                f = primitives(q, gamma)
                nb = f["r"].shape
                s = tuple(slice(b0[a], b0[a] + nb[a]) for a in range(3))
                for k in keys:
                    out[f"avg_{k}"][s[0], s[1]] += f[k].sum(axis=2)
                count[s[0], s[1]] += nb[2]
                if b0[2] <= cut[2] < b0[2] + nb[2]:
                    for k in keys:
                        out[f"xy_{k}"][s[0], s[1]] = f[k][:, :, cut[2] - b0[2]]
                if b0[1] <= cut[1] < b0[1] + nb[1]:
                    for k in keys:
                        out[f"xz_{k}"][s[0], s[2]] = f[k][:, cut[1] - b0[1], :]
                if b0[0] <= cut[0] < b0[0] + nb[0]:
                    for k in keys:
                        out[f"yz_{k}"][s[1], s[2]] = f[k][cut[0] - b0[0], :, :]
    for k in keys:
        out[f"avg_{k}"] /= np.where(count > 0, count, np.nan)
    t = float("nan")  # the checkpoint time is in its XDMF (the HDF5 holds the iteration only)
    for xdmf in work.glob(f"*-{it:09d}.xdmf"):
        m = re.search(r'<Time Value="([^"]+)"', xdmf.read_text()[:4000])
        if m:
            t = float(m.group(1))
    out |= {"x": c[0], "y": c[1], "z": c[2], "cut": np.array([c[a][cut[a]] for a in range(3)]), "time": np.array(t),
            "it": np.array(it), "gamma": np.array(gamma)}
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("work", type=Path)
    ap.add_argument("--it", type=int, default=None, help="checkpoint iteration (default: the last)")
    ap.add_argument("--time", type=float, default=None, help="instead, the checkpoint nearest this time (its XDMF)")
    ap.add_argument("--x", type=float, default=None, help="streamwise position of the y-z plane (default mid-length)")
    ap.add_argument("--y", type=float, default=None, help="height of the x-z plane (default the first cell)")
    ap.add_argument("--z", type=float, default=None, help="spanwise position of the x-y plane (default mid-span)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    it = a.it
    if a.time is not None:
        times = {}
        for xdmf in a.work.glob("*-[0-9]*.xdmf"):
            m = re.search(r'<Time Value="([^"]+)"', xdmf.read_text()[:4000])
            if m and xdmf.stem.split("-")[-1].isdigit() and any(a.work.glob(f"{xdmf.stem}-proc*.h5")):
                times[int(xdmf.stem.split("-")[-1])] = float(m.group(1))
        it = min(times, key=lambda k: abs(times[k] - a.time))
    s = extract(a.work, it, a.x, a.y, a.z)
    out = a.out or a.work / f"slices-{int(s['it']):09d}.npz"
    np.savez_compressed(out, **{k: v.astype(np.float32) if v.ndim == 2 else v for k, v in s.items()})
    print(f"{out}: t {float(s['time']):.4g}, grid {len(s['x'])} x {len(s['y'])} x {len(s['z'])}, cuts x {s['cut'][0]:.4g} y {s['cut'][1]:.4g} z {s['cut'][2]:.4g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
