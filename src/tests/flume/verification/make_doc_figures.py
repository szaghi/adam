#!/usr/bin/env python3
"""Draw the figures of the FLUME documentation (docs/applications/flume/) from verification runs.

Why this script exists: the documentation shows real FLUME results, not sketches. Every figure is drawn from the
checkpoints of a verification input (the case directories of this tree, or the same inputs re-run into `--runs`), with
the exact or reference solution of the corresponding oracle where one exists, so a figure can be regenerated after any
change of the solver and compared with the published one.

Usage:
    make_doc_figures.py --out DIR [--runs DIR] [--only NAME ...]

`--runs` holds the re-runs of inputs whose work directories keep no checkpoint, one sub-directory per case, all written
by `make_doc_runs.sh` (its groups: the MHD gallery, the limited blasts, the Euler near-vacuum problems, the linear-wave
ladders, the rotated shock tube, the quadtree, the runtime regrids, the multi-realm unions, the gravity and
self-gravity cases).
A figure whose data is missing is skipped with a message. Needs numpy, h5py and matplotlib.
"""

from __future__ import annotations

import argparse
import configparser
import math
import re
import sys
import tempfile
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
RF = HERE / "riemann-flux"
TAG = "work-adam_flume_cpu-np2"
DPI = 110
EULER = ("r", "ru", "rv", "rw", "rE")
MHD = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz")


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read a FLUME INI file (`;` comments, no interpolation, case-sensitive keys)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(path)
    return ini


def gamma_of(ini: configparser.ConfigParser) -> float:
    """Return the specific heats ratio of an input."""
    return float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])


def last_files(work: Path, prefix: str = "", it: int | None = None) -> list[Path]:
    """Return the checkpoint files of the last saved iteration of a run (of the realm `prefix`, or of iteration `it`)."""
    files = sorted(p for p in work.glob(f"{prefix}*-proc*.h5") if "restart" not in p.name)
    if not files:
        raise FileNotFoundError(f"no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files) if it is None else it
    return [p for p in files if int(p.name.split("-")[-2]) == last]


def blocks(work: Path, ngc: int, names: tuple[str, ...], prefix: str = "", it: int | None = None) -> list[dict]:
    """Return the blocks of the last checkpoint (of the realm `prefix`, or of iteration `it`): interior fields (x, y,
    z order), cell spacing and interior corner."""
    out = []
    for path in last_files(work, prefix, it):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                # geometry and fields are stored (z, y, x) with ghosts: reverse and strip them
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                fields = {}
                for v in names:
                    if f"{blk}-{v}" in h5:
                        a = h5[f"{blk}-{v}"][()].transpose(2, 1, 0)
                        fields[v] = a[ngc:-ngc, ngc:-ngc, ngc:-ngc] if ngc else a
                out.append({"lo": origin + ngc * dxyz, "d": dxyz, "f": fields})
    return out


def profile(work: Path, ngc: int, names: tuple[str, ...], prefix: str = "") -> tuple[int, np.ndarray, np.ndarray]:
    """Return the active axis, the cell centres and the fields of a 1-D run (first transverse copy)."""
    cells: dict[float, np.ndarray] = {}
    axis = -1
    for b in blocks(work, ngc, names, prefix):
        q = np.stack([b["f"][v] for v in names])
        if axis < 0:
            axis = int(np.argmax(q.shape[1:]))
        n = q.shape[1 + axis]
        col = np.moveaxis(q, 1 + axis, 1).reshape(len(names), n, -1)[:, :, 0]
        for c, v in zip(b["lo"][axis] + (np.arange(n) + 0.5) * b["d"][axis], col.T, strict=True):
            cells[float(c)] = v
    xs = np.array(sorted(cells))
    return axis, xs, np.array([cells[c] for c in xs]).T


def primitives(q: np.ndarray, gamma: float, axis: int = 0) -> dict[str, np.ndarray]:
    """Return density, velocity components along the axis frame, pressure (and field) of conservative profiles."""
    r = q[0]
    out = {"r": r, "u": q[1 + axis] / r, "v": q[1 + (axis + 1) % 3] / r, "w": q[1 + (axis + 2) % 3] / r}
    kin = 0.5 * (q[1] ** 2 + q[2] ** 2 + q[3] ** 2) / r
    mag = 0.5 * (q[5] ** 2 + q[6] ** 2 + q[7] ** 2) if q.shape[0] >= 8 else 0.0
    out["p"] = (gamma - 1.0) * (q[4] - kin - mag)
    if q.shape[0] >= 8:
        out["bn"], out["bt1"], out["bt2"] = q[5 + axis], q[5 + (axis + 1) % 3], q[5 + (axis + 2) % 3]
    return out


def exact_riemann(x: np.ndarray, t: float, x0: float, gamma: float, left: tuple, right: tuple) -> tuple:
    """Return the exact density, velocity and pressure of a 1-D Euler Riemann problem (Toro, chapter 4)."""
    rl, ul, pl = left
    rr, ur, pr = right
    al, ar = math.sqrt(gamma * pl / rl), math.sqrt(gamma * pr / rr)
    g1, g2 = (gamma - 1.0) / (gamma + 1.0), (gamma - 1.0) / (2.0 * gamma)
    ex = 2.0 / (gamma - 1.0)

    def fk(p: float, rk: float, pk: float, ak: float) -> tuple[float, float]:
        if p > pk:
            a_, b_ = 2.0 / ((gamma + 1.0) * rk), g1 * pk
            sq = math.sqrt(a_ / (p + b_))
            return (p - pk) * sq, sq * (1.0 - 0.5 * (p - pk) / (p + b_))
        return ex * ak * ((p / pk) ** g2 - 1.0), (p / pk) ** (-(gamma + 1.0) / (2.0 * gamma)) / (rk * ak)

    p = 0.5 * (pl + pr)
    for _ in range(100):
        fl, dl = fk(p, rl, pl, al)
        fr, dr = fk(p, rr, pr, ar)
        p = max(1e-12, p - (fl + fr + ur - ul) / (dl + dr))
    fl, _ = fk(p, rl, pl, al)
    fr, _ = fk(p, rr, pr, ar)
    us = 0.5 * (ul + ur) + 0.5 * (fr - fl)
    rho, u, pp = np.empty_like(x), np.empty_like(x), np.empty_like(x)
    for n, xn in enumerate(x):
        s = (xn - x0) / t
        if s <= us:
            if p > pl:
                sl = ul - al * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p / pl + g2)
                rho[n], u[n], pp[n] = (rl, ul, pl) if s <= sl else (rl * (p / pl + g1) / (g1 * p / pl + 1.0), us, p)
            elif s <= ul - al:
                rho[n], u[n], pp[n] = rl, ul, pl
            elif s >= us - al * (p / pl) ** g2:
                rho[n], u[n], pp[n] = rl * (p / pl) ** (1.0 / gamma), us, p
            else:
                c = 2.0 / (gamma + 1.0) + g1 / al * (ul - s)
                rho[n], pp[n] = rl * c**ex, pl * c ** (gamma * ex)
                u[n] = 2.0 / (gamma + 1.0) * (al + 0.5 * (gamma - 1.0) * ul + s)
        elif p > pr:
            sr = ur + ar * math.sqrt((gamma + 1.0) / (2.0 * gamma) * p / pr + g2)
            rho[n], u[n], pp[n] = (rr, ur, pr) if s >= sr else (rr * (p / pr + g1) / (g1 * p / pr + 1.0), us, p)
        elif s >= ur + ar:
            rho[n], u[n], pp[n] = rr, ur, pr
        elif s <= us + ar * (p / pr) ** g2:
            rho[n], u[n], pp[n] = rr * (p / pr) ** (1.0 / gamma), us, p
        else:
            c = 2.0 / (gamma + 1.0) - g1 / ar * (s - ur)
            rho[n], pp[n] = rr * c**ex, pr * c ** (gamma * ex)
            u[n] = 2.0 / (gamma + 1.0) * (-ar + 0.5 * (gamma - 1.0) * ur + s)
    return rho, u, pp


def riemann_figure(out: Path, name: str, title: str, runs: list[tuple[str, Path, str]]) -> None:
    """Draw density, velocity and pressure of 1-D Euler Riemann runs against the exact solution."""
    ini = read_ini(next(runs[0][1].glob("*.ini")))
    gamma, t, ngc = gamma_of(ini), float(ini["time"]["time_max"]), int(ini["grid"]["ngc"])
    r1, r2 = ini["initial_conditions_region_1"], ini["initial_conditions_region_2"]
    left, right = (float(r1["r"]), float(r1["u"]), float(r1["p"])), (float(r2["r"]), float(r2["u"]), float(r2["p"]))
    xe = np.linspace(0.0, 1.0, 2001)
    ex = dict(zip(("r", "u", "p"), exact_riemann(xe, t, 0.5, gamma, left, right), strict=True))
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    for ax, key, label in zip(axs, ("r", "u", "p"), (r"$\rho$", r"$u$", r"$p$"), strict=True):
        ax.plot(xe, ex[key], "k-", lw=1.0, label="exact")
        for lab, work, mk in runs:
            axis, xs, q = profile(work, ngc, EULER)
            ax.plot(xs, primitives(q, gamma, axis)[key], mk, ms=3, mfc="none", label=lab)
        ax.set_xlabel("x")
        ax.set_title(label)
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.suptitle(title)
    fig.savefig(out / f"{name}.png", dpi=DPI)
    plt.close(fig)


def fig_sod(out: Path, runs: Path) -> None:
    """V1 / RV-3: Sod, 200 cells, split scheme and weno-riemann HLLC against the exact solution."""
    riemann_figure(out, "sod", "Sod shock tube, 200 cells, t = 0.2",
                   [("weno (split)", HERE / "sod" / f"{TAG}-x", "s"),
                    ("weno-riemann HLLC", RF / f"{TAG}-hllc-6th-weno-characteristic-sod-x", "o")])


def fig_lax(out: Path, runs: Path) -> None:
    """RV-3: Lax, 200 cells, split, LLF and HLLC against the exact solution."""
    riemann_figure(out, "lax", "Lax problem, 200 cells, t = 0.13",
                   [("weno (split)", RF / f"{TAG}-split-6th-weno-characteristic-lax-x", "s"),
                    ("weno-riemann LLF", RF / f"{TAG}-llf-6th-weno-characteristic-lax-x", "^"),
                    ("weno-riemann HLLC", RF / f"{TAG}-hllc-6th-weno-characteristic-lax-x", "o")])


def fig_shu_osher(out: Path, runs: Path) -> None:
    """RV-3: Shu-Osher, 400 cells, split and HLLC against the 16x split reference (cached by shu_osher_oracle)."""
    rho_ref = np.load(Path(tempfile.gettempdir()) / "flume-shu-osher-ref-6400-t1.8.npy")
    xr = -5.0 + (np.arange(rho_ref.size) + 0.5) * 10.0 / rho_ref.size
    fig, axs = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True, gridspec_kw={"width_ratios": [3, 2]})
    for ax, lim in zip(axs, ((-5.0, 5.0), (0.4, 2.6)), strict=True):
        ax.plot(xr, rho_ref, "k-", lw=0.8, label="reference (split, 6400 cells)")
        for lab, tag, mk in (("weno (split)", "split", "s"), ("weno-riemann HLLC", "hllc", "o")):
            _, xs, q = profile(RF / f"{TAG}-{tag}-6th-weno-characteristic-shu-osher-x", 3, EULER)
            ax.plot(xs, q[0], mk, ms=3, mfc="none", label=lab)
        ax.set_xlim(*lim)
        ax.set_xlabel("x")
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(r"$\rho$")
    axs[1].set_ylim(3.0, 4.9)
    axs[1].set_title("entropy waves behind the shock")
    axs[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Shu-Osher shock / entropy-wave interaction, 400 cells, t = 1.8")
    fig.savefig(out / "shu-osher.png", dpi=DPI)
    plt.close(fig)


def field_map(ax, bl: list[dict], values, cmap: str = "viridis", outline: bool = False):  # noqa: ANN001, ANN201
    """Draw a 2-D field (the first z plane of every block) and optionally the block outlines."""
    vals = [values(b) for b in bl]
    lo, hi = min(float(v.min()) for v in vals), max(float(v.max()) for v in vals)
    im = None
    for b, v in zip(bl, vals, strict=True):
        nx, ny = v.shape
        xe = b["lo"][0] + np.arange(nx + 1) * b["d"][0]
        ye = b["lo"][1] + np.arange(ny + 1) * b["d"][1]
        im = ax.pcolormesh(xe, ye, v.T, cmap=cmap, vmin=lo, vmax=hi, shading="flat")
        if outline:
            ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1] - xe[0], ye[-1] - ye[0], fill=False, lw=0.3, ec="w"))
    ax.set_aspect("equal")
    return im


def plane(b: dict, gamma: float, key: str) -> np.ndarray:
    """Return a derived 2-D field of a block (its first z plane): density, pressure, magnetic pressure or Mach."""
    f = {k: v[:, :, 0] for k, v in b["f"].items()}
    r = f["r"]
    kin = 0.5 * (f["ru"] ** 2 + f["rv"] ** 2 + f["rw"] ** 2) / r
    mag = 0.5 * (f["bx"] ** 2 + f["by"] ** 2 + f["bz"] ** 2) if "bx" in f else 0.0
    p = (gamma - 1.0) * (f["rE"] - kin - mag)
    return {"r": r, "p": p, "pb": mag, "mach": np.sqrt(2.0 * kin / r) / np.sqrt(gamma * p / r)}[key]


def first_plane(bl: list[dict]) -> list[dict]:
    """Return the blocks of the lowest z layer (2-D cases are stacked along the null z axis)."""
    zmin = min(b["lo"][2] for b in bl)
    return [b for b in bl if abs(b["lo"][2] - zmin) < 1e-12]


def fig_vortex(out: Path, runs: Path) -> None:
    """V2 / RV-2: isentropic vortex, density field and convergence of four schemes."""
    sys.path.insert(0, str(HERE / "vortex"))
    import vortex_oracle  # noqa: PLC0415

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    work = RF / f"{TAG}-hllc-6th-weno-characteristic-vortex-n128"
    ini = read_ini(next(work.glob("vortex-*.ini")))
    bl = first_plane(blocks(work, int(ini["grid"]["ngc"]), EULER))
    im = field_map(axs[0], bl, lambda b: b["f"]["r"][:, :, 0])
    fig.colorbar(im, ax=axs[0], label=r"$\rho$")
    axs[0].set_title("density at t = 0.1 (a tenth of the period), 128², HLLC")
    ladder = ("064", "128", "256")
    for lab, fmt, dirs in (("weno (split)", "s-", [HERE / "vortex" / f"{TAG}-n{n}" for n in ladder]),
                           ("weno-riemann LLF", "^-", [RF / f"{TAG}-llf-6th-weno-characteristic-vortex-n{n}"
                                                      for n in ladder]),
                           ("weno-riemann HLL", "d-", [RF / f"{TAG}-hll-6th-weno-characteristic-vortex-n{n}"
                                                      for n in ladder]),
                           ("weno-riemann HLLC", "o-", [RF / f"{TAG}-hllc-6th-weno-characteristic-vortex-n{n}"
                                                       for n in ladder])):
        try:
            res = [vortex_oracle.run_error(d) for d in dirs]
        except (StopIteration, SystemExit, FileNotFoundError):
            print(f"vortex: {lab} ladder incomplete, skipped")
            continue
        axs[1].loglog([r[0] for r in res], [r[1] for r in res], fmt, mfc="none", label=lab)
    n = np.array([64, 256])
    axs[1].loglog(n, 1.2e-4 * (64 / n) ** 5, "k:", label="order 5")
    axs[1].set_xlabel("cells per side")
    axs[1].set_ylabel(r"$L_1(\rho)$")
    axs[1].grid(alpha=0.3, which="both")
    axs[1].legend(fontsize=8)
    axs[1].set_title("convergence")
    fig.savefig(out / "vortex.png", dpi=DPI)
    plt.close(fig)


def fig_cylinder(out: Path, runs: Path) -> None:
    """V6: Mach 2 shock over a cylinder, immersed boundary and solid AMR: density with block outlines, Mach."""
    work = RF / f"{TAG}-hllc-6th-weno-characteristic-cylinder"
    ini = read_ini(next(work.glob("*.ini")))
    gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
    sol = ini["solid_1"]
    cx, cy, rad = (float(sol[k]) for k in ("circle_center_x", "circle_center_y", "circle_radius"))
    bl = first_plane(blocks(work, ngc, EULER))

    def masked(b: dict, key: str) -> np.ndarray:
        v = plane(b, gamma, key)
        nx, ny = v.shape
        xc = b["lo"][0] + (np.arange(nx) + 0.5) * b["d"][0]
        yc = b["lo"][1] + (np.arange(ny) + 0.5) * b["d"][1]
        xx, yy = np.meshgrid(xc, yc, indexing="ij")
        return np.ma.masked_where((xx - cx) ** 2 + (yy - cy) ** 2 < rad * rad, v)

    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    for ax, key, lab in zip(axs, ("r", "mach"), (r"$\rho$", "Mach number"), strict=True):
        im = field_map(ax, bl, lambda b, k=key: masked(b, k), cmap="magma" if key == "r" else "coolwarm",
                       outline=(key == "r"))
        ax.add_patch(Circle((cx, cy), rad, fill=True, fc="0.6", ec="k", lw=0.8))
        fig.colorbar(im, ax=ax, label=lab, shrink=0.85)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    axs[0].set_title("density, blocks outlined (refined ring on the surface)")
    axs[1].set_title("Mach number")
    fig.suptitle("Mach 2 shock over a cylinder: immersed boundary + solid AMR, weno-riemann HLLC")
    fig.savefig(out / "shock-cylinder.png", dpi=DPI)
    plt.close(fig)


def fig_step(out: Path, runs: Path) -> None:  # noqa: ARG001
    """EV-step: Woodward-Colella Mach 3 step at t = 4, three realms on quadtrees with static refined boxes: density
    with blocks, Mach number, and a zoom on the Mach stem (on the B-C seam, base grid) and the reflected shock crossing
    C's quadtree 2:1 seam at x = 0.9, against the uniform 1/80 control."""
    work = HERE / "step" / f"{TAG}-full"
    ini = read_ini(work / "step-A.ini")
    gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
    bl = first_plane(blocks(work, ngc, EULER))
    uniform = HERE / "step" / f"{TAG}-full-uniform"
    bu = first_plane(blocks(uniform, ngc, EULER)) if any(uniform.glob("*-proc*.h5")) else None
    fig = plt.figure(figsize=(12, 12.5), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(1.0, 1.0, 1.25))
    axs = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, :]), fig.add_subplot(grid[2, 0]),
           fig.add_subplot(grid[2, 1])]
    panels = ((axs[0], "r", (0.0, 3.0, 0.0, 1.0)), (axs[1], "mach", (0.0, 3.0, 0.0, 1.0)),
              (axs[2], "r", (0.4, 1.4, 0.2, 1.0)), (axs[3], "r" if bu else "p", (0.4, 1.4, 0.2, 1.0)))
    labels = {"r": r"$\rho$", "mach": "Mach number", "p": "p"}
    for ax, key, (x0, x1, y0, y1) in panels:
        cmap = {"r": "magma", "mach": "coolwarm", "p": "viridis"}[key]
        src = bu if (bu and ax is axs[3]) else bl
        im = field_map(ax, src, lambda b, k=key: plane(b, gamma, k), cmap=cmap, outline=(key != "mach"))
        ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2, fill=True, fc="0.6", ec="k", lw=0.8))
        ax.plot([0.0, 0.6], [0.2, 0.2], "w--", lw=0.8)  # seam A-B
        ax.plot([0.6, 0.6], [0.2, 1.0], "w--", lw=0.8)  # seam B-C
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        fig.colorbar(im, ax=ax, label=labels[key], shrink=0.9)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    axs[0].set_title("density, blocks outlined (static refined boxes), inter-realm seams dashed")
    axs[1].set_title("Mach number")
    axs[2].set_title("density, zoom: C refined from x = 0.9 (2:1 seam)")
    axs[3].set_title("density, same zoom, uniform 1/80 (no refined box)" if bu else "pressure, same zoom")
    fig.suptitle("Woodward-Colella Mach 3 step, t = 4: three realms on quadtrees, 1/80 base, 1/160 in static boxes")
    fig.savefig(out / "step.png", dpi=DPI)
    plt.close(fig)


def fig_step_trees(out: Path, runs: Path) -> None:  # noqa: ARG001
    """EV-step trees leg: the step forest at t = 0.5 on quadtrees and on octrees, and their column-wise difference."""
    works = {tree: HERE / "step" / f"{TAG}-trees-{tree}" for tree in ("quad", "oct")}
    ini = read_ini(works["quad"] / "step-A.ini")
    gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
    bl = {tree: blocks(w, ngc, EULER) for tree, w in works.items()}

    def key(b: dict, i: int, j: int) -> tuple[float, float]:
        return (round(float(b["lo"][0] + (i + 0.5) * b["d"][0]), 12),
                round(float(b["lo"][1] + (j + 0.5) * b["d"][1]), 12))

    cols = {}
    for tree in ("quad", "oct"):
        cols[tree] = {}
        for b in bl[tree]:
            q = np.stack([b["f"][v][:, :, 0] for v in EULER])
            for i in range(q.shape[1]):
                for j in range(q.shape[2]):
                    cols[tree].setdefault(key(b, i, j), q[:, i, j])
    scale = np.max(np.abs(np.array(list(cols["oct"].values()))), axis=0)
    scale[1:4] = scale[1:4].max()  # momentum as one vector
    diff = {k: float(np.max(np.abs(cols["quad"][k] - cols["oct"][k]) / scale)) for k in cols["quad"]}
    fig, axs = plt.subplots(3, 1, figsize=(11, 11.5), constrained_layout=True)
    for ax, tree, title in ((axs[0], "quad", "quadtree (nk = 1)"), (axs[1], "oct", "octree, null z (nk = 4, 16 in A)")):
        im = field_map(ax, first_plane(bl[tree]), lambda b: plane(b, gamma, "r"), cmap="magma", outline=True)
        fig.colorbar(im, ax=ax, label=r"$\rho$", shrink=0.9)
        ax.set_title(f"density, {title}, blocks outlined")
    floor = 1.0e-16
    im = None
    for b in bl["quad"]:
        nx, ny = b["f"]["r"].shape[:2]
        v = np.array([[max(diff[key(b, i, j)], floor) for j in range(ny)] for i in range(nx)])
        xe = b["lo"][0] + np.arange(nx + 1) * b["d"][0]
        ye = b["lo"][1] + np.arange(ny + 1) * b["d"][1]
        im = axs[2].pcolormesh(xe, ye, v.T, cmap="viridis", norm=mpl.colors.LogNorm(vmin=floor, vmax=1e-11),
                               shading="flat")
        axs[2].add_patch(Rectangle((xe[0], ye[0]), xe[-1] - xe[0], ye[-1] - ye[0], fill=False, lw=0.3, ec="w"))
    axs[2].set_aspect("equal")
    fig.colorbar(im, ax=axs[2], label="max relative difference (0 shown as 1e-16)", shrink=0.9)
    axs[2].set_title(f"quadtree against octree, per (x, y) column: max {max(diff.values()):.1e}")
    for ax in axs:
        ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2, fill=True, fc="0.6", ec="k", lw=0.8))
        ax.plot([0.0, 0.6], [0.2, 0.2], "c--", lw=0.8)
        ax.plot([0.6, 0.6], [0.2, 1.0], "c--", lw=0.8)
        ax.set_xlim(0.0, 3.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle("Woodward-Colella step, t = 0.5, 1/40 base: the three-realm forest on quadtrees and on octrees")
    fig.savefig(out / "step-trees.png", dpi=DPI)
    plt.close(fig)


def fig_conservation(out: Path, runs: Path) -> None:
    """V3 / RV-4: relative drift of the volume integrals of the periodic AMR box, with and without reflux."""
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    for reflux, ls in (("true", "-"), ("false", "--")):
        hist = next((HERE / "conservation" / f"{TAG}-reflux-{reflux}").glob("*-conservation_history.dat"))
        rows = [ln.split() for ln in hist.read_text().splitlines() if ln.split() and ln.split()[0][0] in "+-0123456789"]
        h = np.array([[float(v) for v in r] for r in rows])
        for c, name in zip(range(2, 7), ("mass", "x momentum", "y momentum", "z momentum", "energy"), strict=True):
            d = np.abs(h[:, c] - h[0, c]) / abs(h[0, c])
            ax.semilogy(h[:, 1], np.maximum(d, 1e-17), ls, label=f"{name}, reflux {reflux}")
    ax.set_xlabel("t")
    ax.set_ylabel("relative drift of the volume integral")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("periodic box with a refined octant: reflux on (solid) and off (dashed)")
    fig.savefig(out / "conservation.png", dpi=DPI)
    plt.close(fig)


def mhd_field_figure(out: Path, name: str, work: Path, panels: tuple, title: str, outline: bool = False) -> None:
    """Draw 2-D MHD fields of a run, one panel per (field key, label, colour map)."""
    ini = read_ini(next(work.glob("*.ini")))
    gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
    bl = first_plane(blocks(work, ngc, (*MHD, "psi")))
    fig, axs = plt.subplots(1, len(panels), figsize=(5.6 * len(panels), 4.8), constrained_layout=True, squeeze=False)
    for ax, (key, lab, cmap) in zip(axs[0], panels, strict=True):
        im = field_map(ax, bl, lambda b, k=key: plane(b, gamma, k), cmap=cmap, outline=outline)
        fig.colorbar(im, ax=ax, label=lab, shrink=0.85)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(lab)
    fig.suptitle(title)
    fig.savefig(out / f"{name}.png", dpi=DPI)
    plt.close(fig)


def fig_orszag_tang(out: Path, runs: Path) -> None:
    """MV-12: Orszag-Tang vortex, 128^2, density and pressure; the AMR variant with its blocks."""
    mhd_field_figure(out, "orszag-tang", runs / "orszag-tang", (("r", r"$\rho$", "magma"), ("p", "$p$", "viridis")),
                     "Orszag-Tang vortex, 128²")
    mhd_field_figure(out, "orszag-tang-amr", runs / "orszag-tang-amr", (("r", r"$\rho$", "magma"),),
                     "Orszag-Tang, 32² with the centre box refined (2:1 seams)", outline=True)


def fig_rotor(out: Path, runs: Path) -> None:
    """MV-13: MHD rotor, density and magnetic pressure."""
    mhd_field_figure(out, "rotor", runs / "rotor", (("r", r"$\rho$", "magma"), ("pb", r"$|B|^2/2$", "viridis")),
                     "MHD rotor")


def fig_field_loop(out: Path, runs: Path) -> None:
    """MV-9: advected field loop, magnetic pressure at the final time."""
    mhd_field_figure(out, "field-loop", runs / "field-loop", (("pb", r"$|B|^2/2$", "viridis"),),
                     "Field loop advected across the periodic box")


def mhd_profile_figure(out: Path, name: str, work: Path, title: str, exact=None,  # noqa: ANN001
                       label: str = "FLUME", extra: tuple[tuple[str, Path, str], ...] = ()) -> None:
    """Draw density, normal velocity, pressure and the tangential field of a 1-D MHD run (and of the `extra` runs)."""
    keys = (("r", r"$\rho$"), ("u", r"$u_n$"), ("p", "$p$"), ("bt1", r"$B_{t1}$"), ("bt2", r"$B_{t2}$"))
    fig, axs = plt.subplots(1, 5, figsize=(16, 3.4), constrained_layout=True)
    for n, (lab_run, w, marker) in enumerate(((label, work, "o"), *extra)):
        ini = read_ini(next(w.glob("*.ini")))
        gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
        axis, xs, q = profile(w, ngc, MHD)
        pr = primitives(q, gamma, axis)
        ex = exact(xs) if exact and n == 0 else None
        for ax, (k, _) in zip(axs, keys, strict=True):
            if ex is not None:
                ax.plot(xs, ex[k], "k-", lw=1.0, label="exact")
            ax.plot(xs, pr[k], marker, ms=2.5, mfc="none", label=lab_run)
    for ax, (_, lab) in zip(axs, keys, strict=True):
        ax.set_title(lab)
        ax.set_xlabel("x")
        ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.suptitle(title)
    fig.savefig(out / f"{name}.png", dpi=DPI)
    plt.close(fig)


def fig_mhd_riemann(out: Path, runs: Path) -> None:
    """MV-4 and the riemann tests: Ryu-Jones 2a against its exact solution, Brio-Wu, Ryu-Jones 4d."""
    sys.path.insert(0, str(HERE / "mhd" / "rj2a"))
    import rj2a_oracle  # noqa: PLC0415

    def rj2a_exact(xs: np.ndarray) -> dict[str, np.ndarray]:
        q = rj2a_oracle.exact(xs)
        kin = 0.5 * (q[1] ** 2 + q[2] ** 2 + q[3] ** 2) / q[0]
        mag = 0.5 * (q[5] ** 2 + q[6] ** 2 + q[7] ** 2)
        return {"r": q[0], "u": q[1] / q[0], "p": (rj2a_oracle.GAMMA - 1.0) * (q[4] - kin - mag), "bt1": q[6],
                "bt2": q[7]}

    mhd_profile_figure(out, "rj2a", runs / "rj2a", "Ryu-Jones 2a (all seven waves), 256 cells", rj2a_exact,
                       label="weno (split)", extra=(("weno-riemann HLLD", runs / "rj2a-hlld", "s"),))
    mhd_profile_figure(out, "brio-wu", runs / "brio-wu", "Brio-Wu shock tube, GLM")
    mhd_profile_figure(out, "rj4d", runs / "rj4d", "Ryu-Jones 4d")


def fig_glm_pulse(out: Path, runs: Path) -> None:
    """MV-3: the (B_x, psi) telegraph system, FLUME against the exact damped solution."""
    sys.path.insert(0, str(HERE / "mhd" / "glm-pulse"))
    import glm_pulse_oracle as g  # noqa: PLC0415

    work = runs / "glm-pulse"
    ini = g.read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), constrained_layout=True)
    for step in ("first", "last"):
        xs, prof = g.load(work, ngc, step)
        t = 0.0 if step == "first" else g.final_time(work)
        bx, psi = g.exact(ini, xs, t)
        for ax, k, e in zip(axs, ("bx", "psi"), (bx, psi), strict=True):
            ax.plot(xs, e, "k-", lw=1.0, label=f"exact, t = {t:g}")
            if k in prof:
                ax.plot(xs, prof[k], "o", ms=2.5, mfc="none", label=f"FLUME, t = {t:g}")
    axs[0].set_title("$B_x$")
    axs[1].set_title(r"$\psi$")
    for ax in axs:
        ax.set_xlabel("x")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle(r"GLM pulse (no damping): the divergence error splits into two pulses travelling at $\pm c_h$")
    fig.savefig(out / "glm-pulse.png", dpi=DPI)
    plt.close(fig)


def fig_blast(out: Path, runs: Path) -> None:
    """PV-1 and PV-2: the limited blasts, density and magnetic pressure on the two schemes."""
    cases = (("blast", r"Balsara-Spicer blast, $\beta = 2.5\cdot10^{-4}$, 128², t = 0.01, EGLM, positivity limiter"),
             ("blast-wushu", r"Wu-Shu blast, $\beta = 2.51\cdot10^{-6}$, 128², t = 0.001, EGLM, positivity limiter"))
    for name, title in cases:
        fig, axs = plt.subplots(2, 2, figsize=(11.2, 9.4), constrained_layout=True)
        for row, (scheme, lab) in zip(axs, (("split", "weno (split)"), ("hlld", "weno-riemann HLLD")), strict=True):
            work = runs / f"{name}-{scheme}"
            ini = read_ini(next(work.glob("*.ini")))
            gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
            bl = first_plane(blocks(work, ngc, (*MHD, "psi")))
            for ax, (key, sym, cmap) in zip(row, (("r", r"$\rho$", "magma"), ("pb", r"$|B|^2/2$", "viridis")),
                                            strict=True):
                im = field_map(ax, bl, lambda b, k=key, g=gamma: plane(b, g, k), cmap=cmap)
                fig.colorbar(im, ax=ax, shrink=0.85)
                ax.set_title(f"{sym}, {lab}")
                ax.set_xlabel("x")
                ax.set_ylabel("y")
        fig.suptitle(title)
        fig.savefig(out / f"{name}.png", dpi=DPI)
        plt.close(fig)


def double_rarefaction_exact(xs: np.ndarray, t: float, gamma: float, state: tuple) -> tuple:
    """Return density, velocity and pressure of two symmetric rarefactions leaving the origin (state: rho, |u|, p)."""
    r0, u0, p0 = state
    a0 = math.sqrt(gamma * p0 / r0)
    xi = np.abs(xs) / t  # the right half; the left one is its mirror image
    a = np.clip(2.0 / (gamma + 1.0) * (a0 - 0.5 * (gamma - 1.0) * (u0 - xi)), 0.0, a0)
    fan = xi < u0 + a0
    r = np.where(fan, r0 * (a / a0) ** (2.0 / (gamma - 1.0)), r0)
    p = np.where(fan, p0 * (a / a0) ** (2.0 * gamma / (gamma - 1.0)), p0)
    u = np.where(fan, xi - a, u0)
    return r, np.sign(xs) * u, p


def fig_vacuum(out: Path, runs: Path) -> None:
    """PV-3: LeBlanc, the double rarefaction and the planar Sedov blast, the two schemes with the limiter."""
    cases = (("leblanc", "LeBlanc shock tube, 400 cells, t = 6", True),
             ("double-rarefaction", "double rarefaction, 400 cells, t = 0.6", False),
             ("sedov", "planar Sedov blast, 800 cells, t = 0.001", True))
    fig, axs = plt.subplots(3, 3, figsize=(13.5, 10.2), constrained_layout=True)
    for row, (name, title, logy) in zip(axs, cases, strict=True):
        for scheme, lab, marker in (("split", "weno (split)", "s"), ("hllc", "weno-riemann HLLC", "o")):
            work = runs / f"{name}-{scheme}"
            ini = read_ini(next(work.glob("*.ini")))
            gamma, ngc = gamma_of(ini), int(ini["grid"]["ngc"])
            axis, xs, q = profile(work, ngc, EULER)
            pr = primitives(q, gamma, axis)
            if scheme == "split" and name != "sedov":
                t = float(ini["time"]["time_max"])
                reg = [ini[f"initial_conditions_region_{i}"] for i in (1, 2)]
                st = [(float(s["r"]), float(s["u"]), float(s["p"])) for s in reg]
                if name == "leblanc":
                    ex = exact_riemann(xs, t, float(reg[0]["emax_x"]), gamma, st[0], st[1])
                else:
                    ex = double_rarefaction_exact(xs, t, gamma, st[1])
                for ax, e in zip(row, ex, strict=True):
                    ax.plot(xs, e, "k-", lw=1.0, label="exact")
            for ax, k in zip(row, ("r", "u", "p"), strict=True):
                ax.plot(xs, pr[k], marker, ms=2.5, mfc="none", label=lab)
        for ax, sym in zip(row, (r"$\rho$", "$u$", "$p$"), strict=True):
            ax.set_title(f"{sym}, {title}", fontsize=9)
            ax.set_xlabel("x")
            ax.grid(alpha=0.3)
        if logy:
            row[0].set_yscale("log")
            row[2].set_yscale("log")
        row[0].legend(fontsize=8)
    fig.suptitle("Euler near vacuum with the positivity limiter (no floors)")
    fig.savefig(out / "near-vacuum.png", dpi=DPI)
    plt.close(fig)


def fig_eglm_energy(out: Path, runs: Path) -> None:  # noqa: ARG001
    """EV-4: drift of the total energy on the Orszag-Tang vortex, GLM and EGLM on the two schemes (history files)."""
    ot = HERE / "mhd" / "orszag-tang"
    curves = (("GLM, weno (split)", TAG, "C0", "-"), ("EGLM, weno (split)", f"{TAG}-eglm", "C1", "-"),
              ("GLM, weno-riemann HLLD", f"{TAG}-hlld-characteristic-6th-weno", "C0", "--"),
              ("EGLM, weno-riemann HLLD", f"{TAG}-hlld-characteristic-6th-weno-eglm", "C1", "--"))
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    for lab, d, colour, style in curves:
        h = np.loadtxt(ot / d / "orszag-tang-conservation_history.dat", skiprows=1)
        ax.plot(h[:, 1], (h[:, 6] - h[0, 6]) / h[0, 6], color=colour, ls=style, label=lab)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$(\int E - \int E_0) / \int E_0$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Orszag-Tang vortex, 128²: drift of the total energy")
    fig.savefig(out / "eglm-energy.png", dpi=DPI)
    plt.close(fig)


def fig_order(out: Path, runs: Path) -> None:
    """MV-5: convergence of the four linear-wave families (1-D), the two schemes."""
    sys.path.insert(0, str(HERE / "mhd" / "linear-wave"))
    import linear_wave_oracle  # noqa: PLC0415

    ladder = (16, 32, 64)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True, sharey=True)
    for ax, (scheme, lab) in zip(axs, (("split", "weno (split)"), ("hlld", "weno-riemann HLLD")), strict=True):
        for wave, marker in (("fast", "o"), ("alfven", "s"), ("slow", "^"), ("entropy", "d")):
            eps = []
            for n in ladder:
                work = runs / f"linear-wave-{scheme}-{wave}-{n}"
                eps.append(linear_wave_oracle.error(work, int(read_ini(next(work.glob("*.ini")))["grid"]["ngc"])))
            orders = " / ".join(f"{math.log2(a / b):.2f}" for a, b in zip(eps[:-1], eps[1:], strict=True))
            ax.loglog(ladder, eps, marker + "-", ms=5, mfc="none", label=f"{wave}, orders {orders}")
        ref = [eps[0] * (ladder[0] / n) ** 5 for n in ladder]
        ax.loglog(ladder, ref, "k:", lw=1.0, label=r"$N^{-5}$")
        ax.set_xticks(ladder, [str(n) for n in ladder])
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        ax.set_xlabel("cells per wavelength")
        ax.set_title(lab)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    axs[0].set_ylabel(r"$\epsilon$ after one period")
    fig.suptitle("MHD linear waves (Stone et al. 2008), amplitude $10^{-7}$: order of accuracy")
    fig.savefig(out / "linear-wave-order.png", dpi=DPI)
    plt.close(fig)


def fig_ghosts(out: Path, runs: Path) -> None:  # noqa: ARG001
    """GP: every face and edge ghost of the step forest on a linear field, its error before and after #65 P0.

    The "before" run is the same probe built from 3f51451d (the inward diagonal copy at the realm edges), kept as
    ghosts/work-before-p0-np2-step; the "after" run is `ghosts/check.sh --cases step --np 2 --keep`.
    """
    sys.path.insert(0, str(HERE / "ghosts"))
    import ghost_probe  # noqa: PLC0415

    cases = (("before #65 P0: realm edges held an inward diagonal copy", HERE / "ghosts" / "work-before-p0-np2-step"),
             ("after: realm edges compose the conditions of the faces", HERE / "ghosts" / f"{TAG}-step"))
    floor = 1.0e-16
    norm = mpl.colors.LogNorm(vmin=floor, vmax=1.0)
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.6), constrained_layout=True, gridspec_kw={"width_ratios": (2.4, 1)})
    sc = None
    for row, (title, work) in enumerate(cases):
        groups = ghost_probe.probe(work)
        pts = np.array([(c[0], c[1], max(err, floor)) for (cls, _), rows in groups.items()
                        for err, _, _, _, c in rows if cls in ("face", "edge")])
        order = np.argsort(pts[:, 2])  # largest errors drawn last, on top
        worst = pts[:, 2].max()
        views = ((-0.05, 3.05, -0.05, 1.05, 2.0), (0.45, 0.75, -0.04, 0.34, 14.0))  # x0, x1, y0, y1, marker size
        for col, (x0, x1, y0, y1, size) in enumerate(views):
            ax = axs[row, col]
            sc = ax.scatter(pts[order, 0], pts[order, 1], c=pts[order, 2], s=size, marker="s", cmap="inferno_r",
                            norm=norm, linewidths=0)
            ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2, fill=True, fc="0.85", ec="k", lw=0.8, zorder=0))
            for x, y in (((0.0, 0.6), (0.2, 0.2)), ((0.6, 0.6), (0.2, 1.0))):
                ax.plot(x, y, "c--", lw=0.8)  # inter-realm seams
            ax.add_patch(Rectangle((0.0, 0.0), 0.6, 1.0, fill=False, ec="k", lw=0.8))
            ax.add_patch(Rectangle((0.6, 0.2), 2.4, 0.8, fill=False, ec="k", lw=0.8))
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_aspect("equal")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        axs[row, 0].set_title(f"{title}: max {worst:.1e}")
        axs[row, 1].set_title("zoom: the step corner of realm A")
    fig.colorbar(sc, ax=axs, label="error of a face or edge ghost, relative to the field scale (0 shown as 1e-16)",
                 shrink=0.8)
    fig.suptitle("Ghost probe GP on the Woodward-Colella step forest (N = 80, refined boxes, 2 ranks): linear field")
    fig.savefig(out / "ghosts.png", dpi=DPI)
    plt.close(fig)


def fig_viscous(out: Path, runs: Path) -> None:  # noqa: ARG001
    """VV-1 to VV-5: convergence of the Navier-Stokes and Ohmic legs, and Becker's shock against the exact profile."""
    sys.path.insert(0, str(HERE / "viscous"))
    import profiles  # noqa: PLC0415
    import waves  # noqa: PLC0415

    vdir = HERE / "viscous"
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 4.6), constrained_layout=True)
    legs = (("VV-1 shear x, order 4", "o-", waves.run_error, "shear-x-o4", (32, 64, 128)),
            ("VV-1 shear x, order 2", "s-", waves.run_error, "shear-x-o2", (32, 64, 128)),
            ("VV-1 shear 45°, order 4", "d-", waves.run_error, "shear-xy-o4", (24, 48, 96)),
            ("VV-2 acoustic, order 4", "^-", waves.run_error, "acoustic-x-o4", (32, 64, 128)),
            ("VV-3 Couette (T), order 4", "v-", profiles.run_error, "couette-o4", (32, 64, 128)),
            ("VV-4 Becker M=2, order 4", "x-", profiles.run_error, "becker-m2-o4", (64, 128, 256)),
            ("VV-4 Becker M=3, order 4", "+-", profiles.run_error, "becker-m3-o4", (64, 128, 256)))
    for lab, fmt, fn, name, ns in legs:
        res = [fn(vdir / f"{TAG}-{name}-n{n}") for n in ns]
        err = [r[2] if name.startswith("couette") else r[1] for r in res]
        axs[0].loglog([r[0] for r in res], np.array(err) / err[0], fmt, mfc="none", label=lab)
    mhd = (("VV-5 Ohmic x, order 4", "o-", "ohmic-x-o4", (32, 64, 128)),
           ("VV-5 Ohmic x, order 2", "s-", "ohmic-x-o2", (32, 64, 128)),
           ("VV-5 Ohmic 45°, order 4", "d-", "ohmic-xy-o4", (24, 48, 96)),
           ("VV-5 Alfvén, $\\mu$ and $\\eta$, order 4", "^-", "alfven-x-o4", (32, 64, 128)),
           ("VV-5 acoustic on MHD, B = 0, order 4", "v-", "acoustic-mhd-x-o4", (32, 64, 128)))
    for lab, fmt, name, ns in mhd:
        res = [waves.run_error(vdir / f"{TAG}-{name}-n{n}") for n in ns]
        err = [r[1] for r in res]
        axs[1].loglog([r[0] for r in res], np.array(err) / err[0], fmt, mfc="none", label=lab)
    n = np.array([24.0, 256.0])
    for ax, title in ((axs[0], "Navier-Stokes convergence (CPU, 2 ranks)"),
                      (axs[1], "Ohmic (MHD) convergence (CPU, 2 ranks)")):
        ax.loglog(n, (24 / n) ** 4, "k:", label="order 4")
        ax.loglog(n, (24 / n) ** 2, "k--", label="order 2")
        ax.set_xlabel("cells along the wave or profile")
        ax.set_ylabel("RMS error / coarsest RMS error")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7)
        ax.set_title(title)
    for mach, color in ((2, "C0"), (3, "C3")):
        for n, mk in ((64, "o"), (256, ".")):
            work = vdir / f"{TAG}-becker-m{mach}-o4-n{n}"
            ini = read_ini(work / "input.ini")
            _, xs, q = profile(work, int(ini["grid"]["ngc"]), EULER)
            u = q[1] / q[0]
            mu = float(ini["physics"]["viscosity"])
            ue, u1, u2 = profiles.becker_u(xs, float(mach), mu)
            axs[2].plot(xs, (u - u2) / (u1 - u2), mk, color=color, ms=4 if n == 64 else 2, mfc="none",
                        label=f"M = {mach}, {n} cells")
        xx = np.linspace(-0.1, 0.1, 801)
        ue, u1, u2 = profiles.becker_u(xx, float(mach), mu)
        axs[2].plot(xx, (ue - u2) / (u1 - u2), "-", color=color, lw=0.9, label=f"M = {mach}, $\\mu$ = {mu}, Becker")
    axs[2].set_xlim(-0.06, 0.06)
    axs[2].set_xlabel("x")
    axs[2].set_ylabel(r"$(u - u_2)/(u_1 - u_2)$")
    axs[2].grid(alpha=0.3)
    axs[2].legend(fontsize=7)
    axs[2].set_title("VV-4: Becker's shock, Pr = 3/4")
    fig.savefig(out / "viscous.png", dpi=DPI)
    plt.close(fig)


def fig_viscous_fields(out: Path, runs: Path) -> None:  # noqa: ARG001
    """VV-1 and VV-8 in 2-D: the diagonal shear wave and its pointwise error on a uniform grid, then the error maps of
    the shear and Ohmic waves on the quadtree refined 2:1 in its centre (where the seam error lives)."""
    sys.path.insert(0, str(HERE / "viscous"))
    import waves  # noqa: PLC0415

    vdir = HERE / "viscous"

    def wave_blocks(name: str) -> tuple[list[dict], callable, callable]:
        work = vdir / f"{TAG}-{name}"
        ini = read_ini(work / "input.ini")
        mode, k, n, coef, _ = waves.linear_solution(ini)
        bl = first_plane(blocks(work, int(ini["grid"]["ngc"]), ("r", "ru", "rv", "bx", "by")))

        def measured(b: dict) -> np.ndarray:
            f = {key: v[:, :, 0] for key, v in b["f"].items()}
            if mode == "magnetic":
                return -n[1] * f["bx"] + n[0] * f["by"]
            return (-n[1] * f["ru"] + n[0] * f["rv"]) / f["r"]

        def error(b: dict) -> np.ndarray:
            nx, ny = b["f"]["r"].shape[:2]
            xc = b["lo"][0] + (np.arange(nx) + 0.5) * b["d"][0]
            yc = b["lo"][1] + (np.arange(ny) + 0.5) * b["d"][1]
            xx, yy = np.meshgrid(xc, yc, indexing="ij")
            exact = np.real(coef * np.exp(1j * k * (n[0] * xx + n[1] * yy)))
            return np.abs(measured(b) - exact)

        return bl, measured, error

    def error_map(ax, bl: list[dict], error, title: str) -> None:  # noqa: ANN001
        vals = [error(b) for b in bl]
        hi = max(float(v.max()) for v in vals)
        norm = mpl.colors.LogNorm(vmin=hi * 1.0e-4, vmax=hi)
        im = None
        for b, v in zip(bl, vals, strict=True):
            nx, ny = v.shape
            xe = b["lo"][0] + np.arange(nx + 1) * b["d"][0]
            ye = b["lo"][1] + np.arange(ny + 1) * b["d"][1]
            im = ax.pcolormesh(xe, ye, np.maximum(v, hi * 1.0e-4).T, cmap="viridis", norm=norm, shading="flat")
            ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1] - xe[0], ye[-1] - ye[0], fill=False, lw=0.4, ec="w"))
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, label="|error|", shrink=0.85)
        ax.set_title(title, fontsize=9)

    fig, axs = plt.subplots(2, 2, figsize=(11.5, 10.5), constrained_layout=True)
    bl, measured, error = wave_blocks("shear-xy-o4-n96")
    im = field_map(axs[0, 0], bl, measured, cmap="RdBu_r", outline=True)
    fig.colorbar(im, ax=axs[0, 0], label=r"$u_t$", shrink=0.85)
    axs[0, 0].set_title(r"VV-1 shear wave at 45$^\circ$, uniform 96$^2$, $\mu$ = 0.01, t = 0.1", fontsize=9)
    error_map(axs[0, 1], bl, error, "VV-1 pointwise error, uniform 96$^2$ (4th order, smooth)")
    bl, _, error = wave_blocks("seam-shear-o4-n48")
    error_map(axs[1, 0], bl, error, "VV-8 shear error, 48$^2$ with the centre refined 2:1 (96$^2$ cells)")
    bl, _, error = wave_blocks("seam-ohmic-o4-n48")
    error_map(axs[1, 1], bl, error, r"VV-8 Ohmic error ($\eta$ = 0.01, mhd-none), same grid")
    for ax in axs.flat:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.savefig(out / "viscous-fields.png", dpi=DPI)
    plt.close(fig)


def fig_couette(out: Path, runs: Path) -> None:  # noqa: ARG001
    """VV-3 in 2-D: the Couette channel (u and T fields, profiles against the exact ones) and where its error lives."""
    sys.path.insert(0, str(HERE / "viscous"))
    import profiles  # noqa: PLC0415

    vdir = HERE / "viscous"
    c = profiles.COUETTE
    k = c["mu"] * profiles.GAMMA / (profiles.GAMMA - 1.0) / c["prandtl"]

    def exact(y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return c["U"] * y / c["H"], c["Tw"] + c["mu"] * c["U"] ** 2 / (k * c["H"]) * (y - y**2 / (2.0 * c["H"]))

    def channel(n: int) -> list[dict]:
        work = vdir / f"{TAG}-couette-o4-n{n}"
        bl = first_plane(blocks(work, int(read_ini(work / "input.ini")["grid"]["ngc"]), ("r", "ru", "rE")))
        for b in bl:
            f = {key: v[:, :, 0] for key, v in b["f"].items()}
            u = f["ru"] / f["r"]
            ny = u.shape[1]
            y = np.broadcast_to(b["lo"][1] + (np.arange(ny) + 0.5) * b["d"][1], u.shape)
            ue, te = exact(y)
            b["u"], b["T"] = u / c["U"], (profiles.GAMMA - 1.0) * (f["rE"] - 0.5 * f["r"] * u**2) / f["r"] / c["Tw"]
            b["eu"], b["eT"], b["y"] = np.abs(b["u"] - ue / c["U"]), np.abs(b["T"] - te / c["Tw"]), y
        return bl

    def log_map(ax, bl: list[dict], key: str, label: str, title: str) -> None:  # noqa: ANN001
        hi = max(float(b[key].max()) for b in bl)
        norm = mpl.colors.LogNorm(vmin=hi * 1.0e-4, vmax=hi)
        im = None
        for b in bl:
            nx, ny = b[key].shape
            xe = b["lo"][0] + np.arange(nx + 1) * b["d"][0]
            ye = b["lo"][1] + np.arange(ny + 1) * b["d"][1]
            im = ax.pcolormesh(xe, ye, np.maximum(b[key], hi * 1.0e-4).T, cmap="viridis", norm=norm, shading="flat")
            ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1] - xe[0], ye[-1] - ye[0], fill=False, lw=0.3, ec="w"))
        fig.colorbar(im, ax=ax, label=label, shrink=0.85)
        ax.set_title(title, fontsize=9)

    fig, axs = plt.subplots(2, 3, figsize=(15.5, 9.5), constrained_layout=True)
    bl = channel(64)
    for ax, key, label, cmap in ((axs[0, 0], "u", "$u/U$", "viridis"), (axs[0, 1], "T", "$T/T_w$", "magma")):
        im = field_map(ax, bl, lambda b, key=key: b[key], cmap=cmap, outline=True)
        fig.colorbar(im, ax=ax, label=label, shrink=0.85)
        ax.set_title(f"{label}, 64 cells across the channel, t = 0.2", fontsize=9)
    log_map(axs[1, 0], bl, "eu", "$|u - u_e|/U$", "pointwise error of $u$ (64 cells)")
    log_map(axs[1, 1], bl, "eT", "$|T - T_e|/T_w$", "pointwise error of $T$ (64 cells)")
    for ax in (axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]):
        ax.set_aspect("equal")
        ax.set_xlabel("x (periodic)")
        ax.set_ylabel("y")
        ax.axhline(0.0, color="c", lw=3)
        ax.axhline(1.0, color="orange", lw=3)
    axs[0, 0].text(0.5, 0.03, "isothermal wall at rest", ha="center", color="c", fontsize=8, transform=axs[0, 0].transAxes)
    axs[0, 0].text(0.5, 0.93, "adiabatic wall moving at U", ha="center", color="orange", fontsize=8,
                   transform=axs[0, 0].transAxes)
    # profiles at 32 cells against the exact solution, and the x-averaged temperature error at every resolution
    yy = np.linspace(0.0, c["H"], 201)
    ue, te = exact(yy)
    ax = axs[0, 2]
    b32 = channel(32)
    ys = np.concatenate([b["y"][0] for b in b32 if b["lo"][0] < 1.0e-12])
    us = np.concatenate([b["u"][0] for b in b32 if b["lo"][0] < 1.0e-12])
    ts = np.concatenate([b["T"][0] for b in b32 if b["lo"][0] < 1.0e-12])
    ax.plot(ue / c["U"], yy, "C0-", lw=1, label="$u/U$ exact")
    ax.plot(us, ys, "C0o", ms=3, mfc="none", label="$u/U$, 32 cells")
    ax2 = ax.twiny()
    ax2.plot(te / c["Tw"], yy, "C3-", lw=1, label="$T/T_w$ exact")
    ax2.plot(ts, ys, "C3s", ms=3, mfc="none", label="$T/T_w$, 32 cells")
    ax.set_xlabel("$u/U$", color="C0")
    ax2.set_xlabel("$T/T_w$", color="C3")
    ax.set_ylabel("y")
    ax.legend(loc="upper left", fontsize=7)
    ax2.legend(loc="lower right", fontsize=7)
    ax = axs[1, 2]
    for n, fmt in ((32, "o-"), (64, "s-"), (128, "^-")):
        rows: dict[float, list[float]] = {}
        for b in channel(n):
            for y, e in zip(b["y"][0], b["eT"].mean(axis=0), strict=True):
                rows.setdefault(round(float(y), 12), []).append(float(e))
        y = np.array(sorted(rows))
        ax.semilogy(y, [np.mean(rows[v]) for v in y], fmt, ms=3, mfc="none", label=f"{n} cells")
    ax.set_xlabel("y")
    ax.set_ylabel(r"$|T - T_e|/T_w$, averaged along x")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    ax.set_title("the error enters at the isothermal wall (2nd-order mirror) and diffuses in;\n"
                 "the moving adiabatic wall's mirror is exact", fontsize=9)
    fig.suptitle("VV-3: compressible Couette flow, $\\mu$ = 0.01, Pr = 0.75, order-4 fluxes (CPU, 2 ranks)")
    fig.savefig(out / "couette.png", dpi=DPI)
    plt.close(fig)


def fig_regrid(out: Path, runs: Path) -> None:  # noqa: ARG001
    """AV (#74 P4): runtime regridding. The tracked vortex (base N = 32, finest N = 128) at t = 0.2 with its blocks and
    its pointwise density error against the exact solution, and the tracked field loop (GLM, N = 64 to 128) at t = 1
    with its blocks: the refined blocks follow the vortex core and the loop edge."""
    sys.path.insert(0, str(HERE / "vortex"))
    import vortex_oracle  # noqa: PLC0415

    vortex = HERE / "regrid" / "work-adam_flume_cpu-av2-tracked"
    loop = HERE / "regrid" / "work-adam_flume_cpu-av4-glm-tracked"
    ini_v, ini_l = read_ini(vortex / "vortex-av.ini"), read_ini(loop / "field-loop.ini")
    bl_v = first_plane(blocks(vortex, int(ini_v["grid"]["ngc"]), EULER))
    bl_l = first_plane(blocks(loop, int(ini_l["grid"]["ngc"]), MHD))
    oracle_ini = vortex_oracle.read_ini(vortex / "vortex-av.ini")

    def error(b: dict) -> np.ndarray:
        r = b["f"]["r"][:, :, 0]
        xc = b["lo"][0] + (np.arange(r.shape[0]) + 0.5) * b["d"][0]
        yc = b["lo"][1] + (np.arange(r.shape[1]) + 0.5) * b["d"][1]
        xx, yy = np.meshgrid(xc, yc, indexing="ij")
        return np.log10(np.abs(r - vortex_oracle.exact_density(xx, yy, oracle_ini)) + 1e-12)

    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    im = field_map(axs[0], bl_v, lambda b: b["f"]["r"][:, :, 0], outline=True)
    fig.colorbar(im, ax=axs[0], label=r"$\rho$", shrink=0.85)
    axs[0].set_title(f"tracked vortex, t = 0.2: density, {len(bl_v)} blocks (uniform 128: 64)")
    im = field_map(axs[1], bl_v, error, cmap="magma", outline=True)
    fig.colorbar(im, ax=axs[1], label=r"$\log_{10}|\rho - \rho_{exact}|$", shrink=0.85)
    axs[1].set_title("its pointwise density error")
    im = field_map(axs[2], bl_l, lambda b: plane(b, gamma_of(ini_l), "pb"), outline=True)
    fig.colorbar(im, ax=axs[2], label=r"$|B|^2/2$", shrink=0.85)
    axs[2].set_title(f"tracked field loop (GLM), t = 1: {len(bl_l)} blocks (uniform 128: 64)")
    for ax in axs:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle("Runtime regridding: the refined blocks follow the vortex and the loop edge (regrid every 5 steps)")
    fig.savefig(out / "regrid.png", dpi=DPI)
    plt.close(fig)


def gravity_oracle():  # noqa: ANN201
    """Import gravity/atmosphere_oracle.py (the GV checks, #73), whose measures the gravity figures reuse."""
    sys.path.insert(0, str(HERE / "gravity"))
    import atmosphere_oracle  # noqa: PLC0415

    return atmosphere_oracle


FLOOR = 1e-18  # an exact zero on a log axis


def fig_gravity_atmospheres(out: Path, runs: Path) -> None:
    """GV (G0): the isothermal and polytropic atmospheres of iso2d and poly2d (g = (0, -1), N = 32, CPU, 2 ranks) after
    one negligible step with poisoned ghosts: interior cells and the `hydrostatic` ghosts against the analytic profile,
    and their relative error. The ghosts beyond a polytropic face follow the isothermal extension of the boundary cell
    (the rule of the `hydrostatic` boundary), so that is their reference."""
    ao = gravity_oracle()
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for case, color in (("iso2d", "C0"), ("poly2d", "C3")):
        work = runs / case
        ini = ao.read_ini(work)
        atm = ao.Atmosphere(ini)
        ngc = int(ini["grid"]["ngc"])
        pts = {"interior": [], "ghost": []}
        for path in last_files(work):
            with h5py.File(path, "r") as h5:
                for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                    origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                    a = {n: h5[f"{blk}-{n}"][()].transpose(2, 1, 0) for n in EULER}
                    nx, ny = a["r"].shape[:2]
                    k = ngc if a["r"].shape[2] > 2 * ngc else 0
                    for i in range(ngc, nx - ngc):
                        for j in range(ny):
                            y = origin[1] + (j + 0.5) * d[1]
                            inside = ngc <= j < ny - ngc
                            if not inside and 0.0 < y < 1.0:
                                continue  # a ghost of a block interface
                            q = {n: float(v[i, j, k]) for n, v in a.items()}
                            rho = q["r"]
                            p = (atm.gamma - 1.0) * (q["rE"] - 0.5 * (q["ru"] ** 2 + q["rv"] ** 2) / rho)
                            c = np.array([origin[0] + (i + 0.5) * d[0], y, 0.0])
                            if inside:
                                re, pe = atm.state(c)
                            else:  # the isothermal extension of the boundary cell
                                c1 = c.copy()
                                c1[1] = 0.5 * d[1] if y < 0.0 else 1.0 - 0.5 * d[1]
                                r1, p1 = atm.state(c1)
                                re = r1 * math.exp(atm.g[1] * (y - c1[1]) * r1 / p1)
                                pe = re * p1 / r1
                            err = max(abs(rho / re - 1.0), abs(p / pe - 1.0))
                            pts["interior" if inside else "ghost"].append((y, rho, p, err))
        ys = np.linspace(0.0, 1.0, 200)
        prof = np.array([atm.state(np.array([0.0, y, 0.0])) for y in ys])
        label = "isothermal" if case == "iso2d" else f"polytropic G = {atm.G:g}"
        axs[0].plot(ys, prof[:, 0], "-", color=color, lw=1, label=rf"$\rho$ exact, {label}")
        axs[0].plot(ys, prof[:, 1], "--", color=color, lw=1, label=rf"$p$ exact, {label}")
        for kind, mk in (("interior", "."), ("ghost", "x")):
            v = np.array(pts[kind])
            axs[0].plot(v[:, 0], v[:, 1], mk, color=color, ms=3)
            axs[0].plot(v[:, 0], v[:, 2], mk, color=color, ms=3)
            axs[1].semilogy(v[:, 0], np.maximum(v[:, 3], FLOOR), mk, color=color, ms=4,
                            label=f"{case}, {kind} ({len(v)} cells), max {v[:, 3].max():.1e}")
    axs[0].axvspan(-0.2, 0.0, color="0.9")
    axs[0].axvspan(1.0, 1.2, color="0.9")
    axs[0].set_xlim(-0.14, 1.14)
    axs[0].set_xlabel("y (shaded: the hydrostatic ghosts)")
    axs[0].set_ylabel(r"$\rho$, $p$ (dots: interior, crosses: ghosts)")
    axs[0].legend(fontsize=7)
    axs[0].set_title("the atmospheres and their hydrostatic ghosts (beyond a polytropic face: the isothermal\n"
                     "extension of the boundary cell, the rule of the boundary)", fontsize=10)
    axs[1].set_xlabel("y")
    axs[1].set_ylabel(r"$\max(|\rho/\rho_e - 1|, |p/p_e - 1|)$")
    axs[1].grid(alpha=0.3, which="both")
    axs[1].legend(fontsize=7)
    axs[1].set_title("relative error against the analytic profile (exact zero drawn at 1e-18)")
    fig.suptitle("GV (G0): stratified atmospheres under g = (0, -1), N = 32, ghosts filled in RK-stage order "
                 "(CPU, 2 ranks)")
    fig.savefig(out / "gravity-atmospheres.png", dpi=DPI)
    plt.close(fig)


def fig_gravity_force(out: Path, runs: Path) -> None:
    """GV (G1, G2a): the force. Free fall in a periodic box (freefall2d: rho v against rho g t, and its error), the
    energy budget of the closed boxes (box2d, box2d-mhd: mass, E + rho Phi and W over time) and the box2d density after
    the gas has fallen onto the bottom wall."""
    ao = gravity_oracle()
    fig, axs = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True)
    work = runs / "freefall2d"
    ini = ao.read_ini(work)
    atm = ao.Atmosphere(ini)
    rows = []
    for _, (t, paths) in sorted(ao.checkpoints(work).items()):
        mom = [q["rv"] for _, _, q in ao.interior(ini, paths)]
        exact = atm.rho1 * atm.g[1] * t
        err = max(abs(m - exact) for m in mom) / abs(exact) if t > 0.0 else 0.0
        rows.append((t, float(np.mean(mom)), err))
    r = np.array(rows)
    ax = axs[0]
    ax.plot(r[:, 0], r[:, 1], "o", ms=5, mfc="none", label=r"FLUME $\rho v$ (mean over the cells)")
    ax.plot(r[:, 0], atm.rho1 * atm.g[1] * r[:, 0], "-", lw=1, label=r"exact $\rho g t$")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\rho v$")
    ax.legend(fontsize=7, loc="lower left")
    ax2 = ax.twinx()
    ax2.semilogy(r[1:, 0], np.maximum(r[1:, 2], FLOOR), "s-", color="C3", ms=3, lw=0.8)
    ax2.set_ylabel(r"$\max |\rho v - \rho g t| / |\rho g t|$", color="C3")
    ax.set_title("GV-1 free fall, every face periodic")
    ax = axs[1]
    for case, ls in (("box2d", "-"), ("box2d-mhd", "--")):
        e = np.array(ao.energy_rows(runs / case))
        t, m, ee, w, b = e[:, 1], e[:, 2], e[:, 3], e[:, 4], e[:, 5]
        tot, wtot = ee + w, ee + w + b
        for val, name, color in ((m, "mass", "C0"), (tot, r"$E + \rho\Phi$", "C1"), (wtot, "W", "C2")):
            ax.semilogy(t, np.maximum(np.abs(val - val[0]) / abs(val[0]), FLOOR), ls, color=color,
                        label=f"{name}, {case}")
    ax.set_xlabel("t")
    ax.set_ylabel("relative drift")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2)
    ax.set_title("GV-3 closed boxes: W is the invariant, E + ρΦ differs by\nthe O(Δx²) wall-row term of the "
                 "fourth-order source")
    work = runs / "box2d"
    ini = ao.read_ini(work)
    bl = first_plane(blocks(work, int(ini["grid"]["ngc"]), EULER))
    im = field_map(axs[2], bl, lambda b: b["f"]["r"][:, :, 0], outline=True)
    fig.colorbar(im, ax=axs[2], label=r"$\rho$", shrink=0.85)
    t_end = max(t for t, _ in ao.checkpoints(work).values())
    axs[2].set_title(f"box2d at t = {t_end:.3f}: the gas has settled into a stratification")
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("y")
    fig.suptitle("GV (G1): uniform gas under g = (0, -1), N = 32, CFL 0.4 (CPU, 2 ranks)")
    fig.savefig(out / "gravity-force.png", dpi=DPI)
    plt.close(fig)


def fig_gravity_order(out: Path, runs: Path) -> None:
    """GV-9 (G2a): the falling density wave of wave2d at N = 32, 64, 128, L1 errors of rho and p against the exact
    solution on `weno` and on `weno-riemann` (HLLC), with slopes 5 and 2 (the order of the G1 face form alone)."""
    ao = gravity_oracle()
    fig, ax = plt.subplots(figsize=(7.5, 5), constrained_layout=True)
    for scheme, label, color in (("weno", "weno", "C0"), ("riemann", "weno-riemann HLLC", "C3")):
        e = np.array([(4 * n, *ao.wave_errors(runs / f"wave-{scheme}-{n}")[1:]) for n in (8, 16, 32)])
        for col, name, ls in ((1, r"$\rho$", "o-"), (2, r"$p/p_1$", "s--")):
            order = math.log2(e[-2, col] / e[-1, col])
            ax.loglog(e[:, 0], e[:, col], ls, color=color, mfc="none",
                      label=f"{name}, {label}: finest order {order:.2f}")
    n = np.array([32.0, 128.0])
    for slope, y0 in ((5, 4e-6), (2, 4e-6)):
        ax.loglog(n, y0 * (n / 32.0) ** -slope, ":", color="0.4", lw=1)
        ax.annotate(f"slope {slope}", (n[-1], y0 * 4.0 ** -slope), fontsize=8, color="0.4")
    ax.set_xlabel("N (cells per side)")
    ax.set_ylabel("L1 error at t = 0.2")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    ax.set_title("GV-9: a density wave falling freely, g = (-0.6, -0.8)\nwith the fourth-order energy source "
                 "(the G1 face form alone gave 2.12 / 2.00)")
    fig.savefig(out / "gravity-order.png", dpi=DPI)
    plt.close(fig)


def fig_gravity_balance(out: Path, runs: Path) -> None:
    """GV (G2b): the well-balanced mode. max |u| / c over time, unbalanced and balanced, on strat2d (isothermal, ten
    scale heights), polystrat2d and strat2d-amr; the spurious Mach number of strat2d-amr at t = 0.5 with its blocks,
    unbalanced and balanced; the response to the 1e-6 and 1e-4 pulses along x = 0.5 at t = 0.1, balanced and not; the
    balanced 1e-6 response field."""
    ao = gravity_oracle()
    fig, axs = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    ax = axs[0, 0]
    for case, color in (("strat2d", "C0"), ("polystrat2d", "C1"), ("strat2d-amr", "C2")):
        for mode, ls in (("unb", "-"), ("bal", "--")):
            work = runs / f"{case}-{mode}"
            ini = ao.read_ini(work)
            m = np.array([(t, ao.max_mach(ini, paths)) for t, paths in
                          (v for _, v in sorted(ao.checkpoints(work).items()))])
            ax.semilogy(m[:, 0], np.maximum(m[:, 1], FLOOR), ls, color=color,
                        label=f"{case}, {'balanced' if mode == 'bal' else 'unbalanced'}: {m[-1, 1]:.1e} at the end")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\max |u| / c$  (exact zero at 1e-18)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    ax.set_title("atmospheres at rest, N = 32")
    for ax, mode in ((axs[0, 1], "unb"), (axs[0, 2], "bal")):
        work = runs / f"strat2d-amr-{mode}"
        ini = ao.read_ini(work)
        gamma = float(ini["physics"]["gamma"])
        bl = first_plane(blocks(work, int(ini["grid"]["ngc"]), EULER))
        im = field_map(ax, bl, lambda b, g=gamma: np.log10(np.maximum(plane(b, g, "mach"), FLOOR)), cmap="magma",
                       outline=True)
        fig.colorbar(im, ax=ax, label=r"$\log_{10} |u| / c$", shrink=0.85)
        ax.set_title(f"strat2d-amr {'balanced' if mode == 'bal' else 'unbalanced'}, t = 0.5: refined box "
                     "[0.25, 0.75]²")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    for ax, mode in ((axs[1, 0], "bal"), (axs[1, 1], "unb")):
        for a, ls in (("6", "-"), ("4", "--")):
            work = runs / f"lin-{mode}-{a}"
            cps = ao.checkpoints(work)
            amp, resp = ao.pulse_response(work, cps[max(cps)][1])
            xs = sorted({k[0] for k in resp})
            x0 = min(xs, key=lambda x: abs(x - 0.5))
            cut = sorted((k[1], v[0]) for k, v in resp.items() if k[0] == x0)
            c = np.array(cut)
            ax.plot(c[:, 0], c[:, 1], ls, label=f"A = {amp:g}")
        ax.set_xlabel(f"y (cut x = {x0:.3f})")
        ax.set_ylabel(r"$(p/p_e - 1)/A$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        ax.set_title("pulse response at t = 0.1, " + ("balanced: the two overlap" if mode == "bal" else
                     "unbalanced: the stationary truncation flow over A"))
    ax = axs[1, 2]
    work = runs / "lin-bal-6"
    cps = ao.checkpoints(work)
    amp, resp = ao.pulse_response(work, cps[max(cps)][1])
    xs = np.array(sorted({k[0] for k in resp}))
    ys = np.array(sorted({k[1] for k in resp}))
    z0 = next(iter(resp))[2]
    grid = np.array([[resp[(x, y, z0)][0] for y in ys] for x in xs])
    pc = ax.pcolormesh(np.r_[xs - 0.5 * (xs[1] - xs[0]), xs[-1] + 0.5 * (xs[1] - xs[0])],
                       np.r_[ys - 0.5 * (ys[1] - ys[0]), ys[-1] + 0.5 * (ys[1] - ys[0])], grid.T, cmap="RdBu_r",
                       shading="flat")
    fig.colorbar(pc, ax=ax, label=r"$(p/p_e - 1)/A$", shrink=0.85)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"balanced, A = {amp:g}, t = 0.1: the acoustic response")
    fig.suptitle("GV (G2b): the well-balanced mode, [gravity] equilibrium = atmosphere (CPU, 2 ranks)")
    fig.savefig(out / "gravity-balance.png", dpi=DPI)
    plt.close(fig)


def fig_gravity_rt_linear(out: Path, runs: Path) -> None:
    """GV-4 (G3): the linear Rayleigh-Taylor stage of rt2d at N = 32 and 64 against the compressible linear theory of
    rt_oracle.py: the mode amplitude, its local growth rate, and its vertical profile against the eigenfunction."""
    ao = gravity_oracle()
    from rt_oracle import column_of, eigenfunction, growth_rate  # noqa: PLC0415

    fig, axs = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True)
    works = [(runs / f"rt-{c}", f"N = {4 * c}", m) for c, m in ((8, "o"), (16, "s"))]
    col, k = column_of(works[0][0] / "input.ini")
    sigma = growth_rate(col, k)
    for work, label, mk in works:
        s = np.array(ao.mode_amplitude(work))
        fit = np.polyfit(s[(s[:, 0] >= 2.0), 0], np.log(s[(s[:, 0] >= 2.0), 1]), 1)[0]
        axs[0].semilogy(s[:, 0], s[:, 1], mk, ms=4, mfc="none", label=f"{label}: fitted rate {fit:.5f}")
        rate = np.diff(np.log(s[:, 1])) / np.diff(s[:, 0])
        axs[1].plot(0.5 * (s[1:, 0] + s[:-1, 0]), rate, mk + "-", ms=3, mfc="none", lw=0.8, label=label)
        ini = ao.read_ini(work)
        cps = ao.checkpoints(work)
        rows, _ = ao.mode_rows(ini, cps[max(cps)][1])
        ys = np.array(list(rows))
        vs = np.abs(np.array(list(rows.values())))
        axs[2].plot(vs / vs.max(), ys, mk, ms=4, mfc="none", label=f"{label}, t = {cps[max(cps)][0]:.1f}")
    t = np.linspace(2.0, 4.0, 2)
    a4 = np.array(ao.mode_amplitude(works[-1][0]))[-1]
    axs[0].semilogy(t, a4[1] * np.exp(sigma * (t - a4[0])), "k-", lw=1, label=f"linear theory {sigma:.5f}")
    axs[0].axvspan(2.0, 4.0, color="0.92")
    axs[0].set_xlabel("t")
    axs[0].set_ylabel(r"mode amplitude $a$ of $v$")
    axs[0].legend(fontsize=7)
    axs[0].set_title("the seeded mode (shaded: the fit window of GV-4)")
    axs[1].axhline(sigma, color="k", lw=1, label="linear theory")
    axs[1].axhspan(sigma * (1 - 1e-3), sigma * (1 + 1e-3), color="0.85", label=r"$\pm 10^{-3}$")
    axs[1].set_ylim(0.97 * sigma, 1.015 * sigma)
    axs[1].set_xlabel("t")
    axs[1].set_ylabel(r"$d \ln a / dt$ between checkpoints")
    axs[1].legend(fontsize=7)
    axs[1].set_title("local growth rate: the solenoidal seed settles on the mode")
    h, v = eigenfunction(col, k, sigma)
    axs[2].plot(np.abs(v), h, "k-", lw=1, label="eigenfunction |v|")
    axs[2].axhline(col.h0, color="0.6", lw=0.6, ls="--")
    axs[2].set_xlabel(r"$|v_k| / \max |v_k|$")
    axs[2].set_ylabel("y")
    axs[2].legend(fontsize=7)
    axs[2].set_title("vertical profile of the mode at t = 4")
    fig.suptitle("GV-4: Rayleigh-Taylor, two isothermal layers 2:1, interface width 0.05, g = 1, wavelength 0.5 "
                 "(CPU, 2 ranks)")
    fig.savefig(out / "gravity-rt-linear.png", dpi=DPI)
    plt.close(fig)


def fig_gravity_rt(out: Path, runs: Path) -> None:
    """G3: the nonlinear single-mode Rayleigh-Taylor (rt2d-nl, interface width 0.01, seed 1e-2, N = 128) without and
    with a horizontal field at half the critical field of the mode (rt2d-mhd, EGLM): density at four times."""
    ao = gravity_oracle()
    times = (2.0, 3.0, 4.0, 5.0)
    fig, axs = plt.subplots(2, len(times), figsize=(4.2 * len(times), 8.6), constrained_layout=True)
    for row, (case, names, label) in enumerate((("rt2d-nl", EULER, "hydro"),
                                                ("rt2d-mhd", MHD, r"MHD, $B_x = 0.1$ (half critical)"))):
        work = runs / case
        ini = ao.read_ini(work)
        ngc = int(ini["grid"]["ngc"])
        cps = ao.checkpoints(work)
        for ax, t in zip(axs[row], times, strict=True):
            step = min(cps, key=lambda s: abs(cps[s][0] - t))
            bl = []
            for path in cps[step][1]:
                with h5py.File(path, "r") as h5:
                    for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                        origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                        a = h5[f"{blk}-r"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, :]
                        bl.append({"lo": origin + ngc * d, "d": d, "f": {"r": a[:, :, :1]}})
            im = field_map(ax, bl, lambda b: b["f"]["r"][:, :, 0], cmap="RdYlBu_r")
            ax.set_title(f"{label}, t = {cps[step][0]:.2f}", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(im, ax=axs[row, -1], label=r"$\rho$", shrink=0.8)
    fig.suptitle("Rayleigh-Taylor single mode: layers 2:1, interface width 0.01, g = 1, wavelength 0.5, seed 1e-2, "
                 "N = 128 (CPU, 2 ranks)")
    fig.savefig(out / "gravity-rt.png", dpi=DPI)
    plt.close(fig)


RT_SYM_CLEANINGS = (("rt2d-nl", "hydro", "k", "-"), ("rt-sym-none", "MHD, no cleaning", "tab:gray", "-"),
                    ("rt-sym-glm", "GLM, $c_h$ = 3", "tab:green", "-"),
                    ("rt-sym-ch1.5", "EGLM, $c_h$ = 1.5", "tab:blue", "--"),
                    ("rt2d-mhd", "EGLM, $c_h$ = 3", "tab:blue", "-"),
                    ("rt-sym-ch6", "EGLM, $c_h$ = 6", "tab:blue", ":"))
RT_SYM_DAMPING = (("rt-sym-a0", r"$\alpha$ = 0", "tab:red"), ("rt2d-mhd", r"$\alpha$ = 0.18 (default)", "tab:blue"),
                  ("rt-sym-a1", r"$\alpha$ = 1", "tab:orange"), ("rt-sym-a5", r"$\alpha$ = 5", "tab:purple"),
                  ("rt-sym-mincell", r"$\alpha$ = 0.18, $L = \Delta x$", "tab:brown"))


def fig_gravity_rt_symmetry(out: Path, runs: Path) -> None:
    """Issue #83: the mirror asymmetry of the magnetised Rayleigh-Taylor run (rt2d-mhd, N = 128) under each divergence
    control (left) and each EGLM damping (centre), and its growth rate over t in [2, 4] against the cleaning speed,
    with the largest growth rate of the linearised cleaned equations on a magnetised shear layer (right)."""
    ao = gravity_oracle()
    from glm_linear import max_growth  # noqa: PLC0415 (gravity/, on the path once gravity_oracle ran)

    series = {case: ao.mirror_asymmetry(runs / case)
              for case in {c for c, *_ in RT_SYM_CLEANINGS} | {c for c, *_ in RT_SYM_DAMPING}}

    def rate(case: str) -> float:
        fit = [(t, np.log(a)) for t, a in series[case] if 2.0 <= t <= 4.0 and a > 0.0]
        return float(np.polyfit(*zip(*fit, strict=True), 1)[0])

    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    for case, label, color, style in RT_SYM_CLEANINGS:
        t, a = zip(*series[case], strict=True)
        axs[0].semilogy(t, np.maximum(a, FLOOR), color=color, ls=style, label=f"{label} ({rate(case):.2f})")
    axs[0].set_title("divergence control (EGLM damping $\\alpha$ = 0.18, $L$ = 1)", fontsize=10)
    t, a = zip(*series["rt-sym-none"], strict=True)
    axs[1].semilogy(t, np.maximum(a, FLOOR), color="tab:gray", label=f"no cleaning ({rate('rt-sym-none'):.2f})")
    for case, label, color in RT_SYM_DAMPING:
        t, a = zip(*series[case], strict=True)
        axs[1].semilogy(t, np.maximum(a, FLOOR), color=color, label=f"{label} ({rate(case):.2f})")
    axs[1].set_title("EGLM, $c_h$ = 3: damping rate $\\alpha c_h / L$", fontsize=10)
    for ax in axs[:2]:
        ax.axvspan(2.0, 4.0, color="0.92", zorder=0)
        ax.set_xlabel("t")
        ax.set_ylabel(r"max $|\rho(x, y) - \rho(\lambda - x, y)|$")
        ax.set_ylim(1e-15, 1e-1)
        ax.legend(fontsize=8, title="rate over t in [2, 4]", title_fontsize=8)
    chs = np.array([0.75, 1.5, 3.0, 6.0])
    measured = [(1.5, "rt-sym-ch1.5"), (3.0, "rt2d-mhd"), (6.0, "rt-sym-ch6")]
    axs[2].plot([c for c, _ in measured], [rate(r) for _, r in measured], "o-", color="tab:blue",
                label="FLUME, EGLM: asymmetry rate")
    axs[2].plot([3.0], [rate("rt-sym-glm")], "s", color="tab:green", label="FLUME, GLM: asymmetry rate")
    axs[2].axhline(rate("rt-sym-none"), color="tab:gray", ls="--", label="FLUME, no cleaning")
    axs[2].axhline(rate("rt2d-nl"), color="k", ls="--", label="FLUME, hydro")
    for model, marker in (("eglm", "^"), ("glm", "v")):
        axs[2].plot(chs, [max_growth(model, 128, 2.0 * np.pi, c) for c in chs], marker, ls=":", color="tab:red",
                    label=f"linearised {model.upper()} PDE, shear layer, $k = 2\\pi$")
    axs[2].set_xlabel(r"cleaning speed $c_h$")
    axs[2].set_ylabel("growth rate")
    axs[2].set_ylim(0.0, None)
    axs[2].legend(fontsize=8)
    axs[2].set_title("the growth rises with $c_h$; the cleaned PDE's does not", fontsize=10)
    fig.suptitle("Magnetised Rayleigh-Taylor ($B_x$ = 0.1, N = 128): the mirror symmetry under each cleaning "
                 "(MHD on FNL, hydro on CPU, 2 ranks)")
    fig.savefig(out / "gravity-rt-symmetry.png", dpi=DPI)
    plt.close(fig)


def fig_ranks(out: Path, runs: Path) -> None:
    """DV-1 (issue #82): the runtime-AMR case of ranks/check.sh on 2 ranks, a checkpoint every step. Left and centre,
    the speed with the block outlines coloured by owner rank before and after the first regrid; right, the blocks of
    each rank against the step, and the np 2 - np 1 difference of the state at every step (bitwise: zero)."""
    ao = gravity_oracle()
    work = runs / "ranks-amr-np2"
    ini = ao.read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    cps, cps1 = ao.checkpoints(work), ao.checkpoints(runs / "ranks-amr-np1")
    colors = ("tab:blue", "tab:orange")

    def blocks(step: int) -> list[dict]:
        bl = []
        for path in cps[step][1]:
            rank = int(path.stem.rsplit("-proc", 1)[1])
            with h5py.File(path, "r") as h5:
                for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                    origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                    f = {n: h5[f"{blk}-{n}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, :1] for n in ("r", "ru", "rv")}
                    bl.append({"lo": origin + ngc * d, "d": d, "f": f, "rank": rank})
        return bl

    steps = sorted(cps)
    regrid = next(s for s in steps if any(n for r, n in _rank_blocks(cps[s][1]).items() if r > 0))
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    for ax, step in zip(axs[:2], (regrid - 1, steps[-1]), strict=True):
        bl = blocks(step)
        im = field_map(ax, bl, lambda b: np.hypot(b["f"]["ru"], b["f"]["rv"])[:, :, 0] / b["f"]["r"][:, :, 0],
                       cmap="Greys")
        for b in bl:
            nx, ny = b["f"]["r"].shape[:2]
            ax.add_patch(Rectangle((b["lo"][0], b["lo"][1]), nx * b["d"][0], ny * b["d"][1], fill=False, lw=1.6,
                                   ec=colors[b["rank"]]))
        count = {r: sum(1 for b in bl if b["rank"] == r) for r in (0, 1)}
        ax.set_title(f"step {step}, t = {cps[step][0]:.4f}: rank 0 {count[0]} blocks, rank 1 {count[1]}", fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, label="|u|", shrink=0.8)
    owned = np.array([[_rank_blocks(cps[s][1]).get(r, 0) for r in (0, 1)] for s in steps])
    for r, (style, width) in zip((0, 1), (("-", 3.5), ("--", 1.8)), strict=True):
        axs[2].step(steps, owned[:, r], where="post", color=colors[r], ls=style, lw=width, label=f"blocks of rank {r}")
    diffs = []
    for s in steps:
        a = {tuple(c.round(12)): q for c, _, q in ao.interior(ini, cps1[s][1])}
        b = {tuple(c.round(12)): q for c, _, q in ao.interior(ini, cps[s][1])}
        diffs.append(max(abs(a[c][n] - b[c][n]) for c in a for n in a[c]))
    axs[2].text(0.97, 0.5, f"max |np 2 - np 1| over every field,\ncell and step ({len(steps)} checkpoints): "
                f"{max(diffs):g}", transform=axs[2].transAxes, ha="right", va="center", fontsize=9,
                bbox={"fc": "w", "ec": "0.6"})
    axs[2].set_xlabel("step")
    axs[2].set_ylabel("blocks")
    axs[2].set_ylim(-0.5, 9.5)
    axs[2].legend(loc="upper left", fontsize=8)
    axs[2].set_title("rank 1 owns no block until the first regrid", fontsize=10)
    fig.suptitle("DV-1: EGLM atmosphere with a pressure pulse, 4 blocks (one sibling family) on 2 ranks, regrids every "
                 "5 steps (FNL, 2 A30 GPUs)")
    fig.savefig(out / "ranks.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_potential(out: Path, runs: Path) -> None:
    """SV-1 (issue #81, M8-P1): the closed box (M7 box2d, g = (0, -1), x periodic, y walls) to t = 0.2 through the
    constant-g path (fourth-order energy source) and through the potential path (`[gravity] potential = linear`, the
    second-order face form of self-gravity). Left, N = 64: the drift of the conserved energy of each path (W for the
    constant one, which adds the O(dx^2) wall term of its fourth-order source; E + rho Phi for the potential one) and of
    E alone, which exchanges with rho Phi. Right: the falling wave (GV-9, exact solution) to t = 0.2, the L1 errors of
    rho and p against N through each path."""
    ao = gravity_oracle()
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    for pot, color in (("uniform", "tab:blue"), ("linear", "tab:red")):
        rows = ao.energy_rows(runs / f"sv1-box2d-{pot}")
        t = np.array([r[1] for r in rows])
        e, w, b = (np.array([r[k] for r in rows]) for k in (3, 4, 5))
        total = e + w + b
        label = "constant g: W = E + rho Phi + wall term" if pot == "uniform" else "potential: E + rho Phi"
        axs[0].semilogy(t[1:], np.maximum(np.abs(total[1:] / total[0] - 1.0), FLOOR), "o-", color=color, ms=3,
                        label=label)
        if pot == "uniform":
            axs[0].semilogy(t[1:], np.maximum(np.abs((e + w)[1:] / (e + w)[0] - 1.0), FLOOR), ":", color=color,
                            label="constant g: E + rho Phi without the wall term")
        axs[0].semilogy(t[1:], np.abs(e[1:] / total[0] - e[0] / total[0]), "--", color=color, alpha=0.6,
                        label=f"{'constant g' if pot == 'uniform' else 'potential'}: E alone")
    axs[0].set_xlabel("t")
    axs[0].set_ylabel("relative drift")
    axs[0].set_ylim(1e-17, 1e-1)
    axs[0].legend(fontsize=8)
    axs[0].set_title("energy budget of the closed box, N = 64", fontsize=10)
    ns = np.array([32, 64, 128], dtype=float)
    for pot, color, name in (("uniform", "tab:blue", "constant g"), ("linear", "tab:red", "potential")):
        errs = np.array([ao.wave_errors(runs / f"sv1-wave-{pot}-{cells}")[1:] for cells in (8, 16, 32)])
        for k, (var, style) in enumerate((("rho", "o-"), ("p", "s--"))):
            order = np.log2(errs[-2, k] / errs[-1, k])
            axs[1].loglog(ns, errs[:, k], style, color=color, label=f"{name}: L1 {var}, order {order:.2f}")
    for slope in (2, 4):
        axs[1].loglog(ns, 1e-5 * (ns[0] / ns) ** slope, ":", color="0.6")
        axs[1].annotate(f"dx^{slope}", (ns[-1], 1e-5 * (ns[0] / ns[-1]) ** slope), fontsize=8, color="0.4")
    axs[1].set_xlabel("N")
    axs[1].set_ylabel("L1 error at t = 0.2")
    axs[1].legend(fontsize=8)
    axs[1].set_title("the falling wave: the potential path is second order", fontsize=10)
    fig.suptitle("SV-1: gravity through the cell potential Phi = -g . x against the constant g (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-potential.png", dpi=DPI)
    plt.close(fig)

def fig_selfgravity_poisson(out: Path, runs: Path) -> None:
    """SV-2 (issue #81, M8-P2): the FFT Poisson solver on the anisotropic periodic box [0,1] x [0,2] x [0,3] for
    rho = rho0 + sin(k1 x) cos(k1 y) cos(k1 z) + 0.5 cos(k2 . x) (test_poisson_fft, CPU, 2 ranks). Left: against N, the
    max error of Phi against the continuum solution (second order, 3-D octree and 2-D quadtree), against the discrete
    solution and the discrete residual (round-off). Middle and right: the 2-D potential at N = 64 and its error against
    the continuum solution."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import poisson_oracle  # noqa: PLC0415

    work = runs / "sv2-poisson"
    rows: dict[str, list[tuple[int, float, float, float]]] = {}
    pattern = re.compile(r"(3-D octree|2-D quadtree|3-D one block),\s+N = (\d+): \|Phi - Phi_h\| (\S+), "
                         r"residual (\S+), \|Phi - Phi_exact\| (\S+),")
    for line in (work / "log.txt").read_text().splitlines():
        m = pattern.search(line)
        if m:
            rows.setdefault(m[1], []).append((int(m[2]), float(m[3]), float(m[4]), float(m[5])))
    fig, axs = plt.subplots(1, 3, figsize=(14.0, 4.8), constrained_layout=True, width_ratios=(1.6, 1.0, 1.0))
    for name, color in (("3-D octree", "tab:blue"), ("2-D quadtree", "tab:red")):
        n, disc, res, cont = (np.array(v) for v in zip(*rows[name], strict=True))
        order = math.log(cont[-2] / cont[-1]) / math.log(n[-1] / n[-2])
        axs[0].loglog(n, cont, "o-", color=color, label=f"{name}: |Phi - Phi_exact|, order {order:.2f}")
        axs[0].loglog(n, disc, "s--", color=color, alpha=0.7, label=f"{name}: |Phi - Phi_h| / max|Phi_h|")
        axs[0].loglog(n, res, "^:", color=color, alpha=0.7, label=f"{name}: residual of L_h, relative")
    ns = np.array([16.0, 64.0])
    axs[0].loglog(ns, 1.5e-3 * (ns[0] / ns) ** 2, "-", color="0.6", lw=0.8)
    axs[0].annotate("dx^2", (32.0, 1.5e-3 * (ns[0] / 32.0) ** 2), xytext=(0, -12), textcoords="offset points",
                    fontsize=8, color="0.4")
    axs[0].set_xlabel("N (cells along x)")
    axs[0].set_ylabel("max error")
    axs[0].set_ylim(1e-16, 1e-1)
    axs[0].legend(fontsize=7)
    axs[0].set_title("second order against the PDE, round-off against L_h", fontsize=10)
    phi = poisson_oracle.load(work / "phi-r4-nk1-l3-p0.bin")[:, :, 0]
    nx, ny = phi.shape
    x = (np.arange(nx) + 0.5) / nx
    y = 2.0 * (np.arange(ny) + 0.5) / ny
    xx, yy = np.meshgrid(x, y, indexing="ij")
    k1 = 2.0 * np.pi * np.array([1.0, 1.0])
    k2 = 2.0 * np.pi * np.array([2.0, 0.5])
    exact = 4.0 * np.pi * (np.sin(k1[0] * xx) * np.cos(k1[1] * yy) / -(k1 @ k1)
                           + 0.5 * np.cos(k2[0] * xx + k2[1] * yy) / -(k2 @ k2))
    extent = (0.0, 1.0, 0.0, 2.0)
    im = axs[1].imshow(phi.T, origin="lower", extent=extent, cmap="RdBu_r")
    fig.colorbar(im, ax=axs[1], shrink=0.85)
    axs[1].set_title("2-D, N = 64: Phi", fontsize=10)
    im = axs[2].imshow((phi - exact).T, origin="lower", extent=extent, cmap="PuOr")
    fig.colorbar(im, ax=axs[2], shrink=0.85)
    axs[2].set_title(f"Phi - Phi_exact (max {np.max(np.abs(phi - exact)):.1e})", fontsize=10)
    for ax in axs[1:]:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    fig.suptitle("SV-2: the FFT Poisson solver, L_h Phi = 4 pi (rho - mean rho) on a periodic box (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-poisson.png", dpi=DPI)
    plt.close(fig)

def fig_selfgravity_isolated(out: Path, runs: Path) -> None:
    """SV-6 (issue #81, M8-P4): the FFT Poisson solver with isolated boundaries on the cube [-1, 1]^3 (test_poisson_fft
    --isolated, CPU, 2 ranks), two spheres of radius 1/2, G = 1. Left: against N, the max error of Phi against the
    exact potential over every cell and ghost (the n = 1 polytrope and the uniform sphere) and the residual of L_h.
    Middle: Phi of the polytrope at N = 64 along the line through the centre, over the whole padded array (the domain
    and its free-space continuation on both sides), against the exact Phi (-G M / r outside). Right: the error on the
    plane through the centre."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import poisson_oracle  # noqa: PLC0415

    work = runs / "sv6-poisson"
    rows: dict[str, list[tuple[int, float, float]]] = {}
    pattern = re.compile(r"(polytrope n=1|uniform sphere),\s+N = (\d+): \|Phi - Phi_h\| \S+, residual (\S+), "
                         r"\|Phi - Phi_exact\| (\S+),")
    for line in (work / "log.txt").read_text().splitlines():
        m = pattern.search(line)
        if m:
            rows.setdefault(m[1], []).append((int(m[2]), float(m[3]), float(m[4])))
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True, width_ratios=(1.3, 1.3, 1.0))
    for name, color in (("polytrope n=1", "tab:blue"), ("uniform sphere", "tab:red")):
        n, res, cont = (np.array(v) for v in zip(*rows[name], strict=True))
        order = -np.polyfit(np.log(n), np.log(cont), 1)[0]
        axs[0].loglog(n, cont, "o-", color=color, label=f"{name}: |Phi - Phi_exact|, least-squares order {order:.2f}")
        axs[0].loglog(n, res, "^:", color=color, alpha=0.7, label=f"{name}: residual of L_h, relative")
    ns = np.array([16.0, 64.0])
    axs[0].loglog(ns, 1.0e-2 * (ns[0] / ns) ** 2, "-", color="0.6", lw=0.8)
    axs[0].annotate("dx^2", (32.0, 1.0e-2 * (ns[0] / 32.0) ** 2), xytext=(0, -12), textcoords="offset points",
                    fontsize=8, color="0.4")
    axs[0].set_xlabel("N (cells per side)")
    axs[0].set_ylabel("max error")
    axs[0].set_ylim(1e-15, 1e-1)
    axs[0].legend(fontsize=7)
    axs[0].set_title("second order against the PDE; the residual is the kernel's", fontsize=10)
    phi = poisson_oracle.load(work / "phi-r8-nk8-l3-p1.bin")
    m, n, h = phi.shape[0], 64, 2.0 / 64
    a, kk = 0.5, math.pi / 0.5
    mass = 4.0 * a**3 / math.pi

    def exact(r: np.ndarray) -> np.ndarray:
        kr = np.maximum(kk * r, 1e-300)
        inside = -mass / a - 4.0 * math.pi / kk**2 * np.sin(kr) / kr
        return np.where(r < a, inside, -mass / np.maximum(r, 1e-300))

    # the padded array holds the free-space potential at the cell offsets -(m - n)/2 .. n - 1 + (m - n)/2 (wrapped)
    t = np.arange(m)
    t = np.where(t >= n + (m - n) // 2, t - m, t)
    order_t = np.argsort(t)
    x = -1.0 + (t[order_t] + 0.5) * h
    line = phi[order_t, n // 2, n // 2]
    yc = -1.0 + (n // 2 + 0.5) * h
    r_line = np.sqrt(x**2 + 2.0 * yc**2)
    axs[1].plot(x, line, ".", ms=3, color="tab:blue", label="Phi, the padded array (N = 64)")
    axs[1].plot(x, exact(r_line), "-", color="0.3", lw=0.8, label="exact, -G M / r outside")
    axs[1].axvspan(-1.0, 1.0, color="0.92", zorder=0, label="the domain")
    axs[1].set_xlabel("x (y = z = h/2)")
    axs[1].set_ylabel("Phi")
    axs[1].legend(fontsize=7)
    axs[1].set_title("the polytrope's potential and its free-space continuation", fontsize=10)
    xs = -1.0 + (np.arange(n) + 0.5) * h
    xx, yy = np.meshgrid(xs, xs, indexing="ij")
    err = phi[:n, :n, n // 2] - exact(np.sqrt(xx**2 + yy**2 + yc**2))
    im = axs[2].imshow(err.T, origin="lower", extent=(-1.0, 1.0, -1.0, 1.0), cmap="PuOr")
    fig.colorbar(im, ax=axs[2], shrink=0.85)
    axs[2].set_title(f"Phi - Phi_exact, z = h/2 (max {np.max(np.abs(err)):.1e})", fontsize=10)
    axs[2].set_xlabel("x")
    axs[2].set_ylabel("y")
    fig.suptitle("SV-6: the FFT Poisson solver with isolated boundaries, zero padding and the lattice Green's function "
                 "(CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-isolated.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_star(out: Path, runs: Path) -> None:
    """SV-7 and SV-8 (issue #81, M8-P4): the n = 1 polytrope star of selfgravity/make_polytrope.py (radius 1/2,
    rho_c = 1, G = 1, gamma = 2) in the closed box [-1, 1]^3 under isolated self-gravity (CPU, 2 ranks). Top: SV-7, the
    equilibrium in a light atmosphere with the positivity limiter, at 32^3 and 64^3: the central density against its
    initial value, and the largest Mach number inside 0.8 a. Bottom left: SV-8, the collapse (pressure halved) at
    32^3, the drift of E + W at CFL 0.4, 0.2, 0.1, with E and W alone. Bottom right: its density on the plane through
    the centre at the densest checkpoint."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415
    import star_oracle  # noqa: PLC0415

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
    for lev, color in ((2, "tab:blue"), (3, "tab:red")):
        rows = star_oracle.star_rows(runs / f"sv7-l{lev}")
        t = np.array([r[0] for r in rows])
        rc = np.array([r[1] for r in rows])
        mach = np.array([r[3] for r in rows])
        n = 8 * 2**lev
        axs[0, 0].plot(t, rc / rc[0] - 1.0, "-", color=color, label=f"{n}^3")
        axs[0, 1].plot(t, mach, "-", color=color, label=f"{n}^3")
    axs[0, 0].axhline(0.0, color="0.6", lw=0.6)
    axs[0, 0].set_xlabel("t")
    axs[0, 0].set_ylabel("central density / initial - 1")
    axs[0, 0].legend(fontsize=8)
    axs[0, 0].set_title("SV-7: the equilibrium star, its core", fontsize=10)
    axs[0, 1].set_xlabel("t")
    axs[0, 1].set_ylabel("largest Mach number inside 0.8 a")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].set_title("SV-7: the spurious motion inside the star", fontsize=10)
    ax = axs[1, 0]
    for cfl, color in (("0.4", "tab:blue"), ("0.2", "tab:orange"), ("0.1", "tab:green")):
        h = jeans_oracle.history(runs / f"sv8-collapse-{cfl}")
        ew0 = h["E_plus_W"][0]
        ax.semilogy(h["time"][1:], np.maximum(np.abs(h["E_plus_W"][1:] - ew0) / abs(ew0), FLOOR), "-", color=color,
                    label=f"E + W, CFL {cfl}")
        if cfl == "0.1":
            ax.semilogy(h["time"][1:], np.abs(h["int_rE"][1:] - h["int_rE"][0]) / abs(ew0), "--", color="0.4",
                        label="E alone")
            ax.semilogy(h["time"][1:], np.abs(h["W"][1:] - h["W"][0]) / abs(ew0), ":", color="0.4", label="W alone")
    ax.set_xlabel("t")
    ax.set_ylabel("relative change")
    ax.legend(fontsize=8)
    ax.set_title("SV-8: the collapse, E + W drifts at the time integration's order", fontsize=10)
    ax = axs[1, 1]
    col = runs / "sv8-collapse-0.1"
    rows = star_oracle.star_rows(col)
    k = int(np.argmax([r[1] for r in rows]))
    ini = star_oracle.read_ini(col)
    cps = star_oracle.checkpoints(col)
    tc, paths = cps[sorted(cps)[k]]
    cells = [(c, q["r"]) for c, dxyz, _, q in star_oracle.interior_steps(ini, paths) if abs(c[2]) < 0.6 * dxyz[2]]
    xs = sorted({round(float(c[0]), 10) for c, _ in cells})
    ys = sorted({round(float(c[1]), 10) for c, _ in cells})
    img = np.zeros((len(xs), len(ys)))
    ix = {x: i for i, x in enumerate(xs)}
    iy = {y: j for j, y in enumerate(ys)}
    for c, r in cells:
        img[ix[round(float(c[0]), 10)], iy[round(float(c[1]), 10)]] = r
    im = ax.imshow(img.T, origin="lower", extent=(-1.0, 1.0, -1.0, 1.0), cmap="magma")
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"SV-8: the density at t = {tc:.2f}, z = h/2 (max {img.max():.3f})", fontsize=10)
    fig.suptitle("SV-7, SV-8: the n = 1 polytrope under isolated self-gravity in a closed box (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-star.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_star_amr(out: Path, runs: Path) -> None:
    """SV-14 (issue #87, P3): the n = 1 polytrope of SV-7 and SV-8 on AMR, isolated self-gravity by the multigrid's
    symmetric James step, the star of radius 0.3 so that the 2:1 seam of the box [-1/2, 1/2]^3 lies in the atmosphere
    (CPU, 2 ranks). Top: the equilibrium with the box refined once over the 32^3 base (the star at the 64^3 spacing)
    against the same star on the uniform 32^3 grid: the central density against its initial value, and the largest
    Mach number inside 0.8 a. Bottom left: the
    collapse on that grid with the per-stage reflux, the drift of E + W at CFL 0.4, 0.2, 0.1, with E and W alone.
    Bottom right: the equilibrium star at the end on the plane z = h/2, block by block, its radius dashed: the 2:1
    seam of the box lies in the atmosphere."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415
    import star_oracle  # noqa: PLC0415

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
    eq = [("sv14-uniform", "uniform 32^3", "tab:blue", "-"),
          ("sv14-box", "AMR, 32^3 + the box [-1/2, 1/2]^3 at 64^3", "tab:green", "-")]
    for name, label, color, ls in eq:
        rows = star_oracle.star_rows(runs / name)
        t = np.array([r[0] for r in rows])
        rc = np.array([r[1] for r in rows])
        axs[0, 0].plot(t, rc / rc[0] - 1.0, ls, color=color, label=label)
        axs[0, 1].plot(t, [r[3] for r in rows], ls, color=color, label=label)
    axs[0, 0].axhline(0.0, color="0.6", lw=0.6)
    axs[0, 0].set_xlabel("t")
    axs[0, 0].set_ylabel("central density / initial - 1")
    axs[0, 0].legend(fontsize=8)
    axs[0, 0].set_title("the equilibrium star on AMR, its core", fontsize=10)
    axs[0, 1].set_xlabel("t")
    axs[0, 1].set_ylabel("largest Mach number inside 0.8 a")
    axs[0, 1].legend(fontsize=8)
    axs[0, 1].set_title("the spurious motion inside the star", fontsize=10)
    ax = axs[1, 0]
    for cfl, color in (("0.4", "tab:blue"), ("0.2", "tab:orange"), ("0.1", "tab:green")):
        h = jeans_oracle.history(runs / f"sv14-collapse-{cfl}")
        ew0 = h["E_plus_W"][0]
        ax.semilogy(h["time"][1:], np.maximum(np.abs(h["E_plus_W"][1:] - ew0) / abs(ew0), FLOOR), "-", color=color,
                    label=f"E + W, CFL {cfl}")
        if cfl == "0.1":
            ax.semilogy(h["time"][1:], np.abs(h["int_rE"][1:] - h["int_rE"][0]) / abs(ew0), "--", color="0.4",
                        label="E alone")
            ax.semilogy(h["time"][1:], np.abs(h["W"][1:] - h["W"][0]) / abs(ew0), ":", color="0.4", label="W alone")
    ax.set_xlabel("t")
    ax.set_ylabel("relative change")
    ax.legend(fontsize=8)
    ax.set_title("the collapse on AMR, E + W drifts at the time integration's order", fontsize=10)
    ax = axs[1, 1]
    work = runs / "sv14-box"
    ini = jeans_oracle.read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    cps = jeans_oracle.checkpoints(work)
    tc, paths = cps[max(cps)]
    rmax, blocks = -math.inf, []
    for path in paths:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                r = h5[f"{blk}-r"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc]
                x0, y0, z0 = (origin + ngc * dxyz)[:3]
                k = int(np.floor((0.0 - z0) / dxyz[2]))  # the cell layer just above z = 0
                if not 0 <= k < r.shape[2]:
                    continue
                rmax = max(rmax, float(r[:, :, k].max()))
                blocks.append((x0, y0, r.shape[0] * dxyz[0], r.shape[1] * dxyz[1], r[:, :, k]))
    for x0, y0, wx, wy, r in blocks:
        im = ax.imshow(r.T, origin="lower", extent=(x0, x0 + wx, y0, y0 + wy), cmap="magma", vmin=0.0, vmax=rmax)
        ax.add_patch(Rectangle((x0, y0), wx, wy, fill=False, lw=0.4, ec="w", alpha=0.6))
    ax.add_patch(Circle((0.0, 0.0), float(ini["initial_conditions"]["star_radius"].split(";")[0]), fill=False, lw=0.8,
                        ec="c", ls="--"))
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.0, 1.0)
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"the AMR star at t = {tc:.2f}, z = h/2: the seam 0.2 off the surface (dashed)", fontsize=10)
    fig.suptitle("SV-14: an n = 1 polytrope (radius 0.3) on AMR, isolated self-gravity by the multigrid's James step "
                 "(CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-star-amr.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_jeans(out: Path, runs: Path) -> None:
    """SV-3 and SV-4 (issue #81, M8-P3): the Jeans modes of selfgravity/make_jeans.py (rho_0 = 1, c = 1,
    k = 2 pi (1, 1), periodic, CPU, 2 ranks). Top left: SV-3, the L1 error of the stable wave (lambda = lambda_J / 2)
    after one period against N, 2-D and 3-D. Top right: SV-4, the amplitude of the unstable mode (lambda = 1.5
    lambda_J, N = 64) against the exact e^(sigma t). Bottom left: the collapse (amplitude 1e-2, N = 64), the drift of
    E + W at CFL 0.4, 0.2, 0.1, with E and W alone. Bottom right: the collapse density at t = 0.4."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415

    fig, axs = plt.subplots(2, 2, figsize=(12.5, 9.0), constrained_layout=True)
    ax = axs[0, 0]
    for label, tag, cells, side, color, mk in (("2-D", "2d", (8, 16, 32), 4, "tab:blue", "o"),
                                               ("3-D", "3d", (8, 16), 2, "tab:red", "s")):
        works = [runs / f"sv3-{tag}-{c}" for c in cells]
        n = np.array([side * c for c in cells], dtype=float)
        err = np.array([jeans_oracle.wave_error(w)[1] for w in works])
        order = math.log2(err[-2] / err[-1])
        ax.loglog(n, err, f"{mk}-", color=color, mfc="none", ms=8,
                  label=f"{label}: L1 rho / (rho_0 A), order {order:.2f}")
    ns = np.array([16.0, 128.0])
    ax.loglog(ns, 2e-3 * (32.0 / ns) ** 2, "-", color="0.6", lw=0.8)
    ax.annotate("dx^2", (64.0, 2e-3 * (32.0 / 64.0) ** 2), xytext=(0, -12), textcoords="offset points", fontsize=8,
                color="0.4")
    ax.set_xlabel("N (cells per side)")
    ax.set_ylabel("L1 error after one period")
    ax.legend(fontsize=8)
    ax.set_title("SV-3: the stable Jeans wave, w = pi sqrt(6)", fontsize=10)
    ax = axs[0, 1]
    lin = runs / "sv4-linear"
    jm = jeans_oracle.Jeans(jeans_oracle.read_ini(lin))
    rows = jeans_oracle.mode_amplitudes(lin)
    t = np.array([r[0] for r in rows])
    a = np.array([r[1] for r in rows])
    sel = (t >= 0.1) & (t <= 0.5)
    sigma = float(np.polyfit(t[sel], np.log(a[sel]), 1)[0])
    ax.semilogy(t, a, "o", ms=3, color="tab:blue", label=f"measured, fit sigma = {sigma:.4f}")
    ax.semilogy(t, jm.amp * np.exp(jm.rate * t), "-", color="0.3", label=f"exact, sigma = pi sqrt(10) = {jm.rate:.4f}")
    ax.set_xlabel("t")
    ax.set_ylabel("mode amplitude")
    ax.legend(fontsize=8)
    ax.set_title("SV-4: the unstable mode, lambda = 1.5 lambda_J", fontsize=10)
    ax = axs[1, 0]
    for cfl, color in (("0.4", "tab:blue"), ("0.2", "tab:orange"), ("0.1", "tab:green")):
        h = jeans_oracle.history(runs / f"sv4-collapse-{cfl}")
        ew0 = h["E_plus_W"][0]
        ax.semilogy(h["time"][1:], np.maximum(np.abs(h["E_plus_W"][1:] - ew0) / abs(ew0), FLOOR), "-", color=color,
                    label=f"E + W, CFL {cfl}")
        if cfl == "0.1":
            ax.semilogy(h["time"][1:], np.abs(h["int_rE"][1:] - h["int_rE"][0]) / abs(ew0), "--", color="0.4",
                        label="E alone")
            ax.semilogy(h["time"][1:], np.abs(h["W"][1:] - h["W"][0]) / abs(ew0), ":", color="0.4", label="W alone")
    ax.set_xlabel("t")
    ax.set_ylabel("relative change")
    ax.legend(fontsize=8)
    ax.set_title("SV-4: the collapse, E + W drifts at the time integration's order", fontsize=10)
    ax = axs[1, 1]
    col = runs / "sv4-collapse-0.1"
    ini = jeans_oracle.read_ini(col)
    cps = jeans_oracle.checkpoints(col)
    tc, paths = cps[max(cps)]
    cells = [(c, q["r"]) for c, _, q in jeans_oracle.interior(ini, paths)]
    xs = sorted({round(float(c[0]), 10) for c, _ in cells})
    ys = sorted({round(float(c[1]), 10) for c, _ in cells})
    img = np.zeros((len(xs), len(ys)))
    ix = {x: i for i, x in enumerate(xs)}
    iy = {y: j for j, y in enumerate(ys)}
    for c, r in cells:
        img[ix[round(float(c[0]), 10)], iy[round(float(c[1]), 10)]] = r
    lo = float(ini["grid"]["emin_x"])
    im = ax.imshow(img.T, origin="lower", extent=(lo, lo + 1.0, lo, lo + 1.0), cmap="magma")
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"the collapse density at t = {tc:.2f} (max {img.max():.3f})", fontsize=10)
    fig.suptitle("SV-3, SV-4: Jeans modes under the self-gravity of the FFT solver (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-jeans.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_amr(out: Path, runs: Path) -> None:
    """SV-9 (issue #81, M8-P5a): the composite Laplacian across 2:1 faces and the multigrid solver on it. Left and
    middle: the 2-D model of selfgravity/seam_operators.py (doubly periodic, coarse x < 1/2, fine x > 1/2, two seams),
    N = 8 to 48, three seam couplings: Tomida & Stone's flux matching (B), its symmetrisation (S1) and the mimetic mean
    flux (M, the one adam_poisson_mg_object assembles); left the max truncation error, middle the seam face gradient
    error (each operator's own seam flux) and the potential error. Right: the multigrid solver on the octree and the
    quadtree of test_poisson_mg (CPU, 2 ranks; the central box refined once, base levels 2 and 3): the potential
    error and the 2:1 face gravity error against the fine cells, the growth of the max truncation error."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import seam_operators  # noqa: PLC0415

    pattern = re.compile(r"SV-9b, AMR (\S+) \S+,\s+base level (\d+), leaves \d+ \(levels \d+-\d+\): iterations (\d+), "
                         r"factor \S+, symmetry (\S+), truncation (\S+), \|Phi - Phi_exact\| (\S+), "
                         r"2:1 face gravity error (\S+),")
    amr: dict[str, list[tuple[int, int, float, float, float, float]]] = {}
    for line in (runs / "sv9-poisson" / "log.txt").read_text(errors="replace").splitlines():
        m = pattern.search(line)
        if m:
            amr.setdefault(m[1], []).append((int(m[2]), int(m[3]), float(m[4]), float(m[5]), float(m[6]),
                                             float(m[7])))
    if not amr:
        raise FileNotFoundError("no SV-9b line in sv9-poisson/log.txt")
    rows = seam_operators.ladder([8, 16, 32, 48])
    style = {"B": ("tab:orange", "s", "B, flux matching (not symmetric)"),
             "S1": ("tab:red", "^", "S1 = (B + B^T) / 2 (symmetric, inconsistent)"),
             "M": ("tab:blue", "o", "M, mimetic mean flux (symmetric, consistent)")}
    fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.8), constrained_layout=True)
    for kind, (color, marker, label) in style.items():
        n, e_phi, tau, e_face, _ = (np.array(v) for v in zip(*rows[kind], strict=True))
        axs[0].loglog(n, tau, marker + "-", color=color, label=label)
        axs[1].loglog(n, e_face, marker + "-", color=color, label=f"{kind}: seam face gradient")
        axs[1].loglog(n, e_phi, marker + ":", color=color, alpha=0.6, label=f"{kind}: potential")
    ns = np.array([8.0, 48.0])
    axs[0].loglog(ns, 30.0 * ns / ns[0], "-", color="0.6", lw=0.8)
    axs[0].annotate("1/h", (24.0, 30.0 * 3.0), xytext=(-18, 2), textcoords="offset points", fontsize=8, color="0.4")
    axs[0].set_xlabel("N (coarse cells per side)")
    axs[0].set_ylabel("max |A Phi_exact / V - Lap Phi_exact|")
    axs[0].legend(fontsize=7)
    axs[0].set_title("truncation: bounded (B, M), growing as 1/h (S1)", fontsize=10)
    for p, y0, lab in ((1, 2.0, "h"), (2, 0.2, "h^2")):
        axs[1].loglog(ns, y0 * (ns[0] / ns) ** p, "-", color="0.6", lw=0.8)
        axs[1].annotate(lab, (20.0, y0 * (8.0 / 20.0) ** p), xytext=(0, -12), textcoords="offset points", fontsize=8,
                        color="0.4")
    axs[1].set_xlabel("N (coarse cells per side)")
    axs[1].set_ylabel("max error")
    axs[1].legend(fontsize=6, ncol=2)
    axs[1].set_title("the seam gradient converges for B and M, not for S1", fontsize=10)
    for dim, color in (("3-D", "tab:blue"), ("2-D", "tab:green")):
        lev, _, _, tau, e_phi, e_g = (np.array(v) for v in zip(*sorted(amr[dim]), strict=True))
        nf = 8.0 * 2.0 ** (lev + 1)  # fine cells along x: 8 per block side, the box refined once

        def order(e: np.ndarray, nf: np.ndarray = nf) -> float:
            return float(np.log(e[0] / e[-1]) / np.log(nf[-1] / nf[0]))

        axs[2].loglog(nf, e_phi, "o-", color=color,
                      label=f"{dim}: potential, order {order(e_phi):.2f} (truncation growth {tau[-1] / tau[0]:.3f})")
        axs[2].loglog(nf, e_g, "s--", color=color, label=f"{dim}: 2:1 face gravity, order {order(e_g):.2f}")
    axs[2].set_xticks([64, 128], ["64", "128"])
    axs[2].minorticks_off()
    axs[2].set_xlim(56, 146)
    axs[2].set_ylim(2e-4, 2e-1)
    axs[2].set_xlabel("fine cells along x")
    axs[2].set_ylabel("max error")
    axs[2].legend(fontsize=7, loc="upper right")
    axs[2].set_title("the multigrid solver on M (box refined once)", fontsize=10)
    fig.suptitle("SV-9: a symmetric, consistent composite Laplacian across 2:1 faces and its multigrid solver "
                 "(CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-amr.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_amr_jeans(out: Path, runs: Path) -> None:
    """SV-10 (issue #81, M8-P5b): self-gravity on a static AMR grid, the Jeans cases of selfgravity/make_jeans.py with
    `--box` (the central box refined once: 2:1 seams the wave crosses; the multigrid solver, the seam gravity work in the
    reflux; CPU, 2 ranks). Left: the L1 error of the stable wave after one period against the base N (the box has twice
    the resolution), and the same wave on the uniform grid of the base N. Middle: the collapse (amplitude 1e-2,
    base 64^2), the drift of E + W at CFL 0.4, 0.2, 0.1 with the end-of-step reflux (solid) and the per-stage one
    (dashed, `reflux_cadence = stage`, P5c), with E and W alone. Right: the collapse density at t = 0.4 on
    the composite grid, the block outlines drawn."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415

    fig, axs = plt.subplots(1, 3, figsize=(16.0, 4.9), constrained_layout=True, width_ratios=(1.0, 1.0, 1.1))
    ax = axs[0]
    for label, tag, color, mk in (("AMR, box refined once", "sv10-2d", "tab:blue", "o"),
                                  ("uniform", "sv10-uniform", "0.5", "s")):
        works = [runs / f"{tag}-{c}" for c in (8, 16, 32)]
        if not all(w.is_dir() for w in works):
            continue
        n = np.array([32.0, 64.0, 128.0])
        err = np.array([jeans_oracle.wave_error(w)[1] for w in works])
        ax.loglog(n, err, f"{mk}-", color=color, mfc="none", ms=8,
                  label=f"{label}: order {math.log2(err[-2] / err[-1]):.2f}")
    ns = np.array([32.0, 128.0])
    ax.loglog(ns, 8e-4 * (32.0 / ns) ** 2, "-", color="0.6", lw=0.8)
    ax.annotate("dx^2", (64.0, 8e-4 * 0.25), xytext=(0, -12), textcoords="offset points", fontsize=8, color="0.4")
    ax.set_xlabel("base N (cells per side)")
    ax.set_ylabel("L1 rho / (rho_0 A) after one period")
    ax.legend(fontsize=8)
    ax.set_title("the stable wave across the 2:1 seams: second order", fontsize=10)
    ax = axs[1]
    for cfl, color in (("0.4", "tab:blue"), ("0.2", "tab:orange"), ("0.1", "tab:green")):
        h = jeans_oracle.history(runs / f"sv10-collapse-{cfl}")
        ew0 = h["E_plus_W"][0]
        ax.semilogy(h["time"][1:], np.maximum(np.abs(h["E_plus_W"][1:] - ew0) / abs(ew0), FLOOR), "-", color=color,
                    label=f"E + W, CFL {cfl}, reflux each step")
        hs = jeans_oracle.history(runs / f"sv10-stage-{cfl}")
        ax.semilogy(hs["time"][1:], np.maximum(np.abs(hs["E_plus_W"][1:] - hs["E_plus_W"][0]) / abs(ew0), FLOOR),
                    "--", color=color, label=f"E + W, CFL {cfl}, reflux each stage")
        if cfl == "0.1":
            ax.semilogy(h["time"][1:], np.abs(h["int_rE"][1:] - h["int_rE"][0]) / abs(ew0), "--", color="0.4",
                        label="E alone")
            ax.semilogy(h["time"][1:], np.abs(h["W"][1:] - h["W"][0]) / abs(ew0), ":", color="0.4", label="W alone (on E alone)")
    ax.set_xlabel("t")
    ax.set_ylabel("relative change")
    ax.legend(fontsize=6)
    ax.set_title("E + W: first order in dt (step reflux), third (stage reflux)", fontsize=10)
    ax = axs[2]
    _draw_composite_density(fig, ax, runs / "sv10-collapse-0.1")
    fig.suptitle("SV-10: Jeans modes under self-gravity on a static AMR grid, the multigrid solver (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-amr-jeans.png", dpi=DPI)
    plt.close(fig)


def _draw_composite_density(fig: plt.Figure, ax: plt.Axes, work: Path, step: int | None = None) -> None:
    """The density of one checkpoint of a 2-D FLUME run (the last, or `step`) on its composite grid, block by block,
    the block outlines drawn."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415

    ini = jeans_oracle.read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    cps = jeans_oracle.checkpoints(work)
    tc, paths = cps[max(cps) if step is None else step]
    rmin, rmax = math.inf, -math.inf
    blocks = []
    for path in paths:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                r = h5[f"{blk}-r"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, 0]
                rmin, rmax = min(rmin, float(r.min())), max(rmax, float(r.max()))
                x0, y0 = origin[0] + ngc * dxyz[0], origin[1] + ngc * dxyz[1]
                blocks.append((x0, y0, r.shape[0] * dxyz[0], r.shape[1] * dxyz[1], r))
    for x0, y0, wx, wy, r in blocks:
        im = ax.imshow(r.T, origin="lower", extent=(x0, x0 + wx, y0, y0 + wy), cmap="magma", vmin=rmin, vmax=rmax)
        ax.add_patch(Rectangle((x0, y0), wx, wy, fill=False, lw=0.4, ec="w", alpha=0.6))
    lo = float(ini["grid"]["emin_x"])
    ax.set_xlim(lo, lo + 1.0)
    ax.set_ylim(lo, lo + 1.0)
    fig.colorbar(im, ax=ax, shrink=0.85)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"the density at t = {tc:.2f}, {len(blocks)} blocks (max {rmax:.3f})", fontsize=10)


def fig_selfgravity_amr_isolated(out: Path, runs: Path) -> None:
    """SV-13 (issue #87, P2): isolated boundaries on AMR, the symmetric James step of adam_poisson_mg_object, from the
    log of test_poisson_mg --isolated (CPU, 2 ranks): SV-6's stars (radius 1/2) on the cube [-1,1]^3, uniform (16^3,
    32^3) or with the star's box refined once over base levels 2 and 3. Left: the max error of the potential on the
    leaves (solid) and of the James ghost values against the exact potential (dashed) against the finest cells per side.
    Right: the two exactness measures against their bounds, the multigrid against the FFT solver's isolated mode
    (uniform) and the symmetry defect of the solution map in situ (the first grid of each leg)."""
    pattern = re.compile(r"(SV-13[ab]), isolated, (uniform|AMR), (polytrope n=1|uniform sphere),\s+base level (\d+), "
                         r"leaves \d+: iterations (\d+) \+ (\d+)(.*?), \|Phi - Phi_exact\| (\S+), "
                         r"boundary ghosts (\S+),")
    legs: dict[str, list[tuple[int, float, float]]] = {}
    exact: list[tuple[str, float, str]] = []
    for line in (runs / "sv13-poisson" / "log.txt").read_text(errors="replace").splitlines():
        m = pattern.search(line)
        if not m:
            continue
        key = f"{m[2]}, {m[3]}"
        legs.setdefault(key, []).append((int(m[4]), float(m[8]), float(m[9])))
        fft = re.search(r"\|Phi_mg - Phi_fft\| / max\|Phi\| ([^,\s]+)", m[7])
        sym = re.search(r"symmetry of V K ([^,\s]+)", m[7])
        tick = f"{m[2]}\n{'polytrope' if m[3].startswith('poly') else 'sphere'}\nlevel {m[4]}"
        if fft:
            exact.append((tick, float(fft[1]), "multigrid - FFT"))
        if sym:
            exact.append((tick, float(sym[1]), "symmetry of V K"))
    if not legs:
        raise FileNotFoundError("no SV-13 line in sv13-poisson/log.txt")
    fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.9), constrained_layout=True, width_ratios=(1.1, 1.0))
    colors = {"uniform, polytrope n=1": "tab:gray", "AMR, polytrope n=1": "tab:blue",
              "AMR, uniform sphere": "tab:orange"}
    for key, rows in legs.items():
        lev, e_phi, e_gh = (np.array(v) for v in zip(*sorted(rows), strict=True))
        nf = 8.0 * 2.0 ** (lev + (1 if key.startswith("AMR") else 0))  # finest cells per side
        order = float(np.log(e_phi[0] / e_phi[-1]) / np.log(nf[-1] / nf[0]))
        axs[0].loglog(nf, e_phi, "o-", color=colors[key], label=f"{key}: potential, order {order:.2f}")
        axs[0].loglog(nf, e_gh, "s--", color=colors[key], alpha=0.7, label=f"{key}: boundary ghosts")
    ns = np.array([16.0, 128.0])
    axs[0].loglog(ns, 2e-2 * (ns[0] / ns) ** 2, "-", color="0.6", lw=0.8)
    axs[0].annotate("h^2", (40.0, 2e-2 * (16.0 / 40.0) ** 2), xytext=(0, -12), textcoords="offset points", fontsize=8,
                    color="0.4")
    axs[0].set_xticks([16, 32, 64, 128], ["16", "32", "64", "128"])
    axs[0].minorticks_off()
    axs[0].set_xlabel("finest cells per side")
    axs[0].set_ylabel("max |Phi - Phi_exact|")
    axs[0].set_ylim(1e-5, 1e-1)
    axs[0].legend(fontsize=7, loc="upper right")
    axs[0].set_title("second order on AMR; the James ghosts converge with the potential", fontsize=10)
    ax = axs[1]
    kinds = {"multigrid - FFT": ("tab:green", 1e-9), "symmetry of V K": ("tab:purple", 1e-9)}
    for i, (lab, val, kind) in enumerate(exact):
        ax.bar(i, val, color=kinds[kind][0], label=kind if kind not in [e[2] for e in exact[:i]] else None)
    ax.axhline(1e-9, color="tab:red", lw=1.0, ls="--", label="bound, 1e-9")
    ax.axhline(1e-12, color="0.5", lw=0.8, ls=":", label="FCG tolerance, 1e-12")
    ax.set_yscale("log")
    ax.set_ylim(1e-16, 1e-7)
    ax.set_xticks(range(len(exact)), [e[0] for e in exact], fontsize=8)
    ax.set_ylabel("relative defect")
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("the free-space solution of the FFT solver; a symmetric map", fontsize=10)
    fig.suptitle("SV-13: isolated self-gravity on AMR, the symmetric James step (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-amr-isolated.png", dpi=DPI)
    plt.close(fig)


def fig_selfgravity_amr_regrid(out: Path, runs: Path) -> None:
    """SV-11 (issue #81, M8-P5c): the collapse of make_jeans.py (base 32^2) under runtime AMR, regridded every 5 steps by
    a Loehner estimator on the density, with the per-stage reflux (CPU, 2 ranks). Left: the change of E + W (relative)
    against t, the regrid steps that changed the grid marked, and the number of blocks of every checkpoint. Right: the
    density at the end on the composite grid."""
    sys.path.insert(0, str(HERE / "selfgravity"))
    import jeans_oracle  # noqa: PLC0415

    work = runs / "sv11-regrid"
    h = jeans_oracle.history(work)
    ew = h["E_plus_W"]
    d = np.diff(ew) / abs(ew[0])
    big = np.abs(d) > 100.0 * np.median(np.abs(d))
    fig, axs = plt.subplots(1, 2, figsize=(13.0, 4.9), constrained_layout=True, width_ratios=(1.2, 1.0))
    ax = axs[0]
    ax.semilogy(h["time"][1:], np.maximum(np.abs(ew[1:] - ew[0]) / abs(ew[0]), FLOOR), "-", color="tab:blue",
                label="|E + W - (E + W)(0)|")
    for i in np.nonzero(big)[0]:
        ax.axvline(h["time"][i + 1], color="tab:red", lw=0.8, ls=":")
        ax.annotate(f"regrid, {d[i]:+.1e}", (h["time"][i + 1], 1e-6), rotation=90, fontsize=7, color="tab:red",
                    xytext=(3, 0), textcoords="offset points")
    ax.set_xlabel("t")
    ax.set_ylabel("relative change of E + W")
    ax.legend(fontsize=8, loc="upper left")
    cps = jeans_oracle.checkpoints(work)
    times = [cps[k][0] for k in sorted(cps)]
    nblk = [jeans_oracle.checkpoint_blocks(cps[k][1]) for k in sorted(cps)]
    ax2 = ax.twinx()
    ax2.step(times, nblk, where="post", color="tab:green", lw=1.0, ls="--", alpha=0.7, label="blocks")
    ax2.set_ylabel("blocks (dashed)", color="tab:green")
    ax2.set_ylim(0, 1.25 * max(nblk))
    ax2.legend(fontsize=8, loc="center left")
    ax.set_title("E + W changes only at the regrids that change the grid", fontsize=10)
    _draw_composite_density(fig, axs[1], work)
    fig.suptitle("SV-11: self-gravity under runtime AMR, the collapse regridded on the density (CPU, 2 ranks)")
    fig.savefig(out / "selfgravity-amr-regrid.png", dpi=DPI)
    plt.close(fig)


def _rank_blocks(paths: list[Path]) -> dict[int, int]:
    """Blocks of each rank in one checkpoint (rank -> count), from the per-rank files."""
    out = {}
    for path in paths:
        with h5py.File(path, "r") as h5:
            out[int(path.stem.rsplit("-proc", 1)[1])] = len({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")})
    return out

def fig_viscous_ladders(out: Path, runs: Path) -> None:  # noqa: ARG001
    """VV-1, VV-2 and VV-5 (issue #85, R1): the absolute RMS error of every wave leg against the resolution, CPU (filled,
    solid) and FNL (open, dotted), and the Ohmic decay: the measured amplitude of the mode over the exact one."""
    sys.path.insert(0, str(HERE / "viscous"))
    import waves  # noqa: PLC0415

    vdir = HERE / "viscous"
    fnl = "work-adam_flume_fnl-np2"
    navier = (("VV-1 shear x, order 4", "o", "shear-x-o4", (32, 64, 128)),
              ("VV-1 shear x, order 2", "s", "shear-x-o2", (32, 64, 128)),
              ("VV-1 shear 45°, order 4", "d", "shear-xy-o4", (24, 48, 96)),
              ("VV-2 acoustic, order 4", "^", "acoustic-x-o4", (32, 64, 128)))
    ohmic = (("VV-5 Ohmic x, order 4", "o", "ohmic-x-o4", (32, 64, 128)),
             ("VV-5 Ohmic x, order 2", "s", "ohmic-x-o2", (32, 64, 128)),
             ("VV-5 Ohmic 45°, order 4", "d", "ohmic-xy-o4", (24, 48, 96)),
             ("VV-5 Alfvén, $\\mu$ and $\\eta$, order 4", "^", "alfven-x-o4", (32, 64, 128)),
             ("VV-5 acoustic on MHD, B = 0, order 4", "v", "acoustic-mhd-x-o4", (32, 64, 128)))
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    for ax, legs, title in ((axs[0], navier, "VV-1, VV-2: Navier-Stokes waves"),
                            (axs[1], ohmic, "VV-5: Ohmic and visco-resistive waves")):
        for i, (lab, mk, name, ns) in enumerate(legs):
            for tag, solid in ((TAG, True), (fnl, False)):
                res = [waves.run_error(vdir / f"{tag}-{name}-n{n}") for n in ns]
                ax.loglog([r[0] for r in res], [r[1] for r in res], mk + ("-" if solid else ":"), color=f"C{i}",
                          mfc=f"C{i}" if solid else "none", ms=5 if solid else 8,
                          label=f"{lab}, {'CPU' if solid else 'FNL'}")
            if name == legs[0][2]:
                e0, n0 = res[0][1], res[0][0]
        n = np.array([24.0, 128.0])
        ax.loglog(n, e0 * 3 * (n0 / n) ** 4, "k:", lw=0.8, label="order 4")
        ax.loglog(n, e0 * 3 * (n0 / n) ** 2, "k--", lw=0.8, label="order 2")
        ax.set_xlabel("cells along the wave")
        ax.set_ylabel("RMS error")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=6.5, ncol=2)
        ax.set_title(title)
    for i, (lab, mk, name, ns) in enumerate(ohmic[:3]):
        for tag, solid in ((TAG, True), (fnl, False)):
            res = [waves.run_error(vdir / f"{tag}-{name}-n{n}") for n in ns]
            axs[2].loglog([r[0] for r in res], [abs(r[2] / r[3] - 1.0) for r in res], mk + ("-" if solid else ":"),
                          color=f"C{i}", mfc=f"C{i}" if solid else "none", ms=5 if solid else 8,
                          label=f"{lab}, {'CPU' if solid else 'FNL'}")
    axs[2].set_xlabel("cells along the wave")
    axs[2].set_ylabel(r"$|A_\mathrm{measured} / A_\mathrm{exact} - 1|$")
    axs[2].grid(alpha=0.3, which="both")
    axs[2].legend(fontsize=6.5)
    axs[2].set_title(r"VV-5: Ohmic decay against the exact rate $e^{-\eta k^2 t}$")
    fig.savefig(out / "viscous-ladders.png", dpi=DPI)
    plt.close(fig)


def proto2d() -> object:
    """Import amr-proto/proto2d.py, the NumPy prototype of #68 Phase B."""
    sys.path.insert(0, str(HERE / "amr-proto"))
    import proto2d as module  # noqa: PLC0415

    return module


def fig_seam_proto_b1(out: Path, runs: Path) -> None:  # noqa: ARG001
    """B1 (#68, issue #85 R1): the 2-D prototype, linear advection across a 2:1 box, N = 32 / 64 / 128. L1 and Linf
    against N per arm (restriction, prolongation, reflux), and where the finest Linf error lives (corner, edge,
    interior of the box seam)."""
    pr = proto2d()
    cells = (32, 64, 128)
    arms = [("uniform", None), ("mean-cubic-on", pr.ARMS["mean-cubic-on"]),
            ("mean-cubic-off", pr.ARMS["mean-cubic-off"]), ("lag6-cubic-off", pr.ARMS["lag6-cubic-off"]),
            ("lag6-quintic-on", pr.ARMS["lag6-quintic-on"]), ("lag6-quintic-off", pr.ARMS["lag6-quintic-off"])]
    res = {name: [pr.run(n, arm, (1.0, 1.0)) for n in cells] for name, arm in arms}
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    for ax, norm, label in ((axs[0], "l1", "$L_1$"), (axs[1], "linf", r"$L_\infty$")):
        for i, (name, _) in enumerate(arms):
            e = [r[norm] for r in res[name]]
            ax.loglog(cells, e, "o-" if name != "uniform" else "k*-", color=None if name == "uniform" else f"C{i}",
                      mfc="none", label=f"{name} ({math.log2(e[-2] / e[-1]):.2f})")
        n = np.array(cells, dtype=float)
        top = res["mean-cubic-on"][0][norm]
        ax.loglog(n, top * (32 / n) ** 2, "k--", lw=0.8, label="order 2")
        ax.loglog(n, top * 0.1 * (32 / n) ** 5, "k:", lw=0.8, label="order 5")
        ax.set_xlabel("coarse cells per side N")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7, title="arm (finest-pair order)", title_fontsize=7)
        ax.set_title(f"{label} error at t = 1")
    names = [name for name, _ in arms[1:]]
    x = np.arange(len(names))
    for j, (region, color) in enumerate((("corner", "C3"), ("edge", "C1"), ("interior", "C2"))):
        axs[2].bar(x + (j - 1) * 0.27, [res[name][-1][f"linf_{region}"] for name in names], 0.27, color=color,
                   label=region)
    axs[2].set_yscale("log")
    axs[2].set_xticks(x, names, rotation=25, fontsize=7)
    axs[2].set_ylabel(r"$L_\infty$ at N = 128")
    axs[2].legend(fontsize=7)
    axs[2].grid(alpha=0.3, axis="y")
    axs[2].set_title("where the error lives (within 3 coarse cells of the seam)")
    fig.suptitle(r"B1: $u_t + u_x + u_y = 0$, the box $[0.25, 0.75]^2$ refined 2:1, 5th-order linear upwind, RK4")
    fig.savefig(out / "seam-proto-b1.png", dpi=DPI)
    plt.close(fig)


def fig_seam_proto_b2(out: Path, runs: Path) -> None:  # noqa: ARG001
    """B2 (#68, issue #85 R1): WENO5-JS on the 2-D prototype, the reflux policies always, troubled, never. The smooth
    order (N = 32 / 64 / 128), and on the Burgers pulses at N = 64 the mean shock position error and the mass defect."""
    pr = proto2d()
    cells = (32, 64, 128)
    arms = [("uniform", None), *pr.B2_SMOOTH.items()]
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    for i, (name, arm) in enumerate(arms):
        e = [pr.run(n, arm, (1.0, 1.0), scheme="weno")["l1"] for n in cells]
        axs[0].loglog(cells, e, "o-" if arm else "k*-", color=f"C{i}" if arm else None, mfc="none",
                      label=f"{name} ({math.log2(e[-2] / e[-1]):.2f})")
    n = np.array(cells, dtype=float)
    axs[0].loglog(n, 2e-3 * (32 / n) ** 5, "k:", lw=0.8, label="order 5")
    axs[0].set_xlabel("coarse cells per side N")
    axs[0].set_ylabel("$L_1$")
    axs[0].grid(alpha=0.3, which="both")
    axs[0].legend(fontsize=7, title="arm (finest-pair order)", title_fontsize=7)
    axs[0].set_title("smooth advection: the order of each policy")
    cases = list(pr.PULSES)
    pols = [("uniform", None), *pr.B2_SMOOTH.items()]
    shock = {(c, name): pr.run_shock(64, arm, c) for c in cases for name, arm in pols}
    x = np.arange(len(cases))
    w = 0.8 / len(pols)
    for j, (name, _) in enumerate(pols):
        axs[1].bar(x + (j - 2) * w, [shock[(c, name)]["shock_mean"] for c in cases], w, label=name,
                   color="k" if name == "uniform" else f"C{j}")
        axs[2].bar(x + (j - 2) * w, [max(abs(shock[(c, name)]["mass"]), 1e-17) for c in cases], w, label=name,
                   color="k" if name == "uniform" else f"C{j}")
    for ax in axs[1:]:
        ax.set_xticks(x, ["c2f (x = 0.25)", "f2c (x = 0.75)", "weak (1e-3)"])
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=7)
    axs[1].axhline(0.0, color="k", lw=0.6)
    axs[1].set_ylabel("mean shock position error (local cells)")
    axs[1].set_title("Burgers shocks crossing the box face, N = 64")
    axs[2].set_yscale("log")
    axs[2].axhline(1e-3, color="k", ls="--", lw=0.8)
    axs[2].set_ylabel("|mass defect| / pulse mass (floor 1e-17)")
    axs[2].set_title("mass defect (dashed: the 1e-3 criterion)")
    fig.suptitle("B2: the reflux policies on WENO5-JS (Lax-Friedrichs splitting), lag6 restriction, quintic fill")
    fig.savefig(out / "seam-proto-b2.png", dpi=DPI)
    plt.close(fig)


def fig_seam_proto_b3(out: Path, runs: Path) -> None:  # noqa: ARG001
    """B3 (#68, issue #85 R1): the footprint of the 6-point restriction against the fine block size, reflux off: Linf
    against N for shift-inward (lag6) at 6 and 8 fine cells per block and for the whole box, and the centred two-pass
    lag6c, with the cubic and the quintic prolongation."""
    pr = proto2d()
    cells = (24, 48, 96)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True, sharey=True)
    for ax, prol in zip(axs, ("cubic", "quintic"), strict=True):
        for i, (kind, block) in enumerate((("lag6", 6), ("lag6", 8), ("lag6", 0), ("lag6c", 0))):
            label = "lag6c (centred, fine ghosts)" if kind == "lag6c" else \
                f"lag6, {block} fine cells per block" if block else "lag6, one block (the whole box)"
            e = [pr.run(n, (kind, prol, False), (1.0, 1.0), block=block)["linf"] for n in cells]
            ax.loglog(cells, e, "o-", color=f"C{i}", mfc="none", label=f"{label} ({math.log2(e[-2] / e[-1]):.2f})")
        n = np.array(cells, dtype=float)
        ax.loglog(n, e[0] * (24 / n) ** (4 if prol == "cubic" else 5), "k:", lw=0.8,
                  label=f"order {4 if prol == 'cubic' else 5}")
        ax.set_xlabel("coarse cells per side N")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7, title="restriction (finest-pair order)", title_fontsize=7)
        ax.set_title(f"{prol} prolongation")
    axs[0].set_ylabel(r"$L_\infty$ at t = 1")
    fig.suptitle("B3: the 6-point restriction footprint against the fine block size (reflux off; 4 fine cells per "
                 "block is infeasible for shift-inward)")
    fig.savefig(out / "seam-proto-b3.png", dpi=DPI)
    plt.close(fig)


# E0 (#68): the plateau density defect behind the shock after it crossed a static 2:1 face, FNL, N = 200, the table of
# docs/applications/flume/verification/seams-studies/capturing-error.md (the runs keep no checkpoint)
E0_DEFECT = {"compatible + mean": (2.32e-2, 1.82e-2, 1.97e-2, 1.83e-2),
             "compatible + lagrange6": (2.32e-2, 1.85e-2, 1.98e-2, 1.84e-2),
             "tricubic + mean": (2.46e-2, 3.10e-2, 1.41e-2, 1.68e-2),
             "tricubic + lagrange6": (2.46e-2, 3.17e-2, 1.41e-2, 1.68e-2)}


def fig_capturing_error(out: Path, runs: Path) -> None:  # noqa: ARG001
    """E0 (#68, issue #85 R1): the plateau density defect against fill x restriction x reflux, Sod and RJ2a."""
    fig, ax = plt.subplots(figsize=(9, 4.4), constrained_layout=True)
    cols = ("Sod, reflux on", "Sod, reflux off", "RJ2a, reflux on", "RJ2a, reflux off")
    x = np.arange(len(cols))
    w = 0.8 / len(E0_DEFECT)
    for j, (arm, vals) in enumerate(E0_DEFECT.items()):
        ax.bar(x + (j - 1.5) * w, vals, w, label=arm, color=f"C{j}")
    ax.set_xticks(x, cols)
    ax.set_ylabel("plateau density defect")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8, title="coarse-to-fine fill + restriction", title_fontsize=8)
    ax.set_title("E0: the defect behind a shock that crossed a 2:1 face (FNL); no operator removes it")
    fig.savefig(out / "capturing-error.png", dpi=DPI)
    plt.close(fig)

def fig_multirealm(out: Path, runs: Path) -> None:
    """Multi-realm legs 1 and 4 (issue #85, R2a): Sod split at the diaphragm into two realms glued by a mirror seam,
    and the refined Sod split at its 2:1 face into a coarse and a fine realm (coupling = refined), each union over its
    single-realm twin; the maximum difference of every cell is in the legend (the check asserts bitwise)."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.4), constrained_layout=True)
    legs = (("sod-2realm", "sod-x", "mirror seam, 1:1, x = 0.5", "realm 1", "realm 2"),
            ("sod-amr-refined", "sod-amr", "refined seam, 2:1, x = 0.75", "realm 1 (coarse)", "realm 2 (fine)"))
    for ax, (multi, single, title, lab1, lab2) in zip(axs, legs, strict=True):
        w, ws = runs / multi, runs / single
        ngc = int(read_ini(next(ws.glob("*.ini")))["grid"]["ngc"])
        _, xs, q = profile(ws, ngc, EULER)
        ax.plot(xs, q[0], "k-", lw=1.0, label=f"single realm ({single})")
        worst = 0.0
        ref = dict(zip(np.round(xs, 12), q[0], strict=True))
        for prefix, lab, mk in ((f"{multi}-r1", lab1, "o"), (f"{multi}-r2", lab2, "s")):
            _, xr, qr = profile(w, ngc, EULER, prefix)
            ax.plot(xr, qr[0], mk, ms=3, mfc="none", label=lab)
            worst = max(worst, max(abs(ref[x] - r) for x, r in zip(np.round(xr, 12), qr[0], strict=True)))
        seam = 0.5 if multi == "sod-2realm" else 0.75
        ax.axvline(seam, color="C3", ls="--", lw=0.8, label="seam")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$\rho$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title=f"max |union - single| = {worst:.1e}", title_fontsize=8)
        ax.set_title(f"Sod at t = 0.2, {title}")
    fig.suptitle("Multi-realm runs: the union of the realms against the single-realm run (CPU, 2 ranks)")
    fig.savefig(out / "multirealm.png", dpi=DPI)
    plt.close(fig)


def fig_rotated_shock_tube(out: Path, runs: Path) -> None:
    """MV-8 (issue #85, R2a): the rotated shock tube at N = 256. Left, the deviation of the normal field from its exact
    uniform value 5/sqrt(4 pi); right, every cell's density and tangential field against s = n.x, over the 1-D
    reference run (the same problem along x, 4 times finer along the normal)."""
    sys.path.insert(0, str(HERE / "mhd" / "rotated-shock-tube"))
    import rotated_shock_tube_oracle as rst  # noqa: PLC0415

    b0 = 5.0 / math.sqrt(4.0 * math.pi)
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.4), constrained_layout=True)
    for row, angle in enumerate((63, 45)):
        w, wr = runs / f"rst-{angle}-256", runs / f"rst-{angle}-ref"
        nx, ny, _, gamma = rst.setup(w)
        nn = math.hypot(nx, ny)
        nx, ny = nx / nn, ny / nn  # setup returns the normal unnormalised, (1, 2) or (1, 1)
        ngc = int(read_ini(next(w.glob("*.ini")))["grid"]["ngc"])
        bl = first_plane(blocks(w, ngc, MHD))
        im = field_map(axs[row, 0], bl, lambda b, nx=nx, ny=ny: (nx * b["f"]["bx"][:, :, 0] +
                                                                  ny * b["f"]["by"][:, :, 0]) / b0 - 1.0,
                       cmap="RdBu_r")
        fig.colorbar(im, ax=axs[row, 0], label=r"$B_\parallel / B_0 - 1$", shrink=0.85)
        axs[row, 0].set_title(f"normal at {math.degrees(math.atan2(ny, nx)):.1f}°: $B_\\parallel$ deviation")
        s2, p2 = rst.primitives(w, ngc)
        sr, pr = rst.primitives(wr, ngc)
        order = np.argsort(sr)
        for col, k, lab in ((1, 0, r"$\rho$"), (2, 5, r"$B_\perp$")):
            axs[row, col].plot(s2, p2[:, k], ".", ms=1, alpha=0.3, color="C0", label="2-D cells, N = 256")
            axs[row, col].plot(sr[order], pr[order, k], "k-", lw=0.9, label="1-D reference")
            axs[row, col].set_xlabel(r"$s = \mathbf{n}\cdot\mathbf{x}$")
            axs[row, col].set_title(lab)
            axs[row, col].grid(alpha=0.3)
            axs[row, col].legend(fontsize=8, markerscale=6)
    fig.suptitle(r"Rotated shock tube (Ryu-Jones 1a), $256^2$, GLM: every 2-D cell against the 1-D reference")
    fig.savefig(out / "rotated-shock-tube.png", dpi=DPI)
    plt.close(fig)


def fig_weno_riemann(out: Path, runs: Path) -> None:
    """RV-5 (issue #85, R2a): Ryu-Jones 2a at 256 cells on the splitting scheme and on weno-riemann HLLD with
    characteristic and with primitive interpolation, over the exact solution."""
    sys.path.insert(0, str(HERE / "mhd" / "rj2a"))
    import rj2a_oracle  # noqa: PLC0415

    def rj2a_exact(xs: np.ndarray) -> dict[str, np.ndarray]:
        q = rj2a_oracle.exact(xs)
        kin = 0.5 * (q[1] ** 2 + q[2] ** 2 + q[3] ** 2) / q[0]
        mag = 0.5 * (q[5] ** 2 + q[6] ** 2 + q[7] ** 2)
        return {"r": q[0], "u": q[1] / q[0], "p": (rj2a_oracle.GAMMA - 1.0) * (q[4] - kin - mag), "bt1": q[6],
                "bt2": q[7]}

    mhd_profile_figure(out, "weno-riemann", runs / "rj2a", "Ryu-Jones 2a, 256 cells: the splitting scheme and HLLD",
                       rj2a_exact, label="weno (split)",
                       extra=(("HLLD, characteristic", runs / "rj2a-hlld", "s"),
                              ("HLLD, primitive", runs / "rj2a-hlld-prim", "^")))


def fig_quadtree(out: Path, runs: Path) -> None:
    """MV-15 (issue #85, R2a): Orszag-Tang 32^2 with [0.25, 0.75]^2 refined 2:1, t = 0.2: the density on the quadtree
    (nk = 1) with its blocks, and the per-cell difference against the octree with a null z axis (first z layer),
    density map and the maximum per variable."""
    def cells(work: Path) -> dict[tuple[float, float], dict[str, float]]:
        ngc = int(read_ini(next(work.glob("*.ini")))["grid"]["ngc"])
        out_: dict[tuple[float, float], dict[str, float]] = {}
        for b in first_plane(blocks(work, ngc, (*MHD, "psi"))):
            nx, ny = b["f"]["r"].shape[:2]
            for i in range(nx):
                for j in range(ny):
                    key = (round(b["lo"][0] + (i + 0.5) * b["d"][0], 10), round(b["lo"][1] + (j + 0.5) * b["d"][1], 10))
                    out_[key] = {v: float(a[i, j, 0]) for v, a in b["f"].items()}
        return out_

    wq, wo = runs / "ot-quad-nk1", runs / "ot-octree"
    cq, co = cells(wq), cells(wo)
    ngc = int(read_ini(next(wq.glob("*.ini")))["grid"]["ngc"])
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    bl = first_plane(blocks(wq, ngc, (*MHD, "psi")))
    gamma = gamma_of(read_ini(next(wq.glob("*.ini"))))
    im = field_map(axs[0], bl, lambda b: plane(b, gamma, "r"), cmap="magma", outline=True)
    fig.colorbar(im, ax=axs[0], label=r"$\rho$", shrink=0.85)
    axs[0].set_title("quadtree (nk = 1): density and blocks")
    keys = sorted(cq)
    xs, ys = np.array([k[0] for k in keys]), np.array([k[1] for k in keys])
    d = np.array([abs(cq[k]["r"] - co[k]["r"]) for k in keys])
    sc = axs[1].scatter(xs, ys, c=np.log10(np.maximum(d, 1e-17)), s=4, marker="s", cmap="viridis")
    fig.colorbar(sc, ax=axs[1], label=r"$\log_{10}|\rho_\mathrm{quad} - \rho_\mathrm{oct}|$ (floor -17)", shrink=0.85)
    axs[1].set_aspect("equal")
    axs[1].set_title("quadtree against octree, density")
    names = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz", "psi")
    worst = [max(abs(cq[k][v] - co[k][v]) / max(max(abs(co[c][v]) for c in keys), 1e-300) for k in keys)
             for v in names]
    axs[2].bar(range(len(names)), np.maximum(worst, 1e-17), color="C0")
    axs[2].set_yscale("log")
    axs[2].set_xticks(range(len(names)), names)
    axs[2].set_ylabel("max |quadtree - octree| / max |octree|")
    axs[2].grid(alpha=0.3, axis="y")
    axs[2].set_title("per variable: round-off, not structure")
    fig.suptitle(r"Quadtree AMR (MV-15): Orszag-Tang, $32^2$ + $[0.25, 0.75]^2$ refined 2:1, t = 0.2, CPU")
    fig.savefig(out / "quadtree.png", dpi=DPI)
    plt.close(fig)


def fig_regrid_runtime(out: Path, runs: Path) -> None:
    """RG rg3 (issue #85, R2a): the isentropic vortex crossing a periodic quadtree, a regrid every 5 steps. The density
    and the blocks at three steps, and the relative drift of the volume integrals per step, conservative against
    linear prolongation (the regrid steps marked)."""
    wc, wl = runs / "regrid-vortex-conservative", runs / "regrid-vortex-linear"
    ini = read_ini(wc / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    its = sorted({int(p.name.split("-")[-2]) for p in wc.glob("*-proc*.h5")})
    picks = (its[0], its[len(its) // 2], its[-1])
    fig, axs = plt.subplots(1, 4, figsize=(20, 4.6), constrained_layout=True)
    for ax, it in zip(axs[:3], picks, strict=True):
        bl = first_plane(blocks(wc, ngc, EULER, it=it))
        im = field_map(ax, bl, lambda b: b["f"]["r"][:, :, 0], cmap="viridis", outline=True)
        ax.set_title(f"step {it}: {len(bl)} blocks")
    fig.colorbar(im, ax=axs[2], label=r"$\rho$", shrink=0.85)
    freq = int(ini["amr"]["frequency"])
    for w, lab, ls in ((wc, "conservative", "-"), (wl, "linear", "--")):
        lines = next(w.glob("*-conservation_history.dat")).read_text().splitlines()
        names = [v.strip('"').removeprefix("int_") for v in lines[0].split("=", 1)[1].split()][2:]
        rows = np.array([[float(v) for v in line.split()] for line in lines[1:] if line.split()])
        for i, v in enumerate(names):
            if v in ("r", "rE"):
                drift = np.abs(rows[:, 2 + i] - rows[0, 2 + i]) / abs(rows[0, 2 + i])
                axs[3].semilogy(rows[:, 0], np.maximum(drift, 1e-17), ls, label=f"{lab}, {v}")
    for k in range(freq, int(rows[-1, 0]) + 1, freq):
        axs[3].axvline(k, color="0.85", lw=0.5, zorder=0)
    axs[3].set_xlabel("step (grey: regrid steps)")
    axs[3].set_ylabel("relative drift (floor 1e-17)")
    axs[3].grid(alpha=0.3)
    axs[3].legend(fontsize=8)
    axs[3].set_title("volume integrals through the regrids")
    fig.suptitle("Runtime regridding (rg3): the vortex on a periodic quadtree, levels 1-3, a regrid every 5 steps")
    fig.savefig(out / "regrid-runtime.png", dpi=DPI)
    plt.close(fig)


def fig_seam_conservation(out: Path, runs: Path) -> None:  # noqa: ARG001
    """VV-7 (issue #85, R2a): the drift of the volume integrals per step with the dissipative fluxes across 2:1 seams
    (reflux on, off, and the ideal twin without dissipation) and across forest seams (mirror 1:1, refined 2:1), from
    the kept histories of viscous/check.sh --leg vv7 (CPU)."""
    vdir = HERE / "viscous"

    def history(path: Path) -> tuple[list[str], np.ndarray]:
        lines = path.read_text().splitlines()
        names = [v.strip('"').removeprefix("int_") for v in lines[0].split("=", 1)[1].split()][2:]
        return names, np.array([[float(v) for v in line.split()] for line in lines[1:] if line.split()])

    def drift(work: Path, var: str) -> tuple[np.ndarray, np.ndarray]:
        """The drift of one integral per step; a forest sums the histories of its realms (one file per realm)."""
        hs = [history(f) for f in sorted(work.glob("*-conservation_history.dat"))]
        names, rows = hs[0][0], sum(h[1] for h in hs)
        i = names.index(var)
        rows[:, 0] /= len(hs)
        scale = np.abs(rows[:, 2 + i]).max() or 1.0
        return rows[:, 0], np.maximum(np.abs(rows[:, 2 + i] - rows[0, 2 + i]) / scale, 1e-17)

    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    for ax, model, var in ((axs[0], "euler", "rE"), (axs[1], "mhd", "rE")):
        for suffix, lab, ls in (("", "reflux on", "-"), ("-off", "reflux off", "--"), ("-ideal", "ideal twin", ":")):
            for v, c in (("r", "C0"), (var, "C3")):
                it, d = drift(vdir / f"{TAG}-conserve-{model}-quadtree{suffix}", v)
                ax.semilogy(it, d, ls, color=c, label=f"{lab}, {v}")
        ax.set_title(f"{'Euler, acoustic' if model == 'euler' else 'MHD-GLM, Alfvén'} wave across a 2:1 box")
    for kind, ls in (("mirror", "-"), ("refined", "--")):
        hf = [history(f) for f in sorted((vdir / f"{TAG}-forest-euler-{kind}-forest").glob("*-conservation_history.dat"))]
        names, rows = hf[0][0], sum(h[1] for h in hf)
        single = history(next((vdir / f"{TAG}-forest-euler-{kind}-single").glob("*-conservation_history.dat")))[1]
        n = min(len(rows), len(single))
        for v, c in (("r", "C0"), ("rE", "C3")):
            i = 2 + names.index(v)
            d = np.abs(rows[:n, i] - single[:n, i]) / np.abs(single[:n, i]).max()
            axs[2].semilogy(single[:n, 0], np.maximum(d, 1e-17), ls, color=c, label=f"{kind} forest, {v}")
    axs[2].set_title("Euler across forest seams: union against the single realm")
    axs[2].set_ylabel("|I_forest - I_single| / max |I_single| (floor 1e-17)")
    for ax in axs[:2]:
        ax.set_ylabel("|I(t) - I(0)| / max |I| (floor 1e-17)")
    for ax in axs:
        ax.set_xlabel("step")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("VV-7: conservation with the dissipative fluxes across seams (CPU, 2 ranks)")
    fig.savefig(out / "seam-conservation.png", dpi=DPI)
    plt.close(fig)

def oracle_module(sub: str, name: str) -> object:
    """Import a verification oracle module from its case directory."""
    sys.path.insert(0, str(HERE / sub))
    return __import__(name)


AO = HERE / "amr-order"
AO_TAG = "work-adam_flume_cpu-np2"
AO_BOX = [0.5, 0.5, 1.0, 1.0]


def ao_ladder(leg: str, uniform: bool = False) -> tuple[list[int], list[dict]]:
    """Return the base cells and the region errors (amr_order_oracle.region_errors) of an AO Euler ladder."""
    ao = oracle_module("amr-order", "amr_order_oracle")
    ns, errs = [], []
    for n in (64, 128, 256):
        work = AO / f"{AO_TAG}-{leg}{'-uniform' if uniform else ''}-n{n}"
        cells, e = ao.region_errors("euler", work, None if uniform else AO_BOX)
        ns.append(cells)
        errs.append(e)
    return ns, errs


def ao_axes(ax_l1, ax_linf, legs: tuple[tuple[str, str, str], ...]) -> None:  # noqa: ANN001
    """Draw the L1 (all cells) and Linf (seam band) ladders of AO legs, with the uniform control."""
    ns, eu = ao_ladder("euler", uniform=True)
    ax_l1.loglog(ns, [e["all"][0] for e in eu], "k*-", label="uniform (coarse)")
    ax_linf.loglog(ns, [e["all"][1] for e in eu], "k*-", label="uniform (coarse)")
    for leg, lab, fmt in legs:
        ns, er = ao_ladder(leg)
        for ax, region, k in ((ax_l1, "all", 0), (ax_linf, "band", 1)):
            v = [e[region][k] for e in er]
            ax.loglog(ns, v, fmt, mfc="none", label=f"{lab} ({math.log2(v[-2] / v[-1]):.2f})")
    n = np.array([64.0, 256.0])
    for ax, top in ((ax_l1, eu[0]["all"][0]), (ax_linf, eu[0]["all"][1])):
        ax.loglog(n, top * 3 * (64 / n) ** 2, "k--", lw=0.8, label="order 2")
        ax.loglog(n, top * (64 / n) ** 5, "k:", lw=0.8, label="order 5")
        ax.set_xlabel("base cells per side N")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7, title="(128 → 256 order)", title_fontsize=7)
    ax_l1.set_ylabel(r"$L_1(\rho)$, all cells")
    ax_linf.set_ylabel(r"$L_\infty(\rho)$, seam band")


def fig_order_across_seams(out: Path, runs: Path) -> None:  # noqa: ARG001
    """AO (#68, issue #85 R2b): the Euler vortex across the static 2:1 quadrant [0.5, 1]^2, CPU, from the kept work
    directories of amr-order/check.sh. L1 and Linf ladders of today's operators (mean + tricubic, reflux on), the
    6-point restriction without reflux and with the quintic fill, against the uniform coarse control; and the density
    error maps at N = 128 of today's operators and of lagrange6 + quintic."""
    ao = oracle_module("amr-order", "amr_order_oracle")
    fig, axs = plt.subplots(1, 4, figsize=(21, 4.8), constrained_layout=True)
    ao_axes(axs[0], axs[1], (("euler", "mean + tricubic, reflux", "o-"), ("euler-l6", "lagrange6, no reflux", "s-"),
                             ("euler-q6", "lagrange6 + quintic, no reflux", "d-")))
    axs[0].set_title("all cells")
    axs[1].set_title("the seam band")
    maps = [(leg, lab) for leg, lab in (("euler", "mean + tricubic, reflux"), ("euler-q6", "lagrange6 + quintic"))]
    errs = []
    for leg, _ in maps:
        work = AO / f"{AO_TAG}-{leg}-n128"
        xy, h, q = ao.last_cells(work, ("r",))
        errs.append((xy, h, np.abs(q[:, 0] - ao.exact_density(xy[:, 0], xy[:, 1], ao.read_ini(work)))))
    lo = np.log10(max(min(float(e[2].max()) for e in errs) * 1e-4, 1e-14))
    hi = np.log10(max(float(e[2].max()) for e in errs))
    for ax, (leg, lab), (xy, h, e) in zip(axs[2:], maps, errs, strict=True):
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=np.log10(np.maximum(e, 10**lo)), s=(h[:, 0] / h.max()) * 6, marker="s",
                        cmap="viridis", vmin=lo, vmax=hi)
        ax.add_patch(Rectangle((0.5, 0.5), 0.5, 0.5, fill=False, ec="w", lw=0.8))
        ax.set_aspect("equal")
        ax.set_title(f"N = 128, {lab}: max {e.max():.1e}")
    fig.colorbar(sc, ax=axs[3], label=r"$\log_{10}|\rho - \rho_\mathrm{exact}|$", shrink=0.85)
    fig.suptitle("AO: the isentropic vortex across the static 2:1 quadrant (white), CPU, 2 ranks")
    fig.savefig(out / "order-across-seams.png", dpi=DPI)
    plt.close(fig)


def sod_columns(work: Path) -> tuple[np.ndarray, np.ndarray, object]:
    """Return the cell centres, the density and the exact Problem of an AS Sod run."""
    ash = oracle_module("amr-shock", "amr_shock_oracle")
    ini = ash.read_ini(work)
    prob = ash.Problem(ini)
    xs, _, q, _ = ash.load(work, ash.EULER, int(ini["grid"]["ngc"]))
    return xs, q[0], prob


def sod_zoom(ax, runs: Path, cases: tuple[tuple[str, str, str], ...], title: str, plateau: bool = True) -> None:  # noqa: ANN001
    """Draw the density behind the Sod shock (the plateau) of several AS runs, the exact value and the face x = 0.75."""
    for case, lab, mk in cases:
        xs, rho, prob = sod_columns(runs / case)
        ax.plot(xs, rho, mk, ms=3, mfc="none", label=lab)
    xe = np.linspace(prob.x_left - 0.05, prob.x_shock + 0.05, 2001)
    ax.plot(xe, prob.density(xe), "k-", lw=0.9, label="exact")
    ax.axvline(0.75, color="C3", ls="--", lw=0.8, label="2:1 face")
    if plateau:
        ax.set_xlim(prob.x_left - 0.03, prob.x_shock + 0.03)
        ax.set_ylim(prob.q2[0] - 0.03, prob.q2[0] + 0.03)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\rho$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    ax.set_title(title)


def fig_shock_crossing(out: Path, runs: Path) -> None:
    """AS (#68, issue #85 R2b): the Sod shock crossing the static 2:1 face at x = 0.75, coarse to fine and fine to
    coarse (today's operators), N = 200, t = 0.2, CPU: the density against the exact solution, and the plateau behind
    the shock with the uniform controls at the coarse and the fine spacing."""
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    for case, lab, mk in (("as-c2f", "c2f (x > 0.75 refined)", "o"), ("as-f2c", "f2c (x < 0.75 refined)", "s")):
        xs, rho, prob = sod_columns(runs / case)
        axs[0].plot(xs, rho, mk, ms=2.5, mfc="none", label=lab)
    xe = np.linspace(0.0, 1.0, 2001)
    axs[0].plot(xe, prob.density(xe), "k-", lw=0.9, label="exact")
    axs[0].axvline(0.75, color="C3", ls="--", lw=0.8, label="2:1 face")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel(r"$\rho$")
    axs[0].grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    axs[0].set_title("Sod at t = 0.2, N = 200")
    sod_zoom(axs[1], runs, (("as-c2f", "c2f", "o"), ("as-uniform-200", "uniform, N = 200", "^"),
                            ("as-uniform-400", "uniform, N = 400", "v")), "the plateau, coarse to fine")
    sod_zoom(axs[2], runs, (("as-f2c", "f2c", "s"), ("as-uniform-200", "uniform, N = 200", "^"),
                            ("as-uniform-400", "uniform, N = 400", "v")), "the plateau, fine to coarse")
    fig.suptitle("AS: a shock crossing a static 2:1 face (mean restriction, tricubic fill, reflux), CPU")
    fig.savefig(out / "shock-crossing.png", dpi=DPI)
    plt.close(fig)


def fig_capturing_error_sod(out: Path, runs: Path) -> None:
    """E0 (#68, issue #85 R2b): the Sod plateau after the c2f crossing, the tricubic and the restriction-compatible
    fills at N = 200 and 400, and the plateau defect (amr_shock_oracle.measure) of every fill x restriction x reflux
    arm of the CPU re-runs."""
    ash = oracle_module("amr-shock", "amr_shock_oracle")
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    sod_zoom(axs[0], runs, (("e0-tricubic-mean-on", "tricubic, N = 200", "o"),
                            ("e0-restriction-compatible-mean-on", "compatible, N = 200", "s"),
                            ("e0-tricubic-400", "tricubic, N = 400", "^"),
                            ("e0-restriction-compatible-400", "compatible, N = 400", "v")),
             "the plateau behind the shock (reflux on, mean restriction)")
    arms = [(f, r) for f in ("restriction-compatible", "tricubic") for r in ("mean", "lagrange6")]
    x = np.arange(len(arms))
    for j, (refl, color) in enumerate((("on", "C0"), ("off", "C1"))):
        d = [ash.measure(runs / f"e0-{f}-{r}-{refl}", 0.75, 0.03)["plateau_r"] for f, r in arms]  # check.sh margin
        axs[1].bar(x + (j - 0.5) * 0.38, d, 0.38, color=color, label=f"reflux {refl}")
    axs[1].set_xticks(x, [f"{'compatible' if f.startswith('restriction') else f} + {r}" for f, r in arms], fontsize=8)
    axs[1].set_ylabel("plateau density defect / jump")
    axs[1].grid(alpha=0.3, axis="y")
    axs[1].legend(fontsize=8)
    axs[1].set_title("Sod, c2f, N = 200 (CPU re-runs)")
    fig.suptitle("E0: the defect behind a shock that crossed a 2:1 face does not converge away and no operator removes it")
    fig.savefig(out / "capturing-error-sod.png", dpi=DPI)
    plt.close(fig)


def fig_quintic_fill(out: Path, runs: Path) -> None:
    """D1 (#68, issue #85 R2b): the ghost probe on a degree-5 polynomial field (box3d, octree, and channel2d, quadtree;
    2 ranks), the largest error per ghost group of today's operators (mean + tricubic: exact to degree 1 and 3) and of
    lagrange6 + quintic (exact to degree 5); and the AO ladder of lagrange6 with the tricubic and the quintic fill."""
    gp = oracle_module("ghosts", "ghost_probe")
    fig, axs = plt.subplots(1, 3, figsize=(20, 6.0), constrained_layout=True, width_ratios=(2.2, 1.1, 1.0))
    for ax, case in zip(axs[:2], ("box3d", "channel2d"), strict=True):
        res = {arm: gp.probe(runs / f"gp5-{arm}-{case}") for arm in ("default", "q6")}
        keys = sorted(set(res["default"]) | set(res["q6"]), key=lambda k: (k[0], k[1]))
        keys = [k for k in keys if k[0] != "solid"]
        x = np.arange(len(keys))
        for j, (arm, lab, color) in enumerate((("default", "mean + tricubic", "C3"), ("q6", "lagrange6 + quintic", "C2"))):
            v = [max((r[0] for r in res[arm].get(k, [])), default=0.0) for k in keys]
            ax.bar(x + (j - 0.5) * 0.4, np.maximum(v, 1e-17), 0.4, color=color, label=lab)
        ax.set_yscale("log")
        ax.set_xticks(x, [f"{c}: {b}" for c, b in keys], rotation=90, fontsize=5.5)
        ax.set_ylabel("largest ghost error / field scale (floor 1e-17)")
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=8)
        ax.set_title(f"ghost probe, degree-5 field, {case}")
    ns, eu = ao_ladder("euler", uniform=True)
    axs[2].loglog(ns, [e["all"][0] for e in eu], "k*-", label="uniform (coarse)")
    for leg, lab, fmt in (("euler-l6", "lagrange6 + tricubic", "s-"), ("euler-q6", "lagrange6 + quintic", "d-")):
        nn, er = ao_ladder(leg)
        v = [e["all"][0] for e in er]
        axs[2].loglog(nn, v, fmt, mfc="none", label=f"{lab} ({math.log2(v[-2] / v[-1]):.2f})")
    axs[2].set_xlabel("base cells per side N")
    axs[2].set_ylabel(r"$L_1(\rho)$, all cells")
    axs[2].grid(alpha=0.3, which="both")
    axs[2].legend(fontsize=8, title="reflux off (128 → 256 order)", title_fontsize=7)
    axs[2].set_title("AO vortex: the quintic fill removes the tricubic error")
    fig.suptitle("D1: the quintic coarse-to-fine fill, in situ (ghost probe) and on the AO vortex (CPU)")
    fig.savefig(out / "quintic-fill.png", dpi=DPI)
    plt.close(fig)


def fig_limited_fills(out: Path, runs: Path) -> None:
    """D2 (#68, issue #85 R2b): the Sod plateau after the crossing with the 6-point fills on lagrange6 (quintic, weno6,
    quintic-mp) against the tricubic fill, coarse to fine and fine to coarse, N = 200, CPU."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for ax, side in zip(axs, ("c2f", "f2c"), strict=True):
        base = (("e0-tricubic-lagrange6-on", "tricubic", "o"),) if side == "c2f" else ()
        sod_zoom(ax, runs, (*base, (f"d-quintic-{side}", "quintic", "s"), (f"d-weno6-{side}", "weno6", "^"),
                            (f"d-quintic-mp-{side}", "quintic-mp", "v")),
                 f"{side}: the plateau behind the shock (lagrange6, reflux on)")
    fig.suptitle("D2: the limited coarse-to-fine fills do not remove the plateau defect (CPU)")
    fig.savefig(out / "limited-fills.png", dpi=DPI)
    plt.close(fig)


def fig_reflux_policy(out: Path, runs: Path) -> None:
    """E2 (#68, issue #85 R2b): reflux_policy = troubled against always. The AO ladder (lagrange6 + quintic) under
    troubled, without reflux and today's operators; and the Sod shock after the crossing, always against troubled."""
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6), constrained_layout=True)
    ns, eu = ao_ladder("euler", uniform=True)
    axs[0].loglog(ns, [e["all"][0] for e in eu], "k*-", label="uniform (coarse)")
    for leg, lab, fmt in (("euler", "mean + tricubic, reflux always", "o-"),
                          ("euler-q6", "lagrange6 + quintic, no reflux", "d-"),
                          ("euler-tq", "lagrange6 + quintic, reflux troubled", "x-")):
        nn, er = ao_ladder(leg)
        v = [e["all"][0] for e in er]
        axs[0].loglog(nn, v, fmt, mfc="none", label=f"{lab} ({math.log2(v[-2] / v[-1]):.2f})")
    axs[0].set_xlabel("base cells per side N")
    axs[0].set_ylabel(r"$L_1(\rho)$, all cells")
    axs[0].grid(alpha=0.3, which="both")
    axs[0].legend(fontsize=7, title="(128 → 256 order)", title_fontsize=7)
    axs[0].set_title("AO vortex: troubled keeps the 6th order")
    for ax, side in zip(axs[1:], ("c2f", "f2c"), strict=True):
        sod_zoom(ax, runs, ((f"d-quintic-{side}", "always", "o"), (f"e2-troubled-{side}", "troubled", "x")),
                 f"Sod {side}: the shock, lagrange6 + quintic", plateau=False)
        _, _, prob = sod_columns(runs / f"d-quintic-{side}")
        ax.set_xlim(prob.x_shock - 0.04, prob.x_shock + 0.03)
        ax.set_ylim(0.1, 0.3)
    fig.suptitle("E2: the reflux policy troubled against always (CPU)")
    fig.savefig(out / "reflux-policy.png", dpi=DPI)
    plt.close(fig)


def fig_integral_drift(out: Path, runs: Path) -> None:  # noqa: ARG001
    """AL (#68, issue #85 R2b): the relative drift of the domain integrals per step, reflux on, off and troubled, on
    the AO vortex and the AS Sod, from the kept work directories of amr-log/check.sh (CPU)."""
    adir = HERE / "amr-log"
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for ax, leg, vars_ in ((axs[0], "vortex", ("r", "rE")), (axs[1], "sod", ("r", "rE"))):
        for run_, lab, ls in (("reflux", "reflux always", "-"), ("noreflux", "reflux off", "--"),
                              ("reflux-troubled", "reflux troubled", ":")):
            path = next((adir / f"{AO_TAG}-{leg}-{run_}").glob("*-conservation_history.dat"))
            lines = path.read_text().splitlines()
            header = lines[0].split("=", 1)[1].replace('"', "").split()
            data = np.array([[float(v) for v in ln.split()] for ln in lines[1:] if ln.strip()])
            for v, c in zip(vars_, ("C0", "C3"), strict=True):
                col = data[:, header.index(f"int_{v}")]
                ax.semilogy(data[:, 0], np.maximum(np.abs(col - col[0]) / abs(col[0]), 1e-17), ls, color=c,
                            label=f"{lab}, {v}")
        ax.set_xlabel("step")
        ax.set_ylabel("relative drift (floor 1e-17)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axs[0].set_title("AO vortex, N = 64, t = 1 (periodic)")
    axs[1].set_title("AS Sod c2f, N = 200, t = 0.2")
    fig.suptitle("AL: the domain integrals across static 2:1 seams (CPU, 2 ranks)")
    fig.savefig(out / "integral-drift.png", dpi=DPI)
    plt.close(fig)

def divb_history(work: Path) -> tuple[list[str], np.ndarray]:
    """Return the column names and the rows of the div(B) history of a run (`[diagnostics] divb_flux`)."""
    lines = next(work.glob("*-divb_history.dat")).read_text().splitlines()
    names = [v.strip('"') for v in lines[0].split("=", 1)[1].split()]
    return names, np.array([[float(v) for v in ln.split()] for ln in lines[1:] if ln.strip()])


def fig_divergence_across_seams(out: Path, runs: Path) -> None:
    """AD (#68, issue #85 R3): the divergence the cleaning acts on (D_flux, from the psi face fluxes) across static 2:1
    seams, FNL: the seam maximum against time on the source-free CPAW (refined N = 64 and 128, the uniform N = 128
    global maximum for scale), and on the A2 field loop under GLM and EGLM (refined N = 64 against uniform N = 128)."""
    adir = HERE / "amr-divb"
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for run_, col, lab, ls in (("cpaw-n64", "seam_max_dflux", "refined N = 64, seam", "-"),
                               ("cpaw-n128", "seam_max_dflux", "refined N = 128, seam", "-"),
                               ("cpaw-uniform-n128", "max_dflux", "uniform N = 128, global", "--")):
        names, rows = divb_history(adir / f"work-adam_flume_fnl-np2-{run_}")
        axs[0].semilogy(rows[1:, 1], rows[1:, names.index(col)], ls, lw=0.8, label=lab)
    axs[0].set_xlabel("t (wave periods)")
    axs[0].set_title("source-free CPAW, GLM: the seam floor stays bounded")
    for cl, c in (("glm", "C0"), ("eglm", "C3")):
        for run_, col, lab, ls in ((f"loop-{cl}-n64", "seam_max_dflux", "refined N = 64, seam", "-"),
                                   (f"loop-{cl}-uniform-n128", "max_dflux", "uniform N = 128, global", ":")):
            names, rows = divb_history(runs / run_)
            axs[1].semilogy(rows[1:, 1], rows[1:, names.index(col)], ls, color=c, lw=0.8, label=f"{cl.upper()}, {lab}")
    axs[1].set_xlabel("t")
    axs[1].set_title("field loop (A2): GLM and EGLM")
    for ax in axs:
        ax.set_ylabel("max |D_flux|")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("AD: the divergence the cleaning acts on, across static 2:1 seams (FNL, 2 ranks)")
    fig.savefig(out / "divergence-across-seams.png", dpi=DPI)
    plt.close(fig)


def fig_glm_damping(out: Path, runs: Path) -> None:
    """G2 (#68, issue #85 R3): the seam D_flux of the source-free CPAW against time for three GLM damping lengths, the
    domain-scale 1.0 (the AD runs), the coarse spacing and min-cell, N = 64 and 128, FNL."""
    adir = HERE / "amr-divb"
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True, sharey=True)
    for ax, n in zip(axs, (64, 128), strict=True):
        for work, lab in ((adir / f"work-adam_flume_fnl-np2-cpaw-n{n}", "1.0 (A2)"),
                          (runs / f"g2-cpaw-coarse-{n}", r"coarse spacing $L_x/N$"),
                          (runs / f"g2-cpaw-mincell-{n}", "min-cell (finest spacing)")):
            names, rows = divb_history(work)
            ax.semilogy(rows[1:, 1], rows[1:, names.index("seam_max_dflux")], lw=0.8, label=lab)
        ax.set_xlabel("t (wave periods)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, title="glm_damping_length", title_fontsize=8)
        ax.set_title(f"refined N = {n}")
    axs[0].set_ylabel("seam max |D_flux|")
    fig.suptitle("G2: the GLM damping length on the source-free CPAW across 2:1 seams (FNL)")
    fig.savefig(out / "glm-damping.png", dpi=DPI)
    plt.close(fig)


def fig_field_loop_seams(out: Path, runs: Path) -> None:
    """G3a and G3b (#68, issue #85 R3): the A2 field loop across the 2:1 box under GLM and EGLM, reflux on and off,
    against the uniform N = 128 run, FNL, today's code: the magnetic energy and <|B_z|>/A0 per checkpoint (every 20
    steps)."""
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for cl, c in (("glm", "C0"), ("eglm", "C3")):
        for suffix, lab, ls in (("n64", "refined, reflux", "-"), ("n64-noreflux", "refined, no reflux", "--"),
                                ("uniform-n128", "uniform N = 128", ":")):
            work = runs / f"loop-{cl}-{suffix}"
            ini = read_ini(next(work.glob("*.ini")))
            ngc = int(ini["grid"]["ngc"])
            its = sorted({int(p.name.split("-")[-2]) for p in work.glob("*-proc*.h5")})
            ts, eb, bz = [], [], []
            for it in its:
                bl = first_plane(blocks(work, ngc, MHD, it=it))
                vol = [float(np.prod(b["d"][:2])) for b in bl]
                eb.append(sum(v * 0.5 * float(np.sum(b["f"]["bx"] ** 2 + b["f"]["by"] ** 2 + b["f"]["bz"] ** 2))
                              for v, b in zip(vol, bl, strict=True)))
                area = sum(v * b["f"]["bz"].size for v, b in zip(vol, bl, strict=True))
                bz.append(sum(v * float(np.sum(np.abs(b["f"]["bz"]))) for v, b in zip(vol, bl, strict=True)) / area)
                ts.append(it)
            a0 = 1.0e-3
            axs[0].plot(ts, np.array(eb) / eb[0], ls, color=c, label=f"{cl.upper()}, {lab}")
            axs[1].semilogy(ts, np.maximum(np.array(bz) / a0, 1e-17), ls, color=c, label=f"{cl.upper()}, {lab}")
    axs[0].set_ylabel("magnetic energy / initial")
    axs[1].set_ylabel(r"$\langle|B_z|\rangle / A_0$")
    for ax in axs:
        ax.set_xlabel("step")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axs[0].set_title("the loop's magnetic energy")
    axs[1].set_title(r"$B_z$ generated by the divergence error ($w = 1$)")
    fig.suptitle("G3a, G3b: the field loop across a 2:1 box, GLM and EGLM, reflux on and off (FNL, today's code)")
    fig.savefig(out / "field-loop-seams.png", dpi=DPI)
    plt.close(fig)


def fig_shock_crossing_rj2a(out: Path, runs: Path) -> None:
    """AS RJ2a (#68, issue #85 R3): the RJ2a fast shock crossing the static 2:1 face at x = 0.75, N = 256, t = 0.2, FNL:
    the density against the exact solution, c2f and f2c, and the plateau behind the shock with the uniform control."""
    ash = oracle_module("amr-shock", "amr_shock_oracle")

    def cols(work: Path) -> tuple[np.ndarray, np.ndarray, object]:
        ini = ash.read_ini(work)
        prob = ash.Problem(ini)
        xs, _, q, _ = ash.load(work, ash.MHD, int(ini["grid"]["ngc"]))
        return xs, q[0], prob

    fig, axs = plt.subplots(1, 2, figsize=(13, 4.6), constrained_layout=True)
    for case, lab, mk in (("as-rj2a-c2f", "c2f", "o"), ("as-rj2a-f2c", "f2c", "s"),
                          ("as-rj2a-uniform", "uniform", "^")):
        xs, rho, prob = cols(runs / case)
        for ax in axs:
            ax.plot(xs, rho, mk, ms=2.5, mfc="none", label=lab)
    xe = np.linspace(0.0, 1.0, 4001)
    for ax in axs:
        ax.plot(xe, prob.density(xe), "k-", lw=0.9, label="exact")
        ax.axvline(0.75, color="C3", ls="--", lw=0.8, label="2:1 face")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$\rho$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    axs[1].set_xlim(prob.x_left - 0.02, prob.x_shock + 0.02)
    axs[1].set_ylim(prob.q2[0] - 0.03, prob.q2[0] + 0.03)
    axs[0].set_title("Ryu-Jones 2a at t = 0.2, N = 256")
    axs[1].set_title("the plateau behind the fast shock")
    fig.suptitle("AS: the RJ2a fast shock crossing a static 2:1 face (mean restriction, tricubic fill, reflux), FNL")
    fig.savefig(out / "shock-crossing-rj2a.png", dpi=DPI)
    plt.close(fig)


def fig_blast_symmetry(out: Path, runs: Path) -> None:
    """#77 (issue #85 R3): the Balsara-Spicer MHD blast (EGLM, positivity limiter), uniform N = 128, FNL: the density at
    t = 0.01, the map of its point-symmetry defect about the centre on the splitting scheme, and the defect of density
    and B_y against the step for the splitting scheme and weno-riemann HLLD."""
    def grid(work: Path, it: int | None, name: str) -> np.ndarray:
        ini = read_ini(next(work.glob("*.ini")))
        ngc = int(ini["grid"]["ngc"])
        bl = first_plane(blocks(work, ngc, MHD, it=it))
        d = min(float(b["d"][0]) for b in bl)
        n = round(1.0 / d)
        g = np.zeros((n, n))
        for b in bl:
            i0, j0 = round(b["lo"][0] / d), round(b["lo"][1] / d)
            a = b["f"][name][:, :, 0]
            g[i0:i0 + a.shape[0], j0:j0 + a.shape[1]] = a
        return g

    def defect(g: np.ndarray, parity: float) -> np.ndarray:
        return np.abs(g - parity * g[::-1, ::-1]) / np.abs(g).max()

    fig, axs = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)
    ws = runs / "blast77-split"
    rho = grid(ws, None, "r")
    im = axs[0].imshow(rho.T, origin="lower", extent=(0, 1, 0, 1), cmap="magma")
    fig.colorbar(im, ax=axs[0], label=r"$\rho$", shrink=0.85)
    axs[0].set_title("splitting scheme, density at t = 0.01")
    dmap = defect(rho, 1.0)
    im = axs[1].imshow(np.log10(np.maximum(dmap, 1e-17)).T, origin="lower", extent=(0, 1, 0, 1), cmap="viridis")
    fig.colorbar(im, ax=axs[1], label=r"$\log_{10}$ point-symmetry defect of $\rho$", shrink=0.85)
    axs[1].set_title(f"splitting: the defect map (max {dmap.max():.1e})")
    for scheme, c in (("split", "C0"), ("hlld", "C3")):
        work = runs / f"blast77-{scheme}"
        its = sorted({int(p.name.split("-")[-2]) for p in work.glob("*-proc*.h5")})
        for name, parity, ls in (("r", 1.0, "-"), ("by", 1.0, "--")):
            dd = [float(defect(grid(work, it, name), parity).max()) for it in its]
            axs[2].semilogy(its, np.maximum(dd, 1e-17), ls, color=c,
                            label=f"{'splitting' if scheme == 'split' else 'HLLD'}, {'density' if name == 'r' else '$B_y$'}")
    axs[2].axhline(1e-10, color="k", ls=":", lw=0.8, label="av3c bound (HLLD)")
    axs[2].set_xlabel("step")
    axs[2].set_ylabel("max point-symmetry defect (floor 1e-17)")
    axs[2].grid(alpha=0.3)
    axs[2].legend(fontsize=8)
    axs[2].set_title("the defect against the step")
    fig.suptitle(r"#77: the MHD blast ($\beta = 2.5\cdot10^{-4}$, EGLM, limiter) loses its point symmetry on the "
                 "flux-splitting path (FNL, N = 128)")
    fig.savefig(out / "mhd-blast-symmetry.png", dpi=DPI)
    plt.close(fig)

def fig_mhd_conservation(out: Path, runs: Path) -> None:  # noqa: ARG001
    """MV-11 (issue #85, R4): MHD conservation across 2:1 AMR faces, the periodic box with the refined octant, from the
    kept conservation histories of mhd/conservation/check.sh (CPU): the relative drift of every volume integral per
    step for GLM with reflux, damped GLM with reflux, GLM without reflux and no cleaning with reflux."""
    cdir = HERE / "mhd" / "conservation"
    variants = (("glm-reflux", "GLM, reflux"), ("glm-damped-reflux", "damped GLM, reflux"),
                ("glm-noreflux", "GLM, no reflux"), ("none-reflux", "no cleaning, reflux"))
    fig, axs = plt.subplots(1, 4, figsize=(20, 4.6), constrained_layout=True, sharey=True)
    for ax, (run_, title) in zip(axs, variants, strict=True):
        path = next((cdir / f"{TAG}-{run_}").glob("*-conservation_history.dat"))
        lines = path.read_text().splitlines()
        names = [v.strip('"').removeprefix("int_") for v in lines[0].split("=", 1)[1].split()][2:]
        rows = np.array([[float(v) for v in ln.split()] for ln in lines[1:] if ln.strip()])
        for i, v in enumerate(names):
            col = rows[:, 2 + i]
            scale = max(abs(col[0]), np.abs(rows[0, 2:]).max() * 1e-12) if v != "psi" else 1.0
            ax.semilogy(rows[:, 0], np.maximum(np.abs(col - col[0]) / scale, 1e-17), label=v)
        ax.set_xlabel("step")
        ax.grid(alpha=0.3)
        ax.set_title(title)
    axs[0].set_ylabel("|I(t) - I(0)| / |I(0)| (psi: absolute; floor 1e-17)")
    axs[0].legend(fontsize=7, ncol=2)
    fig.suptitle("MV-11: the MHD volume integrals across the 2:1 faces of the refined octant (CPU, 2 ranks)")
    fig.savefig(out / "mhd-conservation.png", dpi=DPI)
    plt.close(fig)


FIGURES = {"sod": fig_sod, "lax": fig_lax, "shu-osher": fig_shu_osher, "vortex": fig_vortex,
           "shock-cylinder": fig_cylinder, "step": fig_step, "step-trees": fig_step_trees, "ghosts": fig_ghosts,
           "conservation": fig_conservation,
           "orszag-tang": fig_orszag_tang,
           "rotor": fig_rotor, "field-loop": fig_field_loop, "mhd-riemann": fig_mhd_riemann, "glm-pulse": fig_glm_pulse,
           "blast": fig_blast, "near-vacuum": fig_vacuum, "eglm-energy": fig_eglm_energy, "order": fig_order,
           "viscous": fig_viscous, "viscous-fields": fig_viscous_fields, "couette": fig_couette,
           "regrid": fig_regrid, "gravity-atmospheres": fig_gravity_atmospheres, "gravity-force": fig_gravity_force,
           "gravity-order": fig_gravity_order, "gravity-balance": fig_gravity_balance,
           "gravity-rt-linear": fig_gravity_rt_linear, "gravity-rt": fig_gravity_rt,
           "gravity-rt-symmetry": fig_gravity_rt_symmetry, "ranks": fig_ranks,
           "selfgravity-potential": fig_selfgravity_potential,
           "selfgravity-poisson": fig_selfgravity_poisson,
           "selfgravity-jeans": fig_selfgravity_jeans,
           "selfgravity-isolated": fig_selfgravity_isolated,
           "selfgravity-star": fig_selfgravity_star,
           "selfgravity-amr": fig_selfgravity_amr, "selfgravity-amr-jeans": fig_selfgravity_amr_jeans,
           "selfgravity-amr-regrid": fig_selfgravity_amr_regrid,
           "selfgravity-amr-isolated": fig_selfgravity_amr_isolated,
           "selfgravity-star-amr": fig_selfgravity_star_amr,
           "viscous-ladders": fig_viscous_ladders, "seam-proto-b1": fig_seam_proto_b1,
           "seam-proto-b2": fig_seam_proto_b2, "seam-proto-b3": fig_seam_proto_b3,
           "capturing-error": fig_capturing_error, "multirealm": fig_multirealm,
           "rotated-shock-tube": fig_rotated_shock_tube, "weno-riemann": fig_weno_riemann, "quadtree": fig_quadtree,
           "regrid-runtime": fig_regrid_runtime, "seam-conservation": fig_seam_conservation,
           "order-across-seams": fig_order_across_seams, "shock-crossing": fig_shock_crossing,
           "capturing-error-sod": fig_capturing_error_sod, "quintic-fill": fig_quintic_fill,
           "limited-fills": fig_limited_fills, "reflux-policy": fig_reflux_policy, "integral-drift": fig_integral_drift,
           "divergence-across-seams": fig_divergence_across_seams, "glm-damping": fig_glm_damping,
           "field-loop-seams": fig_field_loop_seams, "shock-crossing-rj2a": fig_shock_crossing_rj2a,
           "mhd-blast-symmetry": fig_blast_symmetry, "mhd-conservation": fig_mhd_conservation}


def main() -> int:
    """Draw the requested figures; report the missing data."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, required=True, help="output directory of the PNG files")
    ap.add_argument("--runs", type=Path, default=HERE, help="directory of the re-run MHD cases")
    ap.add_argument("--only", nargs="*", default=None, help="figures to draw (default: all)")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    status = 0
    for name, fn in FIGURES.items():
        if a.only and name not in a.only:
            continue
        try:
            fn(a.out, a.runs)
            print(f"{name}: done")
        except (FileNotFoundError, StopIteration, KeyError) as err:
            print(f"{name}: skipped ({err!r})")
            status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
