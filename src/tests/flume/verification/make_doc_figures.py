#!/usr/bin/env python3
"""Draw the figures of the FLUME documentation (docs/applications/flume/) from verification runs.

Why this script exists: the documentation shows real FLUME results, not sketches. Every figure is drawn from the
checkpoints of a verification input (the case directories of this tree, or the same inputs re-run into `--runs`), with
the exact or reference solution of the corresponding oracle where one exists, so a figure can be regenerated after any
change of the solver and compared with the published one.

Usage:
    make_doc_figures.py --out DIR [--runs DIR] [--only NAME ...]

`--runs` holds the re-runs of inputs whose work directories keep no checkpoint (the MHD cases), one sub-directory per
case: orszag-tang, orszag-tang-amr, rotor, field-loop, brio-wu, rj4d, glm-pulse, and the cases of `make_doc_runs.sh`
(rj2a, rj2a-hlld, the limited blasts, the Euler near-vacuum problems, the linear-wave ladders). A figure whose data is
missing is skipped with a message. Needs numpy, h5py and matplotlib.
"""

from __future__ import annotations

import argparse
import configparser
import math
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


def last_files(work: Path) -> list[Path]:
    """Return the checkpoint files of the last saved iteration of a run."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        raise FileNotFoundError(f"no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    return [p for p in files if int(p.name.split("-")[-2]) == last]


def blocks(work: Path, ngc: int, names: tuple[str, ...]) -> list[dict]:
    """Return the blocks of the last checkpoint: interior fields (x, y, z order), cell spacing and interior corner."""
    out = []
    for path in last_files(work):
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


def profile(work: Path, ngc: int, names: tuple[str, ...]) -> tuple[int, np.ndarray, np.ndarray]:
    """Return the active axis, the cell centres and the fields of a 1-D run (first transverse copy)."""
    cells: dict[float, np.ndarray] = {}
    axis = -1
    for b in blocks(work, ngc, names):
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
    axs[0].set_title("density after one period, 128², HLLC")
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
                     "Field loop advected across the periodic box, 128²")


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


FIGURES = {"sod": fig_sod, "lax": fig_lax, "shu-osher": fig_shu_osher, "vortex": fig_vortex,
           "shock-cylinder": fig_cylinder, "conservation": fig_conservation, "orszag-tang": fig_orszag_tang,
           "rotor": fig_rotor, "field-loop": fig_field_loop, "mhd-riemann": fig_mhd_riemann, "glm-pulse": fig_glm_pulse,
           "blast": fig_blast, "near-vacuum": fig_vacuum, "eglm-energy": fig_eglm_energy, "order": fig_order}


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
