#!/usr/bin/env python3
"""The wall cells of a wall-modelled channel run of issue #119 T3, against time: the Mach 3 first-cell drain diagnostic.

`cells <work>...`: per checkpoint `channel-<it>-proc*.h5` (time from the `.xdmf`), plane (x, z) averages folded over both walls: the driving force
(= tau_w at a steady state), the Favre mean u, T / T_w and rho of the first cell, the share of the mass in the two first-cell layers, the core rho
and u, and the molecular plus resolved shear stress `mu dU/dy - <rho u'' v''>` at the face between cells 1 and 2 over the driving force.
`budget <work>...`: the total-energy budget in 20-unit windows, `dE/dt = W - Q_w`, the forcing work `W = f int(u dV)` with int(u dV) taken as
int(rho u dV) times their ratio in the last checkpoint (the wall heat flux Q_w is not logged: it is the residual). Per unit wall area.

Checkpoint blocks are (z, y, x) arrays with the ghosts; R = 1, gamma 1.4, the power-law viscosity mu = (T / T_w)^0.7 / Re of the Coleman inputs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np

GAMMA, NGC, OMEGA = 1.4, 3, 0.7


def ini_value(work: Path, key: str) -> float:
    """Return the first `key = value` of the run's input.ini."""
    for line in (work / "input.ini").read_text().splitlines():
        if line.split("=")[0].strip() == key:
            return float(line.split("=")[1])
    raise KeyError(key)


def history(path: Path) -> np.ndarray:
    """Load a history; restarts repeat rows, keep the last of each step."""
    f = np.loadtxt(path, ndmin=2, comments="VARIABLES")
    _, last = np.unique(f[::-1, 0], return_index=True)
    return f[::-1][last]


def checkpoints(work: Path) -> list[tuple[float, str]]:
    """Return (time, iteration tag) of every checkpoint of the run, in time order; the restart file (`channel-restart`, other fields) is not one."""
    out = []
    for x in work.glob("channel-*.xdmf"):
        m = re.search(r'<Time Value="([^"]+)"', x.read_text()[:4000])
        if m and x.stem.split("-")[1].isdigit():
            out.append((float(m.group(1)), x.stem.split("-")[1]))
    return sorted(out)


def layers(work: Path, tag: str) -> tuple[dict[float, list[np.ndarray]], dict[tuple[float, float], list[np.ndarray]], float, float]:
    """Return {distance from the wall: [rho, u, v away from the wall, T]} of every cell layer with both walls folded, the same per wall keyed by
    (side, distance), side -1 the y_min wall and +1 the y_max one, dy, and the volume."""
    acc: dict[tuple[float, float], list[list[np.ndarray]]] = {}
    vol, dy = 0.0, 0.0
    for f5 in sorted(work.glob(f"channel-{tag}-proc*.h5")):
        with h5py.File(f5, "r") as h:
            for b in {k.rsplit("-", 1)[0] for k in h}:
                o, d = h[b + "-origin"][()], h[b + "-dxdydz"][()]
                s = (slice(NGC, -NGC),) * 3
                r, ru, rv, rw, re_ = (h[f"{b}-{v}"][()][s] for v in ("r", "ru", "rv", "rw", "rE"))
                p = (GAMMA - 1.0) * (re_ - 0.5 * (ru**2 + rv**2 + rw**2) / r)
                vol += r.size * float(np.prod(d))
                dy = float(d[1])
                for j in range(r.shape[1]):
                    y = o[1] + (j + NGC + 0.5) * d[1]
                    sgn = 1.0 if y < 0.0 else -1.0
                    z = [r[:, j, :].ravel(), (ru[:, j, :] / r[:, j, :]).ravel(), sgn * (rv[:, j, :] / r[:, j, :]).ravel(), (p[:, j, :] / r[:, j, :]).ravel()]
                    acc.setdefault((-sgn, round(1.0 - abs(y), 10)), []).append(z)
    side = {k: [np.concatenate(c) for c in zip(*v)] for k, v in acc.items()}
    fold: dict[float, list[list[np.ndarray]]] = {}
    for (_, dist), z in side.items():
        fold.setdefault(dist, []).append(z)
    return {k: [np.concatenate(c) for c in zip(*v)] for k, v in fold.items()}, side, dy, vol


def cells(works: list[str]) -> int:
    """Print the wall-cell series of each run."""
    for w in works:
        work = Path(w)
        tw, reynolds = ini_value(work, "wall_temperature"), ini_value(work, "reynolds")
        f = history(work / "channel-forcing_history.dat")
        print(f"== {work.name}  T_w {tw:.5f}  Re {reynolds:g}")
        print("     t   f(+-2.5)    u1      T1/Tw   rho1   mass1  rho_c  u_c    stress12/f |  u1 lo|hi         T1/Tw lo|hi")
        for t, tag in checkpoints(work):
            lay, side, dy, _ = layers(work, tag)
            ks = sorted(lay)
            (r1, u1, v1, t1), (r2, u2, v2, t2) = lay[ks[0]], lay[ks[1]]
            rc, uc = lay[ks[-1]][0], lay[ks[-1]][1]
            fm = f[(f[:, 1] > t - 2.5) & (f[:, 1] <= t + 2.5), 2]
            fm = float(fm.mean()) if fm.size else float("nan")
            fav = lambda r, x: float((r * x).mean() / r.mean())  # noqa: E731
            uf1, uf2, tf1, tf2 = fav(r1, u1), fav(r2, u2), fav(r1, t1), fav(r2, t2)
            mu = (0.5 * (tf1 + tf2) / tw) ** OMEGA / reynolds
            rs = [-float((r * (u - fav(r, u)) * v).mean()) for r, u, v in ((r1, u1, v1), (r2, u2, v2))]
            stress = mu * (uf2 - uf1) / dy + 0.5 * (rs[0] + rs[1])
            mass1 = float(r1.mean()) * dy  # the two first-cell layers (2 dy of the 2-high channel) over the total mass, the mean density being 1
            lo, hi = side[(-1.0, ks[0])], side[(1.0, ks[0])]
            print(
                f" {t:6.1f}  {fm:.3e}  {uf1:+.4f}  {tf1 / tw:.3f}  {r1.mean():6.3f}  {mass1:.3f}  {rc.mean():.3f}  {fav(rc, uc):.3f}  {stress / fm:+.3f}"
                f"     |  {fav(lo[0], lo[1]):+.4f} {fav(hi[0], hi[1]):+.4f}  {fav(lo[0], lo[3]) / tw:.3f} {fav(hi[0], hi[3]) / tw:.3f}"
            )
    return 0


def budget(works: list[str]) -> int:
    """Print the total-energy budget windows of each run."""
    for w in works:
        work = Path(w)
        f, c = history(work / "channel-forcing_history.dat"), history(work / "channel-conservation_history.dat")
        t, tag = checkpoints(work)[-1]
        lay, _, _, vol = layers(work, tag)
        # every block on the same level (the Coleman inputs): unweighted cell sums are volume integrals up to a common factor
        ratio = sum(float(u.sum()) for _, u, _, _ in lay.values()) / sum(float((r * u).sum()) for r, u, _, _ in lay.values())
        area2 = vol  # the two walls, each of area V / 2
        ft = np.interp(c[:, 1], f[:, 1], f[:, 2])
        work_rate = ft * c[:, 3] * ratio
        print(f"== {work.name}  int(u)/int(rho u) at t = {t:g}: {ratio:.4f}")
        print(" window    f          W/2A       dE/dt/2A    Qw/2A      Qw/W")
        for a in range(0, int(c[-1, 1] + 1e-9), 20):
            m = (c[:, 1] > a) & (c[:, 1] <= a + 20.0)
            if m.sum() < 3:
                continue
            de = np.polyfit(c[m, 1], c[m, 6], 1)[0]
            wr = float(work_rate[m].mean())
            print(f" {a:3d}-{a + 20:3d}  {ft[m].mean():.3e}  {wr / area2:.3e}  {de / area2:+.3e}  {(wr - de) / area2:.3e}  {(wr - de) / wr:.3f}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["cells", "budget"])
    ap.add_argument("works", nargs="+")
    a = ap.parse_args()
    return cells(a.works) if a.what == "cells" else budget(a.works)


if __name__ == "__main__":
    sys.exit(main())
