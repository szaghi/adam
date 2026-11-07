#!/usr/bin/env python3
"""Judge the shock-generator runs (issue #121, S3) from their last checkpoint, against the oblique-shock relations.

Subcommands:
    theory --mach M (--deflection T | --angle B)          print the shock angle, deflection and the post-shock ratios
    post <work> --mach M (--deflection T | --angle B) --shock-x X --bound E --angle-bound A
                                                          the uniform region behind the shock against Rankine-Hugoniot: the largest relative
                                                          error of p and rho below E, of the flow angle below A degrees. The region: 16 cells
                                                          behind the captured shock (its oscillations), 4 cells below the streamline from
                                                          the shock's foot (an entropy streak the jump of the ghost state leaves there, printed
                                                          apart) and above y = 0.2 (where the shock leaves through the extrapolated bottom)
    angle <work> --mach M (--deflection T | --angle B) --shock-x X --bound D
                                                          the shock angle fitted from where p crosses (p1 + p2) / 2 along x, y in
                                                          [0.15, 0.85]: within D degrees of the theory, and so its position
    bitwise <a> <b>                                       the last checkpoints of two runs equal bit for bit

The free stream is rho = 1, p = 1 / gamma (c = 1), at Mach M along x (make_shockgen.py); the shock enters through y = 1.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

GAMMA = 1.4
NGC = 3
P_INF = 1.0 / GAMMA


def deflection(m: float, beta: float) -> float:
    """The theta-beta-M relation (radians)."""
    return math.atan(2.0 / math.tan(beta) * (m * m * math.sin(beta) ** 2 - 1.0) / (m * m * (GAMMA + math.cos(2.0 * beta)) + 2.0))


def theory(m: float, theta_deg: float | None, beta_deg: float | None) -> dict:
    """Shock angle, deflection (degrees) and the post-shock ratios; the weak solution when the deflection is given."""
    if beta_deg is not None:
        beta = math.radians(beta_deg)
    else:
        theta = math.radians(theta_deg)
        mu = math.asin(1.0 / m)
        grid = np.linspace(mu + 1e-9, math.pi / 2, 200001)
        th = np.array([deflection(m, b) for b in grid[::100]])
        bmax = grid[::100][int(np.argmax(th))]
        lo, hi = mu, bmax
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if deflection(m, mid) < theta else (lo, mid)
        beta = 0.5 * (lo + hi)
    theta = deflection(m, beta)
    mn1 = m * math.sin(beta)
    p21 = 1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (mn1 * mn1 - 1.0)
    r21 = (GAMMA + 1.0) * mn1 * mn1 / ((GAMMA - 1.0) * mn1 * mn1 + 2.0)
    return {"beta": math.degrees(beta), "theta": math.degrees(theta), "p21": p21, "r21": r21}


def last_checkpoint(work: Path) -> str:
    """Return the iteration tag of the latest checkpoint (by time)."""
    best = None
    for x in work.glob("sg-*.xdmf"):
        m = re.search(r'<Time Value="([^"]+)"', x.read_text()[:4000])
        if m and x.stem.split("-")[1].isdigit() and (best is None or float(m.group(1)) > best[0]):
            best = (float(m.group(1)), x.stem.split("-")[1])
    if best is None:
        raise SystemExit(f"shockgen_oracle: no checkpoint in {work}")
    return best[1]


def cells(work: Path) -> dict:
    """The interior cells of the last checkpoint: x, y, rho, p and the flow angle in the x-y plane (degrees), flattened."""
    import h5py

    tag = last_checkpoint(work)
    out: dict[str, list] = {k: [] for k in ("x", "y", "r", "p", "phi")}
    s = (slice(NGC, -NGC),) * 3
    for f5 in sorted(work.glob(f"sg-{tag}-proc*.h5")):
        with h5py.File(f5, "r") as h:
            for b in sorted({k.rsplit("-", 1)[0] for k in h}):
                o, d = h[f"{b}-origin"][()], h[f"{b}-dxdydz"][()]  # (z, y, x), the origin of the ghost-extended nodes
                r, ru, rv, rw, re_ = (h[f"{b}-{v}"][()][s] for v in ("r", "ru", "rv", "rw", "rE"))
                nz, ny, nx = r.shape
                z, y, x = np.meshgrid(*(o[a] + (np.arange(n) + NGC + 0.5) * d[a] for a, n in zip((0, 1, 2), (nz, ny, nx), strict=True)),
                                      indexing="ij")
                out["x"].append(x.ravel())
                out["y"].append(y.ravel())
                out["r"].append(r.ravel())
                out["p"].append(((GAMMA - 1.0) * (re_ - 0.5 * (ru**2 + rv**2 + rw**2) / r)).ravel())
                out["phi"].append(np.degrees(np.arctan2(rv, ru)).ravel())
                out["dx"] = d[2]
    return {k: (np.concatenate(v) if isinstance(v, list) else v) for k, v in out.items()}


def post(work: Path, th: dict, x0: float, bound: float, angle_bound: float) -> int:
    """The uniform region behind the shock against Rankine-Hugoniot; the entropy streak from the shock's foot printed apart."""
    c = cells(work)
    dx = c["dx"]
    behind = (c["x"] - (x0 + (1.0 - c["y"]) / math.tan(math.radians(th["beta"])))) / dx
    below = (1.0 - (c["x"] - x0) * math.tan(math.radians(th["theta"])) - c["y"]) / dx
    ep = np.abs(c["p"] / P_INF / th["p21"] - 1.0)
    er = np.abs(c["r"] / th["r21"] - 1.0)
    ea = np.abs(c["phi"] + th["theta"])
    m = (behind > 16.0) & (below > 4.0) & (c["y"] > 0.2) & (c["x"] < 2.0 - 4.0 * dx)
    s = (behind > 16.0) & (np.abs(below) < 3.0) & (c["x"] < 2.0 - 4.0 * dx)
    ok = ep[m].max() < bound and er[m].max() < bound and ea[m].max() < angle_bound
    print(f"   {work.name}: {int(m.sum())} cells; |p/p2 - 1| {ep[m].max():.2e}, |rho/rho2 - 1| {er[m].max():.2e} (bound {bound:g}),"
          f" |angle + theta| {ea[m].max():.2e} deg (bound {angle_bound:g}); the streak from the foot |rho/rho2 - 1| {er[s].max():.2e}"
          f" (p2/p1 {th['p21']:.6f}, rho2/rho1 {th['r21']:.6f}, theta {th['theta']:.4f})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def angle(work: Path, th: dict, x0: float, bound: float) -> int:
    """The shock angle fitted from the mid-pressure crossing of every row."""
    c = cells(work)
    pm = P_INF * 0.5 * (1.0 + th["p21"])
    ys, xc = [], []
    for yv in np.unique(np.round(c["y"], 12)):
        if not 0.15 <= yv <= 0.85:
            continue
        m = np.isclose(c["y"], yv)
        order = np.argsort(c["x"][m])
        x, p = c["x"][m][order], c["p"][m][order]
        xu, iu = np.unique(x, return_index=True)  # the z cells of a row share x: average them
        pu = np.array([p[x == v].mean() for v in xu])
        k = np.nonzero((pu[:-1] < pm) & (pu[1:] >= pm))[0]
        if len(k) == 0:
            continue
        k = k[0]
        ys.append(yv)
        xc.append(xu[k] + (pm - pu[k]) / (pu[k + 1] - pu[k]) * (xu[k + 1] - xu[k]))
    slope, icpt = np.polyfit(ys, xc, 1)
    beta = math.degrees(math.atan(-1.0 / slope))
    x_top = slope * 1.0 + icpt
    d = abs(beta - th["beta"])
    ok = d < bound
    print(f"   {work.name}: shock angle {beta:.3f} deg (theory {th['beta']:.3f}, {len(ys)} rows), at y = 1 x = {x_top:.4f} (shock_x {x0:g});"
          f" |error| {d:.3f} deg (bound {bound:g})  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def bitwise(a: Path, b: Path) -> int:
    """The last checkpoints of `a` and `b`, every dataset bit for bit."""
    import h5py

    ta, tb = last_checkpoint(a), last_checkpoint(b)
    n = bad = 0
    for fa in sorted(a.glob(f"sg-{ta}-proc*.h5")):
        fb = b / fa.name.replace(ta, tb)
        with h5py.File(fa, "r") as ha, h5py.File(fb, "r") as hb:
            for k in ha:
                n += 1
                if not np.array_equal(ha[k][()].view(np.uint8), hb[k][()].view(np.uint8)):
                    bad += 1
    ok = n > 0 and bad == 0
    print(f"   {a.name} vs {b.name}: {n} datasets, {bad} differ  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Dispatch the subcommand."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("theory", "post", "angle"):
        q = sub.add_parser(name)
        if name != "theory":
            q.add_argument("work", type=Path)
            q.add_argument("--shock-x", type=float, required=True)
            q.add_argument("--bound", type=float, required=True)
            if name == "post":
                q.add_argument("--angle-bound", type=float, required=True)
        q.add_argument("--mach", type=float, required=True)
        g = q.add_mutually_exclusive_group(required=True)
        g.add_argument("--deflection", type=float)
        g.add_argument("--angle", type=float)
    b = sub.add_parser("bitwise")
    b.add_argument("a", type=Path)
    b.add_argument("b", type=Path)
    a = ap.parse_args()
    if a.cmd == "bitwise":
        return bitwise(a.a, a.b)
    th = theory(a.mach, a.deflection, a.angle)
    if a.cmd == "theory":
        print(f"beta {th['beta']:.6f} deg, theta {th['theta']:.6f} deg, p2/p1 {th['p21']:.8f}, rho2/rho1 {th['r21']:.8f}")
        return 0
    if a.cmd == "post":
        return post(a.work, th, a.shock_x, a.bound, a.angle_bound)
    return angle(a.work, th, a.shock_x, a.bound)


if __name__ == "__main__":
    sys.exit(main())
