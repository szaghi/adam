#!/usr/bin/env python3
"""Oracle of the FLUME stratified atmospheres (issue #73, M7): the hydrostatic state and its ghost cells.

    atmosphere_oracle.py profile <work> [--tol T]
    atmosphere_oracle.py freefall <work> [--tol T]
    atmosphere_oracle.py energy <work> [--tol T] [--report-total]
    atmosphere_oracle.py rest <work> [<work> ...]
    atmosphere_oracle.py agree <workA> <workB> [--tol T]
    atmosphere_oracle.py wave <work> [<work> ...] [--min-order P]

`profile` reads the last checkpoint of a run written by make_atmosphere.py (fields with their ghosts) and checks, with
the constant gravity g and the region-1 state (rho_1, p_1 at the origin) of its input:

- the interior cells against the analytic atmosphere at their centres: isothermal `rho = rho_1 exp(g . x / T_1)`,
  `p = rho T_1`, `T_1 = p_1 / rho_1`; polytropic `rho = rho_1 s^(1/(G-1))`, `p = p_1 s^(G/(G-1))`,
  `s = 1 + (G-1)/G rho_1/p_1 g . x`; the velocity and the field uniform (region 1's);
- the ghosts beyond `hydrostatic` faces against the boundary rule, the isothermal continuation from the first interior
  cell along the normal, `rho = rho_1' exp(g_d (x_d - x1_d) / T_1')` (for the isothermal atmosphere the analytic
  profile itself, so every hydrostatic ghost is checked, the realm edges and corners too; for the polytropic one the
  ghosts beyond one face only), the velocity zero (the normal component odd of a resting state), the field and psi of
  region 1 (even);
- the ghosts across periodic faces against the analytic profile at the wrapped position.

Ghosts inside the domain (block interfaces, 2:1 seams: their restriction and prolongation are not exact on an
exponential) are not checked here (ghosts/ghost_probe.py covers the exchange), nor the boundary ghosts whose projection
on the domain lies outside their block (their donor is such an interface ghost); both are counted. A NaN anywhere checked fails (the
cases poison the ghosts). Exit 0 when every relative error is within `--tol` (default 1e-13).
"""

from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

import h5py
import numpy as np

FIELDS = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz", "psi")
FACES = ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")


def read_ini(work: Path) -> configparser.ConfigParser:
    """Read the case INI of a work directory."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.read(next(work.glob("*.ini")))
    return ini


def gravity(ini: configparser.ConfigParser) -> np.ndarray:
    """Return g from [gravity]: the components, or 1 / froude^2 along the signed axis."""
    sec = ini["gravity"]
    if "froude" in sec:
        axis = sec["axis"].split(";")[0].strip()
        g = np.zeros(3)
        g["xyz".index(axis[1])] = (1.0 if axis[0] == "+" else -1.0) / float(sec["froude"].split(";")[0]) ** 2
        return g
    return np.array([float(sec.get(f"g_{x}", "0").split(";")[0]) for x in "xyz"])


class Atmosphere:
    """The analytic atmosphere of an input."""

    def __init__(self, ini: configparser.ConfigParser) -> None:
        reg = ini["initial_conditions_region_1"]
        self.rho1, self.p1 = float(reg["r"]), float(reg["p"])
        self.vel = np.array([float(reg[k]) for k in ("u", "v", "w")])
        self.mhd = ini["physics"]["physical_model"].split(";")[0].strip() == "mhd-ideal"
        self.b = np.array([float(reg[k]) for k in ("bx", "by", "bz")]) if self.mhd else np.zeros(3)
        self.g = gravity(ini)
        self.gamma = float(ini["physics"]["gamma"].split(";")[0])
        ic = ini["initial_conditions"]
        self.polytropic = ic["type"].split(";")[0].strip() == "polytropic-atmosphere"
        self.G = float(ic["polytropic_gamma"].split(";")[0]) if self.polytropic else 1.0

    def state(self, x: np.ndarray) -> tuple[float, float]:
        """Density and pressure at the point x."""
        gx = float(np.dot(self.g, x))
        if self.polytropic:
            s = 1.0 + (self.G - 1.0) / self.G * self.rho1 / self.p1 * gx
            return self.rho1 * s ** (1.0 / (self.G - 1.0)), self.p1 * s ** (self.G / (self.G - 1.0))
        rho = self.rho1 * math.exp(gx * self.rho1 / self.p1)
        return rho, rho * self.p1 / self.rho1


def primitive(q: dict[str, float], gamma: float, mhd: bool, eglm: bool) -> tuple[float, np.ndarray, float, np.ndarray]:
    """Density, velocity, pressure and field of a conservative state."""
    rho = q["r"]
    u = np.array([q["ru"], q["rv"], q["rw"]]) / rho
    b = np.array([q["bx"], q["by"], q["bz"]]) if mhd else np.zeros(3)
    e_mag = 0.5 * float(np.dot(b, b)) + (0.5 * q["psi"] ** 2 if eglm else 0.0)
    p = (gamma - 1.0) * (q["rE"] - 0.5 * rho * float(np.dot(u, u)) - e_mag)
    return rho, u, p, b


def profile(args: argparse.Namespace) -> int:
    """Check the interior and the boundary ghosts of the last checkpoint against the atmosphere."""
    ini = read_ini(args.work)
    atm = Atmosphere(ini)
    eglm = atm.mhd and ini["mhd"]["divergence_control"].split(";")[0].strip() == "eglm"
    ngc = int(ini["grid"]["ngc"])
    lo = np.array([float(ini["grid"][f"emin_{x}"]) for x in "xyz"])
    hi = np.array([float(ini["grid"][f"emax_{x}"]) for x in "xyz"])
    active = np.array([ini["grid"][f"null_{x}"].strip() != ".true." for x in "xyz"])
    bc = [ini[f"bc_{f}"]["type"].split(";")[0].strip() for f in FACES]
    basename = next(args.work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
    files = sorted(args.work.glob(f"{basename}-*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    err = {"interior rho, p": 0.0, "interior u, B": 0.0, "hydrostatic ghosts rho, p": 0.0,
           "hydrostatic ghosts u, B": 0.0, "periodic ghosts rho, p": 0.0}
    counts = dict.fromkeys(err, 0)
    skipped, nan, beside = 0, 0, 0
    scale_u = math.sqrt(atm.gamma * atm.p1 / atm.rho1) + float(np.linalg.norm(atm.vel))
    scale_b = max(float(np.linalg.norm(atm.b)), 1.0)

    def note(key: str, value: float) -> None:
        err[key] = max(err[key], value)
        counts[key] += 1

    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                arr = {n: h5[f"{blk}-{n}"][()].transpose(2, 1, 0) for n in FIELDS if f"{blk}-{n}" in h5}
                shape = arr["r"].shape
                ks = range(shape[2]) if active[2] else [ngc if shape[2] > 2 * ngc else 0]
                blo = origin + ngc * dxyz  # the block interior
                bhi = origin + (np.array(shape) - ngc) * dxyz
                for i in range(shape[0]):
                    for j in range(shape[1]):
                        for k in ks:
                            c = origin + (np.array([i, j, k]) + 0.5) * dxyz
                            c[~active] = lo[~active]
                            q = {n: float(a[i, j, k]) for n, a in arr.items()}
                            if "psi" not in q:
                                q["psi"] = 0.0
                            out = [d for d in range(3) if active[d] and not (lo[d] < c[d] < hi[d])]
                            inside = all(blo[d] < c[d] < bhi[d] for d in range(3) if active[d])
                            if not out and not inside:
                                continue  # a ghost inside the domain: block interface or seam
                            if out:
                                cc = np.clip(c, lo, hi)
                                if not all(blo[d] < cc[d] < bhi[d] or cc[d] in (lo[d], hi[d])
                                           for d in range(3) if active[d]):
                                    beside += 1
                                    continue  # its donor along the normal is an interface ghost
                            if out and all(bc[2 * d] == "periodic" for d in out):
                                cw = c.copy()
                                for d in out:
                                    cw[d] = lo[d] + (c[d] - lo[d]) % (hi[d] - lo[d])
                                rho, _, p, _ = primitive(q, atm.gamma, atm.mhd, eglm)
                                re, pe = atm.state(cw)
                                note("periodic ghosts rho, p", max(abs(rho / re - 1), abs(p / pe - 1)))
                                continue
                            hyd = [d for d in out if bc[2 * d + (1 if c[d] > hi[d] else 0)] == "hydrostatic"]
                            if out and len(hyd) != len(out):
                                continue  # beyond a non-hydrostatic, non-periodic face (e.g. extrapolation on z)
                            if any(math.isnan(v) for v in q.values()):
                                nan += 1
                                continue
                            rho, u, p, b = primitive(q, atm.gamma, atm.mhd, eglm)
                            if not out:
                                re, pe = atm.state(c)
                                note("interior rho, p", max(abs(rho / re - 1), abs(p / pe - 1)))
                                note("interior u, B", max(float(np.max(np.abs(u - atm.vel))) / scale_u,
                                                          float(np.max(np.abs(b - atm.b))) / scale_b))
                                continue
                            if atm.polytropic and len(out) > 1:
                                skipped += 1
                                continue
                            if atm.polytropic:
                                d = out[0]
                                x1 = c.copy()
                                x1[d] = lo[d] + 0.5 * dxyz[d] if c[d] < lo[d] else hi[d] - 0.5 * dxyz[d]
                                r1, p1 = atm.state(x1)
                                re = r1 * math.exp(atm.g[d] * (c[d] - x1[d]) * r1 / p1)
                                pe = re * p1 / r1
                            else:
                                re, pe = atm.state(c)
                            note("hydrostatic ghosts rho, p", max(abs(rho / re - 1), abs(p / pe - 1)))
                            note("hydrostatic ghosts u, B", max(float(np.max(np.abs(u))) / scale_u,
                                                                float(np.max(np.abs(b - atm.b))) / scale_b))
    status = 0
    for key, value in err.items():
        if counts[key] == 0:
            continue
        ok = value <= args.tol
        status |= 0 if ok else 1
        print(f"   {key:>26}: {counts[key]:7d} cells, max relative error {value:.2e}  {'PASS' if ok else 'FAIL'}")
    if counts["hydrostatic ghosts rho, p"] == 0:
        print("   no hydrostatic ghost checked  FAIL")
        status = 1
    if nan:
        print(f"   {nan} checked cells hold NaN (a ghost read before its fill)  FAIL")
        status = 1
    if skipped:
        print(f"   {skipped} polytropic edge/corner ghosts not checked (the rule is per face)")
    if beside:
        print(f"   {beside} boundary ghosts beside their block not checked (their donor is an interface ghost)")
    print(f"{args.work.name}: {'PASS' if status == 0 else 'FAIL'} (tol {args.tol:g})")
    return status


def checkpoints(work: Path) -> dict[int, tuple[float, list[Path]]]:
    """Every checkpoint of a run: step -> (time, its per-rank files)."""
    basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
    out: dict[int, tuple[float, list[Path]]] = {}
    for f in sorted(work.glob(f"{basename}-*-proc*.h5")):
        step = int(f.name.split("-")[-2])
        if step not in out:
            xdmf = (work / f"{basename}-{step:09d}.xdmf").read_text()
            out[step] = (float(xdmf.split('Time Value="')[1].split('"')[0]), [])
        out[step][1].append(f)
    return out


def interior(ini: configparser.ConfigParser, paths: list[Path]):
    """Yield (centre, volume, conservative state) of every interior cell of one checkpoint."""
    for c, _, vol, q in interior_steps(ini, paths):
        yield c, vol, q


def interior_steps(ini: configparser.ConfigParser, paths: list[Path]):
    """Yield (centre, space steps, volume, conservative state) of every interior cell of one checkpoint."""
    ngc = int(ini["grid"]["ngc"])
    active = np.array([ini["grid"][f"null_{x}"].strip() != ".true." for x in "xyz"])
    for path in paths:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                arr = {n: h5[f"{blk}-{n}"][()].transpose(2, 1, 0) for n in FIELDS if f"{blk}-{n}" in h5}
                shape = arr["r"].shape
                ks = range(ngc, shape[2] - ngc) if active[2] else [ngc if shape[2] > 2 * ngc else 0]
                vol = float(np.prod(dxyz[active]))
                for i in range(ngc, shape[0] - ngc):
                    for j in range(ngc, shape[1] - ngc):
                        for k in ks:
                            c = origin + (np.array([i, j, k]) + 0.5) * dxyz
                            yield c, dxyz, vol, {n: float(a[i, j, k]) for n, a in arr.items()}


def freefall(args: argparse.Namespace) -> int:
    """GV-1: a uniform state in a periodic box falls freely, rho u = rho_1 g t, E = E_1 + rho_1 |g|^2 t^2 / 2, rho and p
    unchanged, at every cell of the last checkpoint."""
    ini = read_ini(args.work)
    atm = Atmosphere(ini)
    cps = checkpoints(args.work)
    t, paths = cps[max(cps)]
    e1 = atm.p1 / (atm.gamma - 1.0)
    err = {"rho": 0.0, "rho u": 0.0, "E": 0.0}
    for _, _, q in interior(ini, paths):
        err["rho"] = max(err["rho"], abs(q["r"] / atm.rho1 - 1.0))
        mom = np.array([q["ru"], q["rv"], q["rw"]])
        err["rho u"] = max(err["rho u"], float(np.max(np.abs(mom - atm.rho1 * atm.g * t))) / (atm.rho1 * np.linalg.norm(atm.g) * t))
        ee = e1 + 0.5 * atm.rho1 * float(np.dot(atm.g, atm.g)) * t * t
        err["E"] = max(err["E"], abs(q["rE"] / ee - 1.0))
    status = 0
    for key, value in err.items():
        ok = value <= args.tol
        status |= 0 if ok else 1
        print(f"   {key:>6} at t = {t:.6e}: max relative error {value:.2e}  {'PASS' if ok else 'FAIL'}")
    print(f"{args.work.name}: {'PASS' if status == 0 else 'FAIL'} (tol {args.tol:g})")
    return status


def energy(args: argparse.Namespace) -> int:
    """GV-3: in a closed box the mass and the total energy are conserved, at every checkpoint, while E itself exchanges
    with rho Phi (Phi = -g . x at the cell centres).

    The conserved total energy of the fourth-order energy source (M7-G2a) is `W = sum (E + rho Phi) V` plus, at each
    `wall-inviscid` face with gravity along its normal, `+- g_n dx_n rho V / 12` over the cells beside it (+ at the
    lower face, - at the upper): the source reads the face beyond the wall, whose mirrored mass flux leaves
    `g dx (F_1 - F_-1) / 24 = d/dt (-+ g dx rho_1 V / 12)` in the budget. W is asserted; the plain sum, which differs by
    that O(dx^2) wall term, is reported."""
    ini = read_ini(args.work)
    atm = Atmosphere(ini)
    lo = np.array([float(ini["grid"][f"emin_{x}"]) for x in "xyz"])
    hi = np.array([float(ini["grid"][f"emax_{x}"]) for x in "xyz"])
    walls = [d for d in range(3) if atm.g[d] != 0.0 and
             all(ini[f"bc_{'xyz'[d]}_{s}"]["type"].split(";")[0].strip() == "wall-inviscid" for s in ("min", "max"))]
    rows = []
    for step, (t, paths) in sorted(checkpoints(args.work).items()):
        m = e = w = b = 0.0
        for c, dxyz, vol, q in interior_steps(ini, paths):
            m += q["r"] * vol
            e += q["rE"] * vol
            w += -q["r"] * float(np.dot(atm.g, c)) * vol
            for d in walls:
                if c[d] - lo[d] < dxyz[d]:
                    b += atm.g[d] * dxyz[d] * q["r"] * vol / 12.0
                if hi[d] - c[d] < dxyz[d]:
                    b -= atm.g[d] * dxyz[d] * q["r"] * vol / 12.0
        rows.append((step, t, m, e, w, b))
    _, _, m0, e0, w0, b0 = rows[0]
    dm = max(abs(r[2] - m0) for r in rows) / abs(m0)
    dw = max(abs(r[3] + r[4] + r[5] - e0 - w0 - b0) for r in rows) / abs(e0 + w0 + b0)
    dt = max(abs(r[3] + r[4] - e0 - w0) for r in rows) / abs(e0 + w0)
    de = max(abs(r[3] - e0) for r in rows) / abs(e0 + w0)
    status = 0
    for key, value, checked in (("mass", dm, True), ("W", dw, not args.report_total),
                                ("E + rho Phi", dt, False), ("E alone", de, False)):
        ok = value <= args.tol
        if checked:
            status |= 0 if ok else 1
        tag = ("PASS" if ok else "FAIL") if checked else "(reported)"
        print(f"   {key:>11}: max relative drift over {len(rows)} checkpoints {value:.2e}  {tag}")
    print(f"{args.work.name}: {'PASS' if status == 0 else 'FAIL'} (tol {args.tol:g})")
    return status


def rest(args: argparse.Namespace) -> int:
    """GV-2 (report): the spurious velocity of an atmosphere at rest, max |u| / c over the interior at the last
    checkpoint of each run, and its order between consecutive runs (the runs given coarse to fine, each twice as fine)."""
    prev = None
    for work in args.works:
        ini = read_ini(work)
        atm = Atmosphere(ini)
        eglm = atm.mhd and ini["mhd"]["divergence_control"].split(";")[0].strip() == "eglm"
        cps = checkpoints(work)
        t, paths = cps[max(cps)]
        mach = 0.0
        for _, _, q in interior(ini, paths):
            rho, u, p, _ = primitive(q, atm.gamma, atm.mhd, eglm)
            mach = max(mach, float(np.linalg.norm(u)) / math.sqrt(atm.gamma * p / rho))
        order = f", order {math.log2(prev / mach):+.2f}" if prev and mach > 0.0 else ""
        print(f"   {work.name}: t = {t:.4f}, max |u| / c = {mach:.3e}{order}")
        prev = mach
    return 0


def agree(args: argparse.Namespace) -> int:
    """GV-7: two runs of the same input (CPU and FNL) at their last common checkpoint, per field max |a - b| over the
    field's scale (a vector component: its vector's largest component)."""
    ca, cb = checkpoints(args.a), checkpoints(args.b)
    step = max(set(ca) & set(cb))
    ini = read_ini(args.a)
    fa = {(tuple(np.round(c, 12))): q for c, _, q in interior(ini, ca[step][1])}
    fb = {(tuple(np.round(c, 12))): q for c, _, q in interior(ini, cb[step][1])}
    worst = 0.0
    group = {"ru": ("ru", "rv", "rw"), "rv": ("ru", "rv", "rw"), "rw": ("ru", "rv", "rw"),
             "bx": ("bx", "by", "bz"), "by": ("bx", "by", "bz"), "bz": ("bx", "by", "bz")}
    for name in FIELDS:
        if name not in next(iter(fa.values())):
            continue
        # a vector component is scaled by its vector's largest component: a component that is zero by symmetry holds
        # round-off only, and its own maximum would make the round-off difference of the backends look O(1)
        scale = max(abs(q[m]) for q in fa.values() for m in group.get(name, (name,)) if m in q) or 1.0
        diff = max(abs(fa[c][name] - fb[c][name]) for c in fa) / scale
        worst = max(worst, diff)
        print(f"   {name:>3}: max |a - b| / max |a| = {diff:.2e}")
    ok = worst <= args.tol
    print(f"{args.a.name} vs {args.b.name} at step {step}: {'PASS' if ok else 'FAIL'} (tol {args.tol:g})")
    return 0 if ok else 1


def wave(args: argparse.Namespace) -> int:
    """GV-9 (M7-G2a): the falling entropy wave of `wave2d`, the L1 errors of rho and p at the last checkpoint against
    the exact solution (the region-1 state, the wave of [initial_conditions] advected by `u_0 t + g t^2 / 2`, p
    uniform), and their orders between consecutive runs (given coarse to fine, each twice as fine). Exit 1 when the
    order of the finest pair is below `--min-order` for either."""
    prev = None
    orders = (0.0, 0.0)
    for work in args.works:
        ini = read_ini(work)
        atm = Atmosphere(ini)
        ic = ini["initial_conditions"]
        angle = math.radians(float(ic["wave_angle"]))
        amp, lam = float(ic["wave_amplitude"]), float(ic["wavelength"])
        n = np.array([math.cos(angle), math.sin(angle), 0.0])
        cps = checkpoints(work)
        t, paths = cps[max(cps)]
        shift = atm.vel * t + 0.5 * atm.g * t * t
        er = ep = vol = 0.0
        for c, v, q in interior(ini, paths):
            rho, _, p, _ = primitive(q, atm.gamma, atm.mhd, False)
            rho_x = atm.rho1 * (1.0 + amp * math.sin(2.0 * math.pi * float(np.dot(n, c - shift)) / lam))
            er += abs(rho - rho_x) * v
            ep += abs(p / atm.p1 - 1.0) * v
            vol += v
        er, ep = er / vol, ep / vol
        if prev:
            orders = (math.log2(prev[0] / er), math.log2(prev[1] / ep))
            tail = f", orders {orders[0]:+.2f} {orders[1]:+.2f}"
        else:
            tail = ""
        print(f"   {work.name}: t = {t:.4f}, L1 rho {er:.3e}, L1 p/p_1 {ep:.3e}{tail}")
        prev = (er, ep)
    ok = min(orders) >= args.min_order
    print(f"GV-9 wave2d: finest orders {orders[0]:+.2f} (rho), {orders[1]:+.2f} (p): {'PASS' if ok else 'FAIL'} "
          f"(min {args.min_order:g})")
    return 0 if ok else 1


def main() -> int:
    """Parse the command line and run the check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="check", required=True)
    for name, tol in (("profile", 1.0e-13), ("freefall", 1.0e-12), ("energy", 1.0e-12)):
        p = sub.add_parser(name)
        p.add_argument("work", type=Path)
        p.add_argument("--tol", type=float, default=tol)
        if name == "energy":
            p.add_argument("--report-total", action="store_true", help="report E + rho Phi, do not assert it")
    p = sub.add_parser("rest")
    p.add_argument("works", type=Path, nargs="+")
    p = sub.add_parser("wave")
    p.add_argument("works", type=Path, nargs="+")
    p.add_argument("--min-order", type=float, default=3.5)
    p = sub.add_parser("agree")
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("--tol", type=float, default=1.0e-10)
    args = parser.parse_args()
    checks = {"profile": profile, "freefall": freefall, "energy": energy, "rest": rest, "agree": agree, "wave": wave}
    return checks[args.check](args)


if __name__ == "__main__":
    raise SystemExit(main())
