#!/usr/bin/env python3
"""Oracle of the FLUME stratified atmospheres (issue #73, M7): the hydrostatic state and its ghost cells.

    atmosphere_oracle.py profile <work> [--tol T]

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


def main() -> int:
    """Parse the command line and run the check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="check", required=True)
    p = sub.add_parser("profile")
    p.add_argument("work", type=Path)
    p.add_argument("--tol", type=float, default=1.0e-13)
    args = parser.parse_args()
    return profile(args)


if __name__ == "__main__":
    raise SystemExit(main())
