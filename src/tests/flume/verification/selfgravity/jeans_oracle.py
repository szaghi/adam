#!/usr/bin/env python3
"""Oracle of the FLUME self-gravity Jeans cases (issue #81, M8-P3, SV-3 and SV-4), written by make_jeans.py.

    jeans_oracle.py wave <work> [<work> ...] [--min-order P]
    jeans_oracle.py growth <work> [--t0 T0] [--t1 T1] [--tol T]
    jeans_oracle.py symmetry <work> [--tol T]
    jeans_oracle.py conservation <work> [<work> ...] [--tol-momentum M] [--min-order P] [--label NAME]
    jeans_oracle.py compare <work_a> <work_b> [--tol T]
    jeans_oracle.py regrid <work> --frequency N [--tol-step T] [--tol-momentum M] [--tol-mass M]

The exact linear solution of a Jeans mode of [initial_conditions] (`sine-wave`, `wave_mode = jeans`) on the uniform
state of region 1 under the self-gravity of [gravity] (`G` or `G_star`): with k = 2 pi / wavelength along
n = (cos a, sin a), the phase theta = k n . x and w^2 = c^2 k^2 - 4 pi G rho_0 (c^2 = gamma p_0 / rho_0),

    stable (w^2 >= 0)    rho = rho_0 (1 + A sin(theta - w t)),
    unstable (w^2 < 0)   rho = rho_0 (1 + A e^(sigma t) sin(theta)),  sigma = sqrt(-w^2).

- `wave` (SV-3): the L1 error of rho at the last checkpoint against it, over rho_0 times the mode amplitude at that
  time, and the orders between consecutive runs (coarse to fine, each twice as fine); exit 1 when the order of the
  finest pair is below `--min-order` (default 1.8).
- `growth` (SV-4): the mode amplitude a(t) = 2 <(rho / rho_0 - 1) sin(theta)> of every checkpoint and sigma from a
  least-squares fit of log a over [t0, t1], against the exact one; exit 1 when the relative error exceeds `--tol`.
- `symmetry` (SV-4): the point reflection about the box corner (make_jeans.py places the density's centre of symmetry
  there): rho even, velocity odd. Per checkpoint the defects max|rho(x) - rho(x')| / (rho_0 a(t)) and
  max|u(x) + u(x')| / max|u|; exit 1 when the last ones exceed `--tol`.
- `conservation` (SV-4): from the conservation history, the drift of E + W (`E_plus_W`, max |(E + W)(t) - (E + W)(0)|
  over |(E + W)(0)|) and the largest |momentum integral|, per run; with several runs of one case at decreasing CFL
  (each half the previous), the order of the E + W drift in dt. Exit 1 when a momentum exceeds `--tol-momentum` or
  the order of the finest pair is below `--min-order`, or (with `--tol-mass`) a relative mass drift exceeds it.
- `regrid` (SV-11, issue #81 P5c): a run under runtime AMR. The grid must change (the blocks of the first and the last
  checkpoint differ); the E + W history is split into jumps (a step change above 100 times the median) and the rest:
  every jump must fall on a regrid step (a multiple of `--frequency`), the largest other step change must stay below
  `--tol-step` (relative to E + W), the momentum below `--tol-momentum` and the mass drift below `--tol-mass`. The jumps
  (the discrete W changes with the grid) are reported, not bounded.
- `compare` (SV-10, issue #81 P5b): two runs of one case on one grid (the FFT and the multigrid solver), their density
  at the last checkpoint, max |rho_a - rho_b| over rho_0 times the mode amplitude; exit 1 above `--tol`.
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gravity"))
from atmosphere_oracle import checkpoints, interior, read_ini  # noqa: E402


def value(sec: configparser.SectionProxy, key: str) -> float:
    """A numeric option, inline comments stripped."""
    return float(sec[key].split(";")[0])


class Jeans:
    """The exact linear Jeans mode of an input."""

    def __init__(self, ini: configparser.ConfigParser) -> None:
        reg, ic, grav = ini["initial_conditions_region_1"], ini["initial_conditions"], ini["gravity"]
        self.rho0, self.p0 = value(reg, "r"), value(reg, "p")
        self.gamma = value(ini["physics"], "gamma")
        big_g = value(grav, "G") if "G" in grav else value(grav, "G_star")
        angle = math.radians(value(ic, "wave_angle"))
        self.amp = value(ic, "wave_amplitude")
        self.k = 2.0 * math.pi / value(ic, "wavelength")
        self.n = np.array([math.cos(angle), math.sin(angle), 0.0])
        self.w2 = self.gamma * self.p0 / self.rho0 * self.k**2 - 4.0 * math.pi * big_g * self.rho0
        self.corner = np.array([value(ini["grid"], f"emin_{x}") for x in "xyz"])
        self.length = np.array([value(ini["grid"], f"emax_{x}") for x in "xyz"]) - self.corner
        self.null = [ini["grid"][f"null_{x}"].strip() == ".true." for x in "xyz"]

    @property
    def stable(self) -> bool:
        """Whether the mode oscillates."""
        return self.w2 >= 0.0

    @property
    def rate(self) -> float:
        """w (stable) or sigma (unstable)."""
        return math.sqrt(abs(self.w2))

    def theta(self, x: np.ndarray) -> float:
        """The phase at x."""
        return self.k * float(np.dot(self.n, x))

    def amplitude(self, t: float) -> float:
        """The relative mode amplitude at t."""
        return self.amp if self.stable else self.amp * math.exp(self.rate * t)

    def rho(self, x: np.ndarray, t: float) -> float:
        """The exact density."""
        if self.stable:
            return self.rho0 * (1.0 + self.amp * math.sin(self.theta(x) - self.rate * t))
        return self.rho0 * (1.0 + self.amplitude(t) * math.sin(self.theta(x)))


def wave_error(work: Path) -> tuple[float, float]:
    """The time and the L1 error of rho over rho_0 times the mode amplitude at the last checkpoint."""
    ini = read_ini(work)
    jm = Jeans(ini)
    cps = checkpoints(work)
    t, paths = cps[max(cps)]
    err = vol = 0.0
    for c, v, q in interior(ini, paths):
        err += abs(q["r"] - jm.rho(c, t)) * v
        vol += v
    return t, err / vol / (jm.rho0 * jm.amplitude(t))


def wave(args: argparse.Namespace) -> int:
    """SV-3: the L1 errors against the exact mode and their orders."""
    prev = None
    order = 0.0
    for work in args.works:
        t, e = wave_error(work)
        tail = ""
        if prev is not None:
            order = math.log2(prev / e)
            tail = f", order {order:+.2f}"
        print(f"   {work.name}: t = {t:.4f}, L1 rho / (rho_0 A) {e:.3e}{tail}")
        prev = e
    ok = order >= args.min_order
    print(f"Jeans wave: finest order {order:+.2f}: {'PASS' if ok else 'FAIL'} (min {args.min_order:g})")
    return 0 if ok else 1


def mode_amplitudes(work: Path) -> list[tuple[float, float]]:
    """Per checkpoint (t, a), a = 2 <(rho / rho_0 - 1) sin(theta)>."""
    ini = read_ini(work)
    jm = Jeans(ini)
    out = []
    for _, (t, paths) in sorted(checkpoints(work).items()):
        s = vol = 0.0
        for c, v, q in interior(ini, paths):
            s += (q["r"] / jm.rho0 - 1.0) * math.sin(jm.theta(c)) * v
            vol += v
        out.append((t, 2.0 * s / vol))
    return out


def growth(args: argparse.Namespace) -> int:
    """SV-4: the growth rate of the unstable mode against the exact one."""
    jm = Jeans(read_ini(args.work))
    rows = [(t, a) for t, a in mode_amplitudes(args.work) if args.t0 <= t <= args.t1 and a > 0.0]
    ts = np.array([r[0] for r in rows])
    la = np.log(np.array([r[1] for r in rows]))
    sigma = float(np.polyfit(ts, la, 1)[0])
    err = abs(sigma / jm.rate - 1.0)
    ok = (not jm.stable) and err <= args.tol
    print(f"   {args.work.name}: {len(rows)} checkpoints in [{args.t0:g}, {args.t1:g}], sigma {sigma:.6f}, exact "
          f"{jm.rate:.6f} (pi sqrt(10) = {math.pi * math.sqrt(10.0):.6f}), relative error {err:.2e}")
    print(f"Jeans growth: {'PASS' if ok else 'FAIL'} (tol {args.tol:g})")
    return 0 if ok else 1


def mirror_defects(work: Path) -> list[tuple[float, float, float]]:
    """Per checkpoint (t, density defect, velocity defect) of the point reflection about the box corner."""
    ini = read_ini(work)
    jm = Jeans(ini)
    amps = dict(mode_amplitudes(work))
    out = []
    for _, (t, paths) in sorted(checkpoints(work).items()):
        cells = {}
        for c, _, q in interior(ini, paths):
            key = tuple(round(float(x), 10) for x in (c - jm.corner) / jm.length)
            cells[key] = (q["r"], np.array([q["ru"], q["rv"], q["rw"]]) / q["r"])
        dr = du = umax = 0.0
        for key, (r, u) in cells.items():
            mirror = tuple(x if null else round((1.0 - x) % 1.0, 10) for x, null in zip(key, jm.null, strict=True))
            r2, u2 = cells[mirror]
            dr = max(dr, abs(r - r2))
            du = max(du, float(np.max(np.abs(u + u2))))
            umax = max(umax, float(np.max(np.abs(u))))
        out.append((t, dr / (jm.rho0 * abs(amps[t])), du / umax if umax > 0.0 else 0.0))
    return out


def symmetry(args: argparse.Namespace) -> int:
    """SV-4: the mirror symmetry of the unstable mode."""
    rows = mirror_defects(args.work)
    t, dr, du = rows[-1]
    worst = max(max(r[1], r[2]) for r in rows)
    ok = max(dr, du) <= args.tol
    print(f"   {args.work.name}: {len(rows)} checkpoints; at t = {t:.4f} density defect {dr:.2e}, velocity defect "
          f"{du:.2e}; largest over the run {worst:.2e}")
    print(f"Jeans symmetry: {'PASS' if ok else 'FAIL'} (tol {args.tol:g})")
    return 0 if ok else 1


def history(work: Path) -> dict[str, np.ndarray]:
    """The conservation history, column name -> values."""
    path = next(work.glob("*-conservation_history.dat"))
    lines = path.read_text().splitlines()
    names = [n.strip('"') for n in lines[0].split("=", 1)[1].split()]
    data = np.array([[float(v) for v in row.split()] for row in lines[1:] if row.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


def conservation(args: argparse.Namespace) -> int:
    """SV-4: the drift of E + W and the momentum, and the order of the drift in dt."""
    drifts = []
    ok = True
    for work in args.works:
        h = history(work)
        ew = h["E_plus_W"]
        drift = float(np.max(np.abs(ew - ew[0])) / abs(ew[0]))
        mom = float(max(np.max(np.abs(h[k])) for k in ("int_ru", "int_rv", "int_rw")))
        w_range = float(np.max(np.abs(h["W"] - h["W"][0])) / abs(ew[0]))
        e_range = float(np.max(np.abs(h["int_rE"] - h["int_rE"][0])) / abs(ew[0]))
        mass = float(np.max(np.abs(h["int_r"] - h["int_r"][0])) / abs(h["int_r"][0]))
        tail = ""
        if drifts:
            tail = f", order {math.log2(drifts[-1] / drift):+.2f}"
        print(f"   {work.name}: {len(ew)} rows, E + W drift {drift:.2e} (E alone {e_range:.2e}, W alone "
              f"{w_range:.2e}), max |momentum| {mom:.2e}, mass drift {mass:.1e}{tail}")
        ok = ok and mom <= args.tol_momentum and (args.tol_mass is None or mass <= args.tol_mass)
        drifts.append(drift)
    if len(drifts) > 1:
        order = math.log2(drifts[-2] / drifts[-1])
        ok = ok and order >= args.min_order
    print(f"{args.label} conservation: {'PASS' if ok else 'FAIL'} (momentum tol {args.tol_momentum:g}, min order "
          f"{args.min_order:g})")
    return 0 if ok else 1


def compare(args: argparse.Namespace) -> int:
    """SV-10: two runs of one case on one grid agree at the last checkpoint."""
    dens = []
    for work in (args.work_a, args.work_b):
        ini = read_ini(work)
        cps = checkpoints(work)
        t, paths = cps[max(cps)]
        dens.append((t, {tuple(round(float(x), 10) for x in c): q["r"] for c, _, q in interior(ini, paths)}))
    jm = Jeans(read_ini(args.work_a))
    (ta, ra), (tb, rb) = dens
    same = ta == tb and ra.keys() == rb.keys()
    diff = max(abs(ra[k] - rb[k]) for k in ra) / (jm.rho0 * jm.amplitude(ta)) if same else math.inf
    ok = diff <= args.tol
    print(f"   {args.work_a.name} against {args.work_b.name}: t = {ta:.4f}, {len(ra)} cells, max |rho_a - rho_b| / "
          f"(rho_0 A) {diff:.2e}")
    print(f"Jeans compare: {'PASS' if ok else 'FAIL'} (tol {args.tol:g})")
    return 0 if ok else 1


def checkpoint_blocks(paths: list[Path]) -> int:
    """The blocks of one checkpoint over its per-rank files."""
    import h5py  # noqa: PLC0415

    n = 0
    for path in paths:
        with h5py.File(path, "r") as h5:
            n += len({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")})
    return n


def regrid(args: argparse.Namespace) -> int:
    """SV-11: E + W, the mass and the momentum of a run under runtime AMR."""
    cps = checkpoints(args.work)
    first, last = checkpoint_blocks(cps[min(cps)][1]), checkpoint_blocks(cps[max(cps)][1])
    h = history(args.work)
    ew = h["E_plus_W"]
    d = np.diff(ew) / abs(ew[0])
    big = np.abs(d) > 100.0 * np.median(np.abs(d))
    its = h["it"][1:].astype(int)
    jumps = [(int(its[i]), float(d[i])) for i in np.nonzero(big)[0]]
    on_regrid = all(it % args.frequency == 0 for it, _ in jumps)
    other = float(np.max(np.abs(d[~big]))) if np.any(~big) else 0.0
    mom = float(max(np.max(np.abs(h[k])) for k in ("int_ru", "int_rv", "int_rw")))
    mass = float(np.max(np.abs(h["int_r"] - h["int_r"][0])) / abs(h["int_r"][0]))
    ok = first != last and on_regrid and other <= args.tol_step and mom <= args.tol_momentum and mass <= args.tol_mass
    print(f"   {args.work.name}: blocks {first} -> {last}, {len(ew)} rows; E + W jumps at regrids "
          f"{', '.join(f'it {it}: {x:+.2e}' for it, x in jumps) or 'none'}; largest other step change {other:.2e}; "
          f"max |momentum| {mom:.2e}; mass drift {mass:.1e}")
    print(f"Jeans regrid: {'PASS' if ok else 'FAIL'} (grid changed {first != last}, jumps on regrid steps {on_regrid}, "
          f"step tol {args.tol_step:g}, momentum tol {args.tol_momentum:g}, mass tol {args.tol_mass:g})")
    return 0 if ok else 1


def main() -> int:
    """Parse the command line and run the check."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="check", required=True)
    p = sub.add_parser("wave")
    p.add_argument("works", type=Path, nargs="+")
    p.add_argument("--min-order", type=float, default=1.8)
    p = sub.add_parser("growth")
    p.add_argument("work", type=Path)
    p.add_argument("--t0", type=float, default=0.1)
    p.add_argument("--t1", type=float, default=0.5)
    p.add_argument("--tol", type=float, default=1.0e-2)
    p = sub.add_parser("symmetry")
    p.add_argument("work", type=Path)
    p.add_argument("--tol", type=float, default=1.0e-10)
    p = sub.add_parser("conservation")
    p.add_argument("works", type=Path, nargs="+")
    p.add_argument("--tol-momentum", type=float, default=1.0e-14)
    p.add_argument("--min-order", type=float, default=2.5)
    p.add_argument("--label", default="Jeans", help="name of the case in the verdict line")
    p.add_argument("--tol-mass", type=float, default=None, help="bound on the relative mass drift (default: none)")
    p = sub.add_parser("regrid")
    p.add_argument("work", type=Path)
    p.add_argument("--frequency", type=int, required=True)
    p.add_argument("--tol-step", type=float, default=1.0e-6)
    p.add_argument("--tol-momentum", type=float, default=1.0e-14)
    p.add_argument("--tol-mass", type=float, default=1.0e-14)
    p = sub.add_parser("compare")
    p.add_argument("work_a", type=Path)
    p.add_argument("work_b", type=Path)
    p.add_argument("--tol", type=float, default=1.0e-6)
    args = parser.parse_args()
    checks = {"wave": wave, "growth": growth, "symmetry": symmetry, "conservation": conservation, "compare": compare,
              "regrid": regrid}
    return checks[args.check](args)


if __name__ == "__main__":
    raise SystemExit(main())
