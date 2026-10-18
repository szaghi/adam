#!/usr/bin/env python3
"""FLUME dual-energy oracle (issue #91, P2): euler-de against euler on the same input.

    de_oracle.py sod <euler work> <euler-de work>   the Euler fields (r, ru, rv, rw, rE) bitwise equal
    de_oracle.py noh <euler work> <euler-de work>   the cold Noh shock at t = 0.6 (shocks at 0.5 +- 0.12): the median
                                                    pressure within 0.08 of the centre (the wall-heating cells aside)
                                                    of both within TOL_P of the exact 1.2, the cold streams beyond 0.2
                                                    within TOL_COLD of p = 1e-6, the densities of the models within
                                                    TOL_RHO

    de_oracle.py agree <cpu euler> <fnl euler> <cpu euler-de> <fnl euler-de>
                                                    CPU and FNL agree (P3): per field, the euler-de difference relative
                                                    to the field's largest magnitude is at most AGREE_FACTOR times the
                                                    euler pair's (the device arithmetic on the same case), or at most
                                                    AGREE_FLOOR; rho s, which euler lacks, against the largest of them

The runs are 1-D along x (the other directions null): each cell is read once, on the line y = z = first cell.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

EULER = ("r", "ru", "rv", "rw", "rE")
TOL_P = 2.0e-2
TOL_COLD = 1.0e-3
TOL_RHO = 2.0e-2
AGREE_FACTOR = 2.0
AGREE_FLOOR = 1.0e-10


def last_line(work: Path) -> tuple[np.ndarray, dict[str, np.ndarray], float]:
    """Return the cell centres along x, the Euler fields along x and gamma of the last checkpoint of a run."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(work / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    phys = ini["physics"]
    gamma = float(phys["gamma"]) if "gamma" in phys else float(phys["cp"]) / float(phys["cv"])
    files = sorted(work.glob("*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    cells: dict[float, list[float]] = {}
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = [h5[f"{blk}-{v}"][()][ngc:-ngc, ngc:-ngc, ngc:-ngc][0, 0, :] for v in EULER]
                for i in range(q[0].size):
                    cells[round(float(origin[0] + (i + ngc + 0.5) * dxyz[0]), 10)] = [float(a[i]) for a in q]
    x = np.array(sorted(cells))
    fields = {v: np.array([cells[c][k] for c in x]) for k, v in enumerate(EULER)}
    return x, fields, gamma


def last_fields(work: Path, names: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Return every interior value of each named field of the last checkpoint, blocks keyed by their origin."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(work / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    files = sorted(work.glob("*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    blocks: dict[tuple[float, ...], dict[str, np.ndarray]] = {}
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                key = tuple(np.round(h5[f"{blk}-origin"][()], 10))
                blocks[key] = {v: h5[f"{blk}-{v}"][()][ngc:-ngc, ngc:-ngc, ngc:-ngc] for v in names if f"{blk}-{v}" in h5}
    keys = sorted(blocks)
    return {v: np.concatenate([blocks[k][v].ravel() for k in keys]) for v in names if v in blocks[keys[0]]}


def relative_gap(a: dict[str, np.ndarray], b: dict[str, np.ndarray]) -> dict[str, float]:
    """Return, per field, the largest difference over the largest magnitude."""
    return {v: float(np.abs(a[v] - b[v]).max() / max(np.abs(a[v]).max(), np.finfo(float).tiny)) for v in a}


def pressure(f: dict[str, np.ndarray], gamma: float) -> np.ndarray:
    """Return the pressure from the total energy."""
    return (gamma - 1.0) * (f["rE"] - 0.5 * (f["ru"] ** 2 + f["rv"] ** 2 + f["rw"] ** 2) / f["r"])


def main() -> int:
    """Run the requested check."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("case", choices=("sod", "noh", "agree"))
    ap.add_argument("works", type=Path, nargs="+", help="sod, noh: euler, euler-de; agree: the four runs")
    a = ap.parse_args()
    if a.case == "agree":
        if len(a.works) != 4:
            ap.error("agree takes four runs: cpu euler, fnl euler, cpu euler-de, fnl euler-de")
        base = relative_gap(last_fields(a.works[0], EULER), last_fields(a.works[1], EULER))
        de = relative_gap(last_fields(a.works[2], (*EULER, "rs")), last_fields(a.works[3], (*EULER, "rs")))
        ok = True
        for v, gap in de.items():
            ref = base.get(v, max(base.values()))
            bound = max(AGREE_FACTOR * ref, AGREE_FLOOR)
            ok &= gap <= bound
            print(f"   {v}: euler-de {gap:.3e}, euler {ref:.3e}{'' if v in base else ' (largest)'}, bound {bound:.3e}")
        print(f"DE agree: {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 1
    if len(a.works) != 2:
        ap.error(f"{a.case} takes two runs: euler, euler-de")
    a.euler, a.de = a.works
    x, fe, gamma = last_line(a.euler)
    xd, fd, _ = last_line(a.de)
    if x.size != xd.size or np.any(x != xd):
        print("   the two runs have different cells  FAIL")
        return 1
    ok = True
    if a.case == "sod":
        for v in EULER:
            same = np.array_equal(fe[v], fd[v])
            ok &= same
            print(f"   {v}: {'bitwise equal' if same else f'differs by {np.abs(fe[v] - fd[v]).max():.3e}'}")
    else:
        drho = float(np.abs(fe["r"] - fd["r"]).max())
        ok &= drho <= TOL_RHO
        print(f"   max |rho_de - rho_euler| = {drho:.3e} (bound {TOL_RHO:g})")
        for label, f in (("euler", fe), ("euler-de", fd)):
            p = pressure(f, gamma)
            post = np.abs(x - 0.5) < 0.08
            ahead = np.abs(x - 0.5) > 0.2
            ep = abs(float(np.median(p[post])) / 1.2 - 1.0)
            ec = float(np.abs(p[ahead] / 1.0e-6 - 1.0).max())
            ok &= ep <= TOL_P and ec <= TOL_COLD
            print(f"   {label}: post-shock p {np.median(p[post]):.5f} (exact 1.2, error {ep:.2e}, bound {TOL_P:g}); cold "
                  f"streams p error {ec:.2e} (bound {TOL_COLD:g})")
    print(f"DE {a.case}: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
