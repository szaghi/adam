#!/usr/bin/env python3
"""FLUME MHD plumbing oracle (issue #41, M2-P1): time step, state preservation, restart, auxiliary fields.

Why this oracle exists: in M2-P1 the MHD model has no face fluxes yet, so its residual is exactly zero. What P1 adds
is the plumbing around the numerics: the input predicate (state width, names), the MHD conversions and fast speed on
host and device, the kernel instances selected by the host, the outputs and the restart. Each check below fails on a
specific plumbing defect:

* --dt: the first step of an unperturbed uniform state must equal CFL / sum_d (|u_d| + c_{f,d}) / dx_d, with the
  fast magnetosonic speed computed here independently of the Fortran (wrong speed formula, wrong B index, wrong aux
  layout on the device, a kernel instance of the wrong model);
* --static: with a zero residual every volume integral is constant to the last bit over the whole run;
* --compare: two runs expected to hold the same state (a long run vs a one-step run, a continuous run vs a restarted
  one) have bitwise identical interior fields, for every saved variable, and the saved variables include the
  model's names (a missing B or psi field, a restart that loses psi);
* --aux: the saved auxiliary fields equal the values recomputed from the saved conservative fields (MHD pressure,
  enthalpy with the magnetic pressure, the B copies);
* --zero (used by the zero-field check, M2-P3): the named variables are exactly zero in the last checkpoint;
* --steady (M2-P3): the last checkpoint equals the first one (step 0), every saved variable, interior cells: bitwise
  by default, or within a tolerance (--steady-tol) relative to max(|variable|, 1).

Usage:
    mhd_plumbing_oracle.py --dt <work-dir> <ini>
    mhd_plumbing_oracle.py --static <work-dir>
    mhd_plumbing_oracle.py --compare <work-dir-a> <work-dir-b> --names r ru rv rw rE bx by bz [psi] [--ngc N]
    mhd_plumbing_oracle.py --aux <work-dir> <ini> [--tol T] [--ngc N]
    mhd_plumbing_oracle.py --zero <work-dir> --names bx by bz [psi] [--ngc N]
    mhd_plumbing_oracle.py --steady <work-dir> [--steady-tol T] [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np


def read_ini(path: Path) -> configparser.ConfigParser:
    """Return the parsed INI file (inline `;` comments stripped)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",))
    ini.read(path)
    return ini


def fast_speed(gamma: float, r: float, p: float, b: np.ndarray, d: int) -> float:
    """Return the fast magnetosonic speed along axis d (textbook form, independent of the Fortran rearrangement)."""
    a2 = gamma * p / r
    b2 = float(b @ b) / r
    bn2 = b[d] ** 2 / r
    return float(np.sqrt(0.5 * (a2 + b2 + np.sqrt((a2 + b2) ** 2 - 4.0 * a2 * bn2))))


def history(work: Path) -> list[list[str]]:
    """Return the data rows of the conservation history of one run, as strings."""
    hist = next(work.glob("*-conservation_history.dat"))
    return [line.split() for line in hist.open() if line.split() and line.split()[0][0] in "+-0123456789"]


def last_fields(work: Path, first: bool = False) -> dict[tuple[float, ...], dict[str, np.ndarray]]:
    """Return every saved variable of the last (or first) checkpoint of one run, keyed by block origin."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"mhd_plumbing_oracle: no checkpoint in {work}")
    steps = [int(p.name.split("-")[-2]) for p in files]
    last = min(steps) if first else max(steps)
    out: dict[tuple[float, ...], dict[str, np.ndarray]] = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            blocks = {k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}
            for blk in blocks:
                key = tuple(float(x) for x in np.round(h5[f"{blk}-origin"][()], 12))
                out[key] = {k[len(blk) + 1 :]: h5[k][()] for k in h5 if k.startswith(f"{blk}-") and h5[k].ndim == 3}
    return out


def interior(a: np.ndarray, ngc: int) -> np.ndarray:
    """Strip ngc ghost layers from a (k, j, i) block array."""
    return a[ngc:-ngc, ngc:-ngc, ngc:-ngc] if ngc > 0 else a


def check_dt(work: Path, ini_path: Path) -> bool:
    """Compare the first step with the dt of the uniform state."""
    ini = read_ini(ini_path)
    cp, cv = ini.getfloat("physics", "cp"), ini.getfloat("physics", "cv")
    gamma = cp / cv
    reg = ini["initial_conditions_region_1"]
    r, p = float(reg["r"]), float(reg["p"])
    u = np.array([float(reg[k]) for k in ("u", "v", "w")])
    b = np.array([float(reg[k]) for k in ("bx", "by", "bz")])
    g, amr = ini["grid"], ini["amr"]
    levels = 2 ** int(amr["iu_ref_levels"])
    dx = [(float(g[f"emax_{a}"]) - float(g[f"emin_{a}"])) / (int(g[f"n{c}"]) * levels)
          for a, c in zip("xyz", "ijk", strict=True)]
    lam = sum((abs(u[d]) + fast_speed(gamma, r, p, b, d)) / dx[d] for d in range(3))
    expected = ini.getfloat("time", "CFL") / lam
    rows = {int(row[0]): float(row[1]) for row in history(work)}
    if 1 not in rows:
        sys.exit(f"mhd_plumbing_oracle: {work} has no history row of step 1")
    rel = abs(rows[1] - expected) / expected
    ok = rel <= 1.0e-13
    print(f"dt of step 1: {rows[1]:.16e} vs fast-speed dt {expected:.16e}, relative difference {rel:.2e}  "
          f"{'PASS' if ok else 'FAIL'} (tol 1e-13)")
    return ok


def check_static(work: Path) -> bool:
    """Every row of the conservation history must hold the integrals of the first row, bit for bit (as printed)."""
    rows = history(work)
    changed = [row[0] for row in rows[1:] if row[2:] != rows[0][2:]]
    ok = len(rows) > 1 and not changed
    print(f"conservation history: {len(rows)} rows, integrals constant to the last printed digit  "
          f"{'PASS' if ok else 'FAIL (rows ' + ', '.join(changed[:5]) + ')'}")
    return ok


def check_compare(a_work: Path, b_work: Path, names: list[str], ngc: int) -> bool:
    """Bitwise comparison of every saved variable of two runs; the model's variables must be present."""
    a, b = last_fields(a_work), last_fields(b_work)
    if a.keys() != b.keys():
        print(f"compare: block origins differ ({len(a)} vs {len(b)} blocks)  FAIL")
        return False
    saved = set(next(iter(a.values())))
    missing = [n for n in names if n not in saved]
    worst, worst_name = 0.0, ""
    for key, fa in a.items():
        for name, arr in fa.items():
            diff = float(np.max(np.abs(interior(arr, ngc) - interior(b[key][name], ngc))))
            if diff > worst:
                worst, worst_name = diff, name
    ok = not missing and worst == 0.0
    print(f"{a_work.name} vs {b_work.name}: {len(a)} blocks, {len(saved)} variables, max |difference| {worst:.3e}"
          f"{' (' + worst_name + ')' if worst_name else ''}; model variables {' '.join(names)}: "
          f"{'present' if not missing else 'MISSING ' + ' '.join(missing)}  {'PASS (bitwise)' if ok else 'FAIL'}")
    return ok


def check_aux(work: Path, ini_path: Path, tol: float, ngc: int) -> bool:
    """The saved auxiliary fields must equal the values recomputed from the saved conservative fields."""
    ini = read_ini(ini_path)
    cp, cv = ini.getfloat("physics", "cp"), ini.getfloat("physics", "cv")
    gamma, rgas = cp / cv, cp - cv
    worst = 0.0
    for fields in last_fields(work).values():
        f = {k: interior(v, ngc) for k, v in fields.items()}
        r = f["r"]
        u, v, w = f["ru"] / r, f["rv"] / r, f["rw"] / r
        pb = 0.5 * (f["bx"] ** 2 + f["by"] ** 2 + f["bz"] ** 2)
        p = (gamma - 1.0) * (f["rE"] - 0.5 * r * (u * u + v * v + w * w) - pb)
        ref = {"rho": r, "u": u, "v": v, "w": w, "p": p, "T": p / (r * rgas), "H": (f["rE"] + p + pb) / r,
               "a": np.sqrt(gamma * p / r), "Bx": f["bx"], "By": f["by"], "Bz": f["bz"]}
        for name, val in ref.items():
            scale = max(float(np.max(np.abs(val))), 1.0e-300)
            worst = max(worst, float(np.max(np.abs(f[name] - val))) / scale)
    ok = worst <= tol
    print(f"auxiliary fields (rho, u, v, w, p, T, H, a, Bx, By, Bz) vs recomputed: max relative difference "
          f"{worst:.3e}  {'PASS' if ok else 'FAIL'} (tol {tol:.1e})")
    return ok


def check_zero(work: Path, names: list[str], ngc: int) -> bool:
    """The named variables must be exactly zero on the interior cells of the last checkpoint."""
    worst = 0.0
    for fields in last_fields(work).values():
        for name in names:
            worst = max(worst, float(np.max(np.abs(interior(fields[name], ngc)))))
    ok = worst == 0.0
    print(f"{work.name}: max |{', '.join(names)}| = {worst:.3e}  {'PASS (exactly zero)' if ok else 'FAIL'}")
    return ok


def check_steady(work: Path, ngc: int, tol: float) -> bool:
    """The last checkpoint must equal the first one (every saved variable, interior cells), bitwise if tol is 0."""
    a, b = last_fields(work, first=True), last_fields(work)
    if a.keys() != b.keys():
        print(f"{work.name}: first and last checkpoints hold different blocks  FAIL")
        return False
    worst, worst_name = 0.0, ""
    for key, fa in a.items():
        for name, arr in fa.items():
            scale = max(float(np.max(np.abs(arr))), 1.0) if tol > 0.0 else 1.0  # floor 1: psi starts at zero
            diff = float(np.max(np.abs(interior(arr, ngc) - interior(b[key][name], ngc)))) / scale
            if diff > worst:
                worst, worst_name = diff, name
    ok = worst <= tol
    kind = "bitwise" if tol == 0.0 else f"relative tol {tol:.1e}"
    print(f"{work.name}: last vs first checkpoint, {len(a)} blocks, max {'relative ' if tol > 0.0 else ''}|difference| "
          f"{worst:.3e}{' (' + worst_name + ')' if worst_name else ''}  {'PASS' if ok else 'FAIL'} ({kind})")
    return ok


def main() -> int:
    """Run the requested check."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dt", type=Path, nargs=2, metavar=("WORK", "INI"))
    parser.add_argument("--static", type=Path, metavar="WORK")
    parser.add_argument("--compare", type=Path, nargs=2, metavar=("WORK_A", "WORK_B"))
    parser.add_argument("--names", nargs="+", default=[])
    parser.add_argument("--aux", type=Path, nargs=2, metavar=("WORK", "INI"))
    parser.add_argument("--zero", type=Path, metavar="WORK")
    parser.add_argument("--steady", type=Path, metavar="WORK")
    parser.add_argument("--steady-tol", type=float, default=0.0, help="--steady relative tolerance (0: bitwise)")
    parser.add_argument("--tol", type=float, default=1.0e-12)
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    if args.dt is not None:
        ok = check_dt(*args.dt)
    elif args.static is not None:
        ok = check_static(args.static)
    elif args.compare is not None:
        ok = check_compare(*args.compare, names=args.names, ngc=args.ngc)
    elif args.aux is not None:
        ok = check_aux(*args.aux, tol=args.tol, ngc=args.ngc)
    elif args.zero is not None:
        ok = check_zero(args.zero, names=args.names, ngc=args.ngc)
    elif args.steady is not None:
        ok = check_steady(args.steady, ngc=args.ngc, tol=args.steady_tol)
    else:
        parser.error("one of --dt, --static, --compare, --aux, --zero, --steady is required")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
