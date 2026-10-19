#!/usr/bin/env python3
"""FLUME sink-particle oracle (issue #94): the replicated sinks and their mesh gravity.

    si_oracle.py rest <work>            a sink at rest in a uniform periodic background: it stays put (displacement
                                        within TOL_REST cells), and gas plus sink momentum is conserved (within
                                        TOL_MOMENTUM of the momentum they exchange)
    si_oracle.py binary <work>          the equal-mass circular binary: the separation within TOL_SEPARATION, the period
                                        (angle of the separation vector, fitted) against Kepler's at the mean
                                        separation within TOL_PERIOD, the orbital energy K - G m1 m2 / d within
                                        TOL_ENERGY, the angular momentum within TOL_ANGULAR, the centre of mass within
                                        TOL_COM
    si_oracle.py same <work a> <work b> two runs bitwise equal: the conservation and sinks histories (the restart leg)
    si_oracle.py accrete <work>         (P2) accretion: gas plus sink mass conserved within TOL_MASS, the sink grown,
                                        gas plus sink momentum within TOL_ACCRETE_MOMENTUM of the momentum exchanged
    si_oracle.py collapse <work>        (P2) the collapsing polytrope formed exactly one sink, within a cell of the
                                        centre, and gas plus sink mass is conserved within TOL_MASS
    si_oracle.py none <work>            (P2) no sink formed (the Jeans-stable box above the threshold)
    si_oracle.py merge <work>           (P2) the two sinks merged: one left, with their total mass and momentum
    si_oracle.py agree <cpu> <fnl>      (P3) CPU and FNL agree: the same steps and sinks, every sink's mass, position
                                        (in cells) and velocity (over the largest speed) within TOL_AGREE;
                                        --mass-only: the masses and the final distance from the centre

The histories are those of the run: `<base>-conservation_history.dat` (the gas integrals, W, E + W, and the sinks'
mass, momentum and kinetic energy) and `<base>-sinks_history.dat` (one row per sink per step: it time id m x y z u v w).
The binary's sinks feel the mesh gravity of the deposit (softened at about two cells), so at a separation of 8 cells its
period and separation are those of Kepler within the mesh-force error, not to round-off (at 8 cells the effective force
is about 0.7 % above Newton's). E + W + K_sinks is no conserved quantity with point masses on a mesh: W includes each
sink's CIC self-energy, about G m^2 / h, which changes as the sink crosses cells (17 % over the orbit); the orbital
energy of the pair is the one checked.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import numpy as np

TOL_REST = 1.0e-3
TOL_MOMENTUM = 1.0e-3
TOL_SEPARATION = 3.0e-2
TOL_ENERGY = 1.0e-2
TOL_PERIOD = 2.0e-2
TOL_ANGULAR = 1.0e-2
TOL_COM = 1.0e-10
TOL_MASS = 1.0e-13
TOL_ACCRETE_MOMENTUM = 5.0e-2
TOL_AGREE = 1.0e-8


def read_history(path: Path) -> dict[str, np.ndarray]:
    """Return the columns of a history file by name."""
    names = path.read_text().splitlines()[0].split("=", 1)[1].replace('"', "").split()
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    return {n: data[:, k] for k, n in enumerate(names)}


def sinks_by_id(work: Path) -> dict[int, dict[str, np.ndarray]]:
    """Return the sinks history split per sink."""
    col = read_history(next(work.glob("*-sinks_history.dat")))
    return {int(i): {k: v[col["id"] == i] for k, v in col.items()} for i in np.unique(col["id"])}


def cell_size(work: Path) -> float:
    """Return the cell spacing of the uniform box."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(work / "input.ini")
    g = ini["grid"]
    blocks = 2 ** int(ini["amr"]["iu_ref_levels"])
    return (float(g["emax_x"]) - float(g["emin_x"])) / (int(g["ni"]) * blocks)


def rest(work: Path) -> bool:
    """Check the sink at rest in the periodic background."""
    s = sinks_by_id(work)[1]
    h = cell_size(work)
    disp = float(np.max(np.sqrt((s["x"] - s["x"][0]) ** 2 + (s["y"] - s["y"][0]) ** 2 + (s["z"] - s["z"][0]) ** 2)))
    ok = True
    status = "PASS" if disp <= TOL_REST * h else "FAIL"
    ok &= status == "PASS"
    print(f"   sink displacement {disp / h:.2e} cells over t = {s['time'][-1]:.3f} (bound {TOL_REST:.0e})  {status}")
    c = read_history(next(work.glob("*-conservation_history.dat")))
    gas = np.stack([c["int_ru"], c["int_rv"], c["int_rw"]], axis=1)
    total = gas + np.stack([c["Px_sinks"], c["Py_sinks"], c["Pz_sinks"]], axis=1)
    exchanged = float(np.max(np.abs(gas)))
    mismatch = float(np.max(np.abs(total)))
    status = "PASS" if mismatch <= TOL_MOMENTUM * max(exchanged, 1.0e-300) else "FAIL"
    ok &= status == "PASS"
    print(f"   momentum: exchanged {exchanged:.2e}, gas + sink {mismatch:.2e} ({mismatch / exchanged:.1e} of it, bound "
          f"{TOL_MOMENTUM:.0e})  {status}")
    e = c["E_plus_W"] + c["K_sinks"]
    print(f"   E + W + K_sinks: drift {abs(e[-1] - e[0]) / abs(e[0]):.2e} (reported)")
    return ok


def binary(work: Path) -> bool:
    """Check the circular binary."""
    sinks = sinks_by_id(work)
    a, b = sinks[1], sinks[2]
    dx, dy, dz = b["x"] - a["x"], b["y"] - a["y"], b["z"] - a["z"]
    d = np.sqrt(dx**2 + dy**2 + dz**2)
    time = a["time"]
    ok = True
    sep = float(np.max(np.abs(d / d[0] - 1.0)))
    status = "PASS" if sep <= TOL_SEPARATION else "FAIL"
    ok &= status == "PASS"
    print(f"   separation: initial {d[0]:.4f}, largest deviation {sep:.2e} (bound {TOL_SEPARATION:.0e})  {status}")
    angle = np.unwrap(np.arctan2(dy, dx))
    omega = float(np.polyfit(time, angle, 1)[0])
    period = 2.0 * np.pi / abs(omega)
    mtot = a["m"][0] + b["m"][0]
    kepler = 2.0 * np.pi * np.sqrt(float(np.mean(d)) ** 3 / mtot)  # G = 1, at the mean separation
    err = abs(period / kepler - 1.0)
    status = "PASS" if err <= TOL_PERIOD else "FAIL"
    ok &= status == "PASS"
    print(f"   period {period:.5f} over {abs(angle[-1] - angle[0]) / (2 * np.pi):.2f} orbits, Kepler at the mean "
          f"separation {kepler:.5f}: "
          f"{err:.2e} (bound {TOL_PERIOD:.0e})  {status}")
    lz = a["m"] * (a["x"] * a["v"] - a["y"] * a["u"]) + b["m"] * (b["x"] * b["v"] - b["y"] * b["u"])
    dl = float(np.max(np.abs(lz / lz[0] - 1.0)))
    status = "PASS" if dl <= TOL_ANGULAR else "FAIL"
    ok &= status == "PASS"
    print(f"   angular momentum: largest deviation {dl:.2e} (bound {TOL_ANGULAR:.0e})  {status}")
    com = np.stack([a["m"] * a[k] + b["m"] * b[k] for k in ("x", "y", "z")], axis=1) / mtot
    dcom = float(np.max(np.abs(com - com[0])))
    status = "PASS" if dcom <= TOL_COM else "FAIL"
    ok &= status == "PASS"
    print(f"   centre of mass: largest displacement {dcom:.2e} (bound {TOL_COM:.0e})  {status}")
    kin = sum(0.5 * s["m"] * (s["u"] ** 2 + s["v"] ** 2 + s["w"] ** 2) for s in (a, b))
    e_orb = kin - a["m"] * b["m"] / d  # G = 1
    de = float(np.max(np.abs(e_orb / e_orb[0] - 1.0)))
    status = "PASS" if de <= TOL_ENERGY else "FAIL"
    ok &= status == "PASS"
    print(f"   orbital energy K - G m1 m2 / d: largest deviation {de:.2e} (bound {TOL_ENERGY:.0e})  {status}")
    return ok


def same(work_a: Path, work_b: Path) -> bool:
    """Check that two runs' conservation and sinks histories are bitwise equal on their common rows."""
    ok = True
    for pattern in ("*-conservation_history.dat", "*-sinks_history.dat"):
        a = read_history(next(work_a.glob(pattern)))
        b = read_history(next(work_b.glob(pattern)))
        rows = min(a["it"].size, b["it"].size)
        equal = set(a) == set(b) and all((a[k][-rows:] == b[k][-rows:]).all() for k in a)
        ok &= equal
        print(f"   {pattern[2:]}: the last {rows} rows {'bitwise equal' if equal else 'DIFFERENT'}")
    print(f"   {'PASS' if ok else 'FAIL'}")
    return ok


def total_mass(work: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return gas plus sink mass per row, and the conservation history."""
    c = read_history(next(work.glob("*-conservation_history.dat")))
    return c["int_r"] + c["M_sinks"], c


def accrete(work: Path) -> bool:
    """Check threshold accretion."""
    mass, c = total_mass(work)
    ok = True
    drift = float(np.max(np.abs(mass / mass[0] - 1.0)))
    status = "PASS" if drift <= TOL_MASS else "FAIL"
    ok &= status == "PASS"
    print(f"   gas + sink mass: drift {drift:.2e} (bound {TOL_MASS:.0e})  {status}")
    grown = c["M_sinks"][-1] - c["M_sinks"][0]
    status = "PASS" if grown > 0.0 else "FAIL"
    ok &= status == "PASS"
    share = grown / c["M_sinks"][0] * 100
    print(f"   accreted {grown:.4e} over {c['it'].size - 1} steps ({share:.1f} % of the sink)  {status}")
    gas = np.stack([c["int_ru"], c["int_rv"], c["int_rw"]], axis=1)
    total = gas + np.stack([c["Px_sinks"], c["Py_sinks"], c["Pz_sinks"]], axis=1)
    exchanged, mismatch = float(np.max(np.abs(gas))), float(np.max(np.abs(total)))
    status = "PASS" if mismatch <= TOL_ACCRETE_MOMENTUM * exchanged else "FAIL"
    ok &= status == "PASS"
    print(f"   momentum: exchanged {exchanged:.2e}, gas + sink {mismatch:.2e} ({mismatch / exchanged:.1e} of it, bound "
          f"{TOL_ACCRETE_MOMENTUM:.0e}: the gravity of the gas, face gradients, against the sink's, CIC)  {status}")
    return ok


def final_sinks(work: Path) -> dict[int, dict[str, float]]:
    """Return the sinks of the last history row."""
    col = read_history(next(work.glob("*-sinks_history.dat")))
    if col["it"].size == 0:
        return {}
    last = col["it"] == col["it"][-1]
    return {int(i): {k: float(v[last][n]) for k, v in col.items()} for n, i in enumerate(col["id"][last])}


def collapse(work: Path) -> bool:
    """Check that the collapsing polytrope formed exactly one central sink."""
    path = next(work.glob("*-sinks_history.dat"))
    rows = len(path.read_text().splitlines()) - 1
    sinks = final_sinks(work) if rows > 0 else {}
    h = cell_size(work)
    ok = len(sinks) == 1
    print(f"   sinks at the end: {len(sinks)} (expected 1)  {'PASS' if ok else 'FAIL'}")
    if ok:
        s = sinks[1]
        r = float(np.sqrt(s["x"] ** 2 + s["y"] ** 2 + s["z"] ** 2))
        status = "PASS" if r <= h else "FAIL"
        ok &= status == "PASS"
        print(f"   the sink at {r / h:.2e} cells from the centre (bound 1), mass {s['m']:.4f}  {status}")
    mass, _ = total_mass(work)
    drift = float(np.max(np.abs(mass / mass[0] - 1.0)))
    status = "PASS" if drift <= TOL_MASS else "FAIL"
    ok &= status == "PASS"
    print(f"   gas + sink mass: drift {drift:.2e} (bound {TOL_MASS:.0e})  {status}")
    return ok


def none(work: Path) -> bool:
    """Check that no sink formed."""
    path = next(work.glob("*-sinks_history.dat"))
    rows = len(path.read_text().splitlines()) - 1
    c = read_history(next(work.glob("*-conservation_history.dat")))
    ok = rows == 0 and float(np.max(np.abs(c["M_sinks"]))) == 0.0
    print(f"   sinks formed: {rows} history rows, largest sink mass {np.max(np.abs(c['M_sinks'])):.1e} over "
          f"{c['it'].size - 1} steps (expected none)  {'PASS' if ok else 'FAIL'}")
    return ok


def merge(work: Path) -> bool:
    """Check that the two sinks merged with their mass and momentum."""
    col = read_history(next(work.glob("*-sinks_history.dat")))
    first = col["it"] == col["it"][0]
    m0 = float(np.sum(col["m"][first]))
    p0 = np.array([float(np.sum(col["m"][first] * col[k][first])) for k in ("u", "v", "w")])
    sinks = final_sinks(work)
    ok = len(sinks) == 1
    verdict = "PASS" if ok else "FAIL"
    print(f"   sinks: {int(first.sum())} at the start, {len(sinks)} at the end (expected 1)  {verdict}")
    if ok:
        s = next(iter(sinks.values()))
        dm = abs(s["m"] / m0 - 1.0)
        dp = float(np.max(np.abs(s["m"] * np.array([s["u"], s["v"], s["w"]]) - p0)))
        status = "PASS" if dm <= TOL_MASS and dp <= TOL_COM else "FAIL"
        ok &= status == "PASS"
        print(f"   merged mass {s['m']:.6f} ({dm:.1e} from the sum), momentum {dp:.1e} from the sum  {status}")
    return ok


def agree(work_cpu: Path, work_fnl: Path, mass_only: bool = False) -> bool:
    """Check that the CPU and FNL sinks histories agree row by row, on physical scales: the positions against the cell
    size, the velocities against the largest speed, the masses against the largest mass (a column of round-off, e.g.
    z in a planar orbit, would make a relative gap meaningless). With mass_only (a sink formed at a symmetric centre:
    the cell that wins the tie is round-off), the masses and the final distance from the centre only."""
    a = read_history(next(work_cpu.glob("*-sinks_history.dat")))
    b = read_history(next(work_fnl.glob("*-sinks_history.dat")))
    if a["it"].size != b["it"].size or not (a["id"] == b["id"]).all() or not (a["it"] == b["it"]).all():
        print(f"   rows: CPU {a['it'].size}, FNL {b['it'].size} (the same steps and sinks expected)  FAIL")
        return False
    h = cell_size(work_cpu)
    gaps = {"time": float(np.max(np.abs(a["time"] - b["time"]))) / max(float(np.max(np.abs(a["time"]))), 1.0e-300),
            "m": float(np.max(np.abs(a["m"] - b["m"]))) / float(np.max(np.abs(a["m"])))}
    if mass_only:
        last = a["it"] == a["it"][-1]
        ra = np.sqrt(a["x"][last] ** 2 + a["y"][last] ** 2 + a["z"][last] ** 2)
        rb = np.sqrt(b["x"][last] ** 2 + b["y"][last] ** 2 + b["z"][last] ** 2)
        gaps["r_final"] = float(np.max(np.abs(ra - rb))) / h
    else:
        gaps["x"] = max(float(np.max(np.abs(a[k] - b[k]))) for k in ("x", "y", "z")) / h
        vmax = max(float(np.max(np.abs(a[k]))) for k in ("u", "v", "w"))
        gaps["v"] = max(float(np.max(np.abs(a[k] - b[k]))) for k in ("u", "v", "w")) / max(vmax, 1.0e-300)
    ok = all(g <= TOL_AGREE for k, g in gaps.items() if k != "r_final") and gaps.get("r_final", 0.0) <= 1.0e-2
    verdict = "PASS" if ok else "FAIL"
    text = ", ".join(f"{k} {g:.1e}" for k, g in gaps.items())
    print(f"   {a['it'].size} rows, gaps (positions in cells, velocities over the largest speed, masses over the "
          f"largest): {text} (bound {TOL_AGREE:.0e}; r_final 1e-2 cells)  {verdict}")
    return ok


def main() -> int:
    """Run the requested check."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("rest", "binary", "same", "accrete", "collapse", "none", "merge", "agree"))
    ap.add_argument("works", type=Path, nargs="+")
    ap.add_argument("--mass-only", action="store_true", help="agree: masses and the final radius only")
    a = ap.parse_args()
    if a.mode == "rest":
        ok = rest(a.works[0])
    elif a.mode == "binary":
        ok = binary(a.works[0])
    elif a.mode == "same":
        ok = same(a.works[0], a.works[1])
    elif a.mode == "agree":
        ok = agree(a.works[0], a.works[1], a.mass_only)
    else:
        ok = {"accrete": accrete, "collapse": collapse, "none": none, "merge": merge}[a.mode](a.works[0])
    print(f"SI {a.mode}: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
