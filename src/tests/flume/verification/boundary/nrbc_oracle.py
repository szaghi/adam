#!/usr/bin/env python3
"""Judge the non-reflecting boundary runs (issue #121, S2) from their last checkpoint.

Subcommands:
    residual <work> --amplitude A [--bound B]   the largest |p - p_inf| left in the box, over gamma p_inf A (the packet's
                                                pressure amplitude): what the faces sent back once the packet has left
    relax <work> --deviation D [--bound B | --above B]  the box-mean |p / p_inf - 1| over the initial deviation D
    ratio <a> <b> --amplitude A --min L --max H the residual of `a` over the one of `b`, within [L, H]
    bitwise <a> <b>                             the last checkpoints of two runs equal bit for bit
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

GAMMA = 1.4
NGC = 3
P_INF = 1.0 / GAMMA


def last_checkpoint(work: Path) -> str:
    """Return the iteration tag of the latest checkpoint (by time)."""
    best = None
    for x in work.glob("nrbc-*.xdmf"):
        m = re.search(r'<Time Value="([^"]+)"', x.read_text()[:4000])
        if m and x.stem.split("-")[1].isdigit() and (best is None or float(m.group(1)) > best[0]):
            best = (float(m.group(1)), x.stem.split("-")[1])
    if best is None:
        raise SystemExit(f"nrbc_oracle: no checkpoint in {work}")
    return best[1]


def pressure(work: Path) -> np.ndarray:
    """The interior pressure of every cell of the last checkpoint, flattened."""
    import h5py

    tag = last_checkpoint(work)
    out = []
    for f5 in sorted(work.glob(f"nrbc-{tag}-proc*.h5")):
        with h5py.File(f5, "r") as h:
            for b in {k.rsplit("-", 1)[0] for k in h}:
                s = (slice(NGC, -NGC),) * 3
                r, ru, rv, rw, re_ = (h[f"{b}-{v}"][()][s] for v in ("r", "ru", "rv", "rw", "rE"))
                out.append(((GAMMA - 1.0) * (re_ - 0.5 * (ru**2 + rv**2 + rw**2) / r)).ravel())
    return np.concatenate(out)


def residual(work: Path, amplitude: float, bound: float | None) -> int:
    """The reflection left in the box, relative to the packet's pressure amplitude."""
    x = np.abs(pressure(work) - P_INF).max() / (GAMMA * P_INF * amplitude)
    ok = bound is None or x < bound
    print(f"   {work.name}: max |p - p_inf| / (gamma p_inf A) = {x:.3e}" + ("" if bound is None else f" (bound {bound:g})  {'PASS' if ok else 'FAIL'}"))
    return 0 if ok else 1


def relax(work: Path, deviation: float, bound: float | None, above: float | None) -> int:
    """The box-mean pressure deviation left, relative to the initial one: below `bound`, or above `above`."""
    x = np.abs(pressure(work) / P_INF - 1.0).mean() / deviation
    ok = (bound is None or x < bound) and (above is None or x > above)
    rule = f" (bound {bound:g})" if bound is not None else f" (min {above:g})" if above is not None else ""
    print(f"   {work.name}: mean |p / p_inf - 1| / initial = {x:.3e}{rule}" + (f"  {'PASS' if ok else 'FAIL'}" if rule else ""))
    return 0 if ok else 1


def ratio(a: Path, b: Path, amplitude: float, lo: float, hi: float) -> int:
    """The residual of `a` over the one of `b`, within [lo, hi]."""
    ra = np.abs(pressure(a) - P_INF).max()
    rb = np.abs(pressure(b) - P_INF).max()
    x = ra / rb
    ok = lo <= x <= hi
    print(f"   {a.name} / {b.name}: residual ratio {x:.3f} (within [{lo:g}, {hi:g}])  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def bitwise(a: Path, b: Path) -> int:
    """The last checkpoints of `a` and `b`, every dataset bit for bit."""
    import h5py

    ta, tb = last_checkpoint(a), last_checkpoint(b)
    n = bad = 0
    for fa in sorted(a.glob(f"nrbc-{ta}-proc*.h5")):
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
    r = sub.add_parser("residual")
    r.add_argument("work", type=Path)
    r.add_argument("--amplitude", type=float, required=True)
    r.add_argument("--bound", type=float)
    x = sub.add_parser("relax")
    x.add_argument("work", type=Path)
    x.add_argument("--deviation", type=float, required=True)
    x.add_argument("--bound", type=float)
    x.add_argument("--above", type=float)
    q = sub.add_parser("ratio")
    q.add_argument("a", type=Path)
    q.add_argument("b", type=Path)
    q.add_argument("--amplitude", type=float, required=True)
    q.add_argument("--min", type=float, required=True)
    q.add_argument("--max", type=float, required=True)
    b = sub.add_parser("bitwise")
    b.add_argument("a", type=Path)
    b.add_argument("b", type=Path)
    a = ap.parse_args()
    if a.cmd == "residual":
        return residual(a.work, a.amplitude, a.bound)
    if a.cmd == "relax":
        return relax(a.work, a.deviation, a.bound, a.above)
    if a.cmd == "ratio":
        return ratio(a.a, a.b, a.amplitude, a.min, a.max)
    return bitwise(a.a, a.b)


if __name__ == "__main__":
    sys.exit(main())
