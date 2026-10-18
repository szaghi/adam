#!/usr/bin/env python3
"""FLUME passive-scalar oracle (issue #92): a run with passive scalars against the same run without them.

    ps_oracle.py <base work> <scalar work> --scalars N [--max-drift TOL]

Asserted: every field of the base run (r, ru, rv, rw, rE, and the entropy rs of euler-de) is bitwise equal in the
scalar run, the scalars do not act on the flow; the integral of each scalar in use is constant within --max-drift
(relative, the conservation history); the scalar slots beyond N are identically zero. Reported: the range of each
specific value s_k = (rho s_k) / rho.

    ps_oracle.py agree <cpu base> <fnl base> <cpu scalars> <fnl scalars>

CPU and FNL agree (P2): per field, the FNL-CPU difference of the scalar runs relative to the field's largest magnitude
is at most AGREE_FACTOR times the base pair's (the device arithmetic on the same case), or AGREE_FLOOR; the scalars,
which the base runs lack, against the largest base gap.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

NS_MAX = 4
AGREE_FACTOR = 2.0
AGREE_FLOOR = 1.0e-10


def last_blocks(work: Path) -> dict[tuple[float, ...], dict[str, np.ndarray]]:
    """Return every dataset of every block of the last checkpoint, keyed by the block origin."""
    files = sorted(work.glob("*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    out: dict[tuple[float, ...], dict[str, np.ndarray]] = {}
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in {k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}:
                key = tuple(np.round(h5[f"{blk}-origin"][()], 10))
                out[key] = {k.rsplit("-", 1)[1]: h5[k][()] for k in h5 if k.startswith(f"{blk}-")}
    return out


def history(work: Path) -> dict[str, np.ndarray]:
    """Return the conservation history, column name -> values."""
    lines = next(work.glob("*-conservation_history.dat")).read_text().splitlines()
    names = [n.strip('"') for n in lines[0].split("=", 1)[1].split()]
    data = np.array([[float(v) for v in row.split()] for row in lines[1:] if row.strip()])
    return {n: data[:, i] for i, n in enumerate(names)}


def agree(works: list[Path]) -> int:
    """CPU against FNL, the scalar runs against the base pair."""
    cb, fb, cs, fs = (last_blocks(w) for w in works)

    def gaps(a: dict, b: dict, names: list[str]) -> dict[str, float]:
        out = {}
        for v in names:
            d = max(float(np.abs(a[k][v] - b[k][v]).max()) for k in a)
            m = max(float(np.abs(a[k][v]).max()) for k in a)
            out[v] = d / m if m > 0.0 else d
        return out

    first = next(iter(cb.values()))
    flow = [v for v in ("r", "ru", "rv", "rw", "rE", "rs") if v in first]
    scal = [v for v in (f"rs{k}" for k in range(1, NS_MAX + 1)) if v in next(iter(cs.values()))]
    base = gaps(cb, fb, flow)
    run = gaps(cs, fs, flow + scal)
    ok = True
    for v, gap in run.items():
        ref = base.get(v, max(base.values()))
        bound = max(AGREE_FACTOR * ref, AGREE_FLOOR)
        ok &= gap <= bound
        print(f"   {v}: scalars run {gap:.3e}, base {ref:.3e}{'' if v in base else ' (largest)'}, bound {bound:.3e}")
    print(f"PS agree: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    """Run the checks."""
    if len(sys.argv) > 1 and sys.argv[1] == "agree":
        return agree([Path(w) for w in sys.argv[2:6]])
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("base", type=Path)
    ap.add_argument("scalars_run", type=Path)
    ap.add_argument("--scalars", type=int, required=True)
    ap.add_argument("--max-drift", type=float, default=1.0e-13)
    a = ap.parse_args()
    base, run = last_blocks(a.base), last_blocks(a.scalars_run)
    ok = sorted(base) == sorted(run)
    fields = [v for v in ("r", "ru", "rv", "rw", "rE", "rs") if v in next(iter(base.values()))]
    for v in fields:
        same = ok and all(np.array_equal(base[k][v], run[k][v]) for k in base)
        ok &= same
        print(f"   {v}: {'bitwise equal' if same else 'differs'} to the run without scalars")
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(a.scalars_run / "input.ini")
    g = int(ini["grid"]["ngc"])
    h = history(a.scalars_run)
    for k in range(1, NS_MAX + 1):
        name = f"rs{k}"
        if k <= a.scalars:
            col = h[f"int_{name}"]
            drift = float(np.max(np.abs(col - col[0])) / max(abs(col[0]), np.finfo(float).tiny))
            vals = np.concatenate([(b[name] / b["r"])[g:-g, g:-g, g:-g].ravel() for b in run.values()])
            ok &= drift <= a.max_drift
            print(f"   s{k}: integral drift {drift:.2e} (bound {a.max_drift:g}); range [{vals.min():.3e}, {vals.max():.6f}]")
        else:
            zero = all(not np.any(b[name]) for b in run.values())
            ok &= zero
            print(f"   rs{k}: {'identically zero' if zero else 'NOT zero'} (unused slot)")
    print(f"PS: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
