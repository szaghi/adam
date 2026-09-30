#!/usr/bin/env python3
"""Oracle of the field loop advection (issue #41, M2-P6, MV-9; Gardiner and Stone 2005, Mignone and Tzeferacos 2010).

Why: with the out-of-plane velocity w the equations give dB_z/dt = w div(B), so the mean |B_z| of the last checkpoint,
relative to the loop amplitude A0, measures the accumulated divergence error: it must stay below --bz-max and decrease
under refinement (--bz-order-min on the finest pair, a positive rate). The magnetic energy E_B = sum |B|^2 / 2 of the
last checkpoint over the first measures the numerical dissipation of the loop: at least --energy-min.

Usage:
    field_loop_oracle.py <work> [<work> ...] [--bz-max E [E ...]] [--bz-order-min P] [--energy-min R [R ...]]
                         [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "linear-wave"))
from linear_wave_oracle import NAMES, profile, volumes  # noqa: E402

IBX, IBY, IBZ = NAMES.index("bx"), NAMES.index("by"), NAMES.index("bz")


def amplitude(work: Path) -> float:
    """Return the loop amplitude A0 of the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return float(ini["initial_conditions"]["loop_amplitude"])


def measures(work: Path, ngc: int) -> tuple[float, float]:
    """Return (<|B_z|> / A0 of the last checkpoint, E_B(last) / E_B(first)), volume-weighted (AMR runs, M2-P7b, have
    cells of two sizes; on a uniform grid the weights are equal)."""
    out = []
    for step in ("first", "last"):
        cells, vol = profile(work, step, ngc), volumes(work, step, ngc)
        keys = sorted(cells)
        out.append((np.array([cells[k] for k in keys]), np.array([vol[k] for k in keys])))
    (first, w_first), (last, w_last) = out
    bz = float(np.sum(np.abs(last[:, IBZ]) * w_last) / np.sum(w_last)) / amplitude(work)
    energy = [float(np.sum((q[:, IBX] ** 2 + q[:, IBY] ** 2 + q[:, IBZ] ** 2) * w))
              for q, w in ((first, w_first), (last, w_last))]
    return bz, energy[1] / energy[0]


def bound(values: list[float] | None, n: int) -> float | None:
    """Return the n-th bound of a per-run list (one value applies to every run)."""
    if values is None:
        return None
    return values[0] if len(values) == 1 else values[n]


def main() -> int:
    """Run the MV-9 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+", help="runs, coarse first, N doubling")
    parser.add_argument("--bz-max", type=float, nargs="+", default=None, help="bound of <|B_z|> / A0 (one or per run)")
    parser.add_argument("--bz-order-min", type=float, default=None, help="minimum refinement rate of <|B_z|>")
    parser.add_argument("--energy-min", type=float, nargs="+", default=None, help="bound of E_B(T) / E_B(0)")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    for name in ("bz_max", "energy_min"):
        values = getattr(args, name)
        if values is not None and len(values) not in (1, len(args.work)):
            sys.exit(f"field_loop_oracle: --{name.replace('_', '-')} takes one bound or one per run")
    ok = True
    bzs = []
    for n, w in enumerate(args.work):
        bz, ratio = measures(w, args.ngc)
        bzs.append(bz)
        line = f"{w.name}: <|B_z|>/A0 {bz:.4e}  E_B(T)/E_B(0) {ratio:.6f}"
        if n > 0:
            line += f"  <|B_z|> rate {math.log2(bzs[-2] / bz):+.2f}"
        print(line)
        b = bound(args.bz_max, n)
        if b is not None:
            good = bz <= b
            ok &= good
            print(f"   <|B_z|>/A0 {bz:.4e} {'<=' if good else '>'} {b:.3e}: {'PASS' if good else 'FAIL'}")
        b = bound(args.energy_min, n)
        if b is not None:
            good = ratio >= b
            ok &= good
            print(f"   E_B(T)/E_B(0) {ratio:.6f} {'>=' if good else '<'} {b:.6f}: {'PASS' if good else 'FAIL'}")
    if args.bz_order_min is not None and len(bzs) > 1:
        p = math.log2(bzs[-2] / bzs[-1])
        good = p >= args.bz_order_min
        ok &= good
        print(f"   <|B_z|> finest-pair rate {p:+.2f} {'>=' if good else '<'} {args.bz_order_min}: "
              f"{'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
