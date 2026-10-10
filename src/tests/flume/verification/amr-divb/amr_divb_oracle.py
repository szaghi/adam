#!/usr/bin/env python3
"""Oracle of the AD cases (#68 A2): D_flux, the divergence the cleaning acts on, across static 2:1 seams.

Reads `<basename>-divb_history.dat` (columns: it, time, max_divb, l1_divb, seam_max_divb, max_dflux, l1_dflux,
seam_max_dflux, max_h_dflux_over_bmax; the D_flux columns need `[diagnostics] divb_flux`).

  growth  the source-free CPAW: per run, the seam D_flux floor (median over the second half of the run), its growth
          (the mean over the last quarter over the mean over the second quarter: about 1 if bounded, > 1 if growing) and
          its peak; across the refined ladder (N doubling) the order of the floor. A uniform run has no seam: its floor
          is that of max_dflux, the reference the seam floor is compared with. Assertions: --ratio-max on the growth of
          every refined run, --floor-order-min on the finest pair.
  loop    the field loop (#78): per run, the final seam and global D_flux and <|B_z|>/A0 (the MV-9 measure of the
          divergence error the out-of-plane velocity turns into B_z); a report, no assertion.

Usage:
    amr_divb_oracle.py growth --refined <work> ... [--uniform <work> ...] [--ratio-max R] [--floor-order-min P]
    amr_divb_oracle.py loop <work> ...
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path

import numpy as np

VERIF = Path(__file__).resolve().parent.parent
COL = {"time": 1, "max_divb": 2, "seam_max_divb": 4, "max_dflux": 5, "l1_dflux": 6, "seam_max_dflux": 7}


def history(work: Path) -> np.ndarray:
    """Return the rows of the div(B) history of a run."""
    files = list(work.glob("*-divb_history.dat"))
    if len(files) != 1:
        sys.exit(f"amr_divb_oracle: expected one div(B) history in {work}, found {len(files)}")
    rows = np.loadtxt(files[0], skiprows=1, ndmin=2)
    if rows.shape[1] < 9:
        sys.exit(f"amr_divb_oracle: {files[0].name} has no D_flux columns ([diagnostics] divb_flux)")
    return rows


def floor_growth(rows: np.ndarray, column: str) -> tuple[float, float, float]:
    """Return (floor, growth, peak) of a D_flux column over time (the step-0 row excluded: no register yet)."""
    t, v = rows[1:, COL["time"]], rows[1:, COL[column]]
    big = t[-1]

    def window(a: float, b: float) -> np.ndarray:
        return v[(t >= a * big) & (t <= b * big)]

    return (float(statistics.median(window(0.5, 1.0))), float(window(0.75, 1.0).mean() / window(0.25, 0.5).mean()),
            float(v.max()))


def growth(args: argparse.Namespace) -> int:
    """The CPAW report and assertions."""
    status, floors = 0, []
    for w in args.refined:
        f, g, p = floor_growth(history(w), "seam_max_dflux")
        floors.append(f)
        ok = g <= args.ratio_max
        status |= 0 if ok else 1
        print(f"{w.name}: seam D_flux floor {f:.4e}, growth {g:.3f} {'<=' if ok else '>'} {args.ratio_max}, "
              f"peak {p:.4e}  {'PASS' if ok else 'FAIL'}")
    orders = [math.log2(a / b) for a, b in zip(floors[:-1], floors[1:], strict=True)]
    if orders:
        print("   seam floor orders: " + " ".join(f"{o:+.2f}" for o in orders))
        if args.floor_order_min is not None:
            ok = orders[-1] >= args.floor_order_min
            status |= 0 if ok else 1
            print(f"   finest-pair floor order {orders[-1]:+.2f} >= {args.floor_order_min}: {'PASS' if ok else 'FAIL'}")
    for k, w in enumerate(args.uniform):
        f, g, _ = floor_growth(history(w), "max_dflux")
        line = f"{w.name}: uniform D_flux floor {f:.4e}, growth {g:.3f}"
        if k < len(floors):
            line += f"; refined seam floor / uniform floor {floors[k] / f:.2f}"
        print(line)
    return status


def loop(args: argparse.Namespace) -> int:
    """The field-loop (#78) report."""
    sys.path.insert(0, str(VERIF / "mhd" / "field-loop"))
    from field_loop_oracle import measures  # noqa: PLC0415

    for w in args.work:
        last = history(w)[-1]
        bz, energy = measures(w, 3)
        print(f"{w.name}: D_flux seam {last[COL['seam_max_dflux']]:.4e} max {last[COL['max_dflux']]:.4e} "
              f"L1 {last[COL['l1_dflux']]:.4e}; centred seam {last[COL['seam_max_divb']]:.4e}; "
              f"<|B_z|>/A0 {bz:.4e}; E_B(T)/E_B(0) {energy:.6f}")
    return 0


def main() -> int:
    """Parse the command line, run the requested report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="case", required=True)
    p = sub.add_parser("growth")
    p.add_argument("--refined", type=Path, nargs="+", required=True, help="refined runs, N doubling")
    p.add_argument("--uniform", type=Path, nargs="*", default=[], help="uniform controls, same N")
    p.add_argument("--ratio-max", type=float, default=1.5, help="bound of the growth of every refined run")
    p.add_argument("--floor-order-min", type=float, default=None, help="minimum finest-pair order of the floor")
    p.set_defaults(func=growth)
    p = sub.add_parser("loop")
    p.add_argument("work", type=Path, nargs="+")
    p.set_defaults(func=loop)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
