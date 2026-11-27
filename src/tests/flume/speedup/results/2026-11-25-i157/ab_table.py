#!/usr/bin/env python3
"""The A/B tables of #157: step costs and AMR ratios of two analyze.py summaries (before and after) side by side.

Usage:
    ab_table.py <summary-A.csv> <summary-B.csv>
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path


def rows(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    """The step-cost rows of a summary.csv, keyed by (case, ranks)."""
    with path.open() as f:
        return {(r["case"], r["np"]): r for r in csv.DictReader(f) if r.get("case")}


def main() -> int:
    """Print the step-cost table and the per-leg ratios, before -> after."""
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    a, b = rows(Path(sys.argv[1])), rows(Path(sys.argv[2]))
    print(
        "| Case | Before (s/step) | After (s/step) | Change | Steps before -> after | Regrid before -> after |"
    )
    print("|---|---|---|---|---|---|")
    for key in sorted(a, key=lambda k: (k[0], int(k[1]))):
        if key not in b:
            continue
        ta, tb = float(a[key]["step_cost_s"]), float(b[key]["step_cost_s"])
        print(
            f"| {key[0]} ({key[1]} GPUs) | {ta:.4f} | {tb:.4f} | {100 * (tb / ta - 1):+.1f} % "
            f"| {float(a[key]['steps_cost_s']):.4f} -> {float(b[key]['steps_cost_s']):.4f} "
            f"| {float(a[key]['regrid_cost_s']):.4f} -> {float(b[key]['regrid_cost_s']):.4f} |"
        )
    print()
    legs = sorted({(r["leg"], r["np"]) for r in a.values()})
    print(
        "| Leg (GPUs) | AMR speedup T(U)/T(R) | Ideal N(U)/N(R) | Grid efficiency | Regrid cost T(R)/T(F) - 1 |"
    )
    print("|---|---|---|---|---|")
    for leg, ranks in legs:
        out = []
        for side in (a, b):
            by_role = {
                r["role"][0].upper(): r  # uniform, runtime, frozen -> U, R, F
                for r in side.values()
                if r["leg"] == leg and r["np"] == ranks
            }
            if not {"U", "R", "F"} <= set(by_role):
                break
            tu, tr, tf = (float(by_role[x]["step_cost_s"]) for x in "URF")
            nu, nr, nf = (float(by_role[x]["cells_mean"]) for x in "URF")
            out.append((tu / tr, nu / nr, tu * nf / nu / tf, tr / tf - 1))
        if len(out) == 2:
            (s0, i0, e0, r0), (s1, _, e1, r1) = out
            print(
                f"| {leg} ({ranks}) | {s0:.2f} -> {s1:.2f} | {i0:.2f} | {e0:.3f} -> {e1:.3f} | {100 * r0:.1f} -> {100 * r1:.1f} % |"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
