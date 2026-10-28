#!/usr/bin/env python3
"""Draw the corner ratios of Boss-Bodenheimer runs (issue #101): the largest speed and the smallest density of the cells
near the edges of the refined cube, over those of the rest of the same radial shell, per checkpoint.

Why: the corner artifact of issue #101 is a transient (steps 140 to 300 of 400), invisible at the end of the run; the
figure shows it, and its absence, over the whole run. The ratios come from `eos_oracle.py corners <work> --csv <file>`,
run where the checkpoints are (they are too large to copy).

Usage:
    corner_figure.py <out.png> <label>=<ratios.csv> [<label>=<ratios.csv> ...] [--max-speed-ratio R] [--min-density-ratio R]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def read(path: Path) -> dict[str, list[float]]:
    """Return the columns of a ratios CSV."""
    with path.open() as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def main() -> int:
    """Draw the figure."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("runs", nargs="+", help="label=ratios.csv")
    ap.add_argument("--max-speed-ratio", type=float, default=2.0)
    ap.add_argument("--min-density-ratio", type=float, default=0.5)
    args = ap.parse_args()
    import matplotlib as mpl  # noqa: PLC0415

    mpl.use("Agg")
    import matplotlib.pyplot as plt  # noqa: PLC0415

    fig, axs = plt.subplots(1, 2, figsize=(12.0, 4.0), constrained_layout=True)
    for spec in args.runs:
        label, _, path = spec.partition("=")
        d = read(Path(path))
        axs[0].semilogy(d["step"], d["speed_ratio"], "o-", ms=3, label=label)
        axs[1].semilogy(d["step"], d["density_ratio"], "o-", ms=3, label=label)
    axs[0].axhline(
        args.max_speed_ratio,
        color="k",
        ls="--",
        lw=1,
        label=f"bound {args.max_speed_ratio:g}",
    )
    axs[1].axhline(
        args.min_density_ratio,
        color="k",
        ls="--",
        lw=1,
        label=f"bound {args.min_density_ratio:g}",
    )
    axs[0].set_ylabel("largest speed: edge cells / rest of the shell")
    axs[1].set_ylabel("smallest density: edge cells / rest of the shell")
    axs[0].set_title("speed (bound: at most the dashed line)", fontsize=9)
    axs[1].set_title("density (bound: at least the dashed line)", fontsize=9)
    for ax in axs:
        ax.set_xlabel("step")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    fig.suptitle(
        "Boss-Bodenheimer cloud, isothermal, 8 ranks: cells near the edges of the refined cube\n"
        "against the rest of the shell 1.55 < r < 1.85 (issue #101)",
        fontsize=9,
    )
    fig.savefig(args.out, dpi=110)
    print(f"corner_figure: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
