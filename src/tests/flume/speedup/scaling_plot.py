#!/usr/bin/env python3
"""Plot the grind-time curves of the FLUME speedup matrix from analyze.py's summary.csv (issue #133 phase 9).

One panel per backend: the grind time g (ns per cell-step, step cost x devices / N) against the cells per device, log x,
one colour per scaling series (the cases of matrix.conf with a `series`: a leg's uniform case and its weak cases) and
grey for the windows of no series. A filled marker is a run on one node, an open one a multi-node run (more devices
than --node-size of its backend); every point is labelled with its device count. Read the multi-node overhead as the
distance of the open markers from the filled curve at the same load.

Needs matplotlib (not stdlib: analyze.py stays stdlib; run this where matplotlib is installed, e.g. WSL after copying
`results/` back).

Usage:
    scaling_plot.py <summary.csv> [--out scaling.png] [--node-size cpu=64,fnl-nvf=4,fnl-amd=8] [--series S,...]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def node_sizes(text: str) -> dict[str, int]:
    """Return the devices per node of each backend from `backend=n,...`.

    Parameters
    ----------
    text : str
        Comma-separated `backend=devices` pairs.

    Returns
    -------
    dict[str, int]
        Backend -> devices per node.
    """
    sizes = {}
    for item in text.split(","):
        name, _, value = item.partition("=")
        sizes[name.strip()] = int(value)
    return sizes


def main() -> int:
    """Parse the command line, read summary.csv, write the figure."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, help="analyze.py's summary.csv")
    parser.add_argument("--out", type=Path, default=Path("scaling.png"), help="the figure (default scaling.png)")
    parser.add_argument("--node-size", default="cpu=64,fnl-nvf=4,fnl-amd=8", help="devices per node per backend")
    parser.add_argument("--series", default="", help="comma-separated series to plot (default: every series)")
    args = parser.parse_args()
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("scaling_plot: matplotlib is not installed (analyze.py's tables need none)", file=sys.stderr)
        return 1
    with args.csv.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if r.get("cells_per_device")]
    if not rows:
        print(f"scaling_plot: no row with cells_per_device in {args.csv} (an analyze.py before phase 9?)",
              file=sys.stderr)
        return 1
    wanted = {s for s in args.series.split(",") if s}
    sizes = node_sizes(args.node_size)
    backends = sorted({r["backend"] for r in rows})
    fig, axes = plt.subplots(1, len(backends), figsize=(5.5 * len(backends), 4.2), squeeze=False)
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    series_all = sorted({r["series"] for r in rows if r["series"]})
    for ax, backend in zip(axes[0], backends, strict=True):
        mine = [r for r in rows if r["backend"] == backend]
        groups = sorted({r["series"] for r in mine}, key=lambda s: (s == "", s))
        for series in groups:
            if wanted and series not in wanted:
                continue
            points = sorted((r for r in mine if r["series"] == series), key=lambda r: float(r["cells_per_device"]))
            colour = colours[series_all.index(series) % len(colours)] if series else "0.6"
            for multi in (False, True):
                sel = [r for r in points if (int(r["np"]) > sizes.get(backend, 1)) == multi]
                if not sel:
                    continue
                ax.plot([float(r["cells_per_device"]) for r in sel], [float(r["ns_per_cell_step"]) for r in sel],
                        "o" if not multi else "s", color=colour, mfc="none" if multi else colour, ls="none",
                        label=f"{series or 'no series'}{' (multi-node)' if multi else ''}")
            for r in points:
                ax.annotate(r["np"], (float(r["cells_per_device"]), float(r["ns_per_cell_step"])),
                            textcoords="offset points", xytext=(4, 4), fontsize=7)
        ax.set_xscale("log")
        ax.set_xlabel("cells per device")
        ax.set_ylabel("g, ns per cell-step (step cost x devices / N)")
        ax.set_title(backend)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"scaling_plot: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
