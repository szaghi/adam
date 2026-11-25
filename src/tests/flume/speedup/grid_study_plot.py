#!/usr/bin/env python3
"""Draw the grid-size study of the FLUME uniform sphere from analyze.py's summary.csv files (issue #133, 2026-11-24).

One figure, `<prefix>.png` in --out (and its silver-ink twin `<prefix>-dark.png` with --theme dark or both,
docs/scripts/figure_theme.py), two panels:

    (a) block size  the grind time g = step cost x devices / N on 1 A30 against the block edge b, for the L3 grids
                    (512 blocks of b^3 cells, `wg-b<b>`) and the same 192^3 grid cut into 64 or 4096 blocks
                    (`wg-l2b48`, `wg-l4b12`), with the least-squares fit g = a + c ((b + 2 ngc) / b)^3 over the L3
                    cases (ngc = 3, WENO-5), drawn over its range only: the 4096-block cut lies 23 ns below its
                    extrapolation, so it describes the L3 series and is no law of b
    (b) GPUs        g against b for the L3 grids on 1, 2 and 4 A30 of one node (the matrix's launch): at one b the grid
                    is the same, the gap is the GPUs'; and, from --bound / --unbound summaries, the A/B of the rank
                    placement on 4 A30 (one node) and 8 A30 (two nodes), diamonds for the grid-study cases and stars
                    for the matrix's 2a (L4, 24^3 blocks): unbound is the matrix's launch, bound pins each rank to its
                    A30's NUMA node (README, "Results"); filled bound, open unbound

Needs matplotlib and numpy (run it on WSL, e.g. exe/.regression-venv/bin/python); analyze.py stays stdlib.

Usage:
    grid_study_plot.py <grid-study summary.csv> --out DIR [--bound CSV ...] [--unbound CSV ...] [--prefix NAME]
                       [--theme light|dark|both]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
mpl.rcParams["figure.facecolor"] = "#f3ede1"  # the washi ground of the docs (docs/BRAND.md)
mpl.rcParams["savefig.facecolor"] = "#f3ede1"
import matplotlib.pyplot as plt  # noqa: E402

DOCS_SCRIPTS = Path(__file__).resolve().parents[4] / "docs" / "scripts"
sys.path.insert(0, str(DOCS_SCRIPTS))
from figure_theme import PAPER_PALETTE, theme, themes  # noqa: E402

NGC = 3  # ghost layers of WENO-5
DPI = 130
NP_COLOUR = {1: PAPER_PALETTE[0], 2: PAPER_PALETTE[2], 4: PAPER_PALETTE[1], 8: PAPER_PALETTE[3]}

Row = dict[str, str]


def read_rows(paths: list[Path]) -> list[Row]:
    """Return the FNL rows of the summary.csv files.

    Parameters
    ----------
    paths
        analyze.py's summary.csv files.

    Returns
    -------
    list[Row]
        Their fnl-nvf rows, in file order.
    """
    rows: list[Row] = []
    for path in paths:
        with path.open(newline="") as f:
            rows += [r for r in csv.DictReader(f) if r["backend"] == "fnl-nvf"]
    return rows


def block_edge(case: str) -> int | None:
    """Return the block edge of a grid-study case (`wg-b<b>`, `wg-l<l>b<b>`), or None for another case."""
    if not case.startswith("wg-"):
        return None
    return int(case.rsplit("b", 1)[1])


def ghost_ratio(b: np.ndarray) -> np.ndarray:
    """Return the stored over interior cells of a b^3 block, ((b + 2 ngc) / b)^3."""
    return ((b + 2 * NGC) / b) ** 3


def fit(rows: list[Row]) -> tuple[float, float]:
    """Return a, c of the least-squares fit g = a + c ghost_ratio(b) over the L3 cases on 1 A30."""
    pts = [(block_edge(r["case"]), float(r["ns_per_cell_step"])) for r in rows
           if r["np"] == "1" and r["case"].startswith("wg-b")]
    b = np.array([p[0] for p in pts], dtype=float)
    g = np.array([p[1] for p in pts])
    design = np.vstack([np.ones_like(b), ghost_ratio(b)]).T
    (a, c), *_ = np.linalg.lstsq(design, g, rcond=None)
    return float(a), float(c)


def draw(grid: list[Row], bound: list[Row], unbound: list[Row], out: Path) -> None:
    """Draw the two panels into `out` (the theme decides the file: figure_theme.theme)."""
    fig, (ax_b, ax_n) = plt.subplots(1, 2, figsize=(9.6, 4.6))
    one = [r for r in grid if r["np"] == "1" and block_edge(r["case"]) is not None]
    l3 = sorted((r for r in one if r["case"].startswith("wg-b")), key=lambda r: block_edge(r["case"]) or 0)
    same = [r for r in one if r["case"] in ("wg-l4b12", "wg-b24", "wg-l2b48")]
    ax_b.plot([block_edge(r["case"]) for r in l3], [float(r["ns_per_cell_step"]) for r in l3], "o",
              color=PAPER_PALETTE[0], label="L3, 512 blocks of $b^3$ (2.1-28.1 M cells)")
    ax_b.plot([block_edge(r["case"]) for r in same], [float(r["ns_per_cell_step"]) for r in same], "s", mfc="none",
              ms=9, color=PAPER_PALETTE[1], label="the same $192^3$ grid in 4096, 512, 64 blocks")
    a, c = fit(grid)
    bs = np.linspace(16, 38, 100)
    ax_b.plot(bs, a + c * ghost_ratio(bs), "-", color="0.45", lw=1,
              label=f"fit over L3: $g = {a:.1f} + {c:.1f}\\,((b+6)/b)^3$ ns")
    ax_b.set_xlabel("block edge $b$ (cells)", fontsize=8)
    ax_b.set_ylabel("g, ns per cell-step (1 A30)", fontsize=8)
    ax_b.set_title("(a) block size, 1 A30", fontsize=9)
    ax_b.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    ax_b.grid(True, alpha=0.3)
    for np_ in (1, 2, 4):
        pts = sorted((r for r in grid if r["np"] == str(np_) and r["case"].startswith("wg-b")),
                     key=lambda r: block_edge(r["case"]) or 0)
        ax_n.plot([block_edge(r["case"]) for r in pts], [float(r["ns_per_cell_step"]) for r in pts], "o-", lw=0.8,
                  color=NP_COLOUR[np_], label=f"{np_} A30, one node (matrix launch)")
    for rows, fill, what in ((unbound, False, "unbound"), (bound, True, "bound")):
        for r in rows:
            np_ = int(r["np"])
            b = block_edge(r["case"]) or 24  # the matrix's 2a: L4, 24^3 blocks
            colour = NP_COLOUR.get(np_, "0.3")
            matrix = not r["case"].startswith("wg-")
            label = None
            if r["case"] == "wg-b38":
                label = f"{np_} A30 {what}"
            elif matrix and np_ == 4 and fill:
                label = "stars: 2a (L4, 24$^3$ blocks)"
            ax_n.plot(b, float(r["ns_per_cell_step"]), "*" if matrix else "D", ms=8 if matrix else 7, color=colour,
                      mfc=colour if fill else "none", label=label)
    ax_n.set_xlabel("block edge $b$ (cells)", fontsize=8)
    ax_n.set_ylabel("g, ns per cell-step", fontsize=8)
    ax_n.set_title("(b) the same grids on 1, 2, 4 A30 (lines), rank placement (marks)", fontsize=9)
    ax_n.legend(fontsize=6.5, frameon=False, ncols=3, loc="upper center", bbox_to_anchor=(0.5, -0.14))
    ax_n.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def main() -> int:
    """Parse the command line, read the summaries, write the figure."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, help="the grid study's summary.csv")
    parser.add_argument("--out", type=Path, required=True, help="the directory of the figure")
    parser.add_argument("--bound", type=Path, nargs="*", default=[], help="summary.csv of bound runs")
    parser.add_argument("--unbound", type=Path, nargs="*", default=[], help="summary.csv of unbound runs")
    parser.add_argument("--prefix", default="perf-grid-study", help="file name (default perf-grid-study)")
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="light",
                        help="light: <name>.png on washi; dark: <name>-dark.png on silver ink; both")
    args = parser.parse_args()
    grid = read_rows([args.csv])
    bound, unbound = read_rows(args.bound), read_rows(args.unbound)
    args.out.mkdir(parents=True, exist_ok=True)
    for t in themes(args.theme):
        with theme(t):
            draw(grid, bound, unbound, args.out / f"{args.prefix}.png")
    a, c = fit(grid)
    print(f"grid_study_plot: {args.prefix}.png in {args.out} ({args.theme}); fit a {a:.2f} c {c:.2f} ns")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
