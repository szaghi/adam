#!/usr/bin/env python3
"""Draw the ghost-fill microbenchmark figures (ghost_fill_bench.F90, README.md, 2026-11-25).

Reads the result files of run.sh (one per GPU: `results/h100-*.txt`, `results/a30-*.txt`) and, for the projection on a
FLUME step, the phase-timer reports of the H100 runs (`src/tests/flume/speedup/results/2026-11-25-spacehpc/phase/
reports/`). Writes into --out, each with its silver-ink twin `<name>-dark.png` under --theme both
(docs/scripts/figure_theme.py):

    ghost-fill-floor.png    the fraction of the GPU's copy floor every variant reaches against the block edge, at
                            ~7 M cells per GPU (lines) and ~56 M (open symbols); H100 and A30 side by side
    ghost-fill-speedup.png  today's fill time over the best structured variant's, and over the 32-bit map's, against
                            the block edge; both loads, both GPUs
    ghost-fill-step.png     a FLUME step on the H100 split into its phases (face fluxes, local fill, BC and blend, MPI
                            exchange, the rest), today, with the local fill at the structured-slab speed, and with
                            that plus a face kernel twice as fast (a scenario, not a measurement); 38^3 blocks, 7 M and
                            56 M cells per GPU, with and without the immersed body

Needs matplotlib and numpy (on WSL: exe/.regression-venv/bin/python).

Usage:
    plot.py [--results DIR] [--phase DIR] --out DIR [--theme light|dark|both]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["figure.facecolor"] = (
    "#f3ede1"  # the washi ground of the docs (docs/BRAND.md)
)
mpl.rcParams["savefig.facecolor"] = "#f3ede1"
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "docs" / "scripts"))
from figure_theme import PAPER_MARKERS, PAPER_PALETTE, theme, themes  # noqa: E402

DPI = 130
VARIANTS = (
    "V0 ADAM map I8P",
    "V1 map I4P",
    "V2 offsets, v inner",
    "V3 slabs, b fastest",
    "V4 slabs, cell fastest",
    "V5 slabs 1 launch, b fast",
    "V6 slabs 1 launch, cell f",
)
LABEL = {
    "V0 ADAM map I8P": "V0 today (per-cell map, I8P)",
    "V1 map I4P": "V1 map in I4P",
    "V2 offsets, v inner": "V2 offsets, variables inside",
    "V3 slabs, b fastest": "V3 slabs, 26 launches",
    "V4 slabs, cell fastest": "V4 slabs, 26 launches, cell-fastest",
    "V5 slabs 1 launch, b fast": "V5 slabs, 1 launch",
    "V6 slabs 1 launch, cell f": "V6 slabs, 1 launch, cell-fastest",
}
STYLE = {v: (PAPER_PALETTE[k], PAPER_MARKERS[k]) for k, v in enumerate(VARIANTS)}
HEAD = re.compile(
    r"# blocks\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+) b\s+(\d+) ngc\s+(\d+) nv\s+(\d+) cells\s+([\d.]+) M"
)
ROW = re.compile(
    r"^\s*(.+?) median\s+([\d.]+) ms\s+min\s+([\d.]+) ms\s+([\d.]+) GB/s\s+of ref\s+([-\d.]+)\s+bitwise (\w)"
)
PHASE = re.compile(
    r"PHASE TIMING: (.{14})\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+) %"
)

Case = dict[str, object]


def read_results(path: Path) -> list[Case]:
    """Return the cases of one run.sh result file.

    Parameters
    ----------
    path
        A result file.

    Returns
    -------
    list of dict
        Per case: b, ngc, cells (M), and per variant (and "ref contiguous copy") the median time (ms), the fraction
        of the copy floor and the effective bandwidth (TB/s, 16 B per ghost value); every bitwise flag must be T.
    """
    cases: list[Case] = []
    for line in path.read_text().splitlines():
        if m := HEAD.search(line):
            cases.append(
                {
                    "b": int(m.group(5)),
                    "ngc": int(m.group(6)),
                    "cells": float(m.group(8)),
                    "t": {},
                    "f": {},
                    "bw": {},
                }
            )
        elif (m := ROW.search(line)) and cases:
            name = m.group(1).strip()
            if m.group(6) != "T":
                sys.exit(
                    f"plot.py: {path.name}: {name} at b {cases[-1]['b']} is not bitwise identical"
                )
            cases[-1]["t"][name] = float(m.group(2))
            cases[-1]["f"][name] = float(m.group(5))
            cases[-1]["bw"][name] = float(m.group(4)) / 1e3  # TB/s
    return cases


def read_phase(path: Path) -> dict[str, float]:
    """Return the mean ms per step of every phase of one phase-timer report."""
    return {
        m.group(1).strip(): float(m.group(3))
        for line in path.read_text().splitlines()
        if (m := PHASE.search(line))
    }


def best(case: Case) -> float:
    """The fastest structured variant's time (V3, V5: ADAM's block-fastest layout)."""
    return min(case["t"]["V3 slabs, b fastest"], case["t"]["V5 slabs 1 launch, b fast"])


def fig_floor(runs: dict[str, list[Case]], out: Path) -> None:
    """Fraction of the copy floor per variant against the block edge, one panel per GPU."""
    fig, axes = plt.subplots(
        1, 2, figsize=(10.5, 4.2), sharey=True, layout="constrained"
    )
    for ax, (gpu, cases) in zip(axes, runs.items(), strict=True):
        small = sorted(
            (c for c in cases if c["cells"] < 20 and c["ngc"] == 3),
            key=lambda c: c["b"],
        )
        large = sorted((c for c in cases if c["cells"] > 20), key=lambda c: c["b"])
        for v in VARIANTS:
            col, mk = STYLE[v]
            ax.plot(
                [c["b"] for c in small],
                [c["f"][v] for c in small],
                marker=mk,
                color=col,
                lw=1.4,
                ms=5,
                label=LABEL[v],
            )
            ax.plot(
                [c["b"] for c in large],
                [c["f"][v] for c in large],
                marker=mk,
                color=col,
                ls="none",
                ms=7,
                mfc="none",
                mew=1.2,
            )
        ax.axhline(1.0, color="0.45", lw=0.8, ls=":")
        ax.text(16, 1.02, "contiguous copy = 1", fontsize=8, color="0.35", va="bottom")
        floor = [c["bw"]["ref contiguous copy"] for c in cases]
        ax.set_title(
            f"{gpu}: copy floor {min(floor):.2f}-{max(floor):.2f} TB/s", fontsize=10
        )
        ax.set_xlabel("block edge b (cells)")
        ax.set_xticks([16, 24, 32, 38, 48])
        ax.set_ylim(0, 1.12)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("fraction of the copy floor")
    axes[0].plot(
        [],
        [],
        color="0.4",
        marker="o",
        ls="none",
        mfc="none",
        label="open: ~56 M cells per GPU",
    )
    axes[0].legend(
        fontsize=7.5, loc="upper left", bbox_to_anchor=(0.0, 0.93), frameon=False
    )
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def fig_speedup(runs: dict[str, list[Case]], out: Path) -> None:
    """Today's time over the best structured variant's (solid) and over the 32-bit map's (dotted)."""
    fig, ax = plt.subplots(figsize=(7.2, 4.2), layout="constrained")
    for k, (gpu, cases) in enumerate(runs.items()):
        for load, sel, mk in (
            ("~7 M", lambda c: c["cells"] < 20 and c["ngc"] == 3, "o"),
            ("~56 M", lambda c: c["cells"] > 20, "s"),
        ):
            cs = sorted((c for c in cases if sel(c)), key=lambda c: c["b"])
            col = PAPER_PALETTE[0] if k == 0 else PAPER_PALETTE[1]
            ax.plot(
                [c["b"] for c in cs],
                [c["t"]["V0 ADAM map I8P"] / best(c) for c in cs],
                marker=mk,
                color=col,
                lw=1.5,
                ms=6,
                mfc=col if load == "~7 M" else "none",
                label=f"{gpu}, {load} cells: slabs",
            )
            ax.plot(
                [c["b"] for c in cs],
                [c["t"]["V0 ADAM map I8P"] / c["t"]["V1 map I4P"] for c in cs],
                marker=mk,
                color=col,
                lw=1.0,
                ls=":",
                ms=4,
                mfc="none",
                label=f"{gpu}, {load} cells: I4P map",
            )
    ax.axhline(1.0, color="0.45", lw=0.8)
    ax.set_xlabel("block edge b (cells)")
    ax.set_xticks([16, 24, 32, 38, 48])
    ax.set_ylabel("speedup of the fill over today's")
    ax.set_ylim(0.8, 6.0)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, ncol=2, frameon=False, loc="upper left")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def fig_step(phase_dir: Path, h100: list[Case], out: Path) -> None:
    """A FLUME step on 4 H100 (38^3 blocks) by phase: today, slab fill, slab fill + face kernel 2x (scenario)."""

    def ratio(cells: float) -> float:
        c = min(
            (c for c in h100 if c["b"] == 38 and c["ngc"] == 3),
            key=lambda c: abs(c["cells"] - cells),
        )
        return best(c) / c["t"]["V0 ADAM map I8P"]

    groups = [
        ("7 M cells/GPU, no body", "n1-noib-wg-b38-np4.txt", 7.0),
        ("7 M cells/GPU, body", "n1-ib-wg-b38-np4.txt", 7.0),
        ("56 M cells/GPU, no body", "n1-noib-wg-l4b38-np4.txt", 56.0),
        ("56 M cells/GPU, body", "n1-ib-wg-l4b38-np4.txt", 56.0),
    ]
    parts = ("face fluxes", "ghost local", "BC + blend", "ghost MPI", "rest")
    colours = dict(
        zip(
            parts,
            (
                PAPER_PALETTE[0],
                PAPER_PALETTE[1],
                PAPER_PALETTE[4],
                PAPER_PALETTE[3],
                "0.65",
            ),
            strict=True,
        )
    )
    scen = ("today", "slab local fill", "slab fill + face kernel 2x")
    fig, axes = plt.subplots(1, 4, figsize=(11.0, 4.3), layout="constrained")
    for ax, (title, fname, cells) in zip(axes, groups, strict=True):
        ph = read_phase(phase_dir / fname)
        step = ph["step (loop)"]
        base = {p: ph[p] for p in parts[:-1]}
        base["rest"] = step - sum(base.values())
        r = ratio(cells)
        bars = [
            dict(base),
            dict(base, **{"ghost local": base["ghost local"] * r}),
            dict(
                base,
                **{
                    "ghost local": base["ghost local"] * r,
                    "face fluxes": base["face fluxes"] / 2,
                },
            ),
        ]
        for x, b in enumerate(bars):
            bottom = 0.0
            for p in parts:
                ax.bar(
                    x,
                    b[p],
                    bottom=bottom,
                    color=colours[p],
                    width=0.7,
                    label=p if x == 0 else None,
                )
                bottom += b[p]
            ghost = b["ghost local"] + b["BC + blend"] + b["ghost MPI"]
            ax.text(
                x,
                bottom,
                f"{bottom:.0f} ms\nghost {100 * ghost / bottom:.0f} %",
                ha="center",
                va="bottom",
                fontsize=7,
            )
        ax.set_xticks(range(3), [s.replace(" + ", "\n+ ") for s in scen], fontsize=7)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(0, step * 1.22)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("ms per step (mean over the 4 H100 ranks)")
    fig.legend(
        *axes[0].get_legend_handles_labels(),
        loc="outside lower center",
        ncol=5,
        fontsize=8,
        frameon=False,
    )
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def main() -> int:
    """Parse the arguments, read the results, draw the three figures in every requested theme."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, default=HERE / "results")
    parser.add_argument(
        "--phase",
        type=Path,
        default=REPO
        / "src/tests/flume/speedup/results/2026-11-25-spacehpc/phase/reports",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="light")
    args = parser.parse_args()
    runs: dict[str, list[Case]] = {}
    for gpu, pattern in (("H100", "h100-*.txt"), ("A30", "a30-*.txt")):
        files = sorted(args.results.glob(pattern))
        if not files:
            sys.exit(f"plot.py: no {pattern} in {args.results}")
        runs[gpu] = read_results(files[-1])
    args.out.mkdir(parents=True, exist_ok=True)
    for t in themes(args.theme):
        with theme(t):
            fig_floor(runs, args.out / "ghost-fill-floor.png")
            fig_speedup(runs, args.out / "ghost-fill-speedup.png")
            fig_step(args.phase, runs["H100"], args.out / "ghost-fill-step.png")
    print(
        f"plot.py: ghost-fill-{{floor,speedup,step}}.png in {args.out} ({args.theme})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
