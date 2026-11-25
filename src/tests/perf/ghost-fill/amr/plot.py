#!/usr/bin/env python3
"""Draw the AMR ghost-fill benchmark figures (test_amr_ghost_bench.F90, README.md, 2026-11-25).

Reads every `results/h100-*.txt` and `results/a30-*.txt` of run.sh (in job order, a later run of a case replacing an
earlier one) and writes into --out, each with its silver-ink
twin `<name>-dark.png` under --theme both (docs/scripts/figure_theme.py). The files of the study (2026-11-25, the
kernel of b8641a59) draw the first three figures; the files of the #157 kernel (they have a V0p line: the density
guard once per fine ghost, with the b8641a59 kernel and the rows sorted by kind measured in the same run) draw the
fourth:

    amr-ghost-kinds.png     the fill on the production-like grid (253 blocks of 24^3, levels 2-3, 2.87 M rows) for every
                            prolongation x restriction pair, split by kind (copy, restriction, prolongation: the V0s
                            launches), with today's single launch (V0) marked; H100 and A30
    amr-ghost-row-cost.png  the cost of one row of each kind (ns per row, all 5 variables) on that grid: the copy, the
                            restriction per regime, the prolongation per regime, with and without the density guard
    amr-ghost-variants.png  the variants on that grid for the default pair (tricubic / mean) and the WENO pair: today
                            (V0), split by kind (V0s), descriptor kernels (D1, D2), guard off (V0ng) and, for the
                            prolongation, the per-anchor kernel (D3) against V0ng's prolongation
    amr-ghost-i157.png      #157 on that grid, every regime pair, H100 and A30, the density guard on and off: the fill
                            of b8641a59 (the study's V0 and V0ng), the #157 fill (V0: the guard once per fine ghost,
                            ADAM's row order; V0ng: the kernel without the guard) and the same kernel on the rows
                            sorted by kind (V0k, V0kng: action A1)

Needs matplotlib and numpy (on WSL: exe/.regression-venv/bin/python).

Usage:
    plot.py [--results DIR] --out DIR [--theme light|dark|both]
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
REPO = HERE.parents[4]
sys.path.insert(0, str(REPO / "docs" / "scripts"))
from figure_theme import PAPER_PALETTE, theme, themes  # noqa: E402

DPI = 130
FILLS = ("tricubic", "quintic", "quintic-mp", "weno6")
RESTRS = ("mean", "lagrange6", "weno6")
KINDS = ("copy", "restriction", "prolongation")
KCOL = {
    "copy": PAPER_PALETTE[0],
    "restriction": PAPER_PALETTE[4],
    "prolongation": PAPER_PALETTE[1],
}
RUN = re.compile(r"^## (\d+) (\d+) ([\d.]+) (\S+) (\S+)")
ROWS = re.compile(
    r"^# rows\s+(\d+) \(copy, restriction, prolongation:\s+(\d+)\s+(\d+)\s+(\d+)\)"
)
V0 = re.compile(r"^V0  .*median\s+([\d.]+) ms")
V0P = re.compile(r"^V0p .*median\s+([\d.]+) ms\s+bitwise vs V0 (\w)")
V0K = re.compile(
    r"^V0k .*median\s+([\d.]+) ms\s+bitwise vs V0 (\w)\s+\(guard off: V0kng\s+([\d.]+) ms\)"
)
SPLIT = re.compile(
    r"^(V0s|D1|D2) .*median\s+([\d.]+) ms\s+\(copy, restr, prol\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\)\s+bitwise (\w)"
)
NG = re.compile(r"^V0ng .*median\s+([\d.]+) ms\s+\(prolongation rows\s+([\d.]+)\)")
D3 = re.compile(r"^D3 .*median\s+([\d.]+) ms\s+bitwise vs V0ng prolongation (\w)")

Run = dict[str, object]


def read(path: Path) -> list[Run]:
    """Return the runs of one result file; stop on any variant that is not bitwise identical to today's fill."""
    runs: list[Run] = []
    for line in path.read_text().splitlines():
        if m := RUN.match(line):
            runs.append(
                {
                    "b": int(m.group(1)),
                    "base": int(m.group(2)),
                    "box": float(m.group(3)),
                    "fill": m.group(4),
                    "restr": m.group(5),
                }
            )
        elif not runs:
            continue
        elif m := ROWS.match(line):
            runs[-1]["rows"] = [int(m.group(k)) for k in (2, 3, 4)]
        elif m := V0.match(line):
            runs[-1]["V0"] = float(m.group(1))
        elif m := V0P.match(line):
            if m.group(2) != "T":
                sys.exit(f"plot.py: {path.name}: V0p not bitwise in {runs[-1]}")
            runs[-1]["V0p"] = float(m.group(1))
        elif m := V0K.match(line):
            if m.group(2) != "T":
                sys.exit(f"plot.py: {path.name}: V0k not bitwise in {runs[-1]}")
            runs[-1]["V0k"] = (float(m.group(1)), float(m.group(3)))
        elif m := SPLIT.match(line):
            if m.group(6) != "T":
                sys.exit(
                    f"plot.py: {path.name}: {m.group(1)} not bitwise in {runs[-1]}"
                )
            runs[-1][m.group(1)] = (
                float(m.group(2)),
                [float(m.group(k)) for k in (3, 4, 5)],
            )
        elif m := NG.match(line):
            runs[-1]["V0ng"] = (float(m.group(1)), float(m.group(2)))
        elif m := D3.match(line):
            if m.group(2) != "T":
                sys.exit(f"plot.py: {path.name}: D3 not bitwise in {runs[-1]}")
            runs[-1]["D3"] = float(m.group(1))
    return [r for r in runs if "V0" in r]


def prod(runs: list[Run], fill: str, restr: str) -> Run:
    """The production-like run (b 24, base 2, box 0.75) of one regime pair."""
    return next(
        r
        for r in runs
        if (r["b"], r["base"], r["box"], r["fill"], r["restr"])
        == (24, 2, 0.75, fill, restr)
    )


def fig_kinds(data: dict[str, list[Run]], out: Path) -> None:
    """Per regime pair, the V0s launches stacked by kind, today's V0 marked."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), layout="constrained")
    short = {
        "tricubic": "tric",
        "quintic": "quin",
        "quintic-mp": "q-mp",
        "weno6": "weno6",
        "mean": "mean",
        "lagrange6": "lag6",
    }
    labels = [f"{short[f]}\n{short[r]}" for f in FILLS for r in RESTRS]
    for ax, (gpu, runs) in zip(axes, data.items(), strict=True):
        for x, (f, r) in enumerate((f, r) for f in FILLS for r in RESTRS):
            run = prod(runs, f, r)
            bottom = 0.0
            for k, kind in enumerate(KINDS):
                h = run["V0s"][1][k]
                ax.bar(
                    x,
                    h,
                    bottom=bottom,
                    color=KCOL[kind],
                    width=0.75,
                    label=kind if x == 0 else None,
                )
                bottom += h
            ax.plot(
                x,
                run["V0"],
                marker="_",
                ms=16,
                mew=2.0,
                color="0.15",
                ls="none",
                label="today: one launch (V0)" if x == 0 else None,
            )
        ax.set_xticks(range(len(labels)), labels, fontsize=6.5)
        ax.set_title(f"{gpu}: 2:1 grid, 253 blocks of 24^3, 2.87 M rows", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.set_xlabel("prolongation / restriction regime")
    axes[0].set_ylabel("ms per fill (rows split by kind, V0s)")
    axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def fig_row_cost(data: dict[str, list[Run]], out: Path) -> None:
    """ns per row of each kind (5 variables) on the production-like grid, guard on (and off for the prolongation)."""
    fig, axes = plt.subplots(
        1, 2, figsize=(11.0, 4.2), sharey=False, layout="constrained"
    )
    for ax, (gpu, runs) in zip(axes, data.items(), strict=True):
        ref = prod(runs, "tricubic", "mean")
        rows = ref["rows"]
        names, vals, cols = ["copy"], [ref["V0s"][1][0] * 1e6 / rows[0]], [KCOL["copy"]]
        for r in RESTRS:
            run = prod(runs, "tricubic", r)
            names.append(f"restr.\n{r}")
            vals.append(run["V0s"][1][1] * 1e6 / rows[1])
            cols.append(KCOL["restriction"])
        for f in FILLS:
            run = prod(runs, f, "mean")
            names.append(f"prol.\n{f}")
            vals.append(run["V0s"][1][2] * 1e6 / rows[2])
            cols.append(KCOL["prolongation"])
        xs = range(len(names))
        ax.bar(xs, vals, color=cols, width=0.7)
        for x, f in zip(range(4, 8), FILLS, strict=True):
            run = prod(runs, f, "mean")
            if "V0ng" in run:
                ng = run["V0ng"][1] * 1e6 / rows[2]
                ax.plot(
                    x,
                    ng,
                    marker="o",
                    color="0.15",
                    ls="none",
                    ms=6,
                    label="prolongation, density guard off" if x == 4 else None,
                )
        for x, v in zip(xs, vals, strict=True):
            ax.text(x, v, f"{v:.1f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(list(xs), names, fontsize=7)
        ax.set_yscale("log")
        ax.set_title(f"{gpu}: cost of one row (5 variables)", fontsize=10)
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("ns per row")
    if any("V0ng" in r for r in data["H100"]):
        axes[0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def fig_variants(data: dict[str, list[Run]], out: Path) -> None:
    """The variants on the production-like grid, the default pair and the WENO pair, both GPUs."""
    pairs = (("tricubic", "mean"), ("weno6", "weno6"))
    variants = ("V0", "V0s", "D1", "D2", "V0ng")
    vcol = dict(
        zip(
            variants,
            (
                PAPER_PALETTE[6],
                PAPER_PALETTE[0],
                PAPER_PALETTE[3],
                PAPER_PALETTE[2],
                PAPER_PALETTE[5],
            ),
            strict=True,
        )
    )
    fig, axes = plt.subplots(1, 4, figsize=(12.0, 4.0), layout="constrained")
    for p, (gpu, runs) in enumerate(data.items()):
        for q, pair in enumerate(pairs):
            ax = axes[2 * p + q]
            run = prod(runs, *pair)
            vals = []
            for v in variants:
                if v == "V0":
                    vals.append(run["V0"])
                elif v == "V0ng":
                    vals.append(run["V0ng"][0] if "V0ng" in run else float("nan"))
                else:
                    vals.append(run[v][0] if v in run else float("nan"))
            ax.bar(
                range(len(variants)), vals, color=[vcol[v] for v in variants], width=0.7
            )
            for x, v in enumerate(vals):
                if v == v:
                    ax.text(x, v, f"{v:.2f}", ha="center", va="bottom", fontsize=7)
            if "D3" in run:
                ax.bar(
                    len(variants),
                    run["V0ng"][1],
                    color=PAPER_PALETTE[5],
                    width=0.7,
                    alpha=0.5,
                )
                ax.bar(len(variants) + 1, run["D3"], color=PAPER_PALETTE[1], width=0.7)
                ax.text(
                    len(variants),
                    run["V0ng"][1],
                    f"{run['V0ng'][1]:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                )
                ax.text(
                    len(variants) + 1,
                    run["D3"],
                    f"{run['D3']:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                )
                ticks = [*variants, "V0ng\nprol.", "D3\nprol."]
            else:
                ticks = list(variants)
            ax.set_xticks(range(len(ticks)), ticks, fontsize=7)
            ax.set_title(f"{gpu}, {pair[0]} / {pair[1]}", fontsize=9)
            ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("ms per fill")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def fig_i157(
    study: dict[str, list[Run]], after: dict[str, list[Run]], out: Path
) -> None:
    """Per regime pair, guard on (top) and off (bottom): b8641a59, #157 (ADAM's row order), the rows sorted by kind."""
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.6), layout="constrained")
    short = {
        "tricubic": "tric",
        "quintic": "quin",
        "quintic-mp": "q-mp",
        "weno6": "weno6",
        "mean": "mean",
        "lagrange6": "lag6",
    }
    labels = [f"{short[f]}\n{short[r]}" for f in FILLS for r in RESTRS]
    names = (
        ("b8641a59 (the study's run)", PAPER_PALETTE[2]),
        ("#157: guard once per ghost, ADAM's row order", PAPER_PALETTE[1]),
        ("#157 on the rows sorted by kind (A1)", PAPER_PALETTE[4]),
    )
    w = 0.27
    for col, gpu in enumerate(after):
        for row, guard in enumerate(("on", "off")):
            ax = axes[row, col]
            for x, (f, r) in enumerate((f, r) for f in FILLS for r in RESTRS):
                old, new = prod(study[gpu], f, r), prod(after[gpu], f, r)
                if guard == "on":
                    vals = (old["V0"], new["V0"], new["V0k"][0])
                else:
                    vals = (old["V0ng"][0], new["V0ng"][0], new["V0k"][1])
                for k, ((name, c), v) in enumerate(zip(names, vals, strict=True)):
                    ax.bar(
                        x + (k - 1) * w,
                        v,
                        width=w,
                        color=c,
                        label=name if x == 0 else None,
                    )
                ax.text(
                    x,
                    vals[1],
                    f"{100 * (vals[1] / vals[0] - 1):+.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=5.5,
                )
            ax.set_xticks(range(len(labels)), labels, fontsize=6.5)
            ax.set_title(f"{gpu}, density guard {guard}", fontsize=10)
            ax.grid(axis="y", alpha=0.3)
        axes[1, col].set_xlabel("prolongation / restriction regime")
    for row in range(2):
        axes[row, 0].set_ylabel("ms per fill (2:1 grid, 253 blocks of 24^3)")
    axes[0, 0].legend(fontsize=8, frameon=False, loc="upper left")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)


def main() -> int:
    """Read the newest result file of each GPU, draw the three figures in every requested theme."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results", type=Path, default=HERE / "results")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="light")
    args = parser.parse_args()
    data: dict[str, list[Run]] = {}
    after: dict[str, list[Run]] = {}
    for gpu, pattern in (("H100", "h100-*.txt"), ("A30", "a30-*.txt")):
        files = sorted(
            args.results.glob(pattern),
            key=lambda p: int(re.sub(r"\D", "", p.stem) or 0),
        )
        if not files:
            sys.exit(f"plot.py: no {pattern} in {args.results}")
        merged: dict[
            tuple, Run
        ] = {}  # every file, in job order: a later run of the same case replaces it
        merged_after: dict[tuple, Run] = {}
        for path in files:
            runs = read(path)
            if any("V0p" in run for run in runs):  # the #157 kernel
                for run in runs:
                    key = (run["b"], run["base"], run["box"], run["fill"], run["restr"])
                    if key in merged_after:
                        # the A30's guarded fill is bimodal between passes (README, "#157"): keep the lower times
                        for v in ("V0", "V0p"):
                            run[v] = min(run[v], merged_after[key][v])
                        run["V0k"] = tuple(
                            map(min, run["V0k"], merged_after[key]["V0k"])
                        )
                    merged_after[key] = run
                continue
            for run in runs:
                key = (run["b"], run["base"], run["box"], run["fill"], run["restr"])
                if key in merged:
                    # V0 on the A30 ran 3-34 % slow in a few passes (README, "Repeatability"): keep its lower time
                    run["V0"] = min(run["V0"], merged[key]["V0"])
                merged[key] = run
        data[gpu] = list(merged.values())
        if merged_after:
            after[gpu] = list(merged_after.values())
    args.out.mkdir(parents=True, exist_ok=True)
    for t in themes(args.theme):
        with theme(t):
            fig_kinds(data, args.out / "amr-ghost-kinds.png")
            fig_row_cost(data, args.out / "amr-ghost-row-cost.png")
            fig_variants(data, args.out / "amr-ghost-variants.png")
            if len(after) == 2:
                fig_i157(data, after, args.out / "amr-ghost-i157.png")
    print(
        f"plot.py: amr-ghost-{{kinds,row-cost,variants}}.png"
        f"{' and amr-ghost-i157.png' if len(after) == 2 else ''} in {args.out} ({args.theme})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
