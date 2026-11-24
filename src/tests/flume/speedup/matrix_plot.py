#!/usr/bin/env python3
"""Draw the summary figures of the FLUME speedup matrix from analyze.py's summary.csv (issue #133, phase 8.1).

Three figures, each `<prefix>-<name>.png` in --out (and its silver-ink twin `<prefix>-<name>-dark.png` with
--theme dark or both, docs/scripts/figure_theme.py):

    backend-speedup     T(cpu node) / T(GPU node) per case: one group per leg, one bar per role (uniform, frozen
                        AMR, runtime AMR); the CPU node is the CPU run with the fewest ranks, the GPU node the FNL run
                        with --gpu-node ranks
    amr-decomposition   per leg and backend node, the step cost of the three roles divided by the uniform case's:
                        the time inside the time steps and the time outside them (the regrid, its marker and its host
                        round trip), with the ideal N(role) / N(uniform) as a tick (the cell ratio: the step cost of
                        an AMR grid with the uniform grid's cost per cell)
    grind               the grind time g = step cost x devices / N against the cells per device, per backend: the
                        scaling series (filled: one node, open: several) and the AMR cases of the sphere

The time outside the steps is `step cost - steps`, both the medians of the window differences: on FNL it is the regrid
(the forest's other phases outside the steps cost nothing measurable there); on the CPU it also holds the wait of the
least loaded ranks, since both timers are maxima over the ranks (README, "Results").

Needs matplotlib (run it on WSL, e.g. exe/.regression-venv/bin/python); analyze.py stays stdlib.

Usage:
    matrix_plot.py <summary.csv> --out DIR [--prefix perf] [--theme light|dark|both] [--cpu-node 64] [--gpu-node 4]
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams["figure.facecolor"] = "#f3ede1"  # the washi ground of the docs (docs/BRAND.md)
mpl.rcParams["savefig.facecolor"] = "#f3ede1"
import matplotlib.pyplot as plt  # noqa: E402

DOCS_SCRIPTS = Path(__file__).resolve().parents[4] / "docs" / "scripts"
sys.path.insert(0, str(DOCS_SCRIPTS))
from figure_theme import PAPER_PALETTE, theme, themes  # noqa: E402

ROLES = ("uniform", "frozen", "runtime")
ROLE_LABEL = {"uniform": "uniform", "frozen": "frozen AMR grid", "runtime": "runtime AMR"}
ROLE_COLOUR = {"uniform": PAPER_PALETTE[0], "frozen": PAPER_PALETTE[4], "runtime": PAPER_PALETTE[2]}
LEG_LABEL = {1: "1 step\n2-D, WENO-7", 2: "2 sphere\n3-D, WENO-5, IB", 3: "3 channel\n3-D, hybrid, viscous"}
SERIES_COLOUR = {"sphere": PAPER_PALETTE[0], "step": PAPER_PALETTE[1]}
DPI = 130

Row = dict[str, str]


def read_rows(path: Path) -> list[Row]:
    """Return the rows of summary.csv.

    Parameters
    ----------
    path : Path
        analyze.py's summary.csv.

    Returns
    -------
    list[Row]
        One dict per window (backend, ranks, case).
    """
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def pick(rows: list[Row], leg: int, role: str, backend: str, np_: int) -> Row | None:
    """Return the row of a matrix case (a leg's role) on a backend and rank count, or None."""
    for r in rows:
        if (r["leg"] == str(leg) and r["role"] == role and r["backend"] == backend and r["np"] == str(np_)
                and r.get("kind", "matrix") == "matrix"):
            return r
    return None


@contextmanager
def figure(path: Path, size: tuple[float, float], ncols: int = 1,
           bottom: float = 0.0) -> Iterator[tuple[plt.Figure, list[plt.Axes]]]:
    """Yield a figure and its axes, then save it to `path` (the theme decides the file: figure_theme.theme); `bottom`
    is the fraction of the height kept free below the axes (a figure legend)."""
    fig, axes = plt.subplots(1, ncols, figsize=size, squeeze=False)
    yield fig, list(axes[0])
    fig.tight_layout(rect=(0.0, bottom, 1.0, 1.0))
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def backend_speedup(rows: list[Row], out: Path, cpu_np: int, gpu_np: int) -> None:
    """Draw T(cpu node) / T(GPU node) per leg and role."""
    with figure(out, (6.4, 3.6)) as (_, (ax,)):
        width = 0.26
        for k, role in enumerate(ROLES):
            xs, ys = [], []
            for leg in (1, 2, 3):
                c, g = pick(rows, leg, role, "cpu", cpu_np), pick(rows, leg, role, "fnl-nvf", gpu_np)
                if c is None or g is None:
                    continue
                xs.append(leg + (k - 1) * width)
                ys.append(float(c["step_cost_s"]) / float(g["step_cost_s"]))
            bars = ax.bar(xs, ys, width=width * 0.92, color=ROLE_COLOUR[role], label=ROLE_LABEL[role])
            ax.bar_label(bars, labels=[f"{y:.1f}" for y in ys], padding=2, fontsize=8, color="0.25")
        ax.set_xticks([1, 2, 3], [LEG_LABEL[leg] for leg in (1, 2, 3)], fontsize=8)
        ax.set_ylabel(f"T(CPU node, {cpu_np} ranks) / T(GPU node, {gpu_np} A30)", fontsize=8)
        ax.set_ylim(0, None)
        ax.margins(y=0.15)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        ax.set_title("GPU node against CPU node, per case", fontsize=9)


def amr_decomposition(rows: list[Row], out: Path, cpu_np: int, gpu_np: int) -> None:
    """Draw the step cost of each role over the uniform case's, split inside / outside the steps, per leg and node."""
    panels = (("fnl-nvf", gpu_np, f"GPU node ({gpu_np} A30)"), ("cpu", cpu_np, f"CPU node ({cpu_np} ranks)"))
    with figure(out, (9.6, 3.8), ncols=2) as (_, axes):
        width = 0.26
        for ax, (backend, np_, title) in zip(axes, panels, strict=True):
            for leg in (1, 2, 3):
                u = pick(rows, leg, "uniform", backend, np_)
                if u is None:
                    continue
                t_u, n_u = float(u["step_cost_s"]), float(u["cells_mean"])
                for k, role in enumerate(ROLES):
                    r = pick(rows, leg, role, backend, np_)
                    if r is None:
                        continue
                    x = leg + (k - 1) * width
                    total = float(r["step_cost_s"]) / t_u
                    inside = min(float(r["steps_cost_s"]) / t_u, total)
                    ax.bar(x, inside, width=width * 0.92, color=ROLE_COLOUR[role],
                           label=ROLE_LABEL[role] if leg == 1 else None)
                    ax.bar(x, total - inside, bottom=inside, width=width * 0.92, color=ROLE_COLOUR[role],
                           hatch="////", edgecolor="0.15", linewidth=0.0, alpha=0.55,
                           label="outside the steps (regrid)" if leg == 1 and role == "runtime" else None)
                    ideal = float(r["cells_mean"]) / n_u
                    ax.plot([x - width * 0.46, x + width * 0.46], [ideal, ideal], color="0.1", lw=1.6,
                            label="ideal N / N(uniform)" if leg == 1 and role == "uniform" else None)
                    ax.annotate(f"{total:.2f}", (x, total), textcoords="offset points", xytext=(0, 2),
                                ha="center", fontsize=7, color="0.25")
            ax.set_xticks([1, 2, 3], [LEG_LABEL[leg] for leg in (1, 2, 3)], fontsize=8)
            ax.set_ylabel("step cost / uniform step cost", fontsize=8)
            ax.set_ylim(0, 1.25)
            ax.grid(axis="y", alpha=0.3)
            ax.set_axisbelow(True)
            ax.set_title(title, fontsize=9)
        axes[0].legend(fontsize=7, frameon=False, loc="upper right", ncols=2)


def grind(rows: list[Row], out: Path, sizes: dict[str, int]) -> None:
    """Draw g = step cost x devices / N against the cells per device, per backend."""
    backends = [b for b in ("fnl-nvf", "cpu") if any(r["backend"] == b for r in rows)]
    unit = {"fnl-nvf": ("ns", 1.0), "cpu": ("us", 1e-3)}
    with figure(out, (9.6, 4.5), ncols=len(backends), bottom=0.17) as (fig, axes):
        for ax, backend in zip(axes, backends, strict=True):
            mine = [r for r in rows if r["backend"] == backend and r.get("cells_per_device")]
            name, scale = unit[backend]
            groups: list[tuple[str, list[Row], str, str]] = []
            for series in ("sphere", "step"):
                pts = [r for r in mine if r.get("series") == series]
                groups.append((f"{series}, uniform", pts, SERIES_COLOUR[series], "o"))
            amr = [r for r in mine if r["leg"] == "2" and r["role"] in ("frozen", "runtime")]
            groups.append(("sphere, frozen AMR", [r for r in amr if r["role"] == "frozen"], PAPER_PALETTE[4], "^"))
            groups.append(("sphere, runtime AMR", [r for r in amr if r["role"] == "runtime"], PAPER_PALETTE[2], "v"))
            for label, pts, colour, marker in groups:
                pts = sorted(pts, key=lambda r: float(r["cells_per_device"]))
                if not pts:
                    continue
                for multi in (False, True):
                    sel = [r for r in pts if (int(r["np"]) > sizes[backend]) == multi]
                    if not sel:
                        continue
                    ax.plot([float(r["cells_per_device"]) for r in sel],
                            [float(r["ns_per_cell_step"]) * scale for r in sel], marker, ls="none", ms=6,
                            color=colour, mfc="none" if multi else colour,
                            label=label if not multi else f"{label}, several nodes")
                for r in pts:
                    ax.annotate(r["np"], (float(r["cells_per_device"]), float(r["ns_per_cell_step"]) * scale),
                                textcoords="offset points", xytext=(4, 3), fontsize=6.5, color="0.35")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("cells per device (GPU or core)", fontsize=8)
            ax.set_ylabel(f"g, {name} per cell-step", fontsize=8)
            what = "A30" if backend == "fnl-nvf" else "core"
            ax.set_title(f"{backend}: grind time (labels: {what}s)", fontsize=9)
            ax.grid(True, which="both", alpha=0.3)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncols=4, fontsize=7, frameon=False)


def main() -> int:
    """Parse the command line, read summary.csv, write the figures."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", type=Path, help="analyze.py's summary.csv")
    parser.add_argument("--out", type=Path, required=True, help="the directory of the figures")
    parser.add_argument("--prefix", default="perf", help="file name prefix (default perf)")
    parser.add_argument("--theme", choices=("light", "dark", "both"), default="light",
                        help="light: <name>.png on washi; dark: <name>-dark.png on silver ink; both")
    parser.add_argument("--cpu-node", type=int, default=64, help="ranks of a CPU node (default 64)")
    parser.add_argument("--gpu-node", type=int, default=4, help="GPUs of a GPU node (default 4)")
    args = parser.parse_args()
    rows = read_rows(args.csv)
    args.out.mkdir(parents=True, exist_ok=True)
    sizes = {"cpu": args.cpu_node, "fnl-nvf": args.gpu_node}
    for t in themes(args.theme):
        with theme(t):
            backend_speedup(rows, args.out / f"{args.prefix}-backend-speedup.png", args.cpu_node, args.gpu_node)
            amr_decomposition(rows, args.out / f"{args.prefix}-amr-decomposition.png", args.cpu_node, args.gpu_node)
            grind(rows, args.out / f"{args.prefix}-grind.png", sizes)
    print(f"matrix_plot: {args.prefix}-{{backend-speedup,amr-decomposition,grind}}.png in {args.out} ({args.theme})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
