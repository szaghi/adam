#!/usr/bin/env python3
"""Draw the PRISM figures of the ADAM paper in the publication style, with their manifest (issue #142 Part D, #145).

Why: the PRISM figures of the paper were digitized from animation frames of runs that predate #104 and #132. These are
drawn from runs of the current tree, in the one publication style of the paper (docs/scripts/figure_theme.py
`paper()`, the block src/tests/flume/verification/paper_style.py applies to the FLUME figures):

* prism-pml-reflection: a Gaussian pulse at a Neumann boundary and at a 20-cell classical PML backed by a perfect
  conductor, incident and reflected B_z(x) (the runs `neumann` and `pml-n20` of pml-reflection/sweep.sh);
* prism-pml-sweep: the reflected amplitude A_r/A_0 against the layer thickness N and against gamma_max at N = 20, with
  the design reflection R = exp(-2 gamma_max L / ((n + 1) c)) (continuum) and its sampled form (the profile as the
  code samples it, gamma_max (1 - m/(N - 1))^n at the cell m from the wall), the CFS layer, the grid-scale wave
  Silver-Mueller and Neumann both send back, and the noise floor (all the runs of sweep.sh).

The rotating-field figures of the paper (prism-rmf-b-phases, prism-rmf-comparison, prism-rmf-divergence) are not
drawn: their case (a 25.6 x 25.6 x 32 cm PEC cavity, 1 mm cells, 99.5 A) is the co-author's and has no input in the
repository (the nearest, src/tests/prism/fnl/fd/BC_tests/Rectangular_array_PEC, is a 38.4 x 38.4 x 48 cm cavity on 2 mm
cells): the manifest lists them as skipped.

Usage:
    make_paper_figures.py --pml-runs <sweep out> --out <dir> [--fnl-runs <sweep out>] [--provenance <json>]
                          [--docs <docs/public/prism>]

--docs also writes the sweep figure for the documentation (pml-sweep.png on the washi ground of the docs, and its
silver-ink twin pml-sweep-dark.png; docs/scripts/figure_theme.py), outside the publication style.

--provenance is a JSON object merged into the manifest (commit, toolchain, device, date, jobs of the runs).
Writes <out>/<name>.pdf and .png at the printed width, and <out>/figures.json.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "docs" / "scripts"))
sys.path.insert(0, str(REPO / "src" / "tests" / "flume" / "verification"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
from figure_theme import (  # noqa: E402
    PAGE_INCLUDES,
    PAGE_WIDTHS,
    PAPER_DPI_PDF,
    PAPER_DPI_PNG,
    PAPER_MARKERS,
    PAPER_PALETTE,
    PAPER_PANEL_LABEL_PT,
    paper,
)
from figure_theme import theme  # noqa: E402
from paper_style import PANEL_LABEL, audit, script_provenance  # noqa: E402

C0 = 299792458.0
DX = 1.0e-3


def design_r(gamma_max: float, cells: int, n: float = 2.0) -> float:
    """The continuum design reflection of a layer of `cells` cells: exp(-2 gamma_max L / ((n + 1) c)), L = cells dx."""
    return math.exp(-2.0 * gamma_max * cells * DX / ((n + 1.0) * C0))


def sampled_r(gamma_max: float, cells: int, n: float = 2.0) -> float:
    """The design reflection of the profile as the code samples it: gamma_max (1 - m/(N - 1))^n at the cell m = 0 ..
    N - 1 from the wall (adam_prism_pml_object.F90 compute_gamma_profile: the outer cell has depth 1, the inner 0)."""
    m = np.arange(cells)
    return math.exp(
        -2.0 * DX * gamma_max * float(((1.0 - m / (cells - 1)) ** n).sum()) / C0
    )


def label(ax, k: int) -> None:  # noqa: ANN001
    """Panel label (a), (b), ... in the left title slot."""
    ax.set_title(PANEL_LABEL % "abcdefgh"[k], loc="left", fontsize=PAPER_PANEL_LABEL_PT)


def save(
    fig,
    out: Path,
    name: str,
    frac: float,
    entries: list[dict],
    runs: list[Path],
    note: str,
) -> None:  # noqa: ANN001
    """Write PDF and PNG at the printed size (frac of \\textwidth), and the manifest entry."""
    left = audit(fig)
    for ext, dpi in (("pdf", PAPER_DPI_PDF), ("png", PAPER_DPI_PNG)):
        fig.savefig(out / f"{name}.{ext}", dpi=dpi)
    w, h = fig.get_size_inches()
    plt.close(fig)
    include = PAGE_INCLUDES["text"] % f"figures/{name}"
    if frac != 1.0:
        include = include.replace(r"width=\textwidth", rf"width={frac:g}\textwidth")
    entries.append(
        {
            "name": name,
            "size_in": [round(float(w), 3), round(float(h), 3)],
            "latex": include,
            "files": [f"{name}.pdf", f"{name}.png"],
            "runs": [str(r) for r in runs],
            "note": note,
            "audit": left,
        }
    )


def fig_reflection(pml: Path, out: Path, entries: list[dict]) -> None:
    """Neumann against the 20-cell PML: incident and reflected B_z(x)."""
    cases = (
        ("neumann", "Neumann boundary"),
        ("pml-n20", "20-cell PML, perfect conductor behind"),
    )
    frac = 0.8
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(frac * PAGE_WIDTHS["text"], 3.0),
        sharex=True,
        constrained_layout=True,
    )
    for row, (case, _title) in enumerate(cases):
        d = np.load(pml / case / "profiles.npz")
        x, t, bz = d["x"] * 100.0, d["t"], d["bz"]
        wall = 0.08 if case == "neumann" else 0.08 + 20 * DX
        t_inc = 0.05 / C0  # incident: pulse centred at x = 5 cm
        # reflected: the PML image pulse centred at x = 0; at the Neumann face, the grid-scale packet it sends back
        # (2.2 c) centred near x = 0, before it reaches the opposite face
        t_ref = (2.0 * wall) / C0 if case != "neumann" else 0.12 / C0
        ki = int(np.argmin(np.abs(t - t_inc)))
        kr = int(np.argmin(np.abs(t - t_ref)))
        a, b = axes[row]
        a.plot(x, bz[ki], color=PAPER_PALETTE[0])
        a.set_ylabel(r"$B_z$ [T]")
        b.plot(x, bz[kr], color=PAPER_PALETTE[1], label="PRISM")
        if case == "pml-n20":
            r = sampled_r(7.8e10, 20)
            b.plot(
                x,
                r * bz[0][::-1],
                color="black",
                ls="--",
                lw=0.7,
                label=rf"incident $\times\, R = {r:.3f}$",
            )
            b.legend(loc="center right", handlelength=1.4, fontsize=7)
        for ax in (a, b):
            if case != "neumann":  # the layers, [0.08, 0.10] m on both sides
                ax.axvspan(8.0, x.max() + 0.05, color="0.88", lw=0)
                ax.axvspan(x.min() - 0.05, -8.0, color="0.88", lw=0)
            ax.set_xlim(x.min() - 0.05, x.max() + 0.05)
            ax.grid(True)
        a.text(
            0.03, 0.92, f"$t = {t[ki] * 1e9:.2f}$ ns", transform=a.transAxes, va="top"
        )
        peak = float(np.abs(bz[kr]).max()) / float(np.abs(bz[0]).max())
        b.text(
            0.03,
            0.92,
            f"$t = {t[kr] * 1e9:.2f}$ ns\npeak {peak:.2g} $A_0$",
            transform=b.transAxes,
            va="top",
        )
        label(a, 2 * row)
        label(b, 2 * row + 1)
    for ax in axes[1]:
        ax.set_xlabel("$x$ [cm]")
    save(
        fig,
        out,
        "prism-pml-reflection",
        frac,
        entries,
        [pml / c for c, _ in cases],
        "rows: Neumann (a, b: the grid-scale wave the face sends back, t = 0.12 m/c), 20-cell classical PML "
        "gamma_max 7.8e10, n 2 (c, d); grey: the layers; dashed: the incident pulse mirrored and scaled by the "
        "sampled design R",
    )


def fig_sweep(
    pml: Path,
    out: Path,
    entries: list[dict],
    fnl: Path | None,
    docs_name: str | None = None,
) -> None:
    """A_r/A_0 against N and against gamma_max, with R(L). With `docs_name` the PNG goes to `out` under that name (the
    docs figure: no PDF, no manifest entry)."""
    res = json.loads((pml / "results.json").read_text())
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(PAGE_WIDTHS["text"], 3.3), constrained_layout=True
    )
    floor = res["far"]["ratio"]
    # (a) thickness at gamma_max = 7.8e10
    ns = [10, 20, 40, 80]
    meas = [res[f"pml-n{n}"]["ratio"] for n in ns]
    nn = np.linspace(8, 84, 200)
    a.semilogy(
        nn,
        [design_r(7.8e10, v) for v in nn],
        color="black",
        lw=0.8,
        label="design $R(L)$",
    )
    a.semilogy(
        ns,
        [sampled_r(7.8e10, v) for v in ns],
        ls="none",
        marker="x",
        color="black",
        ms=4,
        label="design, sampled profile",
    )
    a.semilogy(
        ns,
        meas,
        ls="-",
        marker=PAPER_MARKERS[0],
        color=PAPER_PALETTE[0],
        label=r"PML, $\gamma_{\max} = 7.8\times10^{10}$ s$^{-1}$",
    )
    fixed = {
        10: "pml-n10-g6.24e11",
        20: "pml-n20-g3.12e11",
        40: "pml-n40-g1.56e11",
        80: "pml-n80",
    }
    a.semilogy(
        list(fixed),
        [res[v]["ratio"] for v in fixed.values()],
        ls="-",
        marker=PAPER_MARKERS[1],
        color=PAPER_PALETTE[1],
        label=r"PML, design $R = 9.4\times10^{-7}$",
    )
    a.semilogy(
        20,
        res["cfs-n20-k5"]["ratio"],
        ls="none",
        marker=PAPER_MARKERS[2],
        color=PAPER_PALETTE[2],
        label=r"CFS, $\kappa_{\max} = 5$, $\alpha_{\max} = 2.5\times10^8$ s$^{-1}$",
    )
    if fnl is not None and (fnl / "pml-n20" / "result.json").is_file():
        a.semilogy(
            20,
            json.loads((fnl / "pml-n20" / "result.json").read_text())["ratio"],
            ls="none",
            marker="o",
            mfc="none",
            ms=7,
            color=PAPER_PALETTE[0],
            label="PML, GPU",
        )
    # Silver-Mueller and Neumann send back the same grid-scale wave (the Neumann measure of sweep.sh is the pulse the
    # opposite Neumann face re-emits, not a reflection: pml-reflection.md), so one line for the two
    a.axhline(
        res["silver-muller"]["ratio"],
        color=PAPER_PALETTE[3],
        lw=0.8,
        ls="-.",
        label="Silver–Müller or Neumann, no layer",
    )
    a.axhline(floor, color="0.5", lw=0.8, ls=":", label="noise floor")
    a.set_xlabel("layer thickness $N$ [cells]")
    a.set_ylabel(r"$A_r / A_0$")
    a.set_ylim(3e-7, 2.0)
    a.grid(True, which="major")
    label(a, 0)
    # (b) gamma_max at N = 20
    gs = {
        3.9e10: "pml-n20-g3.9e10",
        7.8e10: "pml-n20",
        1.56e11: "pml-n20-g1.56e11",
        3.12e11: "pml-n20-g3.12e11",
        6.24e11: "pml-n20-g6.24e11",
    }
    gg = np.geomspace(3.0e10, 8.0e11, 200)
    b.loglog(gg, [design_r(v, 20) for v in gg], color="black", lw=0.8)
    b.loglog(
        list(gs),
        [sampled_r(v, 20) for v in gs],
        ls="none",
        marker="x",
        color="black",
        ms=4,
    )
    b.loglog(
        list(gs),
        [res[v]["ratio"] for v in gs.values()],
        ls="-",
        marker=PAPER_MARKERS[0],
        color=PAPER_PALETTE[0],
    )
    b.text(0.97, 0.95, "$N = 20$", transform=b.transAxes, ha="right", va="top")
    b.axhline(floor, color="0.5", lw=0.8, ls=":")
    b.set_xlabel(r"$\gamma_{\max}$ [s$^{-1}$]")
    b.set_ylim(3e-7, 2.0)
    b.grid(True, which="major")
    b.set_xticks([5e10, 1e11, 2e11, 5e11])
    b.set_xticklabels(
        [r"$5\times10^{10}$", r"$10^{11}$", r"$2\times10^{11}$", r"$5\times10^{11}$"]
    )
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    label(b, 1)
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=7)
    runs = sorted(p.parent for p in pml.glob("*/result.json"))
    if docs_name is not None:
        fig.savefig(out / f"{docs_name}.png", dpi=150)
        plt.close(fig)
        return
    save(
        fig,
        out,
        "prism-pml-sweep",
        1.0,
        entries,
        runs,
        "normal incidence, n = 2; A_r/A_0 = peak |B_z| in [-4, 4] cm after the incident pulse left it, over the "
        "initial peak; floor = run without a boundary within reach",
    )


def main() -> int:
    """Draw the figures and the manifest."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pml-runs", type=Path, required=True)
    parser.add_argument("--fnl-runs", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, default=None)
    parser.add_argument("--docs", type=Path, default=None)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    entries: list[dict] = []
    with paper():
        fig_reflection(args.pml_runs, args.out, entries)
        fig_sweep(args.pml_runs, args.out, entries, args.fnl_runs)
    doc = {
        "generated": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": str(Path(__file__)),
        "style": "paper (docs/scripts/figure_theme.py paper())",
        **script_provenance(REPO),
        "matplotlib": matplotlib.__version__,
        "figures": entries,
        "skipped": {
            n: "the co-author's rotating-field cavity has no input in the repository"
            for n in (
                "prism-rmf-b-phases",
                "prism-rmf-comparison",
                "prism-rmf-divergence",
            )
        },
    }
    if args.docs:
        args.docs.mkdir(parents=True, exist_ok=True)
        for name in ("light", "dark"):
            with (
                matplotlib.rc_context(
                    {
                        "figure.facecolor": "#f3ede1",
                        "savefig.facecolor": "#f3ede1",
                        "font.size": 8.0,
                        "legend.fontsize": 6.5,
                    }
                ),
                theme(name),
            ):
                docs_entries: list[dict] = []
                fig_sweep(
                    args.pml_runs,
                    args.docs,
                    docs_entries,
                    args.fnl_runs,
                    docs_name="pml-sweep",
                )
    if args.provenance:
        doc["runs_provenance"] = json.loads(args.provenance.read_text())
    (args.out / "figures.json").write_text(json.dumps(doc, indent=2) + "\n")
    for e in entries:
        print(f"{e['name']}: {', '.join(e['files'])}; audit: {e['audit'] or 'clean'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
