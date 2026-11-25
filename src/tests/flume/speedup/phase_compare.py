#!/usr/bin/env python3
"""Rank the phases of a FLUME step by how much they grow from one node to several (issue #133, 2026-11-24).

Reads the `PHASE TIMING:` report that the measurement-only phase-timer patch prints (results/2026-11-24/phase/,
phase-timer.patch: ADAM_PHASE_TIMING=1) from run directories named `n<nodes>-<backend>-<ib|noib>-bar<1|0>`, each
holding `log.txt`, and writes one markdown table per backend, body and barrier setting: per phase the mean ms per step
over the ranks (and the max) on the smaller and the larger run, the growth in ms per step and the share of the step's
growth, ranked by growth. Equal cells per device on both runs (weak pairs), so a phase that scales perfectly does not
grow. Stdlib only.

Usage:
    phase_compare.py <dir with the run directories or their log.txt copies> [--out FILE]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

LINE = re.compile(r"PHASE TIMING: (.{14})\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+) %")
HEAD = re.compile(r"PHASE TIMING: ranks (\d+), steps timed (\d+), barriers skipped ([TF])")
NAME = re.compile(r"n(\d+)-(cpu|fnl)-(ib|noib)-bar([01])")
SUB = ("  pack", "  post", "  waitall", "  barrier 1", "  unpack", "  barrier 2")
TOTALS = ("step (loop)", "evolve")

Phase = dict[str, float]


def read(log: Path) -> tuple[dict[str, str], dict[str, Phase]] | None:
    """Return the header and the phases of a report, or None when the log has none.

    Parameters
    ----------
    log
        A run's log.txt.

    Returns
    -------
    tuple or None
        The header (ranks, steps, barriers skipped) and, per phase name, min, mean, max (ms per step), calls per step.
    """
    head: dict[str, str] = {}
    phases: dict[str, Phase] = {}
    for line in log.read_text(errors="replace").splitlines():
        if m := HEAD.search(line):
            head = {"ranks": m.group(1), "steps": m.group(2), "skipped": m.group(3)}
        elif m := LINE.search(line):
            phases[m.group(1).rstrip()] = {"min": float(m.group(2)), "mean": float(m.group(3)),
                                           "max": float(m.group(4)), "calls": float(m.group(5))}
    return (head, phases) if phases else None


def table(small: tuple[str, dict], large: tuple[str, dict]) -> list[str]:
    """Return the markdown lines comparing two runs of one backend, body and barrier setting."""
    (ns, (hs, ps)), (nl, (hl, pl)) = small, large
    step_growth = pl["step (loop)"]["mean"] - ps["step (loop)"]["mean"]
    rows = []
    for name, a in ps.items():
        b = pl.get(name)
        if b is None:
            continue
        rows.append((name, a, b, b["mean"] - a["mean"]))
    top = [r for r in rows if r[0] not in TOTALS and r[0] not in SUB]
    top.sort(key=lambda r: -r[3])
    sub = [r for r in rows if r[0] in SUB]
    out = [f"| Phase | {ns} node(s), {hs['ranks']} ranks: mean (max) ms/step | {nl} nodes, {hl['ranks']} ranks: mean (max) |"
           " growth ms/step | share of the step's growth | calls/step |",
           "|---|---|---|---|---|---|"]
    for name, a, b in [(r[0], r[1], r[2]) for r in rows if r[0] in TOTALS]:
        out.append(f"| **{name}** | {a['mean']:.2f} ({a['max']:.2f}) | {b['mean']:.2f} ({b['max']:.2f}) | "
                   f"{b['mean'] - a['mean']:+.2f} | | {b['calls']:.2f} |")
    for name, a, b, g in top:
        out.append(f"| {name} | {a['mean']:.2f} ({a['max']:.2f}) | {b['mean']:.2f} ({b['max']:.2f}) | {g:+.2f} | "
                   f"{100 * g / step_growth if step_growth else 0:.0f} % | {b['calls']:.2f} |")
        if name == "ghost MPI":
            for sname, sa, sb, sg in sub:
                out.append(f"| ... {sname.strip()} | {sa['mean']:.2f} ({sa['max']:.2f}) | {sb['mean']:.2f} "
                           f"({sb['max']:.2f}) | {sg:+.2f} | {100 * sg / step_growth if step_growth else 0:.0f} % | "
                           f"{sb['calls']:.2f} |")
    return out


def main() -> int:
    """Read every run, write the tables."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("root", type=Path, help="directory holding n<N>-<backend>-<ib>-bar<b>/log.txt")
    parser.add_argument("--out", type=Path, help="markdown output (default stdout)")
    args = parser.parse_args()
    runs: dict[tuple[str, str, str], dict[int, tuple]] = {}
    for log in sorted(args.root.glob("n*-*-*-bar*/log.txt")):
        m = NAME.fullmatch(log.parent.name)
        rep = read(log) if m else None
        if not m or rep is None:
            continue
        runs.setdefault((m.group(2), m.group(3), m.group(4)), {})[int(m.group(1))] = rep
    lines = ["# FLUME step, phase by phase, one node against several", ""]
    for (backend, ib, bar), by_n in sorted(runs.items()):
        if len(by_n) < 2:
            continue
        ns, nl = min(by_n), max(by_n)
        what = "with the immersed body" if ib == "ib" else "without the immersed body"
        lines += [f"## {backend}, {what}, ghost-exchange barriers {'on' if bar == '1' else 'off'}", ""]
        lines += table((str(ns), by_n[ns]), (str(nl), by_n[nl])) + [""]
    text = "\n".join(lines) + "\n"
    if args.out:
        args.out.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
