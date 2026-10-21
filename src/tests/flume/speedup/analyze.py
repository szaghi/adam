#!/usr/bin/env python3
"""Analyse the FLUME speedup matrix: the cost of a step of every measured case and leg, and the comparisons.

Why: run.sh records wall times only; this turns them into the cost of one step on the developed flow and into the
ratios the study asks for. Per measurement (``<case>/data/<leg>-np<N>/``):

- the cost of the NWIN-step window of repetition r is D_r = T_long,r - T_short,r (same restart, same start-up); a
  repetition counts only if both its runs exited 0. The step cost is median(D_r) / NWIN, the spread
  (max D_r - min D_r) / median D_r;
- the grid of the window comes from the long run of the first counted repetition: the blocks read from the restart
  (``field blocks number``, summed over the ranks) and the regrids after it (``flume: regrid at step S: blocks A -> B``,
  the grid of the steps after S); a step runs on the grid in force when it starts, so the mean over the window is
  the work it did. Cells = blocks * ni * nj * nk; the throughput is cells / step cost (cell-steps per second, each step
  SSP-33: three stages);
- the window's physical span and mean dt come from the per-step progress lines; the regrids inside the window and, on
  FNL, their host round trip seconds, from the regrid lines.

The comparisons (``results/summary.md``): backend speedup (a reference label against every other, per case), AMR
speedup (noamr against amr: the uniform grid at the finest level against dynamic AMR to the same level), orientation
(xz-oct and yz-oct against xy-oct, same cells), tree tax (xy-oct against xy-quad), and 2-D against 3-D (the cost of
a cell-step of the sphere against the cylinder).

Stdlib only, Python >= 3.7 (medooza's python3 has no numpy).

Usage: analyze.py [--ref LABEL]   LABEL default: the cpu label with the most ranks
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import make_inputs  # noqa: E402

RE_BLOCKS = re.compile(r"\[mpi-(\d+)\].*field blocks number\s+\+?(\d+)")
RE_REGRID = re.compile(r"flume: regrid at step (\d+): blocks (\d+) -> (\d+)(?:.*?, ([-+0-9.Ee]+) s \(host round trip\))?")
RE_STEP = re.compile(r"^\[mpi-0+\]t:\s+\+?(\d+)\s*$")
RE_DT = re.compile(r"^\[mpi-0+\]time step:\s+(\S+)\s*$")
RE_TIME = re.compile(r"^\[mpi-0+\]time:\s+(\S+)\s*$")


def read_info(path: Path) -> dict[str, str]:
    """Return the key=value pairs of a run.info."""
    info = {}
    if path.is_file():
        for line in path.read_text().splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                info[key] = value
    return info


def read_rows(path: Path) -> list[dict[str, str]]:
    """Return the rows of a timings.csv (empty if absent)."""
    if not path.is_file():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def ini_cells(ini: Path) -> int:
    """Return ni * nj * nk of an input."""
    dims = {}
    section = ""
    for line in ini.read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            section = line
        elif section == "[grid]" and "=" in line:
            key, value = (s.strip() for s in line.split("=", 1))
            if key in ("ni", "nj", "nk"):
                dims[key] = int(value)
    return dims["ni"] * dims["nj"] * dims["nk"]


def parse_log(log: Path) -> dict:
    """Return the restart blocks, the regrids and the per-step [dt, time] of a run log (rank 0 progress lines)."""
    blocks: dict[str, int] = {}
    regrids = []
    steps: dict[int, list[float]] = {}
    it = None
    for line in log.read_text(errors="replace").splitlines():
        m = RE_BLOCKS.search(line)
        if m:
            blocks.setdefault(m.group(1), int(m.group(2)))
            continue
        m = RE_REGRID.search(line)
        if m:
            regrids.append((int(m.group(1)), int(m.group(3)), float(m.group(4)) if m.group(4) else 0.0))
            continue
        m = RE_STEP.match(line)
        if m:
            it = int(m.group(1))
            steps.setdefault(it, [0.0, 0.0])
            continue
        m = RE_DT.match(line)
        if m and it is not None:
            steps[it][0] = float(m.group(1))
            continue
        m = RE_TIME.match(line)
        if m and it is not None:
            steps[it][1] = float(m.group(1))
    return {"blocks0": sum(blocks.values()), "regrids": sorted(regrids), "steps": steps}


def measure(d: Path) -> dict | None:
    """Return the metrics of one measurement directory, or None if no repetition counts."""
    info = read_info(d / "run.info")
    rows = read_rows(d / "timings.csv")
    runs = {(row["rep"], row["kind"]): float(row["seconds"]) for row in rows if row["status"] == "0"}
    reps = sorted({r for r, _ in runs if r != "0" and (r, "short") in runs and (r, "long") in runs}, key=int)
    if not reps:
        return None
    nwin, n1, it_w = int(info["nwin"]), int(info["n1"]), int(info["it_warm"])
    diffs = [runs[(r, "long")] - runs[(r, "short")] for r in reps]
    med = statistics.median(diffs)
    long_dir = d / f"rep{reps[0]}-long"
    log = parse_log(long_dir / "log.txt")
    first, last = it_w + n1 + 1, it_w + n1 + nwin
    total = 0
    for k in range(first, last + 1):
        grid = log["blocks0"]
        for s, after, _ in log["regrids"]:
            if s < k:
                grid = after
        total += grid
    in_window = [(s, sec) for s, _, sec in log["regrids"] if first <= s <= last]
    blocks_mean = total / nwin
    cells_mean = blocks_mean * ini_cells(long_dir / "input.ini")
    t_step = med / nwin
    steps = log["steps"]
    dts = [steps[k][0] for k in range(first, last + 1) if k in steps]
    warm = [float(r["seconds"]) for r in rows if r["kind"] == "warmup" and r["status"] == "0"]
    return {
        "label": d.name, "reps": len(reps), "diffs": diffs, "t_step": t_step,
        "spread": (max(diffs) - min(diffs)) / med if med > 0 else float("nan"),
        "blocks_mean": blocks_mean, "cells_mean": cells_mean,
        "throughput": cells_mean / t_step if t_step > 0 else 0.0,
        "regrids": len(in_window), "regrid_seconds": sum(sec for _, sec in in_window),
        "dt_mean": statistics.mean(dts) if dts else float("nan"),
        "t0": steps.get(first - 1, [0.0, float("nan")])[1], "t1": steps.get(last, [0.0, float("nan")])[1],
        "warmup_seconds": warm[0] if warm else float("nan"), "info": info,
    }


def result_md(case: str, m: dict) -> str:
    """Return the result.md of one measurement."""
    info = m["info"]
    n1, nwin, it_w = int(info["n1"]), int(info["nwin"]), int(info["it_warm"])
    gpu = f", GPU {info['gpu']}" if "gpu" in info else ""
    if m["throughput"] > 0:
        throughput = f"{m['throughput'] / 1e6:.3f} Mcell-steps/s ({1e9 / m['throughput']:.2f} ns per cell-step)"
    else:
        throughput = "n/a"
    regrids = str(m["regrids"])
    if m["regrid_seconds"] > 0:
        share = 100 * m["regrid_seconds"] / (m["t_step"] * nwin)
        regrids += f", {m['regrid_seconds']:.3f} s in host round trips ({share:.1f} % of the window)"
    lines = [
        f"# {case}: {m['label']}", "",
        f"Measured {info.get('date', '?')} on `{info.get('host', '?')}` (job {info.get('slurm_job') or '-'}), commit "
        f"`{info.get('commit', '?')[:12]}`, executable sha256 `{info.get('exe_sha256', '?')}`, launcher "
        f"`{info.get('launcher', '?')} {info.get('np', '?')}`{gpu}.", "",
        f"Warmup to t = {info.get('t_warm')} (step {it_w}, {m['warmup_seconds']:.1f} s). Window: steps "
        f"{it_w + n1 + 1}..{it_w + n1 + nwin}, t = {m['t0']:.6g}..{m['t1']:.6g}, mean dt {m['dt_mean']:.4g}.", "",
        "| Quantity | Value |", "|---|---|",
        f"| Step cost (median of {m['reps']} repetitions) | {m['t_step'] * 1e3:.3f} ms |",
        f"| Window costs D_r = T_long - T_short | {', '.join(f'{x:.3f}' for x in m['diffs'])} s |",
        f"| Spread (max - min) / median | {m['spread'] * 100:.1f} % |",
        f"| Mean blocks, cells | {m['blocks_mean']:.1f}, {m['cells_mean']:.4g} |",
        f"| Throughput | {throughput} |",
        f"| Regrids in the window | {regrids} |",
    ]
    if m["spread"] > 0.05:
        lines += ["", f"**Warning:** the spread is {m['spread'] * 100:.1f} % (> 5 %): a noisy node; rerun."]
    return "\n".join(lines) + "\n"


def ratio(a: float | None, b: float | None) -> str:
    """Return a / b formatted, or a dash."""
    return f"{a / b:.2f}" if a and b else "—"


def main() -> None:
    """Collect every measurement, write result.md files, results/summary.{csv,md}, and refresh the READMEs."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ref", default=None, help="reference label of the backend speedups (e.g. cpu-np32)")
    args = parser.parse_args()
    table: dict[tuple[str, str], dict] = {}
    for case in make_inputs.CASES:
        data = ROOT / case / "data"
        for d in sorted(data.glob("*-np*")) if data.is_dir() else []:
            m = measure(d)
            if m is None:
                print(f"analyze: {case}/data/{d.name}: no complete repetition, skipped")
                continue
            (d / "result.md").write_text(result_md(case, m))
            table[(case, d.name)] = m
    labels = sorted({label for _, label in table})
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    keys = ["t_step", "spread", "reps", "blocks_mean", "cells_mean", "throughput", "regrids", "regrid_seconds",
            "dt_mean", "t0", "t1", "warmup_seconds"]
    with (out / "summary.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "label"] + keys)
        for (case, label), m in sorted(table.items()):
            w.writerow([case, label] + [m[k] for k in keys])

    def t(case: str, label: str) -> float | None:
        return table[(case, label)]["t_step"] if (case, label) in table else None

    def ns(case: str, label: str) -> float | None:
        m = table.get((case, label))
        return 1e9 / m["throughput"] if m and m["throughput"] > 0 else None

    cpu = sorted((lb for lb in labels if lb.startswith("cpu-")), key=lambda s: int(s.rsplit("np", 1)[1]))
    ref = args.ref or (cpu[-1] if cpu else (labels[0] if labels else ""))
    head = "| " + " | ".join(labels) + " |"
    rule = "|---|" + "---|" * len(labels)
    md = [
        "# FLUME speedup matrix: results", "",
        "Generated by `analyze.py`; the method is in [../README.md](../README.md).", "",
        "## Step cost", "",
        "| Case | Label | Step [ms] | Spread | Mean cells | ns/cell-step | Regrids |", "|---|---|---|---|---|---|---|",
    ]
    for (case, label), m in sorted(table.items()):
        per_cell = f"{1e9 / m['throughput']:.2f}" if m["throughput"] > 0 else "—"
        md.append(f"| {case} | {label} | {m['t_step'] * 1e3:.3f} | {m['spread'] * 100:.1f} % | {m['cells_mean']:.4g} | "
                  f"{per_cell} | {m['regrids']} |")
    md += ["", f"## Backend speedup: step cost of `{ref}` / step cost of the label", "", "| Case " + head, rule]
    for case in make_inputs.CASES:
        if any((case, lb) in table for lb in labels):
            md.append(f"| {case} | " + " | ".join(ratio(t(case, ref), t(case, lb)) for lb in labels) + " |")
    md += ["", "## AMR speedup: step cost noamr / amr (same finest level)", "", "| Flow " + head, rule]
    for flow in make_inputs.FLOWS:
        md.append(f"| {flow} | " + " | ".join(ratio(t(f"{flow}-noamr", lb), t(f"{flow}-amr", lb)) for lb in labels)
                  + " |")
    md += [
        "", "## Orientation (same cells) and tree tax", "",
        "Step cost of xz-oct and yz-oct relative to xy-oct; xy-oct against xy-quad per step and per cell-step (the "
        "octree carries its null axis: more cells).", "",
        "| Variant | Label | xz-oct / xy-oct | yz-oct / xy-oct | xy-oct / xy-quad (step) | xy-oct / xy-quad (cell-step) |",
        "|---|---|---|---|---|---|",
    ]
    for amr in ("noamr", "amr"):
        xy, quad = f"xy-oct-shock-cylinder-{amr}", f"xy-quad-shock-cylinder-{amr}"
        for lb in labels:
            md.append(f"| {amr} | {lb} | {ratio(t(f'xz-oct-shock-cylinder-{amr}', lb), t(xy, lb))} | "
                      f"{ratio(t(f'yz-oct-shock-cylinder-{amr}', lb), t(xy, lb))} | {ratio(t(xy, lb), t(quad, lb))} | "
                      f"{ratio(ns(xy, lb), ns(quad, lb))} |")
    md += ["", "## 2-D against 3-D: cost of a cell-step, sphere / cylinder", "",
           "| Variant | Label | sphere / xy-quad | sphere / xy-oct |", "|---|---|---|---|"]
    for amr in ("noamr", "amr"):
        sph = f"xyz-shock-sphere-{amr}"
        for lb in labels:
            md.append(f"| {amr} | {lb} | {ratio(ns(sph, lb), ns(f'xy-quad-shock-cylinder-{amr}', lb))} | "
                      f"{ratio(ns(sph, lb), ns(f'xy-oct-shock-cylinder-{amr}', lb))} |")
    full = []
    for case in make_inputs.CASES:
        fdir = ROOT / case / "full"
        for d in sorted(fdir.glob("*-np*")) if fdir.is_dir() else []:
            steps = parse_log(d / "log.txt")["steps"] if (d / "log.txt").is_file() else {}
            for row in read_rows(d / "timings.csv"):
                full.append(f"| {case} | {d.name} | {row['seconds']} | {max(steps, default=0)} | {row['status']} |")
    if full:
        md += ["", "## Full simulations", "", "| Case | Label | Wall [s] | Steps | Exit |", "|---|---|---|---|---|"]
        md += full
    (out / "summary.md").write_text("\n".join(md) + "\n")
    conf = make_inputs.load_conf()
    for case in make_inputs.CASES:
        if (ROOT / case).is_dir():
            (ROOT / case / "README.md").write_text(make_inputs.case_readme(case, conf))
    print(f"analyze: {len(table)} measurements -> {out / 'summary.md'}")


if __name__ == "__main__":
    main()
