#!/usr/bin/env python3
"""Turn the windows of the FLUME speedup matrix into step costs and ratios (issue #133 phase 3, stdlib only).

Reads every window run that run.sh wrote, `<work>/windows/<backend>-np<N>/<case>/r<rep>-n<steps>/` (run.info and
log.txt), and writes `<out>/summary.md` and `<out>/summary.csv`.

Per window (a backend, a rank count, a case), from the forest timing lines that ADAM_SEAM_TIMING=1 prints:

    step cost   median over the repetitions of (T_loop(long) - T_loop(short)) / (long - short): the start-up, the
                restart read and the first-step transients cancel in the difference; the same difference of the steps
                timer, of the regrid timer and of the seam fills (summed over the realms) splits it
    spread      (max - min) / median of the repetitions' differences; above run.info's spread_warning (2 %) the window
                is flagged
    N           the mean cell count over the steps the difference keeps (short, long], from the blocks per realm of the
                window's dump (run.info realm_blocks, cells per block block_cells) and the `flume: regrid at step S:
                blocks A -> B` lines of the long runs (a regrid at step S changes the grid of the steps after S). The
                lines do not name their realm: they are printed in realm order, so each is given to the first realm,
                from the one after the previous line of the same step, whose current count is A
    ns/cell-step  step cost * ranks / N * 1e9 (ranks = devices: cores of a CPU run, GPUs of an FNL run)

Scaling and grind time (phase 9), per window and per series (`SERIES` of the case's
`<work>/cases/<case>/case.conf`, written by make_cases.py from matrix.conf: the uniform case of a leg with a `series`
key and the weak cases of that leg):

    g                  the grind time, step cost * devices / N in ns per cell-step (the column above); cells per
                       device = N / devices
    strong efficiency  per backend and case, g(smallest rank count present) / g(ranks): 1 is ideal
    weak efficiency    per backend and series, the windows grouped by cells per device (a group: loads within
                       LOAD_TOLERANCE, 10 %, of its largest), g(the group's smallest device count) / g(devices): 1 is
                       ideal; a group of one window has none
    The AMR ratios below take the matrix cases only (KIND=matrix): a weak case is a uniform case of the same leg.

The ratios, per leg, backend and rank count, by role (U the uniform case, R the runtime case, F the frozen one; leg 1
calls them a, c, b, legs 2 and 3 a, b, c):

    AMR speedup          T(U) / T(R), its ideal bound N(U) / N(R)
    AMR grid efficiency  [T(U) N(F) / N(U)] / T(F): the frozen AMR grid against the uniform one at constant cell work
                         (ghosts, 2:1 seam fills, reflux); 1 is ideal
    regrid cost          T(R) / T(F) - 1 (marking, regrid, prolongation; the FNL host round trip)
    backend speedup      per case, T(cpu at its largest rank count, the node) / T(backend, ranks)

Runs under `python3 -I` (no third-party package: medooza's python3 has none).

Usage:
    analyze.py [--work DIR] [--out DIR]
"""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import statistics
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOOP = re.compile(r"forest timing \(max over ranks, s\): loop\s+([\d.]+), steps\s+([\d.]+)")
REGRID = re.compile(r"forest regrid timing \(max over ranks, s\): regrid\s+([\d.]+)")
SEAM = re.compile(r"forest seam timing realm (\d+) .*?total\s+([\d.]+)")
REGRID_LINE = re.compile(r"^flume: regrid at step (\d+): blocks (\d+) -> (\d+)", re.M)
ROLE_TAG = {"uniform": "U", "runtime": "R", "frozen": "F"}


class AnalyzeError(Exception):
    """A run directory that cannot be read."""


@dataclass
class Run:
    """One window run: its run.info and the numbers of its log."""

    path: Path
    info: dict[str, str]
    loop: float
    steps: float
    regrid: float
    seam: float
    regrids: list[tuple[int, int, int]]  # (step, blocks before, blocks after), in log order


@dataclass
class Window:
    """The runs of one backend, rank count and case."""

    backend: str
    np: int
    case: str
    leg: int
    role: str
    short: int
    long: int
    spread_warning: float
    runs: dict[tuple[int, int], Run] = field(default_factory=dict)  # (rep, steps) -> run
    notes: list[str] = field(default_factory=list)


def read_info(path: Path) -> dict[str, str]:
    """Return the key=value lines of a run.info."""
    info = {}
    for line in path.read_text().splitlines():
        key, sep, value = line.partition("=")
        if sep:
            info[key] = value
    return info


def read_run(directory: Path) -> Run:
    """Read one run directory (run.info and log.txt)."""
    info = read_info(directory / "run.info")
    log = directory / "log.txt"
    if not log.is_file():
        raise AnalyzeError(f"{directory}: no log.txt")
    text = log.read_text(errors="replace")
    m = LOOP.search(text)
    if not m:
        raise AnalyzeError(f"{directory}: no forest timing line (ADAM_SEAM_TIMING=1 unset, or the run died)")
    r = REGRID.search(text)
    seam = sum(float(s.group(2)) for s in SEAM.finditer(text))
    regrids = [(int(a), int(b), int(c)) for a, b, c in REGRID_LINE.findall(text)]
    return Run(directory, info, float(m.group(1)), float(m.group(2)), float(r.group(1)) if r else 0.0, seam, regrids)


def realm_assignments(blocks: list[int], regrids: list[tuple[int, int, int]], cap: int = 256) -> list[list[int]]:
    """Return the ways of giving each regrid line to a realm that agree with the block counts.

    The lines of one step come in realm order (`forest_object%regrid_realms` ends the realms in order, and a realm whose
    grid did not change prints nothing), so within a step the realms are strictly increasing; a line `A -> B` goes to
    a realm whose current count is A. Several assignments survive when realms share a count; at most `cap` are kept.

    Parameters
    ----------
    blocks : list[int]
        The block count of each realm at the start of the run.
    regrids : list[tuple[int, int, int]]
        The regrid lines (step, blocks before, blocks after), in log order.
    cap : int
        The most assignments returned.

    Returns
    -------
    list[list[int]]
        Each assignment: the realm index of every line.
    """
    found: list[list[int]] = []

    def walk(k: int, current: list[int], last_step: int | None, last_realm: int, chosen: list[int]) -> None:
        if len(found) >= cap:
            return
        if k == len(regrids):
            found.append(list(chosen))
            return
        step, before, after = regrids[k]
        first = last_realm + 1 if step == last_step else 0
        for realm in range(first, len(current)):
            if current[realm] == before:
                saved = current[realm]
                current[realm] = after
                chosen.append(realm)
                walk(k + 1, current, step, realm, chosen)
                chosen.pop()
                current[realm] = saved

    walk(0, list(blocks), None, -1, [])
    return found


def mean_cells(run: Run, short: int, notes: list[str]) -> float:
    """Return the mean cell count of a long run over its steps (s0 + short, s0 + long].

    Parameters
    ----------
    run : Run
        The long run of a window.
    short : int
        The steps of the short run (the first steps, left out of the difference).
    notes : list[str]
        Collects a note when the regrid lines cannot be given to the realms, or not in a single way that matters.

    Returns
    -------
    float
        The step-weighted mean number of cells (the first assignment of the regrid lines to the realms).
    """
    blocks = [int(x) for x in run.info["realm_blocks"].split(",")]
    cells = [int(x) for x in run.info["block_cells"].split(",")]
    s0, steps = int(run.info["start_step"]), int(run.info["steps"])

    def mean_for(assignment: list[int]) -> float:
        total, k_event, current = 0.0, 0, list(blocks)
        for k in range(s0 + 1, s0 + steps + 1):
            while k_event < len(run.regrids) and run.regrids[k_event][0] <= k - 1:  # a regrid at S acts from S + 1
                current[assignment[k_event]] = run.regrids[k_event][2]
                k_event += 1
            if k > s0 + short:
                total += sum(b * c for b, c in zip(current, cells, strict=True))
        return total / (steps - short)

    assignments = realm_assignments(blocks, run.regrids)
    if not assignments:
        notes.append(f"{run.path.parent.name}/{run.path.name}: the regrid lines do not match the block counts "
                     f"{blocks}; N ignores the regrids")
        return float(sum(b * c for b, c in zip(blocks, cells, strict=True)))
    means = [mean_for(a) for a in assignments]
    if max(means) - min(means) > 1.0e-9 * max(means):
        notes.append(f"{run.path.parent.name}/{run.path.name}: realms with equal block counts make N ambiguous, "
                     f"{min(means):,.0f} to {max(means):,.0f} cells (the first kept)")
    return means[0]


def collect(work: Path) -> tuple[dict[tuple[str, int, str], Window], list[str]]:
    """Read every window run under `<work>/windows`; return the windows and the problems met."""
    windows: dict[tuple[str, int, str], Window] = {}
    problems = []
    for info_path in sorted(work.glob("windows/*/*/r*-n*/run.info")):
        try:
            run = read_run(info_path.parent)
        except (AnalyzeError, OSError) as err:
            problems.append(str(err))
            continue
        i = run.info
        key = (i["backend"], int(i["np"]), i["case"])
        if key not in windows:
            windows[key] = Window(i["backend"], int(i["np"]), i["case"], int(i["leg"]), i["role"],
                                  int(i["window_short"]), int(i["window_long"]), float(i["spread_warning"]))
        windows[key].runs[(int(i["rep"]), int(i["steps"]))] = run
    return windows, problems


def read_case_conf(work: Path, case: str) -> dict[str, str]:
    """Return the keys of `<work>/cases/<case>/case.conf` (empty when absent: a case no longer in matrix.conf).

    Parameters
    ----------
    work : Path
        run.sh's run tree.
    case : str
        The case name.

    Returns
    -------
    dict[str, str]
        The shell assignments of the file, unquoted.
    """
    path = work / "cases" / case / "case.conf"
    keys: dict[str, str] = {}
    if not path.is_file():
        return keys
    for line in path.read_text().splitlines():
        name, sep, value = line.partition("=")
        if sep and re.fullmatch(r"[A-Z_0-9]+", name):
            words = shlex.split(value)
            keys[name] = words[0] if words else ""
    return keys


def summarise(w: Window) -> dict[str, object] | None:
    """Return the numbers of one window, or None when no repetition has both runs."""
    reps = sorted({rep for rep, _ in w.runs})
    pairs = [(w.runs[(r, w.short)], w.runs[(r, w.long)]) for r in reps if (r, w.short) in w.runs and (r, w.long) in w.runs]
    missing = len(reps) - len(pairs)
    if missing:
        w.notes.append(f"{missing} repetition(s) without both runs")
    if not pairs:
        return None
    dn = w.long - w.short
    loop = [(b.loop - a.loop) / dn for a, b in pairs]
    cost = statistics.median(loop)
    spread = (max(loop) - min(loop)) / cost if cost > 0 else float("inf")
    cells = statistics.mean(mean_cells(b, w.short, w.notes) for _, b in pairs)
    ranks = w.np
    if any(b.info.get("src_modified") == "yes" for _, b in pairs):
        w.notes.append("built from sources that differ from the recorded commit (src_modified=yes)")
    if len({b.info.get("exe_sha256") for a, b in pairs} | {a.info.get("exe_sha256") for a, b in pairs}) > 1:
        w.notes.append("the runs used different executables")
    if spread > w.spread_warning:
        w.notes.append(f"spread {spread:.1%} above {w.spread_warning:.0%}")
    return {
        "leg": w.leg, "case": w.case, "role": w.role, "backend": w.backend, "np": ranks, "repetitions": len(pairs),
        "step_cost_s": cost, "spread": spread,
        "steps_cost_s": statistics.median((b.steps - a.steps) / dn for a, b in pairs),
        "regrid_cost_s": statistics.median((b.regrid - a.regrid) / dn for a, b in pairs),
        "seam_cost_s": statistics.median((b.seam - a.seam) / dn for a, b in pairs),
        "cells_mean": cells, "ns_per_cell_step": cost * ranks / cells * 1.0e9 if cells > 0 else float("nan"),
        "cells_per_device": cells / ranks,
        "reduced": pairs[0][1].info.get("reduced", "0"), "commit": pairs[0][1].info.get("commit", "unknown"),
        "notes": "; ".join(w.notes),
    }


def ratios(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Return the AMR ratios of every leg, backend and rank count that has its three roles (matrix cases only)."""
    by: dict[tuple[int, str, int], dict[str, dict[str, object]]] = {}
    for row in rows:
        if row["kind"] != "matrix":
            continue
        by.setdefault((row["leg"], row["backend"], row["np"]), {})[str(row["role"])] = row
    out = []
    for (leg, backend, ranks), roles in sorted(by.items()):
        u, r, f = roles.get("uniform"), roles.get("runtime"), roles.get("frozen")
        entry: dict[str, object] = {"leg": leg, "backend": backend, "np": ranks}
        if u and r:
            entry["amr_speedup"] = u["step_cost_s"] / r["step_cost_s"]
            entry["amr_speedup_ideal"] = u["cells_mean"] / r["cells_mean"]
        if u and f:
            entry["grid_efficiency"] = (u["step_cost_s"] * f["cells_mean"] / u["cells_mean"]) / f["step_cost_s"]
        if r and f:
            entry["regrid_cost"] = r["step_cost_s"] / f["step_cost_s"] - 1.0
        if len(entry) > 3:
            out.append(entry)
    return out


def backend_speedups(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Return, per case, T(cpu at its largest rank count) / T(each other backend and rank count)."""
    out = []
    for case in sorted({str(r["case"]) for r in rows}):
        mine = [r for r in rows if r["case"] == case]
        cpu = [r for r in mine if r["backend"] == "cpu"]
        if not cpu:
            continue
        ref = max(cpu, key=lambda r: r["np"])
        for r in mine:
            if r is ref:
                continue
            out.append({"case": case, "reference": f"cpu np {ref['np']}", "backend": r["backend"], "np": r["np"],
                        "speedup": ref["step_cost_s"] / r["step_cost_s"]})
    return out


def scaling(rows: list[dict[str, object]]) -> dict[str, str]:
    """Add the strong and weak efficiencies to the rows; return the label of every series met.

    Every row gains `strong_ref_np`, `strong_efficiency` (per backend and case, against the smallest rank count
    present), `weak_group`, `weak_ref_np` and `weak_efficiency` (per backend and series, the windows grouped by cells
    per device within the series' load tolerance, against the group's smallest device count); None where undefined.

    Parameters
    ----------
    rows : list[dict[str, object]]
        The windows of `summarise`, with their `series`, `series_label` and `load_tolerance`.

    Returns
    -------
    dict[str, str]
        Series name -> its label.
    """
    for row in rows:
        row |= {"strong_ref_np": None, "strong_efficiency": None, "weak_group": None, "weak_ref_np": None,
                "weak_efficiency": None}
    by_case: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        by_case.setdefault((str(row["backend"]), str(row["case"])), []).append(row)
    for members in by_case.values():
        ref = min(members, key=lambda r: r["np"])
        if len(members) > 1:
            for row in members:
                row["strong_ref_np"] = ref["np"]
                row["strong_efficiency"] = ref["ns_per_cell_step"] / row["ns_per_cell_step"]
    labels: dict[str, str] = {}
    by_series: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in rows:
        if row["series"]:
            labels.setdefault(str(row["series"]), str(row["series_label"]))
            by_series.setdefault((str(row["backend"]), str(row["series"])), []).append(row)
    for (_, series), members in by_series.items():
        members.sort(key=lambda r: (-float(r["cells_per_device"]), r["np"]))
        groups: list[list[dict[str, object]]] = []
        for row in members:
            if groups:
                first = float(groups[-1][0]["cells_per_device"])
                if abs(float(row["cells_per_device"]) - first) <= float(row["load_tolerance"]) * first:
                    groups[-1].append(row)
                    continue
            groups.append([row])
        for group in groups:
            load = statistics.median(float(r["cells_per_device"]) for r in group)
            ref = min(group, key=lambda r: r["np"])
            for row in group:
                row["weak_group"] = f"{series} {load / 1.0e6:.3g} M/dev"
                if len(group) > 1:
                    row["weak_ref_np"] = ref["np"]
                    row["weak_efficiency"] = ref["ns_per_cell_step"] / row["ns_per_cell_step"]
    return labels


def fmt(value: object, kind: str) -> str:
    """Format a number for the Markdown tables."""
    if value is None or value == "":
        return "-"
    if kind == "s":
        return f"{value:.4g}"
    if kind == "pct":
        return f"{value:.1%}"
    if kind == "x":
        return f"{value:.3g}"
    if kind == "int":
        return f"{value:,.0f}"
    return str(value)


def write_summary(out: Path, work: Path, rows: list[dict[str, object]], amr: list[dict[str, object]],
                  speed: list[dict[str, object]], labels: dict[str, str], problems: list[str]) -> None:
    """Write summary.csv and summary.md."""
    out.mkdir(parents=True, exist_ok=True)
    columns = ["leg", "case", "role", "backend", "np", "repetitions", "step_cost_s", "spread", "steps_cost_s",
               "regrid_cost_s", "seam_cost_s", "cells_mean", "ns_per_cell_step", "reduced", "commit", "notes",
               "kind", "series", "cells_per_device", "strong_ref_np", "strong_efficiency", "weak_group", "weak_ref_np",
               "weak_efficiency"]
    with (out / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    reduced = any(r["reduced"] == "1" for r in rows)
    lines = [
        "# FLUME speedup matrix: summary", "",
        f"Written by `analyze.py` on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} from `{work}` "
        f"(commits: {', '.join(sorted({str(r['commit'])[:8] for r in rows})) or 'none'}). Protocol and ratios: "
        "`README.md` and [issue #133](https://github.com/szaghi/adam/issues/133).", "",
    ]
    if reduced:
        lines += ["> **Reduced sizes** (`--reduced`): a plumbing test, the timings mean nothing.", ""]
    lines += ["## Step costs", "",
              "| Leg | Case | Role | Backend | Ranks | Step cost (s) | Spread | Steps (s) | Regrid (s/step) | "
              "Seam fills (s/step) | N cells | Cells/device | ns/cell-step |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        flag = " (!)" if "spread" in str(r["notes"]) else ""
        lines.append(f"| {r['leg']} | {r['case']} | {ROLE_TAG[str(r['role'])]} | {r['backend']} | {r['np']} | "
                     f"{fmt(r['step_cost_s'], 's')} | {fmt(r['spread'], 'pct')}{flag} | {fmt(r['steps_cost_s'], 's')} | "
                     f"{fmt(r['regrid_cost_s'], 's')} | {fmt(r['seam_cost_s'], 's')} | {fmt(r['cells_mean'], 'int')} | "
                     f"{fmt(r['cells_per_device'], 'int')} | {fmt(r['ns_per_cell_step'], 'x')} |")
    lines += ["", "## AMR ratios", "",
              "U uniform, R runtime, F frozen. AMR speedup T(U)/T(R) with its ideal N(U)/N(R); grid efficiency "
              "[T(U) N(F)/N(U)]/T(F); regrid cost T(R)/T(F) - 1.", "",
              "| Leg | Backend | Ranks | AMR speedup | ideal | grid efficiency | regrid cost |", "|---|---|---|---|---|---|---|"]
    for a in amr:
        lines.append(f"| {a['leg']} | {a['backend']} | {a['np']} | {fmt(a.get('amr_speedup'), 'x')} | "
                     f"{fmt(a.get('amr_speedup_ideal'), 'x')} | {fmt(a.get('grid_efficiency'), 'x')} | "
                     f"{fmt(a.get('regrid_cost'), 'pct')} |")
    lines += ["", "## Scaling and grind time", "",
              "g = step cost x devices / N, ns per cell-step (devices = ranks: GPUs of an FNL run, cores of a CPU "
              "run). Strong efficiency g(ref)/g per backend and case, ref the smallest rank count present; weak "
              "efficiency g(ref)/g per backend and series at fixed cells per device (the load group), ref the "
              "group's smallest device count. Rows sorted by backend, series and load: read g against cells/device "
              "as the grind curve.", "",
              "| Backend | Series | Case | Devices | N cells | Cells/device | Step cost (s) | g (ns/cell-step) | "
              "Strong eff. (ref np) | Load group | Weak eff. (ref np) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    order = sorted(rows, key=lambda r: (str(r["backend"]), str(r["series"]) or "~",
                                        "" if r["series"] else str(r["case"]), -float(r["cells_per_device"]),
                                        int(r["np"])))
    for r in order:
        strong = f"{fmt(r['strong_efficiency'], 'x')} ({r['strong_ref_np']})" if r["strong_ref_np"] else "-"
        weak = f"{fmt(r['weak_efficiency'], 'x')} ({r['weak_ref_np']})" if r["weak_ref_np"] else "-"
        lines.append(f"| {r['backend']} | {r['series'] or '-'} | {r['case']} | {r['np']} | "
                     f"{fmt(r['cells_mean'], 'int')} | {fmt(r['cells_per_device'], 'int')} | "
                     f"{fmt(r['step_cost_s'], 's')} | {fmt(r['ns_per_cell_step'], 'x')} | {strong} | "
                     f"{r['weak_group'] or '-'} | {weak} |")
    if labels:
        lines += ["", "Series:", "", *[f"- `{name}`: {label}" for name, label in sorted(labels.items())]]
    lines += ["", "## Backend speedups", "", "| Case | Reference | Backend | Ranks | Speedup |", "|---|---|---|---|---|"]
    for s in speed:
        lines.append(f"| {s['case']} | {s['reference']} | {s['backend']} | {s['np']} | {fmt(s['speedup'], 'x')} |")
    notes = [f"- {r['backend']} np {r['np']} {r['case']}: {r['notes']}" for r in rows if r["notes"]]
    notes += [f"- {p}" for p in problems]
    if notes:
        lines += ["", "## Warnings", "", *notes]
    (out / "summary.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    """Parse the command line, read the windows, write the summary."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work", type=Path, default=HERE / "runs" / "full", help="run.sh's run tree (default runs/full)")
    parser.add_argument("--out", type=Path, default=HERE / "results", help="the summary directory (default results/)")
    args = parser.parse_args()
    windows, problems = collect(args.work)
    if not windows:
        print(f"analyze: no window run under {args.work}/windows ({len(problems)} unreadable)", file=sys.stderr)
        for p in problems:
            print(f"   {p}", file=sys.stderr)
        return 1
    rows = []
    for key in sorted(windows, key=lambda k: (windows[k].leg, windows[k].case, k[0], k[1])):
        row = summarise(windows[key])
        if row is None:
            problems.append(f"{key[0]} np {key[1]} {key[2]}: no repetition with both runs")
            continue
        conf = read_case_conf(args.work, key[2])
        row |= {"kind": conf.get("KIND", "matrix"), "series": conf.get("SERIES", ""),
                "series_label": conf.get("SERIES_LABEL", ""),
                "load_tolerance": float(conf.get("LOAD_TOLERANCE", "0.10"))}
        rows.append(row)
    amr, speed = ratios(rows), backend_speedups(rows)
    labels = scaling(rows)
    write_summary(args.out, args.work, rows, amr, speed, labels, problems)
    print(f"analyze: {len(rows)} windows, {len(amr)} AMR ratio rows, {len(speed)} backend speedups -> "
          f"{args.out / 'summary.md'}, {args.out / 'summary.csv'}")
    for r in rows:
        if r["notes"]:
            print(f"   warning: {r['backend']} np {r['np']} {r['case']}: {r['notes']}")
    for p in problems:
        print(f"   warning: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
