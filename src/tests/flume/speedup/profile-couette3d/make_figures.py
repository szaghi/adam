#!/usr/bin/env python3
"""Collect the measurements of the FLUME profiling case and draw the figures of its README.

Why this script exists: the README of the profiling case reports measurements, and a measurement written into a
README by hand drifts from the runs it came from. `collect` reads every run directory that `profile.sh` wrote
(`run.info`, `wall.txt`, the forest timing line of `log.txt`, the Nsight Systems kernel summary) into one CSV row per
run; `figures` draws the README figures from that CSV and from the checkpoints of the correctness runs, so every
number of a figure can be traced to a run directory and redrawn after new runs.

- `results/runs.csv`: one row per run (variant, size, Slurm job, date, commit, executable, wall time, forest loop and
  steps, the mean time per launch of the two flux kernels, the total kernel time).
- `figures/a30-variants.png`: the A30 loop of each build variant split into the face-flux kernel, the dissipative
  kernel, the other kernels and the time with no kernel running (uniform and refined walls).
- `figures/cpu-weno.png`: the CPU runs before and after the WENO hot-path fix (commit f46b56fb).
- `figures/fields.png`: u and T of the correctness runs against the exact Couette solution, their error, and the
  largest difference of every variant from the leg.

Usage:
    make_figures.py collect <runs-dir> [--out results/runs.csv]
    make_figures.py figures <runs-dir> [--results results/runs.csv] [--out figures]

`<runs-dir>` is `runs/` of this directory, or a copy of the one on the machine that ran the profiles. `collect` needs
the kernel summary `kern_cuda_gpu_kern_sum.csv` of an nsys run, and writes it (`nsys stats`) when nsys is on PATH and
the file is missing. `figures` needs numpy, matplotlib and, for the fields, h5py.
"""

from __future__ import annotations

import argparse
import configparser
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KERN = "kern_cuda_gpu_kern_sum.csv"
FACE, DISSIPATIVE = "compute_face_fluxes_dev", "add_dissipative_fluxes_o4_dev"
FIELDS = (
    "variant",
    "size",
    "tool",
    "np",
    "steps",
    "slurm_job",
    "date",
    "host",
    "commit",
    "exe",
    "wall_s",
    "loop_s",
    "steps_s",
    "face_ms",
    "face_s",
    "dissipative_ms",
    "dissipative_s",
    "kernels_s",
)
# the variants of the A30 experiment, in the order of the figures; "leg" is the build of the fnl-nvf leg
VARIANTS = ("leg", "r128", "lto", "lto-r128", "B0", "B0-lto")
# the CPU A/B of the WENO fix (Slurm 1357589): the "before" executable was built in this worktree, the "after" one
# from the uncommitted fix, committed afterwards as f46b56fb (run.info records the HEAD of the time, ff80d8da)
CPU_AB_JOB, CPU_BEFORE = "1357589", "weno-before"
# the runs the variants figure averages: one rank (the kernel summary is rank 0's), the 20 steps of profile.sh
FIG_NP, FIG_STEPS = "1", "20"
# the Couette constants of make_profile.py (and of regression/couette, verification/viscous/profiles.py)
GAMMA, MU, PRANDTL, U, TW, H = 1.4, 0.01, 0.75, 1.0, 1.0, 1.0


def read_info(run: Path) -> dict[str, str]:
    """Return the key=value pairs of a run.info."""
    info = {}
    for line in (run / "run.info").read_text().splitlines():
        key, _, value = line.partition("=")
        info[key] = value
    return info


def variant_of(info: dict[str, str]) -> str:
    """Return the build variant of a run: its tag, or the directory of its executable (runs before --tag existed)."""
    if info.get("tag"):
        return info["tag"]
    exe = info.get("exe") or info["command"].split()[-2]
    parent = Path(exe).parent.name
    if parent.startswith(info["leg"] + "-"):
        return parent[len(info["leg"]) + 1 :]
    return "before" if CPU_BEFORE in exe else "leg"


def kernel_summary(run: Path) -> list[dict[str, str]]:
    """Return the rows of the nsys kernel summary of a run, written first with `nsys stats` if it is missing."""
    path = run / KERN
    rep = run / "nsys-r0.nsys-rep"
    if not path.is_file() and rep.is_file() and shutil.which("nsys"):
        subprocess.run(
            [
                "nsys",
                "stats",
                "--report",
                "cuda_gpu_kern_sum",
                "--format",
                "csv",
                "-o",
                str(run / "kern"),
                str(rep),
            ],
            check=False,
            capture_output=True,
        )
    if not path.is_file():
        return []
    lines = [
        line
        for line in path.read_text().splitlines()
        if line.startswith(("Time", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    ]
    return list(csv.DictReader(lines))


def collect(runs: Path, out: Path) -> None:
    """Write one CSV row per run directory of `runs` (directories without run.info are skipped)."""
    rows = []
    for run in sorted(p for p in runs.iterdir() if (p / "run.info").is_file() and (p / "wall.txt").is_file()):
        info = read_info(run)
        timing = re.search(
            r"forest timing.*?loop\s+([\d.]+), steps\s+([\d.]+)",
            (run / "log.txt").read_text(),
        )
        row = {
            "variant": variant_of(info),
            "size": info["size"],
            "tool": info["tool"],
            "np": info["np"],
            "steps": info.get("steps", ""),
            "slurm_job": info.get("slurm_job", ""),
            "date": info["date"],
            "host": info["host"],
            "commit": info["commit"][:8],
            "exe": info.get("exe") or info["command"].split()[-2],
            "wall_s": (run / "wall.txt").read_text().strip(),
            "loop_s": timing.group(1) if timing else "",
            "steps_s": timing.group(2) if timing else "",
        }
        kernels = kernel_summary(run)
        row["kernels_s"] = f"{sum(float(k['Total Time (ns)']) for k in kernels) / 1e9:.4f}" if kernels else ""
        for key, name in (("face", FACE), ("dissipative", DISSIPATIVE)):
            # every kernel of that name (one per model or per stencil width): total time, and mean per launch
            hit = [k for k in kernels if name in k["Name"]]
            total = sum(float(k["Total Time (ns)"]) for k in hit)
            launches = sum(int(k["Instances"]) for k in hit)
            row[f"{key}_ms"] = f"{total / launches / 1e6:.4f}" if hit else ""
            row[f"{key}_s"] = f"{total / 1e9:.4f}" if hit else ""
        rows.append(row)
    rows.sort(key=lambda r: r["date"])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"collect: {len(rows)} runs -> {out}")


def read_rows(path: Path) -> list[dict[str, str]]:
    """Return the rows of a results CSV."""
    with path.open() as f:
        return list(csv.DictReader(f))


def mean(values: list[float]) -> float:
    """Return the mean of a non-empty list."""
    return sum(values) / len(values)


def fig_variants(rows: list[dict[str, str]], out: Path) -> None:
    """Draw the A30 loop of each variant, split by kernel, for the uniform and the refined-walls runs."""
    import matplotlib.pyplot as plt  # noqa: PLC0415

    fig, axs = plt.subplots(1, 2, figsize=(13.5, 4.6), constrained_layout=True)
    parts = (
        ("face-flux kernel", "C0"),
        ("dissipative kernel", "C1"),
        ("other kernels", "C2"),
        ("no kernel running", "0.75"),
    )
    for ax, size, title in (
        (axs[0], "medium-uniform", "uniform, 16.8 M cells"),
        (axs[1], "medium", "refined walls, 5.8 M cells, a regrid every 10 steps"),
    ):
        names, stacks, loops = [], [], []
        for v in VARIANTS:
            sel = [
                r
                for r in rows
                if r["variant"] == v
                and r["size"] == size
                and r["tool"] == "nsys"
                and r["np"] == FIG_NP
                and r["steps"] == FIG_STEPS
                and r["face_s"]
            ]
            if not sel:
                continue
            if len({r["commit"] for r in sel}) > 1:
                print(f"figures: warning, {v} {size} averages runs of several commits")
            loop = mean([float(r["loop_s"]) for r in sel])
            kern = mean([float(r["kernels_s"]) for r in sel])
            face = mean([float(r["face_s"]) for r in sel])
            diss = mean([float(r["dissipative_s"]) for r in sel])
            names.append(f"{v} ({len(sel)} run{'s' if len(sel) > 1 else ''})")
            stacks.append((face, diss, max(kern - face - diss, 0.0), max(loop - kern, 0.0)))
            loops.append(loop)
        y = range(len(names))
        left = [0.0] * len(names)
        for k, (label, color) in enumerate(parts):
            widths = [s[k] for s in stacks]
            ax.barh(
                list(y),
                widths,
                left=left,
                color=color,
                label=label,
                edgecolor="k",
                lw=0.4,
            )
            left = [a + b for a, b in zip(left, widths, strict=True)]
        for i, loop in enumerate(loops):
            ax.text(
                max(left[i], loop) + 0.6,
                i,
                f"{loop:.1f} s ({100 * (loop / loops[0] - 1):+.0f}%)" if i else f"{loop:.1f} s",
                va="center",
                fontsize=8,
            )
        ax.set_yticks(list(y), names)
        ax.invert_yaxis()
        ax.set_xlim(0, max(loops) * 1.28)
        ax.set_xlabel("forest loop, 20 steps (s)")
        ax.set_title(f"one A30, medium {title}", fontsize=9)
        ax.grid(axis="x", alpha=0.3)
    axs[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(
        "FLUME FNL on one A30: where the loop goes, per build variant (Nsight Systems; mean over the runs)",
        fontsize=10,
    )
    fig.savefig(out / "a30-variants.png", dpi=110)
    plt.close(fig)


def fig_cpu(rows: list[dict[str, str]], out: Path) -> None:
    """Draw the alternated CPU runs before and after the WENO hot-path fix."""
    import matplotlib.pyplot as plt  # noqa: PLC0415

    sel = [r for r in rows if r["slurm_job"] == CPU_AB_JOB and r["size"] == "medium-uniform"]
    if not sel:
        print("figures: no CPU A/B run, cpu-weno.png skipped")
        return
    fig, ax = plt.subplots(figsize=(7.0, 3.6), constrained_layout=True)
    for v, color, label in (
        ("before", "C3", "before (ff80d8da)"),
        ("leg", "C0", "after (the fix, committed as f46b56fb)"),
    ):
        pts = [(i, float(r["wall_s"])) for i, r in enumerate(sel, start=1) if r["variant"] == v]
        ax.plot(
            [p[0] for p in pts],
            [p[1] for p in pts],
            "o",
            color=color,
            ms=8,
            label=label,
        )
        m = mean([p[1] for p in pts])
        ax.axhline(m, color=color, lw=1, ls="--")
        below = v == "before"  # the before line runs along the top of the axes: label it underneath
        ax.text(
            0.6,
            m - 2.0 if below else m + 2.0,
            f"mean {m:.1f} s",
            color=color,
            va="top" if below else "bottom",
            fontsize=8,
        )
    ax.set_xticks(range(1, len(sel) + 1))
    ax.set_xlabel("run, in the order of the job (before and after alternated)")
    ax.set_ylabel("wall time, 20 steps (s)")
    ax.set_xlim(0.5, len(sel) + 0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="center right")
    ax.set_title(
        f"CPU, medium uniform, 16 ranks of a medooza node (Slurm {CPU_AB_JOB}): the WENO hot-path fix",
        fontsize=9,
    )
    fig.savefig(out / "cpu-weno.png", dpi=110)
    plt.close(fig)


def saved_steps(run: Path) -> list[int]:
    """Return the iterations of the checkpoints of a run, in order."""
    return sorted({int(p.name.split("-")[1]) for p in run.glob("couette3d-*-proc*.h5")})


def checkpoint_cells(run: Path, it: int):  # noqa: ANN201
    """Return y and the fields u/U and T/Tw of every interior cell of the checkpoint `it` of a run (np arrays)."""
    import h5py  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(run / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    ys, us, ts = [], [], []
    for path in sorted(run.glob(f"couette3d-{it:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                # stored (z, y, x) with the ghosts: reverse to (x, y, z) and strip them
                origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                f = {
                    v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc]
                    for v in ("r", "ru", "rv", "rw", "rE")
                }
                u = f["ru"] / f["r"]
                kin = 0.5 * (f["ru"] ** 2 + f["rv"] ** 2 + f["rw"] ** 2) / f["r"]
                t = (GAMMA - 1.0) * (f["rE"] - kin) / f["r"]
                y = origin[1] + ngc * d[1] + (np.arange(u.shape[1]) + 0.5) * d[1]
                ys.append(np.broadcast_to(y[None, :, None], u.shape).ravel())
                us.append(u.ravel() / U)
                ts.append(t.ravel() / TW)
    if not ys:
        raise FileNotFoundError(f"no checkpoint {it} in {run}")
    return np.concatenate(ys), np.concatenate(us), np.concatenate(ts)


def exact(y):  # noqa: ANN001, ANN201
    """Return the exact u/U and T/Tw of the compressible Couette flow at heights y."""
    k = MU * GAMMA / (GAMMA - 1.0) / PRANDTL
    return U * y / H / U, (TW + MU * U**2 / (k * H) * (y - y**2 / (2.0 * H))) / TW


def fig_fields(runs: Path, out: Path) -> None:
    """Draw the correctness runs (small, refined walls, `--steps 10 --io`) of every variant against the exact solution."""
    import matplotlib.pyplot as plt  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    base = runs / "fnl-nvf-small-np1-none"
    if not base.is_dir():
        print(f"figures: {base} missing, fields.png skipped")
        return
    first, last = saved_steps(base)[0], saved_steps(base)[-1]
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    y, u, t = checkpoint_cells(base, last)
    yy = np.linspace(0.0, H, 201)
    ue, te = exact(yy)
    levels = np.unique(np.round(y, 12))
    ax = axs[0]
    ax.plot(ue, yy, "C0-", lw=1, label="$u/U$ exact")
    ax.plot(
        [u[np.isclose(y, v)].mean() for v in levels],
        levels,
        "C0o",
        ms=3,
        mfc="none",
        label="$u/U$, leg",
    )
    ax2 = ax.twiny()
    ax2.plot(te, yy, "C3-", lw=1, label="$T/T_w$ exact")
    ax2.plot(
        [t[np.isclose(y, v)].mean() for v in levels],
        levels,
        "C3s",
        ms=3,
        mfc="none",
        label="$T/T_w$, leg",
    )
    ax.set_xlabel("$u/U$", color="C0")
    ax2.set_xlabel("$T/T_w$", color="C3")
    ax.set_ylabel("y")
    ax.legend(loc="upper left", fontsize=7)
    ax2.legend(loc="lower right", fontsize=7)
    ax.set_title(
        f"step {last}, cell values averaged over x and z (the wall layers are one level finer)",
        fontsize=8,
    )
    ax = axs[1]
    for it, fmt in ((first, "C7."), (last, "C0o")):
        yi, ui, ti = checkpoint_cells(base, it)
        eu, et = exact(yi)
        lv = np.unique(np.round(yi, 12))
        ax.semilogy(
            lv,
            [np.abs(ui - eu)[np.isclose(yi, v)].max() for v in lv],
            fmt,
            ms=3,
            mfc="none",
            label=f"$|u - u_e|/U$, step {it}",
        )
        ax.semilogy(
            lv,
            [np.abs(ti - et)[np.isclose(yi, v)].max() for v in lv],
            fmt.replace("o", "s").replace("C0", "C3"),
            ms=3,
            mfc="none",
            label=f"$|T - T_e|/T_w$, step {it}",
        )
    ax.set_xlabel("y")
    ax.set_ylabel("largest error at that height")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    ax.set_title(
        f"the leg against the exact solution (step {first} = the initial condition)",
        fontsize=8,
    )
    ax = axs[2]
    names, diffs = [], []
    for v in VARIANTS[1:]:
        run = runs / f"fnl-nvf-small-np1-none-{v}"
        if not run.is_dir():
            continue
        _, uv, tv = checkpoint_cells(run, last)
        names.append(v)
        diffs.append(max(float(np.abs(uv - u).max()), float(np.abs(tv - t).max())))
    ax.bar(
        names,
        [1.0] * len(names),
        color=["C2" if dv == 0.0 else "C3" for dv in diffs],
        edgecolor="k",
        lw=0.4,
    )
    for i, dv in enumerate(diffs):
        ax.text(
            i,
            0.5,
            "bitwise\nidentical" if dv == 0.0 else f"{dv:.1e}",
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.set_yticks([])
    ax.set_title(f"largest |difference| from the leg, u and T, step {last}, every cell", fontsize=8)
    fig.suptitle(
        f"Correctness of the variants: small refined walls (1.2 M cells), one A30, {last} steps (the regrid at step 10 "
        "leaves the grid unchanged)",
        fontsize=10,
    )
    fig.savefig(out / "fields.png", dpi=110)
    plt.close(fig)


def main() -> int:
    """Run the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("collect")
    p.add_argument("runs", type=Path)
    p.add_argument("--out", type=Path, default=HERE / "results" / "runs.csv")
    p = sub.add_parser("figures")
    p.add_argument("runs", type=Path)
    p.add_argument("--results", type=Path, default=HERE / "results" / "runs.csv")
    p.add_argument("--out", type=Path, default=HERE / "figures")
    args = parser.parse_args()
    if args.cmd == "collect":
        collect(args.runs, args.out)
        return 0
    import matplotlib as mpl  # noqa: PLC0415

    mpl.use("Agg")
    args.out.mkdir(parents=True, exist_ok=True)
    rows = read_rows(args.results)
    fig_variants(rows, args.out)
    fig_cpu(rows, args.out)
    fig_fields(args.runs, args.out)
    print(f"figures: written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
