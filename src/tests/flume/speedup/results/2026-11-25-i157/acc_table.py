#!/usr/bin/env python3
"""Device time per kernel from NVCOMPILER_ACC_TIME logs: the #133 kernel table, before and after #157.

Every rank prints its own report at exit, one block per (function, device); the blocks of ranks that exit together
interleave, which attributes lines to the wrong function on one or two devices of a run (more under mpirun than under
srun). The table therefore takes, for every row, the value most devices agree on within 5 % (the median of that
group), in s of device time per GPU for the 10-step run.

Usage:
    acc_table.py <label>=<log.txt> [<label>=<log.txt> ...]
"""

from __future__ import annotations

import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

HEADER = re.compile(r"^\s{2}(\S+)\s+NVIDIA\s+devicenum=(\d+)")
ELAPSED = re.compile(r"elapsed time\(us\): total=([\d,]+)")
LAUNCHED = re.compile(r"kernel launched (\d+) times")
ROWS = (
    ("all kernels", None),
    ("face fluxes", "compute_face_fluxes_s3_dev"),
    ("intra-rank ghost fill", "update_ghost_local_gpu_dev"),
    ("ghost send buffers", "populate_send_buffer_ghost_gpu_dev"),
    ("boundary conditions", "set_boundary_conditions_dev"),
)


def parse(log: Path) -> tuple[dict[str, dict[int, float]], dict[str, dict[int, int]]]:
    """Return the device time (s) and the launches of every (function, device) of one log."""
    time: dict[str, dict[int, float]] = defaultdict(lambda: defaultdict(float))
    launches: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    fn, dev = None, -1
    for line in log.read_text(errors="replace").splitlines():
        if m := HEADER.match(line):
            fn, dev = m.group(1), int(m.group(2))
        elif fn and (m := ELAPSED.search(line)):
            time[fn][dev] += int(m.group(1).replace(",", "")) * 1e-6
        elif fn and (m := LAUNCHED.search(line)):
            launches[fn][dev] += int(m.group(1))
    return time, launches


def consensus(values: list[float], tol: float = 0.05) -> float:
    """The median of the largest group of values that agree within `tol` (relative) with one of them.

    A zero is a device whose report lost the kernel to the interleaving (every rank runs every kernel): it does not
    vote unless every device shows zero.
    """
    values = [v for v in values if v > 0.0] or values
    best: list[float] = []
    for v in values:
        group = [w for w in values if abs(w - v) <= tol * max(abs(v), 1e-12)]
        if len(group) > len(best):
            best = group
    return statistics.median(best) if best else 0.0


def row_values(time: dict[str, dict[int, float]]) -> dict[str, float]:
    """The consensus over the devices of each table row."""
    devices = sorted({d for per_dev in time.values() for d in per_dev})
    out = {}
    for name, fn in ROWS:
        if fn is None:
            vals = [
                sum(per_dev.get(d, 0.0) for per_dev in time.values()) for d in devices
            ]
        else:
            vals = [time[fn].get(d, 0.0) for d in devices]
        out[name] = consensus(vals)
    return out


def main() -> int:
    """Print the markdown table of the given logs."""
    cols = [arg.split("=", 1) for arg in sys.argv[1:]]
    if not cols:
        print(__doc__)
        return 2
    values = [row_values(parse(Path(path))[0]) for _, path in cols]
    print(
        "| Kernel (per GPU, 10-step run) | "
        + " | ".join(label for label, _ in cols)
        + " |"
    )
    print("|---|" + "---|" * len(cols))
    for name, _ in ROWS:
        print(f"| {name} | " + " | ".join(f"{v[name]:.2f}" for v in values) + " |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
