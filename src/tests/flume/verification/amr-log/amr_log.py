#!/usr/bin/env python3
"""The A4 logs of #68: integral drift across static 2:1 seams, and the WENO sensor on the faces near them.

Why: the reflux policy of #68 E (`[numerics] reflux_policy = always|troubled`) turns the Berger-Colella reflux off at the seam
faces a WENO sensor judges smooth. Two numbers decide whether that is safe and how to tune it, and neither was
logged: how the conserved integrals drift when the reflux is off (bounded, or accumulating?), and what the sensor reads
on the seam faces of a smooth flow (false positives) and of a shock crossing (detections).

* drift: the relative drift (I(t) - I(0)) / |I(0)| of the domain integrals of the conservation history
  ([diagnostics] conservation_history_save), its maximum, its final value, and its growth (the mean |drift| over the
  last quarter of the run over that of the second quarter: about 1 for a bounded drift, 2 or more when it
  accumulates, as a linear drift gives 7/3; "-" when the second quarter holds less than 1 % of the maximum drift,
  as before a shock first crosses a face: the drift is then a jump, read from the maximum);
* sensor: on every checkpoint, the WENO-5 sensor r = min_k w_k / d_k of every face along every non-null axis, the
  minimum over the primitive density, normal velocity and pressure and over the reconstructions of the two cells
  beside the face (the weno-riemann sensor of #47, FLUX_CORRECTION_SENSOR_TAU = 0.2, read on primitive variables).
  The weights are those of adam_weno_object (Jiang-Shu, zeps = 1e-6, exponent S = 3, its smoothness-indicator
  table), recomputed from the checkpoint with its ghosts, which are filled right before the save as the residual
  stencils read them (#65). The faces split into the seam band (within ngc cells of their own spacing from the 2:1
  faces) and the interior; per checkpoint (at most --rows of them, evenly spaced): the minimum r and the fraction of
  faces below tau, per region; then, over all checkpoints, the minimum r and the largest and mean fractions.

The default scheme (`weno`) weights its characteristic split fluxes, not these variables: the map is the sensor a
seam policy would read, not a dump of the scheme's own weights.

Usage:
    amr_log.py drift <work> [<work> ...] [--vars r rE] [--max D]
    amr_log.py sensor <work> [--seam-x X | --box xmin ymin xmax ymax] [--tau 0.2] [--rows 25]
    --max D fails unless the maximum |drift| of every listed variable of every run is at most D
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

VERIF = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(VERIF / "amr-order"))
from amr_order_oracle import box_distance  # noqa: E402

ZEPS = 1.0e-6
WEXP = 3
# optimal weights of the stencils 0, 1, 2 at the left (i - 1/2) and right (i + 1/2) interface (adam_weno_object)
D_LEFT = np.array([0.3, 0.6, 0.1])
D_RIGHT = np.array([0.1, 0.6, 0.3])
# smoothness indicators: IS_s = sum_{a, b} C[s][a][b] v(s - a) v(s - b), the table of initialize_upwind_S3
C = np.zeros((3, 3, 3))
C[0, 0, 0], C[0, 1, 0], C[0, 2, 0], C[0, 1, 1], C[0, 2, 1], C[0, 2, 2] = 10 / 3, -31 / 3, 11 / 3, 25 / 3, -19 / 3, 4 / 3
C[1, 0, 0], C[1, 1, 0], C[1, 2, 0], C[1, 1, 1], C[1, 2, 1], C[1, 2, 2] = 4 / 3, -13 / 3, 5 / 3, 13 / 3, -13 / 3, 4 / 3
C[2, 0, 0], C[2, 1, 0], C[2, 2, 0], C[2, 1, 1], C[2, 2, 1], C[2, 2, 2] = 4 / 3, -19 / 3, 11 / 3, 25 / 3, -31 / 3, 10 / 3


def read_ini(work: Path) -> configparser.ConfigParser:
    """Return the input of a run (its only .ini file)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    return ini


def drift(works: list[Path], names: list[str], bound: float | None) -> int:
    """Print the drift of the listed integrals of every run; return the exit status."""
    status = 0
    print(f"{'run':52s} {'var':>6s} {'max |drift|':>12s} {'final':>12s} {'growth':>8s}")
    for work in works:
        path = next(work.glob("*-conservation_history.dat"))
        lines = path.read_text().splitlines()
        header = lines[0].split("=", 1)[1].replace('"', "").split()
        data = np.array([[float(v) for v in ln.split()] for ln in lines[1:] if ln.strip()])
        n = len(data)
        for name in names:
            col = data[:, header.index(f"int_{name}")]
            rel = (col - col[0]) / abs(col[0])
            q2, q4 = np.abs(rel[n // 4:n // 2]), np.abs(rel[3 * n // 4:])
            peak = float(np.max(np.abs(rel)))
            growth = f"{np.mean(q4) / np.mean(q2):8.2f}" if np.mean(q2) > 0.01 * peak > 0.0 else f"{'-':>8s}"
            verdict = ""
            if bound is not None:
                ok = peak <= bound
                status |= 0 if ok else 1
                verdict = "  PASS" if ok else f"  FAIL (bound {bound:.1e})"
            print(f"{work.name:52s} {name:>6s} {peak:12.3e} {rel[-1]:12.3e} {growth}{verdict}")
    return status


def interface_ratios(v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return min_k w_k / d_k of the left and right reconstructions of every cell of v (axis 0), cells 2 .. n-3."""
    n = v.shape[0]
    cells = np.arange(2, n - 2)
    IS = np.zeros((3, len(cells)) + v.shape[1:])
    for s in range(3):
        for a in range(3):
            for b in range(3):
                if C[s, a, b] != 0.0:
                    IS[s] += C[s, a, b] * v[cells + s - a] * v[cells + s - b]
    inv = 1.0 / (ZEPS + IS) ** WEXP
    ratios = []
    for d in (D_LEFT, D_RIGHT):
        alpha = d.reshape((3,) + (1,) * (IS.ndim - 1)) * inv
        w = alpha / np.sum(alpha, axis=0)
        ratios.append(np.min(w / d.reshape((3,) + (1,) * (IS.ndim - 1)), axis=0))
    return ratios[0], ratios[1]


def face_sensor(prims: list[np.ndarray], axis: int) -> np.ndarray:
    """Return the sensor of the faces between consecutive cells along axis (cells 2 .. n-3 reconstructed): face j sits
    between cells j + 2 and j + 3, the minimum over the variables of the right reconstruction of the cell before and
    the left reconstruction of the cell after."""
    r = None
    for v in prims:
        left, right = interface_ratios(np.moveaxis(v, axis, 0))
        face = np.minimum(right[:-1], left[1:])
        r = face if r is None else np.minimum(r, face)
    return r


def sensor(work: Path, seam_x: float | None, box: list[float] | None, tau: float, rows: int) -> int:
    """Print the sensor per checkpoint and region; return the exit status."""
    ini = read_ini(work)
    ngc = int(ini["grid"]["ngc"])
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    axes = [a for a, x in enumerate("xyz") if ini["grid"][f"null_{x}"].strip().lower() != ".true."]
    lo = np.array([float(ini["grid"][f"emin_{x}"]) for x in "xy"])
    size = np.array([float(ini["grid"][f"emax_{x}"]) for x in "xy"]) - lo
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    steps = sorted({int(p.name.split("-")[-2]) for p in files})
    times = {}
    for line in next(work.glob("*-conservation_history.dat")).read_text().splitlines()[1:]:
        if line.strip():
            it, t = line.split()[:2]
            times[int(it)] = float(t)
    print(f"{work.name}: sensor r = min_k w_k/d_k on (rho, u_n, p), tau = {tau}")
    print(f"{'it':>7s} {'time':>9s} | {'band faces':>10s} {'min r':>9s} {'r < tau':>8s} | {'interior':>9s} "
          f"{'min r':>9s} {'r < tau':>8s}")
    shown = set(steps[:: max(1, -(-len(steps) // rows))]) | {steps[-1]}
    summary: dict[str, list[tuple[float, float]]] = {"band": [], "interior": []}
    for step in steps:
        rs, bands = [], []
        for path in (p for p in files if int(p.name.split("-")[-2]) == step):
            with h5py.File(path, "r") as h5:
                for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                    origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                    q = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in ("r", "ru", "rv", "rw", "rE")}
                    rho = q["r"]
                    vel = [q["ru"] / rho, q["rv"] / rho, q["rw"] / rho]
                    p = (gamma - 1.0) * (q["rE"] - 0.5 * rho * (vel[0] ** 2 + vel[1] ** 2 + vel[2] ** 2))
                    # the interior transverse cells only (a null axis keeps its first interior cell)
                    sl = [slice(ngc, -ngc) if a in axes else slice(ngc, ngc + 1) for a in range(3)]
                    for axis in axes:
                        cut = list(sl)
                        cut[axis] = slice(None)
                        r = face_sensor([rho[tuple(cut)], vel[axis][tuple(cut)], p[tuple(cut)]], axis)
                        # faces from the one between ghost ngc-1 and cell ngc onwards: drop the outer ones
                        n_int = rho.shape[axis] - 2 * ngc
                        r = np.take(r, np.arange(ngc - 3, ngc - 3 + n_int + 1), axis=0)
                        centres = [origin[a] + (np.arange(rho.shape[a]) + 0.5) * dxyz[a] for a in range(3)]
                        fpos = origin[axis] + (np.arange(ngc, ngc + n_int + 1)) * dxyz[axis]
                        grids = [fpos if a == axis else centres[a][sl[a]] for a in range(3)]
                        mesh = np.meshgrid(*grids, indexing="ij")
                        mesh = [np.moveaxis(m, axis, 0) for m in mesh]
                        h = dxyz[axis]
                        if seam_x is not None:
                            band = np.abs(mesh[0] - seam_x) <= ngc * h
                        elif box is not None:
                            xy = np.stack([mesh[0].ravel(), mesh[1].ravel()], axis=1)
                            band = (box_distance(xy, box, lo, size) <= ngc * h).reshape(mesh[0].shape)
                        else:
                            band = np.zeros(r.shape, dtype=bool)
                        rs.append(r.ravel())
                        bands.append(band.ravel())
        r, band = np.concatenate(rs), np.concatenate(bands)
        row = f"{step:7d} {times.get(step, float('nan')):9.4f} |"
        for region, mask in (("band", band), ("interior", ~band)):
            if mask.any():
                low, frac = float(r[mask].min()), float(np.mean(r[mask] < tau))
                summary[region].append((low, frac))
                row += f" {int(mask.sum()):10d} {low:9.3e} {frac:8.3f} |"
            else:
                row += f" {0:10d} {'-':>9s} {'-':>8s} |"
        if step in shown:
            print(row.rstrip(" |"))
    for region, values in summary.items():
        if values:
            lows, fracs = zip(*values, strict=True)
            print(f"summary {region}: min r {min(lows):.3e}, r < tau on at most {max(fracs):.3f} of the faces "
                  f"(mean {sum(fracs) / len(fracs):.3f} over {len(values)} checkpoints)")
    return 0


def main() -> int:
    """Dispatch the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("drift")
    d.add_argument("work", type=Path, nargs="+")
    d.add_argument("--vars", nargs="+", default=["r", "rE"])
    d.add_argument("--max", type=float, default=None)
    s = sub.add_parser("sensor")
    s.add_argument("work", type=Path)
    s.add_argument("--seam-x", type=float, default=None)
    s.add_argument("--box", type=float, nargs=4, default=None)
    s.add_argument("--tau", type=float, default=0.2)
    s.add_argument("--rows", type=int, default=25, help="checkpoints printed (the summary uses all)")
    args = parser.parse_args()
    if args.cmd == "drift":
        return drift(args.work, args.vars, args.max)
    return sensor(args.work, args.seam_x, args.box, args.tau, args.rows)


if __name__ == "__main__":
    sys.exit(main())
