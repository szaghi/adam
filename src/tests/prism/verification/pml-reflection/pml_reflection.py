#!/usr/bin/env python3
"""Measure the amplitude of the pulse reflected by an absorbing boundary (issue #145, item 3).

Why: the reflection coefficient of a PEC-backed PML at normal incidence is, in the continuum,
R = exp(-2 int_0^L gamma dx / c) = exp(-2 gamma_max L / ((n + 1) c)) for the graded profile gamma_max (d/L)^n, the
same for every frequency; the discrete layer adds the reflection of its own grading. The run of
make_pml_reflection.py starts a Gaussian pulse (B_z, centre 0) toward +x; whatever comes back into the window
[-0.04, 0.04] m after the incident pulse has left it is reflection (or the scheme's dispersive trail, which the `far`
input measures). The reflected amplitude is

    A_r / A_0 = max over the checkpoints after t_in of max over the window of |B_z|, over max |B_z| at t = 0,

and, since a reflection that is not a scaled copy of the pulse (the dispersive discrete reflection of a steep grading,
a first-order boundary) has no single amplitude, the energy measure

    (E_r / E_0)^(1/2) = (max over the checkpoints after t_in of the electromagnetic energy in the window, over the
                         energy at t = 0)^(1/2),

which equals A_r / A_0 for a scaled copy. t_in = (0.04 + 5 sigma + 0.01) / c: the incident pulse (cut at 5 sigma) and
the first centimetre of its dispersive trail have left the window, which brings the floor of the measurement (the
`far` input) from 1.2e-6 to 3.2e-7; no reflection reaches the window before (0.08 + 0.04 - 5 sigma) / c, the leading
edge of a reflection at the inner face of a layer (or at a wall without a layer). Every checkpoint is read (the period of make_pml_reflection.py, 4 steps = 1 mm of travel, against sigma = 4 mm: the spatial maximum
over the window catches the peak while the pulse is inside it).

Usage:
    pml_reflection.py <work> [--json out.json] [--profiles out.npz]
    pml_reflection.py --from-profiles <profiles.npz> [--json out.json]   (measure again from a kept profiles.npz)

Prints A_0, A_r, A_r/A_0, the time and position of the reflected maximum; --json writes them, --profiles keeps the
B_z(x) of every checkpoint (for the figure).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

C0 = 299792458.0
T_IN = (0.04 + 5.0 * 0.004 + 0.01) / C0
MU0 = 1.25663706212e-6
EPS0 = 1.0 / (MU0 * C0**2)
WINDOW = 0.04
SIGMA = 0.004
NGC = 3


def profiles(work: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the steps, the cell centres x, B_z(step, x) and the energy density (D^2/eps0 + B^2/mu0)/2 (step, x) on
    the line y = z = 0 (the solution is y, z invariant, checked)."""
    files = sorted(
        (p for p in work.glob("pulse-*-proc*.h5") if "restart" not in p.name),
        key=lambda p: int(p.name.split("-")[1]),
    )
    if not files:
        sys.exit(f"pml_reflection: no checkpoint in {work}")
    steps, rows, dens, x = [], [], [], None
    for path in files:
        with h5py.File(path, "r") as h5:
            blocks = sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")})
            if len(blocks) != 1:
                sys.exit(
                    f"pml_reflection: {path} has {len(blocks)} blocks, the 1-D input has one"
                )
            blk = blocks[0]
            origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
            bz = h5[f"{blk}-Bz"][()][NGC:-NGC, NGC:-NGC, NGC:-NGC]  # (k, j, i)
            prof = bz[bz.shape[0] // 2, bz.shape[1] // 2, :]
            if float(np.abs(bz - prof[None, None, :]).max()) > 1e-12 * max(
                float(np.abs(prof).max()), 1e-300
            ):
                sys.exit(f"pml_reflection: {path} is not y, z invariant")
            em = (
                sum(
                    h5[f"{blk}-{v}"][()][NGC:-NGC, NGC:-NGC, NGC:-NGC] ** 2
                    for v in ("Dx", "Dy", "Dz")
                )
                / EPS0
                + sum(
                    h5[f"{blk}-{v}"][()][NGC:-NGC, NGC:-NGC, NGC:-NGC] ** 2
                    for v in ("Bx", "By", "Bz")
                )
                / MU0
            )
            dens.append(0.5 * em[em.shape[0] // 2, em.shape[1] // 2, :])
            if x is None:
                x = origin[0] + (np.arange(prof.size) + NGC + 0.5) * dxyz[0]
            steps.append(int(path.name.split("-")[1]))
            rows.append(prof)
    return np.array(steps), x, np.array(rows), np.array(dens)


def time_step(log: Path) -> float:
    """Return dt from the run log (the time of the last save over its step, the last step being shorter: second-last
    save over its step)."""
    saves = [
        ln.split("t:")[1].split(",")
        for ln in log.read_text().splitlines()
        if "save HDF5 files t:" in ln
    ]
    pairs = sorted({(int(a), float(b.split("time:")[1])) for a, b in saves})
    it, t = pairs[-2] if len(pairs) > 2 else pairs[-1]
    return t / it


def main() -> int:
    """Measure, print, optionally write."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="?")
    parser.add_argument("--from-profiles", type=Path, default=None)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--profiles", type=Path, default=None)
    args = parser.parse_args()
    if args.from_profiles:
        d = np.load(args.from_profiles)
        steps, t, x, bz = d["steps"], d["t"], d["x"], d["bz"]
        dt, e0, ewin = float(t[1] / steps[1]), 1.0, d["window_energy"]
        args.work = args.from_profiles.parent
    else:
        steps, x, bz, dens = profiles(args.work)
        dt = time_step(args.work / "run.log")
        t = steps * dt
        e0 = float(dens[0].sum())
        ewin = dens[:, np.abs(x) <= WINDOW].sum(axis=1)
    a0 = float(np.abs(bz[0]).max())
    inside = np.abs(x) <= WINDOW
    after = t >= T_IN
    win = np.abs(bz[:, inside]).max(axis=1)
    k = int(np.argmax(np.where(after, win, -1.0)))
    i = int(np.argmax(np.abs(bz[k]) * inside))
    ke = int(np.argmax(np.where(after, ewin, -1.0)))
    out = {
        "work": str(args.work),
        "A0": a0,
        "Ar": float(win[k]),
        "ratio": float(win[k]) / a0,
        "t_max": float(t[k]),
        "x_max": float(x[i]),
        "sign": float(np.sign(bz[k, i])),
        "energy_ratio": float(ewin[ke]) / e0,
        "rho_energy": (float(ewin[ke]) / e0) ** 0.5,
        "t_energy": float(t[ke]),
        "dt": dt,
        "steps": int(steps[-1]),
        "checkpoints": len(steps),
    }
    print(
        f"{args.work.name}: A0 {a0:.6e}  Ar {out['Ar']:.4e}  Ar/A0 {out['ratio']:.4e}  at t {out['t_max']:.4e} s, "
        f"x {out['x_max']:+.4f} m, sign {out['sign']:+.0f};  (Er/E0)^1/2 {out['rho_energy']:.4e} at t "
        f"{out['t_energy']:.4e} s  (dt {dt:.5e} s, {out['steps']} steps)"
    )
    if args.json:
        args.json.write_text(json.dumps(out, indent=1) + "\n")
    if args.profiles and not args.from_profiles:
        np.savez_compressed(
            args.profiles,
            steps=steps,
            t=t,
            x=x,
            bz=bz,
            window_max=win,
            window_energy=ewin / e0,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
