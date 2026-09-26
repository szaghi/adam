#!/usr/bin/env python3
"""Oracle of the GLM pulse (issue #41, M2-P4, MV-3).

Why: in 1-D the (B_x, psi) pair of mixed GLM obeys dB/dt + dpsi/dx = 0, dpsi/dt + c_h^2 dB/dx = -c_d psi
(c_d = c_h^2/c_p^2 = glm_alpha c_h / L), decoupled from the fluid. Each Fourier mode of wavenumber k evolves by the
exponential of [[0, -i k], [-i k c_h^2, -c_d]]: with c_d = 0 this is d'Alembert (two pulses at -+c_h), with c_d > 0 the
telegraph equation (an underdamped mode oscillates in an envelope exp(-c_d t / 2), a long mode is overdamped). The exact
solution of the sampled initial pulse is computed mode by mode through the FFT of the cell values, on the periodic
domain or, for reflecting walls (B_x odd, psi even, D-12), on its odd/even extension of period 2L. Checks, on the last
checkpoint (and the first one with --unchanged):

* 1-D consistency: every transverse copy of a column is bitwise identical, and blocks agree on shared columns;
* --l1-max: L1 of B_x and of psi / c_h against the exact solution, their sum below the bound (one bound, or one per
  run); the relative error (to the L1 norm of the exact solution) is reported;
* --order-min: the observed order of the L1 sum between successive resolutions (several runs) at least the bound;
* --conserved TOL: the integral of B_x (conservation history) never moves more than TOL from its first value;
* --unchanged (divergence_control = none): B_x of the last checkpoint equals the first one BITWISE (its flux is zero).

Usage:
    glm_pulse_oracle.py <work> [<work> ...] [--l1-max L [L ...]] [--order-min P] [--conserved TOL] [--unchanged]
                        [--ngc N]
"""

from __future__ import annotations

import argparse
import configparser
import math
import sys
from pathlib import Path

import h5py
import numpy as np


def read_ini(work: Path) -> configparser.ConfigParser:
    """Return the input of a work directory (its only .ini file)."""
    inis = sorted(work.glob("*.ini"))
    if len(inis) != 1:
        sys.exit(f"glm_pulse_oracle: expected one .ini in {work}, found {len(inis)}")
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.optionxform = str
    ini.read(inis[0])
    return ini


def load(work: Path, ngc: int, step: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Return the cell centres along x and the x profile of bx and psi (if saved) of the first or last checkpoint."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"glm_pulse_oracle: no checkpoint in {work}")
    steps = sorted({int(p.name.split("-")[-2]) for p in files})
    pick = steps[0] if step == "first" else steps[-1]
    cells: dict[float, dict[str, float]] = {}
    for path in (p for p in files if int(p.name.split("-")[-2]) == pick):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                names = [v for v in ("bx", "psi") if f"{blk}-{v}" in h5]
                for v in names:
                    a = h5[f"{blk}-{v}"][()].transpose(2, 1, 0)
                    if ngc:
                        a = a[ngc:-ngc, ngc:-ngc, ngc:-ngc]
                    col = a.reshape(a.shape[0], -1)
                    if not np.all(col == col[:, :1]):
                        sys.exit(f"glm_pulse_oracle: transverse copies of {v} differ in {path.name}:{blk}")
                    centres = origin[0] + (np.arange(a.shape[0]) + ngc + 0.5) * dxyz[0]
                    for c, val in zip(centres, col[:, 0], strict=True):
                        entry = cells.setdefault(float(c), {})
                        if v in entry and entry[v] != val:
                            sys.exit(f"glm_pulse_oracle: two blocks disagree on {v} at {c}")
                        entry[v] = float(val)
    xs = np.array(sorted(cells))
    prof = {v: np.array([cells[c][v] for c in xs]) for v in cells[xs[0]]}
    return xs, prof


def final_time(work: Path) -> float:
    """Return the time of the last conservation history row."""
    hist = next(work.glob("*-conservation_history.dat"))
    return float(hist.read_text().split("\n")[-2].split()[1])


def exact(ini: configparser.ConfigParser, xs: np.ndarray, t: float) -> tuple[np.ndarray, np.ndarray]:
    """Return the exact B_x and psi at the cell centres xs, time t."""
    ic = ini["initial_conditions"]
    centre, width, amp = (float(ic[k]) for k in ("pulse_center", "pulse_width", "pulse_amplitude"))
    b0 = float(ini["initial_conditions_region_1"]["bx"])
    mhd = ini["mhd"]
    ch = float(mhd["glm_ch"])
    dx = xs[1] - xs[0]
    length = mhd["glm_damping_length"].strip()
    length_value = dx if length == "min-cell" else float(length)
    cd = float(mhd["glm_alpha"]) * ch / length_value
    g = b0 + amp * np.exp(-(((xs - centre) / width) ** 2))
    wall = ini["bc_x_min"]["type"].strip() == "wall-inviscid"
    if wall:  # odd (B_x) / even (psi, zero at t = 0) extension about the wall at x = xmin, period 2 L
        if b0 != 0.0:
            sys.exit("glm_pulse_oracle: the wall extension needs a zero background B_x")
        bx0 = np.concatenate([-g[::-1], g])
    else:
        bx0 = g
    n = bx0.size
    k = 2.0 * math.pi * np.fft.fftfreq(n, d=dx)
    bh, ph = np.fft.fft(bx0), np.zeros(n, dtype=complex)
    mats = np.zeros((n, 2, 2), dtype=complex)
    mats[:, 0, 1] = -1j * k
    mats[:, 1, 0] = -1j * k * ch * ch
    mats[:, 1, 1] = -cd
    lam, vec = np.linalg.eig(mats)
    coef = np.linalg.solve(vec, np.stack([bh, ph], axis=1)[..., None])[..., 0]
    sol = np.einsum("nij,nj->ni", vec, coef * np.exp(lam * t))
    bx, psi = np.fft.ifft(sol[:, 0]).real, np.fft.ifft(sol[:, 1]).real
    if wall:
        bx, psi = bx[n // 2 :], psi[n // 2 :]
    return bx, psi


def l1(work: Path, ngc: int) -> tuple[float, float, float]:
    """Return the L1 errors of B_x and psi / c_h at the last checkpoint and the L1 norm of the exact solution."""
    ini = read_ini(work)
    xs, prof = load(work, ngc, "last")
    bx, psi = exact(ini, xs, final_time(work))
    dx = xs[1] - xs[0]
    ch = float(ini["mhd"]["glm_ch"])
    norm = float(np.sum(np.abs(bx)) * dx + np.sum(np.abs(psi)) * dx / ch)
    return float(np.sum(np.abs(prof["bx"] - bx)) * dx), float(np.sum(np.abs(prof["psi"] - psi)) * dx / ch), norm


def main() -> int:
    """Run the MV-3 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="+")
    parser.add_argument("--l1-max", type=float, nargs="+", default=None,
                        help="bound on the L1 sum of B_x and psi / c_h (one, or one per run)")
    parser.add_argument("--order-min", type=float, default=None, help="minimum observed order between resolutions")
    parser.add_argument("--conserved", type=float, default=None, metavar="TOL", help="drift bound of int B_x")
    parser.add_argument("--unchanged", action="store_true", help="B_x of the last checkpoint equals the first one")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    ok = True
    if args.unchanged:
        for w in args.work:
            _, first = load(w, args.ngc, "first")
            _, last = load(w, args.ngc, "last")
            same = np.array_equal(first["bx"], last["bx"])
            ok &= same
            print(f"{w.name}: B_x of the last checkpoint equals the first one: {'PASS (bitwise)' if same else 'FAIL'}")
        return 0 if ok else 1
    bounds = args.l1_max
    if bounds is not None and len(bounds) not in (1, len(args.work)):
        sys.exit("glm_pulse_oracle: --l1-max takes one bound or one per run")
    sums = []
    for n, w in enumerate(args.work):
        eb, ep, norm = l1(w, args.ngc)
        sums.append(eb + ep)
        line = f"{w.name}: L1 B_x {eb:.3e}, psi/c_h {ep:.3e}, sum {eb + ep:.6e} (relative {(eb + ep) / norm:.2e})"
        if bounds is not None:
            bound = bounds[0] if len(bounds) == 1 else bounds[n]
            good = eb + ep <= bound
            ok &= good
            line += f"  {'PASS' if good else 'FAIL'} (max {bound:.3e})"
        print(line)
        if args.conserved is not None:
            hist = next(w.glob("*-conservation_history.dat"))
            rows = hist.read_text().split("\n")
            col = [c.strip('"') for c in rows[0].split("=", 1)[1].split()].index("int_bx")
            vals = np.array([float(r.split()[col]) for r in rows[1:] if r.strip()])
            drift = float(np.max(np.abs(vals - vals[0])))
            good = drift <= args.conserved
            ok &= good
            print(f"   int B_x {vals[0]:.16e}, max drift over {vals.size} rows {drift:.2e}  "
                  f"{'PASS' if good else 'FAIL'} (tol {args.conserved:.1e})")
    if args.order_min is not None and len(sums) > 1:
        for a, b in zip(sums[:-1], sums[1:], strict=True):
            p = math.log2(a / b)
            good = p >= args.order_min
            ok &= good
            print(f"   observed order {p:.2f}  {'PASS' if good else 'FAIL'} (min {args.order_min:.2f})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
