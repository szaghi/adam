#!/usr/bin/env python3
"""Shell-averaged kinetic-energy spectrum of a uniform periodic checkpoint (issue #119, T0.1).

E(k) of the density-weighted velocity w = sqrt(rho) u, so that the sum of E over the shells is the mean of rho |u|^2 / 2
(Parseval); shell k holds the wavenumbers with round(|k|) = k, for a box of side 2 pi L. Writes `k E(k)` as text and checks
that sum against the cells. Refined checkpoints are refused: the FFT needs one cell size.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "multirealm"))
from multirealm_oracle import basenames, last_step, load_cells  # noqa: E402


def main() -> int:
    """Compute and write the spectrum."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--step", type=int, default=None, help="default: the last checkpoint")
    parser.add_argument("--ngc", type=int, default=3)
    parser.add_argument("--out", type=Path, default=None, help="default: <work>/spectrum-<step>.dat")
    parser.add_argument("--peak-share", type=float, default=None, help="check: at least this energy share in --peak-shell")
    parser.add_argument("--peak-shell", type=int, default=2)
    args = parser.parse_args()

    names = basenames(args.work)
    step = last_step(args.work, names) if args.step is None else args.step
    variables, cells = load_cells(args.work, names, step, args.ngc)
    xs = np.array(sorted({c[0] for c in cells}))
    n = len(xs)
    if len(cells) != n**3:
        print(f"spectra: {len(cells)} cells are not a uniform {n}^3 grid  FAIL")
        return 1
    index = {x: i for i, x in enumerate(xs)}
    q = np.zeros((len(variables), n, n, n))
    for (x, y, z), val in cells.items():
        q[:, index[x], index[y], index[z]] = val
    rho = q[variables.index("r")]
    wide = np.stack([q[variables.index(m)] / np.sqrt(rho) for m in ("ru", "rv", "rw")])  # sqrt(rho) u
    k1 = np.fft.fftfreq(n, 1.0 / n)
    kx, ky, kz = np.meshgrid(k1, k1, k1, indexing="ij")
    shell = np.rint(np.sqrt(kx**2 + ky**2 + kz**2)).astype(int)
    density = 0.5 * np.sum(np.abs(np.fft.fftn(wide, axes=(1, 2, 3)) / n**3) ** 2, axis=0)
    spectrum = np.bincount(shell.ravel(), weights=density.ravel())
    mean = float(np.mean(0.5 * np.sum(wide**2, axis=0)))
    parseval = abs(float(spectrum.sum()) - mean) / mean
    out = args.out or args.work / f"spectrum-{step:09d}.dat"
    np.savetxt(out, np.column_stack([np.arange(len(spectrum)), spectrum]), header="k E(k)")
    print(f"{out.name}: {n}^3, mean K {mean:.8e}, Parseval relative {parseval:.2e}")
    ok = parseval < 1e-12
    if args.peak_share is not None:
        share = float(spectrum[args.peak_shell] / spectrum.sum())
        good = share >= args.peak_share
        ok = ok and good
        print(f"energy share in shell {args.peak_shell}: {share:.6f} (min {args.peak_share})  {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
