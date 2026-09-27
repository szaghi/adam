#!/usr/bin/env python3
"""Oracle of the Dedner et al. (2002) div(B) peak (issue #41, M2-P4, MV-10).

Why: the peak in B_x has a non-zero div(B) from t = 0; mixed GLM must remove it, while without cleaning it is only
advected (and the scheme lets it grow). Checks on the div(B) history (`it time max_divb l1_divb seam_max_divb`), the
conservation history and the last checkpoint:

* --decay GLM [NONE]: the GLM integral of |div B| never exceeds its initial value (it is not monotone: the cleaning
  first spreads the peak into outgoing waves, a few-percent transient), ends below --decay-max times the initial
  value and, with a run without cleaning, below --ratio-max times its final value;
* --conserved TOL: the integrals of B_x, B_y, B_z (conservation history) never move more than TOL from the first row;
* --psi-max P: max |psi| of the last checkpoint below P (bounded);
* --seam: the seam-local column is zero at step 0 (the peak is far from the 2:1 seam), never above the global maximum,
  and positive at the end (the cleaning waves and the advected error reach the seam);
* --derived: the saved derived fields pt, beta, bmag equal p + |B|^2/2, 2 p/|B|^2, |B| of the saved auxiliaries, and
  the maximum of the saved divb over the interior equals the last max_divb of the history (host operator vs the
  backend kernel, relative --derived-tol).

Usage:
    divb_peak_oracle.py --decay GLM [NONE] [--decay-max F] [--ratio-max R]
    divb_peak_oracle.py <work> [<work> ...] [--conserved TOL] [--psi-max P] [--seam] [--derived] [--ngc N]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np


def history(work: Path, kind: str) -> dict[str, np.ndarray]:
    """Return the columns of the div(B) or conservation history of a run, by name."""
    path = next(work.glob(f"*-{kind}_history.dat"), None)
    if path is None:
        sys.exit(f"divb_peak_oracle: no {kind} history in {work}")
    rows = path.read_text().split("\n")
    names = [c.strip('"') for c in rows[0].split("=", 1)[1].split()]
    data = np.array([[float(x) for x in r.split()] for r in rows[1:] if r.strip()])
    return {n: data[:, k] for k, n in enumerate(names)}


def last_blocks(work: Path, ngc: int) -> list[dict[str, np.ndarray]]:
    """Return the interior datasets of every block of the last checkpoint."""
    files = sorted(p for p in work.glob("*-proc*.h5") if "restart" not in p.name)
    if not files:
        sys.exit(f"divb_peak_oracle: no checkpoint in {work}")
    last = max(int(p.name.split("-")[-2]) for p in files)
    blocks = []
    for path in (p for p in files if int(p.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                data = {}
                for key in (k for k in h5 if k.startswith(f"{blk}-")):
                    name = key.rsplit("-", 1)[1]
                    if name in ("origin", "dxdydz", "time_iteration"):
                        continue
                    a = h5[key][()]
                    data[name] = a[ngc:-ngc, ngc:-ngc, ngc:-ngc] if ngc else a
                blocks.append(data)
    return blocks


def check_decay(works: list[Path], decay_max: float, ratio_max: float) -> bool:
    """GLM integral of |div B|: bounded by its initial value, decayed, below the run without cleaning."""
    glm = history(works[0], "divb")["l1_divb"]
    ok = bool(np.max(glm) <= glm[0])
    print(f"{works[0].name}: l1 div B initial {glm[0]:.4e}, max {np.max(glm):.4e}, final {glm[-1]:.4e}  "
          f"{'PASS' if ok else 'FAIL'} (never above the initial value)")
    good = glm[-1] <= decay_max * glm[0]
    ok &= good
    print(f"   final / initial {glm[-1] / glm[0]:.3f}  {'PASS' if good else 'FAIL'} (max {decay_max})")
    if len(works) > 1:
        none = history(works[1], "divb")["l1_divb"]
        good = glm[-1] <= ratio_max * none[-1]
        ok &= good
        print(f"   {works[1].name}: l1 div B initial {none[0]:.4e}, final {none[-1]:.4e}; GLM / none "
              f"{glm[-1] / none[-1]:.3f}  {'PASS' if good else 'FAIL'} (max {ratio_max})")
    return ok


def check_conserved(work: Path, tol: float) -> bool:
    """The integrals of B_x, B_y, B_z never move more than tol from the first row."""
    hist = history(work, "conservation")
    ok = True
    parts = []
    for name in ("int_bx", "int_by", "int_bz"):
        drift = float(np.max(np.abs(hist[name] - hist[name][0])))
        ok &= drift <= tol
        parts.append(f"{name} {hist[name][0]:.6e} drift {drift:.2e}")
    print(f"   {'; '.join(parts)}  {'PASS' if ok else 'FAIL'} (tol {tol:.1e})")
    return ok


def check_psi(work: Path, psi_max: float, ngc: int) -> bool:
    """max |psi| of the last checkpoint below the bound."""
    psi = max(float(np.max(np.abs(b["psi"]))) for b in last_blocks(work, ngc))
    ok = psi <= psi_max
    print(f"   max |psi| {psi:.4e}  {'PASS' if ok else 'FAIL'} (max {psi_max:.3e})")
    return ok


def check_seam(work: Path) -> bool:
    """Seam-local column: zero at step 0, never above the global maximum, positive at the end."""
    hist = history(work, "divb")
    seam, top = hist["seam_max_divb"], hist["max_divb"]
    ok = seam[0] == 0.0 and bool(np.all(seam <= top)) and seam[-1] > 0.0
    print(f"   seam max: step 0 {seam[0]:.3e}, rows above the global max {int(np.sum(seam > top))}, final "
          f"{seam[-1]:.4e} (global {top[-1]:.4e})  {'PASS' if ok else 'FAIL'}")
    return ok


def check_derived(work: Path, tol: float, ngc: int) -> bool:
    """Saved pt, beta, bmag against their definitions; the saved divb maximum against the history."""
    err = {"pt": 0.0, "beta": 0.0, "bmag": 0.0}
    divb = 0.0
    for b in last_blocks(work, ngc):
        b2 = b["Bx"] ** 2 + b["By"] ** 2 + b["Bz"] ** 2
        err["pt"] = max(err["pt"], float(np.max(np.abs(b["pt"] - (b["p"] + 0.5 * b2)) / np.abs(b["pt"]))))
        err["beta"] = max(err["beta"], float(np.max(np.abs(b["beta"] - 2.0 * b["p"] / b2) / np.abs(b["beta"]))))
        err["bmag"] = max(err["bmag"], float(np.max(np.abs(b["bmag"] - np.sqrt(b2)) / np.abs(b["bmag"]))))
        divb = max(divb, float(np.max(np.abs(b["divb"]))))
    top = float(history(work, "divb")["max_divb"][-1])
    rel = abs(divb - top) / top
    ok = max(err.values()) <= tol and rel <= tol
    print(f"   derived fields: relative errors pt {err['pt']:.1e}, beta {err['beta']:.1e}, bmag {err['bmag']:.1e}; "
          f"saved max|divb| {divb:.10e} vs history {top:.10e} (relative {rel:.1e})  {'PASS' if ok else 'FAIL'} "
          f"(tol {tol:.0e})")
    return ok


def main() -> int:
    """Run the MV-10 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, nargs="*")
    parser.add_argument("--decay", type=Path, nargs="+", metavar="WORK", help="GLM run [and the run without cleaning]")
    parser.add_argument("--decay-max", type=float, default=0.1, help="bound of final / initial l1 div B (GLM)")
    parser.add_argument("--ratio-max", type=float, default=0.1, help="bound of GLM / none final l1 div B")
    parser.add_argument("--conserved", type=float, default=None, metavar="TOL", help="drift bound of int B")
    parser.add_argument("--psi-max", type=float, default=None, help="bound of max |psi| at the end")
    parser.add_argument("--seam", action="store_true", help="check the seam-local column")
    parser.add_argument("--derived", action="store_true", help="check the saved derived fields")
    parser.add_argument("--derived-tol", type=float, default=1.0e-12)
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    ok = True
    if args.decay:
        ok &= check_decay(args.decay, args.decay_max, args.ratio_max)
    for w in args.work:
        print(f"{w.name}:")
        if args.conserved is not None:
            ok &= check_conserved(w, args.conserved)
        if args.psi_max is not None:
            ok &= check_psi(w, args.psi_max, args.ngc)
        if args.seam:
            ok &= check_seam(w)
        if args.derived:
            ok &= check_derived(w, args.derived_tol, args.ngc)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
