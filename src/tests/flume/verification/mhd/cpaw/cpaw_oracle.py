#!/usr/bin/env python3
"""Oracle of the circularly polarised Alfven wave (issue #41, M2-P5, MV-6; Toth 2000, J. Comput. Phys. 161).

Why: after one period the exact solution is the initial state, so the error of a run is eps of Stone et al. between
its last and its first checkpoint (linear_wave_oracle.error, MV-5). A wrong flux, eigen-decomposition or nonlinear
weight leaves an error that does not converge at the design order (--order-min on the finest pair). The left and right
polarisations are mirror images under z -> -z and rounding is symmetric under a sign change, so their errors agree to
a relative --lr-tol at every N (0 in check.sh: bitwise); a sign error in a transverse component, or a z-dependent term,
breaks the agreement at O(1).

Usage:
    cpaw_oracle.py --right <work> [<work> ...] --left <work> [<work> ...] [--order-min P] [--eps-max E [E ...]]
                   [--lr-tol T] [--ngc N]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "linear-wave"))
from linear_wave_oracle import error  # noqa: E402


def main() -> int:
    """Run the MV-6 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--right", type=Path, nargs="+", required=True, help="right polarisation runs, coarse first")
    parser.add_argument("--left", type=Path, nargs="+", required=True, help="left polarisation runs, coarse first")
    parser.add_argument("--order-min", type=float, default=None, help="minimum observed order of the finest pair")
    parser.add_argument("--eps-max", type=float, nargs="+", default=None, help="bound of eps (one, or one per N)")
    parser.add_argument("--lr-tol", type=float, default=None, help="bound of |eps_right - eps_left| / eps_right")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    if len(args.right) != len(args.left):
        sys.exit("cpaw_oracle: --right and --left take the same number of runs")
    if args.eps_max is not None and len(args.eps_max) not in (1, len(args.right)):
        sys.exit("cpaw_oracle: --eps-max takes one bound or one per N")
    ok = True
    eps: dict[str, list[float]] = {"right": [], "left": []}
    for pol, works in (("right", args.right), ("left", args.left)):
        for n, w in enumerate(works):
            e = error(w, args.ngc)
            eps[pol].append(e)
            line = f"{pol:5s} {w.name}: eps {e:.4e}"
            if n > 0:
                line += f"  order {math.log2(eps[pol][-2] / e):+.2f}"
            if args.eps_max is not None:
                bound = args.eps_max[0] if len(args.eps_max) == 1 else args.eps_max[n]
                good = e <= bound
                ok &= good
                line += f"  {'PASS' if good else 'FAIL'} (max {bound:.3e})"
            print(line)
        if args.order_min is not None and len(works) > 1:
            p = math.log2(eps[pol][-2] / eps[pol][-1])
            good = p >= args.order_min
            ok &= good
            print(f"   {pol} finest-pair order {p:+.2f} {'>=' if good else '<'} {args.order_min}: "
                  f"{'PASS' if good else 'FAIL'}")
    for n, (er, el) in enumerate(zip(eps["right"], eps["left"], strict=True)):
        rel = abs(er - el) / er
        line = f"right vs left, run {n + 1}: relative difference {rel:.2e}"
        if args.lr_tol is not None:
            good = rel <= args.lr_tol
            ok &= good
            line += f"  {'PASS' if good else 'FAIL'} (max {args.lr_tol:.1e})"
        print(line)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
