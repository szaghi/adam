#!/usr/bin/env python3
"""Oracle of the div(B) at the 2:1 AMR seams of an MHD run (issue #41, M2-P7b).

Why: at a 2:1 coarse-fine face the two sides difference B with stencils of different width, so the discrete div(B) the
seam carries is not the one of either level. In PRISM (collocated FD, no damping) that seam source is persistent and
runs away at fixed h (issue #29: 6.9 -> 1002 in 100 steps); FLUME carries mixed GLM, which must keep it bounded. The
div(B) history column seam_max_divb (max |div B| over the cells beside a seam face) is the record, checked by:

* --decay-max R: the final seam value at most R times its peak (no growth after the initial transient; for a smooth
  field whose seam error is set by the initial data, e.g. the field loop);
* --ref WORK --ref-ratio-max F: the seam peak at most F times the peak of max_divb of a reference run (the uniform run
  at the base resolution): the seam must not be a larger div(B) source than the flow features themselves (for a
  shocked flow whose div(B) grows with the shocks, e.g. Orszag-Tang).

Usage:
    seam_divb_oracle.py <work> [--decay-max R] [--ref <work> --ref-ratio-max F]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "divb-peak"))
from divb_peak_oracle import history  # noqa: E402


def main() -> int:
    """Run the requested checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path)
    parser.add_argument("--decay-max", type=float, default=None, help="bound of final / peak seam div(B)")
    parser.add_argument("--ref", type=Path, default=None, help="reference run (uniform, base resolution)")
    parser.add_argument("--ref-ratio-max", type=float, default=None, help="bound of seam peak / reference peak")
    args = parser.parse_args()
    if (args.ref is None) != (args.ref_ratio_max is None):
        sys.exit("seam_divb_oracle: --ref and --ref-ratio-max go together")

    hist = history(args.work, "divb")
    seam, top, it = hist["seam_max_divb"], hist["max_divb"], hist["it"]
    if not np.any(seam > 0.0):
        sys.exit(f"seam_divb_oracle: {args.work} has no seam (seam_max_divb is zero on every row)")
    k = int(np.argmax(seam))
    print(f"{args.work.name}: seam max|div B| peak {seam[k]:.4e} at step {int(it[k])}, final {seam[-1]:.4e}; "
          f"global peak {np.max(top):.4e}, final {top[-1]:.4e}")
    ok = True
    if args.decay_max is not None:
        ratio = seam[-1] / seam[k]
        good = ratio <= args.decay_max
        ok &= good
        print(f"   final / peak {ratio:.3e}  {'PASS' if good else 'FAIL'} (max {args.decay_max:.2e})")
    if args.ref is not None:
        ref = float(np.max(history(args.ref, "divb")["max_divb"]))
        ratio = seam[k] / ref
        good = ratio <= args.ref_ratio_max
        ok &= good
        print(f"   seam peak / {args.ref.name} peak {ref:.4e} = {ratio:.3e}  {'PASS' if good else 'FAIL'} "
              f"(max {args.ref_ratio_max:.2e})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
