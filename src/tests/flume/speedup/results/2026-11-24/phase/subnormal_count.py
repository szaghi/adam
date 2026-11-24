"""Count subnormal doubles in the interior cells of a FLUME restart, per rank and per variable (#133 night)."""

import struct
import sys
from pathlib import Path
import numpy as np

TINY = np.finfo(np.float64).tiny
base = Path(sys.argv[1])
files = sorted(base.parent.glob(base.name + "-proc*.fbd"))
tot = 0
ranks_with = 0
per = []
for f in files:
    with f.open("rb") as h:
        nv, ni, nj, nk, ngc, _a, _b = struct.unpack("=7i", h.read(28))
        (nb,) = struct.unpack("=i", h.read(4))
        n = nv * (ni + 2 * ngc) * (nj + 2 * ngc) * (nk + 2 * ngc)
        rec = 24 + 8 * n
        c = np.zeros(nv, dtype=np.int64)
        for b in range(nb):
            h.seek(32 + b * rec + 24)
            q = np.fromfile(h, dtype="<f8", count=n).reshape(
                (nk + 2 * ngc, nj + 2 * ngc, ni + 2 * ngc, nv)
            )  # Fortran (nv,i,j,k)
            q = q[ngc : ngc + nk, ngc : ngc + nj, ngc : ngc + ni, :]
            s = (q != 0) & (np.abs(q) < TINY)
            c += s.reshape(-1, nv).sum(axis=0)
    per.append((f.name, int(c.sum()), c.tolist()))
    tot += int(c.sum())
    ranks_with += int(c.sum() > 0)
per.sort(key=lambda r: -r[1])
print(f"{base}: {len(files)} ranks, subnormal interior values {tot}, ranks holding any {ranks_with}")
for name, n, c in per[:8]:
    print(f"   {name}: {n} (per variable {c})")
