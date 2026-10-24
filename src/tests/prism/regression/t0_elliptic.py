#!/usr/bin/env python3
"""Assert the t0 elliptic solves of PRISM (issue #84, E3b) on the step-0 checkpoint of a run.

Modes:
  gauss <workdir> <basename>  PIC: div D - rho = -c chi to round-off, chi = (-1)^(i+j+k) over the global cells and
                              c = <rho, chi> / N. PRISM's centred div and grad of order p annihilate chi, so with PEC
                              faces (the potential odd about them) D(G) is singular on chi and the component of rho
                              along it is out of reach of any D: FLAIL drops it, everything else is matched exactly.
  open  <workdir> <basename> margin=<hs>  PIC with exact/open faces (free-space ghosts, #84 E3c): the operator is
                              nonsingular, so div D = rho to the CG tolerance (1e-11: the solve stops at 1e-12 of a
                              right-hand side that carries the ghost term), on the cells farther than `margin` from
                              the domain faces (nearer, div D reads D ghosts that PRISM's own field BC rewrote).
  divb  <workdir> <basename> [margin=<hs>]  max|div B| h / max|B| at round-off after B += curl A (div curl = 0,
                              mimetic), on the cells farther than `margin` (default 0) from the domain faces: with
                              exact/open faces PRISM's field BC copies B into the ghosts, which div B reads nearer.
  same  <dir1> <dir2> <basename> <var>...  the step-0 fields of two runs agree (rank independence).

Reads `div_D`, `rho`, `div_B` and the fields of `<basename>-000000000-proc*.h5` (saved with
save_divergence_fields = .true.), interior cells only. Exit 0 pass, 1 fail, 2 usage or missing data.
Usage: t0_elliptic.py gauss work-np1 pic_pec_t0
"""
import glob
import sys

import h5py
import numpy as np

TOL = 1.0e-12  # relative round-off bound
TOL_CG = 1.0e-11  # free-space ghosts: CG stops at 1e-12 of its right-hand side, which carries the ghost term


def blocks(workdir, basename):
    """Yield (name, origin, dxdydz, getter) for every block of the step-0 checkpoint."""
    files = sorted(glob.glob(f"{workdir}/{basename}-000000000-proc*.h5"))
    if not files:
        sys.exit(f"ERROR: no step-0 checkpoint {workdir}/{basename}-000000000-proc*.h5")
    for fn in files:
        f = h5py.File(fn, "r")
        for b in sorted({k.rsplit("-", 1)[0] for k in f.keys()}):
            # the fields are stored (z, y, x), and so are origin and dxdydz: returned here as (x, y, z)
            yield b, f[f"{b}-origin"][()][::-1], f[f"{b}-dxdydz"][()][::-1], (lambda v, f=f, b=b: f[f"{b}-{v}"][()])


def interior(a, ngc):
    return a[(slice(ngc, a.shape[0] - ngc),) * 3]


def gauss(workdir, basename, ngc):
    data = []
    for _, origin, h, get in blocks(workdir, basename):
        data.append((origin + ngc * h, h, interior(get("rho"), ngc), interior(get("div_D"), ngc)))
    x0 = np.min([d[0] for d in data], axis=0)
    num = den = 0.0
    chis = []
    for start, h, rho, div in data:
        off = np.rint((start - x0) / h).astype(int)  # global cell offsets (x, y, z); the arrays are (z, y, x)
        k, j, i = np.indices(rho.shape)
        chi = (-1.0) ** (i + j + k + off.sum())
        chis.append(chi)
        num += (rho * chi).sum()
        den += (chi * chi).sum()
    c = num / den
    rmax = max(np.abs(d[2]).max() for d in data)
    raw = max(np.abs(d[3] - d[2]).max() for d in data) / rmax
    res = max(np.abs(d[3] - d[2] + c * chi).max() for d, chi in zip(data, chis)) / rmax
    print(f"gauss: max|rho| {rmax:.4e}, c/max|rho| {c / rmax:.3e}, max|div D - rho|/max|rho| {raw:.3e}, "
          f"max|div D - rho + c chi|/max|rho| {res:.3e} (bound {TOL:.0e})")
    return res <= TOL


def gauss_open(workdir, basename, ngc, margin):
    data = []
    for _, origin, h, get in blocks(workdir, basename):
        data.append((origin + ngc * h, h, interior(get("rho"), ngc), interior(get("div_D"), ngc)))
    x0 = np.min([d[0] for d in data], axis=0)
    ncell = np.max([np.rint((d[0] - x0) / d[1]).astype(int) + np.array(d[2].shape[::-1]) for d in data], axis=0)
    rmax = max(np.abs(d[2]).max() for d in data)
    res = 0.0
    for start, h, rho, div in data:
        off = np.rint((start - x0) / h).astype(int)
        k, j, i = np.indices(rho.shape)
        g = np.stack([i + off[0], j + off[1], k + off[2]])  # global 0-based (x, y, z)
        far = np.all((g >= margin) & (g < ncell[:, None, None, None] - margin), axis=0)
        if far.any():
            res = max(res, np.abs(div - rho)[far].max())
    res /= rmax
    print(f"open: max|rho| {rmax:.4e}, max|div D - rho|/max|rho| beyond {margin} cells of the faces {res:.3e} "
          f"(bound {TOL_CG:.0e})")
    return res <= TOL_CG


def divb(workdir, basename, ngc, margin=0):
    data = [(origin + ngc * h, h, get) for _, origin, h, get in blocks(workdir, basename)]
    x0 = np.min([d[0] for d in data], axis=0)
    shape = interior(data[0][2]("div_B"), ngc).shape
    ncell = np.max([np.rint((d[0] - x0) / d[1]).astype(int) for d in data], axis=0) + np.array(shape[::-1])
    bmax = dmax = 0.0
    for start, h, get in data:
        bmax = max(bmax, *(np.abs(interior(get(v), ngc)).max() for v in ("Bx", "By", "Bz")))
        div = interior(get("div_B"), ngc)
        off = np.rint((start - x0) / h).astype(int)
        k, j, i = np.indices(div.shape)
        g = np.stack([i + off[0], j + off[1], k + off[2]])
        far = np.all((g >= margin) & (g < ncell[:, None, None, None] - margin), axis=0)
        if far.any():
            dmax = max(dmax, np.abs(div[far]).max() * h.min())
    res = dmax / bmax
    print(f"divb: max|B| {bmax:.4e}, max|div B| h / max|B| beyond {margin} cells of the faces {res:.3e} "
          f"(bound {TOL:.0e})")
    return bmax > 0.0 and res <= TOL


def same(dir1, dir2, basename, ngc, names):
    def fields(d):
        out = {}
        for _, origin, _, get in blocks(d, basename):
            for v in names:
                out[(tuple(np.round(origin, 12)), v)] = interior(get(v), ngc)
        return out
    a, b = fields(dir1), fields(dir2)
    if a.keys() != b.keys():
        print("same: the two runs have different blocks")
        return False
    scale = max(np.abs(x).max() for x in a.values())
    diff = max(np.abs(a[k] - b[k]).max() for k in a) / scale
    print(f"same: {', '.join(names)}: max difference / max value {diff:.3e} (bound {TOL:.0e})")
    return diff <= TOL


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    mode, args = sys.argv[1], sys.argv[2:]
    opts = dict(a.split("=") for a in args if "=" in a)
    ngc = int(opts.get("ngc", 6))
    args = [a for a in args if "=" not in a]
    if mode == "gauss":
        ok = gauss(args[0], args[1], ngc)
    elif mode == "open":
        ok = gauss_open(args[0], args[1], ngc, int(opts.get("margin", 3)))
    elif mode == "divb":
        ok = divb(args[0], args[1], ngc, int(opts.get("margin", 0)))
    elif mode == "same":
        ok = same(args[0], args[1], args[2], ngc, args[3:])
    else:
        sys.exit(__doc__)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
