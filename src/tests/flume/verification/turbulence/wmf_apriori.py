#!/usr/bin/env python3
"""Issue #121, P3: the convection-fed wall model a priori on the filtered averages a run keeps in its .wmf restart files (per wall column of the y_min
face: the filtered sample u, v, w, p, rho; the filtered pressure gradient; the filtered convection of cells 1-4; rho, rho u_e, rho u_eta,
rho u_e u_eta of cells 1-4). Per column: tau_w of the equilibrium model, of the pressure gradient alone, fed with the convection, with the
resolved stress removed, and the last with the velocity scale of the eddy viscosity sqrt(|tau| / rho) (the forcing kept out of it); and the
corrected model fed with the forcing averaged over the span (the same x), whose difference from the span mean of the per-column tau_w is the
bias of the nonlinear response to a fluctuating forcing.

Without the sample's input filter (wall_model_filter 0) the .wmf holds no filtered sample: its speed is taken from the filtered rho u_e / rho
of the sample cell, its pressure the free stream's (--mach), its temperature a fraction of the wall's (--t-ratio); the variants are compared on
the same samples, so these constants shift them alike.

Usage: wmf_apriori.py <run dir> <dir of sbli_oracle.py> [--dy 1/12] [--nv 26] [--nc 3] [--mach 2.5] [--t-ratio 0.9]
"""

import argparse
import configparser
import sys
from pathlib import Path

import numpy as np


def read_wmf(path: Path, nv: int) -> list[tuple[int, np.ndarray]]:
    """The y_min-side averages of each block of one .wmf: (code, a[nv, ni, nk])."""
    raw = path.read_bytes()
    hdr = np.frombuffer(raw[:20], dtype=np.int32)
    nb, ni, nj, nk = (int(v) for v in hdr[:4])
    off = 20
    out = []
    for _ in range(nb):
        code = int(np.frombuffer(raw[off : off + 8], dtype=np.int64)[0])
        off += 8
        n_x, n_y, n_z = nv * 2 * nj * nk, nv * ni * 2 * nk, nv * ni * nj * 2
        ay = np.frombuffer(raw[off + 8 * n_x : off + 8 * (n_x + n_y)], dtype=np.float64).reshape((nv, ni, 2, nk), order="F")
        off += 8 * (n_x + n_y + n_z)
        out.append((code, ay[:, :, 0, :]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("oracle", type=Path)
    ap.add_argument("--dy", type=float, default=1.0 / 12.0)
    ap.add_argument("--nv", type=int, default=26)
    ap.add_argument("--nc", type=int, default=3)
    ap.add_argument("--mach", type=float, default=2.5)
    ap.add_argument("--t-ratio", type=float, default=0.9)
    a = ap.parse_args()
    sys.path.insert(0, str(a.oracle))
    import sbli_oracle as so

    ini = configparser.ConfigParser(strict=False, inline_comment_prefixes=("!", ";"))
    ini.read(a.run / "input.ini")
    tref = float(ini["physics"]["reference_temperature"])
    gas = {"R": 1.0, "cp": 3.5, "pr": float(ini["physics"]["prandtl"]), "mu0": 1.0 / float(ini["physics"]["reynolds"]), "tref": tref,
           "suth": 0.0, "omega": float(ini["physics"].get("viscosity_exponent", "0.76")), "semilocal": True}
    tw = float(ini["bc_y_min"]["wall_temperature"])
    cols = []
    for f in sorted(a.run.glob("tbl-restart-proc*.wmf")):
        for _, ay in read_wmf(f, a.nv):
            m = ay[10] > 0.0  # the wall columns: the filtered density of cell 1 (the resolved-stress averages) was set
            if m.any():
                cols.append(ay[:, m].T)
    c = np.concatenate(cols)  # [columns, nv]
    print(f"{a.run.name}: {c.shape[0]} wall columns from the .wmf")
    # the sample: without the input filter (wall_model_filter 0) averages 1-5 are not kept; its speed from the filtered rho u_e / rho of the
    # sample cell (cell 3: averages 19, 20), the pressure the free stream's, the temperature a constant (the variants are compared on the
    # same samples, so the constants shift them alike)
    px = c[:, 5]
    us = c[:, 19] / c[:, 18]
    ps = np.full(us.shape, 1.0 / (1.4 * a.mach**2))
    ts = np.full(us.shape, a.t_ratio * float(ini["bc_y_min"]["wall_temperature"]))
    nc = a.nc
    cc = c[:, 6 : 6 + nc]
    rv = c[:, 6 + 4 : 6 + 4 + 4 * 4].reshape(-1, 4, 4)[:, :nc, :]  # [col, cell, (rho, rho ue, rho un, rho ue un)]
    rr = rv[:, :, 3] - rv[:, :, 1] * rv[:, :, 2] / rv[:, :, 0]
    h = np.full(us.shape, 2.5 * a.dy)
    yk = np.concatenate([[0.0], (np.arange(1, nc + 1) - 0.5) * a.dy])

    def make_forcing(f_knots: np.ndarray, r_knots: np.ndarray):
        def forcing(ym: np.ndarray) -> np.ndarray:
            out = np.zeros_like(ym)
            for n in range(ym.shape[0]):
                fine = np.linspace(0.0, float(ym[n, -1]), 801)
                f = np.interp(fine, yk, f_knots[n])
                cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(fine))])
                out[n] = np.interp(ym[n], fine, cum) - np.interp(ym[n], yk, r_knots[n])
            return out
        return forcing

    fk = np.concatenate([np.zeros((us.size, 1)), cc], axis=1) + px[:, None]
    rk = np.concatenate([np.zeros((us.size, 1)), rr], axis=1)
    zero = np.zeros_like(rk)
    res = {
        "equilibrium": so.wall_model(h, us, ts, ps, tw, gas)[0],
        "dp/dx alone": so.wall_model(h, us, ts, ps, tw, gas, px=px, iterations=60)[0],
        "convection-fed": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, zero), iterations=60)[0],
        "  - resolved stress": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60)[0],
        "  - resolved stress, scale sqrt(|tau|)": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60, scale="tau")[0],
    }
    # the corrected model fed with the span-mean forcing (columns grouped by their filtered sample's x is not in the .wmf: use the mean over
    # all columns, a zero-gradient layer is homogeneous enough for the comparison of the means)
    fkm = np.repeat(fk.mean(axis=0, keepdims=True), us.size, axis=0)
    rkm = np.repeat(rk.mean(axis=0, keepdims=True), us.size, axis=0)
    res["  - resolved stress, the mean forcing"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fkm, rkm), iterations=60)[0]
    base = res["equilibrium"].mean()
    for k, v in res.items():
        print(f"  {k:45s} mean tau_w {v.mean():.4e} ({100 * (v.mean() / base - 1):+6.1f} % on the equilibrium); negative in {100 * (v < 0).mean():5.1f} %")
    S_h = make_forcing(fk, rk)(np.tile(np.array([[0.0, 2.5 * a.dy]]), (us.size, 1)))[:, 1]
    print(f"  forcing at h (corrected), S(h): mean {S_h.mean():.3e}, rms {S_h.std():.3e}; |S(h)| mean {np.abs(S_h).mean():.3e}; equilibrium tau_w mean {base:.3e}")
    print(f"  convection C cells 1-{nc}: mean {cc.mean(axis=0)}, rms {cc.std(axis=0)}; resolved stress R: mean {rr.mean(axis=0)}, rms {rr.std(axis=0)}")
    print(f"  px: mean {px.mean():.3e}, rms {px.std():.3e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
