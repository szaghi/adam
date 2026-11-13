#!/usr/bin/env python3
"""Issue #121, P3: the convection-fed wall model a priori on the filtered averages a run keeps in its .wmf restart files (per wall column of the y_min
face: the filtered sample u, v, w, p, rho; the filtered pressure gradient; the filtered convection of cells 1-4; rho, rho u_e, rho u_eta,
rho u_e u_eta of cells 1-4). Per column: tau_w of the equilibrium model, of the pressure gradient alone, fed with the convection, with the
resolved stress removed, and the last with the velocity scale of the eddy viscosity sqrt(|tau| / rho) (the forcing kept out of it); and the
corrected model fed with the forcing averaged over the span (the same x), whose difference from the span mean of the per-column tau_w is the
bias of the nonlinear response to a fluctuating forcing. P4a: the corrected model with two velocity scales that keep the forcing's mean in the
eddy viscosity and its fluctuation out, Duprat et al.'s (2011) u_tau_p, the total stress of the column-mean forcing, and |tau| + |dp/dx| y (column-mean or filtered dp/dx);
P4a': the local modelled stress |tau + S - R| (Park & Moin's local-strain mixing length in the 1-D balance), and the same with kappa matched per
column to the LES's slope at h (cells 2 and 4: Kawai & Larsson's matching of the stress at h). Each column's x comes from its block's Morton code
and its i index, so the statistics are taken over x ranges (--ranges, the developed stretch and the recovery from the inflow); per range, the
mean forcing S, resolved stress R and S - R at the cell centres against the equilibrium tau_w, and with --dns the total shear stress of the
DNS station (Zhang, Duan & Choudhari 2018) at the same y/delta.

Without the sample's input filter (wall_model_filter 0) the .wmf holds no filtered sample: its speed is taken from the filtered rho u_e / rho
of the sample cell, its pressure the free stream's (--mach), its temperature a fraction of the wall's (--t-ratio); the variants are compared on
the same samples, so these constants shift them alike.

Usage: wmf_apriori.py <run dir> <dir of sbli_oracle.py> [--dy 1/12] [--nv 26] [--nc 3] [--mach 2.5] [--t-ratio 0.9]
                      [--ranges 20:40,45:80] [--dns M2p5_Stat.txt --delta 1.925]
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


def block_origin(code: int, ratio: int = 8) -> tuple[int, int, int, int]:
    """(i, j, k, level) of a block from its Morton code (adam_tree_object: child = ratio code + ratio + local, local = i + 2 j + 4 k bits)."""
    lev, c = 1, code
    while c >= ratio:
        c = (c - ratio) // ratio
        lev += 1
    i = j = k = 0
    c, p = code, 0
    while c >= ratio:
        local = c % ratio
        i += (local & 1) << p
        j += ((local >> 1) & 1) << p
        k += ((local >> 2) & 1) << p
        c = (c - ratio) // ratio
        p += 1
    return i, j, k, lev


def dns_total_stress(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """z/delta and the total shear stress / tau_w (mu dU/dz - rho <u'w'>) of a Zhang, Duan & Choudhari (2018) station file (M2p5: U_inf 823.6
    m/s, u_tau 40.6 m/s, T_inf 270 K, T_w 568 K, viscosity ~ T^0.76)."""
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and ln.strip()[0] in "0123456789-."]
    d = np.array([[float(t) for t in r] for r in rows if len(r) == 28])
    zod, zp, u_uinf, t_tinf, r_rhow, upwp = d[:, 1], d[:, 2], d[:, 4], d[:, 6], d[:, 7], d[:, 15]
    up = u_uinf * 823.6 / 40.6
    return zod, (t_tinf * 270.0 / 568.0) ** 0.76 * np.gradient(up, zp) - r_rhow * upwp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("oracle", type=Path)
    ap.add_argument("--dy", type=float, default=1.0 / 12.0)
    ap.add_argument("--nv", type=int, default=26)
    ap.add_argument("--nc", type=int, default=3)
    ap.add_argument("--mach", type=float, default=2.5)
    ap.add_argument("--t-ratio", type=float, default=0.9)
    ap.add_argument("--ranges", default="0:1e9", help="x ranges of the statistics, lo:hi[,lo:hi...] (P4a': 20:40,45:80)")
    ap.add_argument("--dns", type=Path, help="a DNS station file (M2p5_Stat.txt) for the total shear stress at the cell centres")
    ap.add_argument("--only-dynamic", action="store_true", help="the equilibrium model and the dynamic-kappa variant only")
    ap.add_argument("--delta", type=float, default=7.7 / 4.0, help="delta_99 / delta_i at the station (the DNS's 7.7 mm / 4 mm)")
    a = ap.parse_args()
    sys.path.insert(0, str(a.oracle))
    import sbli_oracle as so

    ini = configparser.ConfigParser(strict=False, inline_comment_prefixes=("!", ";"))
    ini.read(a.run / "input.ini")
    tref = float(ini["physics"]["reference_temperature"])
    gas = {"R": 1.0, "cp": 3.5, "pr": float(ini["physics"]["prandtl"]), "mu0": 1.0 / float(ini["physics"]["reynolds"]), "tref": tref,
           "suth": 0.0, "omega": float(ini["physics"].get("viscosity_exponent", "0.76")), "semilocal": True}
    tw = float(ini["bc_y_min"]["wall_temperature"])
    g = ini["grid"]
    x_min, lx = float(g["emin_x"]), float(g["emax_x"]) - float(g["emin_x"])
    cols, xcol = [], []
    for f in sorted(a.run.glob("tbl-restart-proc*.wmf")):
        for code, ay in read_wmf(f, a.nv):
            m = ay[10] > 0.0  # the wall columns: the filtered density of cell 1 (the resolved-stress averages) was set
            if m.any():
                ib, jb, _, lev = block_origin(code)
                if jb != 0:
                    raise ValueError(f"block {code} with wall columns is not on the y_min face (j {jb})")
                ni = ay.shape[1]
                x = x_min + (ib * ni + np.arange(ni) + 0.5) * lx / (2 ** (lev - 1) * ni)
                cols.append(ay[:, m].T)
                xcol.append(np.broadcast_to(x[:, None], m.shape)[m])
    c_all, x_all = np.concatenate(cols), np.concatenate(xcol)  # [columns, nv], [columns]
    print(f"{a.run.name}: {c_all.shape[0]} wall columns from the .wmf, x {x_all.min():.2f}-{x_all.max():.2f}")
    for rng in a.ranges.split(","):
        lo, hi = (float(t) for t in rng.split(":"))
        sel = (x_all >= lo) & (x_all < hi)
        print(f"== x {lo:g}-{hi:g}: {sel.sum()} columns")
        variants(c_all[sel], a, so, gas, ini, tw)
    if a.dns:
        zod, tot = dns_total_stress(a.dns)
        cells = ", ".join(f"y/delta {(n + 0.5) * a.dy / a.delta:.3f}: {np.interp((n + 0.5) * a.dy / a.delta, zod, tot) - 1.0:+.3f}" for n in range(a.nc))
        print(f"  DNS {a.dns.name}: total shear stress - tau_w, per tau_w, at the cell centres: {cells}")
    return 0


def variants(c: np.ndarray, a: argparse.Namespace, so, gas: dict, ini: configparser.ConfigParser, tw: float) -> None:
    """The wall model's variants on the columns `c` [columns, nv]: mean tau_w against the equilibrium model, and the mean forcing."""
    # the sample: without the input filter (wall_model_filter 0) averages 1-5 are not kept; its speed from the filtered rho u_e / rho of the
    # sample cell (cell 3: averages 19, 20), the pressure the free stream's, the temperature a constant (the variants are compared on the
    # same samples, so the constants shift them alike)
    px = c[:, 5]
    us = c[:, 19] / c[:, 18]
    ps = np.full(us.shape, 1.0 / (1.4 * a.mach**2))
    ts = np.full(us.shape, a.t_ratio * float(ini["bc_y_min"]["wall_temperature"]))
    nc = a.nc
    cc = c[:, 6 : 6 + nc]
    rv4 = c[:, 6 + 4 : 6 + 4 + 4 * 4].reshape(-1, 4, 4)  # [col, cell 1-4, (rho, rho ue, rho un, rho ue un)]
    rv = rv4[:, :nc, :]
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
    res = {"equilibrium": so.wall_model(h, us, ts, ps, tw, gas)[0]}
    fkm = np.repeat(fk.mean(axis=0, keepdims=True), us.size, axis=0)
    rkm = np.repeat(rk.mean(axis=0, keepdims=True), us.size, axis=0)
    if not a.only_dynamic:
        res.update({
            "dp/dx alone": so.wall_model(h, us, ts, ps, tw, gas, px=px, iterations=60)[0],
            "convection-fed": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, zero), iterations=60)[0],
            "  - resolved stress": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60)[0],
            "  - resolved stress, scale sqrt(|tau|)": so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60, scale="tau")[0],
        })
        # the corrected model fed with the mean forcing of the x range (in a zero-gradient layer the limit of a long time filter)
        res["  - resolved stress, the mean forcing"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fkm, rkm), iterations=60)[0]
        # #121 P4a: two velocity scales that keep the forcing's mean in the eddy viscosity and its fluctuation out: Duprat et al.'s (2011)
        # u_tau_p from the filtered wall gradient, and the total stress of the slow forcing (here the mean over all columns, the limit of a long
        # filter in a homogeneous layer: the mean over the x range); the forcing in the balance stays the column's own
        res["  - resolved stress, scale Duprat u_tau_p"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60,
                                                                         scale="duprat", px_scale=px)[0]
        res["  - resolved stress, scale |tau + S_mean|"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60,
                                                                         scale="total", s_scale=make_forcing(fkm, rkm))[0]
        res["  - resolved stress, scale |tau| + |px_mean| y"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60,
                                                                             scale="pressure", px_scale=np.full(px.shape, px.mean()))[0]
        res["  - resolved stress, scale |tau| + |px| y (filtered)"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk),
                                                                                   iterations=60, scale="pressure", px_scale=px)[0]
    # #121 P4a': the local modelled stress |tau + S - R| (Park & Moin's local-strain mixing length), the column's own and with kappa matched to
    # the LES's slope at h (the filtered u of cells 2 and 4)
    if not a.only_dynamic:
        res["  - resolved stress, scale |tau + S - R| (local)"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60,
                                                                          scale="total", s_scale=make_forcing(fk, rk))[0]
    st: dict = {}
    ue4 = np.divide(rv4[:, :, 1], rv4[:, :, 0], out=np.full(rv4[:, :, 0].shape, np.nan), where=rv4[:, :, 0] > 0.0)
    slope = np.where(np.isfinite(ue4[:, 3]), (ue4[:, 3] - ue4[:, 1]) / (2.0 * a.dy), (ue4[:, 2] - ue4[:, 1]) / a.dy)  # no cell 4: one-sided
    print(f"  the LES slope at h: central (cells 2, 4) in {100 * np.isfinite(ue4[:, 3]).mean():.1f} % of the columns, one-sided elsewhere")
    res["  - resolved stress, local, dynamic kappa"] = so.wall_model(h, us, ts, ps, tw, gas, forcing=make_forcing(fk, rk), iterations=60,
                                                                    scale="total", s_scale=make_forcing(fk, rk), slope_h=slope, stats=st)[0]
    base = res["equilibrium"].mean()
    for k, v in res.items():
        print(f"  {k:50s} mean tau_w {v.mean():.4e} ({100 * (v.mean() / base - 1):+6.1f} % on the equilibrium); negative in {100 * (v < 0).mean():5.1f} %")
    kf = st["kappa_factor"]
    print(f"  dynamic kappa / 0.41: mean {kf.mean():.2f}, median {np.median(kf):.2f}, rms {kf.std():.2f}, clipped at 0 or 10 in "
          f"{100 * ((kf <= 0.0) | (kf >= 10.0)).mean():.1f} %")
    sk = np.concatenate([np.zeros((us.size, 1)), np.cumsum(0.5 * (fk[:, 1:] + fk[:, :-1]) * np.diff(yk)[None, :], axis=1)], axis=1)[:, 1:]
    for n in range(nc):
        print(f"  cell {n + 1} (y {yk[n + 1]:.4f}): mean S {sk[:, n].mean() / base:+.3f}, R {rr[:, n].mean() / base:+.3f}, S - R "
              f"{(sk[:, n] - rr[:, n]).mean() / base:+.3f} (rms {(sk[:, n] - rr[:, n]).std() / base:.3f}) tau_w; the dp/dx part {px.mean() * yk[n + 1] / base:+.3f}")
    print(f"  px: mean {px.mean():.3e}, rms {px.std():.3e}")


if __name__ == "__main__":
    sys.exit(main())
