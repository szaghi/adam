#!/usr/bin/env python3
"""The shock/turbulent boundary layer interaction of issue #121 S4 (Mach 2.28, flow deflection 8 deg, the IUSTI experiment): its length scales and
wall pressure, from the checkpoints and the wall-model history of a run of make_tbl.py with --top-shock.

`compare <work> --x-imp X --deflection T [--t0 T0] [--json out]`: the checkpoints from T0 on, averaged over the span and in time (tbl_oracle.averages),
and the wall history over the same window:
- the incoming layer: delta(x), delta*, theta and H against the local edge state (tbl_oracle.edge_shape) and the nominal free stream; delta_0, the
  99 % thickness the undisturbed layer would have at x_imp, extrapolated from a linear fit of delta over [x_imp - 15, x_imp - 6] (upstream of the
  interaction's influence); delta* the same way (the Souverein et al. 2013 scaling takes the upstream displacement thickness); Re_theta of both;
- C_f(x) from the wall model's tau_w: the separation (the first C_f < 0 downstream of x_imp - 10 delta_0) and the reattachment (the last), L_sep;
  the same for the time-mean u of the first cell (the resolved flow's separation, which tau_w misses below the sample);
- the mean wall pressure (the first cell), as (p_w - p_inf) / (p_2 - p_inf) with p_2 the pressure behind the incident shock (the normalisation
  of the STREAmS figure), against x* = (x - x_imp) / L;
- the interaction length L_int = x_imp - x_ref, x_ref the apparent origin of the reflected shock: above the layer every row crosses the incident
  shock (p from p_1 to p_2) and then the reflected one (p_2 to p_3, the exact regular reflection), the mid-pressure crossings of the reflected one
  over the rows y in [2.5, 5] delta_0 fitted by a line and extrapolated to the wall; the incident one fitted the same way, its foot against x_imp;
- the Souverein scaling L* = L_int / (delta* g3), g3 = sin(beta - phi) / (sin beta sin phi), and the separation criterion Se = k1 (p_3 - p_1) / q_inf
  with k1 = 3 (low Reynolds number, Volpiani et al. 2018).
Recorded, not judged: the comparison with the digitised references (reference/sbli_m2p28_*) is printed next to each value: the STREAmS DNS and the
IUSTI points of the CPC paper's Fig. 7, whose x* is measured from the reflected shock (shifted here by SHIFT_FIG7 = 1 to x_imp; see `references`), and
the STREAmS DNS at Re_theta 6882 of Bernardini et al. (arXiv:2204.10531, Fig. 4a), x* from x_imp.

`means <work> --x-imp X [--t0 T0] --out <npz>`: the span and time averages near the wall (tbl_oracle.averages, the first 12 rows of the finest grid),
the C_f of the wall-model history, delta_0 and the run's gas and viscosity law, in one file for `apriori` (run it where the checkpoints are).
`fields <work> --x-imp X --deflection T [--t0 T0] --out <npz>`: the span and time averages of every row (tbl_oracle.averages), delta(x), delta_0
and the free-stream temperature, for `compression` (run it where the checkpoints are).
`compression <npz>... [--json out]`: where the interaction's compression starts and how it crosses the layer, per run: the incoming layer's mean
Mach number per cell and sonic height against the DNS profile the inflow is made from (reference/tbl_m2p28_profile.txt), the x where the wall
pressure reaches 5, 10 and 50 % of p_2 - p_inf against the Re_theta 6882 DNS, and the path of the compression through the layer (per height, the
first crossing of p_inf + k (p_2 - p_inf), k = 0.1 and 0.5), and the mean reversed flow of the resolved cells 1-3 (the mean u < 0, its streamwise
extent: the wall model's tau_w, from which `compare` takes L_sep, does not see a reversal below its sample).
`apriori <npz> [--hwm H] [--json out]`: the equilibrium wall model of the run (a port of `wall_model_solve`, adam_flume_wall_model_library.F90) applied
to the mean field at the run's sample (the 3rd cell of the local grid, or the height H between cells as `sample_cells` takes it) per x, with and without
the mean streamwise pressure gradient in the momentum balance, `(mu + mu_t) dU/dy = tau_w + (dp/dx) y`, the eddy viscosity on the velocity scale
`sqrt((|tau_w| + |dp/dx| y) / rho)` (S4d, issue #121): C_f against the history's (the history averages tau of the instantaneous samples, the model
here takes the mean sample), and the C_f < 0 region the pressure gradient opens, if any. A priori: what a run with the term would start from, not
what it would converge to.
`reversal <work> --x-imp X --delta0 D [--t0 T0] [--json out]`: the reversed flow of each checkpoint near the interaction, x in
[x_imp - 4.5 D, x_imp + D], y below 1.5 D, on the finest grid (coarser cells repeated): the area with u < 0 per unit span (the volume over the
span width, in D^2), the height of the highest reversed cell of each (x, z) column (its 50th, 90th and 99th percentiles over the columns that
reverse at all, in D), the fraction of the first and of the third cells (the default sample) reversed, and the streamwise extent of the span-averaged u < 0 in the first cell. For the
instantaneous bubble of the DNS (Bernardini et al. 2022, Fig. 10: in its large state up to y* = y / L_int of about 0.15, 0.5 delta_0, over
x* from about -0.85 to 0.1), which the time means of `compare` do not show (run it where the checkpoints are).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "boundary"))
import tbl_oracle as to  # noqa: E402

GAMMA = 1.4


SHIFT_FIG7 = 1.0  # STREAmS Fig. 7 (CPC 2021) plots x* from the reflected shock's apparent origin (the IUSTI convention): x*(x_imp) = x*(fig) - 1


def references(p21: float) -> dict[str, np.ndarray]:
    """The digitised wall-pressure references as (x*, (p_w - p_inf) / (p_2 - p_inf)), x* = (x - x_imp) / L_int: `dns` and `iusti`, STREAmS Fig. 7
    (CPC 2021) shifted by SHIFT_FIG7 (its origin is the reflected shock, reference/sbli_m2p28_fig7a_*.txt); `dns_hre` and `exp_hre`, the STREAmS
    DNS at Re_theta 6882 and the experiment of Bernardini et al. (arXiv:2204.10531, Fig. 4a), from p_w / p_inf with p2/p1 = `p21`."""
    ref = HERE / "reference"
    out = {}
    for key, name in (("dns", "sbli_m2p28_fig7a_dns.txt"), ("iusti", "sbli_m2p28_fig7a_iusti.txt")):
        a = np.loadtxt(ref / name)
        out[key] = np.column_stack([a[:, 0] - SHIFT_FIG7, a[:, 1]])
    for key, name in (("dns_hre", "sbli_m2p28_bernardini2022_fig4a_dns.txt"), ("exp_hre", "sbli_m2p28_bernardini2022_fig4a_exp.txt")):
        a = np.loadtxt(ref / name)
        out[key] = np.column_stack([a[:, 0], (a[:, 1] - 1.0) / (p21 - 1.0)])
    return out


def shocks(mach: float, theta_deg: float) -> dict:
    """The incident oblique shock and the regular reflection off the wall (shockgen_oracle.reflection, the same gamma)."""
    import shockgen_oracle as so

    return so.reflection(mach, theta_deg)


def crossings(p_row: np.ndarray, x: np.ndarray, level: float, after: float) -> float:
    """The first x > `after` where the row's pressure rises through `level` (linear between cells), NaN if none."""
    k = np.nonzero((p_row[:-1] < level) & (p_row[1:] >= level) & (x[:-1] > after))[0]
    if len(k) == 0:
        return float("nan")
    j = k[0]
    return float(x[j] + (level - p_row[j]) / (p_row[j + 1] - p_row[j]) * (x[j + 1] - x[j]))


def shock_fits(p: np.ndarray, x: np.ndarray, y: np.ndarray, delta0: float, x_imp: float, p_inf: float, p2: float, p3: float,
               beta_deg: float) -> tuple[float, float, float, tuple[list, list, list, list]]:
    """The incident and the reflected shock above the layer: the mid-pressure crossings of the rows 2.5 <= y / delta_0 <= 5 fitted by lines, their
    feet at the wall (the incident one's, the reflected one's apparent origin x_ref), the reflected one's angle to the wall and the crossings."""
    dx = x[1] - x[0]
    pm_i, pm_r = 0.5 * (p_inf + p2), 0.5 * (p2 + p3)
    tb = math.tan(math.radians(beta_deg))
    rows = [j for j in range(len(y)) if 2.5 * delta0 <= y[j] <= 5.0 * delta0]
    xi, xr, yr, yi = [], [], [], []
    for j in rows:
        xa = crossings(p[:, j], x, pm_i, x_imp - (y[j] / tb) - 10.0)
        if np.isfinite(xa):
            xi.append(xa)
            yi.append(y[j])
            xb_ = crossings(p[:, j], x, pm_r, xa + 2.0 * dx)
            if np.isfinite(xb_):
                xr.append(xb_)
                yr.append(y[j])
    ci = np.polyfit(yi, xi, 1) if len(yi) > 3 else [np.nan, np.nan]
    cr = np.polyfit(yr, xr, 1) if len(yr) > 3 else [np.nan, np.nan]
    beta_r = math.degrees(math.atan(1.0 / cr[0])) if len(yr) > 3 else float("nan")
    return float(ci[1]), float(cr[1]), beta_r, (yi, xi, yr, xr)


def compare(work: Path, x_imp: float, theta_deg: float, t0: float, json_out: str | None) -> int:
    """Print the length scales and the wall pressure of an interaction run (and write them as JSON)."""
    cps = [(t, g) for t, g in to.checkpoints(work) if t >= t0]
    if not cps:
        print(f"sbli_oracle: no checkpoint at t >= {t0} in {work}")
        return 1
    mach = 1.0 / math.sqrt(GAMMA * float(to.ini_value(work, "reference_temperature")))
    sh = shocks(mach, theta_deg)
    p_inf = 1.0 / (GAMMA * mach**2)
    p2, p3 = p_inf * sh["incident"]["p21"], p_inf * sh["p31"]
    a = to.averages(work, [g for _, g in cps])
    lx, ly = float(to.ini_value(work, "emax_x")), float(to.ini_value(work, "emax_y"))
    nx, ny = a["u"].shape
    dx, dy = lx / nx, ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    uf = a["ru"] / a["r"]
    p = a["t"] * a["r"]
    delta = to.edge_height(a["u"], y, dy)
    # the nominal integrals stop at 1.5 delta as the edge ones do: above the layer the post-shock region (rho u > 1) would count
    inner = y[None, :] < 1.5 * np.where(np.isfinite(delta), delta, 0.0)[:, None]
    theta = np.where(inner, a["r"] * uf * (1.0 - uf), 0.0).sum(axis=1) * dy
    dstar_nom = np.where(inner, 1.0 - a["r"] * uf, 0.0).sum(axis=1) * dy
    t_inf = float(to.ini_value(work, "reference_temperature"))
    h_edge, rho_e, p_e = to.edge_shape(a, y, delta, t_inf)
    # the edge-normalised displacement thickness: delta* = H theta_e, with theta_e from the same integrals
    dstar_e = np.full(nx, np.nan)
    theta_e = np.full(nx, np.nan)
    for i in range(nx):
        if not np.isfinite(delta[i]):
            continue
        top = (y >= 1.5 * delta[i]) & (y <= 2.5 * delta[i])
        if not top.any():
            continue
        re_, ue = a["r"][i][top].mean(), uf[i][top].mean()
        m = y < 1.5 * delta[i]
        rn, un = a["r"][i][m] / re_, uf[i][m] / ue
        dstar_e[i] = (1.0 - rn * un).sum() * dy
        theta_e[i] = (rn * un * (1.0 - un)).sum() * dy
    up = (x >= x_imp - 15.0) & (x <= x_imp - 6.0) & np.isfinite(delta)
    fit = lambda q: float(np.polyval(np.polyfit(x[up], q[up], 1), x_imp))  # noqa: E731
    delta0, dstar0, dstar0_nom, theta0, theta0_e = fit(delta), fit(dstar_e), fit(dstar_nom), fit(theta), fit(theta_e)
    reynolds = float(to.ini_value(work, "reynolds"))
    re_theta_nom, re_theta_edge = reynolds * theta0, reynolds * theta0_e  # the edge one with rho_e u_e / mu_inf, a first-order correction
    # C_f from the wall model's tau_w over the same window (rho_inf U_inf^2 = 1)
    w = to.wall_rows(work / "tbl-wall_y_min_history.dat")
    w = w[w[:, 1] >= t0]
    nb = (w.shape[1] - 2) // 2
    xb = (np.arange(nb) + 0.5) * lx / nb
    cf = 2.0 * np.interp(x, xb, w[:, 2 : 2 + nb].mean(axis=0))
    neg = np.nonzero((cf < 0.0) & (x > x_imp - 10.0 * delta0))[0]
    x_sep = float(x[neg[0]]) if len(neg) else float("nan")
    x_att = float(x[neg[-1]]) if len(neg) else float("nan")
    l_sep = x_att - x_sep if len(neg) else 0.0
    # the separation of the resolved flow: the time-mean u of the first cell (finest row) below zero, which tau_w misses when the reversed
    # layer lies below the wall model's sample (#121 P2)
    u1 = a["u"][:, 0]
    neg1 = np.nonzero((u1 < 0.0) & (x > x_imp - 10.0 * delta0))[0]
    x_sep1 = float(x[neg1[0]]) if len(neg1) else float("nan")
    x_att1 = float(x[neg1[-1]]) if len(neg1) else float("nan")
    l_sep1 = x_att1 - x_sep1 + (x[1] - x[0]) if len(neg1) else 0.0
    # the reflected and the incident shock, fitted above the layer
    x_inc_wall, x_ref, beta_r, (yi, xi, yr, xr) = shock_fits(p, x, y, delta0, x_imp, p_inf, p2, p3, sh["incident"]["beta"])
    l_int = x_imp - x_ref
    b, ph = math.radians(sh["incident"]["beta"]), math.radians(theta_deg)
    g3 = math.sin(b - ph) / (math.sin(b) * math.sin(ph))
    l_star = l_int / (dstar0 * g3)
    se = 3.0 * (p3 - p_inf) / (0.5 * 1.0 * 1.0)  # k1 = 3, q_inf = rho_inf U_inf^2 / 2 = 1 / 2
    p_w = p[:, 0]
    pw_norm = (p_w - p_inf) / (p2 - p_inf)
    x_star = (x - x_imp) / l_int
    ref = HERE / "reference"
    refs = references(sh["incident"]["p21"])
    dns, iusti, dns_hre = refs["dns"], refs["iusti"], refs["dns_hre"]
    vol = {}
    for line in (ref / "sbli_m2p28_volpiani2018_fig12.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        pan, wall, phi, val = line.split()
        if wall == "adiabatic" and abs(float(phi) - theta_deg) < 0.05:
            vol[pan] = float(val)
    print(f"== {work.name}: {len(cps)} checkpoints from t = {cps[0][0]:.1f} to {cps[-1][0]:.1f}, {nx} x {ny} cells; Mach {mach:.3f}, deflection {theta_deg:g} deg")
    print(f"shocks (exact): beta {sh['incident']['beta']:.3f} deg, p2/p1 {sh['incident']['p21']:.4f}; reflected {sh['angle']:.3f} deg to the wall, p3/p1 {sh['p31']:.4f}")
    print(f"incoming layer at x_imp = {x_imp:g} (fit over [{x_imp - 15:g}, {x_imp - 6:g}]): delta_0 {delta0:.3f}, delta* {dstar0:.4f} (nominal {dstar0_nom:.4f}), "
          f"theta {theta0:.4f}, Re_theta {re_theta_nom:.0f} (against the edge state {re_theta_edge:.0f}); H at x_imp - 10 {h_edge[int((x_imp - 10) / dx)]:.2f}")
    print(f"wall pressure upstream: p_w / p_inf {np.mean(p_w[(x > x_imp - 25.0) & (x < x_imp - 15.0)]) / p_inf:.4f} over [{x_imp - 25:g}, {x_imp - 15:g}]")
    span = f"from x = {x_sep:.2f} to {x_att:.2f}: L_sep / delta_0 {l_sep / delta0:.2f}" if l_sep > 0.0 else "nowhere: L_sep 0"
    print(f"separation: C_f < 0 {span} (Volpiani et al. adiabatic {theta_deg:g} deg, Re_theta 1050: {vol.get('a', float('nan')):.2f})")
    span = f"from x = {x_sep1:.2f} to {x_att1:.2f}: L_sep / delta_0 {l_sep1 / delta0:.2f}" if l_sep1 > 0.0 else "nowhere: L_sep 0"
    print(f"separation of the resolved flow: mean u < 0 in the first cell {span}")
    print(f"incident shock: foot at x = {x_inc_wall:.2f} (x_imp {x_imp:g}); reflected shock: {len(yr)} rows, {beta_r:.2f} deg to the wall, apparent origin x_ref {x_ref:.2f}")
    print(f"L_int = x_imp - x_ref {l_int:.3f}: L_int / delta_0 {l_int / delta0:.2f} (Volpiani et al. {vol.get('b', float('nan')):.2f}); Souverein L* = L_int / (delta* g3) "
          f"{l_star:.3f}, g3 {g3:.3f}, Se {se:.3f}")
    for xs in (-1.0, -0.5, 0.0, 0.5, 1.0, 2.0):
        print(f"  x* {xs:+.1f}: (p_w - p_inf) / (p_2 - p_inf) {np.interp(xs, x_star, pw_norm):.3f} (DNS {np.interp(xs, dns[:, 0], dns[:, 1]):.3f}, "
              f"Re_theta 6882 {np.interp(xs, dns_hre[:, 0], dns_hre[:, 1]):.3f})")
    if json_out:
        out = {"x": x.tolist(), "cf": cf.tolist(), "pw_norm": pw_norm.tolist(), "x_star": x_star.tolist(), "delta": [None if not np.isfinite(v) else float(v) for v in delta],
               "H": [None if not np.isfinite(v) else float(v) for v in h_edge], "p_e": [None if not np.isfinite(v) else float(v) for v in p_e],
               "delta0": delta0, "dstar0": dstar0, "dstar0_nominal": dstar0_nom, "theta0": theta0, "re_theta": re_theta_nom, "re_theta_edge": re_theta_edge,
               "x_sep": x_sep, "x_att": x_att, "l_sep": l_sep, "x_sep_u1": x_sep1, "x_att_u1": x_att1, "l_sep_u1": l_sep1, "u1": u1.tolist(), "x_ref": x_ref, "l_int": l_int, "x_inc_wall": x_inc_wall, "beta_r": beta_r, "l_star": l_star,
               "g3": g3, "se": se, "x_imp": x_imp, "t": [t for t, _ in cps], "volpiani": vol, "iusti": iusti.tolist(), "dns": dns.tolist(),
               "dns_hre": dns_hre.tolist(), "exp_hre": refs["exp_hre"].tolist(),
               "shock_rows": {"yi": yi, "xi": xi, "yr": yr, "xr": xr}}
        Path(json_out).write_text(json.dumps(out))
    return 0


def means(work: Path, x_imp: float, t0: float, out: Path) -> int:
    """Write the near-wall span and time averages of a run, its wall-model C_f, delta_0 and its gas and viscosity law to `out` (an npz)."""
    cps = [(t, g) for t, g in to.checkpoints(work) if t >= t0]
    if not cps:
        print(f"sbli_oracle: no checkpoint at t >= {t0} in {work}")
        return 1
    a = to.averages(work, [g for _, g in cps])
    lx, ly = float(to.ini_value(work, "emax_x")), float(to.ini_value(work, "emax_y"))
    nx, ny = a["u"].shape
    dx, dy = lx / nx, ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    delta = to.edge_height(a["u"], y, dy)
    up = (x >= x_imp - 15.0) & (x <= x_imp - 6.0) & np.isfinite(delta)
    delta0 = float(np.polyval(np.polyfit(x[up], delta[up], 1), x_imp))
    w = to.wall_rows(work / "tbl-wall_y_min_history.dat")
    w = w[w[:, 1] >= t0]
    nb = (w.shape[1] - 2) // 2
    cf = 2.0 * np.interp(x, (np.arange(nb) + 0.5) * lx / nb, w[:, 2 : 2 + nb].mean(axis=0))
    t_inf = float(to.ini_value(work, "reference_temperature"))
    try:
        suth = float(to.ini_value(work, "sutherland_temperature"))
    except (KeyError, ValueError):
        suth = 0.0
    near = slice(0, 12)
    np.savez_compressed(out, x=x, dy=dy, cf=cf, delta0=delta0, x_imp=x_imp, t0=cps[0][0], t1=cps[-1][0], n=len(cps), t_inf=t_inf, suth=suth,
                        reynolds=float(to.ini_value(work, "reynolds")), prandtl=float(to.ini_value(work, "prandtl")),
                        t_wall=float(to.ini_value(work, "wall_temperature")), semilocal=to.ini_value(work, "wall_model_damping") == "semi-local",
                        **{k: a[k][:, near] for k in ("r", "u", "w", "t", "p", "s")})
    print(f"{out}: {len(cps)} checkpoints t {cps[0][0]:.1f}-{cps[-1][0]:.1f}, {nx} x {ny} cells, delta_0 {delta0:.4f}")
    return 0


def fields(work: Path, x_imp: float, theta_deg: float, t0: float, out: Path) -> int:
    """Write the full span and time averages of a run (tbl_oracle.averages, every row), delta_0 and the shocks' pressures to `out` (an npz), for
    `compression` (run it where the checkpoints are)."""
    cps = [(t, g) for t, g in to.checkpoints(work) if t >= t0]
    if not cps:
        print(f"sbli_oracle: no checkpoint at t >= {t0} in {work}")
        return 1
    a = to.averages(work, [g for _, g in cps])
    lx, ly = float(to.ini_value(work, "emax_x")), float(to.ini_value(work, "emax_y"))
    nx, ny = a["u"].shape
    dx, dy = lx / nx, ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    delta = to.edge_height(a["u"], y, dy)
    up = (x >= x_imp - 15.0) & (x <= x_imp - 6.0) & np.isfinite(delta)
    delta0 = float(np.polyval(np.polyfit(x[up], delta[up], 1), x_imp))
    t_inf = float(to.ini_value(work, "reference_temperature"))
    np.savez_compressed(out, x=x, y=y, delta=delta, delta0=delta0, x_imp=x_imp, theta=theta_deg, t_inf=t_inf, t0=cps[0][0], t1=cps[-1][0], n=len(cps),
                        **{k: a[k] for k in a})
    print(f"{out}: {len(cps)} checkpoints t {cps[0][0]:.1f}-{cps[-1][0]:.1f}, {nx} x {ny} cells, delta_0 {delta0:.4f}")
    return 0


def sonic_height(m_col: np.ndarray, y: np.ndarray) -> float:
    """The lowest height where the mean Mach number of a column reaches 1 (linear between cells; the wall, M = 0, below the first cell)."""
    yy, mm = np.concatenate([[0.0], y]), np.concatenate([[0.0], m_col])
    k = np.nonzero(mm >= 1.0)[0]
    if len(k) == 0:
        return float("nan")
    j = k[0]
    return float(yy[j - 1] + (1.0 - mm[j - 1]) / (mm[j] - mm[j - 1]) * (yy[j] - yy[j - 1]))


def compression(npzs: list[Path], json_out: str | None) -> int:
    """Where the compression of the interaction starts and how it crosses the layer, per run (an npz of `fields`): the incoming layer's mean Mach
    profile and sonic height against the DNS profile the inflow is made from, the wall-pressure onset, and the path of the compression through the
    layer (the first crossing of p_inf + k (p_2 - p_inf) per row)."""
    ref = np.loadtxt(HERE / "reference" / "tbl_m2p28_profile.txt")
    res = {}
    for npz in npzs:
        d = np.load(npz)
        x, y, delta0, x_imp = d["x"], d["y"], float(d["delta0"]), float(d["x_imp"])
        dy, dx_ = y[1] - y[0], x[1] - x[0]
        t_inf = float(d["t_inf"])
        mach = 1.0 / math.sqrt(GAMMA * t_inf)
        sh = shocks(mach, float(d["theta"]))
        p_inf = 1.0 / (GAMMA * mach**2)
        p2, p3 = p_inf * sh["incident"]["p21"], p_inf * sh["p31"]
        p = d["p"]
        uf = d["ru"] / d["r"]
        mbar = uf / np.sqrt(GAMMA * d["t"])
        _, x_ref, _, _ = shock_fits(p, x, y, delta0, x_imp, p_inf, p2, p3, sh["incident"]["beta"])
        hre = references(sh["incident"]["p21"])["dns_hre"]
        l_int = x_imp - x_ref
        ys = np.array([sonic_height(mbar[i], y) for i in range(len(x))])
        # the DNS profile: y / delta of its own station, M = U / sqrt(T) M_inf
        m_ref = ref[:, 1] / np.sqrt(ref[:, 2]) * mach
        ys_ref = sonic_height(m_ref[1:], ref[1:, 0])
        up = {}
        for k in (5.0, 10.0):
            i = int(np.argmin(np.abs(x - (x_imp - k * delta0))))
            up[f"x_imp-{k:g}d0"] = {"x": float(x[i]), "y_sonic_d0": ys[i] / delta0, "y_sonic_cells": ys[i] / dy, "delta_d0": float(d["delta"][i]) / delta0,
                                    "M_cells": [float(v) for v in mbar[i, :4]], "u_cells": [float(v) for v in uf[i, :4]]}
        pw = (p[:, 0] - p_inf) / (p2 - p_inf)
        onset = {}
        for lev in (0.05, 0.1, 0.5):
            x0 = crossings(pw, x, lev, x_imp - 20.0)
            onset[f"{lev:g}"] = {"x": x0, "upstream_d0": (x_imp - x0) / delta0, "x_star": (x0 - x_imp) / l_int,
                                 "x_star_dns": crossings(hre[:, 1], hre[:, 0], lev, -3.0)}
        # the path of the compression: per height, the first crossing of p_inf + k (p_2 - p_inf), searched from x_imp - 20
        path = {}
        for lev in (0.1, 0.5):
            lv = p_inf + lev * (p2 - p_inf)
            pts = []
            for yd in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
                j = 0 if yd == 0.0 else int(np.argmin(np.abs(y - yd * delta0)))
                xc = crossings(p[:, j], x, lv, x_imp - 20.0)
                pts.append({"y_d0": float(y[j]) / delta0, "x": xc, "x_star": (xc - x_imp) / l_int, "upstream_d0": (x_imp - xc) / delta0,
                            "M_up": float(mbar[int(np.argmin(np.abs(x - (x_imp - 5.0 * delta0)))), j])})
            path[f"{lev:g}"] = pts
        win = (x > x_imp - 4.0 * delta0) & (x < x_imp + 2.0 * delta0)
        # the mean reversed flow of the resolved cells 1-3 (the mean u, not the wall model's tau_w): its streamwise extent near the interaction
        near = (x > x_imp - 10.0 * delta0) & (x < x_imp + 4.0 * delta0)
        rev = {}
        for j in range(3):
            neg = x[near & (d["u"][:, j] < 0.0)]
            rev[f"cell{j + 1}"] = {"x0": float(neg.min()), "x1": float(neg.max()), "extent_d0": float((neg.max() - neg.min() + dx_) / delta0)} if neg.size else None
        res[npz.stem] = {"t": [float(d["t0"]), float(d["t1"])], "n": int(d["n"]), "delta0": delta0, "dy_d0": dy / delta0, "l_int": l_int,
                         "l_int_d0": l_int / delta0, "upstream": up, "y_sonic_dns_delta": ys_ref, "onset": onset, "path": path,
                         "sonic_line": {"x": x[win].tolist(), "y_d0": (ys[win] / delta0).tolist()}, "mean_reversal": rev}
        print(f"== {npz.stem}: t {float(d['t0']):.0f}-{float(d['t1']):.0f} ({int(d['n'])}), delta_0 {delta0:.3f}, cell {dy / delta0:.4f} delta_0, "
              f"L_int {l_int:.3f} = {l_int / delta0:.2f} delta_0")
        for k, v in up.items():
            print(f"  incoming layer at {k} (x {v['x']:.1f}): sonic height {v['y_sonic_d0']:.4f} delta_0 = {v['y_sonic_cells']:.2f} cells; "
                  f"M of cells 1-4 {', '.join(f'{m:.2f}' for m in v['M_cells'])}")
        print(f"  DNS profile (Zhang et al. M2p5, Re_tau 510, rescaled to Mach {mach:.2f}): sonic height {ys_ref:.4f} delta")
        for k, v in onset.items():
            print(f"  wall pressure reaches {k} of p_2 - p_inf at x {v['x']:.2f}: {v['upstream_d0']:.2f} delta_0 upstream of x_imp, x* {v['x_star']:+.2f} "
                  f"(DNS Re_theta 6882: {v['x_star_dns']:+.2f})")
        for k, pts in path.items():
            print(f"  compression path, p_inf + {k} (p_2 - p_inf): " + "; ".join(f"y {q['y_d0']:.2f}: x* {q['x_star']:+.2f}" for q in pts))
        print("  mean u < 0 of the resolved cells (not the wall model's tau_w): " + "; ".join(
            f"{k} x {v['x0']:.2f}-{v['x1']:.2f}, {v['extent_d0']:.2f} delta_0" if v else f"{k} none" for k, v in rev.items()))
    if json_out:
        Path(json_out).write_text(json.dumps(res))
    return 0


def _law(t: np.ndarray, tref: float, suth: float, omega: float) -> np.ndarray:
    """temperature_law of adam_flume_dissipation_library."""
    r = (t / tref) ** omega
    return r * (tref + suth) / (t + suth) if suth > 0.0 and omega > 0.0 else r


def _log_mean_inverse(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The mean of 1 / c over an interval where c is linear between a and b (log_mean_inverse)."""
    x = (b - a) / (b + a)
    small = np.abs(x) < 1.0e-3
    return np.where(small, 2.0 / (a + b) * (1.0 + x * x / 3.0), np.log(b / a) / np.where(small, 1.0, b - a))


WM_POINTS, WM_ITERATIONS, WM_STRETCH, WM_KAPPA, WM_APLUS, WM_PRT = 40, 16, 1.3, 0.41, 17.0, 0.9


def wall_model(h, us, ts, ps, tw, gas: dict, px=None, iterations: int = 0, ypg: float = 0.0, forcing=None, scale: str = "sum", px_scale=None,
               s_scale=None) -> tuple[np.ndarray, np.ndarray]:
    """tau_w and q_w of the equilibrium model, vectorised over the samples. `px` None: wall_model_solve operation for operation (16 plain
    fixed-point iterations; before #121 S4d 8 with an Aitken step on log tau every second one, which stopped a few samples far from the fixed point). `px` given: the momentum balance with the pressure gradient, `U(h) = tau I0 + px I1`, I0 = int dy / (mu + mu_t), I1 = int y dy /
    (mu + mu_t), solved for tau by a fixed point (tau can be negative; unrelaxed: 16 iterations, the Fortran count, are within 7e-5 of
    |tau| + |px| h over a sweep of the samples, 200 the default here), mu_t on the velocity scale sqrt((|tau| + |px| y) / rho) (with
    px = 0 the same model), the energy flux `cp k dT/dy + tau(y) U` constant across the layer.
    `ypg` > 0 (with `px`): the Bernoulli-consistent variant of Hickel, Touber, Bodart & Larsson (CTR Summer Program 2012, eq. 4.1): the
    convective term cancels the pressure gradient above y_pg = ypg (mu_w^2 / (rho_w |px|))^(1/3) and fades linearly below, so y is replaced
    by g(y) = y - y^2 / (2 y_pg) below y_pg and y_pg / 2 above, in the stress tau + px g(y), in I1 and in the velocity scale.
    `forcing` (with `px`, which it replaces): a callable of the model's heights y [samples, points] returning the integrated forcing
    S(y) = int_0^y (dp/dx + convection) dy' [samples, points]: the stress is tau + S(y), the convection-fed thin-boundary-layer balance.
    `scale`: the velocity scale of the eddy viscosity, "sum" sqrt((|tau| + |px g| or |S|) / rho) (the Fortran models), "tau" sqrt(|tau| / rho)
    (the forcing kept out of the eddy viscosity, entering the balance linearly), "duprat" sqrt(|tau| / rho + u_p^2) with u_p = (mu_w |dp/dx| /
    rho_w^2)^(1/3) of `px_scale` [samples] (the mixed scale of Duprat, Balarac, Metais, Congedo & Brugiere, Phys. Fluids 23, 015101, 2011:
    the scale only, in this model's mixing length and damping), "total" sqrt(|tau + S_bar(y)| / rho) with `s_scale` a callable like
    `forcing` (the total stress of a slow forcing S_bar: its fluctuation kept out of the eddy viscosity, its mean in), "pressure"
    sqrt((|tau| + |px_bar| y) / rho) with px_bar = `px_scale` (a slow wall pressure gradient: its stress in the eddy viscosity, the
    convection's not)."""
    r, cp, pr, mu0, tref, suth, om, sl = (gas[k] for k in ("R", "cp", "pr", "mu0", "tref", "suth", "omega", "semilocal"))
    h, us, ts, ps = (np.asarray(v, dtype=float)[:, None] for v in (h, us, ts, ps))
    dy0 = h * (WM_STRETCH - 1.0) / (WM_STRETCH**WM_POINTS - 1.0)
    y = np.concatenate([np.zeros_like(h), np.cumsum(dy0 * WM_STRETCH ** np.arange(WM_POINTS)[None, :], axis=1)], axis=1)
    t = tw + (ts - tw) * y / h
    rhow = ps / (r * tw)
    muw = mu0 * _law(np.asarray(tw), tref, suth, om)
    tau = muw * us / h
    hp = h * np.sqrt(tau / rhow) * rhow / muw
    for _ in range(3):
        hp = np.maximum(hp, 1.5)
        hp = (us / (np.log(hp) / WM_KAPPA + 5.0)) * h * rhow / muw
    tau = np.where(hp > 11.0, rhow * (hp * muw / (h * rhow)) ** 2, tau)
    dyy = np.diff(y, axis=1)
    p = 0.0 if px is None else np.asarray(px, dtype=float)[:, None]
    g = y
    if forcing is not None:
        p, g = 1.0, forcing(y)
    elif px is not None and ypg > 0.0:
        y_pg = ypg * np.cbrt(muw * muw / (rhow * np.maximum(np.abs(p), 1e-300)))
        g = np.where(y < y_pg, y - y * y / (2.0 * y_pg), 0.5 * y_pg)
    if scale == "duprat":
        up = np.cbrt(muw * np.abs(np.asarray(px_scale, dtype=float)[:, None]) / (rhow * rhow))
        tp = rhow * up * up
    sbar = s_scale(y) if scale == "total" else 0.0
    if scale == "pressure":
        pys = np.abs(np.asarray(px_scale, dtype=float)[:, None]) * y
    for _ in range((iterations or (WM_ITERATIONS if px is None and forcing is None else 200))):
        mu = mu0 * _law(t, tref, suth, om)
        rho = ps / (r * t)
        ts_ = np.abs(tau) + np.abs(p * g) if forcing is not None else np.abs(tau) + np.abs(p) * g
        if scale == "tau":
            ts_ = np.abs(tau) + 0.0 * g
        elif scale == "duprat":
            ts_ = np.abs(tau) + tp + 0.0 * g
        elif scale == "total":
            ts_ = np.abs(tau + sbar)
        elif scale == "pressure":
            ts_ = np.abs(tau) + pys
        damp = 1.0 - np.exp(-y * np.sqrt(rho * ts_) / (mu * WM_APLUS)) if sl else 1.0 - np.exp(-y * np.sqrt(ts_ / rhow) * rhow / (muw * WM_APLUS))
        mt = WM_KAPPA * y * np.sqrt(rho * ts_) * damp * damp
        a = mu / pr + mt / WM_PRT
        mt = mu + mt
        i0 = np.concatenate([np.zeros_like(h), np.cumsum(dyy * _log_mean_inverse(mt[:, :-1], mt[:, 1:]), axis=1)], axis=1)
        if px is None and forcing is None:
            tnew = us / i0[:, -1:]
            u = i0 * tnew
            stress = tnew
        else:
            yt = g / mt
            i1 = np.concatenate([np.zeros_like(h), np.cumsum(0.5 * dyy * (yt[:, :-1] + yt[:, 1:]), axis=1)], axis=1)
            tnew = (us - p * i1[:, -1:]) / i0[:, -1:]
            u = tnew * i0 + p * i1
            stress = tnew + p * g
        da = dyy * _log_mean_inverse(a[:, :-1], a[:, 1:])
        i2 = np.cumsum(da, axis=1)
        su = stress * u
        i3 = np.cumsum(0.5 * da * (su[:, :-1] + su[:, 1:]), axis=1)
        q = (cp * (ts - tw) + i3[:, -1:]) / i2[:, -1:]
        t[:, 1:] = np.maximum(0.5 * t[:, 1:] + 0.5 * (tw + (q * i2 - i3) / cp), 0.1 * np.minimum(tw, ts))
        tau = tnew
    return tau[:, 0], q[:, 0]


def _sample(col: np.ndarray, rep: np.ndarray, dy: float, hwm: float) -> tuple[np.ndarray, np.ndarray]:
    """The sample of `col` [x, rows of the finest grid] as sample_cells takes it on each column's own grid (`rep` finest rows per cell): the 3rd cell
    (hwm 0) or the height hwm between cells m and m + 1; return the values and the heights."""
    val, hh = np.empty(col.shape[0]), np.empty(col.shape[0])
    for i in range(col.shape[0]):
        cells, d = col[i, :: rep[i]], dy * rep[i]
        if hwm > 0.0:
            s = max(hwm / d + 0.5, 1.0)
            m, wgt = int(s), s - int(s)
            val[i] = cells[m - 1] + wgt * (cells[m] - cells[m - 1])
            hh[i] = (m + wgt - 0.5) * d
        else:
            val[i], hh[i] = cells[2], 2.5 * d
    return val, hh


def apriori(npz: Path, hwm: float, json_out: str | None, ypg: float = 4.0) -> int:
    """Apply the equilibrium wall model, with and without the mean pressure gradient, to the mean field of a run (see the module docstring)."""
    d = np.load(npz)
    x, dy, cf_log, delta0, x_imp = d["x"], float(d["dy"]), d["cf"], float(d["delta0"]), float(d["x_imp"])
    gas = {"R": 1.0, "cp": GAMMA / (GAMMA - 1.0), "pr": float(d["prandtl"]), "mu0": 1.0 / float(d["reynolds"]), "tref": float(d["t_inf"]),
           "suth": float(d["suth"]), "omega": 1.5 if float(d["suth"]) > 0.0 else 0.0, "semilocal": bool(d["semilocal"])}
    # a column of a coarser block repeats each value over the finest rows it covers (tbl_oracle.averages): two equal first rows
    rep = np.where(d["u"][:, 0] == d["u"][:, 1], 2, 1)
    us, h = _sample(d["u"], rep, dy, hwm)
    ss, _ = _sample(d["s"], rep, dy, hwm)
    ps, _ = _sample(d["p"], rep, dy, hwm)
    ts = ps / _sample(d["r"], rep, dy, hwm)[0]
    tw = float(d["t_wall"])
    # dp/dx of the mean wall pressure (first cell), central over one cell of each column's own grid
    pw = d["p"][:, 0]
    k = rep
    i = np.arange(len(x))
    lo, hi = np.clip(i - k, 0, len(x) - 1), np.clip(i + k, 0, len(x) - 1)
    px = (pw[hi] - pw[lo]) / (x[hi] - x[lo])
    cf0 = 2.0 * wall_model(h, us, ts, ps, tw, gas)[0]
    cf0s = 2.0 * wall_model(h, ss, ts, ps, tw, gas)[0]
    cfp = 2.0 * wall_model(h, us, ts, ps, tw, gas, px=px)[0]
    cfp_c = 2.0 * wall_model(h, us, ts, ps, tw, gas, px=px, iterations=400)[0]
    cfb = 2.0 * wall_model(h, us, ts, ps, tw, gas, px=px, ypg=ypg)[0]
    print(f"== {npz.name}: {int(d['n'])} checkpoints t {float(d['t0']):.1f}-{float(d['t1']):.1f}; delta_0 {delta0:.3f}; sample "
          f"{'3rd cell' if hwm == 0 else f'hwm {hwm:g}'}, h {h.min():.4f}-{h.max():.4f}")
    print(f"fixed point with dp/dx: max |C_f(200) - C_f(400)| {np.max(np.abs(cfp - cfp_c)):.2e}")
    zpg = (x >= x_imp - 25.0) & (x <= x_imp - 10.0)
    print(f"zero-pressure-gradient window x [{x_imp - 25:g}, {x_imp - 10:g}]: C_f history {cf_log[zpg].mean():.4e}, model on mean u {cf0[zpg].mean():.4e} "
          f"({100 * (cf0[zpg].mean() / cf_log[zpg].mean() - 1):+.1f} %), on mean |u_par| {cf0s[zpg].mean():.4e} "
          f"({100 * (cf0s[zpg].mean() / cf_log[zpg].mean() - 1):+.1f} %), with dp/dx {cfp[zpg].mean():.4e}")
    win = (x >= x_imp - 10.0 * delta0) & (x <= x_imp + 5.0 * delta0)
    for name, c in (("history", cf_log), ("model", cf0), ("model + dp/dx", cfp), (f"Bernoulli y*_pg {ypg:g}", cfb)):
        j = np.argmin(np.where(win, c, np.inf))
        neg = np.nonzero((c < 0.0) & win)[0]
        sep = (f"C_f < 0 from x {x[neg[0]]:.2f} to {x[neg[-1]]:.2f}, L_sep / delta_0 {(x[neg[-1]] - x[neg[0]]) / delta0:.2f}, {len(neg)} columns"
               if len(neg) else "no C_f < 0")
        print(f"  {name:20s}: C_f min {c[j]:.3e} at x {x[j]:.2f}; {sep}")
    print("     x    dp/dx    C_f hist   model    +dp/dx     Bernoulli")
    for xs in np.arange(x_imp - 8.0, x_imp + 6.5, 1.0):
        j = int(np.argmin(np.abs(x - xs)))
        print(f"  {x[j]:6.2f} {px[j]:+.2e} {cf_log[j]:+.3e} {cf0[j]:+.3e} {cfp[j]:+.3e} {cfb[j]:+.3e}")
    if json_out:
        Path(json_out).write_text(json.dumps({"x": x.tolist(), "cf_history": cf_log.tolist(), "cf_model": cf0.tolist(), "cf_model_speed": cf0s.tolist(),
                                              "cf_dpdx": cfp.tolist(), "cf_bernoulli": cfb.tolist(), "ypg": ypg, "dpdx": px.tolist(), "h": h.tolist(), "delta0": delta0, "x_imp": x_imp, "hwm": hwm}))
    return 0


def reversal(work: Path, x_imp: float, delta0: float, t0: float, json_out: str | None) -> int:
    """The reversed flow near the interaction per checkpoint (see the module docstring)."""
    import h5py

    tags = [tag for t, tag in to.checkpoints(work) if t >= t0]
    times = [t for t, _ in to.checkpoints(work) if t >= t0]
    if not tags:
        print(f"no checkpoint after t = {t0} in {work}")
        return 1
    x0, x1, y1 = x_imp - 4.5 * delta0, x_imp + delta0, 1.5 * delta0
    d_min = None
    for f5 in sorted(work.glob(f"tbl-{tags[0]}-proc*.h5")):
        with h5py.File(f5, "r") as h:
            for b in {k.rsplit("-", 1)[0] for k in h if k.endswith("-origin")}:
                d = h[b + "-dxdydz"][()]
                d_min = d if d_min is None else np.minimum(d_min, d)
    dz, dy, dx = (float(v) for v in d_min)
    lz = float(to.ini_value(work, "emax_z"))
    nz, ny, nx = round(lz / dz), round(y1 / dy), round((x1 - x0) / dx)
    rows = []
    for t, tag in zip(times, tags, strict=True):
        u3 = np.full((nz, ny, nx), np.nan)
        for f5 in sorted(work.glob(f"tbl-{tag}-proc*.h5")):
            with h5py.File(f5, "r") as h:
                for b in {k.rsplit("-", 1)[0] for k in h if k.endswith("-origin")}:
                    o, d = h[b + "-origin"][()], h[b + "-dxdydz"][()]
                    lo = o + to.NGC * d
                    shape = np.array(h[b + "-r"].shape) - 2 * to.NGC
                    hi = lo + shape * d
                    if hi[2] <= x0 or lo[2] >= x1 or lo[1] >= y1:
                        continue
                    s = (slice(to.NGC, -to.NGC),) * 3
                    u = h[b + "-ru"][()][s] / h[b + "-r"][()][s]
                    rz, ry, rx = (round(float(d[i]) / float(d_min[i])) for i in range(3))
                    u = np.repeat(np.repeat(np.repeat(u, rz, axis=0), ry, axis=1), rx, axis=2)
                    k0, j0, i0 = round(float(lo[0]) / dz), round(float(lo[1]) / dy), round((float(lo[2]) - x0) / dx)
                    kk, jj, ii = np.arange(k0, k0 + u.shape[0]), np.arange(j0, j0 + u.shape[1]), np.arange(i0, i0 + u.shape[2])
                    mk, mj, mi = (kk >= 0) & (kk < nz), (jj >= 0) & (jj < ny), (ii >= 0) & (ii < nx)
                    u3[np.ix_(kk[mk], jj[mj], ii[mi])] = u[np.ix_(mk, mj, mi)]
        rev = u3 < 0.0
        area = rev.sum() * dx * dy / nz / delta0**2
        yc = (np.arange(ny) + 0.5) * dy
        top = np.where(rev.any(axis=1), yc[ny - 1 - np.argmax(rev[:, ::-1, :], axis=1)], np.nan) / delta0
        cols = top[np.isfinite(top)]
        pct = np.percentile(cols, [50, 90, 99]) if cols.size else np.zeros(3)
        wall = rev[:, 0, :].mean()
        third = rev[:, 2, :].mean()
        um = np.nanmean(u3[:, 0, :], axis=0)
        xs = x0 + (np.arange(nx) + 0.5) * dx
        neg = xs[um < 0.0]
        ext = (neg.max() - neg.min() + dx) / delta0 if neg.size else 0.0
        rows.append([t, area, *pct, wall, ext, third])
        print(f"t {t:8.1f}: reversed area {area:6.3f} delta_0^2; column top p50 {pct[0]:5.3f} p90 {pct[1]:5.3f} p99 {pct[2]:5.3f} delta_0;"
              f" first cells reversed {100 * wall:5.1f} %, third {100 * third:5.1f} %; span-mean first-cell u < 0 over {ext:5.2f} delta_0", flush=True)
    a = np.array(rows)
    print(f"== {work.name}: {len(rows)} checkpoints t {a[0, 0]:.0f}-{a[-1, 0]:.0f}; reversed area mean {a[:, 1].mean():.3f}, max {a[:, 1].max():.3f}"
          f" delta_0^2; column top p99 mean {a[:, 4].mean():.3f}, max {a[:, 4].max():.3f} delta_0; first cells reversed mean {100 * a[:, 5].mean():.1f} %, third cells {100 * a[:, 7].mean():.1f} %;"
          f" span-mean first-cell u < 0 extent mean {a[:, 6].mean():.2f} delta_0")
    if json_out:
        Path(json_out).write_text(json.dumps({"columns": ["t", "area", "top_p50", "top_p90", "top_p99", "wall_fraction", "extent", "third_fraction"], "rows": rows,
                                              "delta0": delta0, "x_imp": x_imp, "window": [x0, x1, y1]}))
    return 0

def convective(npz: Path, json_out: str | None, gas_args: dict) -> int:
    """The wall model fed with the flow's own convection, a priori on the mean field of a run (an npz of `fields`): at each x, the sample of the
    3rd cell, the stress tau + S(y), S = int_0^y F, F = dp/dx + d(rho u u)/dx + d(rho u v)/dy of the mean field in the cells below the sample
    (linear in y between the cell centres, F(0) = dp/dx: no convection at the wall), the eddy viscosity on sqrt((|tau| + |S|) / rho): the
    steady thin-boundary-layer balance of Wang & Moin (2002) or Park & Moin (2014), with the convection taken from the resolved cells instead of
    from its own profile. Against: the equilibrium model, the pressure gradient alone (S4d's), the Bernoulli-consistent one (S4g's), and the mean
    u of the first cell. Three more velocity scales of the eddy viscosity: sqrt(|tau| / rho) (`wall_model_resolved_stress`: on a mean field the
    convection is the mean flow's, without the resolved stress the key removes), Duprat et al.'s (2011) mixed scale, the total stress of
    the forcing, slow here, being a mean, and |tau| + |dp/dx| y. A priori on time means: what the model would return for the mean field, not what a run with it converges to."""
    d = np.load(npz)
    x, y, delta0, x_imp = d["x"], d["y"], float(d["delta0"]), float(d["x_imp"])
    dx, dy = x[1] - x[0], y[1] - y[0]
    t_inf = float(d["t_inf"])
    gas = {"R": 1.0, "cp": GAMMA / (GAMMA - 1.0), "pr": gas_args["prandtl"], "mu0": 1.0 / gas_args["reynolds"], "tref": t_inf,
           "suth": gas_args["suth"], "omega": 1.5, "semilocal": gas_args["semilocal"]}
    r, ru, v, t, p = d["r"], d["ru"], d["v"], d["t"], d["p"]
    u = ru / r
    pw = p[:, 0]
    px = np.gradient(pw, dx)
    # the convection of the mean field in cells 1-3: d(rho u u)/dx (central) + d(rho u v)/dy (the wall's flux zero, central above)
    fx = np.gradient(ru[:, :4] * u[:, :4], dx, axis=0)
    ruv = np.concatenate([np.zeros((len(x), 1)), ru[:, :4] * v[:, :4]], axis=1)
    yy = np.concatenate([[0.0], y[:4]])
    fy = np.gradient(ruv, yy, axis=1)[:, 1:]
    conv = fx + fy  # [x, 4] at the centres of cells 1-4
    near = (x > x_imp - 10.0 * delta0) & (x < x_imp + 4.0 * delta0)
    idx = np.where(near)[0]
    h = np.full(idx.size, 2.5 * dy)
    us, ts, ps = u[idx, 2], t[idx, 2], p[idx, 2]
    tw = gas_args["t_wall"]
    yc = np.concatenate([[0.0], y[:3]])  # the wall and the centres of cells 1-3 (the sample at the 3rd)
    fc = np.concatenate([np.zeros((idx.size, 1)), conv[idx, :3]], axis=1) + px[idx, None]

    def forcing(ym: np.ndarray) -> np.ndarray:
        """S(y) = int_0^y F, F piecewise linear between the wall and the cell centres (constant beyond the last)."""
        out = np.zeros_like(ym)
        for n in range(ym.shape[0]):
            fine = np.linspace(0.0, float(ym[n, -1]), 2001)
            f = np.interp(fine, yc, fc[n])
            cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(fine))])
            out[n] = np.interp(ym[n], fine, cum)
        return out

    cf_eq = 2.0 * wall_model(h, us, ts, ps, tw, gas)[0]
    cf_pg = 2.0 * wall_model(h, us, ts, ps, tw, gas, px=px[idx])[0]
    cf_b = 2.0 * wall_model(h, us, ts, ps, tw, gas, px=px[idx], ypg=4.0)[0]
    cf_c = 2.0 * wall_model(h, us, ts, ps, tw, gas, forcing=forcing)[0]
    cf_ct = 2.0 * wall_model(h, us, ts, ps, tw, gas, forcing=forcing, scale="tau")[0]  # wall_model_resolved_stress's velocity scale
    cf_cd = 2.0 * wall_model(h, us, ts, ps, tw, gas, forcing=forcing, scale="duprat", px_scale=px[idx])[0]  # Duprat et al. (2011)
    cf_cs = 2.0 * wall_model(h, us, ts, ps, tw, gas, forcing=forcing, scale="total", s_scale=forcing)[0]  # the mean field's forcing is slow
    cf_cp = 2.0 * wall_model(h, us, ts, ps, tw, gas, forcing=forcing, scale="pressure", px_scale=px[idx])[0]  # the mean px is slow
    xs = x[idx]

    def sep(cf: np.ndarray) -> str:
        neg = xs[cf < 0.0]
        return f"x {neg.min():.2f}-{neg.max():.2f}, {(neg.max() - neg.min() + dx) / delta0:.2f} delta_0" if neg.size else "none"

    u1 = u[idx, 0]
    print(f"== {npz.stem}: x {xs[0]:.1f}-{xs[-1]:.1f}, sample the 3rd cell (h {2.5 * dy:.3f}), delta_0 {delta0:.3f}")
    print(f"  C_f < 0, equilibrium model: {sep(cf_eq)}; dp/dx alone: {sep(cf_pg)}; Bernoulli-consistent: {sep(cf_b)}; convection-fed: {sep(cf_c)}; and the velocity scale from tau_w: {sep(cf_ct)}")
    print(f"  mean u < 0 in the first cell: {sep(u1)}")
    print(f"  C_f minimum (1e-3): equilibrium {1e3 * cf_eq.min():.3f}, dp/dx {1e3 * cf_pg.min():.3f}, Bernoulli {1e3 * cf_b.min():.3f}, "
          f"convection-fed {1e3 * cf_c.min():.3f}, and the scale from tau_w {1e3 * cf_ct.min():.3f}")
    # the velocity-scale checks (#121 P4a): from x 54 on, the sign changes of C_f and the fastest first cell (against its incoming mean,
    # x 45-52) under a negative C_f; passing: at most 4 changes, and C_f < 0 only under a first cell below half its incoming speed
    u_in = float(u1[(xs >= 45.0) & (xs <= 52.0)].mean())
    down, win = xs >= 54.0, (xs >= 55.0) & (xs <= 64.0)
    for name, cf in (("dp/dx alone", cf_pg), ("convection-fed, |tau| + |S|", cf_c), ("scale from tau_w", cf_ct), ("Duprat u_tau_p", cf_cd),
                     ("total stress of the slow forcing", cf_cs), ("|tau| + |px_bar| y", cf_cp)):
        neg = cf[down] < 0.0
        changes = int(np.count_nonzero(np.diff(neg.astype(int))))
        umax = float((u1[down][neg] / u_in).max()) if neg.any() else 0.0
        verdict = "PASS" if changes <= 4 and umax < 0.5 else "FAIL"
        print(f"  {name:34s} x >= 54: {changes:3d} sign changes, fastest first cell under C_f < 0 {umax:.2f} of the incoming; C_f < 0 over "
              f"{100 * (cf[win] < 0.0).mean():3.0f} % of x 55-64, minimum {1e3 * cf.min():6.3f}: {verdict}")
    if json_out:
        Path(json_out).write_text(json.dumps({"x": xs.tolist(), "cf_eq": cf_eq.tolist(), "cf_pg": cf_pg.tolist(), "cf_bernoulli": cf_b.tolist(),
                                              "cf_convective": cf_c.tolist(), "cf_convective_tau": cf_ct.tolist(),
                                              "cf_convective_duprat": cf_cd.tolist(), "cf_convective_total": cf_cs.tolist(), "cf_convective_pressure": cf_cp.tolist(), "u1": u1.tolist(), "delta0": delta0, "x_imp": x_imp}))
    return 0

def main() -> int:
    """Dispatch the subcommand."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("work", type=Path)
    c.add_argument("--x-imp", type=float, required=True)
    c.add_argument("--deflection", type=float, required=True)
    c.add_argument("--t0", type=float, default=300.0)
    c.add_argument("--json")
    m = sub.add_parser("means")
    m.add_argument("work", type=Path)
    m.add_argument("--x-imp", type=float, required=True)
    m.add_argument("--t0", type=float, default=300.0)
    m.add_argument("--out", type=Path, required=True)
    f = sub.add_parser("fields")
    f.add_argument("work", type=Path)
    f.add_argument("--x-imp", type=float, required=True)
    f.add_argument("--deflection", type=float, required=True)
    f.add_argument("--t0", type=float, default=300.0)
    f.add_argument("--out", type=Path, required=True)
    q = sub.add_parser("compression")
    q.add_argument("npz", type=Path, nargs="+")
    q.add_argument("--json")
    p = sub.add_parser("apriori")
    p.add_argument("npz", type=Path)
    p.add_argument("--hwm", type=float, default=0.0)
    p.add_argument("--ypg", type=float, default=4.0, help="y*_pg of the Bernoulli-consistent model (Hickel et al. 2012)")
    p.add_argument("--json")
    v = sub.add_parser("reversal")
    v.add_argument("work", type=Path)
    v.add_argument("--x-imp", type=float, required=True)
    v.add_argument("--delta0", type=float, required=True)
    v.add_argument("--t0", type=float, default=300.0)
    v.add_argument("--json")
    k = sub.add_parser("convective")
    k.add_argument("npz", type=Path)
    k.add_argument("--reynolds", type=float, default=35500.0)
    k.add_argument("--prandtl", type=float, default=0.71)
    k.add_argument("--suth", type=float, default=0.10662621465945567, help="Sutherland temperature (the S4 runs')")
    k.add_argument("--t-wall", type=float, default=0.26354262850107724)
    k.add_argument("--classical", action="store_true", help="classical damping (default semi-local, the S4 runs')")
    k.add_argument("--json")
    a = ap.parse_args()
    if a.cmd == "convective":
        return convective(a.npz, a.json, {"reynolds": a.reynolds, "prandtl": a.prandtl, "suth": a.suth, "t_wall": a.t_wall, "semilocal": not a.classical})
    if a.cmd == "reversal":
        return reversal(a.work, a.x_imp, a.delta0, a.t0, a.json)
    if a.cmd == "means":
        return means(a.work, a.x_imp, a.t0, a.out)
    if a.cmd == "fields":
        return fields(a.work, a.x_imp, a.deflection, a.t0, a.out)
    if a.cmd == "compression":
        return compression(a.npz, a.json)
    if a.cmd == "apriori":
        return apriori(a.npz, a.hwm, a.json, a.ypg)
    return compare(a.work, a.x_imp, a.deflection, a.t0, a.json)


if __name__ == "__main__":
    sys.exit(main())
