#!/usr/bin/env python3
"""The shock/turbulent boundary layer interaction of issue #121 S4 (Mach 2.28, flow deflection 8 deg, the IUSTI experiment): its length scales and
wall pressure, from the checkpoints and the wall-model history of a run of make_tbl.py with --top-shock.

`compare <work> --x-imp X --deflection T [--t0 T0] [--json out]`: the checkpoints from T0 on, averaged over the span and in time (tbl_oracle.averages),
and the wall history over the same window:
- the incoming layer: delta(x), delta*, theta and H against the local edge state (tbl_oracle.edge_shape) and the nominal free stream; delta_0, the
  99 % thickness the undisturbed layer would have at x_imp, extrapolated from a linear fit of delta over [x_imp - 15, x_imp - 6] (upstream of the
  interaction's influence); delta* the same way (the Souverein et al. 2013 scaling takes the upstream displacement thickness); Re_theta of both;
- C_f(x) from the wall model's tau_w: the separation (the first C_f < 0 downstream of x_imp - 10 delta_0) and the reattachment (the last), L_sep;
- the mean wall pressure (the first cell), as (p_w - p_inf) / (p_2 - p_inf) with p_2 the pressure behind the incident shock (the normalisation
  of the STREAmS figure), against x* = (x - x_imp) / L;
- the interaction length L_int = x_imp - x_ref, x_ref the apparent origin of the reflected shock: above the layer every row crosses the incident
  shock (p from p_1 to p_2) and then the reflected one (p_2 to p_3, the exact regular reflection), the mid-pressure crossings of the reflected one
  over the rows y in [2.5, 5] delta_0 fitted by a line and extrapolated to the wall; the incident one fitted the same way, its foot against x_imp;
- the Souverein scaling L* = L_int / (delta* g3), g3 = sin(beta - phi) / (sin beta sin phi), and the separation criterion Se = k1 (p_3 - p_1) / q_inf
  with k1 = 3 (low Reynolds number, Volpiani et al. 2018).
Recorded, not judged: the comparison with the digitised references (reference/sbli_m2p28_*) is printed next to each value.

`means <work> --x-imp X [--t0 T0] --out <npz>`: the span and time averages near the wall (tbl_oracle.averages, the first 12 rows of the finest grid),
the C_f of the wall-model history, delta_0 and the run's gas and viscosity law, in one file for `apriori` (run it where the checkpoints are).
`apriori <npz> [--hwm H] [--json out]`: the equilibrium wall model of the run (a port of `wall_model_solve`, adam_flume_wall_model_library.F90) applied
to the mean field at the run's sample (the 3rd cell of the local grid, or the height H between cells as `sample_cells` takes it) per x, with and without
the mean streamwise pressure gradient in the momentum balance, `(mu + mu_t) dU/dy = tau_w + (dp/dx) y`, the eddy viscosity on the velocity scale
`sqrt((|tau_w| + |dp/dx| y) / rho)` (S4d, issue #121): C_f against the history's (the history averages tau of the instantaneous samples, the model
here takes the mean sample), and the C_f < 0 region the pressure gradient opens, if any. A priori: what a run with the term would start from, not
what it would converge to.
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
    # the reflected and the incident shock, fitted above the layer
    pm_i, pm_r = 0.5 * (p_inf + p2), 0.5 * (p2 + p3)
    tb = math.tan(math.radians(sh["incident"]["beta"]))
    rows = [j for j in range(ny) if 2.5 * delta0 <= y[j] <= 5.0 * delta0]
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
    x_inc_wall, x_ref = float(ci[1]), float(cr[1])
    l_int = x_imp - x_ref
    beta_r = math.degrees(math.atan(1.0 / cr[0])) if len(yr) > 3 else float("nan")
    b, ph = math.radians(sh["incident"]["beta"]), math.radians(theta_deg)
    g3 = math.sin(b - ph) / (math.sin(b) * math.sin(ph))
    l_star = l_int / (dstar0 * g3)
    se = 3.0 * (p3 - p_inf) / (0.5 * 1.0 * 1.0)  # k1 = 3, q_inf = rho_inf U_inf^2 / 2 = 1 / 2
    p_w = p[:, 0]
    pw_norm = (p_w - p_inf) / (p2 - p_inf)
    x_star = (x - x_imp) / l_int
    ref = HERE / "reference"
    dns = np.loadtxt(ref / "sbli_m2p28_fig7a_dns.txt")
    iusti = np.loadtxt(ref / "sbli_m2p28_fig7a_iusti.txt")
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
    print(f"separation: C_f < 0 from x = {x_sep:.2f} to {x_att:.2f}: L_sep / delta_0 {l_sep / delta0:.2f} (Volpiani et al. adiabatic {theta_deg:g} deg, Re_theta 1050: "
          f"{vol.get('a', float('nan')):.2f})")
    print(f"incident shock: foot at x = {x_inc_wall:.2f} (x_imp {x_imp:g}); reflected shock: {len(yr)} rows, {beta_r:.2f} deg to the wall, apparent origin x_ref {x_ref:.2f}")
    print(f"L_int = x_imp - x_ref {l_int:.3f}: L_int / delta_0 {l_int / delta0:.2f} (Volpiani et al. {vol.get('b', float('nan')):.2f}); Souverein L* = L_int / (delta* g3) "
          f"{l_star:.3f}, g3 {g3:.3f}, Se {se:.3f}")
    for xs in (-1.0, -0.5, 0.0, 0.5, 1.0, 2.0):
        print(f"  x* {xs:+.1f}: (p_w - p_inf) / (p_2 - p_inf) {np.interp(xs, x_star, pw_norm):.3f} (DNS {np.interp(xs, dns[:, 0], dns[:, 1]):.3f})")
    if json_out:
        out = {"x": x.tolist(), "cf": cf.tolist(), "pw_norm": pw_norm.tolist(), "x_star": x_star.tolist(), "delta": [None if not np.isfinite(v) else float(v) for v in delta],
               "H": [None if not np.isfinite(v) else float(v) for v in h_edge], "p_e": [None if not np.isfinite(v) else float(v) for v in p_e],
               "delta0": delta0, "dstar0": dstar0, "dstar0_nominal": dstar0_nom, "theta0": theta0, "re_theta": re_theta_nom, "re_theta_edge": re_theta_edge,
               "x_sep": x_sep, "x_att": x_att, "l_sep": l_sep, "x_ref": x_ref, "l_int": l_int, "x_inc_wall": x_inc_wall, "beta_r": beta_r, "l_star": l_star,
               "g3": g3, "se": se, "x_imp": x_imp, "t": [t for t, _ in cps], "volpiani": vol, "iusti": iusti.tolist(), "dns": dns.tolist(),
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


def wall_model(h, us, ts, ps, tw, gas: dict, px=None, iterations: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """tau_w and q_w of the equilibrium model, vectorised over the samples. `px` None: wall_model_solve operation for operation (16 plain
    fixed-point iterations; before #121 S4d 8 with an Aitken step on log tau every second one, which stopped a few samples far from the fixed point). `px` given: the momentum balance with the pressure gradient, `U(h) = tau I0 + px I1`, I0 = int dy / (mu + mu_t), I1 = int y dy /
    (mu + mu_t), solved for tau by a fixed point (tau can be negative; unrelaxed: 16 iterations, the Fortran count, are within 7e-5 of
    |tau| + |px| h over a sweep of the samples, 200 the default here), mu_t on the velocity scale sqrt((|tau| + |px| y) / rho) (with
    px = 0 the same model), the energy flux `cp k dT/dy + tau(y) U` constant across the layer."""
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
    for _ in range((iterations or (WM_ITERATIONS if px is None else 200))):
        mu = mu0 * _law(t, tref, suth, om)
        rho = ps / (r * t)
        ts_ = np.abs(tau) + np.abs(p) * y
        damp = 1.0 - np.exp(-y * np.sqrt(rho * ts_) / (mu * WM_APLUS)) if sl else 1.0 - np.exp(-y * np.sqrt(ts_ / rhow) * rhow / (muw * WM_APLUS))
        mt = WM_KAPPA * y * np.sqrt(rho * ts_) * damp * damp
        a = mu / pr + mt / WM_PRT
        mt = mu + mt
        i0 = np.concatenate([np.zeros_like(h), np.cumsum(dyy * _log_mean_inverse(mt[:, :-1], mt[:, 1:]), axis=1)], axis=1)
        if px is None:
            tnew = us / i0[:, -1:]
            u = i0 * tnew
            stress = tnew
        else:
            yt = y / mt
            i1 = np.concatenate([np.zeros_like(h), np.cumsum(0.5 * dyy * (yt[:, :-1] + yt[:, 1:]), axis=1)], axis=1)
            tnew = (us - p * i1[:, -1:]) / i0[:, -1:]
            u = tnew * i0 + p * i1
            stress = tnew + p * y
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


def apriori(npz: Path, hwm: float, json_out: str | None) -> int:
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
    print(f"== {npz.name}: {int(d['n'])} checkpoints t {float(d['t0']):.1f}-{float(d['t1']):.1f}; delta_0 {delta0:.3f}; sample "
          f"{'3rd cell' if hwm == 0 else f'hwm {hwm:g}'}, h {h.min():.4f}-{h.max():.4f}")
    print(f"fixed point with dp/dx: max |C_f(200) - C_f(400)| {np.max(np.abs(cfp - cfp_c)):.2e}")
    zpg = (x >= x_imp - 25.0) & (x <= x_imp - 10.0)
    print(f"zero-pressure-gradient window x [{x_imp - 25:g}, {x_imp - 10:g}]: C_f history {cf_log[zpg].mean():.4e}, model on mean u {cf0[zpg].mean():.4e} "
          f"({100 * (cf0[zpg].mean() / cf_log[zpg].mean() - 1):+.1f} %), on mean |u_par| {cf0s[zpg].mean():.4e} "
          f"({100 * (cf0s[zpg].mean() / cf_log[zpg].mean() - 1):+.1f} %), with dp/dx {cfp[zpg].mean():.4e}")
    win = (x >= x_imp - 10.0 * delta0) & (x <= x_imp + 5.0 * delta0)
    for name, c in (("history", cf_log), ("model", cf0), ("model + dp/dx", cfp)):
        j = np.argmin(np.where(win, c, np.inf))
        neg = np.nonzero((c < 0.0) & win)[0]
        sep = (f"C_f < 0 from x {x[neg[0]]:.2f} to {x[neg[-1]]:.2f}, L_sep / delta_0 {(x[neg[-1]] - x[neg[0]]) / delta0:.2f}, {len(neg)} columns"
               if len(neg) else "no C_f < 0")
        print(f"  {name:14s}: C_f min {c[j]:.3e} at x {x[j]:.2f}; {sep}")
    print("     x    dp/dx    C_f hist   model    +dp/dx")
    for xs in np.arange(x_imp - 8.0, x_imp + 6.5, 1.0):
        j = int(np.argmin(np.abs(x - xs)))
        print(f"  {x[j]:6.2f} {px[j]:+.2e} {cf_log[j]:+.3e} {cf0[j]:+.3e} {cfp[j]:+.3e}")
    if json_out:
        Path(json_out).write_text(json.dumps({"x": x.tolist(), "cf_history": cf_log.tolist(), "cf_model": cf0.tolist(), "cf_model_speed": cf0s.tolist(),
                                              "cf_dpdx": cfp.tolist(), "dpdx": px.tolist(), "h": h.tolist(), "delta0": delta0, "x_imp": x_imp, "hwm": hwm}))
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
    p = sub.add_parser("apriori")
    p.add_argument("npz", type=Path)
    p.add_argument("--hwm", type=float, default=0.0)
    p.add_argument("--json")
    a = ap.parse_args()
    if a.cmd == "means":
        return means(a.work, a.x_imp, a.t0, a.out)
    if a.cmd == "apriori":
        return apriori(a.npz, a.hwm, a.json)
    return compare(a.work, a.x_imp, a.deflection, a.t0, a.json)


if __name__ == "__main__":
    sys.exit(main())
