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
    a = ap.parse_args()
    return compare(a.work, a.x_imp, a.deflection, a.t0, a.json)


if __name__ == "__main__":
    sys.exit(main())
