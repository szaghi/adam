#!/usr/bin/env python3
"""The compressible turbulent boundary layer of issue #119 T3.7: reference data, the digital-filter inflow checks, the comparison with the DNS.

`profile <Stat.dat> <out>`: convert a station file of the Zhang, Duan & Choudhari (2018) DNS (NASA Turbulence Modeling Resource, `<case>_Stat.dat`, Tecplot
POINT, z wall-normal, y spanwise) into the profile table of the digital-filter inflow and the `boundary-layer` initial condition: `y / delta`, `U / U_inf`,
`T / T_inf`, `R_uu`, `R_vv`, `R_ww`, `R_uv` in units of `U_inf^2` (v wall-normal, w spanwise), delta the station's 99 % thickness of the file header.
`agree <a> <b>`: two runs on the same step sequence (ranks, or CPU against GPU): the conservation history within 1e-11 of its column scale.
`compare <work> <Stat file> [--t0 T] [--json out]`: the checkpoints from time T on, averaged over the span and in time: per x the 99 % thickness, the momentum
and displacement thicknesses (theta from the Favre mean rho u), H, Re_theta and C_f = 2 d theta / dx (zero pressure gradient; a quadratic fit of theta over
+-3 delta_i); at the station where Re_theta is the DNS's (or the nearest), C_f against the DNS's 2 (u_tau / U_inf)^2 rho_w / rho_inf, H against its H, the
van Driest velocity u+_VD against y+ (with the run's own u_tau) and the rms of u, v (wall-normal), w (spanwise) and <u'v'> over u_tau against y / delta.
Recorded, not judged: the result is the result.
`wall-force <channel work>`: the wall-model history of a channel at a steady state: the mean tau_w of both walls (all bins, the last quarter of the rows)
equals the driving force of the forcing history (per unit wall area, the half height 1) within 1e-3.
`wall-agree <a> <b> [--tol T]`: two runs' wall-model histories (every face file) within T (default 1e-12) of each column's scale.
`inflow <work>`: the ghost layer of the inflow face in the last checkpoint (saved with ghosts): the spanwise mean of u at y = delta / 2 within 10 % of the table,
and fluctuations present (the rms of the spanwise and wall-normal velocity above 1e-3 U_inf below delta).

Checkpoint blocks are (z, y, x) arrays with the ghosts, their origin (the ghost corner) and spacing ordered (z, y, x) too; gamma 1.4."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

NGC, GAMMA = 3, 1.4


def read_stat(path: Path) -> tuple[dict[str, float], np.ndarray]:
    """Return the header scalars (delta_mm, utau, Uinf) and the data rows of a station file."""
    text = path.read_text().splitlines()
    head: dict[str, float] = {}
    for i, line in enumerate(text):
        if line.startswith("# Minf"):
            v = text[i + 1].split()[1:]
            head["Minf"], head["Uinf"], head["Tinf"], head["Tw"] = float(v[0]), float(v[1]), float(v[3]), float(v[4])
        if line.startswith("# x_a/delta_i"):
            v = text[i + 1].split()[1:]
            head["delta_mm"], head["utau"] = float(v[7]), float(v[9])
    start = next(i for i, line in enumerate(text) if "ZONETYPE" in line)
    rows = [list(map(float, line.split())) for line in text[start + 1 :] if line.strip() and not line.lstrip().startswith(("#", "D"))]
    return head, np.array(rows)


def profile(src: str, out: str) -> int:
    """Write the profile table."""
    head, d = read_stat(Path(src))
    ut = head["utau"] / head["Uinf"]
    zod, u, t = d[:, 1], d[:, 4], d[:, 6]
    uu, vv_span, ww_norm, uw = (d[:, 8] * ut) ** 2, (d[:, 9] * ut) ** 2, (d[:, 10] * ut) ** 2, d[:, 15] * ut**2
    lines = [
        f"# From {Path(src).name} (Zhang, Duan & Choudhari 2018, AIAA J 56(11), doi:10.2514/1.J057296; NASA Turbulence Modeling Resource):",
        f"# M_inf {head['Minf']}, T_w / T_inf {head['Tw'] / head['Tinf']:.6f}, u_tau / U_inf {ut:.6f}; y / delta with delta = {head['delta_mm']} mm.",
        "# y/delta  U/U_inf  T/T_inf  R_uu  R_vv(wall-normal)  R_ww(spanwise)  R_uv  [U_inf^2]",
    ]
    lines += [f"{a:.8e} {b:.8e} {c:.8e} {e:.8e} {f:.8e} {g:.8e} {h:.8e}" for a, b, c, e, f, g, h in zip(zod, u, t, uu, ww_norm, vv_span, uw)]
    Path(out).write_text("\n".join(lines) + "\n")
    print(f"profile: {len(zod)} rows to {out}, T_w / T_inf {head['Tw'] / head['Tinf']:.4f}, u_tau / U_inf {ut:.5f}")
    return 0


def history(path: Path) -> np.ndarray:
    """Load a history; restarts repeat rows, keep the last of each step."""
    f = np.loadtxt(path, ndmin=2, comments="VARIABLES")
    _, last = np.unique(f[::-1, 0], return_index=True)
    return f[::-1][last]


def agree(a: str, b: str) -> int:
    """The conservation histories of two runs within 1e-11 of their column scale (floor 1: the transverse momenta are about 0)."""
    x, y = history(Path(a) / "tbl-conservation_history.dat"), history(Path(b) / "tbl-conservation_history.dat")
    if x.shape != y.shape:
        print(f"conservation history: shapes {x.shape} and {y.shape}  FAIL")
        return 1
    d = float(np.max(np.abs(x - y)[:, 2:] / np.maximum(np.abs(x[:, 2:]).max(axis=0), 1.0)))
    ok = d < 1e-11
    print(f"tbl-conservation_history.dat: {d:.2e} (bound 1e-11)  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def ini_value(work: Path, key: str) -> str:
    """Return the first `key = value` of the run's input.ini."""
    for line in (work / "input.ini").read_text().splitlines():
        if line.split("=")[0].strip() == key:
            return line.split("=", 1)[1].strip()
    raise KeyError(key)


def inflow(work_s: str) -> int:
    """The ghost layer of the inflow face (x_min) in the last checkpoint."""
    import h5py

    work = Path(work_s)
    tags = sorted({p.name.split("-")[1] for p in work.glob("tbl-*-proc*.h5") if p.name.split("-")[1].isdigit()})
    tag = tags[-1]
    delta = float(ini_value(work, "df_delta"))
    table = np.loadtxt(work / ini_value(work, "df_profile"), comments="#")
    u_ref = float(np.interp(0.5, table[:, 0], table[:, 1]))
    ys, us, vs, ws = [], [], [], []
    for f5 in sorted(work.glob(f"tbl-{tag}-proc*.h5")):
        with h5py.File(f5, "r") as h:
            for b in {k.rsplit("-", 1)[0] for k in h}:
                o, dd = h[b + "-origin"][()], h[b + "-dxdydz"][()]
                if abs(o[2] + NGC * dd[2]) > 1e-12:  # the block's interior must start at x = 0 (origin and spacing are (z, y, x))
                    continue
                r = h[f"{b}-r"][()][NGC:-NGC, NGC:-NGC, NGC - 1]
                ru, rv, rw = (h[f"{b}-{v}"][()][NGC:-NGC, NGC:-NGC, NGC - 1] for v in ("ru", "rv", "rw"))
                y = o[1] + (np.arange(r.shape[1]) + NGC + 0.5) * dd[1]
                ys.append(np.broadcast_to(y, r.shape).ravel())
                us.append((ru / r).ravel())
                vs.append((rv / r).ravel())
                ws.append((rw / r).ravel())
    if not ys:
        print("inflow: no block at the inflow face  FAIL")
        return 1
    y, u, v, w = (np.concatenate(z) for z in (ys, us, vs, ws))
    near = np.abs(y / delta - 0.5) < 0.1
    below = y < delta
    um = float(u[near].mean())
    vr, wr = float(v[below].std()), float(w[below].std())
    ok_u = abs(um - u_ref) / u_ref < 0.10
    ok_f = vr > 1e-3 and wr > 1e-3
    print(f"inflow ghosts at checkpoint {tag}: <u> at y = delta/2 {um:.4f} against the table {u_ref:.4f} (bound 10 %)  {'PASS' if ok_u else 'FAIL'}")
    print(f"inflow ghosts below delta: rms v {vr:.2e}, rms w {wr:.2e} (min 1e-3)  {'PASS' if ok_f else 'FAIL'}")
    return 0 if ok_u and ok_f else 1


def checkpoints(work: Path) -> list[tuple[float, str]]:
    """Return (time, iteration tag) of every checkpoint of the run in time order (the restart files are not ones)."""
    out = []
    for x in work.glob("tbl-*.xdmf"):
        m = re.search(r'<Time Value="([^"]+)"', x.read_text()[:4000])
        if m and x.stem.split("-")[1].isdigit():
            out.append((float(m.group(1)), x.stem.split("-")[1]))
    return sorted(out)


def averages(work: Path, tags: list[str]) -> dict[str, np.ndarray]:
    """Span and time averages on the (x, y) cell grid: r, ru, u, v, w, T and the second moments uu, vv, ww, uv (Reynolds)."""
    import h5py

    acc: dict[str, np.ndarray] = {}
    count = None
    for tag in tags:
        for f5 in sorted(work.glob(f"tbl-{tag}-proc*.h5")):
            with h5py.File(f5, "r") as h:
                for b in {k.rsplit("-", 1)[0] for k in h}:
                    o, d = h[b + "-origin"][()], h[b + "-dxdydz"][()]
                    s = (slice(NGC, -NGC),) * 3
                    r, ru, rv, rw, re_ = (h[f"{b}-{v}"][()][s] for v in ("r", "ru", "rv", "rw", "rE"))
                    u, v, w = ru / r, rv / r, rw / r
                    t = (GAMMA - 1.0) * (re_ - 0.5 * r * (u**2 + v**2 + w**2)) / r
                    if not acc:
                        nx = int(round(float(ini_value(work, "emax_x")) / d[2]))  # origin and spacing are (z, y, x)
                        ny = int(round(float(ini_value(work, "emax_y")) / d[1]))
                        for k in ("r", "ru", "u", "v", "w", "t", "uu", "vv", "ww", "uv"):
                            acc[k] = np.zeros((nx, ny))
                        count = np.zeros((nx, ny))
                    i0 = int(round((o[2] + NGC * d[2]) / d[2]))
                    j0 = int(round((o[1] + NGC * d[1]) / d[1]))
                    sl = (slice(i0, i0 + r.shape[2]), slice(j0, j0 + r.shape[1]))
                    for k, a in (("r", r), ("ru", ru), ("u", u), ("v", v), ("w", w), ("t", t), ("uu", u * u), ("vv", v * v), ("ww", w * w), ("uv", u * v)):
                        acc[k][sl] += a.sum(axis=0).T
                    count[sl] += r.shape[0]
    return {k: a / count for k, a in acc.items()}


def compare(work_s: str, stat: str, t0: float, json_out: str | None, re_inf: float | None = None) -> int:
    """Print the comparison with the DNS (and write it as JSON)."""
    import json

    work = Path(work_s)
    cps = [(t, g) for t, g in checkpoints(work) if t >= t0]
    if not cps:
        print(f"compare: no checkpoint from t = {t0}")
        return 1
    a = averages(work, [g for _, g in cps])
    reynolds = float(ini_value(work, "reynolds"))
    law = ini_value(work, "viscosity_law")
    re_theta_scale = re_inf if re_inf else reynolds  # Re_theta with the true free-stream viscosity when the power law is matched at the wall
    t_inf, t_w = float(ini_value(work, "reference_temperature")), float(ini_value(work, "wall_temperature"))
    lx, ly = float(ini_value(work, "emax_x")), float(ini_value(work, "emax_y"))
    nx, ny = a["u"].shape
    dx, dy = lx / nx, ly / ny
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    uf = a["ru"] / a["r"]
    edge = np.argmax(a["u"] >= 0.99, axis=1)
    delta = np.where(edge > 0, y[np.maximum(edge, 1) - 1] + dy * (0.99 - a["u"][np.arange(nx), np.maximum(edge, 1) - 1])
                     / (a["u"][np.arange(nx), np.maximum(edge, 1)] - a["u"][np.arange(nx), np.maximum(edge, 1) - 1]), np.nan)
    theta = (a["r"] * uf * (1.0 - uf)).sum(axis=1) * dy
    dstar = (1.0 - a["r"] * uf).sum(axis=1) * dy
    shape = dstar / theta
    re_theta = re_theta_scale * theta
    cf = np.full(nx, np.nan)
    half = int(round(3.0 / dx))
    for i in range(half, nx - half):
        c = np.polyfit(x[i - half : i + half + 1] - x[i], theta[i - half : i + half + 1], 2)
        cf[i] = 2.0 * c[1]
    head, d = read_stat(Path(stat))
    re_dns = {"M2p5": 2835.0, "M6Tw076": 9455.0}.get(Path(stat).name.split("_")[0], float("nan"))
    ut_dns = head["utau"] / head["Uinf"]
    cf_dns = 2.0 * ut_dns**2 * head["Tinf"] / head["Tw"]
    valid = np.isfinite(cf)
    ist = int(np.argmin(np.where(valid, np.abs(re_theta - re_dns), np.inf)))
    rho_w = 1.0 / (t_w / t_inf)  # p_inf / (R T_w) over rho_inf = T_inf / T_w
    cf_mom = cf.copy()
    wall = work / "tbl-wall_y_min_history.dat"
    cf_wall = None
    if wall.exists():  # the wall model's own tau_w (rho_inf U_inf^2 = 1), its bins mapped onto the cells along x, averaged over the same window
        w = wall_rows(wall)
        w = w[w[:, 1] >= t0]
        nb = (w.shape[1] - 2) // 2
        xb = (np.arange(nb) + 0.5) * lx / nb
        cf_wall = 2.0 * np.interp(x, xb, w[:, 2 : 2 + nb].mean(axis=0))
        q_wall = np.interp(x, xb, w[:, 2 + nb : 2 + 2 * nb].mean(axis=0))
        cf = cf_wall
    ut = np.sqrt(0.5 * cf[ist] / rho_w)
    if law == "sutherland":  # (T / T_ref)^(3/2) (T_ref + S) / (T + S), T_ref the free stream
        s_t = float(ini_value(work, "sutherland_temperature"))
        mu_w = (t_w / t_inf) ** 1.5 * (t_inf + s_t) / (t_w + s_t) / reynolds
    else:
        mu_w = (t_w / t_inf) ** float(ini_value(work, "viscosity_exponent")) / reynolds
    yplus = y * ut * rho_w / mu_w
    r_st, u_st = a["r"][ist], a["u"][ist]
    uvd = np.concatenate([[0.0], np.cumsum(np.sqrt(0.5 * (r_st[1:] + r_st[:-1]) / rho_w) * np.diff(u_st))]) + np.sqrt(r_st[0] / rho_w) * u_st[0]
    uvd /= ut
    rms = {k: np.sqrt(np.maximum(a[k + k][ist] - a[k][ist] ** 2, 0.0)) / ut for k in ("u", "v", "w")}
    uv = (a["uv"][ist] - a["u"][ist] * a["v"][ist]) / ut**2
    d_dns = delta[ist]
    print(f"== {work.name}: {len(cps)} checkpoints from t = {cps[0][0]:.1f} to {cps[-1][0]:.1f}, {nx} x {ny} cells")
    print(f"station x = {x[ist]:.2f}: Re_theta {re_theta[ist]:.0f} (DNS {re_dns:.0f}), delta {d_dns:.3f}, theta {theta[ist]:.4f}, H {shape[ist]:.2f} (DNS {float(text_h(stat)):.2f})")
    source = "the wall model's tau_w" if cf_wall is not None else "the momentum integral"
    print(f"C_f ({source}) {cf[ist]:.4e} against the DNS {cf_dns:.4e}: {100.0 * (cf[ist] / cf_dns - 1.0):+.1f} %;  u_tau / U_inf {ut:.5f} (DNS {ut_dns:.5f})")
    b_q = b_q_dns = None
    if cf_wall is not None:
        print(f"C_f (the momentum integral) {cf_mom[ist]:.4e}: {100.0 * (cf_mom[ist] / cf_dns - 1.0):+.1f} %")
        cp = GAMMA / (GAMMA - 1.0)  # R = 1
        b_q = -q_wall[ist] / (rho_w * cp * ut * t_w)  # the DNS's sign: B_q = q_w / (rho_w cp u_tau T_w), q_w into the gas, so negative at a cold wall
        b_q_dns = -float(text_bq(stat))
        print(f"B_q {b_q:+.4f} against the DNS {b_q_dns:+.4f}; wall heat flux {q_wall[ist]:.4e} (into the wall, code units)")
    print(" x      Re_theta  delta   H      C_f       C_f momentum")
    for i in range(half, nx - half, max(1, (nx - 2 * half) // 12)):
        print(f" {x[i]:5.1f}  {re_theta[i]:7.0f}  {delta[i]:6.3f}  {shape[i]:5.2f}  {cf[i]:.4e}  {cf_mom[i]:.4e}")
    zp_dns, uvd_dns = d[:, 2], d[:, 24]
    for yp in (30.0, 100.0, 200.0):
        if yplus[0] <= yp <= yplus[-1]:
            print(f"u+_VD at y+ = {yp:.0f}: {np.interp(yp, yplus, uvd):.2f} (DNS {np.interp(yp, zp_dns, uvd_dns):.2f})")
    for name, col, val in (("u_rms", 8, rms["u"]), ("v_rms (normal)", 10, rms["v"]), ("w_rms (span)", 9, rms["w"])):
        print(f"{name} / u_tau at y / delta = 0.2, 0.5: {np.interp(0.2, y / d_dns, val):.2f}, {np.interp(0.5, y / d_dns, val):.2f} "
              f"(DNS {np.interp(0.2, d[:, 1], d[:, col]):.2f}, {np.interp(0.5, d[:, 1], d[:, col]):.2f})")
    print(f"<u'v'> / u_tau^2 at y / delta = 0.2, 0.5: {np.interp(0.2, y / d_dns, uv):.2f}, {np.interp(0.5, y / d_dns, uv):.2f} "
          f"(DNS {np.interp(0.2, d[:, 1], d[:, 15]):.2f}, {np.interp(0.5, d[:, 1], d[:, 15]):.2f})")
    if json_out:
        out = {"x": x.tolist(), "re_theta": re_theta.tolist(), "cf": [None if not np.isfinite(c) else float(c) for c in cf],
               "cf_mom": [None if not np.isfinite(c) else float(c) for c in cf_mom], "cf_source": source, "H": shape.tolist(),
               "delta": [None if not np.isfinite(c) else float(c) for c in delta], "station": int(ist), "cf_dns": cf_dns, "re_dns": re_dns,
               "y": y.tolist(), "yplus": yplus.tolist(), "uvd": uvd.tolist(), "rms": {k: v.tolist() for k, v in rms.items()}, "uv": uv.tolist(),
               "delta_station": float(d_dns), "t": [t for t, _ in cps], "b_q": b_q, "b_q_dns": b_q_dns, "h_dns": float(text_h(stat))}
        Path(json_out).write_text(json.dumps(out))
    return 0


def text_bq(stat: str) -> float:
    """The -B_q of a station file's header."""
    text = Path(stat).read_text().splitlines()
    i = next(i for i, line in enumerate(text) if line.startswith("# x_a/delta_i"))
    return float(text[i + 1].split()[11])


def text_h(stat: str) -> float:
    """The shape factor H of a station file's header."""
    text = Path(stat).read_text().splitlines()
    i = next(i for i, line in enumerate(text) if line.startswith("# x_a/delta_i"))
    return float(text[i + 1].split()[7])


def wall_rows(path: Path) -> np.ndarray:
    """Load a wall-model history: it, time, the bins' tau, then their q, one row per step. A run resumed from a restart appends again the steps
    after that restart (bitwise the same rows): the last of each step is kept, so a time mean does not weigh those stretches twice."""
    w = np.loadtxt(path, comments=("#", "VARIABLES"), ndmin=2)
    _, last = np.unique(w[::-1, 0], return_index=True)
    return w[np.sort(len(w) - 1 - last)]


def wall_force(work_s: str) -> int:
    """tau_w of both walls against the channel's driving force at a steady state."""
    work = Path(work_s)
    taus = []
    for face in ("y_min", "y_max"):
        w = wall_rows(work / f"channel-wall_{face}_history.dat")
        nb = (w.shape[1] - 2) // 2
        taus.append(w[int(len(w) * 0.75) :, 2 : 2 + nb].mean())
    f = np.loadtxt(work / "channel-forcing_history.dat", ndmin=2, comments="VARIABLES")
    force = float(f[int(len(f) * 0.75) :, 2].mean())
    tau = 0.5 * (taus[0] + taus[1])
    rel = abs(tau - force) / abs(force)
    ok = rel < 1e-3
    print(f"wall history tau_w {tau:.6e} (y_min {taus[0]:.6e}, y_max {taus[1]:.6e}) against the force {force:.6e}: relative {rel:.2e} (bound 1e-3)  {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def wall_agree(a: str, b: str, tol: float) -> int:
    """Every wall-model history of two runs within tol of each column's scale."""
    ok = True
    files = sorted(Path(a).glob("*-wall_*_history.dat"))
    if not files:
        print("wall-agree: no wall history  FAIL")
        return 1
    for fa in files:
        x, y = wall_rows(fa), wall_rows(Path(b) / fa.name)
        if x.shape != y.shape:
            print(f"{fa.name}: shapes {x.shape} and {y.shape}  FAIL")
            ok = False
            continue
        d = float(np.max(np.abs(x - y)[:, 1:] / np.maximum(np.abs(x[:, 1:]).max(axis=0), 1e-300)))
        good = d < tol
        ok = ok and good
        print(f"{fa.name}: {d:.2e} (bound {tol:.0e})  {'PASS' if good else 'FAIL'}")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="what", required=True)
    p = sub.add_parser("profile")
    p.add_argument("src")
    p.add_argument("out")
    p = sub.add_parser("agree")
    p.add_argument("a")
    p.add_argument("b")
    p = sub.add_parser("inflow")
    p.add_argument("work")
    p = sub.add_parser("wall-force")
    p.add_argument("work")
    p = sub.add_parser("wall-agree")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--tol", type=float, default=1e-12)
    p = sub.add_parser("compare")
    p.add_argument("work")
    p.add_argument("stat")
    p.add_argument("--t0", type=float, default=0.0)
    p.add_argument("--json", default=None)
    p.add_argument("--re-inf", type=float, default=None, help="the free-stream Reynolds number per delta_i of Re_theta (default the input's)")
    a = ap.parse_args()
    if a.what == "profile":
        return profile(a.src, a.out)
    if a.what == "agree":
        return agree(a.a, a.b)
    if a.what == "wall-force":
        return wall_force(a.work)
    if a.what == "wall-agree":
        return wall_agree(a.a, a.b, a.tol)
    if a.what == "compare":
        return compare(a.work, a.stat, a.t0, a.json, a.re_inf)
    return inflow(a.work)


if __name__ == "__main__":
    sys.exit(main())
