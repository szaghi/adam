#!/usr/bin/env python3
"""Oracle of the EGLM conservation, EV-4 (issue #47, M3-P4c).

Why: EGLM (Derigs et al. 2018) adds nonconservative sources, -div B (0, B, u.B, u, 0) - (u.grad psi) (0, 0, psi, 0, 1)
in the order (rho, rho u, E, B, psi), so on a periodic box only the density integral is conserved; the momentum, B and
energy integrals move by the time integral of their sources. The oracle bounds each drift by the integral of the
magnitude of its source, measured from the run (owner decision 2026-10-01, the EV-4 row of #47 amended):

* rho: |int rho(t) - int rho(0)| <= --rho-tol max(|int rho(0)|, 1) (round-off);
* rho u (each component): <= K max|B| int_0^t ||div B||_1 dt;
* B (each component):     <= K max|u| int_0^t ||div B||_1 dt;
* E:                      <= K (max|u.B| int_0^t ||div B||_1 dt + int_0^t ||psi (u.grad psi)||_1 dt);

at every row of the conservation history, with K = --factor (default 2). ||div B||_1 = sum |div B| dV is the div(B)
history (every step, the same centred difference of order 2S as the source); the maxima of |B|, |u|, |u.B| are taken
over every saved checkpoint, and ||psi (u.grad psi)||_1 is computed on each checkpoint (second-order centred
differences on the saved ghost cells) and integrated in time by the trapezoidal rule, linearly interpolated to the
history rows: save
checkpoints often (e.g. every 25 steps). The ratio of each drift to its bound is reported.

Usage:
    eglm_conservation_oracle.py <work> [--factor K] [--rho-tol T] [--ngc N]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "divb-peak"))
from divb_peak_oracle import history  # noqa: E402


def checkpoint_terms(path_group: list[Path], ngc: int) -> tuple[float, float, float, float]:
    """Return max|B|, max|u|, max|u.B| and ||psi (u.grad psi)||_1 of one checkpoint (all its rank files)."""
    bmax = umax = ubmax = s_psi = 0.0
    for path in path_group:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                q = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)
                     for v in ("r", "ru", "rv", "rw", "bx", "by", "bz", "psi")}
                d = h5[f"{blk}-dxdydz"][()][::-1]
                inner = tuple(slice(ngc, -ngc) for _ in range(3))
                u = [q[m] / q["r"] for m in ("ru", "rv", "rw")]
                b = [q[m] for m in ("bx", "by", "bz")]
                psi = q["psi"]
                grad = []
                for a in range(3):  # centred differences on the saved ghost cells; a null axis has equal copies
                    g = np.zeros_like(psi)
                    sl_p = [slice(None)] * 3
                    sl_m = [slice(None)] * 3
                    sl_c = [slice(None)] * 3
                    sl_p[a], sl_m[a], sl_c[a] = slice(2, None), slice(None, -2), slice(1, -1)
                    g[tuple(sl_c)] = (psi[tuple(sl_p)] - psi[tuple(sl_m)]) / (2.0 * d[a])
                    grad.append(g)
                bmag = np.sqrt(sum(c**2 for c in b))[inner]
                umag = np.sqrt(sum(c**2 for c in u))[inner]
                ub = sum(uc * bc for uc, bc in zip(u, b, strict=True))[inner]
                upsi = sum(uc * gc for uc, gc in zip(u, grad, strict=True))[inner]
                bmax, umax, ubmax = max(bmax, bmag.max()), max(umax, umag.max()), max(ubmax, np.abs(ub).max())
                s_psi += float(np.sum(np.abs(psi[inner] * upsi)) * d[0] * d[1] * d[2])
    return bmax, umax, ubmax, s_psi


def main() -> int:
    """Run the EV-4 checks, print a report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path)
    parser.add_argument("--factor", type=float, default=2.0, help="safety factor K of the source bounds")
    parser.add_argument("--rho-tol", type=float, default=1.0e-13, help="relative bound of the density drift")
    parser.add_argument("--ngc", type=int, default=3)
    args = parser.parse_args()
    cons = history(args.work, "conservation")
    divb = history(args.work, "divb")
    files = sorted(p for p in args.work.glob("*-proc*.h5") if "restart" not in p.name)
    steps = sorted({int(p.name.split("-")[-2]) for p in files})
    if len(steps) < 3:
        sys.exit(f"eglm_conservation_oracle: {args.work}: {len(steps)} checkpoints, save more often")
    time_of = dict(zip(cons["it"].astype(int), cons["time"], strict=True))
    terms = np.array([checkpoint_terms([p for p in files if int(p.name.split("-")[-2]) == s], args.ngc) for s in steps])
    bmax, umax, ubmax = terms[:, 0].max(), terms[:, 1].max(), terms[:, 2].max()
    t_ck = np.array([time_of[s] for s in steps])
    s_psi = terms[:, 3]
    i_psi_ck = np.concatenate([[0.0], np.cumsum(0.5 * (s_psi[1:] + s_psi[:-1]) * np.diff(t_ck))])
    t = cons["time"]
    l1 = np.interp(t, divb["time"], divb["l1_divb"])
    i_divb = np.concatenate([[0.0], np.cumsum(0.5 * (l1[1:] + l1[:-1]) * np.diff(t))])
    i_psi = np.interp(t, t_ck, i_psi_ck)
    k = args.factor
    bounds = {"r": None,
              "ru": k * bmax * i_divb, "rv": k * bmax * i_divb, "rw": k * bmax * i_divb,
              "bx": k * umax * i_divb, "by": k * umax * i_divb, "bz": k * umax * i_divb,
              "rE": k * (ubmax * i_divb + i_psi)}
    print(f"{args.work.name}: {len(t)} rows, {len(steps)} checkpoints, t = {t[-1]:.6g}; max|B| {bmax:.3e}, max|u| "
          f"{umax:.3e}, max|u.B| {ubmax:.3e}; int ||div B||_1 dt {i_divb[-1]:.3e}, int ||psi u.grad psi||_1 dt "
          f"{i_psi[-1]:.3e}")
    ok = True
    for name, bound in bounds.items():
        col = cons[f"int_{name}"]
        drift = np.abs(col - col[0])
        if bound is None:
            limit = args.rho_tol * max(abs(col[0]), 1.0)
            good = bool(np.all(drift <= limit))
            print(f"   int_{name:2s}: max drift {drift.max():.3e} (bound {limit:.1e}, round-off)  "
                  f"{'PASS' if good else 'FAIL'}")
        else:
            good = bool(np.all(drift <= bound))
            ratio = float(np.max(np.divide(drift, bound, out=np.zeros_like(drift), where=bound > 0.0)))
            if np.any((bound == 0.0) & (drift > 0.0)):
                good, ratio = False, float("inf")
            print(f"   int_{name:2s}: max drift {drift.max():.3e}, final bound {bound[-1]:.3e}, max drift/bound "
                  f"{ratio:.2e}  {'PASS' if good else 'FAIL'}")
        ok = ok and good
    print(f"EV-4 EGLM conservation {'PASS' if ok else 'FAIL'} (K = {k:g})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
