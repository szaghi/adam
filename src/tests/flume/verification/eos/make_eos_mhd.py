#!/usr/bin/env python3
"""Write an isothermal or barotropic MHD input of the P5 verification (issue #100, M13 P5).

    make_eos_mhd.py <out.ini> tube --tube 1|2|3 [--eos isothermal|barotropic] [--cells N] [--scheme ...] [--solver ...]
    make_eos_mhd.py <out.ini> wave --wave fast|alfven|slow [--eos ...] [--cells N] [--angle A]
    make_eos_mhd.py <out.ini> cpaw [--cells N] [--polarisation right|left] [--angle A]

Why: the isothermal and barotropic MHD models (MODEL_MHD_ISO/BARO and their GLM twins) have no reference input; the
cases are built from the regression Sod input (as make_eos.py does for the hydro), with the gas of [eos] and the
magnetic field of the MHD models.

  tube   Mignone (2007, J. Comput. Phys. 225, Table 1) shock tubes 1-3 on [0, 1], the discontinuity at 0.5, a = 1,
         B = H / sqrt(4 pi): 1, fast and slow rarefactions and shocks; 2, rotational discontinuities as well; 3, the
         collision of two Mach-5 streams. The barotropic gas takes POLY of eos_oracle.py (c0 = 1, the first segment
         isothermal, so tubes of density below 2 see the isothermal gas: --rho-scale lifts them into the polytrope).
  wave   a linear wave of one family on the periodic box (the mhd-linear-wave IC with the eigenvectors of Kim et al.
         1999), amplitude 1e-6, one wavelength along the x axis or at --angle degrees in the x-y plane.
  cpaw   the circularly polarised Alfven wave (the mhd-cpaw IC), an exact nonlinear solution for any P(rho).
"""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
from eos_oracle import POLY  # noqa: E402
from make_eos import set_key, unset_key  # noqa: E402

TUBES = {  # Mignone (2007), Table 1: rho, u, v, w, Hx, Hy, Hz (H = sqrt(4 pi) B), t_s
    1: ((1.0, 0.0, 0.0, 0.0, 3.0, 5.0, 0.0), (0.1, 0.0, 0.0, 0.0, 3.0, 2.0, 0.0), 0.1),
    2: ((1.08, 1.2, 0.01, 0.5, 2.0, 3.6, 2.0), (1.0, 0.0, 0.0, 0.0, 2.0, 4.0, 2.0), 0.2),
    3: ((0.1, 5.0, 0.0, 0.0, 0.0, -1.0, -2.0), (0.1, -5.0, 0.0, 0.0, 0.0, 1.0, 2.0), 0.25),
}


def gas(s: str, eos: str, c0: float = 1.0, one_segment: bool = False) -> str:
    """The [eos] section of the gas, no gamma, cp, cv: the isothermal speed c0, or POLY of eos_oracle.py (c0 in place of
    its own), or with `one_segment` the one-segment table exponents = 1 (the isothermal limit)."""
    for key in ("cp", "cv", "gamma"):
        s = unset_key(s, "physics", key)
    s = set_key(s, "eos", "type", eos)
    s = set_key(s, "eos", "sound_speed", repr(c0))
    if eos == "barotropic":
        if one_segment:
            s = set_key(s, "eos", "exponents", "1.0")
        else:
            s = set_key(s, "eos", "densities", " ".join(repr(b) for b in POLY.breaks))
            s = set_key(s, "eos", "exponents", " ".join(repr(g) for g in POLY.exponents))
        s = set_key(s, "eos", "reference_density", repr(POLY.rho_ref))
    return s


def mhd(s: str, dc: str, nv: int) -> str:
    s = set_key(s, "physics", "physical_model", "mhd-ideal")
    s = set_key(s, "mhd", "divergence_control", dc)
    s = set_key(s, "mhd", "divb_tol", "0.0")
    s = set_key(s, "mhd", "divb_error", ".false.")
    s = set_key(s, "mhd", "rho_floor", "0.0")  # no floors with an inert energy slot (they are refused above 0)
    s = set_key(s, "mhd", "p_floor", "0.0")
    if dc == "glm":
        for key, val in (("glm_ch", "8.0"), ("glm_alpha", "0.18"), ("glm_damping_length", "1.0"),
                         ("glm_ch_check", "warning")):
            s = set_key(s, "mhd", key, val)
    return set_key(s, "field", "nv", str(nv + (dc == "glm")))


def numerics(s: str, a: argparse.Namespace) -> str:
    s = set_key(s, "numerics", "scheme_space", a.scheme)
    s = set_key(s, "numerics", "reconstruction_variables", a.recon)
    if a.scheme == "weno-riemann":
        s = set_key(s, "numerics", "riemann_solver", a.solver)
        s = set_key(s, "numerics", "flux_correction", "6th")
        s = set_key(s, "numerics", "flux_correction_sensor", "weno")
    return set_key(s, "time", "CFL", repr(a.cfl))


AXES = {"x": (0, 1, 2), "y": (2, 0, 1), "z": (1, 2, 0)}  # (u, v, w) of the tube state along x -> along the axis


def tube(a: argparse.Namespace) -> str:
    """The tube along --axis on the verification Sod grid of that axis (four blocks along it): the x states permuted
    cyclically (along y: (u, v, w) -> (w, u, v), as the frame of the MHD libraries), so the three runs are the same
    problem in the frames of the libraries (bitwise equal, issue #41 MV-4)."""
    left, right, t_end = TUBES[a.tube]
    base = Path(__file__).resolve().parents[1] / "sod" / f"sod-{a.axis}.ini"
    s = gas(base.read_text(), a.eos, a.c0, a.one_segment)
    s = mhd(s, a.dc, 8)
    s = numerics(s, a)
    s = set_key(s, "grid", {"x": "ni", "y": "nj", "z": "nk"}[a.axis], str(a.cells // 4))
    s = set_key(s, "time", "time_max", repr(a.time_max or t_end))
    h = math.sqrt(4.0 * math.pi)
    perm = AXES[a.axis]
    for k, st in ((1, left), (2, right)):
        sec = f"initial_conditions_region_{k}"
        s = unset_key(s, sec, "p")
        s = set_key(s, sec, "r", repr(st[0] * a.rho_scale))
        vel, bfl = st[1:4], [b / h for b in st[4:7]]
        for key, idx in zip(("u", "v", "w"), perm, strict=True):
            s = set_key(s, sec, key, repr(vel[idx]))
        for key, idx in zip(("bx", "by", "bz"), perm, strict=True):
            s = set_key(s, sec, key, repr(bfl[idx]))
    return s


def cpaw(a: argparse.Namespace) -> str:
    """make_cpaw.py's input (Toth 2000: rho = 1, B_par = 1, B_perp = 0.1, the Alfven speed and the period 1, any P(rho))
    with the gas of [eos]."""
    verif = Path(__file__).resolve().parents[1]
    base = verif / ("sod/sod-x.ini" if a.geometry == "1d" else "vortex/vortex-n064.ini")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "cpaw.ini"
        subprocess.run([sys.executable, str(verif / "mhd/cpaw/make_cpaw.py"), str(base), str(out), "--polarisation",
                        a.polarisation, "--geometry", a.geometry, "--cells", str(a.cells)], check=True)
        s = out.read_text()
    s = gas(s, a.eos)
    s = mhd(s, a.dc, 8)
    s = numerics(s, a)
    s = set_key(s, "time", "CFL", "0.1")
    return unset_key(s, "initial_conditions_region_1", "p")


def wave_speed(eos: str, rho0: float, family: str) -> float:
    """|u_n + c| of the right-going wave on the background of make_linear_wave.py, B = (1, sqrt(2), 1/2) in the wave
    frame, at rest, density rho0 (the mhd_oracle speeds: a^2 = dP/drho)."""
    from eos_oracle import isothermal  # noqa: PLC0415
    from mhd_oracle import kim_eigensystem  # noqa: PLC0415

    gas_law = POLY if eos == "barotropic" else isothermal(1.0)
    w = np.array([rho0, 0.0, 0.0, 0.0, math.sqrt(2.0), 0.5])
    lam, _, _ = kim_eigensystem(gas_law, w, 1.0)
    return float({"fast": lam[5], "alfven": lam[4], "slow": lam[3]}[family])


def wave(a: argparse.Namespace) -> str:
    """make_linear_wave.py's input with the gas of [eos] and the background density rho0 (the period from its speed)."""
    verif = Path(__file__).resolve().parents[1]
    base = verif / ("sod/sod-x.ini" if a.geometry == "1d" else "vortex/vortex-n064.ini")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "lw.ini"
        subprocess.run([sys.executable, str(verif / "mhd/linear-wave/make_linear_wave.py"), str(base), str(out),
                        "--wave", a.wave, "--geometry", a.geometry, "--cells", str(a.cells), "--amplitude",
                        a.amplitude, "--divergence-control", a.dc], check=True)
        s = out.read_text()
    s = gas(s, a.eos)
    s = set_key(s, "mhd", "rho_floor", "0.0")
    s = set_key(s, "mhd", "p_floor", "0.0")
    s = set_key(s, "initial_conditions_region_1", "r", repr(a.rho0))
    s = unset_key(s, "initial_conditions_region_1", "p")
    speed = wave_speed(a.eos, a.rho0, a.wave)
    s = set_key(s, "time", "time_max", repr(a.time_max or 1.0 / speed))
    if a.dc == "glm":
        s = set_key(s, "mhd", "glm_ch", repr(1.1 * wave_speed(a.eos, a.rho0, "fast")))
    return s


def thin(s: str) -> str:
    """One cell per block along every null axis: the base inputs carry four, copies of the same solution that cost 16
    times the work in 1-D (bitwise the same profile, tube 2 checked)."""
    for ax, key in (("x", "ni"), ("y", "nj"), ("z", "nk")):
        if re.search(rf"(?m)^null_{ax}\s*=\s*\.true\.", s):
            s = set_key(s, "grid", key, "1")
    return s


def case_ini(a: argparse.Namespace) -> str:
    s = thin({"tube": tube, "wave": wave, "cpaw": cpaw}[a.case](a))
    s = set_key(s, "IO", "it_save", str(a.it_save))
    for item in a.set:
        sk, value = item.split("=", 1)
        section, key = sk.split(":", 1)
        s = set_key(s, section, key, value)
    return s


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("case", choices=("tube", "wave", "cpaw"))
    ap.add_argument("--tube", type=int, choices=(1, 2, 3), default=2)
    ap.add_argument("--axis", choices=("x", "y", "z"), default="x")
    ap.add_argument("--polarisation", choices=("right", "left"), default="right")
    ap.add_argument("--eos", choices=("isothermal", "barotropic"), default="isothermal")
    ap.add_argument("--c0", type=float, default=1.0, help="tube: the sound speed (POLY's c0 with the barotropic gas)")
    ap.add_argument("--one-segment", action="store_true", help="tube, barotropic: the table exponents = 1")
    ap.add_argument("--rho-scale", type=float, default=1.0)
    ap.add_argument("--wave", choices=("fast", "alfven", "slow"), default="fast")
    ap.add_argument("--geometry", choices=("1d", "2d"), default="1d")
    ap.add_argument("--amplitude", default="1.0e-7")
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--cells", type=int, default=400)
    ap.add_argument("--dc", choices=("none", "glm"), default="none")
    ap.add_argument("--scheme", choices=("weno", "weno-riemann"), default="weno")
    ap.add_argument("--solver", choices=("llf", "hll", "hlld"), default="hlld")
    ap.add_argument("--recon", choices=("characteristic", "conservative", "primitive"), default="characteristic")
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--time-max", type=float, default=None)
    ap.add_argument("--it-save", type=int, default=1000000)
    ap.add_argument("--set", action="append", default=[], help="section:key=value")
    a = ap.parse_args()
    a.out.write_text(case_ini(a))
    return 0


if __name__ == "__main__":
    sys.exit(main())
