#!/usr/bin/env python3
"""Write the inputs of the FLUME equation-of-state verification (issue #100, M13) from the sod-x regression input.

    make_eos.py <out.ini> [--case tube|collide|wave] [--ni N] [--time-max T] [--cfl C] [--c0 C0] [--mach M]
                          [--scheme weno|weno-riemann] [--recon characteristic|conservative|primitive]
                          [--solver llf|hll|hllc] [--rk SCHEME] [--keep-gas] [--null-y] [--level L]
                          [--eos isothermal|barotropic] [--densities "D1 D2 .."] [--exponents "G1 G2 .."]
                          [--rho-ref R] [--rho0 R0] [--left "r u v"] [--right "r u v"]
                          [--set SECTION:KEY=VALUE ...] [--unset SECTION:KEY ...]

The sod-x grid with `--ni` cells along x (default 400; one block, two for the periodic wave, 2^L with `--level L`),
the other directions null, the isothermal gas `[eos] type = isothermal, sound_speed = C0` (default 1) in place of
`cp, cv` (`--keep-gas` keeps them: the contract leg), and no region pressure (P = c0^2 rho is not an input). The cases:

  tube     the isothermal shock tube: rho 1 | 0.125 at x = 0.5, at rest, the transverse velocity v 0.3 | -0.2 (the
           contact carries it), extrapolation boundaries, to t = 0.2 (default): a left rarefaction, the contact, a
           right shock. The y direction is active (nj = 4 cells, the state uniform in y, so its fluxes cancel): a null
           direction freezes its momentum (`null_freeze`), which would freeze v (`--null-y` keeps it null);
  collide  two streams rho = 1, u = +-M (default 5) colliding at x = 0.5, to t = 0.2: two isothermal shocks, the
           dense layer at rest between them, compressed by M_1^2;
  wave     a right-running acoustic wave of amplitude 1e-6 and wavelength 1 on the periodic box [0, 1] (`sine-wave`,
           `wave_mode = acoustic`), to one period t = 1 / c0 (default).

`--eos barotropic` writes the piecewise polytrope `[eos] densities, exponents, reference_density` (P3; without the
lists the code's default, the Mathew & Federrath 2021 table); `--left` / `--right` override the Riemann states
(`r u v`, the leading ones given), `--rho0` the wave's background density.

The space scheme (`--scheme`, default weno), reconstruction (`--recon`, default characteristic), Riemann solver
(`--solver`, weno-riemann only, with the 6th-order flux correction and its WENO sensor) and Runge-Kutta scheme
(`--rk`) are set as given. `--set` writes a key (the section is created when missing), `--unset` removes one: the
contract cases.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2] / "regression" / "sod-x" / "input.ini"


def set_key(s: str, section: str, key: str, value: str) -> str:
    """Set `key = value` in `[section]` (append the section when it is missing)."""
    head = re.search(rf"(?m)^\[{re.escape(section)}\]\s*$", s)
    if head is None:
        return s.rstrip("\n") + f"\n\n[{section}]\n{key} = {value}\n"
    nxt = re.search(r"(?m)^\[", s[head.end() :])
    end = head.end() + nxt.start() if nxt else len(s)
    body = s[head.end() : end]
    if re.search(rf"(?m)^{re.escape(key)}\s*=", body):
        body = re.sub(rf"(?m)^{re.escape(key)}\s*=.*$", f"{key} = {value}", body)
    else:
        body = f"\n{key} = {value}" + body
    return s[: head.end()] + body + s[end:]


def unset_key(s: str, section: str, key: str) -> str:
    """Remove `key` from `[section]`."""
    head = re.search(rf"(?m)^\[{re.escape(section)}\]\s*$", s)
    assert head is not None, section
    nxt = re.search(r"(?m)^\[", s[head.end() :])
    end = head.end() + nxt.start() if nxt else len(s)
    body = re.sub(rf"(?m)^{re.escape(key)}\s*=.*\n?", "", s[head.end() : end])
    return s[: head.end()] + body + s[end:]


def case_ini(a: argparse.Namespace) -> str:
    """Return the INI of a case."""
    s = BASE.read_text()
    if not a.keep_gas:
        s = re.sub(r"(?m)^cp\s*=.*\n^cv\s*=.*\n", "", s)
    s = set_key(s, "eos", "type", a.eos)
    s = set_key(s, "eos", "sound_speed", repr(a.c0))
    if a.eos == "barotropic":  # without --densities/--exponents the code's default table (Mathew & Federrath 2021)
        if a.densities is not None:
            s = set_key(s, "eos", "densities", a.densities)
        if a.exponents is not None:
            s = set_key(s, "eos", "exponents", a.exponents)
        if a.rho_ref is not None:
            s = set_key(s, "eos", "reference_density", repr(a.rho_ref))
    # one block along x; the periodic box two (a single block cannot be its own periodic neighbour); --level L refines
    # the base block L times (2^L blocks along x: a regression case keeps both of its ranks busy)
    level = a.level if a.level is not None else (1 if a.case == "wave" else 0)
    assert a.ni % 2**level == 0, f"--ni {a.ni}: not divisible by the 2^{level} blocks along x"
    s = set_key(s, "grid", "ni", str(a.ni // 2**level))
    s = set_key(s, "amr", "max_level", str(level))
    s = set_key(s, "amr", "iu_ref_levels", str(level))
    s = set_key(s, "time", "time_max", repr(a.time_max))
    s = set_key(s, "time", "CFL", repr(a.cfl))
    s = set_key(s, "numerics", "scheme_space", a.scheme)
    s = set_key(s, "numerics", "reconstruction_variables", a.recon)
    if a.scheme == "weno-riemann":
        s = set_key(s, "numerics", "riemann_solver", a.solver)
        s = set_key(s, "numerics", "flux_correction", "6th")
        s = set_key(s, "numerics", "flux_correction_sensor", "weno")
    if a.rk:
        s = set_key(s, "runge_kutta", "scheme", a.rk)
    for k in (1, 2):
        s = unset_key(s, f"initial_conditions_region_{k}", "p")
    if a.case == "tube":
        s = set_key(s, "initial_conditions_region_2", "r", "0.125")
        s = set_key(s, "initial_conditions_region_1", "v", "0.3")
        s = set_key(s, "initial_conditions_region_2", "v", "-0.2")
        s = set_key(s, "grid", "null_y", ".true." if a.null_y else ".false.")
    elif a.case == "collide":
        s = set_key(s, "initial_conditions_region_2", "r", "1.0")
        s = set_key(s, "initial_conditions_region_1", "u", repr(a.mach * a.c0))
        s = set_key(s, "initial_conditions_region_2", "u", repr(-a.mach * a.c0))
    else:  # the acoustic wave on the periodic box
        s = set_key(s, "initial_conditions", "type", "sine-wave")
        s = set_key(s, "initial_conditions", "regions_number", "1")
        s = set_key(s, "initial_conditions", "wave_mode", "acoustic")
        s = set_key(s, "initial_conditions", "wave_angle", "0.0")
        s = set_key(s, "initial_conditions", "wave_amplitude", "1.0e-6")
        s = set_key(s, "initial_conditions", "wavelength", "1.0")
        s = set_key(s, "initial_conditions_region_1", "emax_x", "1.0")
        s = set_key(s, "initial_conditions_region_1", "r", repr(a.rho0))
        s = set_key(s, "initial_conditions_region_2", "r", repr(a.rho0))
        for face in ("bc_x_min", "bc_x_max"):
            s = set_key(s, face, "type", "periodic")
    if a.case != "wave":  # --left/--right "r u v": the Riemann states (missing values keep the case's)
        for k, st in ((1, a.left), (2, a.right)):
            for key, val in zip(("r", "u", "v"), (st or "").split(), strict=False):
                s = set_key(s, f"initial_conditions_region_{k}", key, val)
    for item in a.set:
        sk, value = item.split("=", 1)
        section, key = sk.split(":", 1)
        s = set_key(s, section, key, value)
    for item in a.unset:
        section, key = item.split(":", 1)
        s = unset_key(s, section, key)
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--case", choices=("tube", "collide", "wave"), default="tube")
    ap.add_argument("--ni", type=int, default=400)
    ap.add_argument("--time-max", type=float, default=None)
    ap.add_argument("--cfl", type=float, default=0.5)
    ap.add_argument("--c0", type=float, default=1.0)
    ap.add_argument("--mach", type=float, default=5.0)
    ap.add_argument("--scheme", choices=("weno", "weno-riemann"), default="weno")
    ap.add_argument("--recon", choices=("characteristic", "conservative", "primitive"), default="characteristic")
    ap.add_argument("--solver", choices=("llf", "hll", "hllc"), default="hllc")
    ap.add_argument("--rk", default=None)
    ap.add_argument("--keep-gas", action="store_true")
    ap.add_argument("--null-y", action="store_true")
    ap.add_argument("--level", type=int, default=None)
    ap.add_argument("--eos", choices=("isothermal", "barotropic"), default="isothermal")
    ap.add_argument("--densities", default=None)
    ap.add_argument("--exponents", default=None)
    ap.add_argument("--rho-ref", type=float, default=None)
    ap.add_argument("--rho0", type=float, default=1.0)
    ap.add_argument("--left", default=None)
    ap.add_argument("--right", default=None)
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--unset", action="append", default=[])
    a = ap.parse_args()
    if a.time_max is None:
        a.time_max = 1.0 / a.c0 if a.case == "wave" else 0.2
    a.out.write_text(case_ini(a))


if __name__ == "__main__":
    main()
