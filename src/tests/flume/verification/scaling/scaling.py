#!/usr/bin/env python3
"""FLUME scaling-covariance oracle (issue #49, N0): rescale an input by powers of two and compare the runs bitwise.

Why: ideal Euler and MHD in FLUME's units (B = B_SI/sqrt(mu0)) carry no dimensionless number, so a run whose input is
rescaled (lengths by 2^j, velocities by 2^k, density by 4^m) is the same run in other units: density scales by 4^m,
momentum by 4^m 2^k, pressure and energy by 4^(k+m), B by 2^(k+m), the mixed-GLM psi by 2^(2k+m), the EGLM psi by
2^(k+m), time by 2^(j-k). Multiplying by a power of two is exact in floating point (away from overflow and underflow), so
a code without absolute constants reproduces the base run BIT FOR BIT after the rescaling; any absolute tolerance (the
WENO zeps, the positivity floor, ...) shows up as a mismatch. This is the measurement N1 of #49 must turn into a pass.

  scaling.py rescale <in.ini> <out.ini> --j J --k K --m M [--weights js|si]
      write the rescaled input; refuses an initial-condition, boundary or [mhd] key it cannot classify, so nothing
      dimensional passes unscaled (an initial condition without parameters, e.g. orszag-tang, cannot be rescaled);
      --weights sets [weno] weights (N1: si, the scale-invariant WENO weights);
  scaling.py compare <base-work> <scaled-work> --j J --k K --m M --ngc N [--psi glm|eglm]
      compare the last common checkpoint: each field of the scaled run against the base field times its exact power
      of two, blocks keyed by their scaled origins; prints the max relative difference per field (0 = bitwise) and
      exits 1 unless every field is bitwise.
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

# dimension of a quantity as the exponents (length, velocity, density) of its scale 2^(a j + b k + c 2m)
LENGTH, VELOCITY, DENSITY = (1, 0, 0), (0, 1, 0), (0, 0, 1)
PRESSURE = (0, 2, 1)       # rho u^2
FIELD = (0, 1, 0.5)        # B = u sqrt(rho)
TIME = (1, -1, 0)
INV_LENGTH = (-1, 0, 0)
NONE = (0, 0, 0)
GRID_KEYS = {f"{e}_{x}": LENGTH for e in ("emin", "emax") for x in "xyz"}
REGION_KEYS = {**GRID_KEYS, "r": DENSITY, "u": VELOCITY, "v": VELOCITY, "w": VELOCITY, "p": PRESSURE,
               "bx": FIELD, "by": FIELD, "bz": FIELD}
IC_KEYS = {"type": NONE, "amr_iterations": NONE, "regions_number": NONE, "axis": NONE, "pulse_axis": NONE,
           "wave": NONE, "polarisation": NONE, "s": NONE,
           "x0": LENGTH, "y0": LENGTH, "radius": LENGTH, "r0": LENGTH, "r1": LENGTH, "strength": VELOCITY,
           "kappa": VELOCITY, "mu": FIELD, "v0": VELOCITY, "rho_in": DENSITY, "pulse_center": LENGTH,
           "pulse_width": LENGTH, "pulse_amplitude": FIELD, "peak_x0": LENGTH, "peak_y0": LENGTH,
           "peak_radius": LENGTH, "peak_amplitude": FIELD, "loop_radius": LENGTH, "loop_amplitude": (1, 1, 0.5),
           "interface": LENGTH, "interface_1": LENGTH, "interface_2": LENGTH, "interface_2_width": LENGTH,
           "period": LENGTH, "normal_x": NONE, "normal_y": NONE, "rho_amplitude": NONE, "rho_wavenumber": INV_LENGTH,
           "wavelength": LENGTH, "wave_angle": NONE, "wave_amplitude": NONE, "b_par": FIELD}
MHD_KEYS = {"divergence_control": NONE, "divb_error": NONE, "glm_alpha": NONE, "glm_ch_check": NONE,
            "glm_ch": VELOCITY, "glm_damping_length": LENGTH, "rho_floor": DENSITY, "p_floor": PRESSURE,
            "divb_tol": (-1, 1, 0.5)}
BOX_KEYS = {f"box_{x}{e}": LENGTH for x in "xyz" for e in ("min", "max")}
NO_IC = ("orszag-tang",)
GEOMETRY = ("origin", "dxdydz", "time_iteration")


def factor(dim: tuple[float, float, float], j: int, k: int, m: int) -> float:
    """Exact power-of-two scale of a quantity of dimension `dim` (m counts powers of 4 of the density)."""
    e = dim[0] * j + dim[1] * k + dim[2] * 2 * m
    if e != int(e):
        raise SystemExit(f"scaling: non-integer power of two {e} (dimension {dim})")
    return 2.0 ** int(e)


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read an INI keeping the key case."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.optionxform = str
    ini.read(path)
    return ini


def rescale(args: argparse.Namespace) -> int:
    """Write the rescaled input."""
    ini = read_ini(args.ini)
    j, k, m = args.j, args.k, args.m
    if ini["initial_conditions"]["type"].strip() in NO_IC:
        raise SystemExit(f"scaling: initial condition '{ini['initial_conditions']['type'].strip()}' has no parameters "
                         "and cannot be rescaled from the input")
    tables = {"grid": GRID_KEYS, "initial_conditions": IC_KEYS, "mhd": MHD_KEYS}
    for section in ini.sections():
        if section.startswith("initial_conditions_region_"):
            table = REGION_KEYS
        elif section.startswith("amr_marker_"):
            if ini[section].get("geo_type", "").strip() != "primitive-box" or int(ini[section].get("mode", "0")) != 1:
                raise SystemExit(f"scaling: [{section}] is not a geometric box marker; not handled")
            table = {**{key: NONE for key in ini[section]}, **BOX_KEYS}
        elif section.startswith("bc_"):
            if any(key != "type" for key in ini[section]):
                raise SystemExit(f"scaling: [{section}] carries a state; boundary states are not handled")
            continue
        elif section in tables:
            table = tables[section]
        else:
            continue
        for key in ini[section]:
            if section == "grid" and key not in GRID_KEYS:
                continue
            if key not in table:
                raise SystemExit(f"scaling: [{section}] {key} is not classified; refusing to leave it unscaled")
            if table[key] == NONE:
                continue
            value = float(ini[section][key].split(";")[0])
            ini[section][key] = repr(value * factor(table[key], j, k, m))
    t = ini["time"]
    t["time_max"] = repr(float(t["time_max"].split(";")[0]) * factor(TIME, j, k, m))
    if args.weights is not None:
        ini["weno"]["weights"] = args.weights
    with open(args.out, "w") as f:
        ini.write(f)
    return 0


def field_scale(name: str, psi: str, j: int, k: int, m: int) -> float:
    """Exact scale of a checkpoint field."""
    dims = {"r": DENSITY, "ru": (0, 1, 1), "rv": (0, 1, 1), "rw": (0, 1, 1), "rE": PRESSURE,
            "bx": FIELD, "by": FIELD, "bz": FIELD}
    if name == "psi":
        return factor((0, 2, 0.5) if psi == "glm" else FIELD, j, k, m)
    if name not in dims:
        raise SystemExit(f"scaling: no scale for field '{name}'")
    return factor(dims[name], j, k, m)


def load(work: Path, step: int, ngc: int) -> tuple[tuple[str, ...], dict]:
    """Interior cells of every checkpoint file of a step: field names and arrays keyed by block origin."""
    blocks = {}
    names: tuple[str, ...] = ()
    for path in sorted(work.glob(f"*-{step:09d}-proc*.h5")):
        with h5py.File(path, "r") as h5:
            if not len(h5):  # a rank that owns no block
                continue
            found = tuple(sorted({key.rsplit("-", 1)[1] for key in h5} - set(GEOMETRY)))
            names = names or found
            for blk in {key.rsplit("-", 1)[0] for key in h5}:
                q = np.stack([h5[f"{blk}-{v}"][()] for v in names])[:, ngc:-ngc, ngc:-ngc, ngc:-ngc]
                blocks[tuple(float(o) for o in h5[f"{blk}-origin"][()])] = q
    return names, blocks


def steps(work: Path) -> set[int]:
    """Checkpoint steps of a work directory (restart files excluded)."""
    return {int(p.name.split("-")[-2]) for p in work.glob("*-[0-9]*-proc*.h5")}


def compare(args: argparse.Namespace) -> int:
    """Compare the scaled run with the base run times the exact scales; 0 if bitwise."""
    j, k, m = args.j, args.k, args.m
    last_base, last_scaled = max(steps(args.base)), max(steps(args.scaled))
    step = max(steps(args.base) & steps(args.scaled))
    names, base = load(args.base, step, args.ngc)
    snames, scaled = load(args.scaled, step, args.ngc)
    lscale = factor(LENGTH, j, k, m)
    keyed = {tuple(o * lscale for o in origin): q for origin, q in base.items()}
    if names != snames or set(keyed) != set(scaled):
        print(f"   step {step}: blocks or fields differ ({len(keyed)} vs {len(scaled)} blocks)  FAIL")
        return 1
    status = int(last_base != last_scaled)
    worst = []
    for v, name in enumerate(names):
        s = field_scale(name, args.psi, j, k, m)
        a = np.concatenate([keyed[o][v].ravel() * s for o in sorted(keyed)])
        b = np.concatenate([scaled[o][v].ravel() for o in sorted(keyed)])
        scale = max(float(np.max(np.abs(a))), float(np.finfo(float).tiny))
        rel = float(np.max(np.abs(a - b))) / scale
        status |= int(rel != 0.0)
        worst.append(f"{name} {rel:.1e}")
    steps_note = "" if last_base == last_scaled else f" (last steps differ: base {last_base}, scaled {last_scaled})"
    print(f"   step {step}{steps_note}: {'BITWISE' if status == 0 else 'not bitwise'}; max relative difference "
          + ", ".join(worst))
    return status


def main() -> int:
    """Dispatch the subcommand."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rescale")
    r.add_argument("ini", type=Path)
    r.add_argument("out", type=Path)
    r.add_argument("--weights", choices=("js", "si"), help="set [weno] weights")
    c = sub.add_parser("compare")
    c.add_argument("base", type=Path)
    c.add_argument("scaled", type=Path)
    c.add_argument("--ngc", type=int, required=True)
    c.add_argument("--psi", choices=("glm", "eglm"), default="glm")
    for p in (r, c):
        p.add_argument("--j", type=int, default=0, help="lengths times 2^j")
        p.add_argument("--k", type=int, default=0, help="velocities times 2^k")
        p.add_argument("--m", type=int, default=0, help="density times 4^m")
    args = parser.parse_args()
    return rescale(args) if args.cmd == "rescale" else compare(args)


if __name__ == "__main__":
    sys.exit(main())
