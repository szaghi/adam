#!/usr/bin/env python3
"""FLUME scaling-covariance oracle (issue #49, N0): rescale an input by powers of two and compare the runs bitwise.

Why: ideal Euler and MHD in FLUME's units (B = B_SI/sqrt(mu0)) carry no dimensionless number, so a run whose input is
rescaled (lengths by 2^j, velocities by 2^k, density by 4^m) is the same run in other units: density scales by 4^m,
momentum by 4^m 2^k, pressure and energy by 4^(k+m), B by 2^(k+m), the mixed-GLM psi by 2^(2k+m), the EGLM psi by
2^(k+m), time by 2^(j-k). Multiplying by a power of two is exact in floating point (away from overflow and underflow), so
a code without absolute constants reproduces the base run BIT FOR BIT after the rescaling; any absolute tolerance (the
WENO zeps, the positivity floor, ...) shows up as a mismatch. This is the measurement N1 of #49 must turn into a pass.

  scaling.py rescale <in.ini> <out.ini> --j J --k K --m M [--weights js|si]
      write the rescaled input; every option FLUME reads is classified (classify(), the table mirroring
      adam_flume_reference_object.F90) and an unclassified one is refused, so nothing dimensional passes unscaled (an
      initial condition without parameters, e.g. orszag-tang, and a forest manifest cannot be rescaled);
      --weights sets [weno] weights (N1: si, the scale-invariant WENO weights);
  scaling.py dimensionalize <in.ini> <out.ini> --j J --k K --m M
      write the same problem in the units L0 = 2^j, u0 = 2^k, rho0 = 4^m with a [reference] section that converts it
      back (issue #49, N2, NV-5): the solver sees the base numbers, so the run must equal the base run bit for bit
      (any weights), and a key the Fortran layer misclassifies shows up as a mismatch. Exit 3 (nothing written) for
      an input the reference layer refuses: orszag-tang, a forest manifest;
  scaling.py check-log <base.ini> <reference-run-log>
      check the conversions the reference layer logged (rank 0): every converted value equals the base value exactly
      (cp, cv: gamma equals cp/cv) and every dimensional option of the base input was converted; exit 1 otherwise. It
      covers the options whose conversion does not show in the solution (e.g. a gradient AMR tolerance far from the
      gradients, an unused one);
  scaling.py compare <base-work> <scaled-work> --j J --k K --m M --ngc N [--psi glm|eglm] [--conservative]
      compare the last common checkpoint: each field of the scaled run against the base field times its exact power
      of two, blocks keyed by their scaled origins; prints the max relative difference per field (0 = bitwise) and
      exits 1 unless every field is bitwise (--conservative: the conservative fields only, e.g. when the temperature
      unit differs).
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
           "peak_radius": LENGTH, "peak_amplitude": FIELD, "loop_radius": LENGTH, "loop_amplitude": FIELD,
           "interface": LENGTH, "interface_1": LENGTH, "interface_2": LENGTH, "interface_2_width": LENGTH,
           "period": LENGTH, "normal_x": NONE, "normal_y": NONE, "rho_amplitude": DENSITY, "rho_wavenumber": INV_LENGTH,
           "wavelength": LENGTH, "wave_angle": NONE, "wave_amplitude": NONE, "b_par": FIELD}
MHD_KEYS = {"divergence_control": NONE, "divb_error": NONE, "glm_alpha": NONE, "glm_ch_check": NONE,
            "glm_ch": VELOCITY, "glm_damping_length": LENGTH, "rho_floor": DENSITY, "p_floor": PRESSURE,
            "divb_tol": (-1, 1, 0.5)}
BOX_KEYS = {f"box_{x}{e}": LENGTH for x in "xyz" for e in ("min", "max")}
NO_IC = ("orszag-tang",)
# dimensions of the reference layer (dimensionalize), beyond the ones above
MOMENTUM = (0, 1, 1)
VELOCITY2 = (0, 2, 0)
PSI_GLM = (0, 2, 0.5)
FIELD_PER_LENGTH = (-1, 1, 0.5)
STATE_KEYS = {"r": DENSITY, "u": VELOCITY, "v": VELOCITY, "w": VELOCITY, "p": PRESSURE, "bx": FIELD, "by": FIELD,
              "bz": FIELD}
DIMENSIONLESS_SECTIONS = ("numerics", "runge_kutta", "weno", "linear-algebra", "fdv", "field", "amr", "solids",
                          "slices", "diagnostics")
FULL_KEYS = {
    "physics": {"physical_model": NONE, "gamma": NONE, "cp": NONE, "cv": NONE},  # cp, cv: the layer makes gamma
    "mhd": {"divergence_control": NONE, "glm_alpha": NONE, "glm_ch_check": NONE, "divb_error": NONE,
            "glm_ch": VELOCITY, "glm_damping_length": LENGTH, "rho_floor": DENSITY, "p_floor": PRESSURE,
            "divb_tol": FIELD_PER_LENGTH},
    "IO": {**{key: NONE for key in ("output_basename", "it_save", "restart", "restart_basename", "restart_save",
                                    "residuals_save", "divergence_history_save", "save_memory_status",
                                    "save_residual_fields", "save_curl_fields", "save_divergence_fields",
                                    "save_gradient_fields", "save_laplacian_fields", "seam_divB_error",
                                    "save_auxiliary_fields")},
           "seam_divB_tol": FIELD_PER_LENGTH},
    "time": {"it_max": NONE, "CFL": NONE, "time_max": TIME},
    "grid": {**{key: NONE for key in ("ni", "nj", "nk", "ngc", "null_x", "null_y", "null_z")}, **GRID_KEYS},
}
BC_SECTIONS = tuple(f"bc_{x}_{e}" for x in "xyz" for e in ("min", "max"))
MARKER_NONE = ("mode", "delta_type", "geo_type", "solid", "target_level", "stl_filename", "field", "var")
SOLID_LENGTHS = tuple(f"{s}_center_{x}" for s in ("sphere", "circle", "rectangle") for x in "xyz") + \
    ("sphere_radius", "circle_radius", "rectangle_edge_1", "rectangle_edge_2")
GEOMETRY = ("origin", "dxdydz", "time_iteration")
CONSERVATIVE = ("r", "ru", "rv", "rw", "rE", "bx", "by", "bz", "psi")


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


def rescale_ini(args: argparse.Namespace) -> configparser.ConfigParser | None:
    """Return the input with every dimensional option rescaled (lengths 2^j, velocities 2^k, density 4^m), or None
    for an input that cannot be rescaled (an initial condition without parameters, a forest manifest)."""
    ini = read_ini(args.ini)
    j, k, m = args.j, args.k, args.m
    if ini.has_section("forest") or ini.get("initial_conditions", "type", fallback="").split(";")[0].strip() in NO_IC:
        return None
    for section in ini.sections():
        for key in ini[section]:
            dim = classify(ini, section, key)
            if dim == NONE:
                continue
            text = ini[section][key].split(";")[0].strip()
            if section == "mhd" and key == "glm_damping_length" and text == "min-cell":
                continue
            ini[section][key] = repr(float(text) * factor(dim, j, k, m))
    return ini


def rescale(args: argparse.Namespace) -> int:
    """Write the rescaled input."""
    ini = rescale_ini(args)
    if ini is None:
        raise SystemExit(f"scaling: {args.ini} cannot be rescaled from the input (forest manifest, or an initial "
                         "condition without parameters, e.g. orszag-tang)")
    if args.weights is not None:
        ini["weno"]["weights"] = args.weights
    with open(args.out, "w") as f:
        ini.write(f)
    return 0

def marker_tol(ini: configparser.ConfigParser, section: str) -> tuple[float, float, float]:
    """Dimension of the gradient tolerance of an AMR marker: the marked variable per length (NONE if unused)."""
    sec = ini[section]
    if int(sec.get("mode", "0").split(";")[0]) != 2:
        return NONE
    field, var = int(sec["field"].split(";")[0]), int(sec["var"].split(";")[0])
    glm = ini.get("mhd", "divergence_control", fallback="none").split(";")[0].strip() == "glm"
    q = {1: DENSITY, 2: MOMENTUM, 3: MOMENTUM, 4: MOMENTUM, 5: PRESSURE, 6: FIELD, 7: FIELD, 8: FIELD,
         9: PSI_GLM if glm else FIELD}
    aux = {1: DENSITY, 2: VELOCITY, 3: VELOCITY, 4: VELOCITY, 5: PRESSURE, 7: VELOCITY2, 8: VELOCITY, 9: FIELD,
           10: FIELD, 11: FIELD}
    dim = (q if field == 1 else aux if field == 2 else {}).get(var)
    if dim is None:
        raise SystemExit(f"scaling: [{section}] tol of field {field} var {var} is not classified")
    return (dim[0] - 1, dim[1], dim[2])


def classify(ini: configparser.ConfigParser, section: str, key: str) -> tuple[float, float, float]:
    """Dimension of an option for the reference layer; refuses an unclassified one."""
    ic = ini["initial_conditions"] if ini.has_section("initial_conditions") else {}
    table: dict = {}
    if section.startswith("initial_conditions_region_"):
        table = {**STATE_KEYS, **GRID_KEYS}
    elif section.startswith("amr_marker_"):
        table = {**{k: NONE for k in MARKER_NONE}, **{k: LENGTH for k in ("delta_fine", "delta_coarse")}, **BOX_KEYS,
                 "tol": marker_tol(ini, section)}
    elif section.startswith("solid_"):
        table = {**{k: NONE for k in ("name", "definition", "bc_type", "circle_axis", "rectangle_axis")},
                 **{k: LENGTH for k in SOLID_LENGTHS}}
    elif section.startswith("slice_"):
        table = {**{k: NONE for k in ("itype", "n_save", "ni", "nj", "nk")}, **GRID_KEYS}
    elif section in ("reference",) + DIMENSIONLESS_SECTIONS:
        return NONE
    elif section in BC_SECTIONS:
        table = {"type": NONE, **STATE_KEYS}
    elif section == "initial_conditions":
        table = dict(IC_KEYS)
        ic_type = ic.get("type", "").split(";")[0].strip()
        wave = ic.get("wave", "").split(";")[0].strip()
        table["wave_amplitude"] = FIELD if ic_type == "mhd-cpaw" else MOMENTUM if wave == "alfven" else DENSITY
    elif section in FULL_KEYS:
        table = FULL_KEYS[section]
    if key not in table:
        raise SystemExit(f"scaling: [{section}] {key} is not classified for the reference layer")
    return table[key]


def dimensionalize(args: argparse.Namespace) -> int:
    """Write the input in the units L0 = 2^j, u0 = 2^k, rho0 = 4^m, with the [reference] section converting it back."""
    ini = rescale_ini(args)
    if ini is None:
        print(f"scaling: {args.ini} is refused by the reference layer (forest manifest or orszag-tang); skipped")
        return 3
    ini["reference"] = {"density": repr(4.0**args.m), "length": repr(2.0**args.j), "velocity": repr(2.0**args.k)}
    with open(args.out, "w") as f:
        ini.write(f)
    return 0

def check_log(args: argparse.Namespace) -> int:
    """Check the logged conversions of a reference run against the base input; 0 if every one is exact and complete."""
    ini = read_ini(args.ini)
    logged: dict[tuple[str, str], str] = {}
    gamma = None
    for line in args.log.read_text(errors="replace").splitlines():
        if not line.startswith("[mpi-00000]  ["):
            continue
        head, _, tail = line[len("[mpi-00000]  ["):].partition("] ")
        if head == "physics" and tail.startswith("cp, cv:"):
            gamma = float(tail.split("-> gamma")[1].split(",")[0])
            continue
        key, sep, values = tail.partition(": ")
        if sep and " -> " in values:
            logged[(head, key)] = values.split(" -> ")[1].strip()
    status = 0
    for section in ini.sections():
        for key in ini[section]:
            text = ini[section][key].split(";")[0].strip()
            if section == "physics" and key in ("cp", "cv"):
                continue
            if classify(ini, section, key) == NONE or (key == "glm_damping_length" and text == "min-cell"):
                continue
            if (section, key) not in logged:
                print(f"   [{section}] {key}: not converted")
                status = 1
            elif float(logged[(section, key)]) != float(text):
                print(f"   [{section}] {key}: converted to {logged[(section, key)]}, base {text}")
                status = 1
    if ini.has_option("physics", "cp"):
        expected = float(ini["physics"]["cp"].split(";")[0]) / float(ini["physics"]["cv"].split(";")[0])
        if gamma != expected:
            print(f"   [physics] cp, cv: gamma {gamma}, expected cp/cv = {expected!r}")
            status = 1
    print(f"   conversions: {len(logged)} logged, {'exact and complete' if status == 0 else 'MISMATCH'}")
    return status


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
        if args.conservative and name not in CONSERVATIVE:
            continue
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
    c.add_argument("--conservative", action="store_true", help="compare the conservative fields only")
    g = sub.add_parser("check-log")
    g.add_argument("ini", type=Path)
    g.add_argument("log", type=Path)
    d = sub.add_parser("dimensionalize")
    d.add_argument("ini", type=Path)
    d.add_argument("out", type=Path)
    for p in (r, c, d):
        p.add_argument("--j", type=int, default=0, help="lengths times 2^j")
        p.add_argument("--k", type=int, default=0, help="velocities times 2^k")
        p.add_argument("--m", type=int, default=0, help="density times 4^m")
    args = parser.parse_args()
    commands = {"rescale": rescale, "compare": compare, "dimensionalize": dimensionalize, "check-log": check_log}
    return commands[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
