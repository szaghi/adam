#!/usr/bin/env python3
"""Write and judge the input-contract cases of the FLUME dissipative terms (issue #65, P1, P2).

Each case is a base input (an Euler or an MHD regression input) with some keys set; `check.sh` runs it and checks:

- `refuse`: the run stops with the expected message (every invalid combination of D-M4-5, D-M4-6, D-M4-7);
- `log`: the run reads the coefficients and logs them: the coefficient implied by each number (`mu = 1/Re`,
  `k = mu cp/Pr`, `eta = 1/Rm`) must be logged exactly (`--check-log`); the run then takes ten steps and logs its
  diffusive time step limit (Euler since P2, MHD since P3);
- `ideal`: zero coefficients and `dissipative_order = 2` change nothing: the run must equal the base run bit for bit;
- `convert`: the case dimensionalised by `../scaling/scaling.py` (a `[reference]` section, powers of two) must log
  every dimensional key converted back exactly (`scaling.py check-log`: the coefficients, the temperatures and the wall
  velocity), then takes ten steps; `convert-refuse`: the dimensionalised case must be refused (a `lundquist` number
  without the Alfvenic preset).

Usage:
    make_contract.py --list                                  # "name kind base expected" per case
    make_contract.py <base.ini> <out.ini> <name>             # write one case
    make_contract.py --check-log <log> <base.ini> <name>     # judge the logged coefficients of a `log` case
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

# the cases that run through stop after a few steps: a viscous sod-x is diffusion-limited, and the `reference` wall,
# 330 times hotter than the gas, makes it so by 2000 times (issue #65, P2)
SHORT = {"it_max": "10"}

# name: (kind, base (euler|mhd), keys {section: {key: value}}, expected message for `refuse`)
CASES: dict[str, tuple[str, str, dict[str, dict[str, str]], str]] = {
    "both-viscosity": ("refuse", "euler", {"physics": {"viscosity": "0.01", "reynolds": "100.0"}},
                       "gives the same term twice (viscosity reynolds)"),
    "both-conductivity": ("refuse", "euler", {"physics": {"viscosity": "0.01", "conductivity": "0.01",
                                                          "prandtl": "0.72"}},
                          "gives the same term twice (conductivity prandtl)"),
    "two-resistivity-numbers": ("refuse", "mhd", {"physics": {"magnetic_reynolds": "100.0", "lundquist": "50.0"}},
                                "gives the same term twice (magnetic_reynolds lundquist)"),
    "prandtl-alone": ("refuse", "euler", {"physics": {"prandtl": "0.72"}}, "(prandtl) needs a viscosity"),
    "resistivity-euler": ("refuse", "euler", {"physics": {"resistivity": "0.01"}},
                          "needs [physics].(physical_model)=mhd-ideal"),
    "negative-viscosity": ("refuse", "euler", {"physics": {"viscosity": "-0.01"}}, "(viscosity) must not be negative"),
    "zero-reynolds": ("refuse", "euler", {"physics": {"reynolds": "0.0"}}, "(reynolds) must be positive"),
    "power-law-inviscid": ("refuse", "euler", {"physics": {"viscosity_law": "power-law"}},
                           "(viscosity_law)=power-law needs a viscosity"),
    "power-law-keys-alone": ("refuse", "euler", {"physics": {"viscosity": "0.01", "viscosity_exponent": "0.7"}},
                             "(viscosity_exponent, reference_temperature) need [physics].(viscosity_law)=power-law"),
    "power-law-no-temperature": ("refuse", "euler", {"physics": {"viscosity": "0.01", "viscosity_law": "power-law",
                                                                 "viscosity_exponent": "0.7"}},
                                 "failed to load [physics].(reference_temperature)"),
    "negative-exponent": ("refuse", "euler", {"physics": {"viscosity": "0.01", "viscosity_law": "power-law",
                                                          "viscosity_exponent": "-0.5", "reference_temperature": "1.0"}},
                          "(viscosity_exponent) must not be negative"),
    "unknown-law": ("refuse", "euler", {"physics": {"viscosity": "0.01", "viscosity_law": "sutherland"}},
                    'unknown [physics].(viscosity_law) "sutherland"'),
    "order-3": ("refuse", "euler", {"numerics": {"dissipative_order": "3"}}, "(dissipative_order) must be 2 or 4"),
    "limiter": ("refuse", "euler", {"physics": {"viscosity": "0.01"}, "numerics": {"positivity_limiter": "cell"}},
                "(positivity_limiter)=cell is refused with dissipative terms"),
    "wall-normal-velocity": ("refuse", "euler", {"bc_x_max": {"type": "wall-noslip", "wall_u": "0.1"}},
                             "(wall_u) is the normal component of the wall velocity and must be 0"),
    "isothermal-no-temperature": ("refuse", "euler", {"bc_x_max": {"type": "wall-isothermal"}},
                                  "failed to load [bc_x_max].(wall_temperature)"),
    "isothermal-negative": ("refuse", "euler", {"bc_x_max": {"type": "wall-isothermal", "wall_temperature": "-1.0"}},
                            "(wall_temperature) must be positive"),
    "order4-ngc2": ("refuse", "euler", {"physics": {"viscosity": "0.01"}, "grid": {"ngc": "2"},
                                        "weno": {"scheme": "weno-u-3"}},
                    "(dissipative_order)=4 needs [grid].(ngc) >= 3"),
    "reynolds": ("log", "euler", {"physics": {"reynolds": "64.0", "prandtl": "0.75"}, "time": SHORT}, ""),
    "coefficients": ("log", "euler", {"physics": {"viscosity": "0.015625", "conductivity": "0.03125"},
                                      "time": SHORT}, ""),
    "magnetic-reynolds": ("log", "mhd", {"physics": {"magnetic_reynolds": "128.0", "reynolds": "32.0"},
                                         "time": SHORT}, ""),
    "zero": ("ideal", "euler", {"physics": {"viscosity": "0.0", "conductivity": "0.0"},
                                "numerics": {"dissipative_order": "2"}}, ""),
    "reference": ("convert", "euler", {"physics": {"viscosity": "0.015625", "conductivity": "0.03125",
                                                   "viscosity_law": "power-law", "viscosity_exponent": "0.7",
                                                   "reference_temperature": "1.5"},
                                       "bc_x_max": {"type": "wall-isothermal", "wall_temperature": "0.9",
                                                    "wall_v": "0.1"}, "time": SHORT}, ""),
    "reference-mhd": ("convert", "mhd", {"physics": {"resistivity": "0.0078125", "viscosity": "0.03125"},
                                         "time": SHORT}, ""),
    "lundquist-reference": ("convert-refuse", "mhd", {"physics": {"lundquist": "50.0"}},
                            "(lundquist) needs [reference].(velocity)=alfvenic"),
}


def read(path: Path) -> configparser.ConfigParser:
    """Return an INI file, keys case-sensitive."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), strict=False, interpolation=None)
    ini.optionxform = str  # type: ignore[assignment,method-assign]
    ini.read(path)
    return ini


def write_case(base: Path, out: Path, name: str) -> None:
    """Write case `name` of the contract from the base input."""
    ini = read(base)
    for section, keys in CASES[name][2].items():
        if not ini.has_section(section):
            ini.add_section(section)
        for key, value in keys.items():
            ini[section][key] = value
    with open(out, "w") as f:
        ini.write(f)


def check_log(log: Path, base: Path, name: str) -> int:
    """Return 0 if the coefficients logged by a `log` case are exactly those its keys imply (mu, k, eta)."""
    physics = read(base)["physics"]
    keys = CASES[name][2]["physics"]
    g = float(physics["gamma"]) if "gamma" in physics else 0.0
    cp = float(physics["cp"]) if "cp" in physics else g * (1.0 / (g - 1.0))  # as the physics object computes it
    mu = float(keys["viscosity"]) if "viscosity" in keys else 1.0 / float(keys["reynolds"]) if "reynolds" in keys \
        else 0.0
    k = float(keys["conductivity"]) if "conductivity" in keys else mu * cp / float(keys["prandtl"]) \
        if "prandtl" in keys else 0.0
    eta = float(keys["resistivity"]) if "resistivity" in keys else 1.0 / float(keys["magnetic_reynolds"]) \
        if "magnetic_reynolds" in keys else 0.0
    logged = {}
    for line in log.read_text(errors="replace").splitlines():
        for label in ("viscosity mu:", "conductivity k:", "resistivity eta:"):
            if line.startswith("[mpi-00000]") and label in line:
                logged[label] = float(line.split(label)[1].split()[0])
    status = 0
    for label, expected in (("viscosity mu:", mu), ("conductivity k:", k), ("resistivity eta:", eta)):
        if label not in logged and expected == 0.0:
            continue  # Euler logs no resistivity
        ok = logged.get(label) == expected
        status |= int(not ok)
        print(f"   {label:18s} logged {logged.get(label)!r}, expected {expected!r}  {'PASS' if ok else 'FAIL'}")
    return status


def main() -> int:
    """List the cases, write one, or judge a logged one."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", help="<base.ini> <out.ini> <name>, or with --check-log <log> <base> <name>")
    parser.add_argument("--list", action="store_true", help="list the cases")
    parser.add_argument("--check-log", action="store_true", help="judge the logged coefficients of a log case")
    args = parser.parse_args()
    if args.list:
        for name, (kind, base, _, expected) in CASES.items():
            print(f"{name}\t{kind}\t{base}\t{expected}")
        return 0
    if len(args.paths) != 3 or args.paths[2] not in CASES:
        parser.error("give three arguments, the last a case name (see --list)")
    if args.check_log:
        return check_log(Path(args.paths[0]), Path(args.paths[1]), args.paths[2])
    write_case(Path(args.paths[0]), Path(args.paths[1]), args.paths[2])
    return 0


if __name__ == "__main__":
    sys.exit(main())
