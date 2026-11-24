#!/usr/bin/env python3
"""Write the cases of the FLUME speedup matrix from matrix.conf (issue #133 phase 3, stdlib only).

Each case gets a directory `<out>/<case>/` (the case names are in matrix.conf, e.g. `1a-step-uniform`) holding:

    the inputs     leg 1: the forest manifest step.ini and the realm inputs step-{A,B,C}.ini (make_step.py --layout
                   aligned); leg 2: input.ini (make_sphere.py); leg 3: input.ini (make_channel.py)
    case.conf      what run.sh needs, as shell assignments: the role, the main input and the realm inputs, the restart
                   basenames, the cells per block of each realm, the tree ratio, the protocol (t_w, settle steps, window
                   lengths, repetitions), the output keys a window switches off, the dump group (DUMPS: `leg<N>` for
                   the cases of a leg, the case name for a weak case), the kind (matrix or weak) and the scaling series
                   of a uniform case (SERIES, SERIES_LABEL, LOAD_TOLERANCE: read by analyze.py)
    README.md      the case: flow, grid, levels, marker, protocol, how run.sh runs it

The inputs are the generators' own, with the run-independent changes made here: the cold-start grid (uniform at the
finest level, or the base level with the runtime marker), the marker and its tolerances, and the frozen twin
(`[amr] frequency = 0`). run.sh makes the per-run changes (restart, step count, outputs) on copies, with ini_set.py.

`--reduced` merges the `[leg.N.reduced]` sections over `[leg.N]` and `[protocol]`: the WSL plumbing test, seconds per
run.

Weak scaling (phase 9): a `[weak.<case>]` section is one more uniform case of the leg named by its `leg` key, at another
size; its keys are those of `[protocol]`, `[leg.N]`, `[weak]` and `[weak.<case>]`, then with `--reduced` those of
`[leg.N.reduced]`, `[weak.reduced]` and `[weak.<case>.reduced]`. It has no runtime or frozen twin, warms up from a cold
start (`warmup_steps`) into its own dump, and its name starts with `w` (run.sh's leg prefixes never select it).

Runs under `python3 -I` (no third-party package: medooza's python3 has none); the generators of legs 1 and 3 are
called as subprocesses from the verification tree.

Usage:
    make_cases.py <out-dir> [--conf matrix.conf] [--reduced] [--legs 1,2,3]
"""

from __future__ import annotations

import argparse
import configparser
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VERIF = HERE.parent / "verification"
MAKE_STEP = VERIF / "step" / "make_step.py"
MAKE_CHANNEL = VERIF / "turbulence" / "make_channel.py"
MAKE_SPHERE = HERE / "make_sphere.py"
ROLES = ("uniform", "runtime", "frozen")
SAVE_SECTIONS = ("IO", "diagnostics")  # the sections whose `*_save` keys a window run sets to 0


class CaseError(Exception):
    """A matrix.conf that cannot be turned into cases."""


def read_conf(path: Path, reduced: bool) -> tuple[dict[int, dict[str, str]], dict[str, dict[str, str]]]:
    """Return the legs and the weak cases of matrix.conf, the reduced overrides merged when asked.

    Parameters
    ----------
    path : Path
        matrix.conf.
    reduced : bool
        Merge the `.reduced` sections (over the full keys: a reduced key always wins).

    Returns
    -------
    tuple[dict[int, dict[str, str]], dict[str, dict[str, str]]]
        Per leg number its keys (the protocol keys included); per weak case name its keys (`[protocol]`, its leg,
        `[weak]`, its section, then the reduced sections of the three), with `uniform` set to the case name and no
        `runtime` or `frozen` key.
    """
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    parser.optionxform = str
    if not parser.read(path):
        raise CaseError(f"no matrix file '{path}'")

    def keys(section: str) -> dict[str, str]:
        """Return a section's keys (empty when absent), each value on one line (continuations joined)."""
        if not parser.has_section(section):
            return {}
        return {k: " ".join(v.split()) for k, v in parser[section].items()}

    protocol = keys("protocol")
    legs: dict[int, dict[str, str]] = {}
    for name in parser.sections():
        m = re.fullmatch(r"leg\.(\d+)", name)
        if not m:
            continue
        leg = protocol | keys(name)
        if reduced:
            leg |= keys(f"{name}.reduced")
        leg["reduced"] = "1" if reduced else "0"
        legs[int(m.group(1))] = leg
    weak: dict[str, dict[str, str]] = {}
    for name in parser.sections():
        m = re.fullmatch(r"weak\.([^.]+)", name)
        if not m or m.group(1) == "reduced":
            continue
        case = m.group(1)
        if not case.startswith("w"):
            raise CaseError(f"[{name}]: a weak case name starts with 'w' (run.sh's leg prefixes must not select it)")
        own = keys(name)
        if "leg" not in own or not own["leg"].isdigit() or f"leg.{own['leg']}" not in parser:
            raise CaseError(f"[{name}]: `leg` must name a [leg.N] section")
        base = f"leg.{own['leg']}"
        merged = protocol | keys(base) | keys("weak") | own
        if reduced:
            merged |= keys(f"{base}.reduced") | keys("weak.reduced") | keys(f"{name}.reduced")
        merged = {k: v for k, v in merged.items() if k not in ("runtime", "frozen")}
        merged |= {"uniform": case, "reduced": "1" if reduced else "0"}
        weak[case] = merged
    return legs, weak


def set_keys(text: str, section: str, keys: dict[str, str]) -> str:
    """Set keys of one INI section, keeping comments; a missing key is added at the section's top.

    Parameters
    ----------
    text : str
        The INI text.
    section : str
        The section name (it must exist).
    keys : dict[str, str]
        The keys and their values.

    Returns
    -------
    str
        The new INI text.
    """
    m = re.search(rf"^\[{re.escape(section)}\][^\n]*\n(.*?)(?=^\[|\Z)", text, flags=re.M | re.S)
    if not m:
        raise CaseError(f"no section [{section}]")
    body = m.group(1)
    for key, value in keys.items():
        body, n = re.subn(rf"^{re.escape(key)}\s*=[^\n]*$", f"{key} = {value}", body, flags=re.M)
        if n == 0:
            body = f"{key} = {value}\n" + body
    return text[: m.start(1)] + body + text[m.end(1):]


def drop_keys(text: str, section: str, keys: tuple[str, ...]) -> str:
    """Remove keys of one INI section (absent keys ignored)."""
    m = re.search(rf"^\[{re.escape(section)}\][^\n]*\n(.*?)(?=^\[|\Z)", text, flags=re.M | re.S)
    if not m:
        raise CaseError(f"no section [{section}]")
    body = re.sub(rf"^(?:{'|'.join(map(re.escape, keys))})\s*=[^\n]*\n", "", m.group(1), flags=re.M)
    return text[: m.start(1)] + body + text[m.end(1):]


def drop_section(text: str, section: str) -> str:
    """Remove one INI section (header and body) when present."""
    return re.sub(rf"^\[{re.escape(section)}\][^\n]*\n.*?(?=^\[|\Z)", "", text, flags=re.M | re.S)


def ini_keys(text: str) -> dict[str, dict[str, str]]:
    """Return the sections and keys of an INI text (inline `;` comments stripped)."""
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    parser.optionxform = str
    parser.read_string(text)
    return {s: dict(parser[s]) for s in parser.sections()}


def run_generator(cmd: list[str]) -> None:
    """Run a generator, stopping with its message when it fails."""
    done = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise CaseError(f"{shlex.join(cmd)} failed:\n{done.stdout}{done.stderr}")


def lohner_or_gradient(leg: dict[str, str], ngc: int, h_base: float, h_fine: float, drho: float) -> dict[str, str]:
    """Return the keys of `[amr_marker_1]` for the density markers of legs 1 and 2 (`marker = lohner | gradient`)."""
    kind = leg["marker"]
    if kind == "lohner":
        buffer = int(leg.get("buffer", "-1"))
        return {"mode": "4", "delta_type": "max", "delta_fine": "0.0", "delta_coarse": "0.0", "field": "1",
                "var": "1", "refine_tol": leg["refine_tol"], "derefine_tol": leg["derefine_tol"], "floor": "0.0",
                "buffer": str(ngc - 1 if buffer < 0 else buffer)}
    if kind == "gradient":
        margin = 1.0 + 1.0e-9
        return {"mode": "2", "delta_type": "max", "delta_fine": repr(h_fine * margin),
                "delta_coarse": repr(h_base * margin), "field": "1", "var": "1",
                "tol": repr(float(leg["tol_factor"]) * drho / h_base)}
    raise CaseError(f"marker '{kind}' is not lohner or gradient")


def marker_section(keys: dict[str, str]) -> str:
    """Return an `[amr_marker_1]` section."""
    return "[amr_marker_1]\n" + "".join(f"{k:12s} = {v}\n" for k, v in keys.items())


def step_case(leg: dict[str, str], role: str, out: Path) -> tuple[str, list[str], list[str]]:
    """Write a leg-1 case (the aligned step forest); return the main input, the realm inputs and the command."""
    weno, m = int(leg["weno"]), int(leg["m"])
    lb, lf = int(leg["base_level"]), int(leg["max_level"])
    ngc = max(3, (weno + 1) // 2)
    cmd = [sys.executable, "-I", str(MAKE_STEP), str(out), "--layout", "aligned", "--dims", "2", "--m", str(m),
           "--weno", str(weno), "--time-max", leg["t_w"], "--cfl", leg["cfl"], "--it-save", leg["huge"],
           "--restart-save", "0", "--residuals-save", "1"]
    if role == "uniform":
        cmd += ["--base-level", str(lf)]
    else:
        cmd += ["--base-level", str(lb), "--max-level", str(lf), "--runtime", "--frequency", leg["regrid_period"]]
    run_generator(cmd)
    inis = [f"step-{r}.ini" for r in ("A", "B", "C")]
    if role != "uniform":
        h_base, h_fine = 0.2 / (m * 2**lb), 0.2 / (m * 2**lf)
        # drho for `marker = gradient`: the jump of a normal shock at the inflow Mach number 3 (rho 1.4 -> 5.4)
        drho = 1.4 * (2.4 * 9.0 / (0.4 * 9.0 + 2.0) - 1.0)
        keys = lohner_or_gradient(leg, ngc, h_base, h_fine, drho=drho)
        for ini in inis:
            path = out / ini
            text = path.read_text()
            text = drop_section(text, "amr_marker_1")
            text = text.replace("[runge_kutta]", marker_section(keys) + "[runge_kutta]", 1)
            if role == "frozen":
                text = set_keys(text, "amr", {"frequency": "0"})
            path.write_text(text)
    return "step.ini", inis, cmd


def sphere_case(leg: dict[str, str], role: str, out: Path) -> tuple[str, list[str], list[str]]:
    """Write a leg-2 case (the sphere); return the main input, the realm inputs and the command."""
    lb, lf = leg["base_level"], leg["max_level"]
    cmd = [sys.executable, "-I", str(MAKE_SPHERE), str(out / "input.ini"), "--block", leg["block"], "--mach",
           leg["mach"], "--shock-x", leg["shock_x"], "--radius", leg["radius"], "--weno", leg["weno"], "--cfl",
           leg["cfl"], "--time-max", leg["t_w"], "--it-save", leg["huge"], "--restart-save", "0",
           "--residuals-save", "1", "--history", "1"]
    if leg["marker"] != "gradient":
        raise CaseError("leg 2: make_sphere.py writes the gradient marker only")
    if role == "uniform":
        cmd += ["--base-level", lf]
    else:
        cmd += ["--base-level", lb, "--max-level", lf, "--runtime", "--frequency", leg["regrid_period"],
                "--tol-factor", leg["tol_factor"]]
    run_generator(cmd)
    if role == "frozen":
        path = out / "input.ini"
        path.write_text(set_keys(path.read_text(), "amr", {"frequency": "0"}))
    return "input.ini", ["input.ini"], cmd


def channel_case(leg: dict[str, str], role: str, out: Path) -> tuple[str, list[str], list[str]]:
    """Write a leg-3 case (the channel); return the main input, the realm inputs and the command."""
    lf = leg["max_level"]
    cmd = [sys.executable, "-I", str(MAKE_CHANNEL), str(out / "input.ini"), "--nx", leg["nx"], "--ny", leg["ny"],
           "--nz", leg["nz"], "--levels", lf, "--mach", leg["mach"], "--reynolds", leg["reynolds"], "--power-law",
           leg["power_law"], "--scheme-space", leg["scheme_space"], "--perturbation", leg["perturbation"], "--cfl",
           leg["cfl"], "--time-max", leg["t_w"], "--it-save", leg["huge"], "--restart-save", "0", "--history", "1"]
    run_generator(cmd)
    path = out / "input.ini"
    # make_channel.py keeps the isentropic-vortex keys of its template, which `[IO] unused_keys = error` (the default
    # since #80 P5) refuses: the channel initial condition reads none of them
    path.write_text(drop_keys(path.read_text(), "initial_conditions", ("x0", "y0", "radius", "strength")))
    if role != "uniform":
        if leg["marker"] != "kolmogorov":
            raise CaseError("leg 3: the runtime marker is the Kolmogorov resolution (marker = kolmogorov)")
        text = path.read_text()
        text = set_keys(text, "amr", {"iu_ref_levels": leg["base_level"], "markers_number": "1",
                                      "frequency": "0" if role == "frozen" else leg["regrid_period"]})
        text = set_keys(text, "numerics", {"reflux": ".true."})
        text = text.rstrip("\n") + "\n\n" + marker_section(
            {"mode": "7             ; local Kolmogorov resolution r = h/eta (AMR_KOLMOGOROV)", "delta_type": "max",
             "delta_fine": "0.0", "delta_coarse": "0.0", "refine_tol": leg["refine_tol"],
             "derefine_tol": leg["derefine_tol"]})
        path.write_text(text)
    return "input.ini", ["input.ini"], cmd


GENERATORS = {"step": step_case, "sphere": sphere_case, "channel": channel_case}


def describe(out: Path, inis: list[str]) -> tuple[list[str], list[int], int, list[str], dict[str, str]]:
    """Read back the realm inputs: restart basenames, cells per block, tree ratio, the `*_save` keys, the [amr] keys."""
    basenames, cells, off = [], [], []
    ratio, amr = 0, {}
    for ini in inis:
        keys = ini_keys((out / ini).read_text())
        basenames.append(keys["IO"]["restart_basename"])
        g = keys["grid"]
        cells.append(int(g["ni"]) * int(g["nj"]) * int(g["nk"]))
        ratio, amr = int(keys["amr"]["ratio"]), keys["amr"]
        for section in SAVE_SECTIONS:
            off += [f"{section}.{k}=0" for k in keys.get(section, {}) if k.endswith("_save")]
    return basenames, cells, ratio, sorted(set(off)), amr


def readme(leg_no: int, leg: dict[str, str], role: str, case: str, conf: dict[str, str], cmd: list[str],
           amr: dict[str, str], out: Path) -> str:
    """Return the README.md of one case."""
    cells = [int(c) for c in conf["BLOCK_CELLS"].split()]
    ratio, lf = int(conf["RATIO"]), int(leg["max_level"])
    uniform_cells = sum(c * ratio**lf for c in cells)
    what = {
        "uniform": f"the uniform grid at the finest level L{lf}: {uniform_cells:,} cells. The warm-up runs it from "
                   f"the cold start to t_w = {leg['t_w']}; its windows restart from the warm-up dump.",
        "runtime": f"runtime AMR from the base level L{leg['base_level']} up to L{lf}, regrid every "
                   f"{leg['regrid_period']} steps ({leg['marker']} marker). The settle restarts it from the warm-up "
                   f"dump for {leg['settle_steps']} steps; its windows restart from the settle dump.",
        "frozen": "the grid of the runtime case frozen ([amr] frequency = 0; D2): its windows restart from the settle "
                  "dump, the same grid and state as the runtime windows.",
    }[role]
    reduced = " **Reduced size** (`make_cases.py --reduced`: plumbing only, timings mean nothing)." \
        if leg["reduced"] == "1" else ""
    if conf["KIND"] == "weak":
        what = (f"the uniform grid at L{lf}: {uniform_cells:,} cells, a weak-scaling case (phase 9). The warm-up runs "
                f"it from the cold start for {leg['warmup_steps']} steps into its own dump; its windows restart from "
                "it. No runtime or frozen twin.")
    series = (f"Scaling series `{conf['SERIES']}`: {conf['SERIES_LABEL']}." if conf["SERIES"] else
              "No scaling series (strong scaling only).")
    lines = [
        f"# {case}", "",
        f"FLUME speedup matrix, leg {leg_no} ([issue #133](https://github.com/szaghi/adam/issues/133)): "
        f"{leg['title']}.{reduced}", "",
        f"**Role: {role}**, {what}", "", series, "",
        "| Item | Value |", "|---|---|",
        f"| main input | `{conf['INPUT']}` (realm inputs: {', '.join(f'`{i}`' for i in conf['INIS'].split())}) |",
        f"| cells per block (per realm) | {', '.join(conf['BLOCK_CELLS'].split())} |",
        f"| tree ratio | {ratio} |",
        f"| levels | base L{amr.get('iu_ref_levels', '?')}, finest L{amr.get('max_level', '?')}, "
        f"[amr] frequency {amr.get('frequency', '0')} |",
        f"| cells at the finest level, uniform | {uniform_cells:,} |",
        f"| warm-up | {'t = ' + leg['t_w'] if int(leg.get('warmup_steps', '0')) <= 0 else leg['warmup_steps'] + ' steps'}"
        " (the leg's uniform case) |",
        f"| settle | {leg['settle_steps']} steps (the leg's runtime case) |" if conf["KIND"] != "weak"
        else "| settle | none (a weak case) |",
        f"| window | {leg['repetitions']} x ({leg['window_short']} + {leg['window_long']} steps), interleaved |",
        "", "Generated by `make_cases.py` from `matrix.conf` (edit there, not here):", "", "```",
        shlex.join([Path(c).name if c.endswith(".py") else c.replace(f"{out}/", "").replace(str(out), ".")
                    for c in cmd[2:]]), "```", "",
        "`run.sh` runs it (`../../README.md`, \"How to run\"); the results go in `results/summary.md`.", "",
    ]
    return "\n".join(lines)


def write_case(leg_no: int, leg: dict[str, str], role: str, root: Path, kind: str = "matrix") -> str:
    """Write one case directory; return its name.

    Parameters
    ----------
    leg_no : int
        The leg (flow) of the case.
    leg : dict[str, str]
        The keys of the leg, or of the weak case.
    role : str
        uniform, runtime or frozen.
    root : Path
        The directory of the case directories.
    kind : str
        matrix (a case of a leg) or weak (a `[weak.<case>]` section).

    Returns
    -------
    str
        The case name.
    """
    case = leg[role]
    out = root / case
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("*.ini"):
        old.unlink()
    generator = GENERATORS.get(leg["generator"])
    if generator is None:
        raise CaseError(f"leg {leg_no}: unknown generator '{leg['generator']}'")
    main_input, inis, cmd = generator(leg, role, out)
    basenames, cells, ratio, off, amr = describe(out, inis)
    for key in ("window_short", "window_long"):
        if int(leg[key]) % int(leg["regrid_period"]):
            raise CaseError(f"leg {leg_no}: {key} {leg[key]} is not a multiple of regrid_period {leg['regrid_period']}")
    conf = {
        "CASE": case, "LEG": str(leg_no), "ROLE": role, "INPUT": main_input, "INIS": " ".join(inis),
        "BASENAMES": " ".join(basenames), "BLOCK_CELLS": " ".join(map(str, cells)), "RATIO": str(ratio),
        "UNIFORM_CASE": leg["uniform"], "RUNTIME_CASE": leg.get("runtime", ""), "FROZEN_CASE": leg.get("frozen", ""),
        "SOURCE": "warmup" if role == "uniform" else "settle",
        "RULE": "units" if role == "runtime" else "siblings",
        "T_W": leg["t_w"], "WARMUP_STEPS": leg.get("warmup_steps", "0"), "SETTLE_STEPS": leg["settle_steps"],
        "REGRID_PERIOD": leg["regrid_period"], "WINDOW_SHORT": leg["window_short"], "WINDOW_LONG": leg["window_long"],
        "REPETITIONS": leg["repetitions"], "SPREAD_WARNING": leg["spread_warning"], "HUGE": leg["huge"],
        "OFF_KEYS": " ".join(off), "REDUCED": leg["reduced"],
        "KIND": kind, "DUMPS": case if kind == "weak" else f"leg{leg_no}",
        "SERIES": leg.get("series", "") if role == "uniform" else "",
        "SERIES_LABEL": leg.get("series_label", "") if role == "uniform" else "",
        "LOAD_TOLERANCE": leg.get("load_tolerance", "0.10"),
    }
    with (out / "case.conf").open("w") as f:
        f.write(f"# {case}: written by make_cases.py from matrix.conf (sourced by run.sh; do not edit)\n")
        f.writelines(f"{k}={shlex.quote(v)}\n" for k, v in conf.items())
    (out / "README.md").write_text(readme(leg_no, leg, role, case, conf, cmd, amr, out))
    return case


def main() -> int:
    """Parse the command line and write the cases."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out", type=Path, help="the directory of the case directories")
    parser.add_argument("--conf", type=Path, default=HERE / "matrix.conf", help="the matrix (default: matrix.conf)")
    parser.add_argument("--reduced", action="store_true", help="the reduced sizes (WSL plumbing test)")
    parser.add_argument("--legs", default="all",
                        help="comma-separated leg numbers, or all (default); a leg brings its weak cases")
    args = parser.parse_args()
    try:
        legs, weak = read_conf(args.conf, args.reduced)
        wanted = sorted(legs) if args.legs == "all" else [int(x) for x in args.legs.split(",")]
        for leg_no in wanted:
            if leg_no not in legs:
                raise CaseError(f"no [leg.{leg_no}] in {args.conf}")
            for role in ROLES:
                print(f"make_cases: {write_case(leg_no, legs[leg_no], role, args.out)}")
        for case, keys in weak.items():
            if int(keys["leg"]) in wanted:
                print(f"make_cases: {write_case(int(keys['leg']), keys, 'uniform', args.out, kind='weak')}")
    except (CaseError, KeyError, ValueError) as err:
        print(f"make_cases: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
