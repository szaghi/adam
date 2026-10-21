#!/usr/bin/env python3
"""Write the inputs and READMEs of the FLUME speedup matrix, and the small helpers its scripts need (stdlib only).

Why this generator exists: the matrix is sixteen cases built from two flows (a Mach 2 shock over a cylinder in 2-D,
over a sphere in 3-D), each without AMR (a uniform grid, a fixed-grid reference) and with dynamic AMR (a
gradient-of-density marker): the cylinder in three orientations (xy, xz, yz) on octrees plus xy on a quadtree, a
heavier quadtree, and the AMD-only sizes of the quadtree and the sphere (FLOWS; each flow reads its grid from
matrix.conf by its own prefix). Writing them by hand would let the rotated copies drift apart; here one description
is rotated and resized, so the cases differ only where they must.

The flow is the V6 verification (verification/shock-cylinder/shock-cylinder.ini): pre-shock rho 1.4, p 1 (sound speed
1), post-shock rho 3.7333, u_n 1.25, p 4.5 upstream of a = 0.2 and at the inflow (a the streamwise axis), a body of
radius 0.13 at the domain centre (Euler wall, immersed boundary), inviscid walls on the cross-stream faces, outflow by
extrapolation; WENO-5 characteristic flux splitting, SSP-33, CFL 0.5, reflux on.

The grid parameters come from matrix.conf; until they are numbers the inputs carry TBD and run.sh refuses them.

Stdlib only: medooza's python3 has no numpy/h5py, and the scripts run there.

Usage:
    make_inputs.py inputs                        write <case>/input.ini and <case>/README.md of every case
    make_inputs.py readmes                       rewrite the READMEs only (build.sh and analyze.py call it)
    make_inputs.py cases                         print the case names, one per line
    make_inputs.py set <ini> <sec.key=value>...  rewrite options of an input in place
    make_inputs.py restart-step <file.time>      print the step and time of a restart (stream: int32, float64)
"""

from __future__ import annotations

import configparser
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TBD = "TBD"
AXES = "xyz"
PRE = {"r": 1.4, "p": 1.0}
POST = {"r": 3.7333333333333334, "un": 1.25, "p": 4.5}
DRHO = POST["r"] - PRE["r"]
SHOCK_AT = 0.2
RADIUS = 0.13
LEGS = ("cpu", "fnl-nvf", "fnl-amd")

# name: (body, axes, tree, grid); axes = streamwise, cross-stream, null (the sphere: x, y, z, all active); grid = the
# prefix of the flow's grid keys in matrix.conf (<grid>_NI, _NOAMR_LEVEL, _BASE_LEVEL, _MAX_LEVEL)
FLOWS = {
    "xy-quad-shock-cylinder": ("cylinder", "xyz", "quad", "CYL"),
    "xy-oct-shock-cylinder": ("cylinder", "xyz", "oct", "CYL"),
    "xz-oct-shock-cylinder": ("cylinder", "xzy", "oct", "CYL"),
    "yz-oct-shock-cylinder": ("cylinder", "yzx", "oct", "CYL"),
    "xyz-shock-sphere": ("sphere", "xyz", "oct", "SPH"),
    "xy-quad-heavy-shock-cylinder": ("cylinder", "xyz", "quad", "QH"),
    "xy-quad-amd-shock-cylinder": ("cylinder", "xyz", "quad", "AQ"),
    "xyz-amd-shock-sphere": ("sphere", "xyz", "oct", "AS"),
}
CASES = [f"{flow}-{amr}" for flow in FLOWS for amr in ("noamr", "amr")]


def load_conf() -> dict[str, str]:
    """Return the KEY=VALUE pairs of matrix.conf (comments and quotes stripped)."""
    conf = {}
    for line in (ROOT / "matrix.conf").read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            key, value = line.split("=", 1)
            conf[key.strip()] = value.strip().strip('"')
    return conf


def grid_of(p: str, conf: dict[str, str]) -> tuple[str, str, str, str]:
    """Return NI, BASE_LEVEL, MAX_LEVEL, NOAMR_LEVEL of grid prefix p (strings, possibly TBD), checked when numbers."""
    ni, base, top, uniform = (conf[f"{p}_{k}"] for k in ("NI", "BASE_LEVEL", "MAX_LEVEL", "NOAMR_LEVEL"))
    if ni != TBD and (int(ni) < 6 or int(ni) % 2):
        raise SystemExit(f"make_inputs: {p}_NI = {ni} must be even and >= 6 (2 ngc along a refined axis)")
    if TBD not in (base, top) and not 0 <= int(base) < int(top):
        raise SystemExit(f"make_inputs: need 0 <= {p}_BASE_LEVEL < {p}_MAX_LEVEL, got {base}, {top}")
    if uniform != TBD and int(uniform) < 0:
        raise SystemExit(f"make_inputs: {p}_NOAMR_LEVEL = {uniform} must be >= 0")
    return ni, base, top, uniform


def uniform_cells(body: str, tree: str, ni: int, level: int, conf: dict[str, str]) -> tuple[int, int]:
    """Return the blocks and cells of a uniform grid at a level (the octree splits its null axis too)."""
    if body == "sphere":
        return 8**level, (ni * 2**level) ** 3
    if tree == "quad":
        return 4**level, (ni * 2**level) ** 2
    return 8**level, (ni * 2**level) ** 2 * int(conf["OCT_NULL_CELLS"]) * 2**level


def spacing(ni: str, level: str, pad: float = 1.0) -> str:
    """Return the cell size 1 / (ni 2^level) (times pad), or TBD."""
    if TBD in (ni, level):
        return TBD
    return repr(pad / (int(ni) * 2 ** int(level)))


def state_lines(a: str, post: bool, pad: str) -> list[str]:
    """Return the r, u, v, w, p lines of the post-shock or pre-shock state, the velocity along axis a."""
    r, un, p = (POST["r"], POST["un"], POST["p"]) if post else (PRE["r"], 0.0, PRE["p"])
    vel = {ax: (un if ax == a else 0.0) for ax in AXES}
    return [f"r{pad} = {r!r}", f"u{pad} = {vel['x']!r}", f"v{pad} = {vel['y']!r}", f"w{pad} = {vel['z']!r}",
            f"p{pad} = {p!r}"]


def case_ini(case: str, conf: dict[str, str]) -> str:
    """Return the input of one case: the full-simulation configuration (run.sh rewrites IO and time per run)."""
    flow, amr = case.rsplit("-", 1)
    body, (a, b, n), tree, grid = FLOWS[flow]
    dynamic = amr == "amr"
    ni, base, top, uniform = grid_of(grid, conf)
    null = n if body == "cylinder" else ""
    cells = {ax: ("1" if tree == "quad" else conf["OCT_NULL_CELLS"]) if ax == null else ni for ax in AXES}
    walls = {b} if body == "cylinder" else {b, n}
    basename = f"shock-{body}"
    sweeps = TBD if TBD in (base, top) else str(int(top) - int(base))
    tol = TBD if TBD in (ni, base) else repr(float(conf["GRAD_TOL_FACTOR"]) * DRHO / float(spacing(ni, base)))
    what = (f"plane {a}{b}, null axis {n}, {tree}tree" if body == "cylinder" else "3-D octree") + (
        f"; dynamic AMR, levels {base}..{top}, gradient-of-density marker" if dynamic else f"; uniform at level {uniform}")
    lines = [
        f"; FLUME speedup matrix, case {case} (make_inputs.py: edit matrix.conf and rerun, not this file)",
        f"; Mach 2 shock over a {body}, streamwise axis {a}, {what}",
        "[IO]", f"output_basename        = {basename}", f"it_save                = {conf['FULL_IT_SAVE']}",
        "restart                = .false.", f"restart_basename       = {basename}-restart", "restart_save           = 0",
        "residuals_save         = 1", "save_memory_status     = .false.", "save_residual_fields   = .false.",
        "save_auxiliary_fields  = .false.", "save_curl_fields       = .false.", "save_divergence_fields = .false.",
        "save_gradient_fields   = .false.", "save_laplacian_fields  = .false.",
        "", "[grid]", f"ni     = {cells['x']}", f"nj     = {cells['y']}", f"nk     = {cells['z']}", "ngc    = 3",
        "emin_x = 0.0", "emin_y = 0.0", "emin_z = 0.0", "emax_x = 1.0", "emax_y = 1.0", "emax_z = 1.0",
    ]
    lines += [f"null_{ax} = {'.true.' if ax == null else '.false.'}" for ax in AXES]
    lines += [
        "", "[amr]", f"max_level      = {top if dynamic else uniform}", f"ratio          = {4 if tree == 'quad' else 8}",
        f"iu_ref_levels  = {base if dynamic else uniform}", "i_prune        = 0", "j_prune        = 0",
        "k_prune        = 0", "l_prune        = -1", f"frequency      = {conf['AMR_FREQUENCY'] if dynamic else 0}",
        f"iters          = {sweeps if dynamic else 1}", f"markers_number = {1 if dynamic else 0}",
    ]
    if dynamic:
        # the admissible spacings are padded by 1e-9 so that a block exactly at the finest (base) spacing is left
        # alone (derefined) whatever the rounding of its dxyz
        lines += [
            "", "[amr_marker_1]", "mode         = 2             ; gradient", "delta_type   = max",
            f"delta_fine   = {spacing(ni, top, 1.0 + 1.0e-9)}   ; finest spacing: blocks where |grad rho| > tol",
            f"delta_coarse = {spacing(ni, base, 1.0 + 1.0e-9)}   ; base spacing: every other block derefines to it",
            "field        = 1             ; conservative variables", "var          = 1             ; density",
            f"tol          = {tol}   ; {conf['GRAD_TOL_FACTOR']} * drho / h(base), drho = {DRHO:.4f}",
        ]
    lines += [
        "", "[field]", "nv = 5", "", "[runge_kutta]", "scheme = runge-kutta-ssp-33",
        "", "[weno]", "scheme              = weno-u-5", "ror_number          = 0", "ror_threshold       = 0.9",
        "ror_vars_number     = 0", "enable_ror_stats    = .false.", "ib_reduction_extent = 0", "ib_reduced_order    = 2",
        "", "[linear-algebra]", "smoothing         = gauss-seidel", "iterations_init   = 3", "iterations_coarse = 10",
        "iterations_fine   = 3", "iterations        = 10", "tolerance         = 1.e-20",
        "", "[fdv]", "fdv_scheme = fd", "fdv_order  = 2",
        "", "[solids]", "number    = 1", "n_eikonal = 2",
        "", "[solid_1]", f"name            = {body}", "bc_type         = 2                 ; Euler wall",
    ]
    if body == "cylinder":
        lines += ["definition      = analytical_circle"]
        lines += [f"circle_center_{ax} = 0.5" for ax in AXES]
        lines += [f"circle_radius   = {RADIUS}", f"circle_axis     = {n}"]
    else:
        lines += ["definition      = analytical_sphere"]
        lines += [f"sphere_center_{ax} = 0.5" for ax in AXES]
        lines += [f"sphere_radius   = {RADIUS}"]
    lines += [
        "", "[slices]", "slices_number = 0",
        "", "[physics]", "physical_model = euler", "cp             = 1040.004", "cv             = 742.86",
        "", "[numerics]", "scheme_space             = weno", "reconstruction_variables = characteristic",
        "reflux                   = .true.", "",
    ]
    for ax in AXES:
        for side in ("min", "max"):
            lines.append(f"[bc_{ax}_{side}]")
            if ax == a and side == "min":
                lines += ["type = inflow"] + state_lines(a, post=True, pad="")
            elif ax in walls:
                lines.append("type = wall-inviscid")
            else:
                lines.append("type = extrapolation")
    lines += [
        "", "[initial_conditions]", "type           = riemann-problem", "regions_number = 2",
        f"amr_iterations = {sweeps if dynamic else 0}",
    ]
    for region, post in ((1, True), (2, False)):
        lines += ["", f"[initial_conditions_region_{region}]"] + state_lines(a, post=post, pad="     ")
        for ax in AXES:
            lo = SHOCK_AT if (ax == a and not post) else 0.0
            hi = SHOCK_AT if (ax == a and post) else 1.0
            lines += [f"emin_{ax} = {lo!r}", f"emax_{ax} = {hi!r}"]
    lines += [
        "", "[time]", "it_max   = -1", f"time_max = {conf['T_END']}", f"CFL      = {conf['CFL']}",
        "", "[diagnostics]", "conservation_history_save = 1",
    ]
    return "\n".join(lines) + "\n"


def grid_size(body: str, tree: str, ni: str, level: str, conf: dict[str, str]) -> str:
    """Return the blocks and cells of the uniform grid at a level, or TBD."""
    if TBD in (ni, level):
        return "TBD (matrix.conf)"
    blocks, cells = uniform_cells(body, tree, int(ni), int(level), conf)
    return f"{blocks} blocks, {cells} cells (h = {spacing(ni, level)})"


def case_readme(case: str, conf: dict[str, str]) -> str:
    """Return the README of one case: description, grid, protocol, layout, the toolchain of every built leg."""
    flow, amr = case.rsplit("-", 1)
    body, (a, b, n), tree, grid = FLOWS[flow]
    dynamic = amr == "amr"
    ni, base, top, uniform = grid_of(grid, conf)
    if body == "cylinder":
        geometry = f"2-D in the {a}{b} plane, {n} null; " + (
            "quadtree (`ratio = 4`), one cell along the null axis" if tree == "quad" else
            f"octree (`ratio = 8`), {conf['OCT_NULL_CELLS']} cells per block along the null axis, which the octree "
            "also splits at every level")
        body_line = f"a cylinder of radius {RADIUS} (axis {n}) through the domain centre"
        walls = b
    else:
        geometry = "3-D octree (`ratio = 8`)"
        body_line = f"a sphere of radius {RADIUS} centred in the unit cube"
        walls = f"{b} and {n}"
    if dynamic:
        amr_row = (f"dynamic: base level {base}, finest {top}; one gradient-of-density marker (mode 2, `field = 1`, "
                   f"`var = 1`), `tol = {conf['GRAD_TOL_FACTOR']}·Δρ/h(base)`; a regrid every "
                   f"{conf['AMR_FREQUENCY']} steps, `iters = max - base`, conservative prolongation")
        grid_row = (f"base level: {grid_size(body, tree, ni, base, conf)}; a uniform grid at the finest level would be "
                    f"{grid_size(body, tree, ni, top, conf)}")
    else:
        amr_row = (f"none: uniform at level {uniform}, no marker, no regrid (a fixed-grid reference, independent of "
                   "the amr levels)")
        grid_row = grid_size(body, tree, ni, uniform, conf)
    has_full = case in conf["FULL_CASES"].split()
    if flow in conf["AMD_ONLY_FLOWS"].split():
        legs_row = "fnl-amd only (sized for the 192 GB AMD GPUs; run.sh skips it on the other legs unless named)"
    else:
        legs_row = "every leg (cpu, fnl-nvf, fnl-amd)"
    out = [
        f"# {case}",
        "",
        "Part of the FLUME speedup matrix ([../README.md](../README.md)). Generated by `make_inputs.py` from",
        "`matrix.conf`: edit those, not this file.",
        "",
        "## Case",
        "",
        f"A Mach 2 planar shock over {body_line}: immersed boundary, Euler wall; the flow of the V6",
        f"verification (`verification/shock-cylinder`). The shock starts at {a} = {SHOCK_AT} and travels along +{a}:",
        f"post-shock (ρ, u_{a}, p) = (3.7333, 1.25, 4.5) upstream and at the inflow, pre-shock (1.4, 0, 1), sound",
        f"speed 1. Inflow at {a}-min, outflow (extrapolation) at {a}-max, inviscid walls on the {walls} faces.",
        f"Compressible Euler, ideal gas γ = 1.4, WENO-5 characteristic flux splitting, SSP-33, CFL {conf['CFL']},",
        "reflux on.",
        "",
        "| Item | Value |",
        "|---|---|",
        f"| Geometry | {geometry} |",
        f"| Cells per block | {ni} per active axis |",
        f"| AMR | {amr_row} |",
        f"| Grid | {grid_row} |",
        f"| Legs | {legs_row} |",
        f"| Warmup | to t = {conf['T_WARM']} (incident shock at mid-body), then a restart |",
        f"| Window | {conf['NWIN']} steps after {conf['N1']}, {conf['REPEATS']} repetitions, all output off |",
        f"| Full simulation | {'to t = ' + conf['T_END'] if has_full else 'no (speedup data only)'} |",
        "",
        "## Layout",
        "",
        "```",
        f"{case}/",
        "  input.ini                      the full-simulation input (run.sh rewrites IO and time per run)",
        "  README.md                      this file",
        "  data/<leg>-np<N>/              the measurement (run.sh <leg> --what data)",
        "    warmup/                      warmup to T_WARM; restart files at its end",
        "    rep<r>-short/ rep<r>-long/   restarted runs of N1 and N1 + NWIN steps, all output off",
        "    timings.csv run.info         wall times; host, executable, launcher, restart step",
        "    result.md                    the cost of this leg (analyze.sh)",
    ]
    if has_full:
        out += ["  full/<leg>-np<N>/              the complete simulation to T_END (run.sh <leg> --what full)"]
    out += ["```", "", "## Toolchain and build", ""]
    for leg in LEGS:
        tc = ROOT / "builds" / leg / "toolchain.md"
        out += [f"### {leg}", ""]
        if tc.is_file():
            out += [tc.read_text().rstrip(), "", f"Full build log: [build.txt](../builds/{leg}/build.txt)."]
        else:
            out += [f"Not built yet (`./build.sh {leg}`)."]
        out += [""]
    out += ["## Results", ""]
    data = ROOT / case / "data"
    results = sorted(data.glob("*/result.md")) if data.is_dir() else []
    out += [f"- [{r.parent.name}](data/{r.parent.name}/result.md)" for r in results] or ["None yet."]
    return "\n".join(out) + "\n"


def set_options(ini: Path, items: list[str]) -> None:
    """Rewrite `section.key=value` options of an input in place (comments are not kept)."""
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    parser.optionxform = str
    parser.read(ini)
    for item in items:
        key, value = item.split("=", 1)
        section, option = key.rsplit(".", 1)
        if not parser.has_section(section):
            parser.add_section(section)
        parser[section][option] = value
    with ini.open("w") as f:
        parser.write(f)


def main() -> None:
    """Dispatch the subcommand."""
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd in ("inputs", "readmes"):
        conf = load_conf()
        for case in CASES:
            (ROOT / case).mkdir(exist_ok=True)
            if cmd == "inputs":
                (ROOT / case / "input.ini").write_text(case_ini(case, conf))
            (ROOT / case / "README.md").write_text(case_readme(case, conf))
    elif cmd == "cases":
        print("\n".join(CASES))
    elif cmd == "set":
        set_options(Path(args[0]), args[1:])
    elif cmd == "restart-step":
        it, time = struct.unpack("=id", Path(args[0]).read_bytes()[:12])
        print(it, repr(time))
    else:
        raise SystemExit(f"make_inputs: unknown subcommand '{cmd}'\n{__doc__}")


if __name__ == "__main__":
    main()
