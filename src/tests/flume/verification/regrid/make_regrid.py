#!/usr/bin/env python3
"""Write the inputs of the FLUME runtime regridding checks (issue #74: check.sh legs rg3-rg5, accuracy.sh legs av2-av7).

Why this generator exists: the runtime regrid legs need cases where the grid really changes during the run, built
from inputs the suite already verifies, so that only the regrid is new:

- ``vortex``: the isentropic vortex in a periodic box (``regression/vortex-periodic``) on a quadtree, base level 1,
  ``max_level`` 3, 8x8 cells per block, a Loehner marker on the density, a regrid every 5 steps, reflux on. The vortex
  travels across blocks, so regrids refine ahead of it and derefine behind it; nothing crosses the boundary, so every
  volume integral must stay constant (rg3) and a restart taken at a regrid step must continue bitwise (rg4).
- ``cylinder``: the shock over the cylinder (``verification/shock-cylinder``) with the solid marker and a Loehner
  marker on the density, a regrid every 5 steps (rg5).
- ``vortex-av``: the accuracy variant of the vortex (AV-2): 16x16 cells per block, a base level, a finest level and
  either runtime tracking (a gradient marker on the density, which also refines the initial grid, a regrid every
  ``--frequency`` steps) or a uniform grid at the base level (``--frequency 0``); the output is named ``vortex-av``
  so that ``vortex/vortex_oracle.py`` reads it.
- ``blast-av``: the Balsara-Spicer MHD blast of the ``blast-limiter`` golden (EGLM, positivity limiter, the disk
  sampled at N = 64), tracked from ``--base-level`` to ``--max-level`` (a gradient marker on the pressure at
  initialisation and the regrids, OR a Loehner marker on it) or uniform at the base level (AV-3b).
- ``loop-av``: the field loop of MV-9 (GLM, or EGLM with ``--eglm``), ``--cells`` cells along x at the base level,
  tracked to ``--max-level`` (a gradient marker on B_x at initialisation and the regrids, OR Loehner markers on B_x and
  B_y) or uniform (AV-4).
- ``cylinder-av``: the shock over the cylinder of V6 to t = 0.25, tracked (the solid marker OR a Loehner marker on the
  density, a regrid every ``--frequency`` steps) or uniform at the finest level of V6 (``--uniform``) (AV-7).
- ``sod-av``: the Sod tube of V1 (``verification/sod/sod-x.ini``) on a quadtree with y null (12 x 4 cells per block,
  so that a level adds blocks along x and y only, not z), tracked from ``--base-level`` to ``--max-level`` (a gradient marker on the density, at initialisation and at
  the regrids, OR a Loehner marker on it) or uniform at the base level (``--frequency 0``) (AV-3a).

Usage:
    make_regrid.py vortex <out.ini> [--prolongation linear|conservative] [--it-max N] [--restart-save N] [--restart]
    make_regrid.py cylinder <out.ini> [--it-max N]
    make_regrid.py vortex-av <out.ini> --base-level L --max-level M [--frequency N] [--time-max T] [--cfl C]
    make_regrid.py sod-av <out.ini> --base-level L --max-level M [--frequency N]
    make_regrid.py blast-av <out.ini> --base-level L --max-level M [--frequency N] [--time-max T]
    make_regrid.py cylinder-av <out.ini> [--uniform] [--frequency N]
    make_regrid.py loop-av <out.ini> --cells N [--max-level M --frequency N] [--eglm] [--time-max T]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
VORTEX = REPO / "src/tests/flume/regression/vortex-periodic/input.ini"
CYLINDER = REPO / "src/tests/flume/verification/shock-cylinder/shock-cylinder.ini"
SOD = REPO / "src/tests/flume/verification/sod/sod-x.ini"
BLAST = REPO / "src/tests/flume/regression/blast-limiter/input.ini"
LOOP = REPO / "src/tests/flume/verification/mhd/field-loop/make_field_loop.py"
LOOP_BASE = REPO / "src/tests/flume/verification/vortex/vortex-n064.ini"

LOHNER = """[amr_marker_{n}]
mode         = 4             ; Loehner estimator (issue #74)
delta_type   = max           ; unused by the Loehner marker; required by the common marker parse
delta_fine   = 0.0
delta_coarse = 0.0
field        = {field}             ; 1 conservative, 2 auxiliary variables
var          = {var}
refine_tol   = {refine}
derefine_tol = {derefine}
floor        = {floor}
buffer       = 2

"""

GRADIENT = """[amr_marker_{n}]
mode         = 2             ; gradient
delta_type   = max
delta_fine   = {fine}
delta_coarse = {coarse}
field        = {field}             ; 1 conservative, 2 auxiliary variables
var          = {var}
tol          = {tol}

"""


def setkey(text: str, key: str, value: str) -> str:
    """Return ``text`` with every ``key = ...`` line set to ``value`` (the key must exist).

    Parameters
    ----------
    text : str
        INI text.
    key : str
        Option name.
    value : str
        New value.

    Returns
    -------
    str
        The edited text.

    Raises
    ------
    KeyError
        If the key is not in the text.
    """
    pattern = re.compile(rf"(?m)^({re.escape(key)}\s*=).*$")
    if not pattern.search(text):
        raise KeyError(key)
    return pattern.sub(lambda m: f"{m.group(1)} {value}", text)


def vortex(args: argparse.Namespace) -> str:
    """Return the input of the travelling vortex case."""
    text = VORTEX.read_text()
    for key, value in (
        ("ni", "8"),
        ("nj", "8"),
        ("iu_ref_levels", "1"),
        ("max_level", "3"),
        ("iters", "2"),
        ("frequency", f"5\nregrid_prolongation = {args.prolongation}"),
        ("markers_number", "1"),
        ("output_basename", "vortex-rg"),
        ("restart_basename", "vortex-rg-restart"),
        ("restart_save", str(args.restart_save)),
        ("restart", ".true." if args.restart else ".false."),
        ("it_max", str(args.it_max)),
        ("time_max", "10.0"),
        ("CFL", "0.4"),
    ):
        text = setkey(text, key, value)
    return text.replace(
        "[field]", LOHNER.format(n=1, field=1, var=1, refine="0.3", derefine="0.1", floor="0.0") + "[field]", 1
    )


def vortex_av(args: argparse.Namespace) -> str:
    """Return the input of the vortex accuracy case: tracked (frequency > 0) or uniform at the base level."""
    text = VORTEX.read_text()
    tracked = args.frequency > 0
    for key, value in (
        ("ni", "16"),
        ("nj", "16"),
        ("iu_ref_levels", str(args.base_level)),
        ("max_level", str(args.max_level)),
        ("iters", "2"),
        ("frequency", str(args.frequency)),
        ("markers_number", "1" if tracked else "0"),
        ("amr_iterations", str(args.max_level - args.base_level) if tracked else "0"),
        ("output_basename", "vortex-av"),
        ("restart_basename", "vortex-av-restart"),
        ("it_max", "-1"),
        ("time_max", args.time_max),
        ("CFL", args.cfl),
    ):
        text = setkey(text, key, value)
    if tracked:
        # the gradient marker acts at initialisation too (the Loehner one only at the runtime regrids), so the vortex
        # starts on the finest level: the finest spacing where |grad rho| > tol, the base spacing elsewhere
        fine, coarse = f"{1.01 / (16 * 2**args.max_level):.6g}", f"{1.01 / (16 * 2**args.base_level):.6g}"
        text = text.replace(
            "[field]", GRADIENT.format(n=1, field=1, var=1, fine=fine, coarse=coarse, tol=args.tol) + "[field]", 1
        )
    return text


def sod_av(args: argparse.Namespace) -> str:
    """Return the input of the Sod accuracy case: tracked (frequency > 0) or uniform at the base level."""
    text = SOD.read_text()
    tracked = args.frequency > 0
    for key, value in (
        ("ni", "12"),
        ("nj", "4"),
        ("nk", "1"),
        ("ratio", "4"),
        ("iu_ref_levels", str(args.base_level)),
        ("max_level", str(args.max_level)),
        ("iters", "2"),
        ("frequency", str(args.frequency)),
        ("markers_number", "2" if tracked else "0"),
        ("amr_iterations", str(args.max_level - args.base_level) if tracked else "0"),
        ("output_basename", "sod-av"),
        ("restart_basename", "sod-av-restart"),
    ):
        text = setkey(text, key, value)
    # the domain [0, length]: on the default [0, 1] the initial jump at x = 0.5 is a block face at every level, the case
    # of issue #76 (the gradient marker, which then read interior cells only, could not see it to refine the initial
    # grid); with 1.2 it is a face at no level (the faces are multiples of 0.3 / 2^(level - 2))
    text = re.sub(r"(?m)^(emax_x\s*=).*$", rf"\g<1> {args.length}", text, count=1)
    head, sep, tail = text.partition("[initial_conditions_region_2]")
    text = head + sep + re.sub(r"(?m)^(emax_x\s*=).*$", rf"\g<1> {args.length}", tail, count=1)
    if tracked:
        # 1 % above the level spacings: the marker compares spacings with ">" and "<=", so an ulp must not decide
        fine = f"{1.01 * args.length / (12 * 2**args.max_level):.6g}"
        coarse = f"{1.01 * args.length / (12 * 2**args.base_level):.6g}"
        markers = GRADIENT.format(n=1, field=1, var=1, fine=fine, coarse=coarse, tol=args.tol)
        markers += LOHNER.format(n=2, field=1, var=1, refine=args.refine, derefine=args.derefine, floor="0.0")
        text = text.replace("[field]", markers + "[field]", 1)
    return text


def blast_av(args: argparse.Namespace) -> str:
    """Return the input of the MHD blast accuracy case: tracked (frequency > 0) or uniform at the base level."""
    text = BLAST.read_text()
    tracked = args.frequency > 0
    for key, value in (
        ("iu_ref_levels", str(args.base_level)),
        ("max_level", str(args.max_level)),
        ("iters", "2"),
        ("frequency", str(args.frequency)),
        ("markers_number", "2" if tracked else "0"),
        ("amr_iterations", str(args.max_level - args.base_level) if tracked else "0"),
        ("output_basename", "blast-av"),
        ("restart_basename", "blast-av-restart"),
        ("it_max", "-1"),
        ("time_max", args.time_max),
    ):
        text = setkey(text, key, value)
    if tracked:
        # the pressure (auxiliary variable 5): the gradient marker refines the disk at initialisation, the Loehner one
        # follows the waves at the regrids
        fine = f"{1.01 / (16 * 2**args.max_level):.6g}"
        coarse = f"{1.01 / (16 * 2**args.base_level):.6g}"
        markers = GRADIENT.format(n=1, field=2, var=5, fine=fine, coarse=coarse, tol=args.tol)
        markers += LOHNER.format(n=2, field=2, var=5, refine=args.refine, derefine=args.derefine, floor="0.0")
        text = text.replace("[field]", markers + "[field]", 1)
    return text


def loop_av(args: argparse.Namespace) -> str:
    """Return the input of the field loop accuracy case: MV-9 (mhd/field-loop/make_field_loop.py) at ``--cells``
    cells along x on its 4x4-block quadtree, tracked to ``--max-level`` (frequency > 0) or uniform."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "loop.ini"
        subprocess.run(
            [
                sys.executable,
                "-I",
                str(LOOP),
                str(LOOP_BASE),
                str(out),
                "--cells",
                str(args.cells),
                "--time-max",
                args.time_max,
            ],
            check=True,
        )
        text = out.read_text()
    tracked = args.frequency > 0
    for key, value in (
        ("max_level", str(args.max_level if tracked else 2)),
        ("iters", "2"),
        ("frequency", str(args.frequency)),
        ("markers_number", "3" if tracked else "0"),
        ("amr_iterations", str(args.max_level - 2) if tracked else "0"),
    ):
        text = setkey(text, key, value)
    if args.eglm:
        text = setkey(text, "divergence_control", "eglm")
    if tracked:
        # the loop edge r = R, where B jumps: a gradient marker on B_x refines it at initialisation (|grad B_x| ~ A0 / h
        # there, ~ A0 / R inside the loop), the Loehner markers on B_x and B_y follow it at the regrids
        h = 2.0 / args.cells
        fine = f"{1.01 * h / 2 ** (args.max_level - 2):.6g}"
        markers = GRADIENT.format(n=1, field=1, var=6, fine=fine, coarse=f"{1.01 * h:.6g}", tol=args.tol)
        # floor: a tenth of the loop amplitude A0 = 1e-3; outside the loop B is GLM residue (~1e-6), which the relative
        # filter alone reads as E ~ 1 (#74 P4)
        markers += LOHNER.format(n=2, field=1, var=6, refine=args.refine, derefine=args.derefine, floor="1e-4")
        markers += LOHNER.format(n=3, field=1, var=7, refine=args.refine, derefine=args.derefine, floor="1e-4")
        text = text.replace("[field]", markers + "[field]", 1)
    return text


def cylinder_av(args: argparse.Namespace) -> str:
    """Return the input of the cylinder accuracy case (AV-7) to t = 0.25: tracked (V6 with the solid marker and a
    Loehner marker on the density, a regrid every ``--frequency`` steps), or uniform at the finest level of V6
    (``--uniform``)."""
    text = CYLINDER.read_text()
    if args.uniform:
        for key, value in (("iu_ref_levels", "3"), ("markers_number", "0"), ("amr_iterations", "0")):
            text = setkey(text, key, value)
        return text
    for key, value in (("frequency", str(args.frequency)), ("markers_number", "2")):
        text = setkey(text, key, value)
    return text.replace(
        "[field]", LOHNER.format(n=2, field=1, var=1, refine="0.6", derefine="0.2", floor="0.0") + "[field]", 1
    )


def cylinder(args: argparse.Namespace) -> str:
    """Return the input of the shock over the cylinder with runtime regridding."""
    text = CYLINDER.read_text()
    for key, value in (
        ("frequency", "5"),
        ("markers_number", "2"),
        ("it_max", str(args.it_max)),
    ):
        text = setkey(text, key, value)
    return text.replace(
        "[field]", LOHNER.format(n=2, field=1, var=1, refine="0.6", derefine="0.2", floor="0.0") + "[field]", 1
    )


def main() -> None:
    """Parse the command line and write the input."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "case", choices=("vortex", "cylinder", "vortex-av", "sod-av", "blast-av", "loop-av", "cylinder-av")
    )
    parser.add_argument("out", type=Path)
    parser.add_argument("--prolongation", choices=("linear", "conservative"), default="conservative")
    parser.add_argument("--it-max", type=int, default=200)
    parser.add_argument("--restart-save", type=int, default=0)
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--base-level", type=int, default=1)
    parser.add_argument("--max-level", type=int, default=3)
    parser.add_argument("--frequency", type=int, default=5)
    parser.add_argument("--time-max", default="0.2")
    parser.add_argument("--cfl", default="0.4")
    parser.add_argument("--refine", default="0.3")
    parser.add_argument("--derefine", default="0.1")
    parser.add_argument("--tol", default="0.5", help="gradient marker tolerance (vortex-av)")
    parser.add_argument("--length", type=float, default=1.0, help="domain length along x (sod-av)")
    parser.add_argument("--cells", type=int, default=32, help="cells along x at the base level (loop-av)")
    parser.add_argument("--eglm", action="store_true", help="EGLM cleaning instead of GLM (loop-av)")
    parser.add_argument("--uniform", action="store_true", help="uniform at the finest level (cylinder-av)")
    args = parser.parse_args()
    text = {
        "vortex": vortex,
        "cylinder": cylinder,
        "vortex-av": vortex_av,
        "sod-av": sod_av,
        "blast-av": blast_av,
        "loop-av": loop_av,
        "cylinder-av": cylinder_av,
    }[args.case](args)
    args.out.write_text(text)


if __name__ == "__main__":
    main()
