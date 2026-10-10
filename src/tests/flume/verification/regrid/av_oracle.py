#!/usr/bin/env python3
"""FLUME runtime AMR verification oracle (issue #74, M5-P4, legs AV-2 to AV-7 of regrid/accuracy.sh).

Why this oracle exists: a runtime-regridded run must be about as accurate as the uniform run at its finest level,
while holding fewer cells, and must not let a regrid spoil what the scheme guarantees (positivity, symmetry,
conservation). Each subcommand compares a tracked run with its uniform references on a quantity with an exact or
reference value, against bounds stated in the #74 P4 plan before any run:

- ``vortex``: the isentropic vortex (AV-2). L1(rho) against the exact solution (vortex/vortex_oracle.py) of the tracked
  run at most ``--fine-ratio`` times the uniform run at the finest level, and at most ``--coarse-ratio`` times the
  uniform run at the base level.
- ``sod``: the Sod tube (AV-3a). L1(rho) against the exact Riemann solution, each cell weighted by its own width
  (the tracked grid mixes levels), at most ``--fine-ratio`` times the uniform-fine run; the contact and the shock
  must sit on finest-level cells at the final time; the copies along the null directions agree within
  ``--spread-max`` (round-off of the seam ghost fill, present without any regrid).
- ``blast``: the Balsara-Spicer MHD blast (AV-3b). The outer shock (the outermost cell where the density departs from
  the ambient 1 by more than ``--threshold``) along the four half-axes through the centre within ``--cells`` finest
  cells of the uniform-fine run, and the tracked run point-symmetric about the centre (rho at (x, y) and at
  (1 - x, 1 - y)) within ``--mirror-tol`` relative to the largest density; the ideal MHD equations are invariant under
  B -> -B, so the reflection maps the blast onto itself (``--no-symmetry`` skips it: the uniform run breaks it too,
  issue #77).
- ``drift``: the volume integrals of a tracked run (conservation histories): the ``--conserved`` ones constant within
  ``--max-drift``, the ``--bounded`` ones drifting at most ``--ratio`` times the uniform-fine run's drift (integrals the
  scheme itself does not conserve, AV-3b: the regrid must add nothing), the rest reported (the momenta of the blast
  start at zero, so their relative drift is noise of the scheme's asymmetry, issue #77).
- ``loop``: the field loop (AV-4). From mhd/field-loop/field_loop_oracle.py (volume-weighted): the magnetic energy
  E_B(T)/E_B(0) of the tracked run at least the uniform-fine run's minus ``--energy-drop`` and its <|B_z|>/A0 (the
  divergence error the out-of-plane velocity turns into B_z) at most ``--bz-ratio`` times the uniform-fine run's; from
  the div(B) histories, max |div B| at the final step at most ``--divb-ratio`` times the uniform-fine run's (GLM or
  EGLM must absorb what the prolongation injects), and the jump of max |div B| at each regrid step reported.
- ``agree``: CPU against FNL (AV-8). Per field, the largest CPU-FNL difference of the tracked run (relative to the
  field's largest magnitude, the metric of conservation_oracle.py --compare) at most ``--ratio`` times the same
  difference of the uniform-fine runs, or below ``--tol``: the regrid must add no divergence between the backends. On
  MHD the backends differ without any regrid (B by ~1e-6 on the field loop, the blast at O(1), issue #77), so a fixed
  tolerance would test that, not the regrid.
- ``cylinder``: the shock over the cylinder (AV-7). The bow shock stand-off on the row next to the symmetry line
  y = 0.5 (the distance from the cylinder front, x = 0.37, to the largest density jump upstream of it) within
  ``--cells`` finest cells of the uniform-fine run. The mirror symmetry and the positivity are checked by
  shock-cylinder/shock_cylinder_oracle.py.

Usage:
    av_oracle.py vortex <tracked> --fine <uniform-fine> --coarse <uniform-coarse> [--fine-ratio R] [--coarse-ratio R]
    av_oracle.py sod <tracked> --fine <uniform-fine> [--fine-ratio R]
    av_oracle.py blast <tracked> --fine <uniform-fine> [--cells C] [--mirror-tol T]
    av_oracle.py drift <tracked> --fine <uniform-fine> --conserved NAME [NAME ...] [--max-drift D] [--ratio R]
    av_oracle.py cylinder <tracked> --fine <uniform-fine> [--cells C] [--report-only]
    av_oracle.py agree <cpu-tracked> <fnl-tracked> <cpu-uniform> <fnl-uniform> [--ratio R] [--tol T]
    av_oracle.py loop <tracked> --fine <uniform-fine> [--energy-drop D] [--bz-ratio R] [--divb-ratio R]
"""

from __future__ import annotations

import argparse
import configparser
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import h5py
import numpy as np

VERIFICATION = Path(__file__).resolve().parents[1]


def load(path: Path) -> ModuleType:
    """Import a sibling oracle by path (the verification directories are not packages).

    Parameters
    ----------
    path : Path
        Oracle source file.

    Returns
    -------
    ModuleType
        The imported module.
    """
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        sys.exit(f"av_oracle: cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verdict(ok: bool) -> str:
    """Return the PASS/FAIL tag."""
    return "PASS" if ok else "FAIL"


def vortex(args: argparse.Namespace) -> int:
    """AV-2: the tracked vortex against the uniform runs at its finest and base levels."""
    oracle = load(VERIFICATION / "vortex" / "vortex_oracle.py")
    l1 = {
        name: oracle.run_error(work)[1]
        for name, work in (("tracked", args.run), ("fine", args.fine), ("coarse", args.coarse))
    }
    for name in ("coarse", "fine", "tracked"):
        print(f"{name:>7}: L1(rho) = {l1[name]:.6e}")
    status = 0
    for ref, bound in (("fine", args.fine_ratio), ("coarse", args.coarse_ratio)):
        ratio = l1["tracked"] / l1[ref]
        ok = ratio <= bound
        status |= 0 if ok else 1
        print(f"tracked / {ref}: {ratio:.3f}  {verdict(ok)} (bound {bound:g})")
    return status


def sod_profile(work: Path, spread_max: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return the cell centres, widths and densities of the last checkpoint of a 1-D run on any AMR grid, and the
    largest spread of its copies along the null directions.

    On a uniform grid the copies are bitwise equal (V1 asserts it); across a coarse-fine seam the ghost fill of a null
    direction differs by round-off between copies (measured 7.5e-14 on a static box-refined grid, with no regrid), so
    the profile takes the mean of the copies and the spread is bounded by ``spread_max`` instead. The widths follow
    from the sorted centres: the first cell starts at the domain minimum, each face is its left face plus twice the
    distance to its centre.
    """
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    ngc = int(ini["grid"]["ngc"])
    basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
    files = sorted(work.glob(f"{basename}-*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    cells: dict[float, float] = {}
    spread = 0.0
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                rho = h5[f"{blk}-r"][()].transpose(2, 1, 0)
                rho = rho[ngc:-ngc, ngc:-ngc, ngc:-ngc] if rho.shape[2] > 2 * ngc else rho[ngc:-ngc, ngc:-ngc, :]
                rho = rho.reshape(rho.shape[0], -1)
                spread = max(spread, float(np.max(np.abs(rho - rho[:, :1]))))
                centres = origin[0] + (np.arange(rho.shape[0]) + ngc + 0.5) * dxyz[0]
                for c, value in zip(centres, rho.mean(axis=1), strict=True):
                    cells.setdefault(float(np.round(c, 12)), float(value))
    xs = np.array(sorted(cells))
    widths = np.empty_like(xs)
    face = float(ini["grid"]["emin_x"])
    for i, c in enumerate(xs):
        widths[i] = 2.0 * (c - face)
        face += widths[i]
    if abs(face - float(ini["grid"]["emax_x"])) > 1e-9 or np.any(widths <= 0.0):
        sys.exit(f"av_oracle: the cells of {work} do not tile the domain")
    if spread > spread_max:
        print(f"   {work.name}: copies along the null directions differ by {spread:.3e}  FAIL (bound {spread_max:g})")
    return xs, widths, np.array([cells[c] for c in xs]), spread


def sod(args: argparse.Namespace) -> int:
    """AV-3a: the tracked Sod tube against the exact solution and the uniform-fine run."""
    oracle = load(VERIFICATION / "sod" / "sod_oracle.py")
    ini = oracle.read_ini(next(args.run.glob("*.ini")))  # physics, time and the two states
    gamma = float(ini["physics"]["cp"]) / float(ini["physics"]["cv"])
    t = float(ini["time"]["time_max"])
    r1, r2 = ini["initial_conditions_region_1"], ini["initial_conditions_region_2"]
    left = (float(r1["r"]), float(r1["u"]), float(r1["p"]))
    right = (float(r2["r"]), float(r2["u"]), float(r2["p"]))
    l1 = {}
    status = 0
    for name, work in (("fine", args.fine), ("tracked", args.run)):
        xs, widths, rho, spread = sod_profile(work, args.spread_max)
        status |= 0 if spread <= args.spread_max else 1
        l1[name] = float(np.sum(np.abs(rho - oracle.exact_riemann_density(xs, t, 0.5, gamma, left, right)) * widths))
        print(f"{name:>7}: {xs.size} cells, L1(rho) = {l1[name]:.6e}, spread of the null copies {spread:.1e}")
        if name == "tracked":
            # the contact and the shock of the exact solution: the jumps of a finely sampled exact density
            sample = np.linspace(0.0, 1.0, 200001)
            jumps = sample[1:][np.abs(np.diff(oracle.exact_riemann_density(sample, t, 0.5, gamma, left, right))) > 1e-3]
            finest = widths.min()
            for where in jumps:
                width = widths[np.argmin(np.abs(xs - where))]
                ok = np.isclose(width, finest)
                status |= 0 if ok else 1
                print(
                    f"   discontinuity at x = {where:.4f}: cell width {width:.4e} (finest {finest:.4e})  {verdict(ok)}"
                )
    ratio = l1["tracked"] / l1["fine"]
    ok = ratio <= args.fine_ratio
    status |= 0 if ok else 1
    print(f"tracked / fine: {ratio:.3f}  {verdict(ok)} (bound {args.fine_ratio:g})")
    return status


def plane_cells(work: Path) -> dict[tuple[float, float], tuple[float, float]]:
    """Return the cells of the last checkpoint of a 2-D run (null z): (x, y) centre -> (width, density)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None)
    ini.read(next(work.glob("*.ini")))
    ngc = int(ini["grid"]["ngc"])
    basename = next(work.glob("*-residuals.dat")).name.removesuffix("-residuals.dat")
    files = sorted(work.glob(f"{basename}-*-proc*.h5"))
    last = max(int(f.name.split("-")[-2]) for f in files)
    cells: dict[tuple[float, float], tuple[float, float]] = {}
    for path in (f for f in files if int(f.name.split("-")[-2]) == last):
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                origin, dxyz = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                rho = h5[f"{blk}-r"][()].transpose(2, 1, 0)
                # the first interior layer along the null z (the checkpoints hold its ghosts when it has them)
                rho = rho[ngc:-ngc, ngc:-ngc, ngc if rho.shape[2] > 2 * ngc else 0]
                xc = origin[0] + (np.arange(rho.shape[0]) + ngc + 0.5) * dxyz[0]
                yc = origin[1] + (np.arange(rho.shape[1]) + ngc + 0.5) * dxyz[1]
                for i, x in enumerate(xc):
                    for j, y in enumerate(yc):
                        cells[(round(float(x), 12), round(float(y), 12))] = (float(dxyz[0]), float(rho[i, j]))
    return cells


def shock_radii(cells: dict[tuple[float, float], tuple[float, float]], threshold: float) -> list[float]:
    """Return the outer shock radius along +x, -x, +y, -y from the centre (0.5, 0.5).

    Along each half-axis the cells are those whose centre lies within half a width of the axis, on the side y > 0.5
    (x > 0.5 for the y half-axes): 0.5 is a block face, so the row just above it is the closest to the axis on every
    grid. The radius is the distance to the outer face of the outermost cell departing from the ambient density.
    """
    radii = []
    for axis, sign in ((0, 1.0), (0, -1.0), (1, 1.0), (1, -1.0)):
        other = 1 - axis
        line = [(c, w, r) for c, (w, r) in cells.items() if 0.5 < c[other] < 0.5 + w and sign * (c[axis] - 0.5) > 0.0]
        out = [abs(c[axis] - 0.5) + 0.5 * w for c, w, r in line if abs(r - 1.0) > threshold]
        radii.append(max(out) if out else 0.0)
    return radii


def blast(args: argparse.Namespace) -> int:
    """AV-3b: the tracked MHD blast against the uniform-fine run, and its point symmetry."""
    status = 0
    tracked, fine = plane_cells(args.run), plane_cells(args.fine)
    finest = min(w for w, _ in tracked.values())
    r_tracked, r_fine = shock_radii(tracked, args.threshold), shock_radii(fine, args.threshold)
    for name, a, b in zip(("+x", "-x", "+y", "-y"), r_tracked, r_fine, strict=True):
        ok = abs(a - b) <= args.cells * finest
        status |= 0 if ok else 1
        print(
            f"shock radius {name}: tracked {a:.5f}, uniform-fine {b:.5f}, |difference| = {abs(a - b) / finest:.2f} "
            f"finest cells  {verdict(ok)} (bound {args.cells:g})"
        )
    print(f"cells: tracked {len(tracked)}, uniform-fine {len(fine)}")
    if args.no_symmetry:
        return status
    rho_max = max(r for _, r in tracked.values())
    unmatched, asym = 0, 0.0
    for (x, y), (_, r) in tracked.items():
        mirror = tracked.get((round(1.0 - x, 12), round(1.0 - y, 12)))
        if mirror is None:
            unmatched += 1
        else:
            asym = max(asym, abs(r - mirror[1]) / rho_max)
    ok = unmatched == 0 and asym <= args.mirror_tol
    status |= 0 if ok else 1
    print(
        f"point symmetry: {len(tracked)} cells, {unmatched} without a mirror cell, max |rho(x) - rho(1 - x)| / "
        f"max rho = {asym:.3e}  {verdict(ok)} (tol {args.mirror_tol:g})"
    )
    return status


def drift(args: argparse.Namespace) -> int:
    """AV-3b: the integrals of the tracked run against the uniform-fine run's drift."""
    oracle = load(VERIFICATION / "conservation" / "conservation_oracle.py")
    names, d_tracked = oracle.drift(args.run)
    names_f, d_fine = oracle.drift(args.fine)
    if names != names_f:
        sys.exit(f"av_oracle: the histories hold different integrals ({names} and {names_f})")
    status = 0
    for name, a, b in zip(names, d_tracked, d_fine, strict=True):
        if name in args.conserved:
            ok = a <= args.max_drift
            text = f"conserved, bound {args.max_drift:g}"
        elif name not in args.bounded or b == 0.0:
            ok, text = True, "reported"
        else:
            ok = a <= args.ratio * b
            text = f"bound {args.ratio:g}x uniform-fine"
        status |= 0 if ok else 1
        print(f"drift {name:>3}: tracked {a:.3e}, uniform-fine {b:.3e}  {verdict(ok)} ({text})")
    return status


def divb_history(work: Path) -> dict[int, float]:
    """Return max |div B| per step from the div(B) history of a run."""
    rows = [line.split() for line in next(work.glob("*-divb_history.dat")).read_text().splitlines()]
    return {int(r[0]): float(r[2]) for r in rows if r and r[0][0] in "0123456789"}


def loop(args: argparse.Namespace) -> int:
    """AV-4: the tracked field loop against the uniform-fine run."""
    oracle = load(VERIFICATION / "mhd" / "field-loop" / "field_loop_oracle.py")
    (bz_t, eb_t), (bz_f, eb_f) = oracle.measures(args.run, 3), oracle.measures(args.fine, 3)
    status = 0
    ok = eb_t >= eb_f - args.energy_drop
    status |= 0 if ok else 1
    print(
        f"E_B(T)/E_B(0): tracked {eb_t:.6f}, uniform-fine {eb_f:.6f}  {verdict(ok)} (at least uniform-fine - "
        f"{args.energy_drop:g})"
    )
    ok = bz_t <= args.bz_ratio * bz_f or args.bz_report
    status |= 0 if ok else 1
    print(
        f"<|B_z|>/A0: tracked {bz_t:.4e}, uniform-fine {bz_f:.4e}, ratio {bz_t / bz_f:.3f}  "
        + ("reported (issue #78)" if args.bz_report else f"{verdict(ok)} (bound {args.bz_ratio:g})")
    )
    h_t, h_f = divb_history(args.run), divb_history(args.fine)
    last_t, last_f = h_t[max(h_t)], h_f[max(h_f)]
    ok = last_t <= args.divb_ratio * last_f
    status |= 0 if ok else 1
    print(
        f"max |div B| at the final step: tracked {last_t:.4e}, uniform-fine {last_f:.4e}, ratio "
        f"{last_t / last_f:.3f}  {verdict(ok)} (bound {args.divb_ratio:g})"
    )
    log = (args.run / "log.txt").read_text(errors="replace").splitlines()
    steps = [int(line.split("step ")[1].split(":")[0]) for line in log if "flume: regrid at step" in line]
    for step in steps:
        before, after = h_t.get(step - 1), h_t.get(step)
        if before is not None and after is not None:
            print(f"   regrid at step {step}: max |div B| {before:.4e} -> {after:.4e} (x{after / before:.2f})")
    return status


def field_differences(a_dir: Path, b_dir: Path, ngc: int) -> tuple[tuple[str, ...], np.ndarray]:
    """Return the field names and, per field, the largest |a - b| over the blocks relative to the field's largest
    magnitude in ``a`` (an identically zero field is compared absolutely)."""
    oracle = load(VERIFICATION / "conservation" / "conservation_oracle.py")
    (names, a), (names_b, b) = oracle.last_fields(a_dir), oracle.last_fields(b_dir)
    if names != names_b or a.keys() != b.keys():
        sys.exit(f"av_oracle: {a_dir.name} and {b_dir.name} hold different fields or blocks")
    g = ngc
    a = {k: v[:, g:-g, g:-g, g:-g] for k, v in a.items()}
    b = {k: v[:, g:-g, g:-g, g:-g] for k, v in b.items()}
    scale = np.max([np.max(np.abs(a[k]), axis=(1, 2, 3)) for k in a], axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    delta = np.max([np.max(np.abs(a[k] - b[k]), axis=(1, 2, 3)) for k in a], axis=0)
    return names, delta / scale


def agree(args: argparse.Namespace) -> int:
    """AV-8: the CPU-FNL difference of a tracked run against the same difference of the uniform-fine runs."""
    names, tracked = field_differences(args.cpu, args.fnl, 3)
    _, uniform = field_differences(args.cpu_uniform, args.fnl_uniform, 3)
    status = 0
    for name, t, u in zip(names, tracked, uniform, strict=True):
        ok = t <= max(args.ratio * u, args.tol)
        status |= 0 if ok else 1
        print(f"   {name:>3}: CPU-FNL tracked {t:.3e}, uniform {u:.3e}  {verdict(ok)}")
    print(
        f"{args.cpu.name} vs {args.fnl.name}: largest {tracked.max():.3e} against uniform {uniform.max():.3e}  "
        f"{verdict(status == 0)} (per field <= {args.ratio:g}x uniform or <= {args.tol:g})"
    )
    return status


def stand_off(cells: dict[tuple[float, float], tuple[float, float]], front: float) -> float:
    """Return the bow shock stand-off on the row just above y = 0.5: the distance from ``front`` to the face between
    the two consecutive cells upstream of it with the largest density jump."""
    row = sorted((c[0], w, r) for c, (w, r) in cells.items() if 0.5 < c[1] < 0.5 + w and c[0] < front)
    jumps = [(abs(b[2] - a[2]), 0.5 * (a[0] + b[0])) for a, b in zip(row[:-1], row[1:], strict=True)]
    return front - max(jumps)[1]


def cylinder(args: argparse.Namespace) -> int:
    """AV-7: the bow shock stand-off of a run of the shock over the cylinder against the uniform-fine run."""
    tracked, fine = plane_cells(args.run), plane_cells(args.fine)
    finest = min(w for w, _ in tracked.values())
    a, b = stand_off(tracked, args.front), stand_off(fine, args.front)
    ok = abs(a - b) <= args.cells * finest
    tag = "reported" if args.report_only else f"{verdict(ok)} (bound {args.cells:g})"
    print(
        f"bow shock stand-off of {args.run.name}: {a:.5f}, uniform-fine {b:.5f}, |difference| = "
        f"{abs(a - b) / finest:.2f} finest cells  {tag}"
    )
    print(f"cells: {args.run.name} {len(tracked)}, uniform-fine {len(fine)}")
    return 0 if ok or args.report_only else 1


def main() -> int:
    """Parse the command line, run the requested check, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="case", required=True)
    p = sub.add_parser("vortex")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--coarse", type=Path, required=True)
    p.add_argument("--fine-ratio", type=float, default=1.5)
    p.add_argument("--coarse-ratio", type=float, default=0.125)
    p.set_defaults(func=vortex)
    p = sub.add_parser("sod")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--fine-ratio", type=float, default=1.3)
    p.add_argument("--spread-max", type=float, default=1e-12)
    p.set_defaults(func=sod)
    p = sub.add_parser("blast")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--cells", type=float, default=2.0)
    p.add_argument("--threshold", type=float, default=0.05)
    p.add_argument("--mirror-tol", type=float, default=1e-10)
    p.add_argument("--no-symmetry", action="store_true")
    p.set_defaults(func=blast)
    p = sub.add_parser("drift")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--conserved", nargs="+", required=True)
    p.add_argument("--max-drift", type=float, default=1e-13)
    p.add_argument("--bounded", nargs="+", default=["rE", "bx", "by", "bz"], help="integrals bounded by --ratio")
    p.add_argument("--ratio", type=float, default=1.1)
    p.set_defaults(func=drift)
    p = sub.add_parser("cylinder")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--cells", type=float, default=2.0)
    p.add_argument("--front", type=float, default=0.37, help="cylinder front (centre 0.5, radius 0.13)")
    p.add_argument("--report-only", action="store_true")
    p.set_defaults(func=cylinder)
    p = sub.add_parser("loop")
    p.add_argument("run", type=Path)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--energy-drop", type=float, default=0.005)
    p.add_argument("--bz-ratio", type=float, default=1.2)
    p.add_argument("--bz-report", action="store_true", help="report <|B_z|> only (EGLM at 2:1 seams, issue #78)")
    p.add_argument("--divb-ratio", type=float, default=2.0)
    p.set_defaults(func=loop)
    p = sub.add_parser("agree")
    p.add_argument("cpu", type=Path)
    p.add_argument("fnl", type=Path)
    p.add_argument("cpu_uniform", type=Path)
    p.add_argument("fnl_uniform", type=Path)
    p.add_argument("--ratio", type=float, default=2.0)
    p.add_argument("--tol", type=float, default=1e-10)
    p.set_defaults(func=agree)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
