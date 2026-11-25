#!/usr/bin/env python3
"""A priori calibration of the Kolmogorov AMR marker (mode 7, r = h/eta) on a FLUME checkpoint (issue #133, step 7.2).

Reads a checkpoint written on a UNIFORM grid (every block at the finest level, e.g. the settled state of the leg-3
channel) and its input, and evaluates the marker exactly as `mark_by_kolmogorov` + `kolmogorov_block_max` do
(`src/app/flume/common/adam_flume_{common_object,kolmogorov_library}.F90`):

    r = h (2 S_d:S_d)^(1/4) / nu^(1/2),   S_d = S - tr(S) I / 3,   S = (grad u + grad u^T) / 2,
    nu = mu0 L(T) / rho,   L(T) = (T / T_ref)^omega [(T_ref + S) / (T + S) under Sutherland],   T = p / (rho R),
    h = (prod_d dx_d)^(1/n) over the n active axes,

the velocity gradient by second-order centred differences across internal and periodic faces (one ghost layer) and
first-order one-sided differences at the first/last cell next to a non-periodic domain face; a block votes on the
LARGEST r over its interior cells: refine above `refine_tol`, derefine below `derefine_tol`. A sibling group (the 8
children of one parent) derefines only when all 8 vote derefine (`tree%sanitize`).

The checkpoint's cells are assembled into one global array (all blocks share one spacing), the ghost layer of a periodic
axis is the wrap of the opposite side: this is what the marker reads after its caller refreshed the ghosts, and does not
depend on the freshness of the ghosts stored in the checkpoint (`--stored-ghosts` cross-checks the per-block evaluation on
the stored ghosts). One level coarser, the field is the 8-cell average of the conservative variables (the restriction of
a derefinement, `adam_field_object`), and r is evaluated on it with the spacing doubled.

`--block i,j,k` re-tiles the same cells into blocks of another shape (a what-if of the block layout: the sibling groups,
hence what can derefine, depend on it).

Usage:
    kolmogorov_apriori.py <dir> <input.ini> [--iteration N] [--block 64,48,32] [--pairs 4:1.9,8:3.5 ...]
                          [--stored-ghosts]
"""

from __future__ import annotations

import argparse
import configparser
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

FIELDS = ("r", "ru", "rv", "rw", "rE")


@dataclass(frozen=True)
class Gas:
    """The gas and its viscosity law, as `physics` and `dissipation%laws` give them to the marker."""

    gamma: float
    rgas: float
    mu0: float
    tref: float
    suth: float
    omega: float

    def nu(self, rho: np.ndarray, temp: np.ndarray) -> np.ndarray:
        """Return the kinematic viscosity `mu0 L(T) / rho` (`temperature_law`).

        Parameters
        ----------
        rho : np.ndarray
            Density.
        temp : np.ndarray
            Temperature.

        Returns
        -------
        np.ndarray
            Kinematic viscosity.
        """
        law = (temp / self.tref) ** self.omega
        if self.suth > 0.0 and self.omega > 0.0:
            law = law * (self.tref + self.suth) / (temp + self.suth)
        return self.mu0 * law / rho


@dataclass(frozen=True)
class Setup:
    """What the marker needs from the input: the gas, the domain, the axes, the block cells and ghosts."""

    gas: Gas
    emin: np.ndarray
    emax: np.ndarray
    active: np.ndarray
    periodic: np.ndarray
    cells: tuple[int, int, int]
    ngc: int
    tolerances: tuple[float, float] | None


def read_ini(path: Path) -> dict[str, dict[str, str]]:
    """Return the sections and keys of an INI file (inline `;` comments stripped).

    Parameters
    ----------
    path : Path
        INI file.

    Returns
    -------
    dict[str, dict[str, str]]
        Sections and keys.
    """
    parser = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    parser.optionxform = str  # type: ignore[assignment,method-assign]
    parser.read_string(path.read_text())
    return {s: dict(parser[s]) for s in parser.sections()}


def is_true(text: str) -> bool:
    """Return the value of a Fortran logical spelled in an INI file.

    Parameters
    ----------
    text : str
        `.true.`, `T`, `true`, ...

    Returns
    -------
    bool
        Its value.
    """
    return text.strip().strip(".").lower() in ("true", "t")


def load_setup(path: Path) -> Setup:
    """Return the marker's inputs from a FLUME input file.

    Parameters
    ----------
    path : Path
        The run's input.ini.

    Returns
    -------
    Setup
        Gas, domain, axes, block layout and the marker's tolerances (None without a mode-7 marker).

    Raises
    ------
    SystemExit
        On an input the marker would refuse or this script does not model.
    """
    ini = read_ini(path)
    phys, grid = ini["physics"], ini["grid"]
    if "gamma" in phys:
        gamma, rgas = float(phys["gamma"]), 1.0
    else:
        cp, cv = float(phys["cp"]), float(phys["cv"])
        gamma, rgas = cp / cv, cp - cv
    if "viscosity" in phys:
        mu0 = float(phys["viscosity"])
    elif "reynolds" in phys:
        mu0 = 1.0 / float(phys["reynolds"])
    else:
        sys.exit(f"{path}: [physics] gives no viscosity (viscosity or reynolds): the marker needs one")
    law = phys.get("viscosity_law", "constant").strip()
    tref, suth, omega = 1.0, 0.0, 0.0
    if law == "power-law":
        tref, omega = float(phys["reference_temperature"]), float(phys["viscosity_exponent"])
    elif law == "sutherland":
        if "sutherland_temperature" not in phys:
            sys.exit(f"{path}: Sutherland from temperature_unit_kelvin is not modelled here: give sutherland_temperature")
        tref, suth, omega = float(phys["reference_temperature"]), float(phys["sutherland_temperature"]), 1.5
    elif law != "constant":
        sys.exit(f"{path}: unknown viscosity_law {law}")
    axes = "xyz"
    emin = np.array([float(grid[f"emin_{a}"]) for a in axes])
    emax = np.array([float(grid[f"emax_{a}"]) for a in axes])
    active = np.array([not is_true(grid.get(f"null_{a}", ".false.")) for a in axes])
    periodic = np.array([ini.get(f"bc_{a}_min", {}).get("type", "").strip() == "periodic" for a in axes])
    cells = (int(grid["ni"]), int(grid["nj"]), int(grid["nk"]))
    tolerances = None
    for name, sec in ini.items():
        if name.startswith("amr_marker_") and sec.get("mode", "").strip() == "7":
            tolerances = (float(sec["refine_tol"]), float(sec["derefine_tol"]))
    gas = Gas(gamma=gamma, rgas=rgas, mu0=mu0, tref=tref, suth=suth, omega=omega)
    return Setup(gas=gas, emin=emin, emax=emax, active=active, periodic=periodic, cells=cells, ngc=int(grid["ngc"]),
                 tolerances=tolerances)


def checkpoint_files(folder: Path, iteration: int | None) -> list[Path]:
    """Return the per-rank files of a checkpoint (`<base>-<iteration>-proc<rank>.h5`), the last one by default.

    Parameters
    ----------
    folder : Path
        Directory holding the checkpoint.
    iteration : int or None
        Iteration of the checkpoint, None for the last one.

    Returns
    -------
    list[Path]
        Its files, by rank.

    Raises
    ------
    SystemExit
        When no checkpoint is found.
    """
    pat = re.compile(r"^(?!.*restart).*-(\d+)-proc(\d+)\.h5$")
    found: dict[int, list[Path]] = {}
    for p in folder.iterdir():
        m = pat.match(p.name)
        if m:
            found.setdefault(int(m.group(1)), []).append(p)
    if not found:
        sys.exit(f"{folder}: no checkpoint <base>-<iteration>-proc<rank>.h5")
    it = max(found) if iteration is None else iteration
    if it not in found:
        sys.exit(f"{folder}: no checkpoint at iteration {it} (found {sorted(found)})")
    return sorted(found[it])


def read_blocks(files: list[Path], ngc: int) -> list[dict]:
    """Return the blocks of a checkpoint: interior corner, spacing and conservative fields WITH their ghosts, (x, y, z).

    Parameters
    ----------
    files : list[Path]
        Per-rank checkpoint files.
    ngc : int
        Ghost cells.

    Returns
    -------
    list[dict]
        Blocks: `lo` (interior corner), `d` (spacings), `q` (5, ni + 2 ngc, nj + 2 ngc, nk + 2 ngc).
    """
    out = []
    for path in files:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in FIELDS])
                out.append({"lo": origin + ngc * d, "d": d, "q": q})
    return out


def auxiliary(q: np.ndarray, gas: Gas) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return density, velocity (3, ...) and temperature of conservative variables (`conservative_to_auxiliary`).

    Parameters
    ----------
    q : np.ndarray
        Conservative variables (5, ...): rho, rho u, rho v, rho w, rho E.
    gas : Gas
        The gas.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        rho, velocity, temperature.
    """
    rho = q[0]
    vel = q[1:4] / rho
    p = (gas.gamma - 1.0) * (q[4] - 0.5 * rho * (vel**2).sum(axis=0))
    return rho, vel, p / (rho * gas.rgas)


def spacing(d: np.ndarray, active: np.ndarray) -> float:
    """Return the estimator's spacing, the geometric mean of the active spacings (`kolmogorov_spacing`).

    Parameters
    ----------
    d : np.ndarray
        Cell spacings.
    active : np.ndarray
        Active axes.

    Returns
    -------
    float
        Spacing (0 without active axis).
    """
    n = int(active.sum())
    return float(np.prod(d[active]) ** (1.0 / n)) if n else 0.0


def cell_ratio(q: np.ndarray, d: np.ndarray, is_open: np.ndarray, setup: Setup) -> np.ndarray:
    """Return r = h / eta on every interior cell of a padded array, the stencil of `kolmogorov_block_max`.

    Parameters
    ----------
    q : np.ndarray
        Conservative variables with ONE ghost layer on every side, (5, n1 + 2, n2 + 2, n3 + 2).
    d : np.ndarray
        Cell spacings.
    is_open : np.ndarray
        (2, 3): the minus/plus face of each axis is internal (its ghost layer is read; one-sided otherwise).
    setup : Setup
        Marker inputs.

    Returns
    -------
    np.ndarray
        r on the interior cells, (n1, n2, n3).
    """
    rho, vel, temp = auxiliary(q, setup.gas)
    n = np.array(q.shape[1:]) - 2
    inner = (slice(1, -1),) * 3
    grad = np.zeros((3, 3, *n))  # grad[a, b] = d u_a / d x_b
    for ax in range(3):
        if not (setup.active[ax] and n[ax] > 1):
            continue

        def shifted(s: int, ax: int = ax) -> np.ndarray:
            idx = list(inner)
            idx[ax] = slice(1 + s, q.shape[1 + ax] - 1 + s)
            return vel[(slice(None), *idx)]

        g = 0.5 * (shifted(1) - shifted(-1)) / d[ax]
        lo, hi = [slice(None)] * 4, [slice(None)] * 4
        lo[1 + ax], hi[1 + ax] = slice(0, 1), slice(-1, None)
        if not is_open[0, ax]:
            g[tuple(lo)] = ((shifted(1) - shifted(0)) / d[ax])[tuple(lo)]
        if not is_open[1, ax]:
            g[tuple(hi)] = ((shifted(0) - shifted(-1)) / d[ax])[tuple(hi)]
        grad[:, ax] = g
    s = 0.5 * (grad + grad.transpose(1, 0, 2, 3, 4))
    tr = s[0, 0] + s[1, 1] + s[2, 2]
    sd2 = np.maximum((s * s).sum(axis=(0, 1)) - tr * tr / 3.0, 0.0)
    nu = setup.gas.nu(rho[inner], temp[inner])
    return spacing(d, setup.active) * np.sqrt(np.sqrt(2.0 * sd2)) / np.sqrt(nu)


def assemble(blocks: list[dict], setup: Setup) -> tuple[np.ndarray, np.ndarray]:
    """Return the global interior array of a uniform checkpoint and its spacing.

    Parameters
    ----------
    blocks : list[dict]
        Blocks (`read_blocks`).
    setup : Setup
        Marker inputs.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Conservative variables (5, Nx, Ny, Nz) and spacing.

    Raises
    ------
    SystemExit
        When the blocks are not all at one spacing (the checkpoint is not uniform).
    """
    d = blocks[0]["d"]
    if any(not np.allclose(b["d"], d, rtol=1e-12) for b in blocks):
        sys.exit("the checkpoint is not uniform (blocks at several levels): this calibration needs a uniform one")
    g = setup.ngc
    ncell = np.array(blocks[0]["q"].shape[1:]) - 2 * g
    big = np.rint((setup.emax - setup.emin) / d).astype(int)
    out = np.full((5, *big), np.nan)
    for b in blocks:
        o = np.rint((b["lo"] - setup.emin) / d).astype(int)
        out[:, o[0]:o[0] + ncell[0], o[1]:o[1] + ncell[1], o[2]:o[2] + ncell[2]] = b["q"][:, g:-g, g:-g, g:-g]
    if np.isnan(out).any():
        sys.exit("the blocks do not tile the domain")
    return out, d


def pad(q: np.ndarray, periodic: np.ndarray) -> np.ndarray:
    """Return a global array padded with one ghost layer: the wrap along periodic axes, the edge copy elsewhere (not
    read: the stencil is one-sided there).

    Parameters
    ----------
    q : np.ndarray
        Conservative variables (5, Nx, Ny, Nz).
    periodic : np.ndarray
        Periodic axes.

    Returns
    -------
    np.ndarray
        (5, Nx + 2, Ny + 2, Nz + 2).
    """
    out = q
    for ax in range(3):
        out = np.pad(out, [(0, 0)] + [(1, 1) if a == ax else (0, 0) for a in range(3)],
                     mode="wrap" if periodic[ax] else "edge")
    return out


def restrict(q: np.ndarray) -> np.ndarray:
    """Return the field one level coarser: the mean of each 2 x 2 x 2 group of cells (the derefinement's restriction).

    Parameters
    ----------
    q : np.ndarray
        Conservative variables (5, Nx, Ny, Nz), even counts.

    Returns
    -------
    np.ndarray
        (5, Nx / 2, Ny / 2, Nz / 2).
    """
    nv, nx, ny, nz = q.shape
    return q.reshape(nv, nx // 2, 2, ny // 2, 2, nz // 2, 2).mean(axis=(2, 4, 6))


def block_reduce(r: np.ndarray, block: tuple[int, int, int], how: str = "max") -> np.ndarray:
    """Return a per-block reduction of a cell field tiled into blocks of `block` cells.

    Parameters
    ----------
    r : np.ndarray
        Cell field (Nx, Ny, Nz).
    block : tuple[int, int, int]
        Cells per block.
    how : str
        `max`, `min` or `mean`.

    Returns
    -------
    np.ndarray
        (Nx / bi, Ny / bj, Nz / bk).
    """
    bi, bj, bk = block
    nx, ny, nz = r.shape
    t = r.reshape(nx // bi, bi, ny // bj, bj, nz // bk, bk)
    return getattr(t, how)(axis=(1, 3, 5))


def stored_ghost_check(blocks: list[dict], setup: Setup, r_fine: np.ndarray) -> float:
    """Return the largest relative difference of the per-block r_max on the STORED ghosts from the assembled one.

    Parameters
    ----------
    blocks : list[dict]
        Blocks (`read_blocks`).
    setup : Setup
        Marker inputs.
    r_fine : np.ndarray
        r on the assembled grid.

    Returns
    -------
    float
        max |r_max(stored) - r_max(assembled)| / r_max(assembled).
    """
    g = setup.ngc
    worst = 0.0
    for b in blocks:
        d, lo = b["d"], b["lo"]
        n = np.array(b["q"].shape[1:]) - 2 * g
        hi = lo + n * d
        is_open = np.array([setup.periodic | (lo > setup.emin + 0.5 * d), setup.periodic | (hi < setup.emax - 0.5 * d)])
        q1 = b["q"][:, g - 1:-(g - 1) or None, g - 1:-(g - 1) or None, g - 1:-(g - 1) or None]
        r_blk = cell_ratio(q1, d, is_open, setup).max()
        o = np.rint((lo - setup.emin) / d).astype(int)
        ref = r_fine[o[0]:o[0] + n[0], o[1]:o[1] + n[1], o[2]:o[2] + n[2]].max()
        worst = max(worst, abs(r_blk - ref) / ref)
    return worst


def wall_units(q: np.ndarray, d: np.ndarray, setup: Setup) -> tuple[float, float]:
    """Return the friction Reynolds number and the viscous length at the y walls (mean over x, z and both walls).

    The wall shear is the mean streamwise velocity of the first cell over its distance to the wall (no slip), the wall
    viscosity the law at the first cell's mean temperature.

    Parameters
    ----------
    q : np.ndarray
        Conservative variables (5, Nx, Ny, Nz).
    d : np.ndarray
        Spacings.
    setup : Setup
        Marker inputs.

    Returns
    -------
    tuple[float, float]
        Re_tau (on the half height) and delta_nu = nu_w / u_tau.
    """
    taus, nus = [], []
    for j in (0, q.shape[2] - 1):
        rho, vel, temp = auxiliary(q[:, :, j, :], setup.gas)
        mu = setup.gas.nu(rho, temp) * rho
        tau = (mu * np.abs(vel[0]) / (0.5 * d[1])).mean()
        taus.append(np.sqrt(tau / rho.mean()))
        nus.append(setup.gas.nu(rho, temp).mean())
    u_tau, nu_w = float(np.mean(taus)), float(np.mean(nus))
    half = 0.5 * (setup.emax[1] - setup.emin[1])
    return u_tau * half / nu_w, nu_w / u_tau


def parse_pairs(text: str) -> list[tuple[float, float]]:
    """Return the (refine_tol, derefine_tol) pairs of `R:D,R:D,...`.

    Parameters
    ----------
    text : str
        Pairs.

    Returns
    -------
    list[tuple[float, float]]
        Pairs.
    """
    return [(float(a), float(b)) for a, b in (p.split(":") for p in text.split(",") if p)]


def predict(rmax_f: np.ndarray, rmax_c: np.ndarray, refine: float, derefine: float) -> dict[str, int | bool]:
    """Return the grid after the first regrids from the uniform finest level with one coarser level allowed.

    Regrid 1: a sibling group (2 x 2 x 2 fine blocks) derefines iff all its 8 children vote derefine (r_max < derefine).
    Regrid 2: a derefined parent votes refine iff its r_max on the restricted field exceeds `refine` (it would then
    oscillate between the levels).

    Parameters
    ----------
    rmax_f : np.ndarray
        Per fine block r_max, (nbx, nby, nbz).
    rmax_c : np.ndarray
        Per coarse block r_max on the restricted field, (nbx / 2, nby / 2, nbz / 2).
    refine : float
        refine_tol.
    derefine : float
        derefine_tol.

    Returns
    -------
    dict[str, int | bool]
        Fine blocks voting derefine, groups derefined, groups that would refine back, fine/coarse blocks after, the
        admissibility of the pair.
    """
    nbx, nby, nbz = rmax_f.shape
    votes = rmax_f < derefine
    groups = votes.reshape(nbx // 2, 2, nby // 2, 2, nbz // 2, 2).all(axis=(1, 3, 5))
    back = groups & (rmax_c > refine)
    stay = groups & ~back
    return {"admissible": 0.0 < derefine < 0.5 * refine, "votes": int(votes.sum()), "groups": int(groups.sum()),
            "back": int(back.sum()), "fine": int(rmax_f.size - 8 * stay.sum()), "coarse": int(stay.sum())}


def report(setup: Setup, q: np.ndarray, d: np.ndarray, block: tuple[int, int, int], pairs: list[tuple[float, float]],
           stored: float | None) -> None:
    """Print the r distribution and the candidate pairs.

    Parameters
    ----------
    setup : Setup
        Marker inputs.
    q : np.ndarray
        Assembled conservative variables.
    d : np.ndarray
        Fine spacing.
    block : tuple[int, int, int]
        Cells per block.
    pairs : list[tuple[float, float]]
        Candidate (refine_tol, derefine_tol).
    stored : float or None
        Result of the stored-ghost cross-check, if run.
    """
    is_open = np.array([setup.periodic, setup.periodic])
    r_f = cell_ratio(pad(q, setup.periodic), d, is_open, setup)
    qc = restrict(q)
    r_c = cell_ratio(pad(qc, setup.periodic), 2.0 * d, is_open, setup)
    rmax_f, rmean_f, rmin_f = (block_reduce(r_f, block, h) for h in ("max", "mean", "min"))
    rmax_c = block_reduce(r_c, block)
    nb = rmax_f.shape
    re_tau, dnu = wall_units(q, d, setup)
    h = spacing(d, setup.active)
    print(f"grid {q.shape[1]} x {q.shape[2]} x {q.shape[3]} cells, d = {d}, h = {h:.5g}; blocks of {block}: "
          f"{nb[0]} x {nb[1]} x {nb[2]} = {rmax_f.size}; one level coarser {rmax_c.size} blocks")
    print(f"gas: gamma {setup.gas.gamma:.4g}, R {setup.gas.rgas:.4g}, mu0 {setup.gas.mu0:.6g}, T_ref {setup.gas.tref:.6g}, "
          f"omega {setup.gas.omega:.3g}, S {setup.gas.suth:.3g}")
    print(f"wall units: Re_tau = {re_tau:.1f}, delta_nu = {dnu:.5g}; dx+ {d[0] / dnu:.2f}, dy+ {d[1] / dnu:.2f}, "
          f"dz+ {d[2] / dnu:.2f}, h+ {h / dnu:.2f} (r -> h+ at the wall, eps+ = 1)")
    if stored is not None:
        print(f"stored-ghost per-block evaluation vs assembled: max relative difference of r_max {stored:.3e}")
    print(f"cell r: min {r_f.min():.3f}, median {np.median(r_f):.3f}, mean {r_f.mean():.3f}, max {r_f.max():.3f}")
    print(f"cell r one level coarser (restricted field, 2h): min {r_c.min():.3f}, median {np.median(r_c):.3f}, "
          f"max {r_c.max():.3f}")

    # by wall distance (both walls folded), cell rows
    ny = q.shape[2]
    print("\nr by wall distance (cell rows from the nearest wall, both walls folded; over x, z):")
    print(f"{'rows':>9} {'y_w':>15} {'y+':>15} {'r mean':>8} {'r max':>8} {'2h: mean':>9} {'2h: max':>8}")
    edges = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
    edges = [e for e in edges if e <= ny // 2]
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        rows = list(range(a, b)) + list(range(ny - b, ny - a))
        crow = sorted({j // 2 for j in rows})
        sel, selc = r_f[:, rows, :], r_c[:, crow, :]
        yw0, yw1 = a * d[1], b * d[1]
        print(f"{a:>4}-{b - 1:<4} {yw0:7.4f}-{yw1:<7.4f} {yw0 / dnu:7.1f}-{yw1 / dnu:<7.1f} {sel.mean():8.3f} "
              f"{sel.max():8.3f} {selc.mean():9.3f} {selc.max():8.3f}")

    # per block row along y
    print("\nper block, by block row along y (r over the block's cells; the marker votes on r_max):")
    print(f"{'row':>3} {'y range':>17} {'blocks':>6} {'r_max min':>9} {'r_max mean':>10} {'r_max max':>9} "
          f"{'r_mean':>7} {'r_min':>7}")
    for jb in range(nb[1]):
        y0 = setup.emin[1] + jb * block[1] * d[1]
        y1 = y0 + block[1] * d[1]
        m = rmax_f[:, jb, :]
        print(f"{jb:>3} {y0:8.4f}-{y1:<8.4f} {m.size:>6} {m.min():9.3f} {m.mean():10.3f} {m.max():9.3f} "
              f"{rmean_f[:, jb, :].mean():7.3f} {rmin_f[:, jb, :].min():7.3f}")
    print("\none level coarser (the parent of each sibling group, r_max on the restricted field), by row along y:")
    print(f"{'row':>3} {'y range':>17} {'blocks':>6} {'r_max min':>9} {'r_max max':>9} {'2 x max child r_max':>19}")
    groups_max = rmax_f.reshape(nb[0] // 2, 2, nb[1] // 2, 2, nb[2] // 2, 2).max(axis=(1, 3, 5))
    for jb in range(rmax_c.shape[1]):
        y0 = setup.emin[1] + jb * 2 * block[1] * d[1]
        y1 = y0 + 2 * block[1] * d[1]
        m = rmax_c[:, jb, :]
        print(f"{jb:>3} {y0:8.4f}-{y1:<8.4f} {m.size:>6} {m.min():9.3f} {m.max():9.3f} "
              f"{2.0 * groups_max[:, jb, :].max():19.3f}")
    print("\nper sibling group: it derefines iff derefine_tol > max child r_max; it stays coarse iff refine_tol > its "
          "coarse r_max")
    print(f"  max child r_max over groups: min {groups_max.min():.3f}, max {groups_max.max():.3f}; "
          f"coarse r_max: min {rmax_c.min():.3f}, max {rmax_c.max():.3f}")

    print("\ncandidate pairs (refine_tol R, derefine_tol D), from the uniform finest grid:")
    print(f"{'R':>6} {'D':>6} {'D<R/2':>5} {'fine votes D':>12} {'groups deref':>12} {'refine back':>11} "
          f"{'fine blocks':>11} {'coarse blocks':>13} {'cells':>12}")
    fine_cells = block[0] * block[1] * block[2]
    for refine, derefine in pairs:
        p = predict(rmax_f, rmax_c, refine, derefine)
        cells = (int(p["fine"]) + int(p["coarse"])) * fine_cells
        print(f"{refine:6.2f} {derefine:6.2f} {'yes' if p['admissible'] else 'NO':>5} {p['votes']:>12} "
              f"{p['groups']:>12} {p['back']:>11} {p['fine']:>11} {p['coarse']:>13} {cells:>12,}")


def main() -> None:
    """Command line."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder", type=Path, help="directory holding the checkpoint <base>-<it>-proc<rank>.h5")
    ap.add_argument("ini", type=Path, help="the run's input.ini")
    ap.add_argument("--iteration", type=int, default=None, help="checkpoint iteration (default: the last)")
    ap.add_argument("--block", default=None, help="cells per block i,j,k (default: the input's ni, nj, nk)")
    ap.add_argument("--pairs", default="4:1.9", help="candidate refine_tol:derefine_tol pairs, comma separated")
    ap.add_argument("--stored-ghosts", action="store_true",
                    help="cross-check the per-block evaluation on the checkpoint's stored ghosts")
    args = ap.parse_args()
    setup = load_setup(args.ini)
    blocks = read_blocks(checkpoint_files(args.folder, args.iteration), setup.ngc)
    q, d = assemble(blocks, setup)
    stored = None
    if args.stored_ghosts:
        stored = stored_ghost_check(blocks, setup, cell_ratio(pad(q, setup.periodic), d,
                                                              np.array([setup.periodic, setup.periodic]), setup))
    del blocks
    block = tuple(int(c) for c in args.block.split(",")) if args.block else setup.cells
    if setup.tolerances:
        print(f"input marker: refine_tol {setup.tolerances[0]}, derefine_tol {setup.tolerances[1]}")
    report(setup, q, d, block, parse_pairs(args.pairs), stored)  # type: ignore[arg-type]


if __name__ == "__main__":
    main()
