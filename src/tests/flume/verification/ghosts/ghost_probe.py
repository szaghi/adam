#!/usr/bin/env python3
"""Check every ghost cell of a FLUME step-0 checkpoint against the linear field it must hold (issue #65, P0).

The run uses `[initial_conditions] type = linear`, so every conservative variable is linear in space,
`q(x) = q_1 (1 + g . x)`, and takes one step of negligible size (`CFL = 1e-30`, `it_max = 0`, `time_max = 0`: the
field moves by ~1e-30 relative). The step-0 checkpoint cannot serve: it is written while each realm initialises, before
the forest connects the realms, so its seam ghosts are unfilled. The checkpoint after the step is written by the
post-step save, which fills the ghosts in the order of a Runge-Kutta stage (inter-realm seams, then the intra-realm
exchange and the boundary conditions) and writes the fields with their ghosts. The same-level copy, the 2:1 restriction and the tricubic
coarse->fine fill are exact on a linear field, so every ghost cell has a known value:

- inside the realm, or across an inter-realm seam, a ghost holds `f(c)` at its own centre `c`;
- along an axis where it lies beyond a physical face of its realm, the boundary condition of that face acts on the
  value the ghost would hold without the face: `wall-inviscid` mirrors the coordinate about the face and negates the
  wall-normal momentum (and the wall-normal field on MHD); `wall-noslip` and `wall-isothermal` mirror it and reflect
  the velocity about the wall velocity (the isothermal wall also sets the ghost temperature `2 T_w - T`);
  `extrapolation` copies the first interior cell along the normal; `periodic` wraps the coordinate; `inflow` holds
  the inflow state. Beyond two physical faces the transforms compose in the backends' order (inflow first, else the
  first face in axis order acts on a donor valued the same way).

This is what a stencil reading the ghost must see. The directional WENO stencils read the face ghosts (outside the
block along one axis); a cross derivative reads the edge ghosts (outside along two axes), such as the tangential
derivatives of the viscous stress at a face (M4, issue #65). Corners (outside along three axes) are reported, not
asserted: no M4 stencil reads them. Neither are the `solid` ghosts, those beyond seams into no realm (the re-entrant
corner of an L-shaped forest: inside the step of the Woodward-Colella tunnel), which hold no fluid value. A null axis
is not an axis here.

The report groups the ghosts by class (face, edge, corner) and by what lies beyond the block along the axes where the
ghost is outside it: `block` (another block of the same realm), `seam` (an inter-realm seam), or the boundary condition
of a physical face (`wall-inviscid`, `extrapolation`, `inflow`, `periodic`). Per group: cells checked and the largest error relative to the field scale `max |q_1| (1 + |g| . L)`.

Usage: ghost_probe.py <work-dir> [--tol 1e-12] [--worst N]
"""

from __future__ import annotations

import argparse
import configparser
import sys
from pathlib import Path

import h5py
import numpy as np

AXES = "xyz"
FACES = {"-x": (0, 0), "+x": (0, 1), "-y": (1, 0), "+y": (1, 1), "-z": (2, 0), "+z": (2, 1)}
MOMENTUM = ("ru", "rv", "rw")
FIELD = ("bx", "by", "bz")


def read_ini(path: Path) -> configparser.ConfigParser:
    """Return an INI file, keys case-sensitive, inline ';' comments stripped."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), strict=False, interpolation=None)
    ini.optionxform = str  # type: ignore[assignment,method-assign]
    ini.read(path)
    return ini


def is_true(value: str) -> bool:
    """Return the value of a Fortran logical."""
    return value.strip().lower() in (".true.", "true", "t")


class Realm:
    """Geometry, boundary conditions, seam faces and step-0 blocks of one realm."""

    def __init__(self, ini_path: Path, work: Path, seams: set[tuple[int, int]]) -> None:
        ini = read_ini(ini_path)
        grid = ini["grid"]
        self.name = ini["IO"]["output_basename"].strip()
        self.ngc = int(grid["ngc"])
        self.lo = np.array([float(grid[f"emin_{a}"]) for a in AXES])
        self.hi = np.array([float(grid[f"emax_{a}"]) for a in AXES])
        self.null = np.array([is_true(grid[f"null_{a}"]) for a in AXES])
        self.seams = seams
        if "gamma" in ini["physics"]:
            self.gamma, self.gas = float(ini["physics"]["gamma"]), 1.0
        else:
            cp, cv = float(ini["physics"]["cp"]), float(ini["physics"]["cv"])
            self.gamma, self.gas = cp / cv, cp - cv
        eglm = ini.has_section("mhd") and ini["mhd"].get("divergence_control", "").strip() == "eglm"
        self.psi_energy = 1.0 if eglm else 0.0
        self.bc: dict[tuple[int, int], str] = {}
        self.inflow: dict[tuple[int, int], dict[str, float]] = {}
        self.wall_velocity: dict[tuple[int, int], np.ndarray] = {}
        self.wall_temperature: dict[tuple[int, int], float] = {}
        for d, a in enumerate(AXES):
            for s, side in enumerate(("min", "max")):
                section = ini[f"bc_{a}_{side}"]
                self.bc[(d, s)] = section["type"].strip()
                if self.bc[(d, s)] == "inflow":
                    self.inflow[(d, s)] = conservative({k: float(v) for k, v in section.items() if k != "type"},
                                                       self.gamma)
                self.wall_velocity[(d, s)] = np.array([float(section.get(f"wall_{c}", "0.0")) for c in "uvw"])
                self.wall_temperature[(d, s)] = float(section.get("wall_temperature", "0.0"))
        ic = ini["initial_conditions"]
        if ic["type"].strip() != "linear":
            sys.exit(f"ghost_probe: {ini_path} has [initial_conditions] type = {ic['type']}, not linear")
        self.gradient = np.array([float(ic[f"gradient_{a}"]) for a in AXES])
        self.blocks = self.load(work)

    def load(self, work: Path) -> list[dict]:
        """Return the blocks of the last checkpoint: interior origin and upper corner, spacing, fields with ghosts."""
        blocks = []
        files = [p for p in work.glob(f"{self.name}-*-proc*.h5") if "-restart-" not in p.name]
        if not files:
            sys.exit(f"ghost_probe: no checkpoint of {self.name} in {work}")
        last = max(int(p.name.split("-")[-2]) for p in files)
        if last == 0:
            sys.exit(f"ghost_probe: {self.name} holds only the step-0 checkpoint (seams unfilled): run one step")
        for path in sorted(p for p in files if int(p.name.split("-")[-2]) == last):
            with h5py.File(path, "r") as h5:
                for blk in sorted({k.rsplit("-", 1)[0] for k in h5}):
                    names = sorted(k.rsplit("-", 1)[1] for k in h5 if k.rsplit("-", 1)[0] == blk)
                    names = [v for v in names if h5[f"{blk}-{v}"].ndim == 3]  # the fields, not the block metadata
                    dx = h5[f"{blk}-dxdydz"][()][::-1]
                    q = np.stack([h5[f"{blk}-{v}"][()].transpose(2, 1, 0) for v in names])
                    lo = h5[f"{blk}-origin"][()][::-1] + self.ngc * dx
                    blocks.append({"lo": lo, "hi": lo + (np.array(q.shape[1:]) - 2 * self.ngc) * dx, "dx": dx,
                                   "q": q, "names": names})
        return blocks


def conservative(prim: dict[str, float], gamma: float) -> dict[str, float]:
    """Return the conservative state of an inflow (Euler or MHD primitive keys), keyed by lower-case variable name."""
    r, u, v, w, p = (prim[k] for k in ("r", "u", "v", "w", "p"))
    b = np.array([prim.get(k, 0.0) for k in ("bx", "by", "bz")])
    state = {"r": r, "ru": r * u, "rv": r * v, "rw": r * w,
             "re": p / (gamma - 1.0) + 0.5 * r * (u * u + v * v + w * w) + 0.5 * float(b @ b)}
    state.update({"bx": b[0], "by": b[1], "bz": b[2], "psi": 0.0})
    return state


def load_realms(work: Path) -> list[Realm]:
    """Return the realms of the run: those of a forest manifest, or the single realm INI of the directory."""
    for path in sorted(work.glob("*.ini")):
        ini = read_ini(path)
        if not ini.has_section("forest"):
            continue
        n = int(ini["forest"]["realms_number"])
        seams: list[set[tuple[int, int]]] = [set() for _ in range(n)]
        for section in ini.sections():
            if section.startswith("forest.topology.face_"):
                face = ini[section]
                seams[int(face["realm_a"]) - 1].add(FACES[face["face_a"].strip()])
                seams[int(face["realm_b"]) - 1].add(FACES[face["face_b"].strip()])
        return [Realm(work / ini[f"realm.{r + 1}"]["ini"], work, seams[r]) for r in range(n)]
    singles = [p for p in sorted(work.glob("*.ini")) if read_ini(p).has_section("grid")]
    if len(singles) != 1:
        sys.exit(f"ghost_probe: {work} holds no forest manifest and {len(singles)} realm INIs")
    return [Realm(singles[0], work, set())]


def noslip_ghost(q: np.ndarray, lower: list[str], d: int, realm: Realm, wall_velocity: np.ndarray,
                 wall_temperature: float | None) -> np.ndarray:
    """Return the no-slip wall ghost of the mirrored state `q` (the rule of `wall_noslip_ghost`, issue #65 P1)."""
    i_r, i_e = lower.index("r"), lower.index("re")
    i_m = [lower.index(n) for n in MOMENTUM]
    i_b = [lower.index(n) for n in FIELD if n in lower]
    u = q[i_m] / q[i_r]
    ug = 2.0 * wall_velocity - u
    e_mag = 0.5 * float(np.sum(q[i_b] ** 2)) + (0.5 * realm.psi_energy * q[lower.index("psi")] ** 2
                                                  if "psi" in lower else 0.0)
    p = (realm.gamma - 1.0) * (q[i_e] - 0.5 * q[i_r] * float(u @ u) - e_mag)
    rho = q[i_r]
    if wall_temperature is not None:
        rho = p / (realm.gas * (2.0 * wall_temperature - p / (q[i_r] * realm.gas)))
    g = q.copy()
    if i_b:
        g[i_b[d]] = -q[i_b[d]]
    g[i_r] = rho
    g[i_m] = rho * ug
    g[i_e] = p / (realm.gamma - 1.0) + 0.5 * rho * float(ug @ ug) + e_mag
    return g


def expected_value(realm: Realm, centre: np.ndarray, blk: dict, q1: np.ndarray, sign_wall: np.ndarray,
                   lower: list[str], realms: list[Realm]) -> tuple[np.ndarray | None, set[str]]:
    """Return the value a ghost centred at `centre` must hold (None outside the forest) and what lies beyond its block.

    The rule of the backends: among the physical realm faces the ghost lies beyond, the first inflow face, else the
    first in axis order, acts on the donor mirrored (walls), clamped (extrapolation) or wrapped (periodic) along its
    axis; the donor is valued the same way, recursively, so the transforms compose in the backends' order.
    """
    beyond: set[str] = set()
    for d in range(3):
        if realm.null[d] or blk["lo"][d] <= centre[d] <= blk["hi"][d]:
            continue
        s = 0 if centre[d] < blk["lo"][d] else 1
        if (s == 0 and centre[d] > realm.lo[d]) or (s == 1 and centre[d] < realm.hi[d]):
            beyond.add("block")
        elif (d, s) in realm.seams:
            beyond.add("seam")
        else:
            beyond.add(realm.bc[(d, s)])

    def value(c: np.ndarray) -> np.ndarray | None:
        faces = [(d, 0 if c[d] < realm.lo[d] else 1) for d in range(3)
                 if not realm.null[d] and (c[d] < realm.lo[d] or c[d] > realm.hi[d])]
        physical = [f for f in faces if f not in realm.seams]
        if not physical:
            if not any(((r.lo <= c) | r.null).all() and ((c <= r.hi) | r.null).all() for r in realms):
                return None  # beyond seams into no realm: a re-entrant corner of the forest, no fluid value
            return q1 * (1.0 + realm.gradient @ c)
        inflow = [f for f in physical if realm.bc[f] == "inflow"]
        d, s = (inflow or physical)[0]
        bc = realm.bc[(d, s)]
        face = realm.lo[d] if s == 0 else realm.hi[d]
        if bc == "inflow":
            return np.array([realm.inflow[(d, s)][v] for v in lower])
        c = c.copy()
        if bc in ("wall-inviscid", "wall-noslip", "wall-isothermal"):
            c[d] = 2.0 * face - c[d]
        elif bc == "extrapolation":
            c[d] = face + (0.5 if s == 0 else -0.5) * blk["dx"][d]
        elif bc == "periodic":
            c[d] += (realm.hi[d] - realm.lo[d]) * (1.0 if s == 0 else -1.0)
        else:
            sys.exit(f"ghost_probe: boundary condition '{bc}' is not modelled")
        donor = value(c)
        if donor is None or bc in ("extrapolation", "periodic"):
            return donor
        if bc == "wall-inviscid":
            return sign_wall[d] * donor
        return noslip_ghost(donor, lower, d, realm, realm.wall_velocity[(d, s)],
                            realm.wall_temperature[(d, s)] if bc == "wall-isothermal" else None)

    return value(centre.copy()), beyond


def probe(work: Path) -> dict[tuple[str, str], list]:
    """Return the ghost cells of the last checkpoint of `work`, grouped by (class, what lies beyond the block).

    Each row is (error relative to the field scale, realm, block lower corner, cell index, cell centre).
    """
    realms = load_realms(work)
    lower = [v.lower() for v in realms[0].blocks[0]["names"]]
    gradient = realms[0].gradient
    if any((r.gradient != gradient).any() for r in realms):
        sys.exit("ghost_probe: the realms must share the gradient")
    # q_1 from one interior cell per block, q = q_1 (1 + g . x)
    samples = []
    for realm in realms:
        for blk in realm.blocks:
            g = realm.ngc
            samples.append(blk["q"][:, g, g, g] / (1.0 + gradient @ (blk["lo"] + 0.5 * blk["dx"])))
    q1 = np.median(np.array(samples), axis=0)
    extent = np.max([np.maximum(np.abs(r.lo), np.abs(r.hi)) for r in realms], axis=0)
    scale = np.max(np.abs(q1)) * (1.0 + np.abs(gradient) @ extent)
    sign_wall = np.ones((3, len(lower)))
    for d in range(3):
        for names in (MOMENTUM, FIELD):
            if names[d] in lower:
                sign_wall[d, lower.index(names[d])] = -1.0

    groups: dict[tuple[str, str], list] = {}
    for realm in realms:
        g = realm.ngc
        for blk in realm.blocks:
            q = blk["q"]
            n = np.array(q.shape[1:]) - 2 * g
            for i, j, k in np.ndindex(*q.shape[1:]):
                idx = np.array([i, j, k]) - g
                if any(realm.null[d] and (idx[d] < 0 or idx[d] >= n[d]) for d in range(3)):
                    continue  # a ghost along a null axis: no stencil differences that axis
                outside = [bool((idx[d] < 0 or idx[d] >= n[d]) and not realm.null[d]) for d in range(3)]
                if not any(outside):
                    continue
                centre = blk["lo"] + (idx + 0.5) * blk["dx"]
                value, beyond = expected_value(realm, centre, blk, q1, sign_wall, lower, realms)
                cls = ("face", "edge", "corner")[sum(outside) - 1]
                if value is None:
                    cls, err = "solid", 0.0
                else:
                    err = float(np.max(np.abs(q[:, i, j, k] - value)) / scale)
                    if not np.isfinite(err):
                        err = np.inf  # a poisoned ghost the fill missed or read before its donor was filled
                groups.setdefault((cls, "+".join(sorted(beyond))), []).append(
                    (err, realm.name, tuple(np.round(blk["lo"], 6)), tuple(idx + 1), tuple(np.round(centre, 6))))
    return groups


def main() -> int:
    """Run the probe, print the report, return the exit status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="work directory of the one-step run")
    parser.add_argument("--tol", type=float, default=1.0e-12, help="error bound relative to the field scale")
    parser.add_argument("--worst", type=int, default=4, help="worst cells listed per failing group")
    args = parser.parse_args()

    groups = probe(args.work)
    status = 0
    print(f"   {'class':6s} {'beyond the block':24s} {'cells':>8s} {'max error':>10s}")
    for (cls, where), rows in sorted(groups.items()):
        worst = max(rows)[0]
        bad = cls in ("face", "edge") and worst > args.tol
        status |= int(bad)
        verdict = "FAIL" if bad else "PASS" if cls in ("face", "edge") else "info"
        print(f"   {cls:6s} {where:24s} {len(rows):8d} {worst:10.1e}  {verdict}")
        if bad:
            for err, name, lo, ijk, c in sorted(rows, reverse=True)[: args.worst]:
                print(f"        {err:9.2e}  {name} block at {lo}, cell {ijk}, centre {c}")
    return status


if __name__ == "__main__":
    sys.exit(main())
