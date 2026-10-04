#!/usr/bin/env python3
"""Draw the figures of the forest cookbook (docs/guide/forest-cookbook.md) as SVG.

Why this script exists: the multi-realm machinery (inter-realm seams, 2:1 jumps, misaligned blocks, the flux register,
the seam cadence, the rank partition) is geometric, and its documentation needs pictures that stay consistent with each
other and can be regenerated when the machinery changes. The figures are schematic but to scale: cells, blocks and
ghost layers are drawn with the counts they stand for. Standard library only (no matplotlib): each figure is a few
primitives (rectangles, grids, arrows, labels) written as SVG text, on a light card that reads in both VitePress themes.

Usage:
    make_forest_figures.py [--out DIR]     (default: docs/public/forest)
"""

from __future__ import annotations

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FONT = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif"
MONO = "'JetBrains Mono', 'Fira Code', Menlo, Consolas, monospace"
INK = "#1F2933"
MUTED = "#52606D"
GRID = "#9AA5B1"
CARD = "#FFFFFF"
CARD_EDGE = "#D9E2EC"
A_EDGE, A_FILL, A_DARK = "#2F6DB5", "#DCE8F7", "#9BBCE3"
B_EDGE, B_FILL, B_DARK = "#D9731F", "#FCE6D3", "#F2B88A"
SEAM = "#C0392B"
REG_EDGE, REG_FILL = "#6B3FA0", "#E7DDF4"
GHOST = "#7B8794"
OK = "#2E8540"


class Svg:
    """A minimal SVG canvas: primitives append elements; `save` wraps them in a card."""

    def __init__(self, width: int, height: int, title: str) -> None:
        self.w, self.h, self.title = width, height, title
        self.items: list[str] = []

    def rect(self, x: float, y: float, w: float, h: float, fill: str = "none", stroke: str = INK, sw: float = 1.0,
             rx: float = 0.0, dash: str = "", opacity: float = 1.0) -> None:
        d = f' stroke-dasharray="{dash}"' if dash else ""
        o = f' fill-opacity="{opacity}"' if opacity < 1.0 else ""
        self.items.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{rx}" fill="{fill}"{o} '
                          f'stroke="{stroke}" stroke-width="{sw}"{d}/>')

    def line(self, x1: float, y1: float, x2: float, y2: float, stroke: str = INK, sw: float = 1.0, dash: str = "",
             arrow: str = "") -> None:
        d = f' stroke-dasharray="{dash}"' if dash else ""
        m = f' marker-end="url(#{arrow})"' if arrow else ""
        self.items.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
                          f'stroke-width="{sw}"{d}{m}/>')

    def path(self, d: str, stroke: str = INK, sw: float = 1.2, fill: str = "none", arrow: str = "",
             dash: str = "") -> None:
        m = f' marker-end="url(#{arrow})"' if arrow else ""
        s = f' stroke-dasharray="{dash}"' if dash else ""
        self.items.append(f'<path d="{d}" stroke="{stroke}" stroke-width="{sw}" fill="{fill}"{m}{s}/>')

    def text(self, x: float, y: float, s: str, size: float = 13, color: str = INK, anchor: str = "start",
             weight: str = "normal", mono: bool = False) -> None:
        f = f' font-family="{MONO}"' if mono else ""
        self.items.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" fill="{color}" text-anchor="{anchor}" '
                          f'font-weight="{weight}" xml:space="preserve"{f}>{s}</text>')

    def grid(self, x: float, y: float, nx: int, ny: int, d: float, stroke: str = GRID, sw: float = 0.6,
             dy: float | None = None) -> None:
        dy = d if dy is None else dy
        for i in range(1, nx):
            self.line(x + i * d, y, x + i * d, y + ny * dy, stroke, sw)
        for j in range(1, ny):
            self.line(x, y + j * dy, x + nx * d, y + j * dy, stroke, sw)

    def cells(self, x: float, y: float, nx: int, ny: int, d: float, fill: str, edge: str, sw: float = 1.6,
              dy: float | None = None, dash: str = "", opacity: float = 1.0) -> None:
        dy = d if dy is None else dy
        self.rect(x, y, nx * d, ny * dy, fill=fill, stroke="none", opacity=opacity)
        self.grid(x, y, nx, ny, d, dy=dy)
        self.rect(x, y, nx * d, ny * dy, stroke=edge, sw=sw, dash=dash)

    def ghost(self, x: float, y: float, w: float, h: float) -> None:
        self.rect(x, y, w, h, fill="url(#hatch)", stroke=GHOST, sw=1.0, dash="4 3")

    def save(self, path: Path) -> None:
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {self.w} {self.h}" width="{self.w}" '
                f'height="{self.h}" font-family="{FONT}" role="img" aria-label="{self.title}">\n'
                f'<title>{self.title}</title>\n<defs>\n'
                '<pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" '
                'patternTransform="rotate(45)"><rect width="6" height="6" fill="#F5F7FA" fill-opacity="0.7"/>'
                f'<line x1="0" y1="0" x2="0" y2="6" stroke="{GHOST}" stroke-width="1.1" stroke-opacity="0.55"/>'
                '</pattern>\n')
        for name, color in (("arr", INK), ("arrA", A_EDGE), ("arrB", B_EDGE), ("arrS", SEAM), ("arrR", REG_EDGE),
                            ("arrG", OK)):
            head += (f'<marker id="{name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
                     f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{color}"/></marker>\n')
        head += '</defs>\n'
        card = (f'<rect x="1" y="1" width="{self.w - 2}" height="{self.h - 2}" rx="12" fill="{CARD}" '
                f'stroke="{CARD_EDGE}" stroke-width="1.5"/>\n')
        path.write_text(head + card + "\n".join(self.items) + "\n</svg>\n")


def fig_anatomy(out: Path) -> None:
    """Two realms glued by a mirror seam: realms, blocks, cells, seam, ghost slab, manifest faces."""
    s = Svg(860, 420, "Anatomy of a two-realm forest")
    s.text(430, 34, "A forest: two realms glued at an inter-realm seam", 17, anchor="middle", weight="bold")
    ox, oy, bw = 70, 80, 70
    s.cells(ox, oy, 4, 4, bw, A_FILL, A_EDGE, sw=2.4)
    s.cells(ox + 4 * bw, oy, 4, 4, bw, B_FILL, B_EDGE, sw=2.4)
    for k in range(1, 4):  # block edges stronger than the grid
        s.line(ox + k * bw, oy, ox + k * bw, oy + 4 * bw, A_EDGE, 1.4)
        s.line(ox, oy + k * bw, ox + 4 * bw, oy + k * bw, A_EDGE, 1.4)
        s.line(ox + (4 + k) * bw, oy, ox + (4 + k) * bw, oy + 4 * bw, B_EDGE, 1.4)
        s.line(ox + 4 * bw, oy + k * bw, ox + 8 * bw, oy + k * bw, B_EDGE, 1.4)
    g = bw / 6
    s.grid(ox + 3 * bw, oy, 6, 6, g, stroke=A_DARK, sw=0.7)
    s.grid(ox + 6 * bw, oy + 3 * bw, 6, 6, g, stroke=B_DARK, sw=0.7)
    s.ghost(ox + 4 * bw, oy, 3 * g, 4 * bw)
    sx = ox + 4 * bw
    s.line(sx, oy - 14, sx, oy + 4 * bw + 14, SEAM, 3.4)
    s.text(sx, oy - 20, "seam", 13, SEAM, "middle", "bold")
    s.text(ox + 2 * bw, oy + 4 * bw + 28, "realm 1  (realm_a, face_a = +x)", 14, A_EDGE, "middle", "bold")
    s.text(ox + 6 * bw, oy + 4 * bw + 28, "realm 2  (realm_b, face_b = -x)", 14, B_EDGE, "middle", "bold")
    s.text(ox + 2 * bw, oy + 4 * bw + 47, "its own INI: grid, numerics, physics", 12, MUTED, "middle")
    s.text(ox + 6 * bw, oy + 4 * bw + 47, "its own INI, its own rank partition", 12, MUTED, "middle")
    s.line(ox + 3.5 * bw, oy + 0.5 * bw, ox + 2.6 * bw, oy - 26, MUTED, 1, arrow="arr")
    s.text(ox + 2.6 * bw - 4, oy - 30, "block: ni x nj x nk cells (one shown)", 12, MUTED, "end")
    tx = ox + 8 * bw + 18
    for k, t in enumerate(["ghost slab of realm 1:", "ngc cells beyond its", "seam face, filled from",
                           "realm 2 interior cells", "(realm 2 has its own,", "over realm 1)"]):
        s.text(tx, oy + 1.0 * bw + 17 * k, t, 12, MUTED)
    s.line(tx - 4, oy + 1.0 * bw + 30, ox + 4 * bw + 3 * g + 3, oy + 1.4 * bw, MUTED, 1, arrow="arr")
    s.save(out / "anatomy.svg")


def fig_mirror_aligned(out: Path) -> None:
    """Zoom on a mirror seam with lined-up blocks: the ghosts of each side copy the peer cells (COPY rows)."""
    s = Svg(860, 400, "Mirror seam, lined-up blocks")
    s.text(430, 34, "coupling = mirror: same cell size, a seam ghost is a copy of the peer cell", 17,
           anchor="middle", weight="bold")
    d, ox, oy = 34, 70, 80
    s.cells(ox, oy, 6, 6, d, A_FILL, A_EDGE, 2.2)
    s.ghost(ox + 6 * d, oy, 3 * d, 6 * d)
    s.line(ox + 6 * d, oy - 10, ox + 6 * d, oy + 6 * d + 10, SEAM, 3.2)
    s.text(ox + 4.5 * d, oy + 6 * d + 26, "realm 1 block + its 3 ghost columns", 13, A_EDGE, "middle", "bold")
    bx = 520
    s.cells(bx, oy, 6, 6, d, B_FILL, B_EDGE, 2.2)
    s.rect(bx, oy, 3 * d, 6 * d, fill=B_DARK, stroke="none", opacity=0.45)
    s.grid(bx, oy, 3, 6, d)
    s.line(bx, oy - 10, bx, oy + 6 * d + 10, SEAM, 3.2)
    s.text(bx + 3 * d, oy + 6 * d + 26, "realm 2 block: its first 3 interior columns", 13, B_EDGE, "middle", "bold")
    for j in (1, 3, 5):
        y = oy + (j + 0.5) * d
        s.path(f"M{bx + 1.5 * d:.1f},{y:.1f} C{bx - 50:.1f},{y - 36:.1f} {ox + 9 * d + 50:.1f},{y - 36:.1f} "
               f"{ox + 7.5 * d:.1f},{y:.1f}", B_EDGE, 1.5, arrow="arrB")
    s.text(430, oy + 6 * d + 58, "the cell centres coincide across the seam: every seam ghost is a COPY row (copied "
           "directly, or packed and sent", 12, MUTED, "middle")
    s.text(430, oy + 6 * d + 76, "when the cell lives on another rank); realm 2 ghosts are filled the same way from "
           "realm 1 interior cells", 12, MUTED, "middle")
    s.save(out / "mirror-aligned.svg")


def fig_mirror_misaligned(out: Path) -> None:
    """Face-on view of a mirror seam whose blocks do not line up (issue #51): overlaps of block faces."""
    s = Svg(880, 510, "Mirror seam, blocks not lined up")
    s.text(440, 34, "Blocks that do not line up (issue #51): the seam face seen along its normal", 17,
           anchor="middle", weight="bold")
    d, oy = 18, 90
    ax = 50
    s.cells(ax, oy, 16, 16, d, REG_FILL, REG_EDGE, 2.0)
    for k in range(1, 4):
        s.line(ax + 4 * k * d, oy, ax + 4 * k * d, oy + 16 * d, REG_EDGE, 1.8)
        s.line(ax, oy + 4 * k * d, ax + 16 * d, oy + 4 * k * d, REG_EDGE, 1.8)
    n = 0
    for j in range(4):
        for i in range(4):
            n += 1
            s.text(ax + (4 * i + 2) * d, oy + (4 * j + 2) * d + 5, str(n), 13, REG_EDGE, "middle", "bold")
    s.text(ax + 8 * d, oy - 14, "realm 1 (realm_a): 16 register faces", 13, REG_EDGE, "middle", "bold")
    s.text(ax + 8 * d, oy + 16 * d + 22, "one register face per realm_a seam block", 12, MUTED, "middle")
    bx = 540
    s.cells(bx, oy, 16, 16, d, B_FILL, B_EDGE, 2.0)
    s.line(bx + 8 * d, oy, bx + 8 * d, oy + 16 * d, B_EDGE, 1.8)
    s.line(bx, oy + 8 * d, bx + 16 * d, oy + 8 * d, B_EDGE, 1.8)
    s.rect(bx, oy, 8 * d, 8 * d, fill=B_DARK, stroke="none", opacity=0.45)
    s.text(bx + 8 * d, oy - 14, "realm 2 (realm_b): 4 block faces", 13, B_EDGE, "middle", "bold")
    s.text(bx + 8 * d, oy + 16 * d + 22, "same cells, blocks twice as large", 12, MUTED, "middle")
    for k, (i, j) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1))):
        s.rect(bx + 4 * i * d + 3, oy + 4 * j * d + 3, 4 * d - 6, 4 * d - 6, stroke=REG_EDGE, sw=1.6, dash="5 3")
        s.text(bx + (4 * i + 2) * d, oy + (4 * j + 2) * d + 5, ("1", "2", "5", "6")[k], 13, REG_EDGE, "middle", "bold")
    s.path(f"M{bx - 10:.1f},{oy + 4 * d:.1f} C{bx - 80:.1f},{oy + 2 * d:.1f} {ax + 16 * d + 60:.1f},"
           f"{oy + 2 * d:.1f} {ax + 16 * d + 8:.1f},{oy + 2 * d:.1f}", REG_EDGE, 1.6, arrow="arrR")
    s.text(440, oy + 2 * d - 26, "scattered over 4 register faces", 12, REG_EDGE, "middle", "bold")
    rows = [("overlap rows of the shaded realm-2 block face (maps%seam_overlap), columns", INK),
            ("[cursor, register offset (inner, outer), block offset (inner, outer), extent (inner, outer), "
             "register inner count]", MUTED),
            ("[1, 0 0, 0 0, 4 4, 4]    [2, 0 0, 4 0, 4 4, 4]    [5, 0 0, 0 4, 4 4, 4]    [6, 0 0, 4 4, 4 4, 4]", MUTED)]
    for k, (r, c) in enumerate(rows):
        s.text(440, oy + 16 * d + 56 + 18 * k, r, 12, c, "middle")
    s.save(out / "mirror-misaligned.svg")


def fig_refined(out: Path) -> None:
    """A 2:1 refined seam (issue #52): interpolate rows (fine ghosts), restrict rows (coarse ghosts)."""
    s = Svg(880, 480, "Refined 2:1 seam")
    s.text(440, 34, "coupling = refined: a 2:1 jump; the owner of the cells evaluates the ghost values", 17,
           anchor="middle", weight="bold")
    D, oy, cx = 44, 80, 60
    s.cells(cx, oy, 6, 6, D, A_FILL, A_EDGE, 2.2)
    sx = cx + 6 * D
    s.cells(sx, oy, 8, 12, D / 2, B_FILL, B_EDGE, 2.2)
    s.line(sx, oy - 12, sx, oy + 6 * D + 12, SEAM, 3.4)
    s.text(cx + 3 * D, oy + 6 * D + 26, "coarse realm (dx)", 13, A_EDGE, "middle", "bold")
    s.text(sx + 2 * D, oy + 6 * D + 26, "fine realm (dx/2)", 13, B_EDGE, "middle", "bold")
    s.rect(cx + 2 * D, oy + 1 * D, 4 * D, 4 * D, fill=A_DARK, stroke=A_EDGE, sw=1.6, opacity=0.55)
    s.grid(cx + 2 * D, oy + 1 * D, 4, 4, D)
    s.rect(sx - D / 2, oy + 2.5 * D, D / 2, D / 2, fill="#FFFFFF", stroke=B_EDGE, sw=2.4)
    s.rect(cx + 2 * D + 6, oy + 1 * D + 6, 168, 20, fill=CARD, stroke=A_EDGE, sw=1.0, rx=4)
    s.text(cx + 2 * D + 90, oy + 1 * D + 21, "INTERPOLATE row: footprint", 12, A_EDGE, "middle", "bold")
    s.text(380, oy + 6 * D + 56, "INTERPOLATE: a fine ghost (white) = interpolant over a 4x4x4 coarse footprint "
           "(tricubic), 3x3x3 (restriction-", 12, MUTED, "middle")
    s.text(380, oy + 6 * D + 72, "compatible) or the anchor cell (injection), the footprint shifted to stay inside "
           "the coarse block interior", 12, MUTED, "middle")
    cgx, cgy = sx, oy + 4 * D
    s.rect(cgx, cgy, D, D, fill=B_DARK, stroke="none", opacity=0.6)
    s.grid(cgx, cgy, 2, 2, D / 2, stroke=B_EDGE, sw=1.0)
    s.rect(cgx, cgy, D, D, stroke=A_EDGE, sw=2.6)
    s.text(sx + 3.6 * D, oy + 5.75 * D, "RESTRICT row", 12, B_EDGE, "middle", "bold")
    s.line(sx + 2.8 * D, oy + 5.6 * D, cgx + D + 4, cgy + D * 0.8, B_EDGE, 1.3, arrow="arrB")
    s.text(380, oy + 6 * D + 96, "RESTRICT: a coarse ghost (blue frame) = mean of the 2x2x2 fine cells under it "
           "(summed, then times 0.125)", 12, MUTED, "middle")
    rx = 700
    s.text(rx + 50, oy + 4, "register face", 13, REG_EDGE, "middle", "bold")
    s.text(rx + 50, oy + 20, "one per coarse seam block", 11, MUTED, "middle")
    s.cells(rx + 10, oy + 36, 4, 4, 20, REG_FILL, REG_EDGE, 2.0)
    for (i, j) in ((0, 0), (1, 0), (0, 1), (1, 1)):
        s.rect(rx + 10 + 40 * i + 2, oy + 36 + 40 * j + 2, 36, 36, stroke=B_EDGE, sw=1.2, dash="4 2")
        s.text(rx + 30 + 40 * i, oy + 61 + 40 * j, f"q{i}{j}", 11, B_EDGE, "middle")
    for k, t in enumerate(["4 fine blocks, each 2:1-", "restricted into its", "quadrant (qi, qj):",
                           "SEAM_KIND_INTER_", "REALM_REFINED"]):
        s.text(rx + 50, oy + 142 + 16 * k, t, 11, MUTED, "middle")
    s.save(out / "refined.svg")


def fig_intra_amr(out: Path) -> None:
    """A single realm refined on a box: intra-realm 2:1 faces."""
    s = Svg(860, 420, "Intra-realm AMR 2:1 faces")
    s.text(430, 34, "Intra-realm AMR: one realm, a refined box, 2:1 faces inside the tree", 17, anchor="middle",
           weight="bold")
    b, ox, oy = 64, 60, 80
    s.cells(ox, oy, 4, 4, b, A_FILL, A_EDGE, 2.4)
    s.rect(ox + 3 * b, oy, b, 4 * b, fill=A_DARK, stroke="none", opacity=0.6)
    s.grid(ox + 3 * b, oy, 2, 8, b / 2, stroke=A_EDGE, sw=1.2)
    s.line(ox + 3 * b, oy, ox + 3 * b, oy + 4 * b, REG_EDGE, 4.0)
    s.text(ox + 2 * b, oy + 4 * b + 26, "4x4 blocks (iu_ref_levels = 2), x > 0.75 one level finer", 13, A_EDGE,
           "middle", "bold")
    s.text(ox + 3 * b, oy - 12, "2:1 faces: SEAM_KIND_INTRA_REALM_AMR", 12, REG_EDGE, "middle", "bold")
    tx = 400
    lines = [("[amr]", INK), ("max_level       = 3", MUTED), ("ratio           = 8      ; octree", MUTED),
             ("iu_ref_levels   = 2", MUTED), ("markers_number  = 1", MUTED),
             ("seam_ghost_fill = tricubic ; default", MUTED), ("", MUTED),
             ("[amr_marker_1]", INK), ("mode         = 1       ; geometric", MUTED),
             ("geo_type     = primitive-box", MUTED), ("box_xmin     = 0.75    ; ... box bounds", MUTED),
             ("target_level = 3", MUTED), ("", MUTED), ("[initial_conditions]", INK),
             ("amr_iterations = 1      ; refine at init", MUTED)]
    s.rect(tx, oy, 420, 15 * 19 + 14, fill="#F5F7FA", stroke=CARD_EDGE, rx=8)
    for k, (t, c) in enumerate(lines):
        s.text(tx + 16, oy + 24 + 19 * k, t, 12, c, mono=True)
    s.save(out / "intra-amr.svg")


def fig_mixed(out: Path) -> None:
    """A forest whose realm 2 is refined inside: both register kinds in one register."""
    s = Svg(860, 420, "Mirror seam plus intra-realm AMR")
    s.text(430, 34, "Both seam families at once: a mirror seam and a 2:1 face inside realm 2", 17, anchor="middle",
           weight="bold")
    b, ox, oy = 60, 190, 80
    s.cells(ox, oy, 4, 4, b, A_FILL, A_EDGE, 2.4)
    s.cells(ox + 4 * b, oy, 4, 4, b, B_FILL, B_EDGE, 2.4)
    s.rect(ox + 7 * b, oy, b, 4 * b, fill=B_DARK, stroke="none", opacity=0.6)
    s.grid(ox + 7 * b, oy, 2, 8, b / 2, stroke=B_EDGE, sw=1.2)
    s.line(ox + 4 * b, oy - 10, ox + 4 * b, oy + 4 * b + 10, SEAM, 3.6)
    s.line(ox + 7 * b, oy, ox + 7 * b, oy + 4 * b, REG_EDGE, 4.0)
    s.text(ox + 4 * b, oy - 16, "inter-realm mirror seam", 12, SEAM, "middle", "bold")
    s.text(ox + 7 * b, oy - 16, "intra-realm 2:1 face", 12, REG_EDGE, "middle", "bold")
    s.text(ox + 2 * b, oy + 4 * b + 26, "realm 1", 14, A_EDGE, "middle", "bold")
    s.text(ox + 6 * b, oy + 4 * b + 26, "realm 2 (refined on x > 0.75)", 14, B_EDGE, "middle", "bold")
    s.text(430, oy + 4 * b + 56, "one flux register: the intra-realm 2:1 faces first (cursors 1..n_intra), "
           "then the inter-realm faces", 12, MUTED, "middle")
    s.text(430, oy + 4 * b + 74, "the seam face itself must keep one cell size per side: refinement may not reach "
           "the seam face", 12, MUTED, "middle")
    s.save(out / "mixed.svg")


def fig_register(out: Path) -> None:
    """The three ways a fine-side skin lands in a register face."""
    s = Svg(880, 340, "Flux register cases")
    s.text(440, 34, "Where the fine-side face fluxes land in the register", 17, anchor="middle", weight="bold")
    d, oy = 26, 110
    for title, x in (("mirror, lined up", 150), ("mirror, blocks not lined up", 440), ("refined 2:1", 730)):
        s.text(x, oy - 34, title, 14, INK, "middle", "bold")
    s.text(40, oy + d - 7, "register", 11, REG_EDGE)
    s.text(40, oy + 70 + d - 7, "fine side", 11, B_EDGE)
    x = 98
    s.cells(x, oy, 4, 1, d, REG_FILL, REG_EDGE, 2)
    s.cells(x, oy + 70, 4, 1, d, B_FILL, B_EDGE, 2)
    s.line(x + 2 * d, oy + 66, x + 2 * d, oy + d + 6, REG_EDGE, 1.6, arrow="arrR")
    s.text(150, oy + 130, "one overlap row, whole skin,", 12, MUTED, "middle")
    s.text(150, oy + 146, "copied 1:1", 12, MUTED, "middle")
    x = 334
    s.cells(x, oy, 4, 1, d, REG_FILL, REG_EDGE, 2)
    s.cells(x + 4 * d + 4, oy, 4, 1, d, REG_FILL, REG_EDGE, 2)
    s.cells(x + 2, oy + 70, 8, 1, d, B_FILL, B_EDGE, 2)
    s.line(x + 2 * d, oy + 66, x + 2 * d, oy + d + 6, REG_EDGE, 1.6, arrow="arrR")
    s.line(x + 6 * d + 4, oy + 66, x + 6 * d + 4, oy + d + 6, REG_EDGE, 1.6, arrow="arrR")
    s.text(440, oy + 130, "one row per register face it overlaps,", 12, MUTED, "middle")
    s.text(440, oy + 146, "rectangles scattered at their offsets", 12, MUTED, "middle")
    x = 678
    s.cells(x, oy, 4, 1, d, REG_FILL, REG_EDGE, 2)
    s.cells(x, oy + 70, 8, 1, d / 2, B_FILL, B_EDGE, 2)
    for k in range(4):
        s.line(x + (k + 0.5) * d, oy + 66, x + (k + 0.5) * d, oy + d + 6, REG_EDGE, 1.2, arrow="arrR")
    s.text(730, oy + 130, "2x2 fine faces averaged into one", 12, MUTED, "middle")
    s.text(730, oy + 146, "coarse face of the block's quadrant", 12, MUTED, "middle")
    s.text(440, oy + 188, "register side (realm_a for mirror, the coarse realm for refined): its skin is F_coarse; "
           "after the step the coarse cells", 12, MUTED, "middle")
    s.text(440, oy + 206, "get the Berger-Colella correction (dt/dx)(F_coarse - sum F_fine); on a mirror seam it is "
           "exactly 0", 12, MUTED, "middle")
    s.save(out / "register.svg")


def fig_cadence(out: Path) -> None:
    """Timelines of the seam fill: beta (every stage), alpha (end of step), alpha with asymmetric K."""
    s = Svg(880, 410, "Seam cadence")
    s.text(440, 34, "When the seam ghosts are filled: coupling_cadence", 17, anchor="middle", weight="bold")
    x0, x1 = 230, 820
    rows = [("beta (stage_coincident)", 5, 5, "every stage, both realms", 110),
            ("alpha (end_of_step)", 5, 5, "once, after the step", 220),
            ("alpha, asymmetric K", 5, 3, "R2 idles its trailing stages", 330)]
    for label, ka, kb, note, y in rows:
        s.text(30, y + 2, label, 13, INK, "start", "bold")
        s.text(30, y + 20, note, 11.5, MUTED)
        for k, (kk, col, dy) in enumerate(((ka, A_EDGE, -16), (kb, B_EDGE, 16))):
            yy = y + dy
            s.line(x0, yy, x1, yy, col, 2.2)
            for st in range(kk):
                xs = x0 + (st + 0.5) * (x1 - x0) / ka
                s.rect(xs - 7, yy - 7, 14, 14, fill=CARD, stroke=col, sw=1.6, rx=3)
            s.text(x0 - 10, yy + 4, "R1" if k == 0 else "R2", 11, col, "end", "bold")
        if label.startswith("beta"):
            for st in range(ka):
                xs = x0 + (st + 0.5) * (x1 - x0) / ka - 18
                s.path(f"M{xs:.1f},{y - 5:.1f} L{xs + 5:.1f},{y:.1f} L{xs:.1f},{y + 5:.1f} L{xs - 5:.1f},{y:.1f} z",
                       SEAM, 1.0, fill=SEAM)
        else:
            s.path(f"M{x1:.1f},{y - 6:.1f} L{x1 + 6:.1f},{y:.1f} L{x1:.1f},{y + 6:.1f} L{x1 - 6:.1f},{y:.1f} z",
                   SEAM, 1.0, fill=SEAM)
        s.line(x1, y - 26, x1, y + 26, MUTED, 1.0, dash="3 3")
        s.line(x0, y - 26, x0, y + 26, MUTED, 1.0, dash="3 3")
    s.text(x0, 70, "step n", 12, MUTED, "middle")
    s.text(x1, 70, "step n+1", 12, MUTED, "middle")
    s.text(440, 392, "squares: RK stages;  red diamonds: seam fills (beta: before each stage's residual; alpha: "
           "after close_step)", 12, MUTED, "middle")
    s.save(out / "cadence.svg")


def fig_ranks(out: Path) -> None:
    """Rank partition: an x split keeps the seam rank-local, a z split puts every seam row across ranks."""
    s = Svg(880, 470, "Seams across ranks")
    s.text(440, 34, "Each realm is partitioned on its own: a seam may join different ranks (issue #40)", 17,
           anchor="middle", weight="bold")
    b = 34
    r0, r1, q0, q1 = "#E4ECF7", "#A9C3E8", "#FDEBDD", "#F4BE90"

    def realm(x: float, y: float, nx: int, ny: int, f0: str, f1: str, edge: str) -> None:
        for j in range(ny):
            for i in range(nx):
                s.rect(x + i * b, y + j * b, b, b, fill=f1 if j < ny // 2 else f0, stroke=GRID, sw=0.8)
        s.rect(x, y, nx * b, ny * b, stroke=edge, sw=2.2)

    x, y = 90, 90
    realm(x, y, 3, 8, r0, r1, A_EDGE)
    realm(x + 3 * b, y, 3, 8, q0, q1, B_EDGE)
    s.line(x + 3 * b, y - 10, x + 3 * b, y + 8 * b + 10, SEAM, 3.2)
    for j in (1, 6):
        s.line(x + 3.5 * b, y + (j + 0.5) * b, x + 2.5 * b + 4, y + (j + 0.5) * b, OK, 1.6, arrow="arrG")
    s.text(x + 3 * b, y + 8 * b + 28, "split along x: each seam pair on one rank", 13, INK, "middle", "bold")
    s.text(x + 3 * b, y + 8 * b + 46, "local rows: a direct copy", 12, OK, "middle")
    x = 520
    realm(x, y + 4 * b, 6, 4, r0, r1, A_EDGE)
    realm(x, y, 6, 4, q0, q1, B_EDGE)
    s.line(x - 10, y + 4 * b, x + 6 * b + 10, y + 4 * b, SEAM, 3.2)
    for i in (1, 3, 5):
        s.line(x + (i + 0.5) * b, y + 4.5 * b, x + (i + 0.5) * b, y + 3.5 * b + 4, SEAM, 1.6, arrow="arrS")
    s.text(x + 3 * b, y + 8 * b + 28, "split along z: rank 1 cells, rank 0 ghosts", 13, INK, "middle", "bold")
    s.text(x + 3 * b, y + 8 * b + 46, "every row received / sent through MPI (tag 4040)", 12, SEAM, "middle")
    s.rect(300, 432, 18, 12, fill=r0, stroke=A_EDGE)
    s.rect(322, 432, 18, 12, fill=r1, stroke=A_EDGE)
    s.text(348, 442, "rank 0 / rank 1 (2 ranks; each realm partitioned along z in Morton order)", 12, MUTED)
    s.save(out / "ranks.svg")


def main() -> None:
    """Write every figure."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "public" / "forest")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for fig in (fig_anatomy, fig_mirror_aligned, fig_mirror_misaligned, fig_refined, fig_intra_amr, fig_mixed,
                fig_register, fig_cadence, fig_ranks):
        fig(args.out)
    print(f"figures written to {args.out}")


if __name__ == "__main__":
    main()
