"""The light (washi) and dark (silver ink) styles of the docs' matplotlib figures (docs/BRAND.md).

Why this module exists: the docs show every figure in both themes. A figure script draws its figures as it always did
(the light style, unchanged); inside `theme("dark")` the same drawing code produces `<name>-dark.png` next to
`<name>.png`, on the silver-ink ground. The docs theme (docs/.vitepress/theme/index.ts) shows the dark twin in the
dark theme when it exists, and the light figure as a washi card otherwise.

How the dark style works, without touching the drawing code of the figures:

* the matplotlib defaults (text, axes, ticks, grid, legend frame, ground) are set to the silver-ink tokens for the
  figures created inside the context;
* just before saving, every colour a figure hard-codes as a neutral grey ('k', 'black', '0.4', '#333', 'w' fills, the
  matplotlib defaults that are neutral) is mapped by `ink`: black becomes the silver ink, white the ground, a grey `g`
  the blend `ground + (1 - g) (ink - ground)`, so the contrast ladder of the light figure is kept. White strokes stay
  white (they are drawn over field images). Coloured strokes and texts too dark on the ground are lightened to a 3:1
  contrast (`legible`); coloured fills are kept;
* field images (images, meshes, filled contours, scatter colour maps, 3-D surfaces) are never changed: they keep their
  colour maps.

The silver-ink tokens are slightly blue, never neutral, so the mapping leaves them alone and is idempotent.
Needs matplotlib.
"""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import matplotlib as mpl
from matplotlib.collections import Collection
from matplotlib.colors import ListedColormap, to_rgba
from matplotlib.contour import ContourSet
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.text import Text

THEMES = ("light", "dark", "both")
WASHI = "#f3ede1"
# silver ink (docs/BRAND.md, "Themes"): ground, bg-alt (the axes panel), text-1, text-2, mesh, border
GROUND, PANEL, INK, TEXT2, MESH, BORDER = (
    "#0e1116",
    "#0a0d11",
    "#e8ecf0",
    "#a9b1bb",
    "#8b949e",
    "#2c333d",
)
DARK_RC = {
    "figure.facecolor": GROUND,
    "savefig.facecolor": GROUND,
    "axes.facecolor": PANEL,
    "axes.edgecolor": TEXT2,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": TEXT2,
    "ytick.color": TEXT2,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "grid.color": MESH,
    "legend.facecolor": "inherit",
    "legend.edgecolor": BORDER,
    "patch.edgecolor": INK,
    "hatch.color": INK,
    "axes3d.xaxis.panecolor": (0.08, 0.094, 0.125, 0.5),
    "axes3d.yaxis.panecolor": (0.08, 0.094, 0.125, 0.5),
    "axes3d.zaxis.panecolor": (0.08, 0.094, 0.125, 0.5),
}
MIN_CONTRAST = 3.0


def _luminance(rgb: tuple[float, ...]) -> float:
    """Return the WCAG relative luminance of an sRGB colour."""
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in rgb[:3]]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def _contrast(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    la, lb = _luminance(a), _luminance(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


def _neutral(rgba: tuple[float, ...]) -> bool:
    """Whether a colour is a neutral grey (black and white included): the colours of the light style to remap."""
    return max(rgba[:3]) - min(rgba[:3]) <= 2.0 / 255.0


def ink(color, face: bool = False) -> tuple[float, float, float, float]:  # noqa: ANN001
    """Return the silver-ink colour of a light-style colour: neutral greys are mapped from the white-to-black ladder
    onto the ground-to-ink ladder; white is kept for strokes (`face` false) and becomes the ground for fills; other
    colours are returned unchanged."""
    rgba = to_rgba(color)
    if rgba[3] == 0.0 or not _neutral(rgba):
        return rgba
    level = sum(rgba[:3]) / 3.0
    if level >= 0.99 and not face:
        return rgba
    g, k = to_rgba(GROUND), to_rgba(INK)
    t = 1.0 - level
    return (*(g[i] + t * (k[i] - g[i]) for i in range(3)), rgba[3])


def legible(color) -> tuple[float, float, float, float]:  # noqa: ANN001
    """Return a stroke colour for the dark ground: neutral greys through `ink`, other colours lightened (blended with
    white) until their contrast with the ground reaches 3:1."""
    rgba = ink(color)
    if rgba[3] == 0.0 or _neutral(rgba):
        return rgba
    g = to_rgba(GROUND)
    out = rgba
    for step in range(1, 21):
        if _contrast(out, g) >= MIN_CONTRAST:
            break
        w = 0.05 * step
        out = (*(c + w * (1.0 - c) for c in rgba[:3]), rgba[3])
    return out


def _restyle_patch(p: Patch) -> None:
    p.set_facecolor(ink(p.get_facecolor(), face=True))
    p.set_edgecolor(legible(p.get_edgecolor()))


def _restyle_collection(c: Collection) -> None:
    if isinstance(c, ContourSet):
        # line contours of explicit colours ('k', 'w'): their colours are a listed colour map of the levels
        if not c.filled and c.colors is not None and isinstance(c.cmap, ListedColormap):
            c.set_cmap(ListedColormap([legible(x) for x in c.cmap.colors]))
        return
    if c.get_array() is not None or type(c).__name__ == "Poly3DCollection":
        return  # colour-mapped data or a 3-D surface: field colours
    c.set_facecolor([ink(x, face=True) for x in c.get_facecolor()])
    c.set_edgecolor([legible(x) for x in c.get_edgecolor()])


def restyle(fig: Figure) -> None:
    """Map the hard-coded light-style colours of a figure onto the silver-ink palette (see the module docstring)."""
    fig.canvas.draw()  # create the ticks of the final layout, so that their colours are mapped too
    for ax in fig.axes:
        if ax.name == "3d":
            ax.patch.set_facecolor(
                GROUND
            )  # a 3-D view has no panel: its panes are its frame
    for a in fig.findobj():
        if isinstance(a, Line2D):
            mec, mfc, mfa = (
                a.get_markeredgecolor(),
                a.get_markerfacecolor(),
                a.get_markerfacecoloralt(),
            )
            a.set_color(legible(a.get_color()))
            a.set_markeredgecolor(legible(mec))
            a.set_markerfacecolor(legible(mfc))
            a.set_markerfacecoloralt(legible(mfa))
        elif isinstance(a, Text):
            a.set_color(legible(a.get_color()))
            if a.get_bbox_patch() is not None:
                _restyle_patch(a.get_bbox_patch())
            if getattr(a, "arrow_patch", None) is not None:
                _restyle_patch(a.arrow_patch)
        elif isinstance(a, Patch) and a is not fig.patch:
            _restyle_patch(a)
        elif isinstance(a, Collection):
            _restyle_collection(a)


def dark_name(fname: str | Path) -> Path:
    """Return the dark twin of a figure file: `<name>-dark.png` next to `<name>.png`."""
    p = Path(fname)
    return p.with_name(f"{p.stem}-dark{p.suffix}")


@contextmanager
def theme(name: str) -> Iterator[None]:
    """Draw the figures of the block in a theme: 'light' changes nothing; 'dark' sets the silver-ink defaults and
    makes every `Figure.savefig` restyle the figure and write its dark twin (`dark_name`) instead."""
    if name not in ("light", "dark"):
        raise ValueError(f"unknown theme {name!r}")
    if name == "light":
        yield
        return
    savefig = Figure.savefig

    def savefig_dark(fig: Figure, fname, *args, **kwargs):  # noqa: ANN001, ANN202
        restyle(fig)
        return savefig(fig, dark_name(fname), *args, **kwargs)

    with mpl.rc_context(DARK_RC):
        Figure.savefig = savefig_dark
        try:
            yield
        finally:
            Figure.savefig = savefig


def themes(choice: str) -> tuple[str, ...]:
    """Return the themes of a `--theme` choice (light, dark or both)."""
    return ("light", "dark") if choice == "both" else (choice,)
