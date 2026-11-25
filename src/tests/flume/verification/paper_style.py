"""Turn the figures of a docs figure script into publication figures, and trace each one to its runs (issue #142).

Why this module exists: the paper's figures must come from the same loaders and drawing code as the documentation
(make_doc_figures.py, 5000+ lines), not from a copy that drifts. A docs figure function draws as it always did; inside
`render(spec, ...)` its `Figure.savefig` is intercepted (the same seam as the dark theme of docs/scripts/figure_theme.py)
and the figure is restyled before it is written:

* the figure title and the panel titles are removed (the caption carries them), or replaced by the spec's short panel
  titles; panels get the labels (a), (b), ... at the same place (the left title slot, so the layout makes room);
* every text goes through `relabel`: the code identifiers of the docs (input keys, leg ids, issue numbers) become
  physical and scheme names, from the one dictionary LABELS; `audit` lists what is left;
* sizes are those of the printed page (figure_theme.PAGE_WIDTHS), every text 8 pt (panel labels 9 pt);
* colour maps: perceptually uniform for unsigned data, diverging (and centred) only for signed data; the 'tab:'
  colours of the docs onto the colour-blind-safe palette; labelled lines without markers get markers (grey scale);
* field maps (meshes, images, filled contours) are rasterized inside the vector PDF at 600 dpi; a PNG is written too.

The drawing itself happens inside figure_theme.paper(), the one publication rcParams block shared with PRISM.
`record_reads` notes every file the drawing opens; `provenance` reads what a run directory says about itself;
`manifest` writes figures.json: per figure its files, source runs and their provenance, and the ADAM commit of the
script run. Needs matplotlib (and h5py, for the record of the HDF5 reads).
"""

from __future__ import annotations

import builtins
import datetime as _dt
import io
import json
import re
import subprocess
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl
from matplotlib.collections import Collection, PathCollection, QuadMesh
from matplotlib.colors import Normalize, TwoSlopeNorm, to_rgba
from matplotlib.contour import ContourSet
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from matplotlib.text import Text

from figure_theme import (
    PAGE_INCLUDES,
    PAGE_WIDTHS,
    PAPER_CMAPS,
    PAPER_DPI_PDF,
    PAPER_DPI_PNG,
    PAPER_FONT_PT,
    PAPER_MARKERS,
    PAPER_PALETTE,
    PAPER_PANEL_LABEL_PT,
    PAPER_TAB,
    paper,
)

# Code identifiers of the docs figures -> names for a reader of the paper. Applied as plain substring replacements,
# longest key first, to every text of a figure (legends, axis and colorbar labels, kept titles, annotations). Add the
# PRISM identifiers here when its figures move to this style (Part D).
LABELS = {
    # schemes (all WENO5 runs: weno-u-5; `weno` = characteristic flux splitting, `weno-riemann` = WENO5 states into a
    # Riemann solver with the sixth-order flux correction)
    "weno (split)": "WENO5, flux splitting",
    "weno-riemann LLF": "WENO5, LLF",
    "weno-riemann HLLC": "WENO5, HLLC",
    "weno-riemann HLLD": "WENO5, HLLD",
    "weno-riemann HLL": "WENO5, HLL",
    "weno-u-5": "WENO5",
    "weno-u-7": "WENO7",
    "weno-u-9": "WENO9",
    # bounds named after oracle legs
    "av3c bound (HLLD)": "bound of the HLLD run",
    # quantities
    "log10 rho": r"$\log_{10}\rho$",
    "alfven": "Alfvén",
    "32^3": "$32^3$",
    "64^3": "$64^3$",
    "128^3": "$128^3$",
    "[-1/2, 1/2]^3": r"$[-\frac{1}{2}, \frac{1}{2}]^3$",
    "central density / initial - 1": r"$\rho_c(t)/\rho_c(0) - 1$",
    "largest Mach number inside 0.8 a": r"max Mach number, $r < 0.8a$",
    "dM / dt over m_0 a^3 / G": r"$\dot M\, G / (m_0 a^3)$",
    "m_0(29.3) = 138.9": r"$m_0(29.3) = 138.9$",
    "E + W": "$E + W$",
    "E alone": "$E$ alone",
    "W alone": "$W$ alone",
    "t / P (Kepler)": r"$t / P_\mathrm{Kepler}$",
    "Kepler circle, d = 0.5": "Kepler circle, $d = 0.5$",
    " (floor 1e-17)": "",
    "step (grey: regrid steps)": "step (grey: regrids)",
}
# the same by pattern: the conserved variables named by their input keys at the end of a legend entry ("linear,
# rE"), and the bare coordinate labels
LABEL_RULES = (
    (re.compile(r"(?<=, )rE$"), r"$\\rho E$"),
    (re.compile(r"(?<=, )r$"), r"$\\rho$"),
    (re.compile(r"^([xyzt])$"), r"$\1$"),  # coordinate and time axis labels in math italics
    (re.compile(r"(?<!flux )\bsplitting(?=, | scheme)"), "flux splitting"),
    (re.compile(r"\bsplitting scheme\b"), "splitting"),
    # reference lines carry the design order
    (re.compile(r"^order (\d+)$"), r"slope \1 (design order)"),
    (re.compile(r"^\$N\^\{-(\d+)\}\$$"), r"$N^{-\1}$ (design order \1)"),
)
_LABEL_RX = re.compile("|".join(re.escape(k) for k in sorted(LABELS, key=len, reverse=True)))
# identifiers that must not reach a paper figure: leg ids (MV-5, SV-7, si11, rg3, av3c), input keys and values
# (weno-u-5, weno-riemann, nk = 1, snake_case), issue numbers, ASCII exponents (32^3) outside math
_SUSPECT = re.compile(
    r"weno-u-\d|weno-riemann|weno \(split\)|\b[A-Z]{1,3}V-\d+\b|\b(?:si|sv|rg|av|co|eo|ps|de|ev|mv|gv)\d+[a-z]?\b"
    r"|#\d+|\bnk = \d|\b[a-z]+_[a-z0-9_]+\b|\w\^\d|\bF\d+\b|\bquadtree \(|\boctree, null"
)
_MATH = re.compile(r"\$[^$]*\$")
_MATH_SPLIT = re.compile(r"(\$[^$]*\$)")
PANEL_LABEL = r"$\mathbf{(%s)}$"  # bold in the math font (serif bold): (a), (b), ...
_SEQUENTIAL_OK = {"viridis", "cividis", "magma", "inferno", "plasma", "gray", "Greys", "binary", "bone"}
_DIVERGING = {"coolwarm", "RdBu", "bwr", "seismic", "PuOr", "BrBG", "RdYlBu", "RdBu_r", "PiYG", "PRGn", "RdGy"}


@dataclass
class Spec:
    """A paper figure: its file name (physical, not the docs name), the figure function that draws it, the docs file
    stem that function saves (its `savefig` is intercepted), the page width ('full', 'text', 'half'), the height in
    inches, the panel titles (None: removed; 'keep': kept through `relabel`; a list: one per panel, '' to remove) and
    the colorbar labels (in colorbar order, None keeps one), legend keywords for every legend (e.g. {"loc": ...}),
    whole-text replacements for this figure only (applied after LABELS, e.g. shorter legend entries), a note for
    the manifest, the fraction of the page width the figure is printed at (`scale`, e.g. 0.72 of the text width: the
    figure is drawn at that size, so its 8 pt stay 8 pt) and the provenance traced by hand (`traced`: the runs, their
    jobs, the commit that recorded the numbers, toolchain, device, dates and the caption facts), for runs whose
    directories carry no run.info (issue #144)."""

    name: str
    draw: Callable
    saves: str
    width: str
    height: float
    titles: list[str] | str | None = None
    cbar_labels: list[str | None] = field(default_factory=list)
    legend: dict = field(default_factory=dict)
    labels: dict = field(default_factory=dict)
    note: str = ""
    scale: float = 1.0
    traced: dict = field(default_factory=dict)

    @property
    def inches(self) -> float:
        """The printed width of the figure in inches."""
        return PAGE_WIDTHS[self.width] * self.scale

    @property
    def latex(self) -> str:
        """The includegraphics line of the figure."""
        include = PAGE_INCLUDES[self.width]
        if self.scale != 1.0:
            include = include.replace(r"width=\textwidth", rf"width={self.scale:g}\textwidth")
        return include % self.name


def relabel(text: str) -> str:
    """Return a text with the code identifiers of LABELS (one pass, longest key first, outside math: never inside a
    replacement) and of LABEL_RULES replaced by their names. Idempotent: relabel(relabel(t)) == relabel(t)."""
    if not text:
        return text
    parts = _MATH_SPLIT.split(text)  # odd entries are the $...$ spans
    out = "".join(p if i % 2 else _LABEL_RX.sub(lambda m: LABELS[m.group(0)], p) for i, p in enumerate(parts))
    for rx, rep in LABEL_RULES:
        out = rx.sub(rep, out)
    return out


def audit(fig: Figure) -> list[str]:
    """Return the texts of a figure that still hold a code identifier (outside math), and the unlabelled colorbars."""
    found = []
    for t in fig.findobj(Text):
        s = t.get_text()
        if s and t.get_visible() and _SUSPECT.search(_MATH.sub("", s)):
            found.append(f"identifier in text: {s!r}")
    for ax in fig.axes:
        if getattr(ax, "_colorbar", None) is not None and not (ax.get_ylabel() or ax.get_xlabel()):
            found.append("colorbar without a label")
    return sorted(set(found))


def _main_axes(fig: Figure) -> list:
    """The panels of a figure in reading order: no colorbars, no twins, no insets."""
    out = []
    for ax in fig.axes:
        if getattr(ax, "_colorbar", None) is not None or hasattr(ax, "_colorbar_info"):
            continue
        twins = [a for a in ax.get_shared_x_axes().get_siblings(ax) if a is not ax and a.bbox.bounds == ax.bbox.bounds]
        if any(a in out for a in twins):
            continue
        if ax.get_axes_locator() is not None:
            continue
        out.append(ax)

    def place(ax) -> tuple:  # noqa: ANN001
        spec = ax.get_subplotspec()
        if spec is not None:  # the grid cell: not moved by equal-aspect shrinking
            return (spec.rowspan.start, spec.colspan.start)
        box = ax.get_position(original=True)
        return (-box.y1, box.x0)

    return sorted(out, key=place)


def _tab_to_palette(color) -> tuple | None:  # noqa: ANN001
    """The palette colour of a hard-coded 'tab:' colour, or None."""
    rgba = to_rgba(color)
    for tab, idx in PAPER_TAB.items():
        if all(abs(a - b) < 1e-3 for a, b in zip(rgba[:3], to_rgba(tab)[:3], strict=True)):
            return (*to_rgba(PAPER_PALETTE[idx])[:3], rgba[3])
    return None


def _cmap_policy(artist) -> None:  # noqa: ANN001
    """Perceptually uniform maps; a diverging map only for signed data, centred on zero."""
    cmap = artist.get_cmap()
    name = cmap.name.removesuffix("_r")
    norm = artist.norm
    signed = norm is not None and norm.vmin is not None and norm.vmin < 0.0 < norm.vmax
    if name in _DIVERGING:
        if signed:
            artist.set_cmap(PAPER_CMAPS["diverging"])
            if type(norm) is Normalize:
                m = max(abs(norm.vmin), abs(norm.vmax))
                artist.set_norm(TwoSlopeNorm(0.0, -m, m))
        else:
            artist.set_cmap(PAPER_CMAPS["alternative"])
    elif name not in _SEQUENTIAL_OK:
        artist.set_cmap(PAPER_CMAPS["sequential"])


def restyle(fig: Figure, spec: Spec) -> list[str]:
    """Restyle a docs figure into the paper figure of `spec` (see the module docstring); return the audit."""
    if fig._suptitle is not None:
        fig._suptitle.set_text("")

    def name(text: str) -> str:
        new = relabel(text)
        return spec.labels.get(new, new)

    panels = _main_axes(fig)
    for i, ax in enumerate(panels):
        title = ""
        if spec.titles == "keep":
            title = name(ax.get_title())
        elif isinstance(spec.titles, list) and i < len(spec.titles):
            title = name(spec.titles[i])
        for loc in ("left", "center", "right"):
            ax.set_title("", loc=loc)
        # the panel label and the panel title share the left slot: one position, no collision
        label = PANEL_LABEL % chr(ord("a") + i) if len(panels) > 1 else ""
        text = f"{label} {title}".strip()
        if text:
            ax.set_title(text, loc="left", fontsize=PAPER_PANEL_LABEL_PT)
    cbars = [ax for ax in fig.axes if getattr(ax, "_colorbar", None) is not None]
    for cax, lab in zip(cbars, spec.cbar_labels, strict=False):
        if lab is not None:
            cax._colorbar.set_label(lab)
    for ax in fig.axes:
        for artist in [*ax.collections, *ax.images]:
            if isinstance(artist, (QuadMesh, AxesImage, ContourSet)) or (
                isinstance(artist, Collection) and artist.get_array() is not None
            ):
                if getattr(ax, "_colorbar", None) is None:
                    _cmap_policy(artist)
                if isinstance(artist, (QuadMesh, AxesImage, ContourSet)) or (
                    isinstance(artist, PathCollection) and len(artist.get_offsets()) > 1000
                ):
                    artist.set_rasterized(True)
        markers = iter(PAPER_MARKERS * 4)
        for line in ax.get_lines():
            for getter, setter in (("get_color", "set_color"), ("get_markerfacecolor", "set_markerfacecolor"),
                                   ("get_markeredgecolor", "set_markeredgecolor")):
                c = getattr(line, getter)()
                if isinstance(c, str) and c == "none":
                    continue
                new = _tab_to_palette(c)
                if new is not None:
                    getattr(line, setter)(new)
            n = len(line.get_xdata())
            rgb = to_rgba(line.get_color())[:3]
            reference = max(rgb) - min(rgb) < 0.02  # black and grey: exact solutions, bounds, slopes
            if (not line.get_label().startswith("_") and line.get_marker() in ("None", "", None, " ")
                    and n > 2 and line.get_linestyle() not in (":", "None") and not reference):
                line.set_marker(next(markers))
                line.set_markevery(max(1, n // 8))
                line.set_markerfacecolor("white")
                line.set_markersize(3.0)
        for coll in ax.collections:
            if isinstance(coll, PathCollection) and coll.get_array() is None:
                fc = coll.get_facecolor()
                if len(fc):
                    new = _tab_to_palette(fc[0])
                    if new is not None:
                        coll.set_facecolor(new)
        ax.tick_params(labelsize=PAPER_FONT_PT, which="both")
    for t in fig.findobj(Text):
        if t.get_text().startswith(PANEL_LABEL[:9]):
            continue  # a panel title: named above, once
        new = name(t.get_text())
        if new != t.get_text():
            t.set_text(new)
        t.set_fontsize(PAPER_FONT_PT)
    for ax in fig.axes:  # rebuilt: the handles carry the markers added above
        leg = ax.get_legend()
        if leg is not None:
            artists = {name(lab): h for h, lab in zip(*ax.get_legend_handles_labels(), strict=True)}
            texts = [t.get_text() for t in leg.get_texts()]
            handles = [artists.get(t, h) for t, h in zip(texts, leg.legend_handles, strict=True)]
            kw = {"loc": leg._loc, "ncols": leg._ncols} | spec.legend
            ax.legend(handles, texts, fontsize=PAPER_FONT_PT, **kw)
    fig.set_size_inches(spec.inches, spec.height)
    engine = fig.get_layout_engine()
    if engine is None:
        try:
            fig.set_layout_engine("constrained")
            engine = fig.get_layout_engine()
        except RuntimeError:
            engine = None
    if engine is not None and hasattr(engine, "set") and type(engine).__name__ == "ConstrainedLayoutEngine":
        engine.set(w_pad=1.5 / 72, h_pad=1.5 / 72, wspace=0.03, hspace=0.03)
    fig.canvas.draw()
    return audit(fig)


@contextmanager
def record_reads() -> Iterator[set[Path]]:
    """Record the files opened for reading in the block (HDF5 through h5py, text through open/io.open and
    numpy.loadtxt, whose C reader bypasses open)."""
    import numpy as np  # noqa: PLC0415

    seen: set[Path] = set()
    real_open, real_io_open, real_loadtxt = builtins.open, io.open, np.loadtxt
    try:
        import h5py  # noqa: PLC0415

        real_h5 = h5py.File
    except ImportError:
        h5py, real_h5 = None, None

    def note(path, mode: str) -> None:  # noqa: ANN001
        if isinstance(path, (str, Path)) and "r" in mode and "+" not in mode:
            seen.add(Path(path).resolve())

    def open_rec(file, mode="r", *args, **kwargs):  # noqa: ANN001, ANN002, ANN003, ANN202
        note(file, mode)
        return real_open(file, mode, *args, **kwargs)

    def loadtxt_rec(fname, *args, **kwargs):  # noqa: ANN001, ANN002, ANN003, ANN202
        note(fname, "r")
        return real_loadtxt(fname, *args, **kwargs)

    builtins.open = io.open = open_rec
    np.loadtxt = loadtxt_rec
    if h5py is not None:

        class FileRec(real_h5):
            def __init__(self, name, mode="r", *args, **kwargs):  # noqa: ANN001, ANN002, ANN003, ANN204
                note(name, mode)
                super().__init__(name, mode, *args, **kwargs)

        h5py.File = FileRec
    try:
        yield seen
    finally:
        builtins.open, io.open = real_open, real_io_open
        np.loadtxt = real_loadtxt
        if h5py is not None:
            h5py.File = real_h5


def _git(repo: Path, *args: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return ""


def script_provenance(repo: Path) -> dict:
    """The ADAM commit of the figure script run, and whether src/ and docs/scripts/ differ from it (the changed files)."""
    head = _git(repo, "rev-parse", "HEAD").strip() or "unknown"
    changed = [ln[3:] for ln in _git(repo, "status", "--porcelain", "--", "src", "docs/scripts").splitlines() if ln]
    return {"adam_commit": head, "src_differs": bool(changed), "changed_files": changed[:200]}


_LOG_KEYS = {
    "commit": re.compile(r"\b(?:commit|git[ _-]?(?:hash|sha|revision))\b\s*[:=]?\s*([0-9a-f]{7,40})", re.IGNORECASE),
    "executable": re.compile(r"\b(?:executable|exe|program)\b\s*[:=]\s*(\S+)", re.IGNORECASE),
    "host": re.compile(r"\b(?:hostname|host|node)\b\s*[:=]\s*(\S+)", re.IGNORECASE),
    "device": re.compile(r"\b(?:device|gpu)\b(?: name)?\s*[:=]\s*(.+)$", re.IGNORECASE),
    "date": re.compile(r"\b(?:date|started)\b\s*[:=]\s*(\S.*)$", re.IGNORECASE),
    "compiler": re.compile(r"\bcompiler(?: version)?\b\s*[:=]\s*(\S.*)$", re.IGNORECASE),
}


def provenance(run: Path) -> dict:
    """What a run directory says about itself: its run.info (key=value) if present, else the header of its log (the
    first 400 lines), else 'unknown'. The executable implied by the work directory name (work-<exe>-...) and the
    modification time of its newest checkpoint are added, labelled as such: they are not the run's own record."""
    info: dict = {"source": "unknown"}
    runinfo = run / "run.info"
    if runinfo.is_file():
        info = {"source": "run.info"}
        for ln in runinfo.read_text(errors="replace").splitlines():
            if "=" in ln:
                k, v = ln.split("=", 1)
                info[k.strip()] = v.strip()
    else:
        logs = sorted([*run.glob("log*.txt"), *run.glob("*.log")])
        if logs:
            head = logs[0].read_text(errors="replace").splitlines()[:400]
            found = {}
            for ln in head:
                for key, rx in _LOG_KEYS.items():
                    m = rx.search(ln)
                    if m and key not in found:
                        found[key] = m.group(1).strip()
            info = {"source": f"log header ({logs[0].name})" if found else f"log header ({logs[0].name}): nothing"}
            info.update(found)
    for key in ("commit", "executable", "host", "device", "date"):
        info.setdefault(key, "unknown")
    m = re.match(r"work-(adam_[a-z0-9_]+?)(?:-np\d+)?-", run.name)
    if m:
        info["executable_from_dirname"] = m.group(1)
    cps = [p for p in run.glob("*.h5")]
    if cps:
        newest = max(p.stat().st_mtime for p in cps)
        info["checkpoint_mtime"] = _dt.datetime.fromtimestamp(newest).astimezone().isoformat(timespec="seconds")
    return info


def runs_of(paths: set[Path]) -> tuple[list[Path], list[Path]]:
    """Split the files read by a figure into run directories (holding checkpoints, a log or histories) and others."""
    skip = (".py", ".pyc", ".ttf", ".otf", ".ttc", ".afm", ".pfb", ".json", ".mplstyle")
    dirs = sorted({p.parent for p in paths if p.suffix.lower() not in skip and "mpl-data" not in p.parts})
    runs, inputs = [], []
    for d in dirs:
        is_run = any(d.glob("*.h5")) or (d / "log.txt").exists() or any(d.glob("*_history.dat"))
        (runs if is_run else inputs).append(d)
    return runs, inputs


def render(spec: Spec, out: Path, runs: Path, usetex: bool = False) -> dict:
    """Draw one paper figure: the spec's function inside the paper style, its savefig intercepted, restyled and
    written as <name>.pdf and <name>.png in `out`. Return its manifest entry (raises what the function raises)."""
    entry: dict = {"name": spec.name, "function": spec.draw.__name__, "docs_name": spec.saves, "width": spec.width,
                   "size_in": [round(spec.inches, 3), spec.height], "latex": spec.latex}
    if spec.scale != 1.0:
        entry["scale"] = spec.scale
    if spec.note:
        entry["note"] = spec.note
    if spec.traced:
        entry["traced"] = spec.traced
    written: list[str] = []
    savefig = Figure.savefig

    def savefig_paper(fig: Figure, fname, *args, **kwargs):  # noqa: ANN001, ANN202, ARG001
        stem = Path(fname).stem
        if stem != spec.saves:
            print(f"{spec.name}: {stem} is not this paper figure ({spec.saves}), not written")
            return None
        entry["audit"] = restyle(fig, spec)
        entry.update(getattr(fig, "paper_notes", {}))
        for ext, dpi in (("pdf", PAPER_DPI_PDF), ("png", PAPER_DPI_PNG)):
            path = out / f"{spec.name}.{ext}"
            savefig(fig, path, dpi=dpi, metadata={"Creator": "ADAM make_doc_figures.py --style paper"}
                    if ext == "pdf" else None)
            written.append(path.name)
        return None

    with record_reads() as seen, paper(usetex):
        Figure.savefig = savefig_paper
        try:
            spec.draw(out, runs)
        finally:
            Figure.savefig = savefig
    if not written:
        raise FileNotFoundError(f"{spec.draw.__name__} saved no {spec.saves}")
    run_dirs, input_dirs = runs_of(seen)
    entry["files"] = written
    entry["runs"] = [{"path": str(d), "provenance": provenance(d)} for d in run_dirs]
    entry["inputs"] = [str(d) for d in input_dirs]
    return entry


def manifest(out: Path, entries: list[dict], repo: Path, script: Path, skipped: dict[str, str]) -> Path:
    """Write figures.json next to the figures: the script's provenance and one entry per figure."""
    doc = {
        "generated": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": str(script),
        "style": "paper",
        **script_provenance(repo),
        "matplotlib": mpl.__version__,
        "page": {k: round(v, 4) for k, v in PAGE_WIDTHS.items()} | {"unit": "inch"},
        "figures": entries,
        "skipped": skipped,
    }
    path = out / "figures.json"
    path.write_text(json.dumps(doc, indent=2) + "\n")
    return path
