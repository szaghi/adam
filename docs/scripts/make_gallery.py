#!/usr/bin/env python3
"""Render the landing-page animations of the ADAM docs from the gallery runs (docs/scripts/make_gallery_runs.sh).

Each animation is drawn twice, on the two grounds of the sumi identity (docs/BRAND.md): washi (light) and silver ink
(dark). The frames carry no title, only a small time label in a strip under the panels, and no margin of another
colour. The colour maps run monotonically in lightness from the ground to the ink (light: washi, seal, ink; dark:
silver ink, seal, amber, paper), so the field reads as ink laid on the page. The colour range is fixed over the whole
animation. Per animation and theme it writes, into --out:

    <name>-<theme>.webm   VP9, constant quality (crf 36), 720 px wide
    <name>-<theme>.mp4    H.264, yuv420p, faststart
    <name>-<theme>.png    the poster: a representative frame

Animations and the run directories they read (sub-directories of --runs):

    orszag-tang        orszag-tang                    density, block outlines
    rayleigh-taylor    rt-hydro | rt-mhd              density, hydro left, MHD right, matched in time
    forward-facing-step step (realms A, B, C)         density, block outlines, inter-realm seams dashed
    prism-faces        pulse-pec | pulse-sm | pulse-pml
                                                      log10 sqrt(|E|^2 + c^2|B|^2), relative to its initial
                                                      peak, in the plane z ~ 0, PEC, Silver-Muller
                                                      and PML faces from left to right, the PML layer outlined

Usage: make_gallery.py --runs DIR [--out DIR] [--only NAME ...] [--fps N] [--width PX] [--keep-frames DIR]
Needs h5py, numpy, matplotlib and ffmpeg on PATH.
"""

from __future__ import annotations

import argparse
import configparser
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
CKPT = re.compile(r"^(?P<prefix>.+)-(?P<it>\d{9})-proc\d{6}\.h5$")
TIME = re.compile(r'<Time Value="([^"]+)"')
MU0, EPS0 = 1.256637e-6, 8.854187e-12  # PRISM's vacuum constants (adam_prism_parameters.F90)
C0 = 1.0 / np.sqrt(MU0 * EPS0)


@dataclass(frozen=True)
class Theme:
    """The colours of one ground of the sumi identity."""

    name: str
    bg: str
    bg_alt: str
    text: str
    mesh: str
    accent: str
    cmap: LinearSegmentedColormap


THEMES = (
    Theme("light", "#f3ede1", "#ebe3d3", "#1d1b18", "#6f6a62", "#b23417",
          LinearSegmentedColormap.from_list("washi", ["#f3ede1", "#e4b48f", "#c73e1d", "#7a2614", "#1d1b18"])),
    Theme("dark", "#0e1116", "#0a0d11", "#e8ecf0", "#8b949e", "#ffb020",
          LinearSegmentedColormap.from_list("ink", ["#0e1116", "#5a1f12", "#c73e1d", "#ffb020", "#fbf3e0"])),
)


@dataclass
class Frame:
    """One time level of one panel: a raster on the finest cells (NaN outside the domain) and the block outlines."""

    time: float
    data: np.ndarray  # (nx, ny), x first
    boxes: list[tuple[float, float, float, float]] = field(default_factory=list)


@dataclass
class Panel:
    """One panel of an animation: its frames, physical extent and decorations."""

    frames: list[Frame]
    extent: tuple[float, float, float, float]  # x0, x1, y0, y1
    outline: bool = False
    decorate: Callable | None = None  # decorate(ax, theme)


# ----------------------------------------------------------------------------------------------------------- reading


def read_ini(path: Path) -> configparser.ConfigParser:
    """Read an ADAM INI file (`;` comments, no interpolation, case-sensitive keys)."""
    ini = configparser.ConfigParser(inline_comment_prefixes=(";",), interpolation=None, strict=False)
    ini.optionxform = str
    ini.read(path)
    return ini


def checkpoints(work: Path, prefix: str | None = None) -> dict[int, tuple[float, list[Path]]]:
    """Return {iteration: (time, files)} of the field checkpoints of a run (of the realm `prefix`)."""
    groups: dict[int, list[Path]] = {}
    for p in work.glob("*-proc*.h5"):
        m = CKPT.match(p.name)
        if m is None or "restart" in p.name or (prefix is not None and m["prefix"] != prefix):
            continue
        groups.setdefault(int(m["it"]), []).append(p)
    out = {}
    for it, files in sorted(groups.items()):
        pre = CKPT.match(files[0].name)["prefix"]
        xdmf = work / f"{pre}-{it:09d}.xdmf"
        hit = TIME.search(xdmf.read_text()) if xdmf.exists() else None
        if hit is None:
            raise RuntimeError(f"no time for {xdmf}")
        out[it] = (float(hit.group(1)), sorted(files))
    if not out:
        raise FileNotFoundError(f"no checkpoint in {work} (prefix {prefix})")
    return out


def read_blocks(files: list[Path], ngc: int, names: tuple[str, ...]) -> list[dict]:
    """Return the blocks of one checkpoint: interior fields (x, y, z order), cell spacing, interior corner."""
    out = []
    for path in files:
        with h5py.File(path, "r") as h5:
            for blk in sorted({k.rsplit("-", 1)[0] for k in h5 if k.endswith("-origin")}):
                # geometry and fields are stored (z, y, x) with ghosts on every axis: reverse and strip them
                origin, d = h5[f"{blk}-origin"][()][::-1], h5[f"{blk}-dxdydz"][()][::-1]
                f = {v: h5[f"{blk}-{v}"][()].transpose(2, 1, 0)[ngc:-ngc, ngc:-ngc, ngc:-ngc] for v in names}
                out.append({"lo": origin + ngc * d, "d": d, "f": f})
    return out


def rasterize(planes: list[tuple[np.ndarray, np.ndarray, np.ndarray]], extent: tuple, fine: float) -> np.ndarray:
    """Paint 2-D block fields (value (nx, ny), lower corner, spacing) on a raster of cell `fine`, NaN where none."""
    x0, x1, y0, y1 = extent
    nx, ny = round((x1 - x0) / fine), round((y1 - y0) / fine)
    out = np.full((nx, ny), np.nan)
    for v, lo, d in planes:
        r = round(d[0] / fine)
        i0, j0 = round((lo[0] - x0) / fine), round((lo[1] - y0) / fine)
        up = np.repeat(np.repeat(v, r, axis=0), r, axis=1)
        out[i0:i0 + up.shape[0], j0:j0 + up.shape[1]] = up
    return out


def density_frames(work: Path, prefixes: list[str | None], extent: tuple, names: tuple[str, ...],
                   value: Callable[[dict], np.ndarray] | None = None) -> list[Frame]:
    """Return the frames of a FLUME run (every realm of `prefixes` on one raster): a 2-D field and block outlines."""
    inis = sorted(work.glob("*.ini"))
    ini = read_ini(next((p for p in inis if p.stem.endswith("-A")), inis[0]))
    ngc = int(ini["grid"]["ngc"])
    value = value or (lambda f: f["r"][:, :, 0])
    cps = {pre: checkpoints(work, pre) for pre in prefixes}
    its = sorted(set.intersection(*(set(c) for c in cps.values())))
    frames = []
    for it in its:
        planes, boxes = [], []
        for pre in prefixes:
            for b in read_blocks(cps[pre][it][1], ngc, names):
                v = value(b["f"])
                planes.append((v, b["lo"], b["d"]))
                boxes.append((b["lo"][0], b["lo"][1], v.shape[0] * b["d"][0], v.shape[1] * b["d"][1]))
        fine = min(p[2][0] for p in planes)
        frames.append(Frame(cps[prefixes[0]][it][0], rasterize(planes, extent, fine), boxes))
    return frames


def pulse_frames(work: Path) -> list[Frame]:
    """Return the frames of a PRISM pulse run in the cell plane just above z = 0: log10 of the field magnitude
    sqrt(|E|^2 + c^2 |B|^2) (the square root of twice the energy density over eps0), relative to its initial peak."""
    ini = read_ini(work / "input.ini")
    ngc = int(ini["grid"]["ngc"])
    g = ini["grid"]
    extent = (float(g["emin_x"]), float(g["emax_x"]), float(g["emin_y"]), float(g["emax_y"]))
    zmid = 0.5 * (float(g["emin_z"]) + float(g["emax_z"]))
    frames, scale = [], None
    for t, files in checkpoints(work).values():
        planes = []
        for b in read_blocks(files, ngc, ("Dx", "Dy", "Dz", "Bx", "By", "Bz")):
            nz = b["f"]["Dx"].shape[2]
            zc = b["lo"][2] + (np.arange(nz) + 0.5) * b["d"][2]
            k = np.nonzero((zc > zmid) & (zc - zmid < b["d"][2]))[0]
            if k.size:
                e2 = sum((b["f"][c][:, :, k[0]] / EPS0) ** 2 for c in ("Dx", "Dy", "Dz"))
                b2 = sum((C0 * b["f"][c][:, :, k[0]]) ** 2 for c in ("Bx", "By", "Bz"))
                mag = np.sqrt(e2 + b2)
                planes.append((mag, b["lo"], b["d"]))
        a = rasterize(planes, extent, min(p[2][0] for p in planes))
        scale = scale or float(np.nanmax(a))
        frames.append(Frame(t, np.log10(np.maximum(a / scale, 1e-30))))
    return frames


# --------------------------------------------------------------------------------------------------------- rendering


def match_times(panels: list[Panel], n_max: int) -> list[list[Frame]]:
    """Return the frame sequence of every panel on one time axis: the times of the first panel (at most n_max,
    evenly picked), the other panels at their checkpoint nearest in time."""
    times = [f.time for f in panels[0].frames]
    pick = np.unique(np.round(np.linspace(0, len(times) - 1, min(n_max, len(times)))).astype(int))
    seqs = []
    for p in panels:
        t_other = np.array([f.time for f in p.frames])
        seqs.append([p.frames[int(np.argmin(np.abs(t_other - times[i])))] for i in pick])
    return seqs


@dataclass
class Animation:
    """An animation: its panels side by side, the colour range and the time label format."""

    name: str
    panels: list[Panel]
    vmin: float | None = None
    vmax: float | None = None
    time_fmt: str = "t = {:.2f}"
    poster: float = 0.7  # the poster frame, as a fraction of the animation
    interp: str = "bilinear"
    gap: int = 6  # pixels between panels


def layout(anim: Animation, width: int, strip: int) -> tuple[int, int, list[tuple[int, int, int, int]]]:
    """Return the figure size in pixels and every panel's (x, y, w, h) in pixels, from the bottom left."""
    aspects = [(p.extent[1] - p.extent[0]) / (p.extent[3] - p.extent[2]) for p in anim.panels]
    avail = width - anim.gap * (len(anim.panels) - 1)
    h = int(round(avail / sum(aspects)))
    h += h % 2
    boxes, x = [], 0
    for i, a in enumerate(aspects):
        w = avail - sum(int(round(h * b)) for b in aspects[:-1]) if i == len(aspects) - 1 else int(round(h * a))
        boxes.append((x, strip, w, h))
        x += w + anim.gap
    height = h + strip
    return width, height + height % 2, boxes


def render(anim: Animation, frames_dir: Path, width: int, n_max: int) -> list[Path]:
    """Draw every frame in both themes into frames_dir/<theme>/NNNN.png; return the theme directories."""
    seqs = match_times(anim.panels, n_max)
    strip = 22
    w_px, h_px, boxes = layout(anim, width, strip)
    data = [f.data for seq in seqs for f in seq]
    vmin = anim.vmin if anim.vmin is not None else float(np.nanpercentile(np.concatenate([d.ravel() for d in data]),
                                                                         0.2))
    vmax = anim.vmax if anim.vmax is not None else float(np.nanpercentile(np.concatenate([d.ravel() for d in data]),
                                                                         99.8))
    dirs = []
    dpi = 100
    for theme in THEMES:
        out = frames_dir / theme.name
        out.mkdir(parents=True, exist_ok=True)
        dirs.append(out)
        cmap = theme.cmap.with_extremes(bad=theme.bg_alt)
        fig = plt.figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi, facecolor=theme.bg)
        axes, images = [], []
        for panel, (x, y, w, h), seq in zip(anim.panels, boxes, seqs, strict=True):
            ax = fig.add_axes((x / w_px, y / h_px, w / w_px, h / h_px))
            ax.set_axis_off()
            ax.set_facecolor(theme.bg)
            x0, x1, y0, y1 = panel.extent
            im = ax.imshow(seq[0].data.T, origin="lower", extent=panel.extent, cmap=cmap, vmin=vmin, vmax=vmax,
                           interpolation=anim.interp, aspect="auto")
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            if panel.decorate is not None:
                panel.decorate(ax, theme)
            # a hairline frame: where the field fades to the ground, it keeps the panels apart
            ax.add_patch(Rectangle((0.0, 0.0), 1.0, 1.0, transform=ax.transAxes, fill=False, ec=theme.mesh, lw=0.8,
                                   alpha=0.7, clip_on=False, zorder=5))
            axes.append(ax)
            images.append(im)
        label = fig.text(1.0 - 8 / w_px, (strip / 2) / h_px, "", ha="right", va="center", color=theme.text,
                         fontsize=9, family="DejaVu Sans")
        for n in range(len(seqs[0])):
            patches = []
            for ax, im, panel, seq in zip(axes, images, anim.panels, seqs, strict=True):
                im.set_data(seq[n].data.T)
                if panel.outline:
                    for bx, by, bw, bh in seq[n].boxes:
                        patches.append(ax.add_patch(Rectangle((bx, by), bw, bh, fill=False, lw=0.3,
                                                              ec=theme.mesh, alpha=0.45)))
            label.set_text(anim.time_fmt.format(seqs[0][n].time))
            fig.savefig(out / f"{n:04d}.png", dpi=dpi, facecolor=theme.bg)
            for p in patches:
                p.remove()
        plt.close(fig)
    return dirs


def encode(frames: Path, stem: Path, fps: int, poster: int) -> None:
    """Encode frames/NNNN.png to stem.webm (VP9) and stem.mp4 (H.264), copy the poster frame to stem.png."""
    src = ["-framerate", str(fps), "-i", str(frames / "%04d.png")]
    common = ["ffmpeg", "-y", "-loglevel", "error"]
    subprocess.run([*common, *src, "-c:v", "libvpx-vp9", "-crf", "36", "-b:v", "0", "-pix_fmt", "yuv420p",
                    "-row-mt", "1", "-deadline", "good", "-cpu-used", "2", "-an", f"{stem}.webm"], check=True)
    subprocess.run([*common, *src, "-c:v", "libx264", "-preset", "slow", "-crf", "26", "-pix_fmt", "yuv420p",
                    "-movflags", "+faststart", "-an", f"{stem}.mp4"], check=True)
    shutil.copyfile(frames / f"{poster:04d}.png", f"{stem}.png")


# ---------------------------------------------------------------------------------------------------- the animations


def step_decorate(ax, theme: Theme) -> None:  # noqa: ANN001
    """The step body and the two inter-realm seams of the forward-facing step."""
    ax.add_patch(Rectangle((0.6, 0.0), 2.4, 0.2, fill=True, fc=theme.bg_alt, ec=theme.mesh, lw=0.8, zorder=3))
    ax.plot([0.0, 0.6], [0.2, 0.2], ls=(0, (3, 2)), color=theme.accent, lw=0.9, zorder=4)  # seam A-B
    ax.plot([0.6, 0.6], [0.2, 1.0], ls=(0, (3, 2)), color=theme.accent, lw=0.9, zorder=4)  # seam B-C


def pml_decorate(width: float, extent: tuple) -> Callable:
    """Outline the PML layers on the x faces."""
    x0, x1, y0, y1 = extent

    def deco(ax, theme: Theme) -> None:  # noqa: ANN001
        for xa in (x0 + width, x1 - width):
            ax.plot([xa, xa], [y0, y1], ls=(0, (3, 2)), color=theme.mesh, lw=0.9)

    return deco


def anim_orszag_tang(runs: Path) -> Animation:
    """FLUME Orszag-Tang vortex, density."""
    names = ("r",)
    frames = density_frames(runs / "orszag-tang", [None], (0.0, 1.0, 0.0, 1.0), names)
    return Animation("orszag-tang", [Panel(frames, (0.0, 1.0, 0.0, 1.0), outline=True)], poster=1.0)


def anim_rayleigh_taylor(runs: Path) -> Animation:
    """FLUME Rayleigh-Taylor single mode, hydro and MHD, density."""
    ext = (0.0, 1.0, 0.0, 1.0)
    panels = [Panel(density_frames(runs / case, [None], ext, ("r",)), ext) for case in ("rt-hydro", "rt-mhd")]
    return Animation("rayleigh-taylor", panels, poster=0.9)


def anim_step(runs: Path) -> Animation:
    """FLUME forward-facing step on the three-realm forest, density."""
    ext = (0.0, 3.0, 0.0, 1.0)
    frames = density_frames(runs / "step", ["step-A", "step-B", "step-C"], ext, ("r",))
    return Animation("forward-facing-step", [Panel(frames, ext, outline=True, decorate=step_decorate)],
                     poster=1.0)


def anim_prism_faces(runs: Path) -> Animation:
    """PRISM Gaussian pulse against PEC, Silver-Muller and PML x faces, log10 of the field magnitude."""
    panels = []
    for face in ("pec", "sm", "pml"):
        work = runs / f"pulse-{face}"
        g = read_ini(work / "input.ini")["grid"]
        ext = (float(g["emin_x"]), float(g["emax_x"]), float(g["emin_y"]), float(g["emax_y"]))
        deco = pml_decorate(0.04, ext) if face == "pml" else None
        panels.append(Panel(pulse_frames(work), ext, decorate=deco))
    return Animation("prism-faces", panels, vmin=-5.0, vmax=0.15, time_fmt="t = {:.2e} s", poster=0.6, gap=8)


ANIMATIONS = {"orszag-tang": anim_orszag_tang, "rayleigh-taylor": anim_rayleigh_taylor,
              "forward-facing-step": anim_step, "prism-faces": anim_prism_faces}


def main() -> int:
    """Render and encode the selected animations."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs", type=Path, required=True, help="the runs directory of make_gallery_runs.sh")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "docs" / "public" / "gallery")
    ap.add_argument("--only", nargs="+", choices=tuple(ANIMATIONS), default=tuple(ANIMATIONS))
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--width", type=int, default=720)
    ap.add_argument("--frames", type=int, default=120, help="at most this many frames per animation")
    ap.add_argument("--keep-frames", type=Path, default=None, help="keep the PNG frames in this directory")
    args = ap.parse_args()
    if shutil.which("ffmpeg") is None:
        print("make_gallery.py: ffmpeg not found on PATH", file=sys.stderr)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    for name in args.only:
        anim = ANIMATIONS[name](args.runs)
        with tempfile.TemporaryDirectory() as tmp:
            base = args.keep_frames / name if args.keep_frames else Path(tmp)
            dirs = render(anim, base, args.width, args.frames)
            n = len(list(dirs[0].glob("*.png")))
            for d in dirs:
                encode(d, args.out / f"{name}-{d.name}", args.fps, min(n - 1, int(round(anim.poster * (n - 1)))))
        sizes = ", ".join(f"{p.name} {p.stat().st_size / 1024:.0f} KiB" for p in sorted(args.out.glob(f"{name}-*")))
        print(f"{name}: {n} frames; {sizes}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
