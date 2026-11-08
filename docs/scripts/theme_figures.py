#!/usr/bin/env python3
"""Put the docs' matplotlib figures on the washi ground of the sumi identity (docs/BRAND.md).

The figures are drawn on white; on the washi page they read as pasted white cards. This replaces the figure
background with the washi colour, and nothing else: the background is the white region reachable from the image
border (a flood fill through near-white pixels), so the axes frames, the curves and the field images bound it and
are never changed. In the dark theme the same image is shown as a washi card (custom.css).

Usage: theme_figures.py [--check] PNG...   (in place; --check lists the files that still have a white border)
Needs numpy and Pillow.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

WASHI = (0xF3, 0xED, 0xE1)


def white_border(img: Image.Image) -> bool:
    a = np.asarray(img.convert("RGB"))
    edge = np.concatenate([a[0], a[-1], a[:, 0], a[:, -1]])
    return bool((edge.min(axis=1) >= 250).mean() > 0.5)


def to_washi(path: Path) -> bool:
    """Recolour the border-connected white background of one PNG; return whether it changed."""
    img = Image.open(path)
    mode = img.mode
    rgb = img.convert("RGB")
    a = np.asarray(rgb)
    near_white = (a.min(axis=2) >= 245).astype(np.uint8) * 255  # pure and anti-aliased paper white
    mask = Image.fromarray(near_white, "L").copy()  # fromarray shares the read-only buffer: floodfill would do nothing
    h, w = near_white.shape
    seeds = [(x, 0) for x in range(0, w, 4)] + [(x, h - 1) for x in range(0, w, 4)]
    seeds += [(0, y) for y in range(0, h, 4)] + [(w - 1, y) for y in range(0, h, 4)]
    for xy in seeds:
        if mask.getpixel(xy) == 255:
            ImageDraw.floodfill(mask, xy, 128)
    bg = np.asarray(mask) == 128
    if not bg.any():
        return False
    out = a.copy()
    out[bg] = WASHI
    res = Image.fromarray(out, "RGB")
    if mode == "RGBA":
        res = res.convert("RGBA")
    res.save(path, optimize=True)
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("png", nargs="+", type=Path)
    ap.add_argument("--check", action="store_true", help="list the files that still have a white border, change nothing")
    a = ap.parse_args()
    bad = 0
    for p in a.png:
        if a.check:
            if white_border(Image.open(p)):
                print(p)
                bad += 1
        elif to_washi(p):
            print(f"washi: {p}")
    return 1 if (a.check and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
