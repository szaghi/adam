# Visual identity: the sumi family (ADAM · FLUME · PRISM)

This is the source of truth for the look of the ADAM docs and its two applications. Read it before touching logos, favicons, colours or the docs theme.

## The idea

Three marks drawn in one hand, in the manner of sumi-e ink painting: a few loaded brush strokes, a seal, and a brushed AMR (adaptive mesh refinement) ground underneath.

The ground is never decoration. It is a balanced quadtree that refines one level deeper exactly under every stroke, so the mesh is finest where the ink is. The finest grid strokes are also the thinnest.

| Mark | Project | Image | What it says |
|---|---|---|---|
| `adam/` | ADAM, the SDK | Bow shock ahead of a body (ensō) | AMR and immersed boundaries, the library's foundations |
| `flume/` | FLUME | Black hole: tilted disk with two rings, lensed arcs doubled, photon ring | Compressible Euler and MHD, gravity, sinks |
| `prism/` | PRISM | RMF/helicon plasma thruster: 4 coils, field lines threading the first coil, nozzle, one central stroke | Maxwell and particle-in-cell plasmas |

All three share the seal (a square carved with a refined quadrant), the brush and the ground.

## Themes

| Token | Light · washi | Dark · silver ink |
|---|---|---|
| Ground (bg) | `#f3ede1` | `#0e1116` |
| bg-alt | `#ebe3d3` | `#0a0d11` |
| Ink / text-1 | `#1d1b18` | `#e8ecf0` |
| text-2 | `#5f5a52` | `#a9b1bb` |
| Mesh | `#6f6a62` | `#8b949e` |
| Links (brand-1) | `#b23417` | `#ffc24d` |
| Seal (brand-2) | `#c73e1d` | `#ffb020` |
| brand-3 | `#8f2a12` | `#c98500` |

Links on washi use `#b23417`, darker than the seal, because the seal colour is under 4.5:1 contrast on the paper. Keep that split.

Typography: body in Inter (as on the other szaghi.github.io sites). Display and hero name in Shippori Mincho 800, letter-spacing 0.04em.

## Files

```
adam/  flume/  prism/
  logo.svg            washi version (light theme), 64×64 viewBox
  logo-dark.svg       silver-ink version (dark theme)
  favicon.svg         simplified: no mesh, strokes thickened, on a washi tile
  favicon-dark.svg    same on a slate tile
  favicon-16.png  favicon-32.png  apple-touch-icon.png (180)
vitepress/
  custom.css          theme tokens for both modes + hero font
  index.ts            loads custom.css into the default theme
  config-snippet.ts   head links (favicons, theme-color) and logo light/dark
```

Install: copy `adam/*` into `docs/public/`, and `flume/` and `prism/` as subfolders of `docs/public/`. Copy `vitepress/custom.css` and `vitepress/index.ts` into `docs/.vitepress/theme/`, then merge `config-snippet.ts` into the VitePress config.

## Rules

- Use `logo.svg` on light grounds and `logo-dark.svg` on dark grounds. Never put the washi logo on a dark ground or the reverse.
- On washi, FLUME's shadow is solid ink. On silver ink it is the ground itself. This difference is intentional; do not "fix" it.
- Below 48 px use the favicon variants; the full logos lose their mesh at small sizes.
- No gradients, glows or drop shadows behind the marks (`--vp-home-hero-image-background-image: none`).
- The SVGs are generated geometry (quadtree plus tapered-brush polygons). Do not hand-edit paths; regenerate if the drawing must change.
- Rejected alternatives (other palettes such as aizome, urushi and sepia; the accretion-disk FLUME; the flux-surface and Hall-thruster PRISM) are archived and are not part of the identity.
