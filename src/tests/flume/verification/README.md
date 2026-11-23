# FLUME verification: figures

`make_doc_figures.py` draws the figures of the FLUME documentation (`docs/applications/flume/`) and of the ADAM paper
from the checkpoints of the verification runs. The runs whose `check.sh` deletes the checkpoints are re-run into a
runs directory by `make_doc_runs.sh <runs-dir> <group>...`; the others are read from the `work-*` directories of
`check.sh --keep` in this tree.

## Documentation figures

```bash
make_doc_figures.py --out <dir> [--runs <runs-dir>] [--only NAME ...] [--theme light|dark|both]
```

One `<name>.png` per figure on the washi ground (`--theme dark`: `<name>-dark.png`, see
`docs/scripts/figure_theme.py`). `--only paper` draws the docs versions of the paper's figures.

## Paper figures (issue #142)

```bash
make_doc_figures.py --out <dir> --style paper --runs <runs-dir> [--only paper | NAME ...] [--usetex]
```

`paper_figures()` lists the paper's figures. Each one has a physical file name (`order-vortex`, `orszag-tang`,
`step-3d-amr`, ...) and is written as `<name>.pdf` and `<name>.png`. The docs figure functions draw them, and
`paper_style.py` restyles them when they are saved:

- white ground; no figure title; panel titles removed, or replaced by short ones set in the spec;
- panel labels (a), (b), ... in the left title slot of every panel;
- the code identifiers mapped to scheme and physical names through the one dictionary `paper_style.LABELS`
  (`LABEL_RULES` covers the patterns). The audit lists any identifier that is left, and any colorbar without a
  label;
- the size of the printed page (`figure_theme.PAGE_WIDTHS`: full = `\textwidth + \extralength` for
  `adjustwidth{-\extralength}` with `\linewidth`, 7.27 in; text = `\textwidth`, 5.46 in), every text at 8 pt and
  the panel labels at 9 pt;
- perceptually uniform colour maps for unsigned data, a centred diverging map only for signed data, the Okabe-Ito
  palette, and markers on labelled lines (but not on reference lines);
- field maps rasterized at 600 dpi inside the vector PDF. With an OpenType-CFF font (TeX Gyre Pagella) the PDF
  fonts are Type 3: matplotlib 3.11 writes Type 42 OpenType-CFF fonts in a form that poppler does not render.

The rcParams block is `figure_theme.PAPER_RC` (`figure_theme.paper()`). It is the only one, and the PRISM figures
will use it too.

`figures.json`, next to the figures, is the manifest. It gives the ADAM commit of the script run and whether `src/`
or `docs/scripts/` differ from it. For each figure it gives its files, its LaTeX include, the run directories it
read, and the provenance of each run: `run.info` if present, otherwise the keys found in the log header, otherwise
`unknown`. The executable named by the work directory and the time of the newest checkpoint are added, each
labelled as such. FLUME logs do not record the commit, the host or the device today. Until the runs write a
`run.info`, those fields stay `unknown`.

Expected run directories (`<runs-dir>/...`): `orszag-tang`, `orszag-tang-amr`, `rotor`, `linear-wave-*`,
`blast77-{split,hlld}`, `regrid-vortex-{conservative,linear}`, `sv7-l{2,3}`, `sv8-collapse-*`, `sv14-*`,
`si-binary`, `si-shu-{32,64}`, `si-collapse-amr`, and `step-3d`. `step-3d` is the 3-D step:
`step/make_step.py <dir> --layout aligned --dims 3 --runtime ...`, with `step-{A,B,C}.ini` and the realms'
checkpoints in the same directory. The vortex and shock-cylinder figures read the `work-*` directories of
`vortex/` and `riemann-flux/`.
