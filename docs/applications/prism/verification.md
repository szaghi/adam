# Verification gallery

PRISM's regression suite (`src/tests/prism/regression/`, [PRISM regression](/tests/prism-regression)) pins every case
to a golden digest; this page shows the cases that also carry an **oracle**: a `check.sh` that asserts a physical or
structural property of the run (an invariant, a symmetry, an agreement between two discretisations of the same
problem) rather than a stored number. The figures are drawn from the runs by
`src/tests/prism/regression/make_doc_figures.py` (each check run with `--keep`); the runs use the CPU backend unless
stated, `PRISM_EXE=exe/adam_prism_fnl` runs a check on the GPU (FNL) backend.

```bash
cd src/tests/prism/regression/rmf-amr-fd-pulse && ./check.sh       # CPU backend
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL backend (GPU)
```

## AMR

### Quadtree 2:1 seams: the source-free pulse

A quadtree (`[amr] ratio = 4`) refines x and y only: its blocks span the domain along z at every level, so a 2:1 seam
is 2:1 in x and y and 1:1 in z. Until [#46](https://github.com/szaghi/adam/issues/46) the seam machinery (ghost
restriction, coarse-fine interpolation, reflux) assumed an octree, and the coarse cells beside a quadtree seam picked up
a spurious z dependence; quadtree AMR was refused at initialisation.

**Case.** `rmf-amr-fd-pulse`: the source-free Gaussian EM pulse (`B0 = 10`, σ = 0.01) travelling along +y across a 2:1
seam at x = 0 (the x < 0 half refined one level), `fd_centered` order 6, SSP-54, 5 steps, 1 rank. The pulse is uniform
in x and z, so the seam is the only thing that can perturb it: it generates By, and the seam div(B) source of
[#29](https://github.com/szaghi/adam/issues/29) (a truncation-order quantity the check pins to a baseline).

**Oracle** (the quadtree leg of `check.sh`, default-on). The same case runs on the octree of the regression case and on
a quadtree with `nk = 16` (z cells of 7.5·10⁻³ against the octree's finest 1.875·10⁻³). PRISM's time step is
`CFL · min(dx, dy, dz) / c`, set on both trees by the finest x and y cells, so the two runs take the same steps; the
pulse is uniform in z, so both stay z-invariant. `quadtree_oracle.py` (shared with FLUME's MV-15) keys the cells by their
(x, y) centre and asserts, per field:

- the z spread of every (x, y) column of each run at most `1e-13`;
- the quadtree against the octree at the same (x, y) at most `1e-13`.

Each field is measured against the largest magnitude of its vector (D, B, J): the pulse carries Dz and Bx, and the seam
adds By; the other components are round-off on the octree (Dx, Dy ~ 3·10⁻¹⁸, Bz ~ 5·10⁻¹⁶ against |Bx| = 10) and exactly
zero on the quadtree, and scaled by their own magnitude they would compare round-off with round-off. The leg also holds
the quadtree's div(D) to the structural round-off invariant of the case.

**Results** (CPU; FNL in parentheses where it differs):

| Quantity | Octree | Quadtree |
|---|---|---|
| z spread of a column, max over D and B | 2.6·10⁻¹⁶ (FNL 3.6·10⁻¹⁶) | 0 |
| against the octree, per (x, y) column, max over D and B | — | 1.8·10⁻¹⁶ (FNL 3.6·10⁻¹⁶) |
| seam max\|div(B)\| at step 5 | 12.65004 | 12.65004 |
| max\|div(D)\| | 3.1·10⁻¹⁵ | 0 |

The quadtree reproduces the octree to round-off, including the seam's div(B) source to every printed digit: on a
z-invariant field the per-axis seam formulas (restriction over 2×2 fine cells, coincident interpolation rows along z)
give the same values as the octree's 2×2×2 ones on two identical z layers. The octree's own z spread and its non-zero
div(D) are round-off of its z interpolation (its tricubic weights depend on the z sub-position); the quadtree has
neither.

![PRISM pulse across a 2:1 seam, octree against quadtree](/prism/quadtree-pulse.png)

Top row: Bx on z = 0 (quadtree, blocks outlined: the refined half x < 0 and the seam at x = 0), and Bx on y = 0 for
both trees: the octree refines its x < 0 blocks in z as well, the quadtree keeps one z layer of blocks. Bottom row: the
quadtree against the octree per (x, y) column (relative to each vector's magnitude); By, the field the seam generates;
and the divergence histories of the two runs, superposed.
