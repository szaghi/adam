# External gravity (M7)

External gravity lands in phases ([#73](https://github.com/szaghi/adam/issues/73), plan in the issue): G0 the inputs,
G1 the force, G2 a well-balanced mode (G2a the order of the energy source, G2b the well-balanced mode itself), G3
Rayleigh–Taylor. Self-gravity is M8 ([#81](https://github.com/szaghi/adam/issues/81), [Self-gravity](../self-gravity/)).
The model is described in [models](../../models#external-gravity), the input in
[input](../../input#gravity-optional-external-gravity-and-self-gravity).

Every case of the group is written by one generator, `src/tests/flume/verification/gravity/make_atmosphere.py`, and run
and judged by `src/tests/flume/verification/gravity/check.sh` with the oracles `atmosphere_oracle.py` and
`rt_oracle.py` beside it. Every input has four ghost cells (`ngc = 4`, the stencil of the fourth-order gravity energy
source with WENO5) and two uniform refinement levels, so that the blocks spread over the ranks.

The figures are drawn from re-runs of the `check.sh` inputs (CPU, 2 ranks), whose checkpoints the checks delete:

```bash
V=src/tests/flume/verification                       # from the repository root
$V/make_doc_runs.sh <dir> gravity                    # G0, G1, G2a, G2b, GV-4 (23 runs, about 27 min)
$V/make_doc_runs.sh <dir> rt                         # the nonlinear Rayleigh-Taylor runs, N = 128 to t = 5 (2 runs)
FLUME_EXE=$PWD/exe/adam_flume_fnl $V/make_doc_runs.sh <dir> rt-sym   # rt2d-mhd under other cleanings (#83, 8 runs)
python3 $V/make_doc_figures.py --out docs/public/flume --runs <dir> --only gravity-atmospheres gravity-force \
   gravity-order gravity-balance gravity-rt-linear gravity-rt gravity-rt-symmetry
```

`gravity-rt` needs the `rt` group; `gravity-rt-symmetry` needs both `rt` and `rt-sym`.

| Case | ID | What it verifies |
|---|---|---|
| [Stratified atmospheres and the `hydrostatic` boundary](./atmospheres) | GV (G0) | the initial atmospheres and the ghosts of the `hydrostatic` faces against the closed form, before any force acts |
| [The gravity force](./force) | GV-1, GV-2, GV-3, GV-6, GV-7 (G1) | exact free fall, mass and total energy of a closed box to round-off, the limiter near vacuum, scaling covariance |
| [The order of the energy source](./energy-order) | GV-9 (G2a) | the fourth-order energy source restores the order of the scheme on a falling density wave |
| [The well-balanced mode](./well-balanced) | GV-2b, GV-5, GV-8, GV-10 (G2b) | an atmosphere at rest stays at rest to round-off, across 2:1 seams too, and small perturbations respond linearly |
| [Rayleigh–Taylor](./rayleigh-taylor) | GV-4 (G3) | the linear growth rate against the compressible linear theory; the nonlinear stage and the mirror symmetry of the MHD run (#83) |
