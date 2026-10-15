# Ideal MHD

The ideal-MHD verification of FLUME: the conservative MHD system with cell-centred $\mathbf{B}$ and mixed GLM cleaning
of M2 ([#41](https://github.com/szaghi/adam/issues/41)), and the Riemann face fluxes (HLL, HLLD), the EGLM cleaning and
the positivity limiter of M3 ([#47](https://github.com/szaghi/adam/issues/47)). Each case has a `check.sh` under
`src/tests/flume/verification/mhd/<case>/` that runs it (CPU by default, `FLUME_EXE=exe/adam_flume_fnl` for the GPU
backend) and an oracle that asserts the physics: an exact solution, a reference run, a symmetry or a conservation
property. All runs use 2 MPI ranks and WENO5 unless a page says otherwise; the models are described in
[models](../../models#ideal-mhd) and the schemes in [numerics](../../numerics).

Seven scripts (`rj2a`, `riemann`, `linear-wave`, `cpaw`, `orszag-tang`, `rotor`, `field-loop`) also run on the
`weno-riemann` scheme through the shared `--numerics` option ([MHD on `weno-riemann`](./weno-riemann)).

| Case | ID | What it verifies |
|---|---|---|
| [Ryu–Jones 2a: all seven waves](./ryu-jones-2a) | MV-4 | the whole MHD path against an exact 7-state solution; x/y/z bitwise; GLM inert in 1-D |
| [Brio–Wu and Ryu–Jones 4d](./brio-wu-rj4d) | — (`mhd/riemann`) | clean runs near the eigensystem degeneracies; the positivity floors engage and stop as they should |
| [GLM divergence cleaning](./glm-cleaning) | MV-3, MV-10 | the GLM telegraph system against its exact solution; div(B) decay of the Dedner peak, uniform and across seams; the `divb_tol` monitor |
| [Order of accuracy](./order-of-accuracy) | MV-5, MV-6, MV-7 | fifth order on linear waves, the circularly polarised Alfvén wave and the magnetised vortex |
| [Orszag–Tang vortex](./orszag-tang) | MV-12 | 180° symmetry, conservation and positivity, uniform and across 2:1 seams |
| [MHD rotor](./rotor) | MV-13 | symmetry under rotation with $\mathbf{B} \to -\mathbf{B}$, positivity |
| [Field loop advection](./field-loop) | MV-9 | the divergence error GLM leaves, the loop dissipation, CPU–FNL agreement on $\mathbf{B}$ |
| [Rotated shock tube](./rotated-shock-tube) | MV-8 | $B_\parallel$ preservation and agreement with the 1-D problem at 63.4° and 45° |
| [MHD on `weno-riemann`](./weno-riemann) | RV-5, RV-6, RV-7 | the HLLD scheme on the M2 verifications |
| [EGLM cleaning](./eglm) | EV-2, EV-3, EV-4 | the energy-consistent cleaning: telegraph solution, EGLM = GLM at $\psi = 0$, drifts bounded by the sources |
| [Positivity limiter](./positivity-limiter) | PV-1 to PV-5 | positivity on magnetised blasts and near vacuum, inactivity on smooth flow, 2:1 seams |
| [Quadtree AMR](./quadtree-amr) | MV-15 | 2:1 seams on quadtrees against the octree reference |

Other MHD legs: the Ryu–Jones 2a split in two realms (MV-14, `mhd/multirealm`) is on the
[multi-realm](../euler/multi-realm) page. `mhd/plumbing` (M2-P1, the MHD model end to end with a zero residual),
`mhd/zero-field` (MV-2, MHD with $\mathbf{B} = 0$ equals Euler, a uniform state across a 2:1 patch) and
`mhd/conservation` (MV-11, MHD conservation across AMR faces with reflux) have no page yet
([#85](https://github.com/szaghi/adam/issues/85)).

## Regenerating the figures

The `check.sh` scripts delete their checkpoints after the oracles, so the figures are drawn from re-runs:

```bash
cd src/tests/flume/verification
./make_doc_runs.sh <runs-dir> blast vacuum rj2a order    # PV-1, PV-2, PV-3, Ryu-Jones 2a, the MV-5 ladders
./make_doc_figures.py --out ../../../../docs/public/flume --runs <runs-dir> \
   --only blast near-vacuum mhd-riemann order glm-pulse orszag-tang rotor field-loop eglm-energy
```

`eglm-energy` reads the conservation histories kept in the `mhd/orszag-tang` work directories. The `--runs`
sub-directories `orszag-tang`, `orszag-tang-amr`, `rotor`, `field-loop`, `brio-wu`, `rj4d` and `glm-pulse` are not
written by `make_doc_runs.sh`: their figures cannot be regenerated from the scripts alone
([#85](https://github.com/szaghi/adam/issues/85)).
