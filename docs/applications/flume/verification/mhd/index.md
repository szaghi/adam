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
| [MHD plumbing](./plumbing) | — (M2-P1) | the MHD model end to end: the fast-speed and $c_h$ time step, the auxiliary fields, a bitwise restart, the refused inputs |
| [MHD with a zero field](./zero-field) | MV-2 | MHD with $\mathbf{B} = 0$ reduces to Euler (Sod along x/y/z, bitwise across axes); a uniform MHD state is steady, bitwise on one level and to round-off across a 2:1 patch |
| [GLM divergence cleaning](./glm-cleaning) | MV-3, MV-10 | the GLM telegraph system against its exact solution; div(B) decay of the Dedner peak, uniform and across seams; the `divb_tol` monitor |
| [Order of accuracy](./order-of-accuracy) | MV-5, MV-6, MV-7 | fifth order on linear waves, the circularly polarised Alfvén wave and the magnetised vortex |
| [MHD conservation across AMR seams](./conservation) | MV-11 | every MHD integral, $\mathbf{B}$ and $\psi$ included, constant to round-off across 2:1 faces with reflux; the no-reflux negative control drifts |
| [Orszag–Tang vortex](./orszag-tang) | MV-12 | 180° symmetry, conservation and positivity, uniform and across 2:1 seams |
| [MHD rotor](./rotor) | MV-13 | symmetry under rotation with $\mathbf{B} \to -\mathbf{B}$, positivity |
| [Field loop advection](./field-loop) | MV-9 | the divergence error GLM leaves, the loop dissipation, CPU–FNL agreement on $\mathbf{B}$ |
| [Rotated shock tube](./rotated-shock-tube) | MV-8 | $B_\parallel$ preservation and agreement with the 1-D problem at 63.4° and 45° |
| [MHD on `weno-riemann`](./weno-riemann) | RV-5, RV-6, RV-7 | the HLLD scheme on the M2 verifications |
| [EGLM cleaning](./eglm) | EV-2, EV-3, EV-4 | the energy-consistent cleaning: telegraph solution, EGLM = GLM at $\psi = 0$, drifts bounded by the sources |
| [Positivity limiter](./positivity-limiter) | PV-1 to PV-5 | positivity on magnetised blasts and near vacuum, inactivity on smooth flow, 2:1 seams |
| [Quadtree AMR](./quadtree-amr) | MV-15 | 2:1 seams on quadtrees against the octree reference |
| [Isothermal and barotropic MHD](./barotropic-mhd) | EO (em) | ideal MHD without an energy equation: Kim's eigensystem and Mignone's isothermal HLLD on $P(\rho)$, fifth-order linear waves, Mignone's tubes, the inert energy slot |

The Ryu–Jones 2a split in two realms (MV-14, `mhd/multirealm`) has no page of its own: it is on the
[multi-realm](../euler/multi-realm) page.

## Regenerating the figures

The `check.sh` scripts delete their checkpoints after the oracles, so the figures are drawn from re-runs:

```bash
cd src/tests/flume/verification
./make_doc_runs.sh <runs-dir> blast vacuum rj2a order mhd rotated quadtree   # every MHD figure run
./make_doc_figures.py --out ../../../../docs/public/flume --runs <runs-dir> \
   --only blast near-vacuum mhd-riemann order glm-pulse orszag-tang rotor field-loop eglm-energy \
          weno-riemann rotated-shock-tube quadtree
./make_doc_runs.sh <runs-dir> eosmhd && ./make_doc_figures.py --out ../../../../docs/public/flume --runs <runs-dir> \
   --only barotropic-mhd   # the isothermal and barotropic MHD (#100 P5)
```

`eglm-energy` reads the conservation histories kept in the `mhd/orszag-tang` work directories. The group `mhd` (issue #85, R2a) writes the runs of the gallery figures: Brio–Wu, RJ4d, the GLM
pulse, Orszag–Tang and its AMR run, the rotor, the field loop, and RJ2a on HLLD with primitive interpolation.
