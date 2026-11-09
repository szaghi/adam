# Ideal MHD

The ideal-MHD verification of FLUME: the conservative MHD system with cell-centred $\mathbf{B}$ and mixed GLM cleaning,
the Riemann face fluxes (HLL, HLLD), the EGLM cleaning, the positivity limiter, the quadtrees and the isothermal and
barotropic MHD. Each case has a `check.sh` under
`src/tests/flume/verification/mhd/<case>/` that runs it (CPU by default, `FLUME_EXE=exe/adam_flume_fnl` for the GPU
backend) and an oracle that asserts the physics: an exact solution, a reference run, a symmetry or a conservation
property. All runs use 2 MPI ranks and WENO5 unless a page says otherwise; the models are described in
[models](../../models#ideal-mhd) and the schemes in [numerics](../../numerics).

Seven scripts (`rj2a`, `riemann`, `linear-wave`, `cpaw`, `orszag-tang`, `rotor`, `field-loop`) also run on the
`weno-riemann` scheme through the shared `--numerics` option ([MHD on `weno-riemann`](./weno-riemann)).

| Case | Type | What it verifies |
|---|---|---|
| [MHD plumbing](./plumbing) | C | the MHD model end to end: the fast-speed and $c_h$ time step, the auxiliary fields, a bitwise restart, the refused inputs |
| [MHD with a zero field](./zero-field) | C | MHD with $\mathbf{B} = 0$ reduces to Euler (Sod along x/y/z, bitwise across axes); a uniform MHD state is steady, bitwise on one level and to round-off across a 2:1 patch |
| [Ryu–Jones 2a: all seven waves](./ryu-jones-2a) | P | the whole MHD path against an exact 7-state solution; x/y/z bitwise; GLM inert in 1-D |
| [Brio–Wu and Ryu–Jones 4d](./brio-wu-rj4d) | P | clean runs near the eigensystem degeneracies; the positivity floors engage and stop as they should |
| [GLM divergence cleaning](./glm-cleaning) | P | the GLM telegraph system against its exact solution; div(B) decay of the Dedner peak, uniform and across seams; the `divb_tol` monitor |
| [Order of accuracy](./order-of-accuracy) | P | fifth order on linear waves, the circularly polarised Alfvén wave and the magnetised vortex |
| [Orszag–Tang vortex](./orszag-tang) | P | 180° symmetry, conservation and positivity, uniform and across 2:1 seams |
| [MHD rotor](./rotor) | P | symmetry under rotation with $\mathbf{B} \to -\mathbf{B}$, positivity |
| [Field loop advection](./field-loop) | P | the divergence error GLM leaves, the loop dissipation, CPU–FNL agreement on $\mathbf{B}$ |
| [Rotated shock tube](./rotated-shock-tube) | P | $B_\parallel$ preservation and agreement with the 1-D problem at 63.4° and 45° |
| [MHD conservation across AMR seams](./conservation) | C | every MHD integral, $\mathbf{B}$ and $\psi$ included, constant to round-off across 2:1 faces with reflux; the no-reflux negative control drifts |
| [MHD on `weno-riemann`](./weno-riemann) | P | the HLLD scheme on the verifications built for the flux splitting |
| [EGLM cleaning](./eglm) | P | the energy-consistent cleaning: telegraph solution, EGLM = GLM at $\psi = 0$, drifts bounded by the sources |
| [Positivity limiter](./positivity-limiter) | P | positivity on magnetised blasts and near vacuum, inactivity on smooth flow, 2:1 seams |
| [Quadtree AMR](./quadtree-amr) | C | 2:1 seams on quadtrees against the octree reference |
| [Isothermal and barotropic MHD](./barotropic-mhd) | P | ideal MHD without an energy equation: Kim's eigensystem and Mignone's isothermal HLLD on $P(\rho)$, fifth-order linear waves, Mignone's tubes, the inert energy slot |

The Ryu–Jones 2a split in two realms is on the [forest seams](../euler/multi-realm) page.

## Regenerating the figures

The `check.sh` scripts delete their checkpoints after the oracles, so the figures are drawn from re-runs:

```bash
cd src/tests/flume/verification
./make_doc_runs.sh <runs-dir> blast vacuum rj2a order mhd rotated quadtree   # every MHD figure run
./make_doc_figures.py --out ../../../../docs/public/flume --runs <runs-dir> \
   --only blast near-vacuum mhd-riemann order glm-pulse orszag-tang rotor field-loop eglm-energy \
          weno-riemann rotated-shock-tube quadtree
./make_doc_runs.sh <runs-dir> eosmhd && ./make_doc_figures.py --out ../../../../docs/public/flume --runs <runs-dir> \
   --only barotropic-mhd   # the isothermal and barotropic MHD (M13)
```

`eglm-energy` reads the conservation histories kept in the `mhd/orszag-tang` work directories. The group `mhd` writes
the runs of the gallery figures: Brio–Wu, RJ4d, the GLM pulse, Orszag–Tang and its AMR run, the rotor, the field loop,
and RJ2a on HLLD with primitive interpolation.

## How it was built

The conservative MHD system with mixed GLM cleaning is milestone M2 ([#41](https://github.com/szaghi/adam/issues/41));
the Riemann face fluxes, the EGLM cleaning and the positivity limiter are M3
([#47](https://github.com/szaghi/adam/issues/47)); the quadtrees are [#46](https://github.com/szaghi/adam/issues/46);
the isothermal and barotropic MHD is M13 ([#100](https://github.com/szaghi/adam/issues/100)). The Ryu–Jones 2a split in
two realms on the forest seams page is MV-14 (M2 · P7c).

| Case | ID | Milestone |
|---|---|---|
| [MHD plumbing](./plumbing) | — (`mhd/plumbing`) | M2 · P1 |
| [MHD with a zero field](./zero-field) | MV-2 | M2 · P3a |
| [Ryu–Jones 2a: all seven waves](./ryu-jones-2a) | MV-4 | M2 · P3b |
| [Brio–Wu and Ryu–Jones 4d](./brio-wu-rj4d) | — (`mhd/riemann`) | M2 · P3c |
| [GLM divergence cleaning](./glm-cleaning) | MV-3, MV-10 | M2 · P4; M6 · G1 |
| [Order of accuracy](./order-of-accuracy) | MV-5, MV-6, MV-7 | M2 · P5 |
| [Orszag–Tang vortex](./orszag-tang) | MV-12 | M2 · P6a |
| [MHD rotor](./rotor) | MV-13 | M2 · P6b |
| [Field loop advection](./field-loop) | MV-9 | M2 · P6c |
| [Rotated shock tube](./rotated-shock-tube) | MV-8 | M2 · P6d |
| [MHD conservation across AMR seams](./conservation) | MV-11 | M2 · P7a |
| [MHD on `weno-riemann`](./weno-riemann) | RV-5, RV-6, RV-7 | M3 · P3c |
| [EGLM cleaning](./eglm) | EV-2, EV-3, EV-4 | M3 · P4 |
| [Positivity limiter](./positivity-limiter) | PV-1 to PV-5 | M3 · P5 |
| [Quadtree AMR](./quadtree-amr) | MV-15 | #46 |
| [Isothermal and barotropic MHD](./barotropic-mhd) | EO (em) | M13 · P5a, P5b |

The [verification catalog](../) places each page in implementation order.
