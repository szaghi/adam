# Runtime AMR

Runtime AMR regrids FLUME during the time loop. `[amr] frequency = 0` keeps the AMR of the initial condition only;
`n > 0` regrids every `n` steps. The library already holds a complete regrid step (`adam_object%amr_update`: the tree
adapts with 2:1 balance, the interior data is prolonged or restricted, the blocks are redistributed, the maps rebuilt).
On top of it FLUME adds a conservative prolongation (`[amr] regrid_prolongation`, default `conservative`), the markers
combined, the Löhner estimator, the solids' distance function and the flux register rebuilt, restarts across regrids
([numerics](../../numerics#adaptive-mesh-refinement)), and the FNL backend regridded by a host round trip. The cases
verify the regrid step in isolation, conservation and bitwise restart through regrids, CPU against FNL, and the accuracy
of tracked runs against uniform-fine runs. The input keys are described in
[setting up runtime AMR](../../input#setting-up-runtime-amr).

| Case | Type | What it verifies |
|---|---|---|
| [Regrid round trip and input contract](./regrid-round-trip) | C | the library regrid in isolation: bookkeeping over 1, 2, 3 ranks, the conservative prolongation exact on the integral and in range, the `frequency` and `regrid_prolongation` contracts |
| [Runtime regrids in FLUME](./runtime-regrids) | C | conservation through regrids, bitwise restart across regrids, the immersed solid, CPU against FNL |
| [Accuracy of runtime AMR](./accuracy) | P | tracked runs against uniform-fine runs on the vortex, Sod, the MHD blast, the field loop and the cylinder; rank and backend invariance |
| [MHD blast symmetry on the flux-splitting path](./mhd-blast-symmetry) | P | why the low-β blast loses its point symmetry on `weno` flux splitting, and the Riemann-flux path that keeps it |
| [Runtime AMR on a forest](./multi-realm) | C | the realms of a forest regridding across their seams: the cross-realm balance, mixed levels along a seam, restart across regrids |

The Jeans-marker contract is on [the polytrope on AMR](../self-gravity/star-amr); the markers FLUME does not implement
are refused at initialisation (leg rg8 of `regrid/check.sh`).

Run every leg:

```bash
cd src/tests/flume/verification/regrid
./check.sh && ./accuracy.sh                                                           # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs rg3,rg4,rg5
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./accuracy.sh                        # FNL
./check.sh --legs rg6 && ./accuracy.sh --legs av8                                     # CPU against FNL
```

The figure (`regrid.png`) is drawn by `make_doc_figures.py` from the CPU work directories `av2-tracked` and
`av4-glm-tracked` that `regrid/accuracy.sh` leaves (no `make_doc_runs.sh` group):

```bash
exe/.regression-venv/bin/python src/tests/flume/verification/make_doc_figures.py --out docs/public/flume --only regrid
```

## How it was built

Runtime AMR is milestone M5 ([#74](https://github.com/szaghi/adam/issues/74)), built in phases:

| Phase | Content |
|---|---|
| P0 | tests the library regrid step on its own and puts the hooks in place |
| P1 | the conservative prolongation FLUME regrids with (`[amr] regrid_prolongation`, default `conservative`) |
| P2 | regrids FLUME on the CPU: the markers combined, the Löhner estimator, the solids' distance function and the flux register rebuilt, restarts across regrids ([numerics](../../numerics#adaptive-mesh-refinement)) |
| P3 | regrids the FNL backend by a host round trip |
| P4 | the accuracy of runtime AMR (AV) |

| Case | ID | Milestone |
|---|---|---|
| [Regrid round trip and input contract](./regrid-round-trip) | RG (rg0 to rg2) | M5 · P0, P1 |
| [Runtime regrids in FLUME](./runtime-regrids) | RG (rg3 to rg6) | M5 · P2, P3 |
| [Accuracy of runtime AMR](./accuracy) | AV (av2 to av8) | M5 · P4 |
| [MHD blast symmetry on the flux-splitting path](./mhd-blast-symmetry) | av3c | [#77](https://github.com/szaghi/adam/issues/77) |
| [Runtime AMR on a forest](./multi-realm) | multirealm legs 7–13 | [#110](https://github.com/szaghi/adam/issues/110) · P1–P5 |

Two later legs of `regrid/check.sh` belong to other issues: rg7, the Jeans-marker contract of M9, and rg8
([#110](https://github.com/szaghi/adam/issues/110) P4.2). The forest page comes from #110, after M13, and the
[verification catalog](../) lists it there.
