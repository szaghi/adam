# Runtime AMR (M5)

Milestone M5 ([#74](https://github.com/szaghi/adam/issues/74)) regrids during the time loop. `[amr] frequency = 0`, the
value of every input so far, keeps the AMR of the initial condition only; `n > 0` regrids every `n` steps. The library
already holds a complete regrid step (`adam_object%amr_update`: the tree adapts with 2:1 balance, the interior data is
prolonged or restricted, the blocks are redistributed, the maps rebuilt).

| Phase | Content |
|---|---|
| P0 | tests the library regrid step on its own and puts the hooks in place |
| P1 | the conservative prolongation FLUME regrids with (`[amr] regrid_prolongation`, default `conservative`) |
| P2 | regrids FLUME on the CPU: the markers combined, the Löhner estimator, the solids' distance function and the flux register rebuilt, restarts across regrids ([numerics](../../numerics#adaptive-mesh-refinement)) |
| P3 | regrids the FNL backend by a host round trip |
| P4 | the accuracy of runtime AMR (AV) |

The input keys are described in [setting up runtime AMR](../../input#setting-up-runtime-amr).

| Case | ID | What it verifies |
|---|---|---|
| [Regrid round trip and input contract](./regrid-round-trip) | RG (rg0 to rg2) | the library regrid in isolation: bookkeeping over 1, 2, 3 ranks, the conservative prolongation exact on the integral and in range, the `frequency` and `regrid_prolongation` contracts |
| [Runtime regrids in FLUME](./runtime-regrids) | RG (rg3 to rg6) | conservation through regrids, bitwise restart across regrids, the immersed solid, CPU against FNL |
| [Accuracy of runtime AMR](./accuracy) | AV (av2 to av8) | tracked runs against uniform-fine runs on the vortex, Sod, the MHD blast, the field loop and the cylinder; rank and backend invariance |
| [MHD blast symmetry on the flux-splitting path](./mhd-blast-symmetry) | #77, av3c | why the low-β blast loses its point symmetry on `weno` flux splitting, and the Riemann-flux path that keeps it |
| [Runtime AMR on a forest](./multi-realm) | multirealm legs 7–13 | the realms of a forest regridding across their seams: the cross-realm balance, mixed levels along a seam, restart across regrids ([#110](https://github.com/szaghi/adam/issues/110)) |

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
