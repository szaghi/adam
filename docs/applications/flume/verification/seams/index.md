# Accuracy and divergence across 2:1 seams

The static 2:1 seams of an AMR grid: what they cost in order of accuracy, in div(B) and at a crossing shock, and the
seam operators that lift the accuracy cap. The cases measure the order of a smooth vortex per region against the uniform
control, the divergence the cleaning acts on, a shock crossing a 2:1 face against the exact solution, and the drift of
the integrals without the reflux; then the 6-point restriction, the high-order and limited coarse-to-fine fills, the
sensor-gated reflux and the EGLM reflux. The transfer operators and the reflux policy change only behind keys; the
defaults are `mean`, `tricubic` and `always`:

| Key | Values (default first) |
|---|---|
| `[amr] seam_restriction` | `mean`, `lagrange6`, `weno6` |
| `[amr] seam_ghost_fill` | `tricubic`, `injection`, `restriction-compatible`, `quintic`, `weno6`, `quintic-mp` |
| `[numerics] reflux_policy` | `always`, `troubled` |
| `[diagnostics] divb_flux` | `.false.`, `.true.` |

One change has no key: under EGLM the reflux leaves $\mathbf{B}$ alone.

Every case runs on 2 ranks from its `check.sh` under `src/tests/flume/verification/amr-*/`; `FLUME_EXE` selects the
backend. Every case has a figure.

| Case | Type | What it verifies |
|---|---|---|
| [Order across static 2:1 seams](./order-across-seams) | P | the order of a smooth vortex crossing a refined quadrant, per region, against the uniform control; the 6-point restriction; the quintic and limited coarse-to-fine fills |
| [Divergence across seams](./divergence-across-seams) | P | D_flux, the divergence the cleaning acts on: bounded and first-order convergent at the seams; the excess $B_z$ of EGLM on a 2:1 grid ([#78](https://github.com/szaghi/adam/issues/78)) is not a divergence difference |
| [Shock crossing a static 2:1 face](./shock-crossing) | P | mass, shock position and the plateau behind a shock crossing a 2:1 face, against the exact solution; the 6-point and WENO restrictions at shocks |
| [Integral drift and the WENO sensor near seams](./integral-drift) | P+C | the drift of the integrals without the reflux, and what the WENO sensor reads at the seams |
| [Reflux policy: troubled against always](./reflux-policy) | P | the sensor-gated reflux recovers the sixth order with shocks unchanged, at a conservation cost |
| [EGLM reflux leaves B alone](./eglm-reflux) | P | the excess $B_z$ of EGLM on a 2:1 grid ([#78](https://github.com/szaghi/adam/issues/78)) is caused by the reflux of $\mathbf{B}$; without it EGLM's $B_z$ drops from 70 to 3.3 times the uniform level |
| [Seam design record](./design-record) | P | appendix: the 2-D prototypes, the restriction footprint, the capturing error behind a crossing shock, the GLM damping length |

Two more seam checks sit on other pages: the 6-point restriction and fill across an inter-realm 2:1 seam on the
[forest seams](../euler/multi-realm) (leg 6) and on the [ghost probe](../infrastructure/ghost-probe), and the `divb_tol`
monitor on [GLM cleaning](../mhd/glm-cleaning).

## How it was built

The seam operators are milestone M6 ([#68](https://github.com/szaghi/adam/issues/68)), which measures, and then lifts,
the accuracy cap of the static 2:1 seams. Phase A adds the instruments (AO, AD, AS, AL); phases C, D, E and G change the
operators. The keys came with phases C1 and C3 (`seam_restriction`), D1 and D2 (`seam_ghost_fill`), E1 and E2
(`reflux_policy`) and A2 (`divb_flux`); the EGLM reflux that leaves $\mathbf{B}$ alone is G3b. Leg 6 of the forest seams
and the ghost probe cover C2/D1; the `divb_tol` monitor is MV-10 (G1).

| Case | ID | Phase |
|---|---|---|
| [Order across static 2:1 seams](./order-across-seams) | AO | A1, C1, D1, D2 |
| [Divergence across seams](./divergence-across-seams) | AD | A2 |
| [Shock crossing a static 2:1 face](./shock-crossing) | AS | A3, C3, C3c |
| [Integral drift and the WENO sensor near seams](./integral-drift) | AL | A4 |
| [Reflux policy: troubled against always](./reflux-policy) | E2 | E1, E2 |
| [EGLM reflux leaves B alone](./eglm-reflux) | — | G3a, G3b |
| [Seam design record](./design-record) | — | B1–B3, E0, G2 |

The [verification catalog](../) places the group in implementation order.
