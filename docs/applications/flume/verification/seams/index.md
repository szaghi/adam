# Seam accuracy and divergence (M6)

Milestone M6 ([#68](https://github.com/szaghi/adam/issues/68)) measures, and then lifts, the accuracy cap of the static
2:1 seams. Phase A adds the instruments (AO, AD, AS, AL); phases C, D, E and G change the operators. The transfer
operators and the reflux policy change only behind new keys; the defaults are `mean`, `tricubic` and `always`:

| Key | Values (default first) | Phase |
|---|---|---|
| `[amr] seam_restriction` | `mean`, `lagrange6`, `weno6` | C1, C3 |
| `[amr] seam_ghost_fill` | `tricubic`, `injection`, `restriction-compatible`, `quintic`, `weno6`, `quintic-mp` | D1, D2 |
| `[numerics] reflux_policy` | `always`, `troubled` | E1, E2 |
| `[diagnostics] divb_flux` | `.false.`, `.true.` | A2 |

One change has no key: under EGLM the reflux leaves $\mathbf{B}$ alone (G3b).

Every case runs on 2 ranks from its `check.sh` under `src/tests/flume/verification/amr-*/`; `FLUME_EXE` selects the
backend. Every case has a figure.

| Case | ID | Phase | Type | What it verifies |
|---|---|---|---|---|
| [Order across static 2:1 seams](./order-across-seams) | AO | A1, C1, D1, D2 | P | the order of a smooth vortex crossing a refined quadrant, per region, against the uniform control; the 6-point restriction; the quintic and limited coarse-to-fine fills |
| [Divergence across seams](./divergence-across-seams) | AD | A2 | P | D_flux, the divergence the cleaning acts on: bounded and first-order convergent at the seams; #78 is not a divergence difference |
| [Shock crossing a static 2:1 face](./shock-crossing) | AS | A3, C3, C3c | P | mass, shock position and the plateau behind a shock crossing a 2:1 face, against the exact solution; the 6-point and WENO restrictions at shocks |
| [Integral drift and the WENO sensor near seams](./integral-drift) | AL | A4 | P+C | the drift of the integrals without the reflux, and what the WENO sensor reads at the seams |
| [Reflux policy: troubled against always](./reflux-policy) | E2 | E1, E2 | P | the sensor-gated reflux recovers the sixth order with shocks unchanged, at a conservation cost |
| [EGLM reflux leaves B alone](./eglm-reflux) | — | G3a, G3b | P | #78 is caused by the reflux of $\mathbf{B}$; without it EGLM's $B_z$ drops from 70 to 3.3 times the uniform level |
| [Seam design record](./design-record) | — | B1–B3, E0, G2 | P | appendix: the 2-D prototypes, the restriction footprint, the capturing error behind a crossing shock, the GLM damping length |

Two M6 legs sit on earlier pages: leg 6 of the [forest seams](../euler/multi-realm) (the 6-point restriction and fill
across an inter-realm 2:1 seam, C2/D1 also on the [ghost probe](../infrastructure/ghost-probe)), and the `divb_tol`
monitor MV-10 (G1) on [GLM cleaning](../mhd/glm-cleaning).
