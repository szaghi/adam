# Seam accuracy and divergence (M6)

Milestone M6 ([issue #68](https://github.com/szaghi/adam/issues/68)) measures, and then lifts, the accuracy cap of the
static 2:1 seams. Phase A adds the instruments (AO, AD, AS, AL); the transfer operators and the reflux policy change
only behind new keys, whose defaults keep today's behaviour:

| Key | Values (default first) | Phase |
|---|---|---|
| `[amr] seam_restriction` | `mean`, `lagrange6`, `weno6` | C1, C3 |
| `[amr] seam_ghost_fill` | `tricubic`, `injection`, `restriction-compatible`, `quintic`, `weno6`, `quintic-mp` | D1, D2 |
| `[numerics] reflux_policy` | `always`, `troubled` | E1, E2 |
| `[diagnostics] divb_flux` | `.false.`, `.true.` | A2 |

The investigations behind the phases, and the decisions they led to, are the [seam studies](../seams-studies/).

Every case runs on 2 ranks from its `check.sh` under `src/tests/flume/verification/amr-*/`; `FLUME_EXE` selects the
backend. No M6 case has a figure yet: they are listed in [#85](https://github.com/szaghi/adam/issues/85).

| Case | ID | What it verifies |
|---|---|---|
| [Order across static 2:1 seams](./order-across-seams) | AO (A1, C1) | the order of a smooth vortex crossing a refined quadrant, per region, against the uniform control; the 6-point restriction |
| [Divergence across seams](./divergence-across-seams) | AD (A2) | D_flux, the divergence the cleaning acts on: bounded and first-order convergent at the seams; #78 is not a divergence difference |
| [Shock crossing a static 2:1 face](./shock-crossing) | AS (A3, C3, C3b, C3c) | mass, shock position and the plateau behind a shock crossing a 2:1 face, against the exact solution; the 6-point and WENO restrictions at shocks |
| [Integral drift and the WENO sensor near seams](./integral-drift) | AL (A4) | the drift of the integrals without the reflux, and what the WENO sensor reads at the seams |
