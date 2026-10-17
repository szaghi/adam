# Seam studies (M6)

These pages are the investigations behind milestone M6 ([issue #68](https://github.com/szaghi/adam/issues/68)): the
prototypes, operator variants and ablations that decided how the static 2:1 seams were changed. The verification cases
they build on, and that freeze the outcome, are in [seam accuracy and divergence](../seams/).

The phases of #68, in the order the evidence led:
- **B (prototypes, NumPy, no ADAM run):** B1 confirms in 2-D that the mean restriction and the reflux each cap the order
  at 2; B3 picks the shift-inward footprint for the 6-point restriction; B2 gives the go for a sensor-gated reflux
  (D-M6-4).
- **C (restriction):** `[amr] seam_restriction = lagrange6`, `weno6`, measured on [AO](../seams/order-across-seams)
  and [AS](../seams/shock-crossing).
- **D (coarse-to-fine fill):** D1 adds the quintic fill, which reaches the scheme's order on smooth flows; D2 adds two
  limited fills and shows that the shock-crossing defect is not the fill's.
- **E (reflux):** E0 shows the shock-crossing defect is a capturing error, no seam operator's; E2 measures
  `[numerics] reflux_policy = troubled` (D-M6-2: the default stays `always`).
- **G (divergence):** G2 keeps the GLM damping length frozen (D-M5-9); G3a finds that the reflux causes #78; G3b fixes it
  by not refluxing $\mathbf{B}$ under EGLM.

The defaults are unchanged by every study except G3b: `mean` restriction, `tricubic` fill, `reflux_policy = always`.
No study has a figure yet: they are listed in [#85](https://github.com/szaghi/adam/issues/85). The prototypes run from
`src/tests/flume/verification/amr-proto/check.sh` (about 5 min; `--legs b1,b2,b3` selects the prototypes, `--long` the
finer ladders).

| Case | ID | What it verifies |
|---|---|---|
| [Quintic coarse-to-fine fill](./quintic-fill) | D1 | the 6-point fill brings smooth flows to sixth order; unlimited, it worsens a shock entering a finer region |
| [Limited coarse-to-fine fills](./limited-fills) | D2 | `weno6` and `quintic-mp` keep the order; no fill controls the shock-crossing defect |
| [Shock-crossing defect is a capturing error](./capturing-error) | E0 | neither the restriction, the fill/restriction consistency nor the reflux causes the defect |
| [Reflux policy: troubled against always](./reflux-policy) | E2 (E1) | the sensor-gated reflux recovers the high order with shocks unchanged, at a conservation cost |
| [#78 is the reflux, not the restricted ghosts](./reflux-not-ghosts) | G3a | EGLM's $B_z$ excess on a 2:1 grid comes from the reflux |
| [GLM damping length on a seam](./glm-damping) | G2 | every damping length keeps the seam divergence bounded; strong damping doubles its floor |
| [EGLM reflux leaves B alone](./eglm-reflux) | G3b | the fix of #78: no reflux of $\mathbf{B}$ under EGLM (70 → 3.3 times uniform) |
| [2:1 caps in a 2-D prototype](./caps-prototype) | B1 | the mean restriction and the reflux each cap the order at 2 in 2-D; the corners are not the cause |
| [Sensor-gated reflux prototype](./sensor-gated-reflux) | B2 | a `troubled` reflux recovers order 5 and keeps the shocks: go for Phase E |
| [Footprint of a high-order restriction](./high-order-restriction) | B3 | shift-inward and two-pass are equivalent with the quintic prolongation; shift-inward chosen |
