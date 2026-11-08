# Inter-realm seams between misaligned blocks

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `rmf-2realm-fv-pulse-blocks` | [#51](https://github.com/szaghi/adam/issues/51), [#54](https://github.com/szaghi/adam/issues/54) | C | CPU, FNL | no | — |

## Case and scope

The 1:1 mirror pulse of [inter-realm 1:1 seams](./inter-realm-mirror), on the finite-volume path, with the seam blocks
of the two realms not lined up: one realm carries the same cells in smaller blocks, or in a single block. Each block
face of the fine side (realm 2, `realm_b`) is scattered into every register face it overlaps. On a mirror seam the
fluxes of the two sides are equal, so every register face must close: the printed mismatch
$\max\lvert F_\text{coarse} - \sum F_\text{fine}\rvert$ must be exactly 0.

What the case catches:
- a fine face scattered into the wrong register face, or into too few
  ([#51](https://github.com/szaghi/adam/issues/51));
- a seam-peer lookup that fails when one realm is a single block: the lookup then walks up from a level-1 code to the
  root, and aborted with error −111 before [#54](https://github.com/szaghi/adam/issues/54).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, no source | [models](../models#models) |
| Space | centred finite volumes of order 6 (`fv_centered`) | [numerics](../numerics#centred-finite-volumes-fv-centered) |
| Time | SSP Runge–Kutta (5,4), CFL 0.5, β cadence | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Register | fine skins scattered over the overlapping register faces (`maps%seam_overlap`) | `accumulate_fine_overlaps`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90); [`adam_maps_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_maps_object.F90) |
| Peer lookup | closest block of the peer tree (`get_closest_block`) | [`adam_tree_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_tree_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Realms | $x \in [-0.06, 0]$ and $[0, 0.06]$, $y, z \in [-0.06, 0.06]$ m, `mirror` seam at $x = 0$, β cadence |
| Default realm grid | 8 blocks of $8\times16\times16$, $\Delta x = 3.75\cdot10^{-3}$ |
| `small1`, `small2` | realm 1 (resp. 2) on 64 blocks of $4\times8\times8$ (`iu_ref_levels = 2`): the same cells |
| `one1`, `one2` | realm 1 (resp. 2) on one block of $16\times32\times32$ (`iu_ref_levels = 0`, `max_level` kept 1) |
| Reference | the single-realm FV pulse: `rmf-amr-fd-pulse` on `fv_centered`, marker off, 8 blocks of $16^3$ |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+y$, `sigma = 0.01`, `B0 = 10` |
| Boundaries | `Neumann` on the outer faces |
| Steps, ranks | 20, 2 |

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/rmf-2realm-fv-pulse-blocks/input.ini`, `realm_{1,2}.ini` | the manifest and the realms |
| `src/tests/prism/regression/rmf-2realm-fv-pulse-blocks/check.sh` | builds the reference and the four legs (`--build`, `--np N`; `PRISM_EXE`) |
| `src/tests/flume/verification/multirealm/multirealm_oracle.py` | the union of the realms against the single realm, cell by cell |

```bash
cd src/tests/prism/regression/rmf-2realm-fv-pulse-blocks && ./check.sh
PRISM_EXE=$PWD/../../../../../exe/adam_prism_fnl ./check.sh        # FNL
```

## Results

Each leg asserts: the run reaches 100 % with no error/abort/NaN; the reflux line is printed and its maximum is exactly
0; the union of the realms matches the single-realm FV pulse within $10^{-13}$ (relative).

| Leg | Register lines (face × step) | max mismatch | Fields against the single realm | Backends |
|---|---|---|---|---|
| `small1` (register side on small blocks: each register face holds $2\times2$ realm-2 faces) | 320 | 0 | $\le 3.7\cdot10^{-15}$ | CPU, FNL |
| `small2` (each realm-2 face scattered into 4 register faces) | 80 | 0 | $\le 3.7\cdot10^{-15}$ | CPU, FNL |
| `one1` | ⟨run⟩ | ⟨run⟩ | ⟨run⟩ | ⟨run⟩ |
| `one2` | ⟨run⟩ | ⟨run⟩ | ⟨run⟩ | ⟨run⟩ |

Related pages: [inter-realm 1:1 seams](./inter-realm-mirror) on aligned blocks; [static AMR](./amr-reflux) for the
register on intra-realm 2:1 faces.

## Findings

- The register closes exactly on every face and step, whichever side carries the small blocks.
- The fields are compared within $10^{-13}$, not bitwise: PRISM FV is not invariant under a change of block layout at
  round-off. The single realm alone on 4 blocks per axis instead of 2 differs by $6.6\cdot10^{-16}$. The register
  check is the exact one.
- A realm made of one block runs: the peer lookup reaches the root without aborting.
