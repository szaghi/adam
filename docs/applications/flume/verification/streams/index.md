# Supersonic boundary layers: boundaries and shock interaction

The compressible inflow and boundary toolchain that a shock/boundary-layer interaction needs, after the STREAmS solver
(Bernardini et al., CPC 2021, 2023): Sutherland's viscosity law, a characteristic non-reflecting outflow and far field,
an oblique-shock generator on a y face, and the wall-modelled interaction itself. The boundaries are verified against
the kinds they replace and against the oblique-shock relations and the exact regular reflection; the interaction
against the STREAmS DNS and the IUSTI experiment. It continues the
[compressible boundary layer](../turbulence/boundary-layer) of the turbulence group.

| Case | Type | What it verifies |
|---|---|---|
| [Non-reflecting boundary](./non-reflecting) | P | `[bc_*] type = non-reflecting` against the kinds it would replace, `extrapolation` and `inflow`: what an outgoing packet leaves in the box, the relaxation of a uniform overpressure, bit-for-bit equality with `extrapolation` where it must reduce to it |
| [Shock generator and regular reflection](./shock-generator) | P | the shock of `type = shock-generator` against the oblique-shock relations, and its regular reflection off a slip wall against the exact reflection |
| [Shock/boundary-layer interaction](./sbli) | P | the wall-modelled interaction at Mach 2.28: the interaction length, the separation from $C_f$ and the wall pressure against the STREAmS DNS and the IUSTI experiment; a refined wall grid under the interaction (in progress) |
| [Sutherland viscosity](#sutherland-viscosity) | P+C | the law against the formula, both wall models under it, its input contract, the cold-wall boundary layer rerun under it |

The non-reflecting and shock-generator legs are in `src/tests/flume/verification/boundary/check.sh`; their `sg` legs are
distinct from the `sg` legs of the [subgrid-scale models](../turbulence/sgs-models).

## Sutherland viscosity

The law $\mu = \mu_\mathrm{ref}(T/T_\mathrm{ref})^{3/2}(T_\mathrm{ref} + S)/(T + S)$ is in [models](../../models#dissipative-terms)
and its keys in [input](../../input#physics). It has no page of its own; it is verified where it is used:

| Where | What it verifies |
|---|---|
| [unit tests](../unit-tests), `test_wall_model` (leg wm2) | the law against the formula ($6.7\cdot10^{-16}$); both wall models under Sutherland against their independent fine-grid references |
| [input contract of the dissipative terms](../dissipative/input-contract) (DC) | the Sutherland and kelvin refusals, $S$ logged exactly from the kelvin anchor, two conversions through the `[reference]` layer |
| [boundary layer, Mach 5.86](../turbulence/boundary-layer#mach-5-86-cold-wall) | the cold-wall layer rerun under Sutherland's law with the DNS's Reynolds number: $C_f$ −3.4 % against −3.8 %, $B_q$ −0.0190 |
| [shock/boundary-layer interaction](./sbli) | the interaction run uses Sutherland's law, $S/T_\infty = 0.776$ |

The power and constant laws are unchanged bit for bit: with $S = 0$ the law takes the operations of the power law.

## How it was built

The toolchain is [#121](https://github.com/szaghi/adam/issues/121), which continues the compressible boundary layer of
[#119](https://github.com/szaghi/adam/issues/119) T3.7.

| Phase | What | Legs | Page |
|---|---|---|---|
| S1 | Sutherland viscosity (`[physics] viscosity_law = sutherland`), the kelvin anchor `temperature_unit_kelvin` | `test_wall_model` (wm2); dissipation contract DC | [Sutherland viscosity](#sutherland-viscosity) |
| S2 | Characteristic non-reflecting outflow and far field (`[bc_*] type = non-reflecting`) | nr1–nr4 | [Non-reflecting boundary](./non-reflecting) |
| S3 | Oblique-shock generator on a y face (`type = shock-generator`) | sg1–sg3 (boundary) | [Shock generator](./shock-generator) |
| S4a | Inviscid regular reflection off a slip wall | sg4 (boundary) | [Shock generator](./shock-generator#regular-reflection-off-a-slip-wall) |
| S4b | Wall-modelled interaction, Mach 2.28, IUSTI | campaign, no leg | [Shock/boundary-layer interaction](./sbli) |
| S4c | S4b with a 2:1 box under the interaction | campaign, in progress | [Shock/boundary-layer interaction](./sbli#a-refined-wall-grid-under-the-interaction) |
| S5 | Recycling–rescaling turbulent inflow | — | not started |

| Page | Phases | Legs | Type |
|---|---|---|---|
| [Non-reflecting boundary](./non-reflecting) | S2 | nr1–nr4 | P |
| [Shock generator](./shock-generator) | S3, S4a | sg1–sg4 (boundary) | P |
| [Shock/boundary-layer interaction](./sbli) | S4b, S4c | campaign | P |
