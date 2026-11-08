# Input validation and refusals

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `refusals.sh`, `rmf-hyperbolic` | [#11](https://github.com/szaghi/adam/issues/11), [#55](https://github.com/szaghi/adam/issues/55), [#56](https://github.com/szaghi/adam/issues/56), [#57](https://github.com/szaghi/adam/issues/57), [#61](https://github.com/szaghi/adam/issues/61), [#84](https://github.com/szaghi/adam/issues/84) (E3b, E3d), [#111](https://github.com/szaghi/adam/issues/111) | C | CPU, FNL | no | `rmf-hyperbolic` (digest + residuals, CPU and FNL) |

## Case and scope

An input PRISM cannot run correctly must stop at initialization with a message that names the option and the value.
Before the issues below, several such inputs ran: from the zero field, with an undefined index, without a boundary
condition, or with a crash at the first step.

`refusals.sh` rewrites a base input (`rmf`, or `pic-pec-t0` for a PIC leg) with one `sed` script per leg and requires
two things: a non-zero exit **and** the refusal's own message in the log. A crash also exits non-zero; only the message
tells a refusal from a crash. A leg whose script changes nothing fails, so a drift of the base input cannot make a leg
pass vacuously.

`rmf-hyperbolic` is the positive counterpart of two refusals: the one case with hyperbolic (GLM) cleaning, on
$\mathbf D$ and $\mathbf B$. It anchors the state layout that [#11](https://github.com/szaghi/adam/issues/11) and
[#61](https://github.com/szaghi/adam/issues/61) fixed: `nv = 11`, $\phi$ and $\psi$ in slots 7 and 8, $\mathbf J$ in
9–11. The full divergence-correction matrix is on the [hazards page](../hazards#divergence-correction).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Parsing | every value checked at initialization; values pass through `strip_control` (CRLF-safe) | [`adam_prism_numerics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_numerics_object.F90), [`adam_prism_physics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_physics_object.F90), [`adam_prism_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_bc_object.F90), [`adam_prism_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_ic_object.F90) |
| `rmf-hyperbolic` | `rmf` with `constrained_transport = DB`, `divergence_correction = hyperbolic`, `chi = 1`, damping `alpha = 1e9` 1/s ($\alpha t \approx 0.1$ over 10 steps) | [models](../models#divergence-correction) |
| Refused combinations | the list and the reasons | [hazards](../hazards#refused-combinations) |

## Initial and boundary conditions

Every leg runs at 1 rank and stops at initialization, in a few seconds. `rmf-hyperbolic` is the `rmf` grid: 64 blocks
of $16^3$ in $[-0.06, 0.06]^3$ m, 4 AC coils, `Neumann` faces, 10 steps, 2 ranks.

## Input files

| File | Role |
|---|---|
| `src/tests/prism/regression/refusals.sh` | the legs (`PRISM_EXE` selects the executable; the backend is its suffix, and a leg naming other backends is skipped) |
| `src/tests/prism/regression/rmf/input.ini`, `pic-pec-t0/input.ini` | the base inputs |
| `src/tests/prism/regression/rmf-hyperbolic/input.ini` | the GLM case |
| `src/tests/lint/check-removed-input-keys.sh` | rejects removed keys (`[fWLayer] C`, `[PIC] elliptic_correction`, `[linear-algebra]`) in every committed input |

```bash
cd src/tests/prism/regression && ./refusals.sh                       # CPU
PRISM_EXE=$PWD/../../../../exe/adam_prism_fnl ./refusals.sh          # FNL
```

## Results

Every leg must exit non-zero with the message below in its log.

| Leg | Backends | Message (`grep -F`) | Refused because | Issue |
|---|---|---|---|---|
| `hyperbolic-no-ct` | cpu, fnl | `hyperbolic needs` | no $\phi$/$\psi$ slot without constrained transport; `var_J*` were left undefined | #61 |
| `poisson` | cpu, fnl | `poisson is not implemented` | no backend implements the Poisson projection | #61 |
| `blanes-moan` | cpu, fnl | `Blanes-Moan time integration` | the integrator object is never initialized; the FNL routines are empty | #55 |
| `cfm` | cpu, fnl | `Commutator-Free Magnus time integration` | as above | #55 |
| `leapfrog` | cpu, fnl | `leapfrog time integration` | as above | #55 |
| `yoshida` | fnl | `Yoshida Runge-Kutta time integration` | `integrate_rk_yoshida_dev` is empty; the CPU integrator works | #55 |
| `ic-riemann-problem` | cpu, fnl | `riemann-problem is not implemented` | the initial condition was never applied: the run started from the zero field | #56 |
| `ic-magnetic-nozzle` | cpu, fnl | `is not implemented as an initial condition` | as above | #56 |
| `ic-unknown` | cpu, fnl | `unknown [initial_conditions].(type)` | an unknown type | #56 |
| `ic-rmf-axis` | cpu, fnl | `unknown [external_fields].(RMF_rotation_axis)` | the component indices were left undefined | #56 |
| `bc-radiative` | cpu, fnl | `radiative is not implemented` | the CPU stopped at the first BC application; the FNL kernel applied none | #57 |
| `scheme-time-typo` | cpu, fnl | `(scheme_time) "runge-kutta"` | a typo fell back silently | #57 |
| `scheme-space-typo` | cpu, fnl | `(scheme_space) "fd_centred"` | as above | #57 |
| `reconstruction-typo` | cpu, fnl | `(reconstruction_variables) "characteristichs"` | as above | #57 |
| `ct-typo` | cpu, fnl | `(constrained_transport) "DD"` | as above | #57 |
| `dc-typo` | cpu, fnl | `(divergence_correction) "hyperbolc"` | as above | #57 |
| `ef-typo` | cpu, fnl | `(external_fields_applied) "rmff"` | as above | #57 |
| `fwl-flag-typo` | cpu, fnl | `(x_minus_layer) "yse"` | as above | #57 |
| `ib-file-off` | cpu, fnl | `(definition) "file.off"` | as above | #57 |
| `ib-axis` | cpu, fnl | `(circle_axis) "w"` | as above | #57 |
| `t0-elliptic-amr` | cpu, fnl | `needs a uniform grid` | a t0 elliptic solve on a refined grid: the wide operator across 2:1 faces is not defined | #84 E3b |
| `elliptic-correction` | cpu, fnl | `elliptic_correction) is removed` | `[PIC] elliptic_correction = .true.`, removed; base input `pic-pec-t0` | #84 E3d |
| `probes-sections` | cpu, fnl | `[probe_N] sections are no longer read` | the former `[probes]` input is refused, not ignored | #111 |
| `probes-number` | cpu, fnl | `(probes_number) is no longer read` | as above | #111 |

| Check | Bound | CPU | FNL |
|---|---|---|---|
| `refusals.sh` | every leg a clean refusal | 23 legs | 24 legs |
| `rmf-hyperbolic` digest and residuals | `rtol = 1e-6`, `atol = 1e-3` | pass | pass |

Related pages: the [input reference](../input) lists every key; [slices](./slices) carries the input validation of
the `[slice_N]` sections; [seam div(B)](./seam-divergence) holds the RK-contract refusal of the staged forest path.

## Findings

- No input in the list runs any more: each stops at initialization and names the option and the value.
- The message requirement separates a refusal from a crash, which also exits non-zero.
- Add a leg whenever an input becomes a refusal.
