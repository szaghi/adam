# Physical models

PRISM solves Maxwell's equations in first-order form for the electric displacement $\mathbf D$ and the magnetic
induction $\mathbf B$, in vacuum ($\mathbf D = \varepsilon_0\mathbf E$, $\mathbf B = \mu_0\mathbf H$), driven by a current
density $\mathbf J$:

$$\frac{\partial \mathbf D}{\partial t} = \nabla\times\frac{\mathbf B}{\mu_0} - \mathbf J, \qquad
\frac{\partial \mathbf B}{\partial t} = -\nabla\times\frac{\mathbf D}{\varepsilon_0}, \qquad
\nabla\cdot\mathbf D = \rho, \qquad \nabla\cdot\mathbf B = 0 .$$

In conservation form $\partial_t\mathbf q + \sum_d \partial_{x_d}\mathbf F_d(\mathbf q) = \mathbf S$, with the fluxes of
$(\mathbf D, \mathbf B)$

$$F_x = \Big(0,\ \tfrac{B_z}{\mu_0},\ -\tfrac{B_y}{\mu_0},\ 0,\ -\tfrac{D_z}{\varepsilon_0},\ \tfrac{D_y}{\varepsilon_0}\Big),$$

and $F_y$, $F_z$ by cyclic permutation; the source is $-\mathbf J$ on $\mathbf D$. The characteristic speeds are
$\pm c_0$, $c_0 = 1/\sqrt{\varepsilon_0\mu_0}$, and two zero speeds. $\mathbf J$ has no flux: it is set at every stage
by the coils and, under PIC, deposited by the particles.

## Models

`[physics] physical_model` selects one of four models ([input](./input#physics); an unknown value is fatal):

| Model | `physical_model` | Current and charge | Units |
|---|---|---|---|
| Electromagnetic | `EM` (and aliases) | $\mathbf J$ from the coils | SI |
| Non-dimensional electromagnetic | `AdimEM` | as EM | $L_0$, $B_0$ (below) |
| Particle-in-cell | `PIC` | $\mathbf J$ from the coils and the particles; $\rho$ from the particles, carried in the state | SI |
| Non-dimensional PIC | `AdimPIC` | as PIC | $L_0$, $B_0$ |

The particles, their deposition and their push are on the [particle-in-cell page](./pic).

## State vector

The state width `nv` follows from the model and the divergence correction, with one source for both the width and the
variable names (`physics%var_phi`, `physics%var_psi`, `var_Jx`; [#61](https://github.com/szaghi/adam/issues/61)):

| Slot | Variable | Present |
|---|---|---|
| 1–3 | $D_x, D_y, D_z$ | always |
| 4–6 | $B_x, B_y, B_z$ | always |
| 7 | $\varphi$ (GLM scalar of $\mathbf D$) or $\psi$ (of $\mathbf B$) | `hyperbolic` with `D` or `B` |
| 7, 8 | $\varphi$, $\psi$ | `hyperbolic` with `DB` |
| next 3 | $J_x, J_y, J_z$ | always |
| last | $\rho$ | PIC models |

So `nv = 9` (EM), 10 (one GLM scalar, or PIC), 11 (`DB`, or PIC with one scalar) and 12 (PIC with `DB`), the maximum
`NV_MAX = 12`; a larger `nv` is fatal at initialization.

## Physical constants

`adam_prism_parameters`:

| Constant | Value |
|---|---|
| $\mu_0$ | `1.256637e-6` H/m |
| $\varepsilon_0$ | `8.854187e-12` F/m |
| $c_0$ | $1/\sqrt{\mu_0\varepsilon_0}$ |
| electron charge | `-1.6e-19` C |
| electron mass | `9.11e-31` kg |
| ion mass (H⁺) | `1.67262192595e-27` kg |
| neutral mass | ion mass + electron mass |
| $k_B$ | `1.380649e-23` J/K |

## Non-dimensional form

Under `AdimEM` and `AdimPIC` the input stays in SI and PRISM scales it by the reference length $L_0$ and field $B_0$
(`[physics] L0`, `B0`, both required and positive):

$$T_0 = \frac{L_0}{c_0}, \qquad D_0 = B_0 c_0\varepsilon_0, \qquad I_0 = \frac{B_0 L_0}{\mu_0}, \qquad
\rho_0 = \frac{D_0}{L_0}, \qquad Q_0 = \rho_0 L_0^3 .$$

The equations become $\partial_t\mathbf D = \nabla\times\mathbf B - \mathbf J$, $\partial_t\mathbf B = -\nabla\times\mathbf D$,
with unit speed of light. Lengths (grid, coil geometry, pulse widths), times, fields and currents are converted at
initialization.

## Divergence correction

The discrete operators keep $\nabla\cdot\mathbf B$ at round-off on a uniform grid ([numerics](./numerics#the-discrete-curl-and-divergence)),
but not across a 2:1 AMR face, a deposit that does not conserve charge, or a source that is not divergence-free.
`[numerics] divergence_correction = hyperbolic` adds the generalized-Lagrange-multiplier (GLM) scalars selected by
`constrained_transport` (`D`: $\varphi$, `B`: $\psi$, `DB`: both):

$$\frac{\partial \mathbf D}{\partial t} = \nabla\times\frac{\mathbf B}{\mu_0} - \nabla\varphi - \mathbf J, \qquad
\frac{\partial \varphi}{\partial t} = -c_h^2\,(\nabla\cdot\mathbf D - \rho) - \alpha\,\varphi,$$

$$\frac{\partial \mathbf B}{\partial t} = -\nabla\times\frac{\mathbf D}{\varepsilon_0} - \nabla\psi, \qquad
\frac{\partial \psi}{\partial t} = -c_h^2\,\nabla\cdot\mathbf B - \alpha\,\psi,$$

with $\rho$ present only under PIC. The cleaning speed is $c_h = \chi c_0$ (`[physics] chi`, required, positive) and
$\alpha$ the damping rate (`[physics] alpha`, default 0; in 1/s under EM and PIC, in units of $1/T_0$ under the non-dimensional models, where it is not converted; with $\alpha > 0$ the time step is also capped at
$1/\alpha$). The fastest signal is then $\max(1, \chi)\,c_0$, which sets the CFL time step.

- `hyperbolic` without `D`, `B` or `DB` is refused at initialization: there would be no scalar to carry.
- `poisson` (a projection) is refused: no backend implements it.
- `no` ignores `constrained_transport` with a warning.

The full matrix of `constrained_transport` × `divergence_correction` is in the
[input hazards](./input#divergence-correction).
