# Shock generator and regular reflection

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| SG | [#121](https://github.com/szaghi/adam/issues/121) · S3, S4a | P | sg1–sg4 (boundary) | CPU, FNL | yes |

## Case and scope

`[bc_y_max] type = shock-generator`
([boundary conditions](../../boundary-conditions#shock-generator-an-oblique-shock-through-a-y-face)) brings an oblique
shock in through the top face. The legs measure the shock it makes against the oblique-shock relations (S3) and its
regular reflection off a slip wall (S4a), with `shockgen_oracle.py`, which solves the relations on its own. These `sg`
legs belong to `src/tests/flume/verification/boundary/check.sh`, not to the
[subgrid-scale models](../turbulence/sgs-models).

## Setup

| Setup | |
|---|---|
| box, cells | $[0, 2]\times[0, 1]\times[0, 0.25]$, $64\times32\times8$, 8 blocks, periodic z ($96\times48\times8$ for sg4) |
| gas | Euler, $\gamma = 1.4$, $\rho = 1$, $p = 1/\gamma$ ($c = 1$), Mach 2 along x |
| boundary | `inflow` at $x = 0$; `extrapolation` at $x = 2$ and $y = 0$ (`wall-inviscid` for sg4); `shock-generator` at $y = 1$, deflection $10^\circ$ from `shock_x` = 0.5 ($\beta = 39.314^\circ$, $p_2/p_1 = 1.70658$, $\rho_2/\rho_1 = 1.45843$) |
| scheme | WENO-U-5 characteristic, SSP-54, CFL 0.5, 2 ranks, to $t = 2.5$ (steady; $t = 5$ for sg4) |
| script | `src/tests/flume/verification/boundary/check.sh` |

## Results

| Leg | Case | Oracle | Measured (CPU; FNL the same to the digits shown) |
|---|---|---|---|
| sg1 | the shock by its deflection | the uniform region behind it (16 cells off the shock, 4 below the streamline from its foot, above $y = 0.2$): $p$ and $\rho$ within 3e-3 of Rankine–Hugoniot, the flow angle within $0.05^\circ$; the shock angle fitted from where $p$ crosses $(p_1 + p_2)/2$ along the rows within $0.3^\circ$ | $p$ 1.47e-3, $\rho$ 1.22e-3, angle $0.026^\circ$; $\beta$ $39.322^\circ$, the foot at $x = 0.496$ |
| sg2 | the same shock by its angle (`shock_angle` 39.313932) | as sg1 | the same digits |
| sg3 | a shock entering beyond the box (`shock_x` 3) | equal to an extrapolated top bit for bit | 0 of 64 datasets differ |
| sg4 | the regular reflection off a slip wall (S4a): `wall-inviscid` at $y = 0$, `shock_x` = 0.1 (the shock meets the wall at $x = 1.321$), $t = 5$ | behind the reflected shock (8 cells off it and off the wall) $p$ and $\rho$ within 4e-3 of the exact $p_3/p_1 = 2.80319$, the flow parallel to the wall within $0.1^\circ$; the reflected shock's angle to the wall within $0.3^\circ$ of $39.384^\circ$ | $p$ 2.10e-3, $\rho$ 1.50e-3, $0.062^\circ$; $39.437^\circ$, its foot at $x = 1.324$ |

## The incident shock

![The shock generator at Mach 2](/flume/shockgen.png)

*Left: $p/p_1$ in the x-y plane of leg sg1; dashed, the shock of the theory from its foot; dotted, the streamline from the foot. Right:
$p/p_1$ and $\rho/\rho_1$ along the row $y = 0.484$ against the step of Rankine–Hugoniot at the theory's shock.*

**What it shows.**
- **The shock is where the relations put it.** Its angle is within $0.01^\circ$ of $\beta$, and it meets the top within a
  cell of `shock_x`; behind it the state is Rankine–Hugoniot's within $1.5\cdot10^{-3}$.
- **What the oracle leaves out, and why.** The captured shock leaves oscillations of a few $10^{-3}$ that decay over 16
  cells behind it (the staircase of an oblique shock on a Cartesian grid). The jump of the ghost state at `shock_x` leaves
  an entropy streak along the streamline from the foot ($\rho$ 7e-3 off, the pressure continuous). Where the shock leaves
  through the extrapolated bottom, the corner downstream sits up to 4 % high: an extrapolated face crossed by an oblique
  shock is not clean. In a shock/boundary-layer interaction the shock meets the wall instead.

## Regular reflection off a slip wall

![The regular reflection off a slip wall](/flume/shockgen-reflect.png)

*Leg sg4. Left: $p/p_1$ in the x-y plane; dashed, the incident shock of the theory; dotted, the reflected one ($M_2 = 1.641$, turned back
by $10^\circ$). Right: $p/p_1$ along the row $y = 0.260$ against the steps of the exact regular reflection.*

**What it shows.**
- **The reflection is the exact one.** Behind the reflected shock the state is within $2\cdot10^{-3}$ of
  $p_3/p_1 = 2.803$, the flow parallel to the wall within $0.06^\circ$, and the reflected shock lies within $0.06^\circ$
  of its angle, its foot within a cell of where the incident shock meets the wall.
- **It needs the steady state.** At $t = 2.5$ the region behind the reflected shock is still 6 % low near the outflow:
  the reflected shock forms at the wall and straightens outward, so the leg runs to $t = 5$.

The viscous interaction with a turbulent boundary layer is [S4b](./sbli).
