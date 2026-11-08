# Non-reflecting boundary

The open boundaries of [#121](https://github.com/szaghi/adam/issues/121): the non-reflecting outflow (S2, legs nr) and the
[shock generator](#shock-generator) (S3, legs sg).

`[bc_*] type = non-reflecting` ([#121](https://github.com/szaghi/adam/issues/121) S2) is a characteristic outflow and far-field
boundary. It extrapolates the outgoing waves, relaxes the incoming acoustic wave toward a target pressure (Rudy & Strikwerda 1980,
Poinsot & Lele 1992), and takes the target's entropy and tangential velocity where flow enters
([boundary conditions](../../boundary-conditions#non-reflecting-a-characteristic-outflow-and-far-field-121-s2)). The legs measure what
it does against the kinds it would replace, `extrapolation` and `inflow`. Script: `src/tests/flume/verification/boundary/check.sh`.

| Setup | |
|---|---|
| box, cells | $[0, 2]\times[0, 0.0625]^2$, $256\times8\times8$ ($64\times8\times8$ for nr3), 8 blocks, periodic y and z |
| gas | Euler, $\gamma = 1.4$, $\rho = 1$, $p_\infty = 1/\gamma$ ($c = 1$), free stream $u$ along x |
| packet | `sine-wave` `acoustic`, wavelength 0.2, amplitude $10^{-3}$, Gaussian envelope of width 0.1 at $x = 1$ |
| boundary | `non-reflecting` with the free stream as target, $L$ = 2 (`nr_length`), $\sigma$ = 0.28 unless stated |
| scheme | WENO-U-5 characteristic, SSP-54, CFL 0.5, 2 ranks |

| Leg | Case | Oracle | Measured (CPU; FNL the same to the digits shown) |
|---|---|---|---|
| nr1 | $\sigma = 0$ at a Mach 0.5 outflow; $\sigma = 0.28$ at a Mach 1.5 outflow (left face `inflow`) | equal to `extrapolation` bit for bit | 0 of 64 datasets differ, both |
| nr2 | the packet leaving through $x = 2$ at rest, $t$ = 1.5 | the largest $\lvert p - p_\infty\rvert$ left, over $\gamma p_\infty A$: below 4e-3 at $\sigma$ = 0.28, 1.1e-3 at 0.07, their ratio in [3, 5] | 2.86e-3, 7.56e-4, ratio 3.78; `extrapolation` 1.15e-4 |
| nr3 | a uniform 2 % overpressure at Mach 0.5, both x faces `non-reflecting`, $t$ = 20 | the box-mean overpressure left, over the initial: below 0.15; with `extrapolation`, above 0.95 | 0.099; 1.000 |
| nr4 | a left-running packet leaving through the subsonic inflow face $x = 0$ (Mach 0.5), $t$ = 3 | what is left: below 6e-3 | 4.28e-3; `inflow` 1.41e-5 |

![The non-reflecting boundary: reflection and pressure anchor](/flume/nrbc.png)

*Left: what is left in the box once the packet of nr2 has gone, against x. Right: the box-mean overpressure of nr3 against time, from
the box integrals; dotted, $e^{-Kt}$.*

**What it shows.**
- **It anchors the pressure; extrapolation does not.** The overpressure falls to 0.099 of itself in 10 acoustic crossings, close to
  $e^{-Kt}$ with $K = \sigma(1 - M^2)c/L$. Extrapolation keeps it at 1.000: a subsonic outflow with `extrapolation` has no pressure of
  its own.
- **It is not less reflecting than the kinds it would replace.** In this characteristic-upwind scheme `extrapolation` already makes the
  incoming characteristic jumps zero, and `inflow` makes them the free stream's: both let linear waves out (1.2e-4, 1.4e-5 left
  behind). The relaxation that holds the pressure reflects in proportion to $\sigma$ (2.9e-3 at 0.28, 7.6e-4 at 0.07, ratio 3.8 for
  a ratio of 4), as $\sigma c/(\omega L)$ predicts.
- **Its limits are exact**: at $\sigma = 0$, and at a supersonic outflow, it is extrapolation bit for bit.

Use it where a subsonic outflow or a far-field face must hold a pressure, with $\sigma$ as small as the transient allows. A supersonic
outflow or a far field that needs no pressure is served as well by `extrapolation`.

## Shock generator

`[bc_y_max] type = shock-generator` ([#121](https://github.com/szaghi/adam/issues/121) S3,
[boundary conditions](../../boundary-conditions#shock-generator-an-oblique-shock-through-a-y-face-121-s3)) brings an oblique shock in through
the top face. The legs measure the shock it makes against the oblique-shock relations (`shockgen_oracle.py`, which solves them on its own).

| Setup | |
|---|---|
| box, cells | $[0, 2]\times[0, 1]\times[0, 0.25]$, $64\times32\times8$, 8 blocks, periodic z |
| gas | Euler, $\gamma = 1.4$, $\rho = 1$, $p = 1/\gamma$ ($c = 1$), Mach 2 along x |
| boundary | `inflow` at $x = 0$; `extrapolation` at $x = 2$ and $y = 0$; `shock-generator` at $y = 1$, deflection $10^\circ$ from `shock_x` = 0.5 ($\beta = 39.314^\circ$, $p_2/p_1 = 1.70658$, $\rho_2/\rho_1 = 1.45843$) |
| scheme | WENO-U-5 characteristic, SSP-54, CFL 0.5, 2 ranks, to $t = 2.5$ (steady) |

| Leg | Case | Oracle | Measured (CPU; FNL the same to the digits shown) |
|---|---|---|---|
| sg1 | the shock by its deflection | the uniform region behind it (16 cells off the shock, 4 below the streamline from its foot, above $y = 0.2$): $p$ and $\rho$ within 3e-3 of Rankine–Hugoniot, the flow angle within $0.05^\circ$; the shock angle fitted from where $p$ crosses $(p_1 + p_2)/2$ along the rows within $0.3^\circ$ | $p$ 1.47e-3, $\rho$ 1.22e-3, angle $0.026^\circ$; $\beta$ $39.322^\circ$, the foot at $x = 0.496$ |
| sg2 | the same shock by its angle (`shock_angle` 39.313932) | as sg1 | the same digits |
| sg3 | a shock entering beyond the box (`shock_x` 3) | equal to an extrapolated top bit for bit | 0 of 64 datasets differ |
| sg4 | the regular reflection off a slip wall ([#121](https://github.com/szaghi/adam/issues/121) S4a): `wall-inviscid` at $y = 0$, $96\times48\times8$ cells, `shock_x` = 0.1 (the shock meets the wall at $x = 1.321$), $t = 5$ | behind the reflected shock (8 cells off it and off the wall) $p$ and $\rho$ within 4e-3 of the exact $p_3/p_1 = 2.80319$, the flow parallel to the wall within $0.1^\circ$; the reflected shock's angle to the wall within $0.3^\circ$ of $39.384^\circ$ | $p$ 2.10e-3, $\rho$ 1.50e-3, $0.062^\circ$; $39.437^\circ$, its foot at $x = 1.324$ |

![The shock generator at Mach 2](/flume/shockgen.png)

*Left: $p/p_1$ in the x-y plane of leg sg1; dashed, the shock of the theory from its foot; dotted, the streamline from the foot. Right:
$p/p_1$ and $\rho/\rho_1$ along the row $y = 0.484$ against the step of Rankine–Hugoniot at the theory's shock.*

**What it shows.**
- **The shock is where the relations put it.** Its angle is within $0.01^\circ$ of $\beta$, and it meets the top within a cell of
  `shock_x`; behind it the state is Rankine–Hugoniot's within $1.5\cdot10^{-3}$.
- **What the oracle leaves out, and why.** The captured shock leaves oscillations of a few $10^{-3}$ that decay over 16 cells behind
  it (the staircase of an oblique shock on a Cartesian grid). The jump of the ghost state at `shock_x` leaves an entropy streak along the
  streamline from the foot ($\rho$ 7e-3 off, the pressure continuous). Where the shock leaves through the extrapolated bottom, the
  corner downstream sits up to 4 % high: an extrapolated face crossed by an oblique shock is not clean. In a shock–boundary-layer
  interaction the shock meets the wall instead.

![The regular reflection off a slip wall](/flume/shockgen-reflect.png)

*Leg sg4. Left: $p/p_1$ in the x-y plane; dashed, the incident shock of the theory; dotted, the reflected one ($M_2 = 1.641$, turned back
by $10^\circ$). Right: $p/p_1$ along the row $y = 0.260$ against the steps of the exact regular reflection.*

- **The reflection is the exact one.** Behind the reflected shock the state is within $2\cdot10^{-3}$ of $p_3/p_1 = 2.803$, the flow
  parallel to the wall within $0.06^\circ$, and the reflected shock lies within $0.06^\circ$ of its angle, its foot within a cell of where the
  incident shock meets the wall.
- **It needs the steady state.** At $t = 2.5$ the region behind the reflected shock is still 6 % low near the outflow: the reflected shock
  forms at the wall and straightens outward, so the leg runs to $t = 5$.
