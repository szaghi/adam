# Non-reflecting boundary

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
