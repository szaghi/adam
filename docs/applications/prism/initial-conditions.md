# Initial conditions and sources

The field at t0 is set by `[initial_conditions] type`; the current density by the coils (`[coils_input]`) and, under
PIC, by the particles; an external field (`[external_fields]`) can act on the particles. The keys are in the
[input reference](./input#initial-conditions); this page gives what each one sets. Under the non-dimensional models
`plane_wave` and the pulses convert their SI inputs with $L_0$, $B_0$, $D_0 = B_0c_0\varepsilon_0$
([models](./models#non-dimensional-form)); `uniform_field` writes its values as given, and `rmf_field` uses its SI
amplitude and frequency with lengths in $L_0$ (no conversion, [#108](https://github.com/szaghi/adam/issues/108)).

After the initial condition, the t0 elliptic solves add the field of the initial charge and current
([numerics](./numerics#t0-elliptic-solves)), and the far-wave layer correction is applied.

## Field initial conditions

| `type` | Field at t0 |
|---|---|
| `vacuum` | every variable zero |
| `plane_wave` | $\mathbf D$, $\mathbf B$ of a plane wave: phase $\cos\!\big(\tfrac{2\pi}{\lambda}(k_x x + k_y y + k_z z)\big)$, $\mathbf D = D_w\,(k_z, k_x, k_y)\cos(\cdot)$, $\mathbf B = B_0\,(k_y, k_z, k_x)\cos(\cdot)$, $D_w = B_0 c_0\varepsilon_0$ |
| `gaussian_pulse` | a linearly polarized Gaussian wave packet (below) |
| `gaussian_derivative_pulse` | a single-cycle pulse with no zero-frequency content (below) |
| `rmf_field` | the rotating magnetic field at t = 0 (below) |
| `uniform_field` | $\mathbf B = (B_x, B_y, B_z)$, $\mathbf D = 0$ |

`riemann-problem`, `magnetic_nozzle` and `rmf_magnetic_nozzle` are refused at load
([#56](https://github.com/szaghi/adam/issues/56)). Every kind leaves the GLM scalars, $\mathbf J$ and $\rho$ at zero.

**`plane_wave` is transverse only along a coordinate axis.** With $(k_x,k_y,k_z) = (1,0,0)$ it gives $D_y$, $B_z$;
for an oblique direction $\mathbf D\cdot\mathbf k = k_xk_z + k_yk_x + k_zk_y \ne 0$ and, with $\mathbf B = B_0(k_y, k_z,
k_x)$, $\mathbf B\cdot\mathbf k = k_yk_x + k_zk_y + k_xk_z \ne 0$: neither initial field is divergence-free, and
without GLM cleaning nothing removes the $\nabla\cdot\mathbf B$.

### Gaussian pulse

Along the axis $a$ of `pulse_direction` ($\pm x$, $\pm y$, $\pm z$; the sign $\sigma_v = \pm1$ is the direction of
travel), with $s = x_a - r_0$:

$$E(s) = B_0 c_0\, e^{-s^2/(2\sigma^2)}\cos(k_0 s), \qquad E = 0 \ \text{for}\ |s| > 5\sigma,$$

the polarization pairs being $(E_y, B_z)$ along $x$, $(E_z, B_x)$ along $y$, $(E_x, B_y)$ along $z$, with
$D = \varepsilon_0 E$ and $B = \sigma_v E/c_0$.

### Gaussian-derivative pulse

The derivative of a Gaussian vector potential $A = B_0\sigma e^{1/2} e^{-s^2/(2\sigma^2)}$, same polarization:

$$E(s) = -B_0 c_0\, e^{1/2}\,\frac{s}{\sigma}\, e^{-s^2/(2\sigma^2)}, \qquad \max|B| = B_0\ \text{at}\ s = \mp\sigma .$$

$A$ vanishes on both sides, so the profile has zero area: no zero-frequency content, no net momentum to free particles,
no quasi-static field left in an absorbing layer. It is not truncated.

### Rotating magnetic field

With the rotation axis $\gamma$ (`[external_fields] RMF_rotation_axis`) and $(x, y)$ the two other coordinates in
cyclic order, $\theta = \operatorname{atan2}(y, x)$, $r = \sqrt{x^2+y^2}$, $\omega = 2\pi f$:

$$B_x = B_{\rm RMF}, \quad B_y = 0, \quad D_\gamma = \varepsilon_0\,\omega\, r\, B_{\rm RMF}\cos\theta,$$

the field at t = 0 of $\mathbf B$ rotating at $\omega$ about $\gamma$, with its induced electric field.

## Coils

A coil is a prescribed current density $\mathbf J$ stamped on the grid, one of four kinds (`coil_type`), with its plane
normal (or axis) among $\pm x$, $\pm y$, $\pm z$:

| Kind | Geometry |
|---|---|
| `rectangular` | a rectangular loop of sides `lx`, `ly` |
| `circular` | a circular loop of radius `r_coil` |
| `solenoid` | `windings` turns of radius `r_coil` over a length `l_solenoid` |
| `helicon` | a polyline antenna through up to 50 points |

**Stamping.** Every kind builds a smooth vector potential $\mathbf A$ of the coil path, with a Gaussian cross-section of
width `sigma`, and sets $\mathbf J = \nabla_h\times\mathbf A$ with the operator of `[fdv] fdv_scheme` (FD or FV). So
$\nabla_h\cdot\mathbf J = 0$ to round-off on a uniform grid, and the coil injects no charge. The first coil's
amplitude is calibrated on its `Amplitude`; **the coils after the first reuse the first coil's calibrated amplitude**,
and **a circular coil's amplitude is never set, so it injects no current** (defects,
[#108](https://github.com/szaghi/adam/issues/108)). The per-coil density is written to the field files as
`coil_NN_j_vec_1/2/3`.

**Waveform.** At time $t$ (the stage time $t + \gamma_s\Delta t$ inside a Runge–Kutta step) a coil adds
$A\,g(s)\cos\theta$ times its stamped density to $\mathbf J$:

$$g(s) = 10s^3 - 15s^4 + 6s^5, \qquad s = \operatorname{clamp}(t/t_d,\ 0,\ 1), \qquad
\theta = \begin{cases}\phi & \texttt{DC\_current}\\ \omega\,(t - t_d) + \phi & \texttt{AC\_current}\end{cases}$$

with $t_d$ = `[coils_input] time_delay`, $\omega = 2\pi\,$`Frequency`, and $\phi$ = `Phase` **in degrees** (0 for DC).
$g$ ramps the current smoothly from zero at t0 to full at $t_d$, with zero first and second derivatives at both ends.
With `time_delay <= 0` the current is full at t0, and its magnetic field is solved for at t0.

In the cells a coil occupies, $\mathbf J$ is reset at every evaluation: a PIC deposit there stops the run.

## External fields

`[external_fields] external_fields_applied = RMF` (a field of amplitude `RMF_B_amplitude` rotating at `RMF_frequency`
about `RMF_rotation_axis`) or `uniform_field` adds a prescribed field **to the fields the particles see**: before every
gather it is added to the state, and subtracted after. It is not part of the evolved state, so in an EM run (no
particles) it does not change the dynamics; it is, however, added to every field file and slice written, which show the
total field. `magnetic_nozzle` and `rmf_and_magnetic_nozzle` are refused at initialization.
