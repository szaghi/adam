# Initial conditions

The initial state is set by the `[initial_conditions]` section: `type` selects one of the sixteen initial conditions
below, and the primitive states of the regions are given in `[initial_conditions_region_N]` sections. The value of
every key is checked: an unknown `type`, a missing key or an inconsistent value stops the run with a message naming the
accepted values.

| `type` | Euler | MHD (`none`) | MHD (`glm`) | Region sections |
|---|:-:|:-:|:-:|---|
| `uniform` | ✓ | ✓ | ✓ | 1 |
| `isentropic-vortex` | ✓ | ✗ | ✗ | 1 |
| `riemann-problem` | ✓ | ✓ | ✓ | `regions_number` |
| `shu-osher` | ✓ | ✗ | ✗ | 2 |
| `rotated-riemann` | ✓ | ✓ | ✓ | 2 |
| `glm-pulse`, `divb-peak`, `mhd-linear-wave`, `mhd-cpaw`, `mhd-vortex`, `mhd-rotor`, `field-loop` | ✗ | ✓ | ✓ | 1 |
| `orszag-tang` | ✗ | ✓ | ✓ | 0 |
| `linear` (verification) | ✓ | ✓ | ✓ | 1 |
| `polytrope` (self-gravity) | ✓ | ✗ | ✗ | 1 (the atmosphere) |
| `singular-isothermal-sphere` (self-gravity) | ✓ | ✗ | ✗ | 1 (the temperature) |
| `bonnor-ebert` (self-gravity, `[eos] type = isothermal`) | ✓ | ✗ | ✗ | 1 (the velocity) |
| `boss-bodenheimer` (self-gravity, `[eos]` isothermal or barotropic) | ✓ | ✗ | ✗ | 1 (the ambient medium) |

"MHD" in the error messages means `[physics] physical_model = mhd-ideal`. `glm` versus `none` is chosen separately by
`[mhd] divergence_control`.

## Conventions

**State vectors.** The conservative variables are
$\mathbf q = (\rho, \rho u, \rho v, \rho w, E\,[, B_x, B_y, B_z\,[, \psi]])$, with
$n_v = 5$ (Euler), $8$ (MHD, `[mhd] divergence_control = none`) or $9$ (MHD with GLM).
The ratio of specific heats is $\gamma = c_p/c_v$ and the gas constant is $R = c_p - c_v$.

**Primitive to conservative conversion** is used by every initial condition and by the inflow boundary condition.
MHD uses rationalised units, so no $4\pi$ appears:

$$
E = \frac{p}{\gamma-1} + \tfrac12\rho\,(u^2+v^2+w^2)\;\Big[+\;\tfrac12\,(B_x^2+B_y^2+B_z^2)\Big]
$$

Sources: Euler at `common/adam_flume_euler_library.F90`, MHD at `common/adam_flume_mhd_library.F90`.
The MHD sums go through `mhd_sum3`, which sorts the three terms before adding them so that a rotated problem gives
bitwise-identical results.

**$\psi$ at $t=0$ is always zero.** The model dispatcher `primitive_state_to_conservative` sets `q(IQ_PSI) = 0` for
`MODEL_MHD_GLM`. Every initial condition goes through this routine,
either directly or through the pre-converted region states `q_region`. The one exception is `mhd-linear-wave`, which
adds its eigenvector only to components 1–8, so $\psi$ also stays 0 there.

**Coordinates.** Every initial condition is a **point value at the cell centre**, not a cell average:
`field%x_cell(i,b)`, `y_cell(j,b)`, `z_cell(k,b)`, and
$x_c = x_{\min} - \Delta x/2 + \text{lin\_space}(i)$.
`set_initial_conditions` fills only the interior cells $1..n_i \times 1..n_j \times 1..n_k$. The ghost cells come from
the `update_ghost` call that follows.

**Keys common to every type.** Both are read from `[initial_conditions]`:

| key | type | meaning |
|---|---|---|
| `type` | string (control characters are blanked by `strip_control`, so CRLF files parse) | one of the 25 names below. Any other value stops the run. |
| `amr_iterations` | integer, **required** | the number of init-time AMR passes (set the IC, refine, repeat). The value is clamped to $\ge 0$. The loop is at `cpu/adam_flume_cpu_object.F90` and `fnl/adam_flume_fnl_object.F90`. When the value is $>0$, every refined non-null axis must have an even block cell count of at least `2 ngc`, otherwise `error_stop`. |

**Region sections** are named `[initial_conditions_region_N]`, with $N = 1..$`regions_number`. They hold the primitive keys `r, u, v, w, p`, plus `bx, by, bz` for either MHD
model. Every one of these keys is required. $\psi$ has no key. The number of region sections
read depends on the type:

| type | regions read | extent keys |
|---|---|---|
| `uniform`, `isentropic-vortex`, `glm-pulse`, `divb-peak`, `mhd-linear-wave`, `mhd-cpaw`, `mhd-vortex`, `mhd-rotor`, `field-loop`, `linear` | 1 | none |
| `shu-osher`, `rotated-riemann` | 2 | none (any extent keys present are ignored) |
| `riemann-problem` | `regions_number` (from the INI) | `emin_x, emin_y, emin_z, emax_x, emax_y, emax_z` (required) |
| `orszag-tang` | 0 | none |

---



## `uniform`: seeded uniform state

*Model:* any. *Keys:* `[initial_conditions] s` (real, required), plus region 1.

The state of region 1 is imposed everywhere, with density and pressure perturbed by a relative amplitude $s$:

$$
\rho = \rho_1\,(1 + s\,h_1),\qquad p = p_1\,(1 + s\,h_2),\qquad h_1,h_2\in[-1,1)
$$

The pair $(h_1, h_2)$ comes from `hash_cell(code, i, j, k)`. The hash is a chained integer mix of the
block Morton code and the cell indexes: xorshift steps (13, 7, 17), then a multiplication by 39083855 modulo $2^{52}$
in 26-bit limbs. It uses no RNG state, so the values do not depend on the rank decomposition or the backend. The
velocity and field are those of region 1.

*Validation:* none on `s`. Setting $|s|\ge 1$ can produce $\rho\le 0$ or $p\le 0$ without any error at this stage.

*Test input:* `src/tests/flume/verification/uniform/input.ini`
```ini
[initial_conditions]
type           = uniform
amr_iterations = 0
s              = 0.0 ; relative amplitude of the seeded density/pressure perturbation

[initial_conditions_region_1]
r = 1.2
u = 10.0
v = 0.0
w = 0.0
p = 101325.0
```
MHD with GLM, from `src/tests/flume/verification/mhd/plumbing/mhd-uniform-glm.ini`:
```ini
[initial_conditions]
type           = uniform
amr_iterations = 0
s              = 0.05 ; relative amplitude of the seeded density/pressure perturbation

[initial_conditions_region_1]
r  = 1.0
u  = 0.3
v  = -0.2
w  = 0.1
p  = 1.0
bx = 0.8
by = 0.5
bz = -0.3
```
The companion file is `mhd-uniform-none.ini`, and `src/tests/flume/regression/amr-periodic-reflux/input.ini` and
`regression/uniform-amr-mhd/input.ini` use this type too.

## `isentropic-vortex`

*Reference:* Shu 1998, ICASE 97-65, §5.1, generalised to any free stream.
*Model:* **Euler only**. Any other model stops the run with `requires [physics].(physical_model) = euler`.
*Keys:* `x0, y0, radius, strength`, plus region 1 (`r,u,v,w,p`, the free stream).
*Validation:* `radius <= 0` stops the run.

Let $\varepsilon$ be `strength`, $R_v$ be `radius`, $\tilde x = (x-x_0)/R_v$, $\tilde y = (y-y_0)/R_v$ and
$e = \exp\!\big(\tfrac12(1-\tilde x^2-\tilde y^2)\big)$. With $T_0 = p_1/\rho_1$:

$$
T = T_0 - \frac{\gamma-1}{\gamma}\,\frac{\varepsilon^2}{8\pi^2}\,e^2,\qquad
\rho = \rho_1\Big(\frac{T}{T_0}\Big)^{1/(\gamma-1)},\qquad p = \rho\,T
$$

$$
u = u_1 - \frac{\varepsilon}{2\pi}\,\tilde y\,e,\qquad v = v_1 + \frac{\varepsilon}{2\pi}\,\tilde x\,e,\qquad w = w_1
$$

Here $T$ stands for $p/\rho$, not the thermodynamic temperature. The flow is isentropic ($p/\rho^\gamma$ is constant)
and it satisfies the radial balance $\partial_r p = \rho v_\theta^2/r$ for any $T_0$ and any $R_v$, so it is an exact
steady solution convected by the free stream. The profile is 2-D; $z$ is not used.

*Test input:* `src/tests/flume/verification/vortex/vortex-n064.ini` (and `-n128`, `-n256`, `regression/vortex-periodic`)
```ini
[initial_conditions]
type           = isentropic-vortex
amr_iterations = 0
x0             = 0.5
y0             = 0.5
radius         = 0.07
strength       = 5.0

[initial_conditions_region_1]
r = 1.0
u = 1.0
v = 1.0
w = 0.0
p = 1.0
```

## `riemann-problem`: piecewise-constant axis-aligned regions

*Model:* any. *Keys:* `[initial_conditions] regions_number` (integer, required, must be $\ge 1$). Each
region section also needs the six extent keys.

A cell with centre $\mathbf c$ belongs to region $r$ if, for every axis $d$, $e_{\min,d}^{(r)} < c_d \le e_{\max,d}^{(r)}$. The region sections are scanned in order and **the first match wins**, so overlapping
regions are allowed and region 1 has priority. A cell covered by no region is **fatal**: `cell center (x, y, z) is covered
by no [initial_conditions_region_*]`. The half-open interval puts a cell centre that sits exactly on an
interface into the lower region.

The $z$ extent must also bracket the cell centres, including on 2-D grids with `null_z`. The tests use
`emin_z = 0, emax_z = 1`.

*Test input (Sod along x):* `src/tests/flume/verification/sod/sod-x.ini`
```ini
[initial_conditions]
type           = riemann-problem
regions_number = 2
amr_iterations = 0

[initial_conditions_region_1]
r      = 1.0
u      = 0.0
v      = 0.0
w      = 0.0
p      = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 0.5
emax_y = 1.0
emax_z = 1.0

[initial_conditions_region_2]
r      = 0.125
u      = 0.0
v      = 0.0
w      = 0.0
p      = 0.1
emin_x = 0.5
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
```
*MHD (Brio-Wu):* `src/tests/flume/regression/brio-wu/input.ini`
```ini
[initial_conditions]
type = riemann-problem
regions_number = 2
amr_iterations = 0

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 1.0
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 0.5
emax_y = 1.0
emax_z = 1.0
bx = 0.75
by = 1.0
bz = 0.0

[initial_conditions_region_2]
r = 0.125
u = 0.0
v = 0.0
w = 0.0
p = 0.1
emin_x = 0.5
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
bx = 0.75
by = -1.0
bz = 0.0
```
Other inputs that use this type: `sod-{y,z}.ini`, `sod-wall-{x,y,z}.ini`, `regression/rj2a-{x,y,z}`,
`regression/sod-{x,y,z}`, and the multi-realm `*-r1.ini`/`*-r2.ini` files.

## `shu-osher`: shock and density-wave interaction

*Reference:* Shu & Osher 1989, J. Comput. Phys. 83, example 8.
*Model:* **Euler only**. Any other model stops the run with `requires [physics].(physical_model) = euler`.
*Keys:* `axis`, `interface`,
`rho_amplitude`, `rho_wavenumber`, plus regions 1 and 2, without extents.
*Validation:* $\rho_2 - |A| \le 0$ stops the run with `rho_amplitude must be smaller than the density of
[initial_conditions_region_2]`.

Let $s$ be the cell-centre coordinate along `axis`:

$$
\mathbf q = \begin{cases}\mathbf q_1 & s \le s_I\\[2pt]
\mathbf q\big(\rho_2 + A\sin(k\,s),\;\mathbf u_2,\;p_2\big) & s > s_I\end{cases}
$$

The phase is $k\,s$ in the **absolute** coordinate, not $k(s-s_I)$.

*Test input:* the RV-3 `shu-osher` leg of `src/tests/flume/verification/riemann-flux/check.sh` (function `set_shu_osher`) rewrites `verification/sod/sod-{x,y,z}.ini`, for example:
```ini
[initial_conditions]
axis           = x
interface      = -4.0
rho_amplitude  = 0.2
rho_wavenumber = 5.0
type           = shu-osher
regions_number = 2
amr_iterations = 0

[initial_conditions_region_1]
r      = 3.857143
u      = 2.629369
v      = 0.0
w      = 0.0
p      = 10.33333
emin_x = 0.0
emin_y = 0.0
emin_z = 0.0
emax_x = 0.5
emax_y = 1.0
emax_z = 1.0

[initial_conditions_region_2]
r      = 1.0
u      = 0.0
v      = 0.0
w      = 0.0
p      = 1.0
emin_x = 0.5
emin_y = 0.0
emin_z = 0.0
emax_x = 1.0
emax_y = 1.0
emax_z = 1.0
```
`regions_number` and the `emin_*`/`emax_*` keys are left over from the Sod input and are ignored for this type.

## `glm-pulse`: Gaussian pulse in one field component

*Reference:* the GLM d'Alembert test (issue #41, MV-3).
*Model:* **MHD** (`mhd-ideal`, with or without GLM). Euler stops the run with `requires [physics].(physical_model) = mhd-ideal`.
*Keys:* `pulse_axis`, `pulse_center`, `pulse_width`, `pulse_amplitude`, plus region 1 (8 keys).
*Validation:* `pulse_width <= 0` stops the run. No check stops a pulse along a null axis.

With $a$ the chosen axis and $s$ the cell-centre coordinate along it:

$$
B_a = B_{a,1} + A\,\exp\!\Big(-\big(\tfrac{s-s_c}{\sigma}\big)^2\Big),
$$

All other primitives are those of region 1, and $\psi = 0$. When $a$ is the axis the pulse varies along, the pulse is
a pure $\nabla\!\cdot\!\mathbf B$ error.

*Test input:* `src/tests/flume/regression/glm-pulse/input.ini`
```ini
[initial_conditions]
type = glm-pulse
amr_iterations = 0
pulse_axis = x
pulse_center = 0.5
pulse_width = 0.1
pulse_amplitude = 0.1

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 1.0
bx = 0.0
by = 0.0
bz = 0.0
```
(The verification variant is generated by `src/tests/flume/verification/mhd/glm-pulse/make_glm_pulse.py`.)

## `divb-peak`: Dedner et al. peak in $B_x$

*Reference:* after Dedner et al. 2002, §5.1 (issue #41, MV-10); the support convention of `peak_radius` (a radius here) is being checked against the paper.
*Model:* **MHD**. *Keys:* `peak_x0, peak_y0, peak_radius, peak_amplitude`, plus region 1.
*Validation:* `peak_radius <= 0` stops the run.

With $\tilde r = \sqrt{(x-x_0)^2+(y-y_0)^2}\,/\,R$:

$$
B_x = B_{x,1} + A\,(1-\tilde r^2)^2 \quad\text{for }\tilde r<1,\qquad B_x = B_{x,1}\ \text{otherwise.}
$$

$\nabla\!\cdot\!\mathbf B\neq 0$ at $t=0$ by construction.

*Test input:* generated by `src/tests/flume/verification/mhd/divb-peak/make_divb_peak.py`, for example:
```ini
[initial_conditions]
type = divb-peak
amr_iterations = 1
peak_x0 = 0.0
peak_y0 = 0.0
peak_radius = 0.125
peak_amplitude = 0.28209479177387814

[initial_conditions_region_1]
r = 1.0
u = 1.0
v = 1.0
w = 0.0
p = 6.0
bx = 0.0
by = 0.0
bz = 0.28209479177387814
```

## `mhd-linear-wave`: small-amplitude MHD eigenmode

*Reference:* Stone et al. 2008, ApJS 178, §8.2 (issue #41, MV-5). The eigenvector
normalisation is that of Stone et al. 2008, appendix B.
*Model:* **MHD**. *Keys:* `wave`, `wave_angle` (degrees, in the
x-y plane), `wave_amplitude`, `wavelength`, plus region 1, which is the background $\mathbf q_0$ in
the **global** frame.
*Validation:* `wavelength <= 0` stops the run.

The background is rotated into the wave frame, $u_n = c\,u + s\,v$, $u_t = -s\,u + c\,v$, and likewise for $B$, with
$c = \cos\alpha$, $s = \sin\alpha$. The right eigenvectors of the $x$-direction system of that frame are
then computed (`mhd_eigenvectors`, `d = 1`; tangents are cyclic, so $t_1 = y'$, $t_2 = z$). The column $k$ is chosen
from the wave, in the order $(u_n-c_f, u_n-c_a, u_n-c_s, u_n, u_n+c_s, u_n+c_a, u_n+c_f, B_n)$:

| `wave` | column | eigenvalue |
|---|---|---|
| `fast` | 7 | $u_n+c_f$ |
| `alfven` | 6 | $u_n+c_a$ |
| `slow` | 5 | $u_n+c_s$ |
| `entropy` | 4 | $u_n$ |

The momentum and field components of $\mathbf r_k$ are rotated back, $r_x = c\,r_n - s\,r_t$, $r_y = s\,r_n + c\,r_t$. Each cell then gets:

$$
\mathbf q(x,y) = \mathbf q_0 + A\,\sin\!\Big(\frac{2\pi\,(x\cos\alpha + y\sin\alpha)}{\lambda}\Big)\,\mathbf r_k
\qquad(\text{components }1..8;\ \psi=0)
$$

The perturbation is added to the **conservative** vector, and $A$ is measured in the normalisation of $\mathbf r_k$.

*Test input:* generated by `src/tests/flume/verification/mhd/linear-wave/make_linear_wave.py`, for example:
```ini
[initial_conditions]
type = mhd-linear-wave
amr_iterations = 0
wave = alfven
wave_angle = 0.0
wave_amplitude = 1.0e-7
wavelength = 1.0

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 0.6
bx = 1.0
by = 1.4142135623730951
bz = 0.5
```

## `mhd-cpaw`: circularly polarised Alfvén wave

*Reference:* Tóth 2000, J. Comput. Phys. 161 (issue #41, MV-6). This is an exact nonlinear
solution of any amplitude.
*Model:* **MHD**. *Keys:* `polarisation` (`right` gives $h=+1$, `left` gives $h=-1$; any other value is
fatal), `wave_angle`, `wave_amplitude`, `wavelength`, `b_par`, plus region 1.
*Validation:* `wavelength <= 0` is fatal and `b_par == 0` is fatal. Region 1 **must** have
`u = v = w = bx = by = bz = 0`, otherwise the run stops. Only $\rho_1$ and $p_1$ are used.

With $\alpha$ = `wave_angle`, $\mathbf n = (\cos\alpha,\sin\alpha,0)$, $\mathbf t_1 = (-\sin\alpha,\cos\alpha,0)$,
$\mathbf t_2 = \hat z$ and $\phi = 2\pi(x\cos\alpha + y\sin\alpha)/\lambda$:

$$
\mathbf B_\perp = A\big(\sin\phi\,\mathbf t_1 + h\cos\phi\,\mathbf t_2\big),\qquad
\mathbf B = b_\parallel\,\mathbf n + \mathbf B_\perp,\qquad
\mathbf u = -\,\mathrm{sign}(b_\parallel)\,\frac{\mathbf B_\perp}{\sqrt{\rho_1}}
$$

$\rho = \rho_1$ and $p = p_1$ everywhere. $|\mathbf B_\perp|$ is uniform, so the total pressure is uniform too. The wave
travels along $+\mathbf n$ at $|b_\parallel|/\sqrt{\rho_1}$, and after one period $\lambda\sqrt{\rho_1}/|b_\parallel|$
the exact solution is the initial state.

*Test input:* generated by `src/tests/flume/verification/mhd/cpaw/make_cpaw.py`, for example:
```ini
[initial_conditions]
type = mhd-cpaw
amr_iterations = 0
polarisation = left
b_par = 1.0
wave_angle = 63.43494882292201
wave_amplitude = 0.1
wavelength = 1.0

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 0.1
bx = 0.0
by = 0.0
bz = 0.0
```

## `mhd-vortex`: magnetised vortex

*Reference:* Balsara 2004, ApJS 151 (issue #41, MV-7).
*Model:* **MHD**. *Keys:* `x0, y0, radius, kappa, mu`, plus region 1 (the free stream).
*Validation:* `radius <= 0` is fatal. Region 1 must have `bx = by = 0`, otherwise `error_stop` with the message
"an equilibrium only without an in-plane free-stream field". `bz` is kept.

With $\tilde x = (x-x_0)/R_v$, $\tilde y = (y-y_0)/R_v$, $\tilde r^2 = \tilde x^2+\tilde y^2$ and
$e = \exp\!\big(\tfrac12(1-\tilde r^2)\big)$:

$$
u = u_1 - \frac{\kappa}{2\pi}\tilde y\,e,\quad v = v_1 + \frac{\kappa}{2\pi}\tilde x\,e,\quad
B_x = -\frac{\mu}{2\pi}\tilde y\,e,\quad B_y = \frac{\mu}{2\pi}\tilde x\,e,
$$
$$
p = p_1 + \frac{\mu^2(1-\tilde r^2) - \rho_1\kappa^2}{8\pi^2}\,e^2,\qquad \rho=\rho_1,\ w=w_1,\ B_z=B_{z,1}.
$$

The field is divergence-free pointwise, and the vortex is an exact steady equilibrium convected by the free stream.

*Test input:* generated by `src/tests/flume/verification/mhd/vortex/make_mhd_vortex.py`, for example:
```ini
[initial_conditions]
type = mhd-vortex
amr_iterations = 0
x0 = 0.0
y0 = 0.0
radius = 1.0
kappa = 1.0
mu = 1.0

[initial_conditions_region_1]
r = 1.0
u = 1.0
v = 1.0
w = 0.0
p = 1.0
bx = 0.0
by = 0.0
bz = 0.0
```

## `orszag-tang`

*Reference:* Stone et al. 2008, §8.4 (issue #41, MV-12).
*Model:* **MHD**. *Keys:* none beyond `type` and `amr_iterations`. `regions_number = 0`, so no region
section is read. The unit period is **assumed**: the domain is not checked against $[0,1]^2$.

With $s = 2x-1$ and $t = 2y-1$:

$$
\rho=\frac{25}{36\pi},\quad p=\frac{5}{12\pi},\quad u=\sin(\pi t),\quad v=-\sin(\pi s),\quad w=0,
$$
$$
B_x=\frac{\sin(\pi t)}{\sqrt{4\pi}},\quad B_y=\frac{\sin(2\pi s)}{\sqrt{4\pi}},\quad B_z=0.
$$

These equal the textbook $u=-\sin 2\pi y$, $v=\sin 2\pi x$, $B_x=-\sin(2\pi y)/\sqrt{4\pi}$, $B_y=\sin(4\pi x)/\sqrt{4\pi}$.
The code evaluates them as odd functions of $s, t$ so that the state is bitwise symmetric under the 180° rotation
about $(1/2,1/2)$ when cell centres are binary fractions.

*Test input:* `src/tests/flume/regression/orszag-tang/input.ini`
```ini
[initial_conditions]
type = orszag-tang
amr_iterations = 0
```

## `mhd-rotor`

*Reference:* Balsara & Spicer 1999 and Tóth 2000, first rotor (issue #41, MV-13).
*Model:* **MHD**. *Keys:* `x0, y0, r0, r1, rho_in, v0`, plus region 1 (ambient).
*Validation:* the run stops unless $0<r_0<r_1$, with the message `needs 0 < r0 < r1`, and stops if
$\rho_{in}\le 0$.

With $\Delta x = x-x_0$, $\Delta y = y-y_0$ and $r = \sqrt{\Delta x^2+\Delta y^2}$:

$$
r<r_0:\quad \rho=\rho_{in},\quad (u,v) = (u_1,v_1) + \frac{v_0}{r_0}(-\Delta y,\ \Delta x)
$$
$$
r_0\le r<r_1:\quad f=\frac{r_1-r}{r_1-r_0},\quad \rho=\rho_1+(\rho_{in}-\rho_1)f,\quad
(u,v) = (u_1,v_1) + \frac{f\,v_0}{r}(-\Delta y,\ \Delta x)
$$

Pressure, $w$ and $\mathbf B$ are those of region 1 everywhere, and the rotation velocity is **added** to the ambient one.

*Test input:* `src/tests/flume/regression/rotor/input.ini`
```ini
[initial_conditions]
type = mhd-rotor
amr_iterations = 0
x0 = 0.5
y0 = 0.5
r0 = 0.1
r1 = 0.115
rho_in = 10.0
v0 = 2.0

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 1.0
bx = 1.4104739588693909
by = 0.0
bz = 0.0
```

## `field-loop`: advected field loop

*Reference:* Gardiner & Stone 2005; Mignone & Tzeferacos 2010, §4.4.1 (issue #41, MV-9).
*Model:* **MHD**. *Keys:* `x0, y0, loop_radius, loop_amplitude`. The centre keys are
plain `x0`/`y0`, not `loop_x0`. Region 1 is also required.
*Validation:* `loop_radius <= 0` is fatal.

The field comes from the vector potential $A_z = A\,(R-r)$ for $r<R$, where $r=|(x-x_0,\,y-y_0)|$. Its curl is
computed analytically at the cell centre:

$$
0<r<R:\quad B_x = B_{x,1} - A\,\frac{y-y_0}{r},\qquad B_y = B_{y,1} + A\,\frac{x-x_0}{r};\qquad\text{otherwise }\mathbf B=\mathbf B_1.
$$

The centre cell itself ($r = 0$ exactly) is excluded. The field is point-sampled and not built from a discrete curl of
$A_z$, so the discrete $\nabla\!\cdot\!\mathbf B$ is non-zero at the loop edge (the discontinuity at $r=R$). With a
uniform $w\ne 0$ this discrete divergence feeds $B_z$.

*Test input:* `src/tests/flume/regression/field-loop/input.ini`
```ini
[initial_conditions]
type = field-loop
amr_iterations = 0
x0 = 0.0
y0 = 0.0
loop_radius = 0.3
loop_amplitude = 1.0e-3

[initial_conditions_region_1]
r = 1.0
u = 2.0
v = 1.0
w = 1.0
p = 1.0
bx = 0.0
by = 0.0
bz = 0.0
```

## `rotated-riemann`: 1-D Riemann problem rotated in x-y, periodic along its normal

*Reference:* Tóth 2000, J. Comput. Phys. 161, §6.3.2 (issue #41, MV-8).
*Model:* **any**. It is not restricted to MHD; for Euler only `u, v` are rotated.
*Keys:* `normal_x, normal_y, interface_1, interface_2, period, interface_2_width`, plus regions 1 and 2.
The region states are given **in the frame of the normal**: `u` and `bx` are the normal components, `v` and `by` the
tangential ones, `w` and `bz` lie along $z$.
*Validation*:
- `normal_x = normal_y = 0` is fatal;
- `period <= 0` is fatal;
- the run stops unless $I_1<I_2$, $w\ge 0$ and $I_2+w<I_1+P$.

**Rotation to x, y**. With $\hat{\mathbf n} = \mathbf n/|\mathbf n|$ and
$\hat{\mathbf t} = (-\hat n_y,\hat n_x)$: $u = u_n\hat n_x - u_t\hat n_y$, $v = u_n\hat n_y + u_t\hat n_x$, and the
same for $B$.

**Region assignment**. Let $s = n_x x + n_y y$ (unnormalised), $f = \mathrm{mod}(s-I_1,\,P)$ and
$d = I_2 - I_1$:

$$
\mathbf q = \begin{cases}
\mathbf q_2 & 0 < f \le d\\
\mathbf q\big(\mathbf P_2 + (\mathbf P_1-\mathbf P_2)\,\tfrac{f-d}{w}\big) & d < f < d+w \quad(\text{linear ramp in the global-frame primitives})\\
\mathbf q_1 & \text{otherwise (including } f=0 \text{ and } f\ge d+w)
\end{cases}
$$

With $w = 0$ the ramp branch is empty, so there is no division by zero. With integer normal components on a square
domain of side $P$, the state is periodic in both $x$ and $y$.

*Test input:* `src/tests/flume/regression/rotated-shock-tube/input.ini`
```ini
[initial_conditions]
type = rotated-riemann
amr_iterations = 0
normal_x = 1.0
normal_y = 2.0
interface_1 = 0.25390625
interface_2 = 0.50390625
interface_2_width = 0.4
period = 1.0

[initial_conditions_region_1]
r = 1.0
u = 10.0
v = 0.0
w = 0.0
p = 20.0
bx = 1.4104739588693909
by = 1.4104739588693909
bz = 0.0

[initial_conditions_region_2]
r = 1.0
u = -10.0
v = 0.0
w = 0.0
p = 1.0
bx = 1.4104739588693909
by = 1.4104739588693909
bz = 0.0
```

## `linear`: a linear field (verification)

Every conservative variable is linear in space,

$$
\mathbf q(\mathbf x) = \mathbf q_1 \,(1 + \mathbf g \cdot \mathbf x),
$$

where $\mathbf q_1$ is the conservative state of `[initial_conditions_region_1]` and $\mathbf g$ =
(`gradient_x`, `gradient_y`, `gradient_z`) a relative gradient (inverse length; the [`[reference]`](./input#reference-optional-dimensional-input)
layer scales it by $L_0$). The three keys are required.

It is not a flow. It is the field on which every ghost fill FLUME uses is exact: the same-level copy, the 2:1
restriction (the mean of the fine cells) and the tricubic coarse-to-fine fill reproduce a linear field to round-off,
and every boundary condition maps it to a known value. The ghost probe GP ([verification](./verification/infrastructure/ghost-probe),
`verification/ghosts/`) runs it for one negligible step (`CFL = 1e-30`) and checks every face and edge ghost against
that value. The state is admissible only while $1 + \mathbf g \cdot \mathbf x > 0$ over the domain; nothing checks it,
since the probe takes no real step. Give the velocity non-zero components, so that every wall sign flip shows.

```ini
[initial_conditions]
type = linear
gradient_x = 0.05
gradient_y = 0.11
gradient_z = 0.07
amr_iterations = 1

[initial_conditions_region_1]
r = 1.0
u = 0.3
v = -0.2
w = 0.1
p = 1.0
```

## `polynomial`: a polynomial field (verification)

The `linear` field raised to a per-direction degree $n$ = `degree` (an integer, at least 1; $n = 1$ is `linear`):

$$
\mathbf q(\mathbf x) = \mathbf q_1 \,\Big(1 + \sum_d g_d\, x_d^{\,n}\Big).
$$

The keys are those of `linear` plus `degree`, all required. The copies and every boundary condition still map it to
a known value (they copy cells), but each seam operator reproduces it only up to its degree: the mean restriction and
the restriction-compatible fill $n = 1$, the tricubic fill $n = 3$, the `lagrange6` restriction and the `quintic` fill
$n = 5$ ([#68](https://github.com/szaghi/adam/issues/68) D1). GP with `--degree N` runs it, so the probe checks the
high-order seam operators in situ: their footprints, shifts and every exchange path. $g_d$ has dimension $L^{-n}$, which
the [`[reference]`](./input#reference-optional-dimensional-input) layer does not convert: a reference with `polynomial`
is fatal.

## `sine-wave`: a small sine wave on a uniform state (verification)

A sine of amplitude $A$ along the unit normal $\mathbf n = (\cos\theta, \sin\theta, 0)$ in the x-y plane, on the state
of `[initial_conditions_region_1]` ($\rho_0$, $\mathbf u_0$, $p_0$), with $s = A\sin(2\pi\,\mathbf n\cdot\mathbf x/\lambda)$:

| `wave_mode` | Perturbation | `wave_amplitude` |
|---|---|---|
| `shear` | $\mathbf u = \mathbf u_0 + s\,\mathbf t$, $\mathbf t = (-\sin\theta, \cos\theta, 0)$ | a velocity |
| `acoustic` | $\rho = \rho_0(1 + s)$, $\mathbf u = \mathbf u_0 + a_0 s\,\mathbf n$, $p = p_0(1 + \gamma s)$: the right-running linear acoustic wave | relative |
| `entropy` | $\rho = \rho_0(1 + s)$: a density wave at uniform velocity and pressure, advected exactly (under a constant gravity it falls freely: the convergence case of the gravity sources, [verification](./verification/gravity/energy-order)) | relative |
| `magnetic` (MHD) | $\mathbf B = \mathbf B_0 + s\,\mathbf t$: a transverse field, which decays by resistivity and, with $\mathbf B_0\cdot\mathbf n \ne 0$, drives a damped Alfvén wave | a field |
| `jeans` (self-gravity) | $\rho = \rho_0(1 + s)$, $p = p_0(1 + \gamma s)$ and, with $\omega^2 = a_0^2k^2 - 4\pi G\rho_0$, $k = 2\pi/\lambda$: for $\omega^2 \ge 0$ the right-running Jeans wave, $\mathbf u = \mathbf u_0 + (\omega/k)\,s\,\mathbf n$; for $\omega^2 < 0$ the growing mode, $\mathbf u = \mathbf u_0 + (\sigma/k)A\cos(2\pi\,\mathbf n\cdot\mathbf x/\lambda)\,\mathbf n$, $\sigma^2 = -\omega^2$ ([#81](https://github.com/szaghi/adam/issues/81), [verification](./verification/self-gravity/)) | relative |

Keys: `wave_mode` (`shear`, `acoustic`, `entropy`, `magnetic` or `jeans`, other → fatal; `magnetic` needs `mhd-ideal`, `jeans` needs `[gravity] self_gravity`), `wave_angle` $\theta$ (degrees), `wave_amplitude` $A$,
`wavelength` $\lambda$ (`> 0`). The values are point values at the cell centres. It is the initial state of the
viscous verifications VV-1 (shear-wave decay) and VV-2 (viscous-thermal acoustic attenuation)
([verification](./verification/dissipative/viscous-waves)) and of the Ohmic ones VV-5 and VV-6
([verification](./verification/dissipative/ohmic-alfven)); the wave is periodic in a box whose sides are multiples
of $\lambda/\cos\theta$ and $\lambda/\sin\theta$.

```ini
[initial_conditions]
type           = sine-wave
amr_iterations = 0
wave_mode      = shear
wave_angle     = 45.0
wave_amplitude = 1.0e-4
wavelength     = 0.7071067811865476
```

## `couette`: exact compressible Couette flow (verification)

The steady compressible Couette flow with constant $\mu$ and $k$ (Euler with a viscosity and a conductivity, constant
law; fatal otherwise): a fixed isothermal wall at $y_0$, an adiabatic wall at $y_0 + H$ moving at $U$ along x, uniform
pressure $p_0$ (`[initial_conditions_region_1] p`; the other region keys are read and ignored), and, with $s = y - y_0$,

$$
u = U\,\frac{s}{H},\qquad T = T_w + \frac{\mu U^2}{k H}\Big(s - \frac{s^2}{2H}\Big),\qquad \rho = \frac{p_0}{R\,T}.
$$

It is an exact solution of the compressible Navier–Stokes equations, so a run started from it keeps it as its exact
solution at every time, and any departure is the error of the scheme (VV-3, with the walls `wall-isothermal` at
$y_0$ and `wall-noslip` with `wall_u = U` at $y_0 + H$). Keys: `wall_y0`, `height` ($H > 0$), `wall_velocity` ($U$),
`wall_temperature` ($T_w > 0$, the code temperature $p/(\rho R)$).

## `becker-shock`: Becker's exact viscous shock (verification)

Becker's (1922) steady shock along x, in the shock frame: constant viscosity and Prandtl number 3/4 (`[physics]
prandtl = 0.75`, fatal otherwise), for which the total enthalpy $H = c_pT + u^2/2$ is uniform and the momentum integral
reduces to $\tfrac43\mu\,u\,u' = m\,\tfrac{\gamma+1}{2\gamma}(u - u_1)(u - u_2)$, $m = \rho_1 u_1$. Its solution is the
implicit

$$
K\,(x - x_0) = F(u) - F\Big(\frac{u_1 + u_2}{2}\Big),\quad F(u) = \frac{u_1\ln(u_1 - u) - u_2\ln(u - u_2)}{u_1 - u_2},\quad
K = \frac{3m(\gamma+1)}{8\gamma\mu},
$$

solved by bisection in each cell ($F$ is monotone). The upstream state is $\rho_1$, $p_1$ of
`[initial_conditions_region_1]` moving at $u_1 = M a_1$ towards $+x$; $u_2$ is the Rankine–Hugoniot downstream
velocity, $\rho = m/u$, $T = (H - u^2/2)/c_p$. Keys: `mach` ($M > 1$), `shock_x0`. Pair it with an `inflow` of the
upstream state at `bc_x_min` and `extrapolation` at `bc_x_max`. It is the initial state of VV-4: like `couette`, an
exact steady solution, so the error needs no relaxation run and no fit of the shock position.

## `isothermal-atmosphere`, `polytropic-atmosphere`: hydrostatic equilibria under gravity

The hydrostatic equilibrium of the constant gravity $\mathbf g$ of `[gravity]` (required, fatal otherwise), with the
state of `[initial_conditions_region_1]` at the origin ([#73](https://github.com/szaghi/adam/issues/73), M7-G0):

$$
\text{isothermal:}\quad \rho = \rho_1\,e^{\mathbf g\cdot\mathbf x/T_1},\quad p = \rho\,T_1,\quad T_1 = p_1/\rho_1;
\qquad
\text{polytropic:}\quad \rho = \rho_1\,s^{1/(G-1)},\quad p = p_1\,s^{G/(G-1)},\quad
s = 1 + \frac{G-1}{G}\,\frac{\rho_1}{p_1}\,\mathbf g\cdot\mathbf x,
$$

with the region's velocity and field uniform. `polytropic_gamma` $G > 1$; the polytrope is fatal where $s \le 0$ (the
top of the atmosphere inside the domain). Optional, both: a relative pressure pulse $p\,(1 + A e^{-|\mathbf x -
\mathbf x_0|^2/w^2})$ (`pulse_amplitude`, `pulse_x`, `pulse_y`, `pulse_z`, `pulse_width`; M7-G2b). Close them with
`hydrostatic` faces along $\mathbf g$; with `[gravity] equilibrium = atmosphere` the atmosphere at rest is an exact
discrete steady state ([verification](./verification/gravity/well-balanced)).

```ini
[initial_conditions]
type             = polytropic-atmosphere
amr_iterations   = 1
polytropic_gamma = 1.2
```

## `rayleigh-taylor`: two isothermal layers, heavy over light

A heavy cold layer over a light hot one under the gravity of `[gravity]`, which must lie along one axis (M7-G3). With
the height $h = -\mathbf g\cdot\mathbf x/|\mathbf g|$ and $z = (h - h_0)/\delta$, the temperature is

$$
T = a + b\tanh z,\qquad a = \frac{T_l + T_h}{2},\quad b = \frac{T_h - T_l}{2},\quad T_l = \frac{p_1}{\rho_1},\quad
T_h = \frac{T_l}{R},
$$

and the pressure is hydrostatic in closed form, $p = p_1\exp\big(-|\mathbf g|\,\delta\,\frac{a z - b\ln(a\cosh z +
b\sinh z) + b\ln a}{a^2 - b^2}\big)$ ($p(h_0) = p_1$), $\rho = p/T$: two isothermal, hence stably stratified, layers
whose interface is the only instability, with density ratio $R$ across it at equal pressure. The state is at rest, with
the field of region 1 (the light layer), and seeded with a divergence-free single mode,

$$
v_h = A\cos(k x_t)\,f(h),\qquad u_t = -\frac{A}{k}\sin(k x_t)\,f'(h),\qquad f = e^{-((h - h_0)/w)^2},\quad
k = \frac{2\pi}{\lambda},
$$

$v_h$ along $-\mathbf g$ and $u_t$ along $x_t$, the first axis other than $\mathbf g$'s. A velocity-only seed that is
not solenoidal also starts sound waves, which beat against the growing mode for several e-foldings; this one does not.
Keys: `rt_height` $h_0$, `rt_width` $\delta$, `rt_density_ratio` $R$ ($R < 1$: light over heavy, the stable
control), `rt_wavelength` $\lambda$, `rt_amplitude` $A$ (0: the unseeded layers), `rt_seed_width` $w$; all required,
$\delta, R, \lambda, w > 0$. The linear growth rate of the mode is the eigenvalue that
`verification/gravity/rt_oracle.py` computes from the same profile (GV-4,
[verification](./verification/gravity/rayleigh-taylor)); `[gravity] equilibrium = atmosphere` keeps the unseeded layers.

```ini
[initial_conditions]
type             = rayleigh-taylor
amr_iterations   = 1
rt_height        = 0.5
rt_width         = 0.05
rt_density_ratio = 2.0
rt_wavelength    = 0.5
rt_amplitude     = 1.0e-5
rt_seed_width    = 0.1
```

## `polytrope`: an n = 1 polytrope star (self-gravity)

A self-gravitating star in equilibrium ([#81](https://github.com/szaghi/adam/issues/81), M8-P4; Euler, `[gravity]
self_gravity` required). The n = 1 polytrope, $p = K\rho^2$, is the Lane–Emden solution
$\rho_s = \rho_c\,\sin(kr)/(kr)$, $k = \pi/a$, of radius $a$, in equilibrium under its own gravity when
$K = 2Ga^2/\pi$; its mass is $M = 4\rho_ca^3/\pi$. Region 1 sets a light atmosphere ($\rho_1$, $p_1$, the velocity):

- inside, $r < a$: $\rho = \rho_1 + \rho_s$, $p = p_1 + fK\rho_s^2$;
- outside: the atmosphere isothermal at $T = p_1/\rho_1$ and hydrostatic in the star's field,
  $\rho = \rho_1\exp\big(GM/T\,(1/r - 1/a)\big)$, $p = \rho T$. A uniform atmosphere would fall onto the star and
  thin at its surface.

Density and pressure are continuous at the surface. `star_pressure_factor` $f = 1$ is the equilibrium of the star
alone (the atmosphere's own weight and gravity perturb it, by the atmosphere's share of the mass); $f < 1$ collapses
it. Keys: `star_radius` $a$, `star_density` $\rho_c$ (both required, $> 0$), `star_x0`, `star_y0`, `star_z0` (the
centre, default 0), `star_pressure_factor` (default 1, $> 0$). With `[gravity] poisson_bc = isolated` the star sits
in a box without its periodic images (SV-7, SV-8, [verification](./verification/self-gravity/isolated-solver)).
A light atmosphere ($\rho_1 \lesssim 10^{-3}\rho_c$) needs `[numerics] positivity_limiter = cell` at the surface.

```ini
[gravity]
self_gravity = .true.
G            = 1.0
poisson_bc   = isolated

[initial_conditions]
type                 = polytrope
amr_iterations       = 0
star_radius          = 0.5
star_density         = 1.0
star_pressure_factor = 1.0

[initial_conditions_region_1]
r = 1.0e-4
u = 0.0
v = 0.0
w = 0.0
p = 1.5915494309189535e-05
```

## `singular-isothermal-sphere`: Shu's isothermal sphere (self-gravity)

The singular isothermal sphere at rest ([#94](https://github.com/szaghi/adam/issues/94), M12-P4; Euler, `[gravity]
self_gravity` required), $\rho = A a^2/(4\pi G r^2)$, $p = \rho a^2$, with the isothermal sound speed
$a^2 = p_1/\rho_1$ of region 1 (only the ratio is used) and region 1's velocity. $A = 2$ is the hydrostatic sphere of
Shu (1977), whose inside-out collapse accretes at $0.975\,a^3/G$; $A > 2$ collapses everywhere, at $m_0(A)\,a^3/G$
(Shu's A-family; $A = 29.3$ is the run of Federrath et al. 2010). Every cell holds the mean of $1/r^2$ over it (the flux
of $\mathbf x/r^2$ through its faces): sampled at the centres, $1/r^2$ would miss the mass of the unresolved core and
the hydrostatic sphere would expand. Keys: `sis_overdensity` $A$ (default 2, $> 0$), `sis_x0`, `sis_y0`, `sis_z0`
(the centre, default 0). A cell centre on the centre is fatal: place it on a cell corner. The gas stays isothermal
with $\gamma$ close to 1 (`gamma = 1.001`, [SI si11](./verification/self-gravity/sink-particles)), or is isothermal
exactly with [`[eos] type = isothermal`](./input#eos-optional-equation-of-state)
([#100](https://github.com/szaghi/adam/issues/100) P4): then $a^2 = c_0^2$ and the regions need no `p`
([EO P4](./verification/euler/barotropic-gravity)).

```ini
[physics]
gamma = 1.001

[initial_conditions]
type            = singular-isothermal-sphere
amr_iterations  = 0
sis_overdensity = 29.3

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
p = 0.25
```

## `bonnor-ebert`: the isothermal Lane–Emden sphere (self-gravity)

The isothermal sphere in hydrostatic equilibrium ([#100](https://github.com/szaghi/adam/issues/100) P4; Euler with
`[eos] type = isothermal`, `[gravity] self_gravity` required), at rest with region 1's velocity:
$\rho = f\rho_c e^{-\psi(\xi)}$, $\xi = r/r_0$, $r_0 = c_0/\sqrt{4\pi G\rho_c}$, where $\psi$ solves the isothermal
Lane–Emden equation $\psi'' + 2\psi'/\xi = e^{-\psi}$, $\psi(0) = \psi'(0) = 0$ (tabulated by RK4 on steps
$\Delta\xi = 10^{-3}$ from the series start, read by cubic Hermite interpolation). The profile fills the domain: in a
closed box the walls confine it, and the sphere inscribed in the box is a Bonnor–Ebert sphere (Bonnor 1956, Ebert
1955), stable when its $\xi$ is below 6.45. The overdensity $f$ pushes it out of equilibrium (a super-critical sphere
left exactly at equilibrium may as well expand). In a cube the profile is not an exact equilibrium: the gas in the
corners pulls through its cubic multipoles ([EO P4](./verification/euler/barotropic-gravity), eg7). Keys:
`be_density` $\rho_c$ (required, $> 0$), `be_overdensity` $f$ (default 1, $> 0$), `be_x0`, `be_y0`, `be_z0` (the
centre, default 0).

```ini
[eos]
type        = isothermal
sound_speed = 1.0

[initial_conditions]
type           = bonnor-ebert
amr_iterations = 0
be_density     = 0.716    ; xi = 3 at r = 1 with G = 1

[initial_conditions_region_1]
r = 1.0
u = 0.0
v = 0.0
w = 0.0
```

## `boss-bodenheimer`: the rotating cloud with an $m = 2$ perturbation (self-gravity)

The standard test of binary fragmentation (Boss & Bodenheimer 1979, as set by Burkert & Bodenheimer 1993;
[#100](https://github.com/szaghi/adam/issues/100) P4; Euler with `[eos]` isothermal or barotropic, `[gravity]
self_gravity` required): a sphere of radius $R$ centred at $\mathbf x_0$, $\rho = \rho_0(1 + A\cos 2\phi)$ ($\phi$ the
azimuth about $z$) in solid-body rotation $\Omega\,\mathbf e_z \times (\mathbf x - \mathbf x_0)$, embedded in region
1's density and velocity (the light ambient medium). Only the density and the velocity are set: the pressure is the
gas's own $P(\rho)$. Keys: `bb_radius` $R$, `bb_density` $\rho_0$ (both required, $> 0$), `bb_omega` $\Omega$
(required), `bb_amplitude` $A$ (default 0.1), `bb_x0`, `bb_y0`, `bb_z0` (default 0). With $G = 1$, $R = 1$ and
$M = 1$, Burkert & Bodenheimer's ratios $\alpha = E_\mathrm{th}/|W| = 0.26$ and $\beta = E_\mathrm{rot}/|W| = 0.16$
give $c_0 = 0.3225$ and $\Omega = 0.6928$ (`eos/make_eos_gravity.py bb`).
