# Physical models

FLUME solves the conservation laws of compressible flow in the form

$$\frac{\partial \mathbf{q}}{\partial t} + \sum_{d=1}^{3} \frac{\partial \mathbf{F}_d(\mathbf{q})}{\partial x_d} = \mathbf{S}(\mathbf{q}),$$

with $\mathbf{q}$ the conservative state, $\mathbf{F}_d$ the flux along direction $d$ and $\mathbf{S}$ the sources (GLM/EGLM
damping, gravity and self-gravity, cooling, sinks, forcing). This page states the equations and the closures; their
discretisation is on the [numerical methods](./numerics) page and the keys are in the [input reference](./input).

`[physics] physical_model` (`euler`, `euler-de`, `mhd-ideal`), `[physics] passive_scalars`, `[mhd] divergence_control`
and `[eos] type` select one of 13 models. The model fixes the state width `nv` and selects, on the host, the compute
kernels compiled for that width; no kernel branches on the model.

| Model | `[physics] physical_model` | `[mhd] divergence_control` | $\mathbf{q}$ | `nv` |
|---|---|---|---|---|
| Compressible Euler | `euler` | — | $(\rho, \rho u, \rho v, \rho w, E)$ | 5 |
| Euler with [dual energy](#dual-energy) | `euler-de` | — | Euler + $\rho s$ | 6 |
| Euler, `euler-de` with [passive scalars](#passive-scalars) | `euler`, `euler-de` + `passive_scalars` > 0 | — | + 4 slots $\rho s_k$ | 9, 10 |
| Ideal MHD with mixed GLM cleaning | `mhd-ideal` | `glm` | $(\rho, \rho u, \rho v, \rho w, E, B_x, B_y, B_z, \psi)$ | 9 |
| Ideal MHD with EGLM cleaning ($E$ includes $\psi^2/2$) | `mhd-ideal` | `eglm` | as GLM | 9 |
| Ideal MHD without divergence control | `mhd-ideal` | `none` | the first 8 of the above | 8 |
| [Isothermal or barotropic](#isothermal-and-barotropic-gases) Euler | `euler` + `[eos] type` | — | as Euler, energy slot inert | 5 |
| [Isothermal or barotropic](#isothermal-and-barotropic-mhd) ideal MHD | `mhd-ideal` + `[eos] type` | `glm`, `none` | as MHD, energy slot inert | 9, 8 |

The gas is calorically perfect (`[eos] type = ideal`, the default) unless `[eos] type` is `isothermal` or `barotropic`,
which apply to `euler` and `mhd-ideal` without passive scalars only. All models are inviscid by default. Viscosity, heat
conduction and Ohmic resistivity ([#65](https://github.com/szaghi/adam/issues/65)) are each given as a coefficient or as
its dimensionless number ([input](./input#dissipative-terms)). With a viscosity or a conductivity each model is its
compressible Navier–Stokes form, and `mhd-ideal` with a resistivity is resistive MHD (Ohm's law
$\mathbf E = -\mathbf u\times\mathbf B + \eta\mathbf J$; [dissipative terms](#dissipative-terms)). The
subgrid-scale closures of large-eddy simulation are on the [turbulence](./turbulence#subgrid-scale-models) page.

## Thermodynamics

The gas is defined by the specific heats at constant pressure and volume, `[physics] cp, cv` (J/(kg K)), or by their
ratio alone, `[physics] gamma`, with the gas constant then 1 (code units):

$$\gamma = \frac{c_p}{c_v}, \qquad R = c_p - c_v, \qquad p = \rho R T, \qquad a = \sqrt{\frac{\gamma p}{\rho}}.$$

Air, the value of most inputs, is `cp = 1040.004`, `cv = 742.86` ($\gamma = 1.4$). The MHD verification problems use
$\gamma = 5/3$ or $2$ (Brio–Wu). $R$ enters only the temperature: the fluxes, eigensystems and Riemann solvers take
$\gamma$ alone.

## Compressible Euler

$$\frac{\partial}{\partial t}\begin{pmatrix}\rho\\ \rho\mathbf{u}\\ E\end{pmatrix}
+ \nabla\cdot\begin{pmatrix}\rho\mathbf{u}\\ \rho\mathbf{u}\otimes\mathbf{u} + p\,\mathbf{I}\\ (E + p)\,\mathbf{u}\end{pmatrix} = 0,
\qquad E = \frac{p}{\gamma-1} + \frac{1}{2}\rho|\mathbf{u}|^2 .$$

The flux along direction $d$ (normal $\mathbf{n} = \mathbf{e}_d$, normal velocity $u_n = \mathbf{u}\cdot\mathbf{n}$) is
$\mathbf{F}_d = \rho u_n (1, \mathbf{u}, H) + p\,(0, \mathbf{n}, 0)$ with the total specific enthalpy
$H = (E + p)/\rho$.

### Eigensystem

The characteristic decomposition used by the WENO schemes is written once for every direction: the tangents are taken
cyclically, $\mathbf{t}_1 = \mathbf{e}_{\mathrm{mod}(d,3)+1}$, $\mathbf{t}_2 = \mathbf{e}_{\mathrm{mod}(d+1,3)+1}$, so a
problem rotated from $x$ to $y$ or $z$ runs the same arithmetic on permuted data. With $b_2 = (\gamma-1)/a^2$ and
$b_1 = b_2 |\mathbf{u}|^2/2$ the right ($\mathbf{r}_k$) and left ($\mathbf{l}_k$) eigenvectors are

| $k$ | $\lambda_k$ | $\mathbf{r}_k$ | $\mathbf{l}_k$ |
|---|---|---|---|
| 1 | $u_n - a$ | $(1,\ \mathbf{u} - a\mathbf{n},\ H - a u_n)$ | $\tfrac12(b_1 + u_n/a,\ -b_2\mathbf{u} - \mathbf{n}/a,\ b_2)$ |
| 2 | $u_n$ | $(1,\ \mathbf{u},\ \lvert\mathbf{u}\rvert^2/2)$ | $(1 - b_1,\ b_2\mathbf{u},\ -b_2)$ |
| 3 | $u_n$ | $(0,\ \mathbf{t}_1,\ \mathbf{u}\cdot\mathbf{t}_1)$ | $(-\mathbf{u}\cdot\mathbf{t}_1,\ \mathbf{t}_1,\ 0)$ |
| 4 | $u_n$ | $(0,\ \mathbf{t}_2,\ \mathbf{u}\cdot\mathbf{t}_2)$ | $(-\mathbf{u}\cdot\mathbf{t}_2,\ \mathbf{t}_2,\ 0)$ |
| 5 | $u_n + a$ | $(1,\ \mathbf{u} + a\mathbf{n},\ H + a u_n)$ | $\tfrac12(b_1 - u_n/a,\ -b_2\mathbf{u} + \mathbf{n}/a,\ b_2)$ |

evaluated at the Roe average of the two cells adjacent to the face
($\sqrt{\rho}$-weighted velocity and enthalpy, $\tilde a^2 = (\gamma-1)(\tilde H - |\tilde{\mathbf{u}}|^2/2)$).
The unit test `test_flume_euler_library` checks $LR = I$, the homogeneity identity $L\,\mathbf{F}(\mathbf{q}) = \Lambda L\,\mathbf{q}$
and the rotation invariance on 10 000 random states in each direction.


## Dual energy

`physical_model = euler-de` ([#91](https://github.com/szaghi/adam/issues/91)) adds the advected entropy $\rho s$,
$s = p/\rho^\gamma$, to the Euler state ($nv = 6$):

$$\partial_t(\rho s) + \nabla\cdot(\rho s\,\mathbf u) = 0,$$

a pure conservation law with no gravity source. Once per step each cell compares its internal energy
$e = E - |\rho\mathbf u|^2/(2\rho)$ with $\eta E$, $\eta$ = `[physics] dual_energy_eta` (default 0.05; the cell-local
test of Bryan et al. 2014):

- in a cold cell, $e < \eta E$, the entropy is a floor: $e = \max(e, e_s)$, $e_s = \rho s\,\rho^{\gamma-1}/(\gamma-1)$,
  and $E$ is rebuilt from it;
- then, everywhere, $\rho s$ follows $e$.

A cell below its adiabat was drained by the truncation error of the total-energy scheme in a cold supersonic flow, and
the floor lifts it back. A cell above its adiabat was heated by a shock, which the advected entropy does not see, and it
keeps its energy. The floor only adds energy, and the log reports it per step. A reset of the energy from the entropy,
as in Cholla (Villasenor et al. 2021), would remove a cold shock's heat ([DE](./verification/euler/dual-energy)). The
runtime-AMR collapse of [#88](https://github.com/szaghi/adam/issues/88) reaches its end with dual energy
([SV-17](./verification/self-gravity/polytrope)). Discretisation: [numerics](./numerics#dual-energy).

## Passive scalars

`[physics] passive_scalars = N` ([#92](https://github.com/szaghi/adam/issues/92), $0 < N \le 4$) selects the scalar
variant of `euler` or `euler-de`: four more conserved variables $\rho s_k$, advected with the flow,

$$\partial_t(\rho s_k) + \nabla\cdot(\rho s_k\,\mathbf u) = 0 .$$

The first $N$ are set by the `s1`..`sN` keys of the initial-condition regions and of the inflows, the others are zero.
The scalars do not act on the flow: the flow fields of a run with scalars are bitwise those of the same run without
them. Discretisation: [numerics](./numerics#passive-scalars).

## Ideal MHD

In code units the magnetic permeability is absorbed into the field, $\mathbf{B} = \mathbf{B}_{SI}/\sqrt{\mu_0}$, so the
magnetic pressure is $|\mathbf{B}|^2/2$ and the fluid variables keep the units of the input
([units and scaling](#units-and-scaling)):

$$\frac{\partial \rho}{\partial t} + \nabla\cdot(\rho\mathbf{u}) = 0,$$

$$\frac{\partial (\rho\mathbf{u})}{\partial t} + \nabla\cdot\Big[\rho\mathbf{u}\otimes\mathbf{u} + p_T\,\mathbf{I} - \mathbf{B}\otimes\mathbf{B}\Big] = 0,
\qquad p_T = p + \tfrac12|\mathbf{B}|^2,$$

$$\frac{\partial E}{\partial t} + \nabla\cdot\Big[(E + p_T)\,\mathbf{u} - (\mathbf{u}\cdot\mathbf{B})\,\mathbf{B}\Big] = 0,
\qquad E = \frac{p}{\gamma-1} + \tfrac12\rho|\mathbf{u}|^2 + \tfrac12|\mathbf{B}|^2,$$

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla\cdot(\mathbf{u}\otimes\mathbf{B} - \mathbf{B}\otimes\mathbf{u}) = 0,
\qquad \nabla\cdot\mathbf{B} = 0 .$$

The flux of $B_j$ along $\mathbf{n}$ is $u_n B_j - B_n u_j$, so the $B_n$ component of the flux vanishes identically:
without a divergence-control term the normal field is transported only through the transverse fluxes. Every input field
(region states, boundary states) is in code units: a field in Tesla must be divided by $\sqrt{\mu_0}$, a Gaussian one by
$\sqrt{4\pi}$, by the user.

**Wave speeds.** Along $\mathbf{n}$ the system has seven waves: $u_n \pm c_f$ (fast), $u_n \pm c_A$ (Alfvén),
$u_n \pm c_s$ (slow) and $u_n$ (entropy), with $b^2 = |\mathbf{B}|^2/\rho$, $b_n^2 = B_n^2/\rho$, $c_A = |b_n|$ and

$$c_{f,s}^2 = \tfrac12\Big(a^2 + b^2 \pm \sqrt{(a^2 + b^2)^2 - 4a^2 b_n^2}\Big).$$

FLUME evaluates the discriminant as $(a^2 - b^2)^2 + 4a^2 b_t^2$ ($b_t^2 = b^2 - b_n^2$), a sum of non-negative terms that
cannot become negative through round-off when $\mathbf{B}$ is aligned with $\mathbf{n}$.

**Eigensystem.** The 7×7 core is the Roe–Balsara-normalised eigensystem of Stone et al. (2008, ApJS 178, appendix B),
evaluated at the arithmetic average of the two face states (a physical state, so the Roe-matrix factors are $X = 0$,
$Y = 1$), with $B_n$ as a parameter. The degenerate states are handled explicitly: a vanishing transverse field
($|B_t| \le 10^{-12}\max(|\mathbf{B}|, \sqrt\rho\,a)$) and the triple umbilic ($c_f^2 - c_s^2 \le 10^{-12} c_f^2$).
Every three-term sum ($|\mathbf{u}|^2$, $|\mathbf{B}|^2$, $\mathbf{u}\cdot\mathbf{B}$) is summed in increasing order of its terms, so
the result does not depend on the order of the components and rotated problems stay bitwise equal. The unit test
`test_flume_mhd_library` checks the eigensystem, including the degenerate states.

### Divergence control: mixed GLM

The induction equation preserves $\nabla\cdot\mathbf{B} = 0$ analytically but not discretely, and a divergence error
acts as a spurious force parallel to $\mathbf{B}$. With `divergence_control = glm` FLUME uses the mixed hyperbolic–parabolic
cleaning of Dedner et al. (2002): a ninth variable $\psi$ is coupled to the induction equation,

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla\cdot(\mathbf{u}\otimes\mathbf{B} - \mathbf{B}\otimes\mathbf{u}) + \nabla\psi = 0,
\qquad \frac{\partial \psi}{\partial t} + c_h^2\,\nabla\cdot\mathbf{B} = -\frac{c_h^2}{c_p^2}\,\psi ,$$

so the $B_n$ flux becomes $\psi$ and the $\psi$ flux $c_h^2 B_n$: along each direction the pair $(B_n, \psi)$ is a
linear system with speeds $\mp c_h$, decoupled from the seven MHD waves. The divergence error travels away at the
constant cleaning speed $c_h$ and is damped at the rate

$$\frac{c_h^2}{c_p^2} = \alpha\,\frac{c_h}{L},$$

with `glm_alpha` $= \alpha$ and `glm_damping_length` $= L$ (a length, or `min-cell` for the minimum cell spacing of the
realm, as Mignone & Tzeferacos 2010; with runtime regridding the spacing of the finest level allowed, so that the damping
rate does not change with the grid). $\psi$ is **not** part of the energy (mixed GLM); the damping source is applied to
$\psi$ only. $c_h$ (`glm_ch`) is constant and uniform, and it bounds the time step (below). `glm_ch_check` warns or stops
when the fastest wave $|u_n| + c_f$ outruns it. Keep $c_h$ just above that speed. Where the discrete div B is not small,
as on the sheared flanks of magnetised Rayleigh–Taylor, the discretised cleaning amplifies small asymmetries at a rate
that grows with $c_h$. The cleaned equations themselves have no such mode
([verification](./verification/gravity/rayleigh-taylor), [#83](https://github.com/szaghi/adam/issues/83)).

The $(B_n, \psi)$ system of one direction is solved exactly at a face by the Riemann solvers of `weno-riemann`
(Dedner et al. 2002, eq. 42):

$$\tilde B_n = \tfrac12(B_{n,L} + B_{n,R}) - \frac{\psi_R - \psi_L}{2c_h}, \qquad
\tilde\psi = \tfrac12(\psi_L + \psi_R) - \frac{c_h}{2}(B_{n,R} - B_{n,L}).$$

`divergence_control = none` drops $\psi$ ($\mathrm{nv} = 8$). It is a diagnostic variant: in multi-dimensional runs the
divergence error grows unchecked (the magnetised vortex reaches a negative pressure), so production runs use `glm`.

### Divergence control: EGLM

In the mixed GLM $\psi$ changes $B_n$ without appearing in the energy, so the change of magnetic energy is taken from the
thermal pressure; at very low plasma $\beta$ this makes the first-order update of high-order stage states
inadmissible (the Balsara–Spicer blast fails, issue #47). `divergence_control = eglm` selects the ideal GLM-MHD of
Derigs et al. (2018, eqs. 3.16–3.18), with $\psi$ in field units and part of the energy:

$$E = \frac{p}{\gamma - 1} + \frac{\rho|\mathbf{u}|^2}{2} + \frac{|\mathbf{B}|^2}{2} + \frac{\psi^2}{2},$$

the $B_n$ flux $c_h\psi$, the $\psi$ flux $c_h B_n$, the energy flux augmented by $c_h\psi B_n$, and the
nonconservative sources

$$-(\nabla\cdot\mathbf{B})\,(0,\ \mathbf{B},\ \mathbf{u}\cdot\mathbf{B},\ \mathbf{u},\ 0)
- (\mathbf{u}\cdot\nabla\psi)\,(0,\ \mathbf{0},\ \psi,\ \mathbf{0},\ 1)$$

in the order $(\rho, \rho\mathbf{u}, E, \mathbf{B}, \psi)$. $\nabla\cdot\mathbf{B}$ and $\nabla\psi$ are centred
differences of order $2S$, the WENO order (second-order sources drop the Alfvén wave to order 2.26). The damping
$-\alpha (c_h/L)\,\psi$ uses the GLM keys and acts on $\psi$ only, so the removed cleaning energy becomes heat. Density
is conserved exactly; momentum, energy and $\mathbf{B}$ up to the $O(\nabla\cdot\mathbf{B})$ sources. In the flux
splitting the $(B_n, \psi)$ block has eigenvectors $\mathbf{r} = (1, \mp 1)$ at speeds $\mp c_h$ and the energy couples
it to the core ($\mathbf{r}(E) = \mp\psi$, $\mathbf{l}_k(\psi) = -\psi\,\mathbf{l}_k(E)$), so with $\psi = 0$ the arithmetic
of the core is that of GLM. With `weno-riemann` the $(B_n, \psi)$ subsystem is solved exactly at the face in field
units, $\tilde B_n = \tfrac12(B_{nL} + B_{nR}) - \tfrac12(\psi_R - \psi_L)$,
$\tilde\psi = \tfrac12(\psi_L + \psi_R) - \tfrac12(B_{nR} - B_{nL})$; the Riemann solver sees the MHD states (energy
without $\psi^2/2$), and the energy flux gains $c_h\tilde\psi\tilde B_n$ plus the cleaning energy carried by the mass
flux, $F_\rho\,\tilde\psi^2/(2\rho_{\text{up}})$ (a passive scalar, upwind density). Both schemes run EGLM.

### Positivity floors

`[mhd] rho_floor, p_floor` (both $\ge 0$; 0 disables) clip density and pressure after each stage, before the ghost
exchange, and count the floored cells (global counters in the log). With the floors disabled, a non-positive state stops
the run. Independently, every step checks the state for non-finite values and stops the run naming the first one.

## Dissipative terms

With a viscosity or a conductivity in `[physics]` ([input](./input#dissipative-terms)), Euler becomes the compressible
Navier–Stokes system ([#65](https://github.com/szaghi/adam/issues/65)), and so does each MHD model, which also takes an
Ohmic resistivity. Each term is the divergence of a flux, added to the inviscid one in direction $d$:

$$
G_{d}^{\rho u_i} = -\tau_{di},\qquad
G_{d}^{E} = -u_i\,\tau_{di} - k\,\frac{\partial T}{\partial x_d},\qquad
\tau = \mu\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \tfrac{2}{3}(\nabla\cdot\mathbf{u})\,I\right),
$$

with no bulk viscosity (Stokes hypothesis), $T = p/(\rho R)$, and the laws $\mu = \mu_0 (T/T_\mathrm{ref})^\omega$,
$k = k_0 (T/T_\mathrm{ref})^{\omega_k}$ ($\omega_k = \omega$ when $k$ comes from the Prandtl number, otherwise 0; the
constant law is $\omega = 0$), or Sutherland's law
$\mu = \mu_0 (T/T_\mathrm{ref})^{3/2}(T_\mathrm{ref} + S)/(T + S)$.

**Ohmic resistivity (MHD).** With the magnetic diffusivity $\eta$ and $\mathbf J = \nabla\times\mathbf B$ (code units),
Ohm's law $\mathbf E = -\mathbf u\times\mathbf B + \eta\mathbf J$ in $\partial_t\mathbf B = -\nabla\times\mathbf E$ and
the Poynting flux $\mathbf E\times\mathbf B$ give

$$
G_{d}^{B_i} = -\eta\left(\frac{\partial B_i}{\partial x_d} - \frac{\partial B_d}{\partial x_i}\right),\qquad
G_{d}^{E} = \eta\,(\mathbf J\times\mathbf B)_d = \eta\,B_j\left(\frac{\partial B_d}{\partial x_j} - \frac{\partial B_j}{\partial x_d}\right),
$$

with $G_d^E = \sum_i B_i\,G_d^{B_i}$, so the field's energy loss is exactly the Ohmic heating $\eta|\mathbf J|^2$ that the
internal energy receives. The field flux is the **curl form**, $\epsilon_{idk}\,\eta J_k$, not $-\eta\,\partial_d B_i$:
the two differ by $\eta\nabla(\nabla\cdot\mathbf B)$, which is not zero under GLM, and the curl form keeps the resistive
term out of the divergence constraint. The decay rate $\eta k^2$ is measured by [VV-5](./verification/dissipative/).

Discretisation: [numerics](./numerics#dissipative-fluxes-navier–stokes).

## Isothermal and barotropic gases

`[eos]` ([#100](https://github.com/szaghi/adam/issues/100)) replaces the ideal gas by a pressure that depends on the
density alone, $P = P(\rho)$. The energy equation drops out: the state keeps its Euler or MHD layout, and the energy
slot $\rho E$ is inert (zero flux, zero source, eigenvalue 0) and stays at its initial 0. The auxiliary fields report
$p = P(\rho)$, $T = P/\rho$ (with $R = 1$), $a = \sqrt{dP/d\rho}$ and $H = 0$.

**Isothermal** (`type = isothermal`): $P = c_0^2\rho$, the sound speed $c_0$ everywhere. In direction $\mathbf n$ the
waves are $u_n - c_0$, the inert slot (0), $u_n$ twice (the transverse velocities) and $u_n + c_0$. The isothermal
system has one wave between the acoustic ones, the contact, and it carries only the transverse velocity.

The isothermal Riemann problem is solved exactly with no energy equation: the shock jump
$u_* = u_K \mp c_0(\rho_*/\rho_K - 1)\sqrt{\rho_K/\rho_*}$ on the compressive side, the fan
$u_* = u_K \mp c_0\ln(\rho_*/\rho_K)$ on the other. Two streams at Mach $M$ compress to $\rho_* = \rho M_1^2$, $M_1$ the
Mach number in the frame of the shock. `src/tests/flume/verification/eos/eos_oracle.py` builds these solutions
independently of the code, and the runs match them ([EO](./verification/euler/isothermal-eos)).

**Barotropic** (`type = barotropic`): a continuous piecewise polytrope,

$$P = K_k\rho^{\Gamma_k}\ \text{on}\ \rho_{k-1} \le \rho < \rho_k,\qquad K_{k+1} = K_k\rho_k^{\Gamma_k - \Gamma_{k+1}},\qquad
c^2 = \frac{dP}{d\rho} = \Gamma_k\frac{P}{\rho},$$

$c^2$ the speed of the Jeans relation. The default table is that of Mathew & Federrath (2021, MNRAS 507, 2448),
$\Gamma = 1, 1.1, 1.4, 1.1, 5/3$, in cgs. Where $\Gamma$ increases upward (a convex break) a rarefaction through the
break splits into two fans with a plateau at the break density. Where it decreases, $\rho c$ drops and the system is
not genuinely nonlinear there, so Riemann solutions crossing it contain composite waves. The schemes (entropy-satisfying
fluxes) do not depend on convexity. The exact solver of the verification refuses such a law, and the default table,
non-convex at its $1.4 \to 1.1$ break, is verified by self-convergence instead
([EO barotropic](./verification/euler/barotropic-eos)).

**Under gravity** both gases take the momentum sources of the gamma-law gas, $\rho\mathbf g$ (constant or from the cell
potential, self-gravity included), and no energy source. Where the ideal gas reads $c_s^2$ or a thermal energy from the
energy slot, these gases take them from the equation of state, with $c_s^2 = dP/d\rho$:

- the Jeans marker's $\lambda_J = \sqrt{\pi c_s^2/(G\rho)}$;
- the sinks' threshold density $\pi c_s^2/(4 G r_\mathrm{acc}^2)$;
- the sinks' Jeans test $|E_\mathrm{grav}| > 2E_\mathrm{th}$, with the virial thermal energy
  $E_\mathrm{th} = \tfrac{3}{2}\int P\,dV$.

The ideal gas keeps its internal energy in the Jeans test; the two agree at $\gamma = 5/3$. At $\gamma \to 1$ the
internal energy $\int P/(\gamma - 1)\,dV$ diverges (at $\gamma = 1.001$ it is a thousand times the virial term), so the
test could not pass, and the singular-isothermal-sphere runs of the gamma-law gas keep creation off.

Discretisation (eigensystem, Riemann fluxes, positivity): [numerics](./numerics#isothermal-and-barotropic-equations-of-state).

### Isothermal and barotropic MHD

With `physical_model = mhd-ideal` the two gases close ideal MHD without an energy equation, with the total pressure
$P(\rho) + |\mathbf B|^2/2$ (Mignone 2007, J. Comput. Phys. 225, 1427, eqs. 1–6). The models keep the MHD layout
$(\rho, \rho\mathbf u, \rho E, \mathbf B[, \psi])$, the energy slot inert. In direction $\mathbf n$ the six waves
$u_n \mp c_f$, $u_n \mp c_a$, $u_n \mp c_s$ act on $(\rho, \rho u_n, \rho u_{t1}, \rho u_{t2}, B_{t1}, B_{t2})$; there is
no entropy mode. The Jacobian holds the gas law only through $a^2 = dP/d\rho$. The normal field is decoupled at speed 0,
or forms the GLM pair $(B_n, \psi)$ at $\mp c_h$.

- **Floors:** the `[mhd]` density and pressure floors are refused, since there is no pressure to floor.
- **Cleaning:** GLM is the mixed GLM of the ideal MHD. EGLM is refused: its Derigs sources are built on the energy
  equation.
- **Gravity:** it moves the momentum only, as for the Euler gases.

Discretisation: [numerics](./numerics#isothermal-and-barotropic-mhd); verification: [EO MHD](./verification/mhd/barotropic-mhd).

## External gravity

A constant gravitational acceleration $\mathbf g$ (`[gravity]`, [input](./input#gravity-optional-external-gravity-and-self-gravity);
[#73](https://github.com/szaghi/adam/issues/73)) adds sources to every model, with the potential
$\Phi = -\mathbf g\cdot\mathbf x$:

$$
\partial_t(\rho\mathbf u) + \dots = \rho\,\mathbf g,\qquad \partial_t E + \dots = \rho\,\mathbf u\cdot\mathbf g,
\qquad\text{so that}\qquad \partial_t(E + \rho\Phi) + \nabla\cdot\big((E + p_T + \rho\Phi)\,\mathbf u - \dots\big) = 0 .
$$

- **Momentum**: $\rho_i\,\mathbf g$ at the cell centre.
- **Energy**: from the scheme's own mass fluxes $F^\rho$, never from $\rho\mathbf u$ at the centre (Mikami et al.
  2008, Mullen et al. 2021): $\sum_d g_d\,(-F^\rho_{i-3/2} + 13F^\rho_{i-1/2} + 13F^\rho_{i+1/2} - F^\rho_{i+3/2})/24$.
  A face's weights sum to one over the cells it touches, so $\sum(E + \rho\Phi)V$ changes only by boundary fluxes, to
  round-off, with periodic boundaries and across block interfaces. The fluxes of a finite-difference scheme are point
  values of the primitive function of $\rho\mathbf u$, so the two-face average $\tfrac12(F^\rho_{i-1/2} + F^\rho_{i+1/2})$
  is a second-order cell value and takes the whole scheme to second order; the four-face form is fourth order
  ([verification](./verification/gravity/energy-order)). It reads the face beyond each block end, so it takes one more
  ghost layer, `ngc >= S + 1`.
- **Walls**: at a `wall-inviscid` wall along $\mathbf g$ the source reads the mirrored face beyond the wall, and the
  exactly conserved total is $\sum(E + \rho\Phi)V \pm \sum g_n\Delta x_n\rho V/12$ over the cells beside the wall (+ at
  the lower wall, − at the upper), which differs from the plain sum by $O(\Delta x^2)$.
- **Positivity limiter**: the backbone carries the time-centred energy $\rho\mathbf g\cdot(\mathbf u + \mathbf g\Delta t/2)$,
  which leaves its internal energy unchanged, and the high-order excess of the source is limited with the fluxes by the
  cell factor $\Lambda$. Gravity never makes the backbone inadmissible (near vacuum, or with $g\Delta t$ above the sound
  speed), and $E + \rho\Phi$ stays exact wherever the limiter does not act.
- **Time step**: no bound of its own. The first step logs $|\mathbf g|\,\Delta t^2/(\mathrm{CFL}\,h_{\min})$, at least
  the velocity gravity adds in a step over the fastest wave speed; near 1 or above, gravity is stiff.
- **Boundary**: `hydrostatic` faces continue the isothermal equilibrium of the boundary cell into the ghosts.
- **Well-balanced mode** (`[gravity] equilibrium = atmosphere`): the residual of the initial atmosphere at rest is
  computed with the run's own space operator once per grid (at the first step, after a restart, after every regrid) and
  subtracted from every residual (the deviation method of Edelmann et al. 2021, §4.3). The atmosphere is then an exact
  discrete steady state, on `weno` and `weno-riemann`, Euler and MHD, across 2:1 seams (the equilibrium's reflux is
  folded into its residual), and the scheme evolves the deviation from it. Without it a stratified atmosphere at rest
  drives a stationary flow of the truncation error's size ([verification](./verification/gravity/well-balanced)). It
  costs one more field and one residual evaluation per grid. The subtracted residual is outside the limiter's factor:
  it is the truncation error of the equilibrium, small against the state wherever the atmosphere is resolved.
- **Through a cell potential** (`[gravity] potential = linear`, [#81](https://github.com/szaghi/adam/issues/81)):
  the same field as the cell potential $\Phi_i = -\mathbf g\cdot\mathbf x_i$, applied by the path self-gravity uses.
  The face gravity is $-(\Phi_{i+1} - \Phi_i)/\Delta x$, the cell gravity the mean of its two faces; the momentum
  gains $\rho_i$ times the cell gravity, the energy $\tfrac12\sum_d (g_{d,i+1/2}F^\rho_{i+1/2} + g_{d,i-1/2}
  F^\rho_{i-1/2})$. That is the face form, second order, and $\sum(E + \rho\Phi)V$ is conserved without a wall term.
  The hydrostatic ghosts take the potential drop, and the limiter the face and cell gravity. It is a test of that path
  ([SV-1](./verification/self-gravity/potential)).

The equilibria FLUME builds are the `isothermal-atmosphere`, `polytropic-atmosphere` and `rayleigh-taylor` initial
conditions ([initial conditions](./initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity)).

### Self-gravity

`[gravity] self_gravity` ([#81](https://github.com/szaghi/adam/issues/81)) adds the gravity of the gas itself:

$$\nabla^2\Phi = 4\pi G(\rho - \bar\rho), \qquad \partial_t(\rho\mathbf u) \mathrel{+}= -\rho\nabla\Phi, \qquad
\partial_t E \mathrel{+}= -\rho\mathbf u\cdot\nabla\Phi,$$

with the mean density removed on a periodic domain (the Jeans swindle), or $\Phi \to 0$ at infinity with isolated
boundaries (`poisson_bc = isolated`). The discretisation is that of Mullen et al. (2021): $\Phi$ solves
$L_h\Phi = 4\pi G(\rho - \bar\rho)$ (periodic) or $L_h\Phi = 4\pi G\rho$ (isolated) with the 7-point Laplacian (5-point
in 2-D, periodic only), and enters through the potential path above (face gravity, cell gravity, the energy from the
face mass fluxes). A constant $\mathbf g$ may be given too; it adds through its own kernels.

**What is conserved, and why.** The solvers are linear and symmetric, $\Phi = 4\pi G\,K\rho$ with $K = K^T$ ($VK$
symmetric on AMR, below): the periodic one inverts $L_h$ on the zero-mean modes, the isolated one convolves with an
even kernel. Summed over the cells, the energy source is
$\sum_i\Phi_i\,(\nabla\cdot F^\rho)_i V = -\sum\Phi\,\dot\rho\,V$, and
$\dot W = \tfrac12\tfrac{d}{dt}(\rho^T K\rho) = \Phi^T\dot\rho$ because $K$ is symmetric: the semi-discrete scheme
conserves $E + W$, $W = \tfrac12\sum\rho\Phi\,V$, exactly (with periodic boundaries or closed walls, where no mass
flux crosses the boundary). The momentum source $-\sum\rho\,\bar g\,V$ vanishes because the centred difference of an
even, translation-invariant kernel is antisymmetric. Neither needs the discrete Poisson equation to hold at round-off
(the isolated residual is about $10^{-11}$); the symmetry is what counts. That holds exactly for the FFT solvers,
which are fixed linear maps. The multigrid solver on AMR stops at a tolerance, so its $K$ is symmetric, and $E + W$
exact, only to that tolerance: the defect measured in situ is $1.7\cdot10^{-11}$ at a tolerance of $10^{-12}$
([SV-13](./verification/self-gravity/multigrid-amr)), the point Mullen et al. make (round-off conservation needs a
round-off Poisson solution). The Runge–Kutta integration (one solve per
stage, nothing time-averaged) leaves a drift of its own order in $\Delta t$; the positivity limiter's gravity energy
correction (above) breaks the exchange wherever it acts.

The coupling is second order in space whatever the WENO order. A linear Jeans mode obeys
$\omega^2 = c^2k^2 - 4\pi G\rho_0$: a wavelength shorter than the Jeans length $\lambda_J = c\sqrt{\pi/(G\rho_0)}$
oscillates, a longer one grows; an n = 1 polytrope star is an isolated equilibrium ([verification](./verification/self-gravity/)).

**Solvers.** Every Runge–Kutta stage solves $\Phi$ from that stage's density, after the ghost update.
`poisson_solver = auto` takes the FFT solver on a uniform grid and the multigrid solver on an AMR grid.

- **FFT** (FFTW on the CPU, cuFFT on FNL; the `*-fft` build modes): distributed over the z-planes of the transform
  array, bitwise the same on every rank count, with $\Phi$ returned with every ghost layer filled. Periodic: each
  Fourier mode is divided by the eigenvalue of $L_h$, so the discrete residual is round-off. Isolated: the density sits
  in a zero-padded array of $2(N + n_{gc})$ cells per direction and is convolved with the lattice Green's function of
  $L_h$ (Hockney–Eastwood), whose residual is about $10^{-13}$ per offset (`scripts/lattice_green.py`); the residual of
  $\Phi$ is that summed over the source, $2\cdot10^{-11}$ of $4\pi G\rho$ at $64^3$, and the ghosts beyond the domain
  hold the free-space potential. On FNL the density and $\Phi$ cross the host every stage.
- **Multigrid** (`adam_poisson_mg_object`): the composite Laplacian of the AMR leaves, flexible CG preconditioned by a
  V-cycle, its coarsest level solved by the FFT solver of the backend; on FNL it runs on the host. On a uniform grid
  `poisson_solver = multigrid` is a check against `fft`; isolated boundaries use the symmetric James step
  ([#87](https://github.com/szaghi/adam/issues/87)).

**On AMR grids** the composite Laplacian carries a mimetic mean flux across each 2:1 face, symmetric in the
volume-weighted form and consistent ([SV-9](./verification/self-gravity/multigrid-amr)). The face gravity there is first
order, the potential second order. The gravity work across each seam rides on the reflux, so $E + W$ stays exact in the
semi-discrete form; the fully discrete drift is first order in $\Delta t$ with `reflux_cadence = step` and of the
Runge–Kutta order with `stage` ([SV-10](./verification/self-gravity/amr-jeans)). A runtime regrid changes $W$ by a
truncation-sized amount ([SV-11](./verification/self-gravity/runtime-amr)). Isolated boundaries on AMR need every
boundary leaf at the coarsest level. The derivation is on the developer page
[ghost fill and AMR gravity](/developer/flume-ghost-fill#self-gravity-on-amr-grids).

**Limiter.** With `positivity_limiter = cell` the limiter's gravity energy correction acts on the potential path too, so
$E + W$ stays exact only where the limiter does not act. A star in a light atmosphere ($10^{-4}$ of its central density)
needs the limiter at its surface and drifts by $2\cdot10^{-3}$ of $E + W$ over a free-fall time at $32^3$, against
$10^{-9}$ without it.

**Log.** The first step logs $G$, the free-fall time $t_{ff} = \sqrt{3\pi/(32G\rho_{\max})}$ of the densest cell,
$\Delta t/t_{ff}$ (a warning above 0.1) and the wall time of one solve. The conservation history gains `W`
$= \tfrac12\sum\rho\Phi\,dV$ and `E_plus_W`, the energy the scheme conserves.

## Radiative cooling

`[cooling]` ([#93](https://github.com/szaghi/adam/issues/93)) adds optically thin cooling to the hydrodynamic models,

$$\partial_t(\rho e) = -n^2\Lambda(T),\qquad n = \rho/m,\qquad T = p/(\rho R),$$

at fixed $\rho$ and $\rho\mathbf u$, with $E$ lowered by the energy removed. $\Lambda$ is the piecewise power law
through the nodes $(T_k, \Lambda_k)$, $\Lambda_k (T/T_k)^{\alpha_k}$ on $[T_k, T_{k+1}]$; above $T_N$ the last segment
continues. $T_1$ is a floor: below it nothing cools, and there is no heating.

With `euler-de` a cooled cell takes the entropy of its cooled energy, $\rho s = (\gamma - 1)e\rho^{1-\gamma}$:
radiation lowers the entropy for real, and an entropy left on the old adiabat would let the next dual-energy floor lift a
cold cooled cell back to it and undo the cooling.

The cooling itself is integrated exactly, so it sets no limit on the step. `[cooling] dt_factor` (default off) limits
the step to `dt_factor` times the shortest cooling time $t_\mathrm{cool} = \rho e/(n^2\Lambda)$, for runs that need the
coupling with the dynamics resolved. The radiated energy is summed over the ranks; the conservation history adds it back
to $E$ (`E_plus_radiated`, constant to round-off). Discretisation: [numerics](./numerics#radiative-cooling).

## Sink particles

`[sinks]` ([#94](https://github.com/szaghi/adam/issues/94)) adds point masses that act on the gas through self-gravity
only. It needs self-gravity and a hydrodynamic model.

**Gravity.** Each sink's mass is added to the density the Poisson solver sees, so the periodic and isolated boundaries
and the AMR composite operator apply to the sinks as to the gas, and $W = \tfrac12\sum\rho_\mathrm{tot}\Phi\,dV$
includes them. A sink feels the gradient of $\Phi$ at its position. Through the mesh the sink gravity is softened at
about two cells.

**Creation** (`creation = .true.`). Once per step a cell of the finest level present (with runtime regridding, of the
finest level allowed) forms a sink when all of these hold, after Federrath et al. (2010):

- its density exceeds $\rho_\mathrm{thr} = \pi c_s^2/(4 G r_\mathrm{acc}^2)$, with $c_s^2 = \gamma p/\rho$ of the cell:
  a Jeans length of $2 r_\mathrm{acc}$, $r_\mathrm{acc}$ = `accretion_radius_cells` cells;
- no sink lies within $r_\mathrm{acc}$;
- the cell is a minimum of the potential over its control volume, the cells within $r_\mathrm{acc}$ (ties allowed: the 8
  cells around a symmetric centre all qualify);
- the flow converges along every principal axis (the symmetric velocity gradient is negative definite);
- the volume is Jeans-unstable, $|E_\mathrm{grav}| > 2E_\mathrm{th}$, with $E_\mathrm{grav}$ summed over the pairs of its
  cells;
- the volume is bound, $E_\mathrm{grav} + E_\mathrm{th} + E_\mathrm{kin} < 0$, the kinetic energy relative to its centre
  of mass.

The candidates of every rank are taken deepest potential first, one per accretion radius. A static box above the
threshold is not enough: at $\gamma = 1.4$ a control volume at the threshold is Jeans-stable,
$|E_\mathrm{grav}|/2E_\mathrm{th} \approx \pi^2\gamma(\gamma-1)/10 = 0.55$. Gas that compresses adiabatically with
$\gamma > 4/3$ never crosses the threshold ($\rho_\mathrm{thr} \propto c_s^2$ grows with $\rho$ at $\gamma = 2$):
collapse needs [cooling](#radiative-cooling).

**Accretion** (`accretion = .true.`). A new sink has no mass. Every sink then takes, once per step, the gas above the
threshold from the cells within its radius that is bound to it and falls in
($\tfrac12|u - v_s|^2 - Gm_s/r < 0$, $(u - v_s)\cdot(x - x_s) < 0$), plus the cell that holds it and every cell of a new
sink without the checks. A cell within several radii gives its gas to the most bound sink. The cell keeps its specific
values: every conservative variable is scaled by $\rho_\mathrm{thr}/\rho$ (velocity, specific energy, passive scalars;
the dual-energy entropy is rebuilt from the energy). The sink gains the mass and the momentum and moves to the centre of
mass, so accretion conserves mass and momentum exactly.

**Mergers** (`merging = .true.`). Two sinks within an accretion radius that are bound and converging
($\tfrac12\mu|\Delta v|^2 - Gm_1m_2/d < 0$, $\Delta v\cdot\Delta x < 0$) merge at their centre of mass, with their mass,
momentum and mass-weighted acceleration (the pair's mutual pull cancels in it).

**Runtime AMR.** With runtime regridding (`[amr] frequency > 0`) the sinks live on the finest level allowed
(`[amr] max_level`), as in Federrath et al. (2010): they form only there, and $r_\mathrm{acc}$ is in cells of that level
for the whole run. At every regrid a block within $r_\mathrm{acc} + \tfrac12\,$`frequency`$ + 2$ cells of a sink (its
accretion radius, the half cell a sink may move per step until the next regrid, its CIC stencil) is refined towards
that level and never derefined, whatever the other markers vote. A sink formed on a coarser level would have its
neighbourhood refined at the next regrid: the infall onto it prolongated and its threshold quadrupled from one step to
the next. A collapse on runtime AMR needs `regrid_prolongation = conservative-state`
([#88](https://github.com/szaghi/adam/issues/88)).

Discretisation (deposit, time stepping, backends): [numerics](./numerics#sink-particles).

## Units and scaling

The equations above carry no dimensionless number. Take three references, a density $\rho_0$, a length $L_0$ and a
velocity $u_0$, and derive the others:

| Quantity | Reference |
|---|---|
| time | $t_0 = L_0/u_0$ |
| pressure, energy density | $p_0 = \rho_0 u_0^2$ |
| magnetic field ($\mathbf{B}_{SI}/\sqrt{\mu_0}$) | $B_0 = u_0\sqrt{\rho_0}$ |
| GLM speed $c_h$ | $u_0$ |
| $\psi$, mixed GLM | $B_0 u_0$ |
| $\psi$, EGLM | $B_0$ |
| temperature | $u_0^2/R$ ($R^* = 1$) |

Dividing every variable by its reference, every term of the Euler and MHD systems, the GLM and EGLM terms and the GLM
damping $\alpha c_h/L$ included, picks up the same factor, which cancels: the equations in the scaled variables are the
same equations (Goedbloed & Poedts 2004, section 4.1.2). A run in SI, in cgs or in units of order 1 is the same run.
PLUTO, Athena++, MPI-AMRVAC, BATS-R-US, FLASH and RAMSES solve the same $\mu_0$-absorbed form and treat units as a
conversion.

What the invariance does not remove:

- **The parameters move into the data.** The Mach number, the Alfvén Mach number and the plasma $\beta$ live in the
  initial and boundary states; one reference speed can set only one of the flow, sound and Alfvén speeds to 1.
- **The temperature needs $R$**, and later a reference temperature, once conduction, cooling or a temperature-dependent
  resistivity exist.
- **Gravity is not free**: an external field brings $g L_0/u_0^2 = 1/\mathrm{Fr}^2$ ([above](#external-gravity)),
  self-gravity $G\rho_0 L_0^2/u_0^2$ ([#81](https://github.com/szaghi/adam/issues/81)); the dissipative terms bring the
  Reynolds, Prandtl, magnetic Reynolds and Lundquist numbers, each accepted in place of its coefficient
  ([input](./input#dissipative-terms)).
- **The discrete scheme must be scale-free too.** FLUME's is, apart from the classic WENO weights: their absolute
  $\varepsilon$ stops the limiting at small magnitudes ([numerics](./numerics#weno-reconstruction)). With
  `[weno] weights = si` a run whose input is scaled by powers of two (lengths, velocities, density) equals the base run
  bit for bit (the [scaling oracle](./verification/infrastructure/scaling-covariance), NV-4); the admissibility floor of the
  positivity limiter and the immersed-boundary cut spacing are relative.

FLUME uses the invariance in two ways. An input in physical units can carry a
[`[reference]`](./input#reference-optional-dimensional-input) section: the values are divided by their references once,
on the loaded file, and the solver runs on numbers of order 1, whatever the weights; `output_units = dimensional`
writes the results back in the input units. Or the input runs as it is with `weights = si`. On Sod at an interstellar
density both match the code-units run to $10^{-13}$, where the classic weights without the layer stop limiting
(NV-6).

## Time-step bound

The time step is the minimum over the realm of

$$\Delta t = \frac{\mathrm{CFL}}{\max_{\text{cells}} \sum_d \dfrac{|u_d| + c_d}{\Delta x_d}},$$

with $c_d = a$ (Euler) or the fast speed along $d$ (MHD), null directions excluded (`[time] CFL`). With GLM the cleaning
waves add the bound $\Delta t \le \mathrm{CFL} / (c_h \max \sum_d 1/\Delta x_d)$. With dissipative terms each cell
sum gains the diffusive part $2\nu\sum_d 1/\Delta x_d^2$
([numerics](./numerics#dissipative-fluxes-navier–stokes)).

## Derived output fields

With `[IO] save_auxiliary_fields = .true.` the checkpoints also hold $\rho, u, v, w, p, T, H, a$ (plus the specific entropy $s$ with dual energy and the passive scalars $s_1..s_4$); MHD adds the copies `Bx, By, Bz`, the total pressure
`pt` $= p_T$, the plasma beta `beta`, the field magnitude `bmag` and the discrete divergence `divb`. Their exact
definitions are on the [input reference](./input#derived-output-fields) page.

## References

- Bryan G. L. et al. (2014), ENZO: an adaptive mesh refinement code for astrophysics, *ApJS* 211, 19,
  doi:10.1088/0067-0049/211/2/19.
- Dedner A. et al. (2002), Hyperbolic divergence cleaning for the MHD equations, *J. Comput. Phys.* 175, 645–673.
- Derigs D. et al. (2018), Ideal GLM-MHD: about the entropy consistent nine-wave magnetic field divergence diminishing
  ideal MHD equations, *J. Comput. Phys.* 364, 420–467.
- Edelmann P. V. F. et al. (2021), Well-balanced treatment of gravity in astrophysical fluid dynamics simulations at low
  Mach numbers, *Astron. Astrophys.* 652, A53, doi:10.1051/0004-6361/202140653.
- Federrath C., Banerjee R., Clark P. C., Klessen R. S. (2010), Modeling collapse and accretion in turbulent gas
  clouds: implementation and comparison of sink particles in AMR and SPH, *ApJ* 713, 269–290,
  doi:10.1088/0004-637X/713/1/269.
- Goedbloed J. P., Poedts S. (2004), *Principles of Magnetohydrodynamics*, Cambridge University Press.
- Hockney R. W., Eastwood J. W. (1988), *Computer Simulation Using Particles*, Adam Hilger, Bristol.
- Mathew S. S., Federrath C. (2021), *MNRAS* 507, 2448.
- Mignone A. (2007), A simple and accurate Riemann solver for isothermal MHD, *J. Comput. Phys.* 225, 1427–1441.
- Mignone A., Tzeferacos P. (2010), A second-order unsplit Godunov scheme for cell-centered MHD: the CTU-GLM scheme,
  *J. Comput. Phys.* 229, 2117–2138.
- Mikami H., Sato Y., Matsumoto T., Hanawa T. (2008), *Astrophys. J.* 683, 357, doi:10.1086/589759.
- Mullen P. D., Hanawa T., Gammie C. F. (2021), An extension of the Athena++ framework for fully conservative
  self-gravitating hydrodynamics, *Astrophys. J. Suppl.* 252, 30, doi:10.3847/1538-4365/abcfbd.
- Stone J. M. et al. (2008), Athena: a new code for astrophysical MHD, *ApJS* 178, 137–177.
- Toro E. F. (2009), *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer.
- Villasenor B., Robertson B., Madau P., Schneider E. (2021), Effects of photoionization and photoheating on
  Lyman-α forest properties from Cholla cosmological simulations, *ApJ* 912, 138, doi:10.3847/1538-4357/abed5a.
