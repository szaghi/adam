# Physical models

FLUME solves hyperbolic conservation laws of compressible, inviscid flow in the form

$$\frac{\partial \mathbf{q}}{\partial t} + \sum_{d=1}^{3} \frac{\partial \mathbf{F}_d(\mathbf{q})}{\partial x_d} = \mathbf{S}(\mathbf{q}),$$

with $\mathbf{q}$ the conservative state, $\mathbf{F}_d$ the flux along direction $d$ and $\mathbf{S}$ a source, present only
in the GLM damping of the MHD model. Four models are implemented; the model fixes the state width `nv` and selects,
on the host, the compute kernels compiled for that width (no branching on the model inside a kernel).

| Model | `[physics] physical_model` | `[mhd] divergence_control` | $\mathbf{q}$ | `nv` |
|---|---|---|---|---|
| Compressible Euler | `euler` | — | $(\rho, \rho u, \rho v, \rho w, E)$ | 5 |
| Ideal MHD with mixed GLM cleaning | `mhd-ideal` | `glm` | $(\rho, \rho u, \rho v, \rho w, E, B_x, B_y, B_z, \psi)$ | 9 |
| Ideal MHD with EGLM cleaning ($E$ includes $\psi^2/2$) | `mhd-ideal` | `eglm` | as GLM | 9 |
| Ideal MHD without divergence control | `mhd-ideal` | `none` | the first 8 of the above | 8 |

All models use a calorically perfect ideal gas and are inviscid by default. Viscosity, heat conduction and Ohmic
resistivity are milestone M4 ([#65](https://github.com/szaghi/adam/issues/65)), each given as a coefficient or as its
number ([input](./input#dissipative-terms-issue-65-m4)): with a viscosity or a conductivity each model is its
compressible Navier–Stokes form, and `mhd-ideal` with a resistivity is resistive MHD (Ohm's law
$\mathbf E = -\mathbf u\times\mathbf B + \eta\mathbf J$; [numerics](./numerics#dissipative-fluxes-navier-stokes)).

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
([verification](./verification#gv-g3-rayleigh-taylor), [#83](https://github.com/szaghi/adam/issues/83)).

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

## External gravity

A constant gravitational acceleration $\mathbf g$ (`[gravity]`, [input](./input#gravity-optional-external-gravity); milestone M7,
[#73](https://github.com/szaghi/adam/issues/73)) adds sources to both models, with the potential $\Phi = -\mathbf
g\cdot\mathbf x$:

$$
\partial_t(\rho\mathbf u) + \dots = \rho\,\mathbf g,\qquad \partial_t E + \dots = \rho\,\mathbf u\cdot\mathbf g,
\qquad\text{so that}\qquad \partial_t(E + \rho\Phi) + \nabla\cdot\big((E + p_T + \rho\Phi)\,\mathbf u - \dots\big) = 0 .
$$

- **Momentum**: $\rho_i\,\mathbf g$ at the cell centre.
- **Energy**: from the scheme's own mass fluxes $F^\rho$, never from $\rho\mathbf u$ at the centre (Mikami et al.
  2008, Mullen et al. 2021): $\sum_d g_d\,(-F^\rho_{i-3/2} + 13F^\rho_{i-1/2} + 13F^\rho_{i+1/2} - F^\rho_{i+3/2})/24$.
  A face's weights sum to one over the cells it touches, so $\sum(E + \rho\Phi)V$ changes only by boundary fluxes,
  to round-off; the fourth-order weights make the source the cell value of $\rho u\,g$ to fourth order (the face
  form $\tfrac12(F_{i-1/2} + F_{i+1/2})$ is second order on a finite-difference scheme, whose fluxes are point values
  of the primitive function of $\rho u$). At a reflecting wall the invariant is $W$, the sum plus an $O(\Delta x^2)$
  term of the wall rows ([verification](./verification#gv-g2a-the-order-of-the-energy-source)). It takes one more
  ghost layer, `ngc >= S + 1`.
- **Positivity limiter**: the backbone carries the time-centred energy $\rho\mathbf g\cdot(\mathbf u + \mathbf g\Delta t/2)$,
  which leaves its internal energy unchanged, and the high-order excess of the source is limited with the fluxes.
- **Time step**: no bound of its own; the regime $|\mathbf g|\Delta t^2/(\mathrm{CFL}\,h)$ is printed at the first step.
- **Boundary**: `hydrostatic` faces continue the isothermal equilibrium of the boundary cell into the ghosts.
- **Well-balanced mode** (`[gravity] equilibrium = atmosphere`): the residual of the initial atmosphere at rest,
  computed with the run's own operator once per grid, is subtracted from every residual (the deviation method of
  Edelmann et al. 2021, §4.3). The atmosphere is then an exact discrete steady state, and the scheme evolves the
  deviation from it; without it a stratified atmosphere at rest drives a stationary flow of the truncation error's size
  ([verification](./verification#gv-g2b-the-well-balanced-mode)).
- **Through a cell potential** (`[gravity] potential = linear`, M8-P1, [#81](https://github.com/szaghi/adam/issues/81)):
  the same field as the cell potential $\Phi_i = -\mathbf g\cdot\mathbf x_i$, applied by the path self-gravity uses.
  The face gravity is $-(\Phi_{i+1} - \Phi_i)/\Delta x$, the cell gravity the mean of its two faces; the momentum
  gains $\rho_i$ times the cell gravity, the energy $\tfrac12\sum_d (g_{d,i+1/2}F^\rho_{i+1/2} + g_{d,i-1/2}
  F^\rho_{i-1/2})$. That is the face form, second order, and $\sum(E + \rho\Phi)V$ is conserved without a wall term.
  The hydrostatic ghosts take the potential drop, and the limiter the face and cell gravity. A test of that path
  ([SV-1](./verification#sv-1-p1-gravity-from-a-cell-potential)); refused with the well-balanced mode.

The equilibria FLUME builds are the `isothermal-atmosphere`, `polytropic-atmosphere` and `rayleigh-taylor` initial
conditions ([initial conditions](./initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity)).

### Self-gravity

`[gravity] self_gravity` ([#81](https://github.com/szaghi/adam/issues/81), M8-P3) adds the gravity of the gas itself:

$$\nabla^2\Phi = 4\pi G(\rho - \bar\rho), \qquad \partial_t(\rho\mathbf u) \mathrel{+}= -\rho\nabla\Phi, \qquad
\partial_t E \mathrel{+}= -\rho\mathbf u\cdot\nabla\Phi,$$

with the mean density removed on a periodic domain (the Jeans swindle), or $\Phi \to 0$ at infinity with isolated
boundaries (`poisson_bc = isolated`, M8-P4). The discretisation is that of Mullen et al. (2021): $\Phi$ solves
$L_h\Phi = 4\pi G(\rho - \bar\rho)$ (periodic) or $L_h\Phi = 4\pi G\rho$ (isolated) with the 7-point Laplacian (5-point
in 2-D, periodic only), and enters through the potential path above (face gravity, cell gravity, the energy from the
face mass fluxes).

**What is conserved, and why.** The solvers are linear and symmetric, $\Phi = 4\pi G\,K\rho$ with $K = K^T$ ($VK$
symmetric on AMR, below): the periodic one inverts $L_h$ on the zero-mean modes, the isolated one convolves with an
even kernel. Summed over the cells, the energy source is
$\sum_i\Phi_i\,(\nabla\cdot F^\rho)_i V = -\sum\Phi\,\dot\rho\,V$, and
$\dot W = \tfrac12\tfrac{d}{dt}(\rho^T K\rho) = \Phi^T\dot\rho$ because $K$ is symmetric: the semi-discrete scheme
conserves $E + W$, $W = \tfrac12\sum\rho\Phi\,V$, exactly (with periodic boundaries or closed walls, where no mass
flux crosses the boundary). The momentum source $-\sum\rho\,\bar g\,V$ vanishes because the centred difference of an
even, translation-invariant kernel is antisymmetric. Neither needs the discrete Poisson equation to hold at round-off
(the isolated residual is about $10^{-11}$); the symmetry is what counts. The Runge–Kutta integration (one solve per
stage, nothing time-averaged) leaves a drift of its own order in $\Delta t$; the positivity limiter's gravity energy
correction (M7-G1) breaks the exchange wherever it acts.

The coupling is second order in space whatever the WENO order. A linear Jeans mode obeys
$\omega^2 = c^2k^2 - 4\pi G\rho_0$: it oscillates above the Jeans length $\lambda_J = c\sqrt{\pi/(G\rho_0)}$ and grows
below it; an n = 1 polytrope star is an isolated equilibrium ([verification](./verification#self-gravity-m8)).

**On AMR grids** (M8-P5b, periodic only; isolated boundaries on AMR: [#87](https://github.com/szaghi/adam/issues/87)),
$L_h$ is the composite Laplacian of the leaves. Across a 2:1 face it carries the **mimetic mean flux**: every fine
face $f$ of a coarse cell C carries $\tfrac23 w\,(\bar\Phi_f - \Phi_C)$, where $\bar\Phi_f$ is the mean of the face's
fine cells and $w$ is the fine face area over the fine spacing, and C receives the sum of these fluxes. The operator is
symmetric in the volume-weighted form, so $\Phi = 4\pi G\,K\rho$ with $VK$ symmetric, and it is consistent (exact for
linear $\Phi$ face by face). [SV-9](./verification#sv-9-p5a-the-poisson-solver-on-amr-grids) explains why the simpler
couplings fail one property or the other.

Both sides of the seam see the same face gravity, $-(\bar\Phi_f - \Phi_C)/(1.5h)$, which is first order there. The
energy sources, however, need the two-point work. For $E + W$ to telescope, the seam's sources must sum to
$\sum_f A_fF^\rho_f(\Phi_C - \Phi_f)$ (with the sign of the mass flux), and the two kernels give
$\sum_f A_fF^\rho_f(\Phi_C - \bar\Phi_f)$ instead. The difference, $\sum_f A_fF^\rho_f(\bar\Phi_f - \Phi_f)$, is small
($\sum_f(\bar\Phi_f - \Phi_f) = 0$) but not zero.

It rides on the reflux. Every stage adds to the flux register's energy entry:
- the coarse side: its own seam source, taken back;
- the fine side: the exact two-point work less the fine cell's own source.

The end-of-step correction then gives the coarse cell exactly the missing work. The semi-discrete $E + W$ is exact
again, the fine cells keep the consistent gravity, and the coarse cell next to the seam takes an $O(h)$ local error.

The fully discrete exchange is first order in $\Delta t$, not of the Runge–Kutta order. The reflux corrects the
coarse cells only at the end of the step, so its stages advance them with their own seam fluxes, an $O(\Delta t)$
inconsistency per step. Measured on the collapsing Jeans mode
([SV-10](./verification#sv-10-p5b-self-gravity-on-a-static-amr-grid)):
- with the seam work: $9\cdot10^{-8}$, $3\cdot10^{-8}$, $1.4\cdot10^{-8}$ at CFL 0.4, 0.2, 0.1;
- without it: $10^{-5}$ at every CFL, a spatial leak.

`[numerics] reflux_cadence = stage` (P5c) restores the Runge–Kutta order: every stage corrects the coarse cells'
residual from that stage's fluxes, and the drift falls at third order in $\Delta t$ ($3.5\cdot10^{-8}$,
$4.3\cdot10^{-9}$, $5.4\cdot10^{-10}$).

**Under runtime AMR** (P5c) every regrid sets the solver up on the new grid and solves $\Phi$ from the regridded
state. The regrid conserves $\rho$, $\rho\mathbf u$ and $E$, but $W = \tfrac12\sum\rho\Phi V$ is the grid's own
discrete energy: a regrid that changes the grid changes $W$ by a truncation-sized amount, and $E + W$ jumps
([SV-11](./verification#sv-11-p5c-self-gravity-under-runtime-amr)). The jump is not a leak of the scheme between
regrids, which stays exact. Making the regrid itself conserve $E + W$ would need a correction of $E$ at each regrid,
which is not done.

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
- **Gravity is not free**: an external field brings $g L_0/u_0^2 = 1/\mathrm{Fr}^2$ (`[gravity]`, milestone M7,
  [#73](https://github.com/szaghi/adam/issues/73), [above](#external-gravity)), self-gravity
  $G\rho_0 L_0^2/u_0^2$ (M8, [#81](https://github.com/szaghi/adam/issues/81)); the
  dissipative terms of M4 bring the Reynolds, Prandtl, magnetic Reynolds and Lundquist numbers, each accepted in
  place of its coefficient ([input](./input#dissipative-terms-issue-65-m4)).
- **The discrete scheme must be scale-free too.** FLUME's is, apart from the classic WENO weights: their absolute
  $\varepsilon$ stops the limiting at small magnitudes ([numerics](./numerics#weno-reconstruction)). With
  `[weno] weights = si` a run whose input is scaled by powers of two (lengths, velocities, density) equals the base run
  bit for bit (the [scaling oracle](./verification#scaling-covariance), NV-4); the admissibility floor of the
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
([numerics](./numerics#dissipative-fluxes-navier-stokes)).

## Derived output fields

With `[IO] save_auxiliary_fields = .true.` the checkpoints also hold $u, v, w, p, H, a$; MHD adds the total pressure
`pt` $= p_T$, the plasma beta `beta`, the field magnitude `bmag` and the discrete divergence `divb`. Their exact
definitions are on the [input reference](./input#derived-output-fields) page.

## References

- Dedner A. et al. (2002), Hyperbolic divergence cleaning for the MHD equations, *J. Comput. Phys.* 175, 645–673.
- Edelmann P. V. F. et al. (2021), Well-balanced treatment of gravity in astrophysical fluid dynamics simulations at low
  Mach numbers, *Astron. Astrophys.* 652, A53, doi:10.1051/0004-6361/202140653.
- Goedbloed J. P., Poedts S. (2004), *Principles of Magnetohydrodynamics*, Cambridge University Press.
- Derigs D. et al. (2018), Ideal GLM-MHD: about the entropy consistent nine-wave magnetic field divergence diminishing
  ideal MHD equations, *J. Comput. Phys.* 364, 420–467.
- Mignone A., Tzeferacos P. (2010), A second-order unsplit Godunov scheme for cell-centered MHD: the CTU-GLM scheme,
  *J. Comput. Phys.* 229, 2117–2138.
- Mikami H., Sato Y., Matsumoto T., Hanawa T. (2008), *Astrophys. J.* 683, 357, doi:10.1086/589759.
- Mullen P. D., Hanawa T., Gammie C. F. (2021), An extension of the Athena++ framework for fully conservative
  self-gravitating hydrodynamics, *Astrophys. J. Suppl.* 252, 30, doi:10.3847/1538-4365/abcfbd.
- Stone J. M. et al. (2008), Athena: a new code for astrophysical MHD, *ApJS* 178, 137–177.
- Toro E. F. (2009), *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer.
