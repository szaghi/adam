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
realm, as Mignone & Tzeferacos 2010). $\psi$ is **not** part of the energy (mixed GLM); the damping source is applied to
$\psi$ only. $c_h$ (`glm_ch`) is constant and uniform, and it bounds the time step (below). `glm_ch_check` warns or stops
when the fastest wave $|u_n| + c_f$ outruns it.

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
- **Gravity is not free** (planned): an external field brings $g L_0/u_0^2$, self-gravity $G\rho_0 L_0^2/u_0^2$; the
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
- Goedbloed J. P., Poedts S. (2004), *Principles of Magnetohydrodynamics*, Cambridge University Press.
- Derigs D. et al. (2018), Ideal GLM-MHD: about the entropy consistent nine-wave magnetic field divergence diminishing
  ideal MHD equations, *J. Comput. Phys.* 364, 420–467.
- Mignone A., Tzeferacos P. (2010), A second-order unsplit Godunov scheme for cell-centered MHD: the CTU-GLM scheme,
  *J. Comput. Phys.* 229, 2117–2138.
- Stone J. M. et al. (2008), Athena: a new code for astrophysical MHD, *ApJS* 178, 137–177.
- Toro E. F. (2009), *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Springer.
