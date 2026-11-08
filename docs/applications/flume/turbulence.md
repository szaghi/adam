# Turbulence

FLUME runs large-eddy simulations of compressible wall-bounded and free turbulence
([#119](https://github.com/szaghi/adam/issues/119)) with these pieces:

| Piece | Selected by | Section |
|---|---|---|
| Algebraic subgrid-scale models (Vreman, WALE, Sigma) | `[turbulence] sgs_model` | [Subgrid-scale models](#subgrid-scale-models) |
| Kinetic-energy-preserving central flux away from shocks | `[numerics] scheme_space = hybrid-kep-weno` | [Hybrid central/WENO](#hybrid-central-weno) |
| Impermeable no-slip walls | `wall-noslip`, `wall-isothermal` faces | [Closed inviscid fluxes](#no-slip-walls-closed-inviscid-fluxes) |
| Equilibrium and inverse-transformation wall models | `wall-modelled` faces | [Wall model](#wall-model) |
| Turbulent boundary-layer inflow | `inflow-digital-filter` faces | [Digital-filter inflow](#digital-filter-inflow) |
| Constant mass flow in a periodic channel | `[forcing] type = mass-flow` | [Mass-flow forcing](#mass-flow-forcing) |
| Kinetic-energy budget and wall histories | `[diagnostics]` | [Diagnostics](#turbulence-diagnostics) |

All of them are Euler only (`euler`, and `euler-de` where stated); MHD is refused. The keys are in the
[input reference](./input); the results are on the [turbulence verification page](./verification/turbulence/).

## Subgrid-scale models

`[turbulence] sgs_model` ([input](./input#turbulence-optional-algebraic-subgrid-scale-model)) adds an eddy viscosity
$\mu_t = \rho\nu_t$ to the [dissipative fluxes](./numerics#dissipative-fluxes-navier–stokes). $\mathbf{g}$ is the velocity
gradient and $\Delta$ the geometric mean of the active steps of the block:

| Model | $\nu_t$ | Vanishes in |
|---|---|---|
| Vreman (2004) | $c\sqrt{B_\beta/(\alpha_{ij}\alpha_{ij})}$, $\alpha_{ij} = \partial_i u_j$, $\beta_{ij} = \sum_m \Delta_m^2 \alpha_{mi}\alpha_{mj}$ | $\alpha_{ij}$ of rank 1, pure shear among them |
| WALE (Nicoud & Ducros 1999) | $(c\Delta)^2 (S^d_{ij}S^d_{ij})^{3/2} / \left((S_{ij}S_{ij})^{5/2} + (S^d_{ij}S^d_{ij})^{5/4}\right)$ | pure shear |
| Sigma (Nicoud et al. 2011) | $(c\Delta)^2 \varsigma_3(\varsigma_1-\varsigma_2)(\varsigma_2-\varsigma_3)/\varsigma_1^2$, $\varsigma_i$ the singular values of $\mathbf{g}$ | pure shear, 2-D flow, solid rotation |

The model constants are taken from other discretisations and are not calibrated on WENO.

**Where it enters.** $\nu_t$ is computed once per stage from the auxiliary variables with 2nd-order central gradients,
on the interior cells and 2 ghost layers, since the fluxes read it up to the second ghost. The dissipative kernels then
add $\tau$ and the heat flux with $\mu = \mu_t$ and $k = \mu_t c_p/Pr_t$ on top of the molecular ones; the flux is linear
in both. The face value is the mean of the two cells, which cannot go negative where a 4-point interpolant could; the
4th-order correction takes the cell value. The isotropic part of the subgrid stress is dropped. The reflux and the
inter-realm register see the total flux, so conservation is unchanged.

**Time step.** The largest $2\max(4/3, \gamma/Pr_t)\,\nu_t\sum_d 1/\Delta x_d^2$ of the domain is added to the limit of
every cell. It is an over-estimate, independent of the partition, taken from the committed state of the previous step;
the growth of $\nu_t$ during a step is not bounded.

**Limits.** $\Delta$ is per block, so $\mu_t$ jumps at a 2:1 face. The first time step of a multi-realm run takes its
limit from seam ghosts that are not filled yet. A model counts as a dissipative term: the positivity limiter and
immersed solids are refused with it.

## Hybrid central/WENO

`scheme_space = hybrid-kep-weno` ([input](./input#numerics)) computes the WENO flux splitting as `weno` does, then
replaces the flux with a non-dissipative one at the faces away from shocks.

**Sensor.** Ducros et al. (1999), per cell from 2nd-order central velocity gradients:
$\Psi = \theta^2 / (\theta^2 + |\boldsymbol\omega|^2 + \varepsilon)$, $\varepsilon = 10^{-8}(a/\Delta x_\min)^2$
(`hybrid_floor`, not calibrated). A cell is also flagged where the relative density jump across it exceeds
`hybrid_jump` (0.1, not calibrated), because the Ducros ratio is zero at a contact. A face takes the central flux when
$\Psi <$ `hybrid_threshold` (0.2) in the 4 cells $-1\ldots2$ along its direction.

**Central flux.** 6th-order, kinetic-energy- and pressure-equilibrium-preserving, in the conservative form of Pirozzoli
(2011) with the two-point flux of Shima et al. (2021), $\langle a\rangle = (a_l+a_r)/2$:

$$
F_{i+1/2} = 2\sum_{l=1}^{3} a_l \sum_{m=0}^{l-1} f(i-m,\,i-m+l),\quad a = (3/4,\,-3/20,\,1/60),
$$
$$
f^\rho = \langle\rho\rangle\langle u_n\rangle,\quad f^{\rho u_i} = \langle\rho\rangle\langle u_n\rangle\langle u_i\rangle + \langle p\rangle\delta_{ni},\quad f^{E} = \langle\rho\rangle\langle u_n\rangle\frac{\mathbf{u}_l\cdot\mathbf{u}_r}{2} + \frac{\langle p\rangle\langle u_n\rangle}{\gamma-1} + \frac{u_{n,l}p_r+u_{n,r}p_l}{2}.
$$

Each face holds one flux, so conservation holds whichever kernel set it. WENO is computed everywhere, so the scheme
costs what `weno` costs plus the overwrite; the two kernel sets are not fused.

**Seams and boundaries.** The central stencil must not read ghost cells that are not a same-level copy. A face along
$d$ whose stencil reaches the ghost cells of a 2:1 face, a physical face or a forest seam normal to $d$ keeps the WENO
flux, so a seam and its reflux stay what they are under `weno`. A flag acts along its own normal only. On such grids
the sensor is exchanged like a state, one more halo exchange per stage: otherwise the two blocks of a same-level
interface can take different fluxes at their common face and lose conservation ($2\cdot10^{-7}$ measured on a 2:1
box).

`hybrid_skin` relaxes this rule for measurement ([input](./input#numerics)). The default stays WENO: on a smooth flow
the central skin is closer to the DNS, but with a shock at a seam there is no reference
([verification](./verification/turbulence/hybrid-scheme#against-the-dns)).

**With a subgrid-scale model.** The eddy viscosity is zero on flagged cells. The rule is per cell: a WENO face next to
a flagged cell, or on a skin, can still carry the model. The sensor and its exchange then run twice per stage.

## No-slip walls: closed inviscid fluxes

At a `wall-noslip` or `wall-isothermal` face of a hydrodynamic model every inviscid face flux except the normal
momentum is set to zero before the dissipative fluxes are added, so the wall is impermeable. The mirrored ghost of an
isothermal wall has another density, and without the closure the split flux carries mass, energy and scalars through
the face at $O(h^3)$: +11 % of the mass over a coarse turbulent channel run. MHD walls keep their fluxes.

## Wall model

A `wall-modelled` face ([input](./input#bc-x-min-bc-x-max-bc-y-min-bc-y-max-bc-z-min-bc-z-max)) mirrors its ghosts as
`wall-inviscid` and replaces the flux through the wall face, after every other dissipative flux:

- no mass flux;
- the parallel momentum flux is $\tau_w$, along the parallel velocity of the sample;
- the energy flux is the wall heat flux $q_w$;
- the normal momentum keeps the mirrored pressure.

The sample is the 3rd cell off the wall, $h = 2.5\,\Delta$, or the height `wall_model_height`, interpolated linearly
between the two cells around it. Kawai & Larsson (2012) place the sample off the first cell. Equilibrium wall models
are usually matched at 0.1–0.3 of $\delta$ with $\Delta/h_{\mathrm{wm}}\lesssim0.3$–0.4: the 3rd cell sits at that bound
and moves with the grid.

### Equilibrium model

The model is the equilibrium thin layer between the wall and the sample, with uniform pressure across it:

$$
\frac{d}{dy}\Big[(\mu+\mu_t)\frac{dU}{dy}\Big]=0,\qquad
\frac{d}{dy}\Big[(\mu+\mu_t)U\frac{dU}{dy}+c_p\Big(\frac{\mu}{Pr}+\frac{\mu_t}{Pr_t}\Big)\frac{dT}{dy}\Big]=0,
$$

with $\mu_t=\kappa y\sqrt{\rho\tau_w}\,\big(1-e^{-y^+/A^+}\big)^2$ and the local density. The damping coordinate is
$y^+=y\sqrt{\tau_w\rho_w}/\mu_w$ (`classical`) or $y^*=y\sqrt{\tau_w\rho}/\mu$ with the local density and viscosity
(`semi-local`, as Yang & Lv 2018 and Chen et al. 2022). $\kappa=0.41$, $A^+=17$, $Pr_t=0.9$; $T=T_w$ and $U=0$ at the
wall, the sample values at $y=h$.

Integrated once, $\tau_w$ and $q_w$ follow from two quadratures for a given temperature profile:

- 40 intervals stretched by 1.3, a rule exact for a coefficient linear on each interval;
- 8 iterations, with an Aitken step on $\ln\tau_w$ every second one (capped at a factor $e$);
- the temperature iterate under-relaxed and kept positive, and the laminar stress if all else fails.

The cost per wall face is fixed, on the device too. Measured by the unit test `test_wall_model`: $\tau_w$ within 0.6 % of
a fine quadrature for $h^+=10^2$–$10^5$ (5.5 % at $10^6$); the Crocco relation to $8\cdot10^{-4}$; a hot compressible
layer against an independent solver to $3\cdot10^{-3}$; no non-finite result in $2\cdot10^5$ random samples up to 4.5
wall sound speeds.

### Inverse-velocity-transformation model

`wall_model = inverse-transformation` follows Griffin, Fu & Moin (2023, J. Fluid Mech. 970, A36). It keeps the
constant-property layer above, at the same semi-local height, $S_t^+=1/(1+\kappa y^*D^2)$, $D=1-e^{-y^*/A^+}$, and maps
its strain to the compressible one by inverting the total-stress-based transformation (their Eq. 2.2):

$$
\frac{dU^+}{dy^*}=\left[\frac1{\mu^+S_t^+}-\frac1{\mu^+}+\sqrt{\rho^+}\left(1+\frac{y^+}{2\rho^+}\frac{d\rho^+}{dy^+}-\frac{y^+}{\mu^+}\frac{d\mu^+}{dy^+}\right)\right]^{-1}.
$$

The temperature comes from the Reynolds analogy anchored at the sample (their Eq. 2.5),
$T=T_w+s\,Pr\,(T_r-T_w)\frac{U}{U_e}\big(1-\frac{U}{U_s}\big)+\big(\frac{U}{U_s}\big)^2(T_s-T_w)$,
$T_r=T_e+Pr^{1/3}U_e^2/(2c_p)$, $s=1.14$, so $q_w=s\,c_p\,\tau_w(T_r-T_w)/U_e$. The edge state $U_e$, $T_e$ is an input
(`wall_model_edge_velocity`, `wall_model_edge_temperature`): the paper locates the edge in the flow, which a face kernel
cannot do.

The property gradients follow from $T(U)$ (uniform pressure, power-law viscosity), so for the physical strain
$g=dU/dy$ Eq. 2.2 is the quadratic $F\beta y\,g^2-(E+F+K\beta y)\,g+K=0$, with $E=\kappa y^*D^2\mu_w/\mu$,
$F=\sqrt{\rho/\rho_w}$, $K=\tau_wF/\mu$, $\beta=T'(U)(\tfrac12+\omega)/T$. FLUME takes its positive root when $\beta<0$
(the other is negative), else the root nearer $K/(E+F)$. That value is the strain without the gradients, the exact
laminar $\tau_w/\mu$; the equation also holds the spurious root $1/(\beta y)$ of $dy^*/dy=0$. The layer is integrated
from the wall on the same 40 stretched intervals, each exact for $1/g$ linear in $y$ and corrected twice, and shot on
$\ln\tau_w$ (secant, 8 shots).

Measured by `test_wall_model`:

- in the constant-property limit it equals the equilibrium model to $2\cdot10^{-9}$;
- against Eq. 2.2 integrated as printed on $2\cdot10^5$ intervals, $3.9\cdot10^{-3}$ on a cold wall and
  $7.8\cdot10^{-3}$ on a hot one, falling to $3\cdot10^{-4}$ and $5\cdot10^{-5}$ with 160 and 640 intervals (the error
  is the grid's, not the algebra's);
- no non-finite result in $2\cdot10^5$ random samples.

On a laminar channel ($h^+\approx1$) it equals the equilibrium model, and on a hot wall it departs from it.

### Input filter

`wall_model_filter` $=T_f$ averages the sampled velocity, pressure and density of each wall face in time,
$\bar a\leftarrow(1-e)\,\bar a+e\,a$, $e=\min(\Delta t/T_f,1)$, as Kawai & Larsson (2012) recommend against the noise
of a single LES sample.

- The average takes the sample once per step, at the first stage (from $q^n$). Every stage uses $(1-e)\bar a+e\,a$
  with its own sample, so $T_f\to0$ is the unfiltered model bit for bit.
- The residual of a balanced gravity equilibrium uses the unfiltered model and leaves the averages alone.
- The averages restart from the next sample at the start and after a regrid.
- A checkpoint carries them (`<restart_basename>-proc<rank>.wmf`, each block's averages after its Morton code, checked
  on load), so a restart continues the run bit for bit. A checkpoint without the file restarts them with a warning.

Near the steady state of a laminar channel $T_f=0.05$ moves the force by $3\cdot10^{-4}$, the lag $T_f\,|df/dt|/f$;
over an early transient $T_f=1$ moves it by $1.4\cdot10^{-2}$.

### Eddy viscosity below the matching height

The face takes $\tau_w$ and $q_w$ from the sample, and nothing in the cells between the wall and $h$ carries them there
unless a subgrid model does. At Mach 3 the first cell then drains: colder than the wall, dense, its flow reversed
([channel at Mach 3](./verification/turbulence/channel#mach-3)).

With `wall_model_eddy_viscosity = .true.` the cells whose centre lies below $h$ take the model's own eddy viscosity,
$\mu_t=\kappa y\sqrt{\rho\tau_w}\,D^2$, $D = 1-e^{-y^+/A^+}$ with the damping of the face (semi-local for the inverse
model). $y$ is the cell centre, $\rho$ and $\mu$ are the cell's, and $\tau_w$ comes from the instantaneous sample, never
the filtered one. The cell takes the larger of it and a subgrid model's, after the hybrid mask. It goes through the
subgrid flux kernels ($k_t=\mu_t c_p/Pr_t$) and the subgrid time-step limit. In a log layer $\mu_t\,dU/dy=\tau_w$, so
the resolved layer carries the modelled stress.

The columns of the two ghost layers along the wall are computed from their ghost cells, so the two blocks of a face
agree on its flux: the same cells split into 2 or 4 blocks per direction agree to round-off, against $2\cdot10^{-6}$
with those columns left out. Without the key nothing changes, bit for bit.

### Limits

- Isothermal walls only.
- The inverse model takes a fixed edge state.
- A restart needs the same ranks, as every FLUME checkpoint.
- The wall faces are rebuilt each step, so after a regrid too, but no test regrids along a wall.

## Digital-filter inflow

An `inflow-digital-filter` face ([input](./input#bc-x-min-bc-x-max-bc-y-min-bc-y-max-bc-z-min-bc-z-max)) feeds a
spatially developing turbulent boundary layer: the mean of a profile table plus fluctuations with its Reynolds
stresses. The method is the digital filter of Klein, Sadiki & Janicka (2003) in the form of Touber & Sandham (2009) and
Xie & Castro (2008):

1. **Plane.** A uniform grid over the face, `df_spacing` apart (adjusted to divide the domain), periodic along a
   periodic span.
2. **Filter.** Uniform random numbers of zero mean and unit variance are filtered along both plane axes with
   $b_k=\exp(-\pi|k|/n)$, normalised to unit sum of squares, $n=L/\Delta$ the integral length in nodes, half-width
   $2n$. The filtered field has unit variance and integral length $2L/\pi$, that of $\exp(-\pi r/2L)$; its correlation
   is $(1+\pi r/n)e^{-\pi r/n}$.
3. **Time.** The streamwise correlation is a temporal one,
   $\psi \leftarrow \psi\,e^{-\pi\Delta t/2T} + \psi_\mathrm{new}\sqrt{1-e^{-\pi\Delta t/T}}$, $T=L_x/U_\infty$, once
   per step. Unit variance is kept, and the correlation is $e^{-\pi t/2T}$.
4. **Stresses.** The signals become velocity fluctuations by the Cholesky factor of the table's stress tensor (Lund,
   Wu & Squires 1998), with $R_{13}=R_{23}=0$.
5. **Thermodynamics.** The strong Reynolds analogy, $T'/\bar T=-(\gamma-1)\bar M^2u'/\bar U$, written
   $T'=-(\gamma-1)\bar U u'\bar\rho/(\gamma p_\infty)$ so that it holds at the wall. The pressure is the free stream's,
   and $\rho'/\bar\rho=-T'/\bar T$.

The plane lives on the host, and every rank builds the same one from one random stream and the forest's global step.
Once per step each rank lists the ghost rows beyond the face and interpolates their states bilinearly in the plane. The
CPU writes them after the face pass of every boundary fill; on FNL they are copied to the device and written by a small
kernel. In the edge and corner passes the face is extrapolated. The plane is saved with the restart
(`<restart_basename>.dff`), so a restart continues bit for bit; without the file a fresh field starts, with a warning.

The `boundary-layer` initial condition ([initial conditions](./initial-conditions#boundary-layer-mean-profile-of-a-boundary-layer))
fills the domain with the mean of the same table.

**Limits.**

- Euler only, one face per realm, the wall on the minimum side of its axis.
- With a periodic span, the ghosts beyond the face and the span are filled by the periodic exchange, so they lag the
  plane by one boundary fill.
- A digital filter needs a recovery length after the inlet, about $40\,\delta_i$ at Mach 2.5 against $20\,\delta_i$
  for recycling (Huang, Duan & Choudhari 2022).

## Mass-flow forcing

`[forcing] type = mass-flow` ([input](./input#forcing-optional-mass-flow-forcing-of-a-periodic-channel)) drives a
periodic channel at a constant mass flow. Once per step a uniform momentum is added along `axis` so that the volume
mean of the momentum is `bulk_momentum`; the total energy gains the change of kinetic energy. Euler only.

The force per unit volume (the added momentum over the step) is written to `<output_basename>-forcing_history.dat`
(header `VARIABLES="it" "time" "f"`, then one row per save) on the conservation-history cadence, in the output units of
`[reference]`. At a statistically steady state its mean times the half height is the wall shear stress. A restart
appends to the file.

## Turbulence diagnostics

### Kinetic-energy budget

`[diagnostics] turbulence_history_save = N` ([input](./input#diagnostics)) writes every $N$ steps the volume integrals
of the kinetic-energy budget of a compressible flow to `<output_basename>-turbulence_history.dat`, to measure how much
the scheme dissipates. It does not change the solution; 0, the default, computes nothing.

$$
K = \int \tfrac12 \rho\,|\mathbf{u}|^2\,dV,\quad
Z = \int \tfrac12 \rho\,|\boldsymbol{\omega}|^2\,dV,\quad
\varepsilon_s = \int \mu\,|\boldsymbol{\omega}|^2\,dV,\quad
\varepsilon_d = \int \tfrac43 \mu\,\theta^2\,dV,\quad
\varepsilon = \int \tau_{ij}\,\frac{\partial u_i}{\partial x_j}\,dV,\quad
p_d = \int p\,\theta\,dV,
$$

with $\theta = \nabla\cdot\mathbf{u}$ and $\tau$, $\mu$ those of the dissipative fluxes. The gradients are 4th-order
central differences of the auxiliary variables, whatever `dissipative_order` is. The integrals are sums, not divided by
the volume. With a subgrid model the file gains a last column, `eps_sgs`.

The exact equation is $dK/dt = p_d - \varepsilon$. What is left of the measured drop of $K$ after the integral of
$\varepsilon - p_d$ is the numerical dissipation of the scheme, plus the truncation error of the integrals. On a periodic
box with constant $\mu$, $\varepsilon = \varepsilon_s + \varepsilon_d$ only in the limit of exact integration.

Each row refreshes the ghosts and the auxiliary variables of the committed state. On a forest of several realms there is
no step-0 row, because the seams are not built yet.

### Wall histories

`[diagnostics] wall_history_save = N` writes `<output_basename>-wall_<face>_history.dat`, one file per `wall-modelled`
face: `it`, `time`, then the wall shear stress along `wall_history_axis` (positive for the drag of a flow along +axis)
and the wall heat flux ($k\,dT/dy$ at the wall, positive when the gas is hotter), each area-averaged over the face cells
in `wall_history_bins` equal bins along the axis. The model is evaluated on the committed state at its sample with the
solvers of the face fluxes, so the history is the stress the solver applies (with `wall_model_filter`, the unfiltered
one).

## References

- Ducros F. et al. (1999), Large-eddy simulation of the shock/turbulence interaction, *J. Comput. Phys.* 152, 517–549.
- Griffin K. P., Fu L., Moin P. (2023), *J. Fluid Mech.* 970, A36.
- Kawai S., Larsson J. (2012), Wall-modeling in large eddy simulation: length scales, grid resolution, and accuracy,
  *Phys. Fluids* 24, 015105.
- Klein M., Sadiki A., Janicka J. (2003), A digital filter based generation of inflow data for spatially developing
  direct numerical or large eddy simulations, *J. Comput. Phys.* 186, 652–665.
- Lund T. S., Wu X., Squires K. D. (1998), Generation of turbulent inflow data for spatially-developing boundary layer
  simulations, *J. Comput. Phys.* 140, 233–258.
- Nicoud F., Ducros F. (1999), Subgrid-scale stress modelling based on the square of the velocity gradient tensor,
  *Flow Turbul. Combust.* 62, 183–200.
- Nicoud F., Baya Toda H., Cabrit O., Bose S., Lee J. (2011), Using singular values to build a subgrid-scale model for
  large eddy simulations, *Phys. Fluids* 23, 085106.
- Shima N., Kuya Y., Tamaki Y., Kawai S. (2021), Preventing spurious pressure oscillations in split convective form
  discretization for compressible flows, *J. Comput. Phys.* 427, 110060.
- Touber E., Sandham N. D. (2009), Large-eddy simulation of low-frequency unsteadiness in a turbulent shock-induced
  separation bubble, *Theor. Comput. Fluid Dyn.* 23, 79–107.
- Vreman A. W. (2004), An eddy-viscosity subgrid-scale model for turbulent shear flow: algebraic theory and
  applications, *Phys. Fluids* 16, 3670–3681.
- Xie Z.-T., Castro I. P. (2008), Efficient generation of inflow conditions for large eddy simulation of street-scale
  flows, *Flow Turbul. Combust.* 81, 449–470.
