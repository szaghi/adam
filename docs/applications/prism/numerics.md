# Numerical methods

PRISM discretizes the equations of the [physical models](./models) by the method of lines: a spatial operator gives the
residual $\mathrm d\mathbf q/\mathrm dt = \mathbf R(\mathbf q)$ at every cell, and a Runge–Kutta scheme advances it. The
spatial scheme is `[numerics] scheme_space` ([input](./input#numerics)).

## Space discretization

### Centred finite differences (`fd_centered`)

The state lives at the cell centres. Every derivative is the centred first derivative of order $2s$,
$s = $ `[fdv] fdv_order`$/2$ (orders 2, 4, 6, 8, 10):

$$\left.\frac{\partial q}{\partial x}\right|_i \approx \frac{1}{\Delta x}\sum_{m=1}^{s} a_m^{(s)}\,(q_{i+m} - q_{i-m}),$$

with $a^{(1)} = \tfrac12$, $a^{(2)} = (8, -1)/12$, $a^{(3)} = (45, -9, 1)/60$, $a^{(4)} = (672, -168, 32, -3)/840$
(`FD1_CC`, `adam_fdv_operators_library`). The residual is evaluated directly from the curls:

$$R_{\mathbf D} = \nabla_h\times\frac{\mathbf B}{\mu_0} - \mathbf J, \qquad R_{\mathbf B} = -\nabla_h\times\frac{\mathbf D}{\varepsilon_0},$$

plus, with GLM cleaning, $-\nabla_h\varphi$, $-\nabla_h\psi$ and the scalars' residuals
$-c_h^2(\nabla_h\cdot\mathbf D - \rho) - \alpha\varphi$, $-c_h^2\nabla_h\cdot\mathbf B - \alpha\psi$. `ngc >= s` is checked at
initialization.

### The discrete curl and divergence

The curl, the divergence and the gradient are all built from the same one-dimensional stencil, so on a uniform grid
the discrete operators commute exactly as their continuous counterparts:

$$\nabla_h\cdot(\nabla_h\times\mathbf F) \equiv 0 .$$

Hence $\partial_t(\nabla_h\cdot\mathbf B) = -\nabla_h\cdot\nabla_h\times(\mathbf D/\varepsilon_0) = 0$: a field that starts
divergence-free stays so to round-off, with no cleaning (a mimetic identity of the centred operators). The same holds
for $\nabla_h\cdot\mathbf D$ with a divergence-free $\mathbf J$. The identity fails across a 2:1 AMR face
([below](#the-2-1-seam-and-div-b)).

### Centred finite volumes (`fv_centered`)

The fluxes $\mathbf F_d(\mathbf q)$ are evaluated at the cell centres (ghosts included), reconstructed at the faces by
the centred FV reconstruction of order $2s$, and differenced:

$$R_v = -\sum_d \frac{F^d_{v,\,i+1/2} - F^d_{v,\,i-1/2}}{\Delta x_d} + S_v .$$

On a uniform grid the centred FV derivative reduces to the FD stencil ([#84](https://github.com/szaghi/adam/issues/84)
E3a measured both against the elliptic solver). The FV path carries the face fluxes, so its 2:1 faces are refluxed
([#38](https://github.com/szaghi/adam/issues/38)). The coils' current density is stamped with the operator of the
chosen family (`[fdv] fdv_scheme`), so an FV run evolves an FV-built source.

### WENO flux splitting (`weno`, CPU only)

A global Lax–Friedrichs splitting in characteristic (`reconstruction_variables = characteristics`) or conservative
variables: with the left eigenvectors $\mathbf L$ of the Maxwell system, $w = \mathbf L\mathbf q$, $g = \mathbf L\mathbf F$,

$$f^{\pm} = \tfrac12\,(g \pm \lambda_{\max} w), \qquad \lambda_{\max} = \max(1,\chi)\,c_0,$$

each split flux is reconstructed at the faces by the library's upwind WENO (`[weno] scheme`, `weights`), projected back
with the right eigenvectors, and differenced. `ngc >= [weno] S` is checked. The FNL backend refuses `weno`.

## Time integration

`[numerics] scheme_time = runge_kutta` is the only scheme that runs; `leapfrog`, `blanes_moan` and
`commutator_free_magnus` are refused at initialization on both backends
([#55](https://github.com/szaghi/adam/issues/55)). The Runge–Kutta scheme is `[runge_kutta] scheme`
([library reference](/reference/library-input#runge-kutta)):

| Scheme | CPU | FNL |
|---|---|---|
| SSP 11, 22, 33, 54 | yes | yes |
| low-storage 1, 2, 3 | yes | yes |
| Yoshida (4-stage symplectic) | yes (EM) | refused |

PIC runs, the PML, the charge-conserving current and the staged path (2:1 faces, forests) need an SSP scheme (PIC:
SSP-22, 33 or 54). At every stage the coils' current is evaluated at the stage time $t + \gamma_s\Delta t$.

**Time step.** $\Delta t = \mathrm{CFL}\,\Delta x_{\min}/\lambda_{\max}$, with $\Delta x_{\min}$ the smallest cell size
over every block and rank and $\lambda_{\max} = \max(1,\chi)\,c_0$, capped by $1/\alpha$ (GLM damping) and by
$C_\gamma/\gamma_{\max}$ (the [PML](./pml)). `[time] it_max > 0` runs that many steps; otherwise the run ends at
`time_max`, the last step shortened to reach it.

## Adaptive mesh refinement

The grid is the library's octree or quadtree of Cartesian blocks ([concepts](/concepts/architecture)). PRISM refines at
initialization only: `[initial_conditions] amr_iterations` passes each set the initial condition and run the markers
of `[amr_marker_N]`:

- mode 1 `primitive-box`: refine the blocks inside a box to `target_level`;
- mode 3 (`AMR_TV`): refine on the total variation of the coils' current density.

The other markers do nothing in PRISM ([input](./input#library-sections-in-prism)). `[amr] frequency` is read and
unused: there is no runtime regrid.

At a 2:1 face the coarse ghosts are restricted from the fine cells and the fine ghosts interpolated from the coarse
ones (the library's seam machinery, octrees and quadtrees, [#46](https://github.com/szaghi/adam/issues/46)); the stages
then run on the staged path of the forest, which needs an SSP scheme.

### The 2:1 seam and div(B)

Across a 2:1 face the two sides difference at different resolutions, the composite operator is no longer a tensor
product, and $\nabla_h\cdot\nabla_h\times \ne 0$ there: the seam injects a $\nabla\cdot\mathbf B$ source in the
grid-oscillation modes. It is truncation-order, $O(h^p)$: it converges under refinement, but it is unbounded in time at
fixed $h$ and independent of $\Delta t$; $\lVert\mathbf B\rVert$ stays bounded, so it is a constraint violation, not an
energy instability. Five local fixes were ruled out ([#29](https://github.com/szaghi/adam/issues/29)); the state is
**accept-truncation**, with an upwinded or dissipative GLM coupling at the seams open as
[#99](https://github.com/szaghi/adam/issues/99) and a staggered-$\mathbf B$ route parked as
[#30](https://github.com/szaghi/adam/issues/30).

- **Monitor.** `[IO] seam_divB_tol` (off by default) warns, on a run with 2:1 faces, when the max
  $\lvert\nabla\cdot\mathbf B\rvert$ outside the absorbing layers exceeds it (a global maximum, not a seam-local one), and `[IO] seam_divB_error = .true.` stops the run instead
  ([library reference](/reference/library-input#io)).
- **Mitigation.** `divergence_correction = hyperbolic` with `constrained_transport = B` or `DB` and a damping
  $\alpha > 0$ delays the growth; it does not cure it.

The details and the measurements are on the [forest page](/concepts/forest#the-2-1-seam-is-not-div-free-accept-truncation).

## Multi-realm runs

A [forest manifest](/reference/library-input#forest-manifest-multi-realm) glues several PRISM realms through seams.
Both backends implement the realm contract of the forest, fast and staged paths and the seam exchange. A seam fills its
ghosts once per step (`coupling_cadence = end_of_step`, α, the default) or at every Runge–Kutta stage
(`stage_coincident`, β). **A 1:1 same-resolution seam needs β to stay divergence-free**: under α the stages after the
first read stale seam ghosts ([#31](https://github.com/szaghi/adam/issues/31)). The divergence diagnostic refills the
inter-realm seam on the committed state before it measures. Worked cases:
`src/tests/prism/regression/rmf-2realm*/`; the cadences are explained on the
[forest page](/concepts/forest#coupling-cadence-α-vs-β).

## t0 elliptic solves

When the initial state has a charge or a current, PRISM solves for the fields they generate before the first step
(`impose_pic_fields_time_zero`), with FLAIL, the library's elliptic solver
([#84](https://github.com/szaghi/adam/issues/84)):

- **Electric field** (PIC models, `[PIC] electric_field_elliptic_initialization`, default on):
  $\mathcal D(\mathcal G)\,\varphi = -\rho/\varepsilon_0$, then $\mathbf D \leftarrow \mathbf D - \varepsilon_0\,\mathcal G\varphi$.
- **Magnetic field** (when $\mathbf J \ne 0$ at t0: PIC currents, or coils with `time_delay <= 0`; under PIC also
  `magnetic_field_elliptic_initialization`): $\mathcal D(\mathcal G)\,A_d = -\mu_0 J_d$ per component, then
  $\mathbf B \leftarrow \mathbf B + \nabla_h\times\mathbf A$.

$\mathcal D(\mathcal G)$ is PRISM's own operator, the centred first derivative of order `fdv_order` composed with
itself on each axis, so the corrected field satisfies $\nabla_h\cdot\mathbf D = \rho$ to round-off with PRISM's own
divergence. Every centred first derivative annihilates the checkerboard mode, so this wide operator has a second null
branch at $\theta = \pi$, which the solver drops exactly. The face kinds come from the
[boundary map](./input#elliptic-boundary-map):

- periodic and PEC faces: a direct FFT solve (FFTW on the CPU, cuFFT on FNL);
- zero-ghost (Dirichlet) faces: a preconditioned CG on rank 0, to $10^{-12}$;
- exact/open faces (Neumann, Silver–Müller): the free-space potential of the whole source in the ghosts, then the
  same CG.

On a singular problem (no Dirichlet face) the mean of the forcing is removed per component; for $\mathbf A$ that mean is
the uniform current, which by periodicity carries no magnetostatic field. A refined grid is refused at the solve (the
wide operator across a 2:1 face is not defined), and so is `fdv_order = 10` (FLAIL accepts 2 to 8).

**Mean fields between PEC walls (`[PIC] elliptic_initialization_closure`).** With one axis $n$ closed by PEC faces and
the others periodic (or null), the transverse mean ($k_\perp = 0$) of the fields is a 1-D problem, and the PEC cavity
leaves two of its invariants free: the potential difference between the walls, $\langle D_n\rangle$, and the trapped
magnetic flux, $\langle\mathbf B_t\rangle$. The solves set both to zero (`grounded`, the default). `isolated` sets them
to the values of the same sources in open space along $n$, where the 1-D field is half the difference of the source
integrated from either side:

$$\langle D_n\rangle = -\frac{1}{V}\sum \rho\,(x_n - x_{mid})\,dV, \qquad
\langle\mathbf B_t\rangle = \mu_0\,\hat n\times\frac{1}{V}\sum \mathbf J\,(x_n - x_{mid})\,dV,$$

with $x_{mid}$ the domain centre and $V$ its volume (a net charge or current included). The mean of the field the solve
adds is shifted to these values; the initial-condition field keeps its own mean. A uniform field has no divergence nor
curl and is even about a PEC face for $D_n$ and $\mathbf B_t$, as PRISM's PEC ghosts reflect them, so the boundary cells
stay consistent with the evolution exactly as under `grounded`. Only the mean mode is open: the modes $k_\perp \ne 0$
keep their images in the walls, off open space by about $e^{-k_\perp d}$ at a distance $d$ from the sources, so the walls
must be far from them for the whole field to match open space. An exact/open face (Neumann, Silver–Müller) is not the
same closure: its ghosts hold the free-space potential of the sources alone, without their periodic images, and do not
match the hyperbolic ghosts of that face.
