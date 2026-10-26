# Particle-in-cell

Under `physical_model = PIC` (or `AdimPIC`) PRISM couples the Maxwell field to charged macro-particles: the particles
deposit their current $\mathbf J$ and charge $\rho$ on the grid, the field is gathered at the particles, and the
particles are pushed by the Lorentz force. The keys are in the [input reference](./input#pic-pic-models-only) (`[PIC]`)
and (`[particle_injection]`).

## Particles

Each particle carries `q_pic(8, n)`: position $(x, y, z)$, velocity $(v_x, v_y, v_z)$, charge $q$ and mass $m$; the
field gathered at it is `pic_fields(6, n)`, $\mathbf D$ and $\mathbf B$.

- **`problem_type = plasma`.** The particles fill `plasma_domain`: the whole domain (`uniform`), a cylinder
  (`cilinder`, note the spelling), a slab or a box. Their number is `nint(plasma_density × volume)` of that region, split
  into neutrals (`neutral_fraction`) and, evenly, ions (H⁺) and electrons. Positions are sampled by
  `space_random_number_generator` (`random` or `layered`), optionally pairing each electron with an ion
  (`space_pairing`); velocities follow, per species, a uniform or a per-axis Maxwellian at the given temperature, or a
  pure drift, plus the drift velocity; `v_av_correction` removes the sampled mean velocity. `random_seed` makes the
  plasma reproducible.
- **`problem_type = single_particle`.** One test particle with the given position, charge, mass and velocity; its
  trajectory and fields go to `single_particle_output.dat` at every step.

## The push

The particles obey the non-relativistic Lorentz equations,

$$\frac{\mathrm d\mathbf x}{\mathrm dt} = \mathbf v, \qquad
\frac{\mathrm d\mathbf v}{\mathrm dt} = \frac{q}{m}\Big(\frac{\mathbf D}{\varepsilon_0} + \mathbf v\times\mathbf B\Big),$$

integrated with the field in the same SSP Runge–Kutta stages (`[PIC] scheme_time = runge_kutta`; SSP-22, 33 or 54).
`leapfrog` is refused at initialization on both backends. An [external field](./initial-conditions#external-fields), when
present, is added to the gathered field.

One step, per stage $s$: locate the particles in the grid; deposit $\mathbf J$ (and fold the deposit that spilled into
ghost cells back onto the cells it belongs to, across blocks and ranks); add the coils' current; evaluate the Maxwell
residual; gather $\mathbf D$, $\mathbf B$ at the particles; assign the field and particle stages. After the last stage:
combine the stages, deposit $\mathbf J$ and $\rho$ of the new positions, wrap the particles across periodic faces,
and migrate the particles that left a rank's blocks.

## Deposition and gather

`particle_weighting_model` (charge), `current_weighting_model` (current) and `field_weighting_model` (gather) choose
the shapes ([input](./input#pic-pic-models-only)):

- **Deposition**: nearest grid point (`NGP`), cloud in cell (`CIC`), triangular shaped cloud (`TSC`), the cubic to
  sextic B-splines, or a Gaussian of width `sigma` cut at `cutoff_sigma` widths.
- **Gather**: interpolations of order 0D to 6D, a Gaussian, or the *coherent* gathers `londrillo`, `totorica`,
  `esirkepov_Q5/Q6/Q8/Q10`, which are the adjoint of the matching charge deposit.

A gather that is the adjoint of the deposit in use cancels the self-force of a particle on itself; PRISM warns when
the chosen gather is not the coherent one of the deposit.

The coils' cells must stay free of plasma: a deposit there stops the run.

### Charge-conserving currents

With `current_weighting_model = conserving` the current is built from the particle trajectories so that the discrete
continuity equation holds,

$$\frac{\rho^{n+1} - \rho^n}{\Delta t} + \sum_s \beta_s\,\nabla_h\cdot\mathbf J_s = 0,$$

with PRISM's own divergence and the stage weights $\beta_s$; the matched solver is `current_conserving_solver`
(`londrillo`, the default; `totorica`, which adapts to any centred FD order; `esirkepov_Q5/Q6/Q8/Q10`, which need
`fdv_order = 6`). It requires `fd_centered`, `runge_kutta` and SSP-22, 33 or 54. **On the CPU backend it runs on one
block and one rank only.**

**With periodic faces** (totorica and esirkepov_Q\<N\>; `londrillo` and `esirkepov_tail_cleanup` stop at initialization)
the matched current is built particle by particle on the box its trajectory spans and folded modulo the periodic
length, using its separable form $J_d = Q_d\,P_{e_1}P_{e_2}\,F^E_d$ (1-D factors only), so the current through the
periodic face, the $k = 0$ mode, is kept; the matched charge is built the same way, and the coherent gather applies
$P^{\mathsf T}$ circularly along periodic axes. Without periodic faces the line solves and block deposit run unchanged.

## Particle boundaries

The particle closure of a face is derived from its field kind, with no key of its own (`adam_prism_pic_bc_object`):

- **`periodic`**: a particle beyond the face is located in the ghost halo of the block on that face (the periodic
  ghosts are images of the opposite interior), so deposits fold back and gathers read the images. Its position is
  wrapped by $\pm L$ at the end of the step only, never inside the Runge–Kutta stages (whose positions combine linearly),
  after the end-of-step deposit and before the migration.
- **every other kind**: the legacy closure, with a warning at initialization: a particle leaving the domain stops the
  run, a conserving current whose support reaches the face stops the run, and the standard deposit drops what spills
  into the boundary ghosts.

## The fields at t0

`[PIC] initialization = coherent` (the only value) builds the initial field from the initial particles: their charge is
deposited and, with `electric_field_elliptic_initialization` (default on), the electric field of that charge is solved
for, so $\nabla_h\cdot\mathbf D = \rho$ to round-off at t0; with `magnetic_field_elliptic_initialization` (default on)
the magnetic field of the initial current is added ([numerics](./numerics#t0-elliptic-solves)). Between PEC walls
along one axis, with the others periodic, `elliptic_initialization_closure = isolated` gives the mean fields of the
plasma in open space along that axis instead of zero (between grounded walls the dipole of the particle positions and
the first moment of the current shift the whole field off open space by a uniform amount). The charge in the null
space of the operator (a net charge without a Dirichlet face, or a checkerboard) cannot be reached by the solve: it is
measured and warned, and a checkerboard under a matched deposit stops the run as a deposit error.

## Charge conservation diagnostic

Every step writes `<basename>-charge_conservation_history.dat`: the total particle and deposited charge and their
difference, the continuity residual $\max\lvert(\rho^{n+1}-\rho^n)/\Delta t + \sum_s\beta_s\nabla_h\cdot\mathbf J_s\rvert$
in the L1, L2 and L∞ norms, the net and absolute charge change, and the Gauss residual
$\nabla_h\cdot\mathbf D - \rho$ in the three norms. With the conserving current the continuity residual is at
round-off; with the other deposits it measures their charge error.
