# Numerical methods

FLUME discretises the [conservation laws](./models) with high-order **finite differences on point values**: the state
$\mathbf{q}_i$ is the value at the centre of cell $i$ of a uniform Cartesian block, and the semi-discrete update is the
conservative flux difference

$$\frac{d\mathbf{q}_i}{dt} = -\sum_{d} \frac{\hat{\mathbf{F}}_{d,i+1/2} - \hat{\mathbf{F}}_{d,i-1/2}}{\Delta x_d} + \mathbf{S}(\mathbf{q}_i),$$

integrated in time by an explicit Runge–Kutta scheme. The numerical flux $\hat{\mathbf{F}}_{i+1/2}$ is built so that the
difference is a high-order approximation of $\partial\mathbf{F}/\partial x$ (Shu & Osher 1989). Two space operators are
available, selected by `[numerics] scheme_space`:

| `scheme_space` | Face flux | Models |
|---|---|---|
| `weno` | WENO reconstruction of the Lax–Friedrichs-split fluxes, per characteristic field | Euler, MHD |
| `weno-riemann` | WENO interpolation of the face states, a Riemann solver and a high-order correction (Chen, Tóth & Gombosi 2016) | Euler (LLF, HLL, HLLC); MHD (LLF, HLL, HLLD) |

Both operators store the face flux at the same place, so seams, AMR reflux and the flux difference are shared.

## Grid, blocks and ghost cells

The domain is a forest of Cartesian blocks of `ni × nj × nk` cells, organised by an octree (or quadtree) with Morton
ordering. Each block carries `ngc` ghost cells per side, filled before every residual evaluation from the neighbouring
blocks (MPI exchange across ranks), from the physical [boundary conditions](./boundary-conditions), or by the coarse–fine
interpolation at an AMR seam. A direction flagged `null_x/y/z = .true.` is inactive: its fluxes are not computed and the
solution is constant along it (1-D and 2-D problems are 3-D runs with null directions). The ghost width must cover the
stencil: `ngc ≥ S` for WENO of order $2S-1$ (`ngc = 3` for WENO5), and `ngc ≥ 2` for `weno-riemann`.

## WENO reconstruction

A WENO scheme of order $2S-1$ (`[weno] scheme = weno-u-{3,5,7,9}`, $S = 2 \ldots 5$) combines $S$ candidate
polynomials $p_k$, each on $S$ points, with nonlinear weights

$$\hat v_{i+1/2} = \sum_{k=0}^{S-1} \omega_k\, p_k, \qquad
\omega_k = \frac{\alpha_k}{\sum_j \alpha_j}, \qquad \alpha_k = \frac{d_k}{(\varepsilon + \beta_k)^{m}},$$

with the optimal (linear) weights $d_k$, the Jiang–Shu smoothness indicators $\beta_k$, $\varepsilon = 10^{-6}$ and the
exponent $m = S$ ($S - 1$ for $S > 4$), the same on the host and on the device. On smooth data $\omega_k \to d_k$ and the
combination reaches order $2S-1$; the stencils crossing a discontinuity get $\omega_k \to 0$. The centred schemes
`weno-c-*` are refused by FLUME.

The absolute $\varepsilon$ ties the weights to the units of the data: $\beta_k$ scales as the square of the field, so
on a field of magnitude $10^{-4}$ the ratio $\beta_k/\varepsilon$ is $10^{8}$ times smaller than on the same field at
magnitude 1, the weights collapse to the linear ones and the scheme stops limiting (measured on Sod at an interstellar
density, $10^{-21}$ kg/m³: the total variation of the density grows by 33%, and by 178% with `weno-riemann`, whose
sensor reads the same weights; [NV-6](./verification#reference-layer-nv-5)). Two cures: the
[reference layer](./models#units-and-scaling), which hands the solver numbers of order 1, or `[weno] weights = si`,
scale-invariant weights after Don, Li, Wang and Wang (2022):

$$\alpha_k = \frac{d_k}{(\varepsilon + \beta_k/\mu^2)^{m}},$$

$\mu$ a magnitude of the reconstructed field with its units, so that $\beta_k/\mu^2$ does not depend on them. Which
magnitude matters. The mean of $|v|$ over the stencil (the descaler of Don et al.) fails on the characteristic fields of
FLUME: projected on the eigenvectors of the face state, a field such as the shear wave $\rho(v - v_\text{Roe})$, or one
carrying only a small wave, is near zero by construction, its $\mu$ is its own variation, $\varepsilon\mu^2 \to 0$ and
the weights turn fully nonlinear on smooth data (issue #49: the isentropic vortex fell from order 5.79 to 4.10). FLUME
takes instead the magnitude of the state that the field projects,

$$\mu_k = \Big\langle \sum_v |l_{kv}|\,\tfrac12\big(|f_v| + \alpha_\text{max}\,\hat q_v\big) \Big\rangle_\text{stencil}$$

for the split fluxes ($\langle \sum_v |l_{kv}|\,\hat q_v \rangle$ for the interpolated fields of `weno-riemann`), $l_k$
the left eigenvector (the identity for conservative variables), $\alpha_\text{max}$ the largest speed of the stencil
(MHD: $\max(|u_n| + c_f)$, never the cleaning speed $c_h$) and $\hat q_v \ge |q_v|$ a reference magnitude of each
component: $|\rho|$, $|E|$, $\|\rho\mathbf{u}\|_1 + \sqrt{\rho p}$ for every momentum component, $\|\mathbf{B}\|_1 +
\sqrt{p}$ for every field component, $|\psi| + c_h \hat q_B$ (GLM) or $|\psi| + \hat q_B$ (EGLM). The reference
magnitudes matter: with $|q_v|$, a field built on components that vanish (the Alfvén rows of a 2-D problem touch only
$\rho w$ and $B_z$) or are exponentially small (the field in the far field of the magnetised vortex, $\sim 10^{-21}$)
gets a descaler of its own size, the weights treat round-off as data, and the vortex went unstable at 256² (issue
#49, MV-7). $\mu_k$ bounds $|f^\pm_k|$ cell by cell, never vanishes on a physical state, and at unit scale acts like
the absolute $\varepsilon$ (on the isentropic vortex, $L_1 = 7.22\cdot10^{-8}$ at 256², order 4.97, against
$8.38\cdot10^{-8}$ with `js`). Every term scales exactly like the field ($\sqrt{\rho p}$ and $\sqrt{p}$ by exact powers
of two), and a power of two is exact in floating point, so multiplying the data by $2^n$ multiplies the reconstruction
by $2^n$ bit for bit (unit tests `test_flume_weno_weights` on the host, `test_flume_weno_interpolation_fnl` on the
device; the whole solver: the [scaling oracle](./verification#scaling-covariance)). $\mu$ is floored at the smallest
normal number. Primitive variables have no projection: a velocity or a field component crossing zero has no magnitude
of its own, so `si` with `weno-riemann` on `primitive` variables is refused. One property of `js` is lost at round-off:
the smoothness indicator of a uniform field is a quadratic form, round-off rather than zero, and the `si` descaler
varies from face to face, so a uniform field in a non-uniform state (the $B_n$ of a 1-D MHD problem) is reconstructed
with face values that differ in the last bits (RJ2a: GLM equals the run without cleaning to $10^{-13}$, not bitwise).

The default `js` is the first formula, bitwise unchanged: the kernel receives $1/\mu$, which is exactly 1 for `js`, so
the choice costs no branch in the kernels.

The WENO kernel takes its coefficient tables as arguments, so one kernel serves two purposes:

- **reconstruction** (`weno`): from point values of a flux, the face value whose difference approximates the derivative;
- **interpolation** (`weno-riemann`): from point values of the state, the state at the face. The interpolation tables
  (for WENO5: candidates exact to degree 2, combination $(3, -20, 90, 60, -5)/128$, linear weights $1/16, 10/16, 5/16$)
  are generated with exact rational arithmetic by `riemann-flux/weno_interpolation_tables.py`.

## `scheme_space = weno`: flux splitting

At each face $i+1/2$ of direction $d$, over the stencil of cells $i-S+1 \ldots i+S$:

1. the eigenvectors are evaluated at the face (Euler: the Roe average of cells $i$, $i+1$; MHD: their arithmetic average);
2. the state and the physical flux of every stencil cell are projected on the left eigenvectors,
   $w_k = \mathbf{l}_k\cdot\mathbf{q}$, $g_k = \mathbf{l}_k\cdot\mathbf{F}$ (`reconstruction_variables = characteristic`), or
   kept component by component (`conservative`);
3. each field is split with a local Lax–Friedrichs speed, $g_k^\pm = \tfrac12(g_k \pm \alpha_k w_k)$, where
   $\alpha_k = \max_{\text{stencil}}|\lambda_k|$ per characteristic field (`conservative`: one speed $\max(|u_n| + c)$,
   and with GLM at least $c_h$);
4. $g^+$ is reconstructed from the left-biased stencil and $g^-$ from the right-biased one, and the face flux is
   projected back, $\hat{\mathbf{F}} = \sum_k \mathbf{r}_k (\hat g_k^+ + \hat g_k^-)$.

For MHD the characteristic decomposition is block-diagonal: the 7×7 wave core plus the $(B_n, \psi)$ pair with speeds
$\mp c_h$ (without cleaning the $B_n$ row has speed 0 and its face flux is exactly zero).

## `scheme_space = weno-riemann`: interpolation, Riemann flux and correction

The face flux of the hybrid scheme (Chen, Tóth & Gombosi 2016) is

$$\hat{\mathbf{F}}_{i+1/2} = \mathbf{F}^{RS} + s\,\Big[c_1\,\mathbf{F}^{RS} + c_2\,(\mathbf{f}_i + \mathbf{f}_{i+1}) + c_3\,(\mathbf{f}_{i-1} + \mathbf{f}_{i+2}) - \mathbf{F}^{RS}\Big],$$

with $\mathbf{F}^{RS}$ the Riemann flux of the two interpolated face states, $\mathbf{f}_m$ the physical flux of cell $m$ and
$s \in \{0, 1\}$ a per-face sensor:

1. **Face states.** The fields of the stencil are interpolated to the face from the left and from the right
   (`reconstruction_variables`: `characteristic`, the recommended default for Euler and MHD, with the eigenvectors of
   the face; or `primitive`, $(\rho, u, v, w, p)$, plus $\mathbf{B}$ and $\psi$ for MHD). A face state with non-positive
   density or pressure is replaced by the adjacent cell's state.
2. **Riemann flux** (`riemann_solver`):
   - `llf`: $\mathbf{F}^{RS} = \tfrac12(\mathbf{f}_L + \mathbf{f}_R) - \tfrac12\alpha(\mathbf{q}_R - \mathbf{q}_L)$,
     $\alpha = \max(|u_{n}| + a)$ of the two states;
   - `hll` (Harten, Lax & van Leer 1983) with the Einfeldt speeds
     $S_L = \min(u_{n,L} - a_L, \tilde u_n - \tilde a)$, $S_R = \max(u_{n,R} + a_R, \tilde u_n + \tilde a)$ (Roe averages):
     $\mathbf{F}^{RS} = \big(S_R\mathbf{f}_L - S_L\mathbf{f}_R + S_L S_R(\mathbf{q}_R - \mathbf{q}_L)\big)/(S_R - S_L)$ for $S_L < 0 < S_R$,
     the upwind physical flux otherwise;
   - `hllc` (Toro, Spruce & Speares 1994) with the same outer speeds and the contact speed of Batten et al. (1997),
     $S_M = \dfrac{p_R - p_L + \rho_L u_{n,L}(S_L - u_{n,L}) - \rho_R u_{n,R}(S_R - u_{n,R})}{\rho_L(S_L - u_{n,L}) - \rho_R(S_R - u_{n,R})}$;
     it resolves an isolated contact exactly.
3. **Correction** (`flux_correction`), coefficients derived with sympy so that the flux difference approximates the
   derivative of point values:

   | `flux_correction` | $c_1$ | $c_2$ | $c_3$ | Order of the flux difference |
   |---|---|---|---|---|
   | `6th` | $64/45$ | $-13/60$ | $1/180$ | 6 (the WENO5 interpolation limits the scheme to 5) |
   | `4th` | $4/3$ | $-1/6$ | $0$ | 4 |
   | `none` | $1$ | $0$ | $0$ | 2 |

4. **Sensor** (`flux_correction_sensor`): with `weno`, $s = 1$ where every interpolated field keeps its nonlinear
   weights close to the linear ones, $\min_k \omega_k/d_k \ge 0.2$ on both sides of the face, and $s = 0$ elsewhere: the
   correction, which amplifies jumps by $64/45$, is switched off at discontinuities. `none` keeps it everywhere.

For MHD (`riemann_solver = llf | hll | hlld`) the solvers of `adam_flume_mhd_riemann_library` work in the frame of the
direction, so the x, y and z runs of a rotated problem stay bitwise equal: LLF, HLL, and HLLD (Miyoshi & Kusano 2005)
with its degenerate cases and a per-face fallback to HLL when a star state is not admissible or the wave speeds are out
of order; the fallbacks are counted and logged per stage (`HLLD fallbacks to HLL: N faces`). With GLM the
$(B_n, \psi)$ subsystem is solved exactly at the face and both states take the resulting $\tilde B_n$; with EGLM the
same in field units, the cleaning energy $\psi^2/2$ kept out of the solver states and added to the energy flux (see
[EGLM](./models#divergence-control-eglm)). On Ryu–Jones 2a
(256 cells) HLLD with characteristic interpolation reaches $L_1 = 3.76\cdot10^{-2}$ (splitting: $3.84\cdot10^{-2}$), with
primitive interpolation $5.13\cdot10^{-2}$; HLL is worse with either ($4.32$, $5.98\cdot10^{-2}$). The MHD verification
suite runs on `weno-riemann` through the `--numerics` option of its scripts (issue #47, M3-P3c): the linear waves and the
circularly polarised Alfvén wave converge at order 4.96–5.00, Brio–Wu and Ryu–Jones 4d need no floor and no HLLD
fallback, Orszag–Tang keeps its symmetry and conserves to $10^{-16}$, the field loop loses as much magnetic energy as
with the splitting (primitive) or less (characteristic). Both interpolations agree on smooth waves (the WENO weights are
linear there); at discontinuities characteristic interpolation is more accurate (the recommended default, issue #47
D-5), primitive keeps the 7×7 eigenvector projection out of the face kernel.

**Measured** (np 2, WENO5, 6th-order correction, CPU and FNL identical to the printed digits; details in the
[verification gallery](./verification)): on the isentropic vortex the HLLC variant reaches order 6.5 between 128² and
256² with $L_1(\rho) = 8.0\cdot10^{-8}$ (splitting: $8.4\cdot10^{-8}$); on Sod and Shu–Osher its $L_1$ error is 9% and 7%
below the splitting scheme's, on Lax 1.6% above.

## Positivity limiter

`[numerics] positivity_limiter = cell` (issue #47, D-9; both space schemes) makes each stage keep the density and the
pressure positive, so the floors are no longer needed on hard problems. Per cell, the **backbone** is the forward-Euler
update with the first-order Lax–Friedrichs fluxes of the cell states (Euler: Rusanov; MHD: the Wu speed, which adds
$|\mathbf{B}_L - \mathbf{B}_R| / (\sqrt{\rho_L} + \sqrt{\rho_R})$ to the fastest wave, and with EGLM at least $c_h$) and,
for EGLM, the second-order nonconservative sources; with the Godunov–Powell (EGLM) sources it is admissible under the
CFL bound (Wu 2018, Wu & Shu 2018). Each face carries the antidiffusive difference $\pm\Delta t/\Delta x\,(F - F^{LF})$
of the high-order flux $F$ (and EGLM the difference of the high- and second-order sources). The cell factor $\Lambda$
is the largest value in $[0, 1]$ for which every corner of the box $[0, \Lambda]^{\text{faces}}$ keeps $\rho$ and the
internal energy above their floors; both are concave in the conservative state, so each corner has a closed-form
factor (Zhang & Shu 2012), and $2^{2D}$ corners (16 in 2-D, 64 in 3-D) bound every combination. The floor of a quantity
is $\kappa$ times its backbone value, $\kappa = 0.1$: positivity alone leaves a state admissible but not usable (on the
planar Sedov blast with HLLC the former absolute floor $10^{-13}$ let one stage drain the centre cell to
$\rho = 10^{-13}$ with its energy kept, the sound speed grew by $10^6$ and the next stage, whose $\Delta t$ is that of
the step's first state, had no admissible backbone), so a stage may not take a cell below a tenth of its first-order
update. The floor is relative only: an absolute one would make the solution depend on the units of the input
(issue #49), and below a backbone value of $10^{-12}$ it was the absolute one that acted. A smooth update differs from the backbone by
$O(\Delta x)$, so the relative floor does not act there. A face takes the smaller
factor of its two cells, $F \leftarrow F^{LF} + \theta(F - F^{LF})$ (Xu 2014; Christlieb et al. 2015), so the update is
conservative and each forward-Euler step of size $\Delta t$ is admissible; an SSP Runge–Kutta stage is a convex
combination of such steps, which is why the limiter requires an SSP scheme. A face whose two cells keep $\Lambda = 1$ is
not touched, so smooth runs are bitwise unchanged. The cell factors are exchanged like a field (intra-realm copies and
MPI); a physical-boundary face takes the interior cell's factor. The log reports, per stage, the limited faces and the
cells whose backbone is inadmissible (zero in every verified case; they would fall back to the floors). A cell with a
non-finite high-order face flux takes $\Lambda = 0$ and its faces the backbone flux itself (a NaN would otherwise pass
every comparison of the corner factors and survive the blend); such cells are logged on a line of their own, and none
occurs in the verified cases.

On the Balsara–Spicer blast (β = 2.5·10⁻⁴, EGLM) the splitting scheme fails at step 3 without the limiter and reaches
$t = 0.01$ with it (at most 4% of the faces limited per stage); `weno-riemann` HLLD passes with or without it.

**At a 2:1 AMR seam** the end-of-step reflux (below) would replace the coarse face flux with the restricted fine one,
a correction that no cell factor bounds: on the blast across a refined box it drives coarse cells beside the seam to
$\rho e < 0$ (issue #50). A limited run therefore synchronises the seam flux at every stage, before the update:

- the fine seam faces take the backbone $F^{LF}(q_C, q_f)$, with the coarse cell $q_C$ as outer state;
- the coarse cell's factor $\Lambda_C$ uses, on the seam face, the means of the four fine backbone and high-order
  fluxes; the coarse backbone update with the mean of $F^{LF}(q_C, q_f)$ is a convex combination of Lax–Friedrichs
  updates of $q_C$, so it is admissible;
- both sides blend with $\theta_s = \min(\Lambda_C, \min \Lambda_f)$, so each cell stays in its corner box and the
  coarse flux is exactly the mean of the fine ones. The stage is conservative, and the reflux only corrects round-off.

The coarse states, the fine means and the factors travel over the ranks in seam skins (four reductions of seam size
per stage); on the GPU the seam cells are gathered to the host, which recomputes their factors and returns them with
the seam face fluxes.

## Time integration

The library Runge–Kutta schemes (`[runge_kutta] scheme`, listed in the [input reference](./input)) integrate the
semi-discrete system: strong-stability-preserving (SSP) and low-storage schemes. Each stage refills the ghost cells,
recomputes the auxiliary variables, applies the positivity floors (MHD) and evaluates the residual. The time step follows
the [CFL bound](./models#time-step-bound); `[time] it_max` and `time_max` stop the run.

## Adaptive mesh refinement

The block tree is refined at initialisation (`[initial_conditions] amr_iterations` passes) by the markers of the `[amr]`
section: a geometric box, a variable gradient, or the surface of an immersed solid. Neighbouring blocks differ by at most
one level (2:1 balance). At a coarse–fine face:

- the **ghost cells** of the fine block are interpolated from the coarse one (the library seam ghost fill), and the
  coarse ghosts are restricted from the fine cells. The interpolant (tricubic by default) works variable by variable, so
  next to a strong shock it can return an inadmissible state: at low β the internal energy is a small residual of
  $E - |\mathbf{B}|^2/2$, and on the Balsara–Spicer blast across a 2:1 seam the fine ghosts reach $\rho e < 0$, which
  makes the high-order flux non-finite (issue #50). After every ghost fill, an inadmissible face ghost is therefore
  blended toward the adjacent interior cell, $g \leftarrow q_\text{in} + t\,(g - q_\text{in})$, with the largest $t$ that
  keeps $\rho$ and $\rho e$ above the limiter's floors of $q_\text{in}$ (a non-finite ghost takes $q_\text{in}$). An
  admissible ghost is never touched, so the blend changes nothing else and costs one scan of the face ghosts; the log
  reports the ghosts blended, per rank;
- **conservation** is restored by the Berger–Colella **reflux**: the fluxes of every Runge–Kutta stage through the seam
  are accumulated in a flux register, weighted by the stage's coefficient in the step, and at the end of the step the
  coarse cells next to the seam are corrected by the difference between the fine fluxes and their own
  (`[numerics] reflux = .true.`; $\mathbf{B}$ and $\psi$ included). With reflux the volume integrals of a periodic box stay
  constant to round-off, without it they drift at $10^{-5}$ ([conservation example](./verification#conservation-across-amr-seams)).

With markers the tree is an octree (`ratio = 8`, and with a null axis `nk ≥ 4`) or a quadtree (`ratio = 4`, any `nk`,
`nk = 1` for a true 2-D run): a quadtree seam is 2:1 in x and y and 1:1 in z, and the ghost restriction, the coarse-fine
interpolation, the reflux and the limiter's seam synchronisation work per axis (issue #46, verified by
[MV-15](./verification#quadtree-amr)). For MHD, a cell-centred $\mathbf{B}$ cannot be kept divergence-free across a 2:1 face (Tóth & Roe 2002): the seams
inject a truncation-level $\nabla\cdot\mathbf{B}$ that GLM transports and damps; the div(B) history reports it
(`seam_max_divb`).

## Immersed boundary (Euler)

Static rigid solids (`[solids]`, `[solid_N]`, e.g. `definition = analytical_circle`) are imposed as inviscid walls:

1. a signed distance function $\phi$ marks the fluid and the solid cells;
2. the fluid state is extended into the solid by `n_eikonal` sweeps of an eikonal extrapolation along the surface
   normal $\nabla\phi$ (Jacobi sweeps: every increment is computed from the old state, so the result does not depend on the
   loop order or the thread count, and the CPU and FNL backends agree);
3. the extended state is inverted into the mirror state of a slip wall (the normal velocity reversed);
4. the fluid cells cut by the surface use the cut spacing in the flux difference, and the solid cells are masked in the
   Runge–Kutta stages. The cut spacing of a fluid cell $c$ with a solid neighbour $s$ is
   $\Delta/2 - \phi_c\,\Delta/(\phi_s - \phi_c)$; the denominator is positive ($\phi_s > 0 > \phi_c$), so it carries no
   guard (the absolute $10^{-12}$ inherited from CHASE was a length, removed with issue #49).

The solid surface can drive the AMR (the refined ring of the [shock-over-a-cylinder example](./verification#shock-over-a-cylinder)).
`mhd-ideal` with solids is refused.

## Multi-realm runs

A forest manifest glues several realms, each with its own input file, through inter-realm seams (mirror coupling of
blocks of the same resolution). The seam ghosts are filled from the peer realm at every Runge–Kutta stage
(`stage_coincident`) or once per step (`end_of_step`); the union of a split problem is bitwise equal to the
single-realm run (verification `multirealm/`).

## Output and diagnostics

Checkpoints are XH5F files (one HDF5 file per rank plus an XDMF description readable by ParaView), written every
`it_save` steps, with optional slices, residual and conservation histories, and restart files. MHD adds the div(B)
history `<basename>-divb_history.dat` (`it time max_divb l1_divb seam_max_divb`) and the `divb_tol` monitor.

## References

- Batten P. et al. (1997), On the choice of wavespeeds for the HLLC Riemann solver, *SIAM J. Sci. Comput.* 18, 1553–1570.
- Berger M. J., Colella P. (1989), Local adaptive mesh refinement for shock hydrodynamics, *J. Comput. Phys.* 82, 64–84.
- Chen Y., Tóth G., Gombosi T. I. (2016), A fifth-order finite difference scheme for hyperbolic equations on
  block-adaptive curvilinear grids, *J. Comput. Phys.* 305, 604–621.
- Christlieb A. J., Liu Y., Tang Q., Xu Z. (2015), Positivity-preserving finite difference weighted ENO schemes with
  constrained transport for ideal magnetohydrodynamic equations, *SIAM J. Sci. Comput.* 37, A1825–A1845,
  doi:10.1137/140971208.
- Einfeldt B. et al. (1991), On Godunov-type methods near low densities, *J. Comput. Phys.* 92, 273–295.
- Don W. S., Li R., Wang B.-S., Wang Y. H. (2022), A novel and robust scale-invariant WENO scheme for hyperbolic
  conservation laws, *J. Comput. Phys.* 448, 110724.
- Jiang G.-S., Shu C.-W. (1996), Efficient implementation of weighted ENO schemes, *J. Comput. Phys.* 126, 202–228.
- Miyoshi T., Kusano K. (2005), A multi-state HLL approximate Riemann solver for ideal MHD, *J. Comput. Phys.* 208, 315–344.
- Shu C.-W., Osher S. (1989), Efficient implementation of essentially non-oscillatory shock-capturing schemes II,
  *J. Comput. Phys.* 83, 32–78.
- Toro E. F., Spruce M., Speares W. (1994), Restoration of the contact surface in the HLL-Riemann solver,
  *Shock Waves* 4, 25–34.
- Tóth G., Roe P. L. (2002), Divergence- and curl-preserving prolongation and restriction formulas,
  *J. Comput. Phys.* 180, 736–750.
- Wu K. (2018), doi:10.1137/18M1168017; Wu K., Shu C.-W. (2018), doi:10.1137/18M1168042: positivity of the first-order
  and high-order schemes for ideal MHD with the Godunov–Powell source.
- Xu Z. (2014), Parametrized maximum principle preserving flux limiters for high order schemes solving hyperbolic
  conservation laws, *Math. Comp.* 83, 2213–2238, doi:10.1090/S0025-5718-2013-02788-3.
- Zhang X., Shu C.-W. (2012), Positivity-preserving high order finite difference WENO schemes for compressible Euler
  equations, *J. Comput. Phys.* 231, 2245–2258, doi:10.1016/j.jcp.2011.11.020.
