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
sensor reads the same weights; [NV-6](./verification/infrastructure/scaling-covariance)). Two cures: the
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
device; the whole solver: the [scaling oracle](./verification/infrastructure/scaling-covariance)). $\mu$ is floored at the smallest
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
$\mp c_h$ (without cleaning the $B_n$ row has speed 0 and its face flux is exactly zero). The core is the Roe–Balsara
set in the normalisation of Stone et al. (2008, appendix B), with $\alpha_f^2 = (a^2 - c_s^2)/(c_f^2 - c_s^2)$ and
$\alpha_s^2 = (c_f^2 - a^2)/(c_f^2 - c_s^2)$. At high β with $\mathbf{B}$ near the face normal one of the two numerators
is a difference of $O(a^2)$ speeds far smaller than them, so FLUME evaluates the smaller one without subtracting,
$2a^2c_t^2/\big((c_f^2 - c_s^2) \pm (v_{A,n}^2 + c_t^2 - a^2)\big)$ with $c_t^2 = B_t^2/\rho$, and the other as its
complement. Subtracted directly, the normalisation carried an absolute error of $\sim\varepsilon^{1/2}$ that the
characteristic WENO turned into a $10^{-6}$ difference of $\mathbf{B}$ between compilers
([#79](https://github.com/szaghi/adam/issues/79)).

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
[verification gallery](./verification/)): on the isentropic vortex the HLLC variant reaches order 6.5 between 128² and
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

With self-gravity the register also carries the seam gravity work ([#81](https://github.com/szaghi/adam/issues/81)
P5b), so the end-of-step correction is not round-off: on a star whose surface a 2:1 seam cuts, in a $10^{-4}$
atmosphere, it drove one coarse atmosphere cell beside the seam to $\rho e < 0$ at step 688, and the run stopped
non-finite ([#88](https://github.com/szaghi/adam/issues/88)). A limited run therefore gives each coarse cell the
largest fraction $\theta \in [0, 1]$ of its end-of-step correction that keeps $\rho$ and $\rho e$ at least $\kappa$
times their values before it (`reflux_positivity_theta`; $\rho e$ is concave in $\theta$, so a bisection finds the
bound). A cell whose full correction is admissible receives it exactly, so a run that never clips is unchanged; a
clipped cell gives up conservation by the part it did not receive, and each step that clips logs
`reflux positivity: +N coarse cells limited`. On that star one cell was limited on five steps, and the run reached
its end. `reflux_cadence = stage` needs no such bound: its correction acts inside the stages, before the limiter.

The coarse states, the fine means and the factors travel over the ranks in seam skins (four reductions of seam size
per stage); on the GPU the seam cells are gathered to the host, which recomputes their factors and returns them with
the seam face fluxes.

## Dual energy

`physical_model = euler-de` ([#91](https://github.com/szaghi/adam/issues/91)) adds the advected entropy $\rho s$,
$s = p/\rho^\gamma$, to the Euler state ($nv = 6$): $\partial_t(\rho s) + \nabla\cdot(\rho s\,\mathbf u) = 0$, a pure
conservation law carried by the same flux machinery (the sixth characteristic field, eigenvalue $u_n$, rides the
acoustic and entropy waves of the Euler system), refluxed and prolonged like the other variables, with no gravity source.
Once per step, after the reflux and any regrid, `sync_dual_energy` compares each cell's internal energy
$e = E - |\rho\mathbf u|^2/(2\rho)$ with $\eta E$, $\eta$ = `[physics] dual_energy_eta` (default 0.05; the
cell-local test of Bryan et al. 2014). In a cold cell, $e < \eta E$, the entropy is a floor: $e = \max(e, e_s)$,
$e_s = \rho s\,\rho^{\gamma-1}/(\gamma-1)$, and $E$ is rebuilt from it. Then everywhere $\rho s$ follows $e$. A cell
below its adiabat was drained by the truncation error of the total-energy scheme in a cold supersonic flow and is lifted
back; a cell above it was heated by a shock, which the advected entropy does not see, and keeps its energy. The floor
only adds energy (logged per step); a reset of the energy from the entropy, as in Cholla (Villasenor et al. 2021), would
remove a cold shock's heat ([DE](./verification/euler/dual-energy)). The runtime-AMR collapse of
[#88](https://github.com/szaghi/adam/issues/88) reaches its end with it ([SV-17](./verification/self-gravity/polytrope)). CPU and FNL run the same
rule (`sync_dual_energy`, `sync_dual_energy_dev`).

## Passive scalars

`[physics] passive_scalars = N` ([#92](https://github.com/szaghi/adam/issues/92), $0 < N \le 4$) selects the scalar
variant of `euler` or `euler-de`: four more conserved variables $\rho s_k$, advected with the flow,
$\partial_t(\rho s_k) + \nabla\cdot(\rho s_k\,\mathbf u) = 0$, by the machinery of the dual-energy entropy (each a
characteristic field of eigenvalue $u_n$, carried by the acoustic and entropy waves, its specific value averaged like
the velocity), refluxed and prolonged like the other variables. They do not act on the flow: the flow fields of a run
with scalars are bitwise those of the same run without them. The capacity is fixed (D-M10-1): the kernels are compiled
for four slots, the first $N$ used and the others zero, so a scalar run pays the projection of the wider state. A scalar
may vanish (a face state needs it non-negative, the entropy positive); WENO keeps it bounded only to its own accuracy:
across a Sod contact a mass fraction moving from 1 to 0 over- and undershoots by $5\cdot10^{-4}$.

## Isothermal and barotropic equations of state

`[eos] type = isothermal` ([#100](https://github.com/szaghi/adam/issues/100), M13 P1) closes the Euler equations
with $P = c_0^2\rho$: mass and momentum only, the sound speed $c_0$ everywhere. The state keeps the Euler layout
(D-M13-1): the energy slot $\rho E$ is inert, with zero flux, zero source and an eigenvalue of 0, and stays at its
initial 0. It is its own model (`MODEL_EULER_ISO`) with its own pointwise library
(`adam_flume_euler_iso_library`) behind the same interface as the Euler one, so the model-agnostic kernel bodies run
it unchanged; the EOS argument the kernels pass as `gamma` is $c_0^2$ for this model.

In direction $\mathbf n$ the waves are $u_n - c_0$, the inert slot (0), $u_n$ twice (the transverse velocities) and
$u_n + c_0$, with

$$\mathbf r_{1,5} = (1,\ \mathbf u \mp c_0\mathbf n,\ 0),\qquad \mathbf l_{1,5} = \Big(\pm\frac{u_n \pm c_0}{2c_0},\
\mp\frac{\mathbf n}{2c_0},\ 0\Big),$$

the transverse waves $\mathbf r = (0, \mathbf t, 0)$, $\mathbf l = (-\mathbf u\cdot\mathbf t, \mathbf t, 0)$. The Roe
average of two isothermal states is the $\sqrt\rho$-weighted velocity with the same $c_0$. The `weno` scheme splits on
this eigensystem; `weno-riemann` interpolates the characteristic fields or the primitives $(\rho, \mathbf u)$, the
inert slot in place of the pressure, and the face state is admissible when $\rho > 0$. `llf` and `hll` (Einfeldt
speeds) are the standard fluxes on the two acoustic speeds. The isothermal HLLC is HLL on the mass and the normal
momentum, the transverse momenta (and any advected scalar) the HLL mass flux times the upwind specific value: the
isothermal system has one wave between the acoustic ones, the contact, and it carries only the transverse velocity, so
the star state needs no pressure. A Riemann flux always returns 0 in the inert slot.

The isothermal Riemann problem is solved exactly with no energy equation: the shock jump
$u_* = u_K \mp c_0(\rho_*/\rho_K - 1)\sqrt{\rho_K/\rho_*}$ on the compressive side, the fan
$u_* = u_K \mp c_0\ln(\rho_*/\rho_K)$ on the other; two streams at Mach $M$ compress to $\rho_* = \rho M_1^2$, $M_1$ the
Mach number in the frame of the shock. `src/tests/flume/verification/eos/eos_oracle.py` builds these solutions
independently of the code, and the runs match them ([EO](./verification/euler/isothermal-eos)). A null direction
freezes its momentum (`null_freeze`) for this model as for every Euler model: a 1-D run cannot carry a transverse
velocity. CPU and FNL (P2) instantiate the same model-agnostic kernel bodies on the same library, so the backends run
the same arithmetic.

`[eos] type = barotropic` (P3) is a continuous piecewise polytrope, $P = K_k\rho^{\Gamma_k}$ on
$\rho_{k-1} \le \rho < \rho_k$, $K_{k+1} = K_k\rho_k^{\Gamma_k - \Gamma_{k+1}}$, with $c^2 = dP/d\rho = \Gamma_k P/\rho$
(D-M13-4: the speed of the Jeans relation). It is its own model, `MODEL_EULER_BARO`, on the same body as the
isothermal gas: the body takes the EOS from four hooks of the including module (the pressure, the temperature, the
sound speed and the speed of the Roe average), and the two libraries differ only there. The isothermal hooks are the
closed forms; the barotropic ones read the realm's table. The kernels receive the EOS through one scalar
(`physics%eos_arg`: $\gamma$, $c_0^2$ or the table id), so the tables live in a fixed-size registry
(`adam_flume_eos_tables_library`): each realm registers its table at initialization, identical tables share an
entry, and the FNL backend copies the registry to the device once.

The Roe average keeps the $\sqrt\rho$-weighted velocity and takes the secant speed $\tilde c^2 = [P]/[\rho]$, the one
for which the Roe matrix satisfies $[F] = \tilde A[U]$ for any $P(\rho)$; $P$ increases with $\rho$, so it is positive,
also across breaks. Below a relative density jump of $10^{-6}$ the secant loses digits to cancellation (about
$\varepsilon$ over the jump) while $c^2$ at the mean density errs by the jump squared, so there the mean is used. Where
$\Gamma$ increases upward (a convex break) a rarefaction through the break splits into two fans with a plateau at the
break density; where it decreases, $\rho c$ drops and the system is not genuinely nonlinear there, so Riemann
solutions crossing it contain composite waves. The schemes (entropy-satisfying fluxes) do not depend on convexity;
the exact solver of the verification refuses such a law, and the default table of Mathew & Federrath (2021), which
is non-convex at its $1.4 \to 1.1$ break, is verified by self-convergence instead
([EO barotropic](./verification/euler/barotropic-eos)).

Under gravity (P4) both gases take the momentum sources of the gamma-law gas, $\rho\mathbf g$ (constant or from the cell
potential, self-gravity included), and no energy source: the residual of the energy slot is cleared after the sources
(`clear_inert_energy`), the fourth-order energy correction of the constant gravity is not computed (nor its `ngc >= 4`
stencil required), and the AMR reflux carries no seam gravity work. Where the ideal gas reads $c_s^2$ or a thermal
energy from the energy slot, these gases take them from the equation of state (D-M13-4): the Jeans marker's
$\lambda_J = \sqrt{\pi c_s^2/(G\rho)}$ and the sinks' threshold density $\pi c_s^2/(4 G r_\mathrm{acc}^2)$ with
$c_s^2 = dP/d\rho$, and the sinks' Jeans test $|E_\mathrm{grav}| > 2E_\mathrm{th}$ with the virial thermal energy
$E_\mathrm{th} = \tfrac{3}{2}\int P\,dV$ (the ideal gas keeps its internal energy, the same at $\gamma = 5/3$; at
$\gamma \to 1$ the internal energy $\int P/(\gamma - 1)\,dV$ diverges: at $\gamma = 1.001$ it is a thousand times the
virial term, so the test could not pass, and the SIS runs of the gamma-law gas keep creation off). The cell positivity
limiter is the density-only one: the models' admissible quantity (`internal_energy` of the kernels) is the density, so a
stage, and the end-of-step reflux, keep every cell's density above `POSITIVITY_LIMITER_KAPPA` times its first-order
backbone value, and an inadmissible seam ghost takes its interior cell ([AMR](#adaptive-mesh-refinement)). A corner
whose density the chord takes to the floor would keep its momentum, so its velocity $\rho\mathbf{u}/\rho$ would be
extrapolated up to $1/\kappa$ times the backbone's; the ideal gas's internal-energy chord rejects that state, and the
barotropic gases take a **kinetic chord** in its place, the same closed form on the concave
$f = \rho W - |\mathbf{m}|^2/(2\rho)$ with $W = (4|\mathbf{u}_b|^2 + c^2)/2$ of the backbone ($c^2 = dP/d\rho$ at its
density; the reflux limiter $c_0^2$): $f \ge \kappa f_b$ bounds the corner's speed by about
$\sqrt{4|\mathbf{u}_b|^2 + c^2}$, twice the backbone's. The factor matters for cold gas: with
$W = (|\mathbf{u}_b|^2 + c^2)/2$ a stage at Mach $M$ that outruns its backbone by $0.45\,c/M$ is already clipped, and
the sinks' cold accretion (eg6) drifted its momentum by 13 % of the exchange, 2.6 % with the factor. Without it the
Boss–Bodenheimer cloud stalled at $1.32\,t_\mathrm{ff}$, the ambient beside the corners of its refined cube drained to
$10^{-5}$ at $|u| = 365$ ([EO P4](./verification/euler/barotropic-gravity)).

### Isothermal and barotropic MHD

With `physical_model = mhd-ideal` (P5) the two gases close ideal MHD without an energy equation, the total pressure
$P(\rho) + |\mathbf B|^2/2$ (Mignone 2007, J. Comput. Phys. 225, 1427, eqs. 1–6). The models are
`MODEL_MHD_ISO`, `MODEL_MHD_BARO` and their mixed-GLM twins `MODEL_MHD_GLM_ISO`, `MODEL_MHD_GLM_BARO`. They keep
the MHD layout $(\rho, \rho\mathbf u, \rho E, \mathbf B[, \psi])$, the energy slot inert as for the Euler gases, and one
library body (`adam_flume_mhd_iso_library_body.INC`) behind the EOS hooks of the isothermal and barotropic modules.
The normal field is a parameter of the 1-D core, as in the adiabatic library: decoupled at speed 0, or the GLM pair
$(B_n, \psi)$ at $\mp c_h$.

**Eigensystem.** In direction $\mathbf n$ the six waves $u_n \mp c_f$, $u_n \mp c_a$, $u_n \mp c_s$ act on
$(\rho, \rho u_n, \rho u_{t1}, \rho u_{t2}, B_{t1}, B_{t2})$; there is no entropy mode. The Jacobian holds the gas law
only through $a^2 = dP/d\rho$, so the eigenvectors of Kim et al. (1999, ApJ 514, 506, eqs. 2-9 to 2-53) serve both
gases with $a^2 = c_0^2$ or the table's $dP/d\rho$. Kim's fast normalisation
$\alpha_f = \sqrt{c_f^2 - b_n^2}/\sqrt{c_f^2 - c_s^2}$ is the Ryu–Jones $\alpha_f$ times $c_f/a$. FLUME computes it from
the cancellation-free pair of the adiabatic core ([#79](https://github.com/szaghi/adam/issues/79)), with that core's
degenerate values and the sign of Kim's eq. 2-53. The characteristic projection takes the eigenvectors of the arithmetic
face average.

**HLLD.** The isothermal HLLD of Mignone (2007, section 3) has three states between the outer speeds $S_{L,R}$ and the
rotational speeds $S^*_{L,R} = u^* \mp |B_n|/\sqrt{\rho^*}$:
- $\rho^*$, $m_n^*$ and the first two flux components are the HLL averages, with $u^* = F_\rho^\mathrm{hll}/\rho^\mathrm{hll}$;
- the transverse components follow eqs. 30–37, and the flux eq. 38.

The outer speeds are those of the adiabatic library (Miyoshi & Kusano 2005, eq. 67). They bound Mignone's Davis
estimates, whose eq. 39 is misprinted with $S_R$ a minimum: it must be the maximum. In the degenerate case,
$(S - S^*_L)(S - S^*_R) \to 0$, the outer transverse state is kept. The flux falls back to HLL when $\rho^* \le 0$ or
when the speeds are out of order. They can be: $u^*$ is the advective velocity, which a strong density jump drives far
off, and on Mignone's own tube 1 ($\rho_L/\rho_R = 10$) the first face has $S^*_R = 3.90 > S_R = 3.37$.

LLF and HLL use each state's own $B_n$, so its jump is dissipated as in the adiabatic library; HLLD uses the average.
With $\mathbf B = 0$ the HLLD flux is the isothermal HLLC's, and an isolated rotational discontinuity is resolved
exactly.

`src/tests/flume/verification/eos/mhd_oracle.py` checks the two papers independently of the code ($LR = I$ and
$R\Lambda L = \partial F/\partial U$ at random and degenerate states, both gases; HLLD consistency, exactness and its
hydrodynamic limit), and the Fortran library agrees with it to $2\cdot10^{-15}$ on 525 cases.

**Positivity, cleaning and refusals.**
- **Floors:** the `[mhd]` density and pressure floors are refused, since there is no pressure to floor.
- **Positivity limiter:** `positivity_limiter = cell` takes the density and the kinetic chord of the Euler gases,
  $c^2 = dP/d\rho$ of the backbone, and an inadmissible seam ghost takes its interior cell (P5b; P5a had left the
  MHD twins on the ideal gas's blend).
- **Cleaning:** GLM is the mixed GLM of the ideal MHD (refused with that limiter, as there). EGLM is refused: its
  Derigs sources are built on the energy equation.
- **Gravity:** it moves the momentum only, as for the Euler gases.
- **Initial conditions:** `orszag-tang`, `mhd-vortex` and the `entropy` wave of `mhd-linear-wave` need a gamma-law
  gas and are refused.
- **Backends:** CPU (P5a) and FNL (P5b), the FNL kernels the device twins of the CPU ones on the same library; the
  barotropic table is read from the device registry.

Verification: [EO P5](./verification/mhd/barotropic-mhd).

## Radiative cooling

`[cooling]` ([#93](https://github.com/szaghi/adam/issues/93), M11) adds optically thin cooling,
$\partial_t(\rho e) = -n^2\Lambda(T)$, $n = \rho/m$, $T = p/(\rho R)$, to the hydrodynamic models. It is operator
split (D-M11-1): once per step, after the dual-energy sync, `apply_cooling` cools every cell isochorically over the step
just advanced, at fixed $\rho$ and $\rho\mathbf u$, and lowers $E$ by the energy it removes. The coupling with the
dynamics is first order in $\Delta t$ (Lie splitting); the cooling itself is exact.

$\Lambda$ is the piecewise power law through the nodes $(T_k, \Lambda_k)$, so Townsend's exact integration (Townsend
2009) applies. Isochorically $dT/dt = -(\gamma - 1)\rho\Lambda(T)/(m^2 R)$; with the temporal evolution function
$Y(T) = (\Lambda_N/T_N)\int_T^{T_N} dT'/\Lambda(T')$, closed form on each segment, $dY/dt$ is the constant
$(\gamma - 1)\rho\Lambda_N/(m^2 R T_N)$. One step adds $\Delta t$ times it to $Y(T)$ and inverts $Y$ on the segment
that holds the result: exact and unconditionally stable for any $\Delta t/t_\mathrm{cool}$, with no cooling limit on
the step. The floor is exact: a cell whose $Y$ passes $Y(T_1)$ stops at $T_1$, and below $T_1$ nothing cools (D-M11-2:
cooling only, no heating). Near a segment exponent $\alpha_k = 1$ both the closed form and its inverse are
$(1 - x^\varepsilon)/\varepsilon$ with $\varepsilon$ small; they are evaluated through `expm1` and `log1p` (written as
series below $10^{-2}$), so every exponent keeps full precision. A cell that does not cool keeps its bits.

With `euler-de` the cooled cells take the entropy of their cooled energy, $\rho s = (\gamma - 1)e\rho^{1-\gamma}$,
computed from the internal energy that the stored $E$ holds: radiation lowers the entropy for real, and an entropy left
on the old adiabat would let the next dual-energy floor lift a cold cooled cell back to it and undo the cooling. Every
cell then leaves the step with its entropy consistent with its energy, as the sync leaves it.

`[cooling] dt_factor` (default off) limits the step to `dt_factor` times the shortest cooling time
$t_\mathrm{cool} = \rho e/(n^2\Lambda)$, for runs that need the coupling resolved. The radiated energy is summed over
the ranks; the conservation history adds it back to $E$ (`E_plus_radiated`, constant to round-off). On a uniform box the
temperature follows the ODE solution to $4\cdot10^{-14}$ relative, at the CFL step as with the limiter, on CPU and FNL
alike ([CO](./verification/euler/radiative-cooling)). The FNL backend runs the same
library arithmetic in its device kernels (`apply_cooling_dev`, `apply_cooling_entropy_dev`), the tables at the fixed
size of 64 nodes.

## Sink particles

`[sinks]` ([#94](https://github.com/szaghi/adam/issues/94), M12) adds point masses that act through self-gravity only.
The list is replicated (D-M12-1): every rank holds every sink and updates it identically from MPI-reduced sums, so no
sink is owned and none migrates; the replicated tree locates any sink's block on any rank.

**Gravity (D-M12-2).** Each sink's mass is spread over its cloud-in-cell stencil (the 8 cell centres of the level of
its block around it) and added to the density the Poisson solver sees, through the packed path of both solvers: the
periodic and isolated boundaries (the James step) and the AMR composite operator apply to the sinks as to the gas, and
`W` $= \tfrac12\sum\rho_\mathrm{tot}\Phi\,dV$ includes them. A sink feels the centred difference of $\Phi$ at those
cells, interpolated with the same weights: the deposit-interpolation pair of Hockney & Eastwood (1988), free of
self-force on a uniform grid. The gas feels $\Phi$ through its face gradients (the potential path of self-gravity), so
gas and sinks exchange momentum equal and opposite to truncation ($4\cdot10^{-5}$ relative on the SI si1 box). Through
the mesh the sink gravity is softened at about two cells.

**Time stepping.** The gas step runs with the sinks at their positions $x^n$ in the Poisson source; then, once per step
on the committed state, kick-drift-kick: $v \mathrel{+}= \tfrac{\Delta t}2 a$, $x \mathrel{+}= \Delta t\,v$, the
potential re-solved, $a$ gathered at the new positions, $v \mathrel{+}= \tfrac{\Delta t}2 a$. Second order for the sinks,
first order in their coupling with the gas. Two time step limits: a sink moves at most half a cell of its level per
step, and gas a cell away from it, pulled by $Gm/h^2$, at most half a cell from rest ($\Delta t \le \tfrac12\sqrt{h^3/(Gm)}$,
the acceleration limit of Federrath et al. 2010: a cold dense medium has a CFL step far longer than its fall onto a sink,
and the self-gravity source is not limited by the CFL condition).

**Creation (D-M12-3).** Once per step, on the committed state (its ghost cells and its potential refreshed), a cell of
the finest level present (with runtime regridding, of the finest level allowed) forms a sink when its density exceeds $\rho_\mathrm{thr} = \pi c_s^2/(4 G r_\mathrm{acc}^2)$
($c_s^2 = \gamma p/\rho$ of the cell, a Jeans length of $2 r_\mathrm{acc}$, $r_\mathrm{acc}$ = `accretion_radius_cells`
cells), no sink lies within $r_\mathrm{acc}$, and its control volume (the cells within $r_\mathrm{acc}$, read through
the ghost cells) passes the checks of Federrath et al. (2010): the cell is a minimum of the potential (ties allowed:
the 8 cells around a symmetric centre all qualify), the flow converges along every principal axis (the symmetric
velocity gradient negative definite), the volume is Jeans-unstable ($|E_\mathrm{grav}| > 2E_\mathrm{th}$, $E_\mathrm{grav}$
summed over the pairs of its cells) and bound ($E_\mathrm{grav} + E_\mathrm{th} + E_\mathrm{kin} < 0$, the kinetic energy
relative to its centre of mass). The candidates of every rank are gathered and taken deepest potential first, one per
accretion radius, so every rank forms the same sinks. A static box above the threshold is not enough: at $\gamma = 1.4$
a control volume at the threshold is Jeans-stable, $|E_\mathrm{grav}|/2E_\mathrm{th} \approx \pi^2\gamma(\gamma-1)/10 =
0.55$. Gas that compresses adiabatically with $\gamma > 4/3$ never crosses the threshold ($\rho_\mathrm{thr} \propto
c_s^2$ grows with $\rho$ at $\gamma = 2$): collapse needs cooling ([radiative cooling](#radiative-cooling)).

**Accretion (D-M12-4).** A new sink has no mass; every sink then takes, once per step, the gas above the threshold from
the cells within its radius that is bound to it and falls in ($\tfrac12|u - v_s|^2 - Gm_s/r < 0$, $(u - v_s)\cdot(x -
x_s) < 0$), the cell that holds it and every cell of a new sink without the checks; a cell within several radii gives
its gas to the most bound sink. The cell keeps its specific values (every conservative variable scaled by
$\rho_\mathrm{thr}/\rho$: velocity, specific energy, passive scalars; the dual-energy entropy rebuilt from the
energy). The sink gains the mass and the momentum and moves to the centre of mass, from MPI-reduced sums: mass and
momentum are conserved exactly by the accretion; the gravity exchange is not, near a point mass the gas's face
gradients and the sink's CIC gather differ by about 1 % ([#94](https://github.com/szaghi/adam/issues/94), leg si5).

**Backends.** CPU and FNL run the same arithmetic: the deposit, the gather and the kick-drift-kick are host work on
both (FNL solves the Poisson equation on the host); accretion runs on the cells within the radii, which FNL gathers from
the device and scatters back (a handful of cells per sink); creation reads the whole committed state, which FNL copies
to the host every step when `creation` is on (correct, not scalable: a device screening of the candidates is a
follow-up). The sink histories of the two backends agree to round-off.

**Mergers.** Two sinks within an accretion radius that are bound and converging ($\tfrac12\mu|\Delta v|^2 - Gm_1m_2/d
< 0$, $\Delta v\cdot\Delta x < 0$) merge at their centre of mass, with their mass, momentum and mass-weighted
acceleration (the pair's mutual pull cancels in it).

**Runtime AMR (D-M12-5).** With runtime regridding (`[amr] frequency > 0`) the sinks live on the finest level allowed
(`[amr] max_level`), as in Federrath et al. (2010): they form only there, $r_\mathrm{acc}$ is in cells of that level for
the whole run, and at every regrid a block within $r_\mathrm{acc} + \tfrac12\,$`frequency`$ + 2$ cells of a sink (its
accretion radius, the half cell a sink may move per step until the next regrid, its CIC stencil) is refined towards that
level and never derefined, whatever the other markers vote. A sink formed on a coarser level would have its
neighbourhood refined at the next regrid: the infall onto it prolongated (supersonic, cold: the children lost their
pressure on the SI collapse) and its threshold quadrupled from one step to the next. A collapse on runtime AMR needs
`regrid_prolongation = conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)). A sink at rest on a composite
grid that is not symmetric about it feels its own potential through the 2:1 seams: a constant truncation force, $3\cdot10^{-2}$
finest cells over $t = 0.5$ on the SI si10 box, against $1.4\cdot10^{-4}$ on the uniform si1 box.

## Dissipative fluxes (Navier–Stokes)

With a viscosity or a conductivity in `[physics]` ([input](./input#dissipative-terms-issue-65-m4)), Euler becomes the
compressible Navier–Stokes system (issue #65, P2), and so does each MHD model (P3), which also takes an Ohmic
resistivity. Each term is the divergence of a flux, added to the inviscid one in direction $d$:

$$
G_{d}^{\rho u_i} = -\tau_{di},\qquad
G_{d}^{E} = -u_i\,\tau_{di} - k\,\frac{\partial T}{\partial x_d},\qquad
\tau = \mu\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T - \tfrac{2}{3}(\nabla\cdot\mathbf{u})\,I\right),
$$

with no bulk viscosity (Stokes hypothesis), $T = p/(\rho R)$, and the laws $\mu = \mu_0 (T/T_\mathrm{ref})^\omega$,
$k = k_0 (T/T_\mathrm{ref})^{\omega_k}$ ($\omega_k = \omega$ when $k$ comes from the Prandtl number, otherwise 0; the
constant law is $\omega = 0$). The pointwise flux is one pure routine shared by both backends
(`adam_flume_dissipation_library`).

**Ohmic resistivity (MHD).** With the magnetic diffusivity $\eta$ and $\mathbf J = \nabla\times\mathbf B$ (code units),
Ohm's law $\mathbf E = -\mathbf u\times\mathbf B + \eta\mathbf J$ in $\partial_t\mathbf B = -\nabla\times\mathbf E$ and
the Poynting flux $\mathbf E\times\mathbf B$ give

$$
G_{d}^{B_i} = -\eta\left(\frac{\partial B_i}{\partial x_d} - \frac{\partial B_d}{\partial x_i}\right),\qquad
G_{d}^{E} = \eta\,(\mathbf J\times\mathbf B)_d = \eta\,B_j\left(\frac{\partial B_d}{\partial x_j} - \frac{\partial B_j}{\partial x_d}\right),
$$

with $G_d^E = \sum_i B_i\,G_d^{B_i}$, so the field's energy loss is exactly the Ohmic heating
$\eta|\mathbf J|^2$ that the internal energy receives. The field flux is the **curl form**, $\epsilon_{idk}\,\eta J_k$,
not $-\eta\,\partial_d B_i$: the two differ by $\eta\nabla(\nabla\cdot\mathbf B)$, which is not zero under GLM, and the
curl form keeps the resistive term out of the divergence constraint. (The plan in #65 printed the field flux with
the opposite sign, an anti-diffusion; VV-5 measures the decay rate $\eta k^2$.) The Ohmic fluxes have their own
kernels, with the same stencils on $W = \mathbf B$, called only when the resistivity is on; the viscous and heat kernels
serve the MHD models unchanged ($(u, v, w, T)$ sit at the same auxiliary indexes). The psi equation and the EGLM
sources are untouched.

**Where it enters.** The dissipative face fluxes are computed by their own kernels, after the inviscid face fluxes and
**before** the seam accumulation and the flux difference, in both space schemes. Two consequences: the AMR reflux and
the inter-realm flux register see the total flux, so conservation across seams needs no new machinery; and the WENO
face kernels are untouched (their FNL register count does not change). A run without coefficients never calls the
kernels and stays bitwise identical to the ideal solver. The kernels read $(u, v, w, T)$ from the auxiliary variables,
whose indexes are common to every model, and the host picks the kernel of the order (`[numerics]
dissipative_order`), never a branch inside a loop.

**Face derivatives.** The scheme stores point values, so a face flux of order $p$ does not by itself give a flux
difference of order $p$: the conservative numerical flux is $h = G - \tfrac{\Delta x^2}{24}G'' + O(\Delta x^4)$
(Shu & Osher 1989). Face $i+\tfrac12$ along $d$, cells $m = -1 \dots 2$ counted from the left one:

| | `dissipative_order = 2` | `dissipative_order = 4` (default) |
|---|---|---|
| face state | $(W_0 + W_1)/2$ | $\big(9(W_0 + W_1) - (W_{-1} + W_2)\big)/16$ |
| normal derivative | $(W_1 - W_0)/\Delta x$ | $\big(27(W_1 - W_0) - (W_2 - W_{-1})\big)/(24\,\Delta x)$ |
| tangential derivative | mean of the 2nd-order central differences of cells 0 and 1 | the 4th-order central differences $D_m$ of cells $-1 \dots 2$, interpolated as the state |
| flux | $G(W_f, \nabla W_f)$ | $G(W_f, \nabla W_f) - (G_2 - G_1 - G_0 + G_{-1})/48$ |

where $G_m$ is the flux of cell $m$ from 2nd-order central gradients, so the correction term is
$\Delta x^2 G''/24$ to second order. The same form serves the linear and the nonlinear terms ($u_i\tau_{di}$,
$\mu(T)$). Across a strong temperature jump the four-point face temperature undershoots, possibly below zero, where a
power law $(T/T_\mathrm{ref})^\omega$ is NaN (a wall 300 times hotter than the gas did it in the input-contract check
DC); where it is not positive, the face takes the mean of its two cells for the laws, which a smooth positive $T$
never triggers. For constant $\mu$ it is a 4th-order $u_{xx}$ on six points, not the compact five-point one; its Fourier
symbol is non-positive at every wavenumber ($-4.67\,\nu/\Delta x^2$ at the grid Nyquist mode against $-5.33$ for the compact
operator), so it damps every mode, and its error is $O(\theta^6)$ in the wavenumber $\theta$. The 4th-order fluxes read
cells $-2 \dots 3$ along $d$ and $\pm 2$ along the tangents, so they need `[grid] ngc` $\ge 3$ (fatal otherwise); the
tangential differences read the edge ghosts, whose fill is exact on every path (the ghost probe GP,
[verification](./verification/infrastructure/ghost-probe)). A null direction has a zero step, so its derivatives are zero exactly and
its ghosts are never read.

**Walls.** `wall-noslip` and `wall-isothermal` ([boundary conditions](./boundary-conditions)) are mirrors, second-order
at the wall: a 4th-order interior with them converges at order 2 overall (VV-3).

**Time step.** The explicit unsplit integration adds the diffusive limit to the cell sum of the CFL bound,

$$
\Delta t = \frac{\mathrm{CFL}}{\max_\text{cells}\Big[\sum_d \frac{|u_d| + a}{\Delta x_d} + 2\,\nu\sum_d\frac{1}{\Delta x_d^2}\Big]},
\qquad \nu = \max\Big(\tfrac{4}{3}\tfrac{\mu}{\rho},\ \gamma\tfrac{k}{\rho c_p},\ \eta\Big),
$$

so under AMR the diffusive part tightens four times per level, against twice for the hyperbolic one. With a power
law and an isothermal wall the laws in $\nu$ are evaluated at $\max(T, T_{w,\max}^2/T)$, the geometric-mirror ghost
temperature of a cell beside the hottest wall: a wall much hotter than the gas gives its face a coefficient above
every cell's, which the cell bound alone misses (a wall 330 times hotter than the gas blew up at CFL 0.5). The first step
logs both maxima, the minimum cell Reynolds number $\sum_d (|u_d|+a)/\Delta x_d \big/ (\nu\sum_d \Delta x_d^{-2})$
($(|u|+a)\Delta x/\nu$ on a uniform 1-D grid), and whether the run is hyperbolic- or diffusion-limited.

**Limits.** Explicit only: a diffusion-limited run pays $\Delta t \propto \Delta x^2$ (super-time-stepping is a
follow-up). The positivity limiter is refused with dissipative terms (D-M4-5), and so are immersed solids, whose walls
are inviscid. A null direction freezes the momentum normal to it
(the CHASE semantics of the flux difference), so a reduced-dimension run cannot carry a velocity along a null axis: a
shear layer or a Couette flow needs its velocity axis active, if thin.

**Across seams.** The dissipative fluxes enter the AMR reflux and the inter-realm register with the inviscid ones, so a
2:1 face conserves to round-off (VV-7). A 2:1 face is second order, for them as for the inviscid fluxes: the mean
restriction of the coarse ghosts and the reflux of point-value fluxes each cap the order at 2 (issue #68 F1); without
the reflux the dissipative flux at the face is first order
([VV-8](./verification/dissipative/seam-accuracy)).

## Time integration

The library Runge–Kutta schemes (`[runge_kutta] scheme`, listed in the [input reference](./input)) integrate the
semi-discrete system: strong-stability-preserving (SSP) and low-storage schemes. Each stage refills the ghost cells,
recomputes the auxiliary variables, applies the positivity floors (MHD) and evaluates the residual. The time step follows
the [CFL bound](./models#time-step-bound), with the diffusive part of the
[dissipative fluxes](#dissipative-fluxes-navier–stokes) when they are on; `[time] it_max` and `time_max` stop the run.

The SSP Butcher tableaux have one source, `ssp_tableau` of `adam_rk_object` (the PRISM boundary and PIC integrators
call it too). SSP(5,4) is the optimal five-stage fourth-order scheme of Spiteri & Ruuth (2002, SIAM J. Numer. Anal.
40, 469), SSP coefficient 1.508: the library stores its Shu–Osher form as printed with 15 digits by Isherwood, Grant &
Gottlieb (2018, arXiv:1810.04800, §4), every coefficient nonnegative, and derives the Butcher tableau from it. The
order conditions then hold to $6\cdot10^{-16}$, and a refinement in time keeps fourth order down to $10^{-13}$. The
tableau the library used before ([#104](https://github.com/szaghi/adam/issues/104)) was accurate to ~$10^{-10}$ only:
its weights summed to $1 - 8.8\cdot10^{-11}$, so a constant force drifted by that relative amount per unit time and
an ODE refinement study stalled at $3\cdot10^{-11}$. The unit test `test_rk_ssp_tableau` pins every SSP tableau
([unit tests](./verification/unit-tests)).

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
  reports the ghosts blended, per rank. The isothermal and barotropic gases take $q_\text{in}$ outright
  ([#100](https://github.com/szaghi/adam/issues/100) P4): a ghost with $\rho \le 0$ makes the segment extrapolate the
  velocity $\rho\mathbf{u}/\rho$, up to $1/\kappa$ times the anchor's at the density floor. The internal-energy chord
  rejects that kinetic energy for the other models, but the density alone does not. On the Boss–Bodenheimer cloud the
  blended ghosts at the corners of the refined cube raised the velocity by a factor per step, until $dt \to 0$.
  Before the blend, the fill itself is bounded ([#101](https://github.com/szaghi/adam/issues/101)), with a
  monotonicity-preserving clamp in the manner of Suresh and Huynh: the interpolated density is clamped to the range of
  the coarse cells that bracket the ghost, widened by the curvature of the data ($\tfrac{3}{32}$ of the M4 minmod of the
  second differences, the largest over the bracketing lines; at the edge of the donor block, the widening of a
  quarter-cell extrapolation), and so is the velocity $\rho\mathbf{u}/\rho$, against the coarse cells' velocities; the
  momentum is that density times that velocity. A ghost whose clamped density is still below a tenth of the coarse
  cell's that holds it (`SEAM_DENSITY_GUARD`) takes, in every variable, that cell's state (the injection value). The
  energy and the other variables are interpolated unchanged. On smooth data the clamp moves nothing (the order across a
  2:1 corner is unchanged, AO `euler-l6`); FLUME turns it on for every model (`maps%seam_density_var`, the density
  index, the momenta following it). Without it, at the corners of a refined region the tricubic overshoots of density
  and momentum across the envelope of the Boss–Bodenheimer cloud gave ghosts whose velocity the blend extrapolated, and
  the corner cells ran away for about 150 steps, to 23 times the sound speed. Clamping the density and the momenta
  apart is not enough (the corner cells still reached 2.9 times the speed around them): the velocity, their ratio, is
  what must be bounded;
- **conservation** is restored by the Berger–Colella **reflux**: the fluxes of every Runge–Kutta stage through the seam
  are accumulated in a flux register, weighted by the stage's coefficient in the step, and at the end of the step the
  coarse cells next to the seam are corrected by the difference between the fine fluxes and their own
  (`[numerics] reflux = .true.`; with `reflux_cadence = stage` the correction is applied at every stage instead). Under
  GLM $\mathbf{B}$ and $\psi$ are refluxed; under EGLM $\mathbf{B}$ is not, since the sources $-\mathbf u\,\nabla\cdot\mathbf B$
  make it a non-conserved integral ([G3b](./verification/seams-studies/eglm-reflux)). With reflux the volume integrals of a periodic box stay
  constant to round-off, without it they drift at $10^{-5}$ ([conservation example](./verification/euler/amr-conservation)).

With markers the tree is an octree (`ratio = 8`, and with a null axis `nk ≥ 4`) or a quadtree (`ratio = 4`, any `nk`,
`nk = 1` for a true 2-D run): a quadtree seam is 2:1 in x and y and 1:1 in z, and the ghost restriction, the coarse-fine
interpolation, the reflux and the limiter's seam synchronisation work per axis (issue #46, verified by
[MV-15](./verification/mhd/quadtree-amr)). For MHD, a cell-centred $\mathbf{B}$ cannot be kept divergence-free across a 2:1 face (Tóth & Roe 2002): the seams
inject a truncation-level $\nabla\cdot\mathbf{B}$ that GLM transports and damps; the div(B) history reports it
(`seam_max_divb`).

**Runtime regridding** ([#74](https://github.com/szaghi/adam/issues/74)). With `[amr]
frequency = n > 0`, after every $n$-th step and before its output, the realm regrids: up to `iters` sweeps, each one
refreshing the ghost cells, combining the markers into one set of flags (a block is refined when any marker asks for it,
coarsened only when all agree, never below the base level), and calling the library regrid (2:1 balance, prolongation,
restriction, redistribution, maps). The distance function of the immersed solids is recomputed on the new grid, the new
state is checked (density and pressure positive, fatal otherwise, with the cell), and the forest rebuilds the 2:1 flux
register, so the reflux follows seams that appear and disappear. The time step stays global (no subcycling). The
FNL backend regrids by a host round trip: each sweep refreshes the ghosts on the device, copies the state to the host,
marks and regrids there with the same code, and copies the state, the block coordinates, the maps and the distance
function back (a map the new grid no longer has is freed on the device); its log line adds the wall time. The state
and distance-function copies move the blocks in use only: the device array leads with the block index, so a kernel
packs those blocks, transposed to the host layout, into a scratch buffer that one contiguous copy moves
([#75](https://github.com/szaghi/adam/issues/75); copying the whole capacity cost about 2.5 s per regrid on the
development box, two ranks of 17195 blocks, whatever the grid, against 0.02 s now on a 52-block grid). The
**Löhner marker** (Löhner 1987, mode 4) is the scale-free second-derivative estimator

$$
E = \sqrt{\frac{\sum_d \left(u_{+} - 2u + u_{-}\right)^2}{\sum_d \left(|u_{+} - u| + |u - u_{-}| +
\epsilon\,(|u_{+}| + 2|u| + |u_{-}|) + \delta\right)^2}} \in [0, 1],
$$

in the dimension-sum form of FLASH and PLUTO ($u_\pm$ the neighbours along $d$): about 0 where the variable is smooth
on the grid, about 1 at a jump. The `epsilon` term filters ripples small relative to the variable; on a variable that
vanishes in quiet regions it filters nothing, and round-off noise reads $E \approx 1$ (the GLM residue of
$|\mathbf{B}| \sim 10^{-6}$ outside a field loop refined the whole box). The absolute `floor` $\delta$ (default 0, the
FLASH/PLUTO form) makes the variations well below it noise. A block takes the largest $E$ over its cells and `buffer` ghost layers; it is refined
above `refine_tol`, coarsened below `derefine_tol`, kept between.
The **Jeans marker** (mode 5, [#91](https://github.com/szaghi/adam/issues/91)) refines where a self-gravitating cell
would be under-resolved: on every cell the Jeans length $\lambda_J = \sqrt{\pi c_s^2/(G\rho)}$, with the thermal
sound speed $c_s^2 = \gamma p/\rho$ (magnetic support not counted, so MHD refines at least as early), is compared with
the block spacing $h$. A block is refined where some cell has $\lambda_J/h$ below `jeans_cells` (default 8; Truelove
et al. 1997 ask for at least 4, against artificial fragmentation), coarsened where every cell is above
`jeans_cells_derefine` (default twice that, the hysteresis band of Federrath et al. 2010), kept between. It reads no
ghost cell, so it also shapes the initial grid: a star refined before its collapse is never prolonged mid-collapse
([SV-16](./verification/self-gravity/polytrope)). It needs self-gravity.

**Keep shocks off the 2:1 seams.** A captured shock that crosses a 2:1 seam into finer cells leaves a defect behind
it (about 2 % of the jump in density, with an entropy defect of the same size) that does not shrink with the resolution and that no
choice of seam operator removes: restriction, coarse-to-fine fill and reflux move it by at most ±30 %
([A3, E0](./verification/seams-studies/capturing-error)). The evidence points to the start-up error of a
captured shock whose numerical structure is squeezed into finer cells. Refine ahead of the shocks instead: a gradient or Löhner
marker on density or pressure with a `buffer` of a few cells and a regrid `frequency` short enough that a shock cannot
outrun the refined region between two regrids, so shocks stay inside the finest level and never cross a seam (the
standard AMR practice; E0 measured the defect, not this remedy).

**Regridding transfers**. When a block is
derefined, the parent takes the mean of its children. When one is refined, `[amr] regrid_prolongation` chooses how the
children are filled. FLUME's default is `conservative`: each child is the parent plus monotonized-central slopes
(van Leer 1977) times a quarter of a cell, one per refined axis, scaled by the largest $\phi \le 1$ that keeps every
child within the parent's and its face neighbours' range (Barth & Jespersen 1989). The children average to the parent,
so a regrid changes no conserved integral and a refine followed by a derefine returns the state exactly; a positive
variable stays positive; linear data are reproduced. The library's `linear` prolongation (tensor interpolation, weights
1/4 and 3/4) is second order too but not conservative ([RG](./verification/runtime-amr/regrid-round-trip)).
The opt-in `conservative-state` ([#88](https://github.com/szaghi/adam/issues/88)) keeps the children's mean on the
parent but builds them from that tensor interpolation, which mixes whole coarse states with the same weights for every
variable, and scales their offsets by one factor $\theta \le 1$ for the whole state, the largest that keeps every
variable of every child within the range of the parent and its $3^d$ neighbours. It is for collapsing self-gravitating
flows: after a regrid of the runtime-AMR collapse of an under-resolved star, with the default's children (limited
variable by variable) the coldest cells keep cooling, a supersonic infall develops and the collapse drains one to a
negative pressure, while with `conservative-state` they recover and the run ends ([SV-15](./verification/self-gravity/polytrope)). It is not the default because it costs accuracy on smooth
flow: an extremum of any one variable flattens the children of every variable, and the tracked vortex of AV misses its
bound with it (L1 1.68 times the uniform fine run's, the bound 1.5). All three read the parent's first ghost layer.

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

The solid surface can drive the AMR (the refined ring of the [shock-over-a-cylinder example](./verification/euler/shock-cylinder)).
`mhd-ideal` with solids is refused.

## Multi-realm runs

A forest manifest glues several realms, each with its own input file, through inter-realm seams (mirror coupling of
blocks of the same resolution). The seam ghosts are filled from the peer realm at every Runge–Kutta stage
(`stage_coincident`) or once per step (`end_of_step`); the union of a split problem is bitwise equal to the
single-realm run (verification `multirealm/`).

## Output and diagnostics

Checkpoints are XH5F files (one HDF5 file per rank plus an XDMF description readable by ParaView), written every
`it_save` steps, with optional slices, residual and conservation histories, and restart files. MHD adds the div(B)
history `<basename>-divb_history.dat` (`it time max_divb l1_divb seam_max_divb`, plus the D_flux columns with
`[diagnostics] divb_flux`) and the `divb_tol` monitor, on the normalised D_flux `max h|D_flux| / max|B|` (#68 G1).

## References

- Barth T. J., Jespersen D. C. (1989), The design and application of upwind schemes on unstructured meshes,
  AIAA paper 89-0366.
- Batten P. et al. (1997), On the choice of wavespeeds for the HLLC Riemann solver, *SIAM J. Sci. Comput.* 18, 1553–1570.
- Berger M. J., Colella P. (1989), Local adaptive mesh refinement for shock hydrodynamics, *J. Comput. Phys.* 82, 64–84.
- Bryan G. L. et al. (2014), ENZO: an adaptive mesh refinement code for astrophysics, *ApJS* 211, 19,
  doi:10.1088/0067-0049/211/2/19.
- Chen Y., Tóth G., Gombosi T. I. (2016), A fifth-order finite difference scheme for hyperbolic equations on
  block-adaptive curvilinear grids, *J. Comput. Phys.* 305, 604–621.
- Christlieb A. J., Liu Y., Tang Q., Xu Z. (2015), Positivity-preserving finite difference weighted ENO schemes with
  constrained transport for ideal magnetohydrodynamic equations, *SIAM J. Sci. Comput.* 37, A1825–A1845,
  doi:10.1137/140971208.
- Einfeldt B. et al. (1991), On Godunov-type methods near low densities, *J. Comput. Phys.* 92, 273–295.
- Federrath C., Banerjee R., Clark P. C., Klessen R. S. (2010), Modeling collapse and accretion in turbulent gas
  clouds: implementation and comparison of sink particles in AMR and SPH, *ApJ* 713, 269–290,
  doi:10.1088/0004-637X/713/1/269.
- Don W. S., Li R., Wang B.-S., Wang Y. H. (2022), A novel and robust scale-invariant WENO scheme for hyperbolic
  conservation laws, *J. Comput. Phys.* 448, 110724.
- Hockney R. W., Eastwood J. W. (1988), *Computer Simulation Using Particles*, Adam Hilger, Bristol.
- Jiang G.-S., Shu C.-W. (1996), Efficient implementation of weighted ENO schemes, *J. Comput. Phys.* 126, 202–228.
- Becker R. (1922), Stoßwelle und Detonation, *Z. Phys.* 8, 321–362.
- Löhner R. (1987), An adaptive finite element scheme for transient problems in CFD, *Comput. Methods Appl. Mech.
  Engrg.* 61, 323–338.
- Miyoshi T., Kusano K. (2005), A multi-state HLL approximate Riemann solver for ideal MHD, *J. Comput. Phys.* 208, 315–344.
- Shu C.-W., Osher S. (1989), Efficient implementation of essentially non-oscillatory shock-capturing schemes II,
  *J. Comput. Phys.* 83, 32–78.
- Toro E. F., Spruce M., Speares W. (1994), Restoration of the contact surface in the HLL-Riemann solver,
  *Shock Waves* 4, 25–34.
- Townsend R. H. D. (2009), An exact integration scheme for radiative cooling in hydrodynamical simulations, *ApJS*
  181, 391–397, doi:10.1088/0067-0049/181/2/391.
- Truelove J. K., Klein R. I., McKee C. F., Holliman J. H. II, Howell L. H., Greenough J. A. (1997), The Jeans
  condition: a new constraint on spatial resolution in simulations of isothermal self-gravitational
  hydrodynamics, *ApJ* 489, L179–L183, doi:10.1086/310975.
- Tóth G., Roe P. L. (2002), Divergence- and curl-preserving prolongation and restriction formulas,
  *J. Comput. Phys.* 180, 736–750.
- Villasenor B., Robertson B., Madau P., Schneider E. (2021), Effects of photoionization and photoheating on
  Lyman-α forest properties from Cholla cosmological simulations, *ApJ* 912, 138, doi:10.3847/1538-4357/abed5a.
- van Leer B. (1977), Towards the ultimate conservative difference scheme. IV. A new approach to numerical
  convection, *J. Comput. Phys.* 23, 276–299.
- Wu K. (2018), doi:10.1137/18M1168017; Wu K., Shu C.-W. (2018), doi:10.1137/18M1168042: positivity of the first-order
  and high-order schemes for ideal MHD with the Godunov–Powell source.
- Xu Z. (2014), Parametrized maximum principle preserving flux limiters for high order schemes solving hyperbolic
  conservation laws, *Math. Comp.* 83, 2213–2238, doi:10.1090/S0025-5718-2013-02788-3.
- Zhang X., Shu C.-W. (2012), Positivity-preserving high order finite difference WENO schemes for compressible Euler
  equations, *J. Comput. Phys.* 231, 2245–2258, doi:10.1016/j.jcp.2011.11.020.
