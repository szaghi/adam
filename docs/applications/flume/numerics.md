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
| `weno-riemann` | WENO interpolation of the face states, a Riemann solver and a high-order correction (Chen, Tóth & Gombosi 2016) | Euler (LLF, HLL, HLLC); MHD in development (issue #47) |

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
   (`reconstruction_variables`: `characteristic`, the Euler default, with the Roe eigenvectors of the face; or
   `primitive`, $(\rho, u, v, w, p)$). A face state with non-positive density or pressure is replaced by the adjacent
   cell's state.
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

The MHD solvers (LLF, HLL, and HLLD of Miyoshi & Kusano 2005, with the $(B_n, \psi)$ subsystem solved exactly and an HLL
fallback) are implemented and unit-tested; their wiring into the scheme is in progress (issue #47, M3-P3).

**Measured** (np 2, WENO5, 6th-order correction, CPU and FNL identical to the printed digits; details in the
[verification gallery](./verification)): on the isentropic vortex the HLLC variant reaches order 6.5 between 128² and
256² with $L_1(\rho) = 8.0\cdot10^{-8}$ (splitting: $8.4\cdot10^{-8}$); on Sod and Shu–Osher its $L_1$ error is 9% and 7%
below the splitting scheme's, on Lax 1.6% above.

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
  coarse ghosts are restricted from the fine cells;
- **conservation** is restored by the Berger–Colella **reflux**: the fluxes of every Runge–Kutta stage through the seam
  are accumulated in a flux register, weighted by the stage's coefficient in the step, and at the end of the step the
  coarse cells next to the seam are corrected by the difference between the fine fluxes and their own
  (`[numerics] reflux = .true.`; $\mathbf{B}$ and $\psi$ included). With reflux the volume integrals of a periodic box stay
  constant to round-off, without it they drift at $10^{-5}$ ([conservation example](./verification#conservation-across-amr-seams)).

With markers the tree must be an octree (`ratio = 8`, and with a null axis `nk ≥ 4`): a quadtree with markers is refused
(issue #46). For MHD, a cell-centred $\mathbf{B}$ cannot be kept divergence-free across a 2:1 face (Tóth & Roe 2002): the seams
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
   Runge–Kutta stages.

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
- Einfeldt B. et al. (1991), On Godunov-type methods near low densities, *J. Comput. Phys.* 92, 273–295.
- Jiang G.-S., Shu C.-W. (1996), Efficient implementation of weighted ENO schemes, *J. Comput. Phys.* 126, 202–228.
- Miyoshi T., Kusano K. (2005), A multi-state HLL approximate Riemann solver for ideal MHD, *J. Comput. Phys.* 208, 315–344.
- Shu C.-W., Osher S. (1989), Efficient implementation of essentially non-oscillatory shock-capturing schemes II,
  *J. Comput. Phys.* 83, 32–78.
- Tóth G., Roe P. L. (2002), Divergence- and curl-preserving prolongation and restriction formulas,
  *J. Comput. Phys.* 180, 736–750.
- Toro E. F., Spruce M., Speares W. (1994), Restoration of the contact surface in the HLL-Riemann solver,
  *Shock Waves* 4, 25–34.
