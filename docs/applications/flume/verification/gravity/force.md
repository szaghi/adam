# The gravity force

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GV-1, GV-2, GV-3, GV-6, GV-7 (G1) | External gravity | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes | — |

## Case and scope

G1 adds the gravity sources: momentum $\rho\mathbf g$, and the energy from the face mass fluxes. Three properties must
hold by construction and are checked to round-off; a fourth is about robustness:

- **free fall is exact** (GV-1): a uniform gas in a periodic box falls freely, $\rho\mathbf u = \rho\mathbf g t$,
  $E = E_0 + \rho g^2 t^2/2$, $\rho$ unchanged. A momentum or energy source with a wrong weight fails it;
- **a closed box conserves its mass and total energy** (GV-3): the gas falls onto the floor and sloshes; mass and
  $E + \rho\Phi$ ($\Phi = -\mathbf g\cdot\mathbf x$; since G2a the invariant $W$, see below) stay constant while $E$
  alone exchanges with $\rho\Phi$. A cell-centred energy source $\rho\mathbf u\cdot\mathbf g$ would fail it;
- **scaling covariance** (GV-6): the run with the scale-invariant WENO weights is the same run in other units, bit for
  bit, as #49 established without gravity;
- **the positivity limiter keeps the run alive** where the gas falls toward vacuum (`box2d-lim`, `box2d-stiff`).

GV-2 reports the spurious flow of the unbalanced scheme on an atmosphere at rest, the baseline of
[the well-balanced mode](./well-balanced); GV-7 compares the backends.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$; `box2d-mhd`: ideal MHD with EGLM (`glm_ch = 3`, `glm_alpha = 0.18`, `glm_damping_length = 1`) | [models](../../models#external-gravity), [EGLM](../../models#divergence-control-eglm) |
| Momentum source | $\rho_i\mathbf g$ at the cell centre | `add_gravity_sources`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `add_gravity_sources_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |
| Energy source | G1: the face form $\tfrac12\sum_d g_d(F^\rho_{d,i-1/2} + F^\rho_{d,i+1/2})$ from the mass fluxes of the continuity equation (Mikami et al. 2008, Mullen et al. 2021); since G2a plus its fourth-order correction ([energy order](./energy-order)) | `add_gravity_sources`, `compute_gravity_correction`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Positivity limiter (option B) | `positivity_limiter = cell` (`box2d-lim`, `box2d-stiff`, `box2d-mhd`): the backbone takes $\rho\mathbf g$ and the time-centred energy $\rho\mathbf g\cdot(\mathbf u + \mathbf g\Delta t/2)$; the face-form excess is limited by $\Lambda$ | `compute_positivity_factors`, [`adam_flume_cpu_positivity_kernels_agnostic.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_positivity_kernels_agnostic.INC); `add_gravity_limited_correction`, `limit_positivity`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#positivity-limiter) |
| Input checks | a gravity component along a null direction and immersed solids are refused; the stiffness $|\mathbf g|\Delta t^2/(\mathrm{CFL}\,h_{min})$ is logged at the first step | `check_gravity`, `log_gravity_regime`, [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`); GV-6 with `weights = si` | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#weno-reconstruction) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | `periodic`; `wall-inviscid` (mirrored ghosts); `hydrostatic` (GV-2) | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); [boundary conditions](../../boundary-conditions#wall-inviscid-slip-wall-by-mirroring-with-a-sign-vector) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90), [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |

## Initial and boundary conditions

Common to every leg: $[0, 1]^2$, quadtree, 4×4 blocks of $8^2$ cells ($N = 32$) unless stated, two uniform levels,
`ngc = 4`, null z (`extrapolation`), np 2, CFL 0.4, $\gamma = 1.4$. The G1 cases start from the
[`uniform`](../../initial-conditions#uniform-seeded-uniform-state) state $\rho = 1$, $p = 1$ at rest under
$\mathbf g = (0, -1)$.

| Leg | Case | State | x faces | y faces | Run |
|---|---|---|---|---|---|
| GV-1 | `freefall2d` | uniform | `periodic` | `periodic` | 50 steps |
| GV-3 | `box2d` | uniform | `periodic` | `wall-inviscid` | 200 steps |
| GV-3 | `box2d-mhd` | uniform, EGLM, $\mathbf B = (0.5, 0, 0)$, limiter on | `periodic` | `wall-inviscid` | 200 steps |
| limiter | `box2d-lim` | $p = 0.05$, limiter on | `periodic` | `wall-inviscid` | 300 steps |
| limiter | `box2d-stiff` | $p = 10^{-3}$, limiter on ($g\Delta t \approx 4.5\,c$) | `periodic` | `wall-inviscid` | 200 steps |
| GV-2 | `iso2d` with the force | [isothermal atmosphere](../../initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity), $\rho_1 = p_1 = 1$ at the origin | `periodic` | `hydrostatic` | to $t = 0.5$, $N = 32, 64, 128$ (8, 16, 32 cells per block) |
| GV-6 | `box2d`, `weights = si` | uniform | `periodic` | `wall-inviscid` | 100 steps, against its input rescaled: lengths by $2^2$, velocities by $2^{-1}$, density by $4^{-2}$ |
| GV-7 | `box2d`, `box2d-mhd` on FNL | as above | | | as above |

Checkpoints every 10 steps (GV-3 asserts every one of them).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes each case (`freefall2d`, `box2d`, `box2d-mhd`, `box2d-lim`, `box2d-stiff`, `iso2d`; `--cells`, `--cfl`, `--time-max`, `--weights`, `--steps`) |
| `src/tests/flume/verification/gravity/check.sh` | runs the legs (`--legs LIST`, e.g. `gv1,gv3,limiter,gv2`; the force legs unless `--no-force`; `--scaling` adds `gv6` (GV-6); `--agree TAG` for GV-7, which needs the legs `gv3`, `gv8`, `gv10`, `gv4c`; `--keep`) |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | `freefall` (GV-1, $10^{-12}$), `energy` (GV-3, $10^{-12}$; `--report-total` for the limiter cases), `rest` (GV-2), `agree` (GV-7, $10^{-10}$) |
| `src/tests/flume/verification/scaling/scaling.py` | GV-6: `rescale --j 2 --k -1 --m -2` and `compare` |
| `src/tests/flume/verification/make_doc_runs.sh` (group `gravity`) | re-runs `freefall2d`, `box2d`, `box2d-mhd` with their checkpoints |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_gravity_force`) | draws the figure |

```bash
cd src/tests/flume/verification/gravity
./check.sh --np 2 --cases iso2d --scaling --keep                         # CPU: G0 iso2d, the force legs, GV-6
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np 2 --cases iso2d --keep \
   --agree adam_flume_cpu                                                # FNL, and GV-7 against the CPU runs
```

## Results

![GV (G1): free fall, the energy budget of the closed boxes and the settled box2d](/flume/gravity-force.png)

*Left, GV-1 free fall; centre, the GV-3 energy budget over the whole run; right, `box2d` settled after 200 steps.*

Drawn with the current scheme (G2a's fourth-order energy source). Left, GV-1: the momentum of the periodic box follows
$\rho g t$, with a pointwise error of one or two units in the last place. Centre, GV-3 over the whole run: the mass and
$W$ stay at round-off, while the plain $E + \rho\Phi$ drifts to $3\cdot10^{-5}$ (`box2d`) through the $O(\Delta x^2)$
wall-row term of G2a; $E$ alone (not drawn) exchanges 3.8 per cent with $\rho\Phi$. Right, `box2d` after 200 steps:
the uniform gas has settled into a stratification, $\rho$ from 1.58 at the floor to 0.42 under the lid.

As measured at G1 (the face-form energy source; CPU, FNL the same by GV-7):

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| GV-1 | $\rho$, $\rho\mathbf u$, $E$ against the exact free fall | $10^{-12}$ | $\rho$ exact, $\rho\mathbf u$ $2.2\cdot10^{-16}$, $E$ $5.6\cdot10^{-16}$ | the same |
| GV-3 `box2d` | mass and $E + \rho\Phi$ constant at every checkpoint | $10^{-12}$ | mass $1.1\cdot10^{-14}$, $E + \rho\Phi$ $1.3\cdot10^{-14}$ ($E$ alone moves $3.8\cdot10^{-2}$) | the same |
| GV-3 `box2d-mhd` | as above | $10^{-12}$ | $8.4\cdot10^{-15}$, $7.5\cdot10^{-15}$ ($E$ $1.8\cdot10^{-2}$) | the same |
| `box2d-lim` | completes, mass conserved; $E + \rho\Phi$ reported | $10^{-12}$ (mass) | 74 304 faces limited, no inadmissible backbone, mass $1.0\cdot10^{-14}$, $E + \rho\Phi$ $4.4\cdot10^{-7}$ | the same |
| `box2d-stiff` | as above | $10^{-12}$ (mass) | 324 424 faces limited, 32 inadmissible backbones, mass $1.2\cdot10^{-14}$, $E + \rho\Phi$ $1.6\cdot10^{-2}$ | the same |
| GV-2 | $\max|\mathbf u|/c$ at $N = 32, 64, 128$, reported | — | $2.9\cdot10^{-12}$, $1.1\cdot10^{-12}$, $2.5\cdot10^{-14}$ | the same, but $2.6\cdot10^{-14}$ at $N = 128$ |
| GV-6 | `weights = si` against the rescaled input | bitwise | bitwise | bitwise |
| GV-7 | FNL against CPU, per field | $10^{-10}$ | — | `box2d` $\le 9.3\cdot10^{-14}$, `box2d-mhd` $\le 9.3\cdot10^{-14}$ |

GV-7 scales a vector component by its vector's largest component: $\rho u$ is zero by symmetry and holds round-off
only. GV-2 runs a mild stratification (scale height 1 on the unit box).

Re-measured with the fourth-order energy source of G2a (the current scheme; CPU and FNL):

| Leg | Asserted | Bound | CPU and FNL |
|---|---|---|---|
| GV-1 | as above | $10^{-12}$ | unchanged (uniform state: the correction is zero) |
| GV-3 `box2d` | mass and $W$ | $10^{-12}$ | mass $7.7\cdot10^{-15}$, $W$ $7.7\cdot10^{-15}$ (the plain $E + \rho\Phi$ $3.3\cdot10^{-5}$, the wall term) |
| GV-3 `box2d-mhd` | mass and $W$ | $10^{-12}$ | $1.2\cdot10^{-14}$, $6.1\cdot10^{-15}$ (plain sum $1.6\cdot10^{-5}$) |
| `box2d-lim` | completes, mass conserved; $W$ reported | $10^{-12}$ (mass) | 107 984 faces limited (the correction is a high-order cell source inside $\Lambda$), no inadmissible backbone, mass $9.4\cdot10^{-15}$, $W$ $1.2\cdot10^{-6}$ (limited cells only) |
| `box2d-stiff` | as above | $10^{-12}$ (mass) | runs, mass $1.2\cdot10^{-14}$ |
| GV-2 | reported | — | unchanged |
| GV-6, NV-5 | bitwise | bitwise | bitwise |
| GV-7 | FNL against CPU, per field | $10^{-10}$ | $\le 7.9\cdot10^{-14}$ |

Related cases: [the order of the energy source](./energy-order) (G2a, the correction and the invariant $W$);
[the well-balanced mode](./well-balanced) removes the GV-2 flow; [positivity limiter](../mhd/positivity-limiter) and
[scaling covariance](../infrastructure/scaling-covariance) without gravity; the same box through the potential path of
self-gravity in [SV-1](../self-gravity/potential).

## Findings

- **The total energy.** With the face-form energy source the conserved total is $E + \rho\Phi$ to round-off at
  every checkpoint, while $E$ alone moves by per cent. Since G2a the exactly conserved total at reflecting walls is
  $W$, $E + \rho\Phi$ plus an $O(\Delta x^2)$ wall-row term ([energy order](./energy-order)); GV-3 asserts $W$ and
  reports the plain sum.
- **The limiter and gravity (owner decision: option B).** The first form put the face-form energy into the backbone.
  A cell beside a wall then gets half of the face form (the wall face carries no mass) while its momentum gains the
  whole $\rho\mathbf g\,\Delta t$: the stage takes about $\tfrac12\rho|\mathbf u||\mathbf g|\Delta t$ out of its
  internal energy, which near vacuum is negative. `box2d-lim` died at step 128 and `box2d-stiff` at step 1 (every
  backbone inadmissible), with and without the limiter. Option B gives the backbone the time-centred energy
  $\rho\mathbf g\cdot(\mathbf u + \mathbf g\Delta t/2)$, which leaves its internal energy unchanged, and limits the
  face-form excess with $\Lambda$; both cases then run. The stiff case's 32 inadmissible backbones come from the
  fluxes, not from gravity.
- In the limited cases $E + \rho\Phi$ (since G2a $W$) is exact only where the limiter does not act; it is reported,
  not asserted.
- `box2d-stiff` is the stiff regime: $g\Delta t$ is about 4.5 times the sound speed. No time-step bound is imposed for
  gravity; the regime is logged at the first step.
- GV-2 is the unbalanced baseline of [the well-balanced mode](./well-balanced): the spurious Mach number of an
  atmosphere at rest falls with $N$ but is not zero.
