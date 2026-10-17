# Rayleigh–Taylor

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GV-4 (G3) | External gravity | M7 ([#73](https://github.com/szaghi/adam/issues/73)); symmetry study [#83](https://github.com/szaghi/adam/issues/83) | CPU, FNL | yes | — |

## Case and scope

Rayleigh–Taylor is the gravity test of the literature, and there it is compared by eye (Liska & Wendroff 2003; Stone
et al. 2008, §8.3; Almgren et al. 2010, §8.3): its nonlinear stage has no converged solution. Its linear stage has
one, the growth rate of the seeded mode, and GV-4 asserts it against the linearised compressible Euler equations of
the exact equilibrium profile (owner decision: B, with the nonlinear stage as figures, C). The incompressible
sharp-interface formula $\sigma^2 = Agk$ is not a reference for that: the finite interface width and compressibility
move $\sigma$ by several per cent.

What the case catches:
- **a wrong equilibrium or `hydrostatic` ghosts** of the two layers (GV-4a), and a spurious flow on the unseeded
  layers (GV-4b);
- **a wrong linear growth rate** (GV-4c): the fitted rate against the oracle's eigenvalue to $10^{-3}$, at two
  resolutions;
- **backend divergence** amplified by the instability (GV-7).

The nonlinear stage, hydro and MHD, is shown as figures without an oracle. The mirror symmetry of the MHD run is the
subject of a study ([#83](https://github.com/szaghi/adam/issues/83)), at the end of the Results: a symmetry that the
equations keep and the scheme loses is either physics (an unstable mode of the flow, fed by round-off) or numerics,
and a user who reads the asymmetry needs to know which.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$; `rt2d-mhd`: ideal MHD, EGLM (`glm_ch = 3`, `glm_alpha = 0.18`, `glm_damping_length = 1`, `glm_ch_check = warning`); the #83 runs also `glm` and `none` | [models](../../models#external-gravity), [EGLM](../../models#divergence-control-eglm), [GLM](../../models#divergence-control-mixed-glm); `set_glm_damping`, [`adam_flume_mhd_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_mhd_object.F90) |
| Equilibrium | two isothermal layers joined by a tanh in temperature at equal pressure, closed-form hydrostatic pressure; divergence-free single-mode seed with a Gaussian envelope | `atmosphere_state`, `set_initial_conditions`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| Gravity sources | momentum $\rho\mathbf g$, fourth-order energy source | `add_gravity_sources`, `compute_gravity_correction`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Well-balanced mode | GV-4b balanced leg only (`equilibrium = atmosphere`) | `refresh_equilibrium_residual`, `subtract_equilibrium_residual`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`); one #83 comparison with `weno-riemann` HLLD | `compute_residuals_weno`, `compute_residuals_riemann`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [`adam_flume_cpu_mhd_eglm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_eglm_kernels.F90), [`adam_flume_cpu_mhd_glm_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_glm_kernels.F90), [`adam_flume_cpu_mhd_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_mhd_kernels.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), CFL 0.4 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| Boundaries | x `periodic`, y `hydrostatic` (exact on isothermal layers) | `hydrostatic_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |
| Oracle | linearised compressible Euler about the closed-form profile, RK4 shooting with bisection on $\sigma$ | `growth_rate`, `src/tests/flume/verification/gravity/rt_oracle.py` |

**The equilibrium** ([`rayleigh-taylor`](../../initial-conditions#rayleigh-taylor-two-isothermal-layers-heavy-over-light)):
two isothermal layers, cold and dense over hot and light, joined by a tanh in temperature at equal pressure. Isothermal
layers are stably stratified, so the interface is the only instability, and the `hydrostatic` boundary is exact on
them. Two layers of uniform density would not do: under gravity a uniform-density layer has entropy falling with
height and is convectively unstable, a second instability competing with the first.

**The oracle** (`gravity/rt_oracle.py`). With perturbations $\propto e^{ikx + \sigma t}$, eliminating $u$ and
$\rho_1$, the adiabatic linear equations reduce to two first-order ODEs in $(v, p_1)$ along the height, with $v = 0$ at
both faces. RK4 shooting from the bottom with bisection on $\sigma$ gives the fundamental mode (no node in $v$). It is
independent of FLUME and checked three ways (`rt_oracle.py selftest`, $g = 0.1$, $k = 4\pi$, density ratio 2):

| Check | Bound | Result |
|---|---|---|
| closed-form pressure against quadrature of $dp/dh = -\rho g$ | $10^{-10}$ | $9.3\cdot10^{-13}$ |
| incompressible ($p_1 = 10^6$), sharp ($\delta$ = 0.004, 0.002) against the two-layer $\sigma^2 = gk\,\Delta\rho/(\rho_h\coth kH + \rho_l\coth kH)$ | $k\delta$ | $2.3\cdot10^{-2}$, $1.1\cdot10^{-2}$: $O(k\delta)$, halving with $\delta$ |
| light over heavy | no growing mode | no growing mode |
| 2000 against 4000 RK4 steps | $10^{-8}$ | $2.4\cdot10^{-11}$ |

## Initial and boundary conditions

| Item | `rt2d` (GV-4) | `rt2d-nl` | `rt2d-mhd` |
|---|---|---|---|
| Domain | $[0, 0.5]\times[0.25, 0.75]$: one wavelength; the mode is down to $e^{-\pi}$ at the y faces, which the oracle includes as walls | $[0, 1]^2$: two wavelengths | as `rt2d-nl` |
| Grid | quadtree, 4×4 blocks of $8^2$ ($N = 32$) and of $16^2$ ($N = 64$), two uniform levels, `ngc = 4` | 4×4 blocks of $32^2$ ($N = 128$) | as `rt2d-nl` |
| Ranks | 2 | 2 | 2 |
| Gravity | $\mathbf g = (0, -1)$ | as `rt2d` | as `rt2d` |
| Layers | density ratio 2:1, light layer $\rho_1 = 1$, $p_1 = 1$ at the interface height $h_0 = 0.5$ | as `rt2d` | as `rt2d` |
| Interface width $\delta$ | 0.05 | 0.01 | 0.01 |
| Seed | wavelength 0.5, amplitude $10^{-5}$, divergence-free, Gaussian envelope of width 0.1 about the interface | amplitude $10^{-2}$ | amplitude $10^{-2}$ |
| Field | — | — | $\mathbf B = (0.1, 0, 0)$, half the critical field of the mode ($\sqrt{\Delta\rho\,g/(2k)} = 0.20$ for a sharp incompressible interface); $\beta \approx 200$ |
| Linear rate (oracle) | $\sigma = 1.523265$ | $\sigma = 1.869$ | — |
| Boundaries | x `periodic`, y `hydrostatic` | as `rt2d` | as `rt2d` |
| CFL | 0.4 | 0.4 | 0.4 |
| Final time | GV-4a one step at CFL $10^{-30}$ (unseeded); GV-4b $t = 1$ (unseeded); GV-4c $t = 4$, checkpoints every $5n$ steps ($n$ cells per block side) | $t = 5$, checkpoint every 100 steps | as `rt2d-nl` |

The #83 runs are `rt2d-mhd` with one cleaning parameter changed at a time: `divergence_control = none`, `glm`; EGLM
with `glm_ch` = 1.5, 6; `glm_alpha` = 0, 1, 5; `glm_damping_length = min-cell`.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes `rt2d`, `rt2d-nl`, `rt2d-mhd` (`--seed`, `--cells`, `--it-save`, `--balanced`, `--divergence-control`, `--glm-ch`, `--glm-alpha`, `--glm-damping-length`) |
| `src/tests/flume/verification/gravity/check.sh` | runs GV-4a, GV-4b, GV-4c with the force legs (`--legs gv4a,gv4b,gv4c`); `--agree TAG` adds GV-7 on `rt-8` |
| `src/tests/flume/verification/gravity/rt_oracle.py` | the linear growth rate (`sigma <input.ini>`) and its self-test (`selftest`) |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | `profile` (GV-4a, $10^{-13}$), `rest` (GV-4b, `--max 1e-12`), `growth` (GV-4c, `--t0 2 --t1 4 --tol 1e-3`), `agree` (GV-7), `symmetry` (#83) |
| `src/tests/flume/verification/gravity/glm_linear.py` | #83: the spectrum of the linearised GLM and EGLM systems about a magnetised double shear layer |
| `src/tests/flume/verification/make_doc_runs.sh` (groups `gravity`, `rt`, `rt-sym`) | `rt-{8,16}`; `rt2d-nl`, `rt2d-mhd`; the eight `rt-sym-*` runs |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_gravity_rt_linear`, `fig_gravity_rt`, `fig_gravity_rt_symmetry`) | draws the three figures |

```bash
cd src/tests/flume/verification/gravity
./rt_oracle.py selftest                                                  # the oracle
./check.sh --np 2 --cases iso2d --keep                                   # CPU (GV-4 runs with the force legs)
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np 2 --cases iso2d --keep \
   --agree adam_flume_cpu                                                # FNL, and GV-7 on rt-8
cd .. && ./make_doc_runs.sh <dir> rt                                     # the nonlinear runs
FLUME_EXE=$PWD/../../../../exe/adam_flume_fnl ./make_doc_runs.sh <dir> rt-sym   # the #83 runs
```

## Results

### Linear stage (GV-4)

![GV-4: the seeded mode, its local growth rate and its vertical profile against the linear theory](/flume/gravity-rt-linear.png)

*The seeded mode against time, its local growth rate, and its vertical profile at $t = 4$ against the oracle's
eigenfunction ($N = 32$, 64).*

Left, the mode grows from $3.6\cdot10^{-6}$ to $7.9\cdot10^{-4}$, at the theory's slope once the transient has passed.
Centre, the local rate between checkpoints: it settles by $t \approx 1.5$ ($\sigma t \approx 2.3$) and then oscillates
around the theory by a few $10^{-3}$, the beat of the seed's stable components. Right, the vertical profile
$|v_k(y)|$ at $t = 4$ against the oracle's eigenfunction: $1.6\cdot10^{-3}$ ($N = 32$) and $1.3\cdot10^{-3}$
($N = 64$) apart at most. FLUME was never given that shape: the seed is a Gaussian.

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| GV-4a unseeded, one negligible step | the layers and their hydrostatic ghosts against the closed form | $10^{-13}$ | interior $1.1\cdot10^{-16}$, ghosts $2.2\cdot10^{-16}$ | the same |
| GV-4b unseeded to $t = 1$ | balanced $\max|\mathbf u|/c$; unbalanced reported | $< 10^{-12}$ | balanced 0; unbalanced $2.2\cdot10^{-6}$ | the same |
| GV-4c seeded, $N = 32$ and 64 to $t = 4$ | the slope of $\ln a$ over $t \in [2, 4]$ against $\sigma = 1.523265$ | $10^{-3}$ | 1.522543 ($4.7\cdot10^{-4}$), 1.522724 ($3.6\cdot10^{-4}$); the two agree to $1.2\cdot10^{-4}$ | the same |
| GV-7 FNL against CPU | `rt-8` ($N = 32$, seeded), per field | $10^{-10}$ | — | $\le 1.3\cdot10^{-12}$ ($\rho$): round-off amplified by the instability, $e^{\sigma t} \approx 400$ over the run |

$a$ is the amplitude of the seeded mode in the vertical velocity,
$a = (\sum_y |\tfrac2{N_x}\sum_x v\,e^{-ikx}|^2\Delta y)^{1/2}$.

### Nonlinear stage (figures, no oracle)

![Rayleigh-Taylor single mode, hydro and MHD, density at t = 2, 3, 4, 5](/flume/gravity-rt.png)

*Density of `rt2d-nl` (top, hydro) and `rt2d-mhd` (bottom, $B_x = 0.1$) at $t = 2, 3, 4, 5$, $N = 128$; each row has its
own colour scale.*

Hydro (top): the spikes roll up through secondary Kelvin–Helmholtz instabilities by $t = 3$ and reach the floor by
$t = 5$. The local growth rate of the seeded mode peaks at 1.77 at $t \approx 1.5$, 5 % under the linear 1.869, before
the mode saturates. This is not a measurement: the 0.01 interface spans 1.3 cells and the $10^{-2}$ seed is nonlinear
almost from the start, which is why GV-4 runs its own resolved case. MHD (bottom): the field slows the mode, with a
peak rate of 1.49, 0.80 of the hydro rate against the sharp-interface estimate $\sqrt{1 - (B/B_c)^2} = 0.87$, and
suppresses the roll-ups, as in Stone et al. (2008), §8.4.

Both runs are periodic with two identical wavelengths, each mirror-symmetric about its own axis. Hydro keeps both
symmetries to round-off, which grows with the flow to $6\cdot10^{-10}$ at $t = 5$. MHD does not. Its mirror asymmetry
stays at round-off to $t \approx 2$, then grows at a rate of 6 to 10, against the mode's 1.5: $2.3\cdot10^{-10}$ at
$t = 3$, $1.5\cdot10^{-7}$ at $t = 4$, $2.8\cdot10^{-3}$ at $t = 5$ (invisible in the figure). `weno-riemann` with HLLD
delays it by about one time unit, $9.8\cdot10^{-12}$ at $t = 3$ and $2.1\cdot10^{-5}$ at $t = 5$, but it grows as fast.
The peak growth rates of the mode agree (1.4931 and 1.4914). It is therefore not the low-$\beta$ flux-splitting
instability of [#77](https://github.com/szaghi/adam/issues/77) ([MHD blast symmetry](../runtime-amr/mhd-blast-symmetry)):
it appears under both flux families and at $\beta \approx 200$, and only with the field on. The study below traces it
to the divergence cleaning.

### The mirror symmetry of the MHD run (#83)

The measure is $a(t) = \max |\rho(x, y) - \rho(\lambda - x, y)|$, which the exact solution keeps at zero
(`atmosphere_oracle.py symmetry`, with its mean growth rate over $t \in [2, 4]$). The runs are `rt2d-mhd` at $N = 128$
with one cleaning parameter changed at a time (`make_doc_runs.sh <dir> rt-sym`, the FNL build, about 2 min each on two
A30 GPUs); the hydro reference ran on the CPU. No pass/fail bound: a report.

![Magnetised Rayleigh-Taylor: mirror asymmetry per cleaning and damping, rate against c_h](/flume/gravity-rt-symmetry.png)

*Left, the mirror asymmetry per divergence control; centre, per EGLM damping; right, its growth rate against $c_h$,
with the largest growth rate of the linearised cleaned equations on a magnetised shear layer.*

| Run (`rt2d-mhd`, $N = 128$) | rate over $t \in [2, 4]$ | $a$ at $t = 4$ | $a$ at $t = 5$ | max $|\nabla\cdot\mathbf B|$ at $t = 4$ |
|---|---|---|---|---|
| hydro (`rt2d-nl`, CPU) | 1.78 | $1.2\cdot10^{-10}$ | $6.3\cdot10^{-10}$ | — |
| `divergence_control = none` | 2.04 | $1.1\cdot10^{-10}$ | $7.0\cdot10^{-6}$ | **8.6** |
| `glm`, `glm_ch` = 3 | 3.27 | $3.5\cdot10^{-10}$ | $4.2\cdot10^{-6}$ | 0.87 |
| `eglm`, `glm_ch` = 1.5 | 4.38 | $9.9\cdot10^{-10}$ | $7.0\cdot10^{-6}$ | 1.00 |
| `eglm`, `glm_ch` = 3 (the case's default) | 6.25 | $1.4\cdot10^{-7}$ | $2.5\cdot10^{-3}$ | 0.84 |
| `eglm`, `glm_ch` = 6 | 10.91 | $2.7\cdot10^{-2}$ | $5.2\cdot10^{-2}$ | 0.81 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 0 | 6.87 | $5.1\cdot10^{-8}$ | $2.8\cdot10^{-3}$ | 0.86 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 1 | 2.91 | $1.4\cdot10^{-9}$ | $2.8\cdot10^{-3}$ | 0.87 |
| `eglm`, `glm_ch` = 3, `glm_alpha` = 5 | 4.85 | $4.0\cdot10^{-9}$ | $2.8\cdot10^{-3}$ | 0.91 |
| `eglm`, `glm_ch` = 3, `glm_damping_length = min-cell` | 4.09 | $5.4\cdot10^{-10}$ | $3.0\cdot10^{-3}$ | 0.76 |

Related cases: [the well-balanced mode](./well-balanced) (the GV-4b balanced leg); [EGLM cleaning](../mhd/eglm) and
[GLM divergence cleaning](../mhd/glm-cleaning) without gravity; [MHD blast symmetry](../runtime-amr/mhd-blast-symmetry)
(#77, the other symmetry loss, of a different cause).

## Findings

- **The seed must be divergence-free.** A seed in $v$ alone starts sound waves that beat against the mode; with one
  the fitted rate was 0.28 % to 0.92 % off depending on the window ($g = 0.1$, $N = 64$).
- **GV-4c runs unbalanced.** The unbalanced scheme's stationary flow (GV-4b) is as large as the seed but uniform in x:
  its content at the seeded wavenumber is $5.8\cdot10^{-23}$.
- **The residual of GV-4c is the fit window, not the discretisation.** The two resolutions agree more closely than
  either agrees with the theory, and the error falls as the window starts later ($1.1\cdot10^{-3}$, $4.7\cdot10^{-4}$,
  $2.6\cdot10^{-4}$ for $t_0$ = 1.5, 2, 2.5 at $N = 32$): the residual is the seed's transient decaying in the window,
  not the discretisation, which is below $1.2\cdot10^{-4}$ with three cells across the interface width.
- **#83: the cleaning drives the MHD asymmetry.** Without cleaning, the MHD run keeps its symmetry about as well as
  hydro does (rate 2.04 against 1.78). GLM raises the rate to 3.27, and EGLM to 4.4, 6.3 and 10.9 at `glm_ch` = 1.5, 3
  and 6. The excess over the uncleaned run grows by about 1.5 `glm_ch`. It is first seen on the shear flanks between
  bubble and spike.
- **It is not a time-stepping artefact.** At CFL 0.2, which takes 2.5 times the steps, the EGLM rate at `glm_ch` = 3 is
  5.31, against 6.25 at CFL 0.4. A fixed gain per step would have doubled the excess.
- **It is not in the cleaned equations.** `gravity/glm_linear.py` linearises the GLM and EGLM systems about an exact,
  divergence-free double shear layer ($V_0 = 0.5$, width 0.05, field 0.3 along the flow) and solves the spectrum by
  Fourier collocation. At $k = 2\pi$ the largest growth rate is 1.008 ($N = 128$) for both systems and for every $c_h$
  from 0.75 to 6 (right panel). Damping ($\alpha c_h / L = c_h$) and a sheared field leave it independent of $c_h$. At
  $k = 4\pi$ it is 0.52 at $N = 128$ and 0.40 at $N = 256$, independent of $c_h$ at both: the critical layer is
  under-resolved there, so the value itself is not converged. The $c_h$-proportional growth is therefore a property of
  the discretised cleaning. It acts where the discrete $\nabla\cdot\mathbf B$ is not small: max
  $|\nabla\cdot\mathbf{B}| \approx 0.8$, so $\Delta x\,|\nabla\cdot\mathbf{B}| / |\mathbf{B}| \approx 6\,\%$ per cell.
  The GLM system without cleaning ($c_h = 0$, conservative MHD without the Powell terms) is the negative control: it is
  not hyperbolic, and its rate rises with $k$ (1.6 at $4\pi$, 3.2 at $8\pi$).
- **The damping does not cure it.** `glm_alpha` from 0 to 5, or $L = \Delta x$, moves $a$ at $t = 4$ within a factor
  of about 40, with no monotone trend. Every `glm_ch` = 3 run ends with the same jump between $t = 4.5$ and 5, to
  $a \approx 3 \cdot 10^{-3}$. The `glm_ch` = 1.5, GLM and uncleaned runs show a smaller late jump instead, to about
  $10^{-5}$: a current sheet on the spike's stem that narrows to 2 cells at $N = 64$. It is a grid-scale effect, common
  to every run.
- **The cleaning still does its job on $\nabla\cdot\mathbf B$:** 0.8–1.0 at $t = 4$ under every cleaning, against 8.6
  without.
- **Guidance.** For a magnetised run whose symmetry matters (mirror-symmetric instabilities, convergence studies of
  the small scales), set `glm_ch` only just above the fastest signal speed. `glm_ch_check = error` stops the run if
  the flow outruns it. Here 1.5 is enough: the check never fired in the `glm_ch` = 1.5 run. Prefer `glm` to `eglm`:
  at the same speed EGLM's excess is about three times GLM's (4.2 against 1.2). `none` keeps the symmetry best but
  lets $\nabla\cdot\mathbf B$ grow tenfold, so use it only as a diagnostic. With EGLM at `glm_ch` = 3, take the
  asymmetry above about $10^{-6}$ after $t \approx 4.3$ as numerical, not physical. Run at a second `glm_ch` to check:
  a physical asymmetry does not depend on it.
