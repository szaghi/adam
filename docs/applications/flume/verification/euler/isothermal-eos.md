# Isothermal equation of state

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EO (eo1 to eo7) | Euler | M13 ([#100](https://github.com/szaghi/adam/issues/100)), P1, P2 | CPU, FNL | yes | `iso-tube`, `iso-collide-hllc` |

## Case and scope

`[eos] type = isothermal` closes the Euler equations with $P = c_0^2\rho$ ([numerics](../../numerics#isothermal-and-barotropic-equations-of-state),
[input](../../input#eos-optional-equation-of-state)). The energy equation is gone, the energy slot `rE` of the state is
inert, and the isothermal Riemann problem has an exact solution built from mass and momentum alone. Every run is
checked against that solution, which `eos_oracle.py` constructs from the jump conditions, independently of the code
(its self-test, P0, checks it against properties it is not built from).

What the case catches:
- **a wrong eigensystem or Roe average:** the characteristic projection would couple the transverse velocity to the
  acoustic waves, or move the shock and the plateau;
- **a wrong isothermal flux:** the plateau, the shock speeds and the compression $\rho_* = \rho M_1^2$ of a collision;
- **an energy slot that is not inert:** any non-zero `rE`, in any checkpoint;
- **a loss of the design order** of the smooth flow, or of mass conservation;
- **a gamma-law path reached by the isothermal model:** every input the isothermal gas does not define yet must be
  refused with its own message, not run on a meaningless pressure.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` with `[eos] type = isothermal`, `sound_speed` $c_0 = 1$: `MODEL_EULER_ISO`, the energy slot inert (D-M13-1) | [`adam_flume_euler_iso_library_body.INC`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_euler_iso_library_body.INC), [`adam_flume_physics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_physics_object.F90) |
| Space | WENO5 characteristic flux splitting (eo1 to eo5); `weno-riemann` with the characteristic interpolation, the 6th-order flux correction and `llf`, `hll`, `hllc` (eo6) | `adam_flume_cpu_euler_iso_kernels.F90`, `adam_flume_cpu_euler_iso_{llf,hll,hllc}_kernels.F90`; FNL `adam_flume_fnl_euler_iso_kernels.F90`, `adam_flume_fnl_euler_iso_{llf,hll,hllc}_kernels.F90` |
| Time | SSP-RK(3,3), CFL 0.5; SSP-RK(5,4) at CFL 0.1 for the order (eo3), so the time error stays below the WENO5 space error | — |

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| eo1 | the shock tube: $\rho = 1 \mid 0.125$ at $x = 0.5$, at rest, the transverse velocity $v = 0.3 \mid -0.2$; 400 cells along $x$, $y$ active with 4 cells (uniform), extrapolation faces, to $t = 0.2$ |
| eo2 | two streams $\rho = 1$, $u = \pm 5$ colliding at $x = 0.5$, 400 cells, to $t = 1$ |
| eo3 | a right-running acoustic wave, $\rho = 1 + 10^{-6}\sin 2\pi x$, $u = c_0(\rho - 1)$, on the periodic $[0, 1]$, $N = 12, 24, 48$ in two blocks on two ranks, one period |
| eo4 | the eo1 input with `[cooling]`, `[sinks]`, `[gravity]`, a viscosity, `positivity_limiter = cell`, `[reference]`, an immersed cylinder, a `wall-noslip` or `hydrostatic` face, `isentropic-vortex`; `mhd-ideal`, passive scalars, `cp`/`cv`, no or a negative `sound_speed`, an unknown `type` |
| eo5 | the eo3 wave at amplitude 0.3, $N = 48$, two ranks, to $t = 1$: it breaks into a shock at $t_b \approx 1/(c_0 A k) = 0.53$ |
| eo6 | eo1 with `weno-riemann` and each solver |
| eo7 | the eo1, eo2, eo3 ($N = 48$) and eo6 `hllc` runs of CPU and FNL |

**The exact solution.** Across a shock the isothermal jump gives $u_* = u_K \mp c_0(\rho_*/\rho_K - 1)
\sqrt{\rho_K/\rho_*}$, across a fan $u_* = u_K \mp c_0\ln(\rho_*/\rho_K)$; $\rho_*$ is the root of the two, and the
contact between them carries only the jump of the transverse velocity. In the tube $\rho_* = 0.34578$,
$u_* = 1.0619$: a left fan, the contact at $x = 0.7124$, a right shock at $0.8326$. The collision compresses the gas to
$\rho_* = M_1^2 = 26.963$, $M_1 = 5.193$ the Mach number in the frame of the shocks, which move out at $\pm 0.1926$.

**A 1-D run cannot carry a transverse velocity.** A null direction freezes its momentum for every Euler model
(`null_freeze`), so in a 1-D tube $\rho v$ stays at its initial profile and the contact cannot move. The tube keeps $y$
active, with four uniform cells whose fluxes cancel exactly.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/eos/make_eos.py` | writes the inputs from the `sod-x` regression input (`--case tube|collide|wave`, the scheme, `--set`/`--unset` for the contract) |
| `src/tests/flume/verification/eos/check.sh` | runs the legs (`--legs eo1,...,eo7`, default all but eo7; `FLUME_EXE` selects the backend; eo7 compares the work directories of both) |
| `src/tests/flume/regression/iso-tube/input.ini`, `iso-collide-hllc/input.ini` | the frozen eo1 tube (WENO5 splitting) and eo2 collision (`weno-riemann` HLLC), 100 cells in four blocks |
| `src/tests/flume/verification/eos/eos_oracle.py` | the exact isothermal (and barotropic) Riemann solutions, the linear wave, and the comparisons with the runs (`riemann`, `wave`, `inert`) |
| `make_doc_runs.sh <runs> eos`, `make_doc_figures.py --only isothermal-eos` | the figure (7 runs, 4 minutes) |

```bash
cd src/tests/flume/verification/eos && ./check.sh                                      # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                           # FNL
./check.sh --legs eo7                                                                 # CPU against FNL
```

## Results

![Isothermal equation of state: the shock tube, the collision, the order and the steepened wave](/flume/isothermal-eos.png)

*Top left: the shock tube against the exact solution; the contact carries the transverse velocity only. Top right: the
Mach-5 collision, the dense layer at $\rho_* = M_1^2$ behind a decaying train of post-shock oscillations. Bottom left:
the L1 error of the acoustic wave after one period, fifth order. Bottom right: the wave at amplitude 0.3, broken
into a shock after $t \approx 0.53$; the total mass is constant to round-off and `rE` stays 0.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| eo1 | the plateau mean against $\rho_*$ | $10^{-3}$ | $8.3\cdot10^{-5}$ (108 cells; largest deviation $1.1\cdot10^{-3}$) | the same |
| eo1 | the right shock; the contact (mid-$v$ crossing) | 1.5 cells; 2 cells | 0.13 cells; 0.11 cells | the same |
| eo1 | `rE`; the L1 errors of $\rho$ and $v$ (reported) | exactly 0; — | 0; $2.1\cdot10^{-3}$, $2.0\cdot10^{-3}$ | the same |
| eo2 | the layer mean against $M_1^2$; both shocks | $10^{-3}$; 1.5 cells | $8.3\cdot10^{-5}$ (144 cells; largest deviation $1.7\cdot10^{-2}$); 0.05 cells each | the same |
| eo3 | the L1 order of the finest pair | $> 4.5$ | 4.84 (4.94 on 12 → 24; L1 $2.8\cdot10^{-12}$ at $N = 48$) | the same |
| eo4 | the 16 invalid inputs | refused | refused | refused |
| eo5 | `rE` in every checkpoint; the total mass | exactly 0; $10^{-13}$ | 0; $1.1\cdot10^{-16}$ | 0; $2.2\cdot10^{-16}$ |
| eo6 | eo1 with `llf`, `hll`, `hllc` | as eo1 | plateau $6.1\cdot10^{-5}$, $3.2\cdot10^{-5}$, $3.2\cdot10^{-5}$; shock 0.15, 0.07, 0.07 cells; contact 0.16 cells each; `rE` 0 | the same |
| eo7 | FNL against CPU, the last checkpoint, per field over its scale (a momentum's at least $\max\rho\,c_0$); `rE` | $10^{-10}$; exactly 0 | (the same comparison) | tube $1.7\cdot10^{-14}$, collision $4.2\cdot10^{-12}$, wave $2.4\cdot10^{-15}$, HLLC tube $1.2\cdot10^{-14}$ (largest field); 0 |

Related cases: [Sod](./sod) is the grid; [radiative cooling](./radiative-cooling) checks the isothermal jump at the
floor temperature of a radiative shock, which an isothermal gas makes exact.

## Findings

- **The isothermal gas needs no pressure in its Riemann problem, and no star pressure in its HLLC.** The contact
  carries only the transverse velocity, so the isothermal HLLC is HLL on mass and normal momentum with the transverse
  momenta upwinded by the mass flux.
- **HLLC's contact advantage does not show here.** On the tube HLL and HLLC agree to the third digit (the
  transverse-velocity L1 $2.139\cdot10^{-3}$ against $2.134\cdot10^{-3}$): the characteristic interpolation and the
  6th-order flux correction already resolve the contact at 400 cells, so the Riemann solver's own contact resolution is
  not what limits it. eo6 checks the solvers are correct, not that HLLC is sharper.
- **A strong isothermal shock leaves a decaying oscillation train.** At Mach 5 the dense layer oscillates about
  $\rho_*$ by 2 % next to each shock, falling to 0.2 % at the centre: the post-shock oscillation of WENO flux splitting
  at a strong shock (Arora & Roe 1997), not an error of the gas. Its mean is within $10^{-3}$ of $M_1^2$.
- **CPU and FNL agree to round-off once a field is measured on its own scale.** On the $10^{-6}$ acoustic wave the
  momentum differs by $2\cdot10^{-15}$ absolute, which is $2\cdot10^{-9}$ of the momentum's own size but $2\cdot10^{-15}$
  of $\rho c_0$. The momentum flux carries the pressure $c_0^2\rho$, so its round-off is that of an $O(1)$ quantity,
  and eo7 measures a momentum against at least $\max\rho\,c_0$. The largest gap of the four runs, $4.2\cdot10^{-12}$
  in the Mach-5 collision's density, is round-off amplified at two strong shocks.
- **The inert slot costs nothing and stays exact.** Zero flux, zero source and a zero eigenvalue keep `rE` at 0 to the
  bit in every checkpoint, also across shocks; the state layout and every model-agnostic kernel are unchanged.
