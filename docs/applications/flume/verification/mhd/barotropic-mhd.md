# Isothermal and barotropic MHD

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| EO (em1 to em10) | Ideal MHD | M13 ([#100](https://github.com/szaghi/adam/issues/100)), P5a (CPU), P5b (FNL) | CPU, FNL | yes | `iso-mhd-tube2`, `baro-mhd-glm-hlld` |

## Case and scope

The [isothermal](../euler/isothermal-eos) and [barotropic](../euler/barotropic-eos) gases gain a magnetic field.
The result is ideal MHD without an energy equation, the total pressure $P(\rho) + |\mathbf B|^2/2$, with the energy
slot inert as for the Euler gases (D-M13-1). Two literature pieces carry it
([numerics](../../numerics#isothermal-and-barotropic-mhd), [input](../../input#eos-optional-equation-of-state)):
- **the six-wave eigensystem of Kim et al. (1999)**, with $a^2 = dP/d\rho$, for the characteristic projection;
- **the isothermal HLLD of Mignone (2007)**, on any $P(\rho)$.

`eos/mhd_oracle.py selftest` checks both papers before any run. It checks $LR = I$ and
$R\Lambda L = \partial F/\partial U$ at random and degenerate states, and the HLLD's consistency, its exactness on a
rotational discontinuity, its hydrodynamic limit and the positivity of its star states. The self-test found eq. 39 of
Mignone misprinted ($S_R$ a minimum) and his speeds out of order on his own tube 1.

What the legs catch:
- **a wrong eigenvector:** a linear wave of one family excites the others; its error stops converging at
  $O(\text{amplitude})$ instead of at fifth order (em1, em2);
- **the sound speed taken as $P/\rho$:** the barotropic waves on $\rho_0 = 5$, where $dP/d\rho = 2.0198$ and
  $P/\rho = 1.4427$, would lag after one period by 12 % of a wavelength (fast) and 4 % (slow), an $O(1)$ error
  (em1);
- **a wrong HLLD:** the Alfvén waves of the circularly polarised wave (em3), Mignone's tubes and their convergence (em4),
  the rotational discontinuities that HLL smears (em6);
- **an axis-dependent frame:** tube 2 along x, y and z (em5);
- **the barotropic path drifting from the isothermal one:** a one-segment table with exponent 1 is the isothermal gas
  (em7);
- **the energy slot written:** every leg checks `rE` exactly 0, gravity included (em10);
- **an input that would silently run wrong:** EGLM, the floors, the gamma-law initial conditions (em8);
- **a backend that disagrees:** one case per model, scheme and solver on CPU and FNL (em9).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Models | `mhd-ideal` with `[eos] type = isothermal` or `barotropic`: `MODEL_MHD_ISO`, `MODEL_MHD_BARO`, and with `divergence_control = glm` `MODEL_MHD_GLM_ISO`, `MODEL_MHD_GLM_BARO` | [`adam_flume_physics_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_physics_object.F90) |
| Library | one body behind the EOS hooks of the isothermal and barotropic modules | `adam_flume_mhd_iso_library_body.INC`, `adam_flume_mhd_{iso,baro}_library.F90` |
| Kernels | the model-agnostic MHD kernel bodies on these libraries (16 CPU and 16 FNL modules; the barotropic table read from the device registry on FNL) | `adam_flume_{cpu,fnl}_mhd[_glm]_{iso,baro}[_{llf,hll,hlld}]_kernels.F90` |
| Space, time | WENO5, characteristic flux splitting (`weno`) or `weno-riemann` with HLLD and the 6th-order correction; SSP-RK(3,3), CFL 0.5 | — |

## Initial and boundary conditions

The barotropic cases use `POLY` of `eos_oracle.py` ($c_0 = 1$, breaks $\rho = 2, 8$, exponents $1, 1.4, 5/3$). The
tubes are Mignone's Table 1 ($a = 1$, $\mathbf B = \mathbf H/\sqrt{4\pi}$) on $[0, 1]$, the discontinuity at 0.5,
extrapolated ends.

| Leg | Setup |
|---|---|
| em1 | one wavelength of a fast, Alfvén or slow wave on the periodic box (the background of Stone et al. 2008, $\mathbf B = (1, \sqrt2, 1/2)$, at rest), amplitude $10^{-7}$, one period: isothermal ($c_f = 2$, $c_a = 1$, $c_s = 1/2$) along x (16, 32 cells) and inclined at $\tan\alpha = 2$ (32², 64²); barotropic on $\rho_0 = 5$ along x; the fast wave with GLM ($c_h = 1.1\,c_f$) |
| em2 | the em1 isothermal waves along x on `weno-riemann` HLLD |
| em3 | the circularly polarised Alfvén wave (Tóth 2000, $B_\perp = 0.1$), isothermal, inclined, right and left polarisation, 32² and 64², CFL 0.1 |
| em4 | tubes 1–3 on `weno-riemann` HLLD, 200, 400, 800 and 3200 cells (one transverse cell per block); tube 3 at 400 cells with `positivity_limiter = cell` |
| em5 | tube 2 along x, y and z (the states permuted cyclically), 200 cells, `weno` and `weno-riemann` |
| em6 | tube 2 at 400 cells, HLLD and HLL, against the 3200-cell HLLD run |
| em7 | tube 2 with $c_0 = 1.3$: the isothermal gas against the barotropic one-segment table of exponent 1 |
| em8 | `eglm`; `p_floor = 1e-6`; `orszag-tang`; the `entropy` wave |
| em9 | the last checkpoints of CPU and FNL: `em1-iso-2d-fast-64`, `em1-baro-1d-slow-32`, `em1-iso-glm-1d-fast-32`, `em2-iso-hlld-1d-alfven-32`, `em3-cpaw-right-64`, `em4-tube2-400`, `em5-weno-y`, `em7-baro`, `em10-barotropic` (not the limited tube 3: the limiter is not bitwise across backends) |
| em10 | the em1 box with amplitude 0 under $g_x = 0.5$, 48 cells, GLM ($c_h = 4$), SSP(5,4), $t = 0.4$: isothermal $\rho_0 = 1$, barotropic $\rho_0 = 5$ |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/eos/make_eos_mhd.py` | writes the cases: `tube` (from the Sod inputs), `wave` (from `mhd/linear-wave/make_linear_wave.py`, the period from the speed of the gas), `cpaw` (from `mhd/cpaw/make_cpaw.py`) |
| `src/tests/flume/verification/eos/mhd.sh` | runs the legs (`--legs em1,...,em10`, default all but em9, which compares the two backends' runs; `FLUME_EXE`, `EM_NP` ranks) |
| `eos/mhd_oracle.py` (`selftest`, `rotation`, `selfconv`, `agree`, `hlld-gain`, `backends`), `mhd/linear-wave/linear_wave_oracle.py`, `mhd/cpaw/cpaw_oracle.py`, `eos/eos_oracle.py` (`freefall`) | the checks |
| `src/tests/flume/regression/iso-mhd-tube2/input.ini`, `baro-mhd-glm-hlld/input.ini` | the frozen tube 2 (isothermal, `weno`) and a barotropic GLM HLLD tube |
| `make_doc_runs.sh <runs> eosmhd`, `make_doc_figures.py --only barotropic-mhd` | the figure |

```bash
cd src/tests/flume/verification/eos && ../../../../../exe/.regression-venv/bin/python mhd_oracle.py selftest
./mhd.sh                                                                                   # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./mhd.sh                                   # FNL
./mhd.sh --legs em9                                                                        # CPU against FNL
```

## Results

![Isothermal and barotropic MHD: linear-wave convergence and Mignone's tubes 1-3](/flume/barotropic-mhd.png)

*Top left: the linear waves along x on 16, 32 and 64 cells, the error after one period over the amplitude, isothermal
and barotropic, the three families; the legend gives the order of the 32/64 pair, 4.97 to 5.05 (slope 5 dotted). Top right: tube 1, $\rho$ at 400 cells against the 3200-cell run. Bottom left: tube
2, $B_y$ at 400 cells with HLLD and with HLL against the reference: the rotational discontinuities. Bottom right:
tube 3, the two Mach-5 streams with $B_n = 0$.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| — | `mhd_oracle.py selftest`: Kim's eigensystem, the embedding, the speeds, the HLLD (consistency, rotational discontinuity, hydro limit, positivity), the CPAW | round-off | 13 checks pass | — |
| em1 | the order of the finest pair: isothermal fast, Alfvén, slow along x; inclined; barotropic; GLM fast | $\ge 4.7$ | 4.96, 4.96, 4.96; 4.97, 4.97, 5.00; 4.97, 4.96, 4.96; 4.96 | the same to the printed digits |
| em2 | the same on `weno-riemann` HLLD: fast, Alfvén, slow | $\ge 4.7$ | 4.96, 4.96, 4.96 | 4.96, 4.96, 4.96 |
| em3 | the CPAW order, right and left; right against left | $\ge 4.7$; 0 | 4.89, 4.89; 0 | 4.89, 4.89; 0 |
| em4 | the ladder order of the L1 distance of $\rho$ to 3200 cells, tubes 1, 2, 3; $\rho > 0$, `rE` | $\ge 0.5$; 0 | 1.15, 1.03, 0.66 (pairs 1.09 and 1.22, 1.03 and 1.03, 0.17 and 1.15); positive, 0 | the same |
| em4 | tube 3 with the limiter: $\rho > 0$, `rE`; the L1 distance (reported) | 0 | positive, 0; $9.0\cdot10^{-3}$ (unlimited $5.4\cdot10^{-4}$) | positive, 0; $9.0\cdot10^{-3}$ |
| em5 | tube 2 along x, y, z, `weno` and `weno-riemann` | 0 | 0, 0 | 0, 0 |
| em6 | the L1 distance of $B_y$ to 3200 cells, HLLD against HLL | HLLD < HLL | $2.35\cdot10^{-3}$ against $2.79\cdot10^{-3}$ (0.84) | the same |
| em7 | the one-segment table against the isothermal gas | $10^{-10}$ | 0 | $6.2\cdot10^{-13}$ |
| em8 | the 4 invalid inputs | refused | refused | refused |
| em9 | CPU against FNL, the 9 cases: the largest gap over the fields | $10^{-10}$ | (the same comparison) | $2.9\cdot10^{-13}$ (tube 2 along y, `weno`); the waves and the CPAW $\le 1.7\cdot10^{-14}$; `rE` 0 |
| em10 | $\rho$ uniform, $\rho u - \rho g t$ relative, `rE`: isothermal, barotropic (SSP(5,4), since [#104](https://github.com/szaghi/adam/issues/104)) | $10^{-12}$; 0 | 0, $1.1\cdot10^{-14}$, 0; 0, $7.0\cdot10^{-15}$, 0 | 0, $1.1\cdot10^{-16}$, 0; 0, $2.2\cdot10^{-15}$, 0 |

Related cases: the [linear waves](./order-of-accuracy) and the [HLLD verification](./weno-riemann) of the ideal gas.

## Findings

- **The literature, checked before the code.** Kim et al. normalise the fast pair with
  $\alpha_f = \sqrt{c_f^2 - b_n^2}/\sqrt{c_f^2 - c_s^2}$, which is the Ryu–Jones $\alpha_f$ times $c_f/a$. The printed
  set passes $LR = I$ and $R\Lambda L = \partial F/\partial U$ at random and degenerate states, both gases. In Mignone's
  section 3, eq. 39 prints $S_R$ as a minimum: it must be the maximum, since it bounds both states' fast speeds. On his
  own tube 1 ($\rho_L/\rho_R = 10$) the first face has $S^*_R = 3.90 > S_R = 3.37$, a star region of negative width that
  the paper does not discuss, because $u^* = F_\rho^\mathrm{hll}/\rho^\mathrm{hll}$ is driven far off by a strong density
  jump. There FLUME's HLLD returns HLL's flux, as the adiabatic HLLD does.
- **The sound speed is $dP/d\rho$, and only there.** The Jacobian holds the gas law only through $a^2 = dP/d\rho$, so
  one eigensystem serves both gases. On $\rho_0 = 5$ the barotropic waves converge at fifth order with
  $c^2 = 2.0198$; with $P/\rho$ the fast and slow waves would lag by 12 % and 4 % of a wavelength per period.
- **HLLD's gain over HLL is modest at this resolution.** On tube 2, 400 cells, its $B_y$ is 16 % closer to the
  reference than HLL's (the leg asserts only that it is closer).
- **A collision with $B_n = 0$ keeps its start-up error.** With no normal field, isothermal MHD has a zero-speed
  tangential mode: a density jump at rest, balanced by the magnetic pressure. Where the two streams of tube 3 meet,
  the start-up error of the collision falls into that mode and stays, the MHD analogue of Noh's wall heating (an
  isothermal hydrodynamic contact cannot carry a density jump, so there the error leaves as sound). At 400 cells the
  centre reaches $\rho = 2.68$ with `weno` (flux splitting) against a plateau of 0.40, two cells wide; `weno-riemann`
  HLLD keeps it at 0.405. The positivity limiter does act there on `weno-riemann` (the limited and unlimited runs differ; the
  likely trigger is its kinetic chord, which bounds a limited cell's speed by about $2|\mathbf u_b| + c$ of the
  first-order backbone, not measured cell by cell), and it raises the centre to 1.42
  (L1 distance $9.0\cdot10^{-3}$, unlimited $5.4\cdot10^{-4}$). The same limiter on the isothermal hydrodynamic Mach-5
  collision lowers the L1 error (from $1.99\cdot10^{-2}$ to $1.43\cdot10^{-2}$, `weno-riemann` HLLC), and on `weno`
  it never acts. For magnetised collisions with $B_n \approx 0$, prefer `weno-riemann` HLLD, and enable the limiter
  only where positivity is at risk.
- **The backends agree to round-off, the one-segment table only to $10^{-12}$ on FNL.** On the CPU the one-segment table
  of exponent 1 and the isothermal gas give the same bits; on FNL they differ by $6.2\cdot10^{-13}$. The table path
  evaluates $K\rho^{\Gamma}$ through a power and the isothermal path $c_0^2\rho$ directly, and the device rounds the
  two differently. Across backends the same input agrees to $2.9\cdot10^{-13}$ at most (em9).
- **The tube-3 criterion was changed after its first run.** Tube 3 has shocks only, so its L1 distance is set by where
  each shock falls within its cell: the pair orders are 0.17 and 1.15. em4 asserts the order over the whole ladder
  (0.66) and prints the pairs.
- **One transverse cell per block.** The Sod inputs the tubes are built from carry four cells per block along each
  null axis: 16 copies of the same 1-D solution. `make_eos_mhd.py` sets one, bitwise the same profile (tube 2 at 400
  cells, checked), and the 3200-cell references take minutes instead of hours.
- **The free fall found a defect of the library's SSP(5,4), now fixed.** Its former weights, printed with 14 digits,
  summed to $1 - 8.8\cdot10^{-11}$, so the free fall on SSP(5,4) gave $\rho u - \rho g t$ of $3.5\cdot10^{-11}$ after
  $t = 0.4$; em10 ran on SSP(3,3) until the tableau was rebuilt from its Shu–Osher form
  ([#104](https://github.com/szaghi/adam/issues/104), [time integration](../../numerics#time-integration)). On SSP(5,4)
  em10 now gives $1.1\cdot10^{-14}$ and $7.0\cdot10^{-15}$, and runs there as the solver-level regression of that fix.
