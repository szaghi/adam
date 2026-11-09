# The hybrid central/WENO scheme

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| HY | [#119](https://github.com/szaghi/adam/issues/119) · T2.1–T2.7 | P+C | hy1–hy13 | CPU, FNL | yes |

## Case and scope

`scheme_space = hybrid-kep-weno` ([turbulence](../../turbulence#hybrid-central-weno)) uses a 6th-order
kinetic-energy- and pressure-equilibrium-preserving central flux at the faces away from shocks and keeps WENO where a
Ducros-type sensor flags a cell (T2.1–T2.3). WENO is kept on the skin of 2:1 seams and physical faces by default
(`hybrid_skin`, T2.4). The legs check the central flux, pressure equilibrium (T2.5), the model masked to central cells
(T2.6), conservation on every grid path and robustness at shocks (T2.7). The campaign reruns the
[ILES baseline](./iles-baseline) cases and compares every scheme with the 512³ DNS.

## Legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | Phase | What it verifies | Measured |
|---|---|---|---|
| hy1 | T2.7 | inviscid Taylor-Green: dissipated $K$, hybrid against WENO | 47 times less on FNL |
| hy2 | T2.3 | one rank against 2 with the hybrid scheme on the random field, where the sensor fires | not recorded |
| hy3 | T2.2 | the central flux (`test_kep`): order, kinetic energy identity, pressure equilibrium | order 5.96–5.99; $4\cdot10^{-16}$; $dp/dt$ $9\cdot10^{-15}$ across 100:1 |
| hy4 | T2.1 | contact of density ratio 8 and 100 | $u$, $p$ within $9\cdot10^{-14}$ |
| hy5 | T2.4 | 2:1 box, conservation; passes with any skin | mass $10^{-16}$ |
| hy6 | T2.4 | Sod against WENO | L1 $8.3\cdot10^{-5}$ |
| hy7 | T2.4 | runtime regridding; a forest with open ends | $10^{-16}$; $5.7\cdot10^{-10}$, the single realm's drift |
| hy8 | T2.5 | density wave, amplitude 0.2; WENO holds it too, so this does not separate them | $p$, $u$ $10^{-14}$ |
| hy9 | T2.6 | model only off flagged cells; one rank against 2; conservation (does not prove the mask on ghost cells) | $\varepsilon_\mathrm{sgs}$ exactly 0 with all flagged; $7\cdot10^{-14}$; $2.3\cdot10^{-16}$ with a model on a 2:1 box; whether a flagged cell faces a free one at an interface is not checked |
| hy10 | T2.4 | the skin relaxed on 2:1 seams (the hy5 case, `central-seams`): what conservation does, the last history rows | a measure with no bound: printed, never fails the script |
| hy11 | T2.4 | relaxing the skin changes the density: seams, open ends (does not catch a shifted band) | $1.6\cdot10^{-4}$, $6.4\cdot10^{-5}$ |
| hy12 | T2.4 | restart on a refined grid, every cell | 0 |
| hy13 | T2.7 | against WENO, L1: a vortex through a shock (2-D), Shu-Osher and a shock through 2:1 seams (1-D, no vorticity, so the Ducros ratio flags any varying flow): robustness, no accuracy claim | $1.0\cdot10^{-4}$, $1.7\cdot10^{-5}$, $3.4\cdot10^{-5}$ |

## Taylor-Green vortex and isotropic turbulence

| Setup | |
|---|---|
| cases | the [ILES baseline](./iles-baseline) inputs with `scheme_space = hybrid-kep-weno` |
| sensor | thresholds 0.2 and 0.1 and the floor at their defaults (not calibrated) |
| model | none, or Vreman acting off flagged cells only (T2.6) |
| run | `turbulence/tgv_campaign.sh`, 4 ranks on A30 |

![The hybrid scheme on the Taylor-Green vortex](/flume/tgv-hybrid.png)

*Top: M 0.1, the hybrid against WENO-5; middle: with Vreman; bottom: $M_0 = 1.25$. Left $-dK/dt$, right the enstrophy; dashed: the 512³ spectral DNS, on the M 0.1 rows.*

| Case, M 0.1 | Peak $-dK/dt$, raw | Over $\Delta t = 0.5$ | at $t$ | Flicker, $7<t<10$ | Peak $Z$ | Numerical share to the peak |
|---|---|---|---|---|---|---|
| WENO-5, $64^3$ | 0.01377 | 0.01357 | 5.09 | 0.3 % | 1.71 | 77 % |
| hybrid, $64^3$ | 0.01591 | 0.01334 | 8.09 | 7.2 % | 3.82 | 50 % |
| hybrid + Vreman, $64^3$ | 0.01347 | 0.01173 | 8.21 | 3.5 % | 2.92 | |
| WENO-5, $128^3$ | 0.01242 | 0.01216 | 8.29 | 0.5 % | 3.40 | 58 % |
| hybrid, $128^3$ | 0.01389 | 0.01317 | 8.48 | 2.0 % | 6.86 | 18 % |
| hybrid + Vreman, $128^3$ | 0.01258 | 0.01226 | 8.58 | 0.8 % | 5.53 | |
| WENO-5, $256^3$ | 0.01272 | 0.01261 | 8.72 | 0.7 % | 6.18 | 27 % |

The flicker is the rms of $-dK/dt$ about its running mean, over the smoothed peak. The share with Vreman is left out: it
does not separate the model from the scheme. Against the same runs with the model everywhere, $Z$ moves under 1 % at
M 0.1, where the sensor rarely fires.

| Case, $M_0 = 1.25$ | Peak $-dK/dt$ | Peak $Z$ | Numerical share |
|---|---|---|---|
| WENO-5, $64^3$ / $128^3$ / $256^3$ | 0.01383 / 0.01081 / 0.01101 | 1.81 / 3.34 / 5.95 | 80 % / 46 % / 20 % |
| hybrid, $64^3$ / $128^3$ / $256^3$ | 0.01266 / 0.01057 / 0.01112 | 2.20 / 4.83 / 7.21 | 73 % / 36 % / 10 % |

| Isotropic turbulence | Peak $Z$, $64^3$ / $128^3$ | $K(t/\tau=4)$, $64^3$ / $128^3$ |
|---|---|---|
| WENO-5 | 34.5 / 51.0 (60.5 at $256^3$) | 0.303 / 0.384 (0.413 at $256^3$) |
| hybrid | 41.0 / 57.4 | 0.367 / 0.408 |
| hybrid + Vreman | 37.8 / 54.1 | |

**What it shows.**
- **At M 0.1 the hybrid at $128^3$ reaches a higher enstrophy peak than WENO-5 at $256^3$**, and at $64^3$ more than
  WENO-5 at $128^3$, with the dissipation peak where the finer runs have it. $Z$ is still below the DNS's 10.3.
- **At $M_0 = 1.25$ the gain is smaller**: the sensor keeps WENO around the shocklets.
- **The raw peak of $-dK/dt$ is not a measure on the hybrid.** It flickers, 7 % at $64^3$, 2 % at $128^3$. The source is
  not established; the sensor switching on and off is the first candidate. Smoothed, the hybrid is 4–6 % above WENO-5 at
  $256^3$.
- **Vreman on the hybrid** removes half of the flicker and lowers $Z$ by 20 %. At $64^3$ it lands below the finer runs
  (0.0117). Against the DNS it makes the hybrid worse at both resolutions (below).

## On a 2:1 box, and what the skin costs

| Setup | |
|---|---|
| grid | M 0.1, base $64^3$, a box of half the side refined once |
| skin | `hybrid_skin = weno` (default) or `central-seams` |

![The hybrid scheme on a 2:1 box](/flume/tgv-hybrid-box.png)

*M 0.1, a box of half the side refined once on $64^3$. Left $-dK/dt$, right the enstrophy.*

| Case | Smoothed peak $-dK/dt$ | Peak $Z$ |
|---|---|---|
| WENO-5, box | 0.01304 | 1.87 |
| hybrid, box, WENO skin (default) | 0.01239 | 3.69 |
| hybrid, box, central skin | 0.01319 | 4.13 |
| hybrid, uniform $64^3$ | 0.01333 | 3.82 |

**What it shows.**
- **With the WENO skin the box ends below the uniform run in $Z$**: the skin takes more than the extra cells give.
- **The central skin is stable to $t = 20$ and above the uniform run.** Against the DNS it is the closer one (table
  below), on this smooth flow; the case with a shock has no reference.
- **A shock through the seams runs with the central skin.** A vortex after Jiang & Shu (the code's profile, same strength
  and radius, unit square) through a Mach 1.1 shock crossing the box seams differs from the WENO skin by at most
  $2.8\cdot10^{-4}$ in density at $t = 0.35$.

## Against the DNS

M 0.1 against the 512³ DNS, from `turbulence/tgv_dns.py`. Peak $-dK/dt$ over a window of 0.5; L2 errors over
$0.3 < t < 19.7$ relative to the DNS, $-dK/dt$ over the same window. This table collects every T0, T1 and T2 run at M 0.1.

| Case | Peak $-dK/dt$ | at $t$ | Peak $Z$ | at $t$ | L2 $-dK/dt$ | L2 $Z$ | L2 $K$ |
|---|---|---|---|---|---|---|---|
| DNS | 0.01286 | 8.97 | 10.29 | 8.97 | | | |
| WENO-5, $64^3$ | 0.01357 | 5.09 | 1.709 | 6.25 | 0.660 | 0.844 | 0.187 |
| WENO-5, $128^3$ | 0.01216 | 8.29 | 3.400 | 8.11 | 0.278 | 0.643 | 0.083 |
| WENO-5, $256^3$ | 0.01260 | 8.72 | 6.182 | 8.56 | 0.090 | 0.331 | 0.020 |
| HLLC, $128^3$ | 0.01235 | 8.32 | 3.585 | 8.11 | 0.260 | 0.623 | 0.070 |
| WENO-7, $128^3$ | 0.01215 | 8.11 | 3.171 | 8.14 | 0.291 | 0.688 | 0.101 |
| WENO-9, $128^3$ | 0.01274 | 8.27 | 4.454 | 8.03 | 0.177 | 0.539 | 0.037 |
| Vreman, $64^3$ | 0.01373 | 5.11 | 1.591 | 6.23 | 0.671 | 0.854 | 0.200 |
| Vreman, $128^3$ | 0.01193 | 8.28 | 3.177 | 8.15 | 0.295 | 0.662 | 0.092 |
| WALE, $128^3$ | 0.01183 | 8.26 | 3.164 | 8.16 | 0.296 | 0.660 | 0.091 |
| Sigma, $128^3$ | 0.01192 | 8.29 | 3.184 | 8.14 | 0.293 | 0.660 | 0.091 |
| WENO-5, 2:1 box | 0.01304 | 6.41 | 1.865 | 6.39 | 0.617 | 0.829 | 0.178 |
| hybrid, $64^3$ | 0.01333 | 8.09 | 3.822 | 7.86 | 0.237 | 0.615 | 0.057 |
| **hybrid, $128^3$** | 0.01317 | 8.47 | 6.856 | 8.25 | 0.067 | 0.271 | 0.014 |
| hybrid + Vreman, $64^3$ | 0.01173 | 8.21 | 2.917 | 8.09 | 0.283 | 0.678 | 0.086 |
| hybrid + Vreman, $128^3$ | 0.01226 | 8.58 | 5.529 | 8.33 | 0.114 | 0.400 | 0.029 |
| hybrid, 2:1 box, WENO skin | 0.01239 | 8.14 | 3.688 | 7.88 | 0.269 | 0.618 | 0.063 |
| hybrid, 2:1 box, central skin | 0.01319 | 8.11 | 4.134 | 7.91 | 0.236 | 0.580 | 0.059 |

**What it shows.**
- **The hybrid at $128^3$ is closer than WENO-5 at $256^3$** on all three errors, with one eighth of the cells. At
  $64^3$ it is closer than WENO-5 at $128^3$.
- **Its peak of $-dK/dt$ is 2.4 % above the DNS** and comes at $t = 8.47$ rather than 8.97; $Z$ at $128^3$ reaches two
  thirds of the DNS peak, so the gradients are not converged.
- **A model makes it worse.** Vreman takes the $-dK/dt$ error from 0.067 to 0.114 and the $Z$ error from 0.27 to 0.40 at
  $128^3$. All three models also worsen WENO-5. At Re 1600 the subgrid scales are not far from resolved, and the
  constants were not tuned; on these integrals no model is recommended.
- **WENO-9** at $128^3$ (0.177) is much closer than WENO-5 (0.278); **WENO-7** is further (0.291), see
  [the ILES baseline](./iles-baseline#taylor-green-vortex).
- **The skin.** On the 2:1 box the WENO skin (0.269) is further than uniform $64^3$ (0.237); the central skin ties it on
  $-dK/dt$ (0.236) and beats it on $Z$ (0.580 against 0.615). Smooth flow only.

## The flow

![Vorticity of the Taylor-Green vortex and of the isotropic turbulence](/flume/tgv-fields.png)

*`hybrid-kep-weno` at $128^3$ (`tgv_cases_vis.txt`). Taylor-Green, Re 1600, Mach 0.1: $\omega_x$ on the
periodic face $x\approx0$ at $t=4.2$ and $9.4$ (the peak of the dissipation is at $t\approx9$; the Workshop's figures show the vorticity norm on
this face) and $\omega_z$ on the symmetry plane $z=\pi$ at $t=9.4$. Isotropic turbulence, $M_t$ 0.6, $Re_\lambda$ 100, $k_0=4$: $\omega_z$ on
$z=\pi$ at $t=2$. Planes from `turbulence/field_slices.py --time`; `make_doc_figures.py --only tgv-fields`.*

- **At $t=4.2$ the vorticity is still in sheets; by $t=9.4$ they have rolled up and broken into the small structures of the transition.**
- **The flow's mirror symmetries survive the breakdown.** On the face $x\approx0$, $\omega_x$ is antisymmetric about $y=\pi$ and about $z=\pi$ to
  $8\cdot10^{-13}$ of its maximum at $t=9.4$ (0 at $t=4.2$), as the initial condition is.
- **On $z=\pi$ the flow is still large-scale at $t=9.4$**: four counter-rotating cells with thin layers at their edges.
- **The isotropic field shows no visible imprint of the blocks** ($2^2$ per direction). One plane, one time: a description, not a measurement.

## Cost

200 steps on $128^3$, two A30, three samples:

| Case | Wall time |
|---|---|
| WENO | 24.1–24.6 s |
| hybrid, every cell flagged | 24.2–24.7 s |
| hybrid, every cell central | 26.7–26.9 s |

The sensor costs under 1 %, the central flux about 10 % on top of WENO, startup included. The scheme runs both flux
kernels; a single kernel with a mask was not tried, so what skipping WENO would save is not measured.

## Not done

- Calibration of the sensor thresholds and floor.
- The source of the flicker of $-dK/dt$.
- A model masked by face rather than by cell; the mask on ghost cells at an interface.
- `central-bcs` on a forest seam (the sensor ghost is computed locally there).
- A shock case with a reference for the relaxed skin.
- A single flux kernel with a mask.
