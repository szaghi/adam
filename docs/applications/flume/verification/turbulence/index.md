# Turbulence (T0, T1, T2, T3)

[#119](https://github.com/szaghi/adam/issues/119) measures what the scheme does to turbulence before anything is built on it. The
diagnostics, the keys and the models are in [input](../../input) and [numerics](../../numerics#subgrid-scale-models).

## Contract legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | What it verifies | Measured |
|---|---|---|
| tt1 | first history row of `taylor-green` against the exact integrals, M 0.1, 1.25 | K $10^{-16}$; Z and $\varepsilon$ order 3.98–4.00; $\varepsilon_d$, $p_d$ $10^{-30}$ |
| tt2 | one rank against 2 | $6\cdot10^{-14}$ |
| tt3 | energy budget to $t = 2$, Re 1600, M 0.1: what the integrals leave of the $K$ drop | 70 %, 8.8 %, 0.28 % at $32^3$, $64^3$, $128^3$ |
| tt6 | pressure work, M 0.8: the budget closes only with $p_d$ | 8.5 % left, 203 % without, 414 % flipped |
| ht1 | `isotropic-turbulence` field: every shell, $K$, dilatation | shells $7\cdot10^{-15}$, $K$ $1.5\cdot10^{-14}$ |
| ht2 | one rank against 2 on the random field | $7\cdot10^{-15}$ |
| sg1 | the formulas against numpy on 200 tensors; the zeros of each model | $3.6\cdot10^{-15}$, $1.2\cdot10^{-15}$, $1.8\cdot10^{-11}$ |
| sg2 | budget with the model, M 0.4: Vreman 36 %, 5.7 %, 0.27 % left at $32^3$, $64^3$, $128^3$ | at $128^3$: 21.6 % without $\varepsilon_\mathrm{sgs}$, 43 % flipped |
| sg3 | no viscosity with each model | $\varepsilon = 0$, mass and energy $3.3\cdot10^{-16}$ |
| sg4 | one rank against 2 with a model | $10^{-14}$ |
| hy1 | inviscid Taylor-Green: dissipated $K$, hybrid against WENO | 47 times less on FNL |
| wm1 | laminar channel, isothermal walls, Re 100: force against $3/Re$; one rank against 2 | order 1.99, $4.8\cdot10^{-4}$ at 64 across; $4\cdot10^{-11}$ |
| wm2 | `test_wall_model`; laminar channel on modelled walls | see [numerics](../../numerics#wall-model); order 1.84, 1.9 % at 64 |
| wm3 | half the wall refined, 2:1: force between the coarse and the fine uniform runs | 0.03 % from their mean, mass $10^{-16}$ |
| wm4 | hot wall, modelled against resolved: energy change; one rank against 2 on modelled walls | 7.8 % (bound 10 %); $4\cdot10^{-11}$ |
| wm5 | matching height 0.3, laminar channel at Re 5 against the exact slip solution $3/Re/(1+3h^2/2)$; semi-local against classical damping | $6.4\cdot10^{-3}$, $2.2\cdot10^{-3}$, $3.0\cdot10^{-5}$ at 16, 32, 64; $4\cdot10^{-7}$ |
| wm6 | semi-local against classical damping on a hot wall, Re 100 | force 0.7 % apart (min 0.01 %) |
| wm7 | input filter: $T_f\to0$ against no filter; $T_f=0.05$ near the steady state; $T_f=1$ over a transient; one rank against 2 | 0 (bit for bit); $3.0\cdot10^{-4}$ (bound $10^{-3}$); $1.4\cdot10^{-2}$ (min $10^{-3}$); $4\cdot10^{-11}$ (CPU), $10^{-13}$ (FNL); CPU against FNL $5\cdot10^{-11}$ |
| wm8 | inverse-transformation model: laminar channel against the equilibrium model; hot wall; one rank against 2 with the filter | $1.0\cdot10^{-4}$, the heat flux of the Reynolds analogy moving the energy by $2.3\cdot10^{-3}$ (bound $10^{-4}$); $9.8\cdot10^{-3}$ (min $10^{-3}$); $4\cdot10^{-11}$ (CPU), 0 (FNL); CPU against FNL $4\cdot10^{-11}$ |
| wm9 | restart with the filter averages (`.wmf`, one per rank): 10 steps + restart to 20 against 20 straight, inverse-transformation model, $T_f=0.01$; control with the `.wmf` removed | fields and histories bit for bit (CPU and FNL); the control differs and warns |
| wm10 | dimensional units on a Vreman channel with the wall model and the forcing, gas `cp, cv`: input bit for bit (NV-5), output against `scaling.py dimensionalize` (NV-9), forcing and turbulence histories included | bit for bit (CPU and FNL); bit for bit, 4 histories |
| hy3 | the central flux: order, kinetic energy identity, pressure equilibrium | order 5.96–5.99; $4\cdot10^{-16}$; $dp/dt$ $9\cdot10^{-15}$ across 100:1 |
| hy4 | contact of density ratio 8 and 100 | $u$, $p$ within $9\cdot10^{-14}$ |
| hy5 | 2:1 box, conservation; passes with any skin | mass $10^{-16}$ |
| hy6 | Sod against WENO | L1 $8.3\cdot10^{-5}$ |
| hy7 | runtime regridding; a forest with open ends | $10^{-16}$; $5.7\cdot10^{-10}$, the single realm's drift |
| hy8 | density wave, amplitude 0.2; WENO holds it too, so this does not separate them | $p$, $u$ $10^{-14}$ |
| hy9 | model only off flagged cells; one rank against 2; conservation (does not prove the mask on ghost cells) | $\varepsilon_\mathrm{sgs}$ exactly 0 with all flagged; $7\cdot10^{-14}$; $2.3\cdot10^{-16}$ with a model on a 2:1 box; whether a flagged cell faces a free one at an interface is not checked |
| hy11 | relaxing the skin changes the density: seams, open ends (does not catch a shifted band) | $1.6\cdot10^{-4}$, $6.4\cdot10^{-5}$ |
| hy12 | restart on a refined grid, every cell | 0 |
| hy13 | against WENO, L1: a vortex through a shock (2-D), Shu-Osher and a shock through 2:1 seams (1-D, no vorticity, so the Ducros ratio flags any varying flow): robustness, no accuracy claim | $1.0\cdot10^{-4}$, $1.7\cdot10^{-5}$, $3.4\cdot10^{-5}$ |

## Taylor-Green vortex: the ILES baseline

Re 1600, $[0, 2\pi]^3$, `runge-kutta-ssp-54`, CFL 0.5, `weights = si`, 4 ranks on A30, from `turbulence/tgv_campaign.sh`. `weno` is the flux splitting, `HLLC`
`weno-riemann` with the 6th-order correction and the WENO sensor. All values per unit volume.

![The ILES baseline on the Taylor-Green vortex](/flume/tgv-iles.png)

*Top: M 0.1, WENO-5 by resolution; middle: M 0.1 at $128^3$ by scheme and order; bottom: $M_0 = 1.25$ by resolution. Left $-dK/dt$, right the enstrophy; dashed: the 512³ spectral DNS, at M 0.1.*

| Case | Peak $-dK/dt$ | at $t$ | Peak $Z$ | at $t$ | Numerical share to the peak | $K(20)$ |
|---|---|---|---|---|---|---|
| M 0.1, WENO-5, $64^3$ | 0.01377 | 5.08 | 1.71 | 6.25 | 77 % | 0.0245 |
| M 0.1, WENO-5, $128^3$ | 0.01242 | 8.24 | 3.40 | 8.11 | 58 % | 0.0219 |
| M 0.1, HLLC, $64^3$ | 0.01305 | 5.42 | 1.83 | 6.56 | 76 % | 0.0284 |
| M 0.1, HLLC, $128^3$ | 0.01249 | 8.28 | 3.59 | 8.11 | 55 % | 0.0227 |
| M 0.1, WENO-7, $128^3$ | 0.01224 | 8.19 | 3.17 | 8.14 | 60 % | 0.0181 |
| M 0.1, WENO-5, $256^3$ | 0.01272 | 8.75 | 6.18 | 8.56 | 27 % | 0.0223 |
| M 0.1, WENO-9, $128^3$ | 0.01299 | 8.28 | 4.45 | 8.03 | 43 % | 0.0213 |
| $M_0$ 1.25, WENO-5, $64^3$ | 0.01383 | 7.12 | 1.81 | 7.62 | 80 % | 0.0206 |
| $M_0$ 1.25, WENO-5, $128^3$ | 0.01081 | 8.62 | 3.34 | 9.09 | 46 % | 0.0243 |
| $M_0$ 1.25, WENO-5, $256^3$ | 0.01101 | 11.12 | 5.95 | 11.94 | 20 % | 0.0245 |
| $M_0$ 1.25, HLLC, $64^3$ | 0.01246 | 7.22 | 2.29 | 7.72 | 69 % | 0.0299 |
| $M_0$ 1.25, HLLC, $128^3$ | 0.01120 | 8.56 | 3.86 | 11.44 | 36 % | 0.0265 |
| $M_0$ 1.25, HLLC, $256^3$ | 0.01122 | 11.11 | 6.31 | 11.99 | 14 % | 0.0243 |

The share is that of the energy lost until the peak of $-dK/dt$ not accounted for by $\varepsilon - p_d$.

**Reference.** The 512³ pseudo-spectral DNS of the High-Order CFD Workshop (case C3.5, `turbulence/reference/spectral_Re1600_512.gdiag`), incompressible, so compared at M 0.1 only: peak $-dK/dt$ 0.01286 and peak $Z$ 10.29, both at $t = 8.97$. The table in [Against the DNS](#against-the-dns) gives the errors. Compressible Mach 1.25 is compared with nothing here.

**What it shows.**
- Numerical dissipation dominates: up to the peak of $-dK/dt$ it is 77 % of the energy lost at $64^3$ and 58 % at $128^3$ with WENO-5 at M 0.1. At $128^3$ the enstrophy peak is a third of the
  reference, so the resolved gradients are far from converged even where the total $-dK/dt$ looks right. At $64^3$ the peak of $-dK/dt$ comes at $t = 5$.
- Scheme and order matter less than resolution at $128^3$: HLLC gains a few percent, WENO-9 reaches $Z$ 4.45. **WENO-7 is below WENO-5** ($Z$ 3.17 against 3.40), which is not explained.
- At $M_0 = 1.25$ the numerical share falls with the resolution, 80 %, 46 %, 20 %.

## Subgrid-scale models (T1) and the go/no-go

![The models and the 2:1 box](/flume/tgv-models.png)

*Top: M 0.1 at $64^3$; middle: at $128^3$; bottom: a 2:1 box on the centre of the $64^3$ grid against the uniform grids.*

| Case, M 0.1 | Peak $-dK/dt$ | at $t$ | Peak $Z$ | Numerical share to the peak | $K(20)$ |
|---|---|---|---|---|---|
| no model, $64^3$ | 0.01377 | 5.08 | 1.71 | 77 % | 0.0245 |
| Vreman, $64^3$ | 0.01385 | 5.05 | 1.59 | 79 % | 0.0236 |
| WALE, $64^3$ | 0.01374 | 5.04 | 1.60 | 79 % | 0.0239 |
| Sigma, $64^3$ | 0.01380 | 5.17 | 1.59 | 79 % | 0.0237 |
| no model, $128^3$ | 0.01242 | 8.24 | 3.40 | 58 % | 0.0219 |
| Vreman, $128^3$ | 0.01210 | 8.29 | 3.18 | 61 % | 0.0213 |
| WALE, $128^3$ | 0.01199 | 8.29 | 3.16 | 61 % | 0.0210 |
| Sigma, $128^3$ | 0.01206 | 8.30 | 3.18 | 60 % | 0.0213 |

Here the share leaves out the model's own dissipation, so it is not a pure numerical share: $\varepsilon_\mathrm{sgs}$ is not subtracted. At $M_0$ 1.25, $128^3$, Vreman gives peak $-dK/dt$ 0.01085 and peak $Z$ 3.10, against 0.01081 and 3.34.

**No model beats no model on this measure.** On a dissipative scheme each one lowers the enstrophy peak further from the reference (1.59 against 1.71 at $64^3$, 3.18 against 3.40 at $128^3$), and the three are nearly the same.
The eddy viscosity adds to a numerical dissipation that already dominates; in the laminar phase, to $t = 2$ at M 0.4, it takes 21–53 % of the energy lost. Under the rule of #119 this makes the hybrid central/WENO scheme (T2) necessary before any model is
recommended. The constants come from other discretizations and were not tuned. A spectral comparison with filtered DNS could still rank the models differently; on the integrals, [the DNS](#against-the-dns) confirms this.

## The effect of a 2:1 box (T0.6)

M 0.1, base $64^3$, the box of side $\pi$ on the centre refined to the $128^3$ spacing, so one eighth of the volume. The first attempt never refined (`amr_iterations = 0`), found because the runs equal the uniform grid to every digit.

| Case | Peak $-dK/dt$ | at $t$ | Peak $Z$ | at $t$ | $K(20)$ | $\max\lvert p_d\rvert$ |
|---|---|---|---|---|---|---|
| uniform $64^3$ | 0.01377 | 5.08 | 1.71 | 6.25 | 0.0245 | 0.00026 |
| box, tricubic fill, mean restriction | 0.01312 | 6.39 | 1.87 | 6.39 | 0.0233 | 0.0070 |
| box, quintic fill, `lagrange6` | 0.01314 | 6.39 | 1.85 | 6.38 | 0.0234 | 0.0115 |
| uniform $128^3$ | 0.01242 | 8.24 | 3.40 | 8.11 | 0.0219 | 0.00007 |

The box delays the peak of $-dK/dt$ from $t = 5.1$ to $6.4$ and lowers it, and the two seam operators give nearly the same integrals: $-dK/dt$ and $Z$ do not tell them apart.

**$p_d$ cannot be read on these grids.** It is 27 to 44 times the uniform value on the box. The likely reason, not verified, is the diagnostic: at M 0.1, $p \approx 71$, and a central difference across a 2:1 face does not telescope, so $\int p\,\theta\,dV$ picks up the mean pressure times a nonzero discrete $\int\theta\,dV$. So the budget and the numerical share are not reported on a 2:1 grid.

## Decaying compressible isotropic turbulence (T0.5)

$M_t = 0.6$, $Re_\lambda = 100$, $k_0 = 4$, $\tau = \lambda_0/u_\mathrm{rms} = 0.5$, power-law viscosity 0.76, Pr 0.7, one seed. The DNS of Johnsen et al. (2010) is not available to this tree: only self-convergence is measured.

![Isotropic turbulence](/flume/hit.png)

*Left: the enstrophy against $t/\tau$; right: the spectrum of $\sqrt\rho\,\mathbf{u}$ at $t/\tau = 4$.*

| Case | Peak $Z$ | at $t/\tau$ | $K(t/\tau = 4)$ | $E(16)$ at $t/\tau = 4$ | $E(30)$ |
|---|---|---|---|---|---|
| no model, $64^3$ | 34.5 | 0.66 | 0.303 | $7.4\cdot10^{-4}$ | $8.3\cdot10^{-6}$ |
| Vreman, $64^3$ | 32.6 | 0.58 | 0.292 | $5.9\cdot10^{-4}$ | $5.8\cdot10^{-6}$ |
| no model, $128^3$ | 51.0 | 1.14 | 0.384 | $3.0\cdot10^{-3}$ | $1.0\cdot10^{-4}$ |
| Vreman, $128^3$ | 48.9 | 1.12 | 0.376 | $2.7\cdot10^{-3}$ | $8.1\cdot10^{-5}$ |
| no model, $256^3$ | 60.5 | 1.36 | 0.413 | $4.2\cdot10^{-3}$ | $2.3\cdot10^{-4}$ |

$K$ keeps changing with the resolution by 7 % from $128^3$ to $256^3$, so $256^3$ with WENO-5 is not converged either. As on the Taylor-Green vortex, Vreman moves away from the finer runs at each resolution.

## The hybrid central/WENO scheme (T2)

`scheme_space = hybrid-kep-weno` ([numerics](../../numerics#hybrid-central-weno)), the same inputs.

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

The flicker is the rms of $-dK/dt$ about its running mean, over the smoothed peak. The share with Vreman is left out: it does not separate the model from the scheme. The model acts off flagged cells only (T2.6); against the same runs with the model everywhere, $Z$ moves under 1 % at M 0.1, where the sensor rarely fires.

| Case, $M_0 = 1.25$ | Peak $-dK/dt$ | Peak $Z$ | Numerical share |
|---|---|---|---|
| WENO-5, $64^3$ / $128^3$ / $256^3$ | 0.01383 / 0.01081 / 0.01101 | 1.81 / 3.34 / 5.95 | 80 % / 46 % / 20 % |
| hybrid, $64^3$ / $128^3$ / $256^3$ | 0.01266 / 0.01057 / 0.01112 | 2.20 / 4.83 / 7.21 | 73 % / 36 % / 10 % |

Isotropic turbulence: peak $Z$ 41.0 and 57.4 at $64^3$, $128^3$ against 34.5 and 51.0 on WENO-5 (60.5 at $256^3$), 37.8 and 54.1 with Vreman; $K(t/\tau=4)$ 0.367 and 0.408 against 0.303 and 0.384 (0.413 at $256^3$).

**What it shows.**
- At M 0.1 the hybrid reaches at $128^3$ a higher enstrophy peak than WENO-5 at $256^3$, and at $64^3$ more than WENO-5 at $128^3$, with the dissipation peak where the finer runs have it. $Z$ is still below the quoted 10.3.
- At $M_0 = 1.25$ the gain is smaller, as the sensor keeps WENO around the shocklets.
- **The raw peak of $-dK/dt$ is not a measure on the hybrid:** it flickers, 7 % at $64^3$, 2 % at $128^3$. The source is not established; the sensor switching on and off is the first candidate. Smoothed, the hybrid is 4–6 % above WENO-5 at $256^3$.
- **Vreman on the hybrid** removes half of the flicker and lowers $Z$ by 20 %. At $64^3$ it lands below the finer runs (0.0117). Against the DNS it makes the hybrid worse at both resolutions (below).
- The thresholds (0.2, 0.1) and the floor are not calibrated.

### On a 2:1 box, and what the skin costs

![The hybrid scheme on a 2:1 box](/flume/tgv-hybrid-box.png)

*M 0.1, a box of half the side refined once on $64^3$. Left $-dK/dt$, right the enstrophy.*

| Case | Smoothed peak $-dK/dt$ | Peak $Z$ |
|---|---|---|
| WENO-5, box | 0.01304 | 1.87 |
| hybrid, box, WENO skin (default) | 0.01239 | 3.69 |
| hybrid, box, central skin | 0.01319 | 4.13 |
| hybrid, uniform $64^3$ | 0.01333 | 3.82 |

With the WENO skin the box ends below the uniform run in $Z$: the skin takes more than the extra cells give. The central skin is stable to $t = 20$ and above it. Against the DNS the central skin is the closer one (table below), on this smooth flow; the case with a shock has no reference. A vortex after Jiang & Shu (the code's profile, same strength and radius, unit square) through a Mach 1.1 shock crossing the box seams also runs with the central skin; it differs from the WENO skin by at most $2.8\cdot10^{-4}$ in density at $t = 0.35$.

### Against the DNS

M 0.1 against the 512³ DNS, from `turbulence/tgv_dns.py`. Peak $-dK/dt$ over a window of 0.5; L2 errors over $0.3 < t < 19.7$ relative to the DNS, $-dK/dt$ over the same window.

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

- **The hybrid at $128^3$ is closer than WENO-5 at $256^3$** on all three errors, with one eighth of the cells. At $64^3$ it is closer than WENO-5 at $128^3$.
- Its peak of $-dK/dt$ is 2.4 % above the DNS and comes at $t = 8.47$ rather than 8.97; $Z$ at $128^3$ reaches two thirds of the DNS peak, so the gradients are not converged.
- **A model makes it worse:** Vreman takes the $-dK/dt$ error from 0.067 to 0.114 and the $Z$ error from 0.27 to 0.40 at $128^3$. All three models also worsen WENO-5. This flow at Re 1600 is not a case where the subgrid scales are far from resolved, and the constants were not tuned; on these integrals no model is recommended here.
- **WENO-9** at $128^3$ (0.177) is much closer than WENO-5 (0.278); **WENO-7** is further (0.291), still unexplained.
- **The skin:** on the 2:1 box the WENO skin (0.269) is further than uniform $64^3$ (0.237); the central skin ties it on $-dK/dt$ (0.236) and beats it on $Z$ (0.580 against 0.615). Smooth flow only.

### Cost

200 steps on $128^3$, two A30, three samples: WENO 24.1–24.6 s, every cell flagged 24.2–24.7 s, every cell central 26.7–26.9 s. The sensor costs under 1 %, the central flux about 10 % on top of WENO, startup included. A single kernel with a mask was not tried, so what skipping WENO would save is not measured.

## Plane channel and the wall model (T3)

Coleman, Kim & Moser (1995): bulk Mach 1.5, bulk Re 3000 on the half height and the wall viscosity, isothermal walls, viscosity $(T/T_w)^{0.7}$,
Pr 0.7; $4\pi\times2\times\frac43\pi$, mass-flow forcing; $Re_\tau$ about 220–222 for the DNS, as cited by Qi et al. (2022), Cao et al. (2024),
Alberti et al. (2025b). The hybrid scheme throughout. From `turbulence/channel_campaign.sh` with `channel_cases_t35.txt`, 4 A30, after the no-slip
leak fix. $\tau_w$ is the mean force per unit volume times $h$; $Re_\tau=3000\sqrt{\tau_w\rho_w}$ with $\rho_w = 1.30$ for every row, so the rows
compare among themselves; the reference run's own wall density, extrapolated from its last checkpoint, is 1.25, which lowers every $Re_\tau$ by 2 %.
$u_c$ and $T_c/T_w$ are the Favre plane means at the centreline of the last checkpoint (one snapshot).

![The plane channel](/flume/channel-wm.png)

*$\tau_w$ against time; dashed: $Re_\tau$ 220 and 222; dotted: the laminar value.*

| Case | Cells | $\tau_w$, $100<t\le150$ | $150<t\le200$ | $Re_\tau$ | $u_c$, $T_c/T_w$ |
|---|---|---|---|---|---|
| resolved walls | $256\times192\times128$, $\Delta y^+\approx2.3$ | 3.82e-3 | – | 211 ($50<t\le100$: 218) | 1.16, 1.26 |
| resolved walls | $64\times32\times32$ | 2.35e-3 | 2.72e-3 | 166, 178 | 1.22, 1.26 |
| equilibrium, 3rd cell | $64\times32\times32$ | 4.37e-3 | 4.49e-3 | 226, 229 | 1.21, 1.29 |
| equilibrium, 3rd cell | $96\times48\times48$ | 4.36e-3 | 4.52e-3 | 226, 230 | 1.19, 1.30 |
| equilibrium, 3rd cell | $128\times64\times64$ | 4.69e-3 | 4.66e-3 | 234, 234 | 1.18, 1.31 |
| equilibrium, 3rd cell, semi-local | $128\times64\times64$ | 4.33e-3 | 4.33e-3 | 225, 225 | 1.18, 1.30 |
| equilibrium, $h=0.15$ | $96\times48\times48$ | 4.71e-3 | 4.77e-3 | 235, 236 | 1.19, 1.30 |
| equilibrium, $h=0.15$ | $128\times64\times64$ | 4.89e-3 | 4.90e-3 | 239, 239 | 1.19, 1.30 |
| equilibrium, $h=0.15$, semi-local | $96\times48\times48$ | 4.15e-3 | 4.09e-3 | 220, 219 | 1.18, 1.29 |
| equilibrium, $h=0.15$, semi-local | $128\times64\times64$ | 4.13e-3 | 4.21e-3 | 220, 222 | 1.17, 1.29 |
| … and the filter $T_f=1$ | $96\times48\times48$ | 4.15e-3 | 4.03e-3 | 220, 217 | 1.18, 1.29 |
| … and the filter $T_f=1$ | $128\times64\times64$ | 4.15e-3 | 4.16e-3 | 220, 221 | 1.19, 1.29 |
| inverse, $h=0.15$, $u_e=1.18$ | $128\times64\times64$ | 4.13e-3 | 4.08e-3 | 220, 218 | 1.16, 1.21 |
| … and the filter $T_f=1$ | $128\times64\times64$ | 4.04e-3 | 4.11e-3 | 217, 219 | 1.17, 1.21 |
| inverse, $h=0.15$, $u_e=1.18$ | $96\times48\times48$ | 2.03e-3 | 1.62e-3 | 154, 138: relaminarises | 1.51, 1.11 |
| inverse, 3rd cell, $u_e=1.18$ | $96\times48\times48$ | 1.61e-3 | 1.40e-3 | 137, 128: relaminarises | 1.44, 1.12 |
| inverse, $h=0.15$, $u_e=1.10$ | $96\times48\times48$ | 3.61e-3 | 4.15e-3 | 206, 220 | 1.18, 1.21 |

The inverse model's edge state is given ($u_e$, $T_e/T_w = 1.29$): the centreline of the $h=0.15$ semi-local runs. Every run conserves mass to $10^{-16}$
(the wall-resolved ones gained 11 % and 1.3 % before the leak fix).

**What it shows.**
- **Matching height 0.15 with semi-local damping puts the equilibrium model in the DNS band on both grids** (219–222). Either change alone is worse:
  the 3rd-cell sample is 226–234 and rises with the grid, $h=0.15$ with the classical damping 235–239, the semi-local damping at the 3rd cell 225.
  The fixed height removes the grid trend; the semi-local damping removes the excess. This is one Mach and Reynolds number, single runs, and the band
  sits within the $\pm2\%$ of the wall density.
- **The input filter does not move the means** (within 1.5 %); on the coarser grid it slows the transient.
- **The inverse model is in the band on the finer grid** with a cooler centreline ($T_c/T_w$ 1.21 against 1.29 for the equilibrium model, 1.26 for
  the resolved run): its heat flux differs. **On $96\times48\times48$ it relaminarises with $u_e = 1.18$**, at $h=0.15$ and at the 3rd cell, but
  not with $u_e=1.10$: the model is sensitive to the given edge velocity, which the flow does not set (the paper takes the edge from the flow; a face
  kernel cannot). Why 1.18 drains the turbulence on that grid is not explained.
- **The wall-resolved reference did not change with the leak fix** (211 against 211 before, $100<t\le150$), so its 4–5 % shortfall from the DNS was
  not the leak: its resolution or the hybrid scheme's dissipation, not tested. The modelled rows in the band agree with the DNS, not with it.
- **The coarse grid with resolved walls** still falls 20 % short: without a model the grid cannot carry the wall layer.
- Only $\tau_w$ and the centreline are compared: no averaged profiles, and the DNS profiles are not in this tree.

### Mach 3

Coleman, Kim & Moser's second case: bulk Mach 3, bulk Re 4880, $Re_\tau$ about 450 for the DNS; the rest as above, the equilibrium model only, to
$t = 200$, `channel_cases_t3_mach3.txt`. The wall density changes with the run at this Mach number, so each row takes its own: $\rho_w$ is the
first-cell mean pressure of the last checkpoint over $T_w$, and $Re_\tau = 4880\sqrt{\tau_w\rho_w}$. $T_1/T_w$ is the first cell (Favre, both walls).
$Q_w/W$ is the wall heat flux over the forcing work for $180<t\le200$, from the total-energy budget (`channel_wall_cells.py budget`; $Q_w$ is the
residual): 1 at a thermal steady state.

![The plane channel at Mach 3](/flume/channel-m3.png)

*Left: $\tau_w$ against time (dotted: the laminar value). Middle and right: the first cell of the two $h=0.15$ runs on $128\times64\times64$ against
time, from a checkpoint every 500 steps (`channel_cases_t3_mach3_drain.txt`, the same runs, bitwise).*

| Case | Cells | $\tau_w$, $100<t\le150$ | $150<t\le200$ | $\rho_w$ | $Re_\tau$ | $u_c$, $T_c/T_w$ | $T_1/T_w$ | $Q_w/W$ |
|---|---|---|---|---|---|---|---|---|
| equilibrium, 3rd cell | $96\times48\times48$ | 1.77e-3 | 1.43e-3 | 1.65 | 264, 238: relaminarises | 1.43, 1.63 | 0.86 | 0.71 |
| equilibrium, $h=0.15$, semi-local | $96\times48\times48$ | 1.89e-3 | 1.59e-3 | 1.54 | 263, 241: relaminarises | 1.50, 1.53 | 0.64 | 0.67 |
| equilibrium, 3rd cell | $128\times64\times64$ | 4.29e-3 | 4.44e-3 | 2.12 | 465, 473 | 1.25, 2.22 | 1.03 | 0.87 |
| equilibrium, 3rd cell, semi-local | $128\times64\times64$ | 3.20e-3 | 3.03e-3 | 2.02 | 392, 381 | 1.20, 2.04 | 1.28 | 0.85 |
| equilibrium, $h=0.15$, semi-local | $128\times64\times64$ | 3.58e-3 | 3.46e-3 | 1.91 | 404, 397 | 1.28, 2.01 | 0.81 | 0.83 |
| equilibrium, $h=0.15$ | $128\times64\times64$ | 6.46e-3 | 7.49e-3 | 2.12 | 571, 615: still rising | 1.48, 2.63 | 0.41 | 0.88 |

**What it shows.**
- **No run is at a thermal steady state at $t=200$**: the total energy still rises in every run ($Q_w/W$ 0.67–0.88), so these $Re_\tau$ are values
  on the way, not means.
- **The Mach 1.5 setting does not carry over**: $h=0.15$ with semi-local damping is 12 % below the DNS; the 3rd-cell sample is 5 % above. The coarse
  grid relaminarises with either setting.
- **The first cell is drained.** The wall face takes the modelled $\tau_w$ and $q_w$ from the sample at $h$; the first cell's own state never enters,
  and nothing below $h$ carries those fluxes (no subgrid model): the molecular plus resolved stress at the face between cells 1 and 2 is 0.54
  $\tau_w$ with semi-local damping and 0.21–0.39 with the classical one (after $t=30$). With $h=0.15$ and the classical damping the first cell falls to $0.26\,T_w$,
  reverses its flow and holds 17 % of the mass at $t\approx80$; $\tau_w$ follows it, also down ($t$ 80–110) and up again: the run does not settle.
  With semi-local damping the drain is slower and $\tau_w$ stays flat, but the first cell is still below $T_w$ from $t\approx125$. A first cell colder
  than an isothermal wall that takes heat out of it is not a boundary layer: the treatment below $h$, not the model's law, is what fails here.

## Not done

- HLLC at $256^3$, M 0.1 was dropped by owner decision.
- Spectra against the DNS: the file holds integrals only. The DNS of Johnsen et al. for the isotropic case.
- Calibration of the sensor; a model masked by face rather than by cell; `central-bcs` on a forest seam.
- T3: Mach 3 to a thermal steady state, with a flux path below $h$ (the drain above); a compressible boundary layer and SBLI (it needs a turbulent inflow); runtime regridding along a modelled wall; mean profiles
  against the DNS; the relaminarisation of the inverse model on $96\times48\times48$ with $u_e = 1.18$; an edge state taken from the flow.
- A model on 2:1 grids beyond conservation; runtime AMR with a model; a wall.
