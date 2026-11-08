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
| tt4 | `spectra.py` on the initial field: Parseval, the energy of the single-mode field in shell 2 | not recorded (bound $10^{-12}$) |
| tt6 | pressure work, M 0.8: the budget closes only with $p_d$ | 8.5 % left, 203 % without, 414 % flipped |
| ht1 | `isotropic-turbulence` field: every shell, $K$, dilatation | shells $7\cdot10^{-15}$, $K$ $1.5\cdot10^{-14}$ |
| ht2 | one rank against 2 on the random field | $7\cdot10^{-15}$ |
| sg1 | the formulas against numpy on 200 tensors; the zeros of each model | $3.6\cdot10^{-15}$, $1.2\cdot10^{-15}$, $1.8\cdot10^{-11}$ |
| sg2 | budget with the model, M 0.4: Vreman 36 %, 5.7 %, 0.27 % left at $32^3$, $64^3$, $128^3$ | at $128^3$: 21.6 % without $\varepsilon_\mathrm{sgs}$, 43 % flipped |
| sg3 | no viscosity with each model | $\varepsilon = 0$, mass and energy $3.3\cdot10^{-16}$ |
| sg4 | one rank against 2 with a model | $10^{-14}$ |
| hy1 | inviscid Taylor-Green: dissipated $K$, hybrid against WENO | 47 times less on FNL |
| hy2 | one rank against 2 with the hybrid scheme on the random field, where the sensor fires | not recorded |
| hy3 | the central flux: order, kinetic energy identity, pressure equilibrium | order 5.96–5.99; $4\cdot10^{-16}$; $dp/dt$ $9\cdot10^{-15}$ across 100:1 |
| hy4 | contact of density ratio 8 and 100 | $u$, $p$ within $9\cdot10^{-14}$ |
| hy5 | 2:1 box, conservation; passes with any skin | mass $10^{-16}$ |
| hy6 | Sod against WENO | L1 $8.3\cdot10^{-5}$ |
| hy7 | runtime regridding; a forest with open ends | $10^{-16}$; $5.7\cdot10^{-10}$, the single realm's drift |
| hy8 | density wave, amplitude 0.2; WENO holds it too, so this does not separate them | $p$, $u$ $10^{-14}$ |
| hy9 | model only off flagged cells; one rank against 2; conservation (does not prove the mask on ghost cells) | $\varepsilon_\mathrm{sgs}$ exactly 0 with all flagged; $7\cdot10^{-14}$; $2.3\cdot10^{-16}$ with a model on a 2:1 box; whether a flagged cell faces a free one at an interface is not checked |
| hy10 | the skin relaxed on 2:1 seams (the hy5 case, `central-seams`): what conservation does, the last history rows | a measure with no bound: printed, never fails the script |
| hy11 | relaxing the skin changes the density: seams, open ends (does not catch a shifted band) | $1.6\cdot10^{-4}$, $6.4\cdot10^{-5}$ |
| hy12 | restart on a refined grid, every cell | 0 |
| hy13 | against WENO, L1: a vortex through a shock (2-D), Shu-Osher and a shock through 2:1 seams (1-D, no vorticity, so the Ducros ratio flags any varying flow): robustness, no accuracy claim | $1.0\cdot10^{-4}$, $1.7\cdot10^{-5}$, $3.4\cdot10^{-5}$ |
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
| wm11 | the wall model's eddy viscosity below the matching height (`wall_model_eddy_viscosity`), perturbed channel at Re 1000, $h = 0.3$: it acts; 2 against 4 blocks per direction; one rank against 2; the same with Vreman | not recorded (bounds: force change above $10^{-3}$; blocks within `channel_oracle agree`'s bounds, $2.4\cdot10^{-6}$ against $10^{-9}$ with the ghost columns left out) |
| tb1 | digital-filter inflow: the generator's stresses, temporal and spanwise correlations, periodic span, plane restart (`test_inflow_df`); a Mach 2.5 box: the inflow ghosts carry the table; one rank against 2; restart with the `.dff` against straight, control without it | stresses 0.4–0.7 % (bound 3 %), correlation exact, integral length 5 % short (bound 10 %); $\langle u\rangle$ at $\delta/2$ 0.877 against 0.891; $2\cdot10^{-16}$; fields and histories bit for bit (CPU and FNL), the control differs; FNL against CPU $10^{-15}$ |
| tb2 | wall-model history: $\tau_w$ of wm5's laminar channel against the driving force; one rank against 2 on the Mach 2.5 box; CPU against GPU | $3.9\cdot10^{-4}$ (bound $10^{-3}$); $4\cdot10^{-16}$; $1.1\cdot10^{-13}$ |
| tt5 | (`--agree TAG`) this executable's tt2 history against the executable TAG: backends | not recorded (bound $10^{-11}$) |

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

### The flow

![Vorticity of the Taylor-Green vortex and of the isotropic turbulence](/flume/tgv-fields.png)

*`hybrid-kep-weno` at $128^3$, re-run for this figure (`tgv_cases_vis.txt`, medooza dev-gn). Taylor-Green, Re 1600, Mach 0.1: $\omega_x$ on the
periodic face $x\approx0$ at $t=4.2$ and $9.4$ (the peak of the dissipation is at $t\approx9$; the Workshop's figures show the vorticity norm on
this face) and $\omega_z$ on the symmetry plane $z=\pi$ at $t=9.4$. Isotropic turbulence, $M_t$ 0.6, $Re_\lambda$ 100, $k_0=4$: $\omega_z$ on
$z=\pi$ at $t=2$. Planes from `turbulence/field_slices.py --time`; `make_doc_figures.py --only tgv-fields`.*

- **At $t=4.2$ the vorticity is still in sheets; by $t=9.4$ they have rolled up and broken into the small structures of the transition.**
- **The flow's mirror symmetries survive the breakdown**: on the face $x\approx0$, $\omega_x$ is antisymmetric about $y=\pi$ and about $z=\pi$ to
  $8\cdot10^{-13}$ of its maximum at $t=9.4$ (0 at $t=4.2$), as the initial condition is.
- **On $z=\pi$ the flow is still large-scale at $t=9.4$**: four counter-rotating cells with thin layers at their edges.
- **The isotropic field shows no visible imprint of the blocks** ($2^2$ per direction). One plane, one time: a description, not a measurement.

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

**Result first:** with the wall model's eddy viscosity below $h$ (`wall_model_eddy_viscosity`, below) and run to a thermal steady state at $t=800$,
$h=0.15$ with semi-local damping gives $Re_\tau = 443$–$445$ against about 451 for the DNS; the classical damping overshoots by 13–21 %. The runs
without it, to $t=200$ only, come first: they show why the key is needed, and their $Re_\tau$ are transient values, not means.

#### Without a flux path below $h$, to $t = 200$

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

#### With the eddy viscosity below $h$, to $t = 800$

The face key `wall_model_eddy_viscosity = .true.` sets $\mu_t = \max(\mu_{t,\mathrm{SGS}},\ \kappa y\sqrt{\rho\tau_w}\,D^2)$ in the cells whose centre is below
$h$ ([numerics](../../numerics.md)): the flux path that the drain above lacked. The three runs on $128\times64\times64$ are the inputs above plus the key, to
$t=200$ (`channel_cases_t3_mach3_ev.txt`), then continued from the restart to $t=800$ (`channel_cases_t3_mach3_ev800.txt`). $\rho_w$ and $T_1/T_w$ are
taken at the last checkpoint as above; $Q_w/W$ over $700<t\le800$.

![The plane channel at Mach 3 with the eddy viscosity below h](/flume/channel-m3-ev.png)

*Left: $Re_\tau$ against time, with the wall density of each checkpoint (dashed: the DNS). Middle: the first cell. Right: the wall heat flux over the
forcing work in 50-unit windows, 1 at a thermal steady state.*

| Case | $\tau_w$, $600<t\le700$ | $700<t\le800$ | $\rho_w$ | $Re_\tau$ | $u_c$, $T_c/T_w$ | $T_1/T_w$ | $Q_w/W$ |
|---|---|---|---|---|---|---|---|
| equilibrium, $h=0.15$, semi-local | 3.52e-3 | 3.48e-3 | 2.36 | **445, 443** | 1.17, 2.44 | 1.79 | 0.99 |
| equilibrium, 3rd cell | 4.51e-3 | 4.35e-3 | 2.51 | 519, 510 | 1.21, 2.65 | 1.57 | 0.99 |
| equilibrium, $h=0.15$ | 5.17e-3 | 5.15e-3 | 2.44 | 548, 547 | 1.21, 2.53 | 1.84 | 1.00 |

**What it shows.**
- **The drain is gone.** The first cell is $1.6$–$1.8\,T_w$ at every time after the transient, the two walls agree within 2 %, and the two first-cell
  layers hold 4–5 % of the mass (3.1 % if uniform), against 16 % at $0.41\,T_w$ without the key.
- **The runs are thermally steady at $t=800$**: $Q_w/W$ is 0.99–1.00 and $dE/dt$ about 60 times smaller than at $t=200$.
- **$h=0.15$ with semi-local damping matches the DNS within 2 %.** With the classical damping $Re_\tau$ is 13 % (3rd cell) and 21 % ($h=0.15$) above.
- **$Re_\tau$ needs the thermal steady state at this Mach number.** The bulk heats until $t\approx600$, $\rho_w$ rises from about 2.1 to 2.4–2.5, and
  $Re_\tau$ rises with it: at $t=200$ the 3rd-cell run gave 438, 448, which looked like a match.
- **Cost:** from the same $t=800$ state of the semi-local run to $t=810$ (medooza, 4 A30), the key on takes 1.69 times the wall time of the key off:
  1.28 times the steps (the eddy viscosity enters the viscous time-step limit) and 1.32 times the time per step.

### The flow

![Instantaneous streamwise velocity in the two channels](/flume/channel-fields.png)

*The streamwise velocity over the bulk velocity of the last checkpoint: the $x$-$y$ plane at mid-span (left), the wall-parallel plane at the
matching height $h=0.15$ above the lower wall (middle, its own colour range) and the cross-flow plane (right). Top: bulk Mach 1.5,
$128\times64\times64$, $h=0.15$ with semi-local damping, $t=200$. Bottom: bulk Mach 3, the same grid and model with the eddy viscosity below $h$,
$t=800$. Planes from `turbulence/field_slices.py`; the figure is `make_doc_figures.py --only channel-fields`.*

- **At the matching height the flow is in streaks** along $x$, but outer-scaled ones. At Mach 3 the first minimum of the spanwise correlation of
  $u$ in that plane puts them $1.05\,h$ apart, about 470 wall units at $Re_\tau\approx445$, against the 100 of the near-wall cycle; at Mach 1.5
  the correlation has no clear minimum within the span. The near-wall cycle lives below the matching height ($y^+\approx33$ at Mach 1.5), in the
  wall model, and the grid ($\Delta z^+\approx14$, $\Delta x^+\approx22$ at Mach 1.5) does not form it above.
- **At Mach 3 the streaks look longer and more regular** than at Mach 1.5. One snapshot each: a description, not a measurement.

## Compressible boundary layer (T3.7)

The zero-pressure-gradient layer of the DNS of Zhang, Duan & Choudhari (2018, AIAA J 56(11)), case M2p5: Mach 2.5, adiabatic wall
($T_w/T_r=1$), $Re_\theta=2835$ and $Re_\tau=510$ at its station. The station profiles come from the NASA Turbulence Modeling Resource
(`turbulence/reference/`). The layer is fed through the face `inflow-digital-filter` ([numerics](../../numerics.md)): the DNS station profile scaled
to the inflow thickness $\delta_i$, with fluctuations that have the DNS Reynolds stresses. $C_f$ is the wall model's own $\tau_w$, written by
`[diagnostics] wall_history_save` (leg tb2).

| Setup | |
|---|---|
| box, cells | $90\times6\times6\,\delta_i$, $1080\times72\times72$, cubic $\delta_i/12$; periodic span, extrapolated top and outflow |
| gas | $Re=U_\infty\delta_i/\nu_\infty=19370$ (the DNS's 4 mm $\delta_i$), $\mu\propto T^{0.76}$ (the DNS uses Sutherland), $Pr=0.71$ |
| wall | `wall-modelled`, $T_w=2.104\,T_\infty$, the 3rd cell, semi-local damping, its eddy viscosity below $h$ (the setting that matched the Mach 3 channel) |
| inflow | integral lengths $0.8$, $0.3$, $0.3\,\delta_i$ (streamwise as a time scale $L_x/U_\infty$, wall-normal, spanwise) |
| scheme | `hybrid-kep-weno`, no subgrid model |
| run | `tbl_cases_t37_l90.txt`, to $t=360\,\delta_i/U_\infty$ (4 flow-throughs); the wall history every 20 steps and 34 checkpoints of $180<t\le360$, averaged over the span and in time (`tbl_oracle.py compare`); 0.40 s per step on 4 A30 |

![The Mach 2.5 boundary layer against the DNS](/flume/tbl-m2p5.png)

*Left: $C_f$ from the wall model's $\tau_w$ (blue) and from the momentum integral $2\,d\theta/dx$ (grey), and $Re_\theta$ (right axis) against
$x$; the square is the DNS at its $Re_\theta$. Middle: $u^+_{VD}$ at the station where $Re_\theta$ is the DNS's, scaled with the run's $u_\tau$. Right:
the rms over $u_\tau$ at the same station; dashed: the DNS.*

| At $Re_\theta=2836$ ($x=55\,\delta_i$) | LES | DNS |
|---|---|---|
| $C_f$, the wall model's $\tau_w$ | **2.21e-3** (the two halves of the window 2.19e-3, 2.22e-3) | 2.31e-3 |
| $C_f$, momentum integral | 1.56e-3 | |
| $H$, against the edge state ([the shape factor](#the-shape-factor-and-the-free-stream)) | 4.08 | 4.10 (its profile: 4.22) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.62, 0.37 |

**What it shows.**
- **$C_f$ within 4.4 % of the DNS** at its $Re_\theta$, with the setting that matched the Mach 3 channel. The two halves of the averaging window
  differ by 1.1 %. Most of the 4.4 % is the inflow's recovery, not the plateau: farther downstream, at $Re_\theta=3000$–$3350$, $C_f$ is within
  1.6–0.4 % of the DNS's value scaled as $Re_\theta^{-1/4}$, and at $\delta_i/16$ it is 0.3 % at the station ([resolution](#resolution-mach-2-5-at-delta-i-16)).
- **The inflow recovers by $x\approx50\,\delta_i$.** The layer reaches the DNS's $Re_\theta$ at $x=55$, against 53 in the DNS. $C_f$ falls from
  $2.8\cdot10^{-3}$ at the inlet to $1.8\cdot10^{-3}$ at $x\approx25$, the adaptation of the synthetic eddies, then rises to a plateau of
  $2.2\cdot10^{-3}$ from $x\approx50$. That matches the $40\,\delta_i$ recovery length of a digital filter at this Mach number (Huang, Duan &
  Choudhari 2022).
- **The momentum integral does not measure $C_f$ here.** Fitted over $\pm3\,\delta_i$, $2\,d\theta/dx$ swings between $0.9$ and $4\cdot10^{-3}$
  along $x$ on the same averaged data, and is $-33$ % at the station. A first run on a $60\,\delta_i$ box had only this estimate, and it reported
  $C_f$ 16–28 % low.
- **The log layer sits 0.9–1.0 above the DNS in $u^+_{VD}$**, and $u_{rms}$ peaks too high, the near-wall streaks of a coarse LES. The resolved
  $-\langle u'v'\rangle$ in the outer layer is above the DNS's. At $\delta_i/16$ the log layer, $u_{rms}$ and $-\langle u'v'\rangle$ move toward
  the DNS. $H$, taken against the local free stream as the DNS takes it, is the DNS's within 1 % (4.08).

### The flow

![Instantaneous fields of the two boundary layers](/flume/tbl-fields.png)

*The last checkpoint of the Mach 2.5 run (top block, $90\,\delta_i$ box, $t=360$) and of the Mach 5.86 cold-wall run (bottom block). Per run: the
streamwise velocity in the $x$-$y$ plane at mid-span over the whole box, with the 99 % line of the spanwise average (black); the wall-parallel plane
at $y=0.29\,\delta_i$ (its own colour range); a $20\,\delta_i$ window of the $x$-$y$ plane and the temperature in the cross-flow plane at its centre.
The figure is `make_doc_figures.py --only tbl-fields`.*

![The Mach 2.5 boundary layer and the interaction in 3-D](/flume/tbl-3d.png)

*The same planes as the faces of a cutaway, $u/U_\infty$, $y$ drawn stretched twice: the Mach 2.5 layer over $15<x<45$, the cross-flow face at
$x=45$ (top), and the Mach 2.28 interaction of [the next section](#shock-boundary-layer-interaction-121-s4b) over $35<x<60$, the cross-flow face at
the impingement $x=60$ (bottom). `make_doc_figures.py --only tbl-3d`.*

- **The streaks at $y=0.29\,\delta_i$ ($0.18\,\delta$) are outer-scaled**: the first minimum of the spanwise correlation of $u$ puts them
  $1.5$–$1.7\,\delta_i$ apart at Mach 2.5 and $1.8$–$2.0\,\delta_i$ at Mach 5.86, about one $\delta$, the same near the inflow ($x<15$) and
  downstream ($30<x<80$). At $\delta_i/12$ the first cell is $0.05\,\delta$, about 26 wall units at $Re_\tau=510$; the near-wall cycle is the
  wall model's, as in the channel.
- **The edge of the layer is corrugated by bulges** of about $\delta$ in length. The cold-wall layer's cross-section shows the hot fluid near the
  wall (up to $5.5\,T_\infty$) lifted in ejections; at Mach 2.5 the near-wall temperature is the adiabatic $2.05\,T_\infty$.

### Resolution: Mach 2.5 at $\delta_i/16$

The Mach 2.5 case on the same $90\,\delta_i$ box at $\delta_i/16$: $1440\times96\times96$ cells, the matching point still the 3rd cell (so $h$
shrinks with the grid, as in Kawai & Larsson 2012), the filter's node spacing the new $\Delta y$; everything else as at $\delta_i/12$. 30
checkpoints of $180<t\le360$ (`tbl_cases_t37_r16.txt`, 35264 steps, 0.80 s per step on 4 A30).

![The Mach 2.5 boundary layer at two resolutions against the DNS](/flume/tbl-m2p5-res.png)

*Left: $C_f$ from the wall model's $\tau_w$ at $\delta_i/12$ (dashed) and $\delta_i/16$ (solid); the square is the DNS at its station. Middle:
$u^+_{VD}$ at the station of each run. Right: the rms and $-\langle u'v'\rangle$ over $u_\tau$; dotted: the DNS.*

| | $\delta_i/12$ | $\delta_i/16$ | DNS |
|---|---|---|---|
| $C_f$ at $Re_\theta=2835$ | 2.208e-3 (halves 2.195e-3, 2.220e-3) | **2.306e-3** (halves 2.299e-3, 2.312e-3) | 2.310e-3 |
| $C_f$ against the DNS scaled as $Re_\theta^{-1/4}$, at $Re_\theta=2600$, 2835, 3000, 3200, 3350 | −12.0, −3.9, −1.6, −0.9, −0.4 % | −3.6, −0.3, +1.0, +0.7, −0.5 % | |
| $H$, against the edge state | 4.08 | 4.23 | 4.10 (its profile: 4.22) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 12.8, 17.3, 19.2 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.79, 0.83, 1.16 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.66, 0.43 | 0.62, 0.37 |

**What it shows.**
- **$C_f$ is grid-converged once the layer has recovered.** Downstream of $x\approx60$ ($Re_\theta\ge3000$) the two grids agree with each other and
  with the DNS's $C_f$ scaled as $Re_\theta^{-1/4}$ within 1–2 %. The scaling is an approximation (the DNS file holds one station), so these
  figures carry its error too.
- **The finer grid shortens the inflow's recovery.** The $C_f$ dip after the inlet ends by $x\approx40$ at $\delta_i/16$, against $x\approx55$–$60$
  at $\delta_i/12$: the $-4.4$ % at the DNS station of the $\delta_i/12$ run is mostly the layer still recovering. At $\delta_i/16$ the station is
  past it, and $C_f$ is $-0.2$ %; the halves of the window differ by 0.6 %.
- **The near-wall and log-layer statistics converge toward the DNS**: the $u_{rms}$ peak at $y/\delta=0.2$ from 2.06 to 1.79, $-\langle
  u'v'\rangle$ from 0.80 to 0.66, the log-layer offset in $u^+_{VD}$ at $y^+=200$ from +1.0 to +0.4.
- **$H$ is the DNS's on both grids** (4.08, 4.23, against 4.10) once taken against the local free stream
  ([the shape factor](#the-shape-factor-and-the-free-stream)).

### The shape factor and the free stream

The DNS's profiles are normalised by its free stream at the station, and are 1 there. The run's free stream is not the nominal one: the layer
thickens downstream and turns the supersonic stream over it, a weak compression of order $\gamma M^2 (d\delta^*/dx)/\sqrt{M^2-1}$ (linear
supersonic theory). `tbl_oracle.py compare` takes $\delta^*$ and $\theta$ against the local edge state: $\rho_e$ and the Favre $u_e$ the means over
$1.5$–$2.5\,\delta$, the integrals to $1.5\,\delta$, the DNS's own profile integrated the same way (4.22 at Mach 2.5, against the header's
rounded 4.1). Earlier versions of this page took them against the nominal $\rho_\infty$, $U_\infty$ and integrated to the top of the box. That
reported $H$ 14 % low at Mach 2.5 and 30 % low at Mach 5.86, and the runs below set out to explain a deficit that was the measurement.

![H against the edge state and against the nominal free stream](/flume/tbl-edge.png)

*Left and right: $H$ of the Mach 2.5 baseline and the Mach 5.86 run against $x$, against the edge state (solid) and the nominal free stream
(dashed); the DNS's header and its profile's at its station. Middle: the edge pressure and density over the nominal ones against $x$.*

- **The free stream is compressed.** At Mach 2.5 $p_e/p_\infty$ rises to 1.027 at $x\approx15$ and is 1.013 at the station, $\rho_e/\rho_\infty$
  1.009; at Mach 5.86 1.08 and 1.042. Integrated to the top of the box ($6\,\delta_i$, more than $3\,\delta$), the denser stream over the layer
  takes 0.5 off $H$ at Mach 2.5 and 4 at Mach 5.86. The dip of the nominal $H$ to 2.8 at $x\approx15$ is the compression's peak, not the
  layer.
- **Against the edge state $H$ is the DNS's**: 4.08 at Mach 2.5 (4.10, its profile 4.22), 13.4 at Mach 5.86 (13.6, its profile 14.9). Taken
  over $2.5$–$3.5$, $3$–$5$ or $4$–$5.5\,\delta_i$, the edge moves it by 0.04 at Mach 2.5.

**Five runs at $\delta_i/12$.** Before the measurement was found, the Mach 2.5 case ran with one setting changed each, to look for the deficit:
the filter's integral lengths raised from $0.8, 0.3, 0.3$ to $1.2, 0.5, 0.5\,\delta_i$ and the span doubled to $12\,\delta_i$
($1080\times72\times144$ cells; `tbl_cases_t37_hdef.txt`), and an algebraic subgrid model, `[turbulence] sgs_model = wale` and `vreman`, their
default constants (`tbl_cases_t37_sgs.txt`; WALE on dev-gn, Vreman on prod-gn). Same averaging, 33–35 checkpoints of $180<t\le360$. They measure
how much the other statistics depend on these settings.

![Runs of the Mach 2.5 layer](/flume/tbl-m2p5-hdef.png)

*Four of the runs (the span of $12\,\delta_i$ in the table only). Left: $H$ against the edge state, against $x$; the square is the DNS at its station. Middle: $u^+_{VD}$ at the station of each run. Right:
$-\langle u'v'\rangle/u_\tau^2$ against $y/\delta$; black, the DNS.*

| At $Re_\theta=2835$ | baseline | lengths $1.2, 0.5, 0.5\,\delta_i$ | span $12\,\delta_i$ | WALE | Vreman | DNS |
|---|---|---|---|---|---|---|
| station $x/\delta_i$ | 55.2 | 51.4 | 54.6 | 54.7 | 57.5 | 53 |
| $H$, against the edge state | 4.08 | 4.08 | 4.11 | 4.05 | 4.11 | 4.10 (its profile: 4.22) |
| $C_f$, the wall model's $\tau_w$ | 2.208e-3 (−4.4 %) | 2.248e-3 (−2.7 %) | 2.205e-3 (−4.6 %) | 2.235e-3 (−3.3 %) | 2.206e-3 (−4.5 %) | 2.310e-3 |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 13.1, 17.5, 19.8 | 13.0, 17.7, 19.8 | 13.0, 17.5, 19.9 | 13.2, 17.4, 19.7 | 13.0, 17.5, 19.7 | 13.6, 16.7, 18.8 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.06, 0.81, 1.22 | 1.97, 0.75, 1.15 | 2.08, 0.79, 1.22 | 1.97, 0.76, 1.19 | 1.94, 0.74, 1.27 | 1.49, 0.93, 1.12 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.80, 0.50 | 0.66, 0.57 | 0.74, 0.56 | 0.72, 0.55 | 0.66, 0.49 | 0.62, 0.37 |

**What it shows.**
- **None of the settings changes the layer much.** $H$ stays within 0.06 of the DNS's header value, $C_f$ within 2 % of the baseline's; longer
  integral lengths shorten the inflow's recovery a little (the layer reaches $Re_\theta=2835$ at $x=51$, $C_f$ −2.7 %), as $\delta_i/16$ did.
- **The subgrid model is active, and does not matter at this grid.** WALE's eddy viscosity, computed from a checkpoint, is about 1.5 times the
  molecular one over $0.3$–$0.75\,\delta$ and 0.4 at the edge.
- **What is left is in the stresses**, not in the mean: the outer $-\langle u'v'\rangle$ is 0.49–0.57 at $0.5\,\delta$ against 0.37 in every run,
  and the near-wall $u_{rms}$ peak 1.9–2.1 against 1.5. At $\delta_i/16$ both move toward the DNS ([resolution](#resolution-mach-2-5-at-delta-i-16)).

### Mach 5.86, cold wall

Case M6Tw076 of the same DNS: Mach 5.86, a cold wall ($T_w/T_r=0.76$, $T_w=5.45\,T_\infty$), $Re_\theta=9455$ and $Re_\tau=453$ at its station
($x=54\,\delta_i$), wall heat-flux parameter $B_q=q_w/(\rho_w c_p u_\tau T_w)=-0.02$ (heat into the wall). The wall history gives $q_w$ as well as
$\tau_w$.

| Setup | |
|---|---|
| box, cells | $120\times6\times6\,\delta_i$, $1440\times72\times72$, cubic $\delta_i/12$: the filter's recovery length grows with the Mach number (Huang, Duan & Choudhari 2022) |
| gas | $\mu\propto T^{0.76}$ with $\mu_w$ matched to the DNS's Sutherland value at the wall, where the model acts ($\mu_w/\mu_\infty=5.13$ against $5.45^{0.76}=3.63$): $Re=97120$, the DNS's $137350$ per $\delta_i$ over $1.414$; $Re_\theta$ is reported with $137350$ (`tbl_oracle.py compare --re-inf`) |
| wall, inflow, scheme | as at Mach 2.5 |
| run | `tbl_cases_t37_m6.txt`, to $t=360$ (21604 steps); 29 checkpoints of $180<t\le360$ and the wall history every 20 steps; 0.54 s per step on 4 A30 |

![The Mach 5.86 cold-wall boundary layer against the DNS](/flume/tbl-m6.png)

*The panels of the Mach 2.5 figure; the title adds $B_q$ at the station.*

| At $Re_\theta=9453$ ($x=55\,\delta_i$) | LES | LES, Sutherland | DNS |
|---|---|---|---|
| $C_f$, the wall model's $\tau_w$ | **9.47e-4** (the two halves of the window 9.70e-4, 9.25e-4); $x=60$–$110$: 9.78e-4 | 9.51e-4; $x=60$–$110$: 9.80e-4 | 9.84e-4 |
| $B_q$, the wall model's $q_w$ | **−0.0189** (halves −0.0191, −0.0187); $x=60$–$110$: −0.0194 | −0.0190 | −0.02 (one significant figure) |
| $\delta/\delta_i$ | 1.61 | 1.61 | 1.72 |
| $H$, against the edge state | 13.41 | 13.71 | 13.6 (its profile: 14.9) |
| $u^+_{VD}$ at $y^+=30$, $100$, $200$ | 12.4, 17.9, 20.8 | 12.4, 17.9, 20.6 | 13.9, 17.2, 19.5 |
| $u_{rms}$, $v_{rms}$, $w_{rms}$ over $u_\tau$ at $y/\delta=0.2$ | 2.21, 0.56, 0.94 | 2.24, 0.55, 0.90 | 1.33, 0.82, 0.99 |
| $-\langle u'v'\rangle/u_\tau^2$ at $y/\delta=0.2$, $0.5$ | 0.64, 0.39 | 0.65, 0.35 | 0.49, 0.25 |

**What it shows.**
- **$C_f$ within 3.8 % of the DNS** at its $Re_\theta$, with no change of the wall-model setting from Mach 2.5; past $x=60$ the plateau is
  0.7 % below. **$B_q$ agrees to the precision of the reference**: the DNS file gives $-B_q=0.02$, one significant figure, and the run gives
  0.0189 at the station, 0.0194 on the plateau.
- **$C_f$ drifts by 4.6 % between the two halves of the window; $q_w/\tau_w$ by 0.4 %.** The layer itself changes slowly over the window, and the
  wall model's heat flux follows its friction. 180 time units are 1.5 flow-throughs of this box, so the $C_f$ figure carries about $\pm2$ %.
- **The inflow recovers by $x\approx55$–$60\,\delta_i$**, farther than at Mach 2.5 (50): $C_f$ falls from $1.4\cdot10^{-3}$ at the inlet to
  $0.8\cdot10^{-3}$ at $x\approx33$, then rises to its plateau. The layer reaches the DNS's $Re_\theta$ at $x=55$, against 54 in the DNS.
- **$H$ is the DNS's** against the edge state (13.4; the header 13.6, its profile 14.9): the free stream over this layer is compressed by 4 %
  ([the shape factor](#the-shape-factor-and-the-free-stream)).
- **The stresses have the Mach 2.5 faults, larger.** $\delta$ is 6 % thin; $u_{rms}$ peaks at 1.7 times the DNS's at $y/\delta=0.2$, and
  $u^+_{VD}$ is 1.5 below the DNS at $y^+=30$ (the first cells under the matching point) and 0.7–1.3 above it at $y^+=100$–$200$. At $\delta_i/12$
  the layer of $Re_\tau=453$ has about 20 cells across $\delta$; at Mach 2.5 a finer grid moves the log layer and the stresses toward the DNS
  ([resolution](#resolution-mach-2-5-at-delta-i-16)).
- **Sutherland's law with the DNS's Reynolds number changes little.** The viscosity, matched at the wall, is 1.41 times the DNS's in the free
  stream. A run under Sutherland's law (`viscosity_law = sutherland`, $S=110.4$ K, the DNS's $Re=137350$ per $\delta_i$; `tbl_cases_t37_hdef.txt`,
  29 checkpoints) removes that: $H$ 13.71 against 13.41, $C_f$ −3.4 % against −3.8 % at the station, the plateau 9.80e-4, $B_q$ −0.0190.
- **The momentum integral is useless here too** ($-62$ % at the station, negative in places).

## Shock/boundary-layer interaction (#121 S4b)

The interaction of the IUSTI experiment (Dupont, Haddad & Debiève 2006, JFM 559) as computed by the STREAmS DNS (Bernardini, Modesti, Salvadore &
Pirozzoli 2021, CPC, Fig. 7): Mach 2.28, an oblique shock of $8^\circ$ flow deflection on an adiabatic wall, $Re_\theta\approx5000$. The layer is
the T3.7 setup; the shock comes in through the top face, the [`shock-generator`](../../boundary-conditions#shock-generator-an-oblique-shock-through-a-y-face-121-s3)
of S3. Nothing is fitted to the interaction: the wall model is the one of the zero-pressure-gradient layers.

| Setup | |
|---|---|
| box, cells | $100\times12\times6\,\delta_i$, $1200\times144\times72$, cubic $\delta_i/12$ |
| shock | $\beta=32.685^\circ$, entering the top at $x=41.30$ so that the inviscid shock reaches the wall at $x_{imp}=60$, past the inflow's recovery; exact regular reflection: $p_2/p_1=1.602$, $p_3/p_1=2.458$, reflected at $29.77^\circ$ |
| gas | $Re=U_\infty\delta_i/\nu_\infty=35500$, Sutherland's law with $S/T_\infty=0.776$ (110.4 K over the experiment's 142 K), $Pr=0.71$ |
| wall, inflow, scheme | `wall-modelled` at $T_w=1.918\,T_\infty$ (the recovery temperature), the 3rd cell, semi-local damping, eddy viscosity below $h$; the M2p5 DNS profile with its temperature excess scaled by $(2.28/2.5)^2$ (`reference/tbl_m2p28_profile.txt`); `hybrid-kep-weno` |
| run | `tbl_cases_s4.txt`, to $t=1200$ on 4 A30 (prod-gn, 0.74 s per step); statistics: the span and the 50 checkpoints of $310\le t\le934$ (`sbli_oracle.py compare --t0 300`) |
| references | the DNS wall pressure and the IUSTI points of STREAmS Fig. 7a and Volpiani, Bernardini & Larsson (2018, PRF 3) Fig. 12, digitised from the vector paths (`reference/sbli_m2p28_*`, error in the headers) |

`sbli_oracle.py` takes $\delta_0$, $\delta^*$ and $\theta$ of the undisturbed layer from a linear fit over $[x_{imp}-15,\,x_{imp}-6]$
extrapolated to $x_{imp}$; $C_f$ from the wall model's $\tau_w$; the reflected shock from the mid-pressure crossings of the rows
$2.5<y/\delta_0<5$, a line extrapolated to the wall: its apparent origin $x_{ref}$ gives $L_{int}=x_{imp}-x_{ref}$, the length STREAmS scales
$x$ with.

![The interaction: wall pressure and skin friction](/flume/sbli.png)

*Left: the mean wall pressure $(p_w-p_\infty)/(p_2-p_\infty)$ against $x^*=(x-x_{imp})/L_{int}$, the DNS (line) and the experiment (triangles).
Right: $C_f$; dotted the apparent origin of the reflected shock ($x^*=-1$), dashed the inviscid impingement. `make_doc_figures.py --only sbli`.*

![The interaction's flow field](/flume/sbli-fields.png)

*One checkpoint ($t=934$). From the top: numerical schlieren of the mid-span plane; the spanwise-averaged pressure with the sonic line (white), the
$\bar u=0$ line (red, absent) and the oracle's shock fits (orange, its crossings as dots); $u$ in the wall-parallel plane $y=0.29$
($0.16\,\delta_0$; dashed the reflected shock's apparent origin, dotted $x_{imp}$); the cross-flow plane at $x_{imp}$ and spanwise-averaged
profiles. `make_doc_figures.py --only sbli-fields`.*

| | FLUME, wall model | Reference |
|---|---|---|
| $Re_\theta$ at $x_{imp}$, nominal / against the edge state | 5247 / 4667 | about 5000 (DNS, experiment) |
| $\delta_0$, $\delta^*_0$, $H$ at $x_{imp}-10$ | $1.807$, $0.531\,\delta_i$, 3.58 | |
| reflected shock: angle to the wall | $28.95^\circ$ | $29.77^\circ$ inviscid |
| $L_{int}/\delta_0$ | **3.15** | 3.22, Volpiani et al. (DNS, $Re_\theta=1050$) |
| separation: $L_{sep}/\delta_0$ | **0**, $C_f$ never negative, minimum $0.22\cdot10^{-3}$ at $x^*\approx-0.25$ | 1.62, Volpiani et al., $Re_\theta=1050$; the DNS and the experiment have a bubble |
| $(p_w-p_\infty)/(p_2-p_\infty)$ at $x^*=-0.5$, $0$, $0.5$, $1$, $2$ | 1.12, 1.72, 2.09, 2.24, 2.37 | DNS: 0.00, 0.64, 1.75, 2.08, 2.29 |
| Souverein $L^*=L_{int}/(\delta^*g_3)$, separation parameter $Se$ | 1.93, 1.20 | |

**What it shows.**
- **The interaction length is right**: $L_{int}/\delta_0=3.15$ against 3.22 in the DNS of Volpiani et al. (2018, at $Re_\theta=1050$), and the reflected
  shock leaves at $28.95^\circ$, $0.8^\circ$ under the inviscid angle. The downstream pressure plateau is
  the DNS's within 4 % ($x^*=2$).
- **The mean flow does not separate.** $C_f$ dips to $0.22\cdot10^{-3}$ and recovers; the DNS, the experiment and Volpiani et al. have a bubble. It comes close: in one checkpoint the
  spanwise-averaged $u$ of the first cell falls to $0.003\,U_\infty$ between $x=55$ and 65, and the wall-parallel plane at $0.16\,\delta_0$ has patches of
  $u<0$ (down to $-0.14\,U_\infty$) between $x=57.7$ and 61.9. The wall model's $\tau_w$ follows the velocity at its matching height, the 3rd cell,
  whose mean stays positive.
- **The wall pressure rises too early and too steeply.** It starts at the reflected shock's apparent origin ($x^*=-1$, where the DNS is still flat)
  and passes $p_2$ at $x^*\approx-0.5$, where the DNS has not risen; the DNS starts at $x^*\approx-0.3$. In the DNS a separation bubble spreads the
  compression into a fan; here the reflected shock reaches the wall almost whole: upstream the sonic line is $0.12\,\delta_i$ ($0.07\,\delta_0$, two cells)
  above the wall, and it thickens to $0.6$–$0.8\,\delta_i$ only past the foot of the shock.
- **The likely cause is the wall model, not yet measured.** The equilibrium model has no pressure-gradient term and computes $\tau_w$ from the
  velocity at the 3rd cell, $0.21\,\delta_i$ or $0.12\,\delta_0$ above the wall: the near-wall deceleration under the adverse gradient reaches
  it late. Testing that needs a non-equilibrium model or a finer wall grid at the interaction (a 2:1 box under it).

## Not done

- HLLC at $256^3$, M 0.1 was dropped by owner decision.
- Spectra against the DNS: the file holds integrals only. The DNS of Johnsen et al. for the isotropic case.
- Calibration of the sensor; a model masked by face rather than by cell; `central-bcs` on a forest seam.
- T3: why the classical damping overshoots at Mach 3; a cheaper flux path below $h$; the boundary layer: the matching height, the outer-layer and near-wall stresses at $\delta_i/12$, a longer window at Mach 5.86; $Re_\theta$ and the station against the edge state (they use the nominal free stream); SBLI: the separation (a wall model with the pressure gradient, or a finer wall grid at the interaction), the wall-pressure fluctuations and the low-frequency unsteadiness; runtime regridding along a modelled wall; mean profiles
  against the DNS; the relaminarisation of the inverse model on $96\times48\times48$ with $u_e = 1.18$; an edge state taken from the flow.
- A model on 2:1 grids beyond conservation; runtime AMR with a model; a wall.
