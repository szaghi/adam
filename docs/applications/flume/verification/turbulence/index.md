# Turbulence (T0, T1, T2)

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

## Taylor-Green vortex: the ILES baseline

Re 1600, $[0, 2\pi]^3$, `runge-kutta-ssp-54`, CFL 0.5, `weights = si`, 4 ranks on A30, from `turbulence/tgv_campaign.sh`. `weno` is the flux splitting, `HLLC`
`weno-riemann` with the 6th-order correction and the WENO sensor. All values per unit volume.

![The ILES baseline on the Taylor-Green vortex](/flume/tgv-iles.png)

*Top: M 0.1, WENO-5 by resolution; middle: M 0.1 at $128^3$ by scheme and order; bottom: $M_0 = 1.25$ by resolution. Left $-dK/dt$, right the enstrophy; dotted: the
enstrophy peak about 10.3 quoted for the incompressible reference by Svärd (2022).*

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

**Reference.** The 512³ pseudo-spectral DNS of the High-Order CFD Workshop (case C3.5, `spectral_Re1600_512.gdiag`) is not available to this tree, so no curve is compared,
only the quoted peak. The incompressible reference has its dissipation peak near $t = 9$. Compressible Mach 1.25 is compared with nothing here.

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
recommended. The constants come from other discretizations and were not tuned. Only peaks are compared: a spectral comparison with filtered DNS could rank the models differently.

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

*Top: M 0.1, the hybrid against WENO-5; middle: with Vreman; bottom: $M_0 = 1.25$. Left $-dK/dt$, right the enstrophy; dotted: the enstrophy peak about 10.3 quoted for the incompressible reference.*

| Case, M 0.1 | Peak $-dK/dt$, raw | Over $\Delta t = 0.5$ | at $t$ | Flicker, $7<t<10$ | Peak $Z$ | Numerical share to the peak |
|---|---|---|---|---|---|---|
| WENO-5, $64^3$ | 0.01377 | 0.01357 | 5.09 | 0.3 % | 1.71 | 77 % |
| hybrid, $64^3$ | 0.01591 | 0.01334 | 8.09 | 7.2 % | 3.82 | 50 % |
| hybrid + Vreman, $64^3$ | 0.01255 | 0.01159 | 8.45 | 3.5 % | 2.94 | |
| WENO-5, $128^3$ | 0.01242 | 0.01216 | 8.29 | 0.5 % | 3.40 | 58 % |
| hybrid, $128^3$ | 0.01389 | 0.01317 | 8.48 | 2.0 % | 6.86 | 18 % |
| hybrid + Vreman, $128^3$ | 0.01271 | 0.01230 | 8.56 | 0.8 % | 5.53 | |
| WENO-5, $256^3$ | 0.01272 | 0.01261 | 8.72 | 0.7 % | 6.18 | 27 % |

The flicker is the rms of $-dK/dt$ about its running mean, over the smoothed peak. The share with Vreman is left out: it does not separate the model from the scheme.

| Case, $M_0 = 1.25$ | Peak $-dK/dt$ | Peak $Z$ | Numerical share |
|---|---|---|---|
| WENO-5, $64^3$ / $128^3$ / $256^3$ | 0.01383 / 0.01081 / 0.01101 | 1.81 / 3.34 / 5.95 | 80 % / 46 % / 20 % |
| hybrid, $64^3$ / $128^3$ / $256^3$ | 0.01266 / 0.01057 / 0.01112 | 2.20 / 4.83 / 7.21 | 73 % / 36 % / 10 % |

Isotropic turbulence: peak $Z$ 41.0 and 57.4 at $64^3$, $128^3$ against 34.5 and 51.0 on WENO-5 (60.5 at $256^3$); $K(t/\tau=4)$ 0.367 and 0.408 against 0.303 and 0.384 (0.413 at $256^3$).

**What it shows.**
- At M 0.1 the hybrid reaches at $128^3$ a higher enstrophy peak than WENO-5 at $256^3$, and at $64^3$ more than WENO-5 at $128^3$, with the dissipation peak where the finer runs have it. $Z$ is still below the quoted 10.3.
- At $M_0 = 1.25$ the gain is smaller, as the sensor keeps WENO around the shocklets.
- **The raw peak of $-dK/dt$ is not a measure on the hybrid:** it flickers, 7 % at $64^3$, 2 % at $128^3$. The source is not established; the sensor switching on and off is the first candidate. Smoothed, the hybrid is 4–6 % above WENO-5 at $256^3$.
- **Vreman on the hybrid** removes half of the flicker and lowers $Z$ by 20 %. At $64^3$ it lands below the finer runs (0.0116). Whether a model helps the hybrid is not settled by these peaks: it needs the DNS curves.
- The thresholds (0.2, 0.1) and the floor are not calibrated.

## Not done

- HLLC at $256^3$, M 0.1 was dropped by owner decision.
- Comparison with DNS curves and spectra, which needs the reference files.
- A model on 2:1 grids beyond conservation; runtime AMR with a model; a wall.
