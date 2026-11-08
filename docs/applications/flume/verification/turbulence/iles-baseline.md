# The ILES baseline (T0)

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| TT, HT | [#119](https://github.com/szaghi/adam/issues/119) · T0.1–T0.6 | P+C | tt1–tt6, ht1–ht2 | CPU, FNL | yes |

## Case and scope

T0 measures what the scheme alone does to turbulence, with no model: the implicit LES (ILES) baseline. The cases are the
Taylor-Green vortex at Mach 0.1 (T0.3) and at $M_0 = 1.25$ (T0.4), the same vortex with a static 2:1 box (T0.6) and
decaying compressible isotropic turbulence (T0.5). The legs check first the diagnostics they rest on (T0.1) and the
initial conditions (T0.2, T0.5). The diagnostics and the keys are on the [turbulence](../../turbulence#turbulence-diagnostics)
page.

## Legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | Phase | What it verifies | Measured |
|---|---|---|---|
| tt1 | T0.1, T0.2 | first history row of `taylor-green` against the exact integrals, M 0.1, 1.25 | K $10^{-16}$; Z and $\varepsilon$ order 3.98–4.00; $\varepsilon_d$, $p_d$ $10^{-30}$ |
| tt2 | T0.1 | one rank against 2 | $6\cdot10^{-14}$ |
| tt3 | T0.3 | energy budget to $t = 2$, Re 1600, M 0.1: what the integrals leave of the $K$ drop | 70 %, 8.8 %, 0.28 % at $32^3$, $64^3$, $128^3$ |
| tt4 | T0.1 | `spectra.py` on the initial field: Parseval, the energy of the single-mode field in shell 2 | not recorded (bound $10^{-12}$) |
| tt5 | T0.1 | (`--agree TAG`) this executable's tt2 history against the executable TAG: backends | not recorded (bound $10^{-11}$) |
| tt6 | T0.1 | pressure work, M 0.8: the budget closes only with $p_d$ | 8.5 % left, 203 % without, 414 % flipped |
| ht1 | T0.5 | `isotropic-turbulence` field: every shell, $K$, dilatation | shells $7\cdot10^{-15}$, $K$ $1.5\cdot10^{-14}$ |
| ht2 | T0.5 | one rank against 2 on the random field | $7\cdot10^{-15}$ |

## Taylor-Green vortex (T0.3, T0.4)

| Setup | |
|---|---|
| box | $[0, 2\pi]^3$, periodic |
| flow | Re 1600; M 0.1 and $M_0 = 1.25$ |
| schemes | `weno` (flux splitting, WENO-5, 7, 9); HLLC: `weno-riemann` with the 6th-order correction and the WENO sensor |
| time | `runge-kutta-ssp-54`, CFL 0.5, `weights = si` |
| run | `turbulence/tgv_campaign.sh`, 4 ranks on A30; all values per unit volume |

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

**Reference.** The 512³ pseudo-spectral DNS of the High-Order CFD Workshop (case C3.5, `turbulence/reference/spectral_Re1600_512.gdiag`), incompressible, so compared at M 0.1 only: peak $-dK/dt$ 0.01286 and peak $Z$ 10.29, both at $t = 8.97$. The L2 errors against it are in [Against the DNS](./hybrid-scheme#against-the-dns). Mach 1.25 has no reference here.

**What it shows.**
- **Numerical dissipation dominates.** Up to the peak of $-dK/dt$ it is 77 % of the energy lost at $64^3$ and 58 % at $128^3$ with WENO-5 at M 0.1.
- **The resolved gradients are far from converged.** At $128^3$ the enstrophy peak is a third of the reference, even where the total $-dK/dt$ looks right. At $64^3$ the peak of $-dK/dt$ comes at $t = 5$, against 9.
- **Resolution matters more than scheme and order at $128^3$.** HLLC gains a few percent, WENO-9 reaches $Z$ 4.45. WENO-7 is below WENO-5 ($Z$ 3.17 against 3.40); the cause is not established.
- **At $M_0 = 1.25$ the numerical share falls with the resolution**: 80 %, 46 %, 20 %.

## A static 2:1 box (T0.6)

| Setup | |
|---|---|
| grid | base $64^3$; the box of side $\pi$ on the centre refined to the $128^3$ spacing (one eighth of the volume) |
| flow, scheme | M 0.1, WENO-5, as above |
| seam operators | tricubic fill with mean restriction; quintic fill with `lagrange6` restriction |

![The models and the 2:1 box](/flume/tgv-models.png)

*Bottom row: the 2:1 box on the centre of the $64^3$ grid against the uniform grids (the upper rows are the [subgrid-scale models](./sgs-models)).*

| Case | Peak $-dK/dt$ | at $t$ | Peak $Z$ | at $t$ | $K(20)$ | $\max\lvert p_d\rvert$ |
|---|---|---|---|---|---|---|
| uniform $64^3$ | 0.01377 | 5.08 | 1.71 | 6.25 | 0.0245 | 0.00026 |
| box, tricubic fill, mean restriction | 0.01312 | 6.39 | 1.87 | 6.39 | 0.0233 | 0.0070 |
| box, quintic fill, `lagrange6` | 0.01314 | 6.39 | 1.85 | 6.38 | 0.0234 | 0.0115 |
| uniform $128^3$ | 0.01242 | 8.24 | 3.40 | 8.11 | 0.0219 | 0.00007 |

**What it shows.**
- **The box delays and lowers the peak of $-dK/dt$**, from $t = 5.1$ to $6.4$.
- **The two seam operators give the same integrals**: $-dK/dt$ and $Z$ do not tell them apart.
- **$p_d$ cannot be read on a 2:1 grid.** It is 27 to 44 times the uniform value. A central difference across a 2:1 face
  does not telescope, so $\int p\,\theta\,dV$ picks up the mean pressure ($p \approx 71$ at M 0.1) times a nonzero discrete
  $\int\theta\,dV$. This explanation is untested. The budget and the numerical share are therefore not reported on a 2:1 grid.

## Decaying compressible isotropic turbulence (T0.5)

| Setup | |
|---|---|
| field | `isotropic-turbulence`, $M_t = 0.6$, $Re_\lambda = 100$, $k_0 = 4$, one seed |
| gas | power-law viscosity 0.76, Pr 0.7 |
| time scale | $\tau = \lambda_0/u_\mathrm{rms} = 0.5$ |
| reference | none: the DNS of Johnsen et al. (2010) is not available to this tree, so only self-convergence is measured |

![Isotropic turbulence](/flume/hit.png)

*Left: the enstrophy against $t/\tau$; right: the spectrum of $\sqrt\rho\,\mathbf{u}$ at $t/\tau = 4$.*

| Case | Peak $Z$ | at $t/\tau$ | $K(t/\tau = 4)$ | $E(16)$ at $t/\tau = 4$ | $E(30)$ |
|---|---|---|---|---|---|
| no model, $64^3$ | 34.5 | 0.66 | 0.303 | $7.4\cdot10^{-4}$ | $8.3\cdot10^{-6}$ |
| Vreman, $64^3$ | 32.6 | 0.58 | 0.292 | $5.9\cdot10^{-4}$ | $5.8\cdot10^{-6}$ |
| no model, $128^3$ | 51.0 | 1.14 | 0.384 | $3.0\cdot10^{-3}$ | $1.0\cdot10^{-4}$ |
| Vreman, $128^3$ | 48.9 | 1.12 | 0.376 | $2.7\cdot10^{-3}$ | $8.1\cdot10^{-5}$ |
| no model, $256^3$ | 60.5 | 1.36 | 0.413 | $4.2\cdot10^{-3}$ | $2.3\cdot10^{-4}$ |

**What it shows.**
- **$256^3$ with WENO-5 is not converged.** $K$ changes by 7 % from $128^3$ to $256^3$.
- **Vreman moves away from the finer runs** at each resolution, as on the Taylor-Green vortex ([subgrid-scale models](./sgs-models)).

The vorticity of both flows on the hybrid scheme is in [The flow](./hybrid-scheme#the-flow).

## Not done

- HLLC at $256^3$, M 0.1.
- Spectra against the DNS: the Workshop file holds integrals only.
- The DNS of Johnsen et al. for the isotropic case.
- The pressure dilatation on a 2:1 grid: the non-telescoping explanation is not tested.
