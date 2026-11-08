# Subgrid-scale models (T1)

| ID | Milestone | Type | Legs | Backends | Figure |
|---|---|---|---|---|---|
| SG | [#119](https://github.com/szaghi/adam/issues/119) · T1 | P+C | sg1–sg4 | CPU, FNL | yes |

## Case and scope

T1 adds the algebraic eddy-viscosity models Vreman, WALE and Sigma (`[turbulence] sgs_model`,
[turbulence](../../turbulence#subgrid-scale-models)) and asks the go/no-go question of #119: does a model improve on the
[ILES baseline](./iles-baseline)? The legs check the formulas, the energy budget with the model and rank invariance; the
campaign reruns the Taylor-Green vortex of T0.3 with each model.

## Legs

`src/tests/flume/verification/turbulence/check.sh`, CPU and FNL.

| Leg | What it verifies | Measured |
|---|---|---|
| sg1 | the formulas against numpy on 200 tensors (`test_sgs`); the zeros of each model | $3.6\cdot10^{-15}$, $1.2\cdot10^{-15}$, $1.8\cdot10^{-11}$ |
| sg2 | budget with the model, M 0.4: Vreman 36 %, 5.7 %, 0.27 % left at $32^3$, $64^3$, $128^3$ | at $128^3$: 21.6 % without $\varepsilon_\mathrm{sgs}$, 43 % flipped |
| sg3 | no viscosity with each model | $\varepsilon = 0$, mass and energy $3.3\cdot10^{-16}$ |
| sg4 | one rank against 2 with a model | $10^{-14}$ |

## Taylor-Green vortex with a model

| Setup | |
|---|---|
| flow, scheme, run | as the [ILES baseline](./iles-baseline#taylor-green-vortex-t0-3-t0-4): Re 1600, M 0.1, WENO-5, `runge-kutta-ssp-54`, CFL 0.5, 4 ranks on A30 |
| models | Vreman, WALE, Sigma, default constants (not tuned) |
| grids | $64^3$, $128^3$; Vreman also at $M_0 = 1.25$, $128^3$ |

![The models and the 2:1 box](/flume/tgv-models.png)

*Top: M 0.1 at $64^3$; middle: at $128^3$; bottom: a [2:1 box](./iles-baseline#a-static-2-1-box-t0-6) on the centre of the $64^3$ grid against the uniform grids.*

| Case, M 0.1 | Peak $-dK/dt$ | at $t$ | Peak $Z$ | Share to the peak | $K(20)$ |
|---|---|---|---|---|---|
| no model, $64^3$ | 0.01377 | 5.08 | 1.71 | 77 % | 0.0245 |
| Vreman, $64^3$ | 0.01385 | 5.05 | 1.59 | 79 % | 0.0236 |
| WALE, $64^3$ | 0.01374 | 5.04 | 1.60 | 79 % | 0.0239 |
| Sigma, $64^3$ | 0.01380 | 5.17 | 1.59 | 79 % | 0.0237 |
| no model, $128^3$ | 0.01242 | 8.24 | 3.40 | 58 % | 0.0219 |
| Vreman, $128^3$ | 0.01210 | 8.29 | 3.18 | 61 % | 0.0213 |
| WALE, $128^3$ | 0.01199 | 8.29 | 3.16 | 61 % | 0.0210 |
| Sigma, $128^3$ | 0.01206 | 8.30 | 3.18 | 60 % | 0.0213 |

The share leaves $\varepsilon_\mathrm{sgs}$ in: it is the energy lost not accounted for by $\varepsilon - p_d$, so it
includes the model's own dissipation. At $M_0$ 1.25, $128^3$, Vreman gives peak $-dK/dt$ 0.01085 and peak $Z$ 3.10,
against 0.01081 and 3.34 without a model.

Against the 512³ DNS (L2 errors over $0.3 < t < 19.7$; the full table is in [Against the DNS](./hybrid-scheme#against-the-dns)):

| Case | L2 $-dK/dt$ | L2 $Z$ | L2 $K$ |
|---|---|---|---|
| WENO-5, $64^3$ | 0.660 | 0.844 | 0.187 |
| Vreman, $64^3$ | 0.671 | 0.854 | 0.200 |
| WENO-5, $128^3$ | 0.278 | 0.643 | 0.083 |
| Vreman, $128^3$ | 0.295 | 0.662 | 0.092 |
| WALE, $128^3$ | 0.296 | 0.660 | 0.091 |
| Sigma, $128^3$ | 0.293 | 0.660 | 0.091 |

**What it shows.**
- **No model beats no model.** On a dissipative scheme each one lowers the enstrophy peak further from the reference
  (1.59 against 1.71 at $64^3$, 3.18 against 3.40 at $128^3$), and all three worsen every L2 error.
- **The three models are nearly the same.**
- **The eddy viscosity adds to a numerical dissipation that already dominates.** In the laminar phase, to $t = 2$ at
  M 0.4, it takes 21–53 % of the energy lost.
- **Go/no-go: no-go for a model on WENO.** Under the rule of #119 this makes the hybrid central/WENO scheme
  ([T2](./hybrid-scheme)) necessary before any model is recommended.

The constants come from other discretizations. The same models in decaying isotropic turbulence are in
[the ILES baseline](./iles-baseline#decaying-compressible-isotropic-turbulence-t0-5); on the hybrid scheme in
[T2](./hybrid-scheme).

## Not done

- A spectral comparison with filtered DNS, which could rank the models differently.
- Tuned constants.
- A model on 2:1 grids beyond conservation (hy9); a model with runtime AMR.
