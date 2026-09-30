# Riemann face fluxes (issue #47, M3)

The M3 scheme (`scheme_space = weno-riemann`, planned): WENO5 **interpolation** of point values to the face, a Riemann
flux F of the two face states, and a high-order correction built from F and the cell fluxes f (Chen, Tóth & Gombosi
2016), so that the flux difference is 5th-order accurate in point-value finite-difference form:

$$\hat F_{i+1/2} = \tfrac{64}{45} F - \tfrac{13}{60}(f_i + f_{i+1}) + \tfrac{1}{180}(f_{i-1} + f_{i+2}) \quad\text{(6th)},\qquad
\hat F_{i+1/2} = \tfrac{4}{3} F - \tfrac{1}{6}(f_i + f_{i+1}) \quad\text{(4th)}.$$

The coefficients match $\hat F = f - \tfrac{h^2}{24} f'' + \tfrac{7 h^4}{5760} f^{(4)}$ at the face (derived with sympy;
the weights sum to 1, leading error $\tfrac{23}{138240} h^6 f^{(6)}$).

## Files

| File | What it does |
|------|--------------|
| `hybrid_proto.py` | NumPy prototype, 1-D Euler: the hybrid scheme (LLF, HLL with Einfeldt speeds, HLLC with Batten speeds; correction 6th, 4th or none; characteristic or primitive interpolation; sensors), and FLUME's current flux-splitting scheme as the baseline (`--scheme split`). |
| `hybrid2d_proto.py` | NumPy prototype, 2-D MHD (periodic): the hybrid scheme with primitive interpolation, HLL or LLF, the (B_n, ψ) subsystem solved exactly at the face, GLM or EGLM (Derigs et al. 2018) with sources of order 2 or 4, and the cell-based positivity limiter; tests `blast`, `cpaw` (oblique Alfvén wave), `brio-wu` (periodic double problem). Imports `lf_proto.py` from `../mhd/positivity-probe`. |

Run with the regression venv, e.g. `exe/.regression-venv/bin/python hybrid_proto.py shu-osher --cells 512 --sensor weno`.

## Measurements (2026-09-30)

Smooth density wave ($\rho = 1 + 0.2 \sin 2\pi x$, $u = p = 1$, periodic, $t = 1$), SSP-RK3 with
$dt \propto h^{5/3}$ (at a fixed CFL the RK3 error, $O(h^3)$, hides the space order from N = 160 on):

| Scheme | L1(ρ), N = 320 | Order 160 → 320 |
|---|---|---|
| split (FLUME today) | 5.71e-10 | 5.06 |
| hybrid, HLLC, 6th | 2.88e-10 | 5.06 |
| hybrid, LLF, 6th | 6.23e-10 | 5.06 |
| hybrid, HLLC, 4th | 3.60e-10 | 4.77 (→ 4) |
| hybrid, HLLC, none | 1.29e-5 | 2.00 |

Shocks, L1(ρ) against an 8× finer split-scheme reference averaged to the grid (N = 256; Shu–Osher N = 512 on
[0, 10]); characteristic interpolation unless stated:

| Scheme | Sod | Lax | Shu–Osher |
|---|---|---|---|
| split (FLUME today) | 2.14e-3 | 6.78e-3 | 0.157 |
| HLLC, 6th, correction always on | 1.82e-3 | 5.89e-3 | 0.114 |
| **HLLC, 6th, WENO-weight sensor** | **1.66e-3** | **5.71e-3** | **0.111** |
| HLLC, 6th, jump sensor | 1.67e-3 | 5.86e-3 | 0.167 |
| HLLC, 4th, always on | 1.77e-3 | 5.78e-3 | 0.109 |
| HLLC, no correction | 1.68e-3 | 5.88e-3 | 0.228 |
| HLL, 6th, WENO-weight sensor | 1.70e-3 | 6.04e-3 | 0.113 |
| LLF, 6th, WENO-weight sensor | 1.99e-3 | 6.67e-3 | 0.146 |
| HLLC, 6th, WENO-weight sensor, primitive | 1.47e-3 | 6.60e-3 (max ρ 1.3087 vs 1.3041) | 0.120 |

No run produced a non-positive state; the maxima of ρ overshoot the reference by at most 3e-4, except the primitive
interpolation on Lax.

## Findings

- The hybrid scheme is 5th order with the 6th-order correction, at half the error of the split scheme (HLLC); the
  4th-order correction degrades toward 4.
- **Sensor: the WENO-weight sensor** (the correction is switched off at a face where, in any characteristic field and
  on either side, $\min_k \omega_k / d_k < 0.2$) is the best on all three problems and inactive on the smooth wave;
  the jump sensor switches the correction off in the smooth entropy waves of Shu–Osher (worse than split).
- The correction is not a source of oscillation on these 1-D Euler problems (always-on is within 10% of the sensor);
  risk R-1 of #47 is to be tested on MHD and 2-D strong shocks.
- Characteristic interpolation is preferred: primitive interpolation overshoots on Lax.
- The reference is the split scheme, which biases the comparison slightly in its favour.

## Measurements, 2-D MHD (2026-09-30, `hybrid2d_proto.py`)

Hybrid scheme: primitive WENO5 interpolation, HLL, 6th-order correction with the WENO-weight sensor, SSP-RK3, CFL 0.4,
unless stated.

**Source order (EGLM).** Circularly polarised Alfvén wave along (1, 2) on [0, 1]² (Tóth 2000: ρ = 1, p = 0.1,
B_par = 1, amplitude 0.1), one period, c_h = 2, dt ∝ h^{5/3}; L1 of B:

| Variant | N = 32 | N = 64 | N = 128 | Order 64 → 128 |
|---|---|---|---|---|
| GLM | 3.455e-4 | 1.140e-5 | 3.524e-7 | 5.02 |
| EGLM, 2nd-order sources | 3.568e-4 | 2.766e-5 | 5.786e-6 | **2.26** |
| EGLM, 4th-order sources | 3.455e-4 | 1.141e-5 | 3.527e-7 | 5.01 |
| EGLM, 4th-order sources, cell limiter on | 3.455e-4 | 1.141e-5 | 3.527e-7 | 5.01 |

Along a diagonal the test does not discriminate: B depends on x + y only and the discrete div B cancels exactly, so the
source vanishes whatever its order.

**Positivity.** Balsara–Spicer blast (β = 2.5e-4, c_h = 60) and its |B| = 1000/√(4π) variant (β = 2.5e-6, c_h = 400),
t = 0.01:

| Variant | Result |
|---|---|
| 64², EGLM, cell limiter, sources of order 4 or 2 | reaches t = 0.01, **0 limited faces**, min p 0.09998 (order 4) and 0.1000 (order 2) |
| 64², EGLM, LLF, cell limiter | reaches t = 0.01, 0 limited faces |
| 64², EGLM, correction always on (no sensor), cell limiter | reaches t = 0.01, 0 limited faces |
| 64², EGLM, **no limiter** (with or without the sensor) | reaches t = 0.01, min p 0.09996–0.09998 |
| 128², EGLM, no limiter / cell limiter | reaches t = 0.01 in 387 steps (identical) |
| 64², β = 2.5e-6, EGLM, no limiter / cell limiter | reaches t = 0.01 in 1450 steps, min p 0.101 (identical) |
| 64², GLM, cell limiter | fails at step 12 (4 backbone-inadmissible cells) |

For comparison, the split scheme with EGLM (`../mhd/positivity-probe/ws_proto.py`) needs the limiter: without it the
blast fails at step 14.

**Oscillations at an MHD shock (risk R-1 of #47).** Periodic double Brio–Wu (γ = 2, [0, 2), N = 400, t = 0.1, no
limiter), against LLF without correction on 3200 cells averaged to the grid:

| Variant | L1(ρ) | L1(B_y) | ρ overshoot | B_y overshoot |
|---|---|---|---|---|
| HLL, 6th, correction always on | 7.17e-3 | 1.10e-2 | +4.0e-5 | +4.8e-5 |
| HLL, 6th, WENO-weight sensor | 6.80e-3 | 1.05e-2 | +2.1e-5 | +2.5e-5 |
| HLL, no correction | 6.77e-3 | 1.05e-2 | +4e-10 | +5e-10 |
| LLF, 6th, WENO-weight sensor | 7.05e-3 | 1.08e-2 | +2.1e-5 | +2.5e-5 |

GLM and EGLM give the same digits here (1-D: ψ = 0 and ∇·B = 0).

## Findings, 2-D MHD

- **EGLM is required** for positivity with either flux: GLM fails the blast with the limiter on.
- **The hybrid flux with primitive interpolation and EGLM passes the blast without the limiter**, at β = 2.5e-4 and
  2.5e-6; the limiter never engages. The split scheme with EGLM needs it. This measurement set the MHD default of
  `weno-riemann` to primitive interpolation (#47, D-5 amended); characteristic interpolation was not prototyped for MHD.
- **The high-order update needs high-order EGLM sources**: 2nd-order sources drop the Alfvén wave to order 2.26.
- The limiter is inactive on the smooth wave (identical errors).
- The correction is not a source of oscillation at the Brio–Wu shocks (overshoot ≤ 5e-5, halved by the sensor).
