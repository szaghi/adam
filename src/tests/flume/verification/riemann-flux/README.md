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
