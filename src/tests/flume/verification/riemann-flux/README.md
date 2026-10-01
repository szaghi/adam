# Riemann face fluxes (issue #47, M3)

The M3 scheme (`scheme_space = weno-riemann`; implemented for Euler with the LLF, HLL and HLLC solvers, MHD to come): WENO5 **interpolation** of point values to the face, a Riemann
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
| `weno_interpolation_tables.py` | Generates (sympy, exact rationals) and checks the WENO interpolation tables of `adam_weno_object%initialize_interpolation`, S = 1 .. 5: candidates exact to degree S−1, linear-weight combination exact to degree 2S−2, weights positive and summing to 1. The unit tests `src/tests/flume/unit/test_flume_weno_interpolation{,_fnl}.F90` check the Fortran tables and the device primitive. |
| `check.sh` | FLUME verification of the implemented scheme on the V1/V2/V3/V6 inputs with the `[numerics]` block rewritten: RV-2 (the V2 isentropic vortex ladder 64/128/256, L1 order of the finest pair ≥ 4.5), RV-3 (Sod, Lax and Shu–Osher along x, y, z: L1(ρ) ≤ bound, y and z bitwise equal to x), RV-4 (the V3 AMR box: integrals constant with reflux, drifting without) and V6 (the shock-cylinder with immersed boundary and solid AMR); options `--leg --solver --correction --sensor --recon --order-min --l1-max --lax-l1-max --so-l1-max`, `--solver split` runs the flux-splitting baseline, `FLUME_EXE` selects the backend. |
| `shu_osher_oracle.py` | RV-3 Shu–Osher oracle: L1(ρ) against the split scheme of `hybrid_proto.py` on 16× finer cells averaged to the grid (computed once, cached in the temporary directory, ~3 min), positivity, x/y/z bitwise. |
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
  The default became characteristic after the M3-P3c measurements (below; owner decision 2026-10-01).
- **The high-order update needs high-order EGLM sources**: 2nd-order sources drop the Alfvén wave to order 2.26.
- The limiter is inactive on the smooth wave (identical errors).
- The correction is not a source of oscillation at the Brio–Wu shocks (overshoot ≤ 5e-5, halved by the sensor).

## Measurements, FLUME (2026-09-30, `check.sh`, np 2)

Euler, LLF, 6th-order correction, characteristic interpolation, `weno-u-5`. CPU and FNL agree to all printed digits.

| Variant | Vortex L1(ρ), N = 64 / 128 / 256 | Order 128 → 256 | Sod L1(ρ) |
|---|---|---|---|
| split (FLUME `weno`, baseline) | 8.218e-5 / 4.621e-6 / 8.379e-8 | 5.79 | 3.244e-3 |
| weno-riemann, LLF, WENO-weight sensor | 1.182e-4 / 7.073e-6 / 1.415e-7 | 5.64 | 3.426e-3 (x, y, z bitwise) |
| weno-riemann, LLF, sensor `none` | 9.113e-5 / 4.748e-6 / 1.415e-7 | 5.07 | — |

- The scheme is 5th order in FLUME; at N = 256 the sensor is inactive (same digits with and without it), at N = 128
  it trips at the vortex extrema (the library weights use the exponent S = 3, the prototype 2) and costs a factor 1.5.
- LLF is 1.7× more dissipative than the per-wave splitting on the vortex and 6% on Sod, as in the prototype
  (LLF row of the 1-D table); HLL/HLLC (M3-P2) are expected to close the gap. The RV-3 bound for LLF is 3.49e-3
  (measured + 2%).

## Measurements, FLUME, HLL and HLLC (2026-09-30, `check.sh`, np 2, M3-P2)

Euler, 6th-order correction, WENO-weight sensor, characteristic interpolation, `weno-u-5`; x, y and z bitwise equal in
every RV-3 run. Lax on the Sod inputs: left `0.445, 0.698, 3.528`, right `0.5, 0, 0.571`, t = 0.13; Shu–Osher: 400 cells
on [−5, 5], t = 1.8, against the 16× split reference.

| L1(ρ) | split (`weno`) | LLF | HLL | HLLC |
|---|---|---|---|---|
| Sod | 3.244e-3 | 3.426e-3 | 3.013e-3 | **2.942e-3** |
| Lax | **1.064e-2** | 1.214e-2 | 1.125e-2 | 1.081e-2 |
| Shu–Osher | 0.2605 (max ρ 4.563) | 0.3009 (4.550) | 0.2490 (4.590) | **0.2425** (4.592; ref 4.672) |
| Vortex, N = 256 (order 128 → 256) | 8.379e-8 (5.79) | 1.415e-7 (5.64) | 8.977e-8 (6.26) | **8.041e-8** (6.49) |

Lax with HLLC, by correction: 6th + sensor 1.081e-2, 6th without sensor 1.146e-2, 4th 1.077e-2, none 1.084e-2.

HLLC: RV-4 drift 0 in all five integrals with reflux, 1–2e-5 without; V6 mirror asymmetry 4.9e-12 (split 1.2e-11),
positive ρ and p.

- **HLLC is the default solver of `check.sh`** and the best of the three on Sod (−9% vs split), Shu–Osher (−7%, the
  entropy waves better resolved) and the vortex; on Lax it is 1.6% worse than split. The correction is not the cause
  (1.8% worse without it): the gap belongs to interpolating the face states and solving the Riemann problem, on a
  problem measured almost entirely at its discontinuities. The sensor is needed (6% worse without it).
- **RV-3 criterion (owner decision, 2026-09-30): L1(ρ) ≤ 1.05 × the split scheme's**, per problem, for the default
  solver: Sod 3.41e-3, Lax 1.12e-2, Shu–Osher 0.274 (the `check.sh` defaults). LLF exceeds them on Lax and Shu–Osher
  and is run with explicit bounds.
- The first Lax measurement was wrong: `sod/sod_oracle.py` took both velocities as zero (Sod only). It now reads the
  normal velocity of the input; the V1 Sod results are unchanged.

## MHD Riemann solvers (2026-09-30, M3-P3a)

`adam_flume_mhd_riemann_library`: LLF, HLL and HLLD (Miyoshi & Kusano 2005) without divergence control and with GLM
(the `(B_n, psi)` subsystem solved exactly at the face), working in the frame of the direction so that rotated problems
run the same arithmetic; face states from primitive or characteristic fields (the MHD default since M3-P3c). Unit tests
`src/tests/flume/unit/test_flume_mhd_riemann{,_fnl}.F90` (RV-0, 5000 random states × 3 directions, γ = 5/3):

| Check | Result |
|---|---|
| RS(q, q) = f(q), six solvers | 1.6e-14 |
| HLLD exact on an isolated contact / tangential / rotational discontinuity | 1.4e-14 / 9.6e-14 / 2.6e-14 (construction Rankine–Hugoniot residual 3.5e-16) |
| cyclic invariance, six solvers | bitwise |
| first-order 1-D update positive, LLF, HLL, HLLD | 0 inadmissible updates |
| device vs host | 5.8e-14, fallback flags identical |

Errors relative to `max(1, |f|, s |q|)`, the round-off scale of the star fluxes. A mutation (the sign of the
tangential-field term of the double-star velocity) fails the rotational check at 0.73.

**HLLD fallback.** The first version tested the star pressure as `p*_T − |B*|²/2`, which differs from the pressure of
the conservative star state in the approximate solver: it fell back on 18% of mild random pairs (ratios within
10^0.25). With the pressure of the conservative star state (from its energy) the rate is 0.015%, 0.4% and 1.8% for
density, pressure and field ratios within 10^0.25, 10 and 10², mostly from `S*_L ≤ S_L` (the speed estimate no longer
bounds the Alfvén wave). The degeneracy threshold 1e-12 against 1e-8 changes nothing measurable.

## MHD on `weno-riemann` (2026-10-01, M3-P3b, M3-P3c)

The MHD verification scripts run on `weno-riemann` through `--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]`
(`mhd/numerics.sh`; recon `characteristic`, correction `6th`, sensor `weno` by default), e.g.
`mhd/rj2a/check.sh --numerics hlld` or `--numerics hlld:primitive`. Without the option they run the flux-splitting
scheme as before. `hlld:primitive` and `hlld:characteristic` (6th, weno) assert bounds measured on them (CPU, plus 2 %;
energy minus 0.2 %), any other spec is a measurement (scheme-independent checks only). np 2, WENO5, SSP; CPU and FNL
agree to every printed digit except the round-off symmetry defects.

| Check (RV) | split | HLLD primitive | HLLD characteristic |
|---|---|---|---|
| RJ2a L1 sum, N = 256 / 512 (RV-6) | 3.842e-2 / 2.117e-2 | 5.133e-2 / 2.926e-2 | 3.764e-2 / 2.123e-2 |
| RJ2a HLL L1 sum, N = 256 (RV-6: HLLD ≤ HLL) | — | 5.983e-2 | 4.316e-2 |
| RJ2a x/y/z (RV-6) | bitwise | bitwise | bitwise |
| RJ2a GLM vs none, CPU / FNL | bitwise | 1.3e-12 / 3.3e-12 | 9.8e-13 |
| Brio–Wu (none, GLM), RJ4d: floors, HLLD fallbacks (RV-6) | 0 | 0, 0 | 0, 0 |
| linear waves, 8 families, finest-pair order (RV-5) | 4.96–5.02 | 4.96–4.98 | 4.96–4.98 |
| CPAW eps N = 128, order; left = right (RV-5) | 1.09e-6, bitwise | 1.153e-6, 5.00, bitwise | 1.151e-6, 4.99, bitwise |
| Orszag–Tang symmetry, conservation (RV-7) | — | 1.4e-11, 1.4e-17 | 9.8e-11, 1.4e-17 |
| Orszag–Tang AMR symmetry, seam div(B) ratio (RV-7) | 9.1e-13, 0.76 | 3.3e-14, 0.585 | 1.4e-13, 0.597 |
| field loop ⟨\|B_z\|⟩/A0, N = 64 / 128 (RV-7) | 1.02e-3 / 5.89e-4 | 9.76e-4 / 5.44e-4 | 1.06e-3 / 5.81e-4 |
| field loop E_B(T)/E_B(0), N = 64 / 128 (RV-7) | 0.872 / 0.937 | 0.869 / 0.936 | 0.887 / 0.944 |
| rotor symmetry, step 100 CPU / FNL (RV-7) | 3.3e-11 | 3.9e-12 / 2.1e-12 | 4.9e-10 / 7.9e-10 |
| rotor symmetry, final CPU / FNL (RV-7) | 1.8e-9 | 5.4e-9 / 2.3e-7 | 1.1e-7 / 6.1e-6 |

- **Interpolated variables (D-5).** On smooth waves the two interpolations agree to 3–4 digits: the WENO weights are the
  linear ones, and a linear interpolation commutes with the frozen eigenvector projection. They differ at
  discontinuities: characteristic is 27 % more accurate on RJ2a and dissipates less magnetic energy in the field loop;
  primitive generates slightly less ⟨|B_z|⟩ and keeps the 7×7 eigenvector projection out of the face kernel, which
  buys no occupancy (below). **Owner decision 2026-10-01: characteristic is the MHD default**; its positivity on the
  Balsara–Spicer blast, where the prototype chose primitive, is settled with EGLM and the limiter (P4, P5). In the
  1-D linear waves the weno-riemann and splitting errors coincide to 4 digits for three families (the time error
  dominates there); the inputs and logs confirm `weno-riemann`.
- **GLM vs no cleaning is no longer bitwise.** In 1-D B_n is uniform and psi starts at zero; the splitting's
  characteristic (B_n, psi) block is inert, but `weno-riemann` interpolates B_n, whose face values differ by round-off,
  which seeds psi (9e-16) and through the shocks the other variables (1e-12). `rj2a/check.sh` passes `--pair-tol 1e-11`
  under `--numerics`.
- **Rotor symmetry.** The final defect grows smoothly from round-off (time series every 50 steps: no jump, so no
  discrete switch flips) and orders by the solver's dissipation, final 3.4e-11 HLL, 1.8e-9 split, 5.4e-9 HLLD, 7.4e-8
  HLLD without the sensor, 1.1e-7 HLLD characteristic; the FNL FMA contraction seeds it 40–55 × larger. The sensor
  does not amplify it (switching it off makes it larger). A symmetry defect of the scheme is at truncation level
  (1e-3) from the first steps, so under `--numerics` the check asserts step 100 within 1e-8 and the end within 1e-4.

**Nsight Compute record (R-2; RTX 4070, cc 8.9, first stage, block 128).** The face kernels are register-limited to two
blocks per SM, MHD and Euler alike: MHD-GLM HLLD (Orszag–Tang 128²) 255 registers, theoretical occupancy 16.7 %, achieved
16.2 %; Euler HLLC (vortex 64²) 174 registers, theoretical 16.7 % (achieved 8.3 %, the grid fills 0.37 waves). The
per-thread local arrays (stencil, eigenvectors, interpolated fields) live in local memory: 528 MB of local load and
store against 29 MB of global load for MHD, 75 MB against 5.2 MB for Euler, about 32 and 18 KB per face, the ratio of the
variables (9 against 5). The 255 registers cost no occupancy relative to Euler, so the kernel split is not needed; the
local-array traffic of the shared face loop (both models) belongs to the performance work. Each face launch is followed
by a one-block kernel (16 registers), the finalisation of the `reduction(+:fallbacks)`. No WSL timing is quoted.
