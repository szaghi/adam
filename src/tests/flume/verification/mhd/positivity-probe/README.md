# Positivity probe: the Balsara–Spicer strong blast (issue #41, M2 stretch S0–S2)

A record, not a check: FLUME cannot run the Balsara & Spicer (1999) strong MHD blast. The M2 measurements (2026-09-30,
issue #41 stretch S0–S2) closed the stretch items "parametrised positivity-preserving flux limiter" and "EGLM / ideal
GLM-MHD switch" without adopting either. The M3 measurements (issue #47) found the cause, the mixed-GLM energy
coupling, and **correct the M2 conclusion on EGLM** (see "Conclusion").

The blast: $\rho = 1$, $p = 0.1$ outside and $p = 1000$ for $r < 0.1$, $\mathbf{B} = 100/\sqrt{4\pi}$ at 45°
(ambient $\beta = 2.5 \cdot 10^{-4}$), $\gamma = 1.4$, periodic $[0, 1]^2$, $t = 0.01$.

## Files

| File | What it does |
|------|--------------|
| `make_blast.py` | Writes the blast as a FLUME input (disk as one `riemann-problem` strip region per cell row). |
| `lf_proto.py` | First-order Rusanov, forward Euler, for four models: FLUME's GLM, ideal GLM-MHD (Derigs et al. 2018, eqs. 3.16–3.18), no cleaning, Powell. |
| `split_check.py` | Counts the one-sided Lax–Friedrichs partial states of the Hu–Adams–Shu limiter that are inadmissible along a first-order run. |
| `weno_proto.py` | WENO5 (component-wise, Lax–Friedrichs split) + SSP-RK3 with no limiter, the one-sided limiter (Hu, Adams & Shu 2013) or the cell-based parametrised limiter (Xu 2014); optional Powell source and a scaled Lax–Friedrichs speed. |
| `ws_proto.py` | M3 gate: WENO5 + SSP-RK3 + cell-based limiter over a first-order LF backbone, for no cleaning, FLUME's GLM and EGLM with its sources; optional Wu (2018) wave speed (local or global) and Powell source; also a random-state test of the backbone alone. |
| `one-sided-limiter-cpu.patch` | The FLUME CPU implementation of the one-sided limiter (`[mhd] positivity_limiter`), reverted; `git apply` restores it. |

Run with the regression venv, e.g. `exe/.regression-venv/bin/python weno_proto.py cell --cells 64`.

## Measurements (2026-09-30)

FLUME CPU, 128², WENO-5 characteristic, SSP-54, GLM ($c_h = 60$), floors off unless stated:

| Case | Result |
|------|--------|
| pressure ratio $10^4$, $\beta = 0.2$ | clean to $t = 0.01$, zero floored cells |
| pressure ratio 100, $\beta = 2.5 \cdot 10^{-4}$, floors armed | floors from step 5, then non-finite |
| full blast (ratio $10^4$, $\beta = 2.5 \cdot 10^{-4}$), floors armed, CFL 0.4 or 0.2 | floors from step 3–5, then non-finite |
| full blast, conservative reconstruction, floors armed | reaches $t = 0.01$ only through 1823 floor events |
| $\beta = 2.5 \cdot 10^{-2}$ or $2.5 \cdot 10^{-3}$, ratio $10^4$, CFL 0.25 | non-positive pressure at step 17 / 6 |
| one-sided limiter, CFL 0.25 | fails at step 5 (4 unresolved faces); CFL 0.1: step 10; no cleaning: step 414 |
| one-sided limiter + floors | stays finite through $4.6 \cdot 10^7$ floor events and a collapsing time step |

Prototypes, 64² (WENO5 + SSP-RK3, CFL 0.25), and first order at 128²:

| Variant | Result |
|---------|--------|
| first-order Rusanov, any of the four models, 128² | clean to $t = 0.01$; GLM and EGLM identical |
| one-sided partial states on first-order states, 128² | never inadmissible (615 steps) |
| WENO5, no limiter | fails at step 23 |
| WENO5, one-sided limiter | fails at step 56 (12 unresolved faces per stage) |
| WENO5, cell-based limiter | fails at step 61 (4 cells with an inadmissible first-order update) |
| WENO5, cell-based + Powell source | fails at step 61 |
| WENO5, cell-based, Lax–Friedrichs speed ×2 or ×4 | fails at step 62–69 |

## Measurements, M3 gate (2026-09-30, `ws_proto.py`)

WENO5 + SSP-RK3, CFL 0.4 (on the backbone speeds), cell-based limiter:

| Model | Backbone | Result |
|---|---|---|
| 64², FLUME's GLM | LF, with or without Powell, standard or Wu speed, local or global | fails at step 38–39; 4 backbone-inadmissible cells |
| 64², no cleaning | same variants | reaches $t = 0.01$ (min p 3–5e-5) |
| 64², EGLM without its sources | LF | stops at step 151 with min p 8.6e-14 (admissibility threshold 1e-13); backbone clean |
| 64², EGLM with the Derigs sources | LF (+ Powell) | reaches $t = 0.01$, min p 3.6e-4, backbone never inadmissible |
| 128², EGLM with the Derigs sources | LF (+ Powell) | reaches $t = 0.01$ in 384 steps, min p 1.3e-4 |
| 128², FLUME's GLM | LF (+ Powell) | fails at step 58 |

The random-state test (`ws_proto.py random`, 32², 20 trials, large B jumps, p down to $e^{-25}$) finds no inadmissible
cell for any variant: random states do not discriminate.

## Conclusion

The failure is caused by FLUME's mixed GLM: $\psi$ changes $B_n$ through the flux but is not in the energy, so the
change of magnetic energy is taken from the thermal pressure, and the first-order Lax–Friedrichs update of the
high-order stage states becomes inadmissible at a few cells; no flux limiter, which only blends toward that update,
can then restore positivity. With EGLM (Derigs et al. 2018: $\psi^2/2$ in the energy, the matching energy flux and the
nonconservative sources) the backbone stays admissible and the cell-based limiter reaches the end of the blast; without
cleaning it does too, but that model is not usable in multi-D. The Wu (2018) wave speed and the Powell source alone
change nothing on this problem.

**Correction of the M2 conclusion.** M2 stated that "the GLM energy term of EGLM is not the cause (GLM and EGLM coincide
at first order)". That test ran first-order updates of first-order states, where GLM is clean as well, so it could not
discriminate; the discriminating test is the backbone applied to high-order stage states, above. EGLM and the
cell-based limiter are planned in issue #47 (M3). The positivity is measured, not proven: no published proof covers
cell-centred GLM or EGLM.

References: Balsara & Spicer 1999 (JCP 149, 270); Hu, Adams & Shu 2013 (JCP 242, 169); Xu 2014 (Math. Comp. 83,
2213); Christlieb et al. 2015 (SIAM J. Sci. Comput. 37, A1825); Wu 2018 (SIAM J. Numer. Anal. 56, 2124); Wu & Shu
2019 (Numer. Math. 142, 995); Derigs et al. 2018 (JCP 364, 420, arXiv:1711.06269).
