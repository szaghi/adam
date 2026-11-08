# Euler

The compressible Euler pages of FLUME. They span several milestones: the Euler core of M1
([#35](https://github.com/szaghi/adam/issues/35)), the Riemann-flux scheme of M3
([#47](https://github.com/szaghi/adam/issues/47)), the forest ([#37](https://github.com/szaghi/adam/issues/37) and the
seam issues after it), the quadtrees ([#46](https://github.com/szaghi/adam/issues/46)), and the models of M9 to M13. The
[timeline](../) places each page in implementation order. The two space schemes, `weno` (WENO5 flux splitting on
characteristic variables) and `weno-riemann` (WENO5 face states, an LLF, HLL or HLLC Riemann flux and a 6th-order flux
correction), are described in [numerics](../../numerics).

| Case | ID | Milestone | Type | What it verifies |
|---|---|---|---|---|
| [Sod shock tube](./sod) | V1, RV-3 | M1, M3 | P | $L_1(\rho)$ against the exact solution; $x$, $y$, $z$ bitwise equal; reflecting walls |
| [Isentropic vortex: order of accuracy](./isentropic-vortex) | V2, RV-2 | M1, M3 | P | the design order on a smooth exact solution, periodic blocks across ranks |
| [Shock over a cylinder](./shock-cylinder) | V6 | M1, M3 | P | the immersed boundary with solid AMR: refinement, mirror symmetry, positivity |
| [Conservation across AMR seams](./amr-conservation) | V3, RV-4 | M1, M3 | C | the reflux at 2:1 faces: the volume integrals constant to round-off |
| [Restart round trip and output channels](./restart) | V7 | M1 | C | a restarted run equals the uninterrupted one bit for bit on four integration paths, an empty rank included; slices and auxiliary fields |
| [Forest seams: the inter-realm contract](./multi-realm) | multirealm legs 1–6, 15; MV-14 | #37 to #71 | C | mirror, misaligned and 2:1 inter-realm seams, across ranks: bitwise equal to one realm |
| [Lax problem](./lax) | RV-3 | M3 | P | the Riemann solvers on a strong contact, bound $1.05\times$ the splitting scheme |
| [Shu–Osher shock / entropy-wave interaction](./shu-osher) | RV-3 | M3 | P | the entropy waves behind a Mach 3 shock against a 16× finer reference |
| [Forward-facing step](./forward-facing-step) | EV-step | #46 | P+C | a three-realm forest on quadtrees against octrees; the Woodward–Colella step at $t = 4$ |
| [Dual energy: the advected entropy of euler-de](./dual-energy) | DE | M9 | P+C | `euler-de` changes nothing where no cell is cold (Sod bitwise) and keeps a cold Noh shock heated |
| [Passive scalars](./passive-scalars) | PS | M10 | P+C | up to four scalars that do not act on the flow, are conserved, and stay uniform through runtime AMR |
| [Radiative cooling](./radiative-cooling) | CO | M11 | P | optically thin cooling by Townsend's exact integration: the isochoric ODE to round-off, the floor, a radiative shock against its steady structure |
| [Isothermal equation of state](./isothermal-eos) | EO (eo) | M13 | P | $P = c_0^2\rho$: the shock tube and a Mach-5 collision against the exact solution, the fifth order on an acoustic wave, the inert energy slot, the contract |
| [Barotropic equation of state](./barotropic-eos) | EO (eb) | M13 | P | a piecewise polytrope: a fan split at a convex break, shocks across segments, the order with $c^2 = P'(\rho_0)$, the Mathew & Federrath table by self-convergence |
| [Barotropic gases under gravity](./barotropic-gravity) | EO (eg) | M13 | P | gravity on the momentum only, Jeans waves with $c^2 = P'(\rho_0)$, the Jeans marker, sinks with $c_s^2 = dP/d\rho$ and the virial thermal energy, Bonnor–Ebert stability, Boss–Bodenheimer |

The figures are drawn by `src/tests/flume/verification/make_doc_figures.py` from the work directories that the
`check.sh` scripts leave behind (the step's with `--keep`):

```bash
cd src/tests/flume/verification
(cd sod && ./check.sh) ; (cd vortex && ./check.sh) ; (cd conservation && ./check.sh)
(cd riemann-flux && ./check.sh                                    # HLLC, every leg
   ./check.sh --solver split --legs lax ; ./check.sh --solver split --legs shu-osher
   ./check.sh --solver llf --legs lax --lax-l1-max 1 ; ./check.sh --solver llf --legs vortex
   ./check.sh --solver hll --legs vortex)                           # LLF exceeds the default Lax bound
(cd step && ./check.sh --legs trees,full --keep)
python3 make_doc_figures.py --out ../../../../docs/public/flume \
   --only sod lax shu-osher vortex shock-cylinder conservation step step-trees
```
