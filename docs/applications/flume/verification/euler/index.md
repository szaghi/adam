# Euler

The compressible Euler cases of FLUME: Riemann problems with exact solutions, a smooth exact solution for the order of
accuracy, the immersed boundary, conservation across AMR seams, and the multi-realm forest. They verify the two space
schemes, `weno` (WENO5 flux splitting on characteristic variables) and `weno-riemann` (WENO5 interpolation of the face
states, an LLF, HLL or HLLC Riemann flux and a 6th-order flux correction), described in [numerics](../../numerics).

Every case runs from its directory under `src/tests/flume/verification/` through `check.sh` (CPU by default,
`FLUME_EXE=exe/adam_flume_fnl` for the FNL backend), with 2 MPI ranks and WENO5 unless the page states otherwise. The
oracles assert the physics (an exact solution, a reference, a symmetry or a conservation property), not a stored
golden; the goldens of the regression suite are listed in each page's header.

| Case | ID | What it verifies |
|---|---|---|
| [Sod shock tube](./sod) | V1, RV-3 | $L_1(\rho)$ against the exact solution; $x$, $y$, $z$ bitwise equal; reflecting walls |
| [Lax problem](./lax) | RV-3 | the Riemann solvers on a strong contact, bound $1.05\times$ the splitting scheme |
| [Shu–Osher shock / entropy-wave interaction](./shu-osher) | RV-3 | the entropy waves behind a Mach 3 shock against a 16× finer reference |
| [Isentropic vortex: order of accuracy](./isentropic-vortex) | V2, RV-2 | the design order on a smooth exact solution, periodic blocks across ranks |
| [Shock over a cylinder](./shock-cylinder) | V6 | the immersed boundary with solid AMR: refinement, mirror symmetry, positivity |
| [Conservation across AMR seams](./amr-conservation) | V3, RV-4 | the reflux at 2:1 faces: the volume integrals constant to round-off |
| [Multi-realm runs](./multi-realm) | multirealm legs 1–6 | mirror, misaligned and 2:1 inter-realm seams, across ranks: bitwise equal to one realm |
| [Forward-facing step](./forward-facing-step) | EV-step | a three-realm forest on quadtrees against octrees; the Woodward–Colella step at $t = 4$ |
| [Dual energy: the advected entropy of euler-de](./dual-energy) | DE | `euler-de` changes nothing where no cell is cold (Sod bitwise) and keeps a cold Noh shock heated |

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
