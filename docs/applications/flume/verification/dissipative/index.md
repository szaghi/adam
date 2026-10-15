# Dissipative terms (M4)

Milestone M4 ([#65](https://github.com/szaghi/adam/issues/65)) adds the dissipative terms to FLUME: viscosity and heat
conduction (Navier–Stokes) on Euler and MHD, resistivity (Ohmic diffusion) on MHD, and the no-slip walls. Each
coefficient is read either as a coefficient or as its dimensionless number (#49 N4); the keys are listed in
[input](../../input#dissipative-terms-issue-65-m4), the fluxes in
[numerics](../../numerics#dissipative-fluxes-navier–stokes), the walls in
[boundary conditions](../../boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls).

The legs grow from the contract (DC) through exact linear and steady solutions on uniform grids (VV-1 to VV-6) to the
2:1 AMR and forest seams (VV-7, VV-8), and end with the run-by-run comparison of the two backends. Every leg runs on
the CPU and on FNL, on 2 ranks.

| Case | ID | What it verifies |
|---|---|---|
| [Input contract](./input-contract) | DC | 25 input cases: invalid combinations refused with a message, implied coefficients logged exactly, zero coefficients bitwise inert, dimensional keys converted back exactly |
| [Viscous waves](./viscous-waves) | VV-1, VV-2 | order 4 and 2 of the viscous flux on a shear wave (along $x$, $y$ and the diagonal), $\mu$ and $k$ together on an acoustic wave, against the exact linear mode |
| [Compressible Couette flow](./couette) | VV-3 | the no-slip walls (isothermal and moving) and the viscous heating, against the exact steady profile |
| [Becker's viscous shock](./becker-shock) | VV-4 | the nonlinear viscous and heat fluxes through a resolved shock, Mach 2 and 3 |
| [Ohmic decay and the visco-resistive Alfvén wave](./ohmic-alfven) | VV-5 | the Ohmic flux (sign, cross terms, coefficient), $\mu$ and $\eta$ together, the viscous fluxes on MHD |
| [Ohmic heating and the divergence constraint](./ohmic-heating) | VV-6 | the Ohmic energy returned to the gas as heat; no divergence added by the curl form (GLM, EGLM) |
| [Conservation across 2:1 seams and forest seams](./seam-conservation) | VV-7 | the dissipative fluxes carried by the reflux and the inter-realm register; forest union equal to the single realm |
| [Accuracy across 2:1 seams](./seam-accuracy) | VV-8 | the order of the composite grid (second), and the reflux as the carrier of the dissipative seam flux |
| [CPU and FNL agree run by run](./cpu-fnl-agreement) | VV `agree` | every one of the 72 VV runs equal on the two backends within $10^{-10}$ of the physical scales |

Run every leg:

```bash
cd src/tests/flume/verification/dissipation && ./check.sh                   # DC, CPU
cd src/tests/flume/verification/viscous && ./check.sh                       # VV-1 to VV-8, CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                 # FNL
./check.sh --leg agree                                                      # CPU against FNL
```

The figures (`viscous.png`, `viscous-fields.png`, `couette.png`) are drawn by `make_doc_figures.py` from the CPU work
directories that `viscous/check.sh` leaves (no `make_doc_runs.sh` group):

```bash
exe/.regression-venv/bin/python src/tests/flume/verification/make_doc_figures.py --out docs/public/flume \
   --only viscous viscous-fields couette
```
