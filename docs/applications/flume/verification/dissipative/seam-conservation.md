# Conservation across 2:1 seams and forest seams

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-7 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | `seam-viscous-amr` |

## Case and scope

A diagonal wave of large amplitude crosses 2:1 coarse-fine faces (one realm with its centre refined) and inter-realm
seams (the same domain split into two realms). The conserved integrals are judged on the refined realm; the fields of
the forest are judged against the single realm.

What the case catches:
- **dissipative fluxes missing from the seam registers:** the dissipative fluxes are added to the inviscid ones before
  the seam accumulation, so the Berger–Colella reflux and the inter-realm register must carry them. At a 2:1 face the
  coarse and the fine side compute different viscous, heat and Ohmic fluxes, and without the register the grid leaks;
- **unfilled edge and corner ghosts:** the dissipative stencils reach further than WENO's: the tangential derivatives
  read the edge and corner ghosts, which no inviscid stencil touches.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Euler with $\mu = 0.01$, $k = 0.02$; ideal MHD with GLM ($c_h = 3$), $B_n = 1$, $\mu = 0.01$, $\eta = 0.01$; $\gamma = 1.4$ | [models](../../models#compressible-euler), [GLM](../../models#divergence-control-mixed-glm) |
| Dissipative fluxes | `dissipative_order = 4`, added to the inviscid fluxes before the seam accumulation | `add_dissipative_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 flux splitting, characteristic variables | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Time | `runge-kutta-ssp-54`, CFL 0.2 | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| 2:1 reflux | Berger–Colella register of the seam fluxes, `[numerics] reflux` on or off | `accumulate_seam_fluxes`, `reflux_stage`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `accumulate_coarse_flux`, `accumulate_fine_flux`, [`adam_flux_register_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_flux_register_object.F90); [numerics](../../numerics#adaptive-mesh-refinement) |
| Coarse-fine ghosts | the library's 2:1 ghost fill | [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Forest seams | `coupling = mirror` (1:1) or `refined` (2:1), `coupling_cadence = stage_coincident`; seam ghost centres wrapped by one period along a tangential axis periodic on both realms | `fill_seam_from_peer_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `wrap_periodic_tangential`, [`adam_forest_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_forest_object.F90); [inter-realm seam](../../boundary-conditions#inter-realm-seam) |
| Conservation history | volume integrals written every step | `compute_conservation`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification) at 45°, $A = 0.01$ (large, so the seams carry a real flux): the acoustic wave on Euler (density, normal velocity and temperature all vary), the magnetic (Alfvén) wave on MHD-GLM; $\rho_0 = p_0 = 1$, wavelength $1/\sqrt2$ |
| Final time | $t = 0.05$ |

| Leg | Grid | Boundaries | Runs |
|---|---|---|---|
| refined box (`conserve-*-quadtree`) | quadtree, $4 \times 4$ blocks of $12^2$ cells (48 per side at the base level), the blocks whose centroid lies in $[0.25, 0.75]^2$ refined 2:1 (8 coarse-fine faces) | `periodic` in $x$ and $y$ | with reflux; without; without reflux and without the dissipative coefficients (the ideal twin) |
| forest, `mirror` | the single realm split at $x = 0.5$ into two realms of $4 \times 4$ blocks of $6 \times 12$ cells, quadtree, glued 1:1 | `extrapolation` (outflow) in $x$, `periodic` in $y$ | single realm; forest |
| forest, `refined` | the single run refines $x > 0.5$; realm 1 the coarse half, realm 2 the fine half at one more uniform level, glued 2:1; octree ($nk = 4$), because the `refined` coupling needs ratio 2 along every axis, null $z$ included | `extrapolation` in $x$, `periodic` in $y$ | single realm; forest |

An octree with null $z$ gave the same drifts as the quadtree to the printed digits at four times the cells, so the
refined-box leg keeps the quadtree; the octree's 2:1 seams are exercised by the refined forest.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/waves.py` | `make --refine ... --reflux`, `split` (the two-realm forest), `conserve` (the drift oracle) |
| `src/tests/flume/verification/multirealm/multirealm_oracle.py` | the forest union against the single realm (`--fields-only`) |
| `src/tests/flume/verification/viscous/check.sh --legs vv7` | runs the 6 refined-box runs and the 8 forest runs (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/make_doc_figures.py --only seam-conservation` | the figure, from the histories that `viscous/check.sh --legs vv7` keeps (CPU) |
| `src/tests/flume/regression/seam-viscous-amr/input.ini` | the frozen Alfvén wave, $A = 0.01$, MHD-GLM with $\mu$ and $\eta$, $48^2$ with the centre refined 2:1, $t = 0.02$ |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --legs vv7             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs vv7       # FNL
```

## Results

![VV-7: drift of the integrals across seams](/flume/seam-conservation.png)

*The drift of the mass and energy integrals per step, CPU. Left and centre: the acoustic wave on Euler and the Alfvén wave on MHD-GLM across a 2:1 box; with the reflux the integrals stay at round-off, without it they drift, and so does the ideal twin without the dissipative coefficients, so the drift is the seam, not the dissipative fluxes. Right: the forests (1:1 mirror and 2:1 refined, outflow in $x$) against their single-realm runs: the summed integrals of the two realms agree to round-off at every step.*

Oracle (`waves.py conserve`). The drift of each conserved integral is divided by $\int|q|\,dV$, the scale of its
round-off (the momentum of a sine wave integrates to zero, so its own value is no scale; $\psi$, sourced only by the
discrete $\nabla\cdot\mathbf B$, is measured against $c_h\int|\mathbf B|\,dV$):
- with reflux, every drift $\le 10^{-13}$;
- without reflux, the largest drift $\ge 10^{-10}$ (the seams are exercised);
- the leak without reflux differs from the ideal twin's by $\ge 10^{-10}$: the dissipative fluxes cross the seams
  unmatched, and the first item shows the register removes that too.

| Case | With reflux (bound $10^{-13}$) | Without ($\ge 10^{-10}$) | Dissipative minus ideal leak ($\ge 10^{-10}$) |
|---|---|---|---|
| Euler | $\le 1.1\cdot10^{-17}$ | $3.0\cdot10^{-5}$ | $2.0\cdot10^{-6}$ |
| MHD-GLM | $\le 1.5\cdot10^{-16}$ | $5.9\cdot10^{-5}$ | $5.1\cdot10^{-6}$ |

CPU and FNL (2 ranks) are identical to the printed digits except the round-off of the conserved runs, at most
$1.1\cdot10^{-16}$ on FNL.

Forests: the union must reproduce the single realm, fields within $10^{-10}$ relative.

| Forest | Bound | CPU | FNL |
|---|---|---|---|
| `mirror` and `refined`, Euler and MHD-GLM | $10^{-10}$ | largest difference $8.7\cdot10^{-12}$ | $4.4\cdot10^{-12}$ |

Related cases: the accuracy across the same seams is [VV-8](./seam-accuracy); the inviscid conservation across AMR
seams is [Euler AMR conservation](../euler/amr-conservation); the forest machinery is [multi-realm runs](../euler/multi-realm).

## Findings

- The forest comparison is not bitwise, unlike MV-14's piecewise-constant states: each realm evaluates the sine from
  its own block origins, one ulp apart at step 0 ($3\cdot10^{-15}$). The largest difference after the run sits at the
  outer $x$ faces, not at the seam.
- **This leg found a library defect.** A seam slab's edge and corner ghosts past a tangential boundary of the realm lie
  outside the peer domain, and were left to the boundary condition; along a periodic axis no condition fills them,
  and they kept their initial values. WENO never reads them, the dissipative tangential derivatives do: the mirror
  split went NaN at the first step (on a single realm of the same half domain the run is clean). The seam enumeration
  now wraps those ghost centres by one period when the axis is periodic on both realms
  ([forest guide](/concepts/forest#seams-across-ranks-issue-40)).
