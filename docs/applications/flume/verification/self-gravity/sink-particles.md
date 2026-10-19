# Sink particles

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| SI (si1 to si11) | Self-gravity | M12 ([#94](https://github.com/szaghi/adam/issues/94)) | CPU, FNL | yes | `sinks-shu` |

## Case and scope

`[sinks]` adds point masses that act through self-gravity only ([numerics](../../numerics#sink-particles)):
- their mass is deposited into the Poisson source by cloud-in-cell;
- the mesh gravity is interpolated back with the same weights;
- they are advanced by kick-drift-kick once per step.

They form where the gas meets the checks of Federrath et al. (2010) and accrete the gas above the threshold density of
their accretion radius. Two bound, converging sinks merge. With runtime regridding they form on the finest level
allowed and keep their neighbourhood there. The list is replicated: every rank holds every sink.

What the case catches:
- **a sink that is not in the Poisson source:** a missing deposit or a wrong normalisation, or a sink that feels the gas
  but pulls none of it;
- **a self-force, or a mesh force far from Newton's:** a sink at rest that moves, or a binary whose period, energy or
  angular momentum drifts;
- **mass that is not conserved:** by the accretion, the creation, the merger or a regrid;
- **a creation that is too eager or too strict:** a Jeans-stable box that forms a sink, or a collapse that forms none or
  several;
- **an accretion rate that is wrong:** against the similarity solution of Shu (1977);
- **a regrid that coarsens a sink's neighbourhood,** or a sink formed on a coarse level;
- **a restart that loses the sinks,** a refused input that runs anyway, or a backend that disagrees.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler` with `[gravity] self_gravity`, isolated (James) or periodic, and `[sinks]` | [`adam_flume_sinks_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_sinks_object.F90) |
| Gravity | the CIC deposit into the packed Poisson source, the centred difference of $\Phi$ gathered with the same weights (Hockney & Eastwood 1988) | `deposit_sinks`, `compute_sink_accelerations` in [`adam_flume_common_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_common_object.F90) |
| Sink step | kick-drift-kick on the committed state once per step; $\Delta t \le \tfrac12 h/|v|$ and $\tfrac12\sqrt{h^3/(Gm)}$ | `kick_sinks`, `drift_sinks`, `sinks_dt_limit` |
| Creation, accretion, mergers | Federrath et al. (2010): the threshold $\rho_\mathrm{thr} = \pi c_s^2/(4 G r_\mathrm{acc}^2)$, the potential minimum, the converging flow, Jeans instability, boundedness; threshold accretion of the bound infalling gas | `create_sinks`, `accrete_sinks`, `merge_sinks`; FNL `accrete_sinks_gpu` |
| Runtime AMR | creation on the finest level allowed, $r_\mathrm{acc}$ in its cells, the neighbourhood kept there at every regrid | `mark_runtime` (`sink_region`), `sink_accretion_radius` |
| Space, time | WENO5 characteristic flux splitting, SSP-RK(3,3), the positivity limiter | as [the polytrope](./polytrope) |

## Initial and boundary conditions

Every case is the [polytrope](./polytrope) box with its numerics: $[-1, 1]^3$, $32^3$ cells (blocks of 16), $G = 1$,
$\gamma = 2$ and $r_\mathrm{acc} = 2.5$ cells, unless the table says otherwise.

| Leg | Setup |
|---|---|
| si1 | one sink of mass 1 at rest at $(0.15, -0.1, 0.05)$, off every symmetry plane, in a uniform periodic background $\rho = p = 10^{-3}$, to $t = 0.5$ |
| si2 | the equal-mass circular binary: masses 1 at $(\mp0.25, 0, 0)$, velocities $(0, \mp1, 0)$ (separation 8 cells), isolated, near vacuum ($\rho = p = 10^{-6}$), one orbit ($t = \pi/2$) |
| si3 | `[sinks]` without self-gravity; with `mhd-ideal`; lists of different lengths; a non-positive mass; a sink outside the domain |
| si4 | the si1 box 20 steps in one run, and 10 + 10 steps across a restart |
| si5 | a sink of mass 0.1 at rest in a dense, cold periodic background ($\rho = 1$, $p = 10^{-3}$, $\rho_\mathrm{thr} \approx 0.06$), to $t = 0.3$ |
| si6 | the polytrope at a tenth of its pressure, kept isothermal at $T = 0.01$ by strong cooling ([radiative cooling](../euler/radiative-cooling)), closed walls, creation on, no initial sink, to $t = 0.6$ |
| si7 | a static periodic box at $\gamma = 1.4$, $\rho = 1$, $p = 0.0185$ ($1.2\,\rho_\mathrm{thr}$), creation on, 20 steps |
| si8 | two sinks of mass 1 at $\mp0.04$ (1.28 cells apart) approaching at $\pm0.1$ |
| si9 | the si2, si5, si6 and si8 runs of CPU and FNL |
| si10 | runtime AMR (a Jeans marker of 8 cells, a regrid every 4 steps, `regrid_prolongation = conservative-state`). The si1 sink on a base of $4^3$ blocks of $8^3$ with one level above it, 40 steps. The si6 collapse at $16^3$ with two levels above the base |
| si11 | the overdense singular isothermal sphere of Federrath et al. (2010, §3.4): $\rho = A a^2/(4\pi G r^2)$ (cell means), $A = 29.3$, $a = 0.5$, $\gamma = 1.001$, at rest, a seed sink of mass $10^{-6}$ at the centre (a cell corner), closed walls, isolated, to $t = 0.25$ |

**The singular isothermal sphere.** The `singular-isothermal-sphere` initial condition sets the cell mean of $1/r^2$,
by the divergence theorem ($\nabla\cdot(\mathbf{x}/r^2) = 1/r^2$) as the flux through the cell's faces with $8\times8$
Gauss–Legendre points per face. Sampled at the cell centres, $1/r^2$ misses the mass of the unresolved core: 2 % of the
$32^3$ box, so that at $A = 2$ the gravity falls short of the pressure and the hydrostatic sphere expands. Federrath et
al. run the overdense sphere, whose collapse is dominated by gravity. Shu's (1977) A-family gives its central mass,
$M = m_0(A)\,a^3 t/G$, from his similarity equations, which `si_oracle.py` integrates: $m_0 = 0.9757$ as $A \to 2^+$
(Shu's 0.975) and $m_0 = 138.94$ at $A = 29.3$, below the pressureless $(2\sqrt2/\pi)A^{3/2} = 142.8$. A shell of
radius $r$ reaches the centre at about $(\pi/2)\,r/(a\sqrt{2A})$, so by $t = 0.25$ the shells within 0.6, ahead of the
walls' rarefaction. $\gamma = 1.001$ keeps the gas isothermal within 1 % over the run's density range. Cooling cannot:
it only removes heat, and expanded gas would stay below its floor.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/sinks/make_si.py` | writes the inputs (`--case rest|binary|accrete|collapse|stable|merge|shu`, the sink lists, `--amr LEVELS`, `--base-levels`, `--overdensity`) |
| `src/tests/flume/verification/sinks/check.sh` | runs the legs (`--legs si1,...,si11`, default every leg but si9; `FLUME_EXE` selects the executable, default `exe/adam_flume_cpu_fft`; si9 compares the work directories of both) |
| `src/tests/flume/verification/sinks/si_oracle.py` | the oracles: Kepler, the orbital energy and angular momentum, the mass and momentum budgets, the AMR levels of the last checkpoint, Shu's similarity equations |
| `src/tests/flume/regression/sinks-shu/input.ini` | si11 at $16^3$, 30 steps |
| `make_doc_runs.sh <runs> sinks`, `make_doc_figures.py --only sink-particles` | the figure (4 runs) |

```bash
cd src/tests/flume/verification/sinks && ./check.sh                                    # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl_fft ./check.sh                       # FNL
./check.sh --legs si9                                                                 # CPU against FNL
```

## Results

![Sink particles: the binary, the Shu accretion rate and the collapse on runtime AMR](/flume/sink-particles.png)

*Top: the circular binary over one orbit, its separation and orbital energy within 2 % (the mesh force is softened at
about two cells). Bottom left: the accretion rate of the overdense isothermal sphere at $32^3$ and $64^3$ against Shu's
$m_0(A)\,a^3/G$. Bottom right: the collapse on runtime AMR, the sink formed on the finest level, its neighbourhood kept
there.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| si1 | the sink's displacement over $t = 0.5$ | $10^{-3}$ cells | $1.36\cdot10^{-4}$ | $1.36\cdot10^{-4}$ |
| si1 | gas plus sink momentum, over the momentum they exchange | $10^{-3}$ | $3.5\cdot10^{-5}$ | $3.5\cdot10^{-5}$ |
| si2 | the separation; the period against Kepler's at the mean separation | 3 %; 2 % | $1.36\cdot10^{-2}$; $7.5\cdot10^{-3}$ | the same |
| si2 | the orbital energy $K - Gm_1m_2/d$; the angular momentum; the centre of mass | 1 %; 1 %; $10^{-10}$ | $1.7\cdot10^{-3}$; $7.2\cdot10^{-4}$; $9.3\cdot10^{-15}$ | $1.7\cdot10^{-3}$; $7.2\cdot10^{-4}$; $7.8\cdot10^{-15}$ |
| si3 | the five invalid inputs | refused | refused | refused |
| si4 | the histories of 20 steps against 10 + 10 across a restart | bitwise | bitwise | bitwise |
| si5 | gas plus sink mass; the sink grows | $10^{-13}$ | exact; 9.7 % of the sink | $2\cdot10^{-16}$; the same |
| si5 | gas plus sink momentum, over the momentum exchanged | 5 % | $1.4\cdot10^{-2}$ | $1.4\cdot10^{-2}$ |
| si6 | sinks formed; distance from the centre; mass | 1; one cell; $10^{-13}$ | 1; $5.0\cdot10^{-3}$ cells; $1.1\cdot10^{-16}$ | 1; the same; $2.2\cdot10^{-16}$ |
| si7 | sinks formed in the Jeans-stable box | none | none | none |
| si8 | sinks left; their mass and momentum against the sum | 1; $10^{-13}$, $10^{-10}$ | 1; exact, $2.4\cdot10^{-16}$ | 1; exact, $3.9\cdot10^{-16}$ |
| si9 | FNL against CPU: positions (cells), velocities, masses, every row | $10^{-8}$ | (the same comparison) | si2 $1.8\cdot10^{-13}$, si5 $3.1\cdot10^{-14}$, si10 $1.2\cdot10^{-12}$; si6, si10 collapse, si11 masses $\le 3.4\cdot10^{-11}$ |
| si10 | the blocks within $r_\mathrm{acc} + 2$ cells of the sink on the finest level; no family refined beyond the reach | all; none | 6 of 6; none (56 coarse blocks beyond) | the same |
| si10 | the sink at rest on the composite grid; gas plus sink mass across the regrids | $10^{-3}$ cells; $10^{-13}$ | $7.9\cdot10^{-5}$; exact | the same |
| si10 | the collapse on AMR: one sink, on the finest level, at the centre; mass | 1; all; one cell; $10^{-13}$ | 1; 8 of 8; $4.2\cdot10^{-5}$ cells; $2.2\cdot10^{-16}$ | the same |
| si11 | the accretion rate over $m_0(29.3)\,a^3/G$, $t \in [0.1, 0.25]$; mass | 6 %; $10^{-13}$ | 0.959; $2.2\cdot10^{-16}$ | 0.959; $6.7\cdot10^{-16}$ |
| si11 at $64^3$ (figure run) | the accretion rate | — | 0.985 | — |

Related cases: [the polytrope](./polytrope) is the box and the collapse of si6 and si10, [radiative
cooling](../euler/radiative-cooling) keeps si6 isothermal, and [self-gravity under runtime AMR](./runtime-amr) is the
regrid that si10 drives.

## Findings

- **The sinks are part of the Poisson source.** A sink at rest stays put within $1.4\cdot10^{-4}$ cells, and gas and sink
  exchange momentum equal and opposite to $3.5\cdot10^{-5}$ of it. The binary follows Kepler within the mesh force's
  error: at 8 cells the force, softened at about two cells, is 0.7 % above Newton's. That is a resolution property, not
  a defect.
- **Mass is exact, momentum is not quite.** Accretion, creation, mergers and regrids conserve the mass to round-off. Near
  a point mass, the gas feels the potential through face gradients and the sink through the CIC gather; the two differ by
  about 1 % of the exchanged momentum (si5).
- **Sinks belong on the finest level allowed.** A sink formed on a coarse level had its neighbourhood refined at the
  next regrid. That interpolated the cold, Mach-8 infall onto the point mass into the new children, and 42 of them lost
  their pressure. It also halved the accretion radius from one step to the next. Forming only on `[amr] max_level`, with
  $r_\mathrm{acc}$ in its cells, is Federrath et al.'s rule and removes both effects. The collapse also needs
  `regrid_prolongation = conservative-state`.
- **A composite grid pulls on a sink.** A sink at rest whose neighbourhood is refined asymmetrically feels its own
  potential through the 2:1 seams. The force is constant, about $3\cdot10^{-2}$ finest cells over $t = 0.5$, against
  $1.4\cdot10^{-4}$ on a uniform grid. It is a truncation effect of the composite Poisson operator, bounded by the 40
  steps of si10.
- **The singular sphere needs cell means.** Sampled at the cell centres, $1/r^2$ loses the mass of the unresolved core:
  2 % of the box, 5 % of the gravity at $r = 0.5$. The hydrostatic sphere then expands instead of collapsing.
- **Shu's rate converges at first order.** The accretion rate is 4.1 % below Shu's at $32^3$ and 1.5 % below at
  $64^3$. Federrath et al. report about 3 % at their resolution. The marginal sphere ($A = 2$) is not a test at these
  resolutions: the threshold accretion zone stores $(\pi^2/3)\,a^2 r_\mathrm{acc}/G$, more than the sink takes before the
  box runs out.
- **CPU and FNL agree to round-off,** within $10^{-11}$, on every compared run. Comparisons of a sink at a symmetric centre
  use its mass and its distance from the centre only: which of 8 tied cells wins, or a velocity of $10^{-11}$, is
  round-off.
