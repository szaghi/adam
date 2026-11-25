# PML reflection at normal incidence

| Cases | Issues | Type | Backends | Figure | Regression |
|---|---|---|---|---|---|
| `pml-reflection` sweep (16 runs) | [#145](https://github.com/szaghi/adam/issues/145) (item 3), [#106](https://github.com/szaghi/adam/issues/106) (open: oblique incidence) | P | CPU, FNL (one point) | yes | no (measurement, not a leg) |

## Case and scope

A plane Gaussian pulse travels along $+x$ in a 1-D domain and meets an absorbing boundary at normal incidence. The
perfectly matched layer is measured against its design reflection: for a layer of thickness $L$ backed by a perfect
conductor, with the damping graded as $\gamma = \gamma_{\max}(d/L)^n$, the continuum reflection coefficient at normal
incidence is

$$R(L) = \exp\Big(-\frac{2}{c}\int_0^L \gamma\,dx\Big) = \exp\Big(-\frac{2\gamma_{\max}L}{(n+1)c}\Big),$$

the same for every frequency, so the reflected pulse is the incident one scaled by $R$. The code samples the profile at
the cell centres, $\gamma_m = \gamma_{\max}(1 - m/(N-1))^n$ for the cell $m = 0 \dots N-1$ from the wall
([profiles](../pml#profiles)), so the design value of the discrete layer is
$R_s = \exp(-2\Delta x\sum_m\gamma_m/c)$; for $N = 20$ it is $2.84\cdot10^{-2}$ against the continuum
$3.11\cdot10^{-2}$ of $L = N\Delta x$. What the layer adds is the reflection of its own grading, which grows with
$\gamma_{\max}\Delta x/c$ and falls as the layer thickens.

The sweep measures:
- the thickness $N = 10, 20, 40, 80$ cells at $\gamma_{\max} = 7.8\cdot10^{10}$ s$^{-1}$, $n = 2$ (the paper's
  absorbing-layer test is the $N = 20$ point);
- $\gamma_{\max}$ from $3.9\cdot10^{10}$ to $6.24\cdot10^{11}$ s$^{-1}$ at $N = 20$;
- $N = 10, 20, 40, 80$ at one design value, $R = 9.4\cdot10^{-7}$ ($\gamma_{\max}L$ fixed): what is left is the
  grading reflection;
- the CFS layer at $N = 20$, with $\kappa_{\max} = 1$, $\alpha_{\max} = 0$ (the classical layer through the CFS path)
  and with the $\kappa_{\max} = 5$, $\alpha_{\max} = 2.5\cdot10^8$ s$^{-1}$ of `rmf-pml-amr`;
- the boundaries without a layer: the first-order Silver–Müller condition, Neumann, PEC;
- the noise floor of the measurement: a wall 0.40 m away, beyond reach in the measured time.

The reflection against the angle of incidence is not measured here ([#106](https://github.com/szaghi/adam/issues/106)).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | Maxwell, SI, no source | [models](../models#models) |
| Space, time | centred finite differences of order 6; SSP Runge–Kutta (5,4), CFL 0.5 | [numerics](../numerics#time-integration) |
| PML | `classic` ($n = 2$ unless stated) or `CFS`, on both $x$ faces, a PEC wall behind; the split variables as SELAH companions; step cap $\Delta t \le 0.1/\gamma_{\max}^{\rm eff}$ | [`adam_prism_pml_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/common/adam_prism_pml_object.F90); [absorbing layers](../pml#perfectly-matched-layer) |
| Silver–Müller | first-order impedance ghosts, normal components solved for zero divergence | [boundary conditions](../boundary-conditions#silver-muller) |
| Neumann | mirror ghosts of every component | [`adam_prism_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/prism/cpu/adam_prism_cpu_object.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | interior $[-0.08, 0.08]$ m, plus the layer of $N$ cells on each $x$ face; $y$, $z$: 6 cells, periodic |
| Grid | one block, $\Delta x = 1$ mm (`ni` = $160 + 2N$) |
| Initial condition | [Gaussian pulse](../initial-conditions#gaussian-pulse) along $+x$: centre 0, $\sigma = 4$ mm, cut at $5\sigma$, $B_z$ peak 10 T ($9.922$ T at the cell centres), $E_y = cB_z$ |
| Time | until the pulse reflected by the wall behind the layer has its centre at $x = -0.04$ m ($2(0.08 + L) + 0.04$ m of travel; a CFS layer adds its $2(\kappa_{\max}-1)L/3$); $\Delta t = 8.34\cdot10^{-13}$ s, lower where the PML cap binds |
| Ranks | 1 |

**Measurement** (`pml_reflection.py`). $A_r/A_0$ is the largest $|B_z|$ in the window $[-0.04, 0.04]$ m over the
checkpoints (every 4 steps, 1 mm of travel) after $t_{in} = 0.07\,\mathrm{m}/c$, when the incident pulse and the first
centimetre of its dispersive trail have left the window, over the initial peak. No reflection reaches the window
before $0.10\,\mathrm{m}/c$. $(E_r/E_0)^{1/2}$ is the same with the electromagnetic energy in the window: it equals
$A_r/A_0$ for a scaled copy of the pulse and differs for a distorted reflection.

## Input files

| File | Role |
|---|---|
| `src/tests/prism/verification/pml-reflection/make_pml_reflection.py` | writes one input (`--boundary pml\|cfs\|silver-muller\|neumann\|pec\|far`, `--cells`, `--gamma-max`, `--exponent`, `--alfa-max`, `--k-max`) |
| `src/tests/prism/verification/pml-reflection/pml_reflection.py` | the measurement; `--profiles` keeps $B_z(x, t)$ for the figure |
| `src/tests/prism/verification/pml-reflection/sweep.sh` | runs the 16 cases (`--jobs`, `--cases`), writes `results.txt`, `results.json` |
| `src/tests/prism/verification/make_paper_figures.py` | the figures from the sweep |

```bash
src/tests/prism/verification/pml-reflection/sweep.sh /tmp/pml --jobs 16                     # CPU, about 15 min on 16 cores
PRISM_EXE=$PWD/exe/adam_prism_fnl src/tests/prism/verification/pml-reflection/sweep.sh /tmp/pml-fnl --cases "pml-n20 silver-muller far"
```

PRISM allocates its block pool from the free memory of the node; set `ADAM_MEMORY_FRACTION` (e.g. `0.02`) when the
runs share a node, or they are killed for memory.

## Results

| Boundary | $N$ | $\gamma_{\max}$ [s$^{-1}$] | $R(L)$ | $R_s$ (sampled) | $A_r/A_0$, CPU | $(E_r/E_0)^{1/2}$, CPU | $A_r/A_0$, FNL |
|---|---|---|---|---|---|---|---|
| classical PML | 10 | $7.8\cdot10^{10}$ | $1.77\cdot10^{-1}$ | $1.60\cdot10^{-1}$ | $1.607\cdot10^{-1}$ | $1.603\cdot10^{-1}$ | |
| classical PML | 20 | $7.8\cdot10^{10}$ | $3.11\cdot10^{-2}$ | $2.84\cdot10^{-2}$ | $2.846\cdot10^{-2}$ | $2.843\cdot10^{-2}$ | $2.846\cdot10^{-2}$ |
| classical PML | 40 | $7.8\cdot10^{10}$ | $9.70\cdot10^{-4}$ | $8.88\cdot10^{-4}$ | $8.878\cdot10^{-4}$ | $8.875\cdot10^{-4}$ | |
| classical PML | 80 | $7.8\cdot10^{10}$ | $9.41\cdot10^{-7}$ | $8.62\cdot10^{-7}$ | $9.6\cdot10^{-7}$ | $1.0\cdot10^{-6}$ | |
| classical PML | 20 | $3.9\cdot10^{10}$ | $1.77\cdot10^{-1}$ | $1.69\cdot10^{-1}$ | $1.688\cdot10^{-1}$ | $1.686\cdot10^{-1}$ | |
| classical PML | 20 | $1.56\cdot10^{11}$ | $9.70\cdot10^{-4}$ | $8.08\cdot10^{-4}$ | $8.15\cdot10^{-4}$ | $8.08\cdot10^{-4}$ | |
| classical PML | 20 | $3.12\cdot10^{11}$ | $9.41\cdot10^{-7}$ | $6.53\cdot10^{-7}$ | $4.69\cdot10^{-5}$ | $6.81\cdot10^{-5}$ | |
| classical PML | 20 | $6.24\cdot10^{11}$ | $8.9\cdot10^{-13}$ | $4.3\cdot10^{-13}$ | $9.42\cdot10^{-5}$ | $1.36\cdot10^{-4}$ | |
| classical PML | 10 | $6.24\cdot10^{11}$ | $9.41\cdot10^{-7}$ | $4.35\cdot10^{-7}$ | $4.39\cdot10^{-4}$ | $6.20\cdot10^{-4}$ | |
| classical PML | 40 | $1.56\cdot10^{11}$ | $9.41\cdot10^{-7}$ | $7.88\cdot10^{-7}$ | $5.63\cdot10^{-6}$ | $8.08\cdot10^{-6}$ | |
| CFS, $\kappa_{\max} = 1$, $\alpha_{\max} = 0$ | 20 | $7.8\cdot10^{10}$ | $3.11\cdot10^{-2}$ | $2.84\cdot10^{-2}$ | $2.846\cdot10^{-2}$ | $2.843\cdot10^{-2}$ | |
| CFS, $\kappa_{\max} = 5$, $\alpha_{\max} = 2.5\cdot10^{8}$ | 20 | $7.8\cdot10^{10}$ | — | — | $2.631\cdot10^{-2}$ | $2.696\cdot10^{-2}$ | |
| Silver–Müller (first order) | — | — | — | — | $3.96\cdot10^{-2}$ | $6.01\cdot10^{-2}$ | $3.96\cdot10^{-2}$ |
| Neumann | — | — | — | — | $3.96\cdot10^{-2}$; then 0.987 (re-emitted by the opposite face, below) | $6.01\cdot10^{-2}$; then 0.999 | |
| PEC | — | — | 1 | 1 | $1.0007$ | $1.0000$ | |
| none within reach (floor) | — | — | — | — | $3.2\cdot10^{-7}$ | $4.4\cdot10^{-7}$ | $3.2\cdot10^{-7}$ |

*Measured: ADAM `0bcdf870` + the night-2026-11-22 tree (HEAD `95ffd728` with the working tree of 2026-11-23), CPU
gfortran 16.2.0 + Open MPI 5.0.11 (`iac_gnu`) on medooza (gnode25), FNL nvfortran 26.9, CUDA 12.9 (`iac_nvf`) on
NVIDIA A30 (gnode23), 2026-11-23 (`pml-reflection/sweep.sh`, jobs 1369621 CPU, 1369626 FNL). FNL reproduces the three CPU values to every printed digit, both measures.*

![PML reflection against the layer thickness and the damping rate](/prism/pml-sweep.png)

*Reflected amplitude against the thickness at $\gamma_{\max} = 7.8\cdot10^{10}$ s$^{-1}$ and at a fixed design
$R = 9.4\cdot10^{-7}$ (left), and against $\gamma_{\max}$ at $N = 20$ (right), with the design reflection (line:
continuum; crosses: the sampled profile); the grid-scale wave of Silver–Müller and Neumann and the noise floor as horizontal lines.*

## Findings

- **The layer reflects what its design says, as long as the grading is resolved.** Up to $R_s \approx 10^{-3}$ the
  measured amplitude is the sampled design value to 1 % ($N = 10, 20, 40$ at $7.8\cdot10^{10}$ s$^{-1}$; $3.9\cdot10^{10}$
  and $1.56\cdot10^{11}$ s$^{-1}$ at $N = 20$), and the reflected pulse is a scaled copy of the incident one (the
  amplitude and the energy measures agree to 0.1 %). The paper's 20-cell layer reflects $2.85\cdot10^{-2}$ of the
  incident amplitude: the sampled profile integrates to 2.6 % more damping than the continuum one, whose $R(L)$ is
  $3.11\cdot10^{-2}$.
- **Beyond, the grading sets the floor.** A steeper profile (larger $\gamma_{\max}$ at fixed $N$) stops paying at
  $R \approx 5\cdot10^{-5}$ for $N = 20$: from $3.12\cdot10^{11}$ to $6.24\cdot10^{11}$ s$^{-1}$ the reflection doubles
  while the design value falls by six decades, and it arrives first (the inner part of the layer), with the opposite
  sign and a distorted shape (the energy measure 1.45 times the amplitude one). At a fixed design
  $R = 9.4\cdot10^{-7}$ the grading reflection falls from $4.4\cdot10^{-4}$ ($N = 10$) to $4.7\cdot10^{-5}$ (20),
  $5.6\cdot10^{-6}$ (40) and $9.6\cdot10^{-7}$ (80, the design value): about $N^{-3}$, the order of a quadratic profile.
  The best measured values are $4.7\cdot10^{-5}$ with 20 cells and $5.6\cdot10^{-6}$ with 40.
- **CFS.** With $\kappa_{\max} = 1$, $\alpha_{\max} = 0$ the CFS path reproduces the classical layer to every printed
  digit. With the `rmf-pml-amr` parameters it reflects $2.63\cdot10^{-2}$: the stretching delays the reflection
  (the pulse crosses the layer at $c/\kappa$) and leaves the amplitude near the classical one.
- **Silver–Müller is not exact at normal incidence on this grid.** The continuum condition reflects nothing at normal
  incidence; the discrete one returns 0.36 % of the energy, $(E_r/E_0)^{1/2} = 6.0\cdot10^{-2}$, as a grid-scale
  (odd–even) wave packet of peak $|B_z| = 3.96\cdot10^{-2}$ of the incident one that travels back at 2.2 times the speed
  of light (the group velocity of the $k\Delta x = \pi$ mode of the sixth-order centred derivative). In the energy
  measure the 20-cell classical layer at the paper's $\gamma_{\max}$ reflects 2.1 times less (4.5 times less energy),
  the 40-cell one 68 times less.
- **Neumann does not reflect the pulse.** The mirror ghosts let 99.6 % of the energy leave through the face, and send
  back the same grid-scale packet as Silver–Müller (0.36 % of the energy; the window energy of the two runs agrees to four digits
  from step 360 to step 560). At the opposite Neumann face that packet turns into a smooth pulse of nearly the
  incident amplitude ($A/A_0 = 0.987$) travelling in $+x$: the energy in the domain grows from 0.36 % to 99.7 % of the
  initial one within 80 steps. A Neumann face is not energy-stable for the centred scheme; the "large reflection" of a
  Neumann boundary in a 1-D run is this re-emission at the opposite face, not a reflection.
- The floor of the measurement, $3.2\cdot10^{-7}$, is the dispersive trail of the pulse (its high-wavenumber content,
  which the centred scheme carries backwards), not round-off.
