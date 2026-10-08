# FLUME

**FLUME** — **F**luid **L**orentz-coupled **U**nsteady **M**agnetohydrodynamic **E**quations — is a compressible magnetohydrodynamics (MHD) application built on the ADAM framework. It solves the equations of an electrically conducting, compressible fluid interacting with its own magnetic field through the Lorentz force.

> **Status: development.** Milestone M1 of [issue #35](https://github.com/szaghi/adam/issues/35) delivered the
> **compressible Euler** equations on both backends; milestone M2 of [issue #41](https://github.com/szaghi/adam/issues/41)
> adds **ideal MHD** with mixed GLM divergence cleaning, verified on both backends including AMR and multi-realm runs;
> milestone M3 of [issue #47](https://github.com/szaghi/adam/issues/47) adds the **`weno-riemann`** scheme (Euler LLF,
> HLL, HLLC; MHD LLF, HLL, HLLD), the energy-consistent **EGLM** cleaning and a **positivity limiter**.
> FLUME supersedes CHASE, which is deprecated.

## Physical Models

| Model | `[physics] physical_model` | `[mhd] divergence_control` | State |
|-------|----------------------------|----------------------------|-------|
| Compressible Euler | `euler` | — | $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E)^\top$, `nv = 5` |
| Ideal MHD, GLM cleaning | `mhd-ideal` | `glm` | $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E, B_x, B_y, B_z, \psi)^\top$, `nv = 9` |
| Ideal MHD, EGLM cleaning | `mhd-ideal` | `eglm` | as GLM, `nv = 9`, with $\psi^2/2$ in $E$ and the Derigs et al. (2018) nonconservative sources |
| Ideal MHD, no control | `mhd-ideal` | `none` | the first 8 of the above, `nv = 8`: diagnostic use only (in multi-D the divergence error grows unchecked: the magnetised vortex reaches a negative pressure) |

All models are inviscid with an ideal gas (`[physics] cp, cv`, J/(kg K)); ideal MHD is a perfectly conducting single
fluid (ideal Ohm's law $\mathbf{E} + \mathbf{u} \times \mathbf{B} = 0$). `mhd-ideal` with immersed solids is refused.
Dissipative effects (viscosity, thermal conduction, Ohmic resistivity) are milestone M4
([#65](https://github.com/szaghi/adam/issues/65)): each coefficient is given as itself or as its number (`[physics]
viscosity` or `reynolds`, `conductivity` or `prandtl`, `resistivity` or `magnetic_reynolds`/`lundquist`, a power-law
viscosity). The viscosity and the conductivity make each model compressible Navier–Stokes, the resistivity makes
`mhd-ideal` resistive MHD (4th-order conservative fluxes, or 2nd order, `[numerics] dissipative_order`; no-slip and
isothermal walls; diffusive time-step limit), on both backends.

## Implemented Capabilities

| Area | What is available | INI |
|------|-------------------|-----|
| Space | WENO flux splitting, per-face Roe eigenvectors, per-wave local Lax-Friedrichs; reconstruction in characteristic or conservative variables; orders `weno-u-3` to `weno-u-9` (the centred `weno-c-*` schemes are refused). MHD: block-diagonal characteristic decomposition, 7×7 Roe–Balsara-normalised eigenvectors (Stone et al. 2008) at the arithmetic mean of the primitive states with $B_n$ as a parameter, plus the $(B_n, \psi)$ block with speeds $\mp c_h$ | `[numerics] scheme_space = weno`, `reconstruction_variables`; `[weno] scheme` |
| Space, `weno-riemann` | WENO interpolation of the face states (characteristic, or primitive) into a Riemann solver: Euler `llf`, `hll`, `hllc`; MHD `llf`, `hll`, `hlld` (HLLD falls back to HLL where its intermediate states are inadmissible, counted in the log); optional 4th/6th-order face-flux correction, switched off by a WENO smoothness sensor | `[numerics] scheme_space = weno-riemann`, `reconstruction_variables`, `riemann_solver`, `flux_correction`, `flux_correction_sensor` |
| Positivity limiter | Cell-based parametrised flux limiter over a first-order Lax–Friedrichs backbone: each stage keeps the density and the pressure positive, and above a tenth of the first-order update; smooth runs bitwise unchanged; at 2:1 AMR seams the seam flux is synchronised at every stage. Euler, MHD without cleaning and EGLM, SSP schemes only; refused with GLM, immersed solids and multi-realm runs | `[numerics] positivity_limiter = none\|cell` |
| Time | Library Runge-Kutta schemes (SSP and low-storage), CFL time step | `[runge_kutta] scheme`; `[time] CFL, it_max, time_max` |
| Boundary conditions | `extrapolation`, `inflow` (primitive state `r, u, v, w, p`; MHD adds `bx, by, bz`, with $\psi = 0$), `wall-inviscid` (MHD: perfectly conducting wall, $u_n$ and $B_n$ odd, the rest and $\psi$ even), `wall-noslip` and `wall-isothermal` (velocity reflected about a tangential wall velocity; temperature mirrored, or set by `wall_temperature`; issue #65 P1), `periodic` (both faces of an axis or neither); the ghosts at a realm edge or corner compose the conditions of its faces (issue #65 P0) | `[bc_{x,y,z}_{min,max}] type`, `wall_u, wall_v, wall_w`, `wall_temperature` |
| Initial conditions | `uniform` (optionally with a seeded perturbation), `isentropic-vortex`, `riemann-problem` (piecewise-constant regions; MHD regions add `bx, by, bz`), `shu-osher` (Euler, along x, y or z); MHD only: `glm-pulse`, `divb-peak`, `mhd-linear-wave`, `mhd-cpaw`, `mhd-vortex`, `orszag-tang`, `mhd-rotor`, `field-loop`, `rotated-riemann`; `linear` (verification of the ghost fills) | `[initial_conditions] type` |
| AMR | Init-time refinement (`amr_iterations` passes) by geometric box, variable gradient or immersed-solid surface; 2:1 coarse-fine faces with stage-weighted conservative reflux (B and $\psi$ included). With markers the tree is an octree (`ratio = 8`; with a null axis, `nk >= 4`) or a quadtree (`ratio = 4`, any `nk`; [#46](https://github.com/szaghi/adam/issues/46), MV-15) | `[amr]`, `[initial_conditions] amr_iterations`, `[numerics] reflux` |
| Immersed boundary | Euler only. Static solids, inviscid wall: distance function, eikonal extrapolation into the solid, cut-cell spacing, solid masks in the Runge-Kutta stages |
| Multi-realm | A forest manifest glues realms through inter-realm seams, across ranks: `mirror` between cells of the same size, `refined` across a 2:1 resolution jump (the AMR seam formulas); per-seam cadence; every model | `[forest]` manifest | `[solids]` |
| Output | XH5F checkpoints (optionally with the auxiliary fields `u, v, w, p, H, a`; MHD adds the derived `pt, beta, bmag, divb`), slices, residuals and conservation histories, restart; MHD adds the div(B) history `<basename>-divb_history.dat` (`it time max_divb l1_divb seam_max_divb`) | `[IO]`, `[slices]` |

Every option value is matched against a fixed list; an unknown value stops the run with a message naming the accepted
spellings.

## Governing Equations

FLUME targets the conservative hyperbolic form of the ideal MHD equations:

$$\frac{\partial \mathbf{q}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{q}) = \mathbf{S}$$

In units where the magnetic permeability is absorbed into $\mathbf{B}$ ($\mathbf{B} \to \mathbf{B}/\sqrt{\mu_0}$), the system reads:

$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = 0$$

$$\frac{\partial (\rho \mathbf{u})}{\partial t} + \nabla \cdot \left[\rho \mathbf{u} \otimes \mathbf{u} + \left(p + \tfrac{1}{2}|\mathbf{B}|^2\right)\mathbf{I} - \mathbf{B} \otimes \mathbf{B}\right] = 0$$

$$\frac{\partial E}{\partial t} + \nabla \cdot \left[\left(E + p + \tfrac{1}{2}|\mathbf{B}|^2\right)\mathbf{u} - (\mathbf{u} \cdot \mathbf{B})\,\mathbf{B}\right] = 0$$

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{B} - \mathbf{B} \otimes \mathbf{u}) = 0$$

with total energy $E = \dfrac{p}{\gamma - 1} + \tfrac{1}{2}\rho|\mathbf{u}|^2 + \tfrac{1}{2}|\mathbf{B}|^2$ (ideal gas closure) and the solenoidal constraint

$$\nabla \cdot \mathbf{B} = 0 .$$

The system is hyperbolic, with seven wave families per direction: two fast magnetosonic, two Alfvén, two slow
magnetosonic and one entropy wave.

**Units.** $\mathbf{B}$ in the state and in every input is $\mathbf{B}_{SI}/\sqrt{\mu_0}$ (magnetic pressure
$|\mathbf{B}|^2/2$); the fluid stays dimensional SI. A field given in Tesla is divided by $\sqrt{\mu_0}$, a Gaussian
one by $\sqrt{4\pi}$, by the user: FLUME applies no conversion.

## Divergence Control

The induction equation preserves $\nabla \cdot \mathbf{B} = 0$ analytically, but not discretely, and a divergence
error feeds back into momentum and energy as a spurious force parallel to $\mathbf{B}$. FLUME uses **mixed GLM
cleaning** (Dedner et al. 2002): a ninth variable $\psi$ couples to the induction equation,

$$\frac{\partial \mathbf{B}}{\partial t} + \nabla \cdot (\mathbf{u} \otimes \mathbf{B} - \mathbf{B} \otimes \mathbf{u}) + \nabla \psi = 0, \qquad
\frac{\partial \psi}{\partial t} + c_h^2 \nabla \cdot \mathbf{B} = -\frac{c_h^2}{c_p^2} \psi ,$$

which transports the error away at the constant speed $c_h$ and damps it; $\psi$ is not part of $E$. It is the
formulation compatible with cell-centred storage, conservative reflux and high-order finite differences; constrained
transport is deferred. **EGLM** (`eglm`, Derigs et al. 2018) adds $\psi^2/2$ to the energy and the nonconservative
sources $-(\nabla\cdot\mathbf{B})(0, \mathbf{B}, \mathbf{u}\cdot\mathbf{B}, \mathbf{u}, 0) - (\mathbf{u}\cdot\nabla\psi)(0, 0, \psi, 0, 1)$, so
that the cleaning is thermodynamically consistent: with it the strongly magnetised blasts stay admissible (with the
positivity limiter on the splitting scheme). The `[mhd]` section (read only with `physical_model = mhd-ideal`):

| Key | Values | Meaning |
|-----|--------|---------|
| `divergence_control` | `glm`, `eglm`, `none` | model variant (state width 9, 9 or 8) |
| `glm_ch` | > 0 | cleaning speed $c_h$, constant and uniform, part of the time-step bound |
| `glm_alpha` | >= 0 | damping $c_h^2/c_p^2 = \alpha\, c_h / L$ |
| `glm_damping_length` | > 0 or `min-cell` | $L$; `min-cell` is the minimum cell spacing of the realm (Mignone & Tzeferacos 2010) |
| `glm_ch_check` | `warning`, `error` | action when the fastest wave outruns $c_h$ |
| `divb_tol`, `divb_error` | >= 0, logical | div(B) monitor: warn (or stop) when `max_divb` exceeds the tolerance; 0 disables it |
| `rho_floor`, `p_floor` | >= 0 | positivity floors with global counters in the log; 0 disables them (a non-positive state then stops the run) |

The GLM keys are required with `glm` and `eglm`; every value is range-checked and an unknown spelling stops the run.

## Source Layout

The layout mirrors PRISM: backend-independent modules in `common/`, one directory per backend.

```
src/app/flume/
├── common/          # Shared across all backends (physics, numerics, BC, IC, I/O, ...)
├── cpu/             # CPU backend (MPI + OpenMP)
└── fnl/             # OpenACC GPU backend (FNL library)
```

Module and file naming follows the ADAM convention: `adam_flume_<name>_object.F90` for types, `adam_flume_<name>_library.F90` for procedure collections, `adam_flume_cpu.F90` / `adam_flume_fnl.F90` for the program entry points.

## Backends

| Backend | Entry point | Key type | Accelerator |
|---------|-------------|----------|-------------|
| CPU | `adam_flume_cpu.F90` | `flume_cpu_object` | MPI + OpenMP |
| FNL | `adam_flume_fnl.F90` | `flume_fnl_object` | MPI + OpenACC (NVIDIA GPU) |

## Building

```bash
fobis build --mode flume-cpu-gnu                         # CPU backend (GNU compiler)
fobis build --mode flume-cpu-gnu-omp                     # CPU backend with OpenMP threads
fobis build --mode flume-fnl-nvf --varset local_nvf      # FNL (OpenACC) backend
```

A run takes the INI file as its argument: `mpirun -np 2 exe/adam_flume_cpu input.ini`.

## Verification and Regression

`src/tests/flume/verification/` holds one `check.sh` per test, each asserting the physics against an oracle:

| Test | Case | Pass criterion |
|------|------|----------------|
| V0 | `unit/`: pointwise Euler library on random states | eigenvector, Jacobian, round-trip and flux-split identities |
| V1 | `sod/`: Sod along x, y, z; reflecting-wall variant | L1 error against the exact solution; the three directions bitwise identical |
| V2 | `vortex/`: isentropic vortex, 64/128/256 cells | observed order of accuracy |
| V3 | `conservation/`: periodic AMR box with reflux | volume integrals constant to round-off; drift without reflux (negative control) |
| V6 | `shock-cylinder/`: Mach 2 shock over a cylinder, IB + solid AMR | refined surface blocks, mirror symmetry, positivity |
| V7 | `io/`: restart round trip, slices, auxiliary fields | bitwise restart; slice and auxiliary values exact |
| — | `multirealm/`: sod-x split in two realms at the diaphragm, mirror seam, beta cadence; the same with one realm refined (2:1 face crossed by the shock); the split along z (every seam row across ranks); sod-amr split at its 2:1 face into a coarse and a fine realm (`coupling = refined`) | union bitwise equal to the single-realm run (issues #37, #40, #52) |

The MHD tests (issue #41) live in `verification/mhd/`:

| Test | Case | Pass criterion |
|------|------|----------------|
| — | `plumbing/`: MHD state end to end | zero residual on a uniform state, restart bitwise, auxiliary fields exact, invalid inputs refused |
| MV-1 | the Euler suite with `physical_model = euler` | every Euler golden byte-identical |
| MV-2 | `zero-field/`: B = 0 through the MHD path (Sod x/y/z); uniform MHD state across a 2:1 patch | Sod L1 within the Euler bound, B and $\psi$ exactly 0; uniform state kept to round-off |
| MV-3 | `glm-pulse/`: a pulse in $B_x$, the $(B_x, \psi)$ telegraph system | L1 against the exact (damped) telegraph solution, observed order; damping, $c_h$ in dt and the wall parity of $\psi$ |
| MV-4 | `rj2a/`: Ryu-Jones 2a along x, y, z | L1 against the exact 7-state solution; x, y, z bitwise in the rotated frame; converges under refinement |
| — | `riemann/`: Brio-Wu, RJ4d; a negative-pressure state | Brio-Wu and RJ4d reach the end with zero floored cells; the negative-pressure state is floored when the floors are armed and stops the run when they are not |
| MV-5 | `linear-wave/`: the 7 linear eigenmodes | design order (~5) |
| MV-6 | `cpaw/`: circularly polarised Alfvén wave, 2-D | design order; left and right polarisations agree bitwise |
| MV-7 | `vortex/`: magnetised vortex | L1 order ~5 |
| MV-8 | `rotated-shock-tube/`: RJ1a at 63.4 and 45 degrees | $B_\parallel$ error against Tóth 2000 |
| MV-9 | `field-loop/`: field loop, $w = 1$; `--amr` across 2:1 seams | $\langle|B_z|\rangle$ at the GLM level and decreasing; magnetic energy decay bounded; the seam div(B) decays |
| MV-10 | `divb-peak/`: Dedner peak | GLM $\|\nabla \cdot \mathbf{B}\|_1$ decays; derived fields exact; monitor warns and stops |
| MV-11 | `conservation/`: V3 box with B and a perturbation | 9 integrals constant to 1e-13 with reflux; negative control without |
| MV-12 | `orszag-tang/`: Orszag-Tang; `--amr` with a symmetric 2:1 box | 180-degree symmetry, conservation, positivity; seam div(B) below the uniform peak |
| MV-13 | `rotor/`: MHD rotor | symmetry under rotation with $\mathbf{B} \to -\mathbf{B}$ |
| MV-14 | `multirealm/`: RJ2a in two realms; restart round trips; RJ2a on the sod-amr cells split at the 2:1 face | union bitwise; restart bitwise with a non-zero $\psi$ |

The M3 tests (issue #47) run the same scripts on the new numerics (`--numerics SOLVER[:RECON[:CORRECTION[:SENSOR]]]` for
the MHD ones, `--divergence-control eglm` for EGLM) plus their own:

| Test | Case | Pass criterion |
|------|------|----------------|
| RV-0 | `unit/`: Euler and MHD Riemann solvers on random states | consistency, exact contact (HLLC, HLLD), positive first-order updates, EGLM = GLM at $\psi = 0$; device = host |
| RV-2..RV-4 | `riemann-flux/`: vortex order, Sod/Lax/Shu–Osher, AMR conservation on `weno-riemann` | order at least 4.8; $L_1(\rho)$ within 1.05 × the splitting scheme's; integrals constant with reflux |
| RV-5..RV-7 | the MHD scripts with `--numerics hlld` | order, symmetry, conservation and recorded bounds as for the splitting scheme |
| EV-2..EV-4 | `glm-pulse/`, `rj2a/`, `orszag-tang/` with EGLM | telegraph solution; EGLM = GLM on RJ2a; $\int\rho$ exact and every other drift within its source bound |
| PV-0 | `unit/test_flume_positivity`: the limiter on random states with perturbed fluxes | every limited update positive and above the relative floor; non-finite fluxes replaced |
| PV-1..PV-4 | `mhd/positivity/`: Balsara–Spicer (2-D, 3-D) and Wu–Shu blasts; LeBlanc, double rarefaction, planar Sedov; smooth cases | final time reached with no inadmissible backbone (no floors); smooth runs bitwise unchanged |

`src/tests/flume/regression/` is the goldened regression suite (a copy of the PRISM harness): `run.sh cpu` runs in CI,
`run-fnl-local.sh` on a GPU workstation. `run-omp-bitwise.sh` (also in CI) requires the OpenMP CPU build to reproduce
the serial one bit for bit on the immersed-boundary case: one unexplained single-ulp divergence is on record (issue #35),
never reproduced since.

Known limitations:

- **Seam div(B) floor.** Cell-centred B cannot be kept divergence-free across a 2:1 coarse-fine face (Tóth & Roe
  2002): the seams inject a truncation-level div(B) source. GLM keeps it bounded (MV-9: it decays after the initial
  transient), unlike PRISM without damping (issue #29); `seam_max_divb` in the div(B) history monitors it.
  Constrained transport, which would remove it, is deferred.
- **Positivity.** Under GLM the density and pressure floors are a heuristic: the mixed-GLM energy coupling ($\psi$
  changes $B_n$ but is not in the energy) makes the first-order update of strongly magnetised states inadmissible, so
  the positivity limiter is refused with GLM. Use EGLM: with the limiter the splitting scheme reaches the end of the
  Balsara–Spicer ($\beta = 2.5 \cdot 10^{-4}$) and Wu–Shu ($\beta = 2.51 \cdot 10^{-6}$) blasts, and `weno-riemann`
  with HLLD runs the first even without it. The limiter guarantees positivity, not accuracy: the limited splitting
  scheme is noisy on those blasts, where HLLD stays clean. Across 2:1 AMR seams the limiter synchronises the seam
  flux at every stage (the coarse flux is the mean of the fine ones, both sides blended with one factor), so the
  reflux corrects round-off only, and inadmissible seam ghosts are blended toward the interior: the blast across
  seams reaches its end on both schemes with the mass constant to round-off
  ([#50](https://github.com/szaghi/adam/issues/50), PV-5).
- **GLM damping and reflux.** With damping, $\int \psi$ is not conserved across 2:1 faces (O(k dt) of the uncorrected
  leak, MV-11); the 8 physical integrals are.
- The FNL backend copies the coarse-fine seam faces to the host at every stage, which dominates its run time on AMR
  cases.
- **Inter-realm seams** ([#40](https://github.com/szaghi/adam/issues/40),
  [#52](https://github.com/szaghi/adam/issues/52)) work across ranks, whatever the partition of each realm. A `mirror`
  seam joins cells of the same size, a `refined` one a 2:1 jump: the fine realm is one level finer on every axis, both
  realms have the same block cells across the seam and the same `[amr] seam_ghost_fill`, and each coarse seam block
  faces 2x2 fine ones. On a `mirror` seam the blocks of the two realms need not line up (different block sizes along the seam,
  [#51](https://github.com/szaghi/adam/issues/51)); the two seam faces must cover each other, the forest stops at
  initialization otherwise. Seams work on octrees and quadtrees, and with single-block realms
  ([#54](https://github.com/szaghi/adam/issues/54)). All realms advance
  with one time step (no subcycling), and the positivity limiter is refused on multi-realm runs.

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
