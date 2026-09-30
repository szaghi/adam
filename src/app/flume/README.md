# FLUME

**FLUME** — **F**luid **L**orentz-coupled **U**nsteady **M**agnetohydrodynamic **E**quations — is a compressible magnetohydrodynamics (MHD) application built on the ADAM framework. It solves the equations of an electrically conducting, compressible fluid interacting with its own magnetic field through the Lorentz force.

> **Status: development.** Milestone M1 of [issue #35](https://github.com/szaghi/adam/issues/35) delivered the
> **compressible Euler** equations on both backends; milestone M2 of [issue #41](https://github.com/szaghi/adam/issues/41)
> adds **ideal MHD** with mixed GLM divergence cleaning, verified on both backends including AMR and multi-realm runs.
> FLUME supersedes CHASE, which is deprecated.

## Physical Models

| Model | `[physics] physical_model` | `[mhd] divergence_control` | State |
|-------|----------------------------|----------------------------|-------|
| Compressible Euler | `euler` | — | $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E)^\top$, `nv = 5` |
| Ideal MHD, GLM cleaning | `mhd-ideal` | `glm` | $\mathbf{q} = (\rho, \rho u, \rho v, \rho w, E, B_x, B_y, B_z, \psi)^\top$, `nv = 9` |
| Ideal MHD, no control | `mhd-ideal` | `none` | the first 8 of the above, `nv = 8`: diagnostic use only (in multi-D the divergence error grows unchecked: the magnetised vortex reaches a negative pressure) |

All models are inviscid with an ideal gas (`[physics] cp, cv`, J/(kg K)); ideal MHD is a perfectly conducting single
fluid (ideal Ohm's law $\mathbf{E} + \mathbf{u} \times \mathbf{B} = 0$). `mhd-ideal` with immersed solids is refused.
Dissipative effects (viscosity, thermal conduction, resistivity) are outside the current scope.

## Implemented Capabilities

| Area | What is available | INI |
|------|-------------------|-----|
| Space | WENO flux splitting, per-face Roe eigenvectors, per-wave local Lax-Friedrichs; reconstruction in characteristic or conservative variables; orders `weno-u-3` to `weno-u-9` (the centred `weno-c-*` schemes are refused). MHD: block-diagonal characteristic decomposition, 7×7 Roe–Balsara-normalised eigenvectors (Stone et al. 2008) at the arithmetic mean of the primitive states with $B_n$ as a parameter, plus the $(B_n, \psi)$ block with speeds $\mp c_h$ | `[numerics] scheme_space = weno`, `reconstruction_variables`; `[weno] scheme` |
| Time | Library Runge-Kutta schemes (SSP and low-storage), CFL time step | `[runge_kutta] scheme`; `[time] CFL, it_max, time_max` |
| Boundary conditions | `extrapolation`, `inflow` (primitive state `r, u, v, w, p`; MHD adds `bx, by, bz`, with $\psi = 0$), `wall-inviscid` (MHD: perfectly conducting wall, $u_n$ and $B_n$ odd, the rest and $\psi$ even), `periodic` (both faces of an axis or neither) | `[bc_{x,y,z}_{min,max}] type` |
| Initial conditions | `uniform` (optionally with a seeded perturbation), `isentropic-vortex`, `riemann-problem` (piecewise-constant regions; MHD regions add `bx, by, bz`); MHD only: `glm-pulse`, `divb-peak`, `mhd-linear-wave`, `mhd-cpaw`, `mhd-vortex`, `orszag-tang`, `mhd-rotor`, `field-loop`, `rotated-riemann` | `[initial_conditions] type` |
| AMR | Init-time refinement (`amr_iterations` passes) by geometric box, variable gradient or immersed-solid surface; 2:1 coarse-fine faces with stage-weighted conservative reflux (B and $\psi$ included). With markers the tree must be an octree (`ratio = 8`; with a null axis, `nk >= 4`): a quadtree with markers is refused ([#46](https://github.com/szaghi/adam/issues/46)) | `[amr]`, `[initial_conditions] amr_iterations`, `[numerics] reflux` |
| Immersed boundary | Euler only. Static solids, inviscid wall: distance function, eikonal extrapolation into the solid, cut-cell spacing, solid masks in the Runge-Kutta stages |
| Multi-realm | A forest manifest glues realms through inter-realm seams (mirror coupling, per-seam cadence); every model | `[forest]` manifest | `[solids]` |
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
transport is deferred. The `[mhd]` section (read only with `physical_model = mhd-ideal`):

| Key | Values | Meaning |
|-----|--------|---------|
| `divergence_control` | `glm`, `none` | model variant (state width 9 or 8) |
| `glm_ch` | > 0 | cleaning speed $c_h$, constant and uniform, part of the time-step bound |
| `glm_alpha` | >= 0 | damping $c_h^2/c_p^2 = \alpha\, c_h / L$ |
| `glm_damping_length` | > 0 or `min-cell` | $L$; `min-cell` is the minimum cell spacing of the realm (Mignone & Tzeferacos 2010) |
| `glm_ch_check` | `warning`, `error` | action when the fastest wave outruns $c_h$ |
| `divb_tol`, `divb_error` | >= 0, logical | div(B) monitor: warn (or stop) when `max_divb` exceeds the tolerance; 0 disables it |
| `rho_floor`, `p_floor` | >= 0 | positivity floors with global counters in the log; 0 disables them (a non-positive state then stops the run) |

The GLM keys are required with `glm` only; every value is range-checked and an unknown spelling stops the run.

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
| — | `multirealm/`: sod-x split in two realms at the diaphragm, mirror seam, beta cadence; the same with one realm refined (2:1 face crossed by the shock) | union bitwise equal to the single-realm run (issue #37) |

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
| MV-14 | `multirealm/`: RJ2a in two realms; restart round trips | union bitwise; restart bitwise with a non-zero $\psi$ |

`src/tests/flume/regression/` is the goldened regression suite (a copy of the PRISM harness): `run.sh cpu` runs in CI,
`run-fnl-local.sh` on a GPU workstation. `run-omp-bitwise.sh` (also in CI) requires the OpenMP CPU build to reproduce
the serial one bit for bit on the immersed-boundary case: one unexplained single-ulp divergence is on record (issue #35),
never reproduced since.

Known limitations:

- **Seam div(B) floor.** Cell-centred B cannot be kept divergence-free across a 2:1 coarse-fine face (Tóth & Roe
  2002): the seams inject a truncation-level div(B) source. GLM keeps it bounded (MV-9: it decays after the initial
  transient), unlike PRISM without damping (issue #29); `seam_max_divb` in the div(B) history monitors it.
  Constrained transport, which would remove it, is deferred.
- **Positivity.** The density and pressure floors are a heuristic under GLM, not a positivity-preserving scheme. Very
  low plasma beta is out of reach: the Balsara–Spicer strong blast ($\beta = 2.5 \cdot 10^{-4}$) fails within a few
  steps. The cause is the mixed-GLM energy coupling ($\psi$ changes $B_n$ but is not in the energy, so the change of
  magnetic energy is taken from the thermal pressure): with it, the first-order Lax–Friedrichs update of the stage
  states is inadmissible and no flux limiter can help. A prototype with EGLM (Derigs et al. 2018) and the cell-based
  limiter reaches the end of the blast (`src/tests/flume/verification/mhd/positivity-probe/`); both are planned in
  [#47](https://github.com/szaghi/adam/issues/47).
- **GLM damping and reflux.** With damping, $\int \psi$ is not conserved across 2:1 faces (O(k dt) of the uncorrected
  leak, MV-11); the 8 physical integrals are.
- **Quadtree AMR** with markers is refused ([#46](https://github.com/szaghi/adam/issues/46)); use an octree.
- The FNL backend copies the coarse-fine seam faces to the host at every stage, which dominates its run time on AMR
  cases; inter-realm seams must join blocks of the same resolution (mirror coupling), and the realms are verified with
  the same block partition on both sides of the seam (the seam fluxes are matched rank-locally).

## License

FLUME is part of the ADAM framework, released under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html) (LGPLv3).

> Copyright (C) Andrea Di Mascio, Federico Negro, Giacomo Rossi, Francesco Salvadore, Stefano Zaghi.
