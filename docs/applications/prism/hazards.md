# Hazards

This page collects the PRISM inputs that are refused, and the ones that run but behave differently than their key
names suggest. Every other PRISM page links here instead of restating these rules. The keys themselves are in the
[input reference](./input).

Two stop points appear below. **Load** is the parsing of the INI file. **Initialization** is the backend set-up after
it, before the first step. A few checks fire later, at the first solve that needs them; the tables say so.

## Divergence correction

The GLM scalars φ (for $\mathbf D$) and ψ (for $\mathbf B$) exist only under `[numerics] divergence_correction =
hyperbolic`. `constrained_transport` selects which of them are carried. `physics` sets their slots and the variable
names from one branch ([#11](https://github.com/szaghi/adam/issues/11), [#61](https://github.com/szaghi/adam/issues/61)):

| `constrained_transport` | `divergence_correction` | Result |
|---|---|---|
| `no` | `no` | `nv = 9` |
| `D`, `B`, `DB` | `no` | constrained transport ignored with a warning, `nv = 9` |
| `D` or `B` | `hyperbolic` | `nv = 10` (φ or ψ in slot 7) |
| `DB` | `hyperbolic` | `nv = 11` (φ in 7, ψ in 8) |
| `no` | `hyperbolic` | **refused** at load: there is no scalar to carry |
| any | `poisson` | **refused** at load: no backend implements the projection |

PIC models add `rho`, one more slot, up to `NV_MAX = 12` ([state vector](./models#state-vector)). The cleaning
equations are on the [models page](./models#divergence-correction).

## Seam div(B)

On a uniform grid the centred curl and divergence commute, so $\nabla_h\cdot\mathbf B$ stays at round-off
([numerics](./numerics#the-discrete-curl-and-divergence)). At a 2:1 AMR face the two sides difference at different
resolutions and the identity fails. The seam injects a $\nabla\cdot\mathbf B$ source in the grid-oscillation modes.

- The source is truncation-order, $O(h^p)$: it converges under refinement.
- At fixed $h$ it grows without bound in time, independent of $\Delta t$. On `rmf-amr-fd-pulse` (source-free pulse,
  1 rank, `ni = 16`) the seam max $\lvert\nabla\cdot\mathbf B\rvert$ grows from 6.9 to 1002 in 100 steps; the
  uniform-grid control stays at 0.
- $\lVert\mathbf B\rVert$ stays bounded (−3 % over the same run): this is a constraint violation, not an energy
  instability.
- The same source appears with `fd_centered` and `fv_centered`, and on octrees and quadtrees.

Five local fixes were ruled out ([#29](https://github.com/szaghi/adam/issues/29)). The resolution is
**accept-truncation**. An upwinded or dissipative GLM coupling at the seams is open as
[#99](https://github.com/szaghi/adam/issues/99); a staggered-$\mathbf B$ route is parked as
[#30](https://github.com/szaghi/adam/issues/30). The measurements are on the
[forest page](/concepts/forest#the-2-1-seam-is-not-div-free-accept-truncation).

| Key (`[IO]`) | Default | Effect |
|---|---|---|
| `seam_divB_tol` | `-1` (off) | When `> 0`, warn when the global max $\lvert\nabla\cdot\mathbf B\rvert$ outside the absorbing layers exceeds it. The maximum is not seam-local; on a grid without 2:1 faces it stays at round-off and the monitor stays silent. |
| `seam_divB_error` | `.false.` | `.true.` stops the run instead of warning. |

**Mitigation.** `divergence_correction = hyperbolic` with `constrained_transport = B` or `DB` and a damping rate
`[physics] alpha > 0` delays the growth. It does not cure it.

**Inter-realm seams.** A 1:1 inter-realm seam is divergence-free only with `coupling_cadence = stage_coincident` (β);
under α the stages after the first read stale seam ghosts ([#31](https://github.com/szaghi/adam/issues/31)). The
cadences are on the [numerics page](./numerics#multi-realm-runs).

## PML constraints

| Backend | Accepted | Checked at |
|---|---|---|
| CPU | `PML_type` `CLASSIC`, `CLASSIC_DIRECT`, `BERMUDEZ`, `CFS`; `scheme_space = fd_centered`; `scheme_time = runge_kutta` with SSP-11/22/33/54; `divergence_correction = no` (hence no constrained transport) | initialization |
| FNL | the same `PML_type` kinds; `fd_centered`; `runge_kutta` with SSP-11/22/33/54; no `hyperbolic` correction, no constrained transport | initialization, before any device allocation |

Any other combination is fatal. The `[PML]` section itself is also checked at load:

- `width > 0` and the six face flags are required; every flag `.false.` is fatal.
- Each `PML_type` has its own required parameters (the [`[PML]` table](./input#pml)). A missing one is fatal, not
  defaulted. An unknown `PML_type` is fatal.
- **`CLASSIC_DIRECT` is not a Runge–Kutta PML.** It skips the auxiliary variables and their integration, and damps the
  tangential components directly. It is a testing instrument, not the production path
  ([absorbing layers](./pml#auxiliary-equations)).

## Elliptic boundary map

The t0 solves ([numerics](./numerics#t0-elliptic-solves)) run on FLAIL's FFT solver with PRISM's own operator
$\mathcal D(\mathcal G)$ of order `fdv_order`. Each face kind maps to an elliptic kind (`build_elliptic_bc_types`), then
to a FLAIL face kind (`solve_time_zero`):

| `[bc_*] type` | Elliptic kind | FLAIL face kind |
|---|---|---|
| `periodic` | periodic | `FLAIL_BC_PERIODIC` |
| `Dirichlet` | Dirichlet | `FLAIL_BC_ZERO_GHOSTS` (every ghost layer zero) |
| `PEC` | PEC | `FLAIL_BC_DIRICHLET` (φ, tangential $\mathbf A$), `FLAIL_BC_NEUMANN` (normal $\mathbf A$) |
| `Neumann`, `Silver_Muller` | exact/open | `FLAIL_BC_FREE_GHOSTS` (the free-space potential of the source in the ghosts, summed over the whole domain) |
| `extrapolation` | — | **fatal**, at the first t0 solve |

The map is total, but it runs **only when a solve first touches the faces**. An input with an `extrapolation` face
passes load and initialization, and stops at t0 if a PIC field or a t0 current needs a solve. A null axis takes the
periodic kind.

## Refused combinations

| Input | Where it stops |
|---|---|
| `[numerics] scheme_time = leapfrog`, `blanes_moan`, `commutator_free_magnus` | initialization, both backends ([#55](https://github.com/szaghi/adam/issues/55)) |
| `[PIC] scheme_time = leapfrog` | initialization, both backends |
| `[numerics] divergence_correction = poisson` | load |
| `divergence_correction = hyperbolic` without constrained transport | load |
| `[bc_*] type = radiative` | load ([#57](https://github.com/szaghi/adam/issues/57)) |
| `[initial_conditions] type = riemann-problem`, `magnetic_nozzle`, `rmf_magnetic_nozzle` | load ([#56](https://github.com/szaghi/adam/issues/56)) |
| `[external_fields] external_fields_applied = magnetic_nozzle`, `rmf_and_magnetic_nozzle` | initialization |
| `[numerics] scheme_space = weno` on FNL | initialization, after the device allocation |
| `runge-kutta-yoshida` on FNL | initialization |
| a non-SSP Runge–Kutta scheme with PIC or the PML | initialization |
| a non-SSP Runge–Kutta scheme on the staged path (2:1 faces, forests) | when the forest first queries the stages |
| `[PIC] current_weighting_model = conserving` on the CPU backend with more than one block or one rank | initialization, at the t0 conserving deposit (the first step on a restart) |
| `[PIC] current_weighting_model = conserving` on the staged path (2:1 faces, forests) | the first staged step, both backends |
| `[PIC] current_weighting_model = conserving` with periodic faces and `current_conserving_solver = londrillo` or `esirkepov_tail_cleanup = .true.` | initialization, both backends ([PIC](./pic)) |
| a refined grid that enters a t0 elliptic solve | the solve ([#84](https://github.com/szaghi/adam/issues/84) E7: the wide operator across 2:1 faces is not defined) |
| `fdv_order = 10` with a t0 elliptic solve | the solve (FLAIL accepts 2, 4, 6, 8) |
| an `extrapolation` face with a t0 elliptic solve | the solve ([above](#elliptic-boundary-map)) |

Every refusal prints a message naming the combination. `src/tests/prism/regression/refusals.sh` keeps them refused,
one leg per refusal ([cookbook](/cookbook/regression#t5-check-that-bad-inputs-stop-at-initialization)).
