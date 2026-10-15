# Ghost probe: every face and edge ghost on a linear field

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GP | Infrastructure | M4 ([#65](https://github.com/szaghi/adam/issues/65) P0, P1), M6 ([#68](https://github.com/szaghi/adam/issues/68) C2, D1) | CPU, FNL | yes | — |

## Case and scope

The directional WENO stencils read only the face ghosts: the cells outside a block along one axis. The dissipative
terms of M4 take cross derivatives, such as the tangential derivatives of the viscous stress or of the current at a
face. Those read the edge ghosts too: the cells outside a block along two axes, in the rows next to a block edge.
Before M4 nothing checked what an edge ghost held. GP checks every one, on every path that fills it: the same-level copy
and the MPI exchange, the 2:1 restriction and tricubic fill, the inter-realm seams (1:1 mirror, misaligned, 2:1
refined), and the physical boundary conditions.

**Oracle.** Each case holds a linear field, $\mathbf q = \mathbf q_1 (1 + \mathbf g \cdot \mathbf x)$, with every
velocity component non-zero so that every wall sign flip shows. The copy, the restriction and the tricubic fill
reproduce such a field exactly. Each case takes one step at `CFL = 1e-30`, which moves the field by $\sim10^{-30}$
relative; the post-step write then fills the ghosts in the order of a Runge–Kutta stage (seams, then the intra-realm
exchange and the boundary conditions) and writes them. `ghost_probe.py` compares every ghost with the value a stencil
must read:
- the field at the ghost's own centre, inside the realm or across a seam;
- along an axis where the ghost lies beyond a physical face, that face's condition applied to the value of its donor:
  the coordinate mirrored with the normal momentum negated for an inviscid wall (and the normal field on MHD); mirrored
  with the velocity reflected about the wall velocity for a no-slip wall (and the temperature $T_w^2/T$ for an
  isothermal one); the first interior cell for extrapolation; the inflow state for inflow; the wrapped coordinate for
  periodic. Beyond two faces the conditions compose in the backends' order (inflow first, else the first face in axis
  order acts on a donor valued the same way).

The no-slip walls are not linear in the state; the probe evaluates them with an independent implementation of the rule
of [Boundary conditions](../../boundary-conditions#wall-noslip-and-wall-isothermal-no-slip-walls), so `walls3d` checks
the plumbing (which donor, which face, which wall parameters, in which order), not the physics of the wall. The physics
comes with the viscous fluxes ([VV-3, Couette](../dissipative/couette)).

**The poison.** The negligible step leaves a stale ghost equal to a fresh one, so the oracle alone cannot see a ghost
read before its donor is written. Every case therefore runs with `[diagnostics] ghost_poison = .true.`: each ghost is set
to NaN before the fill, and a ghost that reads an unfilled donor is written as NaN, which the probe reports as an
infinite error. Its negative control is the defect it found (see Findings).

**The polynomial probe** ([#68](https://github.com/szaghi/adam/issues/68) D1). `check.sh --degree N` starts from the
field $\mathbf q_1 (1 + \sum_d g_d x_d^N)$ instead. The copies and the boundary conditions still reproduce it, but each
seam operator only up to its degree (the mean restriction 1, `tricubic` 3, `lagrange6` and `quintic` 5), so at $N > 1$
the probe checks the high-order operators in situ: their footprints, anchor shifts and every exchange path.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Models | Euler; ideal MHD with GLM (`box3d`) and with EGLM (`walls3d`), $\gamma = 1.4$ | [models](../../models#ideal-mhd) |
| Step | one step at `CFL = 1e-30` (WENO5 splitting, SSP-33; the update is negligible) | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Poison | every ghost NaN before the fill of a field write | `poison_ghosts`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `poison_ghosts_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |
| Intra-realm exchange | same-level copy, MPI exchange, 2:1 restriction and fill (`tricubic` default; `lagrange6`, `quintic` options) | `update_ghost`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `update_ghost_local`, `update_ghost_mpi`, `interp_seam_ghost`, `restrict_seam_high_order`, [`adam_field_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_field_object.F90); [boundary conditions](../../boundary-conditions#_2-1-amr-coarse-fine-ghosts) |
| Inter-realm seams | 1:1 mirror, misaligned, 2:1 refined | `fill_seam_from_peer_forest`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `seam_fill`, [`adam_seam_exchange.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_seam_exchange.F90); [boundary conditions](../../boundary-conditions#inter-realm-seam) |
| Physical boundaries | inflow, extrapolation, periodic, inviscid, no-slip and isothermal walls; edges compose the conditions of their faces | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `realm_edge_face`, `realm_edge_donor`, `wall_noslip_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); `set_boundary_conditions_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90); [boundary conditions](../../boundary-conditions#edges-and-corners-fec-6) |
| Field write | post-step save, ghosts refilled in Runge–Kutta stage order (seams first) | `post_step_forest`, `save_simulation_data`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |

## Initial and boundary conditions

Common to every case: [`linear`](../../initial-conditions#linear-a-linear-field-verification) initial condition with
$\mathbf q_1 = (\rho, u, v, w, p) = (1, 0.3, -0.2, 0.1, 1)$ and $\mathbf B = (0.8, 0.5, -0.3)$ on MHD,
$\mathbf g = (0.05, 0.11, 0.07)$ ($g_z = 0$ on a null $z$ axis); `ngc = 3`; `it_max = 0`, `time_max = 0`,
`CFL = 1e-30`; `ghost_poison = .true.`; MHD cleaning `glm_ch = 3`, `glm_alpha = 0.18`.

| Case (`make_probe.py`) | Grid | Boundaries | Paths |
|---|---|---|---|
| `box3d` (Euler and MHD) | $[0,1]^3$, octree, $2 \times 2 \times 2$ blocks of $8^3$, the block at the origin refined | $x_{\min}$ inflow ($\rho, u, v, w, p = 1.2, 0.4, 0.1, -0.1, 1.1$), $y_{\min}$, $y_{\max}$, $z_{\min}$ inviscid walls, $x_{\max}$, $z_{\max}$ extrapolation | 2:1 seams meeting an inflow face and two walls in 3-D |
| `channel2d` | $[0,2] \times [0,1]$, quadtree (`nk = 1`), blocks of $8^2$, two uniform levels, the box $[0, 0.5] \times [0, 0.25]$ refined | $x$ periodic, $y_{\min}$ inviscid wall, $y_{\max}$ extrapolation | the 2:1 seam crossing the periodic boundary at the wall |
| `mirror3d` | two realms $[0,1]^3$ and $[1,2] \times [0,1]^2$, each $2 \times 2 \times 2$ blocks of $8^3$ | 1:1 mirror seam at $x = 1$; $x_{\min}$ inflow, $y$ and $z_{\min}$ walls, the rest extrapolation | seam edges in every plane |
| `refined3d` | as `mirror3d`, the second realm one level finer | a `coupling = refined` (2:1) seam | the 2:1 inter-realm seam |
| `walls3d` (Euler and MHD with EGLM) | as `box3d` | $y_{\min}$ moving `wall-noslip` (`wall_u = 0.2`, `wall_w = -0.1`), $y_{\max}$ resting `wall-noslip`, $x_{\max}$ `wall-isothermal` (`wall_v = 0.15`, `wall_temperature = 0.9`), $z_{\min}$ inviscid wall, $x_{\min}$ inflow, $z_{\max}$ extrapolation | the no-slip walls of #65 P1 around a refined corner |
| `step` | the [Woodward–Colella forest](../euler/forward-facing-step) at $N = 80$ with its boxes, linearised ($v = 0.5$, $w = 0.3$) | inflow, outflow, walls | misaligned mirror seams, a re-entrant corner |

| Run | Ranks |
|---|---|
| every case, CPU | 1, 2, 3 |
| every case, FNL | 1, 2 |

Faces and edges must hold the expected value to $10^{-12}$ of the field scale. Corners are reported, not asserted, since
no M4 stencil reads them, and neither are the ghosts that lie outside every realm (`solid`, inside the step).

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/ghosts/make_probe.py` | writes the cases; `--linearize` rewrites the step forest into a probe run |
| `src/tests/flume/verification/step/make_step.py` | the step forest (`--cells 80 --refine`) |
| `src/tests/flume/verification/ghosts/check.sh` | runs the cases (`--np "1 2 3"`, `--cases`, `--keep`, `--seam-restriction R`, `--seam-ghost-fill F`, `--degree N`) |
| `src/tests/flume/verification/ghosts/ghost_probe.py` | the expected value of every ghost, the report by class (face, edge, corner) and by what lies beyond, the bound (`--tol 1e-12`) |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_ghosts`) | the figure, from `ghosts/check.sh --cases step --np 2 --keep` and the before-P0 run kept as `ghosts/work-before-p0-np2-step` (built from `3f51451d`) |

```bash
cd src/tests/flume/verification/ghosts && ./check.sh                    # CPU, np 1 2 3, every case
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --np "1 2"  # FNL
./check.sh --cases "box3d channel2d refined3d" --degree 5 --seam-restriction lagrange6 --seam-ghost-fill quintic
```

## Results

![Ghost probe on the step forest, before and after #65 P0](/flume/ghosts.png)

*Each point is a face or edge ghost of the step forest ($N = 80$, refined boxes, 2 ranks), coloured by its error. Top:
before P0, the realm edges at the step corner, the domain corners and the inflow corners carry $O(10^{-1})$ errors;
every other ghost is at round-off. Bottom: after.*

| Probe | Bound | CPU (np 1–3) | FNL (np 1–2) |
|---|---|---|---|
| linear field, every case, every face and edge group | $10^{-12}$ | $\le 1.3\cdot10^{-15}$ | $\le 1.3\cdot10^{-15}$ |
| linear field, `--seam-restriction lagrange6` (`box3d`, `channel2d`, `refined3d`) | $10^{-12}$ | $\le 1.1\cdot10^{-15}$ | $\le 1.1\cdot10^{-15}$ |
| degree 5, `lagrange6` + `quintic` (`box3d`, `channel2d`, `refined3d`) | $10^{-12}$ | $\le 5.8\cdot10^{-15}$ | $\le 5.8\cdot10^{-15}$ |
| degree 5, `tricubic` (negative control) | must fail | misses by $\sim10^{-5}$, seam groups only | |
| degree 2, mean restriction (negative control) | must fail | misses by $\sim4\cdot10^{-5}$, seam groups only | |
| the step forest before P0 | — | realm edges off by up to 0.78 | |

Related cases: [Forward-facing step](../euler/forward-facing-step) and [Multi-realm runs](../euler/multi-realm) run the
seams this probe inspects; [D1, quintic fill](../seams-studies/quintic-fill) uses the polynomial probe.

## Findings

- **What it found** (P0). Before P0, a ghost beyond two realm faces held a copy of its inward diagonal cell. That covers
  wall plus wall at the step corner, wall plus seam where the step face meets the A–B seam, inflow plus wall, and the
  like. On the step forest those ghosts were off by up to 0.78 of the field scale: the wrong position, and the wrong
  sign of the normal momentum. They now compose the conditions of the faces, as
  [Boundary conditions](../../boundary-conditions#edges-and-corners-fec-6) describes. The face ghosts, and the edges
  with another block beyond one of their faces, were already right.
- **What the poison found** (P1, after P0 had been committed). An extrapolation ghost copied the previous ghost of its
  chain along the normal. When the ghost also lies beyond a block interface along another axis, deeper than along the
  normal, that previous ghost is in the same crown, so in the same device launch: on FNL the rows raced. Two identical
  FNL runs of `sod-x` wrote different edge and corner ghosts (the interior identical), and with the poison GP failed
  the `block+extrapolation` edges on every case (NaN read). The CPU walks each row outward and read its donors fresh,
  by the row order alone. The P0 probe could not see it: without the poison the stale value equalled the fresh one. An
  extrapolation ghost now copies the first interior cell along its normal directly (the same value, and a donor the
  exchange has filled), on both backends; with the poison every case passes and two FNL runs are bit for bit equal.
- **A gap in its own instrument.** The step-0 file is written while each realm initialises, before the forest connects
  the realms, so its seam ghosts are unfilled. GP therefore reads the file written after the negligible step, whose
  write refills the seams first.
- The negative controls of the polynomial probe are not part of `check.sh`; they are run by hand with
  `--seam-ghost-fill tricubic --degree 5` and `--seam-restriction mean --degree 2`.
