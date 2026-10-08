# Boundary conditions

Each of the six faces of the domain takes a boundary condition from its section `[bc_x_min]`, `[bc_x_max]`,
`[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`. The condition fills the ghost cells of the blocks that touch the
face before every residual evaluation; interior block faces, inter-realm seams and AMR coarse–fine faces are filled by
the exchange machinery instead (below).


## Input

There are six sections, `[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]` and `[bc_z_max]`, which are
faces 1–6. Each section **requires** `type`. The value goes through
`strip_control`, so CRLF files are accepted. Six values are accepted, and matching is case-sensitive:

| `type` | id | per-face keys |
|---|---|---|
| `extrapolation` | `BC_EXTRAPOLATION = 1` | none |
| `inflow` | `BC_INFLOW = 2` | `r, u, v, w, p` (+ `bx, by, bz` for MHD), all required |
| `wall-inviscid` | `BC_WALL_INVISCID = 3` | none |
| `wall-noslip` | `BC_WALL_NOSLIP = 4` | `wall_u, wall_v, wall_w` (optional, default 0; the normal one must be 0) |
| `wall-isothermal` | `BC_WALL_ISOTHERMAL = 5` | `wall_temperature` (required, $> 0$), `wall_u, wall_v, wall_w` as above |
| `periodic` | library `BC_PERIODIC = -1` | none |

*Validation:*
- an unknown value stops the run and the message lists the six spellings;
- a missing inflow key, or a missing or non-positive `wall_temperature`, is fatal;
- a non-zero normal wall velocity is fatal (the wall does not move through the domain);
- periodicity must be **paired**: `[bc_x_min]` and `[bc_x_max]` must both be periodic or neither, and the same for y and z;
- a physical model other than Euler or MHD is fatal.

## How the ghost cells are filled

`update_ghost` runs three steps in order:
1. the intra-realm local copies;
2. the MPI exchange;
3. `set_boundary_conditions`.

Step 3 walks the boundary crown map `maps%local_map_bc_crown` **crown by crown**, $c = 1..n_{gc}$, from the innermost
crown outwards. Each ghost cell $(i,j,k)$ carries an "inward step" $(\delta_i,\delta_j,\delta_k)$, its BC type and its
`fec` code (face 1–6, edge 7–18, corner 19–26). The CPU implementation is
`cpu/adam_flume_cpu_object.F90`; the FNL twin is `fnl/adam_flume_fnl_kernels.F90`, called crown by
crown from `fnl/adam_flume_fnl_object.F90`. Both apply the same rules.

The grid receives the six BC ids through `grid%set_bc_type`. That call
turns on `is_ijk_periodic` for a periodic axis.

### `extrapolation`: zeroth order (constant)

$$
\mathbf q_{\text{ghost}}(i,j,k) = \mathbf q(i-\delta_i,\,j-\delta_j,\,k-\delta_k)
$$

Each ghost copies its inward neighbour. The crowns are filled innermost first, so every ghost column ends up equal to
the first interior cell: $\mathbf q(1-g) = \mathbf q(1)$ for $g = 1..n_{gc}$. This is order-0 extrapolation, a
zero-gradient outflow, and matches the `!< Zeroth-order extrapolation` comment. All
$n_v$ variables are copied, including $\psi$.

### `inflow`: prescribed state (Dirichlet on the full conservative vector)

$$
\mathbf q_{\text{ghost}} = \mathbf q_{\text{inflow}}^{(f)}\quad\text{on every crown}
$$

$\mathbf q_{\text{inflow}}^{(f)}$ is converted once from the primitive keys. For MHD
with GLM this sets **$\psi = 0$ in the inflow ghosts**.

### `wall-inviscid`: slip wall by mirroring with a sign vector

$$
\mathbf q_{\text{ghost}}(i,j,k) = \mathbf S_d \odot \mathbf q(i_m,j_m,k_m),\qquad d = \lceil f/2\rceil
$$

The mirror index reflects the ghost across the face:

| face | mirror |
|---|---|
| $x_{\min}$ | $i_m = 1-i$ |
| $x_{\max}$ | $i_m = 2n_i+1-i$ |
| $y_{\min}$ | $j_m = 1-j$ |
| $y_{\max}$ | $j_m = 2n_j+1-j$ |
| $z_{\min}$ | $k_m = 1-k$ |
| $z_{\max}$ | $k_m = 2n_k+1-k$ |

So ghost $1-g$ mirrors interior cell $g$, for example $\mathbf q(0) \leftarrow \mathbf S\,\mathbf q(1)$ and
$\mathbf q(-1) \leftarrow \mathbf S\,\mathbf q(2)$.

The sign vector $\mathbf S_d$ is built by `adam_flume_bc_object.F90`:

| variable | Euler | MHD / MHD-GLM |
|---|---|---|
| $\rho$, $E$ | $+1$ (even) | $+1$ |
| $\rho u_n$ (the normal momentum of direction $d$) | $-1$ (odd) | $-1$ |
| $\rho u_{t}$ (tangential momenta) | $+1$ | $+1$ |
| $B_n$ (the normal field of direction $d$) | — | $-1$ (odd), a perfectly conducting reflecting wall (issue #41, D-12) |
| $B_t$ | — | $+1$ |
| $\psi$ | — | $+1$ (**even**) |

$E$ is unchanged by the sign flips because it depends only on squares of the flipped components. The ghost state is
therefore the exact mirror image. As a result, $u_n = 0$ and $B_n = 0$ at the wall face in the sense of the symmetric
average, and $\partial_n\psi = 0$ there.

### `wall-noslip` and `wall-isothermal`: no-slip walls

The two no-slip walls ([#65](https://github.com/szaghi/adam/issues/65), D-M4-7) mirror the cell like `wall-inviscid`
and reflect the **whole** velocity about the wall velocity $\mathbf u_w$ (tangential, `wall_u, wall_v, wall_w`):

$$
\mathbf u_g = 2\,\mathbf u_w - \mathbf u_m, \qquad p_g = p_m,
$$

so the mean of the ghost and its mirror cell is the wall velocity (no slip, no penetration). The pressure is mirrored.

- `wall-noslip` (adiabatic) mirrors the temperature, so $\rho_g = \rho_m$ and $\partial_n T = 0$ at the wall in the
  sense of the symmetric difference.
- `wall-isothermal` sets the ghost temperature to the geometric mirror, $T_g = T_w^2 / T_m$, so the geometric mean
  of the ghost and its mirror cell is the wall temperature, and the density from the mirrored pressure,
  $\rho_g = p_m / (R\,T_g)$, with $T = p/(\rho R)$. The ghost is admissible for any admissible cell. Until #65 P2 the
  mirror was linear, $T_g = 2\,T_w - T_m$, which is negative beside gas hotter than twice the wall: on Sod's right
  state ($T = 0.0027$) against a wall at $T_w = 0.001$ the ghosts reached $\rho = -1.96$, silently, with no non-finite
  value. For $T_m = T_w + \delta$ the two mirrors differ by $\delta^2/T_w$, so both are second-order at the wall.

The total energy is rebuilt from the ghost state, $E_g = p_g/(\gamma-1) + \tfrac12\rho_g|\mathbf u_g|^2 + e_{mag}$, where
the magnetic energy $e_{mag}$ (and $\psi^2/2$ under EGLM) is unchanged by the mirror. The field follows the
`wall-inviscid` rule (normal component odd, a perfectly conducting wall); $\psi$ is even. On a resting adiabatic wall the
rule reduces to a sign vector with every momentum component odd. Both walls are second-order accurate at the wall, by
construction of the mirror.

The rule is the pure procedure `wall_noslip_ghost` in `adam_flume_bc_object.F90`, shared by both backends. It is
model-agnostic: the field components are those present in the state vector (none for Euler), so no kernel branches on
the model. The ghost probe GP checks it on every face and edge, moving and resting walls, Euler and MHD with EGLM
([verification](./verification#ghost-cells)).

The walls are boundary conditions of the ideal solver too: they take effect on any run. Their physical test is the
compressible Couette flow between an isothermal wall and a moving adiabatic one, against its exact profiles
([VV-3](./verification#vv-3-compressible-couette-flow)): the moving-wall mirror is exact for its linear velocity and
symmetric temperature, the isothermal mirror $T_w^2/T$ is second-order, so the run converges at order 2.

### `periodic`

`periodic` gets **no crown rows**. The library tree builds true periodic neighbours, and `update_ghost_local` and
`update_ghost_mpi` fill the ghosts across blocks and ranks. Unlike PRISM's `periodic`, this is the library `BC_PERIODIC`, so it works
with any number of blocks along the axis. The only constraint is the pairing rule in B.1.

### Edges and corners (`fec > 6`)

An edge ghost lies outside its block along two axes, a corner ghost along three. The directional WENO stencils never
read them, but a cross derivative does: the tangential derivatives of the viscous stress or of the current at a face
(M4, [#65](https://github.com/szaghi/adam/issues/65)) read the edge ghosts of the rows next to a block edge. They are
filled so that the ghost layer continues the face ghosts, and the probe GP of the [verification](./verification#ghost-cells)
holds every face and edge ghost to a linear field on every path.

Each crown row carries the boundary its cell lies beyond (`bc_fec`, from the tree), which gives two cases:

- **Beyond one realm face** (`bc_fec` a face). The cell's other directions lead into another block, whose ghosts the
  exchange has filled. The face's kind applies along its normal exactly as for a face ghost, so a wall edge mirrors
  the exchanged ghost and negates its normal momentum.
- **Beyond two or three realm faces** (`bc_fec` an edge or a corner: an edge or a corner of the realm). The kind of
  one physical face among them applies, inflow first (an inflow ghost holds the inflow state whatever else it lies
  beyond), else the first non-seam face, `realm_edge_face`. It acts on the donor `realm_edge_donor`: the cell mirrored
  about that face (any wall kind) or the first interior cell along its normal (extrapolation), with the indexes along the other
  axes kept. The donor therefore lies beyond the other faces only. It is a face ghost (for an edge) or an edge ghost
  (for a corner), already filled by the seam, the exchange or its own face. Two walls compose to a double mirror with
  both normal momenta negated, and a wall beside a seam mirrors the seam ghost.

Beyond seams only, the cell lies outside every realm: the re-entrant corner of an L-shaped forest, inside the step of
the [Woodward–Colella tunnel](./verification#forward-facing-step-a-three-realm-forest-on-quadtrees). It keeps a copy of
its inward diagonal, $\mathbf q(i-\delta_i, j-\delta_j, k-\delta_k)$: finite, and read by no stencil of a fluid cell.

Both backends fill the crown in three passes, the rows beyond one face, then two, then three, each crown by crown, so
that every donor is filled before it is read (`set_boundary_conditions` on the CPU, `set_boundary_conditions_dev` on
FNL). Within the first pass no row reads another row of the pass: a wall row reads the mirrored interior cell, an
extrapolation row the **first interior cell** along its normal (not the previous ghost of the chain, which gives the
same value but is a row of the same launch: on FNL the chain raced, and a ghost beyond a block interface along another
axis could read its donor before it was written; issue #65 P1). The kind of each realm face comes from `[bc_*]`, overridden by `maps%seam_face` for the faces the forest glued.

Before [#65](https://github.com/szaghi/adam/issues/65) P0 every row beyond two realm faces held its inward diagonal
copy. On a linear field that is off by up to 0.78 of the field scale at a wall corner (wrong sign of the normal
momentum, wrong position), measured by GP on the step forest.

### Inter-realm seam

`set_boundary_conditions` leaves seam face rows untouched.
The forest fills them through `fill_seam_from_peer_forest`, which is a **plain injection copy** of the peer's interior
cell into this realm's ghost cell, following the seam ghost map rows. The
buffers used are `q` when `stage_active == 0`, otherwise the active stage of `rk%q_rk`, on each side independently.
`update_ghost` does not fill seam ghosts. Any diagnostic that reads them has to refill them explicitly;
`compute_divb_history` does this. The field write does too: with sibling realms, `save_simulation_data` refills the
seams and then calls `update_ghost`, the order of a Runge-Kutta stage, so the written ghosts (edges included) are those
the stencils read. The step-0 file is the exception: it is written while each realm initialises, before the forest
connects the realms, and its seam ghosts are unfilled.

### 2:1 AMR coarse-fine ghosts

These ghosts are library behaviour, not FLUME-specific. The regime is chosen by the optional key
`[amr] seam_ghost_fill = injection | restriction-compatible | tricubic`. The default when the key is absent is
`tricubic`, and an unknown value is fatal. The fill runs inside
`field%update_ghost_local`.

### Immersed solids

Immersed solids (for example `regression/shock-cylinder-ib`) are not a `[bc_*]` kind. They enter through the level-set
`phi` and the cut spacing in `compute_flux_difference`, and are out of scope here.

## Examples

Walls on x, extrapolation elsewhere, from `src/tests/flume/verification/sod/sod-wall-x.ini`:
```ini
[bc_x_min]
type = wall-inviscid
[bc_x_max]
type = wall-inviscid
[bc_y_min]
type = extrapolation
[bc_y_max]
type = extrapolation
[bc_z_min]
type = extrapolation
[bc_z_max]
type = extrapolation
```
Euler inflow, from `src/tests/flume/verification/shock-cylinder/shock-cylinder.ini`:
```ini
[bc_x_min]
type = inflow
r    = 3.7333333333333334
u    = 1.25
v    = 0.0
w    = 0.0
p    = 4.5
[bc_x_max]
type = extrapolation
[bc_y_min]
type = wall-inviscid
[bc_y_max]
type = wall-inviscid
[bc_z_min]
type = extrapolation
[bc_z_max]
type = extrapolation
```
MHD inflow, from `src/tests/flume/verification/mhd/plumbing/mhd-uniform-glm.ini`:
```ini
[bc_z_min]
type = inflow
r    = 1.0
u    = 0.3
v    = -0.2
w    = 0.1
p    = 1.0
bx   = 0.8
by   = 0.5
bz   = -0.3
```
Periodic faces are written the same way, `type = periodic`, and must be set on both faces of an axis.
