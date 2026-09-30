# Boundary conditions

Each of the six faces of the domain takes a boundary condition from its section `[bc_x_min]`, `[bc_x_max]`,
`[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`. The condition fills the ghost cells of the blocks that touch the
face before every residual evaluation; interior block faces, inter-realm seams and AMR coarse–fine faces are filled by
the exchange machinery instead (below).


## Input

There are six sections, `[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]` and `[bc_z_max]`, which are
faces 1–6. Each section **requires** `type`. The value goes through
`strip_control`, so CRLF files are accepted. Four values are accepted, and matching is case-sensitive:

| `type` | id | per-face keys |
|---|---|---|
| `extrapolation` | `BC_EXTRAPOLATION = 1` | none |
| `inflow` | `BC_INFLOW = 2` | `r, u, v, w, p` (+ `bx, by, bz` for MHD), all required |
| `wall-inviscid` | `BC_WALL_INVISCID = 3` | none |
| `periodic` | library `BC_PERIODIC = -1` | none |

*Validation:*
- an unknown value stops the run and the message lists the four spellings;
- a missing inflow key is fatal;
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

### `periodic`

`periodic` gets **no crown rows**. The library tree builds true periodic neighbours, and `update_ghost_local` and
`update_ghost_mpi` fill the ghosts across blocks and ranks. Unlike PRISM's `periodic`, this is the library `BC_PERIODIC`, so it works
with any number of blocks along the axis. The only constraint is the pairing rule in B.1.

### Edges and corners (`fec > 6`)

For every BC type these ghosts are **extrapolated** from the diagonal inward neighbour,
$\mathbf q(i-\delta_i, j-\delta_j, k-\delta_k)$. This covers wall and inflow edges
too. The directional stencils never read these cells; they are filled only so that the auxiliary-variable pass sees
finite values.

### Inter-realm seam

`set_boundary_conditions` leaves seam face rows untouched.
The forest fills them through `fill_seam_from_peer_forest`, which is a **plain injection copy** of the peer's interior
cell into this realm's ghost cell, following the seam ghost map rows. The
buffers used are `q` when `stage_active == 0`, otherwise the active stage of `rk%q_rk`, on each side independently.
`update_ghost` does not fill seam ghosts. Any diagnostic that reads them has to refill them explicitly;
`compute_divb_history` does this.

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
