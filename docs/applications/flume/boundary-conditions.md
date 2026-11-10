# Boundary conditions

Each of the six faces of the domain takes a boundary condition from its section `[bc_x_min]`, `[bc_x_max]`,
`[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]`, `[bc_z_max]`. The condition fills the ghost cells of the blocks that touch the
face before every residual evaluation; interior block faces, inter-realm seams and AMR coarse–fine faces are filled by
the exchange machinery instead (below).

## Input

There are six sections, `[bc_x_min]`, `[bc_x_max]`, `[bc_y_min]`, `[bc_y_max]`, `[bc_z_min]` and `[bc_z_max]`, which are
faces 1–6. Each section **requires** `type`. The value goes through
`strip_control`, so CRLF files are accepted. Eleven values are accepted, and matching is case-sensitive:

| `type` | id | Meaning | per-face keys |
|---|---|---|---|
| `extrapolation` | `BC_EXTRAPOLATION = 1` | zeroth-order extrapolation (outflow) | none |
| `inflow` | `BC_INFLOW = 2` | prescribed state | `r, u, v, w, p` (+ `bx, by, bz` for MHD), all required |
| `inflow-digital-filter` | `BC_INFLOW_DF = 7` | turbulent boundary-layer inflow: a mean profile table plus digital-filter fluctuations ([turbulence](./turbulence#digital-filter-inflow)) | `r, u, p`, `df_profile`, `df_delta`, `df_spacing`, `df_wall_axis`, `df_length_*` (required), `df_seed`; Euler only, one face per realm ([input](./input#bc-x-min-bc-x-max-bc-y-min-bc-y-max-bc-z-min-bc-z-max), [#119](https://github.com/szaghi/adam/issues/119)) |
| `wall-inviscid` | `BC_WALL_INVISCID = 3` | slip wall | none |
| `wall-modelled` | `BC_WALL_INVISCID = 3` + wall-model flag | isothermal wall whose face flux comes from a wall model; ghosts mirrored as `wall-inviscid` ([turbulence](./turbulence#wall-model)) | `wall_temperature` (required, $> 0$), `wall_model`, `wall_model_height`, `wall_model_damping`, `wall_model_filter`, `wall_model_eddy_viscosity`, `wall_model_pressure_gradient`, `wall_model_pressure_gradient_filter`, `wall_model_edge_velocity`, `wall_model_edge_temperature` ([input](./input#bc-x-min-bc-x-max-bc-y-min-bc-y-max-bc-z-min-bc-z-max), [#119](https://github.com/szaghi/adam/issues/119)) |
| `wall-noslip` | `BC_WALL_NOSLIP = 4` | no-slip adiabatic wall | `wall_u, wall_v, wall_w` (optional, default 0; the normal one must be 0); refused with `[eos]` `isothermal` or `barotropic` |
| `wall-isothermal` | `BC_WALL_ISOTHERMAL = 5` | no-slip wall at a set temperature | `wall_temperature` (required, $> 0$), `wall_u, wall_v, wall_w` as above; refused with `[eos]` `isothermal` or `barotropic` |
| `periodic` | library `BC_PERIODIC = -1` | periodic pair of faces | none |
| `hydrostatic` | `BC_HYDROSTATIC = 6` | hydrostatic ghosts under gravity, mirrored velocity | none (reads `[gravity]`); refused with `[eos]` `isothermal` or `barotropic` |
| `non-reflecting` | `BC_NONREFLECTING = 8` | characteristic outflow and far field | the target `r, u, v, w, p` (required), `nr_length` $L$ (required, $> 0$), `nr_sigma` $\sigma$ (default 0.28); Euler only ([#121](https://github.com/szaghi/adam/issues/121) S2) |
| `shock-generator` | `BC_SHOCK_GENERATOR = 9` | oblique shock entering through a y face | the upstream `r, u, v, w, p`, one of `shock_deflection` and `shock_angle` (degrees), `shock_x`, `shock_x_end` (optional); `bc_y_min` and `bc_y_max` only, Euler only ([#121](https://github.com/szaghi/adam/issues/121) S3) |

*Validation:*
- an unknown value stops the run and the message lists the accepted spellings;
- a missing inflow key, or a missing or non-positive `wall_temperature`, is fatal;
- a non-zero normal wall velocity is fatal (the wall does not move through the domain);
- periodicity must be **paired**: `[bc_x_min]` and `[bc_x_max]` must both be periodic or neither, and the same for y and z;
- a physical model other than Euler or MHD is fatal.

## How the ghost cells are filled

Before every residual evaluation the ghosts are filled in the order of a Runge–Kutta stage: the copies between the
blocks of the realm (with the 2:1 coarse–fine interpolation), the MPI exchange, then the boundary conditions, crown by
crown from the innermost ghost layer outwards. The inter-realm seams are filled by the forest. The sections below give
the rule of each kind; the crown maps, the edge and corner passes and the seam buffers are documented on the developer
page [FLUME ghost fill](/developer/flume-ghost-fill).

### `extrapolation`: zeroth order (constant)

$$
\mathbf q_{\text{ghost}}(i,j,k) = \mathbf q(i-\delta_i,\,j-\delta_j,\,k-\delta_k)
$$

Each ghost copies its inward neighbour. The crowns are filled innermost first, so every ghost column ends up equal to
the first interior cell: $\mathbf q(1-g) = \mathbf q(1)$ for $g = 1..n_{gc}$. This is order-0 extrapolation, a
zero-gradient outflow. All $n_v$ variables are copied, including
$\psi$.

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

The sign vector $\mathbf S_d$ is:

| variable | Euler | MHD / MHD-GLM |
|---|---|---|
| $\rho$, $E$ | $+1$ (even) | $+1$ |
| $\rho u_n$ (the normal momentum of direction $d$) | $-1$ (odd) | $-1$ |
| $\rho u_{t}$ (tangential momenta) | $+1$ | $+1$ |
| $B_n$ (the normal field of direction $d$) | — | $-1$ (odd), a perfectly conducting reflecting wall ([#41](https://github.com/szaghi/adam/issues/41)) |
| $B_t$ | — | $+1$ |
| $\psi$ | — | $+1$ (**even**) |

$E$ is unchanged by the sign flips because it depends only on squares of the flipped components. The ghost state is
therefore the exact mirror image. As a result, $u_n = 0$ and $B_n = 0$ at the wall face in the sense of the symmetric
average, and $\partial_n\psi = 0$ there.

### `wall-noslip` and `wall-isothermal`: no-slip walls

The two no-slip walls ([#65](https://github.com/szaghi/adam/issues/65)) mirror the cell like `wall-inviscid`
and reflect the **whole** velocity about the wall velocity $\mathbf u_w$ (tangential, `wall_u, wall_v, wall_w`):

$$
\mathbf u_g = 2\,\mathbf u_w - \mathbf u_m, \qquad p_g = p_m,
$$

so the mean of the ghost and its mirror cell is the wall velocity (no slip, no penetration). The pressure is mirrored.

- `wall-noslip` (adiabatic) mirrors the temperature, so $\rho_g = \rho_m$ and $\partial_n T = 0$ at the wall in the
  sense of the symmetric difference.
- `wall-isothermal` sets the ghost temperature to the geometric mirror, $T_g = T_w^2 / T_m$, so the geometric mean
  of the ghost and its mirror cell is the wall temperature, and the density from the mirrored pressure,
  $\rho_g = p_m / (R\,T_g)$, with $T = p/(\rho R)$. The ghost is admissible for any admissible cell, which the linear
  mirror $T_g = 2\,T_w - T_m$ is not beside gas hotter than twice the wall. For $T_m = T_w + \delta$ the two mirrors
  differ by $\delta^2/T_w$, so both are second-order at the wall.

The total energy is rebuilt from the ghost state, $E_g = p_g/(\gamma-1) + \tfrac12\rho_g|\mathbf u_g|^2 + e_{mag}$, where
the magnetic energy $e_{mag}$ (and $\psi^2/2$ under EGLM) is unchanged by the mirror. The field follows the
`wall-inviscid` rule (normal component odd, a perfectly conducting wall); $\psi$ is even. On a resting adiabatic wall the
rule reduces to a sign vector with every momentum component odd. Both walls are second-order accurate at the wall, by
construction of the mirror.

The rule is model-agnostic: the field components are those present in the state vector (none for Euler). The ghost
probe GP checks it on every face and edge, moving and resting walls, Euler and MHD with EGLM
([verification](./verification/infrastructure/ghost-probe)).

The walls are boundary conditions of the ideal solver too: they take effect on any run. Their physical test is the
compressible Couette flow between an isothermal wall and a moving adiabatic one, against its exact profiles
([VV-3](./verification/dissipative/couette)): the moving-wall mirror is exact for its linear velocity and
symmetric temperature, the isothermal mirror $T_w^2/T$ is second-order, so the run converges at order 2.

### `hydrostatic`: an atmosphere under gravity ([#73](https://github.com/szaghi/adam/issues/73))

The ghost density and pressure continue the isothermal hydrostatic equilibrium of the first interior cell along the
normal, $\rho_g = \rho_1 \exp(g_d s/T_1)$, $p_g = \rho_g T_1$, with $T_1 = p_1/\rho_1$ and $s$ the signed distance from
that cell to the ghost along the axis $d$: an isothermal atmosphere has exact ghosts. The velocity is mirrored as at
`wall-inviscid` (no mass flux); the magnetic field and $\psi$ are mirrored even. Without gravity along the normal it is a
slip wall that copies the first interior cell's thermodynamic state (Zingale et al. 2002).

<a id="non-reflecting-a-characteristic-outflow-and-far-field-121-s2"></a>

### `non-reflecting`: a characteristic outflow and far field ([#121](https://github.com/szaghi/adam/issues/121))

The ghosts are the first interior cell plus the jumps of the characteristic variables, linearised about it (outward
normal velocity $u_n$, sound speed $c$):

- **Outgoing waves** (speed $> 0$): no jump, extrapolated.
- **The incoming acoustic wave**, $w_1 = p - \rho c\,u_n$ with speed $u_n - c < 0$: its gradient along the normal is set by
  the relaxation of Rudy & Strikwerda (1980) and Poinsot & Lele (1992),
  $\partial w_1/\partial n = K(p - p_t)/(u_n - c)$, $K = \sigma(1 - M^2)c/L$. Ghost $m$ lies $m$ cells off the interior
  cell, and $w_5 = p + \rho c\,u_n$ is unchanged.
- **Flow into the domain** ($u_n < 0$): the entropy $\rho - p/c^2$ and the tangential velocity of the target;
  $u_n + c \le 0$ (supersonic inflow) the whole target.

A supersonic outflow, or $\sigma = 0$ at an outflow, is extrapolation bit for bit. The edge and corner rows beyond the
face are extrapolated.

**When to use it.** In this characteristic-upwind scheme `extrapolation` (outflow) and `inflow` (a prescribed state) already
let linear waves out: they make the incoming characteristic jumps zero, or the free stream's. What they cannot do is hold a
mean pressure at a subsonic outflow, and `extrapolation` keeps whatever pressure the box has. `non-reflecting` pulls it to
$p_t$ on the time scale $1/K$. The price is a reflection of order $\sigma c/(\omega L)$ of a wave of frequency $\omega$
passing through. Measured in the [verification](./verification/streams/non-reflecting):

| | `non-reflecting` | the alternative |
|---|---|---|
| a 2 % overpressure at Mach 0.5, left after $t = 20 = 10\,L/c$ | 0.099 | `extrapolation`: 1.000 |
| a packet leaving at rest ($\lambda = 0.2$, $L = 2$), left behind | 2.9e-3 ($\sigma = 0.28$), 7.6e-4 ($\sigma = 0.07$) | `extrapolation`: 1.2e-4 |
| a packet leaving through a subsonic inflow face, left behind | 4.3e-3 | `inflow`: 1.4e-5 |


<a id="shock-generator-an-oblique-shock-through-a-y-face-121-s3"></a>

### `shock-generator`: an oblique shock through a y face ([#121](https://github.com/szaghi/adam/issues/121))

The face of a shock–boundary-layer interaction that brings the incident shock in (the top of the box). The ghosts with
`shock_x` $\le x <$ `shock_x_end` take the state behind an oblique shock in the upstream state `r, u, v, w, p`; the others are
extrapolated, as are the edge and corner rows beyond the face. The shock is given by its deflection $\theta$
(`shock_deflection`, the weak solution of the $\theta$–$\beta$–$M$ relation) or by its angle $\beta$ from the upstream velocity
(`shock_angle`), and turns the flow toward the opposite face. The post-shock state follows from the normal Mach number
$M\sin\beta$, with $M$ of the velocity in the x-y plane (Rankine–Hugoniot; Anderson, *Modern Compressible Flow*, ch. 4); $w$, tangent
to the shock, is unchanged. The run prints the shock angle, the deflection and $p_2/p_1$.

- **Refused:** a flow that is not supersonic in the x-y plane, an angle outside $(\mu, 90^\circ]$, a deflection beyond detachment,
  a face other than y, a model other than Euler (no dual energy, no passive scalars).
- **The rest of the face is `extrapolation`**: upstream of `shock_x` the free stream leaves or enters at its own state, as at an
  extrapolated top. A reflected shock that reaches the face beyond `shock_x_end` leaves through it.

Measured in the [verification](./verification/streams/shock-generator): at Mach 2 and $10^\circ$ the shock angle
within $0.01^\circ$ of $39.314^\circ$, the state behind it within $1.5\cdot10^{-3}$ of Rankine–Hugoniot; reflected off a `wall-inviscid`
face, the reflected shock within $0.06^\circ$ of the exact regular reflection and the state behind it within $2\cdot10^{-3}$.


### `wall-modelled`

A slip-wall mirror for the ghosts, with the flux through the wall face replaced by the stress and the heat flux of a
wall model sampled off the wall: an isothermal wall for wall-modelled large-eddy simulation. The model, its keys and its
limits are on the [turbulence](./turbulence#wall-model) page.

### `inflow-digital-filter`

A turbulent boundary layer entering through the face: the mean of a profile table plus digital-filter fluctuations with
its Reynolds stresses. The method is on the [turbulence](./turbulence#digital-filter-inflow) page.

### `periodic`

`periodic` is the library periodicity: the tree builds true periodic neighbours, and the block exchange fills the ghosts
across blocks and ranks, with any number of blocks along the axis, one included
([#103](https://github.com/szaghi/adam/issues/103)). The only constraint is the pairing rule of the validation list
above.

### Edges and corners {#edges-and-corners-fec-6}

An edge ghost lies outside its block along two axes, a corner ghost along three. The WENO stencils never read them, but
the tangential derivatives of the viscous stress or of the current at a face do. They are filled so that the ghost layer
continues the face ghosts: beyond one realm face the face's kind applies along its normal; beyond two or three realm
faces the kind of one physical face applies (inflow first), so two walls compose to a double mirror. The ghost probe GP
holds every face and edge ghost to a linear field on every path
([verification](./verification/infrastructure/ghost-probe)). The donor rules and the pass ordering are on the
[developer page](/developer/flume-ghost-fill#edges-and-corners).

### Inter-realm seam

A face glued to another realm by the forest manifest still needs a `type`, but the seam overrides it: the forest fills
the seam ghosts by injection from the peer realm, at every Runge–Kutta stage or once per step (`coupling_cadence`). The
written field files hold the ghosts the stencils read, except the step-0 file, whose seam ghosts are unfilled
([developer page](/developer/flume-ghost-fill#inter-realm-seam)).

### 2:1 AMR coarse-fine ghosts

These ghosts are library behaviour, selected by `[amr] seam_ghost_fill = injection | restriction-compatible | tricubic`
(default `tricubic`; an unknown value is fatal; [library input](/reference/library-input#amr)). FLUME bounds the fill
and blends inadmissible ghosts ([numerics](./numerics#coarse-fine-ghosts)).

### Immersed solids

Immersed solids are not a `[bc_*]` kind. They enter through a level set and the cut spacing of the flux difference
([numerics](./numerics#immersed-boundary-euler)).

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
