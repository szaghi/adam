# FLUME ghost fill and AMR gravity

This page documents two pieces of FLUME internals: how the ghost cells of a realm are filled (crown maps, the
edge and corner passes, the inter-realm seam buffers), and how self-gravity keeps its energy exchange exact on AMR
grids. The user-facing rules are on the [boundary conditions](../applications/flume/boundary-conditions) and
[physical models](../applications/flume/models#self-gravity) pages.

## Ghost fill order

`update_ghost` runs three steps in order:

1. the intra-realm local copies (`field%update_ghost_local`, which also runs the 2:1 coarse–fine fill);
2. the MPI exchange (`update_ghost_mpi`);
3. `set_boundary_conditions`.

The forest fills the inter-realm seam ghosts separately ([below](#inter-realm-seam)).

The grid receives the six BC ids through `grid%set_bc_type`. That call turns on `is_ijk_periodic` for a periodic axis.
`periodic` gets no crown rows: the library tree builds true periodic neighbours, and steps 1 and 2 fill the ghosts
across blocks and ranks, with any number of blocks along the axis. The tree and the maps are rebuilt once the BC types
are known ([#103](https://github.com/szaghi/adam/issues/103)).

## Crown maps

Step 3 walks the boundary crown map `maps%local_map_bc_crown` crown by crown, $c = 1..n_{gc}$, from the innermost crown
outwards. Each ghost cell $(i,j,k)$ carries:

- an inward step $(\delta_i,\delta_j,\delta_k)$;
- its BC type;
- its `fec` code: face 1–6, edge 7–18, corner 19–26;
- the boundary its cell lies beyond (`bc_fec`, from the tree).

The CPU implementation is `cpu/adam_flume_cpu_object.F90` (`set_boundary_conditions`). The FNL twin is
`fnl/adam_flume_fnl_kernels.F90` (`set_boundary_conditions_dev`), called crown by crown from
`fnl/adam_flume_fnl_object.F90`. Both apply the same rules. The per-kind ghost formulas are pure procedures of
`adam_flume_bc_object.F90` (`wall_noslip_ghost`, `hydrostatic_ghost`, `nonreflecting_ghost`, `oblique_shock`) shared by
both backends.

Because the crowns are filled innermost first, an `extrapolation` ghost copying its inward neighbour ends equal to the
first interior cell: $\mathbf q(1-g) = \mathbf q(1)$ for $g = 1..n_{gc}$.

## Edges and corners

An edge ghost lies outside its block along two axes, a corner ghost along three. The directional WENO stencils never
read them, but a cross derivative does: the tangential derivatives of the viscous stress or of the current at a face
read the edge ghosts of the rows next to a block edge. They are filled so that the ghost layer continues the face
ghosts. The ghost probe GP holds every face and edge ghost to a linear field on every path
([verification](../applications/flume/verification/infrastructure/ghost-probe)).

`bc_fec` gives two cases:

- **Beyond one realm face** (`bc_fec` a face). The cell's other directions lead into another block, whose ghosts the
  exchange has filled. The face's kind applies along its normal exactly as for a face ghost, so a wall edge mirrors the
  exchanged ghost and negates its normal momentum.
- **Beyond two or three realm faces** (`bc_fec` an edge or a corner of the realm). The kind of one physical face among
  them applies: inflow first, since an inflow ghost holds the inflow state whatever else it lies beyond; else the first
  non-seam face, `realm_edge_face`. It acts on the donor `realm_edge_donor`: the cell mirrored about that face (any wall
  kind) or the first interior cell along its normal (extrapolation), with the indexes along the other axes kept. The
  donor therefore lies beyond the other faces only. It is a face ghost (for an edge) or an edge ghost (for a corner),
  already filled by the seam, the exchange or its own face. Two walls compose to a double mirror with both normal
  momenta negated, and a wall beside a seam mirrors the seam ghost.

A cell beyond seams only lies outside every realm: the re-entrant corner of an L-shaped forest, inside the step of the
[Woodward–Colella tunnel](../applications/flume/verification/euler/forward-facing-step). It keeps a copy of its inward
diagonal, $\mathbf q(i-\delta_i, j-\delta_j, k-\delta_k)$: finite, and read by no stencil of a fluid cell. Every other
row beyond two realm faces uses a donor; the inward-diagonal copy there is off by up to 0.78 of the field scale at a
wall corner on a linear field (wrong sign of the normal momentum, wrong position), measured by GP on the step forest.

## Pass ordering

Both backends fill the crown in three passes, the rows beyond one face, then two, then three, each crown by crown, so
that every donor is filled before it is read. Within the first pass no row reads another row of the pass:

- a wall row reads the mirrored interior cell;
- an extrapolation row reads the **first interior cell** along its normal, not the previous ghost of the chain. Both
  give the same value, but the previous ghost is a row of the same launch: on FNL a chained read races, and a ghost
  beyond a block interface along another axis can read its donor before it is written.

The kind of each realm face comes from `[bc_*]`, overridden by `maps%seam_face` for the faces the forest glued. On FNL
each launch visits the rows of its own pass only: the rows of every crown are grouped by pass on the host whenever the
device map is copied (`build_bc_pass_rows`), and an empty pass is not launched
([#90](https://github.com/szaghi/adam/issues/90)).

## Inter-realm seam

`set_boundary_conditions` leaves seam face rows untouched. The forest fills them through `fill_seam_from_peer_forest`,
a plain injection copy of the peer's interior cell into this realm's ghost cell, following the seam ghost map rows. The
buffer read is `q` when `stage_active == 0`, otherwise the active stage of `rk%q_rk`, on each side independently.

`update_ghost` does not fill seam ghosts, so any diagnostic that reads them refills them explicitly;
`compute_divb_history` does. The field write does too: with sibling realms, `save_simulation_data` refills the seams and
then calls `update_ghost`, the order of a Runge–Kutta stage, so the written ghosts (edges included) are those the
stencils read. The step-0 file is the exception: it is written while each realm initialises, before the forest connects
the realms, and its seam ghosts are unfilled.

## 2:1 coarse–fine ghosts

The coarse–fine ghosts are library behaviour, filled inside `field%update_ghost_local`. The regime is the optional key
`[amr] seam_ghost_fill = injection | restriction-compatible | tricubic` (default `tricubic`; an unknown value is fatal),
documented in the [library input reference](/reference/library-input#amr). FLUME adds the density/velocity clamp of the
fill and the admissibility blend of the face ghosts
([numerics](../applications/flume/numerics#adaptive-mesh-refinement)).

Immersed solids are not a ghost-fill path: they enter through the level set `phi` and the cut spacing in
`compute_flux_difference` ([numerics](../applications/flume/numerics#immersed-boundary-euler)).

## Self-gravity on AMR grids

The equations, the uniform-grid conservation argument and the solvers are on the
[models page](../applications/flume/models#self-gravity). This section derives the composite operator, the seam work
that keeps $E + W$ exact, and the isolated boundaries on AMR.

### Composite Laplacian

On AMR grids $L_h$ is the composite Laplacian of the leaves. Across a 2:1 face it carries the **mimetic mean flux**:
every fine face $f$ of a coarse cell C carries $\tfrac23 w\,(\bar\Phi_f - \Phi_C)$, where $\bar\Phi_f$ is the mean of
the face's fine cells and $w$ is the fine face area over the fine spacing, and C receives the sum of these fluxes. The
operator is symmetric in the volume-weighted form, so $\Phi = 4\pi G\,K\rho$ with $VK$ symmetric, and it is consistent
(exact for linear $\Phi$ face by face). [SV-9](../applications/flume/verification/self-gravity/multigrid-amr) explains
why the simpler couplings fail one property or the other.

### Seam work on the reflux

Both sides of the seam see the same face gravity, $-(\bar\Phi_f - \Phi_C)/(1.5h)$, which is first order there. The
energy sources need the two-point work. For $E + W$ to telescope, the seam's sources must sum to
$\sum_f A_fF^\rho_f(\Phi_C - \Phi_f)$ (with the sign of the mass flux), and the two kernels give
$\sum_f A_fF^\rho_f(\Phi_C - \bar\Phi_f)$ instead. The difference, $\sum_f A_fF^\rho_f(\bar\Phi_f - \Phi_f)$, is small
($\sum_f(\bar\Phi_f - \Phi_f) = 0$) but not zero.

The difference rides on the reflux. Every stage adds to the flux register's energy entry:

- the coarse side: its own seam source, taken back;
- the fine side: the exact two-point work less the fine cell's own source.

The end-of-step correction then gives the coarse cell exactly the missing work. The semi-discrete $E + W$ is exact
again, the fine cells keep the consistent gravity, and the coarse cell next to the seam takes an $O(h)$ local error.

### Time order

The fully discrete exchange is first order in $\Delta t$, not of the Runge–Kutta order. The reflux corrects the coarse
cells only at the end of the step, so its stages advance them with their own seam fluxes, an $O(\Delta t)$
inconsistency per step. Measured on the collapsing Jeans mode
([SV-10](../applications/flume/verification/self-gravity/amr-jeans)):

- with the seam work: $9\cdot10^{-8}$, $3\cdot10^{-8}$, $1.4\cdot10^{-8}$ at CFL 0.4, 0.2, 0.1;
- without it: $10^{-5}$ at every CFL, a spatial leak.

`[numerics] reflux_cadence = stage` restores the Runge–Kutta order: every stage corrects the coarse cells' residual from
that stage's fluxes, and the drift falls at third order in $\Delta t$ ($3.5\cdot10^{-8}$, $4.3\cdot10^{-9}$,
$5.4\cdot10^{-10}$). `reflux = .false.` or `reflux_policy = troubled` leave the work uncorrected where they leave the
flux.

### Runtime AMR

Every regrid sets the solver up on the new grid (FFT when the grid is uniform, multigrid otherwise) and solves $\Phi$
from the regridded state. The regrid conserves $\rho$, $\rho\mathbf u$ and $E$, but $W = \tfrac12\sum\rho\Phi V$ is the
grid's own discrete energy: a regrid that changes the grid changes $W$ by a truncation-sized amount, and $E + W$ jumps
($2\cdot10^{-4}$ of $E + W$ on the collapse at base $32^2$,
[SV-11](../applications/flume/verification/self-gravity/runtime-amr)). The scheme between regrids stays exact. Making
the regrid itself conserve $E + W$ needs a correction of $E$ at each regrid, which is not done.

### Isolated boundaries on AMR

Isolated boundaries on AMR ([#87](https://github.com/szaghi/adam/issues/87)) use a symmetric form of James's method.
With $A$ the composite operator with zero ghosts at the domain faces and $B$ the map of ghost values onto the boundary
rows, the zero-ghost solve $\Phi_0 = A^{-1}4\pi G V\rho$, extended by zero, satisfies the infinite operator with a
screening charge $s = B^T\Phi_0$ on the ghost cells. The free-space ghost values are $g = -Gs$, $G$ the lattice Green's
function from ghost cell to ghost cell (a zero-padded FFT of the coarsest grid widened by the ghost layer), and the
second solve is $\Phi = A^{-1}(4\pi G V\rho - Bg)$:

$$K = A^{-1} + A^{-1}BGB^TA^{-1},$$

symmetric because $A$ and $G$ are, so the $E + W$ exchange stays exact. On a uniform grid it is the discrete free-space
solution of the FFT solver; on AMR it is second order ($G$ ignores the refinement inside the domain, a truncation-order
error). The textbook James method, with the screening charge on the boundary surface and the potential at the ghost
centres, is not symmetric.

$G$ is the coarse lattice's only if **every leaf on the domain boundary is at the coarsest level**. Anything else is
fatal, at start and at the regrid that would refine a boundary block; an isolated problem keeps its mass away from the
boundary, where refinement is not needed. Each Poisson solve costs two composite solves and one FFT of the coarsest
grid ([SV-13](../applications/flume/verification/self-gravity/multigrid-amr),
[SV-14](../applications/flume/verification/self-gravity/star-amr)).
