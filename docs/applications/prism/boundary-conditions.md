# Boundary conditions

Each domain face takes one kind, `[bc_x_min]` … `[bc_z_max]` `type` ([input](./input#bc-x-min-bc-x-max-bc-y-min-bc-y-max-bc-z-min-bc-z-max)).
The kinds fill the ghost cells of the boundary blocks after every ghost exchange. Waves
leaving the domain are absorbed by the [absorbing layers](./pml), which work on top of a face kind.

## Face kinds

Below, $g$ is a ghost cell and $m$ its mirror about the face (the interior cell at the same distance); $n$ is the face
normal, $t_1, t_2$ the tangential directions.

| `type` | Ghost values |
|---|---|
| `extrapolation` | every variable copied from the adjacent interior cell (zeroth-order extrapolation) |
| `Neumann` | every variable mirrored, $q_g = q_m$ (zero normal derivative) |
| `Dirichlet` | every variable **zero**; there is no key for a prescribed value |
| `PEC` | mirror with signs: $D_n$, $B_{t_1}$, $B_{t_2}$ even, $D_{t_1}$, $D_{t_2}$, $B_n$ odd; $\varphi$ odd, $\psi$ even |
| `Silver_Muller` | first-order absorbing condition (below) |
| `periodic` | the library's periodicity: the ghost exchange fills the ghosts from the opposite side, across blocks and ranks |

- **PEC** is a perfect electric conductor: $\mathbf E_t = 0$ and $B_n = 0$ on the face, imposed by odd reflection.
- **`periodic`** is `BC_PERIODIC = -1` of the library: the tree builds periodic neighbours, so a periodic face is not a
  boundary row. The two faces of an axis must be both periodic or neither (fatal otherwise).
- **Edge and corner ghosts** get no boundary value: only the face ghost rows are set.
- **PIC models**: the boundary ghosts carry no deposited charge or current ($\mathbf J$ and $\rho$ set to zero there).
- `radiative` is refused at load ([#57](https://github.com/szaghi/adam/issues/57)); an unknown spelling is fatal.

## Silver–Müller

The first-order Silver–Müller condition imposes the impedance relation of an outgoing plane wave, $|\mathbf E_t| = c_0|\mathbf B_t|$,
so a normally incident wave leaves the domain. The tangential ghost components come from the mirror cell
with the vacuum admittance $Y = \varepsilon_0 c_0$ (1 non-dimensional), the face sign $s = \pm1$ ($-1$ on a minimum
face) and $(\alpha, \beta, \gamma)$ the two tangential directions and the normal in cyclic order:

$$D^g_\alpha = s\,Y B^m_\beta, \qquad D^g_\beta = -s\,Y B^m_\alpha, \qquad
B^g_\alpha = -s\,D^m_\beta / Y, \qquad B^g_\beta = s\,D^m_\alpha / Y .$$

The normal ghost components are then solved, line by line along the normal, so that the discrete divergence of the
first $h_s = $ `fdv_order`$/2$ interior cells is $\rho$ for $\mathbf D$ (zero without PIC) and zero for $\mathbf B$: the boundary does not
create divergence. This needs `ngc >= fdv_order/2`.

## The elliptic boundary map

The t0 elliptic solves ([numerics](./numerics#t0-elliptic-solves)) see each face through a second, derived kind:
periodic → periodic, `Dirichlet` → zero ghosts, `PEC` → Dirichlet for $\varphi$ and the tangential $\mathbf A$ and
Neumann for the normal $\mathbf A$, `Neumann` and `Silver_Muller` → the free-space (open) potential. `extrapolation` has
no elliptic counterpart and stops the run **at the first t0 solve**, not at load: the table and its consequences are in
the [input hazards](./input#elliptic-boundary-map).

## Particle boundaries

The particles' closure at a face follows its field kind: a periodic face wraps the particles, every other kind keeps
the legacy closure, where a particle leaving the domain stops the run ([PIC](./pic#particle-boundaries)).

When adding a face kind, update the three `select case`s that see it: the field parser (`adam_prism_bc_object`), the
elliptic map (`build_elliptic_bc_types`) and the FLAIL map (`solve_time_zero`); and the particle closure
(`adam_prism_pic_bc_object`), whose map is total.
