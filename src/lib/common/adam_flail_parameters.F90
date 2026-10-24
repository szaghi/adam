!< ADAM, FLAIL elliptic solver library: shared parameters (issue #84).
module adam_flail_parameters
!< ADAM, FLAIL elliptic solver library: shared parameters (issue #84).
!<
!< The boundary-condition kinds of the elliptic solvers, one per domain face, shared by the FFT solver, the multigrid
!< solver and the FLAIL facade (a module of its own: the solvers cannot use the facade, which uses them).
!<   - `FLAIL_BC_PERIODIC`: both faces of the axis; the library's periodicity (`grid%is_ijk_periodic`).
!<   - `FLAIL_BC_DIRICHLET`: the value at the **face** (homogeneous in the solvers: `phi_ghost = -phi_boundary`;
!<     second order: the wall at the face). The FFT solver extends the source odd about the face.
!<   - `FLAIL_BC_NEUMANN`: the normal derivative at the face (homogeneous in the solvers: `phi_ghost = phi_boundary`).
!<     The FFT solver extends the source even about the face. Without a Dirichlet face anywhere the operator is
!<     singular: the solve drops the source mean, as on a periodic grid.
!<   - `FLAIL_BC_DIRICHLET_GHOST`: the value at the centre of the first ghost cell (issue #87: the multigrid's
!<     coarsest level and its James step, where the ghost values are computed, not prescribed). The wall sits half a
!<     cell outside the domain, so as a prescribed condition it is first order; not for clients.
!<   - `FLAIL_BC_ZERO_GHOSTS`: every ghost layer zero (issue #84, E2b): PRISM's elliptic Dirichlet face, which its D
!<     and G read the same way. With the compact operator it is `FLAIL_BC_DIRICHLET_GHOST` (the 7-point stencil reads
!<     one ghost layer); with the wide D(G) it is no reflection, and the FFT solver solves it by preconditioned CG.
!<   - `FLAIL_BC_FREE_GHOSTS`: every ghost layer holds the continuum free-space potential of the source,
!<     `-sum_y factor f(y) dV / (4 pi |x - y|)` (issue #84, E3c): PRISM's exact/open face. The operator is the
!<     zero-ghost one (the two kinds mix on an axis), the ghost values enter the right-hand side; the FFT solver with the
!<     wide D(G) only. Rank 0 sums over the gathered source: the values do not depend on the decomposition.
use :: penf, only : I4P

implicit none
private
public :: FLAIL_BC_PERIODIC
public :: FLAIL_BC_DIRICHLET
public :: FLAIL_BC_NEUMANN
public :: FLAIL_BC_DIRICHLET_GHOST
public :: FLAIL_BC_ZERO_GHOSTS
public :: FLAIL_BC_FREE_GHOSTS

integer(I4P), parameter :: FLAIL_BC_PERIODIC        = 0_I4P !< Periodic (both faces of the axis).
integer(I4P), parameter :: FLAIL_BC_DIRICHLET       = 1_I4P !< Value at the face.
integer(I4P), parameter :: FLAIL_BC_NEUMANN         = 2_I4P !< Normal derivative at the face.
integer(I4P), parameter :: FLAIL_BC_DIRICHLET_GHOST = 3_I4P !< Value at the first ghost centre (issue #87, internal).
integer(I4P), parameter :: FLAIL_BC_ZERO_GHOSTS     = 4_I4P !< Every ghost layer zero (PRISM's Dirichlet, E2b).
integer(I4P), parameter :: FLAIL_BC_FREE_GHOSTS     = 5_I4P !< Ghosts: the free-space potential (PRISM's exact/open, E3c).
endmodule adam_flail_parameters
