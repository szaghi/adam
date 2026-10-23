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
use :: penf, only : I4P

implicit none
private
public :: FLAIL_BC_PERIODIC
public :: FLAIL_BC_DIRICHLET
public :: FLAIL_BC_NEUMANN
public :: FLAIL_BC_DIRICHLET_GHOST

integer(I4P), parameter :: FLAIL_BC_PERIODIC        = 0_I4P !< Periodic (both faces of the axis).
integer(I4P), parameter :: FLAIL_BC_DIRICHLET       = 1_I4P !< Value at the face.
integer(I4P), parameter :: FLAIL_BC_NEUMANN         = 2_I4P !< Normal derivative at the face.
integer(I4P), parameter :: FLAIL_BC_DIRICHLET_GHOST = 3_I4P !< Value at the first ghost centre (issue #87, internal).
endmodule adam_flail_parameters
