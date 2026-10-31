!< ADAM, CPU kernels of the dynamic eddy viscosity (issue #119, T1).

module adam_turbulence_kernels
!< ADAM, CPU kernels of the dynamic eddy viscosity (issue #119, T1).
!<
!< `mu_t = rho nu_t` on the interior cells and on 2 ghost layers per active direction: the 4th-order dissipative fluxes read it
!< at the cells -1 ... 2 along their direction. 2nd-order central velocity gradients, so `ngc >= 3`. The filter width is the
!< geometric mean of the active steps of the block, and so it jumps at a 2:1 face. One kernel per model (the same body), selected
!< on the host; the app passes the auxiliary indexes of the density and of the first velocity component.

! ADAM modules
use :: adam_turbulence_library, only : sgs_sigma, sgs_vreman, sgs_wale
! third party modules
use :: penf,                    only : I4P, R8P

implicit none
private
public :: compute_mu_t_sigma
public :: compute_mu_t_vreman
public :: compute_mu_t_wale

contains
   subroutine compute_mu_t_vreman(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz, is_null, q_aux, mu_t)
   !< Vreman.
#define ADAM_SGS_FORMULA sgs_vreman
#include "adam_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_vreman

   subroutine compute_mu_t_wale(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz, is_null, q_aux, mu_t)
   !< WALE.
#define ADAM_SGS_FORMULA sgs_wale
#include "adam_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_wale

   subroutine compute_mu_t_sigma(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz, is_null, q_aux, mu_t)
   !< Sigma.
#define ADAM_SGS_FORMULA sgs_sigma
#include "adam_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_sigma
endmodule adam_turbulence_kernels
