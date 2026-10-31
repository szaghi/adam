!< ADAM, FNL device kernels of the dynamic eddy viscosity (issue #119, T1).

#include "fundal.H"

module adam_fnl_turbulence_kernels
!< ADAM, FNL device kernels of the dynamic eddy viscosity (issue #119, T1).
!<
!< Device twins of `adam_turbulence_kernels`: see there. Device arrays are transposed, `(b, i, j, k, v)`.

! ADAM modules
use :: adam_turbulence_library, only : sgs_sigma, sgs_vreman, sgs_wale
! third party modules
use :: penf,                    only : I4P, R8P

implicit none
private
public :: compute_mu_t_sigma_dev
public :: compute_mu_t_vreman_dev
public :: compute_mu_t_wale_dev

contains
   subroutine compute_mu_t_vreman_dev(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz_gpu, is_null, q_aux_gpu, mu_t_gpu)
   !< Vreman.
#define ADAM_SGS_FORMULA sgs_vreman
#include "adam_fnl_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_vreman_dev

   subroutine compute_mu_t_wale_dev(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz_gpu, is_null, q_aux_gpu, mu_t_gpu)
   !< WALE.
#define ADAM_SGS_FORMULA sgs_wale
#include "adam_fnl_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_wale_dev

   subroutine compute_mu_t_sigma_dev(ir, iu, ni, nj, nk, ngc, blocks_number, c, dxyz_gpu, is_null, q_aux_gpu, mu_t_gpu)
   !< Sigma.
#define ADAM_SGS_FORMULA sgs_sigma
#include "adam_fnl_turbulence_body.INC"
#undef ADAM_SGS_FORMULA
   endsubroutine compute_mu_t_sigma_dev
endmodule adam_fnl_turbulence_kernels
