!< ADAM, FLUME FNL device kernels of the ideal MHD with GLM cleaning model.

#include "fundal.H"

module adam_flume_fnl_mhd_glm_kernels
!< ADAM, FLUME FNL device kernels of the ideal MHD with GLM cleaning model.
!<
!< The model-agnostic auxiliary-variables kernel (`adam_flume_fnl_aux_kernels_agnostic.INC`) and the MHD dt and
!< conservation kernels (`adam_flume_fnl_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD_GLM`,
!< `NV_AUX_K = NV_AUX_MHD` (issue #41, section 4). No face fluxes yet: the MHD eigensystem and fluxes come with
!< M2-P2/P3. Same kernel rules as `adam_flume_fnl_kernels` (issue #35, D-11/D-12).

! FLUME modules
use :: adam_flume_mhd_library, only : conservative_to_auxiliary=>mhd_conservative_to_auxiliary, mhd_fast_speed
use :: adam_flume_parameters,  only : IA_U, IA_V, IA_W, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_GLM
! third party modules
use :: penf,                   only : I4P, R8P

implicit none
private
public :: compute_conservation_dev
public :: compute_lambda_max_dev
public :: compute_q_aux_dev

contains
   ! public procedures
#include "adam_flume_fnl_aux_kernels_agnostic.INC"

#include "adam_flume_fnl_mhd_kernels_agnostic.INC"
endmodule adam_flume_fnl_mhd_glm_kernels
