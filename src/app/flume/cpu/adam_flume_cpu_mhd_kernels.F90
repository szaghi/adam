!< ADAM, FLUME CPU kernels of the ideal MHD (no divergence control) model.

module adam_flume_cpu_mhd_kernels
!< ADAM, FLUME CPU kernels of the ideal MHD (no divergence control) model.
!<
!< The model-agnostic auxiliary-variables loop (`adam_flume_cpu_aux_kernels_agnostic.INC`) and the MHD dt loop
!< (`adam_flume_cpu_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD`, `NV_AUX_K = NV_AUX_MHD` (issue #41,
!< section 4). No face fluxes yet: the MHD eigensystem and fluxes come with M2-P2/P3.

! FLUME modules
use :: adam_flume_mhd_library, only : conservative_to_auxiliary=>mhd_conservative_to_auxiliary, mhd_fast_speed
use :: adam_flume_parameters,  only : IA_U, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD
! third party modules
use :: penf,                   only : I4P, R8P

implicit none
private
public :: compute_lambda_max
public :: compute_q_aux

contains
   ! public procedures
#include "adam_flume_cpu_aux_kernels_agnostic.INC"

#include "adam_flume_cpu_mhd_kernels_agnostic.INC"
endmodule adam_flume_cpu_mhd_kernels
