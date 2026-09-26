!< ADAM, FLUME FNL device kernels of the ideal MHD with GLM cleaning model.

#include "fundal.H"

module adam_flume_fnl_mhd_glm_kernels
!< ADAM, FLUME FNL device kernels of the ideal MHD with GLM cleaning model.
!<
!< The model-agnostic auxiliary-variables kernel (`adam_flume_fnl_aux_kernels_agnostic.INC`) and the MHD dt and
!< conservation kernels (`adam_flume_fnl_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD_GLM`,
!< `NV_AUX_K = NV_AUX_MHD` (issue #41, section 4). The face-flux kernel is the shared body with the MHD
!< split (M2-P3). Same kernel rules as `adam_flume_fnl_kernels` (issue #35, D-11/D-12).

! ADAM FNL classes, libraries
use :: adam_fnl_weno_kernels,  only : weno_reconstruct_upwind_dev
! FLUME modules
use :: adam_flume_mhd_library, only : compute_face_flux_back_projection=>mhd_glm_face_flux_back_projection, &
                                      conservative_to_auxiliary=>mhd_conservative_to_auxiliary,             &
                                      mhd_fast_speed, mhd_glm_face_split_fluxes, mhd_sum3
use :: adam_flume_parameters,  only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                      IQ_RW, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_GLM, S_MAX
! third party modules
use :: penf,                   only : I4P, R8P

implicit none
private
public :: apply_floors_dev
public :: compute_conservation_dev
public :: compute_face_fluxes_dev
public :: compute_lambda_max_dev
public :: compute_q_aux_dev

contains
   ! public procedures
#include "adam_flume_fnl_face_kernels_agnostic.INC"

#include "adam_flume_fnl_aux_kernels_agnostic.INC"

#include "adam_flume_fnl_mhd_kernels_agnostic.INC"

   ! private procedures
   pure subroutine face_split_fluxes(gamma, ch, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Split adapter of the shared face kernel (issue #41, M2-P3): the MHD split with GLM cleaning.
   real(R8P),    intent(in)  :: gamma                          !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch                             !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d                              !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                              !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic              !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_K,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_K,1-S_MAX:S_MAX)    !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_K) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_K,NV_K)                  !< Right eigenvectors.
   !$acc routine seq
   !$omp declare target

   call mhd_glm_face_split_fluxes(ch=ch, gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                                  fsplit=fsplit, er=er)
   endsubroutine face_split_fluxes
endmodule adam_flume_fnl_mhd_glm_kernels
