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
use :: adam_flume_parameters,  only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                      IQ_RW, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_GLM, S_MAX
! third party modules
use :: penf,                   only : I4P, R8P

implicit none
private
public :: add_glm_damping_dev
public :: apply_floors_dev
public :: compute_conservation_dev
public :: compute_face_fluxes_dev
public :: compute_lambda_max_dev
public :: compute_q_aux_dev
public :: compute_speed_max_dev

contains
   ! public procedures
#include "adam_flume_fnl_face_kernels_agnostic.INC"

#include "adam_flume_fnl_aux_kernels_agnostic.INC"

#include "adam_flume_fnl_mhd_kernels_agnostic.INC"

   subroutine add_glm_damping_dev(ni, nj, nk, ngc, blocks_number, damping, q_gpu, dq_gpu)
   !< Add the GLM damping source to the residuals of the interior cells, `dq(psi) = dq(psi) - (c_h^2 / c_p^2) psi`
   !< (mixed GLM, Dedner et al. 2002; issue #41, section 3.1): a local source, in every Runge-Kutta stage and never in
   !< the fluxes, so the reflux and the conservation of the other variables are untouched.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   real(R8P),    intent(in)    :: damping                            !< Damping rate c_h^2 / c_p^2.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   integer(I4P)                :: b, i, j, k                         !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,damping)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,damping)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dq_gpu(b,i,j,k,IQ_PSI) = dq_gpu(b,i,j,k,IQ_PSI) - damping * q_gpu(b,i,j,k,IQ_PSI)
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_glm_damping_dev

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
