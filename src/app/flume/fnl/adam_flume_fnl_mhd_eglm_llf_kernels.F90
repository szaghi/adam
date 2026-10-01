!< ADAM, FLUME FNL device kernels of the MHD-EGLM model with the LLF Riemann solver (scheme `weno-riemann`).

#include "fundal.H"

module adam_flume_fnl_mhd_eglm_llf_kernels
!< ADAM, FLUME FNL device kernels of the MHD-EGLM model with the LLF Riemann solver (scheme `weno-riemann`).
!<
!< The model- and solver-agnostic face loop (`adam_flume_fnl_riemann_face_kernels_agnostic.INC`, issue #47)
!< instantiated on the MHD-EGLM model and the LLF solver of `adam_flume_mhd_riemann_library`: the local arrays
!< are sized by `NV_K = NV_MHD_EGLM`, `NV_AUX_K = NV_AUX_MHD`, and the four adapters call the MHD libraries.

! ADAM FNL classes, libraries
use :: adam_fnl_weno_kernels,          only : weno_reconstruct_upwind_wratio_dev
! FLUME modules
use :: adam_flume_mhd_library,         only : mhd_eglm_flux
use :: adam_flume_mhd_riemann_library, only : mhd_eglm_face_interpolation_fields, mhd_eglm_face_states, &
                                              mhd_eglm_riemann_llf
use :: adam_flume_parameters,          only : NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_EGLM, S_MAX
! third party modules
use :: penf,                           only : I4P, R8P

implicit none
private
public :: compute_riemann_face_fluxes_dev

contains
   ! public procedures
#include "adam_flume_fnl_riemann_face_kernels_agnostic.INC"

   ! private procedures
   pure subroutine face_interpolation_fields(gamma, ch, d, S, is_characteristic, qs, qas, fint, er)
   !< Interpolation-fields adapter of the shared face kernel: the MHD fields, the EGLM eigenvectors, `ch` unused.
   real(R8P),    intent(in)  :: gamma                        !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch                           !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d                            !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                            !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic            !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_K,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_K,1-S_MAX:S_MAX)  !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_K) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_K,NV_K)                !< Right eigenvectors.
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_face_interpolation_fields(gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, &
                                           qs=qs, qas=qas, fint=fint, er=er)
   endsubroutine face_interpolation_fields

   pure subroutine face_states(gamma, ch, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Face-states adapter of the shared face kernel: the MHD states, `ch` unused.
   real(R8P), intent(in)  :: gamma             !< Specific heats ratio.
   real(R8P), intent(in)  :: ch                !< GLM cleaning speed.
   logical,   intent(in)  :: is_characteristic !< Characteristic (or primitive) variables.
   real(R8P), intent(in)  :: er(NV_K,NV_K)     !< Right eigenvectors.
   real(R8P), intent(in)  :: vr(2,NV_K)        !< Interpolated fields.
   real(R8P), intent(in)  :: q0(NV_K)          !< Conservative variables of cell 0.
   real(R8P), intent(in)  :: q1(NV_K)          !< Conservative variables of cell 1.
   real(R8P), intent(out) :: qL(NV_K)          !< Left state.
   real(R8P), intent(out) :: qR(NV_K)          !< Right state.
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_face_states(gamma=gamma, is_characteristic=is_characteristic, er=er, vr=vr, q0=q0, q1=q1, &
                             qL=qL, qR=qR)
   endsubroutine face_states

   pure subroutine riemann_flux(gamma, ch, d, qL, qR, f, fallback)
   !< Riemann-solver adapter of the shared face kernel: MHD LLF, `ch` the GLM cleaning speed, never a fallback.
   real(R8P),    intent(in)  :: gamma    !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch       !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d        !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K) !< Right state.
   real(R8P),    intent(out) :: f(NV_K)  !< Flux.
   logical,      intent(out) :: fallback !< Fallback flag (always false).
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_riemann_llf(ch=ch, gamma=gamma, d=d, qL=qL, qR=qR, f=f)
   fallback = .false.
   endsubroutine riemann_flux

   pure subroutine cell_flux(gamma, ch, d, q, qa, f)
   !< Cell-flux adapter of the shared face kernel: the MHD physical flux, `gamma` unused.
   real(R8P),    intent(in)  :: gamma        !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch           !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_K)      !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX_K) !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_K)      !< Physical flux.
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_flux(ch=ch, d=d, q=q, qa=qa, f=f)
   endsubroutine cell_flux
endmodule adam_flume_fnl_mhd_eglm_llf_kernels
