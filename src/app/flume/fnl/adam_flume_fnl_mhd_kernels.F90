!< ADAM, FLUME FNL device kernels of the ideal MHD (no divergence control) model.

#include "fundal.H"

module adam_flume_fnl_mhd_kernels
!< ADAM, FLUME FNL device kernels of the ideal MHD (no divergence control) model.
!<
!< The model-agnostic auxiliary-variables kernel (`adam_flume_fnl_aux_kernels_agnostic.INC`) and the MHD dt and
!< conservation kernels (`adam_flume_fnl_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD`,
!< `NV_AUX_K = NV_AUX_MHD` (issue #41, section 4). The face-flux kernel is the shared body with the MHD
!< split (M2-P3). Same kernel rules as `adam_flume_fnl_kernels` (issue #35, D-11/D-12).

! ADAM classes, libraries, parameters
use :: adam_fdv_operators_library, only : compute_derivative1_fd_centered
! ADAM FNL classes, libraries
use :: adam_fnl_weno_kernels,  only : weno_reconstruct_upwind_dev
! FLUME modules
use :: adam_flume_mhd_library, only : compute_face_flux_back_projection=>mhd_face_flux_back_projection, &
                                      conservative_to_auxiliary=>mhd_conservative_to_auxiliary,         &
                                      mhd_fast_speed, mhd_face_split_fluxes, mhd_sum3
use :: adam_flume_mhd_riemann_library, only : mhd_backbone_flux
use :: adam_flume_parameters,  only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                      IQ_RW, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD, POSITIVITY_LIMITER_EPS, &
                                      POSITIVITY_LIMITER_KAPPA, S_MAX
! third party modules
use :: penf,                   only : I4P, I8P, R8P

implicit none
private
public :: apply_floors_dev
public :: blend_positivity_fluxes_dev
public :: compute_conservation_dev
public :: compute_divb_norms_dev
public :: compute_face_fluxes_dev
public :: compute_lambda_max_dev
public :: compute_positivity_factors_dev
public :: compute_q_aux_dev
public :: count_nonfinite_dev
public :: compute_speed_max_dev

contains
   ! public procedures
#include "adam_flume_fnl_face_kernels_agnostic.INC"

#include "adam_flume_fnl_aux_kernels_agnostic.INC"

#include "adam_flume_fnl_mhd_kernels_agnostic.INC"

#include "adam_flume_fnl_positivity_kernels_agnostic.INC"

   ! private procedures
   pure subroutine face_split_fluxes(gamma, ch, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Split adapter of the shared face kernel (issue #41, M2-P3): the MHD split without cleaning, `ch` unused.
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

   call mhd_face_split_fluxes(gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                              fsplit=fsplit, er=er)
   endsubroutine face_split_fluxes

   pure function cleaning_energy(q) result(e)
   !< Cleaning-energy adapter of the floors kernel: zero, `psi` (if any) is not part of the energy of this model.
   real(R8P), intent(in) :: q(NV_K) !< Conservative variables.
   real(R8P)             :: e       !< Cleaning energy.
   !$acc routine seq
   !$omp declare target

   e = 0._R8P
   endfunction cleaning_energy

   pure subroutine backbone_face_flux(gamma, ch, d, qL, qR, f)
   !< Backbone adapter of the positivity limiter: the Lax-Friedrichs flux with the Wu speed, `ch` unused.
   real(R8P),    intent(in)  :: gamma    !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch       !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d        !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K) !< Right state.
   real(R8P),    intent(out) :: f(NV_K)  !< Flux.
   !$acc routine seq
   !$omp declare target

   call mhd_backbone_flux(gamma=gamma, d=d, qL=qL, qR=qR, f=f)
   endsubroutine backbone_face_flux

   pure function internal_energy(q) result(e)
   !< Internal-energy adapter of the positivity limiter: `E` minus the kinetic and magnetic energies, per unit volume.
   real(R8P), intent(in) :: q(NV_K) !< Conservative variables.
   real(R8P)             :: e       !< Internal energy per unit volume.
   !$acc routine seq
   !$omp declare target

   e = q(IQ_RE) - 0.5_R8P * mhd_sum3(q(IQ_RU)**2, q(IQ_RV)**2, q(IQ_RW)**2) / q(IQ_R) - &
       0.5_R8P * mhd_sum3(q(IQ_BX)**2, q(IQ_BY)**2, q(IQ_BZ)**2)
   endfunction internal_energy

   pure subroutine cell_sources(hs, damping, ds, w, qsx, qsy, qsz, qac, s_lo, s_hi)
   !< Sources adapter of the positivity limiter: none (MHD without divergence control).
   integer(I4P), intent(in)  :: hs                     !< Half stencil of the high-order sources.
   real(R8P),    intent(in)  :: damping                !< GLM damping rate.
   real(R8P),    intent(in)  :: ds(3)                  !< Block space steps.
   real(R8P),    intent(in)  :: w(3)                   !< Direction weights: 1 active, 0 null.
   real(R8P),    intent(in)  :: qsx(NV_K,-S_MAX:S_MAX) !< Stencil along x.
   real(R8P),    intent(in)  :: qsy(NV_K,-S_MAX:S_MAX) !< Stencil along y.
   real(R8P),    intent(in)  :: qsz(NV_K,-S_MAX:S_MAX) !< Stencil along z.
   real(R8P),    intent(in)  :: qac(NV_AUX_K)          !< Cell auxiliary variables.
   real(R8P),    intent(out) :: s_lo(NV_K)             !< Backbone sources.
   real(R8P),    intent(out) :: s_hi(NV_K)             !< High-order sources.
   !$acc routine seq
   !$omp declare target

   s_lo = 0._R8P
   s_hi = 0._R8P
   endsubroutine cell_sources
endmodule adam_flume_fnl_mhd_kernels
