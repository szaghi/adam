!< ADAM, FLUME CPU kernels of the Euler model.

module adam_flume_cpu_euler_kernels
!< ADAM, FLUME CPU kernels of the Euler model.
!<
!< The model-agnostic loop bodies (`adam_flume_cpu_{face,aux}_kernels_agnostic.INC`) instantiated on the Euler physics: the
!< local arrays are sized by the Euler constants `NV_K = NV_EULER`, `NV_AUX_K = NV_AUX` (issue #41, section 4), plus
!< the Euler signal-speed loop.

! ADAM classes, libraries, parameters
use :: adam_weno_object,         only : weno_object, weno_reconstruct_upwind
! FLUME modules
use :: adam_flume_dissipation_library, only : dissipative_diffusivity
use :: adam_flume_euler_library, only : compute_face_flux_back_projection, compute_face_split_fluxes,                   &
                                        compute_riemann_llf, conservative_to_auxiliary
use :: adam_flume_parameters,    only : IA_A, IA_R, IA_T, IA_U, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX_K=>NV_AUX, &
                                        NV_K=>NV_EULER, POSITIVITY_LIMITER_KAPPA, S_MAX
! third party modules
use :: penf,                     only : I4P, I8P, R8P

implicit none
private
public :: blend_inadmissible_ghosts
public :: blend_positivity_fluxes
public :: compute_face_fluxes
public :: compute_lambda_max
public :: compute_lambda_max_dissipative
public :: compute_backbone_fluxes
public :: compute_positivity_factors
public :: compute_seam_positivity_factors
public :: compute_q_aux
public :: count_nonfinite

contains
   ! public procedures
#include "adam_flume_cpu_face_kernels_agnostic.INC"

#include "adam_flume_cpu_aux_kernels_agnostic.INC"

#include "adam_flume_cpu_positivity_kernels_agnostic.INC"

#include "adam_flume_cpu_ghost_kernels_agnostic.INC"

   subroutine compute_lambda_max(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz, is_null, q, lambda_max)
   !< Compute `max(sum_d (|u_d| + a) / dx_d)` over the interior cells (null directions excluded).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                         !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                             !< Gas constant.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                   !< Blocks space steps [3, nb].
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                    !< Maximum of sum_d (|u_d| + a) / dx_d.
   real(R8P)                 :: qa(NV_AUX_K)                  !< Auxiliary variables of one cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.
   integer(I4P)              :: d                             !< Direction counter.

   lambda_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q) reduction(max:lambda_max)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               lambda_max = max(lambda_max, sum([((abs(qa(IA_U+d-1)) + qa(IA_A)) / dxyz(d,b), d=1, 3)], mask=.not.is_null))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_lambda_max

   subroutine compute_lambda_max_dissipative(ni, nj, nk, ngc, blocks_number, gamma, R, cp, mu0, k0, eta, t_wall, tref, &
                                             omega_mu, omega_k, dxyz, is_null, q, lambda_max, lambda_hyp, lambda_dif,  &
                                             re_cell_min)
   !< Compute `max(sum_d (|u_d| + a) / dx_d + 2 nu sum_d 1 / dx_d^2)` over the interior cells (issue #65, D-M4-3), `nu`
   !< the largest diffusivity of the cell, with the maxima of the two parts and the minimum cell Reynolds number
   !< `sum_d (|u_d| + a) / dx_d / (nu sum_d 1 / dx_d^2)`, `(|u| + a) dx / nu` on a uniform 1D grid.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma, R, cp                  !< Specific heats ratio, gas constant, cp.
   real(R8P),    intent(in)  :: mu0, k0, eta                  !< Dissipative coefficients.
   real(R8P),    intent(in)  :: t_wall                        !< Hottest isothermal wall temperature, 0 if none.
   real(R8P),    intent(in)  :: tref, omega_mu, omega_k       !< Temperature laws.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                   !< Blocks space steps [3, nb].
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                    !< Maximum of the sum.
   real(R8P),    intent(out) :: lambda_hyp                    !< Maximum of the hyperbolic part.
   real(R8P),    intent(out) :: lambda_dif                    !< Maximum of the diffusive part.
   real(R8P),    intent(out) :: re_cell_min                   !< Minimum cell Reynolds number.
   real(R8P)                 :: qa(NV_AUX_K)                  !< Auxiliary variables of one cell.
   real(R8P)                 :: lh, ld                        !< Hyperbolic and diffusive parts of one cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.
   integer(I4P)              :: d                             !< Direction counter.

   lambda_max = 0._R8P ; lambda_hyp = 0._R8P ; lambda_dif = 0._R8P ; re_cell_min = huge(1._R8P)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q) &
   !$omp& reduction(max:lambda_max,lambda_hyp,lambda_dif) reduction(min:re_cell_min)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               lh = sum([((abs(qa(IA_U+d-1)) + qa(IA_A)) / dxyz(d,b), d=1, 3)], mask=.not.is_null)
               ld = 2._R8P * dissipative_diffusivity(rho=qa(IA_R), T=qa(IA_T), t_wall=t_wall, gamma=gamma, cp=cp, mu0=mu0, k0=k0, &
                                                     eta=eta, tref=tref, omega_mu=omega_mu, omega_k=omega_k) *     &
                    sum([(1._R8P / dxyz(d,b)**2, d=1, 3)], mask=.not.is_null)
               lambda_max = max(lambda_max, lh + ld)
               lambda_hyp = max(lambda_hyp, lh)
               lambda_dif = max(lambda_dif, ld)
               if (ld > 0._R8P) re_cell_min = min(re_cell_min, 2._R8P * lh / ld)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_lambda_max_dissipative

   ! private procedures
   pure subroutine face_split_fluxes(gamma, ch, d, S, is_characteristic, qs, qas, fsplit, er, mu)
   !< Split adapter of the shared face kernel (issue #41, M2-P3): the Euler split, `ch` unused.
   real(R8P),    intent(in)  :: gamma                          !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch                             !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d                              !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                              !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic              !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_K,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_K,1-S_MAX:S_MAX)    !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_K) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_K,NV_K)                  !< Right eigenvectors.
   real(R8P),    intent(out) :: mu(NV_K)                       !< Magnitude of each field (WENO descaler).
   !$acc routine seq
   !$omp declare target

   call compute_face_split_fluxes(gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                                  fsplit=fsplit, er=er, mu=mu)
   endsubroutine face_split_fluxes

   pure subroutine backbone_face_flux(gamma, ch, d, qL, qR, f)
   !< Backbone adapter of the positivity limiter: the Rusanov flux, `ch` unused.
   real(R8P),    intent(in)  :: gamma    !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch       !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d        !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K) !< Right state.
   real(R8P),    intent(out) :: f(NV_K)  !< Flux.
   !$acc routine seq
   !$omp declare target

   call compute_riemann_llf(gamma=gamma, d=d, qL=qL, qR=qR, f=f)
   endsubroutine backbone_face_flux

   pure function internal_energy(q) result(e)
   !< Internal-energy adapter of the positivity limiter: `E` minus the kinetic energy, per unit volume.
   real(R8P), intent(in) :: q(NV_K) !< Conservative variables.
   real(R8P)             :: e       !< Internal energy per unit volume.
   !$acc routine seq
   !$omp declare target

   e = q(IQ_RE) - 0.5_R8P * (q(IQ_RU)**2 + q(IQ_RV)**2 + q(IQ_RW)**2) / q(IQ_R)
   endfunction internal_energy

   pure subroutine cell_sources(ngc, hs, damping, dxyz, is_null, q, q_aux, i, j, k, b, s_lo, s_hi)
   !< Sources adapter of the positivity limiter: none (Euler).
   integer(I4P), intent(in)  :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)  :: hs                                !< Half stencil of the high-order sources.
   real(R8P),    intent(in)  :: damping                           !< GLM damping rate.
   real(R8P),    intent(in)  :: dxyz(3)                           !< Block space steps.
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),    intent(in)  :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P), intent(in)  :: i, j, k, b                        !< Cell indexes.
   real(R8P),    intent(out) :: s_lo(NV_K)                        !< Backbone sources.
   real(R8P),    intent(out) :: s_hi(NV_K)                        !< High-order sources.
   !$acc routine seq
   !$omp declare target

   s_lo = 0._R8P
   s_hi = 0._R8P
   endsubroutine cell_sources
endmodule adam_flume_cpu_euler_kernels
