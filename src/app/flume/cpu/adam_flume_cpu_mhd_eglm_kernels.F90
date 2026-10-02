!< ADAM, FLUME CPU kernels of the ideal MHD with EGLM cleaning model.

module adam_flume_cpu_mhd_eglm_kernels
!< ADAM, FLUME CPU kernels of the ideal MHD with EGLM cleaning model.
!<
!< The model-agnostic auxiliary-variables loop (`adam_flume_cpu_aux_kernels_agnostic.INC`) and the MHD dt loop
!< (`adam_flume_cpu_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD_EGLM`, `NV_AUX_K = NV_AUX_MHD`; the
!< face-flux kernel is the shared body with the EGLM split (Derigs et al. 2018; issue #47, D-8, section 3.3). The
!< energy includes `psi^2 / 2` (the auxiliary variables and the floors read the thermal pressure from the energy
!< without it); the nonconservative sources and the damping are local residual terms.

! ADAM classes, libraries, parameters
use :: adam_fdv_operators_library, only : compute_derivative1_fd_centered
use :: adam_weno_object,           only : weno_object, weno_reconstruct_upwind
! FLUME modules
use :: adam_flume_mhd_library, only : compute_face_flux_back_projection=>mhd_glm_face_flux_back_projection, &
                                      conservative_to_auxiliary=>mhd_eglm_conservative_to_auxiliary,        &
                                      mhd_eglm_face_split_fluxes, mhd_fast_speed, mhd_sum3
use :: adam_flume_mhd_riemann_library, only : mhd_eglm_backbone_flux
use :: adam_flume_parameters,  only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                      IQ_RW, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_EGLM, POSITIVITY_LIMITER_EPS, &
                                      POSITIVITY_LIMITER_KAPPA, S_MAX
! third party modules
use :: penf,                   only : I4P, I8P, R8P

implicit none
private
public :: add_eglm_sources
public :: add_eglm_sources_limited
public :: add_glm_damping
public :: apply_floors
public :: blend_positivity_fluxes
public :: compute_divb_norms
public :: compute_face_fluxes
public :: compute_lambda_max
public :: compute_positivity_factors
public :: compute_q_aux
public :: count_nonfinite
public :: compute_speed_max

contains
   ! public procedures
#include "adam_flume_cpu_face_kernels_agnostic.INC"

#include "adam_flume_cpu_aux_kernels_agnostic.INC"

#include "adam_flume_cpu_mhd_kernels_agnostic.INC"

#include "adam_flume_cpu_positivity_kernels_agnostic.INC"

   subroutine add_eglm_sources(ni, nj, nk, ngc, blocks_number, hs, dxyz, is_null, q, q_aux, dq)
   !< Add the nonconservative EGLM sources to the residuals of the interior cells (Derigs et al. 2018, eqs. 3.16-3.18;
   !< issue #47, section 3.3): `-(div B) (0, B, u.B, u, 0) - (u.grad psi) (0, 0, psi, 0, 1)` in the order
   !< `(rho, rho u, E, B, psi)`. `div B` and `grad psi` are the centred finite differences of the library, half stencil
   !< `hs` (order `2 hs`, the WENO order: second-order sources drop the Alfven wave to order 2.26, P0), null directions
   !< weighted zero. Local sources in every Runge-Kutta stage, never in the fluxes: momentum, energy and B are conserved
   !< up to these O(div B) terms. The 3-term sums are order-independent (`mhd_sum3`). The caller has refreshed the ghost
   !< cells and the auxiliary variables of `q`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   integer(I4P), intent(in)    :: hs                                 !< Finite difference half stencil, <= ngc.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                        !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(inout) :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Residuals.
   real(R8P)                   :: wx, wy, wz                         !< Direction weights: 1 active, 0 null.
   real(R8P)                   :: dbx, dby, dbz                      !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                   :: dpx, dpy, dpz                      !< Derivatives of psi.
   real(R8P)                   :: divb                               !< div B.
   real(R8P)                   :: upsi                               !< u.grad psi.
   real(R8P)                   :: ub                                 !< u.B.
   integer(I4P)                :: b, i, j, k, c                      !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q, q_aux, dq)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(1,b), q=q(IQ_BX,i-hs:i+hs,j,k,b), dq_ds=dbx)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(2,b), q=q(IQ_BY,i,j-hs:j+hs,k,b), dq_ds=dby)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(3,b), q=q(IQ_BZ,i,j,k-hs:k+hs,b), dq_ds=dbz)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(1,b), q=q(IQ_PSI,i-hs:i+hs,j,k,b), dq_ds=dpx)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(2,b), q=q(IQ_PSI,i,j-hs:j+hs,k,b), dq_ds=dpy)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(3,b), q=q(IQ_PSI,i,j,k-hs:k+hs,b), dq_ds=dpz)
               divb = mhd_sum3(wx * dbx, wy * dby, wz * dbz)
               upsi = mhd_sum3(wx * q_aux(IA_U,i,j,k,b) * dpx, wy * q_aux(IA_V,i,j,k,b) * dpy, &
                               wz * q_aux(IA_W,i,j,k,b) * dpz)
               ub   = mhd_sum3(q_aux(IA_U,i,j,k,b) * q(IQ_BX,i,j,k,b), q_aux(IA_V,i,j,k,b) * q(IQ_BY,i,j,k,b), &
                               q_aux(IA_W,i,j,k,b) * q(IQ_BZ,i,j,k,b))
               do c=1, 3
                  dq(IQ_RU+c-1,i,j,k,b) = dq(IQ_RU+c-1,i,j,k,b) - divb * q(IQ_BX+c-1,i,j,k,b)
                  dq(IQ_BX+c-1,i,j,k,b) = dq(IQ_BX+c-1,i,j,k,b) - divb * q_aux(IA_U+c-1,i,j,k,b)
               enddo
               dq(IQ_RE, i,j,k,b) = dq(IQ_RE, i,j,k,b) - divb * ub - upsi * q(IQ_PSI,i,j,k,b)
               dq(IQ_PSI,i,j,k,b) = dq(IQ_PSI,i,j,k,b) - upsi
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_eglm_sources

   subroutine add_eglm_sources_limited(ni, nj, nk, ngc, blocks_number, hs, dxyz, is_null, q, q_aux, lam, dq)
   !< Add the EGLM sources with the positivity limiter (issue #47, D-9): in a cell with factor `Lambda < 1` the
   !< second-order sources plus `Lambda` times the difference to the order-2 hs ones (the sources are linear in `div B`
   !< and `u.grad psi`, so the two scalars are blended); with `Lambda = 1` the arithmetic of `add_eglm_sources`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   integer(I4P), intent(in)    :: hs                                 !< Finite difference half stencil, <= ngc.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                        !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: lam(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Cell factors (component 1).
   real(R8P),    intent(inout) :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Residuals.
   real(R8P)                   :: divb, upsi                         !< div B, u.grad psi (high order or blended).
   real(R8P)                   :: divb_lo, upsi_lo                   !< Second-order div B, u.grad psi.
   real(R8P)                   :: ub                                 !< u.B.
   real(R8P)                   :: lc                                 !< Cell factor.
   integer(I4P)                :: b, i, j, k, c                      !< Counters.

   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q, q_aux, lam, dq)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call eglm_source_scalars(ngc=ngc, hs=hs, dxyz=dxyz(:,b), is_null=is_null, q=q, q_aux=q_aux, i=i, j=j, &
                                        k=k, b=b, divb=divb, upsi=upsi)
               lc = lam(1,i,j,k,b)
               if (lc < 1._R8P) then
                  call eglm_source_scalars(ngc=ngc, hs=1_I4P, dxyz=dxyz(:,b), is_null=is_null, q=q, q_aux=q_aux, i=i, &
                                           j=j, k=k, b=b, divb=divb_lo, upsi=upsi_lo)
                  divb = divb_lo + lc * (divb - divb_lo)
                  upsi = upsi_lo + lc * (upsi - upsi_lo)
               endif
               ub   = mhd_sum3(q_aux(IA_U,i,j,k,b) * q(IQ_BX,i,j,k,b), q_aux(IA_V,i,j,k,b) * q(IQ_BY,i,j,k,b), &
                               q_aux(IA_W,i,j,k,b) * q(IQ_BZ,i,j,k,b))
               do c=1, 3
                  dq(IQ_RU+c-1,i,j,k,b) = dq(IQ_RU+c-1,i,j,k,b) - divb * q(IQ_BX+c-1,i,j,k,b)
                  dq(IQ_BX+c-1,i,j,k,b) = dq(IQ_BX+c-1,i,j,k,b) - divb * q_aux(IA_U+c-1,i,j,k,b)
               enddo
               dq(IQ_RE, i,j,k,b) = dq(IQ_RE, i,j,k,b) - divb * ub - upsi * q(IQ_PSI,i,j,k,b)
               dq(IQ_PSI,i,j,k,b) = dq(IQ_PSI,i,j,k,b) - upsi
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_eglm_sources_limited

   subroutine add_glm_damping(ni, nj, nk, ngc, blocks_number, damping, q, dq)
   !< Add the damping source to the residuals of the interior cells, `dq(psi) = dq(psi) - (alpha c_h / L) psi`, as GLM
   !< (issue #47, section 3.3). Only `psi` is damped: the total energy is unchanged, so the removed cleaning energy
   !< `psi^2 / 2` becomes thermal energy (the pressure grows) and the energy drift stays bounded by the div(B) sources.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                  !< Actual blocks number.
   real(R8P),    intent(in)    :: damping                        !< Damping rate alpha c_h / L.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(inout) :: dq(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   integer(I4P)                :: b, i, j, k                     !< Counters.

   !$omp parallel do collapse(4) default(firstprivate) shared(q, dq)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               dq(IQ_PSI,i,j,k,b) = dq(IQ_PSI,i,j,k,b) - damping * q(IQ_PSI,i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_glm_damping

   ! private procedures
   pure subroutine face_split_fluxes(gamma, ch, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Split adapter of the shared face kernel: the MHD split with EGLM cleaning.
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

   call mhd_eglm_face_split_fluxes(ch=ch, gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                                   fsplit=fsplit, er=er)
   endsubroutine face_split_fluxes

   pure function cleaning_energy(q) result(e)
   !< Cleaning-energy adapter of the floors kernel: `psi^2 / 2`, part of the EGLM energy.
   real(R8P), intent(in) :: q(NV_K) !< Conservative variables.
   real(R8P)             :: e       !< Cleaning energy.
   !$acc routine seq
   !$omp declare target

   e = 0.5_R8P * q(IQ_PSI)**2
   endfunction cleaning_energy

   pure subroutine backbone_face_flux(gamma, ch, d, qL, qR, f)
   !< Backbone adapter of the positivity limiter: the Lax-Friedrichs flux with the Wu speed, at least `ch`.
   real(R8P),    intent(in)  :: gamma    !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch       !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d        !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K) !< Right state.
   real(R8P),    intent(out) :: f(NV_K)  !< Flux.
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_backbone_flux(ch=ch, gamma=gamma, d=d, qL=qL, qR=qR, f=f)
   endsubroutine backbone_face_flux

   pure function internal_energy(q) result(e)
   !< Internal-energy adapter of the positivity limiter: `E` minus the kinetic, magnetic and cleaning energies, per unit volume.
   real(R8P), intent(in) :: q(NV_K) !< Conservative variables.
   real(R8P)             :: e       !< Internal energy per unit volume.
   !$acc routine seq
   !$omp declare target

   e = q(IQ_RE) - 0.5_R8P * mhd_sum3(q(IQ_RU)**2, q(IQ_RV)**2, q(IQ_RW)**2) / q(IQ_R) - &
       0.5_R8P * mhd_sum3(q(IQ_BX)**2, q(IQ_BY)**2, q(IQ_BZ)**2) - 0.5_R8P * q(IQ_PSI)**2
   endfunction internal_energy

   pure subroutine cell_sources(ngc, hs, damping, dxyz, is_null, q, q_aux, i, j, k, b, s_lo, s_hi)
   !< Sources adapter of the positivity limiter (EGLM): the nonconservative sources with second-order (`s_lo`, the
   !< backbone's) and order-2 hs (`s_hi`) centred differences, both with the damping of `psi`.
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
   real(R8P)                 :: divb, upsi                        !< div B, u.grad psi.
   !$acc routine seq
   !$omp declare target

   call eglm_source_scalars(ngc=ngc, hs=1_I4P, dxyz=dxyz, is_null=is_null, q=q, q_aux=q_aux, i=i, j=j, k=k, b=b, &
                            divb=divb, upsi=upsi)
   call eglm_source_vector(q=q(1:NV_K,i,j,k,b), qa=q_aux(1:NV_AUX_K,i,j,k,b), divb=divb, upsi=upsi, s=s_lo)
   call eglm_source_scalars(ngc=ngc, hs=hs, dxyz=dxyz, is_null=is_null, q=q, q_aux=q_aux, i=i, j=j, k=k, b=b, &
                            divb=divb, upsi=upsi)
   call eglm_source_vector(q=q(1:NV_K,i,j,k,b), qa=q_aux(1:NV_AUX_K,i,j,k,b), divb=divb, upsi=upsi, s=s_hi)
   s_lo(IQ_PSI) = s_lo(IQ_PSI) - damping * q(IQ_PSI,i,j,k,b)
   s_hi(IQ_PSI) = s_hi(IQ_PSI) - damping * q(IQ_PSI,i,j,k,b)
   endsubroutine cell_sources

   pure subroutine eglm_source_scalars(ngc, hs, dxyz, is_null, q, q_aux, i, j, k, b, divb, upsi)
   !< Return `div B` and `u.grad psi` of a cell by the centred differences of half stencil `hs` (the arithmetic of
   !< `add_eglm_sources`).
   integer(I4P), intent(in)  :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)  :: hs                                !< Half stencil.
   real(R8P),    intent(in)  :: dxyz(3)                           !< Block space steps.
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),    intent(in)  :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P), intent(in)  :: i, j, k, b                        !< Cell indexes.
   real(R8P),    intent(out) :: divb, upsi                        !< div B, u.grad psi.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: dbx, dby, dbz                     !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                 :: dpx, dpy, dpz                     !< Derivatives of psi.
   !$acc routine seq
   !$omp declare target

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(1), q=q(IQ_BX,i-hs:i+hs,j,k,b), dq_ds=dbx)
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(2), q=q(IQ_BY,i,j-hs:j+hs,k,b), dq_ds=dby)
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(3), q=q(IQ_BZ,i,j,k-hs:k+hs,b), dq_ds=dbz)
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(1), q=q(IQ_PSI,i-hs:i+hs,j,k,b), dq_ds=dpx)
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(2), q=q(IQ_PSI,i,j-hs:j+hs,k,b), dq_ds=dpy)
   call compute_derivative1_fd_centered(s=hs, ds=dxyz(3), q=q(IQ_PSI,i,j,k-hs:k+hs,b), dq_ds=dpz)
   divb = mhd_sum3(wx * dbx, wy * dby, wz * dbz)
   upsi = mhd_sum3(wx * q_aux(IA_U,i,j,k,b) * dpx, wy * q_aux(IA_V,i,j,k,b) * dpy, wz * q_aux(IA_W,i,j,k,b) * dpz)
   endsubroutine eglm_source_scalars

   pure subroutine eglm_source_vector(q, qa, divb, upsi, s)
   !< Return the EGLM source vector of a cell from its `div B` and `u.grad psi` (linear in both).
   real(R8P), intent(in)  :: q(NV_K)      !< Conservative variables.
   real(R8P), intent(in)  :: qa(NV_AUX_K) !< Auxiliary variables.
   real(R8P), intent(in)  :: divb, upsi   !< div B, u.grad psi.
   real(R8P), intent(out) :: s(NV_K)      !< Sources.
   real(R8P)              :: ub           !< u.B.
   integer(I4P)           :: c            !< Counter.
   !$acc routine seq
   !$omp declare target

   ub = mhd_sum3(qa(IA_U) * q(IQ_BX), qa(IA_V) * q(IQ_BY), qa(IA_W) * q(IQ_BZ))
   s = 0._R8P
   do c=1, 3
      s(IQ_RU+c-1) = -divb * q(IQ_BX+c-1)
      s(IQ_BX+c-1) = -divb * qa(IA_U+c-1)
   enddo
   s(IQ_RE)  = -divb * ub - upsi * q(IQ_PSI)
   s(IQ_PSI) = -upsi
   endsubroutine eglm_source_vector
endmodule adam_flume_cpu_mhd_eglm_kernels
