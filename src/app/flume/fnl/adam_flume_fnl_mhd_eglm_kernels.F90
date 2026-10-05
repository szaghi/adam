!< ADAM, FLUME FNL device kernels of the ideal MHD with EGLM cleaning model.

#include "fundal.H"

module adam_flume_fnl_mhd_eglm_kernels
!< ADAM, FLUME FNL device kernels of the ideal MHD with EGLM cleaning model.
!<
!< The model-agnostic auxiliary-variables kernel (`adam_flume_fnl_aux_kernels_agnostic.INC`) and the MHD dt and
!< conservation kernels (`adam_flume_fnl_mhd_kernels_agnostic.INC`) instantiated with `NV_K = NV_MHD_EGLM`,
!< `NV_AUX_K = NV_AUX_MHD`; the face-flux kernel is the shared body with the EGLM split (Derigs et al. 2018; issue #47,
!< D-8, section 3.3). Device twin of `adam_flume_cpu_mhd_eglm_kernels`. Same kernel rules as `adam_flume_fnl_kernels`
!< (issue #35, D-11/D-12).

! ADAM classes, libraries, parameters
use :: adam_fdv_operators_library, only : compute_derivative1_fd_centered
! ADAM FNL classes, libraries
use :: adam_fnl_weno_kernels,  only : weno_reconstruct_upwind_dev
! FLUME modules
use :: adam_flume_mhd_library, only : compute_face_flux_back_projection=>mhd_glm_face_flux_back_projection, &
                                      conservative_to_auxiliary=>mhd_eglm_conservative_to_auxiliary,        &
                                      mhd_eglm_face_split_fluxes, mhd_fast_speed, mhd_sum3
use :: adam_flume_mhd_riemann_library, only : mhd_eglm_backbone_flux
use :: adam_flume_parameters,  only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                      IQ_RW, NV_AUX_K=>NV_AUX_MHD, NV_K=>NV_MHD_EGLM,                         &
                                      POSITIVITY_LIMITER_KAPPA, S_MAX
! third party modules
use :: penf,                   only : I4P, I8P, R8P

implicit none
private
public :: add_eglm_sources_dev
public :: add_eglm_sources_limited_dev
public :: add_glm_damping_dev
public :: apply_floors_dev
public :: blend_inadmissible_ghosts_dev
public :: blend_positivity_fluxes_dev
public :: compute_conservation_dev
public :: compute_divb_norms_dev
public :: compute_face_fluxes_dev
public :: compute_lambda_max_dev
public :: compute_backbone_fluxes_host
public :: compute_positivity_factors_dev
public :: compute_seam_positivity_factors_host
public :: compute_q_aux_dev
public :: count_nonfinite_dev
public :: compute_speed_max_dev

contains
   ! public procedures
#include "adam_flume_fnl_face_kernels_agnostic.INC"

#include "adam_flume_fnl_aux_kernels_agnostic.INC"

#include "adam_flume_fnl_mhd_kernels_agnostic.INC"

#include "adam_flume_fnl_positivity_kernels_agnostic.INC"

#include "adam_flume_fnl_ghost_kernels_agnostic.INC"

   subroutine add_eglm_sources_dev(ni, nj, nk, ngc, blocks_number, hs, dxyz_gpu, is_null, q_gpu, q_aux_gpu, dq_gpu)
   !< Add the nonconservative EGLM sources to the residuals of the interior cells (Derigs et al. 2018, eqs. 3.16-3.18;
   !< issue #47, section 3.3): `-(div B) (0, B, u.B, u, 0) - (u.grad psi) (0, 0, psi, 0, 1)` in the order
   !< `(rho, rho u, E, B, psi)`. `div B` and `grad psi` are the centred finite differences of the library, half stencil
   !< `hs <= S_MAX` (order `2 hs`, the WENO order), null directions weighted zero; the 3-term sums are
   !< order-independent (`mhd_sum3`). Device twin of the CPU `add_eglm_sources`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                        !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                          !< Actual blocks number.
   integer(I4P), intent(in)    :: hs                                     !< Finite difference half stencil.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                        !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                             !< Null directions.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Residuals.
   real(R8P)                   :: sx(1-S_MAX:1+S_MAX)                    !< Private stencil along x.
   real(R8P)                   :: sy(1-S_MAX:1+S_MAX)                    !< Private stencil along y.
   real(R8P)                   :: sz(1-S_MAX:1+S_MAX)                    !< Private stencil along z.
   real(R8P)                   :: wx, wy, wz                             !< Direction weights: 1 active, 0 null.
   real(R8P)                   :: dbx, dby, dbz                          !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                   :: dpx, dpy, dpz                          !< Derivatives of psi.
   real(R8P)                   :: divb                                   !< div B.
   real(R8P)                   :: upsi                                   !< u.grad psi.
   real(R8P)                   :: ub                                     !< u.B.
   integer(I4P)                :: b, i, j, k, m, c                       !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,q_aux_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,hs,wx,wy,wz) private(sx,sy,sz,dbx,dby,dbz,dpx,dpy,dpz,divb,upsi,ub)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,q_aux_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,hs,wx,wy,wz) private(sx,sy,sz,dbx,dby,dbz,dpx,dpy,dpz,divb,upsi,ub)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=-hs, hs
         sx(1+m) = q_gpu(b,i+m,j,k,IQ_BX)
         sy(1+m) = q_gpu(b,i,j+m,k,IQ_BY)
         sz(1+m) = q_gpu(b,i,j,k+m,IQ_BZ)
      enddo
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,1), q=sx(1-hs:1+hs), dq_ds=dbx)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,2), q=sy(1-hs:1+hs), dq_ds=dby)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,3), q=sz(1-hs:1+hs), dq_ds=dbz)
      !$acc loop seq
      do m=-hs, hs
         sx(1+m) = q_gpu(b,i+m,j,k,IQ_PSI)
         sy(1+m) = q_gpu(b,i,j+m,k,IQ_PSI)
         sz(1+m) = q_gpu(b,i,j,k+m,IQ_PSI)
      enddo
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,1), q=sx(1-hs:1+hs), dq_ds=dpx)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,2), q=sy(1-hs:1+hs), dq_ds=dpy)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,3), q=sz(1-hs:1+hs), dq_ds=dpz)
      divb = mhd_sum3(wx * dbx, wy * dby, wz * dbz)
      upsi = mhd_sum3(wx * q_aux_gpu(b,i,j,k,IA_U) * dpx, wy * q_aux_gpu(b,i,j,k,IA_V) * dpy, &
                      wz * q_aux_gpu(b,i,j,k,IA_W) * dpz)
      ub   = mhd_sum3(q_aux_gpu(b,i,j,k,IA_U) * q_gpu(b,i,j,k,IQ_BX), q_aux_gpu(b,i,j,k,IA_V) * q_gpu(b,i,j,k,IQ_BY), &
                      q_aux_gpu(b,i,j,k,IA_W) * q_gpu(b,i,j,k,IQ_BZ))
      !$acc loop seq
      do c=1, 3
         dq_gpu(b,i,j,k,IQ_RU+c-1) = dq_gpu(b,i,j,k,IQ_RU+c-1) - divb * q_gpu(b,i,j,k,IQ_BX+c-1)
         dq_gpu(b,i,j,k,IQ_BX+c-1) = dq_gpu(b,i,j,k,IQ_BX+c-1) - divb * q_aux_gpu(b,i,j,k,IA_U+c-1)
      enddo
      dq_gpu(b,i,j,k,IQ_RE)  = dq_gpu(b,i,j,k,IQ_RE)  - divb * ub - upsi * q_gpu(b,i,j,k,IQ_PSI)
      dq_gpu(b,i,j,k,IQ_PSI) = dq_gpu(b,i,j,k,IQ_PSI) - upsi
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_eglm_sources_dev

   subroutine add_eglm_sources_limited_dev(ni, nj, nk, ngc, blocks_number, hs, dxyz_gpu, is_null, q_gpu, q_aux_gpu, &
                                           lam_gpu, dq_gpu)
   !< Add the EGLM sources with the positivity limiter (issue #47, D-9): in a cell with factor `Lambda < 1` the
   !< second-order sources plus `Lambda` times the difference to the order-2 hs ones, blended through `div B` and
   !< `u.grad psi`; with `Lambda = 1` the arithmetic of `add_eglm_sources_dev` (device twin of the CPU
   !< `add_eglm_sources_limited`).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                        !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                          !< Actual blocks number.
   integer(I4P), intent(in)    :: hs                                     !< Finite difference half stencil.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                        !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                             !< Null directions.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: lam_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Cell factors (component 1).
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Residuals.
   real(R8P)                   :: sx(1-S_MAX:1+S_MAX)                    !< Private stencil along x.
   real(R8P)                   :: sy(1-S_MAX:1+S_MAX)                    !< Private stencil along y.
   real(R8P)                   :: sz(1-S_MAX:1+S_MAX)                    !< Private stencil along z.
   real(R8P)                   :: wx, wy, wz                             !< Direction weights: 1 active, 0 null.
   real(R8P)                   :: dbx, dby, dbz                          !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                   :: dpx, dpy, dpz                          !< Derivatives of psi.
   real(R8P)                   :: lbx, lby, lbz                          !< Second-order derivatives of B.
   real(R8P)                   :: lpx, lpy, lpz                          !< Second-order derivatives of psi.
   real(R8P)                   :: divb, upsi                             !< div B, u.grad psi (high order or blended).
   real(R8P)                   :: divb_lo, upsi_lo                       !< Second-order div B, u.grad psi.
   real(R8P)                   :: ub                                     !< u.B.
   real(R8P)                   :: lc                                     !< Cell factor.
   integer(I4P)                :: b, i, j, k, m, c                       !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu,q_aux_gpu,lam_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,hs,wx,wy,wz)                                                   &
   !$acc& private(sx,sy,sz,dbx,dby,dbz,dpx,dpy,dpz,lbx,lby,lbz,lpx,lpy,lpz,divb,upsi,divb_lo,upsi_lo,ub,lc)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu,q_aux_gpu,lam_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,hs,wx,wy,wz) &
   !$omp& private(sx,sy,sz,dbx,dby,dbz,dpx,dpy,dpz,lbx,lby,lbz,lpx,lpy,lpz,divb,upsi,divb_lo,upsi_lo,ub,lc)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=-hs, hs
         sx(1+m) = q_gpu(b,i+m,j,k,IQ_BX)
         sy(1+m) = q_gpu(b,i,j+m,k,IQ_BY)
         sz(1+m) = q_gpu(b,i,j,k+m,IQ_BZ)
      enddo
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,1), q=sx(1-hs:1+hs), dq_ds=dbx)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,2), q=sy(1-hs:1+hs), dq_ds=dby)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,3), q=sz(1-hs:1+hs), dq_ds=dbz)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,1), q=sx(0:2), dq_ds=lbx)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,2), q=sy(0:2), dq_ds=lby)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,3), q=sz(0:2), dq_ds=lbz)
      !$acc loop seq
      do m=-hs, hs
         sx(1+m) = q_gpu(b,i+m,j,k,IQ_PSI)
         sy(1+m) = q_gpu(b,i,j+m,k,IQ_PSI)
         sz(1+m) = q_gpu(b,i,j,k+m,IQ_PSI)
      enddo
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,1), q=sx(1-hs:1+hs), dq_ds=dpx)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,2), q=sy(1-hs:1+hs), dq_ds=dpy)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,3), q=sz(1-hs:1+hs), dq_ds=dpz)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,1), q=sx(0:2), dq_ds=lpx)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,2), q=sy(0:2), dq_ds=lpy)
      call compute_derivative1_fd_centered(s=1, ds=dxyz_gpu(b,3), q=sz(0:2), dq_ds=lpz)
      divb = mhd_sum3(wx * dbx, wy * dby, wz * dbz)
      upsi = mhd_sum3(wx * q_aux_gpu(b,i,j,k,IA_U) * dpx, wy * q_aux_gpu(b,i,j,k,IA_V) * dpy, &
                      wz * q_aux_gpu(b,i,j,k,IA_W) * dpz)
      lc = lam_gpu(b,i,j,k,1)
      if (lc < 1._R8P) then
         divb_lo = mhd_sum3(wx * lbx, wy * lby, wz * lbz)
         upsi_lo = mhd_sum3(wx * q_aux_gpu(b,i,j,k,IA_U) * lpx, wy * q_aux_gpu(b,i,j,k,IA_V) * lpy, &
                            wz * q_aux_gpu(b,i,j,k,IA_W) * lpz)
         divb = divb_lo + lc * (divb - divb_lo)
         upsi = upsi_lo + lc * (upsi - upsi_lo)
      endif
      ub   = mhd_sum3(q_aux_gpu(b,i,j,k,IA_U) * q_gpu(b,i,j,k,IQ_BX), q_aux_gpu(b,i,j,k,IA_V) * q_gpu(b,i,j,k,IQ_BY), &
                      q_aux_gpu(b,i,j,k,IA_W) * q_gpu(b,i,j,k,IQ_BZ))
      !$acc loop seq
      do c=1, 3
         dq_gpu(b,i,j,k,IQ_RU+c-1) = dq_gpu(b,i,j,k,IQ_RU+c-1) - divb * q_gpu(b,i,j,k,IQ_BX+c-1)
         dq_gpu(b,i,j,k,IQ_BX+c-1) = dq_gpu(b,i,j,k,IQ_BX+c-1) - divb * q_aux_gpu(b,i,j,k,IA_U+c-1)
      enddo
      dq_gpu(b,i,j,k,IQ_RE)  = dq_gpu(b,i,j,k,IQ_RE)  - divb * ub - upsi * q_gpu(b,i,j,k,IQ_PSI)
      dq_gpu(b,i,j,k,IQ_PSI) = dq_gpu(b,i,j,k,IQ_PSI) - upsi
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_eglm_sources_limited_dev

   subroutine add_glm_damping_dev(ni, nj, nk, ngc, blocks_number, damping, q_gpu, dq_gpu)
   !< Add the damping source to the residuals of the interior cells, `dq(psi) = dq(psi) - (alpha c_h / L) psi`, as GLM
   !< (issue #47, section 3.3). Only `psi` is damped: the total energy is unchanged, so the removed cleaning energy
   !< `psi^2 / 2` becomes thermal energy (the pressure grows) and the energy drift stays bounded by the div(B) sources.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   real(R8P),    intent(in)    :: damping                            !< Damping rate alpha c_h / L.
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
   pure subroutine face_split_fluxes(gamma, ch, d, S, is_characteristic, qs, qas, fsplit, er, mu)
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
   real(R8P),    intent(out) :: mu(NV_K)                       !< Magnitude of each field (WENO descaler).
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_face_split_fluxes(ch=ch, gamma=gamma, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                                   fsplit=fsplit, er=er, mu=mu)
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

   pure subroutine cell_sources(hs, damping, ds, w, qsx, qsy, qsz, qac, s_lo, s_hi)
   !< Sources adapter of the positivity limiter (EGLM): the nonconservative sources with second-order (`s_lo`, the
   !< backbone's) and order-2 hs (`s_hi`) centred differences, both with the damping of `psi` (device twin of the CPU
   !< adapter, from the private stencils).
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
   real(R8P)                 :: divb, upsi             !< div B, u.grad psi.
   !$acc routine seq
   !$omp declare target

   call eglm_stencil_scalars(hs=1_I4P, ds=ds, w=w, qsx=qsx, qsy=qsy, qsz=qsz, qac=qac, divb=divb, upsi=upsi)
   call eglm_source_vector(q=qsx(:,0), qa=qac, divb=divb, upsi=upsi, s=s_lo)
   call eglm_stencil_scalars(hs=hs, ds=ds, w=w, qsx=qsx, qsy=qsy, qsz=qsz, qac=qac, divb=divb, upsi=upsi)
   call eglm_source_vector(q=qsx(:,0), qa=qac, divb=divb, upsi=upsi, s=s_hi)
   s_lo(IQ_PSI) = s_lo(IQ_PSI) - damping * qsx(IQ_PSI,0)
   s_hi(IQ_PSI) = s_hi(IQ_PSI) - damping * qsx(IQ_PSI,0)
   endsubroutine cell_sources

   pure subroutine eglm_stencil_scalars(hs, ds, w, qsx, qsy, qsz, qac, divb, upsi)
   !< Return `div B` and `u.grad psi` of a cell from its private stencils, centred differences of half stencil `hs` (the
   !< arithmetic of `add_eglm_sources_dev`).
   integer(I4P), intent(in)  :: hs                     !< Half stencil.
   real(R8P),    intent(in)  :: ds(3)                  !< Block space steps.
   real(R8P),    intent(in)  :: w(3)                   !< Direction weights: 1 active, 0 null.
   real(R8P),    intent(in)  :: qsx(NV_K,-S_MAX:S_MAX) !< Stencil along x.
   real(R8P),    intent(in)  :: qsy(NV_K,-S_MAX:S_MAX) !< Stencil along y.
   real(R8P),    intent(in)  :: qsz(NV_K,-S_MAX:S_MAX) !< Stencil along z.
   real(R8P),    intent(in)  :: qac(NV_AUX_K)          !< Cell auxiliary variables.
   real(R8P),    intent(out) :: divb, upsi             !< div B, u.grad psi.
   real(R8P)                 :: sx(1-S_MAX:1+S_MAX)    !< Stencil of one field along x.
   real(R8P)                 :: sy(1-S_MAX:1+S_MAX)    !< Stencil of one field along y.
   real(R8P)                 :: sz(1-S_MAX:1+S_MAX)    !< Stencil of one field along z.
   real(R8P)                 :: dbx, dby, dbz          !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                 :: dpx, dpy, dpz          !< Derivatives of psi.
   integer(I4P)              :: m                      !< Counter.
   !$acc routine seq
   !$omp declare target

   do m=-hs, hs
      sx(1+m) = qsx(IQ_BX,m)
      sy(1+m) = qsy(IQ_BY,m)
      sz(1+m) = qsz(IQ_BZ,m)
   enddo
   call compute_derivative1_fd_centered(s=hs, ds=ds(1), q=sx(1-hs:1+hs), dq_ds=dbx)
   call compute_derivative1_fd_centered(s=hs, ds=ds(2), q=sy(1-hs:1+hs), dq_ds=dby)
   call compute_derivative1_fd_centered(s=hs, ds=ds(3), q=sz(1-hs:1+hs), dq_ds=dbz)
   do m=-hs, hs
      sx(1+m) = qsx(IQ_PSI,m)
      sy(1+m) = qsy(IQ_PSI,m)
      sz(1+m) = qsz(IQ_PSI,m)
   enddo
   call compute_derivative1_fd_centered(s=hs, ds=ds(1), q=sx(1-hs:1+hs), dq_ds=dpx)
   call compute_derivative1_fd_centered(s=hs, ds=ds(2), q=sy(1-hs:1+hs), dq_ds=dpy)
   call compute_derivative1_fd_centered(s=hs, ds=ds(3), q=sz(1-hs:1+hs), dq_ds=dpz)
   divb = mhd_sum3(w(1) * dbx, w(2) * dby, w(3) * dbz)
   upsi = mhd_sum3(w(1) * qac(IA_U) * dpx, w(2) * qac(IA_V) * dpy, w(3) * qac(IA_W) * dpz)
   endsubroutine eglm_stencil_scalars

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
endmodule adam_flume_fnl_mhd_eglm_kernels
