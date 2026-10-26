!< RIPPLE, the FNL (OpenACC) device kernels.

#include "fundal.H"

module adam_ripple_fnl_kernels
!< RIPPLE, the FNL (OpenACC) device kernels.
!<
!< The device twins of the CPU loops. Device arrays are transposed, `(b, i, j, k, v)`: the block index is stride-1, so
!< it is the innermost loop of every nest. Host scalars are `firstprivate`, device arrays `DEVICEVAR` (FUNDAL's
!< `deviceptr`), private work arrays have constant bounds.

! ADAM classes, libraries, parameters
use :: adam_fnl_weno_kernels, only : weno_reconstruct_upwind_dev
use :: adam_weno_object,      only : S_MAX
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: apply_reflux_face_dev
public :: compute_face_fluxes_dev
public :: compute_flux_difference_dev
public :: pack_face_skin_dev

contains
   subroutine compute_face_fluxes_dev(di, dj, dk, ni, nj, nk, ngc, blocks_number, S, a, weno_a_gpu, weno_p_gpu, &
                                      weno_d_gpu, weno_zeps, q_gpu, fl_gpu)
   !< Face fluxes along the axis of unit step (di, dj, dk), device twin of the CPU `face_fluxes`: the split fluxes
   !< `a+ u` (cells `m = 1-S ... S-1`) and `a- u` (cells `m+1`) of the stencil of face `i+1/2`, reconstructed by WENO
   !< and summed.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Blocks in use.
   integer(I4P), intent(in)    :: S                                 !< WENO half stencil.
   real(R8P),    intent(in)    :: a                                 !< Velocity component.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)              !< WENO optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)           !< WENO polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)           !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                         !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)   !< Face fluxes.
   real(R8P)                   :: v(2,2*S_MAX-1)                    !< Private packed stencil.
   real(R8P)                   :: vr(2)                             !< Private reconstructed split fluxes.
   real(R8P)                   :: ap, am                            !< Split velocities.
   integer(I4P)                :: b, i, j, k, m                     !< Counters.

   ap = max(a, 0._R8P)
   am = min(a, 0._R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,fl_gpu) &
   !$acc& firstprivate(di,dj,dk,ni,nj,nk,blocks_number,S,ap,am,weno_zeps) private(v,vr)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=1-S, S-1
         v(2,m+S) = ap * q_gpu(b,i+m*di,j+m*dj,k+m*dk,1)
         v(1,m+S) = am * q_gpu(b,i+(m+1)*di,j+(m+1)*dj,k+(m+1)*dk,1)
      enddo
      call weno_reconstruct_upwind_dev(S=S, weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu, &
                                       weno_zeps=weno_zeps, weno_rmu=1._R8P, V=v, VR=vr)
      fl_gpu(b,i,j,k,1) = vr(1) + vr(2)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_face_fluxes_dev

   subroutine compute_flux_difference_dev(ni, nj, nk, ngc, blocks_number, w, dxyz_gpu, flx_gpu, fly_gpu, flz_gpu, dq_gpu)
   !< Residuals `dq = -sum_d w_d (F_{i+1/2} - F_{i-1/2}) / dx_d`, `w_d` 1 on an active axis, 0 on a null one.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Blocks in use.
   real(R8P),    intent(in)    :: w(3)                               !< Axis weights.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Block spacings [nb, 3].
   real(R8P),    intent(in)    :: flx_gpu(1:,0:,1:,1:,1:)            !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_gpu(1:,1:,0:,1:,1:)            !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_gpu(1:,1:,1:,0:,1:)            !< Z-face fluxes.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P)                   :: wx, wy, wz                         !< Axis weights.
   integer(I4P)                :: b, i, j, k                         !< Counters.

   wx = w(1) ; wy = w(2) ; wz = w(3)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,flx_gpu,fly_gpu,flz_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,wx,wy,wz)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dq_gpu(b,i,j,k,1) = -(wx * (flx_gpu(b,i,j,k,1) - flx_gpu(b,i-1,j,k,1)) / dxyz_gpu(b,1) + &
                            wy * (fly_gpu(b,i,j,k,1) - fly_gpu(b,i,j-1,k,1)) / dxyz_gpu(b,2) + &
                            wz * (flz_gpu(b,i,j,k,1) - flz_gpu(b,i,j,k-1,1)) / dxyz_gpu(b,3))
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_flux_difference_dev

   subroutine pack_face_skin_dev(fec, b, ni, nj, nk, flx_gpu, fly_gpu, flz_gpu, skin_gpu)
   !< Pack the face fluxes of face `fec` of block `b` into `skin_gpu(1, c)`, `c` running over the two tangential axes,
   !< inner fastest (x faces: j, k; y faces: i, k; z faces: i, j): the flux register's order.
   integer(I4P), intent(in)    :: fec, b                  !< Face (1..6: -x, +x, -y, +y, -z, +z), block.
   integer(I4P), intent(in)    :: ni, nj, nk              !< Grid dimensions.
   real(R8P),    intent(in)    :: flx_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   real(R8P),    intent(inout) :: skin_gpu(1:,1:)         !< Skin (1, cells).
   integer(I4P)                :: f, t1, t2               !< Face index along the normal, tangential counters.

   f = merge(0_I4P, 1_I4P, mod(fec, 2_I4P) == 1_I4P) ! the min face is face 0, the max face is face n
   select case(fec)
   case(1, 2)
      !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(flx_gpu,skin_gpu) firstprivate(b,f,ni,nj,nk)
      do t2=1, nk
      do t1=1, nj
         skin_gpu(1,t1+(t2-1)*nj) = flx_gpu(b,f*ni,t1,t2,1)
      enddo
      enddo
   case(3, 4)
      !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(fly_gpu,skin_gpu) firstprivate(b,f,ni,nj,nk)
      do t2=1, nk
      do t1=1, ni
         skin_gpu(1,t1+(t2-1)*ni) = fly_gpu(b,t1,f*nj,t2,1)
      enddo
      enddo
   case default
      !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(flz_gpu,skin_gpu) firstprivate(b,f,ni,nj,nk)
      do t2=1, nj
      do t1=1, ni
         skin_gpu(1,t1+(t2-1)*ni) = flz_gpu(b,t1,t2,f*nk,1)
      enddo
      enddo
   endselect
   endsubroutine pack_face_skin_dev

   subroutine apply_reflux_face_dev(axis, sgn, b, ni, nj, nk, ngc, nface_cells, scale, delta_gpu, q_gpu)
   !< Add `scale delta(1, c)` to the coarse skin cells of block `b` on the face of normal `axis` and side `sgn`.
   integer(I4P), intent(in)    :: axis, sgn                         !< Face normal axis (1..3) and side (+-1).
   integer(I4P), intent(in)    :: b                                 !< Coarse block.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nface_cells                       !< Skin cells.
   real(R8P),    intent(in)    :: scale                             !< sgn dt / dx_coarse.
   real(R8P),    intent(in)    :: delta_gpu(1:,1:)                  !< Flux mismatch F_coarse - F_fine (1, cells).
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State.
   integer(I4P)                :: c, i, j, k                        !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(delta_gpu,q_gpu) &
   !$acc& firstprivate(axis,sgn,b,ni,nj,nk,nface_cells,scale) private(i,j,k)
   do c=1, nface_cells
      select case(axis)
      case(1) ; i = merge(ni, 1_I4P, sgn > 0) ; j = 1 + mod(c - 1, nj) ; k = 1 + (c - 1) / nj
      case(2) ; i = 1 + mod(c - 1, ni) ; j = merge(nj, 1_I4P, sgn > 0) ; k = 1 + (c - 1) / ni
      case default ; i = 1 + mod(c - 1, ni) ; j = 1 + (c - 1) / ni ; k = merge(nk, 1_I4P, sgn > 0)
      endselect
      q_gpu(b,i,j,k,1) = q_gpu(b,i,j,k,1) + scale * delta_gpu(1,c)
   enddo
   endsubroutine apply_reflux_face_dev
endmodule adam_ripple_fnl_kernels
