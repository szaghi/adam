!< ADAM, FLUME FNL device kernels of the hybrid central/WENO scheme (issue #119, T2).

#include "fundal.H"

module adam_flume_fnl_hybrid_kernels
!< ADAM, FLUME FNL device kernels of the hybrid central/WENO scheme (issue #119, T2).
!<
!< Device twins of `adam_flume_cpu_hybrid_kernels`: see there. Device arrays are transposed, `(b, i, j, k, v)`.

! FLUME modules
use :: adam_flume_hybrid_library, only : compute_ducros, compute_kep_flux
use :: adam_flume_parameters,     only : IA_A, IA_P, IA_R, IA_U, IA_V, IA_W
! third party modules
use :: penf,                      only : I4P, R8P

implicit none
private
public :: compute_hybrid_sensor_dev
public :: mask_mu_t_dev
public :: overwrite_central_fluxes_dev

contains
   subroutine compute_hybrid_sensor_dev(ni, nj, nk, ngc, blocks_number, floor, jump_max, dxyz_gpu, is_null, q_aux_gpu, psi_gpu)
   !< Device twin of `compute_hybrid_sensor`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: floor                                 !< Dimensionless floor.
   real(R8P),    intent(in)    :: jump_max                              !< Bound of the relative density jump.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                       !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                            !< Null directions.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: psi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Sensor, [nb, i, j, k, 1].
   integer(I4P)                :: s1, s2, s3                            !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: w1, w2, w3                            !< Weights along x, y, z: 1 active, 0 null.
   real(R8P)                   :: rdx(3)                                !< Private inverse steps, 0 along null.
   real(R8P)                   :: g(3,3)                                !< Private velocity gradient.
   real(R8P)                   :: rmax                                  !< Private largest inverse step.
   real(R8P)                   :: p                                     !< Private sensor.
   real(R8P)                   :: jump                                  !< Private relative density jump.
   integer(I4P)                :: b, i, j, k, a                         !< Counters.

   s1 = merge(0_I4P, 1_I4P, is_null(1)) ; s2 = merge(0_I4P, 1_I4P, is_null(2)) ; s3 = merge(0_I4P, 1_I4P, is_null(3))
   w1 = real(s1, R8P) ; w2 = real(s2, R8P) ; w3 = real(s3, R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_aux_gpu,psi_gpu)   &
   !$acc& firstprivate(ni,nj,nk,blocks_number,floor,jump_max,s1,s2,s3,w1,w2,w3) private(rdx,g,rmax,p,jump)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_aux_gpu,psi_gpu)                                  &
   !$omp& firstprivate(ni,nj,nk,blocks_number,floor,jump_max,s1,s2,s3,w1,w2,w3) private(rdx,g,rmax,p,jump)
   do k=1-2*s3, nk+2*s3
   do j=1-2*s2, nj+2*s2
   do i=1-2*s1, ni+2*s1
   do b=1, blocks_number
      rdx(1) = w1 / dxyz_gpu(b,1) ; rdx(2) = w2 / dxyz_gpu(b,2) ; rdx(3) = w3 / dxyz_gpu(b,3)
      rmax = max(rdx(1), rdx(2), rdx(3))
      !$acc loop seq
      do a=1, 3
         g(a,1) = 0.5_R8P * rdx(1) * (q_aux_gpu(b,i+s1,j,k,IA_U+a-1) - q_aux_gpu(b,i-s1,j,k,IA_U+a-1))
         g(a,2) = 0.5_R8P * rdx(2) * (q_aux_gpu(b,i,j+s2,k,IA_U+a-1) - q_aux_gpu(b,i,j-s2,k,IA_U+a-1))
         g(a,3) = 0.5_R8P * rdx(3) * (q_aux_gpu(b,i,j,k+s3,IA_U+a-1) - q_aux_gpu(b,i,j,k-s3,IA_U+a-1))
      enddo
      jump = 0._R8P
      jump = max(jump, abs(q_aux_gpu(b,i+s1,j,k,IA_R) - q_aux_gpu(b,i-s1,j,k,IA_R)) / &
                       (q_aux_gpu(b,i+s1,j,k,IA_R) + q_aux_gpu(b,i-s1,j,k,IA_R)))
      jump = max(jump, abs(q_aux_gpu(b,i,j+s2,k,IA_R) - q_aux_gpu(b,i,j-s2,k,IA_R)) / &
                       (q_aux_gpu(b,i,j+s2,k,IA_R) + q_aux_gpu(b,i,j-s2,k,IA_R)))
      jump = max(jump, abs(q_aux_gpu(b,i,j,k+s3,IA_R) - q_aux_gpu(b,i,j,k-s3,IA_R)) / &
                       (q_aux_gpu(b,i,j,k+s3,IA_R) + q_aux_gpu(b,i,j,k-s3,IA_R)))
      call compute_ducros(eps=floor * (q_aux_gpu(b,i,j,k,IA_A) * rmax)**2, jump=jump, jump_max=jump_max, g=g, psi=p)
      psi_gpu(b,i,j,k,1) = p
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_hybrid_sensor_dev

   subroutine mask_mu_t_dev(ni, nj, nk, ngc, blocks_number, threshold, is_null, psi_gpu, mu_t_gpu)
   !< Zero the eddy viscosity of the cells the sensor flags (issue #119, T2.6): there WENO supplies the dissipation.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                      !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                        !< Actual blocks number.
   real(R8P),    intent(in)    :: threshold                            !< Sensor threshold.
   logical,      intent(in)    :: is_null(3)                           !< Null directions.
   real(R8P),    intent(in)    :: psi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Sensor, [nb, i, j, k, 1].
   real(R8P),    intent(inout) :: mu_t_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Eddy viscosity, [nb, i, j, k, 1].
   integer(I4P)                :: s1, s2, s3                           !< Steps along x, y, z: 1 active, 0 null.
   integer(I4P)                :: b, i, j, k                           !< Counters.

   s1 = merge(0_I4P, 1_I4P, is_null(1)) ; s2 = merge(0_I4P, 1_I4P, is_null(2)) ; s3 = merge(0_I4P, 1_I4P, is_null(3))
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(psi_gpu,mu_t_gpu)  &
   !$acc& firstprivate(ni,nj,nk,blocks_number,threshold,s1,s2,s3)
   !$omp OMPLOOP collapse(4) DEVICEPTR(psi_gpu,mu_t_gpu)                                &
   !$omp& firstprivate(ni,nj,nk,blocks_number,threshold,s1,s2,s3)
   do k=1-2*s3, nk+2*s3
   do j=1-2*s2, nj+2*s2
   do i=1-2*s1, ni+2*s1
   do b=1, blocks_number
      if (psi_gpu(b,i,j,k,1) >= threshold) mu_t_gpu(b,i,j,k,1) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine mask_mu_t_dev

   subroutine overwrite_central_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, threshold, faces_gpu, q_aux_gpu, &
                                           psi_gpu, fl_gpu)
   !< Device twin of `overwrite_central_fluxes`.
   integer(I4P), intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: threshold                             !< Sensor threshold.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)                      !< Block faces not to read across [nb, 6].
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: psi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Sensor, [nb, i, j, k, 1].
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   real(R8P)                   :: w(5,-2:3)                             !< Private cell states.
   real(R8P)                   :: f(5)                                  !< Private central flux.
   integer(I4P)                :: nd, l, lo, hi                         !< Cells along `d`, private face index, its flag indexes.
   integer(I4P)                :: b, i, j, k, m, n                      !< Counters.

   nd = ni * di + nj * dj + nk * dk ; lo = 2_I4P * d - 1_I4P ; hi = 2_I4P * d
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,q_aux_gpu,psi_gpu,fl_gpu)  &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,nd,lo,hi,blocks_number,gamma,threshold) private(w,f,l)
   !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,q_aux_gpu,psi_gpu,fl_gpu)                                 &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,nd,lo,hi,blocks_number,gamma,threshold) private(w,f,l)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      l = i * di + j * dj + k * dk
      if (.not.((l <= 2_I4P .and. faces_gpu(b,lo) /= 0_I4P) .or. (l >= nd - 2_I4P .and. faces_gpu(b,hi) /= 0_I4P)) .and. &
          psi_gpu(b,i-di,j-dj,k-dk,1) < threshold .and. psi_gpu(b,i,j,k,1) < threshold .and.                                  &
          psi_gpu(b,i+di,j+dj,k+dk,1) < threshold .and. psi_gpu(b,i+2*di,j+2*dj,k+2*dk,1) < threshold) then
         !$acc loop seq
         do m=-2, 3
            w(1,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,IA_R)
            w(2,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,IA_U)
            w(3,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,IA_V)
            w(4,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,IA_W)
            w(5,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,IA_P)
         enddo
         call compute_kep_flux(d=d, gamma=gamma, w=w, f=f)
         !$acc loop seq
         do n=1, 5
            fl_gpu(b,i,j,k,n) = f(n)
         enddo
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine overwrite_central_fluxes_dev
endmodule adam_flume_fnl_hybrid_kernels
