!< ADAM, FLUME CPU kernels of the hybrid central/WENO scheme (issue #119, T2).

module adam_flume_cpu_hybrid_kernels
!< ADAM, FLUME CPU kernels of the hybrid central/WENO scheme (issue #119, T2).
!<
!< `compute_hybrid_sensor` fills `psi` on the interior and 2 ghost layers per active direction, the WENO face fluxes having been computed
!< as ever; `overwrite_central_fluxes` replaces the flux of direction `d` by the central one at a face whose cells -1 ... 2 along `d` are
!< all below the threshold and whose stencil -2 ... 3 does not reach the ghost cells of a flagged block face along `d`: `ngc >= 3`. The two
!< blocks of a same-level interface see the same stencil and, the ghosts of `psi` being copies, the same mask: one flux per face, conservative.

! FLUME modules
use :: adam_flume_hybrid_library, only : compute_ducros, compute_kep_flux
use :: adam_flume_parameters,     only : IA_A, IA_P, IA_R, IA_U, IA_V, IA_W
! third party modules
use :: penf,                      only : I4P, R8P

implicit none
private
public :: compute_hybrid_sensor
public :: mask_mu_t
public :: overwrite_central_fluxes

contains
   subroutine compute_hybrid_sensor(ni, nj, nk, ngc, blocks_number, floor, jump_max, dxyz, is_null, q_aux, psi)
   !< Compute the sensor from 2nd-order central velocity gradients; the floor of the denominator is `floor (a / dx_min)^2` per cell.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: floor                             !< Dimensionless floor.
   real(R8P),    intent(in)    :: jump_max                          !< Bound of the relative density jump.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: psi(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Sensor, component 1, q-shaped.
   integer(I4P)                :: s(3)                              !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null.
   real(R8P)                   :: g(3,3)                            !< Velocity gradient.
   real(R8P)                   :: rmax                              !< Largest inverse step.
   real(R8P)                   :: jump                              !< Relative density jump of the cell.
   integer(I4P)                :: b, i, j, k, a                     !< Counters.

   s = merge(0_I4P, 1_I4P, is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, psi)
   do b=1, blocks_number
      do k=1-2*s(3), nk+2*s(3)
         do j=1-2*s(2), nj+2*s(2)
            do i=1-2*s(1), ni+2*s(1)
               rdx = real(s, R8P) / dxyz(:,b)
               rmax = maxval(rdx)
               do a=1, 3
                  g(a,1) = 0.5_R8P * rdx(1) * (q_aux(IA_U+a-1,i+s(1),j,k,b) - q_aux(IA_U+a-1,i-s(1),j,k,b))
                  g(a,2) = 0.5_R8P * rdx(2) * (q_aux(IA_U+a-1,i,j+s(2),k,b) - q_aux(IA_U+a-1,i,j-s(2),k,b))
                  g(a,3) = 0.5_R8P * rdx(3) * (q_aux(IA_U+a-1,i,j,k+s(3),b) - q_aux(IA_U+a-1,i,j,k-s(3),b))
               enddo
               jump = 0._R8P
               jump = max(jump, abs(q_aux(IA_R,i+s(1),j,k,b) - q_aux(IA_R,i-s(1),j,k,b)) / &
                                (q_aux(IA_R,i+s(1),j,k,b) + q_aux(IA_R,i-s(1),j,k,b)))
               jump = max(jump, abs(q_aux(IA_R,i,j+s(2),k,b) - q_aux(IA_R,i,j-s(2),k,b)) / &
                                (q_aux(IA_R,i,j+s(2),k,b) + q_aux(IA_R,i,j-s(2),k,b)))
               jump = max(jump, abs(q_aux(IA_R,i,j,k+s(3),b) - q_aux(IA_R,i,j,k-s(3),b)) / &
                                (q_aux(IA_R,i,j,k+s(3),b) + q_aux(IA_R,i,j,k-s(3),b)))
               call compute_ducros(eps=floor * (q_aux(IA_A,i,j,k,b) * rmax)**2, jump=jump, jump_max=jump_max, g=g, psi=psi(1,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_hybrid_sensor

   subroutine mask_mu_t(ni, nj, nk, ngc, blocks_number, threshold, is_null, psi, mu_t)
   !< Zero the eddy viscosity of the cells the sensor flags (issue #119, T2.6): there WENO supplies the dissipation.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: threshold                         !< Sensor threshold.
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: psi(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Sensor, component 1, q-shaped.
   real(R8P),    intent(inout) :: mu_t(1-ngc:,1-ngc:,1-ngc:,1:)     !< Eddy viscosity.
   integer(I4P)                :: s(3)                              !< Steps along x, y, z: 1 active, 0 null.
   integer(I4P)                :: b, i, j, k                        !< Counters.

   s = merge(0_I4P, 1_I4P, is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(psi, mu_t)
   do b=1, blocks_number
      do k=1-2*s(3), nk+2*s(3)
         do j=1-2*s(2), nj+2*s(2)
            do i=1-2*s(1), ni+2*s(1)
               if (psi(1,i,j,k,b) >= threshold) mu_t(i,j,k,b) = 0._R8P
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine mask_mu_t

   subroutine overwrite_central_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, threshold, faces, q_aux, psi, &
                                       fl)
   !< Replace the flux of direction `d` by the central one where the sensor allows it.
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)    :: threshold                         !< Sensor threshold.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(in)    :: psi(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Sensor, component 1, q-shaped.
   integer(I4P), intent(in)    :: faces(1:,1:)                      !< Block faces not to read across [nb, 6].
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   real(R8P)                   :: w(5,-2:3)                         !< Cell states.
   real(R8P)                   :: f(5)                              !< Central flux.
   integer(I4P)                :: n, l                              !< Cells along `d`, index of the face along `d`.
   integer(I4P)                :: b, i, j, k, m                     !< Counters.

   n = ni * di + nj * dj + nk * dk
   !$omp parallel do collapse(4) default(firstprivate) shared(faces, q_aux, psi, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               ! the stencil -2 ... 3 must not reach the ghost cells of a flagged face of its own direction
               l = i * di + j * dj + k * dk
               if ((l <= 2_I4P .and. faces(b,2*d-1) /= 0_I4P) .or. (l >= n - 2_I4P .and. faces(b,2*d) /= 0_I4P)) cycle
               if (psi(1,i-di,j-dj,k-dk,b) >= threshold .or. psi(1,i,j,k,b) >= threshold .or.        &
                   psi(1,i+di,j+dj,k+dk,b) >= threshold .or. psi(1,i+2*di,j+2*dj,k+2*dk,b) >= threshold) cycle
               do m=-2, 3
                  w(1,m) = q_aux(IA_R,i+m*di,j+m*dj,k+m*dk,b)
                  w(2,m) = q_aux(IA_U,i+m*di,j+m*dj,k+m*dk,b)
                  w(3,m) = q_aux(IA_V,i+m*di,j+m*dj,k+m*dk,b)
                  w(4,m) = q_aux(IA_W,i+m*di,j+m*dj,k+m*dk,b)
                  w(5,m) = q_aux(IA_P,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               call compute_kep_flux(d=d, gamma=gamma, w=w, f=f)
               fl(1:5,i,j,k,b) = f
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine overwrite_central_fluxes
endmodule adam_flume_cpu_hybrid_kernels
