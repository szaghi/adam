!< ADAM, IB class FNL kernels (FNL backend of [[ib_object]]).

#include "fundal.H"

module adam_fnl_ib_kernels
!< ADAM, IB class FNL kernels (FNL backend of [[ib_object]]).

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_eikonal_dq_phi_dev
public :: compute_phi_all_solids_dev
public :: compute_phi_analytical_sphere_dev
public :: evolve_eikonal_q_phi_dev
public :: invert_eikonal_q_phi_dev
public :: move_phi_dev
public :: reduce_cell_order_phi_dev

contains
   subroutine compute_eikonal_dq_phi_dev(ib, ni, nj, nk, ngc, nv, blocks_number, dx_gpu, dy_gpu, dz_gpu, phi_gpu, q_gpu, dq_gpu)
   !< Compute eikonal dq inside IB.
   integer(I4P), intent(in)    :: ib                                  !< IB solid index.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in)    :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in)    :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in)    :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in)    :: q_gpu(  1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State variables.
   real(R8P),    intent(inout) :: dq_gpu( 1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State variables variations.
   integer(I4P)                :: i, j, k, b, v                       !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z, n_phi    !< Metrics.

   !$acc parallel loop independent DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, phi_gpu, q_gpu, dq_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dx_gpu, dy_gpu, dz_gpu, phi_gpu, q_gpu, dq_gpu)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            do b=1, blocks_number
               if (phi_gpu(b,i,j,k,ib) > 0._R8P) then
                  n_phi_x = (phi_gpu(b,i+1,j,k,ib) - phi_gpu(b,i-1,j,k,ib) )
                  n_phi_y = (phi_gpu(b,i,j+1,k,ib) - phi_gpu(b,i,j-1,k,ib) )
                  n_phi_z = (phi_gpu(b,i,j,k+1,ib) - phi_gpu(b,i,j,k-1,ib) )
                  n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
                  n_phi = 0.9_R8P / n_phi
                  n_phi_x = n_phi_x * n_phi
                  n_phi_y = n_phi_y * n_phi
                  n_phi_z = n_phi_z * n_phi
                  do v=1, nv
                     dq_gpu(b,i,j,k,v) = 0._R8P
                  enddo
                  if (n_phi_x > 0._R8P) then
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i-1,j,k,v))
                     enddo
                  else
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i+1,j,k,v))
                     enddo
                  endif
                  if (n_phi_y > 0._R8P) then
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j-1,k,v))
                     enddo
                  else
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j+1,k,v))
                     enddo
                  endif
                  if (n_phi_z > 0._R8P) then
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k-1,v))
                     enddo
                  else
                     do v=1, nv
                        dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k+1,v))
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_eikonal_dq_phi_dev

   subroutine compute_phi_all_solids_dev(ni, nj, nk, ngc, blocks_number, phi_gpu)
   !< Compute last phi index, all solids summary.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost grid number.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(inout) :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   integer(I4P)                :: b, i, j, k, s                       !< Counter.
   integer(I4P)                :: all_solids                          !< Last phi index, all solids summary.

   all_solids = ubound(phi_gpu, dim=5)
   !$acc parallel loop independent DEVICEVAR(phi_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            do b=1, blocks_number
               phi_gpu(b,i,j,k,all_solids) = phi_gpu(b,i,j,k,1)
               solids_loop : do s=2, all_solids -1
                  phi_gpu(b,i,j,k,all_solids) = max(phi_gpu(b,i,j,k,all_solids), phi_gpu(b,i,j,k,s))
               enddo solids_loop
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_phi_all_solids_dev

   subroutine compute_phi_analytical_sphere_dev(ib,ni,nj,nk,ngc,blocks_number,sphere,x_cell_gpu,y_cell_gpu,z_cell_gpu,phi_gpu)
   !< Compute phi, distance from ib solid, analytical sphere solid.
   integer(I4P), intent(in)    :: ib                                  !< IB solid index.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost grid number.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in)    :: sphere(1:4)                         !< Sphere center [1:3] and radius [4].
   real(R8P),    intent(in)    :: x_cell_gpu(1:,1-ngc:)               !< Cells x coordinates on GPU.
   real(R8P),    intent(in)    :: y_cell_gpu(1:,1-ngc:)               !< Cells y coordinates on GPU.
   real(R8P),    intent(in)    :: z_cell_gpu(1:,1-ngc:)               !< Cells z coordinates on GPU.
   real(R8P),    intent(inout) :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P)                   :: sphere_center_x                     !< Sphere center x coordinate.
   real(R8P)                   :: sphere_center_y                     !< Sphere center y coordinate.
   real(R8P)                   :: sphere_center_z                     !< Sphere center z coordinate.
   real(R8P)                   :: sphere_radius                       !< Sphere radius.
   integer(I4P)                :: b, i, j, k                          !< Counter.

   sphere_center_x = sphere(1)
   sphere_center_y = sphere(2)
   sphere_center_z = sphere(3)
   sphere_radius   = sphere(4)
   !$acc parallel loop independent DEVICEVAR(x_cell_gpu, y_cell_gpu, z_cell_gpu, phi_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(x_cell_gpu, y_cell_gpu, z_cell_gpu, phi_gpu)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            do b=1, blocks_number
               phi_gpu(b,i,j,k,ib) = -(sqrt((x_cell_gpu(b,i) - sphere_center_x)**2 + &
                                            (y_cell_gpu(b,j) - sphere_center_y)**2 + &
                                            (z_cell_gpu(b,k) - sphere_center_z)**2) - sphere_radius)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_phi_analytical_sphere_dev

   subroutine evolve_eikonal_q_phi_dev(ib, ni, nj, nk, ngc, nv, blocks_number, phi_gpu, dq_gpu, q_gpu)
   !< Evolve eikonal equation over q inside IB.
   integer(I4P), intent(in)    :: ib                                  !< IB solid index.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in)    :: dq_gpu( 1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State variables variation.
   real(R8P),    intent(inout) :: q_gpu(  1:,1-ngc:,1-ngc:,1-ngc:,1:) !< State variables.
   integer(I4P)                :: i, j, k, b, v                       !< Counter.

   !$acc parallel loop independent DEVICEVAR(phi_gpu, dq_gpu, q_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, dq_gpu, q_gpu)
   do k=1, nk
      do j=1, nj
         do i=1,ni
            do b=1, blocks_number
               do v=1, nv
                  if (phi_gpu(b,i,j,k,ib) > 0._R8P) then
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) - dq_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine evolve_eikonal_q_phi_dev

   subroutine invert_eikonal_q_phi_dev(BCS_VISCOUS,BCS_EULER,ib,ni,nj,nk,ngc,nv,blocks_number,bcs_type,phi_gpu,q_gpu)
   !< Invert eikonal equation over q inside IB.
   integer(I4P), intent(in)    :: BCS_VISCOUS                         !< Viscous wall BCS parameter.
   integer(I4P), intent(in)    :: BCS_EULER                           !< Euler wall BCS parameter.
   integer(I4P), intent(in)    :: ib                                  !< IB solid index.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   integer(I4P), intent(in)    :: bcs_type                            !< Immersed boundary type.
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance field.
   real(R8P),    intent(inout) :: q_gpu(  1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   integer(I4P)                :: i, j, k, b                          !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z           !< Distance function normals.
   real(R8P)                   :: n_phi_mod, un_mod                   !< Distance abs normal and normal velocity.

   if     (bcs_type == BCS_VISCOUS) then
      !$acc parallel loop independent DEVICEVAR(phi_gpu, q_gpu)
      !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, q_gpu)
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if (phi_gpu(b,i,j,k,ib) > 0) then
                     q_gpu(b,i,j,k,2) = - q_gpu(b,i,j,k,2)
                     q_gpu(b,i,j,k,3) = - q_gpu(b,i,j,k,3)
                     q_gpu(b,i,j,k,4) = - q_gpu(b,i,j,k,4)
                  endif
               enddo
            enddo
         enddo
      enddo
   elseif (bcs_type == BCS_EULER  ) then
      !$acc parallel loop independent DEVICEVAR(phi_gpu, q_gpu)
      !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, q_gpu)
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if (phi_gpu(b,i,j,k,ib) > 0) then
                     n_phi_x = phi_gpu(b,i+1,j,k,ib) - phi_gpu(b,i-1,j,k,ib)
                     n_phi_y = phi_gpu(b,i,j+1,k,ib) - phi_gpu(b,i,j-1,k,ib)
                     n_phi_z = phi_gpu(b,i,j,k+1,ib) - phi_gpu(b,i,j,k-1,ib)
                     n_phi_mod = sqrt(n_phi_x**2 + n_phi_y**2 + n_phi_z**2)
                     n_phi_x = n_phi_x/n_phi_mod
                     n_phi_y = n_phi_y/n_phi_mod
                     n_phi_z = n_phi_z/n_phi_mod
                     un_mod = q_gpu(b,i,j,k,2)*n_phi_x + q_gpu(b,i,j,k,3)*n_phi_y + q_gpu(b,i,j,k,4)*n_phi_z

                     q_gpu(b,i,j,k,2) = q_gpu(b,i,j,k,2) - 2*un_mod*n_phi_x
                     q_gpu(b,i,j,k,3) = q_gpu(b,i,j,k,3) - 2*un_mod*n_phi_y
                     q_gpu(b,i,j,k,4) = q_gpu(b,i,j,k,4) - 2*un_mod*n_phi_z
                  endif
               enddo
            enddo
         enddo
      enddo
   endif
   endsubroutine invert_eikonal_q_phi_dev

   subroutine move_phi_dev(ni, nj, nk, ngc, blocks_number, velocity, phi_gpu, dphi_gpu)
   !< Move phi and the actual ptree representation.
   integer(I4P), intent(in)    :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                  !< Ghost grid number.
   integer(I4P), intent(in)    :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)    :: velocity(3)                          !< Velocity of the movement.
   real(R8P),    intent(inout) ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(inout) :: dphi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function gradient.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z, n_phi     !< Eikonal direction.
   integer(I4P)                :: b, i, j, k, v                        !< Counter.
   integer(I4P)                :: solids_number                        !< Solids number.

   n_phi_x = velocity(1)
   n_phi_y = velocity(2)
   n_phi_z = velocity(3)
   n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
   n_phi = 0.9_R8P / n_phi
   n_phi_x = n_phi_x * n_phi
   n_phi_y = n_phi_y * n_phi
   n_phi_z = n_phi_z * n_phi

   solids_number = ubound(phi_gpu, dim=5)
   !$acc parallel loop independent DEVICEVAR(phi_gpu, dphi_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, dphi_gpu)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, solids_number
         dphi_gpu(b,i,j,k,v) = 0._R8P
      enddo
      if (n_phi_x > 0._R8P) then
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i-1,j,k,v))
         enddo
      else
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i+1,j,k,v))
         enddo
      endif
      if (n_phi_y > 0._R8P) then
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j-1,k,v))
         enddo
      else
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j+1,k,v))
         enddo
      endif
      if (n_phi_z > 0._R8P) then
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k-1,v))
         enddo
      else
         do v=1, solids_number
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k+1,v))
         enddo
      endif
   enddo
   enddo
   enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(phi_gpu, dphi_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, dphi_gpu)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, solids_number
         phi_gpu(b,i,j,k,v) = phi_gpu(b,i,j,k,v) - dphi_gpu(b,i,j,k,v)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine move_phi_dev

   subroutine reduce_cell_order_phi_dev(ib, ni, nj, nk, ngc, blocks_number, iweno, &
                                        ib_reduced_order, ib_reduction_extent, phi_gpu, cell_scheme_gpu)
   !< Reduce local-cell order of spatial operator close to solids.
   integer(I4P), intent(in)    :: ib                                          !< IB solid index.
   integer(I4P), intent(in)    :: ni                                          !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                          !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                          !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                         !< Ghost grid number.
   integer(I4P), intent(in)    :: blocks_number                               !< Number of blocks.
   integer(I4P), intent(in)    :: iweno                                       !< WENO (half) stencil lenght.
   integer(I4P), intent(in)    :: ib_reduced_order                            !< Reduced order close to IB solids.
   integer(I4P), intent(in)    :: ib_reduction_extent                         !< Extent of reduction close to IBi solids.
   real(R8P),    intent(in)    :: phi_gpu(        1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   integer(I4P), intent(inout) :: cell_scheme_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Modified order close to solids.
   integer(I4P)                :: b, i, j, k, l                               !< Counter.

   !$acc parallel loop independent DEVICEVAR(phi_gpu, cell_scheme_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, cell_scheme_gpu)
   do b=1, blocks_number
      do j=1, nj
         do k=1, nk
            idir: do i=0, ni
               cell_scheme_gpu(b,i,j,k,1) = iweno
               do l=1, ib_reduction_extent
                  if (phi_gpu(b,i+l,j,k,ib)>0.or.phi_gpu(b,i-l+1,j,k,ib)>0) then
                      cell_scheme_gpu(b,i,j,k,1) = ib_reduced_order
                      cycle idir
                  endif
               enddo
            enddo idir
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(phi_gpu, cell_scheme_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, cell_scheme_gpu)
   do b=1, blocks_number
      do i=1, ni
         do k=1, nk
            jdir: do j=0,nj
               cell_scheme_gpu(b,i,j,k,2) = iweno
               do l=1, ib_reduction_extent
                  if (phi_gpu(b,i,j+l,k,ib)>0.or.phi_gpu(b,i,j-l+1,k,ib)>0) then
                      cell_scheme_gpu(b,i,j,k,2) = ib_reduced_order
                      cycle jdir
                  endif
               enddo
            enddo jdir
         enddo
      enddo
   enddo

   !$acc parallel loop independent DEVICEVAR(phi_gpu, cell_scheme_gpu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(phi_gpu, cell_scheme_gpu)
   do b=1, blocks_number
      do j=1, nj
         do i=1, ni
            kdir: do k=0, nk
               cell_scheme_gpu(b,i,j,k,3) = iweno
               do l=1, ib_reduction_extent
                  if (phi_gpu(b,i,j,k+l,ib)>0.or.phi_gpu(b,i,j,k-l+1,ib)>0) then
                      cell_scheme_gpu(b,i,j,k,3) = ib_reduced_order
                      cycle kdir
                  endif
               enddo
            enddo kdir
         enddo
      enddo
   enddo
   endsubroutine reduce_cell_order_phi_dev
endmodule adam_fnl_ib_kernels
