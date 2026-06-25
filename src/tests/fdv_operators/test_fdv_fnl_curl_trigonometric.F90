program test_fdv_fnl_curl_trigonometric
!< Test ADAM fdv FNL operators (namely curl) library with trigonometric functions derivatives, device backend.

!< Test function, vectorial field
!<```
!< q(1:3,x,y,z) = [sin(y), sin(z), sin(x)], (x,y,z) in [0,2PI]x[0,2PI]x[0,2PI].
!<```
!< The exact curl is
!<```
!< curl q(1:3,x,y,z) = [qz,y - qy,z, qx,z - qz,x, qy,x - qx,y]  = [-cos(z), -cos(x), -cos(y)]
!<```
!<
!< Data is initialized on host, copied to device, curl computed on device, result copied back to host.
!< Device memory is managed via FUNDAL dev_alloc/dev_free/dev_memcpy routines.
!< The FNL layout has variables as the last dimension: q(i,j,k,nv) instead of q(nv,i,j,k).

#include "fundal.H"

use :: adam_fnl_fdv_operators_library
use :: fundal
use :: penf

implicit none

integer(I4P), parameter :: n=20_I4P               !< Number of cells in each direction.
integer(I4P), parameter :: gc=5_I4P               !< Number of ghost cells.
real(R8P),    parameter :: PI=4._R8P*atan(1._R8P) !< PI.
real(R8P),    parameter :: h=2_I4P*PI/n           !< Space step.

! host arrays - CPU layout q(nv, i, j, k)
real(R8P) :: q(    1:3, 1-gc:n+gc, 1-gc:n+gc, 1-gc:n+gc) !< Function (vectorial) [sin(y), sin(z), sin(x)].
real(R8P) :: exact(1:3, 1:n,       1:n,        1:n     ) !< Exact curl function.
real(R8P) :: curl( 1:3, 1:n,       1:n,        1:n     ) !< Curl (vectorial) function.
real(R8P) :: error(1:3, 1:n,       1:n,        1:n     ) !< Error function.

! transposed host arrays for device transfer - FNL layout q(i, j, k, nv)
real(R8P) :: q_t(   1-gc:n+gc, 1-gc:n+gc, 1-gc:n+gc, 1:3) !< Transposed field for H2D copy.
real(R8P) :: curl_t(1:n,       1:n,       1:n,       1:3) !< Transposed curl result from D2H copy.

! device pointers (allocated on device via fundal)
real(R8P), pointer :: q_dev(:,:,:,:)    !< Device input field.
real(R8P), pointer :: curl_dev(:,:,:,:) !< Device curl result.

real(R8P)    :: errorL2
integer(I4P) :: s, o, i, j, k, d, ierr

! initialize fundal device environment
myhos   = dev_get_host_num()
devtype = dev_get_device_type()
call dev_set_device_num(0)
mydev = dev_get_device_num()

! initialize host data (same field as CPU reference test)
do k=1-gc, n+gc
   do j=1-gc, n+gc
      do i=1-gc, n+gc
         q(1,i,j,k) = sin(j*h)
         q(2,i,j,k) = sin(k*h)
         q(3,i,j,k) = sin(i*h)
         if (i>=1.and.j>=1.and.k>=1.and.i<=n.and.j<=n.and.k<=n) then
            exact(1,i,j,k) = -cos(k*h)
            exact(2,i,j,k) = -cos(i*h)
            exact(3,i,j,k) = -cos(j*h)
         endif
      enddo
   enddo
enddo

! transpose to FNL layout (variables last) for H2D copy
do k=1-gc, n+gc
   do j=1-gc, n+gc
      do i=1-gc, n+gc
         q_t(i,j,k,1:3) = q(1:3,i,j,k)
      enddo
   enddo
enddo

! allocate device arrays
call dev_alloc(fptr_dev=q_dev,    lbounds=[1-gc,1-gc,1-gc,1], ubounds=[n+gc,n+gc,n+gc,3], ierr=ierr)
call dev_alloc(fptr_dev=curl_dev, lbounds=[1,1,1,1],          ubounds=[n,n,n,3],          ierr=ierr)

! copy input field to device (H2D)
call dev_memcpy_to_device(dst=q_dev, src=q_t)

! test FD curl on device
print '(A)', 'FD schemes (FNL device backend)'
do o=2, 10, 2
   s = o/2
   call compute_curl_fd_dev_kernel(s=s, ni=n, gc_=gc, q_dev=q_dev, curl_dev=curl_dev)
   ! copy curl result back to host (D2H)
   call dev_memcpy_from_device(dst=curl_t, src=curl_dev)
   ! transpose back to CPU layout (variables first)
   do k=1, n
      do j=1, n
         do i=1, n
            curl(1:3,i,j,k) = curl_t(i,j,k,1:3)
         enddo
      enddo
   enddo
   error = abs(curl - exact)
   errorL2 = 0._R8P
   do k=1, n ; do j=1, n ; do i=1, n ; do d=1, 3
      errorL2 = errorL2 + error(d,i,j,k)*error(d,i,j,k)
   enddo ; enddo ; enddo ; enddo
   errorL2 = sqrt(errorL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
enddo

! test FV curl on device (rebuilt finite-volume centered scheme via gathered buffers)
print '(A)', 'FV schemes (FNL device backend)'
do o=2, 10, 2
   s = o/2
   call compute_curl_fv_dev_kernel(s=s, ni=n, gc_=gc, q_dev=q_dev, curl_dev=curl_dev)
   call dev_memcpy_from_device(dst=curl_t, src=curl_dev)
   do k=1, n
      do j=1, n
         do i=1, n
            curl(1:3,i,j,k) = curl_t(i,j,k,1:3)
         enddo
      enddo
   enddo
   error = abs(curl - exact)
   errorL2 = 0._R8P
   do k=1, n ; do j=1, n ; do i=1, n ; do d=1, 3
      errorL2 = errorL2 + error(d,i,j,k)*error(d,i,j,k)
   enddo ; enddo ; enddo ; enddo
   errorL2 = sqrt(errorL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
enddo

! free device arrays
call dev_free(q_dev,    mydev)
call dev_free(curl_dev, mydev)

contains
   subroutine compute_curl_fd_dev_kernel(s, ni, gc_, q_dev, curl_dev)
   !< Compute curl on device using FNL centered FD scheme.
   !< Gathers non-contiguous stencil slices into private local arrays before calling the device routine,
   !< following the same pattern as adam_prism_fnl_object to avoid passing non-contiguous sections to acc routines.
   integer(I4P), intent(in)    :: s                                    !< Half stencil length (order = 2*s).
   integer(I4P), intent(in)    :: ni                                   !< Number of interior cells per direction.
   integer(I4P), intent(in)    :: gc_                                  !< Number of ghost cells.
   real(R8P),    intent(in)    :: q_dev(   1-gc_:, 1-gc_:, 1-gc_:, 1:) !< Device field, FNL layout (i,j,k,nv).
   real(R8P),    intent(inout) :: curl_dev(1:,     1:,     1:,     1:) !< Device curl result (i,j,k,nv).
   integer(I4P) :: i, j, k, ist
   ! private stencil arrays (local automatic, made private per thread by OpenACC)
   real(R8P) :: qsx_y(1-s:1+s) !< Y component over x stencil.
   real(R8P) :: qsx_z(1-s:1+s) !< Z component over x stencil.
   real(R8P) :: qsy_x(1-s:1+s) !< X component over y stencil.
   real(R8P) :: qsy_z(1-s:1+s) !< Z component over y stencil.
   real(R8P) :: qsz_x(1-s:1+s) !< X component over z stencil.
   real(R8P) :: qsz_y(1-s:1+s) !< Y component over z stencil.
   real(R8P) :: dxyz(3) = [h,h,h]
   real(R8P) :: curl_(3)

   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(q_dev, curl_dev)                       &
   !$acc& firstprivate(s, ni, dxyz)                        &
   !$acc& private(qsx_y, qsx_z, qsy_x, qsy_z, qsz_x, qsz_y, curl_)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(q_dev, curl_dev) &
   !$omp& firstprivate(s, ni, dxyz) &
   !$omp& private(qsx_y, qsx_z, qsy_x, qsy_z, qsz_x, qsz_y, curl_)
   do k=1, ni
   do j=1, ni
   do i=1, ni
      !$acc loop seq
      do ist=1-s, 1+s
         qsx_y(ist) = q_dev(i+ist-1, j,       k,       2)
         qsx_z(ist) = q_dev(i+ist-1, j,       k,       3)
         qsy_x(ist) = q_dev(i,       j+ist-1, k,       1)
         qsy_z(ist) = q_dev(i,       j+ist-1, k,       3)
         qsz_x(ist) = q_dev(i,       j,       k+ist-1, 1)
         qsz_y(ist) = q_dev(i,       j,       k+ist-1, 2)
      enddo
      call compute_curl_fd_centered_dev(s=s, dxyz=dxyz,          &
                                        qsx_y=qsx_y, qsx_z=qsx_z,&
                                        qsy_x=qsy_x, qsy_z=qsy_z,&
                                        qsz_x=qsz_x, qsz_y=qsz_y,&
                                        curl=curl_)
      curl_dev(i,j,k,1) = curl_(1)
      curl_dev(i,j,k,2) = curl_(2)
      curl_dev(i,j,k,3) = curl_(3)
   enddo
   enddo
   enddo
   endsubroutine compute_curl_fd_dev_kernel

   subroutine compute_curl_fv_dev_kernel(s, ni, gc_, q_dev, curl_dev)
   !< Compute curl on device using the rebuilt FNL centered FV scheme.
   !< Same gather as the FD kernel (six contiguous directional stencils); calls compute_curl_fv_centered_dev.
   integer(I4P), intent(in)    :: s                                    !< Half stencil length (order = 2*s).
   integer(I4P), intent(in)    :: ni                                   !< Number of interior cells per direction.
   integer(I4P), intent(in)    :: gc_                                  !< Number of ghost cells.
   real(R8P),    intent(in)    :: q_dev(   1-gc_:, 1-gc_:, 1-gc_:, 1:) !< Device field, FNL layout (i,j,k,nv).
   real(R8P),    intent(inout) :: curl_dev(1:,     1:,     1:,     1:) !< Device curl result (i,j,k,nv).
   integer(I4P) :: i, j, k, ist
   real(R8P) :: qsx_y(1-s:1+s) !< Y component over x stencil.
   real(R8P) :: qsx_z(1-s:1+s) !< Z component over x stencil.
   real(R8P) :: qsy_x(1-s:1+s) !< X component over y stencil.
   real(R8P) :: qsy_z(1-s:1+s) !< Z component over y stencil.
   real(R8P) :: qsz_x(1-s:1+s) !< X component over z stencil.
   real(R8P) :: qsz_y(1-s:1+s) !< Y component over z stencil.
   real(R8P) :: dxyz(3) = [h,h,h]
   real(R8P) :: curl_(3)

   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(q_dev, curl_dev)                       &
   !$acc& firstprivate(s, ni, dxyz)                        &
   !$acc& private(qsx_y, qsx_z, qsy_x, qsy_z, qsz_x, qsz_y, curl_)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(q_dev, curl_dev) &
   !$omp& firstprivate(s, ni, dxyz) &
   !$omp& private(qsx_y, qsx_z, qsy_x, qsy_z, qsz_x, qsz_y, curl_)
   do k=1, ni
   do j=1, ni
   do i=1, ni
      !$acc loop seq
      do ist=1-s, 1+s
         qsx_y(ist) = q_dev(i+ist-1, j,       k,       2)
         qsx_z(ist) = q_dev(i+ist-1, j,       k,       3)
         qsy_x(ist) = q_dev(i,       j+ist-1, k,       1)
         qsy_z(ist) = q_dev(i,       j+ist-1, k,       3)
         qsz_x(ist) = q_dev(i,       j,       k+ist-1, 1)
         qsz_y(ist) = q_dev(i,       j,       k+ist-1, 2)
      enddo
      call compute_curl_fv_centered_dev(s=s, dxyz=dxyz,          &
                                        qsx_y=qsx_y, qsx_z=qsx_z,&
                                        qsy_x=qsy_x, qsy_z=qsy_z,&
                                        qsz_x=qsz_x, qsz_y=qsz_y,&
                                        curl=curl_)
      curl_dev(i,j,k,1) = curl_(1)
      curl_dev(i,j,k,2) = curl_(2)
      curl_dev(i,j,k,3) = curl_(3)
   enddo
   enddo
   enddo
   endsubroutine compute_curl_fv_dev_kernel
endprogram test_fdv_fnl_curl_trigonometric
