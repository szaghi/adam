program test_fdv_fnl_operators_vector_trigonometric
!< Test ADAM fdv FNL vector operators (divergence, gradient, laplacian) with trigonometric fields, device backend.
!<
!< Scalar field  f(x,y,z) = sin(x) + sin(y) + sin(z):
!<   grad f      = [cos(x), cos(y), cos(z)]
!<   laplacian f = -sin(x) - sin(y) - sin(z)
!< Vector field  v(x,y,z) = [sin(x), sin(y), sin(z)]:
!<   div v       = cos(x) + cos(y) + cos(z)
!<
!< Each operator is computed on device via the gathered-buffer FV centered _dev routines, exercising the
!< rebuilt FV derivative ladder (d1 for div/grad, d2 for laplacian) on the GPU. FNL layout: variables last.

#include "fundal.H"

use :: adam_fnl_fdv_operators_library
use :: fundal
use :: penf

implicit none

integer(I4P), parameter :: n=20_I4P               !< Number of cells in each direction.
integer(I4P), parameter :: gc=5_I4P               !< Number of ghost cells.
real(R8P),    parameter :: PI=4._R8P*atan(1._R8P) !< PI.
real(R8P),    parameter :: h=2_I4P*PI/n           !< Space step.

! host scalar/vector fields (FNL layout: spatial indices then variables)
real(R8P) :: f_t(1-gc:n+gc, 1-gc:n+gc, 1-gc:n+gc)      !< Scalar field f, device layout.
real(R8P) :: v_t(1-gc:n+gc, 1-gc:n+gc, 1-gc:n+gc, 1:3) !< Vector field v, device layout (i,j,k,nv).

! exact references
real(R8P) :: exact_grad(1:3, 1:n, 1:n, 1:n) !< Exact gradient of f.
real(R8P) :: exact_lap(      1:n, 1:n, 1:n) !< Exact laplacian of f.
real(R8P) :: exact_div(      1:n, 1:n, 1:n) !< Exact divergence of v.

! results (device layout)
real(R8P) :: grad_t(1:n, 1:n, 1:n, 1:3) !< Gradient result (i,j,k,nv).
real(R8P) :: lap_t( 1:n, 1:n, 1:n)      !< Laplacian result.
real(R8P) :: div_t( 1:n, 1:n, 1:n)      !< Divergence result.

! device pointers
real(R8P), pointer :: f_dev(:,:,:)      !< Device scalar field.
real(R8P), pointer :: v_dev(:,:,:,:)    !< Device vector field.
real(R8P), pointer :: grad_dev(:,:,:,:) !< Device gradient result.
real(R8P), pointer :: lap_dev(:,:,:)    !< Device laplacian result.
real(R8P), pointer :: div_dev(:,:,:)    !< Device divergence result.

real(R8P)    :: errL2, errL0, e
integer(I4P) :: s, o, i, j, k, d, ierr

! initialize fundal device environment
myhos   = dev_get_host_num()
devtype = dev_get_device_type()
call dev_set_device_num(0)
mydev = dev_get_device_num()

! initialize host data
do k=1-gc, n+gc
   do j=1-gc, n+gc
      do i=1-gc, n+gc
         f_t(i,j,k)     = sin(i*h) + sin(j*h) + sin(k*h)
         v_t(i,j,k,1)   = sin(i*h)
         v_t(i,j,k,2)   = sin(j*h)
         v_t(i,j,k,3)   = sin(k*h)
         if (i>=1.and.j>=1.and.k>=1.and.i<=n.and.j<=n.and.k<=n) then
            exact_grad(1,i,j,k) = cos(i*h)
            exact_grad(2,i,j,k) = cos(j*h)
            exact_grad(3,i,j,k) = cos(k*h)
            exact_lap(i,j,k)    = -sin(i*h) - sin(j*h) - sin(k*h)
            exact_div(i,j,k)    = cos(i*h) + cos(j*h) + cos(k*h)
         endif
      enddo
   enddo
enddo

! allocate device arrays
call dev_alloc(fptr_dev=f_dev,    lbounds=[1-gc,1-gc,1-gc],   ubounds=[n+gc,n+gc,n+gc],   ierr=ierr)
call dev_alloc(fptr_dev=v_dev,    lbounds=[1-gc,1-gc,1-gc,1], ubounds=[n+gc,n+gc,n+gc,3], ierr=ierr)
call dev_alloc(fptr_dev=grad_dev, lbounds=[1,1,1,1],          ubounds=[n,n,n,3],          ierr=ierr)
call dev_alloc(fptr_dev=lap_dev,  lbounds=[1,1,1],            ubounds=[n,n,n],            ierr=ierr)
call dev_alloc(fptr_dev=div_dev,  lbounds=[1,1,1],            ubounds=[n,n,n],            ierr=ierr)

! H2D
call dev_memcpy_to_device(dst=f_dev, src=f_t)
call dev_memcpy_to_device(dst=v_dev, src=v_t)

print '(A)', 'FV gradient (FNL device backend)'
do o=2, 10, 2
   s = o/2
   call gradient_fv_dev_kernel(s=s, ni=n, gc_=gc, f_dev=f_dev, grad_dev=grad_dev)
   call dev_memcpy_from_device(dst=grad_t, src=grad_dev)
   errL2 = 0._R8P ; errL0 = 0._R8P
   do k=1, n ; do j=1, n ; do i=1, n ; do d=1, 3
      e = abs(grad_t(i,j,k,d) - exact_grad(d,i,j,k))
      errL2 = errL2 + e*e ; errL0 = max(errL0, e)
   enddo ; enddo ; enddo ; enddo
   errL2 = sqrt(errL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(errL0))//', Error L2 '//trim(str(errL2))
enddo

print '(A)', 'FV divergence (FNL device backend)'
do o=2, 10, 2
   s = o/2
   call divergence_fv_dev_kernel(s=s, ni=n, gc_=gc, v_dev=v_dev, div_dev=div_dev)
   call dev_memcpy_from_device(dst=div_t, src=div_dev)
   errL2 = 0._R8P ; errL0 = 0._R8P
   do k=1, n ; do j=1, n ; do i=1, n
      e = abs(div_t(i,j,k) - exact_div(i,j,k))
      errL2 = errL2 + e*e ; errL0 = max(errL0, e)
   enddo ; enddo ; enddo
   errL2 = sqrt(errL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(errL0))//', Error L2 '//trim(str(errL2))
enddo

print '(A)', 'FV laplacian (FNL device backend, rebuilt FV d2)'
do o=2, 10, 2
   s = o/2
   call laplacian_fv_dev_kernel(s=s, ni=n, gc_=gc, f_dev=f_dev, lap_dev=lap_dev)
   call dev_memcpy_from_device(dst=lap_t, src=lap_dev)
   errL2 = 0._R8P ; errL0 = 0._R8P
   do k=1, n ; do j=1, n ; do i=1, n
      e = abs(lap_t(i,j,k) - exact_lap(i,j,k))
      errL2 = errL2 + e*e ; errL0 = max(errL0, e)
   enddo ; enddo ; enddo
   errL2 = sqrt(errL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(errL0))//', Error L2 '//trim(str(errL2))
enddo

! free device arrays
call dev_free(f_dev,    mydev)
call dev_free(v_dev,    mydev)
call dev_free(grad_dev, mydev)
call dev_free(lap_dev,  mydev)
call dev_free(div_dev,  mydev)

contains
   subroutine gradient_fv_dev_kernel(s, ni, gc_, f_dev, grad_dev)
   !< Device kernel: FV centered gradient of a scalar field via gathered directional stencils.
   integer(I4P), intent(in)    :: s                            !< Half stencil length.
   integer(I4P), intent(in)    :: ni                           !< Interior cells per direction.
   integer(I4P), intent(in)    :: gc_                          !< Ghost cells.
   real(R8P),    intent(in)    :: f_dev(   1-gc_:,1-gc_:,1-gc_:) !< Device scalar field.
   real(R8P),    intent(inout) :: grad_dev(1:,    1:,    1:,1:)  !< Device gradient (i,j,k,nv).
   integer(I4P) :: i, j, k, ist
   real(R8P) :: qsx(1-s:1+s), qsy(1-s:1+s), qsz(1-s:1+s) !< Directional stencils.
   real(R8P) :: dxyz(3) = [h,h,h]
   real(R8P) :: grad_(3)
   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(f_dev, grad_dev)                        &
   !$acc& firstprivate(s, ni, dxyz)                         &
   !$acc& private(qsx, qsy, qsz, grad_)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(f_dev, grad_dev) &
   !$omp& firstprivate(s, ni, dxyz) &
   !$omp& private(qsx, qsy, qsz, grad_)
   do k=1, ni
   do j=1, ni
   do i=1, ni
      !$acc loop seq
      do ist=1-s, 1+s
         qsx(ist) = f_dev(i+ist-1, j,        k       )
         qsy(ist) = f_dev(i,        j+ist-1, k       )
         qsz(ist) = f_dev(i,        j,        k+ist-1)
      enddo
      call compute_gradient_fv_centered_dev(s=s, dxyz=dxyz, qsx=qsx, qsy=qsy, qsz=qsz, gradient=grad_)
      grad_dev(i,j,k,1) = grad_(1)
      grad_dev(i,j,k,2) = grad_(2)
      grad_dev(i,j,k,3) = grad_(3)
   enddo
   enddo
   enddo
   endsubroutine gradient_fv_dev_kernel

   subroutine divergence_fv_dev_kernel(s, ni, gc_, v_dev, div_dev)
   !< Device kernel: FV centered divergence of a vector field via gathered component stencils.
   integer(I4P), intent(in)    :: s                              !< Half stencil length.
   integer(I4P), intent(in)    :: ni                             !< Interior cells per direction.
   integer(I4P), intent(in)    :: gc_                            !< Ghost cells.
   real(R8P),    intent(in)    :: v_dev( 1-gc_:,1-gc_:,1-gc_:,1:) !< Device vector field (i,j,k,nv).
   real(R8P),    intent(inout) :: div_dev(1:,    1:,    1:)       !< Device divergence.
   integer(I4P) :: i, j, k, ist
   real(R8P) :: qsx(1-s:1+s), qsy(1-s:1+s), qsz(1-s:1+s) !< Component-x over x, y over y, z over z.
   real(R8P) :: dxyz(3) = [h,h,h]
   real(R8P) :: div_
   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(v_dev, div_dev)                         &
   !$acc& firstprivate(s, ni, dxyz)                         &
   !$acc& private(qsx, qsy, qsz)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(v_dev, div_dev) &
   !$omp& firstprivate(s, ni, dxyz) &
   !$omp& private(qsx, qsy, qsz)
   do k=1, ni
   do j=1, ni
   do i=1, ni
      !$acc loop seq
      do ist=1-s, 1+s
         qsx(ist) = v_dev(i+ist-1, j,        k,        1)
         qsy(ist) = v_dev(i,        j+ist-1, k,        2)
         qsz(ist) = v_dev(i,        j,        k+ist-1, 3)
      enddo
      call compute_divergence_fv_centered_dev(s=s, dxyz=dxyz, qsx=qsx, qsy=qsy, qsz=qsz, divergence=div_)
      div_dev(i,j,k) = div_
   enddo
   enddo
   enddo
   endsubroutine divergence_fv_dev_kernel

   subroutine laplacian_fv_dev_kernel(s, ni, gc_, f_dev, lap_dev)
   !< Device kernel: FV centered laplacian of a scalar field (exercises the rebuilt FV d2).
   integer(I4P), intent(in)    :: s                            !< Half stencil length.
   integer(I4P), intent(in)    :: ni                           !< Interior cells per direction.
   integer(I4P), intent(in)    :: gc_                          !< Ghost cells.
   real(R8P),    intent(in)    :: f_dev(  1-gc_:,1-gc_:,1-gc_:) !< Device scalar field.
   real(R8P),    intent(inout) :: lap_dev(1:,    1:,    1:)     !< Device laplacian.
   integer(I4P) :: i, j, k, ist
   real(R8P) :: qsx(1-s:1+s), qsy(1-s:1+s), qsz(1-s:1+s) !< Directional stencils.
   real(R8P) :: dxyz(3) = [h,h,h]
   real(R8P) :: lap_
   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(f_dev, lap_dev)                         &
   !$acc& firstprivate(s, ni, dxyz)                         &
   !$acc& private(qsx, qsy, qsz)
   !$omp OMPLOOP collapse(3) &
   !$omp& DEVICEPTR(f_dev, lap_dev) &
   !$omp& firstprivate(s, ni, dxyz) &
   !$omp& private(qsx, qsy, qsz)
   do k=1, ni
   do j=1, ni
   do i=1, ni
      !$acc loop seq
      do ist=1-s, 1+s
         qsx(ist) = f_dev(i+ist-1, j,        k       )
         qsy(ist) = f_dev(i,        j+ist-1, k       )
         qsz(ist) = f_dev(i,        j,        k+ist-1)
      enddo
      call compute_laplacian_fv_centered_dev(s=s, dxyz=dxyz, qsx=qsx, qsy=qsy, qsz=qsz, laplacian=lap_)
      lap_dev(i,j,k) = lap_
   enddo
   enddo
   enddo
   endsubroutine laplacian_fv_dev_kernel
endprogram test_fdv_fnl_operators_vector_trigonometric
