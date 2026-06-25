program test_fdv_fnl_operators_trigonometric
!< Test low-level FD centered derivative routines (compute_derivative1_fd_centered, ...)
!< called from OpenACC device loops — FNL device backend equivalent of test_fdv_operators_trigonometric.F90.
!<
!< Test function: sin(x), exact derivatives:
!<   d1 sin = cos(x), d2 sin = -sin(x), d3 sin = -cos(x),
!<   d4 sin = sin(x), d5 sin = cos(x),  d6 sin = -sin(x)
!<
!< Workflow: initialize sin(x) on host → H2D via fundal dev_memcpy_to_device →
!<           OpenACC kernel calls compute_derivativeN_fd_centered (!$acc routine seq) →
!<           D2H via fundal dev_memcpy_from_device → compare with exact.
!<
!< Only FD centered scheme is tested, matching the routines used in the FNL backend.

#include "fundal.H"

use :: adam_fdv_operators_library
use :: fundal
use :: penf

implicit none

integer(I4P), parameter :: n=100_I4P              !< Number of cells.
integer(I4P), parameter :: gc=5_I4P               !< Number of ghost cells (= FDV_S_MAX).
real(R8P),    parameter :: PI=4._R8P*atan(1._R8P) !< PI.
real(R8P),    parameter :: dx=2_I4P*PI/n          !< Space step.

! host arrays
real(R8P) :: x(   1-gc:n+gc)  !< Abscissa.
real(R8P) :: sin0(1-gc:n+gc)  !< sin(x) with ghost cells.
real(R8P) :: sind(1:n)        !< Numerical derivative (one (d,s) pair at a time).
real(R8P) :: exact(1:n,1:6)   !< Exact derivatives d1..d6.
real(R8P) :: error(1:n)       !< Point-wise absolute error.
real(R8P) :: errorL2          !< L2 error norm.

! device pointers (allocated via fundal dev_alloc)
real(R8P), pointer :: sin0_dev(:) => null() !< Device sin(x) field.
real(R8P), pointer :: dsin_dev(:) => null() !< Device derivative result (reused per call).

! stencil half-widths per (derivative, order), same convention as CPU reference test
! s_d(d, o) = minimum stencil half-width for derivative d at FD order o
integer(I4P) :: s_d(1:6, 2:10)
integer(I4P) :: i, o, d, ierr

! initialize fundal device environment
myhos   = dev_get_host_num()
devtype = dev_get_device_type()
call dev_set_device_num(0)
mydev   = dev_get_device_num()

! initialize host data
do i=1-gc, n+gc
   x(i) = i*dx
   sin0(i) = sin(x(i))
enddo
do i=1, n
   exact(i,1) =  cos(x(i))
   exact(i,2) = -sin(x(i))
   exact(i,3) = -cos(x(i))
   exact(i,4) =  sin(x(i))
   exact(i,5) =  cos(x(i))
   exact(i,6) = -sin(x(i))
enddo

! stencil half-widths (FD column from CPU reference test)
do o=2, 10, 2
   s_d(1,o) = o/2     ; s_d(2,o) = o/2
   s_d(3,o) = o/2 + 1 ; s_d(4,o) = o/2 + 1
   s_d(5,o) = o/2 + 2 ; s_d(6,o) = o/2 + 2
enddo

! allocate device arrays
call dev_alloc(fptr_dev=sin0_dev, lbounds=[1-gc], ubounds=[n+gc], ierr=ierr)
call dev_alloc(fptr_dev=dsin_dev, lbounds=[1],    ubounds=[n],    ierr=ierr)

! H2D: copy sin(x) to device (single copy, reused for all (d,o) combinations)
call dev_memcpy_to_device(dst=sin0_dev, src=sin0)

! test FD centered derivatives on device
do d=1, 6
   print '(A)', 'Derivative '//trim(str(d,.true.))
   do o=2, 10, 2
      if (s_d(d,o) > gc) cycle
      select case(d)
         case(1) ; call kernel_d1(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(2) ; call kernel_d2(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(3) ; call kernel_d3(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(4) ; call kernel_d4(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(5) ; call kernel_d5(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(6) ; call kernel_d6(s_d(d,o), n, dx, sin0_dev, dsin_dev)
      endselect
      ! D2H: copy result back
      call dev_memcpy_from_device(dst=sind, src=dsin_dev)
      ! compute error
      error = abs(sind - exact(:,d))
      errorL2 = sqrt(sum(error*error)*dx*dx)
      print '(A)', '  FD CC order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
   enddo
enddo

! test rebuilt FV centered derivative ladder on device (d2..d6 via face-reconstruction + flux difference).
! Same gather as FD (window q(1-s:1+s)); these exercise the new compute_reconstruction_r_fv_centered_d* kernels.
do d=2, 6
   print '(A)', 'FV centered derivative '//trim(str(d,.true.))
   do o=2, 10, 2
      if (s_d(d,o) > gc) cycle
      select case(d)
         case(2) ; call kernel_fv_d2(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(3) ; call kernel_fv_d3(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(4) ; call kernel_fv_d4(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(5) ; call kernel_fv_d5(s_d(d,o), n, dx, sin0_dev, dsin_dev)
         case(6) ; call kernel_fv_d6(s_d(d,o), n, dx, sin0_dev, dsin_dev)
      endselect
      call dev_memcpy_from_device(dst=sind, src=dsin_dev)
      error = abs(sind - exact(:,d))
      errorL2 = sqrt(sum(error*error)*dx*dx)
      print '(A)', '  FV CC order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
   enddo
enddo

! test deconvolution (cell-average -> cell-center point value) on device.
print '(A)', 'FV deconvolution avg->point (reads sin as cell averages, recovers point value)'
do o=2, 10, 2
   if (o/2 > gc) cycle
   call kernel_deconv_a2p(o/2, n, sin0_dev, dsin_dev)
   call dev_memcpy_from_device(dst=sind, src=dsin_dev)
   ! exact point value at cell center x_i is sin(x_i); input sin0 here is treated as the cell-average
   ! sample, so the recovered point value carries the O(h^2) average<->point gap unless sin0 were true
   ! averages -- this checks the device kernel reads DECONV_A2P and runs, not asymptotic order.
   do i=1, n
      error(i) = abs(sind(i) - sin(x(i)))
   enddo
   errorL2 = sqrt(sum(error*error)*dx*dx)
   print '(A)', '  DECONV order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
enddo

! free device arrays
call dev_free(sin0_dev, mydev)
call dev_free(dsin_dev, mydev)

contains
   subroutine kernel_d1(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 1st FD centered derivative at all interior points.
   !< Stencil gathered into private qs to avoid passing sections to acc routine seq.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative1_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), dq_ds=dsin_dev(i))
   enddo
   endsubroutine kernel_d1

   subroutine kernel_d2(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 2nd FD centered derivative at all interior points.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative2_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d2q_ds2=dsin_dev(i))
   enddo
   endsubroutine kernel_d2

   subroutine kernel_d3(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 3rd FD centered derivative at all interior points.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative3_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d3q_ds3=dsin_dev(i))
   enddo
   endsubroutine kernel_d3

   subroutine kernel_d4(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 4th FD centered derivative at all interior points.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative4_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d4q_ds4=dsin_dev(i))
   enddo
   endsubroutine kernel_d4

   subroutine kernel_d5(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 5th FD centered derivative at all interior points.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative5_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d5q_ds5=dsin_dev(i))
   enddo
   endsubroutine kernel_d5

   subroutine kernel_d6(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: compute 6th FD centered derivative at all interior points.
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative6_fd_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d6q_ds6=dsin_dev(i))
   enddo
   endsubroutine kernel_d6

   subroutine kernel_fv_d2(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: rebuilt 2nd FV centered derivative (flux diff of FVR1 face reconstruction).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative2_fv_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d2q_ds2=dsin_dev(i))
   enddo
   endsubroutine kernel_fv_d2

   subroutine kernel_fv_d3(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: rebuilt 3rd FV centered derivative (flux diff of FVR2 face reconstruction).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative3_fv_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d3q_ds3=dsin_dev(i))
   enddo
   endsubroutine kernel_fv_d3

   subroutine kernel_fv_d4(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: rebuilt 4th FV centered derivative (flux diff of FVR3 face reconstruction).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative4_fv_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d4q_ds4=dsin_dev(i))
   enddo
   endsubroutine kernel_fv_d4

   subroutine kernel_fv_d5(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: rebuilt 5th FV centered derivative (flux diff of FVR4 face reconstruction).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative5_fv_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d5q_ds5=dsin_dev(i))
   enddo
   endsubroutine kernel_fv_d5

   subroutine kernel_fv_d6(s_, n_, dx_, sin0_dev, dsin_dev)
   !< OpenACC kernel: rebuilt 6th FV centered derivative (flux diff of FVR5 face reconstruction).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: dx_                         !< Space step.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_, dx_)                &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_, dx_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_derivative6_fv_centered(s=s_, ds=dx_, q=qs(1-s_:1+s_), d6q_ds6=dsin_dev(i))
   enddo
   endsubroutine kernel_fv_d6

   subroutine kernel_deconv_a2p(s_, n_, sin0_dev, dsin_dev)
   !< OpenACC kernel: deconvolution cell-average -> cell-center point value (reads DECONV_A2P on device).
   integer(I4P), intent(in)    :: s_                          !< Stencil half-width.
   integer(I4P), intent(in)    :: n_                          !< Number of interior cells.
   real(R8P),    intent(in)    :: sin0_dev(1-gc:)             !< Device field.
   real(R8P),    intent(inout) :: dsin_dev(1:)                !< Device result.
   real(R8P)                   :: qs(1-FDV_S_MAX:1+FDV_S_MAX) !< Private stencil buffer.
   integer(I4P)                :: i, ist                      !< Counter.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(sin0_dev, dsin_dev)        &
   !$acc& firstprivate(s_)                     &
   !$acc& private(qs)
   !$omp OMPLOOP &
   !$omp& DEVICEPTR(sin0_dev, dsin_dev) &
   !$omp& firstprivate(s_) &
   !$omp& private(qs)
   do i=1, n_
      !$acc loop seq
      do ist=1-s_, 1+s_
         qs(ist) = sin0_dev(i+ist-1)
      enddo
      call compute_avg_to_point(s=s_, q=qs(1-s_:1+s_), qp=dsin_dev(i))
   enddo
   endsubroutine kernel_deconv_a2p
endprogram test_fdv_fnl_operators_trigonometric
