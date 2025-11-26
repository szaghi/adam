program test_fdv_operators_step
!< Test ADAM fdv operators (compute derivatives) library with step function.

!< Step function
!<```
!<  q(x) = exp(bx)*(1 + tanh((x-x0)/a)), x in [0,1], x=0.5, a<1 (0.01)
!<  dq(x)/dx = b*exp(bx)*(1 + tanh((x-x0)/a)) + (exp(bx)/a)*sech^2((x-x0)/a)
!<```

use adam_fdv_operators_library
use penf

implicit none

integer(I4P), parameter :: n=100_I4P           !< Number of cells.
integer(I4P), parameter :: gc=5_I4P            !< Number of ghost cells.
real(R8P),    parameter :: dx=1._R8P/n         !< Space step.
real(R8P),    parameter :: x0=0.5_R8P          !< Pivot.
real(R8P),    parameter ::  a=0.01_R8P         !< Steepness.
real(R8P),    parameter ::  b=5._R8P           !< Transport speed.
real(R8P)               ::    x(1-gc:n+gc)     !< Abscissa.
real(R8P)               ::    q(1-gc:n+gc)     !< Function.
real(R8P)               :: dqdx( 1:n,1:10,1:4) !< Derivative 1 of function for each scheme available.
real(R8P)               :: exact(1:n)          !< Exact derivative.
real(R8P)               :: error(1:n,1:10,1:4) !< Error functions for each scheme available.
real(R8P)               :: errorL2             !< L2 error norm.
real(R8P)               :: y                   !< Buffer.
integer(I4P)            :: s_d(1:10,1:4)       !< Half stencil length of each derivative (depends from order of accuracy).
character(16)           :: msg_head(1:10,1:4)  !< Print message header.
character(23)           :: file_name(1:10,1:4) !< File name.
integer(I4P)            :: i, o, s             !< Counter.

! initialize
do i=1-gc, n+gc
   x(i) = i*dx
   y = (x(i)-x0)/a
   ! q(i) = 1._R8P + tanh(y)
   ! if (i>=1.and.i<=n) exact(i) = (1._R8P/a)*(1._R8P/cosh(y)**2)
   q(i) = exp(b*x(i))*(1._R8P + tanh(y))
   if (i>=1.and.i<=n) exact(i) = b*exp(b*x(i))*(1._R8P + tanh(y)) + (exp(b*x(i))/a)*(1._R8P/cosh(y)**2)
enddo
do o=2, 10, 2
   ! FD                 ! FV centered        ! FV right-upwind   ! FV left-upwind
   s_d(o,1) = o/2 + 0 ; s_d(o,2) = o/2 + 0 ; s_d(o,3) = o + 0  ; s_d(o,4) = o + 0
   msg_head(o,1) = '  FD CC order '//trim(strz(o,2))
   msg_head(o,2) = '  FV CC order '//trim(strz(o,2))
   msg_head(o,3) = '  FV RU order '//trim(strz(o,2))
   msg_head(o,4) = '  FV LU order '//trim(strz(o,2))
   file_name(o,1) = 'FD_CC_order_'//trim(strz(o,2))//'-step.dat'
   file_name(o,2) = 'FV_CC_order_'//trim(strz(o,2))//'-step.dat'
   file_name(o,3) = 'FV_RU_order_'//trim(strz(o,2))//'-step.dat'
   file_name(o,4) = 'FV_LU_order_'//trim(strz(o,2))//'-step.dat'
enddo

! test FDV operators
do o=2, 10, 2
   do i=1, n
      if (s_d(o,1)<=gc) &
      call compute_derivative1_fd_centered(s=s_d(o,1),ds=dx,q=q(i-s_d(o,1):i+s_d(o,1)),dq_ds=dqdx(i,o,1))
      if (s_d(o,2)<=gc) &
      call compute_derivative1_fv_centered(s=s_d(o,2),ds=dx,q=q(i-s_d(o,2):i+s_d(o,2)),dq_ds=dqdx(i,o,2))
      if (s_d(o,3)<=gc) &
      call compute_derivative1_fv_rupwind( s=s_d(o,3),ds=dx,q=q(i         :i+s_d(o,3)),dq_ds=dqdx(i,o,3))
      if (s_d(o,4)<=gc) &
      call compute_derivative1_fv_lupwind( s=s_d(o,4),ds=dx,q=q(i-s_d(o,4):i         ),dq_ds=dqdx(i,o,4))
   enddo
   do s=1, 4
      error(:,o,s) = abs(dqdx(:,o,s) - exact(:))
   enddo
enddo

print '(A)', 'Step function derivative 1'
do o=2, 10, 2
   do s=1, 4
      if (s_d(o,s)<=gc) then
         errorL2 = 0._R8P
         do i=1, n
            errorL2 = errorL2 + error(i,o,s)*error(i,o,s)
         enddo
         errorL2 = sqrt(errorL2*dx*dx)
         print '(A)',msg_head(o,s)//': Error L0 '//trim(str(maxval(error(:,o,s))))//', Error L2 '//trim(str(errorL2))
         call save_file(filename=file_name(o,s),fn=dqdx(:,o,s),fe=exact(:))
      endif
   enddo
enddo

contains
   subroutine save_file(filename, fn, fe)
   !< Save results on file.
   character(*), intent(in) :: filename   !< File name.
   real(R8P),    intent(in) :: fn(1:)     !< Numerical function.
   real(R8P),    intent(in) :: fe(1:)     !< Exact function.
   integer(I4P)             :: j, fu      !< Counter.

   open(newunit=fu, file=trim(adjustl(filename)))
   write(fu, '(A)') 'x numerical_dqdx exact_dqdx'
   do j=1, n
      write(fu, '(A)') trim(str(n=[x(j),fn(j),fe(j)],separator=' '))
   enddo
   close(fu)
   endsubroutine save_file
endprogram test_fdv_operators_step
