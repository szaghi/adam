program test_fdv_operators_trigonometric
!< Test ADAM fdv operators (compute derivatives) library with trigonometric functions derivatives.

use adam_fdv_operators_library
use penf

implicit none

integer(I4P), parameter :: n=100_I4P                !< Number of cells.
integer(I4P), parameter :: gc=5_I4P                 !< Number of ghost cells.
real(R8P),    parameter :: PI=4._R8P*atan(1._R8P)   !< PI.
real(R8P),    parameter :: dx=2_I4P*PI/n            !< Space step.
real(R8P)               ::    x(1-gc:n+gc)          !< Abscissa.
real(R8P)               :: sin0(1-gc:n+gc)          !< Function sin.
real(R8P)               :: sind( 1:n,1:6,1:10,1:4)  !< Derivative d-th of function sin for each scheme available.
real(R8P)               :: exact(1:n,1:6         )  !< Exact derivatives.
real(R8P)               :: error(1:n,1:6,1:10,1:4)  !< Error functions for each scheme available.
real(R8P)               :: errorL2                  !< L2 error norm.
integer(I4P)            :: s_d(1:6,1:10,1:4)        !< Half stencil length of each derivative (depends from order of accuracy).
character(16)           :: msg_head(1:10,1:4)       !< Print message header.
character(12)           :: file_header(1:6)         !< File header.
character(24)           :: file_name(1:6,1:10,1:4)  !< File name.
integer(I4P)            :: i, o, d, s               !< Counter.

! initialize
do i=1-gc, n+gc
   x(i) = i*dx
   sin0(i)    =  sin(x(i))
   if (i>=1.and.i<=n) then
      exact(i,1) =  cos(x(i))
      exact(i,2) = -sin(x(i))
      exact(i,3) = -cos(x(i))
      exact(i,4) =  sin(x(i))
      exact(i,5) =  cos(x(i))
      exact(i,6) = -sin(x(i))
   endif
enddo
do o=2, 10, 2
   ! FD                   ! FV centered          ! FV right-upwind    ! FV left-upwind
   s_d(1,o,1) = o/2 + 0 ; s_d(1,o,2) = o/2 + 0 ; s_d(1,o,3) = o + 0 ; s_d(1,o,4) = o + 0
   s_d(2,o,1) = o/2 + 0 ; s_d(2,o,2) = o/2 + 1 ; s_d(2,o,3) = o + 1 ; s_d(2,o,4) = o + 1
   s_d(3,o,1) = o/2 + 1 ; s_d(3,o,2) = o/2 + 2 ; s_d(3,o,3) = o + 2 ; s_d(3,o,4) = o + 2
   s_d(4,o,1) = o/2 + 1 ; s_d(4,o,2) = o/2 + 3 ; s_d(4,o,3) = o + 3 ; s_d(4,o,4) = o + 3
   s_d(5,o,1) = o/2 + 2 ; s_d(5,o,2) = o/2 + 4 ; s_d(5,o,3) = o + 4 ; s_d(5,o,4) = o + 4
   s_d(6,o,1) = o/2 + 2 ; s_d(6,o,2) = o/2 + 5 ; s_d(6,o,3) = o + 5 ; s_d(6,o,4) = o + 5
   msg_head(o,1) = '  FD CC order '//trim(strz(o,2))
   msg_head(o,2) = '  FV CC order '//trim(strz(o,2))
   msg_head(o,3) = '  FV RU order '//trim(strz(o,2))
   msg_head(o,4) = '  FV LU order '//trim(strz(o,2))
   do d=1, 6
      file_name(d,o,1) = 'FD_CC_order_'//trim(strz(o,2))//'-d'//trim(strz(d,1))//'sin.dat'
      file_name(d,o,2) = 'FV_CC_order_'//trim(strz(o,2))//'-d'//trim(strz(d,1))//'sin.dat'
      file_name(d,o,3) = 'FV_RU_order_'//trim(strz(o,2))//'-d'//trim(strz(d,1))//'sin.dat'
      file_name(d,o,4) = 'FV_LU_order_'//trim(strz(o,2))//'-d'//trim(strz(d,1))//'sin.dat'
   enddo
enddo
file_header(1) = 'x d1sin cos '
file_header(2) = 'x d2sin -sin'
file_header(3) = 'x d3sin -cos'
file_header(4) = 'x d4sin sin '
file_header(5) = 'x d5sin cos '
file_header(6) = 'x d6sin -sin'

! test FDV operators
do o=2, 10, 2
   do i=1, n
      if (s_d(1,o,1)<=gc) &
      call compute_derivative1_fd_centered(s=s_d(1,o,1),ds=dx,q=sin0(i-s_d(1,o,1):i+s_d(1,o,1)),dq_ds  =sind(i,1,o,1))
      if (s_d(2,o,1)<=gc) &
      call compute_derivative2_fd_centered(s=s_d(2,o,1),ds=dx,q=sin0(i-s_d(2,o,1):i+s_d(2,o,1)),d2q_ds2=sind(i,2,o,1))
      if (s_d(3,o,1)<=gc) &
      call compute_derivative3_fd_centered(s=s_d(3,o,1),ds=dx,q=sin0(i-s_d(3,o,1):i+s_d(3,o,1)),d3q_ds3=sind(i,3,o,1))
      if (s_d(4,o,1)<=gc) &
      call compute_derivative4_fd_centered(s=s_d(4,o,1),ds=dx,q=sin0(i-s_d(4,o,1):i+s_d(4,o,1)),d4q_ds4=sind(i,4,o,1))
      if (s_d(5,o,1)<=gc) &
      call compute_derivative5_fd_centered(s=s_d(5,o,1),ds=dx,q=sin0(i-s_d(5,o,1):i+s_d(5,o,1)),d5q_ds5=sind(i,5,o,1))
      if (s_d(6,o,1)<=gc) &
      call compute_derivative6_fd_centered(s=s_d(6,o,1),ds=dx,q=sin0(i-s_d(6,o,1):i+s_d(6,o,1)),d6q_ds6=sind(i,6,o,1))

      if (s_d(1,o,2)<=gc) &
      call compute_derivative1_fv_centered(s=s_d(1,o,2),ds=dx,q=sin0(i-s_d(1,o,2):i+s_d(1,o,2)),dq_ds  =sind(i,1,o,2))
      if (s_d(2,o,2)<=gc) &
      call compute_derivative2_fv_centered(s=s_d(2,o,2),ds=dx,q=sin0(i-s_d(2,o,2):i+s_d(2,o,2)),d2q_ds2=sind(i,2,o,2))
      if (s_d(3,o,2)<=gc) &
      call compute_derivative3_fv_centered(s=s_d(3,o,2),ds=dx,q=sin0(i-s_d(3,o,2):i+s_d(3,o,2)),d3q_ds3=sind(i,3,o,2))
      if (s_d(4,o,2)<=gc) &
      call compute_derivative4_fv_centered(s=s_d(4,o,2),ds=dx,q=sin0(i-s_d(4,o,2):i+s_d(4,o,2)),d4q_ds4=sind(i,4,o,2))
      if (s_d(5,o,2)<=gc) &
      call compute_derivative5_fv_centered(s=s_d(5,o,2),ds=dx,q=sin0(i-s_d(5,o,2):i+s_d(5,o,2)),d5q_ds5=sind(i,5,o,2))
      if (s_d(6,o,2)<=gc) &
      call compute_derivative6_fv_centered(s=s_d(6,o,2),ds=dx,q=sin0(i-s_d(6,o,2):i+s_d(6,o,2)),d6q_ds6=sind(i,6,o,2))

      if (s_d(1,o,3)<=gc) &
      call compute_derivative1_fv_rupwind( s=s_d(1,o,3),ds=dx,q=sin0(i           :i+s_d(1,o,3)),dq_ds  =sind(i,1,o,3))
      if (s_d(2,o,3)<=gc) &
      call compute_derivative2_fv_rupwind( s=s_d(2,o,3),ds=dx,q=sin0(i           :i+s_d(2,o,3)),d2q_ds2=sind(i,2,o,3))
      if (s_d(3,o,3)<=gc) &
      call compute_derivative3_fv_rupwind( s=s_d(3,o,3),ds=dx,q=sin0(i           :i+s_d(3,o,3)),d3q_ds3=sind(i,3,o,3))
      if (s_d(4,o,3)<=gc) &
      call compute_derivative4_fv_rupwind( s=s_d(4,o,3),ds=dx,q=sin0(i           :i+s_d(4,o,3)),d4q_ds4=sind(i,4,o,3))
      if (s_d(5,o,3)<=gc) &
      call compute_derivative5_fv_rupwind( s=s_d(5,o,3),ds=dx,q=sin0(i           :i+s_d(5,o,3)),d5q_ds5=sind(i,5,o,3))
      if (s_d(6,o,3)<=gc) &
      call compute_derivative6_fv_rupwind( s=s_d(6,o,3),ds=dx,q=sin0(i           :i+s_d(6,o,3)),d6q_ds6=sind(i,6,o,3))

      if (s_d(1,o,4)<=gc) &
      call compute_derivative1_fv_lupwind( s=s_d(1,o,4),ds=dx,q=sin0(i-s_d(1,o,4):i           ),dq_ds  =sind(i,1,o,4))
      if (s_d(2,o,4)<=gc) &
      call compute_derivative2_fv_lupwind( s=s_d(2,o,4),ds=dx,q=sin0(i-s_d(2,o,4):i           ),d2q_ds2=sind(i,2,o,4))
      if (s_d(3,o,4)<=gc) &
      call compute_derivative3_fv_lupwind( s=s_d(3,o,4),ds=dx,q=sin0(i-s_d(3,o,4):i           ),d3q_ds3=sind(i,3,o,4))
      if (s_d(4,o,4)<=gc) &
      call compute_derivative4_fv_lupwind( s=s_d(4,o,4),ds=dx,q=sin0(i-s_d(4,o,4):i           ),d4q_ds4=sind(i,4,o,4))
      if (s_d(5,o,4)<=gc) &
      call compute_derivative5_fv_lupwind( s=s_d(5,o,4),ds=dx,q=sin0(i-s_d(5,o,4):i           ),d5q_ds5=sind(i,5,o,4))
      if (s_d(6,o,4)<=gc) &
      call compute_derivative6_fv_lupwind( s=s_d(6,o,4),ds=dx,q=sin0(i-s_d(6,o,4):i           ),d6q_ds6=sind(i,6,o,4))
   enddo
   do s=1, 4
      do d=1, 6
         error(:,d,o,s) = abs(sind(:,d,o,s) - exact(:,d))
      enddo
   enddo
enddo

do d=1, 6
   print '(A)', 'Derivative '//trim(str(d,.true.))
   do o=2, 10, 2
      do s=1, 4
         if (s_d(d,o,s)<=gc) then
            errorL2 = 0._R8P
            do i=1, n
               errorL2 = errorL2 + error(i,d,o,s)*error(i,d,o,s)
            enddo
            errorL2 = sqrt(errorL2*dx*dx)
            print '(A)',msg_head(o,s)//': Error L0 '//trim(str(maxval(error(:,d,o,s))))//', Error L2 '//trim(str(errorL2))
            call save_file(filename=file_name(d,o,s),fileheader=file_header(d),fn=sind(:,d,o,s),fe=exact(:,d))
         endif
      enddo
   enddo
enddo

contains
   subroutine save_file(filename, fileheader, fn, fe)
   !< Save results on file.
   character(*), intent(in) :: filename   !< File name.
   character(*), intent(in) :: fileheader !< Header of file.
   real(R8P),    intent(in) :: fn(1:)     !< Numerical function.
   real(R8P),    intent(in) :: fe(1:)     !< Exact function.
   integer(I4P)             :: j, fu      !< Counter.

   open(newunit=fu, file=trim(adjustl(filename)))
   write(fu, '(A)') trim(adjustl(fileheader))
   do j=1, n
      write(fu, '(A)') trim(str(n=[x(j),fn(j),fe(j)],separator=' '))
   enddo
   close(fu)
   endsubroutine save_file
endprogram test_fdv_operators_trigonometric
