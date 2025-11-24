program test_fdv_operators_trigonometric
!< Test ADM fdv operators library with trigonometric functions derivatives.

use adam_fdv_operators_library
use penf

implicit none

integer(I4P), parameter :: ni=100                 !< Number of cells.
integer(I4P), parameter :: gc=5                   !< Number of ghost cells.
real(R8P),    parameter :: PI=4._R8P*atan(1._R8P) !< PI.
real(R8P),    parameter :: dx=2*PI/ni             !< Space step.
real(R8P)               ::    x(1-gc:ni+gc)       !< Abscissa.
real(R8P)               :: sin0(1-gc:ni+gc)       !< Function sin.
real(R8P)               :: sin1(1   :ni   )       !< Derivative1 of function sin, sin1=cos.
real(R8P)               :: sin2(1   :ni   )       !< Derivative2 of function sin, sin2=-sin.
real(R8P)               :: sin3(1   :ni   )       !< Derivative3 of function sin, sin2=-cos.
real(R8P)               :: sin4(1   :ni   )       !< Derivative4 of function sin, sin4=sin.
real(R8P)               :: sin5(1   :ni   )       !< Derivative5 of function sin, sin5=con.
real(R8P)               :: sin6(1   :ni   )       !< Derivative6 of function sin, sin6=-sin.
integer(I4P)            :: s_d(6)                 !< Half stencil length of each derivative (depends from order of accuracy).
integer(I4P)            :: i, o                   !< Counter.

! initialize
do i=1-gc, ni+gc
   x(i) = i*dx
   sin0(i) = sin(x(i))
enddo

! test FD operators
do o=2, 6, 2
   s_d(1) = o/2 + 0
   s_d(2) = o/2 + 0
   s_d(3) = o/2 + 1
   s_d(4) = o/2 + 1
   s_d(5) = o/2 + 2
   s_d(6) = o/2 + 2
   do i=1, ni
      call compute_derivative1_fd_centered(s=s_d(1),ds=dx,q=sin0(i-s_d(1):i+s_d(1)),dq_ds  =sin1(i))
      call compute_derivative2_fd_centered(s=s_d(2),ds=dx,q=sin0(i-s_d(2):i+s_d(2)),d2q_ds2=sin2(i))
      call compute_derivative3_fd_centered(s=s_d(3),ds=dx,q=sin0(i-s_d(3):i+s_d(3)),d3q_ds3=sin3(i))
      call compute_derivative4_fd_centered(s=s_d(4),ds=dx,q=sin0(i-s_d(4):i+s_d(4)),d4q_ds4=sin4(i))
      call compute_derivative5_fd_centered(s=s_d(5),ds=dx,q=sin0(i-s_d(5):i+s_d(5)),d5q_ds5=sin5(i))
      call compute_derivative6_fd_centered(s=s_d(6),ds=dx,q=sin0(i-s_d(6):i+s_d(6)),d6q_ds6=sin6(i))
   enddo
   call save_file(file_name='fd_d1sin-ord'//trim(strz(o,2))//'.dat',header='x sin1 cos' ,x=x(1:ni),fn=sin1,fe=[( cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d2sin-ord'//trim(strz(o,2))//'.dat',header='x sin2 -sin',x=x(1:ni),fn=sin2,fe=[(-sin(x(i)),i=1,ni)])
   call save_file(file_name='fd_d3sin-ord'//trim(strz(o,2))//'.dat',header='x sin3 -cos',x=x(1:ni),fn=sin3,fe=[(-cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d4sin-ord'//trim(strz(o,2))//'.dat',header='x sin4 sin' ,x=x(1:ni),fn=sin4,fe=[( sin(x(i)),i=1,ni)])
   call save_file(file_name='fd_d5sin-ord'//trim(strz(o,2))//'.dat',header='x sin5 cos' ,x=x(1:ni),fn=sin5,fe=[( cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d6sin-ord'//trim(strz(o,2))//'.dat',header='x sin6 -sin',x=x(1:ni),fn=sin6,fe=[(-sin(x(i)),i=1,ni)])
enddo
do o=8, 8
   s_d(1) = o/2 + 0
   s_d(2) = o/2 + 0
   s_d(3) = o/2 + 1
   s_d(4) = o/2 + 1
   do i=1, ni
      call compute_derivative1_fd_centered(s=s_d(1),ds=dx,q=sin0(i-s_d(1):i+s_d(1)),dq_ds  =sin1(i))
      call compute_derivative2_fd_centered(s=s_d(2),ds=dx,q=sin0(i-s_d(2):i+s_d(2)),d2q_ds2=sin2(i))
      call compute_derivative3_fd_centered(s=s_d(3),ds=dx,q=sin0(i-s_d(3):i+s_d(3)),d3q_ds3=sin3(i))
      call compute_derivative4_fd_centered(s=s_d(4),ds=dx,q=sin0(i-s_d(4):i+s_d(4)),d4q_ds4=sin4(i))
   enddo
   call save_file(file_name='fd_d1sin-ord'//trim(strz(o,2))//'.dat',header='x sin1 cos' ,x=x(1:ni),fn=sin1,fe=[( cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d2sin-ord'//trim(strz(o,2))//'.dat',header='x sin2 -sin',x=x(1:ni),fn=sin2,fe=[(-sin(x(i)),i=1,ni)])
   call save_file(file_name='fd_d3sin-ord'//trim(strz(o,2))//'.dat',header='x sin3 -cos',x=x(1:ni),fn=sin3,fe=[(-cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d4sin-ord'//trim(strz(o,2))//'.dat',header='x sin4 sin' ,x=x(1:ni),fn=sin4,fe=[( sin(x(i)),i=1,ni)])
enddo
do o=10, 10
   s_d(1) = o/2 + 0
   s_d(2) = o/2 + 0
   do i=1, ni
      call compute_derivative1_fd_centered(s=s_d(1),ds=dx,q=sin0(i-s_d(1):i+s_d(1)),dq_ds  =sin1(i))
      call compute_derivative2_fd_centered(s=s_d(2),ds=dx,q=sin0(i-s_d(2):i+s_d(2)),d2q_ds2=sin2(i))
   enddo
   call save_file(file_name='fd_d1sin-ord'//trim(strz(o,2))//'.dat',header='x sin1 cos' ,x=x(1:ni),fn=sin1,fe=[( cos(x(i)),i=1,ni)])
   call save_file(file_name='fd_d2sin-ord'//trim(strz(o,2))//'.dat',header='x sin2 -sin',x=x(1:ni),fn=sin2,fe=[(-sin(x(i)),i=1,ni)])
enddo

! test FV operators
do o=2, 2, 2
   s_d(1) = o/2 + 0
   s_d(2) = o/2 + 1
   s_d(3) = o/2 + 2
   s_d(4) = o/2 + 3
   s_d(5) = o/2 + 4
   do i=1, ni
      call compute_derivative1_fv_centered(s=s_d(1),ds=dx,q=sin0(i-s_d(1):i+s_d(1)),dq_ds  =sin1(i))
      call compute_derivative2_fv_centered(s=s_d(2),ds=dx,q=sin0(i-s_d(2):i+s_d(2)),d2q_ds2=sin2(i))
      call compute_derivative3_fv_centered(s=s_d(3),ds=dx,q=sin0(i-s_d(3):i+s_d(3)),d3q_ds3=sin3(i))
      call compute_derivative4_fv_centered(s=s_d(4),ds=dx,q=sin0(i-s_d(4):i+s_d(4)),d4q_ds4=sin4(i))
      call compute_derivative5_fv_centered(s=s_d(5),ds=dx,q=sin0(i-s_d(5):i+s_d(5)),d5q_ds5=sin5(i))
   enddo
   call save_file(file_name='fv_d1sin-ord'//trim(strz(o,2))//'.dat',header='x sin1 cos' ,x=x(1:ni),fn=sin1,fe=[( cos(x(i)),i=1,ni)])
   call save_file(file_name='fv_d2sin-ord'//trim(strz(o,2))//'.dat',header='x sin2 -sin',x=x(1:ni),fn=sin2,fe=[(-sin(x(i)),i=1,ni)])
   call save_file(file_name='fv_d3sin-ord'//trim(strz(o,2))//'.dat',header='x sin3 -cos',x=x(1:ni),fn=sin3,fe=[(-cos(x(i)),i=1,ni)])
   call save_file(file_name='fv_d4sin-ord'//trim(strz(o,2))//'.dat',header='x sin4 sin' ,x=x(1:ni),fn=sin4,fe=[( sin(x(i)),i=1,ni)])
   call save_file(file_name='fv_d5sin-ord'//trim(strz(o,2))//'.dat',header='x sin5 cos' ,x=x(1:ni),fn=sin5,fe=[( cos(x(i)),i=1,ni)])
enddo
contains
   subroutine save_file(file_name, header, x, fn, fe)
   !< Save results on file.
   character(*), intent(in) :: file_name !< File name.
   character(*), intent(in) :: header    !< Header of file.
   real(R8P),    intent(in) ::  x(1:)    !< Abscissa.
   real(R8P),    intent(in) :: fn(1:)    !< Numerical function.
   real(R8P),    intent(in) :: fe(1:)    !< Exact function.
   integer(I4P)             :: j, fu     !< Counter.

   open(newunit=fu, file=trim(adjustl(file_name)))
   write(fu, '(A)') trim(adjustl(header))
   do j=1, size(x)
      write(fu, '(A)') trim(str(n=[x(j),fn(j),fe(j)],separator=' '))
   enddo
   close(fu)
   endsubroutine save_file
endprogram test_fdv_operators_trigonometric
