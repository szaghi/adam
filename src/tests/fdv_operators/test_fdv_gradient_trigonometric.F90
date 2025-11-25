program test_fdv_gradient_trigonometric
!< Test ADAM fdv operators (namely gradient) library with trigonometric functions derivatives.

!< Test function, scalar field
!<```
!< q(x,y,z) = sin(x)*sin(y)*sin(z), (x,y,z) in [0,2PI]x[0,2PI]x[0,2PI].
!<```
!< The exact gradient is
!<```
!< grad q(x,y,z) = [q,x, q,y, q,z]  = [cos(x)*sin(y)*sin(z), cos(y)*sin(x)*sin(z), cos(z)*sin(x)*sin(y)]
!<```

use adam_fdv_operators_library
use penf

implicit none

integer(I4P), parameter         :: n=100_I4P                                   !< Number of cells in each direction.
integer(I4P), parameter         :: gc=5_I4P                                    !< Number of ghost cells.
real(R8P),    parameter         :: PI=4._R8P*atan(1._R8P)                      !< PI.
real(R8P),    parameter         :: h=2_I4P*PI/n                                !< Space step.
real(R8P), dimension(1-gc:n+gc) :: x, y, z                                     !< Cells centers coordinates.
real(R8P)                       :: nodes(1:3,0:n,0:n,0:n)                      !< Nodes coordinates.
real(R8P)                       :: q(           1-gc:n+gc,1-gc:n+gc,1-gc:n+gc) !< Function sin(x)*sin(y)*sin(z).
real(R8P)                       :: gradient(1:3,1   :n   ,1   :n   ,1   :n   ) !< Gradient (vectorial) function.
real(R8P)                       :: exact(   1:3,1   :n   ,1   :n   ,1   :n   ) !< Exact gradient function.
real(R8P)                       :: error(   1:3,1   :n   ,1   :n   ,1   :n   ) !< Error function.
real(R8P)                       :: errorL2                                     !< L2 error norm.
integer(I4P)                    :: s                                           !< Half stencil length.
integer(I4P)                    :: i, j, k, d, o                               !< Counter.

! initialize
do k=0, n ; do j=0, n ; nodes(1,0:n,j,k) = [(i*h,i=0,n)] ; enddo ; enddo
do k=0, n ; do i=0, n ; nodes(2,i,0:n,k) = [(j*h,j=0,n)] ; enddo ; enddo
do j=0, n ; do i=0, n ; nodes(3,i,j,0:n) = [(k*h,k=0,n)] ; enddo ; enddo
do k=1-gc, n+gc
   z(k) = k*h
   do j=1-gc, n+gc
      y(j) = j*h
      do i=1-gc, n+gc
         x(i) = i*h
         q(i,j,k) = sin(x(i))*sin(y(j))*sin(z(k))
         if (i>=1.and.j>=1.and.k>=1.and.i<=n.and.j<=n.and.k<=n) &
            exact(1:3,i,j,k) = [cos(x(i))*sin(y(j))*sin(z(k)), cos(y(j))*sin(x(i))*sin(z(k)), cos(z(k))*sin(x(i))*sin(y(j))]
      enddo
   enddo
enddo

! test FD gradient
print '(A)', 'FD schemes'
do o=2, 10, 2
   s = o/2
   do k=1, n
   do j=1, n
   do i=1, n
      call compute_gradient_fd_centered(s=s,dxyz=[h,h,h],q=q(i-s:i+s,j-s:j+s,k-s:k+s),gradient=gradient(1:3,i,j,k))
   enddo
   enddo
   enddo
   error = abs(gradient-exact)
   errorL2 = 0._R8P
   do k=1, n
   do j=1, n
   do i=1, n
      do d=1,3
         errorL2 = errorL2 + error(d,i,j,k)*error(d,i,j,k)
      enddo
   enddo
   enddo
   enddo
   errorL2 = sqrt(errorL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
   call save_file(basename='fd_gradient-ord', order=o)
enddo

! test FV gradient
print '(A)', 'FV schemes'
do o=2, 10, 2
   s = o/2
   do k=1, n
   do j=1, n
   do i=1, n
      call compute_gradient_fv_centered(s=s,dxyz=[h,h,h],q=q(i-s:i+s,j-s:j+s,k-s:k+s),gradient=gradient(1:3,i,j,k))
   enddo
   enddo
   enddo
   error = abs(gradient-exact)
   errorL2 = 0._R8P
   do k=1, n
   do j=1, n
   do i=1, n
      do d=1,3
         errorL2 = errorL2 + error(d,i,j,k)*error(d,i,j,k)
      enddo
   enddo
   enddo
   enddo
   errorL2 = sqrt(errorL2*h*h)
   print '(A)', '  Order '//trim(strz(o,2))//': Error L0 '//trim(str(maxval(error)))//', Error L2 '//trim(str(errorL2))
   call save_file(basename='fv_gradient-ord', order=o)
enddo

contains
   subroutine save_file(basename, order)
   !< Save results on file.
   character(*), intent(in) :: basename   !< base file name.
   integer(I4P), intent(in) :: order      !< Order of accuracy.
   integer(I4P)             :: i,j,k,d,fu !< Counter.

   open(newunit=fu, file=trim(adjustl(basename))//trim(strz(order,2))//'.dat')
   write(fu,'(A)') 'TITLE = "Gradient results, order '//trim(strz(order))//'"'
   write(fu,'(A)') 'VARIABLES = "X", "Y", "Z", '//                                          &
                   '"Numerical_Gradient_x" "Numerical_Gradient_y" "Numerical_Gradient_z" '//&
                   '"Exact_Gradient_x" "Exact_Gradient_y" "Exact_Gradient_z" '//            &
                   '"Error_x" "Error_y" "Error_z"'
   write(fu,'(A,I0,A,I0,A,I0,A)')'ZONE I=',n+1,', J=',n+1,', K=',n+1,&
                                 ', DATAPACKING=BLOCK, VARLOCATION=([4,5,6,7,8,9,10,11,12]=CELLCENTERED)'
   do d=1, 3
      do k=0, n
      do j=0, n
      do i=0, n
         write(fu,'(A)') trim(str(nodes(d,i,j,k)))
      enddo
      enddo
      enddo
   enddo
   do d=1, 3
      do k=1, n
      do j=1, n
      do i=1, n
         write(fu,'(A)') trim(str(gradient(d,i,j,k)))
      enddo
      enddo
      enddo
   enddo
   do d=1, 3
      do k=1, n
      do j=1, n
      do i=1, n
         write(fu,'(A)') trim(str(exact(d,i,j,k)))
      enddo
      enddo
      enddo
   enddo
   do d=1, 3
      do k=1, n
      do j=1, n
      do i=1, n
         write(fu,'(A)') trim(str(error(d,i,j,k)))
      enddo
      enddo
      enddo
   enddo
   close(fu)
   endsubroutine save_file
endprogram test_fdv_gradient_trigonometric
