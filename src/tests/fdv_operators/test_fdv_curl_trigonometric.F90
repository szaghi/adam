program test_fdv_curl_trigonometric
!< Test ADAM fdv operators (namely curl) library with trigonometric functions derivatives.

!< Test function, vectorial field
!<```
!< q(1:3,x,y,z) = [sin(y), sin(z), sin(x)], (x,y,z) in [0,2PI]x[0,2PI]x[0,2PI].
!<```
!< The exact curl is
!<```
!< curl q(1:3,x,y,z) = [qz,y - qy,z, qx,z - qz,x, qy,x - qx,y]  = [-cos(z), -cos(x), -cos(y)]
!<```

use adam_fdv_operators_library
use penf

implicit none

integer(I4P), parameter         :: n=20_I4P                                 !< Number of cells in each direction.
integer(I4P), parameter         :: gc=5_I4P                                 !< Number of ghost cells.
real(R8P),    parameter         :: PI=4._R8P*atan(1._R8P)                   !< PI.
real(R8P),    parameter         :: h=2_I4P*PI/n                             !< Space step.
real(R8P), dimension(1-gc:n+gc) :: x, y, z                                  !< Cells centers coordinates.
real(R8P)                       :: nodes(1:3,0:n,0:n,0:n)                   !< Nodes coordinates.
real(R8P)                       :: q(    1:3,1-gc:n+gc,1-gc:n+gc,1-gc:n+gc) !< Function (vectorial) [sin(y), sin(z), sin(x)].
real(R8P)                       :: curl( 1:3,1   :n   ,1   :n   ,1   :n   ) !< Curl (vectorial) function.
real(R8P)                       :: exact(1:3,1   :n   ,1   :n   ,1   :n   ) !< Exact curl function.
real(R8P)                       :: error(1:3,1   :n   ,1   :n   ,1   :n   ) !< Error function.
real(R8P)                       :: errorL2                                  !< L2 error norm.
integer(I4P)                    :: s                                        !< Half stencil length.
integer(I4P)                    :: i, j, k, d, o                            !< Counter.

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
         q(1:3,i,j,k) = [sin(y(j)), sin(z(k)), sin(x(i))]
         if (i>=1.and.j>=1.and.k>=1.and.i<=n.and.j<=n.and.k<=n) exact(1:3,i,j,k) = [-cos(z(k)), -cos(x(i)), -cos(y(j))]
      enddo
   enddo
enddo

! test FD curl
print '(A)', 'FD schemes'
do o=2, 10, 2
   s = o/2
   do k=1, n
   do j=1, n
   do i=1, n
      call compute_curl_fd_centered(s=s,dxyz=[h,h,h],q=q(1:3,i-s:i+s,j-s:j+s,k-s:k+s),curl=curl(1:3,i,j,k))
   enddo
   enddo
   enddo
   error = abs(curl-exact)
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
   call save_file(basename='fd_curl-ord', order=o)
enddo

! test FV curl
print '(A)', 'FV schemes'
do o=2, 10, 2
   s = o/2
   do k=1, n
   do j=1, n
   do i=1, n
      call compute_curl_fv_centered(s=s,dxyz=[h,h,h],q=q(1:3,i-s:i+s,j-s:j+s,k-s:k+s),curl=curl(1:3,i,j,k))
   enddo
   enddo
   enddo
   error = abs(curl-exact)
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
   call save_file(basename='fv_curl-ord', order=o)
enddo

contains
   subroutine save_file(basename, order)
   !< Save results on file.
   character(*), intent(in) :: basename   !< base file name.
   integer(I4P), intent(in) :: order      !< Order of accuracy.
   integer(I4P)             :: i,j,k,d,fu !< Counter.

   open(newunit=fu, file=trim(adjustl(basename))//trim(strz(order,2))//'.dat')
   write(fu,'(A)') 'TITLE = "Curl results, order '//trim(strz(order))//'"'
   write(fu,'(A)') 'VARIABLES = "X", "Y", "Z", '//                              &
                   '"Numerical_Curl_x" "Numerical_Curl_y" "Numerical_Curl_z" '//&
                   '"Exact_Curl_x" "Exact_Curl_y" "Exact_Curl_z" '//            &
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
         write(fu,'(A)') trim(str(curl(d,i,j,k)))
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
endprogram test_fdv_curl_trigonometric
