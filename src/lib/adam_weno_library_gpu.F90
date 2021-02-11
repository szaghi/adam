!< ADAM, WENO reconstructor libray, GPU backend.
module adam_weno_library_gpu
!< ADAM, WENO reconstructor libray, GPU backend.

use PENF
use CUDAFOR
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
private
public :: weno_initialize, weno_l_side, weno_r_side
public :: reconstruct_weno

integer(I4P)                      :: weno_s=1_I4P  !< WENO stencils number/dimension.
real(R8P),    allocatable, device :: weno_a(:,:)   !< Optimal linear weights [1:2,0:s-1].
real(R8P),    allocatable, device :: weno_p(:,:,:) !< Polynomial coefficents [1:2,0:s-1,0:s-1].
real(R8P),    allocatable, device :: weno_d(:,:,:) !< IS         coefficents [0:s-1,0:s-1,0:s-1].
integer(I4P), allocatable, device :: weno_l_side   !< Left side flag.
integer(I4P), allocatable, device :: weno_r_side   !< Left side flag.

contains
   subroutine weno_initialize(weno_stencils)
   !< Initialize WENO library.
   integer(I4P), intent(in), device :: weno_stencils !< WENO stencils number/dimension.

   weno_s = weno_stencils
   if (allocated(weno_a)) deallocate(weno_a)
   if (allocated(weno_p)) deallocate(weno_p)
   if (allocated(weno_d)) deallocate(weno_d)
   weno_l_side = 1_I4P
   weno_r_side = 2_I4P
   if (weno_s>1) then
      allocate(weno_a(1:2,0:weno_s-1                  ))
      allocate(weno_p(1:2,0:weno_s-1,0:weno_s-1       ))
      allocate(weno_d(0:weno_s-1,0:weno_s-1,0:weno_s-1))
   endif
   select case(weno_s)
   case(2)
   case(3)
      ! left reconstruction
      weno_a(1,0) = 0.1_R8P ! stencil 0
      weno_a(1,1) = 0.6_R8P ! stencil 1
      weno_a(1,2) = 0.3_R8P ! stencil 2
      !  cell  0                     ;    cell  1                     ;    cell  2
      weno_p(1,0,0) = 11._R8P/6._R8P ; weno_p(1,1,0) = -7._R8P/6._R8P ; weno_p(1,2,0) =  1._R8P/3._R8P ! stencil 0
      weno_p(1,0,1) =  1._R8P/3._R8P ; weno_p(1,1,1) =  5._R8P/6._R8P ; weno_p(1,2,1) = -1._R8P/6._R8P ! stencil 1
      weno_p(1,0,2) = -1._R8P/6._R8P ; weno_p(1,1,2) =  5._R8P/6._R8P ; weno_p(1,2,2) =  1._R8P/3._R8P ! stencil 2
      ! right reconstruction
      weno_a(2,0) = 0.3_R8P ! stencil 0
      weno_a(2,1) = 0.6_R8P ! stencil 1
      weno_a(2,2) = 0.1_R8P ! stencil 2
      !  cell  0                   ;    cell  1                   ;    cell  2
      weno_p(2,0,0) =  1._R8P/3._R8P ; weno_p(2,1,0) =  5._R8P/6._R8P ; weno_p(2,2,0) = -1._R8P/6._R8P ! stencil 0
      weno_p(2,0,1) = -1._R8P/6._R8P ; weno_p(2,1,1) =  5._R8P/6._R8P ; weno_p(2,2,1) =  1._R8P/3._R8P ! stencil 1
      weno_p(2,0,2) =  1._R8P/3._R8P ; weno_p(2,1,2) = -7._R8P/6._R8P ; weno_p(2,2,2) = 11._R8P/6._R8P ! stencil 2

      ! stencil 0
      !      i*i                      ;       (i-1)*i                   ;       (i-2)*i
      weno_d(0,0,0) =  10._R8P/3._R8P ; weno_d(1,0,0) = -31._R8P/3._R8P ; weno_d(2,0,0) =  11._R8P/3._R8P
      !      /                        ;       (i-1)*(i-1)               ;       (i-2)*(i-1)
      weno_d(0,1,0) =   0._R8P        ; weno_d(1,1,0) =  25._R8P/3._R8P ; weno_d(2,1,0) = -19._R8P/3._R8P
      !      /                        ;        /                        ;       (i-2)*(i-2)
      weno_d(0,2,0) =   0._R8P        ; weno_d(1,2,0) =   0._R8P        ; weno_d(2,2,0) =   4._R8P/3._R8P
      ! stencil 1
      !     (i+1)*(i+1)               ;        i*(i+1)                  ;       (i-1)*(i+1)
      weno_d(0,0,1) =   4._R8P/3._R8P ; weno_d(1,0,1) = -13._R8P/3._R8P ; weno_d(2,0,1) =   5._R8P/3._R8P
      !      /                        ;        i*i                      ;       (i-1)*i
      weno_d(0,1,1) =   0._R8P        ; weno_d(1,1,1) =  13._R8P/3._R8P ; weno_d(2,1,1) = -13._R8P/3._R8P
      !      /                        ;        /                        ;       (i-1)*(i-1)
      weno_d(0,2,1) =   0._R8P        ; weno_d(1,2,1) =   0._R8P        ; weno_d(2,2,1) =   4._R8P/3._R8P
      ! stencil 2
      !     (i+2)*(i+2)               ;       (i+1)*(i+2)               ;        i*(i+2)
      weno_d(0,0,2) =   4._R8P/3._R8P ; weno_d(1,0,2) = -19._R8P/3._R8P ; weno_d(2,0,2) =  11._R8P/3._R8P
      !      /                        ;       (i+1)*(i+1)               ;        i*(i+1)
      weno_d(0,1,2) =   0._R8P        ; weno_d(1,1,2) =  25._R8P/3._R8P ; weno_d(2,1,2) = -31._R8P/3._R8P
      !      /                        ;        /                        ;        i*i
      weno_d(0,2,2) =   0._R8P        ; weno_d(1,2,2) =   0._R8P        ; weno_d(2,2,2) =  10._R8P/3._R8P
   endselect
   endsubroutine weno_initialize

   attributes(device) subroutine reconstruct_weno(side, s, q, qr)
   !< Reconstruct variable with WENO 5 scheme.
   integer(I4P), intent(in)  :: side        !< Side of reconstruction, 1 left, 2 right.
   integer(I4P), intent(in)  :: s           !< Stencils number/dimension.
   real(R8P),    intent(in)  :: q(1-s:-1+s) !< Stencil values.
   real(R8P),    intent(out) :: qr          !< Reconstructed variable.
   real(R8P)                 :: qp(0:s-1)   !< Polynomial reconstructions.
   real(R8P)                 :: IS(0:s-1)   !< Smoothness indicators of the stencils.
   real(R8P)                 :: a(0:s-1)    !< Alpha coefficients for the weights.
   real(R8P)                 :: a_tot       !< Summ of the alpha coefficients.
   real(R8P)                 :: w(0:s-1)    !< Weights of the stencils.
   integer(I4P)              :: s1, s2, s3  !< Counters

   do s1=0,s-1 ! stencil counter
      IS(s1) = 0._R8P
      do s2=0,s-1
         do s3=0,s-1
            IS(s1) = IS(s1) + weno_d(s3,s2,s1) * q(s1-s3) * q(s1-s2)
         enddo
      enddo
   enddo
   a_tot = 0._R8P
   do s1=0,s-1
      a(s1) = weno_a(side,s1) * (1._R8P/(1d-6 + IS(s1))**3)
      a_tot = a_tot + a(s1)
   enddo
   do s1=0,s-1
      w(s1) = a(s1) / a_tot
   enddo
   qp = 0._R8P
   do s1=0,s-1 ! stencil counter
      do s2=0,s-1 ! cell counter counter
         qp(s1) = qp(s1) + weno_p(side,s2,s1) * q(-s2+s1)
      enddo
   enddo
   qr = 0._R8P
   do s1=0,s-1
      qr = qr + w(s1) * qp(s1)
   enddo
   endsubroutine reconstruct_weno
endmodule adam_weno_library_gpu
