!< ADAM, WENO class NVF kernels (NVF backend of [[weno_object]]).
module adam_weno_nvf_kernels
!< ADAM, WENO class NVF kernels (NVF backend of [[weno_object]]).

use penf, only : I4P, R8P
use cudafor

implicit none
private
public :: weno_reconstruct_upwind_kernel

integer(I4P), parameter :: S_max=6    !< Maximum number/dimensions of stencils.
integer(I4P), parameter :: S_max_m1=5 !< Maximum number/dimensions of stencils minus 1.
contains
   ! public procedures
   attributes(device) subroutine weno_reconstruct_upwind_kernel(S, weno_a, weno_p, weno_d, weno_zeps, V, VR)
   !< Reconstruct by WENO upwind method of 2S-1 order, non TBP.
   integer(I4P), intent(in)          :: S                   !< Number of stencils used.
   real(R8P),    intent(in), device  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in), device  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in), device  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)          :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)          :: V (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out)         :: VR(1:2)             !< Left and right (1,2) interface value of reconstructed V.
   real(R8P)                         :: VP(1:2,0:S_max_m1)  !< Polynomial reconstructions.
   real(R8P)                         :: w (1:2,0:S_max_m1)  !< Weights of the stencils.

   call weno_compute_polynomials_kernel(S=S, weno_p=weno_p, V=V(1:2,1-S:-1+S), VP=VP(1:2,0:S-1))
   call weno_compute_weights_kernel(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, V=V(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution_kernel(S=S, VP=VP(1:2,0:S-1), w=w(1:2,0:S-1), VR=VR(1:2))
   endsubroutine weno_reconstruct_upwind_kernel

   ! private procedures
   attributes(device) subroutine weno_compute_convolution_kernel(S, VP, w, VR)
   !< Compute WENO convulution, non TBP.
   integer(I4P), intent(in)  :: S             !< Number of stencils used.
   real(R8P),    intent(in)  :: VP(1:2,0:S-1) !< Polynomial reconstructions.
   real(R8P),    intent(in)  :: w (1:2,0:S-1) !< Weights of the stencils.
   real(R8P),    intent(out) :: VR(1:2      ) !< Left and right (1,2) interface value of reconstructed V.
   integer(I4P)              :: k,f           !< Counter.

   VR = 0._R8P
   do k=0, S-1
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         VR(f) = VR(f) + w(f,k)*VP(f,k)
      enddo
   enddo
   endsubroutine weno_compute_convolution_kernel

   attributes(device) subroutine weno_compute_polynomials_kernel(S, weno_p, V, VP)
   !< Compute WENO polynomials, non TBP.
   integer(I4P), intent(in)         :: S                   !< Number of stencils used.
   real(R8P),    intent(in), device :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)         :: V (1:2,1-S:-1+S)    !< Variable to be reconstructed.
   real(R8P),    intent(out)        :: VP(1:2,0:S-1   )    !< Polynomial reconstructions.
   integer(I4P)                     :: s1,s2,f             !< Counter.

   ! computing the polynomials
   VP = 0._R8P
   do s1=0, S-1 ! stencil counter
      do s2=0, S-1 ! cell counter counter
         do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            VP(f,s1) = VP(f,s1) + weno_p(f,s2,s1,S)*V(f,-s2+s1)
         enddo
      enddo
   enddo
   endsubroutine weno_compute_polynomials_kernel

   attributes(device) subroutine weno_compute_weights_kernel(S, weno_a, weno_d, weno_zeps, V, w)
   !< Compute WENO weights, non TBP.
   integer(I4P), intent(in)         :: S                     !< Number of stencils used.
   real(R8P),    intent(in), device :: weno_a(1:,0:,1:)      !< Optimal weights.
   real(R8P),    intent(in), device :: weno_d(0:,0:,0:,1:)   !< Smoothness indicators coefficients.
   real(R8P),    intent(in)         :: weno_zeps             !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)         :: V    (1:2,1-S:-1+S)   !< Variable to be reconstructed.
   real(R8P),    intent(out)        :: w    (1:2,0:S-1)      !< Weights of the stencils.
   real(R8P)                        :: IS   (1:2,0:S_max_m1) !< Smoothness indicators of the stencils.
   real(R8P)                        :: a    (1:2,0:S_max_m1) !< Alpha coifficients for the weights.
   real(R8P)                        :: a_tot(1:2)            !< Summ of the alpha coefficients.
   integer(I4P)                     :: s1,s2,s3,f            !< Counter.

   ! computing smoothness indicators
   do s1=0,S-1 ! stencil counter
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         IS(f,s1) = 0._R8P
         do s2=0,S-1
            do s3=0,S-1
              IS(f,s1) = IS(f,s1) + weno_d(s3,s2,s1,S)*V(f,s1-s3)*V(f,s1-s2)
            enddo
         enddo
      enddo
   enddo
   ! computing alfa coefficients
   a_tot = 0._R8P
   do s1=0,S-1
      do f=1,2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         a(f,s1) = weno_a(f,s1,S)*(1._R8P/(weno_zeps+IS(f,s1))**S) ; a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   ! computing weights
   do s1=0,S-1
      do f=1,2
         w(f,s1) = a(f,s1)/a_tot(f)
      enddo
   enddo
   endsubroutine weno_compute_weights_kernel
endmodule adam_weno_nvf_kernels
