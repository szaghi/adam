!< ADAM, WENO class FNL kernels (FNL backend of [[weno_fnl_object]]).

module adam_weno_fnl_kernels
!< ADAM, WENO class FNL kernels (FNL backend of [[weno_fnl_object]]).

use adam_weno_object, only : S_max, S_max_m1
use penf, only : I4P, R8P

implicit none
private
public :: weno_reconstruct_upwind_dev
public :: S_max
public :: S_max_m1

contains
   ! public procedures
   subroutine weno_reconstruct_upwind_dev(S, weno_a, weno_p, weno_d, weno_zeps, V, VR)
   !< Reconstruct by WENO upwind method of 2S-1 order, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)    !< Optimal weights.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps           !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variables to be reconstructed.
   real(R8P),    intent(out) :: VR(1:2)             !< Left and right (1,2) interface value of reconstructed V.
   real(R8P)                 :: VP(1:2,0:S_max_m1)  !< Polynomial reconstructions.
   real(R8P)                 :: w (1:2,0:S_max_m1)  !< Weights of the stencils.
   !$acc routine seq

   call weno_compute_polynomials_device(S=S, weno_p=weno_p, V=V(1:2,1-S:-1+S), VP=VP(1:2,0:S-1))
   call weno_compute_weights_device(S=S, weno_a=weno_a, weno_d=weno_d, weno_zeps=weno_zeps, V=V(1:2,1-S:-1+S), w=w(1:2,0:S-1))
   call weno_compute_convolution_device(S=S, VP=VP(1:2,0:S-1), w=w(1:2,0:S-1), VR=VR(1:2))
   endsubroutine weno_reconstruct_upwind_dev

   ! private procedures
   subroutine weno_compute_convolution_device(S, VP, w, VR)
   !< Compute WENO convulution, non TBP.
   integer(I4P), intent(in)  :: S             !< Number of stencils used.
   real(R8P),    intent(in)  :: VP(1:2,0:S-1) !< Polynomial reconstructions.
   real(R8P),    intent(in)  :: w (1:2,0:S-1) !< Weights of the stencils.
   real(R8P),    intent(out) :: VR(1:2      ) !< Left and right (1,2) interface value of reconstructed V.
   integer(I4P)              :: k,f           !< Counter.
   !$acc routine seq

   VR = 0._R8P
   do k=0, S-1
      do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
         VR(f) = VR(f) + w(f,k)*VP(f,k)
      enddo
   enddo
   endsubroutine weno_compute_convolution_device

   subroutine weno_compute_polynomials_device(S, weno_p, V, VP)
   !< Compute WENO polynomials, non TBP.
   integer(I4P), intent(in)  :: S                   !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_p(1:,0:,0:,1:) !< Polinomials coefficients.
   real(R8P),    intent(in)  :: V (1:2,1-S:-1+S)    !< Variable to be reconstructed.
   real(R8P),    intent(out) :: VP(1:2,0:S-1   )    !< Polynomial reconstructions.
   integer(I4P)              :: s1,s2,f             !< Counter.
   !$acc routine seq

   ! computing the polynomials
   VP = 0._R8P
   do s1=0, S-1 ! stencil counter
      do s2=0, S-1 ! cell counter counter
         do f=1, 2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
            VP(f,s1) = VP(f,s1) + weno_p(f,s2,s1,S)*V(f,-s2+s1)
         enddo
      enddo
   enddo
   endsubroutine weno_compute_polynomials_device

   subroutine weno_compute_weights_device(S, weno_a, weno_d, weno_zeps, V, w)
   !< Compute WENO weights, non TBP.
   integer(I4P), intent(in)  :: S                     !< Number of stencils used.
   real(R8P),    intent(in)  :: weno_a(1:,0:,1:)      !< Optimal weights.
   real(R8P),    intent(in)  :: weno_d(0:,0:,0:,1:)   !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: weno_zeps             !< Parameter for avoiding division by zero in computing IS.
   real(R8P),    intent(in)  :: V    (1:2,1-S:-1+S)   !< Variable to be reconstructed.
   real(R8P),    intent(out) :: w    (1:2,0:S-1)      !< Weights of the stencils.
   real(R8P)                 :: IS   (1:2,0:S_max_m1) !< Smoothness indicators of the stencils.
   real(R8P)                 :: a    (1:2,0:S_max_m1) !< Alpha coifficients for the weights.
   real(R8P)                 :: a_tot(1:2)            !< Summ of the alpha coefficients.
   integer(I4P)              :: s1,s2,s3,f            !< Counter.
   !$acc routine seq

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
         a(f,s1) = weno_a(f,s1,S)*(1._R8P/(weno_zeps+IS(f,s1))**(2)) ; a_tot(f) = a_tot(f) + a(f,s1)
      enddo
   enddo
   ! computing weights
   do s1=0,S-1
      do f=1,2
         w(f,s1) = a(f,s1)/a_tot(f)
      enddo
   enddo
   endsubroutine weno_compute_weights_device
endmodule adam_weno_fnl_kernels
