!< ADAM, vector field operators library.
module adam_vector_field_operators_library
!< ADAM, vector field operators library.

use penf

implicit none
save
private
public :: compute_divergence_centered
public :: compute_divergence_upwind
public :: compute_gradient_centered

real(R8P), parameter :: C_CENTERED(5,5) = reshape([                                  &
   ! 2nd order centered (1-point each side)
   0.5_R8P           , 0.0_R8P          , 0.0_R8P         , 0.0_R8P         ,0.0_R8P,&
   ! 4th order centered (2-points each side)
   -1.0_R8P/12.0_R8P , 2.0_R8P/3.0_R8P  , 0.0_R8P         , 0.0_R8P         ,0.0_R8P,&
   ! 6th order centered (3-points each side)
   1.0_R8P/60.0_R8P  ,-3.0_R8P/20.0_R8P , 3.0_R8P/4.0_R8P , 0.0_R8P         ,0.0_R8P,&
   ! 8th order centered (4-points each side)
   -1.0_R8P/280.0_R8P, 4.0_R8P/105.0_R8P,-1.0_R8P/5.0_R8P , 4.0_R8P/5.0_R8P ,0.0_R8P,&
   ! 10th order centered (5-points each side)
   1.0_R8P/1260.0_R8P,-5.0_R8P/504.0_R8P, 5.0_R8P/84.0_R8P,-5.0_R8P/21.0_R8P,5.0_R8P/6.0_R8P],[5,5]) !< Finite difference coeffs
                                                                                                     !< for centered schemes.
real(R8P), parameter :: C_UPWIND(4,5) = reshape([                                 &
    ! Order 1: 1st order upwind
     1.0_R8P        , 0.0_R8P        , 0.0_R8P        , 0.0_R8P         , 0.0_R8P,&
    ! Order 2: 2nd order upwind
    -1.0_R8P/2.0_R8P, 2.0_R8P        ,-3.0_R8P/2.0_R8P, 0.0_R8P         , 0.0_R8P,&
    ! Order 3: 3rd order upwind
     1.0_R8P/3.0_R8P,-3.0_R8P/2.0_R8P, 3.0_R8P        ,-11.0_R8P/6.0_R8P, 0.0_R8P,&
    ! Order 4: 4th order upwind
    -1.0_R8P/4.0_R8P, 4.0_R8P/3.0_R8P,-3.0_R8P        , 4.0_R8P         ,-25.0_R8P/12.0_R8P],[4,5]) !< Finite difference coeffs
                                                                                                    !< for upwind schemes.
contains
   ! public methods
   subroutine compute_divergence_centered(ngc,dxyz,order,q,div)
   !< Compute divergence of q vector field with finite difference centered scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-ngc:i+ngc,j-ngc:j+ngc,k-ngc:k+ngc).
   integer(I4P), intent(in)  :: ngc                        !< Number of ghost cells.
   real(R8P),    intent(in)  :: dxyz(3)                    !< Space steps.
   integer(I4P), intent(in)  :: order                      !< Accuracy order.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:) !< Vector field over the stencil.
   real(R8P),    intent(out) :: div                        !< Divergence of q.
   integer(I4P)              :: m, s                       !< Counter.

   s = order/2
   div = 0.0_R8P
   ! dq(x)/dx
   do m=1, s
      div = div + C_CENTERED(s,m)*(q(1,1+m,1,1) - q(1,1-m,1,1))/dxyz(1)
   enddo
   ! dq(y)/dy
   do m=1, s
      div = div + C_CENTERED(s,m)*(q(2,1,1+m,1) - q(2,1,1-m,1))/dxyz(2)
   enddo
   ! dq(z)/dz
   do m=1, s
      div = div + C_CENTERED(s,m)*(q(3,1,1,1+m) - q(3,1,1,1-m))/dxyz(3)
   enddo
   endsubroutine compute_divergence_centered

   subroutine compute_divergence_upwind(ngc,dxyz,order,q,div,upwind_direction)
   !< Compute divergence of q vector field with finite difference upwind scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-ngc:i+ngc,j-ngc:j+ngc,k-ngc:k+ngc).
   integer(I4P), intent(in)           :: ngc                        !< Number of ghost cells.
   real(R8P),    intent(in)           :: dxyz(3)                    !< Space steps.
   integer(I4P), intent(in)           :: order                      !< Accuracy order.
   real(R8P),    intent(in)           :: q(1:,1-ngc:,1-ngc:,1-ngc:) !< Vector field over the stencil.
   real(R8P),    intent(out)          :: div                        !< Divergence of q.
   integer(I4P), intent(in), optional :: upwind_direction(3)        !< Upwind direction, +1 for forward, -1 for backward.
   integer(I4P)                       :: ud(3)                      !< Upwind direction, local var.
   integer(I4P)                       :: m                          !< Counter.

   ud = -1 ; if (present(upwind_direction)) ud = upwind_direction
   div = 0.0_R8P
   ! dq(x)/dx
   do m=0, order
      div = div + C_UPWIND(order,m+1)*q(1,1+ud(1)*m,1,1)/dxyz(1)
   enddo
   ! dq(y)/dy
   do m=0, order
      div = div + C_UPWIND(order,m+1)*q(2,1,1+ud(2)*m,1)/dxyz(2)
   enddo
   ! dq(z)/dz
   do m=0, order
      div = div + C_UPWIND(order,m+1)*q(3,1,1,1+ud(3)*m)/dxyz(3)
   enddo
   endsubroutine compute_divergence_upwind

   subroutine compute_gradient_centered(ngc,dxyz,order,q,grad)
   !< Compute gradient of q scalar field with finite difference centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-ngc:i+ngc,j-ngc:j+ngc,k-ngc:k+ngc).
   integer(I4P), intent(in)  :: ngc                     !< Number of ghost cells.
   real(R8P),    intent(in)  :: dxyz(3)                 !< Space steps.
   integer(I4P), intent(in)  :: order                   !< Accuracy order.
   real(R8P),    intent(in)  :: q(1-ngc:,1-ngc:,1-ngc:) !< Vector field over the stencil.
   real(R8P),    intent(out) :: grad(3)                 !< Gradient of q.
   integer(I4P)              :: m, s                    !< Counter.

   s = order/2
   grad = 0.0_R8P
   ! dq(ivar)/dx
   do m=1, s
      grad(1) = grad(1) + C_CENTERED(s,m)*(q(1+m,1,1) - q(1-m,1,1))/dxyz(1)
   enddo
   ! dq(ivar)/dy
   do m=1, s
      grad(2) = grad(2) + C_CENTERED(s,m)*(q(1,1+m,1) - q(1,1-m,1))/dxyz(2)
   enddo
   ! dq(ivar)/dz
   do m=1, s
      grad(3) = grad(3) + C_CENTERED(s,m)*(q(1,1,1+m) - q(1,1,1-m))/dxyz(3)
   enddo
   endsubroutine compute_gradient_centered

   ! old prism non TBP, for reference
   !subroutine compute_div(ni, nj, nk, ngc, blocks_number, dxyz, ivar, q, div)
   !!< Compute div(q(ivar). Finite difference central scheme.
   !integer(I4P), intent(in)    :: ni                              !< Grid cells number in I direction.
   !integer(I4P), intent(in)    :: nj                              !< Grid cells number in J direction.
   !integer(I4P), intent(in)    :: nk                              !< Grid cells number in K direction.
   !integer(I4P), intent(in)    :: ngc                             !< Ghost cells number.
   !integer(I4P), intent(in)    :: blocks_number                   !< Number of blocks.
   !real(R8P),    intent(in)    :: dxyz(1:,1:)                     !< Space steps.
   !integer(I4P), intent(in)    :: ivar                            !< Variable (vectorial) of q.
   !real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field variables.
   !real(R8P),    intent(inout) :: div(   1-ngc:,1-ngc:,1-ngc:,1:) !< Divergence of D, B.
   !integer(I4P)                :: i,j,k,b                         !< Counter

   !do b=1, blocks_number
   !do k=1, nk
   !do j=1, nj
   !do i=1, ni
   !   div(i,j,k,b) = 0.5_R8P*((q(ivar  ,i+1,j,k,b) - q(ivar  ,i-1,j,k,b))/dxyz(1,b) + &
   !                           (q(ivar+1,i,j+1,k,b) - q(ivar+1,i,j-1,k,b))/dxyz(2,b) + &
   !                           (q(ivar+2,i,j,k+1,b) - q(ivar+2,i,j,k-1,b))/dxyz(3,b))

   !enddo
   !enddo
   !enddo
   !enddo
   !endsubroutine compute_div
endmodule adam_vector_field_operators_library
