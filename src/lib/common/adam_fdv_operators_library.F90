!< ADAM, finite difference/volume operators approximations library.
module adam_fdv_operators_library
!< ADAM, finite difference/volume operators approximations library.

use penf

implicit none
save
private
! interfaces
public :: compute_divergence_fdv_interface
public :: compute_gradient_fdv_interface
! finite difference
public :: compute_derivative1_fd_centered
public :: compute_divergence_fd_centered
public :: compute_gradient_fd_centered
! finite volume
public :: compute_derivative1_fv_centered
public :: compute_divergence_fv_centered
public :: compute_gradient_fv_centered
public :: compute_reconstruction_r_fv_centered

integer(I4P), parameter :: S_MAX=5_I4P !< Maximum stencil length.
!< Finite Difference (pointwise values at cell centers)
!< Approximate \(\frac{dq}{ds}\) at \(x_i\) as:
!< \[
!< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
!< \]
!< where \(p\) is the order of accuracy (2, 4, 6, 8, 10), \(M = \frac{p}{2}\), and
!< coefficients \(c_m^{(p)}\) are symmetric.
!< Finite Difference Coefficients \(c_m^{(p)}\):
!< | Order \(p\) | Stencil points \(m\)                 | Coefficients \(c_m^{(p)}\)                                     |
!< |-------------|--------------------------------------|----------------------------------------------------------------|
!< | 2nd  (S=1)  |                 -1, 0, 1             | 1/2   *(                    -1,    0, 1                      ) |
!< | 4th  (S=2)  |             -2, -1, 0, 1, 2          | 1/12  *(              1  ,  -8,    0, 8  ,  -1               ) |
!< | 6th  (S=3)  |         -3, -2, -1, 0, 1, 2, 3       | 1/60  *(        -1 ,  9  ,  -45,   0, 45 ,  -9  , 1          ) |
!< | 8th  (S=4)  |     -4, -3, -2, -1, 0, 1, 2, 3, 4    | 1/840 *(    3 , -32,  168,  -672,  0, 672,  -168, 32,  -3    ) |
!< | 10th (S=5)  | -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 | 1/2520*(-3, 33, -210, 1020, -2947, 0, 2947, -1020,210, -33, 3) |
!< Coefficient are antisymmetric respect i, parametrize only half of the stencil coefficients.
                                           !1         2         3        4         5
real(R8P), parameter :: FD_CC_S1(S_MAX)   =[   1._R8P,    0._R8P,  0._R8P,  0._R8P,0._R8P]/2._R8P    !< FD c-coef, S1.
real(R8P), parameter :: FD_CC_S2(S_MAX)   =[   8._R8P,   -1._R8P,  0._R8P,  0._R8P,0._R8P]/12._R8P   !< FD c-coef, S2.
real(R8P), parameter :: FD_CC_S3(S_MAX)   =[  45._R8P,   -9._R8P,  1._R8P,  0._R8P,0._R8P]/60._R8P   !< FD c-coef, S3.
real(R8P), parameter :: FD_CC_S4(S_MAX)   =[ 672._R8P, -168._R8P, 32._R8P, -3._R8P,0._R8P]/840._R8P  !< FD c-coef, S4.
real(R8P), parameter :: FD_CC_S5(S_MAX)   =[2947._R8P,-1020._R8P,210._R8P,-33._R8P,3._R8P]/2520._R8P !< FD c-coef, S5.
real(R8P), parameter :: FD_CC(S_MAX,S_MAX)=reshape([FD_CC_S1, &
                                                    FD_CC_S2, &
                                                    FD_CC_S3, &
                                                    FD_CC_S4, &
                                                    FD_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference centered coefficients.

!< Finite Volume (volumetric averages, derivative from flux differences)
!< Approximate face values \(q_{i+1/2}\) and \(q_{i-1/2}\) as:
!< \[
!< q_{i+1/2} \approx \sum_{m=-N}^{N} a_m^{(p)} q_{i+m}
!< \]
!< \[
!< q_{i-1/2} \approx \sum_{m=-N}^{N} a_m^{(p)} q_{i-1+m}
!< \]
!< with \(N = \frac{p}{2} - 1\).
!< Then,
!< \[
!< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \left( q_{i+1/2} - q_{i-1/2} \right) =
!< \frac{1}{Ds} \sum_{m=-M}^{M} b_m^{(p)} q_{i+m}
!< \]
!< where \(b_m^{(p)} = a_m^{(p)} - a_{m-1}^{(p)}\).
!< | Order \(p\) | Stencil points \(m\)             | Coefficients \(a_m^{(p)}\) for \(q_{i+1/2}\) reconstruction |
!< |-------------|----------------------------------|-------------------------------------------------------------|
!< | 2nd  (S=1)  |                 0, 1             | 1/2   *(                   1   , 1                      )   |
!< | 4th  (S=2)  |             -1, 0, 1, 2          | 1/12  *(             -1  , 7   , 7   , -1               )   |
!< | 6th  (S=3)  |         -2, -1, 0, 1, 2, 3       | 1/60  *(        1  , -8  , 37  , 37  , -8  , 1          )   |
!< | 8th  (S=4)  |     -3, -2, -1, 0, 1, 2, 3, 4    | 1/840 *(   -3 , 29 , -139, 533 , 533 , -139, 29 , -3    )   |
!< | 10th (S=5)  | -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 | 1/2520*(3, -30, 180, -840, 2107, 2107, -840, 180, -30, 3)   |
!< Coefficient are symmetric respect i+1/2, parametrize only half of the stencil coefficients.
                                           !1         2         3        4        5
real(R8P), parameter :: FV_CC_S1(S_MAX)   =[   1._R8P,   0._R8P,  0._R8P,  0._R8P,0._R8P]/2._R8P    !< FV c-coef, S1.
real(R8P), parameter :: FV_CC_S2(S_MAX)   =[   7._R8P,  -1._R8P,  0._R8P,  0._R8P,0._R8P]/12._R8P   !< FV c-coef, S2.
real(R8P), parameter :: FV_CC_S3(S_MAX)   =[  37._R8P,  -8._R8P,  1._R8P,  0._R8P,0._R8P]/60._R8P   !< FV c-coef, S3.
real(R8P), parameter :: FV_CC_S4(S_MAX)   =[ 533._R8P,-139._R8P, 29._R8P, -3._R8P,0._R8P]/840._R8P  !< FV c-coef, S4.
real(R8P), parameter :: FV_CC_S5(S_MAX)   =[2107._R8P,-840._R8P,180._R8P,-30._R8P,3._R8P]/2520._R8P !< FV c-coef, S5.
real(R8P), parameter :: FV_CC(S_MAX,S_MAX)=reshape([FV_CC_S1, &
                                                    FV_CC_S2, &
                                                    FV_CC_S3, &
                                                    FV_CC_S4, &
                                                    FV_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite volume centered coefficients.

interface
   subroutine compute_divergence_fdv_interface(s,dxyz,q,div)
   !< Compute divergence of q vector field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)            !< Space steps.
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: div                  !< Divergence of q.
   endsubroutine compute_divergence_fdv_interface

   subroutine compute_gradient_fdv_interface(s,dxyz,q,grad)
   !< Compute gradient of q scalar field with finite difference centered scheme.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)         !< Space steps.
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: grad(1:3)         !< Gradient of q.
   endsubroutine compute_gradient_fdv_interface
endinterface

contains
   ! public methods
   ! finite difference schemes
   pure subroutine compute_derivative1_fd_centered(s,ds,q,dq_ds)
   !< Compute derivative of order 1 with finite difference centered scheme.
   !< \[
   !< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
   !< \]
   !< The vector field q must be passed with a stencil large enough to computed the derivative with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: dq_ds   !< Derivative of order 1 of q, dq/ds.
   integer(I4P)              :: m       !< Counter.

   dq_ds = 0.0_R8P
   do m=1, s
      dq_ds = dq_ds + FD_CC(s,m)*(q(1+m) - q(1-m))/ds
   enddo
   endsubroutine compute_derivative1_fd_centered

   pure subroutine compute_divergence_fd_centered(s,dxyz,q,div)
   !< Compute divergence of q vector field with finite difference centered scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)            !< Space steps.
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: div                  !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z  !< Divergence components.

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(1,1-s:1+s,1,1),dq_ds=div_x)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(2,1,1-s:1+s,1),dq_ds=div_y)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(3,1,1,1-s:1+s),dq_ds=div_z)
   div = div_x + div_y + div_z
   endsubroutine compute_divergence_fd_centered

   pure subroutine compute_gradient_fd_centered(s,dxyz,q,grad)
   !< Compute gradient of q scalar field with finite difference centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)         !< Space steps.
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: grad(1:3)         !< Gradient of q.

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=grad(1))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=grad(2))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=grad(3))
   endsubroutine compute_gradient_fd_centered

   ! finite volume schemes
   pure subroutine compute_derivative1_fv_centered(s,ds,q,dq_ds)
   !< Compute derivative of order 1 with finite volume centered scheme.
   !< \[
   !< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \left( q_{i+1/2} - q_{i-1/2} \right) =
   !< \frac{1}{Ds} \sum_{m=-M}^{M} b_m^{(p)} q_{i+m}
   !< \]
   !< The vector field q must be passed with a stencil large enough to computed the derivative with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: dq_ds   !< Derivative of order 1 of q, dq/ds.
   real(R8P)                 :: ql,qr   !< Reconstruction of field at left and righ interfaces.

   call compute_reconstruction_r_fv_centered(s=s,q=q(1-s:  s),qr=ql)
   call compute_reconstruction_r_fv_centered(s=s,q=q(2-s:1+s),qr=qr)
   dq_ds = (qr-ql)/ds
   endsubroutine compute_derivative1_fv_centered

   pure subroutine compute_divergence_fv_centered(s,dxyz,q,div)
   !< Compute divergence of q vector field with finite volume centered scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)            !< Space steps.
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: div                  !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z  !< Divergence components.

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1,1-s:1+s,1,1),dq_ds=div_x)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(2,1,1-s:1+s,1),dq_ds=div_y)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(3,1,1,1-s:1+s),dq_ds=div_z)
   div = div_x + div_y + div_z
   endsubroutine compute_divergence_fv_centered

   pure subroutine compute_gradient_fv_centered(s,dxyz,q,grad)
   !< Compute gradient of q scalar field with finite volume centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:3)         !< Space steps.
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: grad(1:3)         !< Gradient of q.

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=grad(1))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=grad(2))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=grad(3))
   endsubroutine compute_gradient_fv_centered

   pure subroutine compute_reconstruction_r_fv_centered(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values. Used for finite volume approach where
   !< first derivative at cell center can be written as fluxes difference at cell interfaces
   !< \[
   !< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \left( q_{i+1/2} - q_{i-1/2} \right) =
   !< \frac{1}{Ds} \sum_{m=-M}^{M} b_m^{(p)} q_{i+m}
   !< \]
   !< with
   !< \[
   !< q_{i+1/2} \approx \sum_{m=-N}^{N} a_m^{(p)} q_{i+m}
   !< \]
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction at right interface of field.
   integer(I4P)              :: m       !< Counter.

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FV_CC(s,m)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered
endmodule adam_fdv_operators_library
