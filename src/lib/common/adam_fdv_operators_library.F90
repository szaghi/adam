!< ADAM, finite difference/volume operators approximations library.
module adam_fdv_operators_library
!< ADAM, finite difference/volume operators approximations library.

use penf

implicit none
save
private
! interfaces
public :: compute_curl_fdv_interface
public :: compute_derivative1_fdv_interface
public :: compute_derivative2_fdv_interface
public :: compute_divergence_fdv_interface
public :: compute_gradient_fdv_interface
public :: compute_laplacian_fdv_interface
! finite difference
public :: compute_curl_fd_centered
public :: compute_derivative1_fd_centered
public :: compute_derivative2_fd_centered
public :: compute_divergence_fd_centered
public :: compute_gradient_fd_centered
public :: compute_laplacian_fd_centered
! finite volume
public :: compute_curl_fv_centered
public :: compute_derivative1_fv_centered
public :: compute_divergence_fv_centered
public :: compute_gradient_fv_centered
public :: compute_reconstruction_r_fv_centered

integer(I4P), parameter :: S_MAX=5_I4P !< Maximum stencil length.

!< Derivative of order 1 Finite Difference (pointwise values at cell centers) centered schemes
!< Approximate \(\frac{dq}{ds}\) at \(x_i\) as:
!< \[
!< \frac{dq}{ds}\bigg|_i \approx \frac{1}{Ds} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
!< \]
!< where \(p\) is the order of accuracy (2, 4, 6, 8, 10), \(M = \frac{p}{2}\), and
!< coefficients \(c_m^{(p)}\) are symmetric.
!< Finite Difference Coefficients \(c_m^{(p)}\):
!< | Order \(p\) | Stencil points \(m\)                 | Coefficients \(c_m^{(p)}\)                                    |
!< |-------------|--------------------------------------|---------------------------------------------------------------|
!< | 2nd  (S=1)  |                 -1, 0, 1             | 1/2   *(                   -1  ,  0, 1                      ) |
!< | 4th  (S=2)  |             -2, -1, 0, 1, 2          | 1/12  *(              1  , -8  ,  0, 8   , -1               ) |
!< | 6th  (S=3)  |         -3, -2, -1, 0, 1, 2, 3       | 1/60  *(        -1  , 9  , -45 ,  0, 45  , -9  ,  1         ) |
!< | 8th  (S=4)  |     -4, -3, -2, -1, 0, 1, 2, 3, 4    | 1/840 *(    3 , -32 , 168, -672,  0, 672 , -168,  32,  -3   ) |
!< | 10th (S=5)  | -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 | 1/2520*(-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2) |
!< Coefficient are antisymmetric respect i, parametrize only half of the stencil coefficients.
                                         !1         2         3        4        5
real(R8P), parameter :: FD1_CC_S1(S_MAX)=[   1._R8P,   0._R8P,  0._R8P,  0._R8P,0._R8P]/2._R8P    !< FD1C, S1.
real(R8P), parameter :: FD1_CC_S2(S_MAX)=[   8._R8P,  -1._R8P,  0._R8P,  0._R8P,0._R8P]/12._R8P   !< FD1C, S2.
real(R8P), parameter :: FD1_CC_S3(S_MAX)=[  45._R8P,  -9._R8P,  1._R8P,  0._R8P,0._R8P]/60._R8P   !< FD1C, S3.
real(R8P), parameter :: FD1_CC_S4(S_MAX)=[ 672._R8P,-168._R8P, 32._R8P, -3._R8P,0._R8P]/840._R8P  !< FD1C, S4.
real(R8P), parameter :: FD1_CC_S5(S_MAX)=[2100._R8P,-600._R8P,150._R8P,-25._R8P,2._R8P]/2520._R8P !< FD1C, S5.
real(R8P), parameter :: FD1_CC(S_MAX,S_MAX)=reshape([FD1_CC_S1, &
                                                     FD1_CC_S2, &
                                                     FD1_CC_S3, &
                                                     FD1_CC_S4, &
                                                     FD1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference derivative 1 centered coefficients.

!< Derivative of order 2 Finite Difference (pointwise values at cell centers) centered schemes
!< Approximate \(\frac{dq}{ds}\) at \(x_i\) as:
!< \[
!< \frac{d^2q}{ds^2}\bigg|_i \approx \frac{1}{Ds} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
!< \]
!< where \(p\) is the order of accuracy (2, 4, 6, 8, 10), \(M = \frac{p}{2}\), and
!< coefficients \(c_m^{(p)}\) are symmetric.
!< Finite Difference Coefficients \(c_m^{(p)}\):
!< | Order \(p\) | Stencil points \(m\)                 | Coefficients \(c_m^{(p)}\)                                       |
!< |-------------|--------------------------------------|------------------------------------------------------------------|
!< | 2nd  (S=1)  |                 -1, 0, 1             | 1      *(                     1 ,    -2, 1                     ) |
!< | 4th  (S=2)  |             -2, -1, 0, 1, 2          | 1/12   *(              -1 ,   16,   -30, 16  , -1              ) |
!< | 6th  (S=3)  |         -3, -2, -1, 0, 1, 2, 3       | 1/180  *(          2,  -27,  270,  -490, 270 , -27 ,  2        ) |
!< | 8th  (S=4)  |     -4, -3, -2, -1, 0, 1, 2, 3, 4    | 1/5040 *(    -9, 128,-1008, 8064,-14350, 8064,-1008, 128,  -9  ) |
!< | 10th (S=5)  | -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5 | 1/25200*(8,-125,1000,-6000,42000,-73766,42000,-6000,1000,-125,8) |
!< Coefficient are antisymmetric respect i, parametrize only half of the stencil coefficients.
                                           !0           1          2          3         4         5
real(R8P), parameter :: FD2_CC_S1(S_MAX+1)=[    -2._R8P,    1._R8P,    0._R8P,   0._R8P,   0._R8P,0._R8P]           !< FD1C, S1.
real(R8P), parameter :: FD2_CC_S2(S_MAX+1)=[   -30._R8P,   16._R8P,   -1._R8P,   0._R8P,   0._R8P,0._R8P]/12._R8P   !< FD2C, S2.
real(R8P), parameter :: FD2_CC_S3(S_MAX+1)=[  -490._R8P,  270._R8P,  -27._R8P,   2._R8P,   0._R8P,0._R8P]/180._R8P  !< FD2C, S3.
real(R8P), parameter :: FD2_CC_S4(S_MAX+1)=[-14350._R8P, 8064._R8P,-1008._R8P, 128._R8P,  -9._R8P,0._R8P]/5040._R8P !< FD2C, S4.
real(R8P), parameter :: FD2_CC_S5(S_MAX+1)=[-73766._R8P,42000._R8P,-6000._R8P,1000._R8P,-125._R8P,8._R8P]/25200._R8P!< FD2C, S5.
real(R8P), parameter :: FD2_CC(S_MAX+1,S_MAX)=reshape([FD2_CC_S1, &
                                                       FD2_CC_S2, &
                                                       FD2_CC_S3, &
                                                       FD2_CC_S4, &
                                                       FD2_CC_S5],&
                                                      [S_MAX+1,S_MAX]) !< Finite difference derivative 1 centered coefficients.

!< Derivative of order 1 Finite Difference Finite Volume (volumetric averages, derivative from flux differences) centered schemes
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
real(R8P), parameter :: FV1_CC_S1(S_MAX)=[   1._R8P,   0._R8P,  0._R8P,  0._R8P,0._R8P]/2._R8P    !< FV1C, S1.
real(R8P), parameter :: FV1_CC_S2(S_MAX)=[   7._R8P,  -1._R8P,  0._R8P,  0._R8P,0._R8P]/12._R8P   !< FV1C, S2.
real(R8P), parameter :: FV1_CC_S3(S_MAX)=[  37._R8P,  -8._R8P,  1._R8P,  0._R8P,0._R8P]/60._R8P   !< FV1C, S3.
real(R8P), parameter :: FV1_CC_S4(S_MAX)=[ 533._R8P,-139._R8P, 29._R8P, -3._R8P,0._R8P]/840._R8P  !< FV1C, S4.
real(R8P), parameter :: FV1_CC_S5(S_MAX)=[2107._R8P,-840._R8P,180._R8P,-30._R8P,3._R8P]/2520._R8P !< FV1C, S5.
real(R8P), parameter :: FV1_CC(S_MAX,S_MAX)=reshape([FV1_CC_S1, &
                                                     FV1_CC_S2, &
                                                     FV1_CC_S3, &
                                                     FV1_CC_S4, &
                                                     FV1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite volume derivative 1 centered coefficients.

interface
   pure subroutine compute_curl_fdv_interface(s,dxyz,q,curl)
   !< Compute curl of q vector field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: curl(1:)             !< Gradient of q [1:3].
   endsubroutine compute_curl_fdv_interface

   pure subroutine compute_derivative1_fdv_interface(s,ds,q,dq_ds)
   !< Compute derivative of order 1 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: dq_ds   !< Derivative of order 1 of q, dq/ds.
   endsubroutine compute_derivative1_fdv_interface

   pure subroutine compute_derivative2_fdv_interface(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   endsubroutine compute_derivative2_fdv_interface

   pure subroutine compute_divergence_fdv_interface(s,dxyz,q,divergence)
   !< Compute divergence of q vector field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: divergence           !< Divergence of q.
   endsubroutine compute_divergence_fdv_interface

   pure subroutine compute_gradient_fdv_interface(s,dxyz,q,gradient)
   !< Compute gradient of q scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].
   endsubroutine compute_gradient_fdv_interface

   pure subroutine compute_laplacian_fdv_interface(s,dxyz,q,laplacian)
   !< Compute laplacian of q scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian         !< Lapliacian of q.
   endsubroutine compute_laplacian_fdv_interface
endinterface

contains
   ! public methods
   ! finite difference schemes
   pure subroutine compute_curl_fd_centered(s,dxyz,q,curl)
   !< Compute curl of q vector field with finite difference centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: curl(1:)             !< Curl of q [1:3].
   real(R8P)                 :: dqx_dy, dqx_dz       !< Derivatives of qx.
   real(R8P)                 :: dqy_dx, dqy_dz       !< Derivatives of qy.
   real(R8P)                 :: dqz_dx, dqz_dy       !< Derivatives of qz.

   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(1,1      ,1-s:1+s,1      ),dq_ds=dqx_dy)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(1,1      ,1      ,1-s:1+s),dq_ds=dqx_dz)

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(2,1-s:1+s,1      ,1      ),dq_ds=dqy_dx)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(2,1      ,1      ,1-s:1+s),dq_ds=dqy_dz)

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(3,1-s:1+s,1      ,1      ),dq_ds=dqz_dx)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(3,1      ,1-s:1+s,1      ),dq_ds=dqz_dy)

   curl(1) = dqz_dy - dqy_dz
   curl(2) = dqx_dz - dqz_dx
   curl(3) = dqy_dx - dqx_dy
   endsubroutine compute_curl_fd_centered

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
      dq_ds = dq_ds + FD1_CC(m,s)*(q(1+m) - q(1-m))/ds
   enddo
   endsubroutine compute_derivative1_fd_centered

   pure subroutine compute_derivative2_fd_centered(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite difference centered scheme.
   !< \[
   !< \frac{d^2q}{ds^2}\bigg|_i \approx \frac{1}{Ds} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
   !< \]
   !< The vector field q must be passed with a stencil large enough to computed the derivative with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   integer(I4P)              :: m       !< Counter.

   d2q_ds2 = FD2_CC(1,s)*q(0)/ds/ds
   do m=1, s
      d2q_ds2 = d2q_ds2 + FD2_CC(m+1,s)*(q(1+m) + q(1-m))/ds/ds
   enddo
   endsubroutine compute_derivative2_fd_centered

   pure subroutine compute_divergence_fd_centered(s,dxyz,q,divergence)
   !< Compute divergence of q vector field with finite difference centered scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: divergence           !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z  !< Divergence components.

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(1,1-s:1+s,1,1),dq_ds=div_x)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(2,1,1-s:1+s,1),dq_ds=div_y)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(3,1,1,1-s:1+s),dq_ds=div_z)
   divergence = div_x + div_y + div_z
   endsubroutine compute_divergence_fd_centered

   pure subroutine compute_gradient_fd_centered(s,dxyz,q,gradient)
   !< Compute gradient of q scalar field with finite difference centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=gradient(1))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=gradient(2))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=gradient(3))
   endsubroutine compute_gradient_fd_centered

   pure subroutine compute_laplacian_fd_centered(s,dxyz,q,laplacian)
   !< Compute laplacian of q scalar field with finite difference centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)                !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:)       !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian               !< Lapliacian of q.
   real(R8P)                 :: d2q_dx2,d2q_dy2,d2q_dz2 !< Lapliacian parts.

   call compute_derivative2_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),d2q_ds2=d2q_dx2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),d2q_ds2=d2q_dy2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),d2q_ds2=d2q_dz2)
   laplacian = d2q_dx2 + d2q_dy2 + d2q_dz2
   endsubroutine compute_laplacian_fd_centered

   ! finite volume schemes
   pure subroutine compute_curl_fv_centered(s,dxyz,q,curl)
   !< Compute curl of q vector field with finite volume centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: curl(1:)             !< Curl of q [1:3].
   real(R8P)                 :: dqx_dy, dqx_dz       !< Derivatives of qx.
   real(R8P)                 :: dqy_dx, dqy_dz       !< Derivatives of qy.
   real(R8P)                 :: dqz_dx, dqz_dy       !< Derivatives of qz.

   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1      ,1-s:1+s,1      ),dq_ds=dqx_dy)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1      ,1      ,1-s:1+s),dq_ds=dqx_dz)

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(2,1-s:1+s,1      ,1      ),dq_ds=dqy_dx)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(2,1      ,1      ,1-s:1+s),dq_ds=dqy_dz)

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(3,1-s:1+s,1      ,1      ),dq_ds=dqz_dx)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(3,1      ,1-s:1+s,1      ),dq_ds=dqz_dy)

   curl(1) = dqz_dy - dqy_dz
   curl(2) = dqx_dz - dqz_dx
   curl(3) = dqy_dx - dqx_dy
   endsubroutine compute_curl_fv_centered

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

   pure subroutine compute_divergence_fv_centered(s,dxyz,q,divergence)
   !< Compute divergence of q vector field with finite volume centered scheme.
   !< The vector field q must be passed with a stencil large enough to computed the divergence with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(1:3,i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1:,1-s:,1-s:,1-s:) !< Vector field over the stencil [1:3,1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: divergence           !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z  !< Divergence components.

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1,1-s:1+s,1,1),dq_ds=div_x)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(2,1,1-s:1+s,1),dq_ds=div_y)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(3,1,1,1-s:1+s),dq_ds=div_z)
   divergence = div_x + div_y + div_z
   endsubroutine compute_divergence_fv_centered

   pure subroutine compute_gradient_fv_centered(s,dxyz,q,gradient)
   !< Compute gradient of q scalar field with finite volume centered scheme.
   !< The scalar field q must be passed with a stencil large enough to computed the gradient with selected order of accuracy and
   !< the stencil must be centered in i,j,k, i.e. q = q(i-order/2:i+order/2,j-order/2:j+order/2,k-order/2:k+order/2).
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=gradient(1))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=gradient(2))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=gradient(3))
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
      qr = qr + FV1_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered
endmodule adam_fdv_operators_library
