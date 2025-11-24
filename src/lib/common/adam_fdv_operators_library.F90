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
public :: compute_derivative3_fdv_interface
public :: compute_derivative4_fdv_interface
public :: compute_derivative5_fdv_interface
public :: compute_derivative6_fdv_interface
public :: compute_divergence_fdv_interface
public :: compute_gradient_fdv_interface
public :: compute_laplacian_fdv_interface
! finite difference
public :: compute_curl_fd_centered
public :: compute_derivative1_fd_centered
public :: compute_derivative2_fd_centered
public :: compute_derivative3_fd_centered
public :: compute_derivative4_fd_centered
public :: compute_derivative5_fd_centered
public :: compute_derivative6_fd_centered
public :: compute_divergence_fd_centered
public :: compute_gradient_fd_centered
public :: compute_laplacian_fd_centered
! finite volume
public :: compute_curl_fv_centered
public :: compute_derivative1_fv_centered
public :: compute_derivative2_fv_centered
public :: compute_derivative3_fv_centered
public :: compute_derivative4_fv_centered
public :: compute_derivative5_fv_centered
public :: compute_derivative6_fv_centered
public :: compute_divergence_fv_centered
public :: compute_gradient_fv_centered
public :: compute_laplacian_fv_centered
public :: compute_reconstruction_r_fv_centered

integer(I4P), parameter :: S_MAX=5_I4P !< Maximum (half) stencil length.

!< Finite Difference method (pointwise values at cell centers) centered schemes.
!< Approximate \(\frac{d^n q}{ds^n}\) at \(x_i\) as:
!< \[
!< \frac{d^n q}{ds^n}\bigg|_i \approx \frac{1}{Ds^n} \sum_{m=-M}^{M} c_m^{(p)} q_{i+m}
!< \]
!< where \(p\) is the order of accuracy (e.g. 2, 4, 6, 8, 10).
!< Derivative of order 1
!< | Order \(p\) | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                    |
!< |-------------|--------------------------|---------------------------------------------------------------|
!< | 2nd  (S=1)  |            -1,0,1        | 1/2   *(                   -1  ,  0, 1                      ) |
!< | 4th  (S=2)  |         -2,-1,0,1,2      | 1/12  *(              1  , -8  ,  0, 8   , -1               ) |
!< | 6th  (S=3)  |      -3,-2,-1,0,1,2,3    | 1/60  *(        -1  , 9  , -45 ,  0, 45  , -9  ,  1         ) |
!< | 8th  (S=4)  |   -4,-3,-2,-1,0,1,2,3,4  | 1/840 *(    3 , -32 , 168, -672,  0, 672 , -168,  32,  -3   ) |
!< | 10th (S=5)  |-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/2520*(-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2) |
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
!< Derivative of order 2
!< | Order \(p\) | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                       |
!< |-------------|--------------------------|------------------------------------------------------------------|
!< | 2nd  (S=1)  |            -1,0,1        | 1      *(                     1 ,    -2, 1                     ) |
!< | 4th  (S=2)  |         -2,-1,0,1,2      | 1/12   *(              -1 ,   16,   -30, 16  , -1              ) |
!< | 6th  (S=3)  |      -3,-2,-1,0,1,2,3    | 1/180  *(          2,  -27,  270,  -490, 270 , -27 ,  2        ) |
!< | 8th  (S=4)  |   -4,-3,-2,-1,0,1,2,3,4  | 1/5040 *(    -9, 128,-1008, 8064,-14350, 8064,-1008, 128,  -9  ) |
!< | 10th (S=5)  |-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/25200*(8,-125,1000,-6000,42000,-73766,42000,-6000,1000,-125,8) |
!< Coefficient are symmetric respect i, parametrize only half of the stencil coefficients.
                                           !0           1          2          3         4         5
real(R8P), parameter :: FD2_CC_S1(S_MAX+1)=[    -2._R8P,    1._R8P,    0._R8P,   0._R8P,   0._R8P,0._R8P]           !< FD2C, S1.
real(R8P), parameter :: FD2_CC_S2(S_MAX+1)=[   -30._R8P,   16._R8P,   -1._R8P,   0._R8P,   0._R8P,0._R8P]/12._R8P   !< FD2C, S2.
real(R8P), parameter :: FD2_CC_S3(S_MAX+1)=[  -490._R8P,  270._R8P,  -27._R8P,   2._R8P,   0._R8P,0._R8P]/180._R8P  !< FD2C, S3.
real(R8P), parameter :: FD2_CC_S4(S_MAX+1)=[-14350._R8P, 8064._R8P,-1008._R8P, 128._R8P,  -9._R8P,0._R8P]/5040._R8P !< FD2C, S4.
real(R8P), parameter :: FD2_CC_S5(S_MAX+1)=[-73766._R8P,42000._R8P,-6000._R8P,1000._R8P,-125._R8P,8._R8P]/25200._R8P!< FD2C, S5.
real(R8P), parameter :: FD2_CC(S_MAX+1,S_MAX)=reshape([FD2_CC_S1, &
                                                       FD2_CC_S2, &
                                                       FD2_CC_S3, &
                                                       FD2_CC_S4, &
                                                       FD2_CC_S5],&
                                                      [S_MAX+1,S_MAX]) !< Finite difference derivative 2 centered coefficients.
!< Derivative of order 3
!< |Order p   | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                                 |
!< |----------|--------------------------|----------------------------------------------------------------------------|
!< |     (S=1)|            -1,0,1        | N.A. (too small stencil)                                                   |
!< |2th  (S=2)|         -2,-1,0,1,2      | 1/2    *(                     -1,    +2,+0,    -2,    +1                  )|
!< |4th  (S=3)|      -3,-2,-1,0,1,2,3    | 1/8    *(              1 ,    -8,   +13,+0,   -13,    +8,    -1           )|
!< |6th  (S=4)|   -4,-3,-2,-1,0,1,2,3,4  | 1/240  *(       -7,   +72,  -338,  +488,+0,  -488,  +338,   -72,   +7     )|
!< |8th  (S=5)|-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/30240*(205,-2522,+14607,-52428,+70098,+0,-70098,+52428,-14607,+2522,-205)|
!< Coefficient are antisymmetric respect i, parametrize only half of the stencil coefficients.
                                         !1           2           3           4          5
real(R8P), parameter :: FD3_CC_S1(S_MAX)=[     0._R8P,     0._R8P,     0._R8P,    0._R8P,   0._R8P]           !< FD3C,S1.
real(R8P), parameter :: FD3_CC_S2(S_MAX)=[    -2._R8P,    +1._R8P,     0._R8P,    0._R8P,   0._R8P]/2._R8P    !< FD3C,S2.
real(R8P), parameter :: FD3_CC_S3(S_MAX)=[   -13._R8P,    +8._R8P,    -1._R8P,    0._R8P,   0._R8P]/8._R8P    !< FD3C,S3.
real(R8P), parameter :: FD3_CC_S4(S_MAX)=[  -488._R8P,  +338._R8P,   -72._R8P,   +7._R8P,   0._R8P]/240._R8P  !< FD3C,S4.
real(R8P), parameter :: FD3_CC_S5(S_MAX)=[-70098._R8P,+52428._R8P,-14607._R8P,+2522._R8P,-205._R8P]/30240._R8P!< FD3C,S5.
real(R8P), parameter :: FD3_CC(S_MAX,S_MAX)=reshape([FD3_CC_S1, &
                                                     FD3_CC_S2, &
                                                     FD3_CC_S3, &
                                                     FD3_CC_S4, &
                                                     FD3_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference derivative 3 centered coefficients.
!< Derivative of order 4
!< |Order p   | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                                 |
!< |----------|--------------------------|----------------------------------------------------------------------------|
!< |     (S=1)|            -1,0,1        | N.A. (too small stencil)                                                   |
!< |2th  (S=2)|         -2,-1,0,1,2      | 1/1    *(                   1,     -4,     6,     -4,    1                 |
!< |4th  (S=3)|      -3,-2,-1,0,1,2,3    | 1/6    *(            -1,   12,    -39,    56,    -39,   12,   -1         ) |
!< |6th  (S=4)|   -4,-3,-2,-1,0,1,2,3,4  | 1/240  *(       7,  -96,  676,  -1952,  2730,  -1952,  676,  -96,   7    ) |
!< |8th  (S=5)|-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/15120*(-82,1261,-9738,52428,-140196,192654,-140196,52428,-9738,1261,-82) |
!< Coefficient are symmetric respect i, parametrize only half of the stencil coefficients.
                                           !0           1            2          3          4         5
real(R8P), parameter :: FD4_CC_S1(S_MAX+1)=[     0._R8P,     0._R8P,     0._R8P,    0._R8P,   0._R8P,  0._R8P]           !< FD4C,S1.
real(R8P), parameter :: FD4_CC_S2(S_MAX+1)=[     6._R8P,     -4._R8P,    1._R8P,    0._R8P,   0._R8P,  0._R8P]           !< FD4C,S2.
real(R8P), parameter :: FD4_CC_S3(S_MAX+1)=[    56._R8P,    -39._R8P,   12._R8P,   -1._R8P,   0._R8P,  0._R8P]/6._R8P    !< FD4C,S3.
real(R8P), parameter :: FD4_CC_S4(S_MAX+1)=[  2730._R8P,  -1952._R8P,  676._R8P,  -96._R8P,   7._R8P,  0._R8P]/240._R8P  !< FD4C,S4.
real(R8P), parameter :: FD4_CC_S5(S_MAX+1)=[192654._R8P,-140196._R8P,52428._R8P,-9738._R8P,1261._R8P,-82._R8P]/15120._R8P!< FD4C,S5.
real(R8P), parameter :: FD4_CC(S_MAX+1,S_MAX)=reshape([FD4_CC_S1, &
                                                       FD4_CC_S2, &
                                                       FD4_CC_S3, &
                                                       FD4_CC_S4, &
                                                       FD4_CC_S5],&
                                                      [S_MAX+1,S_MAX]) !< Finite difference derivative 4 centered coefficients.
!< Derivative of order 5
!< |Order p   | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                     |
!< |----------|--------------------------|----------------------------------------------------------------|
!< |     (S=1)|            -1,0,1        | N.A. (too small stencil)                                       |
!< |     (S=2)|         -2,-1,0,1,2      | N.A. (too small stencil)                                       |
!< |2th  (S=3)|      -3,-2,-1,0,1,2,3    | 1/2  *(           -1,   +4,   -5,+0,   +5,   -4,  +1         ) |
!< |4th  (S=4)|   -4,-3,-2,-1,0,1,2,3,4  | 1/6  *(       1,  -9,  +26,  -29,+0,  +29,  -26,  +9,  -1    ) |
!< |6th  (S=5)|-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/288*(-13,+152,-783,+1872,-1938,+0,+1938,-1872,+783,-152,+13) |
!< Coefficient are antisymmetric respect i, parametrize only half of the stencil coefficients.
                                         !1          2          3         4         5
real(R8P), parameter :: FD5_CC_S2(S_MAX)=[    0._R8P,    0._R8P,   0._R8P,   0._R8P,  0._R8P]         !< FD5C,S1.
real(R8P), parameter :: FD5_CC_S1(S_MAX)=[    0._R8P,    0._R8P,   0._R8P,   0._R8P,  0._R8P]         !< FD5C,S2.
real(R8P), parameter :: FD5_CC_S3(S_MAX)=[   +5._R8P,   -4._R8P,  +1._R8P,   0._R8P,  0._R8P]/2._R8P  !< FD5C,S3.
real(R8P), parameter :: FD5_CC_S4(S_MAX)=[  +29._R8P,  -26._R8P,  +9._R8P,  -1._R8P,  0._R8P]/6._R8P  !< FD5C,S4.
real(R8P), parameter :: FD5_CC_S5(S_MAX)=[+1938._R8P,-1872._R8P,+783._R8P,-152._R8P,+13._R8P]/288._R8P!< FD5C,S5.
real(R8P), parameter :: FD5_CC(S_MAX,S_MAX)=reshape([FD5_CC_S1, &
                                                     FD5_CC_S2, &
                                                     FD5_CC_S3, &
                                                     FD5_CC_S4, &
                                                     FD5_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference derivative 5 centered coefficients.
!< Derivative of order 6
!< |Order p   | Stencil points \(m\)     | Coefficients \(c_m^{(p)}\)                                     |
!< |----------|--------------------------|----------------------------------------------------------------|
!< |     (S=1)|            -1,0,1        | N.A. (too small stencil)                                       |
!< |     (S=2)|         -2,-1,0,1,2      | N.A. (too small stencil)                                       |
!< |2th  (S=3)|      -3,-2,-1,0,1,2,3    | 1/1  *(           1,   -6,  15,   -20,  15,   -6,   1        ) |
!< |4th  (S=4)|   -4,-3,-2,-1,0,1,2,3,4  | 1/4  *(     -1,  12,  -52, 116,  -150, 116,  -52,  12,  -1   ) |
!< |6th  (S=5)|-5,-4,-3,-2,-1,0,1,2,3,4,5| 1/240*(13,-190,1305,-4680,9690,-12276,9690,-4680,1305,-190,13) |
!< Coefficient are symmetric respect i, parametrize only half of the stencil coefficients.
                                           !0           1         2          3         4         5
real(R8P), parameter :: FD6_CC_S2(S_MAX+1)=[    0._R8P,    0._R8P,    0._R8P,   0._R8P,   0._R8P, 0._R8P]         !< FD6C,S1.
real(R8P), parameter :: FD6_CC_S1(S_MAX+1)=[    0._R8P,    0._R8P,    0._R8P,   0._R8P,   0._R8P, 0._R8P]         !< FD6C,S2.
real(R8P), parameter :: FD6_CC_S3(S_MAX+1)=[   -20._R8P,  15._R8P,   -6._R8P,   1._R8P,   0._R8P, 0._R8P]         !< FD6C,S3.
real(R8P), parameter :: FD6_CC_S4(S_MAX+1)=[  -150._R8P, 116._R8P,  -52._R8P,  12._R8P,  -1._R8P, 0._R8P]/4._R8P  !< FD6C,S4.
real(R8P), parameter :: FD6_CC_S5(S_MAX+1)=[-12276._R8P,9690._R8P,-4680._R8P,1305._R8P,-190._R8P,13._R8P]/240._R8P!< FD6C,S5.
real(R8P), parameter :: FD6_CC(S_MAX+1,S_MAX)=reshape([FD6_CC_S1, &
                                                       FD6_CC_S2, &
                                                       FD6_CC_S3, &
                                                       FD6_CC_S4, &
                                                       FD6_CC_S5],&
                                                      [S_MAX+1,S_MAX]) !< Finite difference derivative 6 centered coefficients.

!< Finite Volume (volumetric averages, derivative from flux differences) centered schemes.
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
!< Derivative of order 1
!< | Order \(p\)| Stencil points \(m\)  | Coefficients \(a_m^{(p)}\) for \(q_{i+1/2}\)     |
!< |------------|-----------------------|--------------------------------------------------|
!< | 2nd  (S=1) |            0,1        | 1/2   *(               1   ,1                  ) |
!< | 4th  (S=2) |         -1,0,1,2      | 1/12  *(          -1  ,7   ,7   ,-1            ) |
!< | 6th  (S=3) |      -2,-1,0,1,2,3    | 1/60  *(      1  ,-8  ,37  ,37  ,-8  ,1        ) |
!< | 8th  (S=4) |   -3,-2,-1,0,1,2,3,4  | 1/840 *(  -3 ,29 ,-139,533 ,533 ,-139,29 ,-3   ) |
!< | 10th (S=5) |-4,-3,-2,-1,0,1,2,3,4,5| 1/2520*(3,-30,180,-840,2107,2107,-840,180,-30,3) |
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

   pure subroutine compute_derivative3_fdv_interface(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   endsubroutine compute_derivative3_fdv_interface

   pure subroutine compute_derivative4_fdv_interface(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   endsubroutine compute_derivative4_fdv_interface

   pure subroutine compute_derivative5_fdv_interface(s,ds,q,d5q_ds5)
   !< Compute derivative of order 2 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   endsubroutine compute_derivative5_fdv_interface

   pure subroutine compute_derivative6_fdv_interface(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 of scalar field.
   import :: R8P, I4P
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q, d6q/ds6.
   endsubroutine compute_derivative6_fdv_interface

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
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: dq_ds   !< Derivative of order 1 of q, dq/ds.
   integer(I4P)              :: m       !< Counter.

   dq_ds = 0.0_R8P
   do m=1, s
      dq_ds = dq_ds + FD1_CC(m,s)*(q(1+m) - q(1-m))
   enddo
   dq_ds = dq_ds/ds
   endsubroutine compute_derivative1_fd_centered

   pure subroutine compute_derivative2_fd_centered(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   integer(I4P)              :: m       !< Counter.

   d2q_ds2 = FD2_CC(1,s)*q(1)
   do m=1, s
      d2q_ds2 = d2q_ds2 + FD2_CC(m+1,s)*(q(1+m) + q(1-m))
   enddo
   d2q_ds2 = d2q_ds2/ds/ds
   endsubroutine compute_derivative2_fd_centered

   pure subroutine compute_derivative3_fd_centered(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   integer(I4P)              :: m       !< Counter.

   d3q_ds3 = 0.0_R8P
   do m=1, s
      d3q_ds3 = d3q_ds3 + FD3_CC(m,s)*(q(1+m) - q(1-m))
   enddo
   d3q_ds3 = d3q_ds3/ds/ds/ds
   endsubroutine compute_derivative3_fd_centered

   pure subroutine compute_derivative4_fd_centered(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s-2.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   integer(I4P)              :: m       !< Counter.

   d4q_ds4 = FD4_CC(1,s)*q(1)
   do m=1, s
      d4q_ds4 = d4q_ds4 + FD4_CC(m+1,s)*(q(1+m) + q(1-m))
   enddo
   d4q_ds4 = d4q_ds4/ds/ds/ds/ds
   endsubroutine compute_derivative4_fd_centered

   pure subroutine compute_derivative5_fd_centered(s,ds,q,d5q_ds5)
   !< Compute derivative of order 5 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   integer(I4P)              :: m       !< Counter.

   d5q_ds5 = 0.0_R8P
   do m=1, s
      d5q_ds5 = d5q_ds5 + FD5_CC(m,s)*(q(1+m) - q(1-m))
   enddo
   d5q_ds5 = d5q_ds5/ds/ds/ds/ds/ds
   endsubroutine compute_derivative5_fd_centered

   pure subroutine compute_derivative6_fd_centered(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s-4.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q, d6q/ds6.
   integer(I4P)              :: m       !< Counter.

   d6q_ds6 = FD6_CC(1,s)*q(1)
   do m=1, s
      d6q_ds6 = d6q_ds6 + FD6_CC(m+1,s)*(q(1+m) + q(1-m))
   enddo
   d6q_ds6 = d6q_ds6/ds/ds/ds/ds/ds/ds
   endsubroutine compute_derivative6_fd_centered

   pure subroutine compute_divergence_fd_centered(s,dxyz,q,divergence)
   !< Compute divergence of q vector field with finite difference centered scheme.
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
   integer(I4P), intent(in)  :: s                       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)                !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:)       !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian               !< Laplacian of q.
   real(R8P)                 :: d2q_dx2,d2q_dy2,d2q_dz2 !< Laplacian parts.

   call compute_derivative2_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),d2q_ds2=d2q_dx2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),d2q_ds2=d2q_dy2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),d2q_ds2=d2q_dz2)
   laplacian = d2q_dx2 + d2q_dy2 + d2q_dz2
   endsubroutine compute_laplacian_fd_centered

   ! finite volume schemes
   pure subroutine compute_curl_fv_centered(s,dxyz,q,curl)
   !< Compute curl of q vector field with finite volume centered scheme.
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
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: dq_ds   !< Derivative of order 1 of q, dq/ds.
   real(R8P)                 :: ql,qr   !< Reconstruction of field at left and righ interfaces.

   call compute_reconstruction_r_fv_centered(s=s,q=q(1-s:  s),qr=ql)
   call compute_reconstruction_r_fv_centered(s=s,q=q(2-s:1+s),qr=qr)
   dq_ds = (qr-ql)/ds
   endsubroutine compute_derivative1_fv_centered

   pure subroutine compute_derivative2_fv_centered(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite volume centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   real(R8P)                 :: dql,dqr !< Derivative 1 at left and right cells.

   call compute_derivative1_fv_centered(s=s-1,ds=ds,q=q(1-s  :1+s-2),dq_ds=dql)
   call compute_derivative1_fv_centered(s=s-1,ds=ds,q=q(1-s+2:1+s  ),dq_ds=dqr)
   d2q_ds2 = (dqr - dql)/(2.0_R8P * ds)
   endsubroutine compute_derivative2_fv_centered

   pure subroutine compute_derivative3_fv_centered(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 with finite volume centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   real(R8P)                 :: dql,dqr !< Derivative 2 at left and right cells.

   call compute_derivative2_fv_centered(s=s-1,ds=ds,q=q(1-s  :1+s-2),d2q_ds2=dql)
   call compute_derivative2_fv_centered(s=s-1,ds=ds,q=q(1-s+2:1+s  ),d2q_ds2=dqr)
   d3q_ds3 = (dqr - dql)/(2.0_R8P * ds)
   endsubroutine compute_derivative3_fv_centered

   pure subroutine compute_derivative4_fv_centered(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 with finite volume centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   real(R8P)                 :: dql,dqr !< Derivative 3 at left and right cells.

   call compute_derivative3_fv_centered(s=s-1,ds=ds,q=q(1-s  :1+s-2),d3q_ds3=dql)
   call compute_derivative3_fv_centered(s=s-1,ds=ds,q=q(1-s+2:1+s  ),d3q_ds3=dqr)
   d4q_ds4 = (dqr - dql)/(2.0_R8P * ds)
   endsubroutine compute_derivative4_fv_centered

   pure subroutine compute_derivative5_fv_centered(s,ds,q,d5q_ds5)
   !< Compute derivative of order 5 with finite volume centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   real(R8P)                 :: dql,dqr !< Derivative 4 at left and right cells.

   call compute_derivative4_fv_centered(s=s-1,ds=ds,q=q(1-s  :1+s-2),d4q_ds4=dql)
   call compute_derivative4_fv_centered(s=s-1,ds=ds,q=q(1-s+2:1+s  ),d4q_ds4=dqr)
   d5q_ds5 = (dqr - dql)/(2.0_R8P * ds)
   endsubroutine compute_derivative5_fv_centered

   pure subroutine compute_derivative6_fv_centered(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 with finite volume centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q, d6q/ds6.
   real(R8P)                 :: dql,dqr !< Derivative 5 at left and right cells.

   call compute_derivative5_fv_centered(s=s-1,ds=ds,q=q(1-s  :1+s-2),d5q_ds5=dql)
   call compute_derivative5_fv_centered(s=s-1,ds=ds,q=q(1-s+2:1+s  ),d5q_ds5=dqr)
   d6q_ds6 = (dqr - dql)/(2.0_R8P * ds)
   endsubroutine compute_derivative6_fv_centered

   pure subroutine compute_divergence_fv_centered(s,dxyz,q,divergence)
   !< Compute divergence of q vector field with finite volume centered scheme.
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
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=gradient(1))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=gradient(2))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=gradient(3))
   endsubroutine compute_gradient_fv_centered

   pure subroutine compute_laplacian_fv_centered(s,dxyz,q,laplacian)
   !< Compute laplacian of q scalar field with finite volume centered scheme.
   integer(I4P), intent(in)  :: s                       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)                !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:)       !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian               !< Laplacian of q.
   real(R8P)                 :: d2q_dx2,d2q_dy2,d2q_dz2 !< Laplacian parts.

   call compute_derivative2_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),d2q_ds2=d2q_dx2)
   call compute_derivative2_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),d2q_ds2=d2q_dy2)
   call compute_derivative2_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),d2q_ds2=d2q_dz2)
   laplacian = d2q_dx2 + d2q_dy2 + d2q_dz2
   endsubroutine compute_laplacian_fv_centered

   pure subroutine compute_reconstruction_r_fv_centered(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values. Used for finite volume approach where
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
