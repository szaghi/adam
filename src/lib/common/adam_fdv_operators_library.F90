!< ADAM, finite difference/volume operators approximations library.
module adam_fdv_operators_library
!< ADAM, finite difference/volume operators approximations library.

! third party modules
use :: penf

implicit none
private
public :: FDV_S_MAX
public :: FD0_CC
public :: FD1_CC ! device-readable pair-form first-derivative coefficients (issue #22 F1-bis: buffer-free diagnostics)
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
public :: compute_derivative1_fd_seam ! issue #29 E3: matched-resolution 2:1-seam d/dx
public :: compute_derivative2_fd_centered
public :: compute_derivative3_fd_centered
public :: compute_derivative4_fd_centered
public :: compute_derivative5_fd_centered
public :: compute_derivative6_fd_centered
public :: compute_divergence_fd_centered
public :: compute_gradient_fd_centered
public :: compute_laplacian_fd_centered
public :: compute_reconstruction_r_fd_centered
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
public :: compute_reconstruction_r_fv_centered_d1
public :: compute_reconstruction_r_fv_centered_d2
public :: compute_reconstruction_r_fv_centered_d3
public :: compute_reconstruction_r_fv_centered_d4
public :: compute_reconstruction_r_fv_centered_d5
public :: compute_avg_to_point
public :: compute_point_to_avg
public :: compute_derivative1_fv_rupwind
public :: compute_derivative2_fv_rupwind
public :: compute_derivative3_fv_rupwind
public :: compute_derivative4_fv_rupwind
public :: compute_derivative5_fv_rupwind
public :: compute_derivative6_fv_rupwind
public :: compute_reconstruction_r_fv_rupwind
public :: compute_reconstruction_r_fv_rupwind_d1
public :: compute_reconstruction_r_fv_rupwind_d2
public :: compute_reconstruction_r_fv_rupwind_d3
public :: compute_reconstruction_r_fv_rupwind_d4
public :: compute_derivative1_fv_lupwind
public :: compute_derivative2_fv_lupwind
public :: compute_derivative3_fv_lupwind
public :: compute_derivative4_fv_lupwind
public :: compute_derivative5_fv_lupwind
public :: compute_derivative6_fv_lupwind
public :: compute_reconstruction_r_fv_lupwind
public :: compute_reconstruction_r_fv_lupwind_d1
public :: compute_reconstruction_r_fv_lupwind_d2
public :: compute_reconstruction_r_fv_lupwind_d3
public :: compute_reconstruction_r_fv_lupwind_d4

integer(I4P), parameter :: S_MAX=5_I4P     !< Maximum (half) stencil length.
integer(I4P), parameter :: FDV_S_MAX=S_MAX !< Public name of S_MAX.

!< Finite Difference method (pointwise values at cell centers) centered schemes.
!< Finite Difference (point values, Lagrange interpolation) centered schemes.
!< Pointwise nodal data: q_i = q(x_i).
!< Approximate face values \(q_{i+1/2}\) and \(q_{i-1/2}\) as:
!< \[
!< q_{i+1/2} \approx \sum_{m=-(S-1)}^{S} a_m^{(p)} q_{i+m}
!< \]
!< \[
!< q_{i-1/2} \approx \sum_{m=-(S-1)}^{S} a_m^{(p)} q_{i-1+m}
!< \]
!< with \(p = 2S\).
!<
!< Centered (about i+1/2) Lagrange schemes.
!< | Order \(p\)| Stencil points \(m\)  | Coefficients \(a_m^{(p)}\) for \(q_{i+1/2}\)           |
!< |------------|-----------------------|--------------------------------------------------------|
!< | 2nd  (S=1) |            0,1        | 1/2     *(                 1   ,1                    ) |
!< | 4th  (S=2) |         -1,0,1,2      | 1/16    *(            -1  ,9   ,9   ,-1              ) |
!< | 6th  (S=3) |      -2,-1,0,1,2,3    | 1/256   *(       3  ,-25  ,150 ,150 ,-25 ,3          ) |
!< | 8th  (S=4) |   -3,-2,-1,0,1,2,3,4  | 1/2048  *(  -5 ,49 ,-245,1225,1225,-245,49 ,-5       ) |
!< | 10th (S=5) |-4,-3,-2,-1,0,1,2,3,4,5| 1/65536 *(35,-405,2268,-8820,39690,39690,-8820,2268,-405,35) |
!< Coefficients are symmetric w.r.t. i+1/2, parametrize only half of the stencil coefficients.

                                         !1          2          3         4         5
real(R8P), parameter :: FD0_CC_S1(S_MAX)=[    1._R8P,   0._R8P,   0._R8P,  0._R8P,  0._R8P]/2._R8P       !< FD0C, S1.
real(R8P), parameter :: FD0_CC_S2(S_MAX)=[    9._R8P,  -1._R8P,   0._R8P,  0._R8P,  0._R8P]/16._R8P      !< FD0C, S2.
real(R8P), parameter :: FD0_CC_S3(S_MAX)=[  150._R8P, -25._R8P,   3._R8P,  0._R8P,  0._R8P]/256._R8P     !< FD0C, S3.
real(R8P), parameter :: FD0_CC_S4(S_MAX)=[ 1225._R8P,-245._R8P,  49._R8P, -5._R8P,  0._R8P]/2048._R8P    !< FD0C, S4.
real(R8P), parameter :: FD0_CC_S5(S_MAX)=[39690._R8P,-8820._R8P,2268._R8P,-405._R8P,35._R8P]/65536._R8P  !< FD0C, S5.
real(R8P), parameter :: FD0_CC(S_MAX,S_MAX)=reshape([FD0_CC_S1, &
                                                     FD0_CC_S2, &
                                                     FD0_CC_S3, &
                                                     FD0_CC_S4, &
                                                     FD0_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference derivative 0 centered coefficients.
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

! Issue #29 E3: matched-resolution SEAM first-derivative operators (s=3, order 6).
! At a 2:1 coarse-fine seam the coarse (Δx) and fine (Δx/2) sides use different-
! resolution FD1_CC stencils, which breaks the div_h(curl_h)=0 telescoping (E2 proved
! the source is this OPERATOR mismatch, not a value mismatch). These two matched
! operators restore telescoping to the scheme's full order 6 across the jump (verified
! symbolically: untracked/seam-symbolic-checks/t1_1_*.py). Applied ONLY to the
! seam-NORMAL derivative at seam-adjacent cells; tangential dirs use standard FD1_CC.
!
! FINE-side seam cell: samples at fine-length offsets [0,-1,-2,-3, +3/2,+7/2,+11/2]
! along the outward normal — fine interior (0..-3) then coarse donor centers. Result
! is d/dx in units of the FINE spacing. Order 6, L1=2.49 (well conditioned).
real(R8P), parameter :: FD1_SEAM_FINE(7) = [ 323._R8P/462._R8P, -77._R8P/65._R8P,     &
                                              3._R8P/10._R8P,   -77._R8P/1989._R8P,    &
                                             11._R8P/45._R8P,    -2._R8P/91._R8P,       &
                                             21._R8P/12155._R8P ] !< FINE-side matched seam d/dx weights.
integer(I4P), parameter :: FD1_SEAM_FINE_NODE2(7) = [0,-2,-4,-6, 3,7,11] !< 2*node (integer, fine half-units).

! COARSE-side seam cell: samples at coarse-length offsets [0,-1,-2,-3,-4, +1/4,+3/4]
! — coarse interior (0..-4) then fine ghost centers across the jump. Result is d/dx
! in units of the COARSE spacing. Order 6, L1=7.60 (extrapolates onto near-seam fine
! ghosts — the higher L1 is the quantified E3 risk; PRISM horizon is the real test).
real(R8P), parameter :: FD1_SEAM_COARSE(7) = [ -13._R8P/4._R8P,   -12._R8P/35._R8P,   &
                                                 1._R8P/11._R8P,    -4._R8P/195._R8P,  &
                                                 3._R8P/1292._R8P, 4096._R8P/1105._R8P,&
                                             -4096._R8P/21945._R8P ] !< COARSE-side matched seam d/dx weights.
integer(I4P), parameter :: FD1_SEAM_COARSE_NODE4(7) = [0,-4,-8,-12,-16, 1,3] !< 4*node (integer, coarse quarter-units).

! Issue #29 E3 (C1-ext): tangential Lagrange interpolation of the coarse samples onto
! the fine seam cell's tangential line. Without this the coarse samples are piecewise-
! constant across a coarse cell (their tangential derivative jumps), collapsing the
! div(curl) telescoping to order 1 (the E1-class tangential trap — verified). With
! order-6+ tangential interp, telescoping recovers to order 6 (c3_physical.py). 9-point
! interp from the coarse tangential grid (spacing 2 fine); the fine cell is at one of 2
! sub-positions (offset -/+1/2 from the coarse center). Order 8, L1=1.414 (well cond.).
real(R8P), parameter :: FD_SEAM_TANG(9,0:1) = reshape([ &
   -7293._R8P/8388608._R8P, 9945._R8P/1048576._R8P, -109395._R8P/2097152._R8P,       &
   255255._R8P/1048576._R8P, 3828825._R8P/4194304._R8P, -153153._R8P/1048576._R8P,   &
   85085._R8P/2097152._R8P, -8415._R8P/1048576._R8P, 6435._R8P/8388608._R8P,         & ! sub=0 (offset -1/2)
   6435._R8P/8388608._R8P, -8415._R8P/1048576._R8P, 85085._R8P/2097152._R8P,         &
   -153153._R8P/1048576._R8P, 3828825._R8P/4194304._R8P, 255255._R8P/1048576._R8P,   &
   -109395._R8P/2097152._R8P, 9945._R8P/1048576._R8P, -7293._R8P/8388608._R8P],      & ! sub=1 (offset +1/2)
   [9,2])
   !< Tangential interp weights (node 1:9 at coarse-y offsets [-8,-6,-4,-2,0,2,4,6,8]
   !< fine units, sub 0:1). Coarse sample at the fine seam line = FD_SEAM_TANG(:,sub) .
   !< [9 coarse-tangential cell values].
real(R8P), parameter :: FD1_CC(S_MAX,S_MAX)=reshape([FD1_CC_S1, &
                                                     FD1_CC_S2, &
                                                     FD1_CC_S3, &
                                                     FD1_CC_S4, &
                                                     FD1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite difference derivative 1 centered coefficients.
!$acc declare copyin(FD1_CC_S1,FD1_CC_S2,FD1_CC_S3,FD1_CC_S4,FD1_CC_S5,FD1_CC)

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
!$acc declare copyin(FD2_CC_S1,FD2_CC_S2,FD2_CC_S3,FD2_CC_S4,FD2_CC_S5,FD2_CC)

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
!$acc declare copyin(FD3_CC_S1,FD3_CC_S2,FD3_CC_S3,FD3_CC_S4,FD3_CC_S5,FD3_CC)

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
real(R8P), parameter :: FD4_CC_S1(S_MAX+1)=[     0._R8P,      0._R8P,    0._R8P,    0._R8P,   0._R8P,  0._R8P]           !< FD4C,S1.
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
!$acc declare copyin(FD4_CC_S1,FD4_CC_S2,FD4_CC_S3,FD4_CC_S4,FD4_CC_S5,FD4_CC)

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
!$acc declare copyin(FD5_CC_S1,FD5_CC_S2,FD5_CC_S3,FD5_CC_S4,FD5_CC_S5,FD5_CC)

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
!$acc declare copyin(FD6_CC_S1,FD6_CC_S2,FD6_CC_S3,FD6_CC_S4,FD6_CC_S5,FD6_CC)

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

!< Centered schemes.
!< | Order \(p\)| Stencil points \(m\)  | Coefficients \(a_m^{(p)}\) for \(q_{i+1/2}\)     |
!< |------------|-----------------------|--------------------------------------------------|
!< | 2nd  (S=1) |            0,1        | 1/2   *(               1   ,1                  ) |
!< | 4th  (S=2) |         -1,0,1,2      | 1/12  *(          -1  ,7   ,7   ,-1            ) |
!< | 6th  (S=3) |      -2,-1,0,1,2,3    | 1/60  *(      1  ,-8  ,37  ,37  ,-8  ,1        ) |
!< | 8th  (S=4) |   -3,-2,-1,0,1,2,3,4  | 1/840 *(  -3 ,29 ,-139,533 ,533 ,-139,29 ,-3   ) |
!< | 10th (S=5) |-4,-3,-2,-1,0,1,2,3,4,5| 1/2520*(2,-23,127,-473,1627,1627,-473,127,-23,2) |
!< Coefficient are symmetric respect i+1/2, parametrize only half of the stencil coefficients.
                                         !1         2         3        4        5
real(R8P), parameter :: FV1_CC_S1(S_MAX)=[   1._R8P,   0._R8P,  0._R8P,  0._R8P,0._R8P]/2._R8P    !< FV1C, S1.
real(R8P), parameter :: FV1_CC_S2(S_MAX)=[   7._R8P,  -1._R8P,  0._R8P,  0._R8P,0._R8P]/12._R8P   !< FV1C, S2.
real(R8P), parameter :: FV1_CC_S3(S_MAX)=[  37._R8P,  -8._R8P,  1._R8P,  0._R8P,0._R8P]/60._R8P   !< FV1C, S3.
real(R8P), parameter :: FV1_CC_S4(S_MAX)=[ 533._R8P,-139._R8P, 29._R8P, -3._R8P,0._R8P]/840._R8P  !< FV1C, S4.
real(R8P), parameter :: FV1_CC_S5(S_MAX)=[1627._R8P,-473._R8P,127._R8P,-23._R8P,2._R8P]/2520._R8P !< FV1C, S5.
real(R8P), parameter :: FV1_CC(S_MAX,S_MAX)=reshape([FV1_CC_S1, &
                                                     FV1_CC_S2, &
                                                     FV1_CC_S3, &
                                                     FV1_CC_S4, &
                                                     FV1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< Finite volume centered reconstruction coefficients.
!$acc declare copyin(FV1_CC_S1,FV1_CC_S2,FV1_CC_S3,FV1_CC_S4,FV1_CC_S5,FV1_CC)

!< FV face-reconstruction of the m-th derivative at the face i+1/2 from cell averages (primitive-function
!< construction). FVRm_CC(:,s) reconstructs d^m q/ds^m |_{i+1/2}. The FV centered derivative ladder is built
!< as a conservative flux difference of adjacent face reconstructions: d^n q/ds^n |_i = (g_{i+1/2} - g_{i-1/2})/ds
!< with g = FVR(n-1) reconstruction. m=0 is FV1_CC above. Even m symmetric, odd m antisymmetric about i+1/2.
!< Generated by scripts/fdv_coefficients.py (see scripts/fdv_coefficients.md); do not hand-edit.
!< FV face-reconstruction of derivative 1 (cell averages -> face i+1/2), antisymmetric.
real(R8P), parameter :: FVR1_CC_S1(S_MAX)=[    1._R8P,    0._R8P,0._R8P  ,   0._R8P,0._R8P]
real(R8P), parameter :: FVR1_CC_S2(S_MAX)=[   15._R8P,   -1._R8P,0._R8P  ,   0._R8P,0._R8P]/12._R8P
real(R8P), parameter :: FVR1_CC_S3(S_MAX)=[  245._R8P,  -25._R8P,2._R8P  ,   0._R8P,0._R8P]/180._R8P
real(R8P), parameter :: FVR1_CC_S4(S_MAX)=[ 7175._R8P, -889._R8P,119._R8P,  -9._R8P,0._R8P]/5040._R8P
real(R8P), parameter :: FVR1_CC_S5(S_MAX)=[36883._R8P,-5117._R8P,883._R8P,-117._R8P,8._R8P]/25200._R8P
real(R8P), parameter :: FVR1_CC(S_MAX,S_MAX)=reshape([FVR1_CC_S1, &
                                                      FVR1_CC_S2, &
                                                      FVR1_CC_S3, &
                                                      FVR1_CC_S4, &
                                                      FVR1_CC_S5],&
                                                     [S_MAX,S_MAX]) !< FVR1_CC table.
!$acc declare copyin(FVR1_CC_S1,FVR1_CC_S2,FVR1_CC_S3,FVR1_CC_S4,FVR1_CC_S5,FVR1_CC)

!< FV face-reconstruction of derivative 2 (cell averages -> face i+1/2), symmetric.
real(R8P), parameter :: FVR2_CC_S1(S_MAX)=[     0._R8P,    0._R8P,     0._R8P,   0._R8P,0._R8P   ]
real(R8P), parameter :: FVR2_CC_S2(S_MAX)=[    -1._R8P,    1._R8P,     0._R8P,   0._R8P,0._R8P   ]/2._R8P
real(R8P), parameter :: FVR2_CC_S3(S_MAX)=[    -6._R8P,    7._R8P,    -1._R8P,   0._R8P,0._R8P   ]/8._R8P
real(R8P), parameter :: FVR2_CC_S4(S_MAX)=[  -215._R8P,  273._R8P,   -65._R8P,   7._R8P,0._R8P   ]/240._R8P
real(R8P), parameter :: FVR2_CC_S5(S_MAX)=[-29960._R8P,40138._R8P,-12290._R8P,2317._R8P,-205._R8P]/30240._R8P
real(R8P), parameter :: FVR2_CC(S_MAX,S_MAX)=reshape([FVR2_CC_S1, &
                                                      FVR2_CC_S2, &
                                                      FVR2_CC_S3, &
                                                      FVR2_CC_S4, &
                                                      FVR2_CC_S5],&
                                                     [S_MAX,S_MAX]) !< FVR2_CC table.
!$acc declare copyin(FVR2_CC_S1,FVR2_CC_S2,FVR2_CC_S3,FVR2_CC_S4,FVR2_CC_S5,FVR2_CC)

!< FV face-reconstruction of derivative 3 (cell averages -> face i+1/2), antisymmetric.
real(R8P), parameter :: FVR3_CC_S1(S_MAX)=[     0._R8P,    0._R8P,    0._R8P,   0._R8P,  0._R8P]
real(R8P), parameter :: FVR3_CC_S2(S_MAX)=[    -3._R8P,    1._R8P,    0._R8P,   0._R8P,  0._R8P]
real(R8P), parameter :: FVR3_CC_S3(S_MAX)=[   -28._R8P,   11._R8P,   -1._R8P,   0._R8P,  0._R8P]/6._R8P
real(R8P), parameter :: FVR3_CC_S4(S_MAX)=[ -1365._R8P,  587._R8P,  -89._R8P,   7._R8P,  0._R8P]/240._R8P
real(R8P), parameter :: FVR3_CC_S5(S_MAX)=[-96327._R8P,43869._R8P,-8559._R8P,1179._R8P,-82._R8P]/15120._R8P
real(R8P), parameter :: FVR3_CC(S_MAX,S_MAX)=reshape([FVR3_CC_S1, &
                                                      FVR3_CC_S2, &
                                                      FVR3_CC_S3, &
                                                      FVR3_CC_S4, &
                                                      FVR3_CC_S5],&
                                                     [S_MAX,S_MAX]) !< FVR3_CC table.
!$acc declare copyin(FVR3_CC_S1,FVR3_CC_S2,FVR3_CC_S3,FVR3_CC_S4,FVR3_CC_S5,FVR3_CC)

!< FV face-reconstruction of derivative 4 (cell averages -> face i+1/2), symmetric.
real(R8P), parameter :: FVR4_CC_S1(S_MAX)=[  0._R8P,    0._R8P,  0._R8P,   0._R8P, 0._R8P]
real(R8P), parameter :: FVR4_CC_S2(S_MAX)=[  0._R8P,    0._R8P,  0._R8P,   0._R8P, 0._R8P]
real(R8P), parameter :: FVR4_CC_S3(S_MAX)=[  2._R8P,   -3._R8P,  1._R8P,   0._R8P, 0._R8P]/2._R8P
real(R8P), parameter :: FVR4_CC_S4(S_MAX)=[ 11._R8P,  -18._R8P,  8._R8P,  -1._R8P, 0._R8P]/6._R8P
real(R8P), parameter :: FVR4_CC_S5(S_MAX)=[710._R8P,-1228._R8P,644._R8P,-139._R8P,13._R8P]/288._R8P
real(R8P), parameter :: FVR4_CC(S_MAX,S_MAX)=reshape([FVR4_CC_S1, &
                                                      FVR4_CC_S2, &
                                                      FVR4_CC_S3, &
                                                      FVR4_CC_S4, &
                                                      FVR4_CC_S5],&
                                                     [S_MAX,S_MAX]) !< FVR4_CC table.
!$acc declare copyin(FVR4_CC_S1,FVR4_CC_S2,FVR4_CC_S3,FVR4_CC_S4,FVR4_CC_S5,FVR4_CC)

!< FV face-reconstruction of derivative 5 (cell averages -> face i+1/2), antisymmetric.
real(R8P), parameter :: FVR5_CC_S1(S_MAX)=[   0._R8P,    0._R8P,   0._R8P,   0._R8P, 0._R8P]
real(R8P), parameter :: FVR5_CC_S2(S_MAX)=[   0._R8P,    0._R8P,   0._R8P,   0._R8P, 0._R8P]
real(R8P), parameter :: FVR5_CC_S3(S_MAX)=[  10._R8P,   -5._R8P,   1._R8P,   0._R8P, 0._R8P]
real(R8P), parameter :: FVR5_CC_S4(S_MAX)=[  75._R8P,  -41._R8P,  11._R8P,  -1._R8P, 0._R8P]/4._R8P
real(R8P), parameter :: FVR5_CC_S5(S_MAX)=[6138._R8P,-3552._R8P,1128._R8P,-177._R8P,13._R8P]/240._R8P
real(R8P), parameter :: FVR5_CC(S_MAX,S_MAX)=reshape([FVR5_CC_S1, &
                                                      FVR5_CC_S2, &
                                                      FVR5_CC_S3, &
                                                      FVR5_CC_S4, &
                                                      FVR5_CC_S5],&
                                                     [S_MAX,S_MAX]) !< FVR5_CC table.
!$acc declare copyin(FVR5_CC_S1,FVR5_CC_S2,FVR5_CC_S3,FVR5_CC_S4,FVR5_CC_S5,FVR5_CC)

!< Deconvolution operators (McCorquodale & Colella, JCP 2011; Mignone, JCP 2014 Eq.52): convert between cell
!< averages and cell-center point values to high order, <q>_i = q_i + (h^2/24) d2q + ... . Symmetric, [c0,c1..cs]
!< layout. Order 2s+2. Generated by scripts/fdv_coefficients.py; do not hand-edit.
!< Deconvolution cell-average -> cell-center point value (symmetric).
real(R8P), parameter :: DECONV_A2P_S1(S_MAX+1)=[      26._R8P,     -1._R8P,     0._R8P,     0._R8P,   0._R8P,0._R8P]/24._R8P
real(R8P), parameter :: DECONV_A2P_S2(S_MAX+1)=[    2134._R8P,   -116._R8P,     9._R8P,     0._R8P,   0._R8P,0._R8P]/1920._R8P
real(R8P), parameter :: DECONV_A2P_S3(S_MAX+1)=[  121004._R8P,  -7621._R8P,   954._R8P,   -75._R8P,   0._R8P,0._R8P]/107520._R8P
real(R8P), parameter :: DECONV_A2P_S4(S_MAX+1)=[11702134._R8P,-800216._R8P,125884._R8P,-17000._R8P,1225._R8P,0._R8P]/10321920._R8P
real(R8P), parameter :: DECONV_A2P_S5(S_MAX+1)=[1034788732._R8P,-74586458._R8P,13459192._R8P,-2389025._R8P,306250._R8P,-19845._R8P]&
                                               /908328960._R8P
real(R8P), parameter :: DECONV_A2P(S_MAX+1,S_MAX)=reshape([DECONV_A2P_S1, &
                                                           DECONV_A2P_S2, &
                                                           DECONV_A2P_S3, &
                                                           DECONV_A2P_S4, &
                                                           DECONV_A2P_S5],&
                                                          [S_MAX+1,S_MAX]) !< DECONV_A2P table.
!$acc declare copyin(DECONV_A2P_S1,DECONV_A2P_S2,DECONV_A2P_S3,DECONV_A2P_S4,DECONV_A2P_S5,DECONV_A2P)

!< Deconvolution cell-center point value -> cell-average (symmetric).
real(R8P), parameter :: DECONV_P2A_S1(S_MAX+1)=[    22._R8P,    1._R8P,    0._R8P,  0._R8P,0._R8P,0._R8P]/24._R8P
real(R8P), parameter :: DECONV_P2A_S2(S_MAX+1)=[  5178._R8P,  308._R8P,  -17._R8P,  0._R8P,0._R8P,0._R8P]/5760._R8P
real(R8P), parameter :: DECONV_P2A_S3(S_MAX+1)=[862564._R8P,57249._R8P,-5058._R8P,367._R8P,0._R8P,0._R8P]/967680._R8P
real(R8P), parameter :: DECONV_P2A_S4(S_MAX+1)=[412080590._R8P,29039624._R8P,-3207892._R8P,399032._R8P,-27859._R8P,0._R8P]&
                                               /464486400._R8P
real(R8P), parameter :: DECONV_P2A_S5(S_MAX+1)=[108462733404._R8P,7938579366._R8P,-1002379848._R8P,163655583._R8P,-20312806._R8P,&
                                                1295803._R8P]/122624409600._R8P
real(R8P), parameter :: DECONV_P2A(S_MAX+1,S_MAX)=reshape([DECONV_P2A_S1, &
                                                           DECONV_P2A_S2, &
                                                           DECONV_P2A_S3, &
                                                           DECONV_P2A_S4, &
                                                           DECONV_P2A_S5],&
                                                          [S_MAX+1,S_MAX]) !< DECONV_P2A table.
!$acc declare copyin(DECONV_P2A_S1,DECONV_P2A_S2,DECONV_P2A_S3,DECONV_P2A_S4,DECONV_P2A_S5,DECONV_P2A)

!< Right-upwind schemes (left-upwind are mirrored).
!< | Order \(p\)| Stencil points \(m\)  | Coefficients \(a_m^{(p)}\) for \(q_{i+1/2}\)                   |
!< |------------|-----------------------|----------------------------------------------------------------|
!< | 1st  (S=1) | 0                     | 1/1   *(   1                                                 ) |
!< | 2nd  (S=2) | 0,1                   | 1/2   *(   3,    -1                                          ) |
!< | 3rd  (S=3) | 0,1,2                 | 1/6   *(  11,    -7,    2                                    ) |
!< | 4th  (S=4) | 0,1,2,3               | 1/12  *(  25,   -23,   13,    -3                             ) |
!< | 5th  (S=5) | 0,1,2,3,4             | 1/60  *( 137,  -163,  137,   -63,   12                       ) |
!< | 6th  (S=6) | 0,1,2,3,4,5           | 1/60  *( 147,  -213,  237,  -163,   62,   -10                ) |
!< | 7th  (S=7) | 0,1,2,3,4,5,6         | 1/420 *(1089, -1851, 2559, -2341, 1334,  -430,   60          ) |
!< | 8th  (S=8) | 0,1,2,3,4,5,6,7       | 1/840 *(2283, -4437, 7323, -8357, 6343, -3065,  855, -105    ) |
!< | 9th  (S=9) | 0,1,2,3,4,5,6,7,8     | 1/2520*(7129,-15551,29809,-40751,38629,-24875,10405,-2555,280) |
real(R8P), parameter :: FV1_UR_S1(S_MAX)=[  1._R8P,   0._R8P,  0._R8P,  0._R8P, 0._R8P]         !< FV1UR, S1.
real(R8P), parameter :: FV1_UR_S2(S_MAX)=[  3._R8P,  -1._R8P,  0._R8P,  0._R8P, 0._R8P]/2._R8P  !< FV1UR, S2.
real(R8P), parameter :: FV1_UR_S3(S_MAX)=[ 11._R8P,  -7._R8P,  2._R8P,  0._R8P, 0._R8P]/6._R8P  !< FV1UR, S3.
real(R8P), parameter :: FV1_UR_S4(S_MAX)=[ 25._R8P, -23._R8P, 13._R8P, -3._R8P, 0._R8P]/12._R8P !< FV1UR, S4.
real(R8P), parameter :: FV1_UR_S5(S_MAX)=[137._R8P,-163._R8P,137._R8P,-63._R8P,12._R8P]/60._R8P !< FV1UR, S5.
real(R8P), parameter :: FV1_UR(S_MAX,S_MAX)=reshape([FV1_UR_S1, &
                                                     FV1_UR_S2, &
                                                     FV1_UR_S3, &
                                                     FV1_UR_S4, &
                                                     FV1_UR_S5],&
                                                    [S_MAX,S_MAX]) !< Finite volume right-upwind reconstruction coefficients.
!$acc declare copyin(FV1_UR_S1,FV1_UR_S2,FV1_UR_S3,FV1_UR_S4,FV1_UR_S5,FV1_UR)
real(R8P), parameter :: FV1_UL_S1(S_MAX)=[ 1._R8P,  0._R8P,  0._R8P,   0._R8P,  0._R8P]         !< FV1UL, S1.
real(R8P), parameter :: FV1_UL_S2(S_MAX)=[-1._R8P,  3._R8P,  0._R8P,   0._R8P,  0._R8P]/2._R8P  !< FV1UL, S2.
real(R8P), parameter :: FV1_UL_S3(S_MAX)=[ 2._R8P, -7._R8P, 11._R8P,   0._R8P,  0._R8P]/6._R8P  !< FV1UL, S3.
real(R8P), parameter :: FV1_UL_S4(S_MAX)=[-3._R8P, 13._R8P,-23._R8P,  25._R8P,  0._R8P]/12._R8P !< FV1UL, S4.
real(R8P), parameter :: FV1_UL_S5(S_MAX)=[12._R8P,-63._R8P,137._R8P,-163._R8P,137._R8P]/60._R8P !< FV1UL, S5.
real(R8P), parameter :: FV1_UL(S_MAX,S_MAX)=reshape([FV1_UL_S1, &
                                                     FV1_UL_S2, &
                                                     FV1_UL_S3, &
                                                     FV1_UL_S4, &
                                                     FV1_UL_S5],&
                                                    [S_MAX,S_MAX]) !< Finite volume left-upwind reconstruction coefficients.
!$acc declare copyin(FV1_UL_S1,FV1_UL_S2,FV1_UL_S3,FV1_UL_S4,FV1_UL_S5,FV1_UL)

!< Upwind face-reconstruction of the m-th derivative tables (FVUR right-upwind, FVUL left-upwind), m=1..4.
!< Used by the upwind derivative ladder: dn|_i = (g_{i+1/2}-g_{i-1/2})/ds^n, g = FVU?{n-1} reconstruction.
!< Generated by scripts/fdv_coefficients.py (--emit-fortran-upwind); do not hand-edit.

!< FV right-upwind face-reconstruction of derivative 1 (one-sided, cells <= i).
real(R8P), parameter :: FVUR1_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR1_CC_S2(S_MAX)=[1._R8P,-1._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR1_CC_S3(S_MAX)=[2._R8P,-3._R8P,1._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR1_CC_S4(S_MAX)=[35._R8P,-69._R8P,45._R8P,-11._R8P,0._R8P]/12._R8P
real(R8P), parameter :: FVUR1_CC_S5(S_MAX)=[45._R8P,-109._R8P,105._R8P,-51._R8P,10._R8P]/12._R8P
real(R8P), parameter :: FVUR1_CC(S_MAX,S_MAX)=reshape([FVUR1_CC_S1, &
                                                     FVUR1_CC_S2, &
                                                     FVUR1_CC_S3, &
                                                     FVUR1_CC_S4, &
                                                     FVUR1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUR1_CC table.
!$acc declare copyin(FVUR1_CC_S1,FVUR1_CC_S2,FVUR1_CC_S3,FVUR1_CC_S4,FVUR1_CC_S5,FVUR1_CC)

!< FV right-upwind face-reconstruction of derivative 2 (one-sided, cells <= i).
real(R8P), parameter :: FVUR2_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR2_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR2_CC_S3(S_MAX)=[1._R8P,-2._R8P,1._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR2_CC_S4(S_MAX)=[5._R8P,-13._R8P,11._R8P,-3._R8P,0._R8P]/2._R8P
real(R8P), parameter :: FVUR2_CC_S5(S_MAX)=[17._R8P,-54._R8P,64._R8P,-34._R8P,7._R8P]/4._R8P
real(R8P), parameter :: FVUR2_CC(S_MAX,S_MAX)=reshape([FVUR2_CC_S1, &
                                                     FVUR2_CC_S2, &
                                                     FVUR2_CC_S3, &
                                                     FVUR2_CC_S4, &
                                                     FVUR2_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUR2_CC table.
!$acc declare copyin(FVUR2_CC_S1,FVUR2_CC_S2,FVUR2_CC_S3,FVUR2_CC_S4,FVUR2_CC_S5,FVUR2_CC)

!< FV right-upwind face-reconstruction of derivative 3 (one-sided, cells <= i).
real(R8P), parameter :: FVUR3_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR3_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR3_CC_S3(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR3_CC_S4(S_MAX)=[1._R8P,-3._R8P,3._R8P,-1._R8P,0._R8P]
real(R8P), parameter :: FVUR3_CC_S5(S_MAX)=[3._R8P,-11._R8P,15._R8P,-9._R8P,2._R8P]
real(R8P), parameter :: FVUR3_CC(S_MAX,S_MAX)=reshape([FVUR3_CC_S1, &
                                                     FVUR3_CC_S2, &
                                                     FVUR3_CC_S3, &
                                                     FVUR3_CC_S4, &
                                                     FVUR3_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUR3_CC table.
!$acc declare copyin(FVUR3_CC_S1,FVUR3_CC_S2,FVUR3_CC_S3,FVUR3_CC_S4,FVUR3_CC_S5,FVUR3_CC)

!< FV right-upwind face-reconstruction of derivative 4 (one-sided, cells <= i).
real(R8P), parameter :: FVUR4_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR4_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR4_CC_S3(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR4_CC_S4(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUR4_CC_S5(S_MAX)=[1._R8P,-4._R8P,6._R8P,-4._R8P,1._R8P]
real(R8P), parameter :: FVUR4_CC(S_MAX,S_MAX)=reshape([FVUR4_CC_S1, &
                                                     FVUR4_CC_S2, &
                                                     FVUR4_CC_S3, &
                                                     FVUR4_CC_S4, &
                                                     FVUR4_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUR4_CC table.
!$acc declare copyin(FVUR4_CC_S1,FVUR4_CC_S2,FVUR4_CC_S3,FVUR4_CC_S4,FVUR4_CC_S5,FVUR4_CC)

!< FV left-upwind face-reconstruction of derivative 1 (one-sided, cells >= i+1).
real(R8P), parameter :: FVUL1_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL1_CC_S2(S_MAX)=[1._R8P,-1._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL1_CC_S3(S_MAX)=[-1._R8P,3._R8P,-2._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL1_CC_S4(S_MAX)=[11._R8P,-45._R8P,69._R8P,-35._R8P,0._R8P]/12._R8P
real(R8P), parameter :: FVUL1_CC_S5(S_MAX)=[-10._R8P,51._R8P,-105._R8P,109._R8P,-45._R8P]/12._R8P
real(R8P), parameter :: FVUL1_CC(S_MAX,S_MAX)=reshape([FVUL1_CC_S1, &
                                                     FVUL1_CC_S2, &
                                                     FVUL1_CC_S3, &
                                                     FVUL1_CC_S4, &
                                                     FVUL1_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUL1_CC table.
!$acc declare copyin(FVUL1_CC_S1,FVUL1_CC_S2,FVUL1_CC_S3,FVUL1_CC_S4,FVUL1_CC_S5,FVUL1_CC)

!< FV left-upwind face-reconstruction of derivative 2 (one-sided, cells >= i+1).
real(R8P), parameter :: FVUL2_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL2_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL2_CC_S3(S_MAX)=[1._R8P,-2._R8P,1._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL2_CC_S4(S_MAX)=[-3._R8P,11._R8P,-13._R8P,5._R8P,0._R8P]/2._R8P
real(R8P), parameter :: FVUL2_CC_S5(S_MAX)=[7._R8P,-34._R8P,64._R8P,-54._R8P,17._R8P]/4._R8P
real(R8P), parameter :: FVUL2_CC(S_MAX,S_MAX)=reshape([FVUL2_CC_S1, &
                                                     FVUL2_CC_S2, &
                                                     FVUL2_CC_S3, &
                                                     FVUL2_CC_S4, &
                                                     FVUL2_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUL2_CC table.
!$acc declare copyin(FVUL2_CC_S1,FVUL2_CC_S2,FVUL2_CC_S3,FVUL2_CC_S4,FVUL2_CC_S5,FVUL2_CC)

!< FV left-upwind face-reconstruction of derivative 3 (one-sided, cells >= i+1).
real(R8P), parameter :: FVUL3_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL3_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL3_CC_S3(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL3_CC_S4(S_MAX)=[1._R8P,-3._R8P,3._R8P,-1._R8P,0._R8P]
real(R8P), parameter :: FVUL3_CC_S5(S_MAX)=[-2._R8P,9._R8P,-15._R8P,11._R8P,-3._R8P]
real(R8P), parameter :: FVUL3_CC(S_MAX,S_MAX)=reshape([FVUL3_CC_S1, &
                                                     FVUL3_CC_S2, &
                                                     FVUL3_CC_S3, &
                                                     FVUL3_CC_S4, &
                                                     FVUL3_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUL3_CC table.
!$acc declare copyin(FVUL3_CC_S1,FVUL3_CC_S2,FVUL3_CC_S3,FVUL3_CC_S4,FVUL3_CC_S5,FVUL3_CC)

!< FV left-upwind face-reconstruction of derivative 4 (one-sided, cells >= i+1).
real(R8P), parameter :: FVUL4_CC_S1(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL4_CC_S2(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL4_CC_S3(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL4_CC_S4(S_MAX)=[0._R8P,0._R8P,0._R8P,0._R8P,0._R8P]
real(R8P), parameter :: FVUL4_CC_S5(S_MAX)=[1._R8P,-4._R8P,6._R8P,-4._R8P,1._R8P]
real(R8P), parameter :: FVUL4_CC(S_MAX,S_MAX)=reshape([FVUL4_CC_S1, &
                                                     FVUL4_CC_S2, &
                                                     FVUL4_CC_S3, &
                                                     FVUL4_CC_S4, &
                                                     FVUL4_CC_S5],&
                                                    [S_MAX,S_MAX]) !< FVUL4_CC table.
!$acc declare copyin(FVUL4_CC_S1,FVUL4_CC_S2,FVUL4_CC_S3,FVUL4_CC_S4,FVUL4_CC_S5,FVUL4_CC)

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
   !$acc routine seq

   dq_ds = 0.0_R8P
   do m=1, s
      dq_ds = dq_ds + FD1_CC(m,s)*(q(1+m) - q(1-m))
   enddo
   dq_ds = dq_ds/ds
   endsubroutine compute_derivative1_fd_centered

   pure subroutine compute_derivative1_fd_seam(samples, ds, side, dq_ds)
   !< Issue #29 E3: matched-resolution SEAM first derivative along the normal, for a
   !< seam-adjacent cell at a 2:1 coarse-fine jump. Restores div_h(curl_h)=0
   !< telescoping to order 6 across the jump (E3 T1.1). `samples(1:7)` are the field
   !< values the CALLER gathered along the outward seam-normal, in TABLE ORDER:
   !<   side=+1 (FINE side): [q0, q(-1), q(-2), q(-3), Qc0, Qc1, Qc2] — fine interior
   !<     0..-3 then coarse donor centers; `ds` = FINE spacing.
   !<   side=-1 (COARSE side): [q0, q(-1), q(-2), q(-3), q(-4), Qf0, Qf1] — coarse
   !<     interior 0..-4 then fine ghost centers; `ds` = COARSE spacing.
   !< The weight tables (FD1_SEAM_FINE / FD1_SEAM_COARSE) already carry the geometry;
   !< the caller supplies the samples in order and the local spacing. Same div/curl
   !< primitive both sides (one operator), so the composite telescopes.
   real(R8P),    intent(in)  :: samples(1:7) !< Normal-line samples in table order (see above).
   real(R8P),    intent(in)  :: ds           !< Local spacing (FINE for side=+1, COARSE for side=-1).
   integer(I4P), intent(in)  :: side         !< +1 = fine-side seam cell, -1 = coarse-side seam cell.
   real(R8P),    intent(out) :: dq_ds        !< Matched seam derivative d/ds.
   integer(I4P)              :: n            !< Sample counter.
   !$acc routine seq

   dq_ds = 0.0_R8P
   if (side > 0_I4P) then
      do n=1, 7
         dq_ds = dq_ds + FD1_SEAM_FINE(n) * samples(n)
      enddo
   else
      do n=1, 7
         dq_ds = dq_ds + FD1_SEAM_COARSE(n) * samples(n)
      enddo
   endif
   dq_ds = dq_ds/ds
   endsubroutine compute_derivative1_fd_seam

   pure subroutine compute_derivative2_fd_centered(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite difference centered scheme.
   integer(I4P), intent(in)  :: s       !< Stencil len, order=2*s.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

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
   !$acc routine seq

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
   !$acc routine seq

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
   !$acc routine seq

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
   !$acc routine seq

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
   !$acc routine seq

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

   pure subroutine compute_reconstruction_r_fd_centered(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values. Centered schemes.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction at right interface of field.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FD0_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fd_centered

   ! finite volume schemes
   ! centered
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
   !$acc routine seq

   call compute_reconstruction_r_fv_centered(s=s,q=q(1-s:  s),qr=ql)
   call compute_reconstruction_r_fv_centered(s=s,q=q(2-s:1+s),qr=qr)
   dq_ds = (qr-ql)/ds
   endsubroutine compute_derivative1_fv_centered

   pure subroutine compute_derivative2_fv_centered(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite volume centered scheme.
   !< SOTA: flux difference of the d1 face reconstruction, d2|_i = (g_{i+1/2}-g_{i-1/2})/ds with g=FVR1.
   !< Returns the cell-AVERAGED second derivative (true FV semantics); apply compute_avg_to_point for the
   !< point value at high order.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   real(R8P)                 :: gl,gr   !< d1 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_centered_d1(s=s,q=q(1-s:  s),qr=gl)
   call compute_reconstruction_r_fv_centered_d1(s=s,q=q(2-s:1+s),qr=gr)
   d2q_ds2 = (gr-gl)/(ds*ds)
   endsubroutine compute_derivative2_fv_centered

   pure subroutine compute_derivative3_fv_centered(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 with finite volume centered scheme.
   !< SOTA: flux difference of the d2 face reconstruction (g=FVR2). Cell-averaged third derivative.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   real(R8P)                 :: gl,gr   !< d2 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_centered_d2(s=s,q=q(1-s:  s),qr=gl)
   call compute_reconstruction_r_fv_centered_d2(s=s,q=q(2-s:1+s),qr=gr)
   d3q_ds3 = (gr-gl)/(ds*ds*ds)
   endsubroutine compute_derivative3_fv_centered

   pure subroutine compute_derivative4_fv_centered(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 with finite volume centered scheme.
   !< SOTA: flux difference of the d3 face reconstruction (g=FVR3). Cell-averaged fourth derivative.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   real(R8P)                 :: gl,gr   !< d3 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_centered_d3(s=s,q=q(1-s:  s),qr=gl)
   call compute_reconstruction_r_fv_centered_d3(s=s,q=q(2-s:1+s),qr=gr)
   d4q_ds4 = (gr-gl)/(ds*ds*ds*ds)
   endsubroutine compute_derivative4_fv_centered

   pure subroutine compute_derivative5_fv_centered(s,ds,q,d5q_ds5)
   !< Compute derivative of order 5 with finite volume centered scheme.
   !< SOTA: flux difference of the d4 face reconstruction (g=FVR4). Cell-averaged fifth derivative.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   real(R8P)                 :: gl,gr   !< d4 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_centered_d4(s=s,q=q(1-s:  s),qr=gl)
   call compute_reconstruction_r_fv_centered_d4(s=s,q=q(2-s:1+s),qr=gr)
   d5q_ds5 = (gr-gl)/(ds*ds*ds*ds*ds)
   endsubroutine compute_derivative5_fv_centered

   pure subroutine compute_derivative6_fv_centered(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 with finite volume centered scheme.
   !< SOTA: flux difference of the d5 face reconstruction (g=FVR5). Cell-averaged sixth derivative.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:1+s].
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q, d6q/ds6.
   real(R8P)                 :: gl,gr   !< d5 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_centered_d5(s=s,q=q(1-s:  s),qr=gl)
   call compute_reconstruction_r_fv_centered_d5(s=s,q=q(2-s:1+s),qr=gr)
   d6q_ds6 = (gr-gl)/(ds*ds*ds*ds*ds*ds)
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
   !< Compute reconstruction at right interface from cell center average values. Centered schemes.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction at right interface of field.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FV1_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered

   pure subroutine compute_reconstruction_r_fv_centered_d1(s,q,qr)
   !< Reconstruct d1 q at the right interface i+1/2 from cell averages, centered (antisymmetric, FVR1_CC).
   !< Specialized branch-free kernel (one per derivative order): the FV derivative ladder differences two
   !< adjacent face reconstructions. Avoids per-call dispatch in hot device loops.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction of d1 at right interface.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FVR1_CC(m,s)*(q(m) - q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered_d1

   pure subroutine compute_reconstruction_r_fv_centered_d2(s,q,qr)
   !< Reconstruct d2 q at the right interface i+1/2 from cell averages, centered (symmetric, FVR2_CC).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction of d2 at right interface.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FVR2_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered_d2

   pure subroutine compute_reconstruction_r_fv_centered_d3(s,q,qr)
   !< Reconstruct d3 q at the right interface i+1/2 from cell averages, centered (antisymmetric, FVR3_CC).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction of d3 at right interface.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FVR3_CC(m,s)*(q(m) - q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered_d3

   pure subroutine compute_reconstruction_r_fv_centered_d4(s,q,qr)
   !< Reconstruct d4 q at the right interface i+1/2 from cell averages, centered (symmetric, FVR4_CC).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction of d4 at right interface.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FVR4_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered_d4

   pure subroutine compute_reconstruction_r_fv_centered_d5(s,q,qr)
   !< Reconstruct d5 q at the right interface i+1/2 from cell averages, centered (antisymmetric, FVR5_CC).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr      !< Reconstruction of d5 at right interface.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FVR5_CC(m,s)*(q(m) - q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fv_centered_d5

   pure subroutine compute_avg_to_point(s,q,qp)
   !< Deconvolution: recover the cell-center POINT value from cell AVERAGES, to order 2s+2.
   !< q_i = <q>_i + (h^2/24) d2<q> + ...  (McCorquodale & Colella 2011; Mignone 2014, Eq.52).
   !< Apply after an FV (cell-averaged) derivative to obtain its high-order point value.
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Cell averages over the stencil [1-s:1+s] (centered on cell i at index 1).
   real(R8P),    intent(out) :: qp      !< Cell-center point value at cell i.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qp = DECONV_A2P(1,s)*q(1)
   do m=1, s
      qp = qp + DECONV_A2P(m+1,s)*(q(1+m) + q(1-m))
   enddo
   endsubroutine compute_avg_to_point

   pure subroutine compute_point_to_avg(s,q,qa)
   !< Deconvolution: recover the cell AVERAGE from cell-center POINT values, to order 2s+2.
   !< <q>_i = q_i + (h^2/24) d2 q + ...  (inverse of compute_avg_to_point).
   integer(I4P), intent(in)  :: s       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Cell-center point values over the stencil [1-s:1+s] (cell i at index 1).
   real(R8P),    intent(out) :: qa      !< Cell average at cell i.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qa = DECONV_P2A(1,s)*q(1)
   do m=1, s
      qa = qa + DECONV_P2A(m+1,s)*(q(1+m) + q(1-m))
   enddo
   endsubroutine compute_point_to_avg

   ! right-upwind
   pure subroutine compute_derivative1_fv_rupwind(s,ds,q,dq_ds)
   !< Compute derivative of order 1 with finite volume right-upwind scheme.
   integer(I4P), intent(in)  :: s     !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: ds    !< Space step.
   real(R8P),    intent(in)  :: q(0:) !< Scalar field over the stencil [0:1+s].
   real(R8P),    intent(out) :: dq_ds !< Derivative of order 1 of q, dq/ds.
   real(R8P)                 :: ql,qr !< Reconstruction of field at left and righ interfaces.
   !$acc routine seq

   call compute_reconstruction_r_fv_rupwind(s=s,q=q(0:  s),qr=ql)
   call compute_reconstruction_r_fv_rupwind(s=s,q=q(1:1+s),qr=qr)
   dq_ds = (qr-ql)/ds
   endsubroutine compute_derivative1_fv_rupwind

   pure subroutine compute_derivative2_fv_rupwind(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite volume right-upwind scheme (SOTA flux difference).
   !< d2|_i = (g_{i+1/2}-g_{i-1/2})/ds^2, g = d1 face reconstruction (FVUR1). Cell-averaged d2.
   !< Window q(0:s): cell i at q(1); g_{i+1/2} uses q(1:s) (cells i..i-s+1), g_{i-1/2} uses q(0:s-1).
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< One-sided stencil [1-s:1], cell i at q(1) (reaches cells i-s..i).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   real(R8P)                 :: gl,gr   !< d1 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_rupwind_d1(s=s,q=q(2-s:1),qr=gr)
   call compute_reconstruction_r_fv_rupwind_d1(s=s,q=q(1-s:0),qr=gl)
   d2q_ds2 = (gr-gl)/(ds*ds)
   endsubroutine compute_derivative2_fv_rupwind

   pure subroutine compute_derivative3_fv_rupwind(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 with finite volume right-upwind scheme (flux diff of FVUR2). Cell-avg d3.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< One-sided stencil [1-s:1], cell i at q(1) (reaches cells i-s..i).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   real(R8P)                 :: gl,gr   !< d2 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_rupwind_d2(s=s,q=q(2-s:1),qr=gr)
   call compute_reconstruction_r_fv_rupwind_d2(s=s,q=q(1-s:0),qr=gl)
   d3q_ds3 = (gr-gl)/(ds*ds*ds)
   endsubroutine compute_derivative3_fv_rupwind

   pure subroutine compute_derivative4_fv_rupwind(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 with finite volume right-upwind scheme (flux diff of FVUR3). Cell-avg d4.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< One-sided stencil [1-s:1], cell i at q(1) (reaches cells i-s..i).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   real(R8P)                 :: gl,gr   !< d3 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_rupwind_d3(s=s,q=q(2-s:1),qr=gr)
   call compute_reconstruction_r_fv_rupwind_d3(s=s,q=q(1-s:0),qr=gl)
   d4q_ds4 = (gr-gl)/(ds*ds*ds*ds)
   endsubroutine compute_derivative4_fv_rupwind

   pure subroutine compute_derivative5_fv_rupwind(s,ds,q,d5q_ds5)
   !< Compute derivative of order 5 with finite volume right-upwind scheme (flux diff of FVUR4). Cell-avg d5.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< One-sided stencil [1-s:1], cell i at q(1) (reaches cells i-s..i).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   real(R8P)                 :: gl,gr   !< d4 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_rupwind_d4(s=s,q=q(2-s:1),qr=gr)
   call compute_reconstruction_r_fv_rupwind_d4(s=s,q=q(1-s:0),qr=gl)
   d5q_ds5 = (gr-gl)/(ds*ds*ds*ds*ds)
   endsubroutine compute_derivative5_fv_rupwind

   pure subroutine compute_derivative6_fv_rupwind(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 with finite volume right-upwind scheme.
   !< UNSUPPORTED at FDV_S_MAX=5: the SOTA construction needs the m=5 face reconstruction, a 6-cell
   !< one-sided stencil (FVUR5), which requires S_MAX>=6. Returns 0 as an inert sentinel; raise S_MAX
   !< and add FVUR5_CC to enable. (Not error_stop: this routine is pure / device-callable.)
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< Scalar field over the one-sided stencil.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q (returns 0: unsupported at S_MAX=5).
   !$acc routine seq

   associate(unused_s=>s, unused_q=>q, unused_ds=>ds); end associate
   d6q_ds6 = 0.0_R8P
   endsubroutine compute_derivative6_fv_rupwind

   pure subroutine compute_reconstruction_r_fv_rupwind(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values, left-upwind schemes.
   integer(I4P), intent(in)  :: s     !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:) !< Scalar field over the stencil [0:s].
   real(R8P),    intent(out) :: qr    !< Reconstruction at right interface of field.
   integer(I4P)              :: m     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FV1_UR(m,s)*q(m)
   enddo
   endsubroutine compute_reconstruction_r_fv_rupwind

   pure subroutine compute_reconstruction_r_fv_rupwind_d1(s,q,qr)
   !< Right-upwind reconstruction of d1 q at the face right of the at-face cell (FVUR1_CC).
   !< Specialized branch-free kernel. Window q(1:s) one-sided: q(k) = cell c-(k-1) (c, c-1, ..., c-(s-1)).
   !< The upwind d2 ladder flux-differences two adjacent of these.
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages, one-sided window: q(k)=cell c-(k-1).
   real(R8P),    intent(out) :: qr    !< Reconstruction of d1 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUR1_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_rupwind_d1

   pure subroutine compute_reconstruction_r_fv_rupwind_d2(s,q,qr)
   !< Right-upwind reconstruction of d2 q at the face (FVUR2_CC). Window q(1:s): q(k)=cell c-(k-1).
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d2 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUR2_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_rupwind_d2

   pure subroutine compute_reconstruction_r_fv_rupwind_d3(s,q,qr)
   !< Right-upwind reconstruction of d3 q at the face (FVUR3_CC). Window q(1:s): q(k)=cell c-(k-1).
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d3 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUR3_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_rupwind_d3

   pure subroutine compute_reconstruction_r_fv_rupwind_d4(s,q,qr)
   !< Right-upwind reconstruction of d4 q at the face (FVUR4_CC). Window q(1:s): q(k)=cell c-(k-1).
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d4 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUR4_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_rupwind_d4

   ! left-upwind
   pure subroutine compute_derivative1_fv_lupwind(s,ds,q,dq_ds)
   !< Compute derivative of order 1 with finite volume left-upwind scheme.
   integer(I4P), intent(in)  :: s      !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: ds     !< Space step.
   real(R8P),    intent(in)  :: q(-s:) !< Scalar field over the stencil [-s:0].
   real(R8P),    intent(out) :: dq_ds  !< Derivative of order 1 of q, dq/ds.
   real(R8P)                 :: ql,qr  !< Reconstruction of field at left and righ interfaces.
   !$acc routine seq

   call compute_reconstruction_r_fv_lupwind(s=s,q=q(0-s:-1),qr=ql)
   call compute_reconstruction_r_fv_lupwind(s=s,q=q(1-s:0 ),qr=qr)
   dq_ds = (qr-ql)/ds
   endsubroutine compute_derivative1_fv_lupwind

   pure subroutine compute_derivative2_fv_lupwind(s,ds,q,d2q_ds2)
   !< Compute derivative of order 2 with finite volume left-upwind scheme (SOTA flux difference).
   !< d2|_i = (g_{i+1/2}-g_{i-1/2})/ds^2, g = d1 face reconstruction (FVUL1). Cell-averaged d2.
   !< Window q(0:s): cell i at q(0); g_{i+1/2} uses q(1:s) (cells i+1..i+s), g_{i-1/2} uses q(0:s-1).
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< One-sided stencil [0:s] ahead, cell i at q(0).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d2q_ds2 !< Derivative of order 2 of q, d2q/ds2.
   real(R8P)                 :: gl,gr   !< d1 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_lupwind_d1(s=s,q=q(1:s),qr=gr)
   call compute_reconstruction_r_fv_lupwind_d1(s=s,q=q(0:s-1),qr=gl)
   d2q_ds2 = (gr-gl)/(ds*ds)
   endsubroutine compute_derivative2_fv_lupwind

   pure subroutine compute_derivative3_fv_lupwind(s,ds,q,d3q_ds3)
   !< Compute derivative of order 3 with finite volume left-upwind scheme (flux diff of FVUL2). Cell-avg d3.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< One-sided stencil [0:s] ahead, cell i at q(0).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d3q_ds3 !< Derivative of order 3 of q, d3q/ds3.
   real(R8P)                 :: gl,gr   !< d2 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_lupwind_d2(s=s,q=q(1:s),qr=gr)
   call compute_reconstruction_r_fv_lupwind_d2(s=s,q=q(0:s-1),qr=gl)
   d3q_ds3 = (gr-gl)/(ds*ds*ds)
   endsubroutine compute_derivative3_fv_lupwind

   pure subroutine compute_derivative4_fv_lupwind(s,ds,q,d4q_ds4)
   !< Compute derivative of order 4 with finite volume left-upwind scheme (flux diff of FVUL3). Cell-avg d4.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< One-sided stencil [0:s] ahead, cell i at q(0).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d4q_ds4 !< Derivative of order 4 of q, d4q/ds4.
   real(R8P)                 :: gl,gr   !< d3 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_lupwind_d3(s=s,q=q(1:s),qr=gr)
   call compute_reconstruction_r_fv_lupwind_d3(s=s,q=q(0:s-1),qr=gl)
   d4q_ds4 = (gr-gl)/(ds*ds*ds*ds)
   endsubroutine compute_derivative4_fv_lupwind

   pure subroutine compute_derivative5_fv_lupwind(s,ds,q,d5q_ds5)
   !< Compute derivative of order 5 with finite volume left-upwind scheme (flux diff of FVUL4). Cell-avg d5.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< One-sided stencil [0:s] ahead, cell i at q(0).
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d5q_ds5 !< Derivative of order 5 of q, d5q/ds5.
   real(R8P)                 :: gl,gr   !< d4 reconstructed at left/right faces.
   !$acc routine seq

   call compute_reconstruction_r_fv_lupwind_d4(s=s,q=q(1:s),qr=gr)
   call compute_reconstruction_r_fv_lupwind_d4(s=s,q=q(0:s-1),qr=gl)
   d5q_ds5 = (gr-gl)/(ds*ds*ds*ds*ds)
   endsubroutine compute_derivative5_fv_lupwind

   pure subroutine compute_derivative6_fv_lupwind(s,ds,q,d6q_ds6)
   !< Compute derivative of order 6 with finite volume left-upwind scheme.
   !< UNSUPPORTED at FDV_S_MAX=5: needs the m=5 face reconstruction (FVUL5, 6-cell one-sided), i.e. S_MAX>=6.
   !< Returns 0 as an inert sentinel; raise S_MAX and add FVUL5_CC to enable. (pure / device-callable.)
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(0:)   !< Scalar field over the one-sided stencil.
   real(R8P),    intent(in)  :: ds      !< Space step.
   real(R8P),    intent(out) :: d6q_ds6 !< Derivative of order 6 of q (returns 0: unsupported at S_MAX=5).
   !$acc routine seq

   associate(unused_s=>s, unused_q=>q, unused_ds=>ds); end associate
   d6q_ds6 = 0.0_R8P
   endsubroutine compute_derivative6_fv_lupwind

   pure subroutine compute_reconstruction_r_fv_lupwind(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values, left-upwind schemes.
   integer(I4P), intent(in)  :: s       !< Stencil len, accuracy order.
   real(R8P),    intent(in)  :: q(1-s:) !< Scalar field over the stencil [1-s:0].
   real(R8P),    intent(out) :: qr      !< Reconstruction at right interface of field.
   integer(I4P)              :: m       !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FV1_UL(m,s)*q(1-m)
   enddo
   endsubroutine compute_reconstruction_r_fv_lupwind

   pure subroutine compute_reconstruction_r_fv_lupwind_d1(s,q,qr)
   !< Left-upwind reconstruction of d1 q at the face from one-sided cells ahead {c+1..c+s} (FVUL1_CC).
   !< Specialized branch-free kernel. Ascending window q(1:s): q(1)=cell c+1, q(s)=cell c+s.
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages, one-sided window ahead: q(1)=cell c+1.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d1 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUL1_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_lupwind_d1

   pure subroutine compute_reconstruction_r_fv_lupwind_d2(s,q,qr)
   !< Left-upwind reconstruction of d2 q at the face (FVUL2_CC). Ascending window q(1:s): q(1)=cell c+1.
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window ahead.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d2 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUL2_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_lupwind_d2

   pure subroutine compute_reconstruction_r_fv_lupwind_d3(s,q,qr)
   !< Left-upwind reconstruction of d3 q at the face (FVUL3_CC). Ascending window q(1:s): q(1)=cell c+1.
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window ahead.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d3 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUL3_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_lupwind_d3

   pure subroutine compute_reconstruction_r_fv_lupwind_d4(s,q,qr)
   !< Left-upwind reconstruction of d4 q at the face (FVUL4_CC). Ascending window q(1:s): q(1)=cell c+1.
   integer(I4P), intent(in)  :: s     !< Stencil len.
   real(R8P),    intent(in)  :: q(1:) !< Cell averages over the one-sided window ahead.
   real(R8P),    intent(out) :: qr    !< Reconstruction of d4 at the face.
   integer(I4P)              :: k     !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do k=1, s
      qr = qr + FVUL4_CC(k,s)*q(s-k+1)
   enddo
   endsubroutine compute_reconstruction_r_fv_lupwind_d4
endmodule adam_fdv_operators_library
