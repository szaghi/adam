!< ADAM, pointwise algebraic subgrid-scale models (issue #119, T1), shared by the CPU loops and the FNL device kernels.
module adam_turbulence_library
!< ADAM, pointwise algebraic subgrid-scale models (issue #119, T1), shared by the CPU loops and the FNL device kernels.
!<
!< Every routine returns the kinematic eddy viscosity `nu_t` from the velocity gradient `g(a,b) = d u_a / d x_b`, the constant
!< `c`, the grid steps `dl` and the filter width `delta`. It is `pure`, takes explicit-size dummies and is tagged
!< `!$acc routine seq` + `!$omp declare target`. The formulas come from the literature through the llm-wiki bundle
!< `turbulence-modelling`; their constants were calibrated on other discretizations, not on WENO.
!<
!< * Vreman (2004): `nu_t = c sqrt(B / (a_ij a_ij))`, `a_ij = d u_j / d x_i`, `b_ij = sum_m dl_m^2 a_mi a_mj`,
!<   `B = b_11 b_22 - b_12^2 + b_11 b_33 - b_13^2 + b_22 b_33 - b_23^2`; `c = 2.5 C_s^2`. Zero where `a_ij` has rank 1, pure shear among them.
!< * WALE (Nicoud & Ducros 1999): `nu_t = (c delta)^2 (Sd:Sd)^(3/2) / ((S:S)^(5/2) + (Sd:Sd)^(5/4))`, `S` the strain, `Sd` the
!<   traceless symmetric part of `g g`. Zero in pure shear.
!< * Sigma (Nicoud et al. 2011): `nu_t = (c delta)^2 s_3 (s_1 - s_2) (s_2 - s_3) / s_1^2`, `s_i` the singular values of `g`, from
!<   the invariants of `g^T g` in closed form. Zero in 2-D flow, pure shear and solid rotation.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: sgs_vreman
public :: sgs_wale
public :: sgs_sigma

contains
   pure subroutine sgs_vreman(c, dl, delta, g, nu)
   !< Vreman eddy viscosity.
   real(R8P), intent(in)  :: c        !< Model constant, 2.5 C_s^2.
   real(R8P), intent(in)  :: dl(3)    !< Grid steps.
   real(R8P), intent(in)  :: delta    !< Filter width (unused: the model reads the steps).
   real(R8P), intent(in)  :: g(3,3)   !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R8P), intent(out) :: nu       !< Eddy viscosity.
   real(R8P)              :: b(3,3)   !< Tensor beta.
   real(R8P)              :: bb       !< Invariant B of beta.
   real(R8P)              :: aa       !< a_ij a_ij.
   integer(I4P), parameter :: N = 3_I4P !< Dimensions.
   integer(I4P)            :: i, j, m  !< Counters.
   !$acc routine seq
   !$omp declare target

   aa = 0._R8P
   do j=1, N
      do i=1, N
         b(i,j) = 0._R8P
         aa = aa + g(i,j) * g(i,j)
         do m=1, N
            b(i,j) = b(i,j) + dl(m) * dl(m) * g(i,m) * g(j,m)
         enddo
      enddo
   enddo
   bb = b(1,1) * b(2,2) - b(1,2) * b(1,2) + b(1,1) * b(3,3) - b(1,3) * b(1,3) + b(2,2) * b(3,3) - b(2,3) * b(2,3)
   nu = 0._R8P
   if (aa > 0._R8P .and. bb > 0._R8P) nu = c * sqrt(bb / aa)
   endsubroutine sgs_vreman

   pure subroutine sgs_wale(c, dl, delta, g, nu)
   !< WALE eddy viscosity.
   real(R8P), intent(in)  :: c        !< Model constant C_w.
   real(R8P), intent(in)  :: dl(3)    !< Grid steps (unused).
   real(R8P), intent(in)  :: delta    !< Filter width.
   real(R8P), intent(in)  :: g(3,3)   !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R8P), intent(out) :: nu       !< Eddy viscosity.
   real(R8P)              :: gg(3,3)  !< g g.
   real(R8P)              :: tr       !< Trace of g g over 3.
   real(R8P)              :: ss, sd   !< S:S and Sd:Sd.
   real(R8P)              :: sij, dij !< Components of S and Sd.
   real(R8P)              :: den      !< Denominator.
   integer(I4P), parameter :: N = 3_I4P !< Dimensions.
   integer(I4P)            :: i, j, m  !< Counters.
   !$acc routine seq
   !$omp declare target

   do j=1, N
      do i=1, N
         gg(i,j) = 0._R8P
         do m=1, N
            gg(i,j) = gg(i,j) + g(i,m) * g(m,j)
         enddo
      enddo
   enddo
   tr = (gg(1,1) + gg(2,2) + gg(3,3)) / 3._R8P
   ss = 0._R8P ; sd = 0._R8P
   do j=1, N
      do i=1, N
         sij = 0.5_R8P * (g(i,j) + g(j,i))
         dij = 0.5_R8P * (gg(i,j) + gg(j,i))
         if (i == j) dij = dij - tr
         ss = ss + sij * sij
         sd = sd + dij * dij
      enddo
   enddo
   den = ss**2.5_R8P + sd**1.25_R8P
   nu = 0._R8P
   if (den > 0._R8P) nu = (c * delta)**2 * sd**1.5_R8P / den
   endsubroutine sgs_wale

   pure subroutine sgs_sigma(c, dl, delta, g, nu)
   !< Sigma eddy viscosity.
   real(R8P), intent(in)  :: c                !< Model constant C_sigma.
   real(R8P), intent(in)  :: dl(3)            !< Grid steps (unused).
   real(R8P), intent(in)  :: delta            !< Filter width.
   real(R8P), intent(in)  :: g(3,3)           !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R8P), intent(out) :: nu               !< Eddy viscosity.
   real(R8P), parameter   :: PI3 = 4._R8P * atan(1._R8P) / 3._R8P !< pi / 3.
   real(R8P)              :: t(3,3)           !< g^T g.
   real(R8P)              :: i1, i2, i3       !< Invariants of t.
   real(R8P)              :: a1, a2, a3       !< Auxiliaries of the closed form.
   real(R8P)              :: s1, s2, s3       !< Singular values, s1 >= s2 >= s3.
   integer(I4P), parameter :: N = 3_I4P        !< Dimensions.
   integer(I4P)            :: i, j, m          !< Counters.
   !$acc routine seq
   !$omp declare target

   do j=1, N
      do i=1, N
         t(i,j) = 0._R8P
         do m=1, N
            t(i,j) = t(i,j) + g(m,i) * g(m,j)
         enddo
      enddo
   enddo
   i1 = t(1,1) + t(2,2) + t(3,3)
   i2 = t(1,1) * t(2,2) - t(1,2) * t(1,2) + t(1,1) * t(3,3) - t(1,3) * t(1,3) + t(2,2) * t(3,3) - t(2,3) * t(2,3)
   i3 = t(1,1) * (t(2,2) * t(3,3) - t(2,3) * t(2,3)) - t(1,2) * (t(1,2) * t(3,3) - t(2,3) * t(1,3)) + &
        t(1,3) * (t(1,2) * t(2,3) - t(2,2) * t(1,3))
   nu = 0._R8P
   a1 = i1 * i1 / 9._R8P - i2 / 3._R8P
   ! a1 = 0: every singular value equal (solid rotation or pure dilatation included), where s_1 = s_2 and the model vanishes
   if (.not.(a1 > 0._R8P)) return
   a2 = i1 * i1 * i1 / 27._R8P - i1 * i2 / 6._R8P + i3 / 2._R8P
   a3 = acos(max(-1._R8P, min(1._R8P, a2 / (a1 * sqrt(a1))))) / 3._R8P
   s1 = sqrt(max(0._R8P, i1 / 3._R8P + 2._R8P * sqrt(a1) * cos(a3)))
   s2 = sqrt(max(0._R8P, i1 / 3._R8P - 2._R8P * sqrt(a1) * cos(PI3 + a3)))
   s3 = sqrt(max(0._R8P, i1 / 3._R8P - 2._R8P * sqrt(a1) * cos(PI3 - a3)))
   ! s3^2 from the closed form is round-off when g is singular (2-D flow): s3 stands for 0 below the error of the expression
   if (s3 * s3 < 4._R8P * epsilon(1._R8P) * i1) s3 = 0._R8P
   if (s1 > 0._R8P) nu = max(0._R8P, (c * delta)**2 * s3 * (s1 - s2) * (s2 - s3) / (s1 * s1))
   endsubroutine sgs_sigma
endmodule adam_turbulence_library
