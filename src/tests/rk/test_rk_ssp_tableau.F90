!< Unit test of the SSP Runge-Kutta Butcher tableaux (`ssp_tableau`, issue #104).
program test_rk_ssp_tableau
!< Unit test of the SSP Runge-Kutta Butcher tableaux of `adam_rk_object` (`ssp_tableau`, issue #104).
!<
!< **Why this test exists.** The SSP(5,4) tableau shipped with 14-digit literals accurate to ~1e-10 only: its weights
!< summed to 1 - 8.8e-11 and all eight fourth-order conditions missed by ~1e-11, so a constant force drifted by that
!< relative amount and a refinement study stalled at ~3e-11. A tableau's accuracy is a property of its numbers, not of
!< the solver, so it is pinned here directly, with no grid, field or MPI.
!<
!< **What it pins**, for SSP-11, SSP-22, SSP-33 and SSP-54 (orders 1, 2, 3, 4):
!<   1. the order conditions up to the scheme's order (1, 2, 4 and 8 conditions), evaluated in quadruple precision on
!<      the stored doubles, so the residual is the tableau's and not the arithmetic's: <= 1e-14;
!<   2. the abscissas equal the row sums of A (<= 1e-15), and A is strictly lower triangular;
!<   3. a constant right-hand side y' = 1 integrated over 7 steps with the stage assembly of `rk_object`: exact to 1e-14;
!<   4. y' = -y to T = 1 on 40, 80, 160, 320 steps: the observed order of the finest pair within 0.1 of the scheme's,
!<      and for SSP-54 the error at 320 steps <= 1e-12 (the 1e-10-accurate tableau stalled at 3.2e-11).

use adam_rk_object, only : ssp_tableau, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54
use penf,           only : I4P, R8P, R16P

implicit none

character(len=18), parameter :: SCHEMES(4) = [RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54] !< SSP schemes.
integer(I4P),      parameter :: ORDERS(4)  = [1_I4P, 2_I4P, 3_I4P, 4_I4P]                 !< Their orders.
real(R8P),         parameter :: TOL_ORDER  = 1.e-14_R8P !< Bound on the order-condition residuals.
real(R8P),         parameter :: TOL_ROW    = 1.e-15_R8P !< Bound on |c_i - sum_j a_ij|.
real(R8P),         parameter :: TOL_CONST  = 1.e-14_R8P !< Bound on the constant right-hand side error, relative.
real(R8P),         parameter :: TOL_FLOOR  = 1.e-12_R8P !< Bound on the SSP-54 error at 320 steps.
real(R8P), allocatable       :: alph(:,:)               !< Butcher matrix A.
real(R8P), allocatable       :: beta(:)                 !< Butcher weights b.
real(R8P), allocatable       :: gamm(:)                 !< Butcher abscissas c.
integer(I4P)                 :: nrk                     !< Stages.
integer(I4P)                 :: m                       !< Scheme counter.
logical                      :: ok                      !< Overall result.

ok = .true.
do m=1, size(SCHEMES)
   nrk = 0_I4P
   call ssp_tableau(scheme=SCHEMES(m), nrk=nrk, alph=alph, beta=beta, gamm=gamm)
   print '(A)', 'scheme '//trim(SCHEMES(m))
   ok = check_structure(alph=alph, gamm=gamm) .and. ok
   ok = check_order_conditions(p=ORDERS(m), alph=alph, beta=beta, gamm=gamm) .and. ok
   ok = check_constant(alph=alph, beta=beta) .and. ok
   ok = check_convergence(p=ORDERS(m), alph=alph, beta=beta, floor=(ORDERS(m) == 4_I4P)) .and. ok
enddo
if (.not.ok) error stop 'test_rk_ssp_tableau: FAILED'
print '(A)', 'test_rk_ssp_tableau: PASSED'

contains
   function check_structure(alph, gamm) result(pass)
   !< A strictly lower triangular, c_i equal to the row sums of A.
   real(R8P), intent(in) :: alph(:,:) !< Butcher matrix A.
   real(R8P), intent(in) :: gamm(:)   !< Butcher abscissas c.
   logical               :: pass      !< Result.
   real(R8P)             :: worst     !< Largest |c_i - sum_j a_ij|.
   integer(I4P)          :: i, j      !< Counters.
   logical               :: lower     !< A strictly lower triangular.

   lower = .true.
   do j=1, size(alph, dim=2)
      do i=1, j
         lower = lower .and. alph(i,j) == 0._R8P
      enddo
   enddo
   worst = maxval([(abs(real(gamm(i), R16P) - sum(real(alph(i,:), R16P))), i=1, size(gamm))])
   pass = lower .and. worst <= TOL_ROW
   print '(A,L1,A,ES9.2,A)', '   '//merge('PASS', 'FAIL', pass)//' A strictly lower triangular: ', lower, &
                            ', max |c_i - sum_j a_ij| ', worst, ' (bound 1e-15)'
   endfunction check_structure

   function check_order_conditions(p, alph, beta, gamm) result(pass)
   !< The order conditions up to order p, in quadruple precision on the stored doubles.
   integer(I4P), intent(in) :: p         !< Order of the scheme.
   real(R8P),    intent(in) :: alph(:,:) !< Butcher matrix A.
   real(R8P),    intent(in) :: beta(:)   !< Butcher weights b.
   real(R8P),    intent(in) :: gamm(:)   !< Butcher abscissas c.
   logical                  :: pass      !< Result.
   real(R16P), allocatable  :: a(:,:)    !< A in quadruple precision.
   real(R16P), allocatable  :: b(:)      !< b in quadruple precision.
   real(R16P), allocatable  :: c(:)      !< c in quadruple precision.
   real(R16P), allocatable  :: ac(:)     !< A c.
   real(R16P)               :: res(8)    !< Residuals.
   character(len=14)        :: names(8)  !< Condition names.
   integer(I4P)             :: n         !< Conditions of order <= p.
   integer(I4P)             :: k         !< Counter.
   integer(I4P), parameter  :: CONDITIONS(4) = [1_I4P, 2_I4P, 4_I4P, 8_I4P] !< Conditions up to order 1..4.

   a = real(alph, R16P) ; b = real(beta, R16P) ; c = real(gamm, R16P) ; ac = matmul(a, c)
   names = [character(len=14) :: 'sum b = 1', 'b.c = 1/2', 'b.c^2 = 1/3', 'b.Ac = 1/6', 'b.c^3 = 1/4', &
                                 'b.(c*Ac) = 1/8', 'b.Ac^2 = 1/12', 'b.AAc = 1/24']
   res(1) = sum(b) - 1._R16P
   res(2) = dot_product(b, c) - 1._R16P / 2._R16P
   res(3) = dot_product(b, c**2) - 1._R16P / 3._R16P
   res(4) = dot_product(b, ac) - 1._R16P / 6._R16P
   res(5) = dot_product(b, c**3) - 1._R16P / 4._R16P
   res(6) = dot_product(b, c * ac) - 1._R16P / 8._R16P
   res(7) = dot_product(b, matmul(a, c**2)) - 1._R16P / 12._R16P
   res(8) = dot_product(b, matmul(a, ac)) - 1._R16P / 24._R16P
   n = CONDITIONS(p)
   pass = .true.
   do k=1, n
      pass = pass .and. abs(res(k)) <= TOL_ORDER
      print '(A,ES10.2)', '   '//merge('PASS', 'FAIL', abs(res(k)) <= TOL_ORDER)//' '//names(k)//' residual', &
                          real(res(k), R8P)
   enddo
   endfunction check_order_conditions

   function check_constant(alph, beta) result(pass)
   !< y' = 1, y(0) = 0, 7 steps of 0.1 by the stage assembly of `rk_object`: y(0.7) = 0.7 to round-off.
   real(R8P), intent(in) :: alph(:,:) !< Butcher matrix A.
   real(R8P), intent(in) :: beta(:)   !< Butcher weights b.
   logical               :: pass      !< Result.
   real(R8P)             :: y         !< Solution.
   real(R8P)             :: err       !< Relative error.
   integer(I4P)          :: n         !< Step counter.

   y = 0._R8P
   do n=1, 7
      y = rk_step(alph=alph, beta=beta, lambda=0._R8P, forcing=1._R8P, dt=0.1_R8P, y=y)
   enddo
   err = abs(y - 0.7_R8P) / 0.7_R8P
   pass = err <= TOL_CONST
   print '(A,ES10.2,A)', '   '//merge('PASS', 'FAIL', pass)//" y' = 1 over 7 steps, relative error", err, ' (bound 1e-14)'
   endfunction check_constant

   function check_convergence(p, alph, beta, floor) result(pass)
   !< y' = -y, y(0) = 1, T = 1 on 40, 80, 160, 320 steps: the order of the finest pair within 0.1 of p; with `floor`,
   !< the error at 320 steps <= TOL_FLOOR.
   integer(I4P), intent(in) :: p         !< Order of the scheme.
   real(R8P),    intent(in) :: alph(:,:) !< Butcher matrix A.
   real(R8P),    intent(in) :: beta(:)   !< Butcher weights b.
   logical,      intent(in) :: floor     !< Assert the error floor.
   logical                  :: pass      !< Result.
   integer(I4P), parameter  :: STEPS(4) = [40_I4P, 80_I4P, 160_I4P, 320_I4P] !< Step counts.
   real(R8P)                :: err(4)    !< Errors.
   real(R8P)                :: y         !< Solution.
   real(R8P)                :: order     !< Order of the finest pair.
   integer(I4P)             :: k, n      !< Counters.

   do k=1, size(STEPS)
      y = 1._R8P
      do n=1, STEPS(k)
         y = rk_step(alph=alph, beta=beta, lambda=-1._R8P, forcing=0._R8P, dt=1._R8P / STEPS(k), y=y)
      enddo
      err(k) = abs(y - exp(-1._R8P))
   enddo
   order = log(err(3) / err(4)) / log(2._R8P)
   pass = order >= real(p, R8P) - 0.1_R8P
   print '(A,4ES10.2,A,F6.2)', '   '//merge('PASS', 'FAIL', pass)//" y' = -y, errors on 40..320 steps", err, &
                               ', finest-pair order', order
   if (floor) then
      pass = pass .and. err(4) <= TOL_FLOOR
      print '(A,ES10.2,A)', '   '//merge('PASS', 'FAIL', err(4) <= TOL_FLOOR)//' error at 320 steps', err(4), &
                            ' (bound 1e-12)'
   endif
   endfunction check_convergence

   pure function rk_step(alph, beta, lambda, forcing, dt, y) result(y1)
   !< One step of y' = lambda y + forcing, assembled as `rk_object` does: stage s from the stored stage derivatives
   !< k_ss, ss < s, weighted by alph(s, ss); the update by beta.
   real(R8P), intent(in) :: alph(:,:)       !< Butcher matrix A.
   real(R8P), intent(in) :: beta(:)         !< Butcher weights b.
   real(R8P), intent(in) :: lambda          !< Linear coefficient.
   real(R8P), intent(in) :: forcing         !< Constant forcing.
   real(R8P), intent(in) :: dt              !< Step.
   real(R8P), intent(in) :: y               !< Solution at the start of the step.
   real(R8P)             :: y1              !< Solution at the end of the step.
   real(R8P)             :: k(size(beta))   !< Stage derivatives.
   real(R8P)             :: ys              !< Stage value.
   integer(I4P)          :: s, ss           !< Counters.

   do s=1, size(beta)
      ys = y
      do ss=1, s - 1
         ys = ys + dt * alph(s,ss) * k(ss)
      enddo
      k(s) = lambda * ys + forcing
   enddo
   y1 = y
   do s=1, size(beta)
      y1 = y1 + dt * beta(s) * k(s)
   enddo
   endfunction rk_step
endprogram test_rk_ssp_tableau
