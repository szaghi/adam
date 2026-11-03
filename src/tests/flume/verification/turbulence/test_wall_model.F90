!< The equilibrium wall model of FLUME (issue #119, T3.4), pointwise: leg wm2 of check.sh.

program test_wall_model
!< Laminar limit, log layer against a fine quadrature, and the Crocco relation with Pr = Pr_t.

use :: adam_flume_wall_model_library, only : sample_cells, wall_model_solve
use :: penf,                          only : I4P, I8P, R8P

implicit none
real(R8P), parameter :: CP=3.5_R8P, GAS=1._R8P, KAPPA=0.41_R8P, APLUS=17._R8P
real(R8P)            :: tau, q, us, hplus, err, ref, mu, dy, y1, y2, f1, f2, tauf, hh, uu, tt, pp
real(R8P), parameter :: tw15=1._R8P / (1.4_R8P * 2.25_R8P)
integer(I8P)         :: seed
integer(I4P)         :: nbad, mm
real(R8P)            :: ww
logical              :: ok1
integer(I4P)         :: n
integer(I4P), parameter :: NF=4000000_I4P
logical              :: ok

ok = .true.
! laminar limit: h+ small, so mu_t vanishes: tau = mu U / h, q = cp mu / Pr (Ts - Tw) / h with Couette work
call wall_model_solve(h=1._R8P, us=1.e-3_R8P, ts=1._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1._R8P, tref=1._R8P, omega_mu=0._R8P, tau=tau, q=q)
err = abs(tau - 1.e-3_R8P) / 1.e-3_R8P
print '(A,ES10.2,A)', 'laminar limit, tau relative ', err, ' (bound 1e-6)'
ok = ok .and. err < 1.e-6_R8P
! the log layer: rho = u_tau = 1, mu = 1 / h+, at a temperature so high that the friction heating is negligible (Mach 1e-2): U(h) of the constant-stress layer
   ! by a fine quadrature, then the solver must return tau = 1
do n=1, 4
   hplus = 10._R8P**(n + 1)
   mu = 1._R8P / hplus
   us = 0._R8P
   dy = 1._R8P / real(NF, R8P)
   f1 = 1._R8P / mu
   ! midpoint-free trapezoid on a quadratic map y = s^2, to resolve the wall
   block_i: block
      integer(I4P) :: i
      do i=1, NF
         y1 = (real(i-1, R8P) * dy)**2 ; y2 = (real(i, R8P) * dy)**2
         f2 = 1._R8P / (mu + KAPPA * y2 * (1._R8P - exp(-y2 * hplus / APLUS))**2)
         us = us + 0.5_R8P * (y2 - y1) * (f1 + f2)
         f1 = f2
      enddo
   endblock block_i
   ref = us
   call wall_model_solve(h=1._R8P, us=ref, ts=1.e6_R8P, ps=1.e6_R8P, tw=1.e6_R8P, R=GAS, cp=CP, pr=0.9_R8P, mu0=mu, tref=1._R8P, omega_mu=0._R8P, tau=tau, q=q)
   err = abs(tau - 1._R8P)
   print '(A,ES9.1,A,F8.3,A,ES10.2,A)', 'h+ ', hplus, ', U+ ', ref, ', tau relative ', err, ' (bound 1e-2)'
   ok = ok .and. err < 1.e-2_R8P
enddo
! Crocco with Pr = Pr_t = 0.9: cp (T - Tw) / Pr = Q U / tau - U^2 / 2, so Q = tau (cp (Ts - Tw) / Pr + Us^2 / 2) / Us, with a viscosity law too
call wall_model_solve(h=1._R8P, us=1._R8P, ts=0.6_R8P, ps=1._R8P, tw=0.5_R8P, R=GAS, cp=CP, pr=0.9_R8P, mu0=1.e-4_R8P, tref=0.5_R8P, omega_mu=0.7_R8P, tau=tau, q=q)
ref = tau * (CP * 0.1_R8P / 0.9_R8P + 0.5_R8P)
err = abs(q - ref) / abs(ref)
print '(A,ES10.2,A,ES10.2,A)', 'Crocco, tau ', tau, ', q relative ', err, ' (bound 5e-3)'
ok = ok .and. err < 5.e-3_R8P
! compressible: a hot layer, Pr 0.7, a viscosity law, against the same equations on a fine grid with many iterations: the local density and the heat flux
call reference(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, tau=tauf, q=ref)
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
print '(A,2ES10.2,A)', 'compressible, tau and q relative ', abs(tau - tauf) / tauf, abs(q - ref) / abs(ref), ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P .and. abs(q - ref) / abs(ref) < 2.e-2_R8P
! semi-local damping, the same hot layer, against the reference with the same damping
call reference(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, tau=tauf, q=ref, semilocal=.true.)
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega_mu=0.7_R8P, tau=tau, q=q, &
                      semilocal=.true.)
print '(A,2ES10.2,A)', 'semi-local, tau and q relative ', abs(tau - tauf) / tauf, abs(q - ref) / abs(ref), ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P .and. abs(q - ref) / abs(ref) < 2.e-2_R8P
ref = tau
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
print '(A,ES10.2,A)', 'semi-local against classical on the hot layer, tau relative ', abs(ref - tau) / tau, ' (min 1e-2: the option acts)'
ok = ok .and. abs(ref - tau) / tau > 1.e-2_R8P
! the sample cells: the 3rd cell by default; a height between cells 3 and 4; a height beyond the block, its last cells
call sample_cells(hwm=0._R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = mm == 3 .and. ww == 0._R8P .and. abs(hh - 0.25_R8P) < 1.e-15_R8P
call sample_cells(hwm=0.3_R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = ok1 .and. mm == 3 .and. abs(ww - 0.5_R8P) < 1.e-12_R8P .and. abs(hh - 0.3_R8P) < 1.e-12_R8P
call sample_cells(hwm=2._R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = ok1 .and. mm == 7 .and. ww == 1._R8P .and. abs(hh - 0.75_R8P) < 1.e-12_R8P
print '(A,L1)', 'sample cells: ', ok1
ok = ok .and. ok1
! robustness: 200000 pseudo-random samples about a Mach 1.5 channel wall, parallel speed up to 4.5 wall sound speeds, sample 0.5-3.5 times the wall temperature
nbad = 0 ; seed = 12345_I8P
do n=1, 200000
   hh = 0.01_R8P + 0.3_R8P * lcg(seed) ; uu = 3._R8P * lcg(seed) ; tt = tw15 * (0.5_R8P + 3._R8P * lcg(seed)) ; pp = (0.5_R8P + lcg(seed)) * tw15
   call wall_model_solve(h=hh, us=uu, ts=tt, ps=pp, tw=tw15, R=GAS, cp=CP, pr=0.7_R8P, mu0=1._R8P / 3000._R8P, tref=tw15, omega_mu=0.7_R8P, tau=tau, q=q, &
                         semilocal=mod(n, 2) == 0)
   if (.not.(tau >= 0._R8P .and. tau < huge(1._R8P) .and. abs(q) < huge(1._R8P))) nbad = nbad + 1
enddo
print '(A,I0,A)', 'robustness, non-finite or negative results: ', nbad, ' of 200000 (bound 0)'
ok = ok .and. nbad == 0
if (ok) then
   print '(A)', 'test_wall_model: PASS'
else
   print '(A)', 'test_wall_model: FAIL'
   stop 1
endif
contains
   function lcg(s) result(x)
   !< A linear congruential generator in (0, 1): the same samples everywhere.
   integer(I8P), intent(inout) :: s !< State.
   real(R8P)                   :: x !< Sample.

   s = modulo(s * 48271_I8P, 2147483647_I8P) ! MINSTD: no overflow in 64 bits
   x = real(s, R8P) / 2147483647._R8P
   endfunction lcg

   subroutine reference(h, us, ts, ps, tw, pr, mu0, tref, omega, tau, q, semilocal)
   !< The layer on 200000 intervals uniform in sqrt(y), by trapezoids and 400 under-relaxed iterations: no code shared with the library.
   real(R8P), intent(in)  :: h, us, ts, ps, tw, pr, mu0, tref, omega !< As `wall_model_solve`.
   real(R8P), intent(out) :: tau, q                                   !< Results.
   logical,   intent(in), optional :: semilocal                       !< Semi-local damping.
   real(R8P)              :: dmp
   integer(I4P), parameter :: M=200000_I4P
   real(R8P), allocatable :: y(:), t(:), u(:), mt(:), a(:)
   real(R8P)              :: rhow, muw, d, i1, i2, i3
   integer(I4P)           :: i, it

   allocate(y(0:M), t(0:M), u(0:M), mt(0:M), a(0:M))
   do i=0, M
      y(i) = h * (real(i, R8P) / M)**2
      t(i) = tw + (ts - tw) * y(i) / h
   enddo
   rhow = ps / (GAS * tw) ; muw = mu0 * (tw / tref)**omega
   tau = muw * us / h
   do it=1, 400
      d = 0._R8P
      do i=0, M
         dmp = 1._R8P - exp(-y(i) * sqrt(tau * rhow) / (muw * APLUS))
         if (present(semilocal)) then
            if (semilocal) dmp = 1._R8P - exp(-y(i) * sqrt(tau * ps / (GAS * t(i))) / (mu0 * (t(i) / tref)**omega * APLUS))
         endif
         mt(i) = KAPPA * y(i) * sqrt(ps / (GAS * t(i)) * tau) * dmp**2
         a(i) = mu0 * (t(i) / tref)**omega / pr + mt(i) / 0.9_R8P
         mt(i) = mu0 * (t(i) / tref)**omega + mt(i)
      enddo
      i1 = 0._R8P
      do i=1, M
         i1 = i1 + 0.5_R8P * (y(i) - y(i-1)) * (1._R8P / mt(i) + 1._R8P / mt(i-1))
      enddo
      tau = sqrt(tau * us / i1)
      u(0) = 0._R8P ; i2 = 0._R8P ; i3 = 0._R8P
      do i=1, M
         u(i) = u(i-1) + 0.5_R8P * (y(i) - y(i-1)) * tau * (1._R8P / mt(i) + 1._R8P / mt(i-1))
         i2 = i2 + 0.5_R8P * (y(i) - y(i-1)) * (1._R8P / a(i) + 1._R8P / a(i-1))
         i3 = i3 + 0.5_R8P * (y(i) - y(i-1)) * (u(i) / a(i) + u(i-1) / a(i-1))
      enddo
      q = (CP * (ts - tw) + tau * i3) / i2
      i2 = 0._R8P ; i3 = 0._R8P
      do i=1, M
         i2 = i2 + 0.5_R8P * (y(i) - y(i-1)) * (1._R8P / a(i) + 1._R8P / a(i-1))
         i3 = i3 + 0.5_R8P * (y(i) - y(i-1)) * (u(i) / a(i) + u(i-1) / a(i-1))
         t(i) = 0.5_R8P * t(i) + 0.5_R8P * (tw + (q * i2 - tau * i3) / CP)
      enddo
   enddo
   endsubroutine reference
endprogram test_wall_model
