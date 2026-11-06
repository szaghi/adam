!< The wall models of FLUME (issue #119, T3.4, T3.5), pointwise: leg wm2 of check.sh.

program test_wall_model
!< Laminar limit, log layer against a fine quadrature, and the Crocco relation with Pr = Pr_t; the models under Sutherland's
!< law against the same references, and the law itself.

use :: adam_flume_dissipation_library, only : temperature_law
use :: adam_flume_wall_model_library, only : inverse_strain, sample_cells, wall_model_solve, wall_model_solve_inverse
use :: penf,                          only : I4P, I8P, R8P

implicit none
real(R8P), parameter :: CP=3.5_R8P, GAS=1._R8P, KAPPA=0.41_R8P, APLUS=17._R8P
real(R8P)            :: tau, q, us, hplus, err, ref, mu, dy, y1, y2, f1, f2, tauf, hh, uu, tt, pp
real(R8P), parameter :: tw15=1._R8P / (1.4_R8P * 2.25_R8P)
integer(I8P)         :: seed
integer(I4P)         :: nbad, mm
real(R8P)            :: ww, te0, tauc, qc, c1, g, gl, gs
logical              :: ok1
integer(I4P)         :: n
integer(I4P), parameter :: NF=4000000_I4P
logical              :: ok

ok = .true.
! laminar limit: h+ small, so mu_t vanishes: tau = mu U / h, q = cp mu / Pr (Ts - Tw) / h with Couette work
call wall_model_solve(h=1._R8P, us=1.e-3_R8P, ts=1._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1._R8P, tref=1._R8P, suth=0._R8P, omega_mu=0._R8P, tau=tau, q=q)
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
   call wall_model_solve(h=1._R8P, us=ref, ts=1.e6_R8P, ps=1.e6_R8P, tw=1.e6_R8P, R=GAS, cp=CP, pr=0.9_R8P, mu0=mu, tref=1._R8P, suth=0._R8P, omega_mu=0._R8P, tau=tau, q=q)
   err = abs(tau - 1._R8P)
   print '(A,ES9.1,A,F8.3,A,ES10.2,A)', 'h+ ', hplus, ', U+ ', ref, ', tau relative ', err, ' (bound 1e-2)'
   ok = ok .and. err < 1.e-2_R8P
   ! the inverse model in the same layer, the edge state making Tr = Tw (so T = Tw throughout, uniform properties): the same equation on the same grid as the
   ! equilibrium model, only another iteration (secant shots against Aitken), so the two agree to their convergence
   ! (both at T = 1e10: the friction heating of the equilibrium model, U^2 / cp about 3e-4 of T = 1e6, would vary its density)
   call wall_model_solve(h=1._R8P, us=ref, ts=1.e10_R8P, ps=1.e10_R8P, tw=1.e10_R8P, R=GAS, cp=CP, pr=0.9_R8P, mu0=mu, tref=1._R8P, suth=0._R8P, omega_mu=0._R8P, tau=tauc, q=qc)
   te0 = 1.e10_R8P - 0.9_R8P**(1._R8P / 3._R8P) * ref * ref / (2._R8P * CP)
   call wall_model_solve_inverse(h=1._R8P, us=ref, ts=1.e10_R8P, ps=1.e10_R8P, tw=1.e10_R8P, ue=ref, te=te0, R=GAS, cp=CP, pr=0.9_R8P, mu0=mu, tref=1._R8P, suth=0._R8P, &
                                 omega_mu=0._R8P, tau=tau, q=q)
   print '(A,ES9.1,A,ES10.2,A,ES10.2,A)', '   inverse model: tau relative ', abs(tau - 1._R8P), ' (bound 1e-2), against the equilibrium model ', &
         abs(tau - tauc) / tauc, ' (bound 1e-6), at h+ ', hplus
   ok = ok .and. abs(tau - 1._R8P) < 1.e-2_R8P .and. abs(tau - tauc) / tauc < 1.e-6_R8P
enddo
! Crocco with Pr = Pr_t = 0.9: cp (T - Tw) / Pr = Q U / tau - U^2 / 2, so Q = tau (cp (Ts - Tw) / Pr + Us^2 / 2) / Us, with a viscosity law too
call wall_model_solve(h=1._R8P, us=1._R8P, ts=0.6_R8P, ps=1._R8P, tw=0.5_R8P, R=GAS, cp=CP, pr=0.9_R8P, mu0=1.e-4_R8P, tref=0.5_R8P, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
ref = tau * (CP * 0.1_R8P / 0.9_R8P + 0.5_R8P)
err = abs(q - ref) / abs(ref)
print '(A,ES10.2,A,ES10.2,A)', 'Crocco, tau ', tau, ', q relative ', err, ' (bound 5e-3)'
ok = ok .and. err < 5.e-3_R8P
! compressible: a hot layer, Pr 0.7, a viscosity law, against the same equations on a fine grid with many iterations: the local density and the heat flux
call reference(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, tau=tauf, q=ref)
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
print '(A,2ES10.2,A)', 'compressible, tau and q relative ', abs(tau - tauf) / tauf, abs(q - ref) / abs(ref), ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P .and. abs(q - ref) / abs(ref) < 2.e-2_R8P
! semi-local damping, the same hot layer, against the reference with the same damping
call reference(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, tau=tauf, q=ref, semilocal=.true.)
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q, &
                      semilocal=.true.)
print '(A,2ES10.2,A)', 'semi-local, tau and q relative ', abs(tau - tauf) / tauf, abs(q - ref) / abs(ref), ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P .and. abs(q - ref) / abs(ref) < 2.e-2_R8P
ref = tau
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
print '(A,ES10.2,A)', 'semi-local against classical on the hot layer, tau relative ', abs(ref - tau) / tau, ' (min 1e-2: the option acts)'
ok = ok .and. abs(ref - tau) / tau > 1.e-2_R8P
! the inverse model (T3.5) on a cold wall, the gas twice as hot, a viscosity law, against Eq. 2.2 of Griffin et al. (2023) integrated as printed on a fine grid
call reference_inverse(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, ue=1.2_R8P, te=2.1_R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, &
                       tau=tauf)
call wall_model_solve_inverse(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, ue=1.2_R8P, te=2.1_R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, &
                              omega_mu=0.7_R8P, tau=tau, q=q)
print '(A,ES10.2,A)', 'inverse model, cold wall, tau relative to Eq. 2.2 on a fine grid ', abs(tau - tauf) / tauf, ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P
! and on a hot wall, the gas half as hot (beta < 0: the roots of the strain have opposite signs)
call reference_inverse(h=0.3_R8P, us=1._R8P, ts=1._R8P, ps=1._R8P, tw=2._R8P, ue=1.2_R8P, te=0.9_R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=0.7_R8P, &
                       tau=tauf)
call wall_model_solve_inverse(h=0.3_R8P, us=1._R8P, ts=1._R8P, ps=1._R8P, tw=2._R8P, ue=1.2_R8P, te=0.9_R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, &
                              omega_mu=0.7_R8P, tau=tauc, q=qc)
print '(A,ES10.2,A)', 'inverse model, hot wall, tau relative to Eq. 2.2 on a fine grid ', abs(tauc - tauf) / tauf, ' (bound 2e-2)'
ok = ok .and. abs(tauc - tauf) / tauf < 2.e-2_R8P
! its heat flux is the Reynolds analogy, q = s cp tau (Tr - Tw) / Ue (the cold wall)
ref = 1.14_R8P * CP * tau * (2.1_R8P + 0.7_R8P**(1._R8P / 3._R8P) * 1.44_R8P / (2._R8P * CP) - 1._R8P) / 1.2_R8P
print '(A,ES10.2,A)', 'inverse model, Reynolds analogy, q relative ', abs(q - ref) / abs(ref), ' (bound 1e-12)'
ok = ok .and. abs(q - ref) / abs(ref) < 1.e-12_R8P
! and it is not the equilibrium model there
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0._R8P, omega_mu=0.7_R8P, tau=tauc, q=qc)
print '(A,ES10.2,A)', 'inverse against equilibrium model on the cold wall, tau relative ', abs(tau - tauc) / tauc, ' (min 1e-2: the option acts)'
ok = ok .and. abs(tau - tauc) / tauc > 1.e-2_R8P
! the root of the strain: in a laminar layer (y* 0.04, so S_t = 1 to 1e-6) cooled so hard that beta y g > 1, the strain is tau / mu, not the spurious 1 / (beta y)
c1 = 400._R8P ! dT/dU at the wall
call strain_point(g, gl, gs)
print '(A,ES10.2,A,ES10.2,A)', 'inverse strain beyond dy*/dy = 0: relative to tau / mu ', abs(g - gl) / gl, ' (bound 2e-3), the other root ', abs(gs - gl) / gl, &
      ' (min 0.1)'
ok = ok .and. abs(g - gl) / gl < 2.e-3_R8P .and. abs(gs - gl) / gl > 0.1_R8P
! the sample cells: the 3rd cell by default; a height between cells 3 and 4; a height beyond the block, its last cells
call sample_cells(hwm=0._R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = mm == 3 .and. ww == 0._R8P .and. abs(hh - 0.25_R8P) < 1.e-15_R8P
call sample_cells(hwm=0.3_R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = ok1 .and. mm == 3 .and. abs(ww - 0.5_R8P) < 1.e-12_R8P .and. abs(hh - 0.3_R8P) < 1.e-12_R8P
call sample_cells(hwm=2._R8P, dx=0.1_R8P, n=8_I4P, m=mm, w=ww, h=hh)
ok1 = ok1 .and. mm == 7 .and. ww == 1._R8P .and. abs(hh - 0.75_R8P) < 1.e-12_R8P
print '(A,L1)', 'sample cells: ', ok1
ok = ok .and. ok1
! Sutherland's law: the library's temperature_law against the formula; the power law and a zero exponent unchanged bit for bit
err = abs(temperature_law(t=2.3_R8P, tref=1.1_R8P, suth=0.4_R8P, omega=1.5_R8P) - sqrt(2.3_R8P / 1.1_R8P)**3 * 1.5_R8P / 2.7_R8P)
ok1 = temperature_law(t=2.3_R8P, tref=1.1_R8P, suth=0._R8P, omega=0.7_R8P) == (2.3_R8P / 1.1_R8P)**0.7_R8P .and. &
      temperature_law(t=2.3_R8P, tref=1.1_R8P, suth=0.4_R8P, omega=0._R8P) == 1._R8P
print '(A,ES10.2,A,L1,A)', 'Sutherland law against the formula ', err, ' (bound 1e-14); the power law and a zero exponent bitwise: ', ok1
ok = ok .and. err < 1.e-14_R8P .and. ok1
! the equilibrium model under Sutherland (S = 0.4 T_ref, air at 276 K), the hot layer, against the reference with the same law
call reference(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=1.5_R8P, tau=tauf, q=ref, suth=0.4_R8P)
call wall_model_solve(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, suth=0.4_R8P, omega_mu=1.5_R8P, &
                      tau=tau, q=q)
print '(A,2ES10.2,A)', 'Sutherland, equilibrium, tau and q relative ', abs(tau - tauf) / tauf, abs(q - ref) / abs(ref), ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P .and. abs(q - ref) / abs(ref) < 2.e-2_R8P
! the inverse model under Sutherland on the cold wall: the property gradients of the reference are differences, the library's the analytic slope 3/2 - T / (T + S)
call reference_inverse(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, ue=1.2_R8P, te=2.1_R8P, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, omega=1.5_R8P, &
                       tau=tauf, suth=0.4_R8P)
call wall_model_solve_inverse(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, ue=1.2_R8P, te=2.1_R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, &
                              suth=0.4_R8P, omega_mu=1.5_R8P, tau=tau, q=q)
print '(A,ES10.2,A)', 'Sutherland, inverse model, cold wall, tau relative to Eq. 2.2 on a fine grid ', abs(tau - tauf) / tauf, ' (bound 2e-2)'
ok = ok .and. abs(tau - tauf) / tauf < 2.e-2_R8P
! and the law acts: the power law of the same exponent gives another stress
call wall_model_solve_inverse(h=0.3_R8P, us=1._R8P, ts=2._R8P, ps=1._R8P, tw=1._R8P, ue=1.2_R8P, te=2.1_R8P, R=GAS, cp=CP, pr=0.7_R8P, mu0=1.e-4_R8P, tref=1._R8P, &
                              suth=0._R8P, omega_mu=1.5_R8P, tau=tauc, q=qc)
print '(A,ES10.2,A)', 'Sutherland against the power law 3/2, inverse model, tau relative ', abs(tau - tauc) / tauc, ' (min 1e-2: the law acts)'
ok = ok .and. abs(tau - tauc) / tauc > 1.e-2_R8P
! robustness: 200000 pseudo-random samples about a Mach 1.5 channel wall, parallel speed up to 4.5 wall sound speeds, sample 0.5-3.5 times the wall temperature
nbad = 0 ; seed = 12345_I8P
do n=1, 200000
   hh = 0.01_R8P + 0.3_R8P * lcg(seed) ; uu = 3._R8P * lcg(seed) ; tt = tw15 * (0.5_R8P + 3._R8P * lcg(seed)) ; pp = (0.5_R8P + lcg(seed)) * tw15
   call wall_model_solve(h=hh, us=uu, ts=tt, ps=pp, tw=tw15, R=GAS, cp=CP, pr=0.7_R8P, mu0=1._R8P / 3000._R8P, tref=tw15, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q, &
                         semilocal=mod(n, 2) == 0)
   if (.not.(tau >= 0._R8P .and. tau < huge(1._R8P) .and. abs(q) < huge(1._R8P))) nbad = nbad + 1
   ! the inverse model on the same sample, an edge 1 to 1.5 times faster and 0.8 to 1.3 times as hot
   call wall_model_solve_inverse(h=hh, us=uu, ts=tt, ps=pp, tw=tw15, ue=uu * (1._R8P + 0.5_R8P * lcg(seed)) + 1.e-3_R8P, te=tt * (0.8_R8P + 0.5_R8P * lcg(seed)), &
                                 R=GAS, cp=CP, pr=0.7_R8P, mu0=1._R8P / 3000._R8P, tref=tw15, suth=0._R8P, omega_mu=0.7_R8P, tau=tau, q=q)
   if (.not.(tau >= 0._R8P .and. tau < huge(1._R8P) .and. abs(q) < huge(1._R8P))) nbad = nbad + 1
enddo
print '(A,I0,A)', 'robustness, non-finite or negative results: ', nbad, ' of 400000 (bound 0)'
ok = ok .and. nbad == 0
if (ok) then
   print '(A)', 'test_wall_model: PASS'
else
   print '(A)', 'test_wall_model: FAIL'
   stop 1
endif
contains
   pure function law(t, tref, omega, s) result(r)
   !< The viscosity ratio of the test, written apart from the library: the power law, or Sutherland's for s > 0.
   real(R8P), intent(in) :: t, tref, omega, s
   real(R8P)             :: r

   if (s > 0._R8P) then
      r = sqrt(t / tref)**3 * (tref + s) / (t + s)
   else
      r = (t / tref)**omega
   endif
   endfunction law

   function lcg(s) result(x)
   !< A linear congruential generator in (0, 1): the same samples everywhere.
   integer(I8P), intent(inout) :: s !< State.
   real(R8P)                   :: x !< Sample.

   s = modulo(s * 48271_I8P, 2147483647_I8P) ! MINSTD: no overflow in 64 bits
   x = real(s, R8P) / 2147483647._R8P
   endfunction lcg

   subroutine strain_point(g, gl, gs)
   !< The strain of the inverse model at one point of a laminar layer (y = tau = mu_w = 1, y* about 0.04) where beta y g > 1: `g` the model's, `gl` the
   !< laminar tau / mu, `gs` the root 1 / (beta y) of dy*/dy = 0; u = 0.02 of us = 1, Tw = 1, Ts = 1, omega 1.
   real(R8P), intent(out) :: g, gl, gs !< Strains.
   real(R8P)              :: t, dtdu  !< Temperature, its slope.

   g = inverse_strain(y=1._R8P, u=0.02_R8P, tau=1._R8P, us=1._R8P, ts=1._R8P, ps=1._R8P, tw=1._R8P, c1=c1, R=GAS, mu0=1._R8P, tref=1._R8P, suth=0._R8P, omega_mu=1._R8P)
   t = 1._R8P + c1 * 0.02_R8P * 0.98_R8P
   dtdu = c1 * (1._R8P - 0.04_R8P)
   gl = 1._R8P / t        ! tau / mu, mu = T
   gs = t / (1.5_R8P * dtdu) ! 1 / (beta y)
   endsubroutine strain_point

   subroutine reference_inverse(h, us, ts, ps, tw, ue, te, pr, mu0, tref, omega, tau, suth)
   !< The inverse model by Eq. 2.2 of Griffin et al. (2023) as printed, in y+ and y*, the property gradients by differences of the previous points, explicit steps on
   !< 200000 intervals uniform in sqrt(y), bisection on ln(tau): no code shared with the library.
   real(R8P), intent(in)  :: h, us, ts, ps, tw, ue, te, pr, mu0, tref, omega !< As `wall_model_solve_inverse`.
   real(R8P), intent(out) :: tau                                              !< Wall stress.
   real(R8P), intent(in), optional :: suth                                   !< Sutherland constant (default 0, the power law).
   integer(I4P), parameter :: M=200000_I4P
   real(R8P)              :: rhow, muw, cc, lo, hi, x, u, t, rp, mp, yp, ys, rp0, mp0, yp0, drp, dmp, st, inv, dys, ut, uh, y0, y1, s
   integer(I4P)           :: i, it

   s = 0._R8P ; if (present(suth)) s = suth
   rhow = ps / (GAS * tw) ; muw = mu0 * law(tw, tref, omega, s)
   cc = 1.14_R8P * pr * (te + pr**(1._R8P / 3._R8P) * ue * ue / (2._R8P * CP) - tw) / ue
   lo = log(muw * us / h) - 1._R8P ; hi = lo + 20._R8P
   do it=1, 80
      x = 0.5_R8P * (lo + hi)
      tau = exp(x) ; ut = sqrt(tau / rhow)
      u = 0._R8P ; rp0 = 1._R8P ; mp0 = 1._R8P ; yp0 = 0._R8P ; drp = 0._R8P ; dmp = 0._R8P
      do i=1, M
         y0 = h * (real(i - 1, R8P) / M)**2 ; y1 = h * (real(i, R8P) / M)**2
         t = tw + cc * u * (1._R8P - u / us) + (u / us)**2 * (ts - tw)
         rp = tw / t ; mp = law(t, tref, omega, s) / law(tw, tref, omega, s) ; yp = y0 * sqrt(tau * rhow) / muw ; ys = yp * sqrt(rp) / mp
         if (i > 2) then
            drp = (rp - rp0) / (yp - yp0) ; dmp = (mp - mp0) / (yp - yp0)
         endif
         st = 1._R8P / (1._R8P + KAPPA * ys * (1._R8P - exp(-ys / APLUS))**2)
         dys = sqrt(rp) * (1._R8P + yp / (2._R8P * rp) * drp - yp / mp * dmp) / mp                   ! dy*/dy+
         inv = 1._R8P / (mp * st) - 1._R8P / mp + sqrt(rp) * (1._R8P + yp / (2._R8P * rp) * drp - yp / mp * dmp) ! (dU+/dy*)^-1, Eq. 2.2
         u = u + ut * dys / inv * (y1 - y0) * sqrt(tau * rhow) / muw
         rp0 = rp ; mp0 = mp ; yp0 = yp
         if (u > us) exit ! tau too large already (and T(U) is a quadratic fitted on [0, us] only)
      enddo
      uh = u
      if (uh > us) then
         hi = x
      else
         lo = x
      endif
   enddo
   tau = exp(0.5_R8P * (lo + hi))
   endsubroutine reference_inverse

   subroutine reference(h, us, ts, ps, tw, pr, mu0, tref, omega, tau, q, semilocal, suth)
   !< The layer on 200000 intervals uniform in sqrt(y), by trapezoids and 400 under-relaxed iterations: no code shared with the library.
   real(R8P), intent(in)  :: h, us, ts, ps, tw, pr, mu0, tref, omega !< As `wall_model_solve`.
   real(R8P), intent(out) :: tau, q                                   !< Results.
   logical,   intent(in), optional :: semilocal                       !< Semi-local damping.
   real(R8P), intent(in), optional :: suth                            !< Sutherland constant (default 0, the power law).
   real(R8P)              :: dmp
   integer(I4P), parameter :: M=200000_I4P
   real(R8P), allocatable :: y(:), t(:), u(:), mt(:), a(:)
   real(R8P)              :: rhow, muw, d, i1, i2, i3, s
   integer(I4P)           :: i, it

   allocate(y(0:M), t(0:M), u(0:M), mt(0:M), a(0:M))
   do i=0, M
      y(i) = h * (real(i, R8P) / M)**2
      t(i) = tw + (ts - tw) * y(i) / h
   enddo
   s = 0._R8P ; if (present(suth)) s = suth
   rhow = ps / (GAS * tw) ; muw = mu0 * law(tw, tref, omega, s)
   tau = muw * us / h
   do it=1, 400
      d = 0._R8P
      do i=0, M
         dmp = 1._R8P - exp(-y(i) * sqrt(tau * rhow) / (muw * APLUS))
         if (present(semilocal)) then
            if (semilocal) dmp = 1._R8P - exp(-y(i) * sqrt(tau * ps / (GAS * t(i))) / (mu0 * law(t(i), tref, omega, s) * APLUS))
         endif
         mt(i) = KAPPA * y(i) * sqrt(ps / (GAS * t(i)) * tau) * dmp**2
         a(i) = mu0 * law(t(i), tref, omega, s) / pr + mt(i) / 0.9_R8P
         mt(i) = mu0 * law(t(i), tref, omega, s) + mt(i)
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
