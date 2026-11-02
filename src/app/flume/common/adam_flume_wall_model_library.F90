!< ADAM, FLUME equilibrium wall model (issue #119, T3.4).

#include "fundal.H"

module adam_flume_wall_model_library
!< ADAM, FLUME equilibrium wall model (issue #119, T3.4).
!<
!< The thin-layer equations between the wall and a sample at height `h`, parallel velocity `us`, temperature `ts`, pressure `ps` (uniform across the layer):
!< `d/dy[(mu + mu_t) dU/dy] = 0`, `d/dy[(mu + mu_t) U dU/dy + cp (mu/Pr + mu_t/Pr_t) dT/dy] = 0`, with
!< `mu_t = kappa y sqrt(rho tau_w) (1 - exp(-y+ / A+))^2`, with the local density, `y+ = y sqrt(tau_w rho_w) / mu_w`, kappa 0.41, A+ 17, Pr_t 0.9 (Kawai & Larsson 2012).
!< Both integrate once: `(mu + mu_t) dU/dy = tau_w`, `cp (mu/Pr + mu_t/Pr_t) dT/dy + tau_w U = Q_w`, `Q_w = k_w dT/dy|_w`. Given T(y), the two quadratures on a fixed stretched
!< grid give `tau_w` and `Q_w` from the sample values: a fixed number of iterations, so a fixed cost per wall face on the device, with no tridiagonal solve.

use :: penf, only : I4P, R8P

implicit none
private
public :: wall_model_solve
public :: WM_ITERATIONS, WM_POINTS

integer(I4P), parameter :: WM_POINTS=40_I4P     !< Intervals of the layer; tau within 1 % up to h+ = 1e5, 5.5 % off at 1e6.
integer(I4P), parameter :: WM_ITERATIONS=8_I4P  !< Iterations, Aitken every second one.
real(R8P),    parameter :: WM_STRETCH=1.3_R8P   !< Ratio of consecutive intervals.
real(R8P),    parameter :: WM_KAPPA=0.41_R8P    !< von Karman constant.
real(R8P),    parameter :: WM_APLUS=17._R8P     !< Damping constant.
real(R8P),    parameter :: WM_PRT=0.9_R8P       !< Turbulent Prandtl number.

contains
   pure subroutine wall_model_solve(h, us, ts, ps, tw, R, cp, pr, mu0, tref, omega_mu, tau, q)
   !< Return the wall shear stress `tau` and the heat flux `q` (`k dT/dy` at the wall, positive when the gas is hotter than the wall) of an isothermal wall.
   real(R8P), intent(in)  :: h, us, ts, ps       !< Height, parallel speed, temperature and pressure of the sample.
   real(R8P), intent(in)  :: tw                  !< Wall temperature.
   real(R8P), intent(in)  :: R, cp, pr           !< Gas constant, specific heat, Prandtl number.
   real(R8P), intent(in)  :: mu0, tref, omega_mu !< Viscosity law, as `compute_dissipative_flux`.
   real(R8P), intent(out) :: tau, q              !< Wall shear stress, wall heat flux.
   real(R8P)              :: y(0:WM_POINTS)      !< Heights.
   real(R8P)              :: t(0:WM_POINTS)      !< Temperature.
   real(R8P)              :: u(0:WM_POINTS)      !< Velocity.
   real(R8P)              :: mt(0:WM_POINTS)     !< mu + mu_t.
   real(R8P)              :: a(0:WM_POINTS)      !< mu / Pr + mu_t / Pr_t.
   real(R8P)              :: dy0, mu, rho, rhow, muw, ut, damp, i2, i3, dy, hp, tnew, x0, x1, x2 !< Work.
   integer(I4P)           :: n, it               !< Counters.
   !$acc routine seq
   !$omp declare target

   dy0 = h * (WM_STRETCH - 1._R8P) / (WM_STRETCH**WM_POINTS - 1._R8P)
   y(0) = 0._R8P
   t(0) = tw
   do n=1, WM_POINTS
      y(n) = y(n-1) + dy0 * WM_STRETCH**(n-1)
      t(n) = tw + (ts - tw) * y(n) / h
   enddo
   rhow = ps / (R * tw)
   muw = mu0 * (tw / tref)**omega_mu
   q = muw * cp / pr * (ts - tw) / h
   tau = 0._R8P
   if (us <= 0._R8P) return
   tau = muw * us / h
   ! the log law as a first guess, when above the laminar one: tau grows as (U kappa / ln(h+))^2
   hp = h * sqrt(tau / rhow) * rhow / muw
   do it=1, 3
      hp = max(hp, 1.5_R8P)
      hp = (us / (log(hp) / WM_KAPPA + 5._R8P)) * h * rhow / muw
   enddo
   if (hp > 11._R8P) tau = rhow * (hp * muw / (h * rhow))**2
   x0 = log(tau) ; x1 = x0
   do it=1, WM_ITERATIONS
      ut = sqrt(tau / rhow)
      do n=0, WM_POINTS
         mu = mu0 * (t(n) / tref)**omega_mu
         rho = ps / (R * t(n))
         damp = 1._R8P - exp(-y(n) * ut * rhow / (muw * WM_APLUS))
         mt(n) = WM_KAPPA * y(n) * sqrt(rho * tau) * damp * damp
         a(n) = mu / pr + mt(n) / WM_PRT
         mt(n) = mu + mt(n)
      enddo
      ! exact for a coefficient linear in y on each interval, so for the log layer
      u(0) = 0._R8P ; i2 = 0._R8P ; i3 = 0._R8P
      do n=1, WM_POINTS
         u(n) = u(n-1) + (y(n) - y(n-1)) * log_mean_inverse(mt(n-1), mt(n))
      enddo
      tnew = us / u(WM_POINTS)
      do n=0, WM_POINTS
         u(n) = u(n) * tnew
      enddo
      do n=1, WM_POINTS
         dy = (y(n) - y(n-1)) * log_mean_inverse(a(n-1), a(n))
         i2 = i2 + dy
         i3 = i3 + 0.5_R8P * dy * (u(n) + u(n-1))
      enddo
      q = (cp * (ts - tw) + tnew * i3) / i2
      i2 = 0._R8P ; i3 = 0._R8P
      do n=1, WM_POINTS
         dy = (y(n) - y(n-1)) * log_mean_inverse(a(n-1), a(n))
         i2 = i2 + dy
         i3 = i3 + 0.5_R8P * dy * (u(n) + u(n-1))
         ! under-relaxed and kept positive: with a poor tau the frictional heating term can drive the iterate below zero (NaN in a sweep at Mach 2 samples)
         t(n) = max(0.5_R8P * t(n) + 0.5_R8P * (tw + (q * i2 - tnew * i3) / cp), 0.1_R8P * min(tw, ts))
      enddo
      ! Aitken on log(tau) every second iteration: tau = g(tau) is a power law in the log layer
      x2 = log(tnew)
      if (mod(it, 2_I4P) == 0_I4P .and. abs(x2 - 2._R8P * x1 + x0) > 1.e-12_R8P) then
         x2 = min(max(x0 - (x1 - x0)**2 / (x2 - 2._R8P * x1 + x0), x1 - 1._R8P), x1 + 1._R8P) ! at most a factor e per extrapolation
      endif
      x0 = x1 ; x1 = x2
      if (mod(it, 2_I4P) == 0_I4P) x0 = x2
      tau = exp(x2)
   enddo
   if (.not.(tau < huge(1._R8P) .and. tau >= 0._R8P .and. abs(q) < huge(1._R8P))) then ! a last resort, never seen after the guards above
      tau = muw * us / h
      q = muw * cp / pr * (ts - tw) / h
   endif
   endsubroutine wall_model_solve

   pure function log_mean_inverse(a, b) result(r)
   !< The mean of 1 / c(y) over an interval where c is linear between `a` and `b`.
   real(R8P), intent(in) :: a, b !< End values.
   real(R8P)             :: r    !< Mean.
   real(R8P)             :: x    !< Relative difference.
   !$acc routine seq
   !$omp declare target

   x = (b - a) / (b + a)
   if (abs(x) < 1.e-3_R8P) then
      r = 2._R8P / (a + b) * (1._R8P + x * x / 3._R8P)
   else
      r = log(b / a) / (b - a)
   endif
   endfunction log_mean_inverse
endmodule adam_flume_wall_model_library
