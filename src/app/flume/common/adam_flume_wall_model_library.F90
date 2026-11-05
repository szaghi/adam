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
public :: sample_cells
public :: wall_eddy_viscosity
public :: wall_model_solve
public :: wall_model_solve_inverse
public :: inverse_strain
public :: WM_ITERATIONS, WM_POINTS

integer(I4P), parameter :: WM_POINTS=40_I4P     !< Intervals of the layer; tau within 1 % up to h+ = 1e5, 5.5 % off at 1e6.
integer(I4P), parameter :: WM_ITERATIONS=8_I4P  !< Iterations, Aitken every second one.
real(R8P),    parameter :: WM_STRETCH=1.3_R8P   !< Ratio of consecutive intervals.
real(R8P),    parameter :: WM_KAPPA=0.41_R8P    !< von Karman constant.
real(R8P),    parameter :: WM_APLUS=17._R8P     !< Damping constant.
real(R8P),    parameter :: WM_PRT=0.9_R8P       !< Turbulent Prandtl number.
real(R8P),    parameter :: WM_SRA=1.14_R8P      !< Reynolds analogy factor of the inverse model (Griffin et al. 2023, after Zhang et al. 2014).
integer(I4P), parameter :: WM_SHOTS=8_I4P       !< Shots of the inverse model, secant on ln(tau) after the first.
integer(I4P), parameter :: WM_SUBITERATIONS=2_I4P !< Corrections of each step of the inverse model.

contains
   pure function wall_eddy_viscosity(y, rho, mu, tau, rhow, muw, semilocal) result(mt)
   !< The eddy viscosity of the layer, `kappa y sqrt(rho tau_w) D^2`, `D = 1 - exp(-y+ / A+)` with `y+ = y sqrt(tau_w rho_w) / mu_w` (wall units) or, `semilocal`,
   !< `y* = y sqrt(rho tau_w) / mu` (local density and viscosity): the `mu_t` of `wall_model_solve`, at the height `y` of a cell below the matching height with that
   !< cell's density and viscosity (issue #119, `wall_model_eddy_viscosity`).
   real(R8P), intent(in) :: y, rho, mu   !< Height, density and viscosity of the cell.
   real(R8P), intent(in) :: tau          !< Wall shear stress.
   real(R8P), intent(in) :: rhow, muw    !< Density and viscosity at the wall.
   logical,   intent(in) :: semilocal    !< Semi-local damping.
   real(R8P)             :: mt           !< Eddy viscosity.
   real(R8P)             :: damp         !< Damping.
   !$acc routine seq
   !$omp declare target

   if (semilocal) then
      damp = 1._R8P - exp(-y * sqrt(rho * tau) / (mu * WM_APLUS))
   else
      damp = 1._R8P - exp(-y * sqrt(tau * rhow) / (muw * WM_APLUS))
   endif
   mt = WM_KAPPA * y * sqrt(rho * tau) * damp * damp
   endfunction wall_eddy_viscosity

   pure subroutine sample_cells(hwm, dx, n, m, w, h)
   !< The cells of a wall-model sample: from the matching height `hwm` (0: the 3rd cell), the cell m counted from the wall and the weight w of cell
   !< m + 1, both within the `n` cells of the block, and the height of the sample actually taken.
   real(R8P),    intent(in)  :: hwm, dx !< Matching height, spacing along the normal.
   integer(I4P), intent(in)  :: n       !< Cells of the block along the normal.
   integer(I4P), intent(out) :: m       !< Cell counted from the wall.
   real(R8P),    intent(out) :: w       !< Weight of cell m + 1.
   real(R8P),    intent(out) :: h       !< Height of the sample.
   real(R8P)                 :: s       !< Height in cells, from the wall.
   !$acc routine seq
   !$omp declare target

   if (hwm > 0._R8P) then
      s = max(hwm / dx + 0.5_R8P, 1._R8P)
      m = int(s, I4P) ; w = s - real(m, R8P)
      if (m >= n) then
         m = n - 1_I4P ; w = 1._R8P ! the block holds no farther cell: its last one
      endif
      h = (real(m, R8P) + w - 0.5_R8P) * dx
   else
      m = 3_I4P ; w = 0._R8P ; h = 2.5_R8P * dx
   endif
   endsubroutine sample_cells

   pure subroutine wall_model_solve(h, us, ts, ps, tw, R, cp, pr, mu0, tref, omega_mu, tau, q, semilocal)
   !< Return the wall shear stress `tau` and the heat flux `q` (`k dT/dy` at the wall, positive when the gas is hotter than the wall) of an isothermal wall.
   real(R8P), intent(in)  :: h, us, ts, ps       !< Height, parallel speed, temperature and pressure of the sample.
   real(R8P), intent(in)  :: tw                  !< Wall temperature.
   real(R8P), intent(in)  :: R, cp, pr           !< Gas constant, specific heat, Prandtl number.
   real(R8P), intent(in)  :: mu0, tref, omega_mu !< Viscosity law, as `compute_dissipative_flux`.
   real(R8P), intent(out) :: tau, q              !< Wall shear stress, wall heat flux.
   logical,   intent(in), optional :: semilocal  !< Semi-local damping, y* = y sqrt(rho tau) / mu with the local rho, mu (else wall units).
   logical                :: sl                  !< Local copy.
   real(R8P)              :: y(0:WM_POINTS)      !< Heights.
   real(R8P)              :: t(0:WM_POINTS)      !< Temperature.
   real(R8P)              :: u(0:WM_POINTS)      !< Velocity.
   real(R8P)              :: mt(0:WM_POINTS)     !< mu + mu_t.
   real(R8P)              :: a(0:WM_POINTS)      !< mu / Pr + mu_t / Pr_t.
   real(R8P)              :: dy0, mu, rho, rhow, muw, ut, damp, i2, i3, dy, hp, tnew, x0, x1, x2 !< Work.
   integer(I4P)           :: n, it               !< Counters.
   !$acc routine seq
   !$omp declare target

   sl = .false. ; if (present(semilocal)) sl = semilocal
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
         if (sl) then
            damp = 1._R8P - exp(-y(n) * sqrt(rho * tau) / (mu * WM_APLUS))
         else
            damp = 1._R8P - exp(-y(n) * ut * rhow / (muw * WM_APLUS))
         endif
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

   pure subroutine wall_model_solve_inverse(h, us, ts, ps, tw, ue, te, R, cp, pr, mu0, tref, omega_mu, tau, q)
   !< Return the wall shear stress `tau` and the heat flux `q` (`k dT/dy` at the wall, positive when the gas is hotter than the wall) of an isothermal wall by the model of
   !< Griffin, Fu & Moin (2023, JFM 970 A36): the strain of the constant-property layer at the same semi-local height, `S_t = 1 / (1 + kappa y* D^2)`,
   !< `D = 1 - exp(-y* / A+)`, `y* = y sqrt(rho tau) / mu`, mapped to the compressible one by the inverse total-stress-based transformation (their Eq. 2.2),
   !< with the temperature from the matching-anchored Reynolds analogy (Eq. 2.5)
   !< `T = Tw + s Pr (Tr - Tw) (U / Ue) (1 - U / Us) + (U / Us)^2 (Ts - Tw)`, `Tr = Te + Pr^(1/3) Ue^2 / (2 cp)`, s 1.14; so `q = s cp tau (Tr - Tw) / Ue`.
   !< `Ue`, `Te`: the state at the edge of the layer (the channel centre), given.
   !<
   !< Eq. 2.2 written for the physical strain g = dU/dy, its property gradients through T(U): with `E = (1 / S_t - 1) mu_w / mu`, `F = sqrt(rho / rho_w)`,
   !< `K = tau F / mu`, `beta = (dT/dU) (1/2 + omega_mu) / T` (ideal gas at uniform pressure, power-law viscosity), `F beta y g^2 - (E + F + K beta y) g + K = 0`;
   !< its positive root when `beta` < 0 (the other is negative), else the one nearer `K / (E + F)`, the strain without the gradients (the exact `tau / mu` of the
   !< laminar layer, where the equation also holds the spurious root `1 / (beta y)` of `dy*/dy = 0`). The layer is integrated from the wall on the stretched grid of `wall_model_solve`, each interval
   !< exact for 1 / g linear in y (the viscous and the log layer), corrected `WM_SUBITERATIONS` times, then shot on ln(tau) to the sample speed.
   real(R8P), intent(in)  :: h, us, ts, ps       !< Height, parallel speed, temperature and pressure of the sample.
   real(R8P), intent(in)  :: tw                  !< Wall temperature.
   real(R8P), intent(in)  :: ue, te              !< Speed and temperature at the edge of the layer.
   real(R8P), intent(in)  :: R, cp, pr           !< Gas constant, specific heat, Prandtl number.
   real(R8P), intent(in)  :: mu0, tref, omega_mu !< Viscosity law, as `compute_dissipative_flux`.
   real(R8P), intent(out) :: tau, q              !< Wall shear stress, wall heat flux.
   real(R8P)              :: rhow, muw, c1, hp, x0, x1, f0, f1, slope, dx !< Work.
   integer(I4P)           :: it                  !< Counter.
   !$acc routine seq
   !$omp declare target

   rhow = ps / (R * tw)
   muw = mu0 * (tw / tref)**omega_mu
   ! dT/dU at the wall: s Pr (Tr - Tw) / Ue
   c1 = WM_SRA * pr * (te + pr**(1._R8P / 3._R8P) * ue * ue / (2._R8P * cp) - tw) / ue
   q = muw * cp / pr * (ts - tw) / h
   tau = 0._R8P
   if (us <= 0._R8P) return
   ! the first guess of wall_model_solve: laminar, or the log law when above it
   tau = muw * us / h
   hp = h * sqrt(tau / rhow) * rhow / muw
   do it=1, 3
      hp = max(hp, 1.5_R8P)
      hp = (us / (log(hp) / WM_KAPPA + 5._R8P)) * h * rhow / muw
   enddo
   slope = 1._R8P ! d ln U(h) / d ln tau: 1 in the viscous layer, about 1/2 in the log layer
   if (hp > 11._R8P) then
      tau = rhow * (hp * muw / (h * rhow))**2
      slope = 0.6_R8P
   endif
   x1 = log(tau)
   f1 = log(inverse_layer_speed(h, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu) / us)
   do it=1, WM_SHOTS
      x0 = x1 ; f0 = f1
      dx = min(max(-f0 / slope, -1._R8P), 1._R8P) ! at most a factor e per shot
      x1 = x0 + dx
      f1 = log(inverse_layer_speed(h, exp(x1), us, ts, ps, tw, c1, R, mu0, tref, omega_mu) / us)
      if (abs(x1 - x0) > 1.e-14_R8P) slope = min(max((f1 - f0) / (x1 - x0), 0.2_R8P), 1.5_R8P)
   enddo
   tau = exp(x1 + min(max(-f1 / slope, -1._R8P), 1._R8P)) ! a last secant step
   q = cp * c1 * tau / pr
   if (.not.(tau < huge(1._R8P) .and. tau >= 0._R8P .and. abs(q) < huge(1._R8P))) then ! a last resort
      tau = muw * us / h
      q = muw * cp / pr * (ts - tw) / h
   endif
   endsubroutine wall_model_solve_inverse

   pure function inverse_layer_speed(h, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu) result(u)
   !< The speed at the height `h` of the layer of wall stress `tau` of the inverse model (`wall_model_solve_inverse`).
   real(R8P), intent(in) :: h, tau, us, ts, ps, tw  !< Height, wall stress, sample speed, temperature, pressure, wall temperature.
   real(R8P), intent(in) :: c1                      !< dT/dU at the wall.
   real(R8P), intent(in) :: R, mu0, tref, omega_mu  !< Gas constant, viscosity law.
   real(R8P)             :: u                       !< Speed at h.
   real(R8P)             :: dy0, y0, y1, g0, g1, u1 !< Work.
   integer(I4P)          :: n, k                    !< Counters.
   !$acc routine seq
   !$omp declare target

   dy0 = h * (WM_STRETCH - 1._R8P) / (WM_STRETCH**WM_POINTS - 1._R8P)
   y0 = 0._R8P ; u = 0._R8P
   g0 = inverse_strain(y0, u, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu)
   do n=1, WM_POINTS
      y1 = y0 + dy0 * WM_STRETCH**(n-1)
      u1 = u + (y1 - y0) * g0
      do k=1, WM_SUBITERATIONS
         g1 = inverse_strain(y1, u1, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu)
         u1 = u + (y1 - y0) * log_mean_inverse(1._R8P / g0, 1._R8P / g1)
      enddo
      g0 = inverse_strain(y1, u1, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu)
      u = u1 ; y0 = y1
   enddo
   endfunction inverse_layer_speed

   pure function inverse_strain(y, u, tau, us, ts, ps, tw, c1, R, mu0, tref, omega_mu) result(g)
   !< The strain dU/dy at the height `y` and speed `u` of the layer of wall stress `tau` of the inverse model (`wall_model_solve_inverse`).
   real(R8P), intent(in) :: y, u, tau, us, ts, ps, tw !< Height, speed, wall stress, sample speed, temperature, pressure, wall temperature.
   real(R8P), intent(in) :: c1                        !< dT/dU at the wall.
   real(R8P), intent(in) :: R, mu0, tref, omega_mu    !< Gas constant, viscosity law.
   real(R8P)             :: g                         !< Strain.
   real(R8P)             :: t, dtdu, rho, mu, ys, d, e, f, k, b, p, a, dsc, glin, g2 !< Work.
   !$acc routine seq
   !$omp declare target

   t = max(tw + c1 * u * (1._R8P - u / us) + (u / us)**2 * (ts - tw), 0.1_R8P * min(tw, ts))
   dtdu = c1 * (1._R8P - 2._R8P * u / us) + 2._R8P * u * (ts - tw) / (us * us)
   rho = ps / (R * t)
   mu = mu0 * (t / tref)**omega_mu
   ys = y * sqrt(rho * tau) / mu
   d = 1._R8P - exp(-ys / WM_APLUS)
   e = WM_KAPPA * ys * d * d * (mu0 * (tw / tref)**omega_mu) / mu ! (1 / S_t - 1) / mu+
   f = sqrt(tw / t)                                                ! sqrt(rho+), rho+ = Tw / T at uniform pressure
   k = tau * f / mu
   b = dtdu * (0.5_R8P + omega_mu) / t * y
   p = e + f + k * b
   a = f * b
   dsc = p * p - 4._R8P * a * k
   glin = k / (e + f)
   ! beta < 0 (the gas cooler away from the wall): the roots have opposite signs (product K / a < 0), the positive one is 2K / (p + sqrt(dsc)), p + sqrt(dsc) > 0;
   ! beta > 0: both positive (p > 0), the one nearer the strain without the gradients; no real root: that strain
   if (dsc < 0._R8P) then
      g = glin
   else
      g = 2._R8P * k / (p + sqrt(dsc))
      if (a > 0._R8P) then
         g2 = (p + sqrt(dsc)) / (2._R8P * a)
         if (abs(g2 - glin) < abs(g - glin)) g = g2
      endif
   endif
   endfunction inverse_strain

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
