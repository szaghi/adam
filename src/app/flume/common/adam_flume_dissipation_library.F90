!< ADAM, FLUME pointwise dissipative physics shared by the CPU and FNL backends (issue #65, M4).
module adam_flume_dissipation_library
!< ADAM, FLUME pointwise dissipative physics shared by the CPU and FNL backends (issue #65, M4).
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`, as
!< in `adam_flume_euler_library`: the CPU loops and the FNL device kernels call the same source.
!<
!< The dissipative flux of direction `d` (issue #65, section 2), added to the inviscid one, is
!<```
!< momentum u_i: -tau_{di},  tau = mu (grad u + grad u^T - 2/3 div(u) I)   (Stokes hypothesis)
!< energy:       -u_i tau_{di} - k dT/dx_d
!<```
!< with `mu = mu0 (T / T_ref)^omega_mu`, `k = k0 (T / T_ref)^omega_k`: a constant law has `T_ref = 1` and zero
!< exponents, and `x**0 = 1` exactly, so it costs flops, never a branch. The gradient is stored `g(a,b) = d W_a / d x_b`
!< with `W = (u, v, w, T)`.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_dissipative_flux
public :: dissipative_diffusivity

contains
   ! public procedures
   pure subroutine compute_dissipative_flux(d, mu0, k0, tref, omega_mu, omega_k, w, g, f)
   !< Compute the dissipative flux of direction `d` from the state `w = (u, v, w, T)` and its gradient `g`: `f(1:3)` the
   !< momentum components, `f(4)` the energy.
   integer(I4P), intent(in)  :: d                  !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: mu0, k0            !< Viscosity and conductivity at the reference temperature.
   real(R8P),    intent(in)  :: tref               !< Reference temperature of the laws.
   real(R8P),    intent(in)  :: omega_mu, omega_k  !< Exponents of the laws.
   real(R8P),    intent(in)  :: w(4)               !< Velocity and temperature.
   real(R8P),    intent(in)  :: g(4,3)             !< Gradient, g(a,b) = d w_a / d x_b.
   real(R8P),    intent(out) :: f(4)               !< Flux: momentum (3), energy.
   real(R8P)                 :: mu, k              !< Transport coefficients at the temperature.
   real(R8P)                 :: div                !< Velocity divergence.
   real(R8P)                 :: tau(3)             !< Stress row tau_{d,:}.
   integer(I4P)              :: i                  !< Counter.
   !$acc routine seq
   !$omp declare target

   mu  = mu0 * (w(4) / tref)**omega_mu
   k   = k0  * (w(4) / tref)**omega_k
   div = g(1,1) + g(2,2) + g(3,3)
   do i=1, 3
      tau(i) = mu * (g(i,d) + g(d,i))
   enddo
   tau(d) = tau(d) - 2._R8P / 3._R8P * mu * div
   f(1:3) = -tau
   f(4)   = -(w(1) * tau(1) + w(2) * tau(2) + w(3) * tau(3)) - k * g(4,d)
   endsubroutine compute_dissipative_flux

   pure function dissipative_diffusivity(rho, T, t_wall, gamma, cp, mu0, k0, eta, tref, omega_mu, omega_k) result(nu)
   !< Return the largest diffusivity of a cell, `max(4/3 mu / rho, gamma k / (rho cp), eta)`, the one of the diffusive
   !< time step limit (issue #65, D-M4-3).
   !<
   !< The laws are evaluated at `max(T, t_wall^2 / T)`: `t_wall` is the hottest isothermal wall temperature (0 without
   !< one, or with constant laws), and `t_wall^2 / T` the geometric-mirror ghost temperature of the cell beside it,
   !< which bounds the temperature of the wall face. Without it a power law at a wall much hotter than the gas gives
   !< the wall face a coefficient far above every cell's, and the cell bound is not a bound (a wall 330 times hotter than
   !< the gas blew up at CFL 0.5 and ran at 0.01). The exponents are non-negative, so the larger temperature is the bound.
   real(R8P), intent(in) :: rho, T             !< Density and temperature.
   real(R8P), intent(in) :: t_wall             !< Hottest isothermal wall temperature, 0 if none.
   real(R8P), intent(in) :: gamma, cp          !< Specific heats ratio and specific heat at constant pressure.
   real(R8P), intent(in) :: mu0, k0, eta       !< Viscosity, conductivity (at the reference temperature), resistivity.
   real(R8P), intent(in) :: tref               !< Reference temperature of the laws.
   real(R8P), intent(in) :: omega_mu, omega_k  !< Exponents of the laws.
   real(R8P)             :: nu                 !< Largest diffusivity.
   real(R8P)             :: te                 !< Temperature of the laws.
   !$acc routine seq
   !$omp declare target

   te = max(T, t_wall**2 / T)
   nu = max(4._R8P / 3._R8P * mu0 * (te / tref)**omega_mu / rho, gamma * k0 * (te / tref)**omega_k / (rho * cp), eta)
   endfunction dissipative_diffusivity
endmodule adam_flume_dissipation_library
