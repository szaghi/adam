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
!< with `mu = mu0 L(T, omega_mu)`, `k = k0 L(T, omega_k)` (`temperature_law`): `L = (T / T_ref)^omega`, times the
!< Sutherland factor `(T_ref + S) / (T + S)` when `S > 0` and `omega > 0` (Sutherland's law is `omega = 3/2` with it). A
!< constant law has `T_ref = 1`, zero exponents and `S = 0`, and `x**0 = 1` exactly, so it costs flops; the Sutherland
!< test is on scalars, the same in every thread. The gradient is stored `g(a,b) = d W_a / d x_b`
!< with `W = (u, v, w, T)`.
!<
!< The Ohmic flux (MHD, issue #65 P3), with the magnetic diffusivity `eta` and `J = curl B` (code units), follows from
!< `dB/dt = -curl(E)`, `E = -u x B + eta J`, and the Poynting flux `E x B`:
!<```
!< field B_i: -eta (dB_i/dx_d - dB_d/dx_i)                      (the curl form: eps_idk eta J_k)
!< energy:    eta (J x B)_d = eta (B_j dB_d/dx_j - B_j dB_j/dx_d)
!<```
!< The curl form differs from `-eta grad B` by `eta grad(div B)`, which is not zero under GLM: the curl form keeps the
!< resistive term out of the divergence constraint. (#65 section 2 printed the field flux with the opposite sign, an
!< anti-diffusion; the sign here is the derived one, and VV-5 measures the decay rate.) The gradient is stored
!< `g(a,b) = d B_a / d x_b`.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_dissipative_flux
public :: compute_resistive_flux
public :: compute_turbulence_integrands
public :: dissipative_diffusivity
public :: temperature_law

contains
   ! public procedures
   pure subroutine compute_dissipative_flux(d, mu0, k0, tref, suth, omega_mu, omega_k, w, g, f)
   !< Compute the dissipative flux of direction `d` from the state `w = (u, v, w, T)` and its gradient `g`: `f(1:3)` the
   !< momentum components, `f(4)` the energy.
   integer(I4P), intent(in)  :: d                  !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: mu0, k0            !< Viscosity and conductivity at the reference temperature.
   real(R8P),    intent(in)  :: tref, suth         !< Reference temperature, Sutherland constant (0: power law).
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

   mu  = mu0 * temperature_law(t=w(4), tref=tref, suth=suth, omega=omega_mu)
   k   = k0  * temperature_law(t=w(4), tref=tref, suth=suth, omega=omega_k)
   div = g(1,1) + g(2,2) + g(3,3)
   do i=1, 3
      tau(i) = mu * (g(i,d) + g(d,i))
   enddo
   tau(d) = tau(d) - 2._R8P / 3._R8P * mu * div
   f(1:3) = -tau
   f(4)   = -(w(1) * tau(1) + w(2) * tau(2) + w(3) * tau(3)) - k * g(4,d)
   endsubroutine compute_dissipative_flux

   pure subroutine compute_resistive_flux(d, eta, b, g, f)
   !< Compute the Ohmic flux of direction `d` from the field `b` and its gradient `g`: `f(1:3)` the field components,
   !< `f(4)` the energy.
   integer(I4P), intent(in)  :: d      !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: eta    !< Magnetic diffusivity.
   real(R8P),    intent(in)  :: b(3)   !< Magnetic field.
   real(R8P),    intent(in)  :: g(3,3) !< Gradient, g(a,b) = d B_a / d x_b.
   real(R8P),    intent(out) :: f(4)   !< Flux: field (3), energy.
   integer(I4P)              :: i      !< Counter.
   !$acc routine seq
   !$omp declare target

   do i=1, 3
      f(i) = -eta * (g(i,d) - g(d,i))
   enddo
   f(4) = eta * (b(1) * (g(d,1) - g(1,d)) + b(2) * (g(d,2) - g(2,d)) + b(3) * (g(d,3) - g(3,d)))
   endsubroutine compute_resistive_flux

   pure subroutine compute_turbulence_integrands(mu0, tref, suth, omega_mu, rho, p, T, u, g, f)
   !< Return the integrands of the turbulence history at a cell (issue #119, T0.1), per unit volume: `f(1)` kinetic energy
   !< `rho |u|^2 / 2`, `f(2)` enstrophy `rho |omega|^2 / 2`, `f(3)` solenoidal dissipation `mu |omega|^2`, `f(4)` dilatational
   !< dissipation `4/3 mu theta^2`, `f(5)` the viscous work `tau_ij du_i/dx_j`, `f(6)` the pressure dilatation `p theta`.
   real(R8P), intent(in)  :: mu0, tref, suth, omega_mu !< Viscosity law, as in `compute_dissipative_flux`.
   real(R8P), intent(in)  :: rho, p, T           !< Density, pressure, temperature.
   real(R8P), intent(in)  :: u(3)                !< Velocity.
   real(R8P), intent(in)  :: g(3,3)              !< Velocity gradient, g(a,b) = d u_a / d x_b.
   real(R8P), intent(out) :: f(6)                !< Integrands.
   real(R8P)              :: mu                  !< Viscosity at the temperature.
   real(R8P)              :: theta               !< Velocity divergence.
   real(R8P)              :: w1, w2, w3          !< Vorticity.
   real(R8P)              :: w_sq                !< Vorticity squared.
   !$acc routine seq
   !$omp declare target

   mu    = mu0 * temperature_law(t=T, tref=tref, suth=suth, omega=omega_mu)
   theta = g(1,1) + g(2,2) + g(3,3)
   w1    = g(3,2) - g(2,3)
   w2    = g(1,3) - g(3,1)
   w3    = g(2,1) - g(1,2)
   w_sq  = w1 * w1 + w2 * w2 + w3 * w3
   f(1)  = 0.5_R8P * rho * (u(1) * u(1) + u(2) * u(2) + u(3) * u(3))
   f(2)  = 0.5_R8P * rho * w_sq
   f(3)  = mu * w_sq
   f(4)  = 4._R8P / 3._R8P * mu * theta * theta
   f(5)  = mu * ((g(1,2) + g(2,1)) * (g(1,2) + g(2,1)) + (g(1,3) + g(3,1)) * (g(1,3) + g(3,1)) +    &
                 (g(2,3) + g(3,2)) * (g(2,3) + g(3,2)) +                                             &
                 2._R8P * (g(1,1) * g(1,1) + g(2,2) * g(2,2) + g(3,3) * g(3,3))) -                   &
           2._R8P / 3._R8P * mu * theta * theta
   f(6)  = p * theta
   endsubroutine compute_turbulence_integrands

   pure function dissipative_diffusivity(rho, T, t_wall, gamma, cp, mu0, k0, eta, tref, suth, omega_mu, omega_k) result(nu)
   !< Return the largest diffusivity of a cell, `max(4/3 mu / rho, gamma k / (rho cp), eta)`, the one of the diffusive
   !< time step limit (issue #65, D-M4-3).
   !<
   !< The laws are evaluated at `max(T, t_wall^2 / T)`: `t_wall` is the hottest isothermal wall temperature (0 without
   !< one, or with constant laws), and `t_wall^2 / T` the geometric-mirror ghost temperature of the cell beside it,
   !< which bounds the temperature of the wall face (both laws increase with T, Sutherland's too). Without it a power law at a wall much hotter than the gas gives
   !< the wall face a coefficient far above every cell's, and the cell bound is not a bound (a wall 330 times hotter than
   !< the gas blew up at CFL 0.5 and ran at 0.01). The exponents are non-negative, so the larger temperature is the bound.
   real(R8P), intent(in) :: rho, T             !< Density and temperature.
   real(R8P), intent(in) :: t_wall             !< Hottest isothermal wall temperature, 0 if none.
   real(R8P), intent(in) :: gamma, cp          !< Specific heats ratio and specific heat at constant pressure.
   real(R8P), intent(in) :: mu0, k0, eta       !< Viscosity, conductivity (at the reference temperature), resistivity.
   real(R8P), intent(in) :: tref, suth         !< Reference temperature, Sutherland constant (0: power law).
   real(R8P), intent(in) :: omega_mu, omega_k  !< Exponents of the laws.
   real(R8P)             :: nu                 !< Largest diffusivity.
   real(R8P)             :: te                 !< Temperature of the laws.
   !$acc routine seq
   !$omp declare target

   te = max(T, t_wall**2 / T)
   nu = max(4._R8P / 3._R8P * mu0 * temperature_law(t=te, tref=tref, suth=suth, omega=omega_mu) / rho, &
            gamma * k0 * temperature_law(t=te, tref=tref, suth=suth, omega=omega_k) / (rho * cp), eta)
   endfunction dissipative_diffusivity

   pure function temperature_law(t, tref, suth, omega) result(r)
   !< Return the ratio of a transport coefficient to its reference value, `(t / tref)^omega`, times the Sutherland factor
   !< `(tref + suth) / (t + suth)` when `suth > 0` and `omega > 0`: the power law for `suth = 0` (the same operations
   !< as before Sutherland, so the power and constant laws are unchanged bit for bit), Sutherland's law for `omega = 3/2`.
   !< A zero exponent is a constant coefficient under either law (the conductivity not given by the Prandtl number).
   real(R8P), intent(in) :: t     !< Temperature.
   real(R8P), intent(in) :: tref  !< Reference temperature.
   real(R8P), intent(in) :: suth  !< Sutherland constant (code temperature units), 0 for the power law.
   real(R8P), intent(in) :: omega !< Exponent.
   real(R8P)             :: r     !< Ratio.
   !$acc routine seq
   !$omp declare target

   r = (t / tref)**omega
   if (suth > 0._R8P .and. omega > 0._R8P) r = r * (tref + suth) / (t + suth)
   endfunction temperature_law
endmodule adam_flume_dissipation_library
