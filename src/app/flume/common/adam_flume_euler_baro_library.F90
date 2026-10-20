!< ADAM, FLUME pointwise barotropic Euler physics shared by the CPU and FNL backends (issue #100, M13 P3).
module adam_flume_euler_baro_library
!< ADAM, FLUME pointwise barotropic Euler physics shared by the CPU and FNL backends (issue #100, M13 P3).
!<
!< The body of adam_flume_euler_iso_library on a piecewise polytrope, `[eos] type = barotropic`: the energy slot is
!< inert (D-M13-1) and the EOS argument named `gamma` carries the id of the realm's table in
!< adam_flume_eos_tables_library (as a real). Only the four EOS hooks below differ from the isothermal library; the
!< eigensystem, the Roe average and the LLF/HLL/HLLC fluxes are the body's, with c = sqrt(dP/drho) of each state.

! FLUME modules
use :: adam_flume_eos_tables_library, only : barotropic_c2, barotropic_pressure
use :: adam_flume_parameters, only : IA_A, IA_H, IA_P, IA_R, IA_S, IA_T, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                     IQ_RW, NV_AUX_K=>NV_AUX, NV_AUX_BASE=>NV_AUX, NV_BASE=>NV_EULER, NV_K=>NV_EULER, &
                                     NV_POSITIVE=>NV_EULER, S_MAX
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: compute_eigenvalues
public :: compute_eigenvectors
public :: compute_face_flux_back_projection
public :: compute_face_interpolation_fields
public :: compute_face_split_fluxes
public :: compute_face_states
public :: compute_flux
public :: compute_riemann_hll
public :: compute_riemann_hllc
public :: compute_riemann_llf
public :: compute_roe_average
public :: conservative_to_auxiliary
public :: primitive_to_conservative

real(R8P), parameter :: ROE_SECANT_MIN = 1.0e-6_R8P !< Relative density gap below which the Roe speed is c(rho mean).

contains
#include "adam_flume_euler_iso_library_body.INC"

   ! EOS hooks of the body: a piecewise polytrope, `gamma` = the table id
   pure function eos_pressure(gamma, rho) result(p)
   !< P(rho) of the table.
   real(R8P), intent(in) :: gamma !< Table id.
   real(R8P), intent(in) :: rho   !< Density.
   real(R8P)             :: p     !< Pressure.
   !$acc routine seq
   !$omp declare target

   p = barotropic_pressure(id=nint(gamma, I4P), rho=rho)
   endfunction eos_pressure

   pure function eos_temperature(gamma, R, rho) result(t)
   !< T = P / (rho R).
   real(R8P), intent(in) :: gamma !< Table id.
   real(R8P), intent(in) :: R     !< Gas constant.
   real(R8P), intent(in) :: rho   !< Density.
   real(R8P)             :: t     !< Temperature.
   !$acc routine seq
   !$omp declare target

   t = barotropic_pressure(id=nint(gamma, I4P), rho=rho) / (rho * R)
   endfunction eos_temperature

   pure function eos_sound_speed(gamma, rho) result(a)
   !< a = sqrt(dP/drho).
   real(R8P), intent(in) :: gamma !< Table id.
   real(R8P), intent(in) :: rho   !< Density.
   real(R8P)             :: a     !< Speed of sound.
   !$acc routine seq
   !$omp declare target

   a = sqrt(barotropic_c2(id=nint(gamma, I4P), rho=rho))
   endfunction eos_sound_speed

   pure function eos_roe_sound_speed(gamma, rhoL, rhoR, pL, pR) result(a)
   !< The speed of the Roe average, the secant sqrt([P] / [rho]): with it the Roe matrix satisfies [F] = A~ [U] for
   !< any P(rho). P increases with rho, so the secant is positive, between the two one-sided speeds of a segment and
   !< across breaks too. Below a relative density gap of ROE_SECANT_MIN the secant loses digits to cancellation (about
   !< eps / gap) while the mid-point speed errs by O(gap^2): there the speed is c(rho mean), both errors below 2e-10.
   real(R8P), intent(in) :: gamma      !< Table id.
   real(R8P), intent(in) :: rhoL, rhoR !< Densities.
   real(R8P), intent(in) :: pL, pR     !< Pressures.
   real(R8P)             :: a          !< Speed of sound of the average.
   !$acc routine seq
   !$omp declare target

   if (abs(rhoR - rhoL) > ROE_SECANT_MIN * max(rhoL, rhoR)) then
      a = sqrt((pR - pL) / (rhoR - rhoL))
   else
      a = sqrt(barotropic_c2(id=nint(gamma, I4P), rho=0.5_R8P * (rhoL + rhoR)))
   endif
   endfunction eos_roe_sound_speed
endmodule adam_flume_euler_baro_library
