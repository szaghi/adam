!< ADAM, FLUME pointwise isothermal Euler physics shared by the CPU and FNL backends (issue #100, M13).
module adam_flume_euler_iso_library
!< ADAM, FLUME pointwise isothermal Euler physics shared by the CPU and FNL backends (issue #100, M13).
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`,
!< like adam_flume_euler_library, whose public interface it repeats routine by routine (the model-agnostic kernel
!< bodies call either one). Isothermal gas, `p = c0^2 rho`, `a = c0`: the energy slot is inert (D-M13-1), and the EOS
!< argument named `gamma` carries c0^2. The eigensystem, the Roe average and the HLLC flux are documented in the body.

! FLUME modules
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

contains
#include "adam_flume_euler_iso_library_body.INC"

   ! EOS hooks of the body: the isothermal gas, `gamma` = c0^2
   pure function eos_pressure(gamma, rho) result(p)
   !< P = c0^2 rho.
   real(R8P), intent(in) :: gamma !< c0^2.
   real(R8P), intent(in) :: rho   !< Density.
   real(R8P)             :: p     !< Pressure.
   !$acc routine seq
   !$omp declare target

   p = gamma * rho
   endfunction eos_pressure

   pure function eos_temperature(gamma, R, rho) result(t)
   !< T = P / (rho R) = c0^2 / R.
   real(R8P), intent(in) :: gamma !< c0^2.
   real(R8P), intent(in) :: R     !< Gas constant.
   real(R8P), intent(in) :: rho   !< Density (unused: T is uniform).
   real(R8P)             :: t     !< Temperature.
   !$acc routine seq
   !$omp declare target

   t = gamma / R
   endfunction eos_temperature

   pure function eos_sound_speed(gamma, rho) result(a)
   !< a = c0.
   real(R8P), intent(in) :: gamma !< c0^2.
   real(R8P), intent(in) :: rho   !< Density (unused: a is uniform).
   real(R8P)             :: a     !< Speed of sound.
   !$acc routine seq
   !$omp declare target

   a = sqrt(gamma)
   endfunction eos_sound_speed

   pure function eos_roe_sound_speed(gamma, rhoL, rhoR, pL, pR) result(a)
   !< The speed of the Roe average: c0 (the secant [P] / [rho] is c0^2 exactly).
   real(R8P), intent(in) :: gamma      !< c0^2.
   real(R8P), intent(in) :: rhoL, rhoR !< Densities (unused).
   real(R8P), intent(in) :: pL, pR     !< Pressures (unused).
   real(R8P)             :: a          !< Speed of sound of the average.
   !$acc routine seq
   !$omp declare target

   a = sqrt(gamma)
   endfunction eos_roe_sound_speed
endmodule adam_flume_euler_iso_library
