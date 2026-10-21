!< ADAM, FLUME pointwise isothermal MHD physics and Riemann solvers shared by the CPU and FNL backends (issue #100, M13 P5).
module adam_flume_mhd_iso_library
!< ADAM, FLUME pointwise isothermal MHD physics and Riemann solvers shared by the CPU and FNL backends (issue #100, M13 P5).
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`. The
!< public interface repeats, routine by routine, the one adam_flume_mhd_library and adam_flume_mhd_riemann_library give
!< the MHD kernels and objects without divergence control and with the mixed GLM (EGLM is refused for these gases).
!< Isothermal gas, `p = c0^2 rho`, `a = c0`: the energy slot is inert (D-M13-1), and the EOS argument named `gamma`
!< carries c0^2. The body (adam_flume_mhd_iso_library_body.INC) documents the eigensystem (Kim et al. 1999) and the HLLD
!< solver (Mignone 2007).

! FLUME modules
use :: adam_flume_parameters, only : IA_A, IA_BX, IA_BY, IA_BZ, IA_H, IA_P, IA_R, IA_T, IA_U, IA_V, IA_W, IQ_BX, &
                                     IQ_BY, IQ_BZ, IQ_PSI, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX_MHD, NV_MHD,      &
                                     NV_MHD_GLM, S_MAX
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: mhd_backbone_flux
public :: mhd_conservative_to_auxiliary
public :: mhd_eigenvalues
public :: mhd_eigenvectors
public :: mhd_face_average
public :: mhd_face_flux_back_projection
public :: mhd_face_interpolation_fields
public :: mhd_face_split_fluxes
public :: mhd_face_states
public :: mhd_fast_speed
public :: mhd_flux
public :: mhd_frame_indexes
public :: mhd_glm_eigenvalues
public :: mhd_glm_eigenvectors
public :: mhd_glm_face_flux_back_projection
public :: mhd_glm_face_interpolation_fields
public :: mhd_glm_face_split_fluxes
public :: mhd_glm_face_states
public :: mhd_glm_flux
public :: mhd_glm_riemann_hll
public :: mhd_glm_riemann_hlld
public :: mhd_glm_riemann_llf
public :: mhd_primitive_to_conservative
public :: mhd_riemann_hll
public :: mhd_riemann_hlld
public :: mhd_riemann_llf
public :: mhd_state_magnitudes
public :: mhd_sum3

real(R8P), parameter :: EPS_BT=1.e-12_R8P   !< Degenerate transverse field: |B_t| <= EPS_BT max(|B|, sqrt(rho) a).
real(R8P), parameter :: EPS_FS=1.e-12_R8P   !< Degenerate fast-slow pair (triple umbilic): c_f^2 - c_s^2 <= EPS_FS c_f^2.
real(R8P), parameter :: EPS_HLLD=1.e-12_R8P !< Degenerate HLLD star: rho* |(S-S*_L)(S-S*_R)| <= EPS_HLLD (rho (S-u)^2 + B_n^2).

contains
#include "adam_flume_mhd_iso_library_body.INC"

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
endmodule adam_flume_mhd_iso_library
