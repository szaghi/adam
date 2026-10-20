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
endmodule adam_flume_euler_iso_library
