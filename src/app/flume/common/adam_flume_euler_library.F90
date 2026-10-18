!< ADAM, FLUME pointwise Euler physics shared by the CPU and FNL backends.
module adam_flume_euler_library
!< ADAM, FLUME pointwise Euler physics shared by the CPU and FNL backends.
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`,
!< so CPU loops and FNL device kernels call the SAME source: CPU/FNL agreement is then a property of the loops only.
!< Calorically perfect ideal gas: `p = (gamma-1) (rE - rho |u|^2 / 2)`, `T = p / (rho R)`, `a = sqrt(gamma p / rho)`.
!<
!< Eigensystem in direction `d` (issue #35, section 3.3): normal `n = e_d`, tangents taken cyclically,
!< `t1 = e_{mod(d,3)+1}`, `t2 = e_{mod(d+1,3)+1}`; with `b2 = (gamma-1)/a^2`, `b1 = b2 |u|^2 / 2`, `u_n = u.n`:
!<```
!< r1 = (1, u - a n, H - a u_n)        l1 = (b1 + u_n/a, -b2 u - n/a, b2) / 2
!< r2 = (1, u, |u|^2 / 2)              l2 = (1 - b1, b2 u, -b2)
!< r3 = (0, t1, u.t1)                  l3 = (-u.t1, t1, 0)
!< r4 = (0, t2, u.t2)                  l4 = (-u.t2, t2, 0)
!< r5 = (1, u + a n, H + a u_n)        l5 = (b1 - u_n/a, -b2 u + n/a, b2) / 2
!<```
!< eigenvalues `(u_n - a, u_n, u_n, u_n, u_n + a)`. Storage, stated once: `el(k,:)` is the row `l_k`, `er(:,k)` is the
!< column `r_k`; projection `w_k = sum_v el(k,v) q_v`, back-projection `q_v = sum_k er(v,k) w_k`.
!<
!< The WENO reconstruction primitive differs between host (`adam_weno_object`) and device (`adam_fnl_weno_kernels`),
!< so the per-face flux is split in two physics halves around it (issue #35, section 7.1): `compute_face_split_fluxes`
!< projects and splits the stencil, the backend reconstructs each field, `compute_face_flux_back_projection` returns
!< to conservative variables.

! FLUME modules
use :: adam_flume_parameters, only : IA_A, IA_H, IA_P, IA_R, IA_S, IA_T, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, &
                                     IQ_RW, NV_AUX_K=>NV_AUX, NV_AUX_BASE=>NV_AUX, NV_BASE=>NV_EULER, NV_K=>NV_EULER, S_MAX
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
#include "adam_flume_euler_library_body.INC"
endmodule adam_flume_euler_library
