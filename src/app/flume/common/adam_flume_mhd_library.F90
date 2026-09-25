!< ADAM, FLUME pointwise ideal MHD physics shared by the CPU and FNL backends.
module adam_flume_mhd_library
!< ADAM, FLUME pointwise ideal MHD physics shared by the CPU and FNL backends.
!<
!< Same contract as `adam_flume_euler_library`: every routine is `pure`, takes explicit-size dummies and is tagged
!< `!$acc routine seq` + `!$omp declare target`, and never branches on the model (issue #41, D-5). Code units
!< (issue #41, D-2): `B` is rationalised (`B_SI / sqrt(mu0)`), so the magnetic pressure is `|B|^2 / 2` and
!<```
!< E = p/(gamma-1) + rho |u|^2 / 2 + |B|^2 / 2,     p = (gamma-1) (E - rho |u|^2 / 2 - |B|^2 / 2),     T = p / (rho R)
!<```
!< The routines work on the first `NV_MHD` conservative variables: the GLM variant (`NV_MHD_GLM`) passes its whole
!< state, whose leading `NV_MHD` entries are the same (sequence association); `psi` is not part of the energy (D-3).
!<
!< Scope in M2-P1: conversions and the fast magnetosonic speed (the dt and the auxiliary variables of a zero-residual
!< run); the eigensystem and the fluxes come with M2-P2.

! FLUME modules
use :: adam_flume_parameters, only : IA_A, IA_BX, IA_BY, IA_BZ, IA_H, IA_P, IA_R, IA_T, IA_U, IA_V, IA_W, &
                                     IQ_BX, IQ_BY, IQ_BZ, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX_MHD, NV_MHD
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: mhd_conservative_to_auxiliary
public :: mhd_fast_speed
public :: mhd_primitive_to_conservative

contains
   ! public procedures
   pure subroutine mhd_conservative_to_auxiliary(gamma, R, q, qa)
   !< Compute the auxiliary (primitive and derived) variables of a cell from its conservative variables.
   !<
   !< `qa(IA_H)` is the total specific enthalpy including the magnetic pressure, `(E + p + |B|^2 / 2) / rho`; the
   !< magnetic field is copied into `qa(IA_BX:IA_BZ)` so the eigen-routines read one array.
   real(R8P), intent(in)  :: gamma          !< Specific heats ratio.
   real(R8P), intent(in)  :: R              !< Gas constant.
   real(R8P), intent(in)  :: q(NV_MHD)      !< Conservative variables.
   real(R8P), intent(out) :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P)              :: pb             !< Magnetic pressure, |B|^2 / 2.
   !$acc routine seq
   !$omp declare target

   pb = 0.5_R8P * (q(IQ_BX)**2 + q(IQ_BY)**2 + q(IQ_BZ)**2)
   qa(IA_R)  = q(IQ_R)
   qa(IA_U)  = q(IQ_RU) / q(IQ_R)
   qa(IA_V)  = q(IQ_RV) / q(IQ_R)
   qa(IA_W)  = q(IQ_RW) / q(IQ_R)
   qa(IA_P)  = (gamma - 1._R8P) * (q(IQ_RE) - 0.5_R8P * q(IQ_R) * (qa(IA_U)**2 + qa(IA_V)**2 + qa(IA_W)**2) - pb)
   qa(IA_T)  = qa(IA_P) / (q(IQ_R) * R)
   qa(IA_H)  = (q(IQ_RE) + qa(IA_P) + pb) / q(IQ_R)
   qa(IA_A)  = sqrt(gamma * qa(IA_P) / q(IQ_R))
   qa(IA_BX) = q(IQ_BX)
   qa(IA_BY) = q(IQ_BY)
   qa(IA_BZ) = q(IQ_BZ)
   endsubroutine mhd_conservative_to_auxiliary

   pure function mhd_fast_speed(d, qa) result(cf)
   !< Return the fast magnetosonic speed of a cell along direction `d`.
   !<
   !< `cf^2 = (a^2 + b^2 + sqrt((a^2 - b^2)^2 + 4 a^2 bt^2)) / 2` with `b^2 = |B|^2 / rho` and `bt^2` its part transverse
   !< to `d`: the discriminant `(a^2 + b^2)^2 - 4 a^2 bn^2` written as a sum of non-negative terms, so round-off cannot
   !< make it negative when `B` is aligned with `d`.
   integer(I4P), intent(in) :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in) :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P)                :: cf             !< Fast magnetosonic speed.
   real(R8P)                :: a2             !< Squared sound speed.
   real(R8P)                :: b2             !< Squared Alfven speed, |B|^2 / rho.
   real(R8P)                :: bt2            !< Transverse part of b2.
   !$acc routine seq
   !$omp declare target

   a2  = qa(IA_A)**2
   b2  = (qa(IA_BX)**2 + qa(IA_BY)**2 + qa(IA_BZ)**2) / qa(IA_R)
   bt2 = max(b2 - qa(IA_BX+d-1)**2 / qa(IA_R), 0._R8P)
   cf  = sqrt(0.5_R8P * (a2 + b2 + sqrt((a2 - b2)**2 + 4._R8P * a2 * bt2)))
   endfunction mhd_fast_speed

   pure subroutine mhd_primitive_to_conservative(gamma, r, u, v, w, p, bx, by, bz, q)
   !< Compute the conservative variables of a cell from its primitive variables (`psi`, if any, is set by the caller).
   real(R8P), intent(in)  :: gamma      !< Specific heats ratio.
   real(R8P), intent(in)  :: r          !< Density.
   real(R8P), intent(in)  :: u, v, w    !< Velocity components.
   real(R8P), intent(in)  :: p          !< Pressure.
   real(R8P), intent(in)  :: bx, by, bz !< Magnetic field components.
   real(R8P), intent(out) :: q(NV_MHD)  !< Conservative variables.
   !$acc routine seq
   !$omp declare target

   q(IQ_R)  = r
   q(IQ_RU) = r * u
   q(IQ_RV) = r * v
   q(IQ_RW) = r * w
   q(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * r * (u**2 + v**2 + w**2) + 0.5_R8P * (bx**2 + by**2 + bz**2)
   q(IQ_BX) = bx
   q(IQ_BY) = by
   q(IQ_BZ) = bz
   endsubroutine mhd_primitive_to_conservative
endmodule adam_flume_mhd_library
