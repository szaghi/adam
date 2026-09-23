!< ADAM, FLUME pointwise Euler physics shared by the CPU and FNL backends.
module adam_flume_euler_library
!< ADAM, FLUME pointwise Euler physics shared by the CPU and FNL backends.
!<
!< Every routine is `pure`, takes explicit-size dummies and is tagged `!$acc routine seq` + `!$omp declare target`,
!< so CPU loops and FNL device kernels call the SAME source: CPU/FNL agreement is then a property of the loops only.
!< Calorically perfect ideal gas: `p = (gamma-1) (rE - rho |u|^2 / 2)`, `T = p / (rho R)`, `a = sqrt(gamma p / rho)`.

! FLUME modules
use :: adam_flume_parameters, only : IA_A, IA_H, IA_P, IA_R, IA_T, IA_U, IA_V, IA_W, &
                                     IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX, NV_EULER
! third party modules
use :: penf,                  only : R8P

implicit none
private
public :: conservative_to_auxiliary
public :: primitive_to_conservative

contains
   ! public procedures
   pure subroutine conservative_to_auxiliary(gamma, R, q, qa)
   !< Compute the auxiliary (primitive and derived) variables of a cell from its conservative variables.
   real(R8P), intent(in)  :: gamma       !< Specific heats ratio.
   real(R8P), intent(in)  :: R           !< Gas constant.
   real(R8P), intent(in)  :: q(NV_EULER) !< Conservative variables.
   real(R8P), intent(out) :: qa(NV_AUX)  !< Auxiliary variables.
   !$acc routine seq
   !$omp declare target

   qa(IA_R) = q(IQ_R)
   qa(IA_U) = q(IQ_RU) / q(IQ_R)
   qa(IA_V) = q(IQ_RV) / q(IQ_R)
   qa(IA_W) = q(IQ_RW) / q(IQ_R)
   qa(IA_P) = (gamma - 1._R8P) * (q(IQ_RE) - 0.5_R8P * q(IQ_R) * (qa(IA_U)**2 + qa(IA_V)**2 + qa(IA_W)**2))
   qa(IA_T) = qa(IA_P) / (q(IQ_R) * R)
   qa(IA_H) = (q(IQ_RE) + qa(IA_P)) / q(IQ_R)
   qa(IA_A) = sqrt(gamma * qa(IA_P) / q(IQ_R))
   endsubroutine conservative_to_auxiliary

   pure subroutine primitive_to_conservative(gamma, r, u, v, w, p, q)
   !< Compute the conservative variables of a cell from its primitive variables.
   real(R8P), intent(in)  :: gamma       !< Specific heats ratio.
   real(R8P), intent(in)  :: r           !< Density.
   real(R8P), intent(in)  :: u, v, w     !< Velocity components.
   real(R8P), intent(in)  :: p           !< Pressure.
   real(R8P), intent(out) :: q(NV_EULER) !< Conservative variables.
   !$acc routine seq
   !$omp declare target

   q(IQ_R)  = r
   q(IQ_RU) = r * u
   q(IQ_RV) = r * v
   q(IQ_RW) = r * w
   q(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * r * (u**2 + v**2 + w**2)
   endsubroutine primitive_to_conservative
endmodule adam_flume_euler_library
