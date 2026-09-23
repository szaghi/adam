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
use :: adam_flume_parameters, only : IA_A, IA_H, IA_P, IA_R, IA_T, IA_U, IA_V, IA_W, &
                                     IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX, NV_EULER, S_MAX
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: compute_eigenvalues
public :: compute_eigenvectors
public :: compute_face_flux_back_projection
public :: compute_face_split_fluxes
public :: compute_flux
public :: compute_roe_average
public :: conservative_to_auxiliary
public :: primitive_to_conservative

contains
   ! public procedures
   pure subroutine compute_eigenvalues(d, qa, lambda)
   !< Compute the eigenvalues of the flux Jacobian in direction `d`: `(u_n - a, u_n, u_n, u_n, u_n + a)`.
   integer(I4P), intent(in)  :: d                !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX)       !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_EULER) !< Eigenvalues.
   real(R8P)                 :: un               !< Normal velocity.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   lambda(1) = un - qa(IA_A)
   lambda(2) = un
   lambda(3) = un
   lambda(4) = un
   lambda(5) = un + qa(IA_A)
   endsubroutine compute_eigenvalues

   pure subroutine compute_eigenvectors(gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors of the flux Jacobian in direction `d` at state `qa`.
   real(R8P),    intent(in)  :: gamma                  !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                      !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX)             !< Auxiliary variables (e.g. a Roe average).
   real(R8P),    intent(out) :: el(NV_EULER,NV_EULER)  !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_EULER,NV_EULER)  !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: u(3)                   !< Velocity.
   real(R8P)                 :: a, H                   !< Speed of sound, total specific enthalpy.
   real(R8P)                 :: b1, b2                 !< Eigenvectors coefficients.
   integer(I4P)              :: d1, d2                 !< Tangential directions.
   integer(I4P)              :: c                      !< Counter.
   !$acc routine seq
   !$omp declare target

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   u  = [qa(IA_U), qa(IA_V), qa(IA_W)]
   a  = qa(IA_A)
   H  = qa(IA_H)
   b2 = (gamma - 1._R8P) / (a * a)
   b1 = 0.5_R8P * b2 * (u(1)**2 + u(2)**2 + u(3)**2)
   er = 0._R8P
   el = 0._R8P
   ! right eigenvectors (columns)
   er(1,1) = 1._R8P
   er(1,2) = 1._R8P
   er(1,5) = 1._R8P
   do c=1, 3
      er(1+c,1) = u(c)
      er(1+c,2) = u(c)
      er(1+c,5) = u(c)
   enddo
   er(1+d,1)  = u(d) - a
   er(1+d,5)  = u(d) + a
   er(1+d1,3) = 1._R8P
   er(1+d2,4) = 1._R8P
   er(5,1)    = H - a * u(d)
   er(5,2)    = 0.5_R8P * (u(1)**2 + u(2)**2 + u(3)**2)
   er(5,3)    = u(d1)
   er(5,4)    = u(d2)
   er(5,5)    = H + a * u(d)
   ! left eigenvectors (rows)
   el(1,1) = 0.5_R8P * (b1 + u(d) / a)
   el(2,1) = 1._R8P - b1
   el(3,1) = -u(d1)
   el(4,1) = -u(d2)
   el(5,1) = 0.5_R8P * (b1 - u(d) / a)
   do c=1, 3
      el(1,1+c) = -0.5_R8P * b2 * u(c)
      el(2,1+c) = b2 * u(c)
      el(5,1+c) = -0.5_R8P * b2 * u(c)
   enddo
   el(1,1+d)  = el(1,1+d) - 0.5_R8P / a
   el(5,1+d)  = el(5,1+d) + 0.5_R8P / a
   el(3,1+d1) = 1._R8P
   el(4,1+d2) = 1._R8P
   el(1,5)    = 0.5_R8P * b2
   el(2,5)    = -b2
   el(5,5)    = 0.5_R8P * b2
   endsubroutine compute_eigenvectors

   pure subroutine compute_face_flux_back_projection(is_characteristic, er, vr, flux)
   !< Return the face flux in conservative variables from the reconstructed split fields, `F = R (f+ + f-)`.
   logical,   intent(in)  :: is_characteristic     !< Reconstruction variables: characteristic or conservative.
   real(R8P), intent(in)  :: er(NV_EULER,NV_EULER) !< Right eigenvectors, er(:,k) = r_k (unused if conservative).
   real(R8P), intent(in)  :: vr(2,NV_EULER)        !< Reconstructed split fields at the face.
   real(R8P), intent(out) :: flux(NV_EULER)        !< Face flux.
   integer(I4P)           :: k, v                  !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_EULER
         flux(v) = 0._R8P
         do k=1, NV_EULER
            flux(v) = flux(v) + er(v,k) * (vr(1,k) + vr(2,k))
         enddo
      enddo
   else
      do v=1, NV_EULER
         flux(v) = vr(1,v) + vr(2,v)
      enddo
   endif
   endsubroutine compute_face_flux_back_projection

   pure subroutine compute_face_split_fluxes(gamma, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2`, ready for the WENO upwind reconstruction.
   !<
   !< The stencil holds cells `m = 1-S ... S` relative to cell `i` (cell `i+1` is `m = 1`). Characteristic variant:
   !< Roe eigenvectors of cells 0 and 1, per-wave speeds `alpha_k = max_m |lambda_k(m)|`, `w = L q`, `g = L f_d(q)`;
   !< conservative variant: no projection, one speed `alpha = max_m (|u_n| + a)`. Split `f+ = (g + alpha w) / 2`,
   !< `f- = g - f+`, stored in the layout of the WENO upwind primitive: `fsplit(2,m,k)` = f+ of cell `m` for
   !< `m = 1-S ... S-1`, `fsplit(1,m,k)` = f- of cell `m+1`. `er` returns the right eigenvectors for the back-projection.
   real(R8P),    intent(in)  :: gamma                                !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                                    !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                    !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                    !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_EULER,1-S_MAX:S_MAX)           !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX,1-S_MAX:S_MAX)            !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_EULER)   !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_EULER,NV_EULER)                !< Right eigenvectors (identity if conservative).
   real(R8P)                 :: el(NV_EULER,NV_EULER)                !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: roe(NV_AUX)                          !< Roe average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_EULER)                     !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_EULER)                      !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_EULER)                          !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                             !< Projected state, projected flux, split flux.
   integer(I4P)              :: k, m, v                              !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      call compute_roe_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), roe=roe)
      call compute_eigenvectors(gamma=gamma, d=d, qa=roe, el=el, er=er)
      alpha = 0._R8P
      do m=1-S, S
         call compute_eigenvalues(d=d, qa=qas(:,m), lambda=lambda)
         do k=1, NV_EULER
            alpha(k) = max(alpha(k), abs(lambda(k)))
         enddo
      enddo
   else
      el = 0._R8P
      er = 0._R8P
      do k=1, NV_EULER
         el(k,k) = 1._R8P
         er(k,k) = 1._R8P
      enddo
      alpha = 0._R8P
      do m=1-S, S
         alpha(1) = max(alpha(1), abs(qas(IA_U+d-1,m)) + qas(IA_A,m))
      enddo
      alpha = alpha(1)
   endif
   do m=1-S, S
      call compute_flux(d=d, q=qs(:,m), qa=qas(:,m), f=f)
      do k=1, NV_EULER
         w = 0._R8P
         g = 0._R8P
         do v=1, NV_EULER
            w = w + el(k,v) * qs(v,m)
            g = g + el(k,v) * f(v)
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
      enddo
   enddo
   endsubroutine compute_face_split_fluxes

   pure subroutine compute_flux(d, q, qa, f)
   !< Compute the physical flux in direction `d`: `rho u_d (1, u, v, w, H) + p (0, e_d, 0)`.
   integer(I4P), intent(in)  :: d           !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_EULER) !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX)  !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_EULER) !< Physical flux.
   real(R8P)                 :: un          !< Normal velocity.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   f(IQ_R)  = q(IQ_R)  * un
   f(IQ_RU) = q(IQ_RU) * un
   f(IQ_RV) = q(IQ_RV) * un
   f(IQ_RW) = q(IQ_RW) * un
   f(IQ_RU+d-1) = f(IQ_RU+d-1) + qa(IA_P)
   f(IQ_RE) = (q(IQ_RE) + qa(IA_P)) * un
   endsubroutine compute_flux

   pure subroutine compute_roe_average(gamma, qaL, qaR, roe)
   !< Compute the Roe average of two states: `sqrt(rho)`-weighted velocity and total enthalpy.
   !<
   !< The average holds density, velocity, total enthalpy, speed of sound and the matching pressure
   !< `p = rho a^2 / gamma`; temperature is not defined for the Roe state and is set to zero.
   real(R8P), intent(in)  :: gamma       !< Specific heats ratio.
   real(R8P), intent(in)  :: qaL(NV_AUX) !< Left state auxiliary variables.
   real(R8P), intent(in)  :: qaR(NV_AUX) !< Right state auxiliary variables.
   real(R8P), intent(out) :: roe(NV_AUX) !< Roe average.
   real(R8P)              :: wL, wR      !< Square roots of the densities.
   real(R8P)              :: cL, cR      !< Normalized weights.
   !$acc routine seq
   !$omp declare target

   wL = sqrt(qaL(IA_R))
   wR = sqrt(qaR(IA_R))
   cL = wL / (wL + wR)
   cR = wR / (wL + wR)
   roe(IA_R) = wL * wR
   roe(IA_U) = cL * qaL(IA_U) + cR * qaR(IA_U)
   roe(IA_V) = cL * qaL(IA_V) + cR * qaR(IA_V)
   roe(IA_W) = cL * qaL(IA_W) + cR * qaR(IA_W)
   roe(IA_H) = cL * qaL(IA_H) + cR * qaR(IA_H)
   roe(IA_A) = sqrt((gamma - 1._R8P) * (roe(IA_H) - 0.5_R8P * (roe(IA_U)**2 + roe(IA_V)**2 + roe(IA_W)**2)))
   roe(IA_P) = roe(IA_R) * roe(IA_A)**2 / gamma
   roe(IA_T) = 0._R8P
   endsubroutine compute_roe_average

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
