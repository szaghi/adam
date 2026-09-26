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
!< Two model variants (issue #41, D-1/D-14): without divergence control (`mhd_*`, `NV_MHD`) and with mixed GLM cleaning
!< (`mhd_glm_*`, `NV_MHD_GLM`, the cleaning speed `ch` an explicit argument). Their eigensystems share the 7x7
!< Roe-Balsara core of Stone et al. (2008, ApJS 178, appendix B) evaluated at the arithmetic face average, a physical
!< state (so the Roe-matrix factors are `X = 0`, `Y = 1`): the `B_n` row is decoupled (speed 0) without cleaning, the
!< `(B_n, psi)` pair is the 2x2 block of speeds `-+c_h` with GLM (block-diagonal eigenvectors). The face split and the
!< back-projection follow the Euler library (`compute_face_split_fluxes`), each variant at its own compile-time size.

! FLUME modules
use :: adam_flume_parameters, only : IA_A, IA_BX, IA_BY, IA_BZ, IA_H, IA_P, IA_R, IA_T, IA_U, IA_V, IA_W,  &
                                     IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX_MHD, &
                                     NV_MHD, NV_MHD_GLM, S_MAX
! third party modules
use :: penf,                  only : I4P, R8P

implicit none
private
public :: EPS_BT
public :: EPS_FS
public :: mhd_conservative_to_auxiliary
public :: mhd_eigenvalues
public :: mhd_eigenvectors
public :: mhd_face_average
public :: mhd_face_flux_back_projection
public :: mhd_face_split_fluxes
public :: mhd_fast_speed
public :: mhd_flux
public :: mhd_glm_eigenvalues
public :: mhd_glm_eigenvectors
public :: mhd_glm_face_flux_back_projection
public :: mhd_glm_face_split_fluxes
public :: mhd_glm_flux
public :: mhd_primitive_to_conservative
public :: mhd_sum3

real(R8P), parameter :: EPS_BT=1.e-12_R8P !< Degenerate transverse field: |B_t| <= EPS_BT max(|B|, sqrt(rho) a).
real(R8P), parameter :: EPS_FS=1.e-12_R8P !< Degenerate fast-slow pair (triple umbilic): c_f^2 - c_s^2 <= EPS_FS c_f^2.

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

   pb = 0.5_R8P * mhd_sum3(q(IQ_BX)**2, q(IQ_BY)**2, q(IQ_BZ)**2)
   qa(IA_R)  = q(IQ_R)
   qa(IA_U)  = q(IQ_RU) / q(IQ_R)
   qa(IA_V)  = q(IQ_RV) / q(IQ_R)
   qa(IA_W)  = q(IQ_RW) / q(IQ_R)
   qa(IA_P)  = (gamma - 1._R8P) * (q(IQ_RE) - 0.5_R8P * q(IQ_R) * mhd_sum3(qa(IA_U)**2, qa(IA_V)**2, qa(IA_W)**2) - pb)
   qa(IA_T)  = qa(IA_P) / (q(IQ_R) * R)
   qa(IA_H)  = (q(IQ_RE) + qa(IA_P) + pb) / q(IQ_R)
   qa(IA_A)  = sqrt(gamma * qa(IA_P) / q(IQ_R))
   qa(IA_BX) = q(IQ_BX)
   qa(IA_BY) = q(IQ_BY)
   qa(IA_BZ) = q(IQ_BZ)
   endsubroutine mhd_conservative_to_auxiliary

   pure subroutine mhd_eigenvalues(d, qa, lambda)
   !< Compute the eigenvalues of the MHD system without divergence control in direction `d`:
   !< `(u_n - c_f, u_n - c_a, u_n - c_s, u_n, u_n + c_s, u_n + c_a, u_n + c_f, 0)`, the last one of the decoupled `B_n`.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_MHD) !< Eigenvalues.
   !$acc routine seq
   !$omp declare target

   call mhd_wave_eigenvalues(d=d, qa=qa, lambda=lambda)
   lambda(8) = 0._R8P
   endsubroutine mhd_eigenvalues

   pure subroutine mhd_eigenvectors(gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors of the MHD system without divergence control in
   !< direction `d` at state `qa`: the 7x7 Roe-Balsara core on `(rho, rho u_n, rho u_t1, rho u_t2, E, B_t1, B_t2)` and
   !< the decoupled `B_n` (identity row and column).
   real(R8P),    intent(in)  :: gamma             !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                 !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)    !< Auxiliary variables (e.g. a face average).
   real(R8P),    intent(out) :: el(NV_MHD,NV_MHD) !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD) !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: l7(7,7)           !< Core left eigenvectors.
   real(R8P)                 :: r7(7,7)           !< Core right eigenvectors.
   integer(I4P)              :: mp(7)             !< Core variables in the full state.
   integer(I4P)              :: i, k              !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_eigenvectors_core(gamma=gamma, d=d, qa=qa, l7=l7, r7=r7, mp=mp)
   er = 0._R8P
   el = 0._R8P
   do k=1, 7
      do i=1, 7
         er(mp(i),k) = r7(i,k)
         el(k,mp(i)) = l7(k,i)
      enddo
   enddo
   er(IQ_BX+d-1,8) = 1._R8P
   el(8,IQ_BX+d-1) = 1._R8P
   endsubroutine mhd_eigenvectors

   pure subroutine mhd_face_average(gamma, qaL, qaR, avg)
   !< Compute the arithmetic average of two states for the face eigenvectors (issue #41, section 3.3): density,
   !< velocity, pressure and magnetic field averaged, speed of sound and total enthalpy recomputed from them, so the
   !< average is a physical state and its eigensystem is exact (no Roe-matrix `X`, `Y` factors).
   real(R8P), intent(in)  :: gamma           !< Specific heats ratio.
   real(R8P), intent(in)  :: qaL(NV_AUX_MHD) !< Left state auxiliary variables.
   real(R8P), intent(in)  :: qaR(NV_AUX_MHD) !< Right state auxiliary variables.
   real(R8P), intent(out) :: avg(NV_AUX_MHD) !< Average state auxiliary variables.
   !$acc routine seq
   !$omp declare target

   avg(IA_R)  = 0.5_R8P * (qaL(IA_R)  + qaR(IA_R))
   avg(IA_U)  = 0.5_R8P * (qaL(IA_U)  + qaR(IA_U))
   avg(IA_V)  = 0.5_R8P * (qaL(IA_V)  + qaR(IA_V))
   avg(IA_W)  = 0.5_R8P * (qaL(IA_W)  + qaR(IA_W))
   avg(IA_P)  = 0.5_R8P * (qaL(IA_P)  + qaR(IA_P))
   avg(IA_T)  = 0.5_R8P * (qaL(IA_T)  + qaR(IA_T))
   avg(IA_BX) = 0.5_R8P * (qaL(IA_BX) + qaR(IA_BX))
   avg(IA_BY) = 0.5_R8P * (qaL(IA_BY) + qaR(IA_BY))
   avg(IA_BZ) = 0.5_R8P * (qaL(IA_BZ) + qaR(IA_BZ))
   avg(IA_A)  = sqrt(gamma * avg(IA_P) / avg(IA_R))
   avg(IA_H)  = gamma / (gamma - 1._R8P) * avg(IA_P) / avg(IA_R)                        + &
                0.5_R8P * mhd_sum3(avg(IA_U)**2, avg(IA_V)**2, avg(IA_W)**2)            + &
                mhd_sum3(avg(IA_BX)**2, avg(IA_BY)**2, avg(IA_BZ)**2) / avg(IA_R)
   endsubroutine mhd_face_average

   pure subroutine mhd_face_flux_back_projection(is_characteristic, er, vr, flux)
   !< Return the face flux in conservative variables from the reconstructed split fields, `F = R (f+ + f-)`.
   logical,   intent(in)  :: is_characteristic !< Reconstruction variables: characteristic or conservative.
   real(R8P), intent(in)  :: er(NV_MHD,NV_MHD) !< Right eigenvectors, er(:,k) = r_k (unused if conservative).
   real(R8P), intent(in)  :: vr(2,NV_MHD)      !< Reconstructed split fields at the face.
   real(R8P), intent(out) :: flux(NV_MHD)      !< Face flux.
   integer(I4P)           :: k, v              !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_MHD
         flux(v) = 0._R8P
         do k=1, NV_MHD
            flux(v) = flux(v) + er(v,k) * (vr(1,k) + vr(2,k))
         enddo
      enddo
   else
      do v=1, NV_MHD
         flux(v) = vr(1,v) + vr(2,v)
      enddo
   endif
   endsubroutine mhd_face_flux_back_projection

   pure subroutine mhd_face_split_fluxes(gamma, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2` (MHD without divergence control), ready for the WENO
   !< upwind reconstruction: the algorithm and layout of the Euler `compute_face_split_fluxes`, with the eigenvectors of
   !< the arithmetic face average. The `B_n` row has speed 0 in both variants, so its face flux is exactly zero.
   real(R8P),    intent(in)  :: gamma                            !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                                !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)    !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_MHD) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD)                !< Right eigenvectors (identity if conservative).
   real(R8P)                 :: el(NV_MHD,NV_MHD)                !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: avg(NV_AUX_MHD)                  !< Face average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_MHD)                   !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_MHD)                    !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_MHD)                        !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                         !< Projected state, projected flux, split flux.
   integer(I4P)              :: pv(NV_MHD)                       !< State in the frame order of direction `d`.
   integer(I4P)              :: k, m, v                          !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   if (is_characteristic) then
      call mhd_face_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), avg=avg)
      call mhd_eigenvectors(gamma=gamma, d=d, qa=avg, el=el, er=er)
      alpha = 0._R8P
      do m=1-S, S
         call mhd_eigenvalues(d=d, qa=qas(:,m), lambda=lambda)
         do k=1, NV_MHD
            alpha(k) = max(alpha(k), abs(lambda(k)))
         enddo
      enddo
   else
      el = 0._R8P
      er = 0._R8P
      do k=1, NV_MHD
         el(k,k) = 1._R8P
         er(k,k) = 1._R8P
      enddo
      alpha = 0._R8P
      do m=1-S, S
         alpha(1) = max(alpha(1), abs(qas(IA_U+d-1,m)) + mhd_fast_speed(d=d, qa=qas(:,m)))
      enddo
      alpha = alpha(1)
      alpha(IQ_BX+d-1) = 0._R8P
   endif
   do m=1-S, S
      call mhd_flux(d=d, q=qs(:,m), qa=qas(:,m), f=f)
      do k=1, NV_MHD
         w = 0._R8P
         g = 0._R8P
         do v=1, NV_MHD
            w = w + el(k,pv(v)) * qs(pv(v),m)
            g = g + el(k,pv(v)) * f(pv(v))
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
      enddo
   enddo
   endsubroutine mhd_face_split_fluxes

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
   b2  = mhd_sum3(qa(IA_BX)**2, qa(IA_BY)**2, qa(IA_BZ)**2) / qa(IA_R)
   bt2 = max(b2 - qa(IA_BX+d-1)**2 / qa(IA_R), 0._R8P)
   cf  = sqrt(0.5_R8P * (a2 + b2 + sqrt((a2 - b2)**2 + 4._R8P * a2 * bt2)))
   endfunction mhd_fast_speed

   pure subroutine mhd_flux(d, q, qa, f)
   !< Compute the physical flux of the MHD system without divergence control in direction `d` (issue #41, section 3.1):
   !< `(rho u_n, rho u u_n + p_t n - B B_n, (E + p_t) u_n - (u.B) B_n, B u_n - u B_n)`, the `B_n` component zero
!< (induction `dB/dt = curl(u x B) = -div(u B - B u)`: the flux of `B_j` along `n` is `u_n B_j - B_n u_j`).
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_MHD)      !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_MHD)      !< Physical flux.
   real(R8P)                 :: un, bn         !< Normal velocity and magnetic field.
   real(R8P)                 :: pt             !< Total pressure, p + |B|^2 / 2.
   real(R8P)                 :: ub             !< u.B.
   integer(I4P)              :: c              !< Counter.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   bn = qa(IA_BX+d-1)
   pt = qa(IA_P) + 0.5_R8P * mhd_sum3(qa(IA_BX)**2, qa(IA_BY)**2, qa(IA_BZ)**2)
   ub = mhd_sum3(qa(IA_U) * qa(IA_BX), qa(IA_V) * qa(IA_BY), qa(IA_W) * qa(IA_BZ))
   f(IQ_R) = q(IQ_R) * un
   do c=1, 3
      f(IQ_RU+c-1) = q(IQ_RU+c-1) * un - bn * qa(IA_BX+c-1)
      f(IQ_BX+c-1) = qa(IA_BX+c-1) * un - qa(IA_U+c-1) * bn
   enddo
   f(IQ_RU+d-1) = f(IQ_RU+d-1) + pt
   f(IQ_BX+d-1) = 0._R8P
   f(IQ_RE) = (q(IQ_RE) + pt) * un - ub * bn
   endsubroutine mhd_flux

   pure subroutine mhd_glm_eigenvalues(ch, d, qa, lambda)
   !< Compute the eigenvalues of the MHD system with GLM cleaning in direction `d`: the 7 MHD waves and the
   !< `(B_n, psi)` pair `-c_h, +c_h`.
   real(R8P),    intent(in)  :: ch                 !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d                  !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)     !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_MHD_GLM) !< Eigenvalues.
   !$acc routine seq
   !$omp declare target

   call mhd_wave_eigenvalues(d=d, qa=qa, lambda=lambda)
   lambda(8) = -ch
   lambda(9) =  ch
   endsubroutine mhd_glm_eigenvalues

   pure subroutine mhd_glm_eigenvectors(ch, gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors of the MHD system with GLM cleaning in direction `d`:
   !< block-diagonal, the 7x7 Roe-Balsara core and the `(B_n, psi)` pair of the flux `(psi, c_h^2 B_n)`,
   !< `r = (1, -+c_h)`, `l = (1, -+1/c_h) / 2` (issue #41, section 3.3).
   real(R8P),    intent(in)  :: ch                        !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                     !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                         !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)            !< Auxiliary variables (e.g. a face average).
   real(R8P),    intent(out) :: el(NV_MHD_GLM,NV_MHD_GLM) !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM) !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: l7(7,7)                   !< Core left eigenvectors.
   real(R8P)                 :: r7(7,7)                   !< Core right eigenvectors.
   integer(I4P)              :: mp(7)                     !< Core variables in the full state.
   integer(I4P)              :: i, k                      !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_eigenvectors_core(gamma=gamma, d=d, qa=qa, l7=l7, r7=r7, mp=mp)
   er = 0._R8P
   el = 0._R8P
   do k=1, 7
      do i=1, 7
         er(mp(i),k) = r7(i,k)
         el(k,mp(i)) = l7(k,i)
      enddo
   enddo
   er(IQ_BX+d-1,8) = 1._R8P
   er(IQ_PSI,8)    = -ch
   er(IQ_BX+d-1,9) = 1._R8P
   er(IQ_PSI,9)    =  ch
   el(8,IQ_BX+d-1) =  0.5_R8P
   el(8,IQ_PSI)    = -0.5_R8P / ch
   el(9,IQ_BX+d-1) =  0.5_R8P
   el(9,IQ_PSI)    =  0.5_R8P / ch
   endsubroutine mhd_glm_eigenvectors

   pure subroutine mhd_glm_face_flux_back_projection(is_characteristic, er, vr, flux)
   !< Return the face flux in conservative variables from the reconstructed split fields, `F = R (f+ + f-)` (GLM).
   logical,   intent(in)  :: is_characteristic         !< Reconstruction variables: characteristic or conservative.
   real(R8P), intent(in)  :: er(NV_MHD_GLM,NV_MHD_GLM) !< Right eigenvectors, er(:,k) = r_k (unused if conservative).
   real(R8P), intent(in)  :: vr(2,NV_MHD_GLM)          !< Reconstructed split fields at the face.
   real(R8P), intent(out) :: flux(NV_MHD_GLM)          !< Face flux.
   integer(I4P)           :: k, v                      !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_MHD_GLM
         flux(v) = 0._R8P
         do k=1, NV_MHD_GLM
            flux(v) = flux(v) + er(v,k) * (vr(1,k) + vr(2,k))
         enddo
      enddo
   else
      do v=1, NV_MHD_GLM
         flux(v) = vr(1,v) + vr(2,v)
      enddo
   endif
   endsubroutine mhd_glm_face_flux_back_projection

   pure subroutine mhd_glm_face_split_fluxes(ch, gamma, d, S, is_characteristic, qs, qas, fsplit, er)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2` (MHD with GLM cleaning): as `mhd_face_split_fluxes`,
   !< with the `(B_n, psi)` pair at speeds `-+c_h` (conservative variant: one speed, `max(|u_n| + c_f, c_h)`).
   real(R8P),    intent(in)  :: ch                                   !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                                !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                                    !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                    !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                    !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD_GLM,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)        !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_MHD_GLM) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM)            !< Right eigenvectors (identity if conservative).
   real(R8P)                 :: el(NV_MHD_GLM,NV_MHD_GLM)            !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: avg(NV_AUX_MHD)                      !< Face average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_MHD_GLM)                   !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_MHD_GLM)                    !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_MHD_GLM)                        !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                             !< Projected state, projected flux, split flux.
   integer(I4P)              :: pv(NV_MHD_GLM)                       !< State in the frame order of direction `d`.
   integer(I4P)              :: k, m, v                              !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   pv(NV_MHD_GLM) = IQ_PSI
   if (is_characteristic) then
      call mhd_face_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), avg=avg)
      call mhd_glm_eigenvectors(ch=ch, gamma=gamma, d=d, qa=avg, el=el, er=er)
      alpha = 0._R8P
      do m=1-S, S
         call mhd_glm_eigenvalues(ch=ch, d=d, qa=qas(:,m), lambda=lambda)
         do k=1, NV_MHD_GLM
            alpha(k) = max(alpha(k), abs(lambda(k)))
         enddo
      enddo
   else
      el = 0._R8P
      er = 0._R8P
      do k=1, NV_MHD_GLM
         el(k,k) = 1._R8P
         er(k,k) = 1._R8P
      enddo
      alpha = ch
      do m=1-S, S
         alpha(1) = max(alpha(1), abs(qas(IA_U+d-1,m)) + mhd_fast_speed(d=d, qa=qas(:,m)))
      enddo
      alpha = alpha(1)
   endif
   do m=1-S, S
      call mhd_glm_flux(ch=ch, d=d, q=qs(:,m), qa=qas(:,m), f=f)
      do k=1, NV_MHD_GLM
         w = 0._R8P
         g = 0._R8P
         do v=1, NV_MHD_GLM
            w = w + el(k,pv(v)) * qs(pv(v),m)
            g = g + el(k,pv(v)) * f(pv(v))
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
      enddo
   enddo
   endsubroutine mhd_glm_face_split_fluxes

   pure subroutine mhd_glm_flux(ch, d, q, qa, f)
   !< Compute the physical flux of the MHD system with GLM cleaning in direction `d`: the MHD flux with the `B_n`
   !< component `psi` and the `psi` component `c_h^2 B_n` (mixed GLM, issue #41, section 3.1).
   real(R8P),    intent(in)  :: ch              !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_MHD_GLM)   !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)  !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)   !< Physical flux.
   !$acc routine seq
   !$omp declare target

   call mhd_flux(d=d, q=q, qa=qa, f=f)
   f(IQ_BX+d-1) = q(IQ_PSI)
   f(IQ_PSI)    = ch * ch * q(IQ_BX+d-1)
   endsubroutine mhd_glm_flux

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
   q(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * r * mhd_sum3(u**2, v**2, w**2) + 0.5_R8P * mhd_sum3(bx**2, by**2, bz**2)
   q(IQ_BX) = bx
   q(IQ_BY) = by
   q(IQ_BZ) = bz
   endsubroutine mhd_primitive_to_conservative

   pure function mhd_sum3(a, b, c) result(s)
   !< Return `a + b + c` summed in increasing order, so the result does not depend on the order of the arguments.
   !<
   !< Why: `|u|^2`, `|B|^2`, `u.B` summed in the fixed x, y, z order round differently when the components are permuted,
   !< so a problem rotated from x to y or z (cyclic tangents) drifts at round-off (issue #41, MV-4). Sorting the three
   !< terms (branchless `min`/`max`) makes the sum a function of the set of terms: rotations are bitwise invariant.
   real(R8P), intent(in) :: a, b, c !< Terms.
   real(R8P)             :: s       !< Sum.
   !$acc routine seq
   !$omp declare target

   s = (min(a, b, c) + max(min(a, b), min(max(a, b), c))) + max(a, b, c)
   endfunction mhd_sum3

   ! private procedures
   pure subroutine mhd_eigenvectors_core(gamma, d, qa, l7, r7, mp)
   !< Compute the 7x7 Roe-Balsara eigenvectors of ideal MHD in direction `d` at the physical state `qa` (Stone et al.
   !< 2008, ApJS 178, appendix B, conservative variables, with the Roe-matrix factors `X = 0`, `Y = 1` of a physical
   !< state). Core variables, in the full-state indexes returned in `mp`: `(rho, rho u_n, rho u_t1, rho u_t2, E, B_t1,
   !< B_t2)`, tangents taken cyclically as in the Euler library; `B_n` enters as a parameter. Waves ordered
   !< `(u_n - c_f, u_n - c_a, u_n - c_s, u_n, u_n + c_s, u_n + c_a, u_n + c_f)`; `l7(k,:)` is the row `l_k`, `r7(:,k)`
   !< the column `r_k`, `l7 r7 = I`.
   !<
   !< Degeneracies (named thresholds, issue #41, section 3.3): a transverse field below `EPS_BT` of the field (or of
   !< `sqrt(rho) a`) takes the direction `(beta_t1, beta_t2) = (1, 0)`; a fast-slow separation `c_f^2 - c_s^2` below
   !< `EPS_FS c_f^2` (triple umbilic) takes `(alpha_f, alpha_s) = (1, 0)`; `alpha_f`, `alpha_s` are clamped to [0, 1].
   real(R8P),    intent(in)  :: gamma          !< Specific heats ratio.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: l7(7,7)        !< Left eigenvectors, rows.
   real(R8P),    intent(out) :: r7(7,7)        !< Right eigenvectors, columns.
   integer(I4P), intent(out) :: mp(7)          !< Full-state index of each core variable.
   real(R8P)                 :: rho, di        !< Density and its inverse.
   real(R8P)                 :: v1, v2, v3     !< Velocity: normal, tangential 1 and 2.
   real(R8P)                 :: b1, b2, b3     !< Magnetic field: normal, tangential 1 and 2.
   real(R8P)                 :: vsq, btsq, bt  !< |u|^2, |B_t|^2, |B_t|.
   real(R8P)                 :: vaxsq, vax     !< Squared normal Alfven speed and the speed.
   real(R8P)                 :: hp             !< Total enthalpy without the magnetic part, h - |B|^2 / rho.
   real(R8P)                 :: asq, a         !< Squared sound speed and the speed.
   real(R8P)                 :: ct2, tsum      !< Squared transverse Alfven speed, sum of the squared speeds.
   real(R8P)                 :: tdif, cf2_cs2  !< Difference of the squared speeds, c_f^2 - c_s^2.
   real(R8P)                 :: cfsq, cf       !< Squared fast speed and the speed.
   real(R8P)                 :: cssq, cs       !< Squared slow speed and the speed.
   real(R8P)                 :: bet2, bet3     !< Transverse field direction.
   real(R8P)                 :: vbet           !< u_t . beta.
   real(R8P)                 :: alf, als       !< Fast and slow normalisation coefficients.
   real(R8P)                 :: sqrtd, isqrtd  !< sqrt(rho) and its inverse.
   real(R8P)                 :: sgn            !< Sign of B_n (+1 if zero).
   real(R8P)                 :: qf, qs         !< c_f alpha_f sgn, c_s alpha_s sgn.
   real(R8P)                 :: afp, asp       !< a alpha_f / sqrt(rho), a alpha_s / sqrt(rho).
   real(R8P)                 :: afpbb, aspbb   !< afp |B_t|, asp |B_t|.
   real(R8P)                 :: qa_, qb, qc, qd !< Row temporaries.
   real(R8P)                 :: nrm            !< Left normalisation, 1 / (2 a^2), then (gamma-1) / (2 a^2).
   real(R8P)                 :: cff, css       !< Normalised c_f alpha_f, c_s alpha_s.
   real(R8P)                 :: qfn, qsn       !< Normalised qf, qs.
   real(R8P)                 :: afr, asr       !< Normalised afp rho, asp rho.
   real(R8P)                 :: afpb, aspb     !< Normalised afp |B_t|, asp |B_t|.
   real(R8P)                 :: afb, asb       !< alpha_f, alpha_s normalised by (gamma-1) / (2 a^2).
   real(R8P)                 :: vqstr          !< u_t . beta.
   real(R8P)                 :: nrm2           !< (gamma-1) / a^2.
   integer(I4P)              :: d1, d2         !< Tangential directions.
   !$acc routine seq
   !$omp declare target

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   mp(1) = IQ_R
   mp(2) = IQ_RU + d  - 1
   mp(3) = IQ_RU + d1 - 1
   mp(4) = IQ_RU + d2 - 1
   mp(5) = IQ_RE
   mp(6) = IQ_BX + d1 - 1
   mp(7) = IQ_BX + d2 - 1
   rho = qa(IA_R)
   di  = 1._R8P / rho
   v1  = qa(IA_U  + d  - 1)
   v2  = qa(IA_U  + d1 - 1)
   v3  = qa(IA_U  + d2 - 1)
   b1  = qa(IA_BX + d  - 1)
   b2  = qa(IA_BX + d1 - 1)
   b3  = qa(IA_BX + d2 - 1)
   vsq   = v1 * v1 + v2 * v2 + v3 * v3
   btsq  = b2 * b2 + b3 * b3
   vaxsq = b1 * b1 * di
   hp    = qa(IA_H) - (vaxsq + btsq * di)
   asq   = max((gamma - 1._R8P) * (hp - 0.5_R8P * vsq), tiny(1._R8P))
   a     = sqrt(asq)
   ! fast and slow speeds
   ct2     = btsq * di
   tsum    = vaxsq + ct2 + asq
   tdif    = vaxsq + ct2 - asq
   cf2_cs2 = sqrt(tdif * tdif + 4._R8P * asq * ct2)
   cfsq    = 0.5_R8P * (tsum + cf2_cs2)
   cf      = sqrt(cfsq)
   cssq    = asq * vaxsq / cfsq
   cs      = sqrt(cssq)
   ! transverse field direction
   bt = sqrt(btsq)
   if (bt <= EPS_BT * max(sqrt(b1 * b1 + btsq), sqrt(rho) * a)) then
      bet2 = 1._R8P
      bet3 = 0._R8P
   else
      bet2 = b2 / bt
      bet3 = b3 / bt
   endif
   vbet = v2 * bet2 + v3 * bet3
   ! fast and slow normalisation
   if (cfsq - cssq <= EPS_FS * cfsq) then
      alf = 1._R8P
      als = 0._R8P
   elseif (asq - cssq <= 0._R8P) then
      alf = 0._R8P
      als = 1._R8P
   elseif (cfsq - asq <= 0._R8P) then
      alf = 1._R8P
      als = 0._R8P
   else
      alf = min(1._R8P, sqrt((asq - cssq) / (cfsq - cssq)))
      als = min(1._R8P, sqrt((cfsq - asq) / (cfsq - cssq)))
   endif
   sqrtd  = sqrt(rho)
   isqrtd = 1._R8P / sqrtd
   sgn    = 1._R8P ; if (b1 < 0._R8P) sgn = -1._R8P
   qf     = cf * alf * sgn
   qs     = cs * als * sgn
   afp    = a * alf * isqrtd
   asp    = a * als * isqrtd
   afpbb  = afp * bt
   aspbb  = asp * bt
   vax    = sqrt(vaxsq)
   ! right eigenvectors, columns
   r7(1,1) = alf
   r7(1,2) = 0._R8P
   r7(1,3) = als
   r7(1,4) = 1._R8P
   r7(1,5) = als
   r7(1,6) = 0._R8P
   r7(1,7) = alf
   r7(2,1) = alf * (v1 - cf)
   r7(2,2) = 0._R8P
   r7(2,3) = als * (v1 - cs)
   r7(2,4) = v1
   r7(2,5) = als * (v1 + cs)
   r7(2,6) = 0._R8P
   r7(2,7) = alf * (v1 + cf)
   qa_ = alf * v2
   qb  = als * v2
   qc  = qs * bet2
   qd  = qf * bet2
   r7(3,1) = qa_ + qc
   r7(3,2) = -bet3
   r7(3,3) = qb - qd
   r7(3,4) = v2
   r7(3,5) = qb + qd
   r7(3,6) = bet3
   r7(3,7) = qa_ - qc
   qa_ = alf * v3
   qb  = als * v3
   qc  = qs * bet3
   qd  = qf * bet3
   r7(4,1) = qa_ + qc
   r7(4,2) = bet2
   r7(4,3) = qb - qd
   r7(4,4) = v3
   r7(4,5) = qb + qd
   r7(4,6) = -bet2
   r7(4,7) = qa_ - qc
   r7(5,1) = alf * (hp - v1 * cf) + qs * vbet + aspbb
   r7(5,2) = -(v2 * bet3 - v3 * bet2)
   r7(5,3) = als * (hp - v1 * cs) - qf * vbet - afpbb
   r7(5,4) = 0.5_R8P * vsq
   r7(5,5) = als * (hp + v1 * cs) + qf * vbet - afpbb
   r7(5,6) = -r7(5,2)
   r7(5,7) = alf * (hp + v1 * cf) - qs * vbet + aspbb
   r7(6,1) = asp * bet2
   r7(6,2) = -bet3 * sgn * isqrtd
   r7(6,3) = -afp * bet2
   r7(6,4) = 0._R8P
   r7(6,5) = r7(6,3)
   r7(6,6) = r7(6,2)
   r7(6,7) = r7(6,1)
   r7(7,1) = asp * bet3
   r7(7,2) = bet2 * sgn * isqrtd
   r7(7,3) = -afp * bet3
   r7(7,4) = 0._R8P
   r7(7,5) = r7(7,3)
   r7(7,6) = r7(7,2)
   r7(7,7) = r7(7,1)
   ! left eigenvectors, rows
   nrm   = 0.5_R8P / asq
   cff   = nrm * alf * cf
   css   = nrm * als * cs
   qfn   = nrm * qf
   qsn   = nrm * qs
   afr   = nrm * afp * rho
   asr   = nrm * asp * rho
   afpb  = nrm * afp * bt
   aspb  = nrm * asp * bt
   nrm   = nrm * (gamma - 1._R8P)
   afb   = alf * nrm
   asb   = als * nrm
   vqstr = v2 * bet2 + v3 * bet3
   nrm2  = 2._R8P * nrm
   l7(1,1) = afb * (vsq - hp) + cff * (cf + v1) - qsn * vqstr - aspb
   l7(1,2) = -afb * v1 - cff
   l7(1,3) = -afb * v2 + qsn * bet2
   l7(1,4) = -afb * v3 + qsn * bet3
   l7(1,5) = afb
   l7(1,6) = asr * bet2 - afb * b2
   l7(1,7) = asr * bet3 - afb * b3
   l7(2,1) = 0.5_R8P * (v2 * bet3 - v3 * bet2)
   l7(2,2) = 0._R8P
   l7(2,3) = -0.5_R8P * bet3
   l7(2,4) = 0.5_R8P * bet2
   l7(2,5) = 0._R8P
   l7(2,6) = -0.5_R8P * sqrtd * bet3 * sgn
   l7(2,7) = 0.5_R8P * sqrtd * bet2 * sgn
   l7(3,1) = asb * (vsq - hp) + css * (cs + v1) + qfn * vqstr + afpb
   l7(3,2) = -asb * v1 - css
   l7(3,3) = -asb * v2 - qfn * bet2
   l7(3,4) = -asb * v3 - qfn * bet3
   l7(3,5) = asb
   l7(3,6) = -afr * bet2 - asb * b2
   l7(3,7) = -afr * bet3 - asb * b3
   l7(4,1) = 1._R8P - nrm2 * 0.5_R8P * vsq
   l7(4,2) = nrm2 * v1
   l7(4,3) = nrm2 * v2
   l7(4,4) = nrm2 * v3
   l7(4,5) = -nrm2
   l7(4,6) = nrm2 * b2
   l7(4,7) = nrm2 * b3
   l7(5,1) = asb * (vsq - hp) + css * (cs - v1) - qfn * vqstr + afpb
   l7(5,2) = -asb * v1 + css
   l7(5,3) = -asb * v2 + qfn * bet2
   l7(5,4) = -asb * v3 + qfn * bet3
   l7(5,5) = asb
   l7(5,6) = l7(3,6)
   l7(5,7) = l7(3,7)
   l7(6,1) = -l7(2,1)
   l7(6,2) = 0._R8P
   l7(6,3) = -l7(2,3)
   l7(6,4) = -l7(2,4)
   l7(6,5) = 0._R8P
   l7(6,6) = l7(2,6)
   l7(6,7) = l7(2,7)
   l7(7,1) = afb * (vsq - hp) + cff * (cf - v1) + qsn * vqstr - aspb
   l7(7,2) = -afb * v1 + cff
   l7(7,3) = -afb * v2 - qsn * bet2
   l7(7,4) = -afb * v3 - qsn * bet3
   l7(7,5) = afb
   l7(7,6) = l7(1,6)
   l7(7,7) = l7(1,7)
   endsubroutine mhd_eigenvectors_core

   pure subroutine mhd_frame_indexes(d, pv)
   !< Return the conservative variables in the frame order of direction `d`, `(rho, rho u_n, rho u_t1, rho u_t2, E, B_n,
   !< B_t1, B_t2)`, tangents cyclic as in `mhd_eigenvectors_core`: the projections sum in this order, so a problem
   !< rotated from x to y or z sums the same terms in the same order (bitwise invariant, issue #41, MV-4).
   integer(I4P), intent(in)  :: d          !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(out) :: pv(NV_MHD) !< Full-state index of each frame variable.
   integer(I4P)              :: d1, d2     !< Tangential directions.
   !$acc routine seq
   !$omp declare target

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   pv(1) = IQ_R
   pv(2) = IQ_RU + d  - 1
   pv(3) = IQ_RU + d1 - 1
   pv(4) = IQ_RU + d2 - 1
   pv(5) = IQ_RE
   pv(6) = IQ_BX + d  - 1
   pv(7) = IQ_BX + d1 - 1
   pv(8) = IQ_BX + d2 - 1
   endsubroutine mhd_frame_indexes

   pure subroutine mhd_wave_eigenvalues(d, qa, lambda)
   !< Compute the 7 MHD wave speeds in direction `d`, `(u_n - c_f, u_n - c_a, u_n - c_s, u_n, u_n + c_s, u_n + c_a,
   !< u_n + c_f)`, the fast speed as `mhd_fast_speed`, the slow one as `a^2 b_n^2 / (rho c_f^2)` (no cancellation).
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(7)      !< Wave speeds.
   real(R8P)                 :: un             !< Normal velocity.
   real(R8P)                 :: cf, ca, cs     !< Fast, Alfven and slow speeds.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   cf = mhd_fast_speed(d=d, qa=qa)
   ca = abs(qa(IA_BX+d-1)) / sqrt(qa(IA_R))
   cs = sqrt(qa(IA_A)**2 * ca**2 / max(cf**2, tiny(1._R8P)))
   lambda(1) = un - cf
   lambda(2) = un - ca
   lambda(3) = un - cs
   lambda(4) = un
   lambda(5) = un + cs
   lambda(6) = un + ca
   lambda(7) = un + cf
   endsubroutine mhd_wave_eigenvalues
endmodule adam_flume_mhd_library
