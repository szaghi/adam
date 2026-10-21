   ! Body of the isothermal and barotropic MHD libraries (issue #100, M13 P5): adam_flume_mhd_iso_library and
   ! adam_flume_mhd_baro_library contain it and differ only in their EOS hooks (`eos_pressure`, `eos_temperature`,
   ! `eos_sound_speed` of the argument named `gamma`: c0^2, or the id of the barotropic table).
   !
   ! The public interface repeats, routine by routine, the one adam_flume_mhd_library and adam_flume_mhd_riemann_library
   ! give the MHD kernels and objects (the model-agnostic kernel bodies call either). The system is ideal MHD without
   ! an energy equation (Mignone 2007, eqs. 1-6), the total pressure `P(rho) + |B|^2 / 2`; the energy slot is inert
   ! (D-M13-1): zero flux, zero state, its eigenvector the unit vector at speed 0. The normal field is a parameter of
   ! the 1-D core (decoupled at speed 0, or the GLM pair at -+c_h), as in the adiabatic library.
   !
   ! **Eigensystem**: the six waves of Kim et al. (1999, ApJ 514, 506, eqs. 2-9..2-53) on `(rho, rho u_n, rho u_t1,
   ! rho u_t2, B_t1, B_t2)`, with `a^2 = dP/drho` (the Jacobian holds the gas law only there, so the isothermal vectors
   ! are the barotropic ones). Kim's `alpha_f = sqrt(c_f^2 - b_n^2) / sqrt(c_f^2 - c_s^2)` is the Ryu-Jones `alpha_f`
   ! times `c_f / a`, `alpha_s` the same: both come from the cancellation-free pair of the adiabatic core (issue #79),
   ! with its degenerate values (`beta = (1, 0)` at `|B_t| = 0`, `(1, 0)` for the pair at the triple umbilic) and the
   ! continuity sign of Kim's eq. 2-53. eos/mhd_oracle.py checks this form (L R = I, R Lambda L = dF/dU) at random and
   ! degenerate states, both gases.
   !
   ! **HLLD** (Mignone 2007, J. Comput. Phys. 225, 1427, section 3): three states between the outer speeds (those of
   ! the adiabatic library, Miyoshi & Kusano 2005 eq. 67) and the rotational speeds `S*_{L,R} = u* -+ |B_n| /
   ! sqrt(rho*)`; `rho*`, `m_n*` and the first two flux components the HLL averages, `u* = F_rho^hll / rho^hll`, the
   ! transverse components from eqs. 30-37, the flux by eq. 38. `(S - S*_L)(S - S*_R) -> 0` keeps the outer transverse
   ! state (Mignone's degenerate case). HLL is used, and `fallback` set, when `rho* <= 0` or the speeds are out of order
   ! (`S_L < S*_L <= S*_R < S_R`). An isolated rotational discontinuity is resolved exactly, and with `B = 0` the flux is
   ! the isothermal HLLC's (eos/mhd_oracle.py).

   ! public procedures
   pure subroutine mhd_conservative_to_auxiliary(gamma, R, q, qa)
   !< Compute the auxiliary variables of a cell: `qa(IA_P) = P(rho)`, `qa(IA_A) = sqrt(dP/drho)`, `qa(IA_H)` the total
   !< pressure over the density (no enthalpy without an energy equation), the field copied into `qa(IA_BX:IA_BZ)`.
   real(R8P), intent(in)  :: gamma          !< EOS parameter.
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
   qa(IA_P)  = eos_pressure(gamma, q(IQ_R))
   qa(IA_T)  = eos_temperature(gamma, R, q(IQ_R))
   qa(IA_H)  = (qa(IA_P) + pb) / q(IQ_R)
   qa(IA_A)  = eos_sound_speed(gamma, q(IQ_R))
   qa(IA_BX) = q(IQ_BX)
   qa(IA_BY) = q(IQ_BY)
   qa(IA_BZ) = q(IQ_BZ)
   endsubroutine mhd_conservative_to_auxiliary

   pure subroutine mhd_eigenvalues(d, qa, lambda)
   !< Compute the eigenvalues in direction `d`: `(u_n - c_f, u_n - c_a, u_n - c_s, u_n + c_s, u_n + c_a, u_n + c_f, 0, 0)`,
   !< the last two of the inert energy slot and of the decoupled `B_n`.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_MHD) !< Eigenvalues.
   !$acc routine seq
   !$omp declare target

   call mhd_wave_eigenvalues(d=d, qa=qa, lambda=lambda)
   lambda(7) = 0._R8P
   lambda(8) = 0._R8P
   endsubroutine mhd_eigenvalues

   pure subroutine mhd_eigenvectors(gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors in direction `d` at state `qa`: the 6x6 core of Kim et al.
   !< on `(rho, rho u_n, rho u_t1, rho u_t2, B_t1, B_t2)`, the inert energy slot and the decoupled `B_n` (unit rows).
   real(R8P),    intent(in)  :: gamma             !< EOS parameter (unused: the core reads `qa(IA_A)`).
   integer(I4P), intent(in)  :: d                 !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)    !< Auxiliary variables (e.g. a face average).
   real(R8P),    intent(out) :: el(NV_MHD,NV_MHD) !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD) !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: l6(6,6)           !< Core left eigenvectors.
   real(R8P)                 :: r6(6,6)           !< Core right eigenvectors.
   integer(I4P)              :: mp(6)             !< Core variables in the full state.
   integer(I4P)              :: i, k              !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_eigenvectors_core(d=d, qa=qa, l6=l6, r6=r6, mp=mp)
   er = 0._R8P
   el = 0._R8P
   do k=1, 6
      do i=1, 6
         er(mp(i),k) = r6(i,k)
         el(k,mp(i)) = l6(k,i)
      enddo
   enddo
   er(IQ_RE,7) = 1._R8P
   el(7,IQ_RE) = 1._R8P
   er(IQ_BX+d-1,8) = 1._R8P
   el(8,IQ_BX+d-1) = 1._R8P
   endsubroutine mhd_eigenvectors

   pure subroutine mhd_face_average(gamma, qaL, qaR, avg)
   !< Compute the arithmetic average of two states for the face eigenvectors: density, velocity and field averaged, the
   !< pressure and the speed of sound of the averaged density (a physical state, so its eigensystem is exact).
   real(R8P), intent(in)  :: gamma           !< EOS parameter.
   real(R8P), intent(in)  :: qaL(NV_AUX_MHD) !< Left state auxiliary variables.
   real(R8P), intent(in)  :: qaR(NV_AUX_MHD) !< Right state auxiliary variables.
   real(R8P), intent(out) :: avg(NV_AUX_MHD) !< Average state auxiliary variables.
   !$acc routine seq
   !$omp declare target

   avg(IA_R)  = 0.5_R8P * (qaL(IA_R)  + qaR(IA_R))
   avg(IA_U)  = 0.5_R8P * (qaL(IA_U)  + qaR(IA_U))
   avg(IA_V)  = 0.5_R8P * (qaL(IA_V)  + qaR(IA_V))
   avg(IA_W)  = 0.5_R8P * (qaL(IA_W)  + qaR(IA_W))
   avg(IA_T)  = 0.5_R8P * (qaL(IA_T)  + qaR(IA_T))
   avg(IA_BX) = 0.5_R8P * (qaL(IA_BX) + qaR(IA_BX))
   avg(IA_BY) = 0.5_R8P * (qaL(IA_BY) + qaR(IA_BY))
   avg(IA_BZ) = 0.5_R8P * (qaL(IA_BZ) + qaR(IA_BZ))
   avg(IA_P)  = eos_pressure(gamma, avg(IA_R))
   avg(IA_A)  = eos_sound_speed(gamma, avg(IA_R))
   avg(IA_H)  = (avg(IA_P) + 0.5_R8P * mhd_sum3(avg(IA_BX)**2, avg(IA_BY)**2, avg(IA_BZ)**2)) / avg(IA_R)
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

   pure subroutine mhd_face_split_fluxes(gamma, d, S, is_characteristic, qs, qas, fsplit, er, mu)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2` (no divergence control), as the adiabatic library:
   !< the eigenvectors of the arithmetic face average; the energy slot and the `B_n` row have speed 0 and zero flux.
   real(R8P),    intent(in)  :: gamma                            !< EOS parameter.
   integer(I4P), intent(in)  :: d                                !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)    !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_MHD) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD)                !< Right eigenvectors (identity if conservative).
   real(R8P),    intent(out) :: mu(NV_MHD)                       !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_MHD,NV_MHD)                !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: avg(NV_AUX_MHD)                  !< Face average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_MHD)                   !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_MHD)                    !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_MHD)                        !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                         !< Projected state, projected flux, split flux.
   real(R8P)                 :: a                                !< Absolute projection of one cell (descaler).
   real(R8P)                 :: qm(NV_MHD)                       !< Reference magnitudes of one cell (descaler).
   real(R8P)                 :: amax                             !< Fast speed bound of the stencil (descaler).
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
   amax = 0._R8P
   do m=1-S, S
      amax = max(amax, abs(qas(IA_U+d-1,m)) + mhd_fast_speed(d=d, qa=qas(:,m)))
   enddo
   mu = 0._R8P
   do m=1-S, S
      call mhd_flux(d=d, q=qs(:,m), qa=qas(:,m), f=f)
      call mhd_state_magnitudes(q=qs(:,m), qa=qas(:,m), qm=qm)
      do k=1, NV_MHD
         w = 0._R8P
         g = 0._R8P
         a = 0._R8P
         do v=1, NV_MHD
            w = w + el(k,pv(v)) * qs(pv(v),m)
            g = g + el(k,pv(v)) * f(pv(v))
            a = a + abs(el(k,pv(v))) * (abs(f(pv(v))) + amax * qm(pv(v)))
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
         mu(k) = mu(k) + 0.5_R8P * a
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine mhd_face_split_fluxes

   pure function mhd_fast_speed(d, qa) result(cf)
   !< Return the fast magnetosonic speed of a cell along direction `d` (as the adiabatic library, `a = qa(IA_A)`).
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
   !< Compute the physical flux in direction `d`: `(rho u_n, rho u u_n + p_t n - B B_n, 0, B u_n - u B_n)`, `p_t = P(rho)
   !< + |B|^2 / 2`, the energy slot and the `B_n` component zero.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_MHD)      !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_MHD)      !< Physical flux.
   real(R8P)                 :: un, bn         !< Normal velocity and magnetic field.
   real(R8P)                 :: pt             !< Total pressure, P + |B|^2 / 2.
   integer(I4P)              :: c              !< Counter.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   bn = qa(IA_BX+d-1)
   pt = qa(IA_P) + 0.5_R8P * mhd_sum3(qa(IA_BX)**2, qa(IA_BY)**2, qa(IA_BZ)**2)
   f(IQ_R) = q(IQ_R) * un
   do c=1, 3
      f(IQ_RU+c-1) = q(IQ_RU+c-1) * un - bn * qa(IA_BX+c-1)
      f(IQ_BX+c-1) = qa(IA_BX+c-1) * un - qa(IA_U+c-1) * bn
   enddo
   f(IQ_RU+d-1) = f(IQ_RU+d-1) + pt
   f(IQ_BX+d-1) = 0._R8P
   f(IQ_RE) = 0._R8P
   endsubroutine mhd_flux

   pure subroutine mhd_frame_indexes(d, pv)
   !< Return the conservative variables in the frame order of direction `d`, `(rho, rho u_n, rho u_t1, rho u_t2, slot5,
   !< B_n, B_t1, B_t2)`, tangents cyclic (bitwise rotation invariance, issue #41, MV-4).
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

   pure subroutine mhd_glm_eigenvalues(ch, d, qa, lambda)
   !< Compute the eigenvalues with GLM cleaning in direction `d`: the six waves, the energy slot (0) and the `(B_n, psi)`
   !< pair `-c_h, +c_h`.
   real(R8P),    intent(in)  :: ch                 !< GLM cleaning speed.
   integer(I4P), intent(in)  :: d                  !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)     !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_MHD_GLM) !< Eigenvalues.
   !$acc routine seq
   !$omp declare target

   call mhd_wave_eigenvalues(d=d, qa=qa, lambda=lambda)
   lambda(7) = 0._R8P
   lambda(8) = -ch
   lambda(9) =  ch
   endsubroutine mhd_glm_eigenvalues

   pure subroutine mhd_glm_eigenvectors(ch, gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors with GLM cleaning in direction `d`: block-diagonal, the
   !< core, the energy slot and the `(B_n, psi)` pair of the flux `(psi, c_h^2 B_n)`, `r = (1, -+c_h)`,
   !< `l = (1, -+1/c_h) / 2`.
   real(R8P),    intent(in)  :: ch                        !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                     !< EOS parameter (unused: the core reads `qa(IA_A)`).
   integer(I4P), intent(in)  :: d                         !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD)            !< Auxiliary variables (e.g. a face average).
   real(R8P),    intent(out) :: el(NV_MHD_GLM,NV_MHD_GLM) !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM) !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: l6(6,6)                   !< Core left eigenvectors.
   real(R8P)                 :: r6(6,6)                   !< Core right eigenvectors.
   integer(I4P)              :: mp(6)                     !< Core variables in the full state.
   integer(I4P)              :: i, k                      !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_eigenvectors_core(d=d, qa=qa, l6=l6, r6=r6, mp=mp)
   er = 0._R8P
   el = 0._R8P
   do k=1, 6
      do i=1, 6
         er(mp(i),k) = r6(i,k)
         el(k,mp(i)) = l6(k,i)
      enddo
   enddo
   er(IQ_RE,7)     = 1._R8P
   el(7,IQ_RE)     = 1._R8P
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

   pure subroutine mhd_glm_face_split_fluxes(ch, gamma, d, S, is_characteristic, qs, qas, fsplit, er, mu)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2` (GLM cleaning): as `mhd_face_split_fluxes`, with the
   !< `(B_n, psi)` pair at speeds `-+c_h` (conservative variant: one speed, `max(|u_n| + c_f, c_h)`).
   real(R8P),    intent(in)  :: ch                                   !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                                !< EOS parameter.
   integer(I4P), intent(in)  :: d                                    !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                    !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                    !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD_GLM,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)        !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_MHD_GLM) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM)            !< Right eigenvectors (identity if conservative).
   real(R8P),    intent(out) :: mu(NV_MHD_GLM)                       !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_MHD_GLM,NV_MHD_GLM)            !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: avg(NV_AUX_MHD)                      !< Face average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_MHD_GLM)                   !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_MHD_GLM)                    !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_MHD_GLM)                        !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                             !< Projected state, projected flux, split flux.
   real(R8P)                 :: a                                    !< Absolute projection of one cell (descaler).
   real(R8P)                 :: qm(NV_MHD_GLM)                       !< Reference magnitudes of one cell (descaler).
   real(R8P)                 :: amax                                 !< Fast speed bound of the stencil (descaler).
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
   amax = 0._R8P
   do m=1-S, S
      amax = max(amax, abs(qas(IA_U+d-1,m)) + mhd_fast_speed(d=d, qa=qas(:,m)))
   enddo
   mu = 0._R8P
   do m=1-S, S
      call mhd_glm_flux(ch=ch, d=d, q=qs(:,m), qa=qas(:,m), f=f)
      call mhd_state_magnitudes(q=qs(1:NV_MHD,m), qa=qas(:,m), qm=qm(1:NV_MHD))
      qm(IQ_PSI) = abs(qs(IQ_PSI,m)) + ch * qm(IQ_BX)
      do k=1, NV_MHD_GLM
         w = 0._R8P
         g = 0._R8P
         a = 0._R8P
         do v=1, NV_MHD_GLM
            w = w + el(k,pv(v)) * qs(pv(v),m)
            g = g + el(k,pv(v)) * f(pv(v))
            a = a + abs(el(k,pv(v))) * (abs(f(pv(v))) + amax * qm(pv(v)))
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
         mu(k) = mu(k) + 0.5_R8P * a
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine mhd_glm_face_split_fluxes

   pure subroutine mhd_glm_flux(ch, d, q, qa, f)
   !< Compute the physical flux with GLM cleaning in direction `d`: the flux with the `B_n` component `psi` and the
   !< `psi` component `c_h^2 B_n` (mixed GLM).
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
   !< Compute the conservative variables of a cell from its primitive variables: the energy slot is inert (0), the
   !< pressure `p` is that of the gas law and is not read (`psi`, if any, is set by the caller).
   real(R8P), intent(in)  :: gamma      !< EOS parameter (unused).
   real(R8P), intent(in)  :: r          !< Density.
   real(R8P), intent(in)  :: u, v, w    !< Velocity components.
   real(R8P), intent(in)  :: p          !< Pressure (unused: P(rho)).
   real(R8P), intent(in)  :: bx, by, bz !< Magnetic field components.
   real(R8P), intent(out) :: q(NV_MHD)  !< Conservative variables.
   !$acc routine seq
   !$omp declare target

   q(IQ_R)  = r
   q(IQ_RU) = r * u
   q(IQ_RV) = r * v
   q(IQ_RW) = r * w
   q(IQ_RE) = 0._R8P
   q(IQ_BX) = bx
   q(IQ_BY) = by
   q(IQ_BZ) = bz
   endsubroutine mhd_primitive_to_conservative

   pure function mhd_sum3(a, b, c) result(s)
   !< Return `a + b + c` summed in increasing order (rotation-invariant sums, issue #41, MV-4).
   real(R8P), intent(in) :: a, b, c !< Terms.
   real(R8P)             :: s       !< Sum.
   !$acc routine seq
   !$omp declare target

   s = (min(a, b, c) + max(min(a, b), min(max(a, b), c))) + max(a, b, c)
   endfunction mhd_sum3

   pure subroutine mhd_state_magnitudes(q, qa, qm)
   !< Return the reference magnitude of each conservative variable of a cell (the WENO descaler scales, as the adiabatic
   !< library: momentum `|rho u|_1 + sqrt(rho p)`, field `|B|_1 + sqrt(p)`, `p = P(rho)`; the inert slot its own value).
   real(R8P), intent(in)  :: q(NV_MHD)      !< Conservative variables.
   real(R8P), intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P), intent(out) :: qm(NV_MHD)     !< Reference magnitudes.
   real(R8P)              :: pm, bm         !< Momentum and field magnitudes.
   !$acc routine seq
   !$omp declare target

   pm = mhd_sum3(abs(q(IQ_RU)), abs(q(IQ_RV)), abs(q(IQ_RW))) + sqrt(abs(qa(IA_R) * qa(IA_P)))
   bm = mhd_sum3(abs(q(IQ_BX)), abs(q(IQ_BY)), abs(q(IQ_BZ))) + sqrt(abs(qa(IA_P)))
   qm(IQ_R)  = abs(q(IQ_R))
   qm(IQ_RU) = pm
   qm(IQ_RV) = pm
   qm(IQ_RW) = pm
   qm(IQ_RE) = abs(q(IQ_RE))
   qm(IQ_BX) = bm
   qm(IQ_BY) = bm
   qm(IQ_BZ) = bm
   endsubroutine mhd_state_magnitudes

   ! public procedures, Riemann solvers and face states of `weno-riemann`
   pure subroutine mhd_backbone_flux(gamma, d, qL, qR, f)
   !< Compute the positivity backbone flux of two states in direction `d` (no divergence control): the Lax-Friedrichs
   !< flux with the Wu speed, as the adiabatic library.
   real(R8P),    intent(in)  :: gamma           !< EOS parameter.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD)      !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD)      !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD)       !< Flux.
   real(R8P)                 :: qaL(NV_AUX_MHD) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_MHD) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_MHD)      !< Left physical flux.
   real(R8P)                 :: fR(NV_MHD)      !< Right physical flux.
   real(R8P)                 :: sigma           !< Lax-Friedrichs speed.
   integer(I4P)              :: v               !< Counter.
   !$acc routine seq
   !$omp declare target

   call mhd_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qL, qa=qaL)
   call mhd_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qR, qa=qaR)
   call mhd_flux(d=d, q=qL, qa=qaL, f=fL)
   call mhd_flux(d=d, q=qR, qa=qaR, f=fR)
   sigma = wu_speed(d=d, qaL=qaL, qaR=qaR)
   do v=1, NV_MHD
      f(v) = 0.5_R8P * (fL(v) + fR(v)) - 0.5_R8P * sigma * (qR(v) - qL(v))
   enddo
   endsubroutine mhd_backbone_flux

   pure subroutine mhd_face_interpolation_fields(gamma, d, S, is_characteristic, qs, qas, fint, er, mu)
   !< Compute the fields of the stencil of face `i+1/2` to interpolate (no divergence control), in the WENO upwind
   !< layout: characteristic, or primitive `(rho, u, v, w, 0, B_x, B_y, B_z)` (the inert slot in place of `p`).
   real(R8P),    intent(in)  :: gamma                          !< EOS parameter.
   integer(I4P), intent(in)  :: d                              !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                              !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic              !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)  !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_MHD) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD)              !< Right eigenvectors (unused if primitive).
   real(R8P),    intent(out) :: mu(NV_MHD)                     !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_MHD,NV_MHD)              !< Left eigenvectors.
   real(R8P)                 :: avg(NV_AUX_MHD)                !< Face average of cells 0 and 1.
   real(R8P)                 :: w(NV_MHD)                      !< Fields of one cell.
   real(R8P)                 :: qm(NV_MHD)                     !< Reference magnitudes of one cell (descaler).
   integer(I4P)              :: pv(NV_MHD)                     !< State in the frame order of direction `d`.
   integer(I4P)              :: k, m, v                        !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   if (is_characteristic) then
      call mhd_face_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), avg=avg)
      call mhd_eigenvectors(gamma=gamma, d=d, qa=avg, el=el, er=er)
   else
      er = 0._R8P
   endif
   mu = 0._R8P
   do m=1-S, S
      if (is_characteristic) then
         call mhd_state_magnitudes(q=qs(:,m), qa=qas(:,m), qm=qm)
         do k=1, NV_MHD
            w(k) = 0._R8P
            do v=1, NV_MHD
               w(k) = w(k) + el(k,pv(v)) * qs(pv(v),m)
               mu(k) = mu(k) + abs(el(k,pv(v))) * qm(pv(v))
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = 0._R8P
         w(6) = qas(IA_BX,m)
         w(7) = qas(IA_BY,m)
         w(8) = qas(IA_BZ,m)
      endif
      do k=1, NV_MHD
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
         if (.not.is_characteristic) mu(k) = mu(k) + abs(w(k))
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine mhd_face_interpolation_fields

   pure subroutine mhd_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Compute the two face states from the interpolated fields (no divergence control): the left state from `vr(2,:)`,
   !< the right one from `vr(1,:)`; a state with non-positive density falls back to its cell.
   real(R8P), intent(in)  :: gamma             !< EOS parameter.
   logical,   intent(in)  :: is_characteristic !< Characteristic (or primitive) variables.
   real(R8P), intent(in)  :: er(NV_MHD,NV_MHD) !< Right eigenvectors (unused if primitive).
   real(R8P), intent(in)  :: vr(2,NV_MHD)      !< Interpolated fields.
   real(R8P), intent(in)  :: q0(NV_MHD)        !< Conservative variables of cell 0.
   real(R8P), intent(in)  :: q1(NV_MHD)        !< Conservative variables of cell 1.
   real(R8P), intent(out) :: qL(NV_MHD)        !< Left state.
   real(R8P), intent(out) :: qR(NV_MHD)        !< Right state.
   integer(I4P)           :: k, v              !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_MHD
         qL(v) = 0._R8P
         qR(v) = 0._R8P
         do k=1, NV_MHD
            qL(v) = qL(v) + er(v,k) * vr(2,k)
            qR(v) = qR(v) + er(v,k) * vr(1,k)
         enddo
      enddo
   else
      call mhd_primitive_to_conservative(gamma=gamma, r=vr(2,1), u=vr(2,2), v=vr(2,3), w=vr(2,4), p=vr(2,5), &
                                         bx=vr(2,6), by=vr(2,7), bz=vr(2,8), q=qL)
      call mhd_primitive_to_conservative(gamma=gamma, r=vr(1,1), u=vr(1,2), v=vr(1,3), w=vr(1,4), p=vr(1,5), &
                                         bx=vr(1,6), by=vr(1,7), bz=vr(1,8), q=qR)
   endif
   if (.not.is_admissible(q=qL)) qL = q0
   if (.not.is_admissible(q=qR)) qR = q1
   endsubroutine mhd_face_states

   pure subroutine mhd_glm_face_interpolation_fields(ch, gamma, d, S, is_characteristic, qs, qas, fint, er, mu)
   !< Compute the fields of the stencil of face `i+1/2` to interpolate (GLM): as `mhd_face_interpolation_fields`, `psi`
   !< the ninth field (primitive) or the GLM eigenvectors (characteristic).
   real(R8P),    intent(in)  :: ch                                 !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                              !< EOS parameter.
   integer(I4P), intent(in)  :: d                                  !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                  !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                  !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD_GLM,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)      !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_MHD_GLM) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM)          !< Right eigenvectors (unused if primitive).
   real(R8P),    intent(out) :: mu(NV_MHD_GLM)                     !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_MHD_GLM,NV_MHD_GLM)          !< Left eigenvectors.
   real(R8P)                 :: avg(NV_AUX_MHD)                    !< Face average of cells 0 and 1.
   real(R8P)                 :: w(NV_MHD_GLM)                      !< Fields of one cell.
   real(R8P)                 :: qm(NV_MHD_GLM)                     !< Reference magnitudes of one cell (descaler).
   integer(I4P)              :: pv(NV_MHD_GLM)                     !< State in the frame order of direction `d`.
   integer(I4P)              :: k, m, v                            !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   pv(NV_MHD_GLM) = IQ_PSI
   if (is_characteristic) then
      call mhd_face_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), avg=avg)
      call mhd_glm_eigenvectors(ch=ch, gamma=gamma, d=d, qa=avg, el=el, er=er)
   else
      er = 0._R8P
   endif
   mu = 0._R8P
   do m=1-S, S
      if (is_characteristic) then
         call mhd_state_magnitudes(q=qs(1:NV_MHD,m), qa=qas(:,m), qm=qm(1:NV_MHD))
         qm(IQ_PSI) = abs(qs(IQ_PSI,m)) + ch * qm(IQ_BX)
         do k=1, NV_MHD_GLM
            w(k) = 0._R8P
            do v=1, NV_MHD_GLM
               w(k) = w(k) + el(k,pv(v)) * qs(pv(v),m)
               mu(k) = mu(k) + abs(el(k,pv(v))) * qm(pv(v))
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = 0._R8P
         w(6) = qas(IA_BX,m)
         w(7) = qas(IA_BY,m)
         w(8) = qas(IA_BZ,m)
         w(9) = qs(IQ_PSI,m)
      endif
      do k=1, NV_MHD_GLM
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
         if (.not.is_characteristic) mu(k) = mu(k) + abs(w(k))
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine mhd_glm_face_interpolation_fields

   pure subroutine mhd_glm_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Compute the two face states from the interpolated fields (GLM): as `mhd_face_states`, `psi` included.
   real(R8P), intent(in)  :: gamma                     !< EOS parameter.
   logical,   intent(in)  :: is_characteristic         !< Characteristic (or primitive) variables.
   real(R8P), intent(in)  :: er(NV_MHD_GLM,NV_MHD_GLM) !< Right eigenvectors (unused if primitive).
   real(R8P), intent(in)  :: vr(2,NV_MHD_GLM)          !< Interpolated fields.
   real(R8P), intent(in)  :: q0(NV_MHD_GLM)            !< Conservative variables of cell 0.
   real(R8P), intent(in)  :: q1(NV_MHD_GLM)            !< Conservative variables of cell 1.
   real(R8P), intent(out) :: qL(NV_MHD_GLM)            !< Left state.
   real(R8P), intent(out) :: qR(NV_MHD_GLM)            !< Right state.
   integer(I4P)           :: k, v                      !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_MHD_GLM
         qL(v) = 0._R8P
         qR(v) = 0._R8P
         do k=1, NV_MHD_GLM
            qL(v) = qL(v) + er(v,k) * vr(2,k)
            qR(v) = qR(v) + er(v,k) * vr(1,k)
         enddo
      enddo
   else
      call mhd_primitive_to_conservative(gamma=gamma, r=vr(2,1), u=vr(2,2), v=vr(2,3), w=vr(2,4), p=vr(2,5), &
                                         bx=vr(2,6), by=vr(2,7), bz=vr(2,8), q=qL)
      call mhd_primitive_to_conservative(gamma=gamma, r=vr(1,1), u=vr(1,2), v=vr(1,3), w=vr(1,4), p=vr(1,5), &
                                         bx=vr(1,6), by=vr(1,7), bz=vr(1,8), q=qR)
      qL(IQ_PSI) = vr(2,9)
      qR(IQ_PSI) = vr(1,9)
   endif
   if (.not.is_admissible(q=qL)) qL = q0
   if (.not.is_admissible(q=qR)) qR = q1
   endsubroutine mhd_glm_face_states

   pure subroutine mhd_glm_riemann_hll(ch, gamma, d, qL, qR, f)
   !< Compute the HLL flux of two states in direction `d` (GLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< EOS parameter.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_GLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_GLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)  !< Flux.
   real(R8P)                 :: wL(NV_MHD)     !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)     !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)     !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)     !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)     !< Frame flux.
   real(R8P)                 :: bn, psi        !< Face normal field and GLM scalar.
   integer(I4P)              :: pv(NV_MHD)     !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call glm_face_subsystem(ch=ch, bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_hll(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   call scatter_glm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, f=f)
   endsubroutine mhd_glm_riemann_hll

   pure subroutine mhd_glm_riemann_hlld(ch, gamma, d, qL, qR, f, fallback)
   !< Compute the HLLD flux of two states in direction `d` (GLM, the `(B_n, psi)` subsystem solved exactly); `fallback`
   !< is set when HLL replaced HLLD.
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< EOS parameter.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_GLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_GLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)  !< Flux.
   logical,      intent(out) :: fallback       !< HLL used instead of HLLD.
   real(R8P)                 :: wL(NV_MHD)     !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)     !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)     !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)     !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)     !< Frame flux.
   real(R8P)                 :: bn, psi        !< Face normal field and GLM scalar.
   integer(I4P)              :: pv(NV_MHD)     !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call glm_face_subsystem(ch=ch, bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_hlld(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr, fallback=fallback)
   call scatter_glm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, f=f)
   endsubroutine mhd_glm_riemann_hlld

   pure subroutine mhd_glm_riemann_llf(ch, gamma, d, qL, qR, f)
   !< Compute the LLF flux of two states in direction `d` (GLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< EOS parameter.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_GLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_GLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)  !< Flux.
   real(R8P)                 :: wL(NV_MHD)     !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)     !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)     !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)     !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)     !< Frame flux.
   real(R8P)                 :: bn, psi        !< Face normal field and GLM scalar.
   integer(I4P)              :: pv(NV_MHD)     !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call glm_face_subsystem(ch=ch, bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_llf(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   call scatter_glm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, f=f)
   endsubroutine mhd_glm_riemann_llf

   pure subroutine mhd_riemann_hll(gamma, d, qL, qR, f)
   !< Compute the HLL flux of two states in direction `d` (no divergence control, each state's own `B_n`).
   real(R8P),    intent(in)  :: gamma      !< EOS parameter.
   integer(I4P), intent(in)  :: d          !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD)  !< Flux.
   real(R8P)                 :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD) !< Frame flux.
   integer(I4P)              :: pv(NV_MHD) !< State in the frame order of direction `d`.
   integer(I4P)              :: k          !< Counter.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=qL(pv(6)), w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=qR(pv(6)), w=wR, u=uR)
   call frame_hll(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   do k=1, NV_MHD
      f(pv(k)) = fr(k)
   enddo
   endsubroutine mhd_riemann_hll

   pure subroutine mhd_riemann_hlld(gamma, d, qL, qR, f, fallback)
   !< Compute the HLLD flux of two states in direction `d` (no divergence control, `B_n` the average of the two states,
   !< its flux zero); `fallback` is set when HLL replaced HLLD.
   real(R8P),    intent(in)  :: gamma      !< EOS parameter.
   integer(I4P), intent(in)  :: d          !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD)  !< Flux.
   logical,      intent(out) :: fallback   !< HLL used instead of HLLD.
   real(R8P)                 :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD) !< Frame flux.
   real(R8P)                 :: bn         !< Face normal field.
   integer(I4P)              :: pv(NV_MHD) !< State in the frame order of direction `d`.
   integer(I4P)              :: k          !< Counter.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   bn = 0.5_R8P * (qL(pv(6)) + qR(pv(6)))
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_hlld(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr, fallback=fallback)
   do k=1, NV_MHD
      f(pv(k)) = fr(k)
   enddo
   endsubroutine mhd_riemann_hlld

   pure subroutine mhd_riemann_llf(gamma, d, qL, qR, f)
   !< Compute the LLF flux of two states in direction `d` (no divergence control, each state's own `B_n`).
   real(R8P),    intent(in)  :: gamma      !< EOS parameter.
   integer(I4P), intent(in)  :: d          !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD)  !< Flux.
   real(R8P)                 :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD) !< Frame flux.
   integer(I4P)              :: pv(NV_MHD) !< State in the frame order of direction `d`.
   integer(I4P)              :: k          !< Counter.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call frame_state(gamma=gamma, pv=pv, q=qL, bn=qL(pv(6)), w=wL, u=uL)
   call frame_state(gamma=gamma, pv=pv, q=qR, bn=qR(pv(6)), w=wR, u=uR)
   call frame_llf(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   do k=1, NV_MHD
      f(pv(k)) = fr(k)
   enddo
   endsubroutine mhd_riemann_llf

   ! private procedures
   pure subroutine mhd_eigenvectors_core(d, qa, l6, r6, mp)
   !< Compute the 6x6 eigenvectors of Kim et al. (1999) in direction `d` at the physical state `qa` (`a = qa(IA_A)`),
   !< core variables `(rho, rho u_n, rho u_t1, rho u_t2, B_t1, B_t2)` in the full-state indexes `mp`, waves ordered
   !< `(u_n - c_f, u_n - c_a, u_n - c_s, u_n + c_s, u_n + c_a, u_n + c_f)`; `l6(k,:)` the row `l_k`, `r6(:,k)` the column
   !< `r_k`, `l6 r6 = I`. Kim's `b = B / sqrt(rho)`; `alpha_f = alf c_f / a`, `alpha_s = als` with the cancellation-free
   !< Ryu-Jones pair `(alf, als)` and the degenerate values of the adiabatic core (see the module); Kim's eq. 2-53 flips
   !< the slow pair (`a^2 > c_a^2`) or the fast pair (`a^2 < c_a^2`) with the sign of the transverse field.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: l6(6,6)        !< Left eigenvectors, rows.
   real(R8P),    intent(out) :: r6(6,6)        !< Right eigenvectors, columns.
   integer(I4P), intent(out) :: mp(6)          !< Full-state index of each core variable.
   real(R8P)                 :: rho, sqrtd     !< Density and its square root.
   real(R8P)                 :: v1, v2, v3     !< Velocity: normal, tangential 1 and 2.
   real(R8P)                 :: b1, b2, b3     !< Magnetic field: normal, tangential 1 and 2.
   real(R8P)                 :: bxs            !< Kim's b_x = B_n / sqrt(rho), signed.
   real(R8P)                 :: btsq, bt       !< |B_t|^2, |B_t|.
   real(R8P)                 :: asq, a         !< Squared sound speed and the speed.
   real(R8P)                 :: vaxsq, ca      !< Squared normal Alfven speed and the speed.
   real(R8P)                 :: ct2, tsum      !< Squared transverse Alfven speed, sum of the squared speeds.
   real(R8P)                 :: tdif, cf2_cs2  !< Difference of the squared speeds, c_f^2 - c_s^2.
   real(R8P)                 :: asq_m_cssq     !< a^2 - c_s^2, cancellation-free.
   real(R8P)                 :: cfsq_m_asq     !< c_f^2 - a^2, cancellation-free.
   real(R8P)                 :: cfsq, cf       !< Squared fast speed and the speed.
   real(R8P)                 :: cs             !< Slow speed.
   real(R8P)                 :: bet2, bet3     !< Transverse field direction.
   real(R8P)                 :: alf, als       !< Ryu-Jones normalisation pair.
   real(R8P)                 :: af, as         !< Kim's alpha_f, alpha_s.
   real(R8P)                 :: sbx, sbt       !< sign(B_n) (+1 if zero), Kim's sign(b_T).
   real(R8P)                 :: th1, th2       !< Kim's theta_1, theta_2.
   real(R8P)                 :: vt             !< beta . u_t.
   real(R8P)                 :: sg             !< -1, +1 for the two members of a pair.
   integer(I4P)              :: d1, d2         !< Tangential directions.
   integer(I4P)              :: k, c           !< Pair member and column.
   !$acc routine seq
   !$omp declare target

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   mp(1) = IQ_R
   mp(2) = IQ_RU + d  - 1
   mp(3) = IQ_RU + d1 - 1
   mp(4) = IQ_RU + d2 - 1
   mp(5) = IQ_BX + d1 - 1
   mp(6) = IQ_BX + d2 - 1
   rho   = qa(IA_R)
   sqrtd = sqrt(rho)
   v1    = qa(IA_U  + d  - 1)
   v2    = qa(IA_U  + d1 - 1)
   v3    = qa(IA_U  + d2 - 1)
   b1    = qa(IA_BX + d  - 1)
   b2    = qa(IA_BX + d1 - 1)
   b3    = qa(IA_BX + d2 - 1)
   a     = qa(IA_A)
   asq   = a * a
   bxs   = b1 / sqrtd
   btsq  = b2 * b2 + b3 * b3
   vaxsq = b1 * b1 / rho
   ca    = sqrt(vaxsq)
   ! fast and slow speeds
   ct2     = btsq / rho
   tsum    = vaxsq + ct2 + asq
   tdif    = vaxsq + ct2 - asq
   cf2_cs2 = sqrt(tdif * tdif + 4._R8P * asq * ct2)
   cfsq    = 0.5_R8P * (tsum + cf2_cs2)
   cf      = sqrt(cfsq)
   cs      = sqrt(asq * vaxsq / cfsq)
   ! transverse field direction
   bt = sqrt(btsq)
   if (bt <= EPS_BT * max(sqrt(b1 * b1 + btsq), sqrtd * a)) then
      bet2 = 1._R8P
      bet3 = 0._R8P
   else
      bet2 = b2 / bt
      bet3 = b3 / bt
   endif
   ! the Ryu-Jones pair without cancellation (as the adiabatic core), then Kim's alpha_f = alf c_f / a
   if (cf2_cs2 <= EPS_FS * cfsq) then
      alf = 1._R8P
      als = 0._R8P
   else
      if (tdif >= 0._R8P) then
         asq_m_cssq = 2._R8P * asq * ct2 / (cf2_cs2 + tdif)
         cfsq_m_asq = cf2_cs2 - asq_m_cssq
      else
         cfsq_m_asq = 2._R8P * asq * ct2 / (cf2_cs2 - tdif)
         asq_m_cssq = cf2_cs2 - cfsq_m_asq
      endif
      alf = min(1._R8P, sqrt(asq_m_cssq / cf2_cs2))
      als = min(1._R8P, sqrt(cfsq_m_asq / cf2_cs2))
   endif
   af  = alf * cf / a
   as  = als
   sbx = 1._R8P ; if (b1 < 0._R8P) sbx = -1._R8P
   sbt = 1._R8P ; if (b2 < 0._R8P .or. (b2 == 0._R8P .and. b3 < 0._R8P)) sbt = -1._R8P
   th1 = 2._R8P * (af * af * asq + as * as * cfsq)
   th2 = 2._R8P * (af * af * cf * a + as * as * ca * cs)
   vt  = bet2 * v2 + bet3 * v3
   do k=1, 2
      sg = real(2 * k - 3, R8P)
      ! fast pair, Kim eqs. 2-21 and 2-31..2-36: columns 1 (u - c_f) and 6 (u + c_f)
      c = merge(1_I4P, 6_I4P, k == 1)
      r6(1,c) = af
      r6(2,c) = af * (v1 + sg * cf)
      r6(3,c) = af * v2 - sg * as * bet2 * bxs
      r6(4,c) = af * v3 - sg * as * bet3 * bxs
      r6(5,c) = as * bet2 * cf / sqrtd
      r6(6,c) = as * bet3 * cf / sqrtd
      l6(c,1) = af * asq / th1 + sg * (-af * a * v1 + as * cs * vt * sbx) / th2
      l6(c,2) = sg * af * a / th2
      l6(c,3) = -sg * as * bet2 * cs * sbx / th2
      l6(c,4) = -sg * as * bet3 * cs * sbx / th2
      l6(c,5) = as * bet2 * cf * sqrtd / th1
      l6(c,6) = as * bet3 * cf * sqrtd / th1
      ! Alfven pair, eqs. 2-22 and 2-38..2-43: columns 2 (u - c_a) and 5 (u + c_a)
      c = merge(2_I4P, 5_I4P, k == 1)
      r6(1,c) = 0._R8P
      r6(2,c) = 0._R8P
      r6(3,c) = -sg * bet3 * sbx
      r6(4,c) =  sg * bet2 * sbx
      r6(5,c) =  bet3 / sqrtd
      r6(6,c) = -bet2 / sqrtd
      l6(c,1) = sg * 0.5_R8P * (bet3 * v2 - bet2 * v3) * sbx
      l6(c,2) = 0._R8P
      l6(c,3) = -sg * 0.5_R8P * bet3 * sbx
      l6(c,4) =  sg * 0.5_R8P * bet2 * sbx
      l6(c,5) =  0.5_R8P * bet3 * sqrtd
      l6(c,6) = -0.5_R8P * bet2 * sqrtd
      ! slow pair, eqs. 2-23 and 2-45..2-50: columns 3 (u - c_s) and 4 (u + c_s)
      c = merge(3_I4P, 4_I4P, k == 1)
      r6(1,c) = as
      r6(2,c) = as * (v1 + sg * cs)
      r6(3,c) = as * v2 + sg * af * bet2 * a * sbx
      r6(4,c) = as * v3 + sg * af * bet3 * a * sbx
      r6(5,c) = -af * bet2 * asq / (cf * sqrtd)
      r6(6,c) = -af * bet3 * asq / (cf * sqrtd)
      l6(c,1) = as * cfsq / th1 - sg * (as * ca * v1 + af * cf * vt * sbx) / th2
      l6(c,2) = sg * as * ca / th2
      l6(c,3) = sg * af * bet2 * cf * sbx / th2
      l6(c,4) = sg * af * bet3 * cf * sbx / th2
      l6(c,5) = -af * bet2 * cf * sqrtd / th1
      l6(c,6) = -af * bet3 * cf * sqrtd / th1
   enddo
   ! Kim eq. 2-53: continuity across a sign change of the transverse field
   if (asq > vaxsq) then
      r6(:,3:4) = sbt * r6(:,3:4)
      l6(3:4,:) = sbt * l6(3:4,:)
   elseif (asq < vaxsq) then
      r6(:,1) = sbt * r6(:,1)
      r6(:,6) = sbt * r6(:,6)
      l6(1,:) = sbt * l6(1,:)
      l6(6,:) = sbt * l6(6,:)
   endif
   endsubroutine mhd_eigenvectors_core

   pure subroutine mhd_wave_eigenvalues(d, qa, lambda)
   !< Compute the six wave speeds in direction `d`, `(u_n - c_f, u_n - c_a, u_n - c_s, u_n + c_s, u_n + c_a, u_n + c_f)`.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_MHD) !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(6)      !< Wave speeds.
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
   lambda(4) = un + cs
   lambda(5) = un + ca
   lambda(6) = un + cf
   endsubroutine mhd_wave_eigenvalues

   pure function frame_fast_speed(gamma, w) result(cf)
   !< Return the fast magnetosonic speed along the frame normal, `a^2 = dP/drho` of the state.
   real(R8P), intent(in) :: gamma     !< EOS parameter.
   real(R8P), intent(in) :: w(NV_MHD) !< Frame primitive variables.
   real(R8P)             :: cf        !< Fast magnetosonic speed.
   real(R8P)             :: a2        !< Squared sound speed.
   real(R8P)             :: b2        !< Squared Alfven speed, |B|^2 / rho.
   real(R8P)             :: bt2       !< Transverse part of b2.
   !$acc routine seq
   !$omp declare target

   a2  = eos_sound_speed(gamma, w(1))**2
   b2  = mhd_sum3(w(6)**2, w(7)**2, w(8)**2) / w(1)
   bt2 = max(b2 - w(6)**2 / w(1), 0._R8P)
   cf  = sqrt(0.5_R8P * (a2 + b2 + sqrt((a2 - b2)**2 + 4._R8P * a2 * bt2)))
   endfunction frame_fast_speed

   pure subroutine frame_flux(w, u, f)
   !< Compute the physical flux along the frame normal, the energy slot and the `B_n` component zero.
   real(R8P), intent(in)  :: w(NV_MHD) !< Frame primitive variables (`w(5) = P(rho)`).
   real(R8P), intent(in)  :: u(NV_MHD) !< Frame conservative variables.
   real(R8P), intent(out) :: f(NV_MHD) !< Frame flux.
   real(R8P)              :: pt        !< Total pressure.
   !$acc routine seq
   !$omp declare target

   pt = w(5) + 0.5_R8P * mhd_sum3(w(6)**2, w(7)**2, w(8)**2)
   f(1) = u(1) * w(2)
   f(2) = u(2) * w(2) - w(6) * w(6) + pt
   f(3) = u(3) * w(2) - w(6) * w(7)
   f(4) = u(4) * w(2) - w(6) * w(8)
   f(5) = 0._R8P
   f(6) = 0._R8P
   f(7) = w(7) * w(2) - w(3) * w(6)
   f(8) = w(8) * w(2) - w(4) * w(6)
   endsubroutine frame_flux

   pure subroutine frame_hll(gamma, wL, uL, wR, uR, f)
   !< Compute the HLL flux in the frame.
   real(R8P), intent(in)  :: gamma      !< EOS parameter.
   real(R8P), intent(in)  :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P), intent(in)  :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P), intent(in)  :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P), intent(in)  :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P), intent(out) :: f(NV_MHD)  !< Frame flux.
   real(R8P)              :: fL(NV_MHD) !< Left physical flux.
   real(R8P)              :: fR(NV_MHD) !< Right physical flux.
   real(R8P)              :: sL, sR     !< Outer wave speeds.
   !$acc routine seq
   !$omp declare target

   call frame_flux(w=wL, u=uL, f=fL)
   call frame_flux(w=wR, u=uR, f=fR)
   call frame_speeds(gamma=gamma, wL=wL, wR=wR, sL=sL, sR=sR)
   call hll_select(sL=sL, sR=sR, uL=uL, uR=uR, fL=fL, fR=fR, f=f)
   endsubroutine frame_hll

   pure subroutine frame_hlld(gamma, wL, uL, wR, uR, f, fallback)
   !< Compute the isothermal HLLD flux in the frame (Mignone 2007, section 3), `B_n = wL(6) = wR(6)`: the HLL averages,
   !< the advective velocity `u* = F_rho^hll / rho^hll`, the rotational speeds `u* -+ |B_n| / sqrt(rho*)`, the outer
   !< star states (eqs. 30-33) and the centre state (eqs. 34-37), sampled by eq. 38.
   real(R8P), intent(in)  :: gamma          !< EOS parameter.
   real(R8P), intent(in)  :: wL(NV_MHD)     !< Left frame primitive variables.
   real(R8P), intent(in)  :: uL(NV_MHD)     !< Left frame conservative variables.
   real(R8P), intent(in)  :: wR(NV_MHD)     !< Right frame primitive variables.
   real(R8P), intent(in)  :: uR(NV_MHD)     !< Right frame conservative variables.
   real(R8P), intent(out) :: f(NV_MHD)      !< Frame flux.
   logical,   intent(out) :: fallback       !< HLL used instead of HLLD.
   real(R8P)              :: fL(NV_MHD)     !< Left physical flux.
   real(R8P)              :: fR(NV_MHD)     !< Right physical flux.
   real(R8P)              :: usL(NV_MHD)    !< Left star state.
   real(R8P)              :: usR(NV_MHD)    !< Right star state.
   real(R8P)              :: sL, sR         !< Outer wave speeds.
   real(R8P)              :: ssL, ssR       !< Rotational wave speeds.
   real(R8P)              :: rs, ms         !< HLL density and normal momentum.
   real(R8P)              :: frs, fms       !< HLL mass and normal-momentum fluxes.
   real(R8P)              :: us             !< Advective velocity, F_rho^hll / rho^hll.
   real(R8P)              :: sq, x          !< sqrt(rho*), sqrt(rho*) sign(B_n).
   real(R8P)              :: ct(2), cb(2)   !< Centre transverse momentum and field.
   integer(I4P)           :: k              !< Counter.
   !$acc routine seq
   !$omp declare target

   call frame_flux(w=wL, u=uL, f=fL)
   call frame_flux(w=wR, u=uR, f=fR)
   call frame_speeds(gamma=gamma, wL=wL, wR=wR, sL=sL, sR=sR)
   fallback = .false.
   if (sL >= 0._R8P) then
      f = fL
      return
   elseif (sR <= 0._R8P) then
      f = fR
      return
   endif
   rs  = (sR * uR(1) - sL * uL(1) - fR(1) + fL(1)) / (sR - sL)
   ms  = (sR * uR(2) - sL * uL(2) - fR(2) + fL(2)) / (sR - sL)
   frs = (sR * fL(1) - sL * fR(1) + sL * sR * (uR(1) - uL(1))) / (sR - sL)
   fms = (sR * fL(2) - sL * fR(2) + sL * sR * (uR(2) - uL(2))) / (sR - sL)
   if (.not.(rs > 0._R8P)) then
      fallback = .true.
      call hll_select(sL=sL, sR=sR, uL=uL, uR=uR, fL=fL, fR=fR, f=f)
      return
   endif
   us  = frs / rs
   sq  = sqrt(rs)
   ssL = us - abs(wL(6)) / sq
   ssR = us + abs(wL(6)) / sq
   if (.not.(sL < ssL .and. ssR < sR)) then
      fallback = .true.
      call hll_select(sL=sL, sR=sR, uL=uL, uR=uR, fL=fL, fR=fR, f=f)
      return
   endif
   call hlld_star(s=sL, ssL=ssL, ssR=ssR, rs=rs, ms=ms, us=us, w=wL, ust=usL)
   call hlld_star(s=sR, ssL=ssL, ssR=ssR, rs=rs, ms=ms, us=us, w=wR, ust=usR)
   if (ssL >= 0._R8P) then
      do k=1, NV_MHD
         f(k) = fL(k) + sL * (usL(k) - uL(k))
      enddo
   elseif (ssR <= 0._R8P) then
      do k=1, NV_MHD
         f(k) = fR(k) + sR * (usR(k) - uR(k))
      enddo
   else
      ! the centre state (reached only with B_n /= 0: B_n = 0 gives ssL = ssR = u*)
      x = sq * sign(1._R8P, wL(6))
      do k=1, 2
         ct(k) = 0.5_R8P * (usL(2+k) + usR(2+k)) + 0.5_R8P * x * (usR(6+k) - usL(6+k))
         cb(k) = 0.5_R8P * (usL(6+k) + usR(6+k)) + 0.5_R8P * (usR(2+k) - usL(2+k)) / x
      enddo
      f(1) = frs
      f(2) = fms
      f(3) = ct(1) * us - wL(6) * cb(1)
      f(4) = ct(2) * us - wL(6) * cb(2)
      f(5) = 0._R8P
      f(6) = 0._R8P
      f(7) = cb(1) * us - wL(6) * ct(1) / rs
      f(8) = cb(2) * us - wL(6) * ct(2) / rs
   endif
   endsubroutine frame_hlld

   pure subroutine frame_llf(gamma, wL, uL, wR, uR, f)
   !< Compute the LLF (Rusanov) flux in the frame, speed `max(|u_n| + c_f)` of the two states.
   real(R8P), intent(in)  :: gamma      !< EOS parameter.
   real(R8P), intent(in)  :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P), intent(in)  :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P), intent(in)  :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P), intent(in)  :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P), intent(out) :: f(NV_MHD)  !< Frame flux.
   real(R8P)              :: fL(NV_MHD) !< Left physical flux.
   real(R8P)              :: fR(NV_MHD) !< Right physical flux.
   real(R8P)              :: alpha      !< Lax-Friedrichs speed.
   integer(I4P)           :: k          !< Counter.
   !$acc routine seq
   !$omp declare target

   call frame_flux(w=wL, u=uL, f=fL)
   call frame_flux(w=wR, u=uR, f=fR)
   alpha = max(abs(wL(2)) + frame_fast_speed(gamma=gamma, w=wL), abs(wR(2)) + frame_fast_speed(gamma=gamma, w=wR))
   do k=1, NV_MHD
      f(k) = 0.5_R8P * (fL(k) + fR(k)) - 0.5_R8P * alpha * (uR(k) - uL(k))
   enddo
   endsubroutine frame_llf

   pure subroutine frame_speeds(gamma, wL, wR, sL, sR)
   !< Compute the outer wave speeds `S_L = min(u_nL, u_nR) - max(c_fL, c_fR)`, `S_R = max(u_nL, u_nR) + max(c_fL, c_fR)`.
   real(R8P), intent(in)  :: gamma      !< EOS parameter.
   real(R8P), intent(in)  :: wL(NV_MHD) !< Left frame primitive variables.
   real(R8P), intent(in)  :: wR(NV_MHD) !< Right frame primitive variables.
   real(R8P), intent(out) :: sL, sR     !< Outer wave speeds.
   real(R8P)              :: cf         !< Largest fast speed.
   !$acc routine seq
   !$omp declare target

   cf = max(frame_fast_speed(gamma=gamma, w=wL), frame_fast_speed(gamma=gamma, w=wR))
   sL = min(wL(2), wR(2)) - cf
   sR = max(wL(2), wR(2)) + cf
   endsubroutine frame_speeds

   pure subroutine frame_state(gamma, pv, q, bn, w, u)
   !< Compute the frame primitive and conservative variables of a state whose normal field is set to `bn`: `w(5) = P(rho)`,
   !< the energy slot of `u` zero (`q` may be longer than `NV_MHD`: its leading entries).
   real(R8P),    intent(in)  :: gamma      !< EOS parameter.
   integer(I4P), intent(in)  :: pv(NV_MHD) !< State in the frame order.
   real(R8P),    intent(in)  :: q(NV_MHD)  !< Conservative variables.
   real(R8P),    intent(in)  :: bn         !< Normal field of the frame state.
   real(R8P),    intent(out) :: w(NV_MHD)  !< Frame primitive variables.
   real(R8P),    intent(out) :: u(NV_MHD)  !< Frame conservative variables.
   !$acc routine seq
   !$omp declare target

   w(1) = q(IQ_R)
   w(2) = q(pv(2)) / q(IQ_R)
   w(3) = q(pv(3)) / q(IQ_R)
   w(4) = q(pv(4)) / q(IQ_R)
   w(5) = eos_pressure(gamma, q(IQ_R))
   w(6) = bn
   w(7) = q(pv(7))
   w(8) = q(pv(8))
   u(1) = q(IQ_R)
   u(2) = q(pv(2))
   u(3) = q(pv(3))
   u(4) = q(pv(4))
   u(5) = 0._R8P
   u(6) = bn
   u(7) = w(7)
   u(8) = w(8)
   endsubroutine frame_state

   pure subroutine glm_face_subsystem(ch, bnL, bnR, psiL, psiR, bn, psi)
   !< Solve the linear `(B_n, psi)` subsystem of the mixed GLM exactly at the face (Dedner et al. 2002, eq. 42).
   real(R8P), intent(in)  :: ch         !< GLM cleaning speed.
   real(R8P), intent(in)  :: bnL, bnR   !< Normal fields of the two states.
   real(R8P), intent(in)  :: psiL, psiR !< GLM scalars of the two states.
   real(R8P), intent(out) :: bn, psi    !< Face normal field and GLM scalar.
   !$acc routine seq
   !$omp declare target

   bn  = 0.5_R8P * (bnL + bnR) - 0.5_R8P * (psiR - psiL) / ch
   psi = 0.5_R8P * (psiL + psiR) - 0.5_R8P * ch * (bnR - bnL)
   endsubroutine glm_face_subsystem

   pure subroutine hll_select(sL, sR, uL, uR, fL, fR, f)
   !< Return the HLL flux of two frame states from their outer speeds and physical fluxes.
   real(R8P), intent(in)  :: sL, sR     !< Outer wave speeds.
   real(R8P), intent(in)  :: uL(NV_MHD) !< Left frame conservative variables.
   real(R8P), intent(in)  :: uR(NV_MHD) !< Right frame conservative variables.
   real(R8P), intent(in)  :: fL(NV_MHD) !< Left physical flux.
   real(R8P), intent(in)  :: fR(NV_MHD) !< Right physical flux.
   real(R8P), intent(out) :: f(NV_MHD)  !< Flux.
   integer(I4P)           :: k          !< Counter.
   !$acc routine seq
   !$omp declare target

   if (sL >= 0._R8P) then
      f = fL
   elseif (sR <= 0._R8P) then
      f = fR
   else
      do k=1, NV_MHD
         f(k) = (sR * fL(k) - sL * fR(k) + sL * sR * (uR(k) - uL(k))) / (sR - sL)
      enddo
   endif
   endsubroutine hll_select

   pure subroutine hlld_star(s, ssL, ssR, rs, ms, us, w, ust)
   !< Build the outer star state of a side (Mignone 2007, eqs. 20-21 and 30-33): the HLL density and normal momentum,
   !< the transverse momentum and field from the jump across the outer wave. Near the degeneracy `(s - S*_L)(s - S*_R)
   !< -> 0` (no transverse field, `B_n^2 > a^2 rho`) the outer transverse state is kept.
   real(R8P), intent(in)  :: s           !< Outer wave speed of the side.
   real(R8P), intent(in)  :: ssL, ssR    !< Rotational wave speeds.
   real(R8P), intent(in)  :: rs, ms      !< HLL density and normal momentum.
   real(R8P), intent(in)  :: us          !< Advective velocity.
   real(R8P), intent(in)  :: w(NV_MHD)   !< Frame primitive variables of the side.
   real(R8P), intent(out) :: ust(NV_MHD) !< Star state.
   real(R8P)              :: den         !< (s - S*_L)(s - S*_R).
   real(R8P)              :: rdd         !< rho (s - u_n)^2.
   integer(I4P)           :: k           !< Counter.
   !$acc routine seq
   !$omp declare target

   den = (s - ssL) * (s - ssR)
   rdd = w(1) * (s - w(2))**2
   ust(1) = rs
   ust(2) = ms
   ust(5) = 0._R8P
   ust(6) = w(6)
   if (abs(den) * rs <= EPS_HLLD * (rdd + w(6)**2)) then
      do k=1, 2
         ust(2+k) = rs * w(2+k)
         ust(6+k) = w(6+k)
      enddo
   else
      do k=1, 2
         ust(2+k) = rs * w(2+k) - w(6) * w(6+k) * (us - w(2)) / den
         ust(6+k) = w(6+k) / rs * (rdd - w(6)**2) / den
      enddo
   endif
   endsubroutine hlld_star

   pure function is_admissible(q) result(ok)
   !< Return true when a state has positive density (the one positive quantity without an energy equation).
   real(R8P), intent(in) :: q(NV_MHD) !< Conservative variables.
   logical               :: ok        !< Admissible state.
   !$acc routine seq
   !$omp declare target

   ok = q(IQ_R) > 0._R8P
   endfunction is_admissible

   pure subroutine scatter_glm_flux(ch, pv, fr, bn, psi, f)
   !< Return the GLM face flux in the global order: the frame flux scattered, the `B_n` flux `psi~`, the `psi` flux
   !< `c_h^2 B~_n`.
   real(R8P),    intent(in)  :: ch            !< GLM cleaning speed.
   integer(I4P), intent(in)  :: pv(NV_MHD)    !< State in the frame order.
   real(R8P),    intent(in)  :: fr(NV_MHD)    !< Frame flux.
   real(R8P),    intent(in)  :: bn, psi       !< Face normal field and GLM scalar.
   real(R8P),    intent(out) :: f(NV_MHD_GLM) !< Flux.
   integer(I4P)              :: k             !< Counter.
   !$acc routine seq
   !$omp declare target

   do k=1, NV_MHD
      f(pv(k)) = fr(k)
   enddo
   f(pv(6))  = psi
   f(IQ_PSI) = ch * ch * bn
   endsubroutine scatter_glm_flux

   pure function wu_speed(d, qaL, qaR) result(sigma)
   !< Return the Wu (2018) Lax-Friedrichs speed of two states in direction `d` (as the adiabatic library).
   integer(I4P), intent(in) :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in) :: qaL(NV_AUX_MHD) !< Left auxiliary variables.
   real(R8P),    intent(in) :: qaR(NV_AUX_MHD) !< Right auxiliary variables.
   real(R8P)                :: sigma           !< Speed.
   real(R8P)                :: srL, srR        !< sqrt(rho) of the two states.
   real(R8P)                :: cfL, cfR        !< Fast speeds.
   real(R8P)                :: un              !< sqrt(rho)-weighted normal velocity.
   !$acc routine seq
   !$omp declare target

   srL = sqrt(qaL(IA_R))
   srR = sqrt(qaR(IA_R))
   cfL = mhd_fast_speed(d=d, qa=qaL)
   cfR = mhd_fast_speed(d=d, qa=qaR)
   un  = (srL * qaL(IA_U+d-1) + srR * qaR(IA_U+d-1)) / (srL + srR)
   sigma = max(abs(qaL(IA_U+d-1)) + cfL, abs(qaR(IA_U+d-1)) + cfR, abs(un) + max(cfL, cfR)) + &
           sqrt(mhd_sum3((qaL(IA_BX) - qaR(IA_BX))**2, (qaL(IA_BY) - qaR(IA_BY))**2, (qaL(IA_BZ) - qaR(IA_BZ))**2)) / &
           (srL + srR)
   endfunction wu_speed
