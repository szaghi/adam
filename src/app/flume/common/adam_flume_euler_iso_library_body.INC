!< ADAM, FLUME pointwise isothermal Euler physics, the body of adam_flume_euler_iso_library (issue #100, M13).
!<
!< The including module renames NV_K, NV_AUX_K (its widths) and NV_BASE = NV_EULER, NV_AUX_BASE = NV_AUX (the Euler
!< widths): the components NV_BASE+1 ... NV_K are advected scalars, as in the Euler body; with NV_K = NV_BASE every
!< scalar loop is empty. The public interface is the Euler one, routine by routine, so the model-agnostic kernel bodies
!< call either library unchanged; the EOS parameter keeps the Euler keyword `gamma`, and HERE IT IS THE ISOTHERMAL
!< SOUND SPEED SQUARED c0^2 (the host passes `physics%cs2`).
!<
!< Isothermal gas, P = c0^2 rho (D-M13-1): the energy slot IQ_RE is inert. Its flux is zero, every Riemann flux forces
!< it to zero, its eigenvalue is zero and its eigenvector the unit vector, so it keeps its initial value, 0. The
!< auxiliary pressure is c0^2 rho, the speed of sound c0, the temperature c0^2 / R, the total enthalpy undefined (0).
!<
!< Eigensystem in direction `d`: normal `n = e_d`, tangents taken cyclically `t1 = e_{mod(d,3)+1}`, `t2 =
!< e_{mod(d+1,3)+1}`, waves in the Euler order with the inert slot in the place of the entropy wave:
!<```
!< r1 = (1, u - a n, 0)       l1 = ((u_n + a) / (2 a), -n / (2 a), 0)       lambda1 = u_n - a
!< r2 = (0, 0, 1)             l2 = (0, 0, 1)                               lambda2 = 0      (the inert slot)
!< r3 = (0, t1, 0)            l3 = (-u.t1, t1, 0)                          lambda3 = u_n
!< r4 = (0, t2, 0)            l4 = (-u.t2, t2, 0)                          lambda4 = u_n
!< r5 = (1, u + a n, 0)       l5 = (-(u_n - a) / (2 a), n / (2 a), 0)      lambda5 = u_n + a
!<```
!< checked against `dF/dU` in src/tests/flume/verification/eos/eos_oracle.py (eigen). The Roe average of two isothermal
!< states is the sqrt(rho)-weighted velocity with the same speed of sound.
   ! public procedures
   pure subroutine compute_eigenvalues(d, qa, lambda)
   !< Compute the eigenvalues of the flux Jacobian in direction `d`: `(u_n - a, 0, u_n, u_n, u_n + a)`.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_K) !< Auxiliary variables.
   real(R8P),    intent(out) :: lambda(NV_K) !< Eigenvalues.
   real(R8P)                 :: un           !< Normal velocity.
   integer(I4P)              :: v            !< Counter.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   lambda(1) = un - qa(IA_A)
   lambda(2) = 0._R8P
   lambda(3) = un
   lambda(4) = un
   lambda(5) = un + qa(IA_A)
   do v=NV_BASE+1, NV_K ! advected scalars
      lambda(v) = un
   enddo
   endsubroutine compute_eigenvalues

   pure subroutine compute_eigenvectors(gamma, d, qa, el, er)
   !< Compute the left (rows) and right (columns) eigenvectors of the flux Jacobian in direction `d` at state `qa`.
   real(R8P),    intent(in)  :: gamma         !< Isothermal sound speed squared c0^2 (unused: a is in qa).
   integer(I4P), intent(in)  :: d             !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qa(NV_AUX_K)  !< Auxiliary variables (e.g. a Roe average).
   real(R8P),    intent(out) :: el(NV_K,NV_K) !< Left eigenvectors, el(k,:) = l_k.
   real(R8P),    intent(out) :: er(NV_K,NV_K) !< Right eigenvectors, er(:,k) = r_k.
   real(R8P)                 :: u(3)          !< Velocity.
   real(R8P)                 :: a             !< Speed of sound.
   integer(I4P)              :: d1, d2        !< Tangential directions.
   integer(I4P)              :: c, v          !< Counters.
   real(R8P)                 :: sc            !< Specific advected scalar.
   !$acc routine seq
   !$omp declare target

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   u  = [qa(IA_U), qa(IA_V), qa(IA_W)]
   a  = qa(IA_A)
   er = 0._R8P
   el = 0._R8P
   ! right eigenvectors (columns)
   er(1,1) = 1._R8P
   er(1,5) = 1._R8P
   do c=1, 3
      er(1+c,1) = u(c)
      er(1+c,5) = u(c)
   enddo
   er(1+d,1)  = u(d) - a
   er(1+d,5)  = u(d) + a
   er(IQ_RE,2) = 1._R8P
   er(1+d1,3) = 1._R8P
   er(1+d2,4) = 1._R8P
   ! left eigenvectors (rows)
   el(1,1)    = 0.5_R8P * (u(d) + a) / a
   el(1,1+d)  = -0.5_R8P / a
   el(2,IQ_RE) = 1._R8P
   el(3,1)    = -u(d1)
   el(3,1+d1) = 1._R8P
   el(4,1)    = -u(d2)
   el(4,1+d2) = 1._R8P
   el(5,1)    = -0.5_R8P * (u(d) - a) / a
   el(5,1+d)  = 0.5_R8P / a
   do v=NV_BASE+1, NV_K ! advected scalars: s rides the acoustic waves, its own wave is e_v
      sc = qa(IA_S+v-NV_BASE-1)
      er(v,1) = sc
      er(v,5) = sc
      er(v,v) = 1._R8P
      el(v,1) = -sc
      el(v,v) = 1._R8P
   enddo
   endsubroutine compute_eigenvectors

   pure subroutine compute_face_flux_back_projection(is_characteristic, er, vr, flux)
   !< Return the face flux in conservative variables from the reconstructed split fields, `F = R (f+ + f-)`.
   logical,   intent(in)  :: is_characteristic !< Reconstruction variables: characteristic or conservative.
   real(R8P), intent(in)  :: er(NV_K,NV_K)     !< Right eigenvectors, er(:,k) = r_k (unused if conservative).
   real(R8P), intent(in)  :: vr(2,NV_K)        !< Reconstructed split fields at the face.
   real(R8P), intent(out) :: flux(NV_K)        !< Face flux.
   integer(I4P)           :: k, v              !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_K
         flux(v) = 0._R8P
         do k=1, NV_K
            flux(v) = flux(v) + er(v,k) * (vr(1,k) + vr(2,k))
         enddo
      enddo
   else
      do v=1, NV_K
         flux(v) = vr(1,v) + vr(2,v)
      enddo
   endif
   endsubroutine compute_face_flux_back_projection

   pure subroutine compute_face_interpolation_fields(gamma, d, S, is_characteristic, qs, qas, fint, er, mu)
   !< Return the fields of the stencil of face `i+1/2` to be WENO-interpolated to the face (scheme `weno-riemann`), in
   !< the layout of the Euler body: characteristic `w = L q` on the Roe average of cells 0 and 1, or primitive
   !< `(rho, u, v, w, slot)`, the inert slot in the place of the pressure (P = c0^2 rho needs no separate field).
   real(R8P),    intent(in)  :: gamma                         !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d                             !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                             !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic             !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_K,1-S_MAX:S_MAX)        !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_K,1-S_MAX:S_MAX)   !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_K)  !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_K,NV_K)                 !< Right eigenvectors (unused if primitive).
   real(R8P),    intent(out) :: mu(NV_K)                      !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_K,NV_K)                 !< Left eigenvectors.
   real(R8P)                 :: roe(NV_AUX_K)                 !< Roe average of cells 0 and 1.
   real(R8P)                 :: w(NV_K)                       !< Fields of one cell.
   real(R8P)                 :: qm(NV_K)                      !< Reference magnitudes of one cell (descaler).
   integer(I4P)              :: k, m, v                       !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      call compute_roe_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), roe=roe)
      call compute_eigenvectors(gamma=gamma, d=d, qa=roe, el=el, er=er)
   else
      er = 0._R8P
   endif
   mu = 0._R8P
   do m=1-S, S
      if (is_characteristic) then
         call compute_state_magnitudes(q=qs(:,m), qa=qas(:,m), qm=qm)
         do k=1, NV_K
            w(k) = 0._R8P
            do v=1, NV_K
               w(k) = w(k) + el(k,v) * qs(v,m)
               mu(k) = mu(k) + abs(el(k,v)) * qm(v)
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = qs(IQ_RE,m)
         do v=NV_BASE+1, NV_K ! advected scalars: the specific value
            w(v) = qas(IA_S+v-NV_BASE-1,m)
         enddo
      endif
      do k=1, NV_K
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
         if (.not.is_characteristic) mu(k) = mu(k) + abs(w(k))
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine compute_face_interpolation_fields

   pure subroutine compute_face_split_fluxes(gamma, d, S, is_characteristic, qs, qas, fsplit, er, mu)
   !< Project and Lax-Friedrichs-split the stencil of face `i+1/2`, ready for the WENO upwind reconstruction: the Euler
   !< body's algorithm on the isothermal eigensystem (the inert field has speed 0, state 0 and flux 0: it splits to 0).
   real(R8P),    intent(in)  :: gamma                          !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d                              !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                              !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic              !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)  :: qs(NV_K,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_K,1-S_MAX:S_MAX)    !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fsplit(2,1-S_MAX:S_MAX-1,NV_K) !< Split fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_K,NV_K)                  !< Right eigenvectors (identity if conservative).
   real(R8P),    intent(out) :: mu(NV_K)                       !< Magnitude of each field (WENO descaler).
   real(R8P)                 :: el(NV_K,NV_K)                  !< Left eigenvectors (identity if conservative).
   real(R8P)                 :: roe(NV_AUX_K)                  !< Roe average of cells 0 and 1.
   real(R8P)                 :: lambda(NV_K)                   !< Eigenvalues of one cell.
   real(R8P)                 :: alpha(NV_K)                    !< Lax-Friedrichs speeds.
   real(R8P)                 :: f(NV_K)                        !< Physical flux of one cell.
   real(R8P)                 :: w, g, fp                       !< Projected state, projected flux, split flux.
   real(R8P)                 :: a                              !< Absolute projection of one cell (descaler).
   real(R8P)                 :: qm(NV_K)                       !< Reference magnitudes of one cell (descaler).
   real(R8P)                 :: amax                           !< Largest speed of the stencil (descaler).
   integer(I4P)              :: k, m, v                        !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      call compute_roe_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), roe=roe)
      call compute_eigenvectors(gamma=gamma, d=d, qa=roe, el=el, er=er)
      alpha = 0._R8P
      do m=1-S, S
         call compute_eigenvalues(d=d, qa=qas(:,m), lambda=lambda)
         do k=1, NV_K
            alpha(k) = max(alpha(k), abs(lambda(k)))
         enddo
      enddo
   else
      el = 0._R8P
      er = 0._R8P
      do k=1, NV_K
         el(k,k) = 1._R8P
         er(k,k) = 1._R8P
      enddo
      alpha = 0._R8P
      do m=1-S, S
         alpha(1) = max(alpha(1), abs(qas(IA_U+d-1,m)) + qas(IA_A,m))
      enddo
      alpha = alpha(1)
   endif
   amax = maxval(alpha)
   mu = 0._R8P
   do m=1-S, S
      call compute_flux(d=d, q=qs(:,m), qa=qas(:,m), f=f)
      call compute_state_magnitudes(q=qs(:,m), qa=qas(:,m), qm=qm)
      do k=1, NV_K
         w = 0._R8P
         g = 0._R8P
         a = 0._R8P
         do v=1, NV_K
            w = w + el(k,v) * qs(v,m)
            g = g + el(k,v) * f(v)
            a = a + abs(el(k,v)) * (abs(f(v)) + amax * qm(v))
         enddo
         fp = 0.5_R8P * (g + alpha(k) * w)
         if (m < S)     fsplit(2,m,k)   = fp
         if (m > 1 - S) fsplit(1,m-1,k) = g - fp
         mu(k) = mu(k) + 0.5_R8P * a
      enddo
   enddo
   mu = max(mu / real(2*S, R8P), tiny(1._R8P))
   endsubroutine compute_face_split_fluxes

   pure subroutine compute_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Return the conservative left and right states of face `i+1/2` from the WENO-interpolated fields. A state with a
   !< non-positive density (or an inadmissible scalar) is replaced by the state of its cell: the face falls back to
   !< first order. The pressure follows from the density, so the density is the only thermodynamic condition.
   real(R8P), intent(in)  :: gamma             !< Isothermal sound speed squared c0^2.
   logical,   intent(in)  :: is_characteristic !< Characteristic (or primitive) variables.
   real(R8P), intent(in)  :: er(NV_K,NV_K)     !< Right eigenvectors (unused if primitive).
   real(R8P), intent(in)  :: vr(2,NV_K)        !< Interpolated fields.
   real(R8P), intent(in)  :: q0(NV_K)          !< Conservative variables of cell 0 (left of the face).
   real(R8P), intent(in)  :: q1(NV_K)          !< Conservative variables of cell 1 (right of the face).
   real(R8P), intent(out) :: qL(NV_K)          !< Left state.
   real(R8P), intent(out) :: qR(NV_K)          !< Right state.
   integer(I4P)           :: k, v              !< Counters.
   logical                :: ok                !< Admissibility of a state.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_K
         qL(v) = 0._R8P
         qR(v) = 0._R8P
         do k=1, NV_K
            qL(v) = qL(v) + er(v,k) * vr(2,k)
            qR(v) = qR(v) + er(v,k) * vr(1,k)
         enddo
      enddo
   else
      call primitive_to_conservative(gamma=gamma, r=vr(2,1), u=vr(2,2), v=vr(2,3), w=vr(2,4), p=0._R8P, q=qL)
      call primitive_to_conservative(gamma=gamma, r=vr(1,1), u=vr(1,2), v=vr(1,3), w=vr(1,4), p=0._R8P, q=qR)
      do v=NV_BASE+1, NV_K ! advected scalars
         qL(v) = vr(2,1) * vr(2,v)
         qR(v) = vr(1,1) * vr(1,v)
      enddo
   endif
   ! positive density (false for NaN too); no internal procedure: device routine
   ok = qL(IQ_R) > 0._R8P
   do v=NV_BASE+1, NV_K ! advected scalars: a passive scalar non-negative (issue #92, D-M10-2)
      if (ok) ok = qL(v) > 0._R8P .or. (v > NV_POSITIVE .and. qL(v) >= 0._R8P)
   enddo
   if (.not.ok) qL = q0
   ok = qR(IQ_R) > 0._R8P
   do v=NV_BASE+1, NV_K
      if (ok) ok = qR(v) > 0._R8P .or. (v > NV_POSITIVE .and. qR(v) >= 0._R8P)
   enddo
   if (.not.ok) qR = q1
   endsubroutine compute_face_states

   pure subroutine compute_flux(d, q, qa, f)
   !< Compute the physical flux in direction `d`: `rho u_d (1, u, v, w, 0) + p (0, e_d, 0)`, p = c0^2 rho.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: q(NV_K)      !< Conservative variables.
   real(R8P),    intent(in)  :: qa(NV_AUX_K) !< Auxiliary variables.
   real(R8P),    intent(out) :: f(NV_K)      !< Physical flux.
   real(R8P)                 :: un           !< Normal velocity.
   integer(I4P)              :: v            !< Counter.
   !$acc routine seq
   !$omp declare target

   un = qa(IA_U+d-1)
   f(IQ_R)  = q(IQ_R)  * un
   f(IQ_RU) = q(IQ_RU) * un
   f(IQ_RV) = q(IQ_RV) * un
   f(IQ_RW) = q(IQ_RW) * un
   f(IQ_RU+d-1) = f(IQ_RU+d-1) + qa(IA_P)
   f(IQ_RE) = 0._R8P
   do v=NV_BASE+1, NV_K ! advected scalars
      f(v) = q(v) * un
   enddo
   endsubroutine compute_flux

   pure subroutine compute_roe_average(gamma, qaL, qaR, roe)
   !< Compute the Roe average of two isothermal states: the `sqrt(rho)`-weighted velocity, the density `sqrt(rho_L
   !< rho_R)`, the speed of sound c0 and the matching pressure `c0^2 rho`; enthalpy and temperature are set to zero.
   real(R8P), intent(in)  :: gamma         !< Isothermal sound speed squared c0^2.
   real(R8P), intent(in)  :: qaL(NV_AUX_K) !< Left state auxiliary variables.
   real(R8P), intent(in)  :: qaR(NV_AUX_K) !< Right state auxiliary variables.
   real(R8P), intent(out) :: roe(NV_AUX_K) !< Roe average.
   real(R8P)              :: wL, wR        !< Square roots of the densities.
   real(R8P)              :: cL, cR        !< Normalized weights.
   integer(I4P)           :: v             !< Counter.
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
   roe(IA_H) = 0._R8P
   roe(IA_A) = sqrt(gamma)
   roe(IA_P) = gamma * roe(IA_R)
   roe(IA_T) = 0._R8P
   do v=NV_AUX_BASE+1, NV_AUX_K ! advected scalars, sqrt(rho)-weighted like the velocity
      roe(v) = cL * qaL(v) + cR * qaR(v)
   enddo
   endsubroutine compute_roe_average

   pure subroutine compute_riemann_hll(gamma, d, qL, qR, f)
   !< Compute the HLL flux of two isothermal states in direction `d` with the Einfeldt speeds `S_L = min(u_nL - c0,
   !< u~_n - c0)`, `S_R = max(u_nR + c0, u~_n + c0)` (`~` the Roe average); the inert slot's flux is zero.
   real(R8P),    intent(in)  :: gamma         !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d             !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K)      !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K)      !< Right state.
   real(R8P),    intent(out) :: f(NV_K)       !< Flux.
   real(R8P)                 :: qaL(NV_AUX_K) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_K) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_K)      !< Left physical flux.
   real(R8P)                 :: fR(NV_K)      !< Right physical flux.
   real(R8P)                 :: sL, sR        !< Left and right wave speeds.
   !$acc routine seq
   !$omp declare target

   call hll_flux(gamma=gamma, d=d, qL=qL, qR=qR, qaL=qaL, qaR=qaR, fL=fL, fR=fR, sL=sL, sR=sR, f=f)
   endsubroutine compute_riemann_hll

   pure subroutine compute_riemann_hllc(gamma, d, qL, qR, f)
   !< Compute the isothermal HLLC flux of two states in direction `d`: HLL for the mass and the normal momentum, the
   !< transverse momenta and the advected scalars the HLL mass flux times the upwind specific value (the contact, the
   !< only wave the isothermal system has between the acoustic ones, carries the tangential velocity and the scalars
   !< and nothing else: it needs no star pressure). Checked against eos_oracle.py (flux: consistency, isolated contact).
   real(R8P),    intent(in)  :: gamma         !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d             !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K)      !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K)      !< Right state.
   real(R8P),    intent(out) :: f(NV_K)       !< Flux.
   real(R8P)                 :: qaL(NV_AUX_K) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_K) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_K)      !< Left physical flux.
   real(R8P)                 :: fR(NV_K)      !< Right physical flux.
   real(R8P)                 :: sL, sR        !< Left and right wave speeds.
   integer(I4P)              :: d1, d2        !< Tangential directions.
   integer(I4P)              :: v             !< Counter.
   !$acc routine seq
   !$omp declare target

   call hll_flux(gamma=gamma, d=d, qL=qL, qR=qR, qaL=qaL, qaR=qaR, fL=fL, fR=fR, sL=sL, sR=sR, f=f)
   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   if (f(IQ_R) >= 0._R8P) then
      f(IQ_RU+d1-1) = f(IQ_R) * qaL(IA_U+d1-1)
      f(IQ_RU+d2-1) = f(IQ_R) * qaL(IA_U+d2-1)
      do v=NV_BASE+1, NV_K
         f(v) = f(IQ_R) * qaL(IA_S+v-NV_BASE-1)
      enddo
   else
      f(IQ_RU+d1-1) = f(IQ_R) * qaR(IA_U+d1-1)
      f(IQ_RU+d2-1) = f(IQ_R) * qaR(IA_U+d2-1)
      do v=NV_BASE+1, NV_K
         f(v) = f(IQ_R) * qaR(IA_S+v-NV_BASE-1)
      enddo
   endif
   endsubroutine compute_riemann_hllc

   pure subroutine compute_riemann_llf(gamma, d, qL, qR, f)
   !< Compute the local Lax-Friedrichs (Rusanov) flux of two states in direction `d`,
   !< `F = (f(qL) + f(qR)) / 2 - alpha (qR - qL) / 2`, `alpha = max(|u_n|) + c0`; the inert slot's flux is zero.
   real(R8P),    intent(in)  :: gamma         !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d             !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K)      !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K)      !< Right state.
   real(R8P),    intent(out) :: f(NV_K)       !< Flux.
   real(R8P)                 :: qaL(NV_AUX_K) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_K) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_K)      !< Left physical flux.
   real(R8P)                 :: fR(NV_K)      !< Right physical flux.
   real(R8P)                 :: alpha         !< Lax-Friedrichs speed.
   integer(I4P)              :: v             !< Counter.
   !$acc routine seq
   !$omp declare target

   ! the gas constant only sets the temperature, which the flux does not use
   call conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qL, qa=qaL)
   call conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qR, qa=qaR)
   call compute_flux(d=d, q=qL, qa=qaL, f=fL)
   call compute_flux(d=d, q=qR, qa=qaR, f=fR)
   alpha = max(abs(qaL(IA_U+d-1)) + qaL(IA_A), abs(qaR(IA_U+d-1)) + qaR(IA_A))
   do v=1, NV_K
      f(v) = 0.5_R8P * (fL(v) + fR(v)) - 0.5_R8P * alpha * (qR(v) - qL(v))
   enddo
   f(IQ_RE) = 0._R8P
   endsubroutine compute_riemann_llf

   pure subroutine conservative_to_auxiliary(gamma, R, q, qa)
   !< Compute the auxiliary variables of a cell from its conservative variables: `p = c0^2 rho`, `a = c0`,
   !< `T = c0^2 / R`, `H = 0` (no energy equation).
   real(R8P), intent(in)  :: gamma        !< Isothermal sound speed squared c0^2.
   real(R8P), intent(in)  :: R            !< Gas constant.
   real(R8P), intent(in)  :: q(NV_K)      !< Conservative variables.
   real(R8P), intent(out) :: qa(NV_AUX_K) !< Auxiliary variables.
   integer(I4P)           :: v            !< Counter.
   !$acc routine seq
   !$omp declare target

   qa(IA_R) = q(IQ_R)
   qa(IA_U) = q(IQ_RU) / q(IQ_R)
   qa(IA_V) = q(IQ_RV) / q(IQ_R)
   qa(IA_W) = q(IQ_RW) / q(IQ_R)
   qa(IA_P) = gamma * q(IQ_R)
   qa(IA_T) = gamma / R
   qa(IA_H) = 0._R8P
   qa(IA_A) = sqrt(gamma)
   do v=NV_BASE+1, NV_K ! advected scalars, specific
      qa(IA_S+v-NV_BASE-1) = q(v) / q(IQ_R)
   enddo
   endsubroutine conservative_to_auxiliary

   pure subroutine primitive_to_conservative(gamma, r, u, v, w, p, q)
   !< Compute the conservative variables of a cell from its primitive variables: the pressure is not a state variable
   !< of the isothermal gas (P = c0^2 rho), so `p` is ignored and the inert slot set to 0 (the advected scalars, if any,
   !< are the caller's).
   real(R8P), intent(in)  :: gamma   !< Isothermal sound speed squared c0^2 (unused).
   real(R8P), intent(in)  :: r       !< Density.
   real(R8P), intent(in)  :: u, v, w !< Velocity components.
   real(R8P), intent(in)  :: p       !< Pressure (ignored).
   real(R8P), intent(out) :: q(NV_K) !< Conservative variables.
   !$acc routine seq
   !$omp declare target

   q(IQ_R)  = r
   q(IQ_RU) = r * u
   q(IQ_RV) = r * v
   q(IQ_RW) = r * w
   q(IQ_RE) = 0._R8P
   endsubroutine primitive_to_conservative

   ! private procedures
   pure subroutine hll_flux(gamma, d, qL, qR, qaL, qaR, fL, fR, sL, sR, f)
   !< The HLL flux with the Einfeldt speeds, its auxiliary states, physical fluxes and speeds (shared by HLL and HLLC).
   real(R8P),    intent(in)  :: gamma         !< Isothermal sound speed squared c0^2.
   integer(I4P), intent(in)  :: d             !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_K)      !< Left state.
   real(R8P),    intent(in)  :: qR(NV_K)      !< Right state.
   real(R8P),    intent(out) :: qaL(NV_AUX_K) !< Left auxiliary variables.
   real(R8P),    intent(out) :: qaR(NV_AUX_K) !< Right auxiliary variables.
   real(R8P),    intent(out) :: fL(NV_K)      !< Left physical flux.
   real(R8P),    intent(out) :: fR(NV_K)      !< Right physical flux.
   real(R8P),    intent(out) :: sL, sR        !< Left and right wave speeds.
   real(R8P),    intent(out) :: f(NV_K)       !< Flux.
   real(R8P)                 :: roe(NV_AUX_K) !< Roe average.
   integer(I4P)              :: v             !< Counter.
   !$acc routine seq
   !$omp declare target

   ! the gas constant only sets the temperature, which the flux does not use
   call conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qL, qa=qaL)
   call conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qR, qa=qaR)
   call compute_flux(d=d, q=qL, qa=qaL, f=fL)
   call compute_flux(d=d, q=qR, qa=qaR, f=fR)
   call compute_roe_average(gamma=gamma, qaL=qaL, qaR=qaR, roe=roe)
   sL = min(qaL(IA_U+d-1) - qaL(IA_A), roe(IA_U+d-1) - roe(IA_A))
   sR = max(qaR(IA_U+d-1) + qaR(IA_A), roe(IA_U+d-1) + roe(IA_A))
   if (sL >= 0._R8P) then
      f = fL
   elseif (sR <= 0._R8P) then
      f = fR
   else
      do v=1, NV_K
         f(v) = (sR * fL(v) - sL * fR(v) + sL * sR * (qR(v) - qL(v))) / (sR - sL)
      enddo
   endif
   f(IQ_RE) = 0._R8P
   endsubroutine hll_flux

   pure subroutine compute_state_magnitudes(q, qa, qm)
   !< Return the reference magnitude of each conservative variable of a cell, the WENO descaler scale (the Euler body's
   !< definition: `|rho|`, the inert slot `|q_5|`, every momentum component `|rho u|_1 + sqrt(rho p)`).
   real(R8P), intent(in)  :: q(NV_K)      !< Conservative variables.
   real(R8P), intent(in)  :: qa(NV_AUX_K) !< Auxiliary variables.
   real(R8P), intent(out) :: qm(NV_K)     !< Reference magnitudes.
   real(R8P)              :: a, b, c      !< Absolute momentum components.
   real(R8P)              :: pm           !< Momentum magnitude.
   integer(I4P)           :: v            !< Counter.
   !$acc routine seq
   !$omp declare target

   a = abs(q(IQ_RU))
   b = abs(q(IQ_RV))
   c = abs(q(IQ_RW))
   pm = ((min(a, b, c) + max(min(a, b), min(max(a, b), c))) + max(a, b, c)) + sqrt(abs(qa(IA_R) * qa(IA_P)))
   qm(IQ_R)  = abs(q(IQ_R))
   qm(IQ_RU) = pm
   qm(IQ_RV) = pm
   qm(IQ_RW) = pm
   qm(IQ_RE) = abs(q(IQ_RE))
   do v=NV_BASE+1, NV_K ! advected scalars
      qm(v) = abs(q(v))
   enddo
   endsubroutine compute_state_magnitudes
