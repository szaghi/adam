!< ADAM, FLUME pointwise ideal MHD Riemann solvers and face states of the scheme `weno-riemann` (issue #47).
module adam_flume_mhd_riemann_library
!< ADAM, FLUME pointwise ideal MHD Riemann solvers and face states of the scheme `weno-riemann` (issue #47).
!<
!< Same contract as `adam_flume_mhd_library`: `pure`, explicit-size dummies, `!$acc routine seq` + `!$omp declare target`,
!< no branching on the model; one set of public routines per model variant (`mhd_*`: `NV_MHD`, no divergence control;
!< `mhd_glm_*`: `NV_MHD_GLM`, mixed GLM cleaning with speed `ch`; `mhd_eglm_*`: `NV_MHD_EGLM`, EGLM cleaning, `psi` in B
!< units and `psi^2 / 2` in the energy).
!<
!< **Frame.** Every solver works in the frame of direction `d`, primitive `w = (rho, u_n, u_t1, u_t2, p, B_n, B_t1, B_t2)`
!< and conservative `u = (rho, rho u_n, rho u_t1, rho u_t2, E, B_n, B_t1, B_t2)`, tangents cyclic
!< (`mhd_frame_indexes`), every three-term sum through `mhd_sum3`: a problem rotated from x to y or z runs the same
!< arithmetic, so the rotated runs are bitwise equal (issue #41, MV-4).
!<
!< **Solvers** (issue #47, section 3.2), wave speeds `S_L = min(u_nL, u_nR) - max(c_fL, c_fR)`,
!< `S_R = max(u_nL, u_nR) + max(c_fL, c_fR)` (Miyoshi & Kusano 2005, eq. 67):
!< * LLF (Rusanov), speed `max(|u_n| + c_f)`;
!< * HLL (Harten, Lax & van Leer 1983);
!< * HLLD (Miyoshi & Kusano 2005, J. Comput. Phys. 208), one normal field `B_n` for both states. Degenerate cases: when
!<   `rho_a (S_a - u_na)(S_a - S_M) - B_n^2` is negligible the star state takes the outer transverse velocity and field;
!<   `B_n = 0` gives `S*_L = S*_R = S_M` and no double-star region. HLL is used instead, and `fallback` is set, when a
!<   star state has non-positive density or pressure (the pressure of the conservative star state, from its energy)
!<   or the speeds are out of order (`S_L < S*_L <= S_M <= S*_R < S_R`); measured on random pairs whose density,
!<   pressure and field ratios are within 10, 10^2: 0.4%, 1.8% of the faces, mostly `S*_L <= S_L` (the speed
!<   estimate no longer bounds the Alfven wave).
!<
!< **Normal field.** Without divergence control LLF and HLL take each state's own `B_n` (the `B_n` flux is zero, the
!< jump is dissipated), HLLD the average. With GLM the linear `(B_n, psi)` subsystem is solved exactly at the face
!< (Dedner et al. 2002, eq. 42): `B~_n = (B_nL + B_nR)/2 - (psi_R - psi_L)/(2 c_h)`,
!< `psi~ = (psi_L + psi_R)/2 - c_h (B_nR - B_nL)/2`; both states take `B~_n` (pressure kept, energy recomputed: `psi`
!< is not in the energy of the mixed GLM), the `B_n` flux is `psi~` and the `psi` flux `c_h^2 B~_n`. With EGLM (issue
!< #47, section 3.2) the subsystem in B units, `B~_n = (B_nL + B_nR)/2 - (psi_R - psi_L)/2`,
!< `psi~ = (psi_L + psi_R)/2 - (B_nR - B_nL)/2`; the solvers see the MHD states (energy without `psi^2 / 2`); the `B_n`
!< flux is `c_h psi~`, the `psi` flux `c_h B~_n`, and the energy flux gains `c_h psi~ B~_n` plus the advected cleaning
!< energy, a passive scalar carried by the mass flux, `F_rho psi~^2 / (2 rho_up)` with the upwind density of the sign of
!< `F_rho` (Larrouturou 1991): the flux of two equal states is the physical flux.
!<
!< **Positivity backbone** (issue #47, D-9): the first-order Lax-Friedrichs flux of the cell states, `(f(qL) + f(qR)) /
!< 2 - sigma (qR - qL) / 2`, each state with its own `B_n`, `sigma` the Wu (2018, doi:10.1137/18M1168017) speed
!< `max(s_L, s_R, |u~_n| + max(c_fL, c_fR)) + |B_L - B_R| / (sqrt(rho_L) + sqrt(rho_R))`, `s = |u_n| + c_f`, `u~_n` the
!< sqrt(rho)-weighted normal velocity (with EGLM at least `c_h`): with the Godunov-Powell (EGLM) sources it keeps the
!< first-order update admissible (Wu & Shu 2018, doi:10.1137/18M1168042), the backbone of the M3 limiter (prototype
!< M3-P0).
!<
!< **Face states.** Interpolated fields: characteristic (the eigenvectors of the arithmetic face average, projected in the
!< frame order; the MHD default, #47 D-5) or primitive `(rho, u, v, w, p, B_x, B_y, B_z[, psi])`, `p` the thermal
!< pressure; a face state with non-positive density or pressure is replaced by the adjacent cell's state.

! FLUME modules
use :: adam_flume_mhd_library, only : mhd_conservative_to_auxiliary, mhd_eglm_conservative_to_auxiliary,         &
                                      mhd_eglm_eigenvectors, mhd_eglm_flux, mhd_eigenvectors, mhd_face_average,     &
                                      mhd_fast_speed, mhd_flux, mhd_frame_indexes, mhd_glm_eigenvectors,            &
                                      mhd_primitive_to_conservative, mhd_sum3
use :: adam_flume_parameters,  only : IA_BX, IA_BY, IA_BZ, IA_P, IA_R, IA_U, IA_V, IA_W, IQ_BX, IQ_BY, IQ_BZ, IQ_PSI, &
                                      IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW, NV_AUX_MHD, NV_MHD, NV_MHD_EGLM, NV_MHD_GLM,  &
                                      S_MAX
! third party modules
use :: penf,                   only : I4P, R8P

implicit none
private
public :: EPS_HLLD
public :: mhd_backbone_flux
public :: mhd_eglm_backbone_flux
public :: mhd_eglm_face_interpolation_fields
public :: mhd_eglm_face_states
public :: mhd_eglm_riemann_hll
public :: mhd_eglm_riemann_hlld
public :: mhd_eglm_riemann_llf
public :: mhd_face_interpolation_fields
public :: mhd_face_states
public :: mhd_glm_face_interpolation_fields
public :: mhd_glm_face_states
public :: mhd_glm_riemann_hll
public :: mhd_glm_riemann_hlld
public :: mhd_glm_riemann_llf
public :: mhd_riemann_hll
public :: mhd_riemann_hlld
public :: mhd_riemann_llf

real(R8P), parameter :: EPS_HLLD=1.e-12_R8P !< Degenerate HLLD star: |rho d (S-S_M) - B_n^2| <= EPS_HLLD (|rho d (S-S_M)| + B_n^2).

contains
   ! public procedures
   pure subroutine mhd_backbone_flux(gamma, d, qL, qR, f)
   !< Compute the positivity backbone flux of two states in direction `d` (MHD without divergence control): the
   !< Lax-Friedrichs flux with the Wu speed.
   real(R8P),    intent(in)  :: gamma          !< Specific heats ratio.
   integer(I4P), intent(in)  :: d              !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD)     !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD)     !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD)      !< Flux.
   real(R8P)                 :: qaL(NV_AUX_MHD) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_MHD) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_MHD)     !< Left physical flux.
   real(R8P)                 :: fR(NV_MHD)     !< Right physical flux.
   real(R8P)                 :: sigma          !< Lax-Friedrichs speed.
   integer(I4P)              :: v              !< Counter.
   !$acc routine seq
   !$omp declare target

   ! the gas constant only sets the temperature, which the flux does not use
   call mhd_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qL, qa=qaL)
   call mhd_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qR, qa=qaR)
   call mhd_flux(d=d, q=qL, qa=qaL, f=fL)
   call mhd_flux(d=d, q=qR, qa=qaR, f=fR)
   sigma = wu_speed(d=d, qaL=qaL, qaR=qaR)
   do v=1, NV_MHD
      f(v) = 0.5_R8P * (fL(v) + fR(v)) - 0.5_R8P * sigma * (qR(v) - qL(v))
   enddo
   endsubroutine mhd_backbone_flux

   pure subroutine mhd_eglm_backbone_flux(ch, gamma, d, qL, qR, f)
   !< Compute the positivity backbone flux of two states in direction `d` (MHD with EGLM): the Lax-Friedrichs flux of the
   !< EGLM physical fluxes with the Wu speed, at least `c_h`.
   real(R8P),    intent(in)  :: ch              !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma           !< Specific heats ratio.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_EGLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_EGLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_EGLM)  !< Flux.
   real(R8P)                 :: qaL(NV_AUX_MHD) !< Left auxiliary variables.
   real(R8P)                 :: qaR(NV_AUX_MHD) !< Right auxiliary variables.
   real(R8P)                 :: fL(NV_MHD_EGLM) !< Left physical flux.
   real(R8P)                 :: fR(NV_MHD_EGLM) !< Right physical flux.
   real(R8P)                 :: sigma           !< Lax-Friedrichs speed.
   integer(I4P)              :: v               !< Counter.
   !$acc routine seq
   !$omp declare target

   call mhd_eglm_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qL, qa=qaL)
   call mhd_eglm_conservative_to_auxiliary(gamma=gamma, R=1._R8P, q=qR, qa=qaR)
   call mhd_eglm_flux(ch=ch, d=d, q=qL, qa=qaL, f=fL)
   call mhd_eglm_flux(ch=ch, d=d, q=qR, qa=qaR, f=fR)
   sigma = max(wu_speed(d=d, qaL=qaL, qaR=qaR), ch)
   do v=1, NV_MHD_EGLM
      f(v) = 0.5_R8P * (fL(v) + fR(v)) - 0.5_R8P * sigma * (qR(v) - qL(v))
   enddo
   endsubroutine mhd_eglm_backbone_flux

   pure subroutine mhd_eglm_face_interpolation_fields(gamma, d, S, is_characteristic, qs, qas, fint, er)
   !< Compute the fields of the stencil of face `i+1/2` to interpolate (MHD with EGLM): as
   !< `mhd_glm_face_interpolation_fields`, with the EGLM eigenvectors (`psi` the average of cells 0 and 1) or the
   !< primitive fields with the thermal pressure.
   real(R8P),    intent(in)  :: gamma                               !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                                   !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                   !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                   !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD_EGLM,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)       !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_MHD_EGLM) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD_EGLM,NV_MHD_EGLM)          !< Right eigenvectors (unused if primitive).
   real(R8P)                 :: el(NV_MHD_EGLM,NV_MHD_EGLM)          !< Left eigenvectors.
   real(R8P)                 :: avg(NV_AUX_MHD)                     !< Face average of cells 0 and 1.
   real(R8P)                 :: w(NV_MHD_EGLM)                      !< Fields of one cell.
   integer(I4P)              :: pv(NV_MHD_EGLM)                     !< State in the frame order of direction `d`.
   integer(I4P)              :: k, m, v                             !< Counters.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   pv(NV_MHD_EGLM) = IQ_PSI
   if (is_characteristic) then
      call mhd_face_average(gamma=gamma, qaL=qas(:,0), qaR=qas(:,1), avg=avg)
      call mhd_eglm_eigenvectors(gamma=gamma, d=d, qa=avg, psi=0.5_R8P * (qs(IQ_PSI,0) + qs(IQ_PSI,1)), el=el, er=er)
   else
      er = 0._R8P
   endif
   do m=1-S, S
      if (is_characteristic) then
         do k=1, NV_MHD_EGLM
            w(k) = 0._R8P
            do v=1, NV_MHD_EGLM
               w(k) = w(k) + el(k,pv(v)) * qs(pv(v),m)
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = qas(IA_P,m)
         w(6) = qas(IA_BX,m)
         w(7) = qas(IA_BY,m)
         w(8) = qas(IA_BZ,m)
         w(9) = qs(IQ_PSI,m)
      endif
      do k=1, NV_MHD_EGLM
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
      enddo
   enddo
   endsubroutine mhd_eglm_face_interpolation_fields

   pure subroutine mhd_eglm_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Compute the two face states from the interpolated fields (MHD with EGLM): as `mhd_glm_face_states`, the primitive
   !< states with `psi^2 / 2` in the energy; the admissibility is that of the thermal pressure.
   real(R8P), intent(in)  :: gamma                       !< Specific heats ratio.
   logical,   intent(in)  :: is_characteristic           !< Characteristic (or primitive) variables.
   real(R8P), intent(in)  :: er(NV_MHD_EGLM,NV_MHD_EGLM) !< Right eigenvectors (unused if primitive).
   real(R8P), intent(in)  :: vr(2,NV_MHD_EGLM)          !< Interpolated fields.
   real(R8P), intent(in)  :: q0(NV_MHD_EGLM)            !< Conservative variables of cell 0.
   real(R8P), intent(in)  :: q1(NV_MHD_EGLM)            !< Conservative variables of cell 1.
   real(R8P), intent(out) :: qL(NV_MHD_EGLM)            !< Left state.
   real(R8P), intent(out) :: qR(NV_MHD_EGLM)            !< Right state.
   integer(I4P)           :: k, v                        !< Counters.
   !$acc routine seq
   !$omp declare target

   if (is_characteristic) then
      do v=1, NV_MHD_EGLM
         qL(v) = 0._R8P
         qR(v) = 0._R8P
         do k=1, NV_MHD_EGLM
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
      qL(IQ_RE)  = qL(IQ_RE) + 0.5_R8P * qL(IQ_PSI)**2
      qR(IQ_RE)  = qR(IQ_RE) + 0.5_R8P * qR(IQ_PSI)**2
   endif
   if (.not.eglm_is_admissible(q=qL)) qL = q0
   if (.not.eglm_is_admissible(q=qR)) qR = q1
   endsubroutine mhd_eglm_face_states

   pure subroutine mhd_eglm_riemann_hll(ch, gamma, d, qL, qR, f)
   !< Compute the HLL flux of two states in direction `d` (MHD with EGLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch              !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma           !< Specific heats ratio.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_EGLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_EGLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_EGLM)  !< Flux.
   real(R8P)                 :: wL(NV_MHD)      !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)      !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)      !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)      !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)      !< Frame flux.
   real(R8P)                 :: bn, psi         !< Face normal field and EGLM scalar.
   integer(I4P)              :: pv(NV_MHD)      !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call eglm_face_subsystem(bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_hll(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   call scatter_eglm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, rhoL=qL(IQ_R), rhoR=qR(IQ_R), f=f)
   endsubroutine mhd_eglm_riemann_hll

   pure subroutine mhd_eglm_riemann_hlld(ch, gamma, d, qL, qR, f, fallback)
   !< Compute the HLLD flux of two states in direction `d` (MHD with EGLM, the `(B_n, psi)` subsystem solved exactly);
   !< `fallback` is set when HLL replaced HLLD.
   real(R8P),    intent(in)  :: ch              !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma           !< Specific heats ratio.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_EGLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_EGLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_EGLM)  !< Flux.
   logical,      intent(out) :: fallback        !< HLL used instead of HLLD.
   real(R8P)                 :: wL(NV_MHD)      !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)      !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)      !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)      !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)      !< Frame flux.
   real(R8P)                 :: bn, psi         !< Face normal field and EGLM scalar.
   integer(I4P)              :: pv(NV_MHD)      !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call eglm_face_subsystem(bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_hlld(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr, fallback=fallback)
   call scatter_eglm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, rhoL=qL(IQ_R), rhoR=qR(IQ_R), f=f)
   endsubroutine mhd_eglm_riemann_hlld

   pure subroutine mhd_eglm_riemann_llf(ch, gamma, d, qL, qR, f)
   !< Compute the LLF flux of two states in direction `d` (MHD with EGLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch              !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma           !< Specific heats ratio.
   integer(I4P), intent(in)  :: d               !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: qL(NV_MHD_EGLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_EGLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_EGLM)  !< Flux.
   real(R8P)                 :: wL(NV_MHD)      !< Left frame primitive variables.
   real(R8P)                 :: wR(NV_MHD)      !< Right frame primitive variables.
   real(R8P)                 :: uL(NV_MHD)      !< Left frame conservative variables.
   real(R8P)                 :: uR(NV_MHD)      !< Right frame conservative variables.
   real(R8P)                 :: fr(NV_MHD)      !< Frame flux.
   real(R8P)                 :: bn, psi         !< Face normal field and EGLM scalar.
   integer(I4P)              :: pv(NV_MHD)      !< State in the frame order of direction `d`.
   !$acc routine seq
   !$omp declare target

   call mhd_frame_indexes(d=d, pv=pv)
   call eglm_face_subsystem(bnL=qL(pv(6)), bnR=qR(pv(6)), psiL=qL(IQ_PSI), psiR=qR(IQ_PSI), bn=bn, psi=psi)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qL, bn=bn, w=wL, u=uL)
   call eglm_frame_state(gamma=gamma, pv=pv, q=qR, bn=bn, w=wR, u=uR)
   call frame_llf(gamma=gamma, wL=wL, uL=uL, wR=wR, uR=uR, f=fr)
   call scatter_eglm_flux(ch=ch, pv=pv, fr=fr, bn=bn, psi=psi, rhoL=qL(IQ_R), rhoR=qR(IQ_R), f=f)
   endsubroutine mhd_eglm_riemann_llf

   pure subroutine mhd_face_interpolation_fields(gamma, d, S, is_characteristic, qs, qas, fint, er)
   !< Compute the fields of the stencil of face `i+1/2` to interpolate (MHD without divergence control), in the WENO
   !< upwind layout of the Euler `compute_face_interpolation_fields`: `fint(2,m,k)` and `fint(1,m-1,k)` hold field `k`
   !< of cell `m`.
   real(R8P),    intent(in)  :: gamma                          !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                              !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                              !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic              !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)  !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_MHD) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD,NV_MHD)              !< Right eigenvectors (unused if primitive).
   real(R8P)                 :: el(NV_MHD,NV_MHD)              !< Left eigenvectors.
   real(R8P)                 :: avg(NV_AUX_MHD)                !< Face average of cells 0 and 1.
   real(R8P)                 :: w(NV_MHD)                      !< Fields of one cell.
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
   do m=1-S, S
      if (is_characteristic) then
         do k=1, NV_MHD
            w(k) = 0._R8P
            do v=1, NV_MHD
               w(k) = w(k) + el(k,pv(v)) * qs(pv(v),m)
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = qas(IA_P,m)
         w(6) = qas(IA_BX,m)
         w(7) = qas(IA_BY,m)
         w(8) = qas(IA_BZ,m)
      endif
      do k=1, NV_MHD
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
      enddo
   enddo
   endsubroutine mhd_face_interpolation_fields

   pure subroutine mhd_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Compute the two face states from the interpolated fields (MHD without divergence control): the left state from
   !< `vr(2,:)`, the right one from `vr(1,:)`; a state with non-positive density or pressure falls back to its cell.
   real(R8P), intent(in)  :: gamma             !< Specific heats ratio.
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

   pure subroutine mhd_glm_face_interpolation_fields(ch, gamma, d, S, is_characteristic, qs, qas, fint, er)
   !< Compute the fields of the stencil of face `i+1/2` to interpolate (MHD with GLM): as
   !< `mhd_face_interpolation_fields`, `psi` the ninth field (primitive) or the GLM eigenvectors (characteristic).
   real(R8P),    intent(in)  :: ch                                 !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma                              !< Specific heats ratio.
   integer(I4P), intent(in)  :: d                                  !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)  :: S                                  !< WENO stencil half-width, S <= S_MAX.
   logical,      intent(in)  :: is_characteristic                  !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)  :: qs(NV_MHD_GLM,1-S_MAX:S_MAX)       !< Stencil conservative variables.
   real(R8P),    intent(in)  :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)      !< Stencil auxiliary variables.
   real(R8P),    intent(out) :: fint(2,1-S_MAX:S_MAX-1,NV_MHD_GLM) !< Fields in the WENO upwind layout.
   real(R8P),    intent(out) :: er(NV_MHD_GLM,NV_MHD_GLM)          !< Right eigenvectors (unused if primitive).
   real(R8P)                 :: el(NV_MHD_GLM,NV_MHD_GLM)          !< Left eigenvectors.
   real(R8P)                 :: avg(NV_AUX_MHD)                    !< Face average of cells 0 and 1.
   real(R8P)                 :: w(NV_MHD_GLM)                      !< Fields of one cell.
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
   do m=1-S, S
      if (is_characteristic) then
         do k=1, NV_MHD_GLM
            w(k) = 0._R8P
            do v=1, NV_MHD_GLM
               w(k) = w(k) + el(k,pv(v)) * qs(pv(v),m)
            enddo
         enddo
      else
         w(1) = qas(IA_R,m)
         w(2) = qas(IA_U,m)
         w(3) = qas(IA_V,m)
         w(4) = qas(IA_W,m)
         w(5) = qas(IA_P,m)
         w(6) = qas(IA_BX,m)
         w(7) = qas(IA_BY,m)
         w(8) = qas(IA_BZ,m)
         w(9) = qs(IQ_PSI,m)
      endif
      do k=1, NV_MHD_GLM
         if (m < S)     fint(2,m,k)   = w(k)
         if (m > 1 - S) fint(1,m-1,k) = w(k)
      enddo
   enddo
   endsubroutine mhd_glm_face_interpolation_fields

   pure subroutine mhd_glm_face_states(gamma, is_characteristic, er, vr, q0, q1, qL, qR)
   !< Compute the two face states from the interpolated fields (MHD with GLM): as `mhd_face_states`, `psi` included.
   real(R8P), intent(in)  :: gamma                     !< Specific heats ratio.
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
   !< Compute the HLL flux of two states in direction `d` (MHD with GLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< Specific heats ratio.
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
   !< Compute the HLLD flux of two states in direction `d` (MHD with GLM, the `(B_n, psi)` subsystem solved exactly);
   !< `fallback` is set when HLL replaced HLLD.
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< Specific heats ratio.
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
   !< Compute the LLF flux of two states in direction `d` (MHD with GLM, the `(B_n, psi)` subsystem solved exactly).
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   real(R8P),    intent(in)  :: gamma          !< Specific heats ratio.
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
   !< Compute the HLL flux of two states in direction `d` (MHD without divergence control, each state's own `B_n`).
   real(R8P),    intent(in)  :: gamma      !< Specific heats ratio.
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
   !< Compute the HLLD flux of two states in direction `d` (MHD without divergence control, `B_n` the average of the two
   !< states, its flux zero); `fallback` is set when HLL replaced HLLD.
   real(R8P),    intent(in)  :: gamma      !< Specific heats ratio.
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
   !< Compute the LLF flux of two states in direction `d` (MHD without divergence control, each state's own `B_n`).
   real(R8P),    intent(in)  :: gamma      !< Specific heats ratio.
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
   pure subroutine double_star(us, vss, bss, e, uss)
   !< Build the HLLD double-star state of a side from its star state (density, normal velocity and `B_n` unchanged).
   real(R8P), intent(in)  :: us(NV_MHD)  !< Star state of the side.
   real(R8P), intent(in)  :: vss(2)      !< Double-star tangential velocity.
   real(R8P), intent(in)  :: bss(2)      !< Double-star tangential field.
   real(R8P), intent(in)  :: e           !< Double-star energy.
   real(R8P), intent(out) :: uss(NV_MHD) !< Double-star state.
   !$acc routine seq
   !$omp declare target

   uss(1) = us(1)
   uss(2) = us(2)
   uss(3) = us(1) * vss(1)
   uss(4) = us(1) * vss(2)
   uss(5) = e
   uss(6) = us(6)
   uss(7) = bss(1)
   uss(8) = bss(2)
   endsubroutine double_star

   pure subroutine eglm_face_subsystem(bnL, bnR, psiL, psiR, bn, psi)
   !< Solve the linear `(B_n, psi)` subsystem of EGLM (flux `c_h (psi, B_n)`, `psi` in B units) exactly at the face:
   !< the characteristic fields `B_n + psi` (speed `+c_h`, from the left) and `B_n - psi` (`-c_h`, from the right).
   real(R8P), intent(in)  :: bnL, bnR   !< Normal fields of the two states.
   real(R8P), intent(in)  :: psiL, psiR !< EGLM scalars of the two states.
   real(R8P), intent(out) :: bn, psi    !< Face normal field and EGLM scalar.
   !$acc routine seq
   !$omp declare target

   bn  = 0.5_R8P * (bnL + bnR) - 0.5_R8P * (psiR - psiL)
   psi = 0.5_R8P * (psiL + psiR) - 0.5_R8P * (bnR - bnL)
   endsubroutine eglm_face_subsystem

   pure subroutine eglm_frame_state(gamma, pv, q, bn, w, u)
   !< Compute the frame primitive and conservative variables of an EGLM state whose normal field is set to `bn`: the
   !< MHD state, the energy without the cleaning energy `psi^2 / 2` (as `frame_state`).
   real(R8P),    intent(in)  :: gamma           !< Specific heats ratio.
   integer(I4P), intent(in)  :: pv(NV_MHD)      !< State in the frame order.
   real(R8P),    intent(in)  :: q(NV_MHD_EGLM)  !< Conservative variables.
   real(R8P),    intent(in)  :: bn              !< Normal field of the frame state.
   real(R8P),    intent(out) :: w(NV_MHD)       !< Frame primitive variables.
   real(R8P),    intent(out) :: u(NV_MHD)       !< Frame conservative variables.
   real(R8P)                 :: q8(NV_MHD)      !< MHD state.
   !$acc routine seq
   !$omp declare target

   q8 = q(1:NV_MHD)
   q8(IQ_RE) = q(IQ_RE) - 0.5_R8P * q(IQ_PSI)**2
   call frame_state(gamma=gamma, pv=pv, q=q8, bn=bn, w=w, u=u)
   endsubroutine eglm_frame_state

   pure function eglm_is_admissible(q) result(ok)
   !< Return true when an EGLM state has positive density and thermal pressure (energy without `psi^2 / 2`).
   real(R8P), intent(in) :: q(NV_MHD_EGLM) !< Conservative variables.
   logical               :: ok             !< Admissible state.
   real(R8P)             :: q8(NV_MHD)     !< MHD state.
   !$acc routine seq
   !$omp declare target

   q8 = q(1:NV_MHD)
   q8(IQ_RE) = q(IQ_RE) - 0.5_R8P * q(IQ_PSI)**2
   ok = is_admissible(q=q8)
   endfunction eglm_is_admissible

   pure function frame_fast_speed(gamma, w) result(cf)
   !< Return the fast magnetosonic speed along the frame normal (as `mhd_fast_speed`, from the frame primitive state).
   real(R8P), intent(in) :: gamma     !< Specific heats ratio.
   real(R8P), intent(in) :: w(NV_MHD) !< Frame primitive variables.
   real(R8P)             :: cf        !< Fast magnetosonic speed.
   real(R8P)             :: a2        !< Squared sound speed.
   real(R8P)             :: b2        !< Squared Alfven speed, |B|^2 / rho.
   real(R8P)             :: bt2       !< Transverse part of b2.
   !$acc routine seq
   !$omp declare target

   a2  = gamma * w(5) / w(1)
   b2  = mhd_sum3(w(6)**2, w(7)**2, w(8)**2) / w(1)
   bt2 = max(b2 - w(6)**2 / w(1), 0._R8P)
   cf  = sqrt(0.5_R8P * (a2 + b2 + sqrt((a2 - b2)**2 + 4._R8P * a2 * bt2)))
   endfunction frame_fast_speed

   pure subroutine frame_flux(w, u, f)
   !< Compute the physical flux along the frame normal: `(rho u_n, rho u u_n + p_t n - B B_n, (E + p_t) u_n - (u.B) B_n,
   !< B u_n - u B_n)`, the `B_n` component zero.
   real(R8P), intent(in)  :: w(NV_MHD) !< Frame primitive variables.
   real(R8P), intent(in)  :: u(NV_MHD) !< Frame conservative variables.
   real(R8P), intent(out) :: f(NV_MHD) !< Frame flux.
   real(R8P)              :: pt        !< Total pressure.
   real(R8P)              :: ub        !< u.B.
   !$acc routine seq
   !$omp declare target

   pt = w(5) + 0.5_R8P * mhd_sum3(w(6)**2, w(7)**2, w(8)**2)
   ub = mhd_sum3(w(2) * w(6), w(3) * w(7), w(4) * w(8))
   f(1) = u(1) * w(2)
   f(2) = u(2) * w(2) - w(6) * w(6) + pt
   f(3) = u(3) * w(2) - w(6) * w(7)
   f(4) = u(4) * w(2) - w(6) * w(8)
   f(5) = (u(5) + pt) * w(2) - ub * w(6)
   f(6) = 0._R8P
   f(7) = w(7) * w(2) - w(3) * w(6)
   f(8) = w(8) * w(2) - w(4) * w(6)
   endsubroutine frame_flux

   pure subroutine frame_hll(gamma, wL, uL, wR, uR, f)
   !< Compute the HLL flux in the frame.
   real(R8P), intent(in)  :: gamma      !< Specific heats ratio.
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
   !< Compute the HLLD flux in the frame (Miyoshi & Kusano 2005, section 5), `B_n = wL(6) = wR(6)`.
   real(R8P), intent(in)  :: gamma          !< Specific heats ratio.
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
   real(R8P)              :: uss(NV_MHD)    !< Double-star state of the side of the face.
   real(R8P)              :: sL, sR, sM     !< Outer and contact wave speeds.
   real(R8P)              :: ssL, ssR       !< Alfven wave speeds.
   real(R8P)              :: dL, dR         !< S - u_n of the two states.
   real(R8P)              :: den            !< Denominator of the contact speed and star pressure.
   real(R8P)              :: ptL, ptR       !< Total pressures.
   real(R8P)              :: pts            !< Star total pressure.
   real(R8P)              :: sqL, sqR       !< Square roots of the star densities.
   real(R8P)              :: sg             !< sign(B_n).
   real(R8P)              :: vss(2), bss(2) !< Double-star tangential velocity and field.
   real(R8P)              :: ubs, ubss      !< u.B of a star and of the double-star state.
   logical                :: okL, okR       !< Admissible star states.
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
   ptL = wL(5) + 0.5_R8P * mhd_sum3(wL(6)**2, wL(7)**2, wL(8)**2)
   ptR = wR(5) + 0.5_R8P * mhd_sum3(wR(6)**2, wR(7)**2, wR(8)**2)
   dL  = sL - wL(2)
   dR  = sR - wR(2)
   den = dR * wR(1) - dL * wL(1)
   sM  = (dR * wR(1) * wR(2) - dL * wL(1) * wL(2) - ptR + ptL) / den
   pts = (dR * wR(1) * ptL - dL * wL(1) * ptR + wL(1) * wR(1) * dR * dL * (wR(2) - wL(2))) / den
   call hlld_star(s=sL, sM=sM, pts=pts, w=wL, u=uL, us=usL, ok=okL)
   call hlld_star(s=sR, sM=sM, pts=pts, w=wR, u=uR, us=usR, ok=okR)
   if (okL .and. okR) then
      ssL = sM - abs(wL(6)) / sqrt(usL(1))
      ssR = sM + abs(wL(6)) / sqrt(usR(1))
      if (.not.(sL < ssL .and. ssL <= sM .and. sM <= ssR .and. ssR < sR)) okL = .false.
   endif
   if (.not.(okL .and. okR)) then
      fallback = .true.
      call hll_select(sL=sL, sR=sR, uL=uL, uR=uR, fL=fL, fR=fR, f=f)
      return
   endif
   if (ssL >= 0._R8P) then
      do k=1, NV_MHD
         f(k) = fL(k) + sL * (usL(k) - uL(k))
      enddo
   elseif (ssR <= 0._R8P) then
      do k=1, NV_MHD
         f(k) = fR(k) + sR * (usR(k) - uR(k))
      enddo
   else
      ! double-star states (reached only with B_n /= 0: B_n = 0 gives ssL = sM = ssR)
      sqL = sqrt(usL(1))
      sqR = sqrt(usR(1))
      sg  = sign(1._R8P, wL(6))
      do k=1, 2
         vss(k) = (sqL * usL(2+k) / usL(1) + sqR * usR(2+k) / usR(1) + (usR(6+k) - usL(6+k)) * sg) / (sqL + sqR)
         bss(k) = (sqL * usR(6+k) + sqR * usL(6+k) + sqL * sqR * (usR(2+k) / usR(1) - usL(2+k) / usL(1)) * sg) / &
                  (sqL + sqR)
      enddo
      ubss = mhd_sum3(sM * wL(6), vss(1) * bss(1), vss(2) * bss(2))
      if (sM >= 0._R8P) then
         ubs = mhd_sum3(sM * wL(6), usL(3) / usL(1) * usL(7), usL(4) / usL(1) * usL(8))
         call double_star(us=usL, vss=vss, bss=bss, e=usL(5) - sqL * (ubs - ubss) * sg, uss=uss)
         do k=1, NV_MHD
            f(k) = fL(k) + sL * (usL(k) - uL(k)) + ssL * (uss(k) - usL(k))
         enddo
      else
         ubs = mhd_sum3(sM * wL(6), usR(3) / usR(1) * usR(7), usR(4) / usR(1) * usR(8))
         call double_star(us=usR, vss=vss, bss=bss, e=usR(5) + sqR * (ubs - ubss) * sg, uss=uss)
         do k=1, NV_MHD
            f(k) = fR(k) + sR * (usR(k) - uR(k)) + ssR * (uss(k) - usR(k))
         enddo
      endif
   endif
   endsubroutine frame_hlld

   pure subroutine frame_llf(gamma, wL, uL, wR, uR, f)
   !< Compute the LLF (Rusanov) flux in the frame, speed `max(|u_n| + c_f)` of the two states.
   real(R8P), intent(in)  :: gamma      !< Specific heats ratio.
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
   real(R8P), intent(in)  :: gamma      !< Specific heats ratio.
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
   !< Compute the frame primitive and conservative variables of a state whose normal field is set to `bn`: the pressure
   !< of the state is kept, the energy recomputed with `bn` (`q` may be longer than `NV_MHD`: its leading entries).
   real(R8P),    intent(in)  :: gamma      !< Specific heats ratio.
   integer(I4P), intent(in)  :: pv(NV_MHD) !< State in the frame order.
   real(R8P),    intent(in)  :: q(NV_MHD)  !< Conservative variables.
   real(R8P),    intent(in)  :: bn         !< Normal field of the frame state.
   real(R8P),    intent(out) :: w(NV_MHD)  !< Frame primitive variables.
   real(R8P),    intent(out) :: u(NV_MHD)  !< Frame conservative variables.
   real(R8P)                 :: ek         !< Squared velocity, |u|^2.
   !$acc routine seq
   !$omp declare target

   w(1) = q(IQ_R)
   w(2) = q(pv(2)) / q(IQ_R)
   w(3) = q(pv(3)) / q(IQ_R)
   w(4) = q(pv(4)) / q(IQ_R)
   ek   = mhd_sum3(w(2)**2, w(3)**2, w(4)**2)
   w(5) = (gamma - 1._R8P) * (q(IQ_RE) - 0.5_R8P * q(IQ_R) * ek - 0.5_R8P * mhd_sum3(q(pv(6))**2, q(pv(7))**2, q(pv(8))**2))
   w(6) = bn
   w(7) = q(pv(7))
   w(8) = q(pv(8))
   u(1) = q(IQ_R)
   u(2) = q(pv(2))
   u(3) = q(pv(3))
   u(4) = q(pv(4))
   u(5) = w(5) / (gamma - 1._R8P) + 0.5_R8P * q(IQ_R) * ek + 0.5_R8P * mhd_sum3(bn**2, w(7)**2, w(8)**2)
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

   pure subroutine hlld_star(s, sM, pts, w, u, us, ok)
   !< Build the HLLD star state of a side (Miyoshi & Kusano 2005, eqs. 43-48); `ok` is false when its density or
   !< pressure (from its energy) is not positive. Near the degeneracy `rho d (s - S_M) = B_n^2` the outer transverse state is kept.
   real(R8P), intent(in)  :: s            !< Outer wave speed of the side.
   real(R8P), intent(in)  :: sM           !< Contact wave speed.
   real(R8P), intent(in)  :: pts          !< Star total pressure.
   real(R8P), intent(in)  :: w(NV_MHD)    !< Frame primitive variables of the side.
   real(R8P), intent(in)  :: u(NV_MHD)    !< Frame conservative variables of the side.
   real(R8P), intent(out) :: us(NV_MHD)   !< Star state.
   logical,   intent(out) :: ok           !< Admissible star state.
   real(R8P)              :: d            !< s - u_n.
   real(R8P)              :: rs           !< Star density.
   real(R8P)              :: e            !< rho d (s - S_M) - B_n^2.
   real(R8P)              :: rdd          !< rho d (s - S_M).
   real(R8P)              :: vs(2), bs(2) !< Star tangential velocity and field.
   real(R8P)              :: pt           !< Total pressure of the side.
   real(R8P)              :: ub, ubs      !< u.B of the side and of the star state.
   integer(I4P)           :: k            !< Counter.
   !$acc routine seq
   !$omp declare target

   d   = s - w(2)
   rs  = w(1) * d / (s - sM)
   rdd = w(1) * d * (s - sM)
   e   = rdd - w(6)**2
   if (abs(e) <= EPS_HLLD * (abs(rdd) + w(6)**2)) then
      do k=1, 2
         vs(k) = w(2+k)
         bs(k) = w(6+k)
      enddo
   else
      do k=1, 2
         vs(k) = w(2+k) - w(6) * w(6+k) * (sM - w(2)) / e
         bs(k) = w(6+k) * (w(1) * d**2 - w(6)**2) / e
      enddo
   endif
   pt  = w(5) + 0.5_R8P * mhd_sum3(w(6)**2, w(7)**2, w(8)**2)
   ub  = mhd_sum3(w(2) * w(6), w(3) * w(7), w(4) * w(8))
   ubs = mhd_sum3(sM * w(6), vs(1) * bs(1), vs(2) * bs(2))
   us(1) = rs
   us(2) = rs * sM
   us(3) = rs * vs(1)
   us(4) = rs * vs(2)
   us(5) = (d * u(5) - pt * w(2) + pts * sM + w(6) * (ub - ubs)) / (s - sM)
   us(6) = w(6)
   us(7) = bs(1)
   us(8) = bs(2)
   ! the pressure of the conservative star state (from its energy), not pts - |B*|^2/2: the two differ in the
   ! approximate solver, and the admissibility of the flux F + S (U* - U) depends on U*
   ok = rs > 0._R8P
   if (ok) ok = us(5) - 0.5_R8P * rs * mhd_sum3(sM**2, vs(1)**2, vs(2)**2) - &
                0.5_R8P * mhd_sum3(w(6)**2, bs(1)**2, bs(2)**2) > 0._R8P
   endsubroutine hlld_star

   pure function is_admissible(q) result(ok)
   !< Return true when a state (its leading `NV_MHD` entries) has positive density and pressure (the pressure up to the
   !< positive factor `gamma - 1`).
   real(R8P), intent(in) :: q(NV_MHD) !< Conservative variables.
   logical               :: ok        !< Admissible state.
   !$acc routine seq
   !$omp declare target

   ok = q(IQ_R) > 0._R8P
   if (ok) ok = q(IQ_RE) - 0.5_R8P * mhd_sum3(q(IQ_RU)**2, q(IQ_RV)**2, q(IQ_RW)**2) / q(IQ_R) - &
                0.5_R8P * mhd_sum3(q(IQ_BX)**2, q(IQ_BY)**2, q(IQ_BZ)**2) > 0._R8P
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
   !< Return the Wu (2018) Lax-Friedrichs speed of two states in direction `d` (see the module header).
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

   pure subroutine scatter_eglm_flux(ch, pv, fr, bn, psi, rhoL, rhoR, f)
   !< Return the EGLM face flux in the global order: the frame flux scattered, the `B_n` flux `c_h psi~`, the `psi`
   !< flux `c_h B~_n`, the energy flux plus `c_h psi~ B~_n` and the cleaning energy carried by the mass flux,
   !< `F_rho psi~^2 / (2 rho_up)`, `rho_up` the density of the side the mass flows from.
   real(R8P),    intent(in)  :: ch             !< GLM cleaning speed.
   integer(I4P), intent(in)  :: pv(NV_MHD)     !< State in the frame order.
   real(R8P),    intent(in)  :: fr(NV_MHD)     !< Frame flux.
   real(R8P),    intent(in)  :: bn, psi        !< Face normal field and EGLM scalar.
   real(R8P),    intent(in)  :: rhoL, rhoR     !< Densities of the two states.
   real(R8P),    intent(out) :: f(NV_MHD_EGLM) !< Flux.
   real(R8P)                 :: rho_up         !< Upwind density.
   integer(I4P)              :: k              !< Counter.
   !$acc routine seq
   !$omp declare target

   do k=1, NV_MHD
      f(pv(k)) = fr(k)
   enddo
   rho_up = merge(rhoL, rhoR, fr(1) >= 0._R8P)
   f(pv(6))  = ch * psi
   f(IQ_PSI) = ch * bn
   f(IQ_RE)  = f(IQ_RE) + ch * psi * bn + fr(1) * (0.5_R8P * psi**2 / rho_up)
   endsubroutine scatter_eglm_flux
endmodule adam_flume_mhd_riemann_library
