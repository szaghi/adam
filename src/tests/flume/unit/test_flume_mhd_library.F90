!< Unit test MV-0 of the FLUME MHD library: eigensystem, fluxes, conversions, face average, split consistency.
program test_flume_mhd_library
!< Unit test MV-0 of the FLUME MHD library: eigensystem, fluxes, conversions, face average, split consistency.
!<
!< **Why this test exists** (issue #41, sections 3.3 and 9, MV-0). The MHD eigensystem is the most error-prone piece of
!< M2: a sign or a normalisation slip in one of the 49 entries of the Roe-Balsara core is invisible in 1-D x-direction
!< runs and in smooth flows, and its degeneracies (field along the normal, no normal field, the triple umbilic) are
!< where codes break. It is pinned here, for both variants (without cleaning, `nv = 8`; GLM, `nv = 9`) and in all
!< three directions, before any solver uses it.
!<
!< **What it pins**, on N deterministic random admissible states (plasma beta from ~0.01 to ~100) and each direction d:
!< 1. `L R = I` of the full systems (8x8, 9x9);
!< 2. the 7x7 core against a central-difference flux Jacobian taken with `B_n` (and `psi`) held fixed:
!<    `L7 A7 R7 = Lambda7` (truncation-limited tolerance), which checks every entry of the core and its eigenvalues;
!< 3. the GLM `(B_n, psi)` block: `L2 A2 R2 = diag(-c_h, c_h)` exactly and the block structure (zeros outside it);
!<    without cleaning the `B_n` flux is exactly zero; the GLM flux equals the plain one on the core variables;
!< 4. the round trip primitive -> conservative -> auxiliary, and the fast speed of the eigenvalues against
!<    `mhd_fast_speed`;
!< 5. the face average of a state with itself returns the state;
!< 6. cyclic invariance: flux and eigenvectors in direction d of a state equal those in direction 1 of the cyclically
!<    permuted state (velocity and field permuted), components permuted back;
!< 7. split + back-projection consistency: a uniform stencil returns the physical flux, S = 1..S_MAX, both variants,
!<    characteristic and conservative;
!< 8. the degenerate states: B = 0 (hydrodynamic limit), B along one axis (B_t = 0 in that direction, B_n = 0 in the
!<    others), the triple umbilic, and the transverse-field threshold `EPS_BT` straddled by one ulp, through checks 1-3.

use :: adam_flume_mhd_library, only : EPS_BT, mhd_conservative_to_auxiliary, mhd_eigenvalues, mhd_eigenvectors,    &
                                      mhd_face_average, mhd_face_flux_back_projection, mhd_face_split_fluxes,      &
                                      mhd_fast_speed, mhd_flux, mhd_glm_eigenvalues, mhd_glm_eigenvectors,         &
                                      mhd_glm_face_flux_back_projection, mhd_glm_face_split_fluxes, mhd_glm_flux, &
                                      mhd_primitive_to_conservative
use :: adam_flume_parameters,  only : IA_A, IA_BX, IA_H, IA_P, IA_R, IA_U, IQ_BX, IQ_PSI, IQ_R, IQ_RE, IQ_RU, &
                                      NV_AUX_MHD, NV_MHD, NV_MHD_GLM, S_MAX
use :: penf,                   only : I4P, R8P, str

implicit none

integer(I4P), parameter :: N=5000_I4P           !< Random states number.
integer(I4P), parameter :: NC=8_I4P             !< Checks number.
real(R8P),    parameter :: GAMMA=5._R8P/3._R8P  !< Specific heats ratio.
real(R8P),    parameter :: R=1._R8P             !< Gas constant.
real(R8P),    parameter :: TOL_EXACT=1.e-11_R8P !< Tolerance of the exact identities (relative).
real(R8P),    parameter :: TOL_FD=1.e-6_R8P     !< Tolerance of the finite-difference Jacobian checks (relative).
real(R8P)               :: err(NC)              !< Maximum error of each check.
real(R8P)               :: prim(8)              !< Primitive state (r, u, v, w, p, bx, by, bz).
real(R8P)               :: psi, ch              !< GLM scalar and cleaning speed.
real(R8P)               :: a                    !< Sound speed.
integer(I4P)            :: seed(64)             !< Random generator seed.
integer(I4P)            :: n_, c, ns            !< Counters, seed size.
logical                 :: test_passed          !< Aggregate pass flag.
character(len=48)       :: check_name(NC)       !< Checks names.

check_name = ['L R = I (8x8, 9x9)                              ', &
              'L7 A7_fd R7 = Lambda7 (7x7 core, fin. diff.)    ', &
              'GLM (B_n, psi) block, B_n flux without cleaning ', &
              'primitive -> conservative -> auxiliary, c_f     ', &
              'face average of a state with itself             ', &
              'cyclic invariance of fluxes, eigenvectors       ', &
              'uniform stencil split -> physical flux          ', &
              'degenerate states (checks 1-3 at degeneracies)  ']
err = 0._R8P
call random_seed(size=ns)
if (ns > size(seed)) error stop 'random seed larger than expected'
seed = [(20260925_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns))
do n_=1, N
   prim(1) = 0.1_R8P + 9.9_R8P * uniform()
   prim(5) = 0.1_R8P + 9.9_R8P * uniform()
   a = sqrt(GAMMA * prim(5) / prim(1))
   do c=2, 4
      prim(c) = 3._R8P * a * (2._R8P * uniform() - 1._R8P)
   enddo
   do c=6, 8
      prim(c) = 10._R8P**(2._R8P * uniform() - 1._R8P) * sqrt(prim(1)) * a * (2._R8P * uniform() - 1._R8P)
   enddo
   psi = 2._R8P * uniform() - 1._R8P
   ch  = 0.5_R8P + 4.5_R8P * uniform()
   call check_state(prim=prim, psi=psi, ch=ch, e=err(1:3))
   call check_conversions(prim=prim, e=err(4))
   call check_average(prim=prim, e=err(5))
   call check_cyclic(prim=prim, psi=psi, ch=ch, e=err(6))
   call check_split(prim=prim, psi=psi, ch=ch, e=err(7))
enddo
call check_degenerate(e=err(8))

test_passed = .true.
do c=1, NC
   if (((c == 2 .or. c == 8) .and. err(c) > TOL_FD) .or. (c /= 2 .and. c /= 8 .and. err(c) > TOL_EXACT)) then
      print '(A)', 'FAIL: '//check_name(c)//' max error '//trim(str(err(c)))
      test_passed = .false.
   else
      print '(A)', 'PASS: '//check_name(c)//' max error '//trim(str(err(c)))
   endif
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: flume mhd library ('//trim(str(N))//' states x 3 directions x 2 variants)'
else
   print '(A)', 'TEST FAILED: flume mhd library'
   error stop 1
endif

contains
   function uniform() result(x)
   !< Return a pseudo-random number in [0, 1) from the intrinsic generator (seeded once, deterministic).
   real(R8P) :: x !< Random number.

   call random_number(x)
   endfunction uniform

   subroutine state(prim, psi, q8, q9, qa)
   !< Return the conservative states of both variants and the auxiliary variables of a primitive state.
   real(R8P), intent(in)  :: prim(8)        !< Primitive state.
   real(R8P), intent(in)  :: psi            !< GLM scalar.
   real(R8P), intent(out) :: q8(NV_MHD)     !< Conservative variables, no cleaning.
   real(R8P), intent(out) :: q9(NV_MHD_GLM) !< Conservative variables, GLM.
   real(R8P), intent(out) :: qa(NV_AUX_MHD) !< Auxiliary variables.

   call mhd_primitive_to_conservative(gamma=GAMMA, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), &
                                      bx=prim(6), by=prim(7), bz=prim(8), q=q8)
   q9(1:NV_MHD) = q8
   q9(IQ_PSI)   = psi
   call mhd_conservative_to_auxiliary(gamma=GAMMA, R=R, q=q8, qa=qa)
   endsubroutine state

   subroutine core_map(d, mp)
   !< Return the full-state indexes of the 7 core variables in direction d (independent restatement of the layout).
   integer(I4P), intent(in)  :: d      !< Direction.
   integer(I4P), intent(out) :: mp(7)  !< Core variables indexes.
   integer(I4P)              :: d1, d2 !< Tangential directions.

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   mp = [IQ_R, IQ_RU+d-1, IQ_RU+d1-1, IQ_RU+d2-1, IQ_RE, IQ_BX+d1-1, IQ_BX+d2-1]
   endsubroutine core_map

   function identity_error(m, nv) result(e)
   !< Return `max|m - I|`.
   integer(I4P), intent(in) :: nv       !< Size.
   real(R8P),    intent(in) :: m(nv,nv) !< Matrix.
   real(R8P)                :: e        !< Error.
   integer(I4P)             :: i, j     !< Counters.

   e = 0._R8P
   do j=1, nv
      do i=1, nv
         e = max(e, abs(m(i,j) - merge(1._R8P, 0._R8P, i == j)))
      enddo
   enddo
   endfunction identity_error

   subroutine check_state(prim, psi, ch, e)
   !< Checks 1-3 on one state, all directions, both variants; errors accumulated into e(1:3).
   real(R8P), intent(in)    :: prim(8)                                          !< Primitive state.
   real(R8P), intent(in)    :: psi, ch                                          !< GLM scalar and cleaning speed.
   real(R8P), intent(inout) :: e(3)                                             !< Maximum errors of checks 1-3.
   real(R8P)                :: q8(NV_MHD), q9(NV_MHD_GLM)                       !< Conservative variables.
   real(R8P)                :: qa(NV_AUX_MHD)                                   !< Auxiliary variables.
   real(R8P)                :: el8(NV_MHD,NV_MHD), er8(NV_MHD,NV_MHD)           !< Eigenvectors, no cleaning.
   real(R8P)                :: el9(NV_MHD_GLM,NV_MHD_GLM)                       !< Left eigenvectors, GLM.
   real(R8P)                :: er9(NV_MHD_GLM,NV_MHD_GLM)                       !< Right eigenvectors, GLM.
   real(R8P)                :: lam8(NV_MHD), lam9(NV_MHD_GLM)                   !< Eigenvalues.
   real(R8P)                :: f8(NV_MHD), f9(NV_MHD_GLM)                       !< Fluxes.
   real(R8P)                :: A7(7,7), D7(7,7)                                 !< Core Jacobian, L7 A7 R7 - Lambda7.
   real(R8P)                :: A2(2,2), D2(2,2)                                 !< GLM block Jacobian, residual.
   real(R8P)                :: qp(NV_MHD_GLM), qm(NV_MHD_GLM)                   !< Perturbed states.
   real(R8P)                :: qap(NV_AUX_MHD), qam(NV_AUX_MHD)                 !< Perturbed auxiliary states.
   real(R8P)                :: fp(NV_MHD_GLM), fm(NV_MHD_GLM)                   !< Perturbed fluxes.
   real(R8P)                :: h                                                !< Perturbation.
   integer(I4P)             :: mp(7), bp(2)                                     !< Core and GLM block indexes.
   integer(I4P)             :: d, c, i                                          !< Counters.

   call state(prim=prim, psi=psi, q8=q8, q9=q9, qa=qa)
   do d=1, 3
      call core_map(d=d, mp=mp)
      bp = [IQ_BX+d-1, IQ_PSI]
      call mhd_eigenvectors(gamma=GAMMA, d=d, qa=qa, el=el8, er=er8)
      call mhd_glm_eigenvectors(ch=ch, gamma=GAMMA, d=d, qa=qa, el=el9, er=er9)
      call mhd_eigenvalues(d=d, qa=qa, lambda=lam8)
      call mhd_glm_eigenvalues(ch=ch, d=d, qa=qa, lambda=lam9)
      ! 1. L R = I
      e(1) = max(e(1), identity_error(matmul(el8, er8), NV_MHD), identity_error(matmul(el9, er9), NV_MHD_GLM))
      ! 2. core against the finite-difference Jacobian, B_n and psi held fixed
      do c=1, 7
         h = 1.e-6_R8P * max(1._R8P, abs(q9(mp(c))))
         qp = q9 ; qp(mp(c)) = qp(mp(c)) + h
         qm = q9 ; qm(mp(c)) = qm(mp(c)) - h
         call mhd_conservative_to_auxiliary(gamma=GAMMA, R=R, q=qp, qa=qap)
         call mhd_conservative_to_auxiliary(gamma=GAMMA, R=R, q=qm, qa=qam)
         call mhd_glm_flux(ch=ch, d=d, q=qp, qa=qap, f=fp)
         call mhd_glm_flux(ch=ch, d=d, q=qm, qa=qam, f=fm)
         A7(:,c) = (fp(mp) - fm(mp)) / (2._R8P * h)
      enddo
      D7 = matmul(el9(1:7,mp), matmul(A7, er9(mp,1:7)))
      do i=1, 7
         D7(i,i) = D7(i,i) - lam9(i)
      enddo
      e(2) = max(e(2), maxval(abs(D7)) / maxval(abs(lam9(1:7))))
      D7 = matmul(el8(1:7,mp), matmul(A7, er8(mp,1:7)))
      do i=1, 7
         D7(i,i) = D7(i,i) - lam8(i)
      enddo
      e(2) = max(e(2), maxval(abs(D7)) / maxval(abs(lam8(1:7))))
      ! 3. GLM block (linear: exact) and its structure; B_n flux without cleaning
      A2 = reshape([0._R8P, ch * ch, 1._R8P, 0._R8P], [2, 2])
      D2 = matmul(el9(8:9,bp), matmul(A2, er9(bp,8:9)))
      D2(1,1) = D2(1,1) + ch
      D2(2,2) = D2(2,2) - ch
      e(3) = max(e(3), maxval(abs(D2)) / ch, maxval(abs(el9(8:9,mp))), maxval(abs(er9(mp,8:9))), &
                 maxval(abs(el9(1:7,bp))), maxval(abs(er9(bp,1:7))), abs(lam9(8) + ch), abs(lam9(9) - ch))
      call mhd_flux(d=d, q=q8, qa=qa, f=f8)
      call mhd_glm_flux(ch=ch, d=d, q=q9, qa=qa, f=f9)
      e(3) = max(e(3), abs(f8(IQ_BX+d-1)), abs(lam8(8)), abs(f9(IQ_BX+d-1) - psi),                   &
                 abs(f9(IQ_PSI) - ch * ch * q9(IQ_BX+d-1)) / (ch * ch), maxval(abs(f9(mp) - f8(mp))) / &
                 maxval(abs(f8)))
   enddo
   endsubroutine check_state

   subroutine check_conversions(prim, e)
   !< Check the round trip primitive -> conservative -> auxiliary and the fast speed of the eigenvalues.
   real(R8P), intent(in)    :: prim(8)                    !< Primitive state.
   real(R8P), intent(inout) :: e                          !< Maximum error.
   real(R8P)                :: q8(NV_MHD), q9(NV_MHD_GLM) !< Conservative variables.
   real(R8P)                :: qa(NV_AUX_MHD)             !< Auxiliary variables.
   real(R8P)                :: lam(NV_MHD)                !< Eigenvalues.
   real(R8P)                :: back(8)                    !< Recovered primitive state.
   real(R8P)                :: cf                         !< Fast speed.
   integer(I4P)             :: d                          !< Direction.

   call state(prim=prim, psi=0._R8P, q8=q8, q9=q9, qa=qa)
   back = [qa(IA_R), qa(IA_U), qa(IA_U+1), qa(IA_U+2), qa(IA_P), qa(IA_BX), qa(IA_BX+1), qa(IA_BX+2)]
   e = max(e, maxval(abs(back - prim) / max(1._R8P, abs(prim))))
   do d=1, 3
      call mhd_eigenvalues(d=d, qa=qa, lambda=lam)
      cf = mhd_fast_speed(d=d, qa=qa)
      e = max(e, abs((lam(7) - qa(IA_U+d-1)) - cf) / cf, abs((qa(IA_U+d-1) - lam(1)) - cf) / cf)
   enddo
   endsubroutine check_conversions

   subroutine check_average(prim, e)
   !< Check that the face average of a state with itself returns the state.
   real(R8P), intent(in)    :: prim(8)                    !< Primitive state.
   real(R8P), intent(inout) :: e                          !< Maximum error.
   real(R8P)                :: q8(NV_MHD), q9(NV_MHD_GLM) !< Conservative variables.
   real(R8P)                :: qa(NV_AUX_MHD)             !< Auxiliary variables.
   real(R8P)                :: avg(NV_AUX_MHD)            !< Average.
   integer(I4P)             :: idx(9)                     !< Compared auxiliary variables.

   call state(prim=prim, psi=0._R8P, q8=q8, q9=q9, qa=qa)
   idx = [IA_R, IA_U, IA_U+1, IA_U+2, IA_P, IA_H, IA_A, IA_BX, IA_BX+2]
   call mhd_face_average(gamma=GAMMA, qaL=qa, qaR=qa, avg=avg)
   e = max(e, maxval(abs(avg(idx) - qa(idx)) / max(1._R8P, abs(qa(idx)))))
   endsubroutine check_average

   subroutine check_cyclic(prim, psi, ch, e)
   !< Check that direction d of a state equals direction 1 of the cyclically permuted state, permuted back.
   real(R8P), intent(in)    :: prim(8)                                              !< Primitive state.
   real(R8P), intent(in)    :: psi, ch                                              !< GLM scalar, cleaning speed.
   real(R8P), intent(inout) :: e                                                    !< Maximum error.
   real(R8P)                :: pp(8)                                                !< Permuted primitive state.
   real(R8P)                :: q8(NV_MHD), q9(NV_MHD_GLM)                           !< Original states.
   real(R8P)                :: q8p(NV_MHD), q9p(NV_MHD_GLM)                         !< Permuted states.
   real(R8P)                :: qa(NV_AUX_MHD), qap(NV_AUX_MHD)                      !< Auxiliary variables.
   real(R8P)                :: el(NV_MHD_GLM,NV_MHD_GLM), er(NV_MHD_GLM,NV_MHD_GLM)   !< Original eigenvectors.
   real(R8P)                :: elp(NV_MHD_GLM,NV_MHD_GLM), erp(NV_MHD_GLM,NV_MHD_GLM) !< Permuted eigenvectors.
   real(R8P)                :: f(NV_MHD_GLM), fp(NV_MHD_GLM)                        !< Fluxes.
   integer(I4P)             :: perm(NV_MHD_GLM)                                     !< Permuted state components.
   integer(I4P)             :: d                                                    !< Direction.

   call state(prim=prim, psi=psi, q8=q8, q9=q9, qa=qa)
   do d=1, 3
      pp = [prim(1), cshift(prim(2:4), d - 1), prim(5), cshift(prim(6:8), d - 1)]
      call state(prim=pp, psi=psi, q8=q8p, q9=q9p, qa=qap)
      perm = [1, 1+mod(d-1,3)+1, 1+mod(d,3)+1, 1+mod(d+1,3)+1, 5, 5+mod(d-1,3)+1, 5+mod(d,3)+1, 5+mod(d+1,3)+1, 9]
      call mhd_glm_flux(ch=ch, d=d, q=q9, qa=qa, f=f)
      call mhd_glm_flux(ch=ch, d=1, q=q9p, qa=qap, f=fp)
      call mhd_glm_eigenvectors(ch=ch, gamma=GAMMA, d=d, qa=qa, el=el, er=er)
      call mhd_glm_eigenvectors(ch=ch, gamma=GAMMA, d=1, qa=qap, el=elp, er=erp)
      e = max(e, maxval(abs(f(perm) - fp)) / maxval(abs(fp)))
      e = max(e, maxval(abs(er(perm,:) - erp)) / maxval(abs(erp)))
      e = max(e, maxval(abs(el(:,perm) - elp)) / maxval(abs(elp)))
   enddo
   endsubroutine check_cyclic

   subroutine check_split(prim, psi, ch, e)
   !< Check that a uniform stencil, split and back-projected, returns the physical flux (both variants, S = 1..S_MAX).
   real(R8P), intent(in)    :: prim(8)                                        !< Primitive state.
   real(R8P), intent(in)    :: psi, ch                                        !< GLM scalar and cleaning speed.
   real(R8P), intent(inout) :: e                                              !< Maximum error.
   real(R8P)                :: q8(NV_MHD), q9(NV_MHD_GLM)                     !< Conservative variables.
   real(R8P)                :: qa(NV_AUX_MHD)                                 !< Auxiliary variables.
   real(R8P)                :: qs8(NV_MHD,1-S_MAX:S_MAX)                      !< Stencil, no cleaning.
   real(R8P)                :: qs9(NV_MHD_GLM,1-S_MAX:S_MAX)                  !< Stencil, GLM.
   real(R8P)                :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)                  !< Stencil auxiliary variables.
   real(R8P)                :: fs8(2,1-S_MAX:S_MAX-1,NV_MHD)                  !< Split fields, no cleaning.
   real(R8P)                :: fs9(2,1-S_MAX:S_MAX-1,NV_MHD_GLM)              !< Split fields, GLM.
   real(R8P)                :: er8(NV_MHD,NV_MHD), er9(NV_MHD_GLM,NV_MHD_GLM) !< Right eigenvectors.
   real(R8P)                :: fl8(NV_MHD), fl9(NV_MHD_GLM)                   !< Face fluxes.
   real(R8P)                :: f8(NV_MHD), f9(NV_MHD_GLM)                     !< Physical fluxes.
   integer(I4P)             :: d, S, m, variant                               !< Counters.

   call state(prim=prim, psi=psi, q8=q8, q9=q9, qa=qa)
   do m=1-S_MAX, S_MAX
      qs8(:,m) = q8
      qs9(:,m) = q9
      qas(:,m) = qa
   enddo
   do d=1, 3
      call mhd_flux(d=d, q=q8, qa=qa, f=f8)
      call mhd_glm_flux(ch=ch, d=d, q=q9, qa=qa, f=f9)
      do variant=1, 2
         do S=1, S_MAX
            call mhd_face_split_fluxes(gamma=GAMMA, d=d, S=S, is_characteristic=(variant == 1), qs=qs8, qas=qas, &
                                       fsplit=fs8, er=er8)
            call mhd_face_flux_back_projection(is_characteristic=(variant == 1), er=er8, vr=fs8(:,0,:), flux=fl8)
            call mhd_glm_face_split_fluxes(ch=ch, gamma=GAMMA, d=d, S=S, is_characteristic=(variant == 1), qs=qs9, &
                                           qas=qas, fsplit=fs9, er=er9)
            call mhd_glm_face_flux_back_projection(is_characteristic=(variant == 1), er=er9, vr=fs9(:,0,:), flux=fl9)
            e = max(e, maxval(abs(fl8 - f8)) / maxval(abs(f8)), maxval(abs(fl9 - f9)) / maxval(abs(f9)))
         enddo
      enddo
   enddo
   endsubroutine check_split

   subroutine check_degenerate(e)
   !< Checks 1-3 at the degenerate states, and across the transverse-field threshold (one ulp either side).
   real(R8P), intent(inout) :: e      !< Maximum error (over checks 1-3).
   real(R8P)                :: ed(3)  !< Errors of checks 1-3.
   real(R8P)                :: p(8)   !< Primitive state.
   real(R8P)                :: a_     !< Sound speed of the reference state (rho = p = 1).
   real(R8P)                :: t      !< Threshold of the transverse field.
   integer(I4P)             :: k      !< Counter.

   ed = 0._R8P
   a_ = sqrt(GAMMA)
   ! B = 0: hydrodynamic limit
   p = [1._R8P, 0.3_R8P, -0.2_R8P, 0.1_R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P]
   call check_state(prim=p, psi=0.1_R8P, ch=2._R8P, e=ed)
   ! field along x only (B_t = 0 in x, B_n = 0 in y and z), weak and strong
   do k=1, 2
      p = [1._R8P, 0.3_R8P, -0.2_R8P, 0.1_R8P, 1._R8P, merge(0.5_R8P, 3._R8P, k == 1), 0._R8P, 0._R8P]
      call check_state(prim=p, psi=0.1_R8P, ch=2._R8P, e=ed)
   enddo
   ! triple umbilic in x: B along x with b_n^2 / rho = a^2
   p = [1._R8P, 0.3_R8P, -0.2_R8P, 0.1_R8P, 1._R8P, a_, 0._R8P, 0._R8P]
   call check_state(prim=p, psi=0.1_R8P, ch=2._R8P, e=ed)
   ! transverse field straddling the EPS_BT threshold by one ulp (direction x, B_t = by; |B| = 2 to the last bit)
   t = EPS_BT * 2._R8P
   do k=-1, 1, 2
      p = [1._R8P, 0.3_R8P, -0.2_R8P, 0.1_R8P, 1._R8P, 2._R8P, nearest(t, real(k, R8P)), 0._R8P]
      call check_state(prim=p, psi=0.1_R8P, ch=2._R8P, e=ed)
   enddo
   do k=1, 3
      if (ed(k) /= ed(k)) ed(k) = huge(1._R8P)
   enddo
   e = max(e, maxval(ed))
   print '(A)', 'degenerate states: L R = I '//trim(str(ed(1)))//', core fin. diff. '//trim(str(ed(2)))// &
                ', GLM block '//trim(str(ed(3)))
   endsubroutine check_degenerate
endprogram test_flume_mhd_library
