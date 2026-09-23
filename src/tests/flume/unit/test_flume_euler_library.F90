!< Unit test V0 of the FLUME Euler library: eigensystem, flux, conversions, Roe average, split consistency.
program test_flume_euler_library
!< Unit test V0 of the FLUME Euler library: eigensystem, flux, conversions, Roe average, split consistency.
!<
!< **Why this test exists** (issue #35, sections 1.3 and 11). The Euler eigensystem of CHASE is applied transposed and
!< is singular in y/z; NASTO's FNL y-direction eigenvectors are suspect (#36). Both went unnoticed because only 1-D
!< x-direction problems were ever run. FLUME's eigensystem is therefore validated here, in all three directions,
!< before any solver uses it.
!<
!< **What it pins**, on N deterministic random admissible states and each direction d:
!< 1. `L R = I`;
!< 2. the exact homogeneity identity `L F(q) = Lambda L q` (the Euler flux is homogeneous of degree one, `F = A q`,
!<    and `A = R Lambda L`), which checks the whole eigen-decomposition to round-off;
!< 3. a finite-difference Jacobian check `L A_fd R = Lambda` (truncation-limited tolerance);
!< 4. the round trip primitive -> conservative -> auxiliary;
!< 5. the Roe average of a state with itself returns the state;
!< 6. cyclic invariance: flux and eigenvectors in direction d of a state equal those in direction 1 of the cyclically
!<    permuted state, components permuted back (the check that exposes transposed or mis-rotated eigenvectors);
!< 7. split + back-projection consistency: a uniform stencil returns the physical flux, for S = 1..S_MAX and for both
!<    characteristic and conservative variables.

use :: adam_flume_euler_library, only : compute_eigenvalues, compute_eigenvectors, compute_face_flux_back_projection, &
                                        compute_face_split_fluxes, compute_flux, compute_roe_average,                 &
                                        conservative_to_auxiliary, primitive_to_conservative
use :: adam_flume_parameters,    only : IA_A, IA_H, IA_P, IA_R, IA_U, IA_V, IA_W, NV_AUX, NV_EULER, S_MAX
use :: penf,                     only : I4P, R8P, str

implicit none

integer(I4P), parameter :: N=10000_I4P            !< Random states number.
real(R8P),    parameter :: GAMMA=1.4_R8P          !< Specific heats ratio.
real(R8P),    parameter :: R=287.05_R8P           !< Gas constant.
real(R8P),    parameter :: TOL_EXACT=1.e-12_R8P   !< Tolerance of the exact identities (relative).
real(R8P),    parameter :: TOL_FD=1.e-6_R8P       !< Tolerance of the finite-difference Jacobian check (relative).
real(R8P)               :: err(7)                 !< Maximum error of each check.
real(R8P)               :: prim(5)                !< Primitive state (r, u, v, w, p).
real(R8P)               :: q(NV_EULER)            !< Conservative variables.
real(R8P)               :: qa(NV_AUX)             !< Auxiliary variables.
real(R8P)               :: el(NV_EULER,NV_EULER)  !< Left eigenvectors.
real(R8P)               :: er(NV_EULER,NV_EULER)  !< Right eigenvectors.
real(R8P)               :: lambda(NV_EULER)       !< Eigenvalues.
integer(I4P)            :: seed(64)               !< Random generator seed.
integer(I4P)            :: n_, d, ns              !< Counters, seed size.
logical                 :: test_passed            !< Aggregate pass flag.
character(len=40)       :: check_name(7)          !< Checks names.

check_name = ['L R = I                                 ', &
              'L F(q) = Lambda L q (homogeneity)       ', &
              'L A_fd R = Lambda (finite differences)  ', &
              'primitive -> conservative -> auxiliary  ', &
              'Roe average of a state with itself      ', &
              'cyclic invariance of flux, eigenvectors ', &
              'uniform stencil split -> physical flux  ']
err = 0._R8P
call random_seed(size=ns)
if (ns > size(seed)) error stop 'random seed larger than expected'
seed = [(20260923_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns))
do n_=1, N
   prim(1) = 0.1_R8P + 9.9_R8P * uniform()
   prim(5) = 0.1_R8P + 9.9_R8P * uniform()
   prim(2) = 3._R8P * sqrt(GAMMA * prim(5) / prim(1)) * (2._R8P * uniform() - 1._R8P)
   prim(3) = 3._R8P * sqrt(GAMMA * prim(5) / prim(1)) * (2._R8P * uniform() - 1._R8P)
   prim(4) = 3._R8P * sqrt(GAMMA * prim(5) / prim(1)) * (2._R8P * uniform() - 1._R8P)
   call primitive_to_conservative(gamma=GAMMA, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), q=q)
   call conservative_to_auxiliary(gamma=GAMMA, R=R, q=q, qa=qa)
   err(4) = max(err(4), maxval(abs([qa(IA_R), qa(IA_U), qa(IA_V), qa(IA_W), qa(IA_P)] - prim) / &
                        max(1._R8P, abs(prim))))
   call check_roe(qa=qa, e=err(5))
   do d=1, 3
      call compute_eigenvectors(gamma=GAMMA, d=d, qa=qa, el=el, er=er)
      call compute_eigenvalues(d=d, qa=qa, lambda=lambda)
      err(1) = max(err(1), identity_error(matmul(el, er)))
      call check_homogeneity(d=d, q=q, qa=qa, el=el, lambda=lambda, e=err(2))
      call check_fd_jacobian(d=d, q=q, el=el, er=er, lambda=lambda, e=err(3))
      call check_cyclic(d=d, prim=prim, e=err(6))
      call check_split(d=d, q=q, qa=qa, e=err(7))
   enddo
enddo

test_passed = .true.
do n_=1, 7
   if ((n_ == 3 .and. err(n_) > TOL_FD) .or. (n_ /= 3 .and. err(n_) > TOL_EXACT)) then
      print '(A)', 'FAIL: '//check_name(n_)//' max error '//trim(str(err(n_)))
      test_passed = .false.
   else
      print '(A)', 'PASS: '//check_name(n_)//' max error '//trim(str(err(n_)))
   endif
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: flume euler library ('//trim(str(N))//' states x 3 directions)'
else
   print '(A)', 'TEST FAILED: flume euler library'
   error stop 1
endif

contains
   function uniform() result(x)
   !< Return a pseudo-random number in [0, 1) from the intrinsic generator (seeded once, deterministic).
   real(R8P) :: x !< Random number.

   call random_number(x)
   endfunction uniform

   function identity_error(m) result(e)
   !< Return `max|m - I|`.
   real(R8P), intent(in) :: m(NV_EULER,NV_EULER)  !< Matrix.
   real(R8P)             :: e                     !< Error.
   real(R8P)             :: id(NV_EULER,NV_EULER) !< Identity.
   integer(I4P)          :: i                     !< Counter.

   id = 0._R8P
   do i=1, NV_EULER
      id(i,i) = 1._R8P
   enddo
   e = maxval(abs(m - id))
   endfunction identity_error

   subroutine check_homogeneity(d, q, qa, el, lambda, e)
   !< Check `L F(q) = Lambda L q`, relative to `max|Lambda L q|`.
   integer(I4P), intent(in)    :: d                     !< Direction.
   real(R8P),    intent(in)    :: q(NV_EULER)           !< Conservative variables.
   real(R8P),    intent(in)    :: qa(NV_AUX)            !< Auxiliary variables.
   real(R8P),    intent(in)    :: el(NV_EULER,NV_EULER) !< Left eigenvectors.
   real(R8P),    intent(in)    :: lambda(NV_EULER)      !< Eigenvalues.
   real(R8P),    intent(inout) :: e                     !< Maximum error.
   real(R8P)                   :: f(NV_EULER)           !< Physical flux.
   real(R8P)                   :: rhs(NV_EULER)         !< Lambda L q.

   call compute_flux(d=d, q=q, qa=qa, f=f)
   rhs = lambda * matmul(el, q)
   e = max(e, maxval(abs(matmul(el, f) - rhs)) / maxval(abs(rhs)))
   endsubroutine check_homogeneity

   subroutine check_fd_jacobian(d, q, el, er, lambda, e)
   !< Check `L A_fd R = Lambda` with a central-difference flux Jacobian, relative to `max|lambda|`.
   integer(I4P), intent(in)    :: d                     !< Direction.
   real(R8P),    intent(in)    :: q(NV_EULER)           !< Conservative variables.
   real(R8P),    intent(in)    :: el(NV_EULER,NV_EULER) !< Left eigenvectors.
   real(R8P),    intent(in)    :: er(NV_EULER,NV_EULER) !< Right eigenvectors.
   real(R8P),    intent(in)    :: lambda(NV_EULER)      !< Eigenvalues.
   real(R8P),    intent(inout) :: e                     !< Maximum error.
   real(R8P)                   :: A(NV_EULER,NV_EULER)  !< Finite-difference Jacobian.
   real(R8P)                   :: D_(NV_EULER,NV_EULER) !< L A R - Lambda.
   real(R8P)                   :: qp(NV_EULER)          !< Perturbed state, plus.
   real(R8P)                   :: qm(NV_EULER)          !< Perturbed state, minus.
   real(R8P)                   :: qap(NV_AUX)           !< Perturbed auxiliary state, plus.
   real(R8P)                   :: qam(NV_AUX)           !< Perturbed auxiliary state, minus.
   real(R8P)                   :: fp(NV_EULER)          !< Perturbed flux, plus.
   real(R8P)                   :: fm(NV_EULER)          !< Perturbed flux, minus.
   real(R8P)                   :: h                     !< Perturbation.
   integer(I4P)                :: c, i                  !< Counters.

   do c=1, NV_EULER
      h = 1.e-6_R8P * max(1._R8P, abs(q(c)))
      qp = q ; qp(c) = qp(c) + h
      qm = q ; qm(c) = qm(c) - h
      call conservative_to_auxiliary(gamma=GAMMA, R=R, q=qp, qa=qap)
      call conservative_to_auxiliary(gamma=GAMMA, R=R, q=qm, qa=qam)
      call compute_flux(d=d, q=qp, qa=qap, f=fp)
      call compute_flux(d=d, q=qm, qa=qam, f=fm)
      A(:,c) = (fp - fm) / (2._R8P * h)
   enddo
   D_ = matmul(el, matmul(A, er))
   do i=1, NV_EULER
      D_(i,i) = D_(i,i) - lambda(i)
   enddo
   e = max(e, maxval(abs(D_)) / maxval(abs(lambda)))
   endsubroutine check_fd_jacobian

   subroutine check_roe(qa, e)
   !< Check that the Roe average of a state with itself returns the state.
   real(R8P), intent(in)    :: qa(NV_AUX)  !< Auxiliary variables.
   real(R8P), intent(inout) :: e           !< Maximum error.
   real(R8P)                :: roe(NV_AUX) !< Roe average.
   integer(I4P)             :: idx(7)      !< Compared auxiliary variables.

   idx = [IA_R, IA_U, IA_V, IA_W, IA_H, IA_A, IA_P]
   call compute_roe_average(gamma=GAMMA, qaL=qa, qaR=qa, roe=roe)
   e = max(e, maxval(abs(roe(idx) - qa(idx)) / max(1._R8P, abs(qa(idx)))))
   endsubroutine check_roe

   subroutine check_cyclic(d, prim, e)
   !< Check that direction `d` of a state equals direction 1 of the cyclically permuted state, permuted back.
   integer(I4P), intent(in)    :: d                      !< Direction.
   real(R8P),    intent(in)    :: prim(5)                !< Primitive state.
   real(R8P),    intent(inout) :: e                      !< Maximum error.
   real(R8P)                   :: q(NV_EULER)            !< Original state.
   real(R8P)                   :: qa(NV_AUX)             !< Original auxiliary state.
   real(R8P)                   :: qp(NV_EULER)           !< Permuted state.
   real(R8P)                   :: qap(NV_AUX)            !< Permuted auxiliary state.
   real(R8P)                   :: el(NV_EULER,NV_EULER)  !< Original left eigenvectors.
   real(R8P)                   :: er(NV_EULER,NV_EULER)  !< Original right eigenvectors.
   real(R8P)                   :: elp(NV_EULER,NV_EULER) !< Permuted left eigenvectors.
   real(R8P)                   :: erp(NV_EULER,NV_EULER) !< Permuted right eigenvectors.
   real(R8P)                   :: f(NV_EULER)            !< Original flux.
   real(R8P)                   :: fp(NV_EULER)           !< Permuted flux.
   real(R8P)                   :: vel(3)                 !< Velocity.
   integer(I4P)                :: perm(NV_EULER)         !< Original components of the permuted state.

   vel = prim(2:4)
   call primitive_to_conservative(gamma=GAMMA, r=prim(1), u=vel(1), v=vel(2), w=vel(3), p=prim(5), q=q)
   call conservative_to_auxiliary(gamma=GAMMA, R=R, q=q, qa=qa)
   vel = cshift(prim(2:4), d - 1)
   call primitive_to_conservative(gamma=GAMMA, r=prim(1), u=vel(1), v=vel(2), w=vel(3), p=prim(5), q=qp)
   call conservative_to_auxiliary(gamma=GAMMA, R=R, q=qp, qa=qap)
   perm = [1, 1+mod(d-1,3)+1, 1+mod(d,3)+1, 1+mod(d+1,3)+1, 5]
   call compute_flux(d=d, q=q, qa=qa, f=f)
   call compute_flux(d=1, q=qp, qa=qap, f=fp)
   call compute_eigenvectors(gamma=GAMMA, d=d, qa=qa, el=el, er=er)
   call compute_eigenvectors(gamma=GAMMA, d=1, qa=qap, el=elp, er=erp)
   e = max(e, maxval(abs(f(perm) - fp)) / maxval(abs(fp)))
   e = max(e, maxval(abs(er(perm,:) - erp)) / maxval(abs(erp)))
   e = max(e, maxval(abs(el(:,perm) - elp)) / maxval(abs(elp)))
   endsubroutine check_cyclic

   subroutine check_split(d, q, qa, e)
   !< Check that a uniform stencil, split and back-projected, returns the physical flux (S = 1..S_MAX, both variants).
   integer(I4P), intent(in)    :: d                                  !< Direction.
   real(R8P),    intent(in)    :: q(NV_EULER)                        !< Conservative variables.
   real(R8P),    intent(in)    :: qa(NV_AUX)                         !< Auxiliary variables.
   real(R8P),    intent(inout) :: e                                  !< Maximum error.
   real(R8P)                   :: qs(NV_EULER,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P)                   :: qas(NV_AUX,1-S_MAX:S_MAX)          !< Stencil auxiliary variables.
   real(R8P)                   :: fsplit(2,1-S_MAX:S_MAX-1,NV_EULER) !< Split fields.
   real(R8P)                   :: er(NV_EULER,NV_EULER)              !< Right eigenvectors.
   real(R8P)                   :: vr(2,NV_EULER)                     !< Reconstructed split fields.
   real(R8P)                   :: flux(NV_EULER)                     !< Face flux.
   real(R8P)                   :: f(NV_EULER)                        !< Physical flux.
   integer(I4P)                :: S, m, variant                      !< Counters.

   do m=1-S_MAX, S_MAX
      qs(:,m)  = q
      qas(:,m) = qa
   enddo
   call compute_flux(d=d, q=q, qa=qa, f=f)
   do variant=1, 2
      do S=1, S_MAX
         call compute_face_split_fluxes(gamma=GAMMA, d=d, S=S, is_characteristic=(variant == 1), qs=qs, qas=qas, &
                                        fsplit=fsplit, er=er)
         vr = fsplit(:,0,:)
         call compute_face_flux_back_projection(is_characteristic=(variant == 1), er=er, vr=vr, flux=flux)
         e = max(e, maxval(abs(flux - f)) / maxval(abs(f)))
      enddo
   enddo
   endsubroutine check_split
endprogram test_flume_euler_library
