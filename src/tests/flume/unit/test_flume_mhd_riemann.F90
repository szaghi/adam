!< Unit test RV-0 (MHD) of the FLUME MHD Riemann solvers of the scheme `weno-riemann` (issue #47).
program test_flume_mhd_riemann
!< Unit test RV-0 (MHD) of the FLUME MHD Riemann solvers of the scheme `weno-riemann` (issue #47).
!<
!< **Why this test exists** (issue #47, section 8, RV-0). The HLLD solver is exact on isolated contacts, tangential and
!< rotational discontinuities, and a wrong star state (a sign, a missing `B_n` term) usually still gives a consistent,
!< dissipative flux that runs: only these exactness identities expose it. Bitwise cyclic invariance is what makes the
!< x, y and z runs of a problem bitwise equal (issue #41, MV-4).
!<
!< **What it pins**, on N deterministic random admissible states and each direction d, for LLF, HLL, HLLD without
!< divergence control and with GLM (errors relative to `max(1, |f|, s |q|)`, the round-off scale of the star fluxes,
!< `s` the largest `|u_n| + c_f`):
!< 1. consistency `RS(q, q) = f(q)` (all six solvers);
!< 2. HLLD exact on an isolated contact (`B_n /= 0`, only the density jumps): the flux of the upwind state;
!< 3. HLLD exact on an isolated tangential discontinuity (`B_n = 0`; density, tangential velocity and field jump, total
!<    pressure continuous);
!< 4. HLLD exact on an isolated rotational discontinuity (Alfven wave, speed `u_n +- |B_n|/sqrt(rho)`, the tangential
!<    field rotated, the velocity jump `-(S - u_n) [B_t] / B_n`); the construction is checked first by its
!<    Rankine-Hugoniot residual `|F_R - F_L - S (U_R - U_L)|`, printed and required below the same tolerance;
!< 5. cyclic invariance, BITWISE: the flux in direction d of a pair equals the direction-1 flux of the cyclically
!<    permuted pair, components permuted back (all six solvers);
!< 6. the first-order 1-D update `q_i - dt/dx (F(q_i, q_i+1) - F(q_i-1, q_i))` keeps density and pressure positive
!<    for LLF and HLL at `dt/dx = 1/(2 s)` (count of inadmissible updates; one `B_n` per triple, 1-D div B = 0);
!< 7. the same update with HLLD (count, reported; not part of the issue's criterion);
!< 8. HLLD fallbacks to HLL on random pairs (count, reported).

use :: adam_flume_mhd_library,         only : mhd_conservative_to_auxiliary, mhd_fast_speed, mhd_flux, mhd_glm_flux, &
                                              mhd_primitive_to_conservative
use :: adam_flume_mhd_riemann_library, only : mhd_glm_riemann_hll, mhd_glm_riemann_hlld, mhd_glm_riemann_llf,        &
                                              mhd_riemann_hll, mhd_riemann_hlld, mhd_riemann_llf
use :: adam_flume_parameters,          only : IA_U, IQ_BX, IQ_PSI, IQ_R, IQ_RE, IQ_RU, NV_AUX_MHD, NV_MHD, NV_MHD_GLM
use :: penf,                           only : I4P, R8P, str

implicit none

integer(I4P), parameter :: N=5000_I4P           !< Random states number.
integer(I4P), parameter :: NC=8_I4P             !< Checks number.
real(R8P),    parameter :: GAMMA=5._R8P/3._R8P  !< Specific heats ratio.
real(R8P),    parameter :: R=1._R8P             !< Gas constant.
real(R8P),    parameter :: CH=2._R8P            !< GLM cleaning speed.
real(R8P),    parameter :: TOL_EXACT=1.e-12_R8P !< Tolerance of the exact identities (relative).
real(R8P)               :: err(NC)              !< Maximum error (or count) of each check.
real(R8P)               :: rh                   !< Largest Rankine-Hugoniot residual of the rotational discontinuities.
integer(I4P)            :: seed(64)             !< Random generator seed.
integer(I4P)            :: n_, c, d, ns         !< Counters, seed size.
logical                 :: test_passed          !< Aggregate pass flag.
character(len=48)       :: check_name(NC)       !< Checks names.

check_name = ['RS(q, q) = f(q), 6 solvers                      ', &
              'HLLD exact on an isolated contact               ', &
              'HLLD exact on a tangential discontinuity        ', &
              'HLLD exact on a rotational discontinuity        ', &
              'cyclic invariance, bitwise, 6 solvers           ', &
              '1-D update positive, LLF + HLL (count)          ', &
              '1-D update positive, HLLD (count, info)         ', &
              'HLLD fallbacks on random pairs (count, info)    ']
err = 0._R8P
rh = 0._R8P
call random_seed(size=ns)
if (ns > size(seed)) error stop 'random seed larger than expected'
seed = [(20260930_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns))
do n_=1, N
   do d=1, 3
      call check_consistency(d=d, e=err(1))
      call check_contact(d=d, e=err(2))
      call check_tangential(d=d, e=err(3))
      call check_rotational(d=d, e=err(4), rh=rh)
      call check_cyclic(d=d, e=err(5))
      call check_positivity(d=d, e=err(6:7))
      call check_fallback(d=d, e=err(8))
   enddo
enddo

test_passed = rh <= TOL_EXACT
print '(A)', 'rotational discontinuity construction: max Rankine-Hugoniot residual '//trim(str(rh))
do c=1, NC
   if ((c <= 4 .and. err(c) > TOL_EXACT) .or. ((c == 5 .or. c == 6) .and. err(c) > 0._R8P)) then
      print '(A)', 'FAIL: '//check_name(c)//' max error '//trim(str(err(c)))
      test_passed = .false.
   else
      print '(A)', 'PASS: '//check_name(c)//' max error '//trim(str(err(c)))
   endif
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: flume mhd riemann ('//trim(str(N))//' states x 3 directions)'
else
   print '(A)', 'TEST FAILED: flume mhd riemann'
   error stop 1
endif

contains
   function uniform() result(x)
   !< Return a pseudo-random number in [0, 1) from the intrinsic generator (seeded once, deterministic).
   real(R8P) :: x !< Random number.

   call random_number(x)
   endfunction uniform

   subroutine random_prim(prim)
   !< Return a random admissible primitive state `(r, u, v, w, p, bx, by, bz)`: density and pressure log-uniform in
   !< [1e-3, 1e1], velocity components within +-2 a, field components of magnitude 0.1 to 10 times sqrt(rho) a.
   real(R8P), intent(out) :: prim(8) !< Primitive state.
   real(R8P)              :: a       !< Sound speed.
   integer(I4P)           :: c       !< Counter.

   prim(1) = 10._R8P**(-3._R8P + 4._R8P * uniform())
   prim(5) = 10._R8P**(-3._R8P + 4._R8P * uniform())
   a = sqrt(GAMMA * prim(5) / prim(1))
   do c=2, 4
      prim(c) = 2._R8P * a * (2._R8P * uniform() - 1._R8P)
   enddo
   do c=6, 8
      prim(c) = 10._R8P**(2._R8P * uniform() - 1._R8P) * sqrt(prim(1)) * a * (2._R8P * uniform() - 1._R8P)
   enddo
   endsubroutine random_prim

   subroutine to_conservative(prim, psi, q)
   !< Return the GLM-sized conservative state of a primitive state and a GLM scalar.
   real(R8P), intent(in)  :: prim(8)       !< Primitive state.
   real(R8P), intent(in)  :: psi           !< GLM scalar.
   real(R8P), intent(out) :: q(NV_MHD_GLM) !< Conservative variables.

   call mhd_primitive_to_conservative(gamma=GAMMA, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), bx=prim(6), &
                                      by=prim(7), bz=prim(8), q=q)
   q(IQ_PSI) = psi
   endsubroutine to_conservative

   subroutine solver(s, d, qL, qR, f, fallback)
   !< Return the flux of solver `s`: 1-3 LLF, HLL, HLLD without divergence control, 4-6 the same with GLM (1-3 set the
   !< first `NV_MHD` components of `f`, the last one zero).
   integer(I4P), intent(in)  :: s              !< Solver.
   integer(I4P), intent(in)  :: d              !< Direction.
   real(R8P),    intent(in)  :: qL(NV_MHD_GLM) !< Left state.
   real(R8P),    intent(in)  :: qR(NV_MHD_GLM) !< Right state.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)  !< Flux.
   logical,      intent(out) :: fallback       !< HLLD fallback flag.

   f = 0._R8P
   fallback = .false.
   select case(s)
   case(1)
      call mhd_riemann_llf(gamma=GAMMA, d=d, qL=qL(1:NV_MHD), qR=qR(1:NV_MHD), f=f(1:NV_MHD))
   case(2)
      call mhd_riemann_hll(gamma=GAMMA, d=d, qL=qL(1:NV_MHD), qR=qR(1:NV_MHD), f=f(1:NV_MHD))
   case(3)
      call mhd_riemann_hlld(gamma=GAMMA, d=d, qL=qL(1:NV_MHD), qR=qR(1:NV_MHD), f=f(1:NV_MHD), fallback=fallback)
   case(4)
      call mhd_glm_riemann_llf(ch=CH, gamma=GAMMA, d=d, qL=qL, qR=qR, f=f)
   case(5)
      call mhd_glm_riemann_hll(ch=CH, gamma=GAMMA, d=d, qL=qL, qR=qR, f=f)
   case(6)
      call mhd_glm_riemann_hlld(ch=CH, gamma=GAMMA, d=d, qL=qL, qR=qR, f=f, fallback=fallback)
   endselect
   endsubroutine solver

   subroutine physical_flux(glm, d, q, f, s)
   !< Return the physical flux of a state (without divergence control or with GLM) and its `|u_n| + c_f`.
   logical,      intent(in)  :: glm            !< GLM flux.
   integer(I4P), intent(in)  :: d              !< Direction.
   real(R8P),    intent(in)  :: q(NV_MHD_GLM)  !< Conservative variables.
   real(R8P),    intent(out) :: f(NV_MHD_GLM)  !< Physical flux.
   real(R8P),    intent(out) :: s              !< |u_n| + c_f.
   real(R8P)                 :: qa(NV_AUX_MHD) !< Auxiliary variables.

   call mhd_conservative_to_auxiliary(gamma=GAMMA, R=R, q=q(1:NV_MHD), qa=qa)
   f = 0._R8P
   if (glm) then
      call mhd_glm_flux(ch=CH, d=d, q=q, qa=qa, f=f)
   else
      call mhd_flux(d=d, q=q(1:NV_MHD), qa=qa, f=f(1:NV_MHD))
   endif
   s = abs(qa(IA_U+d-1)) + mhd_fast_speed(d=d, qa=qa)
   endsubroutine physical_flux

   function scaled_error(fr, f, s, qL, qR, nv) result(e)
   !< Return `max|fr - f| / max(1, max|f|, s max(|qL|, |qR|))` over the first `nv` components.
   real(R8P),    intent(in) :: fr(NV_MHD_GLM) !< Solver flux.
   real(R8P),    intent(in) :: f(NV_MHD_GLM)  !< Exact flux.
   real(R8P),    intent(in) :: s              !< Wave speed scale.
   real(R8P),    intent(in) :: qL(NV_MHD_GLM) !< Left state.
   real(R8P),    intent(in) :: qR(NV_MHD_GLM) !< Right state.
   integer(I4P), intent(in) :: nv             !< Components compared.
   real(R8P)                :: e              !< Error.

   e = maxval(abs(fr(1:nv) - f(1:nv))) / max(1._R8P, maxval(abs(f(1:nv))), &
                                             s * max(maxval(abs(qL(1:nv))), maxval(abs(qR(1:nv)))))
   endfunction scaled_error

   subroutine frame_to_prim(d, fp, prim)
   !< Return the global primitive state of a frame primitive state `(rho, u_n, u_t1, u_t2, p, B_n, B_t1, B_t2)` of
   !< direction `d` (tangents cyclic).
   integer(I4P), intent(in)  :: d       !< Direction.
   real(R8P),    intent(in)  :: fp(8)   !< Frame primitive state.
   real(R8P),    intent(out) :: prim(8) !< Global primitive state.
   integer(I4P)              :: d1, d2  !< Tangential directions.

   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   prim(1)    = fp(1)
   prim(5)    = fp(5)
   prim(1+d)  = fp(2)
   prim(1+d1) = fp(3)
   prim(1+d2) = fp(4)
   prim(5+d)  = fp(6)
   prim(5+d1) = fp(7)
   prim(5+d2) = fp(8)
   endsubroutine frame_to_prim

   subroutine check_consistency(d, e)
   !< Update `e` with the consistency error `RS(q, q) - f(q)` of the six solvers on a random state.
   integer(I4P), intent(in)    :: d              !< Direction.
   real(R8P),    intent(inout) :: e              !< Maximum error.
   real(R8P)                   :: prim(8)        !< Primitive state.
   real(R8P)                   :: q(NV_MHD_GLM)  !< Conservative variables.
   real(R8P)                   :: f(NV_MHD_GLM)  !< Physical flux.
   real(R8P)                   :: fr(NV_MHD_GLM) !< Solver flux.
   real(R8P)                   :: s              !< Wave speed scale.
   logical                     :: fb             !< Fallback flag.
   integer(I4P)                :: k              !< Counter.

   call random_prim(prim=prim)
   call to_conservative(prim=prim, psi=2._R8P * uniform() - 1._R8P, q=q)
   do k=1, 6
      call physical_flux(glm=(k > 3), d=d, q=q, f=f, s=s)
      call solver(s=k, d=d, qL=q, qR=q, f=fr, fallback=fb)
      e = max(e, scaled_error(fr=fr, f=f, s=max(s, CH), qL=q, qR=q, nv=merge(NV_MHD_GLM, NV_MHD, k > 3)))
   enddo
   endsubroutine check_consistency

   subroutine check_contact(d, e)
   !< Update `e` with the HLLD error on an isolated contact: only the density jumps, `B_n /= 0`.
   integer(I4P), intent(in)    :: d              !< Direction.
   real(R8P),    intent(inout) :: e              !< Maximum error.
   real(R8P)                   :: prim(8)        !< Left primitive state.
   real(R8P)                   :: primR(8)       !< Right primitive state.
   real(R8P)                   :: qL(NV_MHD_GLM) !< Left state.
   real(R8P)                   :: qR(NV_MHD_GLM) !< Right state.
   real(R8P)                   :: f(NV_MHD_GLM)  !< Exact flux.
   real(R8P)                   :: fr(NV_MHD_GLM) !< HLLD flux.
   real(R8P)                   :: s, sR          !< Wave speed scales.
   logical                     :: fb             !< Fallback flag.
   integer(I4P)                :: k              !< Counter.

   call random_prim(prim=prim)
   primR = prim
   primR(1) = 10._R8P**(-3._R8P + 4._R8P * uniform())
   call to_conservative(prim=prim, psi=0._R8P, q=qL)
   call to_conservative(prim=primR, psi=0._R8P, q=qR)
   do k=3, 6, 3
      call physical_flux(glm=(k > 3), d=d, q=qR, f=f, s=sR)
      s = sR
      if (prim(1+d) >= 0._R8P) call physical_flux(glm=(k > 3), d=d, q=qL, f=f, s=s)
      call solver(s=k, d=d, qL=qL, qR=qR, f=fr, fallback=fb)
      e = max(e, scaled_error(fr=fr, f=f, s=max(s, sR), qL=qL, qR=qR, nv=NV_MHD))
   enddo
   endsubroutine check_contact

   subroutine check_tangential(d, e)
   !< Update `e` with the HLLD error on an isolated tangential discontinuity: `B_n = 0`, the density, the tangential
   !< velocity and field jump, the normal velocity and the total pressure are continuous.
   integer(I4P), intent(in)    :: d              !< Direction.
   real(R8P),    intent(inout) :: e              !< Maximum error.
   real(R8P)                   :: fp(8), fpR(8)  !< Frame primitive states.
   real(R8P)                   :: prim(8)        !< Global primitive state.
   real(R8P)                   :: qL(NV_MHD_GLM) !< Left state.
   real(R8P)                   :: qR(NV_MHD_GLM) !< Right state.
   real(R8P)                   :: f(NV_MHD_GLM)  !< Exact flux.
   real(R8P)                   :: fr(NV_MHD_GLM) !< HLLD flux.
   real(R8P)                   :: s, sR, pt      !< Wave speed scales, total pressure.
   logical                     :: fb             !< Fallback flag.
   integer(I4P)                :: k              !< Counter.

   call random_prim(prim=fp)
   fp(6) = 0._R8P
   pt = fp(5) + 0.5_R8P * (fp(7)**2 + fp(8)**2)
   call random_prim(prim=fpR)
   fpR(2) = fp(2)
   fpR(6) = 0._R8P
   ! the right magnetic pressure is at most 0.9 pt, the gas pressure takes the rest
   fpR(7:8) = fpR(7:8) * min(1._R8P, sqrt(1.8_R8P * pt / max(fpR(7)**2 + fpR(8)**2, tiny(1._R8P))))
   fpR(5) = pt - 0.5_R8P * (fpR(7)**2 + fpR(8)**2)
   call frame_to_prim(d=d, fp=fp, prim=prim)
   call to_conservative(prim=prim, psi=0._R8P, q=qL)
   call frame_to_prim(d=d, fp=fpR, prim=prim)
   call to_conservative(prim=prim, psi=0._R8P, q=qR)
   do k=3, 6, 3
      call physical_flux(glm=(k > 3), d=d, q=qR, f=f, s=sR)
      s = sR
      if (fp(2) >= 0._R8P) call physical_flux(glm=(k > 3), d=d, q=qL, f=f, s=s)
      call solver(s=k, d=d, qL=qL, qR=qR, f=fr, fallback=fb)
      e = max(e, scaled_error(fr=fr, f=f, s=max(s, sR), qL=qL, qR=qR, nv=NV_MHD))
   enddo
   endsubroutine check_tangential

   subroutine check_rotational(d, e, rh)
   !< Update `e` with the HLLD error on an isolated rotational discontinuity and `rh` with the Rankine-Hugoniot
   !< residual of its construction.
   integer(I4P), intent(in)    :: d               !< Direction.
   real(R8P),    intent(inout) :: e               !< Maximum error.
   real(R8P),    intent(inout) :: rh              !< Maximum Rankine-Hugoniot residual.
   real(R8P)                   :: fp(8), fpR(8)   !< Frame primitive states.
   real(R8P)                   :: prim(8)         !< Global primitive state.
   real(R8P)                   :: qL(NV_MHD_GLM)  !< Left state.
   real(R8P)                   :: qR(NV_MHD_GLM)  !< Right state.
   real(R8P)                   :: fL(NV_MHD_GLM)  !< Left physical flux.
   real(R8P)                   :: fR(NV_MHD_GLM)  !< Right physical flux.
   real(R8P)                   :: fr_(NV_MHD_GLM) !< HLLD flux.
   real(R8P)                   :: s, sL_, sR_     !< Wave speed, speed scales.
   real(R8P)                   :: th, ct, st      !< Rotation angle, its cosine and sine.
   logical                     :: fb              !< Fallback flag.
   integer(I4P)                :: k               !< Counter.

   call random_prim(prim=fp)
   if (abs(fp(6)) < 1.e-2_R8P * sqrt(GAMMA * fp(5))) fp(6) = sign(sqrt(GAMMA * fp(5)), fp(6))
   th = acos(-1._R8P) * (2._R8P * uniform() - 1._R8P)
   ct = cos(th)
   st = sin(th)
   fpR = fp
   fpR(7) = ct * fp(7) - st * fp(8)
   fpR(8) = st * fp(7) + ct * fp(8)
   ! wave speed S = u_n + sigma |B_n| / sqrt(rho): velocity jump -(S - u_n) [B_t] / B_n
   s = fp(2) + sign(1._R8P, uniform() - 0.5_R8P) * abs(fp(6)) / sqrt(fp(1))
   fpR(3) = fp(3) - (s - fp(2)) * (fpR(7) - fp(7)) / fp(6)
   fpR(4) = fp(4) - (s - fp(2)) * (fpR(8) - fp(8)) / fp(6)
   call frame_to_prim(d=d, fp=fp, prim=prim)
   call to_conservative(prim=prim, psi=0._R8P, q=qL)
   call frame_to_prim(d=d, fp=fpR, prim=prim)
   call to_conservative(prim=prim, psi=0._R8P, q=qR)
   call physical_flux(glm=.false., d=d, q=qL, f=fL, s=sL_)
   call physical_flux(glm=.false., d=d, q=qR, f=fR, s=sR_)
   rh = max(rh, maxval(abs(fR(1:NV_MHD) - fL(1:NV_MHD) - s * (qR(1:NV_MHD) - qL(1:NV_MHD)))) / &
                max(1._R8P, maxval(abs(fL(1:NV_MHD))), max(sL_, sR_) * max(maxval(abs(qL(1:NV_MHD))), &
                                                                          maxval(abs(qR(1:NV_MHD))))))
   do k=3, 6, 3
      call physical_flux(glm=(k > 3), d=d, q=qL, f=fL, s=sL_)
      call physical_flux(glm=(k > 3), d=d, q=qR, f=fR, s=sR_)
      call solver(s=k, d=d, qL=qL, qR=qR, f=fr_, fallback=fb)
      if (s >= 0._R8P) then
         e = max(e, scaled_error(fr=fr_, f=fL, s=max(sL_, sR_), qL=qL, qR=qR, nv=NV_MHD))
      else
         e = max(e, scaled_error(fr=fr_, f=fR, s=max(sL_, sR_), qL=qL, qR=qR, nv=NV_MHD))
      endif
   enddo
   endsubroutine check_rotational

   subroutine check_cyclic(d, e)
   !< Update `e` with the largest difference between the direction-`d` flux of a random pair and the direction-1 flux of
   !< the cyclically permuted pair, components permuted back (must be exactly zero).
   integer(I4P), intent(in)    :: d               !< Direction.
   real(R8P),    intent(inout) :: e               !< Maximum difference.
   real(R8P)                   :: pL(8), pR(8)    !< Primitive states.
   real(R8P)                   :: fpL(8), fpR(8)  !< Frame primitive states of direction d.
   real(R8P)                   :: gL(8), gR(8)    !< Permuted primitive states (direction 1).
   real(R8P)                   :: qL(NV_MHD_GLM)  !< Left state.
   real(R8P)                   :: qR(NV_MHD_GLM)  !< Right state.
   real(R8P)                   :: q1L(NV_MHD_GLM) !< Permuted left state.
   real(R8P)                   :: q1R(NV_MHD_GLM) !< Permuted right state.
   real(R8P)                   :: f(NV_MHD_GLM)   !< Flux in direction d.
   real(R8P)                   :: f1(NV_MHD_GLM)  !< Flux of the permuted pair in direction 1.
   real(R8P)                   :: psiL, psiR      !< GLM scalars.
   logical                     :: fb              !< Fallback flag.
   integer(I4P)                :: k, d1, d2       !< Counter, tangential directions.

   call random_prim(prim=pL)
   call random_prim(prim=pR)
   psiL = 2._R8P * uniform() - 1._R8P
   psiR = 2._R8P * uniform() - 1._R8P
   d1 = mod(d, 3) + 1
   d2 = mod(d + 1, 3) + 1
   fpL = [pL(1), pL(1+d), pL(1+d1), pL(1+d2), pL(5), pL(5+d), pL(5+d1), pL(5+d2)]
   fpR = [pR(1), pR(1+d), pR(1+d1), pR(1+d2), pR(5), pR(5+d), pR(5+d1), pR(5+d2)]
   call frame_to_prim(d=1, fp=fpL, prim=gL)
   call frame_to_prim(d=1, fp=fpR, prim=gR)
   call to_conservative(prim=pL, psi=psiL, q=qL)
   call to_conservative(prim=pR, psi=psiR, q=qR)
   call to_conservative(prim=gL, psi=psiL, q=q1L)
   call to_conservative(prim=gR, psi=psiR, q=q1R)
   do k=1, 6
      call solver(s=k, d=d, qL=qL, qR=qR, f=f, fallback=fb)
      call solver(s=k, d=1, qL=q1L, qR=q1R, f=f1, fallback=fb)
      e = max(e, abs(f(IQ_R) - f1(IQ_R)), abs(f(IQ_RE) - f1(IQ_RE)), abs(f(IQ_PSI) - f1(IQ_PSI)),                  &
              abs(f(IQ_RU+d-1) - f1(IQ_RU)), abs(f(IQ_RU+d1-1) - f1(IQ_RU+1)), abs(f(IQ_RU+d2-1) - f1(IQ_RU+2)), &
              abs(f(IQ_BX+d-1) - f1(IQ_BX)), abs(f(IQ_BX+d1-1) - f1(IQ_BX+1)), abs(f(IQ_BX+d2-1) - f1(IQ_BX+2)))
   enddo
   endsubroutine check_cyclic

   subroutine check_positivity(d, e)
   !< Update the counts `e(1)` (LLF, HLL) and `e(2)` (HLLD) of inadmissible first-order 1-D updates on a random triple
   !< of states sharing one normal field.
   integer(I4P), intent(in)    :: d                  !< Direction.
   real(R8P),    intent(inout) :: e(2)               !< Inadmissible updates counts.
   real(R8P)                   :: prim(8)            !< Primitive state.
   real(R8P)                   :: q(NV_MHD_GLM,-1:1) !< Conservative variables of the triple.
   real(R8P)                   :: f(NV_MHD_GLM)      !< Physical flux (unused).
   real(R8P)                   :: fm(NV_MHD_GLM)     !< Flux of the face i-1/2.
   real(R8P)                   :: fp(NV_MHD_GLM)     !< Flux of the face i+1/2.
   real(R8P)                   :: qn(NV_MHD)         !< Updated state.
   real(R8P)                   :: s, smax, bn, pn    !< Wave speed, largest one, normal field, updated pressure.
   logical                     :: fb                 !< Fallback flag.
   integer(I4P)                :: m, k               !< Counters.

   smax = 0._R8P
   bn = 0._R8P
   do m=-1, 1
      call random_prim(prim=prim)
      if (m == -1) bn = prim(5+d)
      prim(5+d) = bn
      call to_conservative(prim=prim, psi=0._R8P, q=q(:,m))
      call physical_flux(glm=.false., d=d, q=q(:,m), f=f, s=s)
      smax = max(smax, s)
   enddo
   do k=1, 3
      call solver(s=k, d=d, qL=q(:,-1), qR=q(:,0), f=fm, fallback=fb)
      call solver(s=k, d=d, qL=q(:,0), qR=q(:,1), f=fp, fallback=fb)
      qn = q(1:NV_MHD,0) - 0.5_R8P / smax * (fp(1:NV_MHD) - fm(1:NV_MHD))
      pn = -1._R8P
      if (qn(IQ_R) > 0._R8P) pn = qn(IQ_RE) - 0.5_R8P * sum(qn(IQ_RU:IQ_RU+2)**2) / qn(IQ_R) - &
                                  0.5_R8P * sum(qn(IQ_BX:IQ_BX+2)**2)
      if (.not.(pn > 0._R8P)) then
         if (k < 3) then
            e(1) = e(1) + 1._R8P
         else
            e(2) = e(2) + 1._R8P
         endif
      endif
   enddo
   endsubroutine check_positivity

   subroutine check_fallback(d, e)
   !< Update the count `e` of HLLD fallbacks to HLL on a random pair (both variants).
   integer(I4P), intent(in)    :: d              !< Direction.
   real(R8P),    intent(inout) :: e              !< Fallbacks count.
   real(R8P)                   :: pL(8), pR(8)   !< Primitive states.
   real(R8P)                   :: qL(NV_MHD_GLM) !< Left state.
   real(R8P)                   :: qR(NV_MHD_GLM) !< Right state.
   real(R8P)                   :: f(NV_MHD_GLM)  !< Flux.
   logical                     :: fb             !< Fallback flag.
   integer(I4P)                :: k              !< Counter.

   call random_prim(prim=pL)
   call random_prim(prim=pR)
   call to_conservative(prim=pL, psi=0._R8P, q=qL)
   call to_conservative(prim=pR, psi=0._R8P, q=qR)
   do k=3, 6, 3
      call solver(s=k, d=d, qL=qL, qR=qR, f=f, fallback=fb)
      if (fb) e = e + 1._R8P
   enddo
   endsubroutine check_fallback
endprogram test_flume_mhd_riemann
