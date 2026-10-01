!< Unit test RV-0 (MHD, device half) of the FLUME MHD Riemann solvers: device results against host results.
module test_flume_mhd_riemann_fnl_kernels
!< Unit test RV-0 (MHD, device half) of the FLUME MHD Riemann solvers: the evaluation routine run on host and device.
!<
!< Module-level (not a `contains`-internal procedure): device routines must be module procedures.

use :: adam_flume_mhd_library,         only : mhd_primitive_to_conservative
use :: adam_flume_mhd_riemann_library, only : mhd_eglm_riemann_hll, mhd_eglm_riemann_hlld, mhd_eglm_riemann_llf, &
                                              mhd_glm_riemann_hll, mhd_glm_riemann_hlld, mhd_glm_riemann_llf,    &
                                              mhd_riemann_hll, mhd_riemann_hlld, mhd_riemann_llf
use :: adam_flume_parameters,          only : IQ_PSI, IQ_RE, NV_MHD, NV_MHD_EGLM, NV_MHD_GLM
use :: penf,                           only : I4P, R8P

implicit none
private
public :: evaluate
public :: CH
public :: GAMMA
public :: IFB1
public :: IFB2
public :: IFB3
public :: NR

integer(I4P), parameter :: NR=3*NV_MHD+3*NV_MHD_GLM+3*NV_MHD_EGLM+3 !< Results per pair and direction: nine fluxes, three
                                                                   !< fallback flags.
integer(I4P), parameter :: IFB1=3*NV_MHD+1                  !< Result index of the HLLD fallback flag (no divergence control).
integer(I4P), parameter :: IFB2=3*NV_MHD+3*NV_MHD_GLM+2     !< Result index of the HLLD fallback flag (GLM).
integer(I4P), parameter :: IFB3=NR                          !< Result index of the HLLD fallback flag (EGLM).
real(R8P),    parameter :: GAMMA=5._R8P/3._R8P        !< Specific heats ratio.
real(R8P),    parameter :: CH=2._R8P                  !< GLM cleaning speed.

contains
   subroutine evaluate(d, prim, res)
   !< Evaluate the nine solvers on one pair of primitive states `prim(1:9,1:2)` (`psi` the ninth), direction `d`; the
   !< EGLM states are the GLM ones with `psi^2 / 2` added to the energy.
   integer(I4P), intent(in)  :: d               !< Direction.
   real(R8P),    intent(in)  :: prim(9,2)       !< Primitive states of the pair.
   real(R8P),    intent(out) :: res(NR)         !< Results.
   real(R8P)                 :: q(NV_MHD_GLM,2) !< Conservative states.
   real(R8P)                 :: qe(NV_MHD_EGLM,2) !< EGLM conservative states.
   real(R8P)                 :: f(NV_MHD)       !< Flux, no divergence control.
   real(R8P)                 :: fg(NV_MHD_GLM)  !< Flux, GLM.
   logical                   :: fb              !< Fallback flag.
   integer(I4P)              :: m, c, k         !< Counters.
   !$acc routine seq
   !$omp declare target

   do m=1, 2
      call mhd_primitive_to_conservative(gamma=GAMMA, r=prim(1,m), u=prim(2,m), v=prim(3,m), w=prim(4,m), p=prim(5,m), &
                                         bx=prim(6,m), by=prim(7,m), bz=prim(8,m), q=q(:,m))
      q(IQ_PSI,m) = prim(9,m)
   enddo
   c = 0
   call mhd_riemann_llf(gamma=GAMMA, d=d, qL=q(1:NV_MHD,1), qR=q(1:NV_MHD,2), f=f)
   do k=1, NV_MHD
      c = c + 1 ; res(c) = f(k)
   enddo
   call mhd_riemann_hll(gamma=GAMMA, d=d, qL=q(1:NV_MHD,1), qR=q(1:NV_MHD,2), f=f)
   do k=1, NV_MHD
      c = c + 1 ; res(c) = f(k)
   enddo
   call mhd_riemann_hlld(gamma=GAMMA, d=d, qL=q(1:NV_MHD,1), qR=q(1:NV_MHD,2), f=f, fallback=fb)
   do k=1, NV_MHD
      c = c + 1 ; res(c) = f(k)
   enddo
   c = c + 1 ; res(c) = merge(1._R8P, 0._R8P, fb)
   call mhd_glm_riemann_llf(ch=CH, gamma=GAMMA, d=d, qL=q(:,1), qR=q(:,2), f=fg)
   do k=1, NV_MHD_GLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   call mhd_glm_riemann_hll(ch=CH, gamma=GAMMA, d=d, qL=q(:,1), qR=q(:,2), f=fg)
   do k=1, NV_MHD_GLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   call mhd_glm_riemann_hlld(ch=CH, gamma=GAMMA, d=d, qL=q(:,1), qR=q(:,2), f=fg, fallback=fb)
   do k=1, NV_MHD_GLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   c = c + 1 ; res(c) = merge(1._R8P, 0._R8P, fb)
   do m=1, 2
      qe(:,m) = q(:,m)
      qe(IQ_RE,m) = q(IQ_RE,m) + 0.5_R8P * q(IQ_PSI,m)**2
   enddo
   call mhd_eglm_riemann_llf(ch=CH, gamma=GAMMA, d=d, qL=qe(:,1), qR=qe(:,2), f=fg)
   do k=1, NV_MHD_EGLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   call mhd_eglm_riemann_hll(ch=CH, gamma=GAMMA, d=d, qL=qe(:,1), qR=qe(:,2), f=fg)
   do k=1, NV_MHD_EGLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   call mhd_eglm_riemann_hlld(ch=CH, gamma=GAMMA, d=d, qL=qe(:,1), qR=qe(:,2), f=fg, fallback=fb)
   do k=1, NV_MHD_EGLM
      c = c + 1 ; res(c) = fg(k)
   enddo
   c = c + 1 ; res(c) = merge(1._R8P, 0._R8P, fb)
   endsubroutine evaluate
endmodule test_flume_mhd_riemann_fnl_kernels

program test_flume_mhd_riemann_fnl
!< Unit test RV-0 (MHD, device half) of the FLUME MHD Riemann solvers: device results against host results.
!<
!< **Why this test exists** (issue #47, RV-0). The MHD Riemann solvers are shared by the CPU loops and the FNL kernels
!< (`!$acc routine seq` + `!$omp declare target`); `test_flume_mhd_riemann` pins their correctness on the host. This
!< test pins that the device build of the SAME source returns the host results: for N random admissible pairs and every
!< direction, the LLF, HLL and HLLD fluxes without divergence control, with GLM and with EGLM, and the three HLLD fallback
!< flags,
!< compared relative to each quantity's magnitude (host and device may contract differently into FMAs, so not bitwise;
!< the fallback flags must agree exactly, reported as a count of mismatches).

use :: penf,                               only : I4P, R8P, str
use :: test_flume_mhd_riemann_fnl_kernels, only : evaluate, GAMMA, IFB1, IFB2, IFB3, NR

implicit none

integer(I4P), parameter :: N=4096_I4P       !< Random pairs number.
real(R8P),    parameter :: TOL=1.e-12_R8P   !< Relative tolerance.
real(R8P)               :: prim(9,2,N)      !< Pairs of primitive states.
real(R8P)               :: res_dev(NR,3,N)  !< Results, device.
real(R8P)               :: res_host(NR,3,N) !< Results, host.
real(R8P)               :: err              !< Maximum relative difference.
real(R8P)               :: x, a             !< Random number, sound speed.
integer(I4P)            :: seed(64)         !< Random generator seed.
integer(I4P)            :: ns_, n_, m, d, c !< Counters, seed size.
integer(I4P)            :: flags            !< Fallback flag mismatches.
logical                 :: test_passed      !< Aggregate pass flag.

call random_seed(size=ns_)
if (ns_ > size(seed)) error stop 'random seed larger than expected'
seed = [(20260930_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns_))
do n_=1, N
   do m=1, 2
      call random_number(x) ; prim(1,m,n_) = 10._R8P**(-1._R8P + 2._R8P * x)
      call random_number(x) ; prim(5,m,n_) = 10._R8P**(-1._R8P + 2._R8P * x)
      a = sqrt(GAMMA * prim(5,m,n_) / prim(1,m,n_))
      do c=2, 4
         call random_number(x) ; prim(c,m,n_) = a * (2._R8P * x - 1._R8P)
      enddo
      do c=6, 8
         call random_number(x) ; prim(c,m,n_) = sqrt(prim(1,m,n_)) * a * (2._R8P * x - 1._R8P)
      enddo
      call random_number(x) ; prim(9,m,n_) = 2._R8P * x - 1._R8P
   enddo
enddo

!$acc parallel loop gang vector collapse(2) copyin(prim) copyout(res_dev)
!$omp target teams distribute parallel do collapse(2) map(to:prim) map(from:res_dev)
do n_=1, N
do d=1, 3
   call evaluate(d=d, prim=prim(:,:,n_), res=res_dev(:,d,n_))
enddo
enddo
do n_=1, N
   do d=1, 3
      call evaluate(d=d, prim=prim(:,:,n_), res=res_host(:,d,n_))
   enddo
enddo

err = maxval(abs(res_dev - res_host) / max(1._R8P, abs(res_host)))
flags = count(res_dev(IFB1,:,:) /= res_host(IFB1,:,:)) + count(res_dev(IFB2,:,:) /= res_host(IFB2,:,:)) + &
        count(res_dev(IFB3,:,:) /= res_host(IFB3,:,:))
test_passed = (err <= TOL) .and. (flags == 0)
print '(A)', 'device vs host, nine MHD Riemann fluxes: max relative difference '//trim(str(err))// &
             ', fallback flag mismatches '//trim(str(flags))
if (test_passed) then
   print '(A)', 'TEST PASSED: flume mhd riemann, device vs host ('//trim(str(N))//' pairs x 3 directions)'
else
   print '(A)', 'TEST FAILED: flume mhd riemann, device vs host'
   error stop 1
endif
endprogram test_flume_mhd_riemann_fnl
