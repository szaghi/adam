!< Unit test V0 (device half) of the FLUME Euler library: device results against host results.
module test_flume_euler_library_fnl_kernels
!< Unit test V0 (device half) of the FLUME Euler library: the evaluation routine run on both host and device.
!<
!< Module-level (not a `contains`-internal procedure): device routines must be module procedures.

use :: adam_flume_euler_library, only : compute_eigenvalues, compute_eigenvectors, compute_face_split_fluxes, compute_flux, &
                                        conservative_to_auxiliary, primitive_to_conservative
use :: adam_flume_parameters,    only : NV_AUX, NV_EULER, S_MAX
use :: penf,                     only : I4P, R8P

implicit none
private
public :: evaluate
public :: GAMMA
public :: NR
public :: NS
public :: R
public :: S

integer(I4P), parameter :: S=3_I4P                           !< WENO stencil half-width of the split test.
integer(I4P), parameter :: NR=NV_AUX+NV_EULER*(2+2*NV_EULER) !< Pointwise results per state and direction.
integer(I4P), parameter :: NS=2*(2*S_MAX)*NV_EULER           !< Split results per state and direction.
real(R8P),    parameter :: GAMMA=1.4_R8P                     !< Specific heats ratio.
real(R8P),    parameter :: R=287.05_R8P                      !< Gas constant.

contains
   subroutine evaluate(d, prim, res, split)
   !< Evaluate the library on one stencil of primitive states, direction `d`: pointwise results of cell 0 and the
   !< split fields of the whole stencil.
   integer(I4P), intent(in)  :: d                                  !< Direction.
   real(R8P),    intent(in)  :: prim(5,1-S_MAX:S_MAX)              !< Stencil of primitive states.
   real(R8P),    intent(out) :: res(NR)                            !< Pointwise results.
   real(R8P),    intent(out) :: split(NS)                          !< Split results.
   real(R8P)                 :: qs(NV_EULER,1-S_MAX:S_MAX)         !< Stencil conservative variables.
   real(R8P)                 :: qas(NV_AUX,1-S_MAX:S_MAX)          !< Stencil auxiliary variables.
   real(R8P)                 :: el(NV_EULER,NV_EULER)              !< Left eigenvectors.
   real(R8P)                 :: er(NV_EULER,NV_EULER)              !< Right eigenvectors.
   real(R8P)                 :: lambda(NV_EULER)                   !< Eigenvalues.
   real(R8P)                 :: f(NV_EULER)                        !< Physical flux.
   real(R8P)                 :: fsplit(2,1-S_MAX:S_MAX-1,NV_EULER) !< Split fields.
   integer(I4P)              :: m, i, j, k, c                      !< Counters.
   !$acc routine seq
   !$omp declare target

   do m=1-S_MAX, S_MAX
      call primitive_to_conservative(gamma=GAMMA, r=prim(1,m), u=prim(2,m), v=prim(3,m), w=prim(4,m), p=prim(5,m), &
                                     q=qs(:,m))
      call conservative_to_auxiliary(gamma=GAMMA, R=R, q=qs(:,m), qa=qas(:,m))
   enddo
   call compute_eigenvalues(d=d, qa=qas(:,0), lambda=lambda)
   call compute_flux(d=d, q=qs(:,0), qa=qas(:,0), f=f)
   call compute_eigenvectors(gamma=GAMMA, d=d, qa=qas(:,0), el=el, er=er)
   c = 0
   do i=1, NV_AUX
      c = c + 1 ; res(c) = qas(i,0)
   enddo
   do i=1, NV_EULER
      c = c + 1 ; res(c) = lambda(i)
      c = c + 1 ; res(c) = f(i)
   enddo
   do j=1, NV_EULER
      do i=1, NV_EULER
         c = c + 1 ; res(c) = el(i,j)
         c = c + 1 ; res(c) = er(i,j)
      enddo
   enddo
   call compute_face_split_fluxes(gamma=GAMMA, d=d, S=S, is_characteristic=.true., qs=qs, qas=qas, fsplit=fsplit, er=er)
   split = 0._R8P
   c = 0
   do k=1, NV_EULER
      do m=1-S_MAX, S_MAX-1
         do i=1, 2
            c = c + 1
            if (m >= 1-S .and. m <= S-1) split(c) = fsplit(i,m,k)
         enddo
      enddo
   enddo
   endsubroutine evaluate
endmodule test_flume_euler_library_fnl_kernels

program test_flume_euler_library_fnl
!< Unit test V0 (device half) of the FLUME Euler library: device results against host results.
!<
!< **Why this test exists** (issue #35, section 11, V0). The Euler library is shared by the CPU loops and the FNL
!< kernels (`!$acc routine seq` + `!$omp declare target`); its correctness is pinned on the host by
!< `test_flume_euler_library`. This test pins that the device build of the SAME source returns the host results:
!< for N random admissible stencils and every direction, auxiliary variables, eigenvalues, flux and eigenvectors of the
!< central cell, and the split fields of the non-uniform stencil (S = 3, exercising the Roe average and the per-wave
!< speeds), are computed in an offloaded loop and compared with the host evaluation, relative to each quantity's
!< magnitude. Host and device may contract differently into FMAs, so the comparison is not required to be bitwise.

use :: adam_flume_parameters,                only : S_MAX
use :: penf,                                 only : I4P, R8P, str
use :: test_flume_euler_library_fnl_kernels, only : evaluate, GAMMA, NR, NS

implicit none

integer(I4P), parameter :: N=4096_I4P              !< Random stencils number.
real(R8P),    parameter :: TOL=1.e-13_R8P          !< Relative tolerance.
real(R8P)               :: prim(5,1-S_MAX:S_MAX,N) !< Stencils of primitive states.
real(R8P)               :: res_dev(NR,3,N)         !< Pointwise results, device.
real(R8P)               :: res_host(NR,3,N)        !< Pointwise results, host.
real(R8P)               :: split_dev(NS,3,N)       !< Split results, device.
real(R8P)               :: split_host(NS,3,N)      !< Split results, host.
real(R8P)               :: err_res                 !< Maximum pointwise relative difference.
real(R8P)               :: err_split               !< Maximum split relative difference.
real(R8P)               :: x                       !< Random number.
integer(I4P)            :: seed(64)                !< Random generator seed.
integer(I4P)            :: ns_, n_, m, d, c        !< Counters, seed size.
logical                 :: test_passed             !< Aggregate pass flag.

call random_seed(size=ns_)
if (ns_ > size(seed)) error stop 'random seed larger than expected'
seed = [(20260923_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns_))
do n_=1, N
   do m=1-S_MAX, S_MAX
      call random_number(x) ; prim(1,m,n_) = 0.1_R8P + 9.9_R8P * x
      call random_number(x) ; prim(5,m,n_) = 0.1_R8P + 9.9_R8P * x
      do c=2, 4
         call random_number(x) ; prim(c,m,n_) = 3._R8P * sqrt(GAMMA * prim(5,m,n_) / prim(1,m,n_)) * (2._R8P * x - 1._R8P)
      enddo
   enddo
enddo

!$acc parallel loop gang vector collapse(2) copyin(prim) copyout(res_dev, split_dev)
!$omp target teams distribute parallel do collapse(2) map(to:prim) map(from:res_dev, split_dev)
do n_=1, N
do d=1, 3
   call evaluate(d=d, prim=prim(:,:,n_), res=res_dev(:,d,n_), split=split_dev(:,d,n_))
enddo
enddo
do n_=1, N
   do d=1, 3
      call evaluate(d=d, prim=prim(:,:,n_), res=res_host(:,d,n_), split=split_host(:,d,n_))
   enddo
enddo

err_res   = maxval(abs(res_dev   - res_host  ) / max(1._R8P, abs(res_host  )))
err_split = maxval(abs(split_dev - split_host) / max(1._R8P, abs(split_host)))
test_passed = (err_res <= TOL) .and. (err_split <= TOL)
print '(A)', 'device vs host, pointwise (aux, eigenvalues, flux, eigenvectors): max relative difference '// &
             trim(str(err_res))
print '(A)', 'device vs host, split fields of random stencils (S=3):            max relative difference '// &
             trim(str(err_split))
if (test_passed) then
   print '(A)', 'TEST PASSED: flume euler library, device vs host ('//trim(str(N))//' stencils x 3 directions)'
else
   print '(A)', 'TEST FAILED: flume euler library, device vs host'
   error stop 1
endif
endprogram test_flume_euler_library_fnl
