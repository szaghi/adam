!< Unit test MV-0 (device half) of the FLUME MHD library: device results against host results.
module test_flume_mhd_library_fnl_kernels
!< Unit test MV-0 (device half) of the FLUME MHD library: the evaluation routine run on both host and device.
!<
!< Module-level (not a `contains`-internal procedure): device routines must be module procedures.

use :: adam_flume_mhd_library, only : mhd_conservative_to_auxiliary, mhd_eigenvalues, mhd_eigenvectors,               &
                                      mhd_face_split_fluxes, mhd_flux, mhd_glm_eigenvalues, mhd_glm_eigenvectors,     &
                                      mhd_glm_face_split_fluxes, mhd_glm_flux, mhd_primitive_to_conservative
use :: adam_flume_parameters,  only : IQ_PSI, NV_AUX_MHD, NV_MHD, NV_MHD_GLM, S_MAX
use :: penf,                   only : I4P, R8P

implicit none
private
public :: evaluate
public :: GAMMA
public :: NR
public :: NS

integer(I4P), parameter :: S=3_I4P                                        !< WENO stencil half-width of the split test.
integer(I4P), parameter :: NR=NV_AUX_MHD+2*(NV_MHD+NV_MHD_GLM)+ &
                              2*(NV_MHD**2+NV_MHD_GLM**2)                  !< Pointwise results per state and direction.
integer(I4P), parameter :: NS=2*(2*S_MAX)*(NV_MHD+NV_MHD_GLM)             !< Split results per state and direction.
real(R8P),    parameter :: GAMMA=5._R8P/3._R8P                            !< Specific heats ratio.
real(R8P),    parameter :: R=1._R8P                                       !< Gas constant.

contains
   subroutine evaluate(d, ch, prim, res, split)
   !< Evaluate the library on one stencil of primitive states (both variants), direction `d`: pointwise results of
   !< cell 0 and the characteristic split fields of the whole stencil.
   integer(I4P), intent(in)  :: d                                      !< Direction.
   real(R8P),    intent(in)  :: ch                                     !< GLM cleaning speed.
   real(R8P),    intent(in)  :: prim(9,1-S_MAX:S_MAX)                  !< Stencil of primitive states and psi.
   real(R8P),    intent(out) :: res(NR)                                !< Pointwise results.
   real(R8P),    intent(out) :: split(NS)                              !< Split results.
   real(R8P)                 :: qs8(NV_MHD,1-S_MAX:S_MAX)              !< Stencil, no cleaning.
   real(R8P)                 :: qs9(NV_MHD_GLM,1-S_MAX:S_MAX)          !< Stencil, GLM.
   real(R8P)                 :: qas(NV_AUX_MHD,1-S_MAX:S_MAX)          !< Stencil auxiliary variables.
   real(R8P)                 :: el8(NV_MHD,NV_MHD), er8(NV_MHD,NV_MHD) !< Eigenvectors, no cleaning.
   real(R8P)                 :: el9(NV_MHD_GLM,NV_MHD_GLM)             !< Left eigenvectors, GLM.
   real(R8P)                 :: er9(NV_MHD_GLM,NV_MHD_GLM)             !< Right eigenvectors, GLM.
   real(R8P)                 :: lam8(NV_MHD), lam9(NV_MHD_GLM)         !< Eigenvalues.
   real(R8P)                 :: f8(NV_MHD), f9(NV_MHD_GLM)             !< Fluxes.
   real(R8P)                 :: fs8(2,1-S_MAX:S_MAX-1,NV_MHD)          !< Split fields, no cleaning.
   real(R8P)                 :: fs9(2,1-S_MAX:S_MAX-1,NV_MHD_GLM)      !< Split fields, GLM.
   integer(I4P)              :: m, i, j, k, c                          !< Counters.
   !$acc routine seq
   !$omp declare target

   do m=1-S_MAX, S_MAX
      call mhd_primitive_to_conservative(gamma=GAMMA, r=prim(1,m), u=prim(2,m), v=prim(3,m), w=prim(4,m), p=prim(5,m), &
                                         bx=prim(6,m), by=prim(7,m), bz=prim(8,m), q=qs8(:,m))
      do i=1, NV_MHD
         qs9(i,m) = qs8(i,m)
      enddo
      qs9(IQ_PSI,m) = prim(9,m)
      call mhd_conservative_to_auxiliary(gamma=GAMMA, R=R, q=qs8(:,m), qa=qas(:,m))
   enddo
   call mhd_eigenvalues(d=d, qa=qas(:,0), lambda=lam8)
   call mhd_glm_eigenvalues(ch=ch, d=d, qa=qas(:,0), lambda=lam9)
   call mhd_flux(d=d, q=qs8(:,0), qa=qas(:,0), f=f8)
   call mhd_glm_flux(ch=ch, d=d, q=qs9(:,0), qa=qas(:,0), f=f9)
   call mhd_eigenvectors(gamma=GAMMA, d=d, qa=qas(:,0), el=el8, er=er8)
   call mhd_glm_eigenvectors(ch=ch, gamma=GAMMA, d=d, qa=qas(:,0), el=el9, er=er9)
   c = 0
   do i=1, NV_AUX_MHD
      c = c + 1 ; res(c) = qas(i,0)
   enddo
   do i=1, NV_MHD
      c = c + 1 ; res(c) = lam8(i)
      c = c + 1 ; res(c) = f8(i)
   enddo
   do i=1, NV_MHD_GLM
      c = c + 1 ; res(c) = lam9(i)
      c = c + 1 ; res(c) = f9(i)
   enddo
   do j=1, NV_MHD
      do i=1, NV_MHD
         c = c + 1 ; res(c) = el8(i,j)
         c = c + 1 ; res(c) = er8(i,j)
      enddo
   enddo
   do j=1, NV_MHD_GLM
      do i=1, NV_MHD_GLM
         c = c + 1 ; res(c) = el9(i,j)
         c = c + 1 ; res(c) = er9(i,j)
      enddo
   enddo
   call mhd_face_split_fluxes(gamma=GAMMA, d=d, S=S, is_characteristic=.true., qs=qs8, qas=qas, fsplit=fs8, er=er8)
   call mhd_glm_face_split_fluxes(ch=ch, gamma=GAMMA, d=d, S=S, is_characteristic=.true., qs=qs9, qas=qas, &
                                  fsplit=fs9, er=er9)
   split = 0._R8P
   c = 0
   do k=1, NV_MHD
      do m=1-S_MAX, S_MAX-1
         do i=1, 2
            c = c + 1
            if (m >= 1-S .and. m <= S-1) split(c) = fs8(i,m,k)
         enddo
      enddo
   enddo
   do k=1, NV_MHD_GLM
      do m=1-S_MAX, S_MAX-1
         do i=1, 2
            c = c + 1
            if (m >= 1-S .and. m <= S-1) split(c) = fs9(i,m,k)
         enddo
      enddo
   enddo
   endsubroutine evaluate
endmodule test_flume_mhd_library_fnl_kernels

program test_flume_mhd_library_fnl
!< Unit test MV-0 (device half) of the FLUME MHD library: device results against host results.
!<
!< **Why this test exists** (issue #41, section 9, MV-0). The MHD library is shared by the CPU loops and the FNL kernels
!< (`!$acc routine seq` + `!$omp declare target`); its correctness is pinned on the host by `test_flume_mhd_library`.
!< This test pins that the device build of the SAME source returns the host results: for N random admissible stencils,
!< every direction and both variants (no cleaning, GLM), auxiliary variables, eigenvalues, fluxes and eigenvectors of
!< the central cell, and the characteristic split fields of the non-uniform stencil (S = 3, exercising the face average
!< and the per-wave speeds), computed in an offloaded loop and compared with the host evaluation, relative to each
!< quantity's magnitude. Host and device may contract differently into FMAs: the comparison is not required bitwise.

use :: adam_flume_parameters,              only : S_MAX
use :: penf,                               only : I4P, R8P, str
use :: test_flume_mhd_library_fnl_kernels, only : evaluate, GAMMA, NR, NS

implicit none

integer(I4P), parameter :: N=4096_I4P              !< Random stencils number.
real(R8P),    parameter :: TOL=1.e-12_R8P          !< Relative tolerance.
real(R8P)               :: prim(9,1-S_MAX:S_MAX,N) !< Stencils of primitive states and psi.
real(R8P)               :: ch(N)                   !< GLM cleaning speeds.
real(R8P)               :: res_dev(NR,3,N)         !< Pointwise results, device.
real(R8P)               :: res_host(NR,3,N)        !< Pointwise results, host.
real(R8P)               :: split_dev(NS,3,N)       !< Split results, device.
real(R8P)               :: split_host(NS,3,N)      !< Split results, host.
real(R8P)               :: err_res                 !< Maximum pointwise relative difference.
real(R8P)               :: err_split               !< Maximum split relative difference.
real(R8P)               :: x, a                    !< Random number, sound speed.
integer(I4P)            :: seed(64)                !< Random generator seed.
integer(I4P)            :: ns_, n_, m, d, c        !< Counters, seed size.
logical                 :: test_passed             !< Aggregate pass flag.

call random_seed(size=ns_)
if (ns_ > size(seed)) error stop 'random seed larger than expected'
seed = [(20260925_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns_))
do n_=1, N
   call random_number(x) ; ch(n_) = 0.5_R8P + 4.5_R8P * x
   do m=1-S_MAX, S_MAX
      call random_number(x) ; prim(1,m,n_) = 0.1_R8P + 9.9_R8P * x
      call random_number(x) ; prim(5,m,n_) = 0.1_R8P + 9.9_R8P * x
      a = sqrt(GAMMA * prim(5,m,n_) / prim(1,m,n_))
      do c=2, 4
         call random_number(x) ; prim(c,m,n_) = 3._R8P * a * (2._R8P * x - 1._R8P)
      enddo
      do c=6, 8
         call random_number(x) ; prim(c,m,n_) = sqrt(prim(1,m,n_)) * a * 3._R8P * (2._R8P * x - 1._R8P)
      enddo
      call random_number(x) ; prim(9,m,n_) = 2._R8P * x - 1._R8P
   enddo
enddo

!$acc parallel loop gang vector collapse(2) copyin(prim, ch) copyout(res_dev, split_dev)
!$omp target teams distribute parallel do collapse(2) map(to:prim, ch) map(from:res_dev, split_dev)
do n_=1, N
do d=1, 3
   call evaluate(d=d, ch=ch(n_), prim=prim(:,:,n_), res=res_dev(:,d,n_), split=split_dev(:,d,n_))
enddo
enddo
do n_=1, N
   do d=1, 3
      call evaluate(d=d, ch=ch(n_), prim=prim(:,:,n_), res=res_host(:,d,n_), split=split_host(:,d,n_))
   enddo
enddo

err_res   = maxval(abs(res_dev   - res_host  ) / max(1._R8P, abs(res_host  )))
err_split = maxval(abs(split_dev - split_host) / max(1._R8P, abs(split_host)))
test_passed = (err_res <= TOL) .and. (err_split <= TOL)
print '(A)', 'device vs host, pointwise (aux, eigenvalues, fluxes, eigenvectors): max relative difference '// &
             trim(str(err_res))
print '(A)', 'device vs host, split fields of random stencils (S=3):              max relative difference '// &
             trim(str(err_split))
if (test_passed) then
   print '(A)', 'TEST PASSED: flume mhd library, device vs host ('//trim(str(N))//' stencils x 3 directions x 2 variants)'
else
   print '(A)', 'TEST FAILED: flume mhd library, device vs host'
   error stop 1
endif
endprogram test_flume_mhd_library_fnl
