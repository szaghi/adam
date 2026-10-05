!< Unit test RV-1 (device half) of the WENO interpolation tables: device kernel against host kernel.
module test_flume_weno_interpolation_fnl_kernels
!< Unit test RV-1 (device half) of the WENO interpolation tables: the device WENO primitive on one packed stencil.
!<
!< Module-level (not a `contains`-internal procedure): device routines must be module procedures.

use :: adam_fnl_weno_kernels, only : weno_reconstruct_upwind_dev
use :: adam_weno_object,      only : weno_descaling
use :: penf,                  only : I4P, R8P

implicit none
private
public :: interpolate

contains
   subroutine interpolate(S, a, p, d, zeps, sigma, v, vr)
   !< Interpolate one packed stencil with the device WENO primitive.
   integer(I4P), intent(in)  :: S              !< Number of stencils used.
   real(R8P),    intent(in)  :: a(1:,0:,1:)    !< Interpolation optimal weights.
   real(R8P),    intent(in)  :: p(1:,0:,0:,1:) !< Interpolation polynomials coefficients.
   real(R8P),    intent(in)  :: d(0:,0:,0:,1:) !< Smoothness indicators coefficients.
   real(R8P),    intent(in)  :: zeps           !< Parameter avoiding division by zero.
   real(R8P),    intent(in)  :: sigma          !< Descaler switch: 0 Jiang-Shu, 1 scale-invariant.
   real(R8P),    intent(in)  :: v(1:2,1-S:S-1) !< Packed stencil.
   real(R8P),    intent(out) :: vr(1:2)        !< Interface values.
   !$acc routine seq
   !$omp declare target

   call weno_reconstruct_upwind_dev(S=S, weno_a=a, weno_p=p, weno_d=d, weno_zeps=zeps,                      &
                                    weno_rmu=weno_descaling(S=S, weno_sigma=sigma, v=v), V=v, VR=vr)
   endsubroutine interpolate
endmodule test_flume_weno_interpolation_fnl_kernels

program test_flume_weno_interpolation_fnl
!< Unit test RV-1 (device half) of the WENO interpolation tables: device kernel against host kernel.
!<
!< **Why this test exists** (issue #47, M3-P1a). The hybrid face flux of M3 calls the device WENO primitive
!< `weno_reconstruct_upwind_dev` with the interpolation tables `a_interp`, `p_interp` in place of the reconstruction
!< ones. The host half (`test_flume_weno_interpolation`) pins the tables; this test pins that the device primitive fed
!< with them returns the host primitive's values: for every upwind scheme and N random stencils of point values (smooth
!< and rough), both interfaces, relative to each value's magnitude. Host and device may contract differently into FMAs,
!< so the comparison is not required to be bitwise. Both weights are checked, Jiang-Shu and scale-invariant (issue #49):
!< for the latter the device must also satisfy WENO(2**n v) = 2**n WENO(v) bit for bit (NV-1, device half).

use :: adam_globals,                              only : mpih
use :: adam_weno_object,                          only : weno_descaling, weno_object, weno_reconstruct_upwind
use :: penf,                                      only : I4P, R8P, str
use :: test_flume_weno_interpolation_fnl_kernels, only : interpolate

implicit none

character(8), parameter :: SCHEMES(5)=['weno-u-1', 'weno-u-3', 'weno-u-5', 'weno-u-7', 'weno-u-9'] !< Schemes.
integer(I4P), parameter :: N=4096_I4P       !< Random stencils number.
real(R8P),    parameter :: TOL=1.e-13_R8P   !< Relative tolerance.
type(weno_object)       :: weno(5)          !< WENO objects, one per scheme (`initialize` allocates once).
real(R8P), allocatable  :: a(:,:,:)         !< Interpolation optimal weights (plain copy for the data clauses).
real(R8P), allocatable  :: p(:,:,:,:)       !< Interpolation polynomials coefficients (plain copy).
real(R8P), allocatable  :: d(:,:,:,:)       !< Smoothness indicators coefficients (plain copy).
real(R8P), allocatable  :: v(:,:,:)         !< Packed stencils.
real(R8P), allocatable  :: vr_dev(:,:)      !< Interface values, device.
real(R8P), allocatable  :: vr_host(:,:)     !< Interface values, host.
real(R8P), allocatable  :: v_scal(:,:,:)    !< Packed stencils times LAMBDA.
real(R8P), allocatable  :: vr_scal(:,:)     !< Interface values of the scaled stencils, device.
real(R8P)               :: sigma            !< Descaler switch.
real(R8P), parameter    :: LAMBDA=2._R8P**(-40) !< Scaling of the device covariance check.
real(R8P)               :: zeps             !< Parameter avoiding division by zero.
real(R8P)               :: err              !< Maximum relative difference.
real(R8P)               :: x                !< Random number.
integer(I4P)            :: seed(64)         !< Random generator seed.
integer(I4P)            :: ns_, n_, m, c, S !< Counters, seed size.
integer(I4P)            :: wt               !< Weights counter: 1 Jiang-Shu, 2 scale-invariant.
logical                 :: test_passed      !< Aggregate pass flag.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
call random_seed(size=ns_)
if (ns_ > size(seed)) error stop 'random seed larger than expected'
seed = [(20260930_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns_))
test_passed = .true.
do c=1, size(SCHEMES)
   call weno(c)%initialize(scheme=SCHEMES(c), nb=1, ngc=5, ni=1, nj=1, nk=1)
   call weno(c)%initialize_interpolation
   S = weno(c)%S
   a = weno(c)%a_interp ; p = weno(c)%p_interp ; d = weno(c)%d ; zeps = weno(c)%zeps
   if (allocated(v)) deallocate(v, v_scal, vr_dev, vr_host, vr_scal)
   allocate(v(2,1-S:S-1,N), v_scal(2,1-S:S-1,N), vr_dev(2,N), vr_host(2,N), vr_scal(2,N))
   do n_=1, N
      call random_number(x)
      if (x < 0.5_R8P) then ! smooth: samples of a random quadratic
         call random_number(x)
         do m=1-S, S-1
            v(:,m,n_) = 1._R8P + x * m + 0.1_R8P * x * m * m
         enddo
      else                  ! rough: independent random values
         do m=1-S, S-1
            call random_number(x) ; v(1,m,n_) = 10._R8P * x - 5._R8P
            call random_number(x) ; v(2,m,n_) = 10._R8P * x - 5._R8P
         enddo
      endif
   enddo
   v_scal = LAMBDA * v
   do wt=1, 2
   sigma = merge(0._R8P, 1._R8P, wt == 1)
   !$acc parallel loop gang vector copyin(a, p, d, v, v_scal) copyout(vr_dev, vr_scal)
   !$omp target teams distribute parallel do map(to:a, p, d, v, v_scal) map(from:vr_dev, vr_scal)
   do n_=1, N
      call interpolate(S=S, a=a, p=p, d=d, zeps=zeps, sigma=sigma, v=v(:,:,n_), vr=vr_dev(:,n_))
      call interpolate(S=S, a=a, p=p, d=d, zeps=zeps, sigma=sigma, v=v_scal(:,:,n_), vr=vr_scal(:,n_))
   enddo
   do n_=1, N
      call weno_reconstruct_upwind(S=S, weno_a=a, weno_p=p, weno_d=d, weno_zeps=zeps,                            &
                                   weno_rmu=weno_descaling(S=S, weno_sigma=sigma, v=v(:,:,n_)), v=v(:,:,n_), &
                                   vr=vr_host(:,n_))
   enddo
   err = maxval(abs(vr_dev - vr_host) / max(1._R8P, abs(vr_host)))
   print '(A)', SCHEMES(c)//merge(' Jiang-Shu      ', ' scale-invariant', wt == 1)//': device vs host, max relative '// &
                'difference '//trim(str(n=err))//', stencils breaking WENO(2^-40 v) = 2^-40 WENO(v) on device '//      &
                trim(str(n=count(any(vr_scal /= LAMBDA * vr_dev, dim=1)), no_sign=.true.))
   if (err > TOL) test_passed = .false.
   if (wt == 2 .and. any(vr_scal /= LAMBDA * vr_dev)) then
      print '(A)', SCHEMES(c)//': FAILED, the scale-invariant weights are not scale-covariant on device'
      test_passed = .false.
   endif
   enddo
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: WENO interpolation, device vs host ('//trim(str(n=N, no_sign=.true.))// &
                ' stencils x 5 schemes x 2 weights), scale-invariant weights scale-covariant on device'
else
   print '(A)', 'TEST FAILED: WENO interpolation, device vs host'
endif
call mpih%finalize
if (.not.test_passed) error stop 1
endprogram test_flume_weno_interpolation_fnl
