!< Unit test NV-1 of the scale-invariant WENO weights (issue #49, N1).
program test_flume_weno_weights
!< Unit test NV-1 of the scale-invariant WENO weights (issue #49, N1).
!<
!< **Why this test exists.** The Jiang-Shu weights `d / (zeps + IS)**wexp` compare the smoothness indicators with an
!< absolute `zeps`, so the weights depend on the units of the data: at small magnitudes they collapse to the linear
!< weights and the scheme stops limiting. `[weno] weights = si` divides the indicators by the square of the descaler
!< `mu = mean |v|` of the global stencil (`seth_weno_descaling`, the generic descaler; FLUME passes the magnitude of the
!< projected state instead, which scales the same way), so the reconstruction is homogeneous of degree one in the data.
!< Scaling by a power of two is exact in floating point, so the property is checked bit for bit. Since #139 (S6) the
!< kernels are SETH's, called through the adapter's dispatchers (`weno_reconstruct`, `weno_reconstruct_interp_wr`): the
!< test checks what FLUME relies on. For every upwind scheme (`weno-u-1` to `weno-u-9`), reconstruction and interpolation, N random
!< stencils (smooth and rough, magnitudes from 1e-12 to 1e12) and scalings 2**n, n in NS:
!<
!< 1. scale-invariant weights: WENO(2**n v) = 2**n WENO(v) bitwise, every stencil, both interfaces;
!< 2. Jiang-Shu weights (negative control): the same identity fails on some stencil, otherwise check 1 proves nothing;
!< 3. on quadratic data, which every candidate of S >= 3 reproduces, both weights give the same value to round-off at
!<    any magnitude: the scale-invariant weights stay finite and normalised where the descaler is tiny;
!< 4. data identically zero (descaler floored at `tiny`) reconstruct to zero, not to a NaN.

use :: adam_globals,     only : mpih
use :: adam_weno_object, only : weno_object, weno_reconstruct, weno_reconstruct_interp_wr, WENO_WEIGHTS_JS
use :: penf,             only : I4P, R8P, str
use :: seth,             only : seth_weno_descaling

implicit none

character(8), parameter :: SCHEMES(5)=['weno-u-1', 'weno-u-3', 'weno-u-5', 'weno-u-7', 'weno-u-9'] !< Schemes.
integer(I4P), parameter :: NS(6)=[-60_I4P, -20_I4P, -1_I4P, 1_I4P, 20_I4P, 60_I4P] !< Exponents of the scalings.
integer(I4P), parameter :: N=2048_I4P      !< Random stencils number.
real(R8P),    parameter :: TOL_QUAD=1.e-12_R8P !< Agreement of the two weights on quadratic data, relative.
type(weno_object)       :: weno(5)         !< WENO objects, one per scheme.
real(R8P)               :: v(2,-4:4)       !< Packed stencil.
real(R8P)               :: vr(2)           !< Interface values.
real(R8P)               :: vs(2)           !< Interface values of the scaled stencil.
real(R8P)               :: vj(2)           !< Interface values, Jiang-Shu weights.
real(R8P)               :: x, mag          !< Random number, magnitude.
real(R8P)               :: lambda          !< Scaling.
real(R8P)               :: sigma           !< Descaler switch.
real(R8P)               :: dev             !< Maximum relative difference of the two weights on quadratic data.
integer(I4P)            :: seed(64)        !< Random generator seed.
integer(I4P)            :: broken(2)       !< Stencils breaking the identity: scale-invariant, Jiang-Shu.
integer(I4P)            :: ns_, n_, m, c, S, t, l, wt !< Counters, seed size.
logical                 :: smooth          !< Smooth (quadratic) stencil flag.
logical                 :: test_passed     !< Aggregate pass flag.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
call random_seed(size=ns_)
if (ns_ > size(seed)) error stop 'random seed larger than expected'
seed = [(20261004_I4P + n_, n_=1, size(seed))]
call random_seed(put=seed(1:ns_))
test_passed = .true.
do c=1, size(SCHEMES)
   call weno(c)%initialize(scheme=SCHEMES(c))
   if (weno(c)%weights /= WENO_WEIGHTS_JS .or. weno(c)%sigma /= 0._R8P) then
      print '(A)', SCHEMES(c)//': FAILED, the default weights are not Jiang-Shu'
      test_passed = .false.
   endif
   S = weno(c)%S
   do t=1, 2 ! 1 reconstruction tables, 2 interpolation tables
      broken = 0_I4P
      dev = 0._R8P
      do n_=1, N
         call random_number(x)
         smooth = x < 0.5_R8P
         call random_number(x)
         mag = 10._R8P**(24._R8P * x - 12._R8P)
         call random_number(x)
         do m=1-S, S-1
            if (smooth) then ! samples of a random quadratic, offset or crossing zero
               v(:,m) = mag * (x - 0.5_R8P + 0.1_R8P * x * m + 0.01_R8P * x * m * m)
            else             ! independent random values of both signs
               call random_number(x) ; v(1,m) = mag * (2._R8P * x - 1._R8P)
               call random_number(x) ; v(2,m) = mag * (2._R8P * x - 1._R8P)
            endif
         enddo
         do wt=1, 2 ! 1 scale-invariant, 2 Jiang-Shu
            sigma = merge(1._R8P, 0._R8P, wt == 1)
            call reconstruct(v=v(:,1-S:S-1), vr=vr)
            do l=1, size(NS)
               lambda = 2._R8P**NS(l)
               call reconstruct(v=lambda * v(:,1-S:S-1), vr=vs)
               if (any(vs /= lambda * vr)) then
                  broken(wt) = broken(wt) + 1_I4P
                  exit
               endif
            enddo
         enddo
         if (smooth .and. S >= 3) then
            call reconstruct_si_js(v=v(:,1-S:S-1), vr=vr, vj=vj)
            dev = max(dev, maxval(abs(vr - vj)) / maxval(abs(v(:,1-S:S-1))))
         endif
      enddo
      print '(A)', SCHEMES(c)//merge(': reconstruction', ': interpolation ', t == 1)//' tables, stencils breaking '// &
                   'WENO(2^n v) = 2^n WENO(v): scale-invariant '//trim(str(n=broken(1), no_sign=.true.))//           &
                   ', Jiang-Shu '//trim(str(n=broken(2), no_sign=.true.))//' of '//trim(str(n=N, no_sign=.true.))//   &
                   '; quadratic data, max relative difference of the two weights '//trim(str(n=dev))
      ! 1. scale-invariant weights covariant bit for bit
      if (broken(1) > 0_I4P) then
         print '(A)', SCHEMES(c)//': FAILED, the scale-invariant weights are not scale-covariant'
         test_passed = .false.
      endif
      ! 2. negative control: the Jiang-Shu weights are not (weno-u-1 has one weight and is linear)
      if (S > 1 .and. broken(2) == 0_I4P) then
         print '(A)', SCHEMES(c)//': FAILED, the Jiang-Shu weights pass too: the check proves nothing'
         test_passed = .false.
      endif
      ! 3. both weights reproduce quadratic data
      if (dev > TOL_QUAD) then
         print '(A)', SCHEMES(c)//': FAILED, the two weights disagree on quadratic data'
         test_passed = .false.
      endif
   enddo
   ! 4. zero data
   v = 0._R8P
   sigma = 1._R8P
   do t=1, 2
      call reconstruct(v=v(:,1-S:S-1), vr=vr)
      if (any(vr /= 0._R8P)) then
         print '(A)', SCHEMES(c)//': FAILED, zero data do not reconstruct to zero'
         test_passed = .false.
      endif
   enddo
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: scale-invariant WENO weights (NV-1 host, '//trim(str(n=N, no_sign=.true.))// &
                ' stencils x 5 schemes x 2 tables)'
else
   print '(A)', 'TEST FAILED: scale-invariant WENO weights'
endif
call mpih%finalize
if (.not.test_passed) error stop 1
contains
   subroutine reconstruct(v, vr)
   !< Reconstruct (`t` 1) or interpolate (`t` 2) one stencil of scheme `c` with the weights `sigma`, SETH's kernels
   !< through the adapter; the packed rows map to SETH's stencils as `vp(m) = v(2,m)`, `vm(m+1) = v(1,m)`.
   real(R8P), intent(in)  :: v(1:2,1-S:S-1) !< Packed stencil.
   real(R8P), intent(out) :: vr(1:2)        !< Interface values.
   real(R8P)              :: vp(1-S:S), vm(1-S:S) !< SETH's stencils.
   real(R8P)              :: rmu, wrp, wrm     !< Inverse descaler, weights ratios (unused).
   integer(I4P)           :: k                 !< Counter.

   vp = 0._R8P ; vm = 0._R8P
   do k=1-S, S-1
      vp(k)   = v(2,k)
      vm(k+1) = v(1,k)
   enddo
   rmu = seth_weno_descaling(S=S, sigma=sigma, vp=vp, vm=vm)
   if (t == 1) then
      call weno_reconstruct(S=S, vp=vp, vm=vm, eps=weno(c)%eps, rmu=rmu, rp=vr(2), rm=vr(1))
   else
      call weno_reconstruct_interp_wr(S=S, vp=vp, vm=vm, eps=weno(c)%eps, rmu=rmu, rp=vr(2), rm=vr(1), wrp=wrp, wrm=wrm)
   endif
   endsubroutine reconstruct

   subroutine reconstruct_si_js(v, vr, vj)
   !< Reconstruct one stencil with both weights.
   real(R8P), intent(in)  :: v(1:2,1-S:S-1) !< Packed stencil.
   real(R8P), intent(out) :: vr(1:2)        !< Interface values, scale-invariant weights.
   real(R8P), intent(out) :: vj(1:2)        !< Interface values, Jiang-Shu weights.

   sigma = 1._R8P ; call reconstruct(v=v, vr=vr)
   sigma = 0._R8P ; call reconstruct(v=v, vr=vj)
   endsubroutine reconstruct_si_js
endprogram test_flume_weno_weights
