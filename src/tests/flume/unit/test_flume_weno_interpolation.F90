!< Unit test RV-1 of the WENO interpolation tables (issue #47, M3-P1a).
program test_flume_weno_interpolation
!< Unit test RV-1 of the WENO interpolation tables (issue #47, M3-P1a).
!<
!< **Why this test exists.** `weno_object%initialize_interpolation` adds the tables of the WENO interpolation of point
!< values (Jiang, Shu & Zhang 2013) next to the reconstruction tables of cell averages, for the M3 hybrid face flux.
!< The addition must leave the old library untouched and the new tables exact. For every upwind scheme (`weno-u-1` to
!< `weno-u-9`) and every S' <= S stored in the tables:
!<
!< 1. the reconstruction tables `a`, `p`, `d` are bitwise the ones before the call;
!< 2. each candidate polynomial reproduces the monomials of degree <= S'-1 at both interfaces, x = -1/2 and +1/2;
!< 3. the linear-weight combination reproduces the monomials of degree <= 2S'-2, and the weights are positive and sum
!<    to 1;
!< 4. `weno_reconstruct_upwind` fed with the interpolation tables converges on point values of a smooth function at
!<    order >= ORDER_MIN(S) over the spacings h/2, h/4 (printed for every scheme). The error is the maximum over NX
!<    stencil centres spread over a period: at a single point the nonlinear weights can cancel the error by chance.
!<    The bounds allow for the order lost by the Jiang-Shu weights at critical points (Henrick, Aslam & Powers 2005),
!<    a property of the kernel shared with the reconstruction, not of the tables: measured 1.00, 2.61, 7.12, 6.71, 8.79
!<    (weno-u-1 .. weno-u-9, 2026-09-30). This check pins the kernel-table wiring; checks 1-3 pin the tables.

use :: adam_globals,     only : mpih
use :: adam_weno_object, only : weno_descaling, weno_object, weno_reconstruct_upwind
use :: penf,             only : I4P, R8P, str

implicit none

character(8), parameter :: SCHEMES(5)=['weno-u-1', 'weno-u-3', 'weno-u-5', 'weno-u-7', 'weno-u-9'] !< Schemes.
real(R8P),    parameter :: ORDER_MIN(5)=[0.9_R8P, 2.0_R8P, 4.5_R8P, 6.0_R8P, 8.0_R8P] !< Minimum orders (see note 4).
real(R8P),    parameter :: TOL=1.e-14_R8P  !< Relative tolerance of the exactness checks.
real(R8P),    parameter :: H0=0.2_R8P      !< Coarsest spacing of the convergence check.
integer(I4P), parameter :: NX=64_I4P       !< Stencil centres of the convergence check.
type(weno_object)       :: weno(5)         !< WENO objects, one per scheme (`initialize` allocates once).
real(R8P), allocatable  :: a0(:,:,:)       !< Reconstruction weights before the call.
real(R8P), allocatable  :: p0(:,:,:,:)     !< Reconstruction polynomials before the call.
real(R8P), allocatable  :: d0(:,:,:,:)     !< Smoothness coefficients before the call.
real(R8P)               :: x(2)            !< Interfaces abscissae.
real(R8P)               :: comb, cand      !< Combination and candidate values.
real(R8P)               :: scale           !< Scale of a combination.
real(R8P)               :: err(3,2)        !< Convergence errors per spacing and interface.
real(R8P)               :: order(2)        !< Measured orders per interface.
real(R8P)               :: v(2,-4:4)       !< Packed stencil.
real(R8P)               :: vr(2)           !< Interpolated values.
real(R8P)               :: h               !< Spacing.
real(R8P)               :: x0              !< Centre of a convergence stencil.
integer(I4P)            :: c, S, Sk, f, deg, s1, s2, m, l, n !< Counters.
logical                 :: test_passed     !< Aggregate pass flag.

call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
test_passed = .true.
x = [-0.5_R8P, 0.5_R8P]
do c=1, size(SCHEMES)
   associate(w=>weno(c))
   call w%initialize(scheme=SCHEMES(c), nb=1, ngc=5, ni=1, nj=1, nk=1)
   S = w%S
   a0 = w%a ; p0 = w%p ; d0 = w%d
   call w%initialize_interpolation
   ! 1. old tables untouched
   if (any(w%a /= a0) .or. any(w%p /= p0) .or. any(w%d /= d0)) then
      print '(A)', SCHEMES(c)//': FAILED, the reconstruction tables changed'
      test_passed = .false.
   endif
   do Sk=1, S
      do f=1, 2
         ! 3. weights positive, summing to 1
         if (any(w%a_interp(f,0:Sk-1,Sk) <= 0._R8P) .or. abs(sum(w%a_interp(f,0:Sk-1,Sk)) - 1._R8P) > TOL) then
            print '(A)', SCHEMES(c)//': FAILED, weights of S='//trim(str(n=Sk, no_sign=.true.))//' interface '// &
                         trim(str(n=f, no_sign=.true.))
            test_passed = .false.
         endif
         do deg=0, 2*Sk-2
            comb = 0._R8P ; scale = 0._R8P
            do s1=0, Sk-1
               cand = 0._R8P
               do s2=0, Sk-1
                  cand = cand + w%p_interp(f,s2,s1,Sk) * real(s1-s2, R8P)**deg
                  scale = scale + abs(w%a_interp(f,s1,Sk) * w%p_interp(f,s2,s1,Sk) * real(s1-s2, R8P)**deg)
               enddo
               comb = comb + w%a_interp(f,s1,Sk) * cand
               ! 2. candidates exact to degree Sk-1
               if (deg <= Sk-1 .and. abs(cand - x(f)**deg) > TOL * max(1._R8P, scale)) then
                  print '(A)', SCHEMES(c)//': FAILED, candidate '//trim(str(n=s1, no_sign=.true.))//' of S='// &
                               trim(str(n=Sk, no_sign=.true.))//' at degree '//trim(str(n=deg, no_sign=.true.))
                  test_passed = .false.
               endif
            enddo
            ! 3. combination exact to degree 2Sk-2
            if (abs(comb - x(f)**deg) > TOL * max(1._R8P, scale)) then
               print '(A)', SCHEMES(c)//': FAILED, combination of S='//trim(str(n=Sk, no_sign=.true.))//' at degree '// &
                            trim(str(n=deg, no_sign=.true.))
               test_passed = .false.
            endif
         enddo
      enddo
   enddo
   ! 4. convergence of the kernel with the interpolation tables on g(x) = sin(x) + 0.3 cos(2 x)
   err = 0._R8P
   do l=1, 3
      h = H0 / 2._R8P**(l-1)
      do n=1, NX
         x0 = 2._R8P * acos(-1._R8P) * (n - 1) / NX
         do m=1-S, S-1
            v(:,m) = sin(x0 + m * h) + 0.3_R8P * cos(2._R8P * (x0 + m * h))
         enddo
         call weno_reconstruct_upwind(S=S, weno_a=w%a_interp, weno_p=w%p_interp, weno_d=w%d, weno_zeps=w%zeps, &
                                      weno_rmu=weno_descaling(S=S, weno_sigma=w%sigma, v=v(:,1-S:S-1)),      &
                                      v=v(:,1-S:S-1), vr=vr)
         do f=1, 2
            err(l,f) = max(err(l,f), abs(vr(f) - (sin(x0 + x(f) * h) + 0.3_R8P * cos(2._R8P * (x0 + x(f) * h)))))
         enddo
      enddo
   enddo
   do f=1, 2
      order(f) = log(err(2,f) / err(3,f)) / log(2._R8P)
   enddo
   print '(A)', SCHEMES(c)//': errors (h, h/2, h/4) left '//trim(str(n=err(1,1)))//' '//trim(str(n=err(2,1)))//' '// &
                trim(str(n=err(3,1)))//' order '//trim(str(n=order(1)))
   print '(A)', SCHEMES(c)//': errors (h, h/2, h/4) right '//trim(str(n=err(1,2)))//' '//trim(str(n=err(2,2)))//' '// &
                trim(str(n=err(3,2)))//' order '//trim(str(n=order(2)))
   if (any(order < ORDER_MIN(c))) then
      print '(A)', SCHEMES(c)//': FAILED, order below '//trim(str(n=ORDER_MIN(c)))
      test_passed = .false.
   endif
   endassociate
enddo
if (test_passed) then
   print '(A)', 'TEST PASSED: WENO interpolation tables (old tables untouched, exactness, weights, convergence)'
else
   print '(A)', 'TEST FAILED: WENO interpolation tables'
endif
call mpih%finalize
if (.not.test_passed) error stop 1
endprogram test_flume_weno_interpolation
