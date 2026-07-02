!< Unit test of the seam coarse->fine ghost-fill weight tables (issue #21 N0).
program test_seam_interpolation
!< Unit test of `adam_seam_interpolation_library` — the 1D weight tables and
!< tensor evaluators for the 2:1 coarse->fine seam ghost fill (issue #21, plan
!< amendment v2, milestone N0).
!<
!< **Why this test exists.** The E0' gate (N3) can only attribute a seam
!< div(B) slope change to the FILL MECHANISM if the weights themselves are
!< beyond suspicion. This test pins them with no maps, no MPI, no INI:
!<   1. Partition of unity: every regime/state reproduces constants (exactly
!<      representable weights => near-0-ulp).
!<   2. Tricubic exactness: every per-direction monomial x^a y^b z^c with
!<      a,b,c <= 3, at all 8 octants and all 4^3 anchor (shift) states => q=4.
!<   3. Compatible-arm exactness on linears, all 3^3 anchor states => q=2.
!<   4. Compatible arm, THE DEFINING PROPERTY: exact restriction round trip
!<      R o P = I on RANDOM data — the 8-octant mean of the interpolants
!<      equals the anchor value to round-off, for every anchor state.
!<   5. THE INCOMPATIBILITY THEOREM, executably pinned (G1 amendment,
!<      #21 issuecomment-4869444520): the tricubic fill of f = x^2 gives every
!<      octant exactly 1/16 (correct interpolation!), whose mean minus the
!<      anchor value f(0)=0 is 1/16 =/= 0 — restriction-compatibility and
!<      q>2 are mutually exclusive for point values. This asserts the defect
!<      EQUALS H^2/16: it is expected behavior, not a failure.
!<   6. Quadratic safety set: exact on per-direction degree <= 2, centered.

use adam_seam_interpolation_library, only : SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC, &
                                            seam_interpolate_tricubic, seam_interpolate_compatible, &
                                            seam_interpolate_quadratic, &
                                            seam_tricubic_centered_pos, seam_compatible_centered_pos
use penf, only : I4P, R8P

implicit none

real(R8P), parameter :: tol_exact = 1.e-14_R8P !< Tolerance for identities exact up to accumulation round-off.
real(R8P), parameter :: eta(1:2) = [-0.25_R8P, +0.25_R8P] !< Octant sub-position coordinates.
logical              :: test_passed             !< Aggregate pass flag.

test_passed = .true.

partition_of_unity: block
   real(R8P)    :: fp4(1:4,1:4,1:4), fp3(1:3,1:3,1:3)
   real(R8P)    :: got
   integer(I4P) :: sx, sy, sz, px, py, pz
   logical      :: ok

   ok = .true.
   fp4 = 1._R8P
   fp3 = 1._R8P
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      do pz=1,4 ; do py=1,4 ; do px=1,4
         got = seam_interpolate_tricubic(fp4, [sx,sy,sz], [px,py,pz])
         if (abs(got - 1._R8P) > tol_exact) ok = .false.
      enddo ; enddo ; enddo
      do pz=1,3 ; do py=1,3 ; do px=1,3
         got = seam_interpolate_compatible(fp3, [sx,sy,sz], [px,py,pz])
         if (abs(got - 1._R8P) > tol_exact) ok = .false.
      enddo ; enddo ; enddo
      got = seam_interpolate_quadratic(fp3, [sx,sy,sz])
      if (abs(got - 1._R8P) > tol_exact) ok = .false.
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: partition of unity (constants reproduced by every regime/state)'
   else
      write(*,'(A)') 'FAIL: a weight set does not sum to 1'
      test_passed = .false.
   endif
endblock partition_of_unity

tricubic_monomials: block
   real(R8P)    :: fp4(1:4,1:4,1:4)
   real(R8P)    :: got, want, x, y, z
   integer(I4P) :: a, b, c, sx, sy, sz, px, py, pz, i, j, k
   logical      :: ok

   ok = .true.
   do a=0,3 ; do b=0,3 ; do c=0,3
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do pz=1,4 ; do py=1,4 ; do px=1,4
            do k=1,4 ; do j=1,4 ; do i=1,4
               x = real(i - px, R8P) ; y = real(j - py, R8P) ; z = real(k - pz, R8P)
               fp4(i,j,k) = x**a * y**b * z**c
            enddo ; enddo ; enddo
            want = eta(sx)**a * eta(sy)**b * eta(sz)**c
            got  = seam_interpolate_tricubic(fp4, [sx,sy,sz], [px,py,pz])
            if (abs(got - want) > tol_exact*max(1._R8P, abs(want))) ok = .false.
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: tricubic exact on all per-direction monomials deg <= 3, all octants, all shift states (q=4)'
   else
      write(*,'(A)') 'FAIL: tricubic missed a degree <= 3 monomial'
      test_passed = .false.
   endif
endblock tricubic_monomials

compatible_linears: block
   real(R8P)    :: fp3(1:3,1:3,1:3)
   real(R8P)    :: got, want, x, y, z
   integer(I4P) :: a, b, c, sx, sy, sz, px, py, pz, i, j, k
   logical      :: ok

   ok = .true.
   do a=0,1 ; do b=0,1 ; do c=0,1
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do pz=1,3 ; do py=1,3 ; do px=1,3
            do k=1,3 ; do j=1,3 ; do i=1,3
               x = real(i - px, R8P) ; y = real(j - py, R8P) ; z = real(k - pz, R8P)
               fp3(i,j,k) = x**a * y**b * z**c
            enddo ; enddo ; enddo
            want = eta(sx)**a * eta(sy)**b * eta(sz)**c
            got  = seam_interpolate_compatible(fp3, [sx,sy,sz], [px,py,pz])
            if (abs(got - want) > tol_exact*max(1._R8P, abs(want))) ok = .false.
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: compatible arm exact on linears, all octants, all shift states (q=2)'
   else
      write(*,'(A)') 'FAIL: compatible arm missed a linear monomial'
      test_passed = .false.
   endif
endblock compatible_linears

compatible_round_trip: block
   real(R8P)    :: fp3(1:3,1:3,1:3)
   real(R8P)    :: mean8, anchor
   integer(I4P) :: sx, sy, sz, px, py, pz
   logical      :: ok

   ok = .true.
   call random_seed()
   do pz=1,3 ; do py=1,3 ; do px=1,3
      call random_number(fp3)
      anchor = fp3(px,py,pz)
      mean8  = 0._R8P
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         mean8 = mean8 + seam_interpolate_compatible(fp3, [sx,sy,sz], [px,py,pz])
      enddo ; enddo ; enddo
      mean8 = mean8 * 0.125_R8P
      if (abs(mean8 - anchor) > tol_exact) ok = .false.
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: compatible arm satisfies exact R o P = I on random data (8-octant mean = anchor)'
   else
      write(*,'(A)') 'FAIL: compatible arm broke the restriction round trip'
      test_passed = .false.
   endif
endblock compatible_round_trip

incompatibility_theorem: block
   real(R8P)    :: fp4(1:4,1:4,1:4)
   real(R8P)    :: got, mean8, x
   integer(I4P) :: sx, sy, sz, px, py, pz, i, j, k
   logical      :: ok

   ! f = x^2, centered per-sub anchor states (the production kernel convention).
   ok    = .true.
   mean8 = 0._R8P
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      px = seam_tricubic_centered_pos(sx)
      py = seam_tricubic_centered_pos(sy)
      pz = seam_tricubic_centered_pos(sz)
      do k=1,4 ; do j=1,4 ; do i=1,4
         x = real(i - px, R8P)
         fp4(i,j,k) = x*x
      enddo ; enddo ; enddo
      got = seam_interpolate_tricubic(fp4, [sx,sy,sz], [px,py,pz])
      ! correct interpolation: each octant must get exactly eta^2 = 1/16
      if (abs(got - 0.0625_R8P) > tol_exact) ok = .false.
      mean8 = mean8 + got
   enddo ; enddo ; enddo
   mean8 = mean8 * 0.125_R8P
   ! the theorem: mean - f(anchor) = H^2/16 exactly (anchor value is 0)
   if (abs(mean8 - 0.0625_R8P) > tol_exact) ok = .false.
   if (ok) then
      write(*,'(A)') 'PASS: incompatibility theorem pinned (tricubic octants of x^2 exact = 1/16; mean defect = H^2/16)'
   else
      write(*,'(A)') 'FAIL: tricubic x^2 octant values or mean defect differ from the theorem'
      test_passed = .false.
   endif
endblock incompatibility_theorem

quadratic_safety: block
   real(R8P)    :: fp3(1:3,1:3,1:3)
   real(R8P)    :: got, want, x, y, z
   integer(I4P) :: a, b, c, sx, sy, sz, i, j, k
   logical      :: ok

   ok = .true.
   do a=0,2 ; do b=0,2 ; do c=0,2
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do k=1,3 ; do j=1,3 ; do i=1,3
            x = real(i - 2, R8P) ; y = real(j - 2, R8P) ; z = real(k - 2, R8P)
            fp3(i,j,k) = x**a * y**b * z**c
         enddo ; enddo ; enddo
         want = eta(sx)**a * eta(sy)**b * eta(sz)**c
         got  = seam_interpolate_quadratic(fp3, [sx,sy,sz])
         if (abs(got - want) > tol_exact*max(1._R8P, abs(want))) ok = .false.
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: quadratic safety set exact on per-direction monomials deg <= 2, centered (q=3)'
   else
      write(*,'(A)') 'FAIL: quadratic safety set missed a degree <= 2 monomial'
      test_passed = .false.
   endif
endblock quadratic_safety

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: seam interpolation weight tables'
else
   write(*,'(A)') 'TEST FAILED: seam interpolation weight tables'
endif
if (.not. test_passed) error stop 1
endprogram test_seam_interpolation
