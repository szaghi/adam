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
!<   6b. Quintic fill (issue #68 D1): exact on per-direction degree <= 5, all octants and 6^3 shift states; an unrefined
!<      direction reads the anchor only.
!<   7. The fine->coarse restrictions (issue #68): the 6-point Lagrange (C1) and the WENO with its troubled-cell
!<      fallback to the mean (C3).

use adam_seam_interpolation_library, only : SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC, &
                                            seam_interpolate_tricubic, seam_interpolate_compatible, &
                                            seam_interpolate_quadratic, &
                                            seam_tricubic_centered_pos, seam_compatible_centered_pos, &
                                            seam_shift_anchor_pos, seam_meta_pack, seam_meta_unpack, &
                                            SEAM_W_LAGRANGE6, SEAM_RESTRICT_IDENTITY, seam_restriction_code, &
                                            seam_restriction_meta_pack, seam_restriction_meta_unpack, &
                                            SEAM_RESTRICT_LAGRANGE6, SEAM_RESTRICT_WENO6, &
                                            seam_restriction_weno6_1d, seam_restriction_weno6, &
                                            seam_interpolate_quintic, seam_quintic_centered_pos, SEAM_SUB_COINCIDENT
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
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

quintic_monomials: block
   ! #68 D1: the 6-point coarse->fine fill, exact for per-direction degree <= 5 at every octant and shift state; along an
   ! unrefined direction (sub = 3, eta = 0) only the anchor is read: the other nodes hold NaN and must not leak
   real(R8P)    :: fp6(1:6,1:6,1:6)
   real(R8P)    :: got, want, x, y, z, e(1:3)
   integer(I4P) :: a, b, c, sx, sy, sz, px, py, pz, i, j, k
   logical      :: ok

   ok = .true.
   do a=0,5 ; do b=0,5 ; do c=0,5
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do pz=1,6 ; do py=1,6 ; do px=1,6
            do k=1,6 ; do j=1,6 ; do i=1,6
               x = real(i - px, R8P) ; y = real(j - py, R8P) ; z = real(k - pz, R8P)
               fp6(i,j,k) = x**a * y**b * z**c
            enddo ; enddo ; enddo
            want = eta(sx)**a * eta(sy)**b * eta(sz)**c
            got  = seam_interpolate_quintic(fp6, [sx,sy,sz], [px,py,pz])
            ! relative to the largest footprint term: x^5 y^5 z^5 of offsets up to 5 rounds in the 3-D sum
            if (abs(got - want) > tol_exact*max(1._R8P, maxval(abs(fp6)))) ok = .false.
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   ! a quadtree row (z unrefined): every z node but the anchor is NaN; degree <= 5 in x and y, degree 0 in z
   do a=0,5 ; do b=0,5
      do sy=1,2 ; do sx=1,2
         do py=1,6 ; do px=1,6
            pz = seam_quintic_centered_pos(SEAM_SUB_COINCIDENT)
            fp6 = ieee_value(1._R8P, ieee_quiet_nan)
            do j=1,6 ; do i=1,6
               fp6(i,j,pz) = real(i - px, R8P)**a * real(j - py, R8P)**b
            enddo ; enddo
            e = [eta(sx), eta(sy), 0._R8P]
            want = e(1)**a * e(2)**b
            got  = seam_interpolate_quintic(fp6, [sx,sy,SEAM_SUB_COINCIDENT], [px,py,pz])
            if (.not.(abs(got - want) <= tol_exact*max(1._R8P, 5._R8P**(a+b)))) ok = .false.
         enddo ; enddo
      enddo ; enddo
   enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: quintic exact on all per-direction monomials deg <= 5, all octants and shift states (q=6), '// &
                     'unrefined direction reads the anchor only'
   else
      write(*,'(A)') 'FAIL: quintic missed a degree <= 5 monomial or read a node of an unrefined direction'
      test_passed = .false.
   endif
endblock quintic_monomials

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

metadata_round_trip: block
   integer(I4P) :: sub(1:3), p4(1:3), p3(1:3), p6(1:3)
   integer(I4P) :: usub(1:3), up4(1:3), up3(1:3), up6(1:3)
   integer(I4P) :: s1, s2, s3, a1, a2, a3, b1, b2, b3, c1, c2, c3, meta
   logical      :: ok

   ok = .true.
   do s3=1,3 ; do s2=1,3 ; do s1=1,3 ! sub = 3: an unrefined direction (issue #46)
      do a3=1,4 ; do a2=1,4 ; do a1=1,4
         do b3=1,3 ; do b2=1,3 ; do b1=1,3
            do c3=1,6,5 ; do c2=1,6,5 ; do c1=1,6 ! the p6 extremes per direction, every value along x
               sub = [s1,s2,s3] ; p4 = [a1,a2,a3] ; p3 = [b1,b2,b3] ; p6 = [c1,c2,c3]
               meta = seam_meta_pack(sub=sub, p4=p4, p3=p3, p6=p6)
               call seam_meta_unpack(meta=meta, sub=usub, p4=up4, p3=up3, p6=up6)
               if (any(usub /= sub) .or. any(up4 /= p4) .or. any(up3 /= p3) .or. any(up6 /= p6)) ok = .false.
               if (meta < 0) ok = .false. ! 27 bits: never the sign bit
            enddo ; enddo ; enddo
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: metadata pack/unpack round trip over all sub x p4 x p3 x p6 combinations'
   else
      write(*,'(A)') 'FAIL: metadata pack/unpack round trip corrupted a field'
      test_passed = .false.
   endif
endblock metadata_round_trip

shift_clamp: block
   integer(I4P) :: a, n, sub, p, pc
   logical      :: ok

   ok = .true.
   do n=4,16 ! every block size from the smallest legal one to the E0 size
      do a=1,n
         do sub=1,2
            ! tricubic footprint (4 nodes)
            pc = seam_tricubic_centered_pos(sub)
            p  = seam_shift_anchor_pos(anchor=a, n_cells=n, p_centered=pc, footprint_n=4_I4P)
            if (p < 1 .or. p > 4) ok = .false.
            if (a + 1 - p < 1 .or. a + 4 - p > n) ok = .false.               ! footprint inside real cells
            if (a + 1 - pc >= 1 .and. a + 4 - pc <= n .and. p /= pc) ok = .false. ! centered when it fits
            ! compatible footprint (3 nodes)
            pc = seam_compatible_centered_pos(sub)
            p  = seam_shift_anchor_pos(anchor=a, n_cells=n, p_centered=pc, footprint_n=3_I4P)
            if (p < 1 .or. p > 3) ok = .false.
            if (a + 1 - p < 1 .or. a + 3 - p > n) ok = .false.
            if (a + 1 - pc >= 1 .and. a + 3 - pc <= n .and. p /= pc) ok = .false.
            ! quintic footprint (6 nodes, #68 D1), on blocks of at least 6 cells
            if (n >= 6) then
               pc = seam_quintic_centered_pos(sub)
               p  = seam_shift_anchor_pos(anchor=a, n_cells=n, p_centered=pc, footprint_n=6_I4P)
               if (p < 1 .or. p > 6) ok = .false.
               if (a + 1 - p < 1 .or. a + 6 - p > n) ok = .false.
               if (a + 1 - pc >= 1 .and. a + 6 - pc <= n .and. p /= pc) ok = .false.
            endif
         enddo
      enddo
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: shift-inward clamp keeps every footprint inside real cells and stays centered when possible'
   else
      write(*,'(A)') 'FAIL: shift-inward clamp produced an out-of-block footprint or shifted needlessly'
      test_passed = .false.
   endif
endblock shift_clamp

restriction_lagrange6: block
   ! #68 C1: the 6-point restriction to the coarse centre, which sits midway between the fine pair (first, first + 1),
   ! the nodes first - code .. first - code + 5 (shift-inward), code 5 the identity of an unrefined direction
   real(R8P)    :: x(1:6, 0:5), got, want, d(1:6), wt, scale
   integer(I4P) :: c(1:3), a, b, e, m1, m2, m3, cx, cy, cz, meta, code(1:3), f, n
   logical      :: ok
   real(R8P)    :: rnd(1:6)

   ! node positions in fine units, the coarse centre at 0 (the fine pair at -1/2, +1/2)
   do cx=0, 4
      do m1=1, 6
         x(m1,cx) = -0.5_R8P - real(cx, R8P) + real(m1 - 1, R8P)
      enddo
   enddo
   x(:,SEAM_RESTRICT_IDENTITY) = 0._R8P ! the identity reads the cell itself (node 1), at the coarse centre

   ! 1. exactness: every monomial x^a y^b z^c, a, b, c <= 5, for every code triple (identity included), restricts to its
   !    value at the coarse centre (1 for the constant, 0 otherwise; the identity direction carries degree 0 only)
   ok = .true.
   do cz=0, 5 ; do cy=0, 5 ; do cx=0, 5
      c = [cx, cy, cz]
      do e=0, 5 ; do b=0, 5 ; do a=0, 5
         if (c(1) == SEAM_RESTRICT_IDENTITY .and. a > 0) cycle
         if (c(2) == SEAM_RESTRICT_IDENTITY .and. b > 0) cycle
         if (c(3) == SEAM_RESTRICT_IDENTITY .and. e > 0) cycle
         got = 0._R8P
         scale = 0._R8P
         do m3=1, 6 ; do m2=1, 6 ; do m1=1, 6
            wt = SEAM_W_LAGRANGE6(m1,c(1)) * SEAM_W_LAGRANGE6(m2,c(2)) * SEAM_W_LAGRANGE6(m3,c(3))
            if (wt /= 0._R8P) then
               got = got + wt * x(m1,c(1))**a * x(m2,c(2))**b * x(m3,c(3))**e
               scale = scale + abs(wt * x(m1,c(1))**a * x(m2,c(2))**b * x(m3,c(3))**e)
            endif
         enddo ; enddo ; enddo
         want = merge(1._R8P, 0._R8P, a + b + e == 0)
         ! relative to the sum of |terms|: x^15 of half-integers up to 11/2 exceeds 53 bits, so the 3-D sum rounds
         if (abs(got - want) > 1.e-14_R8P * scale) ok = .false.
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   if (ok) then
      write(*,'(A)') 'PASS: 6-point restriction exact on per-direction degree <= 5, every shift and identity code'
   else
      write(*,'(A)') 'FAIL: 6-point restriction not exact on a degree <= 5 monomial'
      test_passed = .false.
   endif

   ! 2. symmetry: the weights of code o are those of code 4 - o reversed (exactly), so on random data the mirrored
   !    footprint restricts to the same value
   ok = .true.
   do cx=0, 4
      if (any(SEAM_W_LAGRANGE6(:,cx) /= SEAM_W_LAGRANGE6(6:1:-1,4-cx))) ok = .false.
      if (abs(sum(SEAM_W_LAGRANGE6(:,cx)) - 1._R8P) > 0._R8P) ok = .false. ! partition of unity, exact (k/256)
   enddo
   call random_seed()
   do n=1, 100
      call random_number(rnd)
      d = rnd
      do cx=0, 4
         got  = sum(SEAM_W_LAGRANGE6(:,cx) * d)
         want = sum(SEAM_W_LAGRANGE6(:,4-cx) * d(6:1:-1))
         if (abs(got - want) > tol_exact) ok = .false.
      enddo
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: 6-point restriction weights mirror-symmetric (exact) and a partition of unity'
   else
      write(*,'(A)') 'FAIL: 6-point restriction weights not mirror-symmetric or not a partition of unity'
      test_passed = .false.
   endif

   ! 3. meta round trip on all 216 code triples of both regimes, never 0 (0 is the rows that keep the mean)
   ok = .true.
   do f=SEAM_RESTRICT_LAGRANGE6, SEAM_RESTRICT_WENO6
      do cz=0, 5 ; do cy=0, 5 ; do cx=0, 5
         meta = seam_restriction_meta_pack(code=[cx, cy, cz], regime=f)
         call seam_restriction_meta_unpack(meta=meta, code=code, regime=n)
         if (meta == 0 .or. any(code /= [cx, cy, cz]) .or. n /= f) ok = .false.
      enddo ; enddo ; enddo
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: restriction meta pack/unpack round trip, 0 never packed'
   else
      write(*,'(A)') 'FAIL: restriction meta pack/unpack'
      test_passed = .false.
   endif

   ! 4. the codes: identity on an unrefined direction, -1 outside the donor block, shift-inward inside it (the footprint
   !    first - code .. first - code + 5 always within 1..n, centred, code 2, when it fits)
   ok = .true.
   if (seam_restriction_code(first=1, n_cells=4, ratio=1) /= SEAM_RESTRICT_IDENTITY) ok = .false.
   do n=6, 12
      if (seam_restriction_code(first=0, n_cells=n, ratio=2) /= -1) ok = .false.
      if (seam_restriction_code(first=n, n_cells=n, ratio=2) /= -1) ok = .false.
      do f=1, n - 1
         cx = seam_restriction_code(first=f, n_cells=n, ratio=2)
         if (cx < 0 .or. cx > 4) ok = .false.
         if (f - cx < 1 .or. f - cx + 5 > n) ok = .false.
         if (f - 2 >= 1 .and. f + 3 <= n .and. cx /= 2) ok = .false.
      enddo
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: restriction codes keep the footprint inside the donor block, centred when it fits'
   else
      write(*,'(A)') 'FAIL: restriction codes'
      test_passed = .false.
   endif
endblock restriction_lagrange6

restriction_weno6: block
   ! #68 C3: the WENO restriction, per direction a convex combination of the cubic interpolants of the jump-free 4-point
   ! sub-stencils containing the coarse centre; no candidate left (a jump inside the fine pair) -> the mean of the pair
   real(R8P)    :: u(1:6), x(1:6), got, want, err(1:2), h, x0, amp, off, rnd(1:6), cube(1:6,1:6,1:6)
   real(R8P)    :: px(1:6), py(1:6), pz(1:6)
   integer(I4P) :: code, m, l, n, jmp, i, j, k
   logical      :: ok
   real(R8P), parameter :: order_min(0:4) = [3.7_R8P, 4.7_R8P, 5.7_R8P, 4.7_R8P, 3.7_R8P] !< Design order - 0.3.

   ! 1. cubic exactness, scale and offset invariant: f = x^3 + 6 x + 1 (no jump: f' > 0 and smooth), every code, under
   !    amplitudes and offsets of 1e-6..1e5; the identity code returns the cell itself
   ok = .true.
   do code=0, 4
      do m=1, 6
         x(m) = -0.5_R8P - real(code, R8P) + real(m - 1, R8P)
      enddo
      do l=-6, 5, 1
         amp = 10._R8P**l
         off = 10._R8P**(-l)
         u = amp * (x**3 + 6._R8P * x + 1._R8P) + off
         got = seam_restriction_weno6_1d(u=u, code=code)
         want = amp + off
         if (abs(got - want) > 1.e-13_R8P * (amp * maxval(abs(x**3 + 6._R8P * x + 1._R8P)) + off)) ok = .false.
      enddo
   enddo
   u = [7._R8P, -1._R8P, 3._R8P, 2._R8P, 5._R8P, 0._R8P]
   if (seam_restriction_weno6_1d(u=u, code=SEAM_RESTRICT_IDENTITY) /= 7._R8P) ok = .false.
   ! the tensor sweep, separable cubics p(x) q(y) r(z): exact (each line is a scaled cubic), identity axes included
   do code=0, 4
      do m=1, 6
         x(m) = -0.5_R8P - real(code, R8P) + real(m - 1, R8P)
      enddo
      px = x**3 + 6._R8P * x + 1._R8P
      py = 2._R8P * x**3 + x + 3._R8P
      pz = -x**3 - 4._R8P * x + 2._R8P
      do k=1, 6 ; do j=1, 6 ; do i=1, 6
         cube(i,j,k) = px(i) * py(j) * pz(k)
      enddo ; enddo ; enddo
      got = seam_restriction_weno6(code=[code, code, code], u=cube)
      if (abs(got - 6._R8P) > 1.e-12_R8P * maxval(abs(cube))) ok = .false.
      got = seam_restriction_weno6(code=[code, SEAM_RESTRICT_IDENTITY, code], u=cube)
      if (abs(got - 2._R8P * py(1)) > 1.e-12_R8P * maxval(abs(cube))) ok = .false.
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: WENO restriction exact on cubics (1-D every code, scale-invariant; 3-D tensor, identity axes)'
   else
      write(*,'(A)') 'FAIL: WENO restriction not exact on a cubic'
      test_passed = .false.
   endif

   ! 2. design order on smooth data, f = sin(2 pi x) + 0.3 cos(6 x), max error over 101 positions, h = 0.05 and 0.025:
   !    6 centred (code 2), 5 one cell from the block edge (codes 1, 3), 4 at the edge (codes 0, 4)
   ok = .true.
   do code=0, 4
      do l=1, 2
         h = 0.05_R8P / real(l, R8P)
         err(l) = 0._R8P
         do n=0, 100
            x0 = real(n, R8P) / 100._R8P
            do m=1, 6
               x(m) = x0 + h * (real(m - 1 - code, R8P) - 0.5_R8P)
            enddo
            u = sin(2._R8P * acos(-1._R8P) * x) + 0.3_R8P * cos(6._R8P * x)
            got = seam_restriction_weno6_1d(u=u, code=code)
            err(l) = max(err(l), abs(got - (sin(2._R8P * acos(-1._R8P) * x0) + 0.3_R8P * cos(6._R8P * x0))))
         enddo
      enddo
      if (log(err(1) / err(2)) / log(2._R8P) < order_min(code)) then
         ok = .false.
         write(*,'(A,I2,A,F6.2)') '   code', code, ': order ', log(err(1) / err(2)) / log(2._R8P)
      endif
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: WENO restriction design order on smooth data (6 centred, 5 and 4 towards the block edge)'
   else
      write(*,'(A)') 'FAIL: WENO restriction below its design order on smooth data'
      test_passed = .false.
   endif

   ! 3. steps: a jump in any interval but the fine pair's restricts to the fine pair's own state (exactly, up to the
   !    weights' normalisation), a jump inside the fine pair to their mean (the troubled-cell fallback); also with a
   !    density-like 1e-3 / 1 contrast. The tensor sweep troubled on every line returns the 8-cell mean.
   ok = .true.
   do code=0, 4
      do jmp=1, 5
         do m=1, 6
            u(m) = merge(1._R8P, 1.e-3_R8P, m > jmp)
         enddo
         got = seam_restriction_weno6_1d(u=u, code=code)
         if (jmp == code + 1) then
            want = 0.5_R8P * (u(code+1) + u(code+2))
         else
            want = u(code+1)
         endif
         if (abs(got - want) > tol_exact) ok = .false.
      enddo
   enddo
   px = [1._R8P, 1._R8P, 1._R8P, 2._R8P, 2._R8P, 2._R8P] ! a jump inside the centred (code 2) fine pair, nodes 3, 4
   do k=1, 6 ; do j=1, 6 ; do i=1, 6
      cube(i,j,k) = px(i) * px(j) * px(k)
   enddo ; enddo ; enddo
   if (abs(seam_restriction_weno6(code=[2, 2, 2], u=cube) - 1.5_R8P**3) > tol_exact) ok = .false.
   if (ok) then
      write(*,'(A)') 'PASS: WENO restriction takes the fine pair''s side of a jump, the mean of a jump inside the pair'
   else
      write(*,'(A)') 'FAIL: WENO restriction on a step'
      test_passed = .false.
   endif

   ! 4. mirror symmetry on random data (jumps included): code o on u equals code 4 - o on u reversed; constants exact
   ok = .true.
   do n=1, 1000
      call random_number(rnd)
      if (mod(n, 2) == 0) rnd = rnd**8 ! sharper, jump-like data
      do code=0, 4
         got  = seam_restriction_weno6_1d(u=rnd, code=code)
         u = rnd(6:1:-1)
         want = seam_restriction_weno6_1d(u=u, code=4-code)
         if (abs(got - want) > tol_exact) ok = .false.
      enddo
   enddo
   do code=0, 4
      u = 0.7_R8P
      if (abs(seam_restriction_weno6_1d(u=u, code=code) - 0.7_R8P) > tol_exact) ok = .false.
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: WENO restriction mirror-symmetric on random data, exact on constants'
   else
      write(*,'(A)') 'FAIL: WENO restriction not mirror-symmetric or not exact on constants'
      test_passed = .false.
   endif
endblock restriction_weno6

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: seam interpolation weight tables'
else
   write(*,'(A)') 'TEST FAILED: seam interpolation weight tables'
endif
if (.not. test_passed) error stop 1
endprogram test_seam_interpolation
