!< Unit test of the 2:1 coarse-fine flux restriction (`restrict_fine_face_to_quadrant`).
program test_amr_reflux_restrict
!< Unit test of `prism_cpu_object`'s pure 2:1 flux-restriction helper
!< `restrict_fine_face_to_quadrant` — the arithmetic that maps a fine block's
!< face flux onto its quadrant of the coarse-face skin for Berger-Colella reflux
!< at an intra-realm AMR jump (issue #13 §7.5, milestone M3).
!<
!< **Why this test exists.** M3's new and error-prone code is the restriction:
!< placing each of the 4 fine blocks into the correct disjoint quadrant of the
!< coarse face and conservatively averaging the 2x2 fine cells under each coarse
!< cell. A wrong quadrant offset or a wrong linear index silently mis-accumulates
!< F_fine_sum and corrupts the reflux correction. This test pins the helper
!< directly with synthetic fluxes, with no solve, no INI, no IO.
!<
!< **What it pins:**
!<   1. Conservative averaging: a coarse cell gets the arithmetic mean of the
!<      2x2 fine cells under it (a constant fine field restricts to that constant).
!<   2. Quadrant placement: each (ioff,joff) writes exactly its inner_n/2 * outer_n/2
!<      block of coarse cells, at the right linear-index offset, and nowhere else.
!<   3. Disjoint tiling: the 4 quadrants together cover every coarse cell exactly
!<      once (no gap, no overlap).

use adam_prism_cpu_object, only : restrict_fine_face_to_quadrant
use adam_globals, only : mpih
use penf, only : I4P, R8P

implicit none

integer(I4P), parameter :: nv      = 2_I4P            !< State-vector width.
integer(I4P), parameter :: inner_n = 4_I4P            !< Coarse-face inner cell count (so quadrant = 2).
integer(I4P), parameter :: outer_n = 4_I4P            !< Coarse-face outer cell count.
integer(I4P), parameter :: nface   = inner_n*outer_n  !< Coarse-skin cell count (16).
real(R8P),    parameter :: tol     = 1.0e-13_R8P      !< Round-off tolerance.

real(R8P) :: fine_face(1:nv, 1:inner_n, 1:outer_n) !< One fine block's face flux.
real(R8P) :: slab(1:nv, 1:nface)                   !< Coarse-skin slab (filled by the helper).
real(R8P) :: cover(1:nface)                        !< How many quadrants wrote each coarse cell.
integer(I4P) :: ioff, joff, fi, fo, c, ic, oc, v
logical :: test_passed

test_passed = .true.
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)

! --- Test 1: constant fine field restricts to the same constant on its quadrant. ---
! A fine face that is uniformly 7.0 must average to 7.0 in every covered coarse cell.
fine_face = 7.0_R8P
slab = -1.0_R8P  ! poison: untouched cells must remain at the poison value
call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, &
                                    ioff=0_I4P, joff=0_I4P, slab=slab)
block
   logical :: ok
   ok = .true.
   ! Quadrant (0,0) covers coarse cells (ic=1..2, oc=1..2): linear c=(oc-1)*inner_n+ic.
   do oc = 1_I4P, outer_n/2_I4P
      do ic = 1_I4P, inner_n/2_I4P
         c = (oc - 1_I4P)*inner_n + ic
         do v = 1_I4P, nv
            if (abs(slab(v, c) - 7.0_R8P) > tol) ok = .false.
         enddo
      enddo
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: constant fine field restricts to the same constant (conservative mean)'
   else
      write(*,'(A)') 'FAIL: constant-field restriction did not reproduce the constant'
      test_passed = .false.
   endif
endblock

! --- Test 2: exact 2x2 average of a non-constant field, in the correct cell. ---
! Set distinct values in the 2x2 fine block under coarse cell (ic=1,oc=1) of quadrant (1,0).
fine_face = 0.0_R8P
fine_face(1, 1, 1) = 1.0_R8P ; fine_face(1, 2, 1) = 2.0_R8P
fine_face(1, 1, 2) = 3.0_R8P ; fine_face(1, 2, 2) = 4.0_R8P  ! mean = 2.5
slab = 0.0_R8P
call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, &
                                    ioff=1_I4P, joff=0_I4P, slab=slab)
block
   integer(I4P) :: c_expected
   ! Quadrant (1,0): inner offset = ioff*inner_n/2 = 2, outer offset = 0.
   ! Coarse cell (ic=1,oc=1) within the quadrant → global (ic=3, oc=1).
   c_expected = (0_I4P + 1_I4P - 1_I4P)*inner_n + (2_I4P + 1_I4P)  ! = 3
   if (abs(slab(1, c_expected) - 2.5_R8P) <= tol) then
      write(*,'(A)') 'PASS: 2x2 conservative mean = 2.5 lands in the correct quadrant cell'
   else
      write(*,'(A,F8.4,A,I0)') 'FAIL: expected 2.5 at c=', 2.5_R8P, ' got slab(1,c) at c=', c_expected
      write(*,'(A,F12.8)') '      slab value = ', slab(1, c_expected)
      test_passed = .false.
   endif
endblock

! --- Test 3: the 4 quadrants tile the coarse face exactly once (no gap/overlap). ---
cover = 0.0_R8P
do joff = 0_I4P, 1_I4P
   do ioff = 0_I4P, 1_I4P
      fine_face = 1.0_R8P
      slab = 0.0_R8P
      call restrict_fine_face_to_quadrant(fine_face=fine_face, inner_n=inner_n, outer_n=outer_n, &
                                          ioff=ioff, joff=joff, slab=slab)
      ! Each quadrant writes 1.0 into its cells; accumulate a coverage count.
      do c = 1_I4P, nface
         if (abs(slab(1, c) - 1.0_R8P) <= tol) cover(c) = cover(c) + 1.0_R8P
      enddo
   enddo
enddo
block
   logical :: ok
   ok = .true.
   do c = 1_I4P, nface
      if (abs(cover(c) - 1.0_R8P) > tol) ok = .false.   ! exactly one quadrant per cell
   enddo
   if (ok) then
      write(*,'(A)') 'PASS: the 4 quadrants tile every coarse cell exactly once (disjoint, complete)'
   else
      write(*,'(A)') 'FAIL: quadrant tiling has a gap or overlap'
      test_passed = .false.
   endif
endblock

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: amr reflux restrict'
else
   write(*,'(A)') 'TEST FAILED: amr reflux restrict'
endif

call mpih%finalize
if (.not. test_passed) error stop 1
endprogram test_amr_reflux_restrict
