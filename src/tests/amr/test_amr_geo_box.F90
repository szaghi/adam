!< Unit test of the PRISM-CPU AMR_GEO primitive-box marker (`mark_by_geometry`).
program test_amr_geo_box
!< Unit test of `prism_cpu_object%mark_by_geometry` — the AMR_GEO + AMR_GEO_PRIMITIVE_BOX
!< marker that flags blocks for refinement by an axis-aligned box, deterministically
!< and independently of the solution (issue #13 §7.5, milestone M0).
!<
!< **Why this test exists.** The only working AMR marker today is the total-variation
!< marker (`mark_by_j_vec_total_variation`), which is solution-dependent and therefore a
!< poor fixture: the set of refined blocks depends on the field state. M0 adds a
!< geometric box marker so a test can demand a *deterministic* 2:1 coarse-fine jump at
!< known coordinates — the fixture the whole intra-realm AMR seam reflux path (M1-M4)
!< rests on. This test pins the marking rule directly, with no INI, no solve, no IO.
!<
!< **What it pins** (verbatim from `mark_by_geometry`):
!<```
!< refinements_needed(b) = TO_BE_REFINED   iff  centroid(b) in [box_emin, box_emax]
!<                                              AND  tree%level(code(b)) < target_level
!<```
!< plus: blocks outside the box are not refined; blocks inside the box but already at or
!< above `target_level` are not refined. The marker is *additive* — its `do_init` default
!< seeds the whole query to `TO_NOT_TOUCH`, so both an out-of-box block and an
!< in-box-but-deep block read `TO_NOT_TOUCH` (it does not force coarsening outside the box,
!< unlike the total-variation marker).
!<
!< The realm's `field` (emin/emax/code/blocks_number) and `tree` members are populated by
!< hand to the minimum the TBP reads.

use adam_prism_cpu_object, only : prism_cpu_object
use adam_globals, only : mpih
use adam_parameters, only : TO_BE_REFINED, TO_NOT_TOUCH
use penf, only : I4P, I8P, R8P

implicit none

integer(I4P), parameter :: nb           = 4_I4P                       !< Blocks.
integer(I4P), parameter :: target_level = 2_I4P                       !< Refine blocks below level 2.
real(R8P),    parameter :: box_emin(3)  = [0.0_R8P, 0.0_R8P, 0.0_R8P] !< Box minimum corner.
real(R8P),    parameter :: box_emax(3)  = [1.0_R8P, 1.0_R8P, 1.0_R8P] !< Box maximum corner.

type(prism_cpu_object) :: realm          !< Realm under test.
integer(I4P)           :: expected(1:nb) !< Expected refinements_needed per block.
integer(I4P)           :: b              !< Counter.
integer(I4P), target   :: nb_target = nb !< Pointer target for realm%blocks_number (wired by realm in production).
logical                :: test_passed    !< Aggregate pass flag.

test_passed = .true.
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)

! Build the realm's minimal field/tree state by hand.
realm%adam%tree%ratio = 8_I4P       ! 3D octree refinement ratio (level() reads it)
realm%adam%field%blocks_number = nb
realm%blocks_number => nb_target    ! production wires this to field%blocks_number in adam_realm_object

allocate(realm%adam%field%emin(1:3, 1:nb))
allocate(realm%adam%field%emax(1:3, 1:nb))
allocate(realm%adam%field%code(1:nb))
allocate(realm%adam%field%refinements_needed(1:nb))

! Block 1: centroid (0.5,0.5,0.5) IN box, code 0 -> level 1 (< 2)  => TO_BE_REFINED
realm%adam%field%emin(:,1) = [0.0_R8P, 0.0_R8P, 0.0_R8P]
realm%adam%field%emax(:,1) = [1.0_R8P, 1.0_R8P, 1.0_R8P]
realm%adam%field%code(1)   = 0_I8P

! Block 2: centroid (2.5,0.5,0.5) OUTSIDE box (x>1), code 0 -> level 1  => TO_NOT_TOUCH
! (a refine-region marker is additive: outside the box it leaves blocks untouched, it does
! not force coarsening — that is the TV marker's policy, not this one's).
realm%adam%field%emin(:,2) = [2.0_R8P, 0.0_R8P, 0.0_R8P]
realm%adam%field%emax(:,2) = [3.0_R8P, 1.0_R8P, 1.0_R8P]
realm%adam%field%code(2)   = 0_I8P

! Block 3: centroid (0.25,0.25,0.25) IN box, code 8 -> level 2 (NOT < 2) => TO_NOT_TOUCH
realm%adam%field%emin(:,3) = [0.0_R8P, 0.0_R8P, 0.0_R8P]
realm%adam%field%emax(:,3) = [0.5_R8P, 0.5_R8P, 0.5_R8P]
realm%adam%field%code(3)   = 8_I8P

! Block 4: centroid (0.75,0.75,0.25) IN box, code 1 -> level 1 (< 2)  => TO_BE_REFINED
realm%adam%field%emin(:,4) = [0.5_R8P, 0.5_R8P, 0.0_R8P]
realm%adam%field%emax(:,4) = [1.0_R8P, 1.0_R8P, 0.5_R8P]
realm%adam%field%code(4)   = 1_I8P

expected = [TO_BE_REFINED, TO_NOT_TOUCH, TO_NOT_TOUCH, TO_BE_REFINED]

! Sanity: the level construction is what the test assumes.
block
   logical :: levels_ok
   levels_ok = (realm%adam%tree%level(realm%adam%field%code(1)) == 1_I4P) .and. &
               (realm%adam%tree%level(realm%adam%field%code(3)) == 2_I4P)
   if (levels_ok) then
      write(*,'(A)') 'PASS: tree%level() of the seed codes is as assumed (1 and 2)'
   else
      write(*,'(A)') 'FAIL: tree%level() of the seed codes is NOT as assumed — fixture invalid'
      write(*,'(A,I0,A,I0)') '      level(code1)=', realm%adam%tree%level(realm%adam%field%code(1)), &
                             '  level(code3)=', realm%adam%tree%level(realm%adam%field%code(3))
      test_passed = .false.
   endif
endblock

! Apply the marker.
call realm%mark_by_geometry(box_emin=box_emin, box_emax=box_emax, target_level=target_level)

! Assertion: every block matches the expected refinement flag exactly.
block
   logical :: marks_ok
   marks_ok = .true.
   do b = 1_I4P, nb
      if (realm%adam%field%refinements_needed(b) /= expected(b)) then
         marks_ok = .false.
         write(*,'(A,I0,A,I0,A,I0)') 'FAIL: block ', b, ' flag=', realm%adam%field%refinements_needed(b), &
                                     ' expected=', expected(b)
      endif
   enddo
   if (marks_ok) then
      write(*,'(A)') 'PASS: in-box blocks below target_level refine; out-of-box / deep blocks do not'
   else
      test_passed = .false.
   endif
endblock

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: prism amr geo box'
else
   write(*,'(A)') 'TEST FAILED: prism amr geo box'
endif

call mpih%finalize
if (.not. test_passed) error stop 1
endprogram test_amr_geo_box
