!< Unit test of the block partition over the ranks: `contiguous_min_max_split` of `adam_tree_object` (issue #123).
program test_tree_balance
!< Unit test of the block partition over the ranks: `contiguous_min_max_split` of `adam_tree_object` (issue #123).
!<
!< **Why this test exists.** `tree%mpi_redistribute` cuts the Morton-sorted blocks into one contiguous run per rank, a complete
!< set of sibling leaves (the `ratio` children of one parent, which derefinement reads on one rank) being a unit the cut never
!< splits. The greedy walk it replaces put the 120 blocks of S4c (#121, a 2:1 box) 20/13/19/17/13/19/19/0 on 8 ranks; the split of
!< the units is now `contiguous_min_max_split`, a pure function this test checks without a tree or MPI ranks.
!<
!< **What it pins.**
!< - The S4c unit sequence (36 units: 8 sets of 8 refined children and 28 single blocks, in the Morton order of the tree): the
!<   busiest rank holds 30, 18 and 10 blocks on 4, 8 and 16 ranks (the optimum keeping the sets whole, issue #123), no rank empty.
!< - A uniform tree (8 sets of 8): 8 per rank on 8 ranks; on 16 ranks 8 ranks hold 8 and 8 are empty (fewer units than ranks).
!< - Every split is contiguous (the part is non-decreasing along the units, by steps of 0 or 1), uses only parts 0 to P-1, and
!<   leaves no part empty when there are at least as many units as parts.
!< - Optimality: on 2000 random sequences (1 to 40 units of weight 1 or 8, 1 to 12 parts) the largest part equals the exact
!<   minimum of a dynamic program over all contiguous splits.
!<
!< Usage: `mpirun -np 1 exe/test_tree_balance` (the library module links MPI; the test calls none).

use :: adam_tree_object, only : contiguous_min_max_split
use :: penf,             only : I4P, I8P, R8P, str

implicit none

integer(I8P), parameter   :: S4C(36) = [1_I8P, 8_I8P, 1_I8P, 1_I8P, 1_I8P, 8_I8P, 1_I8P, 1_I8P, 8_I8P, 1_I8P, 1_I8P, 1_I8P, &
                                        8_I8P, 1_I8P, 1_I8P, 1_I8P, 8_I8P, 8_I8P, 1_I8P, 8_I8P, 1_I8P, 1_I8P, 1_I8P, 8_I8P, &
                                        1_I8P, 1_I8P, 8_I8P, 1_I8P, 1_I8P, 1_I8P, 8_I8P, 1_I8P, 1_I8P, 1_I8P, 8_I8P, 8_I8P] !< S4c units.
integer(I8P)              :: uniform(8)          !< Uniform tree units.
integer(I8P), allocatable :: w(:)                !< Random units.
integer(I4P)              :: parts               !< Number of parts.
integer(I4P)              :: trial               !< Random trial counter.
integer(I4P)              :: n                   !< Number of random units.
integer(I4P)              :: u                   !< Counter.
integer(I4P)              :: failures            !< Failed checks.
real(R8P)                 :: r                   !< Random number.
integer                   :: seed_size           !< Size of the random seed.
integer,      allocatable :: seed(:)             !< Random seed.

failures = 0_I4P
call check_case('S4c, 4 ranks', S4C, 4_I4P, 30_I8P, .false.)
call check_case('S4c, 8 ranks', S4C, 8_I4P, 18_I8P, .false.)
call check_case('S4c, 16 ranks', S4C, 16_I4P, 10_I8P, .false.)
uniform = 8_I8P
call check_case('uniform, 8 ranks', uniform, 8_I4P, 8_I8P, .false.)
call check_case('uniform, 16 ranks', uniform, 16_I4P, 8_I8P, .true.)
call check_case('S4c, 1 rank', S4C, 1_I4P, 120_I8P, .false.)
call random_seed(size=seed_size)
allocate(seed(seed_size))
seed = 20261110
call random_seed(put=seed)
do trial=1, 2000
   call random_number(r) ; n = 1_I4P + int(r * 40._R8P, I4P)
   call random_number(r) ; parts = 1_I4P + int(r * 12._R8P, I4P)
   if (allocated(w)) deallocate(w)
   allocate(w(n))
   do u=1, n
      call random_number(r)
      w(u) = merge(8_I8P, 1_I8P, r < 0.3_R8P)
   enddo
   call check_case('random '//trim(str(trial, .true.)), w, parts, optimum(w, parts), n < parts, quiet=.true.)
enddo
if (failures == 0_I4P) then
   print '(A)', 'test_tree_balance: PASSED'
else
   print '(A)', 'test_tree_balance: FAILED ('//trim(str(failures, .true.))//' checks)'
   error stop 1
endif

contains
   subroutine check_case(label, weights, parts, expected_max, may_be_empty, quiet)
   !< Split `weights` into `parts` and check contiguity, the part range, no empty part (unless allowed) and the largest part.
   character(*), intent(in)           :: label        !< Case label.
   integer(I8P), intent(in)           :: weights(:)   !< Units.
   integer(I4P), intent(in)           :: parts        !< Number of parts.
   integer(I8P), intent(in)           :: expected_max !< Expected largest part.
   logical,      intent(in)           :: may_be_empty !< Empty parts are allowed (fewer units than parts).
   logical,      intent(in), optional :: quiet        !< Print only a failure.
   integer(I4P)                       :: part(size(weights)) !< Part of each unit.
   integer(I8P)                       :: load(0:parts-1)     !< Load of each part.
   logical                            :: ok                  !< Case result.
   integer(I4P)                       :: v                   !< Counter.
   character(:), allocatable          :: loads               !< Loads, for the report.

   part = contiguous_min_max_split(weights=weights, parts=parts)
   ok = all(part >= 0_I4P .and. part < parts)
   if (ok) ok = part(1) == 0_I4P
   do v=2, size(weights)
      if (part(v) - part(v-1) < 0_I4P .or. part(v) - part(v-1) > 1_I4P) ok = .false.
   enddo
   load = 0_I8P
   if (ok) then
      do v=1, size(weights)
         load(part(v)) = load(part(v)) + weights(v)
      enddo
   endif
   if (.not.may_be_empty .and. any(load == 0_I8P)) ok = .false.
   if (maxval(load) /= expected_max) ok = .false.
   if (.not.ok) failures = failures + 1_I4P
   if (.not.ok .or. .not.present(quiet)) then
      loads = ''
      do v=0, parts - 1
         loads = loads//' '//trim(str(load(v), .true.))
      enddo
      print '(A)', '  '//label//': loads'//loads//', largest '//trim(str(maxval(load), .true.))//' (expected '// &
                   trim(str(expected_max, .true.))//')'//merge('  ok  ', '  FAIL', ok)
   endif
   endsubroutine check_case

   function optimum(weights, parts) result(best)
   !< Return the minimum over all contiguous splits of `weights` into at most `parts` parts of the largest part (dynamic program).
   integer(I8P), intent(in) :: weights(:) !< Units.
   integer(I4P), intent(in) :: parts      !< Number of parts.
   integer(I8P)             :: best       !< Minimum largest part.
   integer(I8P)             :: f(0:size(weights), 0:parts) !< f(i, k): best for the first i units in k parts.
   integer(I8P)             :: prefix(0:size(weights))     !< Prefix sums.
   integer(I4P)             :: i, j, k                     !< Counters.

   prefix(0) = 0_I8P
   do i=1, size(weights)
      prefix(i) = prefix(i-1) + weights(i)
   enddo
   f = huge(1_I8P)
   f(0,:) = 0_I8P
   do k=1, parts
      do i=1, size(weights)
         do j=0, i - 1
            if (f(j,k-1) < huge(1_I8P)) f(i,k) = min(f(i,k), max(f(j,k-1), prefix(i) - prefix(j)))
         enddo
      enddo
   enddo
   best = minval(f(size(weights),1:parts))
   endfunction optimum
endprogram test_tree_balance
