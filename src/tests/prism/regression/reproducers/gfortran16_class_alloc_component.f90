! Reproducer: gfortran-16-experimental emits a FALSE `-fcheck=bounds` failure for
! a two-levels-deep allocatable array `self%face(i)%a` accessed through a
! `class()` `pass(self)` TBP. The computation itself is correct; only the
! bounds-check INSTRUMENTATION is wrong.
!
! WHY THIS FILE EXISTS
! --------------------
! `test_prism_reflux_apply` SIGSEGVs/aborts at runtime ONLY in the `-gnu-debug`
! build mode (which uses `-fcheck=all`) when built with the gfortran-16
! experimental trunk this box's `mpif90` wraps:
!
!     Fortran runtime error: Index '1' of dimension 1 of array 'self%face'
!     above upper bound of 0
!
! The crash is NOT in ADAM's code, and it is NOT a real descriptor corruption.
! The same test built in the OPTIMIZED `-gnu` mode (no `-fcheck`) PASSES all
! assertions, and this reproducer runs to SURVIVED without `-fcheck=bounds`.
!
! EXACTLY WHAT HAPPENS (statements split so the fault is unambiguous; a one-liner
! `allocate(...); ... = 0.0` hides which sub-statement faults)
! --------------------------------------------------------------------------
!   A  allocate(self%face(1:n))            -> SUCCEEDS; ubound(self%face,1)=1 reads correctly
!   B  allocate(self%face(1)%a(3))         -> SUCCEEDS
!   X  ... ubound(self%face(1)%a, 1) ...   -> with -fcheck=bounds: ABORTS here (false positive)
!                                             without it: reads 3 correctly, runs to SURVIVED
!
! So the allocations and the actual value are fine. The first expression that
! DOUBLE-INDEXES the component (`self%face(1)%a`, two levels deep through the
! `class` passed-object) makes the gfortran-16 bounds-check instrumentation
! re-fetch the OUTER `self%face` descriptor as garbage (ubound 0) and trip a
! spurious out-of-bounds. A single-level `ubound(self%face,1)` (statement A's
! print) is instrumented correctly. In the original test this surfaced at
! `self%face(1)%a = 0.0` (an array assignment whose bounds check reads the
! nested-array descriptor), which is why the backtrace pointed one statement
! AFTER the allocation.
!
! TRIGGER ISOLATION
! -----------------
!   -fcheck=all       -> CRASH
!   -fcheck=bounds    -> CRASH        (this is the specific culprit)
!   -finit-derived    -> SURVIVED     (not involved)
!   (no -fcheck)      -> SURVIVED     (incl. -O2)
!
! BISECTION
! ---------
!   gfortran 13/14/15 (stable), even with -fcheck=all : SURVIVED (X reads ubound = 3)
!   gfortran 16.0.1 20260315 (experimental) [trunk r16-8100-g3aca3bae8ee]
!                                       with -fcheck=bounds : false abort at X
!
! It is specifically the POLYMORPHIC (`class`) passed-object: replacing `init`
! with a module subroutine taking a `type(reg_t)` dummy makes the instrumented
! access succeed on gfortran-16-experimental too. (flux_register is never
! extended, so the `class` is incidental — but TBPs require it.)
!
! BUILD / RUN
! -----------
!   <gf16-mpif90> -O0 -g                 ... && ./repro   # SURVIVED  (no bounds check)
!   <gf16-mpif90> -O0 -g -fcheck=bounds  ... && ./repro   # false abort at the X line
!   gfortran-15   -O0 -g -fcheck=all     ... && ./repro   # SURVIVED  (stable compiler)
!
! IMPACT / WORKAROUND for ADAM
! ----------------------------
! The solver and the OPTIMIZED tests are unaffected (correct results). Only the
! `-gnu-debug` test mode false-aborts under gfortran-16-experimental. Run the
! debug-mode CPU tests with a STABLE gfortran (13/14/15), or run those tests in
! the optimized `-gnu` mode, until the compiler is fixed. File upstream against
! gcc Fortran if untracked.

module m
   implicit none
   type :: face_t
      real, allocatable :: a(:)            !< element type has an allocatable component
   end type
   type :: reg_t
      integer :: n = 0
      type(face_t), allocatable :: face(:) !< allocatable array of the above
   contains
      procedure, pass(self) :: init        !< class() pass(self) — the trigger
   end type
contains
   subroutine init(self, n)
      class(reg_t), intent(inout) :: self
      integer, intent(in) :: n

      self%n = n

      ! A — allocate the outer array. SUCCEEDS; single-level bound reads correctly.
      allocate(self%face(1:n))
      print *, 'A ok: ubound(self%face,1) =', ubound(self%face, 1)

      ! B — allocate the nested array. SUCCEEDS.
      allocate(self%face(1)%a(3))
      print *, 'B ok: allocated(self%face(1)%a) =', allocated(self%face(1)%a)

      ! X — first DOUBLE-INDEX read of the nested array. With -fcheck=bounds on
      ! gfortran-16-experimental this false-aborts ("self%face above upper bound
      ! of 0"); without it (or on stable gfortran) it reads 3 and continues.
      print *, 'X    : ubound(self%face(1)%a,1) =', ubound(self%face(1)%a, 1)

      self%face(1)%a = 0.0
      print *, 'assignment ok'
   end subroutine
end module

program p
   use m
   type(reg_t) :: r
   call r%init(1)
   print *, 'SURVIVED'
end program
