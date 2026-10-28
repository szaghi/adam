!< A refinement criterion of one's own, for RIPPLE (issue #110, P4.5): the worked example of the build-your-own tutorial.
module ripple_crest_criterion
!< A refinement criterion of one's own, for RIPPLE (issue #110, P4.5): the worked example of the build-your-own tutorial.
!<
!< `crest_criterion` refines the blocks where the advected profile rises above `threshold`, until they reach
!< `target_level`, and leaves the others alone. It reads the app's state, so its vote selects the app's type.
use :: adam_amr_criterion_object, only : amr_criterion_object
use :: adam_parameters,           only : TO_BE_REFINED, TO_NOT_TOUCH
use :: adam_realm_object,         only : realm_object
use :: adam_ripple_common_object, only : ripple_common_object
use :: penf,                      only : I4P, R8P

implicit none
private
public :: crest_criterion

type, extends(amr_criterion_object) :: crest_criterion
   !< Refine where u > threshold.
   real(R8P)    :: threshold=0.9_R8P  !< Refine the blocks where some cell holds u above it.
   integer(I4P) :: target_level=3_I4P !< Refine up to this level.
   contains
      procedure, pass(self) :: vote
endtype crest_criterion

contains
   function vote(self, realm, b) result(mark)
   !< Refine block `b` when one of its cells holds u above the threshold and the block is below the target level;
   !< otherwise leave it as it is.
   class(crest_criterion), intent(inout) :: self  !< The criterion.
   class(realm_object),    intent(in)    :: realm !< The realm.
   integer(I4P),           intent(in)    :: b     !< Local block index.
   integer(I4P)                          :: mark  !< The vote.

   mark = TO_NOT_TOUCH
   select type(realm)
   class is(ripple_common_object)
      if (realm%adam%tree%level(realm%adam%field%code(b)) >= self%target_level) return
      if (any(realm%q(1, 1:realm%ni, 1:realm%nj, 1:realm%nk, b) > self%threshold)) mark = TO_BE_REFINED
   endselect
   endfunction vote
endmodule ripple_crest_criterion

!< RIPPLE CPU driver with a registered refinement criterion (issue #110, P4.5).
program test_ripple_crest_criterion
!< RIPPLE CPU driver with a registered refinement criterion (issue #110, P4.5).
!<
!< The RIPPLE CPU driver (adam_ripple_cpu.F90) plus the registration of `crest_criterion`: an input selects it by
!< `[amr_marker_N] mode = 6, name = crest` (src/tests/ripple/crest.ini). RIPPLE refines at initialization only, so the
!< criterion acts there (`runtime_only = .false.`).
use :: adam_forest_object,     only : forest_object
use :: adam_realm_object,      only : realm_pointer
use :: adam_ripple_cpu_object, only : ripple_cpu_object
use :: ripple_crest_criterion, only : crest_criterion

implicit none

type(ripple_cpu_object), target :: realm(1)          !< The realm.
type(realm_pointer)             :: realm_pointers(1) !< A pointer to it: what the forest takes.
type(forest_object)             :: forest            !< Orchestrator.
type(crest_criterion),   target :: crest             !< The criterion (the program owns it).
character(999)                  :: input_file_name   !< Input file name.

input_file_name = 'input.ini'
if (command_argument_count() > 0) call get_command_argument(1, input_file_name)
crest%name = 'crest'
crest%runtime_only = .false. ! RIPPLE refines at initialization only
call realm(1)%register_amr_criterion(crest)
realm_pointers(1)%ptr => realm(1)
call forest%simulate(realm=realm_pointers, filename=trim(adjustl(input_file_name)))
endprogram test_ripple_crest_criterion
