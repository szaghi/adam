!< ADAM AMR criterion: a refinement criterion an application or a user defines (issue #110, P4).
module adam_amr_criterion_object
!< ADAM AMR criterion: a refinement criterion an application or a user defines (issue #110, P4).
!<
!< The built-in markers of `adam_amr_object` (geometry, gradient, Loehner, Jeans) are selected by a numeric mode in the
!< `[amr_marker_N]` sections. A criterion of this kind is code: a type extending `amr_criterion_object` whose `vote`
!< returns, for one block of a realm, `TO_BE_REFINED`, `TO_NOT_TOUCH` or `TO_BE_DEREFINED`. The program registers an
!< instance with its realm before the simulation starts (the app's `register_amr_criterion`), and an input selects it
!< by name, `[amr_marker_N] mode = 6` (`AMR_USER`) and `name = <its name>`. The app combines its votes with those of the
!< other markers as it combines theirs. `vote` receives the realm as `class(realm_object)`: a criterion selects the
!< app's concrete type (`select type`) to read its fields. It is called on the host, for the blocks this rank owns.
!<
!< `mark` gives the votes of all the blocks of the rank at once; its default asks `vote` block by block. A criterion that
!< is not local to a block (it compares neighbours, dilates its marks, refreshes ghosts first) overrides it: the
!< built-in markers of an app are criteria of this kind, so built-in and registered criteria share one list. `mark` is
!< called with `initial = .true.` by the AMR of the initial condition, unless `runtime_only` (the default: a
!< registered criterion acts at the runtime regrids only).
use :: adam_realm_object, only : realm_object
use :: penf

implicit none
private
public :: amr_criterion_object
public :: amr_criterion_pointer

type, abstract :: amr_criterion_object
   !< A refinement criterion: one vote per block.
   character(:), allocatable :: name                !< The name an input selects it by (`[amr_marker_N] name`).
   logical                   :: runtime_only=.true. !< Skipped by the AMR of the initial condition.
   contains
      procedure(vote_interface), pass(self), deferred :: vote !< The vote of one block.
      procedure,                 pass(self)           :: mark !< The votes of every block of the rank.
endtype amr_criterion_object

type :: amr_criterion_pointer
   !< A registered criterion (the program owns the instance, which must outlive the simulation).
   class(amr_criterion_object), pointer :: ptr => null() !< The criterion.
endtype amr_criterion_pointer

abstract interface
   function vote_interface(self, realm, b) result(mark)
   !< The vote of block `b` (a local block index of this rank) of `realm`: TO_BE_REFINED, TO_NOT_TOUCH or
   !< TO_BE_DEREFINED (adam_parameters).
   import :: amr_criterion_object, realm_object, I4P
   class(amr_criterion_object), intent(inout) :: self  !< The criterion.
   class(realm_object),         intent(in)    :: realm !< The realm the block belongs to.
   integer(I4P),                intent(in)    :: b     !< Local block index.
   integer(I4P)                               :: mark  !< The vote.
   endfunction vote_interface
endinterface

contains
   subroutine mark(self, realm, votes, initial)
   !< The votes of every block of this rank (`votes(b)` for local block `b`), by default `vote` block by block. An
   !< override may read and update the realm (ghost refreshes, auxiliary fields), and may act differently for the AMR of
   !< the initial condition (`initial`).
   class(amr_criterion_object), intent(inout) :: self     !< The criterion.
   class(realm_object),         intent(inout) :: realm    !< The realm.
   integer(I4P),                intent(out)   :: votes(:) !< One vote per local block.
   logical,                     intent(in)    :: initial  !< The AMR of the initial condition (the default votes alike).
   integer(I4P)                               :: b        !< Counter.

   do b=1, size(votes, kind=I4P)
      votes(b) = self%vote(realm=realm, b=b)
   enddo
   endsubroutine mark
endmodule adam_amr_criterion_object
