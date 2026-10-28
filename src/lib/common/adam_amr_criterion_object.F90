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
use :: adam_realm_object, only : realm_object
use :: penf

implicit none
private
public :: amr_criterion_object
public :: amr_criterion_pointer

type, abstract :: amr_criterion_object
   !< A refinement criterion: one vote per block.
   character(:), allocatable :: name !< The name an input selects it by (`[amr_marker_N] name`).
   contains
      procedure(vote_interface), pass(self), deferred :: vote !< The vote of one block.
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
endmodule adam_amr_criterion_object
