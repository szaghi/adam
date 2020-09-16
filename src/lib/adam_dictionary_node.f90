!< ADAM, dictionary node class definition.
module adam_dictionary_node_object
!< ADAM, dictionary node class definition.

use PENF, only : I4P

implicit none
private
public :: destroy_dictionary_node
public :: dictionary_node_object

integer(I4P), parameter :: KEY_LEN = 50 !< Length of dictionary node's key.

type :: dictionary_node_object
!< Dictionary node class definition.
   character(len=KEY_LEN),                public :: key=''           !< The key.
   integer(I4P),                          public :: content=0_I4P    !< The content.
   type(dictionary_node_object), pointer, public :: next=>null()     !< The next node in the dictionary.
   type(dictionary_node_object), pointer, public :: previous=>null() !< The previous node in the dictionary.
   contains
      ! public methods
      procedure, pass(self) :: destroy !< Destroy dictionary node.
      ! operators
      generic :: assignment(=) => dictionary_node_assign_dictionary_node      !< Overload `=`.
      procedure, pass(lhs), private :: dictionary_node_assign_dictionary_node !< Operator `=`.
endtype dictionary_node_object

contains
  ! public non TBP
  recursive subroutine destroy_dictionary_node(node)
  !< Destroy dictionary node and its subsequent ones.
  type(dictionary_node_object), pointer, intent(inout) :: node !< The node.

  if (associated(node)) then
    call node%destroy
    call destroy_dictionary_node(node=node%next)
    node%previous => null()
    deallocate(node)
    node => null()
  endif
  endsubroutine destroy_dictionary_node

   ! public methods
   elemental subroutine destroy(self)
   !< Destroy dictionary node.
   class(dictionary_node_object), intent(inout) :: self  !< Dictionary node.
   type(dictionary_node_object)                 :: fresh !< Fresh instance of dictionary node.

   self = fresh
   endsubroutine destroy

   ! operators
   ! =
   pure subroutine dictionary_node_assign_dictionary_node(lhs, rhs)
   !< Operator `=`.
   class(dictionary_node_object), intent(inout) :: lhs !< Left hand side.
   type(dictionary_node_object),  intent(in)    :: rhs !< Right hand side.

   lhs%key = rhs%key
   lhs%content = rhs%content
   endsubroutine dictionary_node_assign_dictionary_node
endmodule adam_dictionary_node_object
