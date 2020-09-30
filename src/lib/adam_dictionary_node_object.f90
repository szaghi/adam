!< ADAM, dictionary node class definition.
module adam_dictionary_node_object
!< ADAM, dictionary node class definition.

use adam_tree_topology
! use MORTIF, only : morton3D, morton2D
use PENF, only : I4P, I8P, str, cton

implicit none
private
public :: KEY_LEN
public :: destroy_dictionary_node
public :: key_str, key_int
public :: dictionary_node_object

integer(I4P), parameter :: KEY_LEN = 49 !< Length of dictionary node's key.

type :: dictionary_node_object
   !< Dictionary node class definition.
   character(len=KEY_LEN),                public :: key=''           !< The key.
   integer(I8P),                          public :: code=-1_I8P      !< The Morton code.
   integer(I8P),                          public :: content=0_I8P    !< The content.
   type(dictionary_node_object), pointer, public :: next=>null()     !< The next node in the dictionary.
   type(dictionary_node_object), pointer, public :: previous=>null() !< The previous node in the dictionary.
   contains
      ! public methods
      procedure, pass(self) :: compute_morton_code !< Compute Morton code.
      procedure, pass(self) :: destroy             !< Destroy dictionary node.
      procedure, pass(self) :: initialize          !< Initialize dictionary node with key/content pair.
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

   pure function key_str(l, tijk, bijk) result(key)
   !< Return key in string format.
   integer(I4P), intent(in) :: l       !< Refinementl level.
   integer(I4P), intent(in) :: tijk(3) !< Tree coordinates into the forest.
   integer(I4P), intent(in) :: bijk(3) !< Block coordinates into the tree.
   character(len=KEY_LEN)   :: key     !< The key.

   key = trim(str('(I7)',l))//trim(str('(I7)',tijk(1)))//trim(str('(I7)',tijk(2)))//trim(str('(I7)',tijk(3)))//&
                              trim(str('(I7)',bijk(1)))//trim(str('(I7)',bijk(2)))//trim(str('(I7)',bijk(3)))
   endfunction key_str

   subroutine key_int(key, l, tijk, bijk)
   !< Return key in string format.
   character(len=*), intent(in)            :: key     !< The key.
   integer(I4P),     intent(out), optional :: l       !< Refinementl level.
   integer(I4P),     intent(out), optional :: tijk(3) !< Tree coordinates into the forest.
   integer(I4P),     intent(out), optional :: bijk(3) !< Block coordinates into the tree.

   if (present(l)) l = cton(key(1:7), 1_I4P)
   if (present(tijk)) then
      tijk(1) = cton(key(8 :14), 1_I4P)
      tijk(2) = cton(key(15:21), 1_I4P)
      tijk(3) = cton(key(22:28), 1_I4P)
   endif
   if (present(bijk)) then
      bijk(1) = cton(key(29:35), 1_I4P)
      bijk(2) = cton(key(36:42), 1_I4P)
      bijk(3) = cton(key(43:49), 1_I4P)
   endif
   endsubroutine key_int

   ! public methods
   subroutine compute_morton_code(self, ratio)
   !< Compute Morton code.
   class(dictionary_node_object), intent(inout)        :: self    !< Dictionary node.
   integer(I4P),                  intent(in), optional :: ratio   !< Refinement ratio.
   integer(I4P)                                        :: ijkl(4) !< IJKL coordinates.
   integer(I4P)                                        :: ratio_  !< Refinement ratio, local variable.

   ratio_ = 8_I4P ; if (present(ratio)) ratio_ = ratio
   call key_int(key=self%key, l=ijkl(4), bijk=ijkl(1:3))
   if (   ratio_==8_I8P) then
      self%code = coordinates_to_morton(i=ijkl(1), j=ijkl(2), k=ijkl(3), l=ijkl(4), ratio=ratio_)
   elseif(ratio_==4_I8P) then
      self%code = coordinates_to_morton(i=ijkl(1), j=ijkl(2), l=ijkl(4), ratio=ratio_)
   endif
   endsubroutine compute_morton_code

   elemental subroutine destroy(self)
   !< Destroy dictionary node.
   class(dictionary_node_object), intent(inout) :: self  !< Dictionary node.
   type(dictionary_node_object)                 :: fresh !< Fresh instance of dictionary node.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, content, key, code, ratio)
   !< Initialize dictionary node with key/content pair.
   class(dictionary_node_object), intent(inout)        :: self    !< Dictionary node.
   integer(I8P),                  intent(in)           :: content !< The content.
   character(len=*),              intent(in), optional :: key     !< The key.
   integer(I8P),                  intent(in), optional :: code    !< The Morton code.
   integer(I4P),                  intent(in), optional :: ratio   !< Refinement ratio.

   self%content = content
   if (present(key)) then
      self%key = key
      call self%compute_morton_code(ratio=ratio)
   elseif (present(code)) then
      ! @TODO to be implemented
   endif
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine dictionary_node_assign_dictionary_node(lhs, rhs)
   !< Operator `=`.
   class(dictionary_node_object), intent(inout) :: lhs !< Left hand side.
   type(dictionary_node_object),  intent(in)    :: rhs !< Right hand side.

   lhs%key = rhs%key
   lhs%code = rhs%code
   lhs%content = rhs%content
   endsubroutine dictionary_node_assign_dictionary_node
endmodule adam_dictionary_node_object
