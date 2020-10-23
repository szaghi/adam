!< ADAM, tree node class definition.
module adam_tree_node_object
!< ADAM, tree node class definition.

use PENF, only : I4P, I8P!, str, cton

implicit none
private
public :: destroy_tree_node
public :: tree_node_object
public :: NODE_TO_BE_REFINED, NODE_TO_BE_DEREFINED, NODE_TO_NOT_TOUCH

! integer(I4P), parameter :: KEY_LEN = 49 !< Length of tree node's key.
integer(I4P), parameter :: NODE_TO_BE_REFINED=1_I4P    !< Flag for node to be refined.
integer(I4P), parameter :: NODE_TO_BE_DEREFINED=-1_I4P !< Flag for node to be derefined.
integer(I4P), parameter :: NODE_TO_NOT_TOUCH=0_I4P     !< Flag for node to be untouched.

type :: tree_node_object
   !< Tree node class definition.
   integer(I8P),                    public :: code=-2_I8P             !< The Morton code.
   integer(I4P),                    public :: refinement_needed=0_I4P !< Flag for refinement/derefinement algorithm.
   integer(I4P),                    public :: myrank=0_I4P            !< MPI rank process.
   integer(I8P),                    public :: block_index=1_I8P       !< Block index in the field array.
   type(tree_node_object), pointer, public :: next=>null()            !< The next node in the tree.
   type(tree_node_object), pointer, public :: previous=>null()        !< The previous node in the tree.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy tree node.
      procedure, pass(self) :: initialize !< Initialize tree node.
      ! operators
      generic :: assignment(=) => tree_node_assign_tree_node      !< Overload `=`.
      procedure, pass(lhs), private :: tree_node_assign_tree_node !< Operator `=`.
endtype tree_node_object

contains
   ! public non TBP
   recursive subroutine destroy_tree_node(node)
   !< Destroy tree node and its subsequent ones.
   type(tree_node_object), pointer, intent(inout) :: node !< The node.

   if (associated(node)) then
      call node%destroy
      call destroy_tree_node(node=node%next)
      node%previous => null()
      deallocate(node)
      node => null()
   endif
   endsubroutine destroy_tree_node

   ! public methods
   elemental subroutine destroy(self)
   !< Destroy tree node.
   class(tree_node_object), intent(inout) :: self  !< Tree node.
   type(tree_node_object)                 :: fresh !< Fresh instance of tree node.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, code, refinement_needed, myrank, block_index)
   !< Initialize tree node.
   class(tree_node_object), intent(inout)        :: self              !< Tree node.
   integer(I8P),            intent(in)           :: code              !< The Morton code.
   integer(I4P),            intent(in), optional :: refinement_needed !< Flag for refinement/derefinement algorithm.
   integer(I4P),            intent(in), optional :: myrank            !< MPI rank process.
   integer(I8P),            intent(in), optional :: block_index       !< Block index in the field array.

   call self%destroy
   self%code = code
   if (present(refinement_needed)) self%refinement_needed = refinement_needed
   if (present(myrank           )) self%myrank            = myrank
   if (present(block_index      )) self%block_index       = block_index
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine tree_node_assign_tree_node(lhs, rhs)
   !< Operator `=`.
   class(tree_node_object), intent(inout) :: lhs !< Left hand side.
   type(tree_node_object),  intent(in)    :: rhs !< Right hand side.

   lhs%code = rhs%code
   lhs%refinement_needed = rhs%refinement_needed
   lhs%myrank = rhs%myrank
   lhs%block_index = rhs%block_index
   endsubroutine tree_node_assign_tree_node
endmodule adam_tree_node_object
