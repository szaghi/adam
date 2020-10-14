!< ADAM, tree node class definition.
module adam_tree_node_object
!< ADAM, tree node class definition.

! use adam_tree_topology
! use MORTIF, only : morton3D, morton2D
use PENF, only : I4P, I8P!, str, cton

implicit none
private
public :: destroy_tree_node
public :: tree_node_object
public :: TO_BE_REFINED, TO_BE_DEREFINED, TO_NOT_TOUCH

! integer(I4P), parameter :: KEY_LEN = 49 !< Length of tree node's key.
integer(I4P), parameter :: TO_BE_REFINED=1_I4P    !< Flag for node to be refined.
integer(I4P), parameter :: TO_BE_DEREFINED=-1_I4P !< Flag for node to be derefined.
integer(I4P), parameter :: TO_NOT_TOUCH=0_I4P     !< Flag for node to be untouched.

type :: tree_node_object
   !< Tree node class definition.
   ! character(len=KEY_LEN),          public :: key=''                  !< The key.
   integer(I8P),                    public :: code=-2_I8P             !< The Morton code.
   integer(I8P),                    public :: finest_code=0_I8P       !< The Morton code.
   integer(I8P),                    public :: content=0_I8P           !< The content.
   integer(I4P),                    public :: refinement_needed=0_I4P !< Flag for refinement/derefinement algorithm.
   integer(I4P),                    public :: myrank=0_I4P            !< MPI rank process.
   integer(I8P),                    public :: block_index=1_I8P       !< Block index in the field array.
   type(tree_node_object), pointer, public :: next=>null()            !< The next node in the tree.
   type(tree_node_object), pointer, public :: previous=>null()        !< The previous node in the tree.
   contains
      ! public methods
      ! procedure, pass(self) :: compute_morton_code !< Compute Morton code.
      procedure, pass(self) :: compute_finest_code !< Compute finest Morton code, the latest child of code given the level.
      procedure, pass(self) :: destroy             !< Destroy tree node.
      procedure, pass(self) :: initialize          !< Initialize tree node with key/content pair.
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

   !pure function key_str(l, tijk, bijk) result(key)
   !!< Return key in string format.
   !integer(I4P), intent(in) :: l       !< Refinementl level.
   !integer(I4P), intent(in) :: tijk(3) !< Tree coordinates into the forest.
   !integer(I4P), intent(in) :: bijk(3) !< Block coordinates into the tree.
   !character(len=KEY_LEN)   :: key     !< The key.

   !key = trim(str('(I7)',l))//trim(str('(I7)',tijk(1)))//trim(str('(I7)',tijk(2)))//trim(str('(I7)',tijk(3)))//&
   !                           trim(str('(I7)',bijk(1)))//trim(str('(I7)',bijk(2)))//trim(str('(I7)',bijk(3)))
   !endfunction key_str

   !subroutine key_int(key, l, tijk, bijk)
   !!< Return key in string format.
   !character(len=*), intent(in)            :: key     !< The key.
   !integer(I4P),     intent(out), optional :: l       !< Refinementl level.
   !integer(I4P),     intent(out), optional :: tijk(3) !< Tree coordinates into the forest.
   !integer(I4P),     intent(out), optional :: bijk(3) !< Block coordinates into the tree.

   !if (present(l)) l = cton(key(1:7), 1_I4P)
   !if (present(tijk)) then
   !   tijk(1) = cton(key(8 :14), 1_I4P)
   !   tijk(2) = cton(key(15:21), 1_I4P)
   !   tijk(3) = cton(key(22:28), 1_I4P)
   !endif
   !if (present(bijk)) then
   !   bijk(1) = cton(key(29:35), 1_I4P)
   !   bijk(2) = cton(key(36:42), 1_I4P)
   !   bijk(3) = cton(key(43:49), 1_I4P)
   !endif
   !endsubroutine key_int

   ! public methods
   !subroutine compute_morton_code(self, ratio)
   !!< Compute Morton code.
   !class(tree_node_object), intent(inout)        :: self    !< Tree node.
   !integer(I4P),                  intent(in), optional :: ratio   !< Refinement ratio.
   !integer(I4P)                                        :: ijkl(4) !< IJKL coordinates.
   !integer(I4P)                                        :: ratio_  !< Refinement ratio, local variable.

   !ratio_ = 8_I4P ; if (present(ratio)) ratio_ = ratio
   !call key_int(key=self%key, l=ijkl(4), bijk=ijkl(1:3))
   !if (   ratio_==8_I8P) then
   !   self%code = coordinates_to_morton(i=ijkl(1), j=ijkl(2), k=ijkl(3), l=ijkl(4), ratio=ratio_)
   !elseif(ratio_==4_I8P) then
   !   self%code = coordinates_to_morton(i=ijkl(1), j=ijkl(2), l=ijkl(4), ratio=ratio_)
   !endif
   !endsubroutine compute_morton_code

   elemental subroutine compute_finest_code(self, max_level)
   !< Compute finest Morton code, the latest child of code given the level.
   class(tree_node_object), intent(inout) :: self      !< Tree node.
   integer(I4P),            intent(in)    :: max_level !< Maximum refinement level considered for computing the finest code.
   integer(I4P)                           :: level_    !< Node level.
   integer(I4P)                           :: l         !< Counter.

   self%finest_code = self%code
   ! level_ = level(code=self%code)
   if (level_ < max_level) then
      do l=level_+1, max_level
         !< @TODO to be implemented
      enddo
   endif
   endsubroutine compute_finest_code

   elemental subroutine destroy(self)
   !< Destroy tree node.
   class(tree_node_object), intent(inout) :: self  !< Tree node.
   type(tree_node_object)                 :: fresh !< Fresh instance of tree node.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, code, content, finest_code, refinement_needed, myrank, block_index)
   !< Initialize tree node with key/content pair.
   class(tree_node_object), intent(inout)        :: self              !< Tree node.
   integer(I8P),            intent(in)           :: code              !< The Morton code.
   integer(I8P),            intent(in)           :: content           !< The content.
   integer(I8P),            intent(in), optional :: finest_code       !< The finest Morton code.
   integer(I4P),            intent(in), optional :: refinement_needed !< Flag for refinement/derefinement algorithm.
   integer(I4P),            intent(in), optional :: myrank            !< MPI rank process.
   integer(I8P),            intent(in), optional :: block_index       !< Block index in the field array.

   call self%destroy
   self%code = code
   self%content = content
   if (present(finest_code      )) self%finest_code       = finest_code
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
   lhs%finest_code = rhs%finest_code
   lhs%content = rhs%content
   lhs%refinement_needed = rhs%refinement_needed
   lhs%myrank = rhs%myrank
   lhs%block_index = rhs%block_index
   endsubroutine tree_node_assign_tree_node
endmodule adam_tree_node_object
