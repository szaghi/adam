!< ADAM, tree bucket class definition.
module adam_tree_bucket_object
!< ADAM, tree bucket class definition.
!< The bucket is implemented as a dictionary based on a double linked list.

use adam_tree_node_object, only : destroy_tree_node, tree_node_object
use PENF, only : I4P, I8P, str

implicit none
private
public :: tree_bucket_object
public :: iterator_interface
public :: len

type :: tree_bucket_object
   !< tree bucket class definition.
   type(tree_node_object), pointer :: head=>null()    !< The first node in the tree bucket.
   type(tree_node_object), pointer :: tail=>null()    !< The last node in the tree bucket.
   integer(I4P)                    :: nodes_number=0  !< Number of nodes in the tree bucket.
   integer(I8P)                    :: code(1:2)=[0,0] !< Minimum and maximum unique code values actually stored.
   contains
      ! public methods
      procedure, pass(self) :: add_node    !< Add a node pointer to the tree bucket.
      procedure, pass(self) :: destroy     !< Destroy the tree bucket.
      procedure, pass(self) :: has_code    !< Check if the code is present in the tree bucket.
      procedure, pass(self) :: loop        !< Sentinel while-loop on nodes returning the code/content pair.
      procedure, pass(self) :: node        !< Return a pointer to a node.
      procedure, pass(self) :: remove_node !< Remove a node from the tree bucket, given the code.
      procedure, pass(self) :: traverse    !< Traverse tree bucket from head to tail calling the iterator procedure.
      ! private methods
      procedure, pass(self), private :: add_code          !< Add a code to codes list.
      procedure, pass(self), private :: remove_by_pointer !< Remove node from tree bucket, given pointer to it.
      procedure, pass(self), private :: remove_code       !< Remove a code from codes list.
endtype tree_bucket_object

abstract interface
   subroutine iterator_interface(node, done)
   !< Iterator procedure for traversing all nodes in a tree bucket.
   import :: tree_node_object
   type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the tree bucket.
   logical,                         intent(out) :: done !< Flag to set to true to stop traversing.
   endsubroutine iterator_interface
endinterface

interface len
  !< Overload `len` builtin for accepting a [[tree_bucket_object]].
  module procedure tree_bucket_len
endinterface

contains
   ! public non TBP
   elemental function tree_bucket_len(self) result(length)
   !< Return the number of nodes of the tree bucket, namely the tree bucket length.
   type(tree_bucket_object), intent(in) :: self   !< The tree bucket.
   integer(I4P)                         :: length !< The tree bucket length.

   length = self%nodes_number
   endfunction tree_bucket_len

   ! public methods
   subroutine add_node(self, code, content, finest_code, refinement_needed, myrank, block_index)
   !< Add a node pointer to the tree bucket.
   !<
   !< @note If a node with the same code is already in the tree bucket, it is removed and the new one will replace it.
   class(tree_bucket_object), intent(inout)        :: self              !< The tree bucket.
   integer(I8P),              intent(in)           :: code              !< The Morton code.
   integer(I8P),              intent(in)           :: content           !< The content.
   integer(I8P),              intent(in), optional :: finest_code       !< The finest Morton code.
   integer(I4P),              intent(in), optional :: refinement_needed !< Flag for refinement/derefinement algorithm.
   integer(I4P),              intent(in), optional :: myrank            !< MPI rank process.
   integer(I8P),              intent(in), optional :: block_index       !< Block index in the field array.
   type(tree_node_object), pointer                 :: p                 !< Pointer to scan the tree bucket.

   ! if the node is already there, then remove it
   p => self%node(code=code)
   if (associated(p)) call self%remove_by_pointer(p)

   ! update next/previous pointers
   if (associated(self%tail)) then ! insert new node at the end
      allocate(self%tail%next)
      p => self%tail%next
      p%previous => self%tail
   else
      allocate(self%head) ! insert new node as first node
      p => self%head
   end if
   self%tail => p

   call p%initialize(code=code, content=content, finest_code=finest_code, refinement_needed=refinement_needed, &
                     myrank=myrank, block_index=block_index)

   call self%add_code(code=p%code)

   self%nodes_number = self%nodes_number + 1
   endsubroutine add_node

   subroutine destroy(self)
   !< Destroy the tree bucket.
   class(tree_bucket_object), intent(inout) :: self !< The tree bucket.

   if (associated(self%head)) call destroy_tree_node(node=self%head)
   self%head => null()
   self%tail => null()
   self%nodes_number = 0
   endsubroutine destroy

   function has_code(self, code)
   !< Check if the code is present in the tree bucket.
   class(tree_bucket_object), intent(in)  :: self     !< The tree bucket.
   integer(I8P),              intent(in)  :: code     !< The Morton code.
   logical                                :: has_code !< Check result.

   has_code = .false.
   call self%traverse(iterator=code_iterator_search)
   contains
      subroutine code_iterator_search(node, done)
      !< Iterator procedure for searching a code.
      type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the tree bucket.
      logical,                         intent(out) :: done !< Flag to set to true to stop traversing.

      has_code = node%code==code
      done = has_code
      endsubroutine code_iterator_search
   endfunction has_code

   function loop(self, code, content) result(again)
   !< Sentinel while-loop on nodes returning the code/content pair (for tree bucket looping).
   class(tree_bucket_object), intent(in)            :: self      !< The tree bucket.
   integer(I8P),              intent(out), optional :: code      !< The Morton code.
   integer(I8P),              intent(out), optional :: content   !< The content.
   logical                                          :: again     !< Sentinel flag to contine the loop.
   type(tree_node_object), pointer, save            :: p=>null() !< Pointer to current node.

   again = .false.
   if (present(content)) content = 0
   if (self%nodes_number>0) then
      if (.not.associated(p)) then
         p => self%head
         if (present(code)) code = p%code
         if (present(content)) content = p%content
         again = .true.
      elseif (associated(p%next)) then
         p => p%next
         if (present(code)) code = p%code
         if (present(content)) content = p%content
         again = .true.
      else
         p => null()
         again = .false.
      endif
   endif
   endfunction loop

   function node(self, code) result(p)
   !< Return a pointer to a node in the tree bucket.
   class(tree_bucket_object), intent(in)  :: self !< The tree bucket.
   integer(I8P),              intent(in)  :: code !< The Morton code.
   type(tree_node_object), pointer        :: p    !< Pointer to node queried.

   p => null()
   call self%traverse(iterator=code_iterator_search)
   contains
      subroutine code_iterator_search(node, done)
      !< Iterator procedure for searching a code.
      type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the tree bucket.
      logical,                         intent(out) :: done !< Flag to set to true to stop traversing.

      done = node%code==code
      if (done) p => node
      endsubroutine code_iterator_search
   endfunction node

   subroutine remove_node(self, code)
   !< Remove a node from the tree bucket, given the code.
   class(tree_bucket_object), intent(inout) :: self !< The tree bucket.
   integer(I8P),              intent(in)    :: code !< The Morton code.
   type(tree_node_object), pointer          :: p    !< Pointer to scan the tree bucket.

   p => self%node(code=code)
   if (associated(p)) call self%remove_by_pointer(p=p)
   endsubroutine remove_node

   subroutine traverse(self, iterator)
   !< Traverse tree bucket from head to tail calling the iterator procedure.
   class(tree_bucket_object), intent(in) :: self     !< The tree bucket.
   procedure(iterator_interface)         :: iterator !< The iterator procedure to call for each node.
   type(tree_node_object), pointer       :: p        !< Pointer to scan the tree bucket.
   logical                               :: done     !< Flag to set to true to stop traversing.

   done = .false.
   p => self%head
   do
     if (associated(p)) then
       call iterator(node=p, done=done)
       if (done) exit
       p => p%next
     else
       exit
     endif
   enddo
   endsubroutine traverse

   ! private methods
   pure subroutine add_code(self, code)
   !< Add a code to minimum and maximum unique code values.
   class(tree_bucket_object), intent(inout) :: self !< The tree bucket.
   integer(I8P),              intent(in)    :: code !< The Morton code.

   self%code(1) = min(self%code(1), code) ; if (self%code(1)==0) self%code(1) = code
   self%code(2) = max(self%code(2), code)
   endsubroutine add_code

   subroutine remove_by_pointer(self, p)
   !< Remove node from tree bucket, given pointer to it.
   class(tree_bucket_object),       intent(inout) :: self         !< The tree bucket.
   type(tree_node_object), pointer, intent(inout) :: p            !< Pointer to the node to remove.
   logical                                        :: has_next     !< Check if tree bucket node has a next item.
   logical                                        :: has_previous !< Check if tree bucket node has a previous item.

   if (associated(p)) then
     call self%remove_code(code=p%code)
     call p%destroy ! destroy the node contents
     has_next     = associated(p%next)
     has_previous = associated(p%previous)
     if (has_next.and.has_previous) then ! neither first nor last in tree bucket
       p%previous%next => p%next
       p%next%previous => p%previous
     elseif (has_next.and.(.not.has_previous)) then ! first one in tree bucket
       self%head          => p%next
       self%head%previous => null()
     elseif (has_previous.and.(.not.has_next)) then ! last one in tree bucket
       self%tail      => p%previous
       self%tail%next => null()
     elseif ((.not.has_previous).and.(.not.has_next)) then ! only one in the tree bucket
       self%head => null()
       self%tail => null()
     endif
     deallocate(p)
     p => null()
     self%nodes_number = self%nodes_number - 1
   endif
   endsubroutine remove_by_pointer

   subroutine remove_code(self, code)
   !< Remove a code to minimum and maximum unique code values.
   class(tree_bucket_object), intent(inout) :: self !< The tree bucket.
   integer(I8P),              intent(in)    :: code !< The Morton code.
   type(tree_node_object), pointer          :: p    !< Pointer to scan the tree bucket.

   if (self%nodes_number==1) then
      self%code = 0
   elseif (self%nodes_number>=2) then
      p => null()
      if (self%code(1)==code) then
         call self%traverse(iterator=code_iterator_search)
         if (associated(p)) then
            if (associated(p%next)) self%code(1) = p%next%code
         endif
      elseif (self%code(2)==code) then
         call self%traverse(iterator=code_iterator_search)
         if (associated(p)) then
            if (associated(p%previous)) self%code(2) = p%previous%code
         endif
      endif
   endif
   contains
     subroutine code_iterator_search(node, done)
     !< Iterator procedure for searching a code.
     type(tree_node_object), pointer, intent(in)  :: node !< Actual node pointer in the tree bucket.
     logical,                         intent(out) :: done !< Flag to set to true to stop traversing.

     done = node%code==code
     if (done) p => node
     endsubroutine code_iterator_search
   endsubroutine remove_code
endmodule adam_tree_bucket_object
