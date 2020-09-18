!< ADAM, dictionary class definition.
module adam_dictionary_object
!< ADAM, dictionary class definition.

use adam_dictionary_node_object, only : KEY_LEN, destroy_dictionary_node, dictionary_node_object
use PENF, only : I4P, I8P

implicit none
private
public :: dictionary_object
public :: iterator_interface
public :: len

type :: dictionary_object
   !< Dictionary class definition.
   type(dictionary_node_object), pointer :: head=>null()    !< The first node in the dictionary.
   type(dictionary_node_object), pointer :: tail=>null()    !< The last node in the dictionary.
   integer(I4P)                          :: nodes_number=0  !< Number of nodes in the dictionary.
   character(len=KEY_LEN)                :: keys(2)=['',''] !< Minimum and maximum key values actually stored.
   contains
      ! public methods
      procedure, pass(self) :: add_node    !< Add a node pointer to the dictionary.
      procedure, pass(self) :: destroy     !< Destroy the dictionary.
      procedure, pass(self) :: has_key     !< Check if the key is present in the dictionary.
      procedure, pass(self) :: loop        !< Sentinel while-loop on nodes returning the key/content pair (for dictionary looping).
      procedure, pass(self) :: node        !< Return a pointer to a node.
      procedure, pass(self) :: remove_node !< Remove a node from the dictionary, given the key.
      procedure, pass(self) :: traverse    !< Traverse dictionary from head to tail calling the iterator procedure.
      ! private methods
      procedure, pass(self), private :: add_key           !< Add key to keys list.
      procedure, pass(self), private :: remove_by_pointer !< Remove node from dictionary, given pointer to it.
      procedure, pass(self), private :: remove_key        !< Remove key to keys list.
      ! operators
      ! generic :: assignment(=) => dictionary_assign_dictionary      !< Overload `=`.
      ! procedure, pass(lhs), private :: dictionary_assign_dictionary !< Operator `=`.
endtype dictionary_object

abstract interface
   subroutine iterator_interface(node, done)
   !< Iterator procedure for traversing all nodes in a dictionary.
   import :: dictionary_node_object
   type(dictionary_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
   logical,                               intent(out) :: done !< Flag to set to true to stop traversing.
   endsubroutine iterator_interface
endinterface

interface len
  !< Overload `len` builtin for accepting a [[dictionary]].
  module procedure dictionary_len
endinterface

contains
   ! public non TBP
   elemental function dictionary_len(self) result(length)
   !< Return the number of nodes of the dictionary, namely the dictionary length.
   type(dictionary_object), intent(in) :: self   !< The dictionary.
   integer(I4P)                        :: length !< The dictionary length.

   length = self%nodes_number
   endfunction dictionary_len

   ! public methods
   subroutine add_node(self, key, content)
   !< Add a node pointer to the dictionary.
   !<
   !< @note If a node with the same key is already in the dictionary, it is removed and the new one will replace it.
   class(dictionary_object), intent(inout) :: self    !< The dictionary.
   character(len=*),         intent(in)    :: key     !< The key.
   integer(I8P),             intent(in)    :: content !< The content.
   type(dictionary_node_object), pointer   :: p       !< Pointer to scan the dictionary.

   ! if the node is already there, then remove it
   p => self%node(key=key)
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

   call p%initialize(key=key, content=content) ! fill the new node with provided contents

   call self%add_key(key=p%key)

   self%nodes_number = self%nodes_number + 1
   endsubroutine add_node

   subroutine destroy(self)
   !< Destroy the dictionary.
   class(dictionary_object), intent(inout) :: self !< The dictionary.

   if (associated(self%head)) call destroy_dictionary_node(node=self%head)
   self%head => null()
   self%tail => null()
   self%nodes_number = 0
   self%keys = ['','']
   endsubroutine destroy

   function has_key(self, key)
   !< Check if the key is present in the dictionary.
   class(dictionary_object), intent(in)  :: self    !< The dictionary.
   character(len=*),         intent(in)  :: key     !< The key.
   logical                               :: has_key !< Check result.

   has_key = .false.
   call self%traverse(iterator=key_iterator_search)
   contains
      subroutine key_iterator_search(node, done)
      !< Iterator procedure for searching a key.
      type(dictionary_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
      logical,                               intent(out) :: done !< Flag to set to true to stop traversing.

      has_key = node%key==key
      done = has_key
      endsubroutine key_iterator_search
   endfunction has_key

   function loop(self, key, content) result(again)
   !< Sentinel while-loop on nodes returning the key/content pair (for dictionary looping).
   class(dictionary_object), intent(in)            :: self      !< The dictionary.
   character(len=*),         intent(out), optional :: key       !< The key.
   integer(I8P),             intent(out), optional :: content   !< The content.
   logical                                         :: again     !< Sentinel flag to contine the loop.
   type(dictionary_node_object), pointer, save     :: p=>null() !< Pointer to current node.

   again = .false.
   if (present(content)) content = 0
   if (self%nodes_number>0) then
      if (.not.associated(p)) then
         p => self%head
         if (present(key)) key = p%key
         if (present(content)) content = p%content
         again = .true.
      elseif (associated(p%next)) then
         p => p%next
         if (present(key)) key = p%key
         if (present(content)) content = p%content
         again = .true.
      else
         p => null()
         again = .false.
      endif
   endif
   endfunction loop

   function node(self, key) result(p)
   !< Return a pointer to a node in the dictionary.
   class(dictionary_object), intent(in)  :: self !< The dictionary.
   character(len=*),         intent(in)  :: key  !< The key.
   type(dictionary_node_object), pointer :: p    !< Pointer to node queried.

   p => null()
   call self%traverse(iterator=key_iterator_search)
   contains
      subroutine key_iterator_search(node, done)
      !< Iterator procedure for searching a key.
      type(dictionary_node_object), pointer, intent(in)  :: node !< Actual node pointer in the dictionary.
      logical,                               intent(out) :: done !< Flag to set to true to stop traversing.

      done = node%key==key
      if (done) p => node
      endsubroutine key_iterator_search
   endfunction node

   subroutine remove_node(self, key)
   !< Remove a node from the dictionary, given the key.
   class(dictionary_object), intent(inout) :: self !< The dictionary.
   character(len=*),         intent(in)    :: key  !< The key.
   type(dictionary_node_object), pointer   :: p    !< Pointer to scan the dictionary.

   p => self%node(key=key)
   if (associated(p)) call self%remove_by_pointer(p=p)
   endsubroutine remove_node

   subroutine traverse(self, iterator)
   !< Traverse dictionary from head to tail calling the iterator procedure.
   class(dictionary_object), intent(in)  :: self     !< The dictionary.
   procedure(iterator_interface)         :: iterator !< The iterator procedure to call for each node.
   type(dictionary_node_object), pointer :: p        !< Pointer to scan the dictionary.
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
   pure subroutine add_key(self, key)
   !< Add key to keys list.
   class(dictionary_object), intent(inout) :: self !< The dictionary.
   character(len=*),         intent(in)    :: key  !< The key.

   !< @TODO to be implemented
   endsubroutine add_key

   subroutine remove_by_pointer(self, p)
   !< Remove node from dictionary, given pointer to it.
   class(dictionary_object),              intent(inout) :: self         !< The dictionary.
   type(dictionary_node_object), pointer, intent(inout) :: p            !< Pointer to the node to remove.
   logical                                              :: has_next     !< Check if dictionary node has a next item.
   logical                                              :: has_previous !< Check if dictionary node has a previous item.

   if (associated(p)) then
     call self%remove_key(key=p%key)
     call p%destroy ! destroy the node contents
     has_next     = associated(p%next)
     has_previous = associated(p%previous)
     if (has_next.and.has_previous) then ! neither first nor last in dictionary
       p%previous%next => p%next
       p%next%previous => p%previous
     elseif (has_next.and.(.not.has_previous)) then ! first one in dictionary
       self%head          => p%next
       self%head%previous => null()
     elseif (has_previous.and.(.not.has_next)) then ! last one in dictionary
       self%tail      => p%previous
       self%tail%next => null()
     elseif ((.not.has_previous).and.(.not.has_next)) then ! only one in the dictionary
       self%head => null()
       self%tail => null()
     endif
     deallocate(p)
     p => null()
     self%nodes_number = self%nodes_number - 1
   endif
   endsubroutine remove_by_pointer

   subroutine remove_key(self, key)
   !< Remove key to keys list.
   class(dictionary_object), intent(inout) :: self !< The dictionary.
   character(len=*),         intent(in)    :: key  !< The key.

   !< @TODO to be implemented
   endsubroutine remove_key
endmodule adam_dictionary_object
