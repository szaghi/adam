!< ADAM, hash table class definition.
module adam_hash_table_object
!< ADAM, hash table class definition.

use adam_dictionary_node_object, only : KEY_LEN, destroy_dictionary_node, dictionary_node_object
use adam_dictionary_object, only : dictionary_object, iterator_interface, len
use PENF, only : I4P, I8P

implicit none
private
public :: hash_table_object

integer(I4P), parameter :: HT_BUCKETS_NUMBER_DEF = 9973_I4P !< Default number of buckets of hash table.

type :: hash_table_object
   !< Hash table class definition.
   type(dictionary_object), allocatable :: bucket(:)               !< Hash table buckets.
   integer(I4P)                         :: buckets_number=0_I4P    !< Number of buckets used.
   integer(I4P)                         :: nodes_number=0_I4P      !< Number of nodes actually stored, namely the hash table length.
   character(len=KEY_LEN), allocatable  :: keys(:,:)               !< Minimum and maximum key values actually stored.
   logical                              :: is_initialized_=.false. !< Initialization status.
   contains
      ! public methods
      procedure, pass(self) :: add_node     !< Add a node pointer to the hash table.
      procedure, pass(self) :: destroy      !< Destroy the hash table.
      procedure, pass(self) :: hash         !< Hash the key.
      procedure, pass(self) :: has_key      !< Check if the key is present in the hash table.
      procedure, pass(self) :: initialize   !< Initialize the hash table.
      procedure, pass(self) :: node_content !< Return node's content, given the key.
      procedure, pass(self) :: remove_node  !< Remove a node from the hash table, given the key.
      ! procedure, pass(self) :: traverse    !< Traverse hash table calling the iterator procedure.
endtype hash_table_object

contains
   ! public methods
   subroutine add_node(self, key, content)
   !< Add a node pointer to the hash table.
   !<
   !< @note If a node with the same key is already in the hash table, it is removed and the new one will replace it.
   class(hash_table_object), intent(inout) :: self          !< The hash table.
   character(len=*),         intent(in)    :: key           !< The key.
   integer(I8P),             intent(in)    :: content       !< The content.
   integer(I4P)                            :: b             !< Bucket index, namely hashed key.
   integer(I4P)                            :: bucket_length !< Length (nodes number) of the bucket where the node is placed.

   if (.not.self%is_initialized_) call self%initialize ! initialize the table with default options
   b = self%hash(key=key)
   bucket_length = len(self%bucket(b))
   call self%bucket(b)%add_node(key=key, content=content)
   self%nodes_number = self%nodes_number + (len(self%bucket(b)) - bucket_length)
   ! self%keys(1:2, b) = self%bucket(b)%keys()
   endsubroutine add_node

   subroutine destroy(self)
   !< Destroy the hash table.
   class(hash_table_object), intent(inout) :: self !< The hash table.
   integer(I4P)                            :: b    !< Counter.

   if (allocated(self%bucket)) then
      do b=1, size(self%bucket, dim=1)
        call self%bucket(b)%destroy
      enddo
      deallocate(self%bucket)
   endif
   self%buckets_number = 0_I4P
   self%nodes_number = 0_I4P
   if (allocated(self%keys)) deallocate(self%keys)
   endsubroutine destroy

   function has_key(self, key)
   !< Check if the key is present in the hash table.
   class(hash_table_object), intent(in) :: self    !< The hash table.
   character(len=*),         intent(in) :: key     !< Key to hash.
   logical                              :: has_key !< Check result.

   has_key = .false.
   if (self%is_initialized_) has_key = self%bucket(self%hash(key=key))%has_key(key=key)
   endfunction has_key

   elemental function hash(self, key) result(bucket)
   !< Hash the key.
   class(hash_table_object), intent(in) :: self   !< The hash table.
   character(len=*),         intent(in) :: key    !< Key to hash.
   integer(I4P)                         :: bucket !< Bucket index corresponding to the key.

   bucket = 0
   !< @TODO implement a real hash function
   ! if (self%is_initialized_) bucket = key%hash(buckets_number=self%buckets_number)
   endfunction hash

   subroutine initialize(self, buckets_number)
   !< Initialize the hash table.
   class(hash_table_object), intent(inout)        :: self           !< The hash table.
   integer(I4P),             intent(in), optional :: buckets_number !< Number of buckets for initialize the hash table.

   call self%destroy
   self%buckets_number = HT_BUCKETS_NUMBER_DEF ; if (present(buckets_number)) self%buckets_number = buckets_number
   allocate(self%bucket(1:self%buckets_number))
   allocate(self%keys(1:2,1:self%buckets_number))
   endsubroutine initialize

   function node_content(self, key) result(content)
   !< Return node's content, given the key.
   class(hash_table_object), intent(in)  :: self    !< The hash table.
   character(len=*),         intent(in)  :: key     !< The key.
   integer(I8P)                          :: content !< Content pointer of the queried node.
   type(dictionary_node_object), pointer :: p       !< Pointer to current node.

   content = 0
   if (self%is_initialized_) then
      p => self%bucket(self%hash(key=key))%node(key=key)
      if (associated(p)) content = p%content
   endif
   endfunction node_content

   subroutine remove_node(self, key)
   !< Remove a node from the hash table, given the key.
   class(hash_table_object), intent(inout) :: self !< The hash table.
   character(len=*),         intent(in)    :: key  !< The key.
   integer(I4P)                            :: b    !< Bucket index, namely hashed key.

   if (self%is_initialized_) then
     b = self%hash(key=key)
     call self%bucket(b)%remove_node(key=key)
     ! self%keys(1:2, b) = self%bucket(b)%keys()
   endif
   endsubroutine remove_node

   subroutine traverse(self, iterator)
   !< Traverse hash table calling the iterator procedure.
   class(hash_table_object), intent(in) :: self     !< The hash_table.
   procedure(iterator_interface)        :: iterator !< The (key) iterator procedure to call for each node.
   integer(I4P)                         :: b        !< Counter.

   if (self%is_initialized_) then
      do b=1, self%buckets_number
         call self%bucket(b)%traverse(iterator)
      enddo
   endif
   endsubroutine traverse
endmodule adam_hash_table_object
