!< ADAM, hash table class definition.
module adam_hash_table_object
!< ADAM, hash table class definition.

use adam_dictionary_node_object, only : KEY_LEN, destroy_dictionary_node, dictionary_node_object
use adam_dictionary_object, only : dictionary_object, iterator_interface, len
use PENF, only : I1P, I4P, I8P, R8P

implicit none
private
public :: hash_table_object

integer(I4P), parameter :: HT_BUCKETS_NUMBER_DEF = 9973_I4P !< Default number of buckets of hash table.
real(R8P),    parameter :: HT_MAX_LOAD = 0.9_R8P            !< Maximum load of hash table buckets.

type :: hash_table_object
   !< Hash table class definition.
   type(dictionary_object), allocatable :: bucket(:)               !< Hash table buckets.
   integer(I8P), allocatable            :: keys(:,:)               !< Minimum and maximum key values actually stored.
   integer(I4P)                         :: buckets_number=0_I4P    !< Number of buckets used.
   integer(I4P)                         :: nodes_number=0_I4P      !< Number of nodes actually stored, namely the hash table length.
   real(R8P)                            :: max_load=HT_MAX_LOAD    !< Maximum load of hash table buckets.
   logical                              :: is_initialized_=.false. !< Initialization status.
   contains
      ! public methods
      procedure, pass(self) :: add_node             !< Add a node pointer to the hash table.
      procedure, pass(self) :: prime_buckets_number !< Return the buckets number as the nearest prime number given nodes number.
      procedure, pass(self) :: destroy              !< Destroy the hash table.
      procedure, pass(self) :: hash                 !< Hash the key.
      procedure, pass(self) :: has_key              !< Check if the key is present in the hash table.
      procedure, pass(self) :: initialize           !< Initialize the hash table.
      procedure, pass(self) :: node                 !< Return a pointer to a node.
      procedure, pass(self) :: node_content         !< Return node's content, given the key.
      procedure, pass(self) :: remove_node          !< Remove a node from the hash table, given the key.
      procedure, pass(self) :: resize               !< Resize the hash table.
      procedure, pass(self) :: traverse             !< Traverse hash table calling the iterator procedure.
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

   if (.not.self%is_initialized_) call self%initialize ! initialize the table with default options
   b = self%hash(key=key)
   call self%bucket(b)%add_node(key=key, content=content)
   self%nodes_number = self%nodes_number + 1
   ! self%keys(1:2, b) = self%bucket(b)%keys()
   endsubroutine add_node

   elemental function prime_buckets_number(self, nodes_number) result(buckets_number)
   !< Return the buckets number as the nearest prime number given nodes number.
   !<
   !< @note The balanced buckets number is computing considering the hash table load defined in `self` and using the
   !< Sieve of Eratoshenes for findining the nearest prime number.
   class(hash_table_object), intent(in) :: self           !< The hash_table.
   integer(I4P),             intent(in) :: nodes_number   !< Nodes number to be stored in the hash table.
   integer(I4P)                         :: buckets_number !< Well balanced, prime buckets number.
   logical, allocatable                 :: is_prime(:)    !< List of prime numbers up to buckets number.
   integer(I4P)                         :: b              !< Counter.

   buckets_number = int((1._R8P / self%max_load) * nodes_number, I4P)
   allocate(is_prime(buckets_number))
   is_prime = .true.
   is_prime(1) = .false.
   do b=2, int(sqrt(real(buckets_number, R8P)), I4P)
      if (is_prime(b)) is_prime(b*b:buckets_number:b) = .false.
   enddo
   b = buckets_number
   do while(.not.is_prime(b))
      b = b - 1
   enddo
   buckets_number = b
   endfunction prime_buckets_number

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
   if (allocated(self%keys)) deallocate(self%keys)
   self%buckets_number = 0_I4P
   self%nodes_number = 0_I4P
   self%max_load = HT_MAX_LOAD
   self%is_initialized_ = .false.
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
   if (self%is_initialized_) bucket = abs(mod(murmurhash3(key=key), self%buckets_number)) + 1
   endfunction hash

   subroutine initialize(self, max_load, nodes_number, buckets_number)
   !< Initialize the hash table.
   class(hash_table_object), intent(inout)        :: self           !< The hash table.
   real(R8P),                intent(in), optional :: max_load       !< Maximum load of hash table buckets.
   integer(I4P),             intent(in), optional :: nodes_number   !< Nodes number to be stored in the hash table.
   integer(I4P),             intent(in), optional :: buckets_number !< Number of buckets for initialize the hash table.

   call self%destroy
   if (present(max_load)) self%max_load = max_load
   if (present(nodes_number)) then
      self%buckets_number = self%prime_buckets_number(nodes_number=nodes_number)
   else
      self%buckets_number = HT_BUCKETS_NUMBER_DEF ; if (present(buckets_number)) self%buckets_number = buckets_number
   endif
   allocate(self%bucket(1:self%buckets_number))
   allocate(self%keys(1:2,1:self%buckets_number))
   self%is_initialized_ = .true.
   endsubroutine initialize

   function node(self, key) result(p)
   !< Return a pointer to a node in the dictionary.
   class(hash_table_object), intent(in)  :: self !< The hash table.
   character(len=*),         intent(in)  :: key  !< The key.
   type(dictionary_node_object), pointer :: p    !< Pointer to node queried.

   p => null()
   if (self%is_initialized_) p => self%bucket(self%hash(key=key))%node(key=key)
   endfunction node

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
      self%nodes_number = self%nodes_number - 1
      ! self%keys(1:2, b) = self%bucket(b)%keys()
   endif
   endsubroutine remove_node

   subroutine resize(self, nodes_number, max_load)
   !< Resize the hash table.
   class(hash_table_object), intent(inout)        :: self         !< The hash table.
   integer(I4P),             intent(in)           :: nodes_number !< Nodes number to be stored in the hash table.
   real(R8P),                intent(in), optional :: max_load     !< Maximum load of hash table buckets.
   type(hash_table_object)                        :: swap         !< Temporary (swap) hash table.
   character(len=KEY_LEN)                         :: key          !< Hash table node key.
   integer(I8P)                                   :: content      !< Hash table node content.
   integer(I4P)                                   :: b            !< Counter.

   if (self%is_initialized_) then
      if (present(max_load)) self%max_load = max_load
      if (self%nodes_number > int((1._R8P/self%max_load)*nodes_number, I4P)) return ! new size too small, cannot previous nodes
      call swap%initialize(max_load=self%max_load, nodes_number=nodes_number)
      do b=1, self%buckets_number
         do while(self%bucket(b)%loop(key=key, content=content))
            call swap%add_node(key=key, content=content)
         enddo
      enddo
      call move_alloc(from=swap%bucket, to=self%bucket)
      call move_alloc(from=swap%keys, to=self%keys)
      self%buckets_number  = swap%buckets_number
      self%nodes_number    = swap%nodes_number
      self%max_load        = swap%max_load
      self%is_initialized_ = swap%is_initialized_
   else
      call self%initialize(nodes_number=nodes_number)
   endif
   endsubroutine resize

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

   ! private non TBP
   pure function murmurhash3(key) result(hash)
   !< MurMurHash v3, implementation taken by https://github.com/jannisteunissen/murmur3-fortran/blob/master/m_murmur3.f90.
   character(len=*), intent(in) :: key                    !< The key.
   integer(I4P)                 :: hash                   !< The hash.
   integer(I4P)                 :: klen                   !< The key length.
   integer(I4P), parameter      :: seed=42                !< Rondomizing seed.
   integer(I4P)                 :: i, i0, n, nblocks      !< Counters.
   integer(I4P)                 :: h1, k1                 !< Counters.
   integer(I4P), parameter      :: c1        = -862048943 !< 0xcc9e2d51.
   integer(I4P), parameter      :: c2        = 461845907  !< 0x1b873593.
   integer(I4P), parameter      :: shifts(3) = [0, 8, 16] !< Shift offsets.

   klen=len(key)
   h1      = seed
   ! nblocks = shiftr(klen, 2)    ! nblocks/4
   nblocks = ishft(klen, -2)    ! nblocks/4

   ! body
   do i = 1, nblocks
      k1 = transfer(key(i*4-3:i*4), k1)

      k1 = k1 * c1
      k1 = rotl32(k1,15_I1P)
      k1 = k1 * c2

      h1 = ieor(h1, k1)
      h1 = rotl32(h1,13_I1P)
      h1 = h1 * 5 - 430675100  ! 0xe6546b64
   end do

   ! tail
   k1 = 0
   i  = iand(klen, 3)
   i0 = 4 * nblocks

   do n = i, 1, -1
      ! k1 = ieor(k1, shiftl(iachar(key(i0+n:i0+n)), shifts(n)))
      k1 = ieor(k1, ishft(iachar(key(i0+n:i0+n)), shifts(n)))
   end do

   ! Check if the above loop was executed
   if (i >= 1) then
      k1 = k1 * c1
      k1 = rotl32(k1,15_I1P)
      k1 = k1 * c2
      h1 = ieor(h1, k1)
   end if

   ! finalization
   h1 = ieor(h1, klen)
   h1 = fmix32(h1)
   hash = h1
   contains
      pure function rotl32(x, r)
      integer(I4P), intent(in) :: x
      integer(I1P), intent(in) :: r
      integer(I4P)             :: rotl32
      ! rotl32 = ior(shiftl(x, r), shiftr(x, (32 - r)))
      rotl32 = ior(ishft(x, r), ishft(x, -(32 - r)))
      endfunction rotl32

      pure function fmix32(h_in) result(h)
      !< Finalization mix - force all bits of a hash block to avalanche.
      integer(I4P), intent(in) :: h_in
      integer(I4P)             :: h
      h = h_in
      ! h = ieor(h, shiftr(h, 16))
      h = ieor(h, ishft(h, -16))
      h = h * (-2048144789) !0x85ebca6b
      ! h = ieor(h, shiftr(h, 13))
      h = ieor(h, ishft(h, -13))
      h = h * (-1028477387) !0xc2b2ae35
      ! h = ieor(h, shiftr(h, 16))
      h = ieor(h, ishft(h, -16))
      endfunction fmix32
   endfunction murmurhash3
endmodule adam_hash_table_object
