!< ADAM, tree class definition.
module adam_tree_object
!< ADAM, tree class definition.
!< The tree data structure is organized as a hash table and can be arranged as octree or quadtree simply adopting the proper
!< refinement ratio, namely 8 or 4 ratio respectively.

!< The tree exploits the Morton order for the linearization of the tree. A prototype of quadtree is represented below.
!<
!<            j^
!<    L=0      |
!<      L=1    |
!<        L=2  |
!<          L=3|
!<             |------------------------------------------------
!<          7  |  10   11  |  14   15  |  58   59  |  62   63  |
!<        3----|    10     |    11     |    14     |    15     |
!<          6  |  8    9   |  12   13  |  56   57  |  60   61  |
!<      1------|---------- 2 ----------|---------- 3 ----------|
!<          5  |  34   45  |  38   39  |  50   51  |  54   55  |
!<        2----|     8     |     9     |    12     |    13     |
!<          4  |  32   33  |  36   37  |  48   49  |  52   53  |
!<    0--------|-----------|--------- -1 ----------|-----------|
!<          3  |  10   11  |  14   15  |  26   27  |  30   31  |
!<        1----|     2     |     3     |     6     |     7     |
!<          2  |  8    9   |  12   13  |  24   25  |  28   29  |
!<      0------|---------- 0 ----------|---------- 1 ----------|
!<          1  |  2    3   |  6    7   |  18   19  |  22   23  |
!<        0----|     0     |     1     |     4     |     5     |
!<          0  |  0    1   |  4    5   |  16   17  |  20   21  |
!<             O------------------------------------------------------>
!<                0 |  1   |  2 |  3   |  4 |  5   |  6 |  7     L=3  i
!<                  0      |    1      |    2      |    3        L=2
!<                         0           |           1             L=1
!<                                     0                         L=0
!<
!< In the above representation each refinement level is represented *alone*, namely without taking into account the existence of
!< previous levels. However, the real hierarchy can be taken into account by simply adding the offset derived from the previous
!< levels, thus the numbering becomes as below.
!<
!<            j^
!<    L=0      |
!<      L=1    |
!<        L=2  |
!<          L=3|
!<             |------------------------------------------------
!<          7  |  62   63  |  66   67  |  78   79  |  82   63  |
!<        3----|    14     |    15     |    18     |    19     |
!<          6  |  60   61  |  64   65  |  76   77  |  80   81  |
!<      1------|---------- 2 ----------|---------- 3 ----------|
!<          5  |  54   55  |  58   59  |  70   71  |  74   75  |
!<        2----|    12     |    13     |    16     |    17     |
!<          4  |  52   53  |  56   57  |  68   69  |  72   73  |
!<    0--------|-----------|--------- -1 ----------|-----------|
!<          3  |  30   31  |  34   35  |  46   47  |  50   51  |
!<        1----|     6     |     7     |    10     |    11     |
!<          2  |  28   29  |  32   33  |  44   45  |  48   49  |
!<      0------|---------- 0 ----------|---------- 1 ----------|
!<          1  |  22   23  |  26   27  |  38   39  |  42   43  |
!<        0----|     4     |     5     |     8     |     9     |
!<          0  |  20   21  |  24   25  |  36   37  |  40   41  |
!<             O------------------------------------------------------>
!<                0 |  1   |  2 |  3   |  4 |  5   |  6 |  7     L=3  i
!<                  0      |    1      |    2      |    3        L=2
!<                         0           |           1             L=1
!<                                     0                         L=0
!<
!< This last numbering is the complete Morton order where the Morton code ideintifing a node entails all the spatial information,
!< the refinement level L and the spatial coordinatates IJ.
!< Sometimes it is convenient to use the representation where level L is not encoded into the Morton order, but the conversion is as
!< as adding the offset of previous level, namely using the funcion `first_at_level`.
!< Note that the ancestor node has Morton code -1, it is the ancestor of all nodes, it is **Adam**.
!< For the octree case (3D case) the ordering is equivalent, simply the local children fall in [0,7] instead of [0,3], namely the
!< local numbering is always [0, ratio-1], as represented below.
!<
!<  ^
!< /|\Z
!<  |
!<  |                            *----------*----------*
!<  |                           /|         /|         /|
!<  |                          / |        / |        / |
!<  |                         /  |       /  |       /  |
!<  |                        /   |      /   |      /   |
!<  |                       /    *-----/----*-----/----*
!<  |            6 <-------/--+ /|    /    /|    /   +---------> 7
!<  |                     *----------*----------*    / |
!<  |                    /|   /  |  /|   /  |  /|   /  |
!<  |                   / |  /   | / |  /   | / |  /   |
!<  |                  /  | /    */--|-/----*/--|-/----*
!<  |         2 <-----/---|/--+ //   |/    //   |/  +---------> 3
!<  |                /    *-----/----*-----/----*    /
!<  |    4 <--------/--+ /|   //    /|   //   +------------> 5
!<  |              *----------*----------*    / |  /
!<  |              |   /  | / |   /  | / |   /  | /
!<  |              |  /   |/  |  /   |/  |  /   |/
!<  |              | /    *---|-/----*---|-/----*
!<  |      0 <-----|/--+ /    |/    /    |/  +------------> 1
!<  |              *----------*----------*    /
!<  |              |   /      |   /      |   /
!<  |              |  /       |  /       |  /
!<  |   _ Y        | /        | /        | /
!<  |   /|         |/         |/         |/
!<  |  /           *----------*----------*
!<  | /
!<  |/                                                    X
!<  o------------------------------------------------------------------->

use adam_tree_node_object, only : destroy_tree_node, tree_node_object, TO_BE_REFINED, TO_BE_DEREFINED, TO_NOT_TOUCH
use adam_tree_bucket_object, only : tree_bucket_object, iterator_interface, len
use MORTIF, only : morton2D, morton3D, demorton2D, demorton3D
use PENF, only : I1P, I4P, I8P, R8P, str

implicit none
private
public :: tree_object

! tree defaults
integer(I4P), parameter :: TREE_BUCKETS_NUMBER_DEF = 9973_I4P    !< Default number of buckets of hash table.
real(R8P),    parameter :: TREE_MAX_LOAD = 0.9_R8P               !< Maximum load of hash table buckets.
integer(I4P), parameter :: TREE_MAX_SANITIZE_ITERATIONS = 10_I4P !< Default number of tree sanitize iterations.
! nodes types
integer(I4P), parameter :: STANDARD_NODE = 0_I4P           !< Standard node type.
integer(I4P), parameter :: MORE_REFINED_NODE = 1_I4P       !< More refined node type.
integer(I4P), parameter :: BOUNDARY_CONDITION_NODE = 2_I4P !< Boundary condition node type.

type :: tree_object
   !< Tree class definition.
   type(tree_bucket_object), allocatable :: bucket(:)               !< Tree buckets.
   integer(I8P), allocatable             :: code(:,:)               !< Min and max code values actually stored [2,buckets_number].
   integer(I4P)                          :: buckets_number=0_I4P    !< Number of buckets used.
   integer(I4P)                          :: nodes_number=0_I4P      !< Number of nodes actually stored, namely the tree length.
   real(R8P)                             :: max_load=TREE_MAX_LOAD  !< Maximum load of tree buckets.
   integer(I4P)                          :: ratio=8_I4P             !< Refinement ratio.
   integer(I4P)                          :: max_level=12_I4P        !< Maximum refinement level.
   logical                               :: is_initialized_=.false. !< Initialization status.
   integer(I8P), allocatable             :: to_refine(:)            !< List of nodes to be refined.
   ! type(tree_bucket_object)              :: to_derefine             !< List of node to be derefined.
   contains
      ! public methods
      procedure, pass(self) :: add_node             !< Add a node pointer to the tree.
      procedure, pass(self) :: codes                !< Return the list of (sorted) codes actually stored in the tree.
      procedure, pass(self) :: destroy              !< Destroy the tree.
      procedure, pass(self) :: loop                 !< Sentinel while-loop on nodes returning the code/content pair.
      procedure, pass(self) :: hash                 !< Hash the key.
      procedure, pass(self) :: has_code             !< Check if the code is present in the tree.
      procedure, pass(self) :: initialize           !< Initialize the tree.
      procedure, pass(self) :: node                 !< Return a pointer to a node.
      procedure, pass(self) :: node_content         !< Return node's content, given the key.
      procedure, pass(self) :: prime_buckets_number !< Return the buckets number as the nearest prime number given nodes number.
      procedure, pass(self) :: refine               !< Refine nodes.
      procedure, pass(self) :: remove_node          !< Remove a node from the tree, given the key.
      procedure, pass(self) :: resize               !< Resize the tree.
      procedure, pass(self) :: sanitize             !< Sanitize the tree.
      procedure, pass(self) :: traverse             !< Traverse tree calling the iterator procedure.
      ! Morton ordering methods
      generic               :: coordinates_to_morton => &
                               coordinates3D_to_morton, &
                               coordinates2D_to_morton !< Return the Morton code given space-level coordinates.
      generic               :: morton_to_coordinates => &
                               morton_to_coordinates3D, &
                               morton_to_coordinates2D !< Return the space-level coordinates given Morton code.
      procedure, pass(self) :: child                   !< Return the i-th child given Morton code.
      procedure, pass(self) :: child_local             !< Return the child index in the local numbering.
      procedure, pass(self) :: children                !< Return the children list given Morton code.
      procedure, pass(self) :: finest_at_level         !< Return the finest node code at given level.
      procedure, pass(self) :: first_at_level          !< Return the first node code at given level.
      procedure, pass(self) :: first_common_parent     !< Return the first common parent given two codes.
      procedure, pass(self) :: last_at_level           !< Return the last node code at given level.
      procedure, pass(self) :: level                   !< Return the refinement level given the code.
      procedure, pass(self) :: lower                   !< Return true if code is lower than other.
      procedure, pass(self) :: get_neighbor            !< Return the neighbor in a given face of given Morton code.
      procedure, pass(self) :: greater                 !< Return true if code is greater than other.
      procedure, pass(self) :: parent                  !< Return the parent given Morton code.
      procedure, pass(self) :: parent_at_level         !< Return the parent given Morton code at given level.
      procedure, pass(self) :: path                    !< Return the path codes, the list of codes from given node to root.
      procedure, pass(self) :: print_code_topology     !< Print all code topology data.
      procedure, pass(self) :: siblings                !< Return the siblings Morton code given Morton code.
      ! private methods
      procedure, pass(self), private :: coordinates3D_to_morton !< Return the Morton code given ijkl coordinates.
      procedure, pass(self), private :: coordinates2D_to_morton !< Return the Morton code given ijl coordinates.
      procedure, pass(self), private :: morton_to_coordinates3D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: morton_to_coordinates2D !< Return the ijkl coordinates given Morton code.
      procedure, pass(self), private :: update_to_refine        !< Update list of nodes to be refined.
endtype tree_object

contains
   ! public methods
   subroutine add_node(self, code, content, finest_code, refinement_needed, max_load, nodes_number, buckets_number, ratio)
   !< Add a node pointer to the tree.
   !<
   !< @note If a node with the same key is already in the tree, it is removed and the new one will replace it.
   class(tree_object), intent(inout)        :: self              !< The tree.
   integer(I8P),       intent(in)           :: code              !< The Morton code.
   integer(I8P),       intent(in)           :: content           !< The content.
   integer(I8P),       intent(in), optional :: finest_code       !< The finest Morton code.
   integer(I4P),       intent(in), optional :: refinement_needed !< Flag for refinement/derefinement algorithm.
   real(R8P),          intent(in), optional :: max_load          !< Maximum load of tree buckets.
   integer(I4P),       intent(in), optional :: nodes_number      !< Nodes number to be stored in the tree.
   integer(I4P),       intent(in), optional :: buckets_number    !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in), optional :: ratio             !< Refinement ratio.
   integer(I4P)                             :: b                 !< Bucket index, namely hashed key.

   if (.not.self%is_initialized_) &
      call self%initialize(max_load=max_load, nodes_number=nodes_number, buckets_number=buckets_number, ratio=ratio)
   b = self%hash(code=code)
   call self%bucket(b)%add_node(code=code, content=content, finest_code=finest_code, refinement_needed=refinement_needed)
   self%nodes_number = self%nodes_number + 1
   self%code(1:2, b) = self%bucket(b)%code
   endsubroutine add_node

   function codes(self)
   !< Return the list of (sorted) codes actually stored in the tree.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P), allocatable      :: codes(:) !< List of codes.
   integer(I8P)                   :: code     !< Counter.
   integer(I8P)                   :: c        !< Counter.
   integer(I8P), allocatable      :: work(:)  !< Working memory for sorting codes list.

   allocate(codes(self%nodes_number))
   allocate(work((self%nodes_number+1)/2))
   c = 0
   do while(self%loop(code=code))
      c = c + 1
      codes(c) = code
   enddo
   call mergesort(array=codes)

   contains
      recursive subroutine mergesort(array)
      !< Sort input array by means of mergesort algorithm.
      integer(I8P), intent(inout) :: array(:) !< Array to be sorted.
      integer(I8P)                :: half     !< Half size counter.

      half = (size(array) + 1) / 2
      if (size(array) < 2) then
         continue
      else if (size(array) == 2) then
         ! if (array(1) > array(2)) call swap_element(array(1), array(2))
         if (self%greater(array(1), array(2))) call swap_element(array(1), array(2))
      else
         call mergesort(array( : half))
         call mergesort(array(half + 1 :))
         ! if (array(half) > array(half + 1)) then
         if (self%greater(array(half), array(half + 1))) then
            work(1 : half) = array(1 : half)
            call merge_array(work(1 : half), array(half + 1:), array)
         endif
      endif
      endsubroutine mergesort

      subroutine merge_array(A, B, C)
      !< Merge arrays A/B in C.
      integer(I8P), target, intent(in)    :: A(:), B(:) !< Input arrays.
      integer(I8P), target, intent(inout) :: C(:)       !< Output array.
      integer(I8P)                        :: i, j, k    !< Counter.

      if (size(A) + size(B) > size(C)) stop

      i = 1 ; j = 1
      do k = 1, size(C)
         if (i <= size(A) .and. j <= size(B)) then
            ! if (A(i) <= B(j)) then
            if (self%lower(A(i), B(j))) then
               C(k) = A(i)
               i = i + 1
            else
               C(k) = B(j)
               j = j + 1
            endif
         else if (i <= size(A)) then
            C(k) = A(i)
            i = i + 1
         else if (j <= size(B)) then
            C(k) = B(j)
            j = j + 1
         endif
      enddo
      endsubroutine merge_array

      subroutine swap_element(x, y)
      !< Swap array element.
      integer(I8P), intent(inout) :: x, y !< Array element to be swaped.
      integer(I8P)                :: tmp  !< Temporary memory.

      tmp = x ; x = y ; y = tmp
      endsubroutine swap_element
   endfunction codes

   subroutine destroy(self)
   !< Destroy the tree.
   class(tree_object), intent(inout) :: self !< The tree.
   integer(I4P)                      :: b    !< Counter.

   if (allocated(self%bucket)) then
      do b=1, size(self%bucket, dim=1)
        call self%bucket(b)%destroy
      enddo
      deallocate(self%bucket)
   endif
   if (allocated(self%code)) deallocate(self%code)
   self%buckets_number = 0_I4P
   self%nodes_number = 0_I4P
   self%max_load = TREE_MAX_LOAD
   self%ratio = 8_I4P
   self%max_level = 12_I4P
   self%is_initialized_ = .false.
   if (allocated(self%to_refine)) deallocate(self%to_refine)
   endsubroutine destroy

   function loop(self, code, content, node) result(again)
   !< Sentinel while-loop on nodes returning the code/content pair (for tree looping).
   class(tree_object),     intent(in)                      :: self      !< The tree bucket.
   integer(I8P),           intent(out), optional           :: code      !< The Morton code.
   integer(I8P),           intent(out), optional           :: content   !< The content.
   type(tree_node_object), intent(out), optional, pointer  :: node      !< Pointer to current node.
   logical                                                 :: again     !< Sentinel flag to contine the loop.
   integer(I4P), save                                      :: b=1_I4P   !< Bucket counter.
   type(tree_node_object), pointer, save                   :: p=>null() !< Pointer to current node.

   again = .false.
   if (self%nodes_number>0) then
      do
         if (b>self%buckets_number) then
            b = 1
            p => null()
            again = .false.
            return
         else
            if (self%bucket(b)%nodes_number>0) then
               if (.not.associated(p)) then
                  p => self%bucket(b)%head
                  if (present(code)) code = p%code
                  if (present(content)) content = p%content
                  if (present(node)) node => p
                  again = .true.
                  return
               elseif (associated(p%next)) then
                  p => p%next
                  if (present(code)) code = p%code
                  if (present(content)) content = p%content
                  if (present(node)) node => p
                  again = .true.
                  return
               else
                  b = b + 1
                  p => null()
               endif
            else
               b = b + 1
               p => null()
            endif
         endif
      enddo
   endif
   endfunction loop

   function has_code(self, code)
   !< Check if the key is present in the tree.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P),       intent(in) :: code     !< The Morton code.
   logical                        :: has_code !< Check result.

   has_code = .false.
   if (self%is_initialized_) has_code = self%bucket(self%hash(code=code))%has_code(code=code)
   endfunction has_code

   elemental function hash(self, code) result(bucket)
   !< Hash the key.
   class(tree_object), intent(in) :: self   !< The tree.
   integer(I8P),       intent(in) :: code   !< The Morton code.
   integer(I4P)                   :: bucket !< Bucket index corresponding to the key.

   bucket = 0
   ! if (self%is_initialized_) bucket = abs(mod(murmurhash3(key=key), self%buckets_number)) + 1
   if (self%is_initialized_) bucket = int(mod(code, int(self%buckets_number, I8P)), I4P) + 2
   endfunction hash

   subroutine initialize(self, max_load, nodes_number, buckets_number, ratio, max_level)
   !< Initialize the tree.
   class(tree_object), intent(inout)        :: self           !< The tree.
   real(R8P),          intent(in), optional :: max_load       !< Maximum load of tree buckets.
   integer(I4P),       intent(in), optional :: nodes_number   !< Nodes number to be stored in the tree.
   integer(I4P),       intent(in), optional :: buckets_number !< Number of buckets for initialize the tree.
   integer(I4P),       intent(in), optional :: ratio          !< Refinement ratio.
   integer(I4P),       intent(in), optional :: max_level      !< Maximum refinement level.

   call self%destroy
   if (present(max_load)) self%max_load = max_load
   if (present(nodes_number)) then
      self%buckets_number = self%prime_buckets_number(nodes_number=nodes_number)
   else
      self%buckets_number = TREE_BUCKETS_NUMBER_DEF ; if (present(buckets_number)) self%buckets_number = buckets_number
   endif
   allocate(self%bucket(1:self%buckets_number))
   allocate(self%code(1:2,1:self%buckets_number))
   self%code = 0_I8P
   if (present(ratio)) self%ratio = ratio
   if (present(max_level)) self%max_level = max_level
   self%is_initialized_ = .true.
   endsubroutine initialize

   function node(self, code) result(p)
   !< Return a pointer to a node in the tree.
   class(tree_object), intent(in)  :: self !< The tree.
   integer(I8P),       intent(in)  :: code !< The Morton code.
   type(tree_node_object), pointer :: p    !< Pointer to node queried.

   p => null()
   if (self%is_initialized_) p => self%bucket(self%hash(code=code))%node(code=code)
   endfunction node

   function node_content(self, code) result(content)
   !< Return node's content, given the code.
   class(tree_object), intent(in)  :: self    !< The tree.
   integer(I8P),       intent(in)  :: code    !< The Morton code.
   integer(I8P)                    :: content !< Content pointer of the queried node.
   type(tree_node_object), pointer :: p       !< Pointer to current node.

   content = 0
   if (self%is_initialized_) then
      p => self%bucket(self%hash(code=code))%node(code=code)
      if (associated(p)) content = p%content
   endif
   endfunction node_content

   elemental function prime_buckets_number(self, nodes_number) result(buckets_number)
   !< Return the buckets number as the nearest prime number given nodes number.
   !<
   !< @note The balanced buckets number is computing considering the tree load defined in `self` and using the
   !< Sieve of Eratoshenes for findining the nearest prime number.
   class(tree_object), intent(in) :: self           !< The tree.
   integer(I4P),       intent(in) :: nodes_number   !< Nodes number to be stored in the tree.
   integer(I4P)                   :: buckets_number !< Well balanced, prime buckets number.
   logical, allocatable           :: is_prime(:)    !< List of prime numbers up to buckets number.
   integer(I4P)                   :: b              !< Counter.

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

   subroutine refine(self, force_all)
   !< Refine nodes.
   class(tree_object), intent(inout)        :: self           !< The tree.
   logical,            intent(in), optional :: force_all      !< Force all nodes to be refined.
   integer(I8P)                             :: refined_number !< Number of nodes to be refined.
   integer(I8P)                             :: n              !< Counter.
   integer(I4P)                             :: i              !< Counter.

   call self%update_to_refine(refined_number=refined_number, force_all=force_all)
   do n=1, refined_number
      do i=0, self%ratio-1
         call self%add_node(code=self%child(code=self%to_refine(n), i=i), content=n)
      enddo
      call self%remove_node(code=self%to_refine(n))
   enddo
   endsubroutine refine

   subroutine remove_node(self, code)
   !< Remove a node from the tree, given the code.
   class(tree_object), intent(inout) :: self !< The tree.
   integer(I8P),       intent(in)    :: code !< The Morton code.
   integer(I4P)                      :: b    !< Bucket index, namely hashed key.

   if (self%is_initialized_) then
      b = self%hash(code=code)
      call self%bucket(b)%remove_node(code=code)
      self%nodes_number = self%nodes_number - 1
      self%code(1:2, b) = self%bucket(b)%code
   endif
   endsubroutine remove_node

   subroutine resize(self, nodes_number, max_load, ratio)
   !< Resize the tree.
   class(tree_object), intent(inout)        :: self         !< The tree.
   integer(I4P),       intent(in)           :: nodes_number !< Nodes number to be stored in the tree.
   real(R8P),          intent(in), optional :: max_load     !< Maximum load of tree buckets.
   integer(I4P),       intent(in), optional :: ratio        !< Refinement ratio.
   type(tree_object)                        :: swap         !< Temporary (swap) tree.
   integer(I8P)                             :: code         !< The Morton code.
   integer(I8P)                             :: content      !< Tree node content.
   integer(I4P)                             :: b            !< Counter.

   if (self%is_initialized_) then
      if (present(max_load)) self%max_load = max_load
      if (self%nodes_number > int((1._R8P/self%max_load)*nodes_number, I4P)) return ! new size too small, cannot previous nodes
      call swap%initialize(max_load=self%max_load, nodes_number=nodes_number)
      do b=1, self%buckets_number
         do while(self%bucket(b)%loop(code=code, content=content))
            call swap%add_node(code=code, content=content)
         enddo
      enddo
      call move_alloc(from=swap%bucket, to=self%bucket)
      call move_alloc(from=swap%code, to=self%code)
      self%buckets_number  = swap%buckets_number
      self%nodes_number    = swap%nodes_number
      self%is_initialized_ = swap%is_initialized_
   else
      call self%initialize(nodes_number=nodes_number, max_load=max_load, ratio=ratio)
   endif
   endsubroutine resize

   subroutine sanitize(self, to_derefine, iterations_number)
   !< Sanitize the tree.
   class(tree_object),        intent(inout)        :: self                 !< The tree.
   integer(I8P), allocatable, intent(out)          :: to_derefine(:)       !< List of nodes to be derefined.
   integer(I4P),              intent(in), optional :: iterations_number    !< Sanitazie iterations number.
   integer(I4P)                                    :: iterations_number_   !< Sanitazie iterations number.
   type(tree_node_object), pointer                 :: node                 !< Pointer to node.
   type(tree_node_object), pointer                 :: sibling              !< Pointer to node sibling.
   integer(I8P)                                    :: code                 !< Code.
   integer(I8P), allocatable                       :: siblings(:)          !< List of code siblings.
   integer(I8P), allocatable                       :: neighbor(:)          !< List of code neighbors.
   type(tree_node_object), pointer                 :: neigh                !< Pointer to node neighbor.
   integer(I4P)                                    :: neighbor_type        !< Neighbors type.
   logical                                         :: is_sanitize_complete !< Flag for finishing sanitize.
   logical                                         :: can_be_derefined     !< Flag for checking derefinement possibility.
   integer(I8P), allocatable                       :: codes_analyzed(:)    !< List of codes analyzed.
   integer(I4P)                                    :: new_level            !< New level counter.
   integer(I4P)                                    :: new_level_n          !< Neighbor new level counter.
   integer(I4P)                                    :: s, sib, f, n         !< Counter.

   iterations_number_ = TREE_MAX_SANITIZE_ITERATIONS ; if (present(iterations_number)) iterations_number_ = iterations_number

   !< @TODO remove check on max level in update_to_refine
   min_max_check_loop : do while(self%loop(node=node))
      new_level = self%level(code=node%code) + node%refinement_needed
      if ((new_level > self%max_level).or.(new_level < 0)) then
         node%refinement_needed = TO_NOT_TOUCH
      endif
   enddo min_max_check_loop

   sanitize_loop : do s=1, iterations_number
      is_sanitize_complete = .true.

      ! check for the sanity of derefinement
      to_derefine = [-2_I8P] ! initialize list with no-sense code
      codes_analyzed = [-2_I8P] ! initialize list with no-sense code
      derefine_loop : do while(self%loop(node=node))
         ! check if I want to be derefined and I have not been analyzed yet
         if ((node%refinement_needed == TO_BE_DEREFINED).and.(findloc(codes_analyzed, node%code, dim=1)==0)) then
            ! check sibling for derefinement
            can_be_derefined = .true.
            code = node%code
            siblings = self%siblings(code=code)
            sibs_check_loop : do sib=1, self%ratio -1
               if (.not.self%has_code(code=siblings(sib))) then
                  can_be_derefined = .false.
                  exit sibs_check_loop
               endif
               sibling => self%node(code=siblings(sib))
               if (sibling%refinement_needed /= TO_BE_DEREFINED) then
                  can_be_derefined = .false.
                  exit sibs_check_loop
               endif
            enddo sibs_check_loop
            if (can_be_derefined) then
               to_derefine = [to_derefine, [code], siblings]
               codes_analyzed = [codes_analyzed, [code], siblings]
            else
               is_sanitize_complete = .false.
               node%refinement_needed = TO_NOT_TOUCH
               codes_analyzed = [codes_analyzed, [code]] ! check with Francesco
               do sib=1, self%ratio -1
                  if (self%has_code(code=siblings(sib))) then
                     codes_analyzed = [codes_analyzed, [siblings(sib)]] ! check with Francesco
                     sibling => self%node(code=siblings(sib))
                     if (sibling%refinement_needed == TO_BE_DEREFINED) then
                        ! due some of your siblings you cannot be derefined, you need to be altered
                        sibling%refinement_needed = TO_NOT_TOUCH
                     endif
                  endif
               enddo
            endif
         endif
      enddo derefine_loop

     ! check for the sanity of refinement
     do while(self%loop(node=node))
        new_level = self%level(code=node%code) + node%refinement_needed
        do f=1, 6
           call self%get_neighbor(code=node%code, face=f, neighbor=neighbor, neighbor_type=neighbor_type)
           if (neighbor_type /= BOUNDARY_CONDITION_NODE) then
              if (neighbor_type == MORE_REFINED_NODE) then
                 do n=1, size(neighbor, dim=1)
                    ! check level
                    neigh => self%node(code=neighbor(n))
                    new_level_n = self%level(code=neighbor(n)) + neigh%refinement_needed
                    if (new_level_n > new_level + 1) then
                       ! a neighbour want to be refined 2 levels more than me, I have to refine more too
                       is_sanitize_complete = .false.
                       node%refinement_needed = min(node%refinement_needed+1, 1)
                       exit
                    endif
                 enddo
              endif
           endif
        enddo

        if (node%refinement_needed > 1) then
           print '(A)',  'CANNOT REFINE TWICE IN A ROW. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           stop
        endif

        new_level = self%level(code=node%code) + node%refinement_needed
        if (new_level > self%max_level) then
           print '(A)',  'CANNOT REFINE MORE. SOMETHING WENT TERRIBLY WRONG. EXIT!'
           stop
        endif
     enddo
     if (is_sanitize_complete) exit sanitize_loop
   enddo sanitize_loop

   if (.not.is_sanitize_complete) then
      print '(A)',  'SANITZE CANNOT BE COMPLETED. SOMETHING WENT TERRIBLY WRONG. EXIT!'
      stop
   endif
   endsubroutine sanitize

   subroutine traverse(self, iterator)
   !< Traverse tree calling the iterator procedure.
   class(tree_object), intent(in) :: self     !< The tree.
   procedure(iterator_interface)  :: iterator !< The (key) iterator procedure to call for each node.
   integer(I4P)                   :: b        !< Counter.

   if (self%is_initialized_) then
      do b=1, self%buckets_number
         call self%bucket(b)%traverse(iterator)
      enddo
   endif
   endsubroutine traverse

   ! Morton ordering methods
   elemental function child(self, code, i)
   !< Return the i-th child given Morton code.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I4P),       intent(in) :: i     !< Child index [0, ratio-1].
   integer(I8P)                   :: child !< Child Morton code.

   child = self%ratio * code + self%ratio + i
   endfunction child

   elemental function child_local(self, code) result(child)
   !< Return the child index in the local numbering.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I8P)                   :: child !< Child Morton code.

   child = 0
   if (code>0) child = int(code + self%ratio - ((code + self%ratio)/self%ratio)*self%ratio, I4P)
   endfunction child_local

   pure function children(self, code)
   !< Return the children given Morton code.
   class(tree_object), intent(in) :: self        !< The tree.
   integer(I8P),       intent(in) :: code        !< Morton code.
   integer(I8P), allocatable      :: children(:) !< Children Morton code.

   select case(self%ratio)
   case(4_I4P)
      children = [self%child(code=code, i=0), &
                  self%child(code=code, i=1), &
                  self%child(code=code, i=2), &
                  self%child(code=code, i=3)]
   case(8_I4P)
      children = [self%child(code=code, i=0), &
                  self%child(code=code, i=1), &
                  self%child(code=code, i=2), &
                  self%child(code=code, i=3), &
                  self%child(code=code, i=4), &
                  self%child(code=code, i=5), &
                  self%child(code=code, i=6), &
                  self%child(code=code, i=7)]
   endselect
   endfunction children

   elemental function finest_at_level(self, code, level) result(finest)
   !< Return the inest node code at given level, namely the last child at a given level.
   class(tree_object), intent(in) :: self       !< The tree.
   integer(I8P),       intent(in) :: code       !< Morton code.
   integer(I4P),       intent(in) :: level      !< Refinement level.
   integer(I8P)                   :: finest     !< Finest Morton code.
   integer(I4P)                   :: code_level !< Level of the given code.
   integer(I4P)                   :: l          !< Counter.

   finest = code
   code_level = self%level(code=code)
   if (code_level<level) then
      finest = code
      do l=code_level+1, level
         finest = self%child(code=finest, i=self%ratio-1)
      enddo
   endif
   endfunction finest_at_level

   elemental function first_at_level(self, level) result(code)
   !< Return the first node code at given level.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: level !< Refinement level.
   integer(I8P)                   :: code  !< Morton code.
   integer(I4P)                   :: l     !< Counter.

   code = 0
   if (level>1) then
      do l=2, level
         code = self%child(code=code, i=0_I4P)
      enddo
   endif
   endfunction first_at_level

   elemental function first_common_parent(self, code1, code2) result(fc_parent)
   !< Return the first common parent given two codes.
   class(tree_object), intent(in) :: self         !< The tree.
   integer(I8P),       intent(in) :: code1, code2 !< Morton codes.
   integer(I8P)                   :: fc_parent    !< First common parent.
   integer(I4P)                   :: level(2)     !< Levels of two codes.
   integer(I8P)                   :: parent(2)    !< Parents of two codes.
   logical                        :: is_found     !< Flag for searching exit.

   parent(1) = self%parent(code=code1)
   parent(2) = self%parent(code=code2)
   level(1) = self%level(code=parent(1))
   level(2) = self%level(code=parent(2))
   if     (level(1)>level(2)) then
      do while(self%level(code=parent(1))==level(2))
         parent(1) = self%parent(code=parent(1))
      enddo
   elseif (level(1)<level(2)) then
      do while(self%level(code=parent(2))==level(1))
         parent(2) = self%parent(code=parent(2))
      enddo
   endif
   is_found = parent(1) == parent(2)
   do while(.not.is_found)
      parent(1) = self%parent(code=parent(1))
      parent(2) = self%parent(code=parent(2))
      is_found = parent(1) == parent(2)
   enddo
   fc_parent = parent(1)
   endfunction first_common_parent

   elemental function last_at_level(self, level) result(code)
   !< Return the last node code at given level.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: level !< Refinement level.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=level) + self%ratio**level - 1
   endfunction last_at_level

   elemental function level(self, code)
   !< Return the refinement level given the code.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I8P),       intent(in) :: code  !< Morton code.
   integer(I4P)                   :: level !< Refinement level.
   integer(I8P)                   :: c     !< Counter.

   if (code==-1) then
      level = 0 ! ancestor of all has level 0
   else
      level = 1
      c = code
      do while(c>=self%ratio)
         c = (c - self%ratio) / self%ratio
         if (c>=0) level = level + 1
      enddo
   endif
   endfunction level

   elemental function lower(self, lhs, rhs) result(res)
   !< Return true if code is lower than other.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I8P),       intent(in) :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in) :: rhs  !< Right hand side of code comparison.
   logical                        :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) < self%finest_at_level(code=rhs, level=self%max_level)
   endfunction lower

   subroutine get_neighbor(self, code, face, neighbor, neighbor_type)
   !< Return the neighbor in a given face of given Morton code.
   !<
   !< We define *direct neighbor* the neighbor of given code in the given face at the same level of the given code
   !< either if it exists or not.
   class(tree_object), intent(in)               :: self                   !< The tree.
   integer(I8P),       intent(in)               :: code                   !< Morton code.
   integer(I4P),       intent(in)               :: face                   !< Face queried.
   integer(I8P),       intent(out), allocatable :: neighbor(:)            !< Neighbors codes list, [1] or [ratio/2].
   integer(I4P),       intent(out)              :: neighbor_type          !< Type of neighbor.
   integer(I8P)                                 :: direct_neighbor        !< Morton code of direct neighbor.
   integer(I8P)                                 :: direct_neighbor_parent !< Morton code of direct neighbor parent.
   integer(I4P)                                 :: i_dn(4)                !< I coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: j_dn(4)                !< J coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: k_dn(4)                !< K coordinate of direct neighbor, or its 4 children.
   integer(I4P)                                 :: l_dn                   !< L coordinate of direct neighbor.
   integer(I4P)                                 :: i_dn_offset            !< I coordinate offset of direct neighbor l+1.
   integer(I4P)                                 :: j_dn_offset            !< J coordinate offset of direct neighbor l+1.
   integer(I4P)                                 :: k_dn_offset            !< K coordinate offset of direct neighbor l+1.

   ! compute coordinates of code
   select case(self%ratio)
   case(4_I4P)
      call self%morton_to_coordinates(code=code, i=i_dn(1), j=j_dn(1), l=l_dn)
   case(8_I4P)
      call self%morton_to_coordinates(code=code, i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn)
   endselect

   ! compute coordinates of direct neighbor and check if it falls outside the ancestor, in case
   ! it is a boundary condition node
   select case(face)
   case(1_I4P)
      i_dn(1) = i_dn(1) - 1
      if (i_dn(1) < 0) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   case(2_I4P)
      i_dn(1) = i_dn(1) + 1
      if (i_dn(1) > 2**l_dn - 1) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   case(3_I4P)
      j_dn(1) = j_dn(1) - 1
      if (j_dn(1) < 0) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   case(4_I4P)
      j_dn(1) = j_dn(1) + 1
      if (j_dn(1) > 2**l_dn - 1) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   case(5_I4P)
      k_dn(1) = k_dn(1) - 1
      if (k_dn(1) < 0) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   case(6_I4P)
      k_dn(1) = k_dn(1) + 1
      if (k_dn(1) > 2**l_dn - 1) then
         neighbor_type = BOUNDARY_CONDITION_NODE
         return
      endif
   endselect

   ! compute direct neighbor code
   select case(self%ratio)
   case(4_I4P)
      direct_neighbor = self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), l=l_dn)
   case(8_I4P)
      direct_neighbor = self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn)
   endselect

   ! check if direct neighbor is a sibling
   if (findloc(array=self%siblings(code=code), value=direct_neighbor, dim=1)>0) then
      ! direct neighbor is a sibling, thus surely exist and it is unique, return it
      neighbor = [direct_neighbor]
      neighbor_type = STANDARD_NODE
      return
   endif

   ! direct neighbor is not a sibling, check if it exists
   if (self%has_code(code=direct_neighbor)) then
      neighbor = [direct_neighbor]
      neighbor_type = STANDARD_NODE
      return
   endif

   ! direct neighbor does not exist, check if its parent exists
   direct_neighbor_parent = self%parent(code=direct_neighbor)
   if (self%has_code(code=direct_neighbor_parent)) then
      neighbor = [direct_neighbor_parent]
      neighbor_type = STANDARD_NODE
      return
   endif

   ! direct neighbor parent does not exist, its children must exist or 2-1 rule has been broken
   ! using ijk coordinates at level l one can find the ijk coordinates of neighbor at level l+1
   l_dn = l_dn + 1
   i_dn_offset = (i_dn(1) - 1) * 2 + 1
   j_dn_offset = (j_dn(1) - 1) * 2 + 1
   k_dn_offset = (k_dn(1) - 1) * 2 + 1
   select case(face)
   case(1_I4P)
      i_dn(1) = i_dn_offset + 2
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 2
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(2_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 1
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 1
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(3_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 2
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 2
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(4_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 1
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 1
      k_dn(4) = k_dn_offset + 2
   case(5_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 2

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 2

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 2

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 2
   case(6_I4P)
      i_dn(1) = i_dn_offset + 1
      j_dn(1) = j_dn_offset + 1
      k_dn(1) = k_dn_offset + 1

      i_dn(2) = i_dn_offset + 2
      j_dn(2) = j_dn_offset + 1
      k_dn(2) = k_dn_offset + 1

      i_dn(3) = i_dn_offset + 1
      j_dn(3) = j_dn_offset + 2
      k_dn(3) = k_dn_offset + 1

      i_dn(4) = i_dn_offset + 2
      j_dn(4) = j_dn_offset + 2
      k_dn(4) = k_dn_offset + 1
   endselect
   select case(self%ratio)
   case(4_I4P)
      neighbor = [self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(2), j=j_dn(2), l=l_dn)]
      neighbor_type = MORE_REFINED_NODE
   case(8_I4P)
      neighbor = [self%coordinates_to_morton(i=i_dn(1), j=j_dn(1), k=k_dn(1), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(2), j=j_dn(2), k=k_dn(2), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(3), j=j_dn(3), k=k_dn(3), l=l_dn), &
                  self%coordinates_to_morton(i=i_dn(4), j=j_dn(4), k=k_dn(4), l=l_dn)]
      neighbor_type = MORE_REFINED_NODE
   endselect
   endsubroutine get_neighbor

   elemental function greater(self, lhs, rhs) result(res)
   !< Return true if code is greater than other.
   class(tree_object), intent(in) :: self !< The tree.
   integer(I8P),       intent(in) :: lhs  !< Left hand side of code comparison.
   integer(I8P),       intent(in) :: rhs  !< Right hand side of code comparison.
   logical                        :: res  !< Comparison result.

   res = self%finest_at_level(code=lhs, level=self%max_level) > self%finest_at_level(code=rhs, level=self%max_level)
   endfunction greater

   elemental function parent(self, code)
   !< Return the parent given Morton code.
   class(tree_object), intent(in) :: self   !< The tree.
   integer(I8P),       intent(in) :: code   !< Morton code.
   integer(I8P)                   :: parent !< Parent Morton code.

   parent = -1 ! ancestor of all has not parent
   if (code>self%ratio-1) parent = int(real(code - self%ratio) / self%ratio, kind=I8P)
   endfunction parent

   elemental function parent_at_level(self, code, level) result(parent)
   !< Return the parent given Morton code at a given level.
   class(tree_object), intent(in) :: self     !< The tree.
   integer(I8P),       intent(in) :: code     !< Morton code.
   integer(I4P),       intent(in) :: level    !< Refinement level.
   integer(I8P)                   :: parent   !< Parent Morton code.
   integer(I8P), allocatable      :: path_(:) !< Parent Morton code.

   parent = -1_I8P ! ancestor of all nodes
   if (level>=1) then
      path_ = self%path(code=code)
      parent = path_(size(path_) - level + 1)
   endif
   endfunction parent_at_level

   pure function path(self, code)
   !< Return the path codes, the list of codes from given node to root.
   class(tree_object), intent(in) :: self    !< The tree.
   integer(I8P),       intent(in) :: code    !< Morton code.
   integer(I8P), allocatable      :: path(:) !< Path codes from node to root.
   integer(I8P), allocatable      :: temp(:) !< Temporary path list.
   integer(I8P)                   :: c       !< Counter.

   path = [code]
   c = code
   do while(self%level(code=c)>1)
      allocate(temp(1:size(path)+1))
      temp(1:size(path)) = path
      temp(size(path)+1) = self%parent(code=c)
      call move_alloc(from=temp,to=path)
      c = self%parent(code=c)
   enddo
   endfunction path

   subroutine print_code_topology(self, code,  &
                                  level,       &
                                  parent,      &
                                  parents,     &
                                  path,        &
                                  child,       &
                                  child_local, &
                                  finest,      &
                                  siblings,    &
                                  neighbor,    &
                                  whole)
   !< Print all code topology data.
   class(tree_object), intent(in)           :: self          !< The tree.
   integer(I8P),       intent(in)           :: code          !< Morton code.
   logical,            intent(in), optional :: level         !< Level of node.
   logical,            intent(in), optional :: parent        !< Parent of code.
   logical,            intent(in), optional :: parents       !< Parents list.
   logical,            intent(in), optional :: path          !< Path from node to parent of first level.
   logical,            intent(in), optional :: child         !< (First) Child of code.
   logical,            intent(in), optional :: child_local   !< Local child-numbering of code.
   logical,            intent(in), optional :: finest        !< Finest Morton code.
   logical,            intent(in), optional :: siblings      !< Siblings of code.
   logical,            intent(in), optional :: neighbor      !< Neighbor of code.
   logical,            intent(in), optional :: whole         !< Whole topology data.
   logical                                  :: level_        !< Level of node, local var.
   logical                                  :: parent_       !< Parent of code, local var.
   logical                                  :: child_        !< (First) Child of code, local var.
   logical                                  :: child_local_  !< Local child-numbering of code, local var.
   logical                                  :: finest_       !< Finest Morton code, local var.
   logical                                  :: siblings_     !< Siblings of code, local var.
   logical                                  :: path_         !< Path from node to parent of first level, local var.
   logical                                  :: parents_      !< Parents list, local var.
   logical                                  :: neighbor_     !< Neighbor of code, local var.
   logical                                  :: whole_        !< Whole topology data, local var.
   character(:), allocatable                :: topology      !< Topology string.
   character(:), allocatable                :: parents_str   !< Parents string list.
   character(:), allocatable                :: neighbors_str !< Neighbors string.
   integer(I8P), allocatable                :: neighbors(:)  !< Neighbor of code.
   integer(I4P)                             :: neighbor_type !< Type of neighbor.
   integer(I4P)                             :: f, l          !< Counter.

   level_       = .false. ; if (present(level      ))  level_       = level
   parent_      = .false. ; if (present(parent     ))  parent_      = parent
   child_       = .false. ; if (present(child      ))  child_       = child
   child_local_ = .false. ; if (present(child_local))  child_local_ = child_local
   finest_      = .false. ; if (present(finest     ))  finest_      = finest
   siblings_    = .false. ; if (present(siblings   ))  siblings_    = siblings
   path_        = .false. ; if (present(path       ))  path_        = path
   parents_     = .false. ; if (present(parents    ))  parents_     = parents
   neighbor_    = .false. ; if (present(neighbor   ))  neighbor_    = neighbor
   whole_       = .false. ; if (present(whole      ))  whole_       = whole

   parents_str = ''
   do l=self%level(code=code)-1, 0, -1
      parents_str = parents_str//' pal_'//trim(str(l,.true.))//' '//trim(str(self%parent_at_level(code=code, level=l),.true.))
   enddo

   neighbors_str = ''
   do f=1, 6
      if (self%ratio==4_I4P.and.f>4) exit
      call self%get_neighbor(code=code, neighbor=neighbors, face=f, neighbor_type=neighbor_type)
      if (allocated(neighbors)) then
         neighbors_str = neighbors_str//' f_'//trim(str(f,.true.))//' '//trim(str(neighbors,.true.))//&
                         ' type-'//trim(str(neighbor_type,.true.))
      else
         neighbors_str = neighbors_str//' f-'//trim(str(f,.true.))//' type_'//trim(str(neighbor_type,.true.))
      endif
   enddo

   topology = ' code: '//trim(str(code))
   if (level_.or.whole_      ) topology = topology//' level: '//trim(str(self%level(code=code),.true.))
   if (parent_.or.whole_     ) topology = topology//' parent: '//trim(str(self%parent(code=code),.true.))
   if (parents_.or.whole_    ) topology = topology//' parents: '//parents_str
   if (path_.or.whole_       ) topology = topology//' path: '//trim(str(self%path(code=code),.true.))
   if (child_.or.whole_      ) topology = topology//' child: '//trim(str(self%child(code=code, i=0),.true.))
   if (child_local_.or.whole_) topology = topology//' child_local: '//trim(str(self%child_local(code=code),.true.))
   if (finest_.or.whole_     ) topology = topology//' finest: '//trim(str(self%finest_at_level(code=code, &
                                                                                               level=self%max_level),.true.))
   if (siblings_.or.whole_   ) topology = topology//' siblings: '//trim(str(self%siblings(code=code),.true.))
   if (neighbor_.or.whole_   ) topology = topology//' neighbor: '//neighbors_str

   print '(A)', topology
   endsubroutine print_code_topology

   pure function siblings(self, code)
   !< Return the siblings Morton code given Morton code.
   class(tree_object), intent(in) :: self                     !< The tree.
   integer(I8P),       intent(in) :: code                     !< Morton code.
   integer(I8P)                   :: siblings(1:self%ratio-1) !< Siblings Morton codes [1:ratio-1].
   integer(I4P)                   :: local                    !< Local child code [0,ratio].
   integer(I4P)                   :: start                    !<
   integer(I4P)                   :: l, s                     !< Counter.

   if (code==-1) then
      siblings = -1_I8P ! anceestor of all has not siblings
   else
      local = self%child_local(code=code)
      start = code - local + 1
      s = 0
      do l=0, self%ratio - 1
         if (l/=local) then
            s = s + 1
            siblings(s) = start + l - 1
         endif
      enddo
   endif
   endfunction siblings

   ! private methods
   subroutine update_to_refine(self, refined_number, force_all)
   !< List of nodes to be refined.
   class(tree_object), intent(inout)        :: self           !< The tree.
   integer(I8P),       intent(out)          :: refined_number !< Number of nodes to be refined.
   logical,            intent(in), optional :: force_all      !< Force all nodes to be refined.
   logical                                  :: force_all_     !< Force all nodes to be refined, local var.
   type(tree_node_object), pointer          :: node           !< Pointer to current node.
   integer(I8P)                             :: n              !< Counter.

   force_all_ = .false. ; if (present(force_all)) force_all_ = force_all
   if (allocated(self%to_refine)) deallocate(self%to_refine)
   allocate(self%to_refine(1:self%nodes_number))
   self%to_refine = -1_I8P
   n = 0_I8P
   do while(self%loop(node=node))
      if (self%level(code=node%code)+1<=self%max_level) then
         if (node%refinement_needed==TO_BE_REFINED.or.force_all_) then
            n = n + 1
            self%to_refine(n) = node%code
         endif
      endif
   enddo
   refined_number = n
   endsubroutine update_to_refine

   subroutine morton_to_coordinates2D(self, code, i, j, l)
   !< Return the ijkl coordinates given Morton code.
   class(tree_object), intent(in)  :: self     !< The tree.
   integer(I8P),       intent(in)  :: code     !< Morton code.
   integer(I4P),       intent(out) :: i        !< I coordinate.
   integer(I4P),       intent(out) :: j        !< J coordinate.
   integer(I4P),       intent(out) :: l        !< L coordinate.
   integer(I4P)                    :: ij(2)    !< IJ local coordinates.
   integer(I8P), allocatable       :: path_(:) !< Path from code to root.
   integer(I4P)                    :: p        !< Counter.
   integer(I8P)                    :: c        !< Counter.

   i = 0_I4P
   j = 0_I4P
   l = self%level(code=code)
   path_ = self%path(code=code)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton2D(code=self%child_local(code=c), i=ij(1), j=ij(2))
      i = i + ij(1) * 2**(p-1)
      j = j + ij(2) * 2**(p-1)
   enddo
   endsubroutine morton_to_coordinates2D

   subroutine morton_to_coordinates3D(self, code, i, j, k, l)
   !< Return the ijkl coordinates given Morton code.
   class(tree_object), intent(in)  :: self     !< The tree.
   integer(I8P),       intent(in)  :: code     !< Morton code.
   integer(I4P),       intent(out) :: i        !< I coordinate.
   integer(I4P),       intent(out) :: j        !< J coordinate.
   integer(I4P),       intent(out) :: k        !< K coordinate.
   integer(I4P),       intent(out) :: l        !< L coordinate.
   integer(I4P)                    :: ijk(3)   !< IJK local coordinates.
   integer(I8P), allocatable       :: path_(:) !< Path from code to root.
   integer(I4P)                    :: p        !< Counter.
   integer(I8P)                    :: c        !< Counter.

   i = 0_I4P
   j = 0_I4P
   k = 0_I4P
   l = self%level(code=code)
   path_ = self%path(code=code)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton3D(code=self%child_local(code=c), i=ijk(1), j=ijk(2), k=ijk(3))
      i = i + ijk(1) * 2**(p-1)
      j = j + ijk(2) * 2**(p-1)
      k = k + ijk(3) * 2**(p-1)
   enddo
   endsubroutine morton_to_coordinates3D

   function coordinates2D_to_morton(self, i, j, l) result(code)
   !< Return the Morton code given ijl coordinates.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: i     !< I coordinate.
   integer(I4P),       intent(in) :: j     !< J coordinate.
   integer(I4P),       intent(in) :: l     !< L coordinate.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=l) + morton2D(i=i, j=j)
   endfunction coordinates2D_to_morton

   function coordinates3D_to_morton(self, i, j, k, l) result(code)
   !< Return the Morton code given ijkl coordinates.
   class(tree_object), intent(in) :: self  !< The tree.
   integer(I4P),       intent(in) :: i     !< I coordinate.
   integer(I4P),       intent(in) :: j     !< J coordinate.
   integer(I4P),       intent(in) :: k     !< K coordinate.
   integer(I4P),       intent(in) :: l     !< L coordinate.
   integer(I8P)                   :: code  !< Morton code.

   code = self%first_at_level(level=l) + morton3D(i=i, j=j, k=k)
   endfunction coordinates3D_to_morton

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
endmodule adam_tree_object
