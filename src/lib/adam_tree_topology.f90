!< ADAM, tree topology library.
module adam_tree_topology
!< ADAM, tree topology library: collection of procedures for handling trees topology.
!< The tree is organized as Octree/Quadtree where Morton order is exploited for the linearization of the tree.
!< A prototype of quadtree is represented below.
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
!<    0--------|-----------|-----------O-----------|-----------|
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
!<    0--------|-----------|-----------O-----------|-----------|
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

use MORTIF, only : morton2D, morton3D, demorton2D, demorton3D
use PENF, only : I1P, I4P, I8P, R8P

implicit none
private
public :: morton_to_coordinates2D, morton_to_coordinates3D
public :: coordinates2D_to_morton, coordinates3D_to_morton
public :: child
public :: child_local
public :: first_at_level
public :: last_at_level
public :: level
public :: parent
public :: path
public :: siblings

contains
   ! public procedures
   function morton_to_coordinates2D(code, ratio) result(ijl)
   !< Return the ijkl coordinates given Morton code.
   integer(I8P), intent(in)  :: code     !< Morton code.
   integer(I4P), intent(in)  :: ratio    !< Refinement ratio.
   integer(I4P)              :: ijl(3)   !< IJL coordinates.
   integer(I4P)              :: ij(2)    !< IJ local coordinates.
   integer(I8P), allocatable :: path_(:) !< Path from code to root.
   integer(I4P)              :: p        !< Counter.
   integer(I8P)              :: c        !< Counter.

   ijl = 0_I4P
   ijl(3) = level(code=code, ratio=ratio)
   path_ = path(code=code, ratio=ratio)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton2D(code=child_local(code=c, ratio=ratio), i=ij(1), j=ij(2))
      ijl(1) = ijl(1) + ij(1) * 2**(p-1)
      ijl(2) = ijl(2) + ij(2) * 2**(p-1)
   enddo
   endfunction morton_to_coordinates2D

   function morton_to_coordinates3D(code, ratio) result(ijkl)
   !< Return the ijkl coordinates given Morton code.
   integer(I8P), intent(in)  :: code     !< Morton code.
   integer(I4P), intent(in)  :: ratio    !< Refinement ratio.
   integer(I4P)              :: ijkl(4)  !< IJKL coordinates.
   integer(I4P)              :: ijk(3)   !< IJK local coordinates.
   integer(I8P), allocatable :: path_(:) !< Path from code to root.
   integer(I4P)              :: p        !< Counter.
   integer(I8P)              :: c        !< Counter.

   ijkl = 0_I4P
   ijkl(4) = level(code=code, ratio=ratio)
   path_ = path(code=code, ratio=ratio)
   do p=1, size(path_, dim=1)
      c = path_(p)
      call demorton3D(code=child_local(code=c, ratio=ratio), i=ijk(1), j=ijk(2), k=ijk(3))
      ijkl(1) = ijkl(1) + ijk(1) * 2**(p-1)
      ijkl(2) = ijkl(2) + ijk(2) * 2**(p-1)
      ijkl(3) = ijkl(3) + ijk(3) * 2**(p-1)
   enddo
   endfunction morton_to_coordinates3D

   function coordinates2D_to_morton(ijl, ratio) result(code)
   !< Return the Morton code given ijkl coordinates.
   integer(I4P), intent(in)  :: ijl(3) !< IJL coordinates.
   integer(I4P), intent(in)  :: ratio  !< Refinement ratio.
   integer(I8P)              :: code   !< Morton code.

   code = first_at_level(level=ijl(3), ratio=ratio) + morton2D(i=ijl(1), j=ijl(2))
   endfunction coordinates2D_to_morton

   function coordinates3D_to_morton(ijkl, ratio) result(code)
   !< Return the Morton code given ijkl coordinates.
   integer(I4P), intent(in)  :: ijkl(4) !< IJKL coordinates.
   integer(I4P), intent(in)  :: ratio   !< Refinement ratio.
   integer(I8P)              :: code    !< Morton code.

   code = first_at_level(level=ijkl(4), ratio=ratio) + morton3D(i=ijkl(1), j=ijkl(2), k=ijkl(3))
   endfunction coordinates3D_to_morton

   ! pure/elemental procedures
   elemental function child(code, i, ratio)
   !< Return the i-th child given Morton code.
   integer(I8P), intent(in) :: code  !< Morton code.
   integer(I4P), intent(in) :: i     !< Child index [0, ratio-1].
   integer(I4P), intent(in) :: ratio !< Refinement ratio.
   integer(I8P)             :: child !< Child Morton code.

   child = ratio * code + ratio + i
   endfunction child

   elemental function child_local(code, ratio) result(child)
   !< Return the child index in the local numbering.
   integer(I8P), intent(in) :: code  !< Morton code.
   integer(I4P), intent(in) :: ratio !< Refinement ratio.
   integer(I8P)             :: child !< Child Morton code.

   child = 0
   if (code>0) child = int(code + ratio - ((code + ratio)/ratio)*ratio, I4P)
   endfunction child_local

   elemental function first_at_level(level, ratio) result(code)
   !< Return the first node code at given level.
   integer(I4P), intent(in) :: level !< Refinement level.
   integer(I4P), intent(in) :: ratio !< Refinement ratio.
   integer(I8P)             :: code  !< Morton code.
   integer(I4P)             :: l     !< Counter.

   code = 0
   if (level>1) then
      do l=2, level
         code = child(code=code, i=0_I4P, ratio=ratio)
      enddo
   endif
   endfunction first_at_level

   elemental function last_at_level(level, ratio) result(code)
   !< Return the last node code at given level.
   integer(I4P), intent(in) :: level !< Refinement level.
   integer(I4P), intent(in) :: ratio !< Refinement ratio.
   integer(I8P)             :: code  !< Morton code.

   code = first_at_level(level=level, ratio=ratio) + ratio**level - 1
   endfunction last_at_level

   elemental function level(code, ratio)
   !< Return the refinement level given the code.
   integer(I8P), intent(in) :: code  !< Morton code.
   integer(I4P), intent(in) :: ratio !< Refinement ratio.
   integer(I4P)             :: level !< Refinement level.
   integer(I8P)             :: c     !< Counter.

   level = 1
   c = code
   do while(c>=ratio)
      c = (c - ratio) / ratio
      if (c>=0) level = level + 1
   enddo
   endfunction level

   elemental function parent(code, ratio)
   !< Return the parent given Morton code.
   integer(I8P), intent(in) :: code   !< Morton code.
   integer(I4P), intent(in) :: ratio  !< Refinement ratio.
   integer(I8P)             :: parent !< Parent Morton code.

   parent = 0
   if (code>ratio) parent = int(real(code - ratio) / ratio, kind=I8P)
   endfunction parent

   pure function path(code, ratio)
   integer(I8P), intent(in)  :: code    !< Morton code.
   integer(I4P), intent(in)  :: ratio   !< Refinement ratio.
   integer(I8P), allocatable :: path(:) !< Path codes from node to root.
   integer(I8P), allocatable :: temp(:) !< Temporary path list.
   integer(I8P)              :: c       !< Counter.

   path = [code]
   c = code
   do while(level(code=c, ratio=ratio)>1)
      allocate(temp(1:size(path)+1))
      temp(1:size(path)) = path
      temp(size(path)+1) = parent(code=c, ratio=ratio)
      call move_alloc(from=temp,to=path)
      c = parent(code=c, ratio=ratio)
   enddo
   endfunction path

   pure function siblings(code, ratio)
   !< Return the parent given Morton code.
   integer(I8P), intent(in) :: code                !< Morton code.
   integer(I4P), intent(in) :: ratio               !< Refinement ratio.
   integer(I8P)             :: siblings(1:ratio-1) !< Siblings Morton codes [1:ratio-1].
   integer(I4P)             :: local               !< Local child code [0,ratio].
   integer(I4P)             :: start               !<
   integer(I4P)             :: l, s                !< Counter.

   local = child_local(code=code, ratio=ratio)
   start = code - local + 1
   s = 0
   do l=0, ratio - 1
      if (l/=local) then
         s = s + 1
         siblings(s) = start + l - 1
      endif
   enddo
   endfunction siblings
endmodule adam_tree_topology
