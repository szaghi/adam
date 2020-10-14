!< ADAM, field class definition.
module adam_field_object
!< ADAM, field class definition.

!< A structured block is composed of hexahedron finite volumes with quadrilateral faces using the
!< following internal numeration for nodes and faces:
!<```
!< /|\Z
!<  |                            F(4)         _ F(6)
!<  |                            /|\          /!
!<  |                        7    |          /    8
!<  |                         *------------------*
!<  |                        /|   |        /    /|
!<  |                       / |   |       /    / |
!<  |                      /  |   |      /    /  |
!<  |                     /   |   |     /    /   |
!<  |                    /    |   |    +    /    |
!<  |                   /     |   |        /     |
!<  |                  /      |   +       /      |
!<  |                 /      3|          /       |4
!<  |                /        * --------/--------*
!<  |      F(1)<----/----+   /         /        /
!<  |              *------------------*    +-------->F(2)
!<  |             5|       /          |6      /
!<  |              |      /           |      /
!<  |              |     /        +   |     /
!<  |              |    /         |   |    /
!<  |              |   /      +   |   |   /
!<  |              |  /      /    |   |  /
!<  |              | /      /     |   | /
!<  |              |/      /      |   |/
!<  |              *------------------*
!<  |             1      /        |    2
!<  |                   /        \|/
!<  |   _ Y           |/_       F(3)
!<  |   /|         F(5)
!<  |  /
!<  | /
!<  |/                                                    X
!<  O----------------------------------------------------->
!<```
!< Each hexadron cells is faces-connected to its neighboring, thus the cells build a structured block with implicit
!< connectivity, e.g. in 2D space a FriVolous block could be as the following:
!<```
!<                 _ J
!<                 /|                          _____
!<               5+ ...*----*----*----*----*...     |
!<               /    /    /    /    /    /         |
!<              /    /    /    /    /    /          |
!<            4+ ...*----*----*----*----*...        |
!<            /    /    /    /    /    /            |
!<           /    /    /    /    /    /             |
!<         3+ ...*----*----*----*----*...           |  Structured block of 4x4 Finite Volumes
!<         /    /    / FV /    /    /               |
!<        /    /    /    /    /    /                |
!<      2+ ...*----*----*----*----*...              |
!<      /    /    /    /    /    /                  |
!<     /    /    /    /    /    /                   |
!<   1+ ...*----*----*----*----*...                 |
!<   /     .    .    .    .    .                    |
!<  /      .    .    .    .    .               _____
!< O-------+----+----+----+----+-------------> I
!<         1    2    3    4    5
!<```

use PENF, only : I4P, R8P
! use VecFor, only : vector

implicit none
private
public :: field_object

type :: field_object
   integer(I4P) :: ni=32_I4P            !< Number of cells in X direction.
   integer(I4P) :: nj=32_I4P            !< Number of cells in Y direction.
   integer(I4P) :: nk=32_I4P            !< Number of cells in Z direction.
   integer(I4P) :: gc1=4_I4P            !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P) :: gc2=4_I4P            !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P) :: gc3=4_I4P            !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P) :: gc4=4_I4P            !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P) :: gc5=4_I4P            !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P) :: gc6=4_I4P            !< Number of ghost cells in k+ direction for boundary conditions.
   integer(I4P) :: nb=0_I4P             !< Number of blocks.
   ! type(vector) :: emin                 !< Coordinates of minimum abscissa of the block.
   ! type(vector) :: emax                 !< Coordinates of maximum abscissa of the block.
   real(R8P), allocatable :: emin(:,:)  ! [    3,b]
   real(R8P), allocatable :: emax(:,:)  ! [    3,b]
   real(R8P), allocatable :: u(:,:,:,:) ! [i,j,k,b]
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy the field.
      procedure, pass(self) :: initialize !< Initialize the field.
endtype field_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy field.
   class(field_object), intent(inout) :: self !< The field.

   self%ni  = 32_I4P
   self%nj  = 32_I4P
   self%nk  = 32_I4P
   self%gc1 = 4_I4P
   self%gc2 = 4_I4P
   self%gc3 = 4_I4P
   self%gc4 = 4_I4P
   self%gc5 = 4_I4P
   self%gc6 = 4_I4P
   self%nb  = 0_I4P
   if (allocated(self%emin)) deallocate(self%emin)
   if (allocated(self%emax)) deallocate(self%emax)
   if (allocated(self%u   )) deallocate(self%u   )
   endsubroutine destroy

   pure subroutine initialize(self, ni, nj, nk, gc, nb, emin, emax)
   !< Initialize field.
   class(field_object), intent(inout)        :: self      !< The field.
   integer(I4P),        intent(in), optional :: ni        !< Number of cells in X direction.
   integer(I4P),        intent(in), optional :: nj        !< Number of cells in Y direction.
   integer(I4P),        intent(in), optional :: nk        !< Number of cells in Z direction.
   integer(I4P),        intent(in), optional :: gc(6)     !< Number of ghost cells in each direction.
   integer(I4P),        intent(in), optional :: nb        !< Number of blocks.
   real(R8P),           intent(in), optional :: emin(:,:)
   real(R8P),           intent(in), optional :: emax(:,:)

   call self%destroy
   if (present(ni  )) self%ni   = ni
   if (present(nj  )) self%nj   = nj
   if (present(nk  )) self%nk   = nk
   if (present(gc  )) self%gc1  = gc(1)
   if (present(gc  )) self%gc2  = gc(2)
   if (present(gc  )) self%gc3  = gc(3)
   if (present(gc  )) self%gc4  = gc(4)
   if (present(gc  )) self%gc5  = gc(5)
   if (present(gc  )) self%gc6  = gc(6)
   if (present(nb  )) self%nb   = nb
   if (present(emin)) self%emin = emin
   if (present(emax)) self%emax = emax
   if (nb>0) then
      allocate(self%u(1-self%gc1:self%ni+self%gc2, 1-self%gc3:self%nj+self%gc4, 1-self%gc5:self%nk+self%gc6, 1:self%nb))
      self%u = 1._R8P
   endif
   endsubroutine initialize
endmodule adam_field_object
