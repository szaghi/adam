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

use PENF, only : I8P, I4P, R8P, str
use vtk_fortran, only : vtK_file, vtm_file

implicit none
private
public :: field_object

type :: field_object
   integer(I4P)           :: ni=32_I4P           !< Number of cells in i direction.
   integer(I4P)           :: nj=32_I4P           !< Number of cells in j direction.
   integer(I4P)           :: nk=32_I4P           !< Number of cells in k direction.
   integer(I4P)           :: gc1=4_I4P           !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P)           :: gc2=4_I4P           !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P)           :: gc3=4_I4P           !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P)           :: gc4=4_I4P           !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P)           :: gc5=4_I4P           !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P)           :: gc6=4_I4P           !< Number of ghost cells in k+ direction for boundary conditions.
   integer(I4P)           :: nb=0_I4P            !< Number of all blocks that can be stored.
   integer(I4P)           :: blocks_number=0_I4P !< Number of blocks actually stored.
   real(R8P), allocatable :: emin(:,:)           !< Coordinates of minimum abscissa of each block [3,nb]
   real(R8P), allocatable :: emax(:,:)           !< Coordinates of maximum abscissa of each block [3,nb]
   real(R8P), allocatable :: u(:,:,:,:)          !< Field cell centered variables [ni,nj,nk,nb]
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy the field.
      procedure, pass(self) :: initialize !< Initialize the field.
      procedure, pass(self) :: refine     !< Refine blocks.
      procedure, pass(self) :: save_vtk   !< Save field in VTK files.
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
   class(field_object), intent(inout)        :: self    !< The field.
   integer(I4P),        intent(in), optional :: ni      !< Number of cells in X direction.
   integer(I4P),        intent(in), optional :: nj      !< Number of cells in Y direction.
   integer(I4P),        intent(in), optional :: nk      !< Number of cells in Z direction.
   integer(I4P),        intent(in), optional :: gc(6)   !< Number of ghost cells in each direction.
   integer(I4P),        intent(in), optional :: nb      !< Number of all blocks that can be stored.
   real(R8P),           intent(in), optional :: emin(3) !< Coordinates of minium abscissa.
   real(R8P),           intent(in), optional :: emax(3) !< Coordinates of maxium abscissa.

   call self%destroy
   if (present(ni)) self%ni  = ni
   if (present(nj)) self%nj  = nj
   if (present(nk)) self%nk  = nk
   if (present(gc)) self%gc1 = gc(1)
   if (present(gc)) self%gc2 = gc(2)
   if (present(gc)) self%gc3 = gc(3)
   if (present(gc)) self%gc4 = gc(4)
   if (present(gc)) self%gc5 = gc(5)
   if (present(gc)) self%gc6 = gc(6)
   if (present(nb)) self%nb  = nb
   if (nb>0) then
      allocate(self%emin(3,nb))
      allocate(self%emax(3,nb))
      if (present(emin)) then
         self%emin(:,1) = emin
      else
         self%emin(:,1) = 0._R8P
      endif
      if (present(emax)) then
         self%emax(:,1) = emax
      else
         self%emax(:,1) = 1._R8P
      endif

      allocate(self%u(1-self%gc1:self%ni+self%gc2, &
                      1-self%gc3:self%nj+self%gc4, &
                      1-self%gc5:self%nk+self%gc6, 1:self%nb))
      self%u = 1._R8P
   endif
   endsubroutine initialize

   pure subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
   class(field_object),       intent(inout) :: self               !< The field.
   integer(I4P),              intent(in)    :: ratio              !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:)   !< List of refined blocks.
   real(R8P)                                :: dx, dy, dz         !< Space deltas.
   integer(I4P)                             :: b, i, j, k         !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii         !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4 !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8 !< Counter.

   self%blocks_number = maxval(block_refined, dim=1)
   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)

      ic1 = block_refined((b-1)*ratio+1) ; self%u(:,:,:,ic1) = self%u(:,:,:,ib)
      ic2 = block_refined((b-1)*ratio+2) ; self%u(:,:,:,ic2) = self%u(:,:,:,ib)
      ic3 = block_refined((b-1)*ratio+3) ; self%u(:,:,:,ic3) = self%u(:,:,:,ib)
      ic4 = block_refined((b-1)*ratio+4) ; self%u(:,:,:,ic4) = self%u(:,:,:,ib)
      ic5 = block_refined((b-1)*ratio+5) ; self%u(:,:,:,ic5) = self%u(:,:,:,ib)
      ic6 = block_refined((b-1)*ratio+6) ; self%u(:,:,:,ic6) = self%u(:,:,:,ib)
      ic7 = block_refined((b-1)*ratio+7) ; self%u(:,:,:,ic7) = self%u(:,:,:,ib)
      ic8 = block_refined((b-1)*ratio+8) ; self%u(:,:,:,ic8) = self%u(:,:,:,ib)
   enddo

   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)

      dx = self%emax(1,ib) - self%emin(1,ib)
      dy = self%emax(2,ib) - self%emin(2,ib)
      dz = self%emax(3,ib) - self%emin(3,ib)

      ii = 1
      do k=0, 1
         do j=0, 1
            do i=0, 1
               ic = block_refined((b-1)*ratio + ii)

               self%emin(1,ic) = self%emin(1,ib) + i * dx/2
               self%emin(2,ic) = self%emin(2,ib) + j * dy/2
               self%emin(3,ic) = self%emin(3,ib) + k * dz/2

               self%emax(1,ic) = self%emin(1,ic) + dx/2
               self%emax(2,ic) = self%emin(2,ic) + dy/2
               self%emax(3,ic) = self%emin(3,ic) + dz/2

               ii = ii + 1
            enddo
         enddo
      enddo
   enddo
   endsubroutine refine

   subroutine save_vtk(self, basename)
   !< Save field in VTK files.
   class(field_object), intent(in) :: self      !< The field.
   character(*),        intent(in) :: basename  !< Base name of output files.
   integer(I4P)                    :: error     !< Error trapping flag.
   character(:), allocatable       :: filenames !< File names list.
   type(vtk_file)                  :: vtk       !< VTK file handler.
   type(vtm_file)                  :: vtm       !< VTM file handler.
   integer(I4P)                    :: b         !< Counter.

   filenames = ''
   do b=1, self%blocks_number
      filenames = filenames//trim(basename)//'-block-'//trim(str(b,.true.))//'.vtr '
      error = vtk%initialize(format='raw', filename=trim(basename)//'-block-'//trim(str(b,.true.))//'.vtr', &
                             mesh_topology='RectilinearGrid',                                              &
                             nx1=0, nx2=1, ny1=0, ny2=1, nz1=0, nz2=1)
      error = vtk%xml_writer%write_piece(nx1=0, nx2=1, ny1=0, ny2=1, nz1=0, nz2=1)
      error = vtk%xml_writer%write_geo(x=[self%emin(1,b), self%emax(1,b)], &
                                       y=[self%emin(2,b), self%emax(2,b)], &
                                       z=[self%emin(3,b), self%emax(3,b)])
      error = vtk%xml_writer%write_dataarray(location='cell', action='open')
      error = vtk%xml_writer%write_dataarray(data_name='morton', x=[self%u(1,1,1,b)])
      error = vtk%xml_writer%write_dataarray(location='cell', action='close')
      error = vtk%xml_writer%write_piece()
      error = vtk%finalize()
   enddo

   error = vtm%initialize(filename=trim(basename)//'.vtm')
   error = vtm%write_block(filenames=trim(filenames), name='adam')
   error = vtm%finalize()
   endsubroutine save_vtk
endmodule adam_field_object
