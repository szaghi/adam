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

implicit none
private
public :: field_object

type :: field_object
   integer(I4P)              :: ni=32_I4P           !< Number of cells in i direction.
   integer(I4P)              :: nj=32_I4P           !< Number of cells in j direction.
   integer(I4P)              :: nk=32_I4P           !< Number of cells in k direction.
   integer(I4P)              :: gc1=4_I4P           !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P)              :: gc2=4_I4P           !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P)              :: gc3=4_I4P           !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P)              :: gc4=4_I4P           !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P)              :: gc5=4_I4P           !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P)              :: gc6=4_I4P           !< Number of ghost cells in k+ direction for boundary conditions.
   integer(I4P)              :: nb=0_I4P            !< Number of all blocks that can be stored.
   integer(I4P)              :: nv=1_I4P            !< Number of field variables.
   integer(I4P)              :: blocks_number=0_I4P !< Number of blocks actually stored.
   integer(I8P), allocatable :: code(:)             !< Morton codes [nb].
   real(R8P),    allocatable :: emin(:,:)           !< Coordinates of minimum abscissa of each block [3,nb]
   real(R8P),    allocatable :: emax(:,:)           !< Coordinates of maximum abscissa of each block [3,nb]
   real(R8P),    allocatable :: x(:,:)              !< Grids nodes x abscissa [ni+gc1+gc2+1,nb]
   real(R8P),    allocatable :: y(:,:)              !< Grids nodes y abscissa [nj+gc3+gc4+1,nb]
   real(R8P),    allocatable :: z(:,:)              !< Grids nodes z abscissa [nk+gc5+gc6+1,nb]
   real(R8P),    allocatable :: u(:,:,:,:)          !< Field cell centered variables [ni+gc12,nj+gc34,nk+gc56,nv,nb]
   contains
      ! public methods
      procedure, pass(self) :: adapt       !< Adapt field accordingly to refine/derefine necessity.
      procedure, pass(self) :: compute_xyz !< Compute grids coordinates from grids extents emin/emax.
      procedure, pass(self) :: destroy     !< Destroy the field.
      procedure, pass(self) :: initialize  !< Initialize the field.
      ! private methods
      procedure, pass(self), private :: derefine !< Derefine blocks.
      procedure, pass(self), private :: refine   !< Refine blocks.
endtype field_object

contains
   ! public methods
   subroutine adapt(self, ratio, block_to_refine, block_refined, block_to_derefine, block_derefined)
   !< Adapt field accordingly to refine/derefine necessity.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:)   !< List of field blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:)   !< List of field refined blocks with Morton code.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of field blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of field derefined blocks with Morton code.

   call self%refine(  ratio=ratio, block_to_refine=block_to_refine,     block_refined=block_refined    )
   call self%derefine(ratio=ratio, block_to_derefine=block_to_derefine, block_derefined=block_derefined)
   endsubroutine adapt

   elemental subroutine compute_xyz(self, b)
   !< Compute grids coordinates from grids extents emin/emax of b-th block.
   class(field_object), intent(inout) :: self       !< The field.
   integer(I4P),        intent(in)    :: b          !< Block index.
   real(R8P)                          :: dx, dy, dz !< Space deltas.
   integer(I4P)                       :: i, j, k    !< Counter.

   associate(emin=>self%emin, emax=>self%emax, ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             gc1=>self%gc1, gc2=>self%gc2, gc3=>self%gc3,  gc4=>self%gc4, gc5=>self%gc5, gc6=>self%gc6)
      dx = (emax(1,b) - emin(1,b)) / ni
      dy = (emax(2,b) - emin(2,b)) / nj
      dz = (emax(3,b) - emin(3,b)) / nk
      do i=0-gc1, ni+gc2
         self%x(i,b) = emin(1,b) + i * dx
      enddo
      do j=0-gc3, nj+gc4
         self%y(j,b) = emin(2,b) + j * dy
      enddo
      do k=0-gc5, nk+gc6
         self%z(k,b) = emin(3,b) + k * dz
      enddo
   endassociate
   endsubroutine compute_xyz

   pure subroutine derefine(self, ratio, block_to_derefine, block_derefined)
   !< Derefine blocks.
   class(field_object),       intent(inout) :: self                 !< The field.
   integer(I4P),              intent(in)    :: ratio                !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_derefine(:) !< List of blocks to be derefined.
   integer(I8P), allocatable, intent(in)    :: block_derefined(:,:) !< List of derefined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz           !< Space deltas.
   integer(I4P)                             :: b, i, j, k           !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii           !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4   !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8   !< Counter.

   do b=1, size(block_derefined, dim=2)
      ib = block_derefined(2,b)

      ic1 = block_to_derefine((b-1)*ratio+1)
      ic2 = block_to_derefine((b-1)*ratio+2)
      ic3 = block_to_derefine((b-1)*ratio+3)
      ic4 = block_to_derefine((b-1)*ratio+4)
      ic5 = block_to_derefine((b-1)*ratio+5)
      ic6 = block_to_derefine((b-1)*ratio+6)
      ic7 = block_to_derefine((b-1)*ratio+7)
      ic8 = block_to_derefine((b-1)*ratio+8)

      self%u(:,:,:,ib) = block_derefined(1,b) ; self%code(ib) = block_derefined(1,b)

   enddo

   do b=1, size(block_derefined, dim=2)
      ib = block_derefined(2,b)

      ic1 = block_to_derefine((b-1)*ratio+1)

      dx = self%emax(1,ic1) - self%emin(1,ic1)
      dy = self%emax(2,ic1) - self%emin(2,ic1)
      dz = self%emax(3,ic1) - self%emin(3,ic1)

      self%emin(1,ib) = self%emin(1,ic1)
      self%emin(2,ib) = self%emin(2,ic1)
      self%emin(3,ib) = self%emin(3,ic1)

      self%emax(1,ib) = self%emin(1,ic1) + 2 * dx
      self%emax(2,ib) = self%emin(2,ic1) + 2 * dy
      self%emax(3,ib) = self%emin(3,ic1) + 2 * dz

      call self%compute_xyz(b=ib)

   enddo
   endsubroutine derefine

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
   self%nv  = 1_I4P
   if (allocated(self%code)) deallocate(self%code)
   if (allocated(self%emin)) deallocate(self%emin)
   if (allocated(self%emax)) deallocate(self%emax)
   if (allocated(self%x   )) deallocate(self%x   )
   if (allocated(self%y   )) deallocate(self%y   )
   if (allocated(self%z   )) deallocate(self%z   )
   if (allocated(self%u   )) deallocate(self%u   )
   endsubroutine destroy

   pure subroutine initialize(self, ni, nj, nk, gc, nv, nb, emin, emax)
   !< Initialize field.
   class(field_object), intent(inout)        :: self    !< The field.
   integer(I4P),        intent(in), optional :: ni      !< Number of cells in X direction.
   integer(I4P),        intent(in), optional :: nj      !< Number of cells in Y direction.
   integer(I4P),        intent(in), optional :: nk      !< Number of cells in Z direction.
   integer(I4P),        intent(in), optional :: gc(6)   !< Number of ghost cells in each direction.
   integer(I4P),        intent(in), optional :: nv      !< Number of field variables.
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
   if (present(nv)) self%nv  = nv
   if (present(nb)) self%nb  = nb
   if (nb>0) then
      allocate(self%code(nb))
      self%code    = -2_I8P
      self%code(1) = -1_I8P ! first block is assumed to be ADAM

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

      allocate(self%x(0-self%gc1:self%ni+self%gc2, 1:self%nb))
      allocate(self%y(0-self%gc3:self%nj+self%gc4, 1:self%nb))
      allocate(self%z(0-self%gc5:self%nk+self%gc6, 1:self%nb))
      self%x = 0._R8P
      self%y = 0._R8P
      self%z = 0._R8P

      allocate(self%u(1-self%gc1:self%ni+self%gc2, & ! TODO add nv rank
                      1-self%gc3:self%nj+self%gc4, &
                      1-self%gc5:self%nk+self%gc6, 1:self%nb))
      self%u = 0._R8P

      call self%compute_xyz(b=1) ! compute grid of ADAM
   endif
   endsubroutine initialize

   pure subroutine refine(self, ratio, block_to_refine, block_refined)
   !< Refine blocks.
   class(field_object),       intent(inout) :: self               !< The field.
   integer(I4P),              intent(in)    :: ratio              !< Refinement ratio.
   integer(I8P), allocatable, intent(in)    :: block_to_refine(:) !< List of blocks to be refined.
   integer(I8P), allocatable, intent(in)    :: block_refined(:,:) !< List of refined blocks with Morton code.
   real(R8P)                                :: dx, dy, dz         !< Space deltas.
   integer(I4P)                             :: b, i, j, k         !< Spatial counter.
   integer(I4P)                             :: ib, ic, ii         !< Counter.
   integer(I4P)                             :: ic1, ic2, ic3, ic4 !< Counter.
   integer(I4P)                             :: ic5, ic6, ic7, ic8 !< Counter.

   self%blocks_number = maxval(block_refined(2,:))
   do b=1, size(block_to_refine, dim=1)
      ib = block_to_refine(b)

      ic1 = block_refined(2,(b-1)*ratio+1)
      ic2 = block_refined(2,(b-1)*ratio+2)
      ic3 = block_refined(2,(b-1)*ratio+3)
      ic4 = block_refined(2,(b-1)*ratio+4)
      ic5 = block_refined(2,(b-1)*ratio+5)
      ic6 = block_refined(2,(b-1)*ratio+6)
      ic7 = block_refined(2,(b-1)*ratio+7)
      ic8 = block_refined(2,(b-1)*ratio+8)

      self%u(:,:,:,ic1) = block_refined(1,(b-1)*ratio+1) ; self%code(ic1) = block_refined(1,(b-1)*ratio+1)
      self%u(:,:,:,ic2) = block_refined(1,(b-1)*ratio+2) ; self%code(ic2) = block_refined(1,(b-1)*ratio+2)
      self%u(:,:,:,ic3) = block_refined(1,(b-1)*ratio+3) ; self%code(ic3) = block_refined(1,(b-1)*ratio+3)
      self%u(:,:,:,ic4) = block_refined(1,(b-1)*ratio+4) ; self%code(ic4) = block_refined(1,(b-1)*ratio+4)
      self%u(:,:,:,ic5) = block_refined(1,(b-1)*ratio+5) ; self%code(ic5) = block_refined(1,(b-1)*ratio+5)
      self%u(:,:,:,ic6) = block_refined(1,(b-1)*ratio+6) ; self%code(ic6) = block_refined(1,(b-1)*ratio+6)
      self%u(:,:,:,ic7) = block_refined(1,(b-1)*ratio+7) ; self%code(ic7) = block_refined(1,(b-1)*ratio+7)
      self%u(:,:,:,ic8) = block_refined(1,(b-1)*ratio+8) ; self%code(ic8) = block_refined(1,(b-1)*ratio+8)
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
               ic = block_refined(2,(b-1)*ratio + ii)

               self%emin(1,ic) = self%emin(1,ib) + i * dx/2
               self%emin(2,ic) = self%emin(2,ib) + j * dy/2
               self%emin(3,ic) = self%emin(3,ib) + k * dz/2

               self%emax(1,ic) = self%emin(1,ic) + dx/2
               self%emax(2,ic) = self%emin(2,ic) + dy/2
               self%emax(3,ic) = self%emin(3,ic) + dz/2

               call self%compute_xyz(b=ic)

               ii = ii + 1
            enddo
         enddo
      enddo
   enddo
   endsubroutine refine
endmodule adam_field_object
