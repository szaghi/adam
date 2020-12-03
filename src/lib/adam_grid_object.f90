!< ADAM, grid class definition.
module adam_grid_object
!< ADAM, grid class definition.

use adam_parameters
use PENF

implicit none
private
public :: grid_object

type :: grid_object
   !< Grid class definition.
   real(R8P)    :: domain_emin(3)=[0._R8P,0._R8P,0._R8P] !< Coordinates of minimum abscissa of whole domain.
   real(R8P)    :: domain_emax(3)=[1._R8P,1._R8P,1._R8P] !< Coordinates of maximum abscissa of whole domain.
   integer(I4P) :: ni=16_I4P                             !< Number of cells in i direction.
   integer(I4P) :: nj=16_I4P                             !< Number of cells in j direction.
   integer(I4P) :: nk=16_I4P                             !< Number of cells in k direction.
   integer(I4P) :: gci=2_I4P                             !< Number of ghost cells in i direction for boundary conditions.
   integer(I4P) :: gcj=2_I4P                             !< Number of ghost cells in j direction for boundary conditions.
   integer(I4P) :: gck=2_I4P                             !< Number of ghost cells in k direction for boundary conditions.
   integer(I4P) :: weight_neighbor(26)=0_I4P             !< Weight of neighbors (cells number).
   contains
      ! public methods
      procedure, pass(self) :: compute_metrics         !< Compute metrics of a block.
      procedure, pass(self) :: compute_weight_neighbor !< Compute weight of neighbors.
      procedure, pass(self) :: destroy                 !< Destroy the field.
      procedure, pass(self) :: initialize              !< Initialize the field.
      ! operators
      generic :: assignment(=) => grid_assign_grid      !< Overload `=`.
      procedure, pass(lhs), private :: grid_assign_grid !< Operator `=`.
endtype grid_object

contains
   ! public methods
   subroutine compute_metrics(self, coordinates,      &
                              dx, dy, dz,             &
                              emin, emax,             &
                              x_node, y_node, z_node, &
                              x_cell, y_cell, z_cell)
   !< Compute metrics of a block.
   class(grid_object), intent(in)            :: self                                 !< The grid.
   integer(I4P),       intent(in)            :: coordinates(4)                       !< Block coordinates.
   real(R8P),          intent(out), optional :: dx, dy, dz                           !< Space steps.
   real(R8P),          intent(out), optional :: emin(3), emax(3)                     !< Min/max abscissa of block.
   real(R8P),          intent(out), optional :: x_node(0-self%gci:self%ni+self%gci)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%gcj:self%nj+self%gcj)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%gck:self%nk+self%gck)  !< Z coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%gci:self%ni+self%gci)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%gcj:self%nj+self%gcj)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%gck:self%nk+self%gck)  !< Z coordinates.
   real(R8P)                                 :: emin_(3), emax_(3)                   !< Min/max abscissa of block, local var.
   real(R8P)                                 :: x_node_(0-self%gci:self%ni+self%gci) !< X coordinates, local var.
   real(R8P)                                 :: y_node_(0-self%gcj:self%nj+self%gcj) !< Y coordinates, local var.
   real(R8P)                                 :: z_node_(0-self%gck:self%nk+self%gck) !< Z coordinates, local var.
   real(R8P)                                 :: x_cell_(1-self%gci:self%ni+self%gci) !< X coordinates, local var.
   real(R8P)                                 :: y_cell_(1-self%gcj:self%nj+self%gcj) !< Y coordinates, local var.
   real(R8P)                                 :: z_cell_(1-self%gck:self%nk+self%gck) !< Z coordinates, local var.
   real(R8P)                                 :: dx_, dy_, dz_                        !< Space steps, local var.
   integer(I4P)                              :: i, j, k, l                           !< Counter.

   i = coordinates(1)
   j = coordinates(2)
   k = coordinates(3)
   l = coordinates(4)
   dx_ = (self%domain_emax(1) - self%domain_emin(1)) / 2**l
   dy_ = (self%domain_emax(2) - self%domain_emin(2)) / 2**l
   dz_ = (self%domain_emax(3) - self%domain_emin(3)) / 2**l
   emin_(1) = i * dx_ ; emax_(1) = emin_(1) + dx_
   emin_(2) = j * dy_ ; emax_(2) = emin_(2) + dy_
   emin_(3) = k * dz_ ; emax_(3) = emin_(3) + dz_
   dx_ = dx_ / self%ni
   dy_ = dy_ / self%nj
   dz_ = dz_ / self%nk
   do i=0-self%gci, self%ni+self%gci
      x_node_(i) = emin_(1) + i * dx_
   enddo
   do j=0-self%gcj, self%nj+self%gcj
      y_node_(j) = emin_(2) + j * dy_
   enddo
   do k=0-self%gck, self%nk+self%gck
      z_node_(k) = emin_(3) + k * dz_
   enddo
   do i=1-self%gci, self%ni+self%gci
      x_cell_(i) = x_node_(i-1) + dx_ * 0.5_R8P
   enddo
   do j=1-self%gcj, self%nj+self%gcj
      y_cell_(j) = y_node_(j-1) + dy_ * 0.5_R8P
   enddo
   do k=1-self%gck, self%nk+self%gck
      z_cell_(k) = z_node_(k-1) + dz_ * 0.5_R8P
   enddo
   if (present(dx)) dx = dx_
   if (present(dy)) dy = dy_
   if (present(dz)) dz = dz_
   if (present(emin)) emin = emin_
   if (present(emax)) emax = emax_
   if (present(x_node)) x_node = x_node_
   if (present(y_node)) y_node = y_node_
   if (present(z_node)) z_node = z_node_
   if (present(x_cell)) x_cell = x_cell_
   if (present(y_cell)) y_cell = y_cell_
   if (present(z_cell)) z_cell = z_cell_
   endsubroutine compute_metrics

   subroutine compute_weight_neighbor(self)
   !< Compute weight of neighbors.
   class(grid_object), intent(inout) :: self    !< The grid.
   integer(I4P)                      :: nijk(3) !< Ni, nj, nk in array.
   integer(I4P)                      :: gc(3)   !< Ghost cell in array.
   integer(I4P)                      :: fec     !< Counter.
   integer(I4P)                      :: i       !< Counter.

   self%weight_neighbor = 1_I4P
   nijk = [self%ni, self%nj, self%nk]
   gc   = [self%gci, self%gcj, self%gck]
   do fec=1, 26
      do i=1, 3
         self%weight_neighbor(fec) = self%weight_neighbor(fec) * (  abs(delta_neighbor(i,fec))  * gc(i) + &
                                                                 (1-abs(delta_neighbor(i,fec))) * nijk(i))
      enddo
   enddo
   endsubroutine compute_weight_neighbor

   elemental subroutine destroy(self)
   !< Destroy field.
   class(grid_object), intent(inout) :: self  !< The grid.
   type(grid_object)                 :: fresh !< Fresh grid.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, ni, nj, nk, gc, emin, emax)
   !< Initialize field.
   class(grid_object), intent(inout)        :: self    !< The grid.
   integer(I4P),       intent(in), optional :: ni      !< Number of cells in X direction.
   integer(I4P),       intent(in), optional :: nj      !< Number of cells in Y direction.
   integer(I4P),       intent(in), optional :: nk      !< Number of cells in Z direction.
   integer(I4P),       intent(in), optional :: gc(3)   !< Number of ghost cells in each direction.
   real(R8P),          intent(in), optional :: emin(3) !< Coordinates of minium abscissa.
   real(R8P),          intent(in), optional :: emax(3) !< Coordinates of maxium abscissa.

   call self%destroy
   ! grid data
   if (present(emin)) then
      self%domain_emin = emin
   else
      self%domain_emin = 0._R8P
   endif
   if (present(emax)) then
      self%domain_emax = emax
   else
      self%domain_emax = 1._R8P
   endif
   if (present(ni)) self%ni  = ni
   if (present(nj)) self%nj  = nj
   if (present(nk)) self%nk  = nk
   if (present(gc)) self%gci = gc(1)
   if (present(gc)) self%gcj = gc(2)
   if (present(gc)) self%gck = gc(3)
   call self%compute_weight_neighbor
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine grid_assign_grid(lhs, rhs)
   !< Operator `=`.
   class(grid_object), intent(inout) :: lhs !< Left hand side.
   type(grid_object),  intent(in)    :: rhs !< Right hand side.

   lhs%domain_emin     = rhs%domain_emin
   lhs%domain_emax     = rhs%domain_emax
   lhs%ni              = rhs%ni
   lhs%nj              = rhs%nj
   lhs%nk              = rhs%nk
   lhs%gci             = rhs%gci
   lhs%gcj             = rhs%gcj
   lhs%gck             = rhs%gck
   lhs%weight_neighbor = rhs%weight_neighbor
   endsubroutine grid_assign_grid
endmodule adam_grid_object
