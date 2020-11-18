!< ADAM, grid class definition.
module adam_grid_object
!< ADAM, grid class definition.

use PENF

implicit none
private
public :: grid_object

type :: grid_object
   !< Grid class definition.
   real(R8P)    :: domain_emin(3) !< Coordinates of minimum abscissa of whole domain.
   real(R8P)    :: domain_emax(3) !< Coordinates of maximum abscissa of whole domain.
   integer(I4P) :: ni=64_I4P      !< Number of cells in i direction.
   integer(I4P) :: nj=64_I4P      !< Number of cells in j direction.
   integer(I4P) :: nk=64_I4P      !< Number of cells in k direction.
   integer(I4P) :: gc1=0_I4P      !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P) :: gc2=0_I4P      !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P) :: gc3=0_I4P      !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P) :: gc4=0_I4P      !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P) :: gc5=0_I4P      !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P) :: gc6=0_I4P      !< Number of ghost cells in k+ direction for boundary conditions.
   contains
      ! public methods
      procedure, pass(self) :: compute_metrics !< Compute metrics of a block.
      procedure, pass(self) :: destroy         !< Destroy the field.
      procedure, pass(self) :: initialize      !< Initialize the field.
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
   class(grid_object), intent(in)            :: self                                 !< The field.
   integer(I4P),       intent(in)            :: coordinates(4)                       !< Block coordinates.
   real(R8P),          intent(out), optional :: dx, dy, dz                           !< Space steps.
   real(R8P),          intent(out), optional :: emin(3), emax(3)                     !< Min/max abscissa of block.
   real(R8P),          intent(out), optional :: x_node(0-self%gc1:self%ni+self%gc2)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%gc3:self%nj+self%gc4)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%gc5:self%nk+self%gc6)  !< Z coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%gc1:self%ni+self%gc2)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%gc3:self%nj+self%gc4)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%gc5:self%nk+self%gc6)  !< Z coordinates.
   real(R8P)                                 :: emin_(3), emax_(3)                   !< Min/max abscissa of block, local var.
   real(R8P)                                 :: x_node_(0-self%gc1:self%ni+self%gc2) !< X coordinates, local var.
   real(R8P)                                 :: y_node_(0-self%gc3:self%nj+self%gc4) !< Y coordinates, local var.
   real(R8P)                                 :: z_node_(0-self%gc5:self%nk+self%gc6) !< Z coordinates, local var.
   real(R8P)                                 :: x_cell_(1-self%gc1:self%ni+self%gc2) !< X coordinates, local var.
   real(R8P)                                 :: y_cell_(1-self%gc3:self%nj+self%gc4) !< Y coordinates, local var.
   real(R8P)                                 :: z_cell_(1-self%gc5:self%nk+self%gc6) !< Z coordinates, local var.
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
   do i=0-self%gc1, self%ni+self%gc2
      x_node_(i) = emin_(1) + i * dx_
   enddo
   do j=0-self%gc3, self%nj+self%gc4
      y_node_(j) = emin_(2) + j * dy_
   enddo
   do k=0-self%gc5, self%nk+self%gc6
      z_node_(k) = emin_(3) + k * dz_
   enddo
   do i=1-self%gc1, self%ni+self%gc2
      x_cell_(i) = x_node_(i-1) + dx_ * 0.5_R8P
   enddo
   do j=1-self%gc3, self%nj+self%gc4
      y_cell_(j) = y_node_(j-1) + dy_ * 0.5_R8P
   enddo
   do k=1-self%gc5, self%nk+self%gc6
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

   elemental subroutine destroy(self)
   !< Destroy field.
   class(grid_object), intent(inout) :: self  !< The field.
   type(grid_object)                 :: fresh !< Fresh field.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, ni, nj, nk, gc, emin, emax)
   !< Initialize field.
   class(grid_object), intent(inout)        :: self    !< The field.
   integer(I4P),       intent(in), optional :: ni      !< Number of cells in X direction.
   integer(I4P),       intent(in), optional :: nj      !< Number of cells in Y direction.
   integer(I4P),       intent(in), optional :: nk      !< Number of cells in Z direction.
   integer(I4P),       intent(in), optional :: gc(6)   !< Number of ghost cells in each direction.
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
   if (present(gc)) self%gc1 = gc(1)
   if (present(gc)) self%gc2 = gc(2)
   if (present(gc)) self%gc3 = gc(3)
   if (present(gc)) self%gc4 = gc(4)
   if (present(gc)) self%gc5 = gc(5)
   if (present(gc)) self%gc6 = gc(6)
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine grid_assign_grid(lhs, rhs)
   !< Operator `=`.
   class(grid_object), intent(inout) :: lhs !< Left hand side.
   type(grid_object),  intent(in)    :: rhs !< Right hand side.

   lhs%domain_emin  = rhs%domain_emin
   lhs%domain_emax  = rhs%domain_emax
   lhs%ni           = rhs%ni
   lhs%nj           = rhs%nj
   lhs%nk           = rhs%nk
   lhs%gc1          = rhs%gc1
   lhs%gc2          = rhs%gc2
   lhs%gc3          = rhs%gc3
   lhs%gc4          = rhs%gc4
   lhs%gc5          = rhs%gc5
   lhs%gc6          = rhs%gc6
   endsubroutine grid_assign_grid
endmodule adam_grid_object
