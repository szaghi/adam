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
   integer(I4P) :: ngc=2_I4P                             !< Number of ghost cells.
   integer(I4P) :: weight_neighbor(26)=0_I4P             !< Weight of neighbors (cells number).
   integer(I4P) :: bc_type(6)=0_I4P                      !< Type of boundary conditions in the 6 faces of grid.
   logical      :: is_ijk_periodic(3)=.false.            !< Flag to indicate if the direction i, j or k is periodic.
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
   real(R8P),          intent(out), optional :: x_node(0-self%ngc:self%ni+self%ngc)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_node(0-self%ngc:self%nj+self%ngc)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_node(0-self%ngc:self%nk+self%ngc)  !< Z coordinates.
   real(R8P),          intent(out), optional :: x_cell(1-self%ngc:self%ni+self%ngc)  !< X coordinates.
   real(R8P),          intent(out), optional :: y_cell(1-self%ngc:self%nj+self%ngc)  !< Y coordinates.
   real(R8P),          intent(out), optional :: z_cell(1-self%ngc:self%nk+self%ngc)  !< Z coordinates.
   real(R8P)                                 :: emin_(3), emax_(3)                   !< Min/max abscissa of block, local var.
   real(R8P)                                 :: x_node_(0-self%ngc:self%ni+self%ngc) !< X coordinates, local var.
   real(R8P)                                 :: y_node_(0-self%ngc:self%nj+self%ngc) !< Y coordinates, local var.
   real(R8P)                                 :: z_node_(0-self%ngc:self%nk+self%ngc) !< Z coordinates, local var.
   real(R8P)                                 :: x_cell_(1-self%ngc:self%ni+self%ngc) !< X coordinates, local var.
   real(R8P)                                 :: y_cell_(1-self%ngc:self%nj+self%ngc) !< Y coordinates, local var.
   real(R8P)                                 :: z_cell_(1-self%ngc:self%nk+self%ngc) !< Z coordinates, local var.
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
   do i=0-self%ngc, self%ni+self%ngc
      x_node_(i) = emin_(1) + i * dx_
   enddo
   do j=0-self%ngc, self%nj+self%ngc
      y_node_(j) = emin_(2) + j * dy_
   enddo
   do k=0-self%ngc, self%nk+self%ngc
      z_node_(k) = emin_(3) + k * dz_
   enddo
   do i=1-self%ngc, self%ni+self%ngc
      x_cell_(i) = x_node_(i-1) + dx_ * 0.5_R8P
   enddo
   do j=1-self%ngc, self%nj+self%ngc
      y_cell_(j) = y_node_(j-1) + dy_ * 0.5_R8P
   enddo
   do k=1-self%ngc, self%nk+self%ngc
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
   integer(I4P)                      :: ngc     !< Ghost cells number.
   integer(I4P)                      :: fec     !< Counter.
   integer(I4P)                      :: i       !< Counter.

   self%weight_neighbor = 1_I4P
   nijk = [self%ni, self%nj, self%nk]
   ngc = self%ngc
   do fec=1, 26
      do i=1, 3
         self%weight_neighbor(fec) = self%weight_neighbor(fec) * (  abs(fec_to_delta(i,fec))  * ngc + &
                                                                 (1-abs(fec_to_delta(i,fec))) * nijk(i))
      enddo
   enddo
   endsubroutine compute_weight_neighbor

   elemental subroutine destroy(self)
   !< Destroy field.
   class(grid_object), intent(inout) :: self  !< The grid.
   type(grid_object)                 :: fresh !< Fresh grid.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, ni, nj, nk, ngc, emin, emax, bc_type)
   !< Initialize field.
   class(grid_object), intent(inout)        :: self               !< The grid.
   integer(I4P),       intent(in), optional :: ni                 !< Number of cells in X direction.
   integer(I4P),       intent(in), optional :: nj                 !< Number of cells in Y direction.
   integer(I4P),       intent(in), optional :: nk                 !< Number of cells in Z direction.
   integer(I4P),       intent(in), optional :: ngc                !< Number of ghost cells.
   real(R8P),          intent(in), optional :: emin(3)            !< Coordinates of minium abscissa.
   real(R8P),          intent(in), optional :: emax(3)            !< Coordinates of maxium abscissa.
   integer(I4P),       intent(in), optional :: bc_type(6)         !< Type of boundary conditions in the 6 faces of grid.

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
   if (present(ni))  self%ni = ni
   if (present(nj))  self%nj = nj
   if (present(nk))  self%nk = nk
   if (present(ngc)) self%ngc = ngc
   call self%compute_weight_neighbor
   if (present(bc_type)) self%bc_type = bc_type
   if (any(self%bc_type(1:2)==BC_PERIODIC)) self%is_ijk_periodic(1) = .true.
   if (any(self%bc_type(3:4)==BC_PERIODIC)) self%is_ijk_periodic(2) = .true.
   if (any(self%bc_type(5:6)==BC_PERIODIC)) self%is_ijk_periodic(3) = .true.
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
   lhs%ngc             = rhs%ngc
   lhs%weight_neighbor = rhs%weight_neighbor
   lhs%bc_type         = rhs%bc_type
   lhs%is_ijk_periodic = rhs%is_ijk_periodic
   endsubroutine grid_assign_grid
endmodule adam_grid_object
