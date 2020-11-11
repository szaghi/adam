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
   integer(I4P) :: ni=4_I4P       !< Number of cells in i direction.
   integer(I4P) :: nj=4_I4P       !< Number of cells in j direction.
   integer(I4P) :: nk=4_I4P       !< Number of cells in k direction.
   integer(I4P) :: gc1=0_I4P      !< Number of ghost cells in i- direction for boundary conditions.
   integer(I4P) :: gc2=0_I4P      !< Number of ghost cells in i+ direction for boundary conditions.
   integer(I4P) :: gc3=0_I4P      !< Number of ghost cells in j- direction for boundary conditions.
   integer(I4P) :: gc4=0_I4P      !< Number of ghost cells in j+ direction for boundary conditions.
   integer(I4P) :: gc5=0_I4P      !< Number of ghost cells in k- direction for boundary conditions.
   integer(I4P) :: gc6=0_I4P      !< Number of ghost cells in k+ direction for boundary conditions.
   contains
      ! public methods
      procedure, pass(self) :: compute_emin_emax !< Compute emin/emax of a block.
      procedure, pass(self) :: dxyz              !< Return space deltas.
      procedure, pass(self) :: xyz               !< Return grids coordinates.
      procedure, pass(self) :: destroy           !< Destroy the field.
      procedure, pass(self) :: initialize        !< Initialize the field.
      ! operators
      generic :: assignment(=) => grid_assign_grid      !< Overload `=`.
      procedure, pass(lhs), private :: grid_assign_grid !< Operator `=`.
endtype grid_object

contains
   ! public methods
   subroutine compute_emin_emax(self, coordinates, emin, emax)
   !< Compute emin/emax of a block.
   class(grid_object), intent(in)  :: self           !< The field.
   integer(I4P),       intent(in)  :: coordinates(4) !< Block coordinates.
   real(R8P),          intent(out) :: emin(3)        !< Coordinates of minimum abscissa of block.
   real(R8P),          intent(out) :: emax(3)        !< Coordinates of maximum abscissa of block.
   real(R8P)                       :: dx, dy, dz     !< Domain delta space.
   real(R8P)                       :: dxl, dyl, dzl  !< Local delta space.
   integer(I4P)                    :: i, j, k, l     !< Counter.

   dx = self%domain_emax(1) - self%domain_emin(1)
   dy = self%domain_emax(2) - self%domain_emin(2)
   dz = self%domain_emax(3) - self%domain_emin(3)
   i = coordinates(1)
   j = coordinates(2)
   k = coordinates(3)
   l = coordinates(4)
   dxl = dx / 2**l
   dyl = dy / 2**l
   dzl = dz / 2**l
   emin(1) = i * dxl ; emax(1) = emin(1) + dxl
   emin(2) = j * dyl ; emax(2) = emin(2) + dyl
   emin(3) = k * dzl ; emax(3) = emin(3) + dzl
   endsubroutine compute_emin_emax

   pure function dxyz(self, emin, emax, axis)
   !< Return space deltas.
   class(grid_object), intent(in) :: self    !< The field.
   real(R8P),          intent(in) :: emin(3) !< Coordinates of minimum abscissa of block.
   real(R8P),          intent(in) :: emax(3) !< Coordinates of maximum abscissa of block.
   character(1),       intent(in) :: axis    !< Axis direction queried ['x','y','z'].
   real(R8P), allocatable         :: dxyz    !< Grid delta space.

   select case(axis)
   case('x')
      dxyz = (emax(1) - emin(1)) / self%ni
   case('y')
      dxyz = (emax(2) - emin(2)) / self%nj
   case('z')
      dxyz = (emax(3) - emin(3)) / self%nk
   endselect
   endfunction dxyz

   pure function xyz(self, emin, emax, axis)
      !< Return grids coordinates
   class(grid_object), intent(in) :: self    !< The field.
   real(R8P),          intent(in) :: emin(3) !< Coordinates of minimum abscissa of block.
   real(R8P),          intent(in) :: emax(3) !< Coordinates of maximum abscissa of block.
   character(1),       intent(in) :: axis    !< Axis direction queried ['x','y','z'].
   real(R8P), allocatable         :: xyz(:)  !< Grid coordinates.
   real(R8P)                      :: dxyz    !< Space delta.
   integer(I4P)                   :: i       !< Counter.

   dxyz = self%dxyz(emin=emin, emax=emax, axis=axis)
   select case(axis)
   case('x')
      allocate(xyz(0:self%ni))
      do i=0, self%ni
         xyz(i) = emin(1) + i * dxyz
      enddo
   case('y')
      allocate(xyz(0:self%nj))
      do i=0, self%nj
         xyz(i) = emin(2) + i * dxyz
      enddo
   case('z')
      allocate(xyz(0:self%nk))
      do i=0, self%nk
         xyz(i) = emin(3) + i * dxyz
      enddo
   endselect
   endfunction xyz

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
