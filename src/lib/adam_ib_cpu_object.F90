!< ADAM, IB class definition, CPU backend.
module adam_ib_cpu_object
!< ADAM, IB class definition, CPU backend.

use adam_field_object, only : field_object
use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: ib_cpu_object
public :: BCS_VISCOUS
public :: BCS_EULER

integer(I4P), parameter :: BCS_VISCOUS = 1_I4P !< Visous wall.
integer(I4P), parameter :: BCS_EULER   = 2_I4P !< Inviscid wall.

type :: ib_cpu_object
   !< IB class definition, CPU backend.
   type(mpih_object)         :: mpih            !< MPI handler.
   integer(I4P)              :: solids_number=0 !< Number of solids (only 1 supported now).
   integer(I4P), allocatable :: bcs_type(:)     !< Immersed boundary condition type.
   real(R8P),    allocatable :: bcs_vars(:, :)  !< Variables' array for immersed boundary conditions.
   ! Pointers to ADAM data for easy handling.
   type(field_object), pointer :: field=>null() !< The field.
   type(grid_object),  pointer :: grid =>null() !< The grid.
   ! Large arrays.
   real(R8P), allocatable ::  phi(:,:,:,:,:) !< IB distance function.
   contains
      ! public methods
      procedure, pass(self) :: description      !< Return pretty-printed object description.
      procedure, pass(self) :: evolve_eikonal_q !< Evolve eikonal q.
      procedure, pass(self) :: initialize       !< Initialize IB.
      procedure, pass(self) :: move_phi         !< Move phi and the actual ptree representation.
      procedure, pass(self) :: update_phi       !< Update IB distance function.
endtype ib_cpu_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(ib_cpu_object), intent(in) :: self             !< IB.
   character(len=:), allocatable    :: desc             !< Description.
   character(len=1), parameter      :: NL=new_line('a') !< New line character.
   integer(I4P)                     :: s                !< Counter.

   desc =       self%mpih%myrankstr//'IB main data'//NL
   desc = desc//self%mpih%myrankstr//'  solids number: '//trim(str(self%solids_number))
   do s=1, self%solids_number
      desc = desc//NL//self%mpih%myrankstr//'  BC type('//trim(str(s,.true.))//'):    '//trim(str(self%bcs_type(:)  ))
   enddo
   endfunction description

   subroutine initialize(self, grid, field, file_parameters)
   !< Initialize the equation.
   class(ib_cpu_object), intent(inout)      :: self            !< IB.
   type(grid_object),    intent(in), target :: grid            !< The grid.
   type(field_object),   intent(in), target :: field           !< The field.
   type(file_ini),       intent(inout)      :: file_parameters !< INI file handler.

   call self%mpih%initialize

   print '(A)', self%mpih%myrankstr//'ib_cpu_object%initialize start'

   ! associate ADAM main data
   self%field         => field
   self%grid          => grid

   call load_solids_from_ini_file

   ! allocate large arrays
   associate(ngc=>self%grid%ngc, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, nb=>self%field%nb, &
             solids_number=>self%solids_number)
   if (solids_number > 0) then
      allocate(self%phi(1:solids_number, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%phi = -1._R8P
   endif
   endassociate

   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'ib_cpu_object%initialize finish'
   contains
      subroutine load_solids_from_ini_file
      !< Parse immersed boundary solids setting from input file.
      character(999) :: buf_CHAR !< String buffer.
      integer(I4P)   :: buf_I4   !< I4 buffer.
      character(999) :: sname    !< Section name.
      integer(I4P)   :: i_solid  !< Counter.

      call file_parameters%get(section_name='solids', option_name='solids_number', val=buf_I4)
      self%solids_number = buf_I4
      if (self%solids_number > 0) then
         allocate(self%bcs_type(self%solids_number))
         allocate(self%bcs_vars(6,self%solids_number))
         do i_solid=1, self%solids_number
            sname = 'solid_'//trim(str(i_solid,.true.))
            call file_parameters%get(section_name=sname, option_name='name', val=buf_CHAR)
            call file_parameters%get(section_name=sname, option_name='bcs_type', val=buf_I4)
            self%bcs_type(i_solid) = buf_I4
         enddo
      endif
      endsubroutine load_solids_from_ini_file
   endsubroutine initialize

   subroutine move_phi(self, velocity, s)
   !< Move phi.
   class(ib_cpu_object), intent(inout) :: self                             !< IB.
   real(R8P),            intent(in)    :: velocity(3)                      !< Velocity of the movement.
   integer(I4P),         intent(in)    :: s                                !< Solid index.
   ! real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi !< Eikonal direction.
   ! integer(I4P)                        :: b, i, j, k                       !< Counter.

   ! associate (ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, blocks_number=>self%field%blocks_number, phi=>self%phi)
   ! n_phi_x = velocity(1)
   ! n_phi_y = velocity(2)
   ! n_phi_z = velocity(3)
   ! n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
   ! n_phi = 0.9_R8P / n_phi
   ! n_phi_x = n_phi_x * n_phi
   ! n_phi_y = n_phi_y * n_phi
   ! n_phi_z = n_phi_z * n_phi

   ! do b=1, blocks_number
   ! do k=1, nk
   ! do j=1, nj
   ! do i=1, ni
   !    dphi(i,j,k,b) = 0._R8P
   !    if (n_phi_x > 0._R8P) then
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_x) * (phi(s,i,j,k,b) - phi(s,i-1,j,k,b))
   !    else
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_x) * (phi(s,i,j,k,b) - phi(s,i+1,j,k,b))
   !    endif
   !    if (n_phi_y > 0._R8P) then
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_y) * (phi(s,i,j,k,b) - phi(s,i,j-1,k,b))
   !    else
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_y) * (phi(s,i,j,k,b) - phi(s,i,j+1,k,b))
   !    endif
   !    if (n_phi_z > 0._R8P) then
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_z) * (phi(s,i,j,k,b) - phi(s,i,j,k-1,b))
   !    else
   !       dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_z) * (phi(s,i,j,k,b) - phi(s,i,j,k+1,b))
   !    endif
   ! enddo
   ! enddo
   ! enddo
   ! enddo

   ! do b=1, blocks_number
   ! do k=1, nk
   ! do j=1, nj
   ! do i=1, ni
   !    phi(s,i,j,k,b) = phi(s,i,j,k,b) - dphi(i,j,k,b)
   ! enddo
   ! enddo
   ! enddo
   ! enddo
   ! endassociate
   endsubroutine move_phi

   subroutine update_phi(self)
   !< Update phi.
   class(ib_cpu_object), intent(inout) :: self           !< IB.
   integer(I4P)                        :: b, i, j, k, ib !< Counter.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,                                   &
             phi=>self%phi, solids_number=>self%solids_number)
   if (solids_number > 0) then
      print '(A)', self%mpih%myrankstr//'ib_cpu_object%update IB distance start'
      do b=1,blocks_number
         do i=1-ngc,ni+ngc
            do j=1-ngc,nj+ngc
               do k=1-ngc,nk+ngc
                  do ib=1, solids_number
                     phi(ib,i,j,k,b) = - (sqrt((x_cell(i,b)-10._R8P)**2+(y_cell(j,b)-10._R8P)**2+(z_cell(k,b)-10._R8P)**2)-1.0_R8P)
                  enddo
               enddo
            enddo
         enddo
      enddo
      print '(A)', self%mpih%myrankstr//'ib_cpu_object%update IB distance finish'
   endif
   endassociate
   endsubroutine update_phi

   subroutine evolve_eikonal_q(self, q)
   !< Evolve eikonal q.
   class(ib_cpu_object), intent(in)    :: self                             !< IB.
   real(R8P),            intent(inout) ::  q(1:,               &
                                             1-self%grid%ngc:, &
                                             1-self%grid%ngc:, &
                                             1-self%grid%ngc:, &
                                             1:)                           !< Conservative variables.
   real(R8P)                           :: dq(1:self%field%nv)              !< Conservative variables differences.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi !< Eikonal directions.
   integer(I4P)                        :: i, j, k, b, s                    !< Counter.

   associate(blocks_number=>self%field%blocks_number, ni=>self%grid%ni, nj=>self%grid%nj, nk=>self%grid%nk, ngc=>self%grid%ngc, &
             nv=>self%field%nv, phi=>self%phi, solids_number=>self%solids_number)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1,ni
               solids_loop : do s=1, size(phi, dim=1)
                  if (phi(s,i,j,k,b) > 0._R8P) then
                     ! compute dq
                     n_phi_x = (phi(s,i+1,j,k,b) - phi(s,i-1,j,k,b))
                     n_phi_y = (phi(s,i,j+1,k,b) - phi(s,i,j-1,k,b))
                     n_phi_z = (phi(s,i,j,k+1,b) - phi(s,i,j,k-1,b))
                     n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
                     n_phi = 0.9_R8P / n_phi
                     n_phi_x = n_phi_x * n_phi
                     n_phi_y = n_phi_y * n_phi
                     n_phi_z = n_phi_z * n_phi
                     dq(:) = 0._R8P
                     if (n_phi_x > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_x) * (q(:,i,j,k,b) - q(:,i-1,j,k,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_x) * (q(:,i,j,k,b) - q(:,i+1,j,k,b))
                     endif
                     if (n_phi_y > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_y) * (q(:,i,j,k,b) - q(:,i,j-1,k,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_y) * (q(:,i,j,k,b) - q(:,i,j+1,k,b))
                     endif
                     if (n_phi_z > 0._R8P) then
                        dq(:) = dq(:) + abs(n_phi_z) * (q(:,i,j,k,b) - q(:,i,j,k-1,b))
                     else
                        dq(:) = dq(:) + abs(n_phi_z) * (q(:,i,j,k,b) - q(:,i,j,k+1,b))
                     endif
                     ! evolve q
                     q(:,i,j,k,b) = q(:,i,j,k,b) - dq(:)
                     exit solids_loop
                  endif
               enddo solids_loop
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine evolve_eikonal_q
endmodule adam_ib_cpu_object
