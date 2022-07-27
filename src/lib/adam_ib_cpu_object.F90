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

type :: ib_cpu_object
   !< IB class definition, CPU backend.
   type(mpih_object)         :: mpih            !< MPI handler.
   integer(I4P)              :: solids_number=0 !< Number of solids (only 1 supported now).
   integer(I4P), allocatable :: bcs_type(:)     !< Immersed boundary condition type.
   ! Pointers to ADAM data for easy handling.
   type(field_object), pointer :: field        =>null() !< The field.
   type(grid_object),  pointer :: grid         =>null() !< The grid.
   integer(I4P),       pointer :: nb           =>null() !< Total blocks number for MPI.
   integer(I4P),       pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P),       pointer :: nv           =>null() !< Number of field variables.
   integer(I4P),       pointer :: ngc          =>null() !< Number of ghost cells.
   integer(I4P),       pointer :: ni           =>null() !< Number of cells in i direction.
   integer(I4P),       pointer :: nj           =>null() !< Number of cells in j direction.
   integer(I4P),       pointer :: nk           =>null() !< Number of cells in k direction.
   ! Large arrays.
   real(R8P), allocatable ::  phi(:,:,:,:,:) !< IB distance function.
   real(R8P), allocatable :: dphi(  :,:,:,:) !< IB distance function delta.
   contains
      ! public methods
      procedure, pass(self) :: evolve_eikonal_q !< Evolve eikonal q.
      procedure, pass(self) :: initialize       !< Initialize IB.
      procedure, pass(self) :: move_phi         !< Move phi and the actual ptree representation.
      procedure, pass(self) :: update_phi       !< Update IB distance function.
endtype ib_cpu_object

contains
   ! public methods
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
   self%nb            => field%nb
   self%blocks_number => field%blocks_number
   self%nv            => field%nv
   self%ngc           => grid%ngc
   self%ni            => grid%ni
   self%nj            => grid%nj
   self%nk            => grid%nk

   call load_solids_from_ini_file

   ! allocate large arrays
   associate(ngc=>self%ngc, ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb, solids_number=>self%solids_number)
   if (solids_number > 0) then
      allocate(self%phi(1:solids_number, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      allocate(self%dphi(                1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, 1:nb))
      self%phi = -1._R8P
   endif
   endassociate

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
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi !< Eikonal direction.
   integer(I4P)                        :: b, i, j, k                       !< Counter.

   associate (ni=>self%ni, nj=>self%nj, nk=>self%nk, blocks_number=>self%blocks_number, dphi=>self%dphi, phi=>self%phi)
   n_phi_x = velocity(1)
   n_phi_y = velocity(2)
   n_phi_z = velocity(3)
   n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
   n_phi = 0.9_R8P / n_phi
   n_phi_x = n_phi_x * n_phi
   n_phi_y = n_phi_y * n_phi
   n_phi_z = n_phi_z * n_phi

   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      dphi(i,j,k,b) = 0._R8P
      if (n_phi_x > 0._R8P) then
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_x) * (phi(s,i,j,k,b) - phi(s,i-1,j,k,b))
      else
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_x) * (phi(s,i,j,k,b) - phi(s,i+1,j,k,b))
      endif
      if (n_phi_y > 0._R8P) then
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_y) * (phi(s,i,j,k,b) - phi(s,i,j-1,k,b))
      else
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_y) * (phi(s,i,j,k,b) - phi(s,i,j+1,k,b))
      endif
      if (n_phi_z > 0._R8P) then
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_z) * (phi(s,i,j,k,b) - phi(s,i,j,k-1,b))
      else
         dphi(i,j,k,b) = dphi(i,j,k,b) + abs(n_phi_z) * (phi(s,i,j,k,b) - phi(s,i,j,k+1,b))
      endif
   enddo
   enddo
   enddo
   enddo

   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      phi(s,i,j,k,b) = phi(s,i,j,k,b) - dphi(i,j,k,b)
   enddo
   enddo
   enddo
   enddo
   endassociate
   endsubroutine move_phi

   subroutine update_phi(self)
   !< Update phi.
   class(ib_cpu_object), intent(inout) :: self           !< IB.
   integer(I4P)                        :: b, i, j, k, ib !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, &
             x_cell=>self%field%x_cell, y_cell=>self%field%y_cell, z_cell=>self%field%z_cell,         &
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

   subroutine evolve_eikonal_q(self, dq, q)
   !< Evolve eikonal q.
   class(ib_cpu_object), intent(in)    :: self                                          !< IB.
   real(R8P),            intent(inout) :: dq(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Conservative variables differences.
   real(R8P),            intent(inout) ::  q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, v                                 !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, &
             phi=>self%phi, solids_number=>self%solids_number)
   if ((solids_number > 0).and.(blocks_number > 0)) then
      call compute_eikonal_dq(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, phi=phi, dq=dq, q=q)
   endif

   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1,ni
               do v=1, nv
                  if (phi(1,i,j,k,b) > 0._R8P) then
                     q(v,i,j,k,b) = q(v,i,j,k,b) - dq(v,i,j,k,b)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine evolve_eikonal_q

   ! non TBP methods
   subroutine compute_eikonal_dq(ni, nj, nk, ngc, nv, blocks_number, phi, q, dq)
   !< Compute eikonal q-differences.
   integer(I4P), intent(in)    :: ni                               !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                               !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                               !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                              !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                               !< Number of conservative varibales.
   integer(I4P), intent(in)    :: blocks_number                    !< Number of blocks.
   real(R8P),    intent(in)    :: phi(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Distance function.
   real(R8P),    intent(in)    ::   q(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(inout) ::  dq(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables differences.
   integer(I4P)                :: i, j, k, b, v                    !< Counter.
   real(R8P)                   :: n_phi_x, n_phi_y, n_phi_z, n_phi !< Eikonal directions.

   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               if (phi(1,i,j,k,b) > 0._R8P) then
                  n_phi_x = (phi(1,i+1,j,k,b) - phi(1,i-1,j,k,b))
                  n_phi_y = (phi(1,i,j+1,k,b) - phi(1,i,j-1,k,b))
                  n_phi_z = (phi(1,i,j,k+1,b) - phi(1,i,j,k-1,b))
                  n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
                  n_phi = 0.9_R8P / n_phi
                  n_phi_x = n_phi_x * n_phi
                  n_phi_y = n_phi_y * n_phi
                  n_phi_z = n_phi_z * n_phi
                  do v=1, nv
                     dq(v,i,j,k,b) = 0._R8P
                  enddo
                  if (n_phi_x > 0._R8P) then
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_x) * (q(v,i,j,k,b) - q(v,i-1,j,k,b))
                     enddo
                  else
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_x) * (q(v,i,j,k,b) - q(v,i+1,j,k,b))
                     enddo
                  endif
                  if (n_phi_y > 0._R8P) then
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_y) * (q(v,i,j,k,b) - q(v,i,j-1,k,b))
                     enddo
                  else
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_y) * (q(v,i,j,k,b) - q(v,i,j+1,k,b))
                     enddo
                  endif
                  if (n_phi_z > 0._R8P) then
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_z) * (q(v,i,j,k,b) - q(v,i,j,k-1,b))
                     enddo
                  else
                     do v=1, nv
                        dq(v,i,j,k,b) = dq(v,i,j,k,b) + abs(n_phi_z) * (q(v,i,j,k,b) - q(v,i,j,k+1,b))
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_eikonal_dq
endmodule adam_ib_cpu_object
