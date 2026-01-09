!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
module adam_prism_pic_object
!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
! ADAM modules
use :: adam_mpih_object, only : mpih_object
use adam_field_object, only : field_object
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: INI_SECTION_NAME
public :: prism_pic_object
public :: particle_weighting_interface
public :: CIC_WEIGHTING_MODEL
public :: NGP_WEIGHTING_MODEL
public :: TSC_WEIGHTING_MODEL
public :: CIC_weighting
public :: NGP_weighting
public :: TSC_weighting

character(len=3), parameter :: INI_SECTION_NAME      = 'PIC' !< INI file section name for PIC configuration.
character(len=3), parameter :: CIC_WEIGHTING_MODEL   = 'CIC' !< CIC weighting model.
character(len=3), parameter :: NGP_WEIGHTING_MODEL   = 'NGP' !< NGP weighting model.
character(len=3), parameter :: TSC_WEIGHTING_MODEL   = 'TSC' !< TSC weighting model.
! PIC variables layout in q_pic array:
!q_pic(1) = x
!q_pic(2) = y
!q_pic(3) = z
!q_pic(4) = vx
!q_pic(5) = vy
!q_pic(6) = vz
!q_pic(7) = charge

type :: prism_pic_object
   type(mpih_object)         :: mpih                      !< MPI handler.
   integer(I4P)              :: particle_number = 0_I4P   !< Total number of particles.
   character(len=99)         :: particle_weighting_model  !< Particle weighting model.
   integer(I4P), allocatable :: neighbour_list(:,:)       !< Particle grid positions array.
contains
   procedure, pass(self) :: description                   !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                    !< Initialize IC.
   procedure, pass(self) :: load_from_file                !< Load config from file.
   procedure, pass(self) :: particle_cartesian_grid_index !< Compute the grid index corresponding to a particle position.
   procedure, pass(self) :: CIC_weighting                 !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_weighting                 !< Nearest Grid Point weighting of particle quantities to the grid.
   !procedure, pass(self) :: TSC_weighting    !< Triangular Shaped Cloud weighting of particle quantities to the grid.
endtype prism_pic_object

interface
   subroutine particle_weighting_interface(self, field, q, q_pic, nv)
   import :: prism_pic_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                                                                          !< External fields.
   type(field_object),                  intent(inout) :: field                                                             !< The field.
   real(R8P),                           intent(inout) :: q(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)   !< Field variables.
   real(R8P),                           intent(in)    :: q_pic(1:,1:)   !< PIC variables.
   integer(I4P),                        intent(in)    :: nv                                                                !< Number of variables.
   endsubroutine particle_weighting_interface
endinterface

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_pic_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable                   :: desc             !< Description.
   character(len=1), parameter                     :: NL=new_line('a') !< New line character.
   desc =       self%mpih%myrankstr//'PIC object description:'
   desc = desc//NL//self%mpih%myrankstr//'    Number of particles: '//trim(str(self%particle_number))
   desc = desc//NL//self%mpih%myrankstr//'    Particle weighting model: '//trim(self%particle_weighting_model)
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize PIC.
   class(prism_pic_object), intent(inout) :: self            !< External fields.
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_pic_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()

   allocate(self%neighbour_list(self%particle_number, 4))

   print '(A)', self%mpih%myrankstr//'prism_pic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load PIC configuration from file.
	class(prism_pic_object), intent(inout)   			 :: self             !< PIC object.
	type(file_ini),          intent(in)		  			 :: file_parameters  !< File handler.
   logical,                 intent(in), optional    :: go_on_fail      	!< Go on if load fails.
   logical                                          :: go_on_fail_     	!< Go on if load fails.
   integer(I4P)                                     :: error           	!< Error status.
   character(99)                                    :: buff       		!< Option character buffer.

	go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%particle_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%particle_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%particle_weighting_model = TSC_WEIGHTING_MODEL
	case default
		call self%mpih%error_stop(msg=': invalid particle weighting model ['//trim(adjustl(buff))//'] in  & 
      ['//INI_SECTION_NAME//'].(weighting_model)')
	endselect

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='particle_number', &
   val=self%particle_number, error=error)
   if (.not.go_on_fail_.and.error>0) & 
   call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(particle_number)')
   endsubroutine load_from_file

   subroutine particle_cartesian_grid_index(self, field, q_PIC)
   !< Compute the grid index corresponding to a particle position. Good for cartesian grids only.
   class(prism_pic_object), intent(inout) :: self                                                            !< External fields.
   type(field_object),      intent(in)    :: field                                                           !< The field.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)                                                    !< PIC variables.
   real(R8P)                              :: n                                                               !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p                                              !< Particle grid indices

   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, &
            nk=>field%grid%nk, ngc=>field%grid%ngc, dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:), &
            np => self%particle_number, e_min => field%grid%domain_emin, e_max => field%grid%domain_emax, &
            neighbour_list => self%neighbour_list)

   !Va completato considerando la presenza di più blocchi, questo funziona per un blocco solo
   do n = 1, np
      i_p = (q_PIC(n,1) - e_min(1)) / dx(1)
      j_p = (q_PIC(n,2) - e_min(2)) / dy(1)
      k_p = (q_PIC(n,3) - e_min(3)) / dz(1)
      b_p = 1 ! Single block only for now

      neighbour_list(n,1) = ceiling(b_p)
      neighbour_list(n,2) = ceiling(i_p)
      neighbour_list(n,3) = ceiling(j_p)
      neighbour_list(n,4) = ceiling(k_p)
   enddo
   endassociate
   endsubroutine particle_cartesian_grid_index

   subroutine NGP_weighting(self, field, q, q_PIC, nv)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                                                            !< External fields.
   type(field_object),      intent(inout) :: field                                                           !< The field.
   real(R8P),               intent(inout) :: q(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)                                                    !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                                              !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                                   !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p                                              !< Particle grid indices
   real(R8P)                              :: x_p, y_p, z_p                                                   !< Particle position scalar
   real(R8P)                              :: dx, dy, dz                                                      !< Grid spacing
   real(R8P)                              :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                             y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                             z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)
   real(R8P)                              :: wx, wy, wz                                                      !< Weighting factors
   real(R8P)                              :: cell_coord(3)                                                   !< Cell coordinates 

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(n,1)
      i_p = self%neighbour_list(n,2)
      j_p = self%neighbour_list(n,3)
      k_p = self%neighbour_list(n,4)

      ! Qua va capito come gestire la questione dei blocchi multipli
      call field%grid%cell_xyz(coordinates = field%coordinates(:,b_p), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera
      q(nv, i_p, j_p, k_p, b_p) = q(nv, i_p, j_p, k_p, b_p) + q_PIC(7,n)
      !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
      !se scritta in questo modo
   enddo
   endsubroutine NGP_weighting

   subroutine CIC_weighting(self, field, q, q_PIC, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                                                            !< External fields.
   type(field_object),      intent(inout) :: field                                                           !< The field.
   real(R8P),               intent(inout) :: q(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)                                                    !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                                              !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                                   !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p                                              !< Particle grid indices
   real(R8P)                              :: x_p, y_p, z_p                                                   !< Particle position scalar
   real(R8P)                              :: dx, dy, dz                                                      !< Grid spacing
   real(R8P)                              :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                             y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                             z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)
   real(R8P)                              :: wx, wy, wz                                                      !< Weighting factors
   real(R8P)                              :: cell_coord(3)                                                   !< Cell coordinates 

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(n,1)
      i_p = self%neighbour_list(n,2)
      j_p = self%neighbour_list(n,3)
      k_p = self%neighbour_list(n,4)

      ! Qua va capito come gestire la questione dei blocchi multipli
      call field%grid%cell_xyz(coordinates = field%coordinates(:,b_p), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
               if (abs((q_PIC(1,n) - cell_coord(1))/dx) <= 1.0_R8P) then
                  Wx = 1.0_R8P - abs((q_PIC(1,n) - cell_coord(1))/dx)
               else
                  Wx = 0.0_R8P
               end if   
               if (abs((q_PIC(2,n) - cell_coord(2))/dy) <= 1.0_R8P) then
                  Wy = 1.0_R8P - abs((q_PIC(2,n) - cell_coord(2))/dy)
               else
                  Wy = 0.0_R8P
               end if
               if (abs((q_PIC(3,n) - cell_coord(3))/dz) <= 1.0_R8P) then
                  Wz = 1.0_R8P - abs((q_PIC(3,n) - cell_coord(3))/dz)
               else
                  Wz = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_PIC(7,n) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endsubroutine CIC_weighting

   subroutine TSC_weighting(self, field, q, q_PIC, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                                                            !< External fields.
   type(field_object),      intent(inout) :: field                                                           !< The field.
   real(R8P),               intent(inout) :: q(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)                                                    !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                                              !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                                   !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p                                              !< Particle grid indices
   real(R8P)                              :: x_p, y_p, z_p                                                   !< Particle position scalar
   real(R8P)                              :: dx, dy, dz                                                      !< Grid spacing
   real(R8P)                              :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                             y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                             z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc)
   real(R8P)                              :: wx, wy, wz                                                      !< Weighting factors
   real(R8P)                              :: cell_coord(3)                                                   !< Cell coordinates 

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(n,1)
      i_p = self%neighbour_list(n,2)
      j_p = self%neighbour_list(n,3)
      k_p = self%neighbour_list(n,4)

      ! Qua va capito come gestire la questione dei blocchi multipli
      call field%grid%cell_xyz(coordinates = field%coordinates(:,b_p), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
               if (abs((q_PIC(1,n) - cell_coord(1))/dx) <= 0.5_R8P) then
                  Wx = 0.75_R8P - ((q_PIC(1,n) - cell_coord(1))/dx)**2
               elseif (abs((q_PIC(1,n) - cell_coord(1))/dx) <= 1.5_R8P .and. abs((q_PIC(1,n) - cell_coord(1))/dx) > 0.5_R8P) then
                  Wx = 0.5_R8P * (1.5_R8P - abs((q_PIC(1,n) - cell_coord(1))/dx))**2
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_PIC(2,n) - cell_coord(2))/dy) <= 0.5_R8P) then
                  Wy = 0.75_R8P - ((q_PIC(2,n) - cell_coord(2))/dy)**2
               elseif (abs((q_PIC(2,n) - cell_coord(2))/dy) <= 1.5_R8P .and. abs((q_PIC(2,n) - cell_coord(2))/dy) > 0.5_R8P) then
                  Wy = 0.5_R8P * (1.5_R8P - abs((q_PIC(2,n) - cell_coord(2))/dy))**2
               else
                  Wy = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_PIC(7,n) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endsubroutine TSC_weighting

endmodule adam_prism_pic_object