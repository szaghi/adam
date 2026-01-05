!< ADAM, PRISM Initial Conditions class definition, CPU backend.
module adam_prism_ic_object
!< ADAM, PRISM Initial Conditions class definition, CPU backend.

! ADAM modules
use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
! PRISM modules
use adam_prism_physics_object, only : prism_physics_object
! third party modules
use finer
use penf
use adam_prism_parameters

implicit none
private
public :: IC_TYPE_RP
public :: IC_TYPE_VACUUM
public :: IC_TYPE_PLANE_WAVE
public :: IC_TYPE_RMF
public :: IC_TYPE_MAGNETIC_NOZZLE
public :: IC_TYPE_RMF_NOZZLE
public :: prism_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions"     !< INI (config) file section name containing IC configs.
character(len=6),  parameter :: IC_TYPE_VACUUM="vacuum"                   !< Vacuum IC TYPE parameter.
character(len=15), parameter :: IC_TYPE_RP="riemann-problem"              !< Riemann Problem IC TYPE parameter.   
character(len=10), parameter :: IC_TYPE_PLANE_WAVE="plane_wave"           !< Riemann Problem IC TYPE parameter.
character(len=10), parameter :: IC_TYPE_RMF="rmf_field"             !< Rotating Magnetic Field IC TYPE parameter.
character(len=15), parameter :: IC_TYPE_MAGNETIC_NOZZLE="magnetic_nozzle" !< Nozzle IC TYPE parameter.
character(len=15), parameter :: IC_TYPE_RMF_NOZZLE="rmf_magnetic_nozzle"  !< Rotating Magnetic Field Nozzle IC TYPE parameter.

type :: prism_ic_object
   !< Initial Conditions class definition, CPU backend.
   type(mpih_object)         :: mpih                 !< MPI handler.
   integer(I4P)              :: amr_iterations=1_I4P !< Number of AMR iterations imposing IC.
   character(:), allocatable :: ic_type              !< IC type.
   integer(I4P)              :: regions_number=1_I4P !< Number of IC regions.
   real(R8P), allocatable    :: q(:,:)               !< Primitive variables (Dx,Dy,Dz,Bx,By,Bz,Jx,Jy,Jz).
   real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   real(R8P)                 :: kx=0.0_R8P           !< Plane wave number in x direction.
   real(R8P)                 :: ky=0.0_R8P           !< Plane wave number in y direction.
   real(R8P)                 :: kz=0.0_R8P           !< Plane wave number in z direction.
   real(R8P)                 :: lambda=0.0_R8P       !< Plane wave wavelength.
   real(R8P)                 :: B0=0.0_R8P           !< Plane wave background magnetic field amplitude.
   real(R8P)                 :: RMF_frequency        !< Rotating magnetic field frequency.
   real(R8P)                 :: RMF_B_amplitude      !< Rotating magnetic field amplitude.
	character(len=99)         :: RMF_rotation_axis 	  !< Rotating magnetic field rotation axis (X, Y, Z).
	integer(I4P)              :: alpha                !< RMF rotation axis coordinate 1
	integer(I4P)              :: beta                 !< RMF rotation axis coordinate 2
	integer(I4P)              :: gamma                !< RMF rotation axis coordinate 3
   contains
      ! public methods
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize IC.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_initial_conditions !< Set initial conditions on PRISM fields.
endtype prism_ic_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_ic_object), intent(in) :: self             !< IC.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: r, v             !< Counter.

   desc =       self%mpih%myrankstr//'IC main data'//NL
   desc = desc//self%mpih%myrankstr//'  regions number: '//trim(str(self%regions_number))
   do r=1, self%regions_number
      desc = desc//NL//self%mpih%myrankstr//'  region('//trim(str(r,.true.))//')'
      do v=1, size(self%q, dim=1)
         desc = desc//NL//self%mpih%myrankstr//'    q('//trim(str(v,.true.))//'): '//trim(str(self%q(v,r)))
      enddo
         desc = desc//NL//self%mpih%myrankstr//'    emin: '//trim(str(self%emin(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    emax: '//trim(str(self%emax(:,r)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(prism_ic_object), intent(inout) :: self            !< IC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_ic_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_ic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_ic_object), intent(inout)        :: self            !< IC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                    :: sname           !< Section name.
   integer(I4P)                                 :: i               !< Counter.
   integer(I4P)                                 :: error           !< Error status.
   character(99)                                :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='amr_iterations', val=self%amr_iterations, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(amr_iterations)')
   self%amr_iterations = max(0_I4P, self%amr_iterations)
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='type', val=buff_char, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(type)')
   self%ic_type = trim(adjustl(buff_char))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')
   if (self%regions_number>=1) then
      allocate(   self%q(1:9, 1:self%regions_number))
      allocate(self%emin(1:3, 1:self%regions_number))
      allocate(self%emax(1:3, 1:self%regions_number))
      do i=1, self%regions_number
         sname = INI_SECTION_NAME//'_region_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='Dx', val=self%q(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dx)')
         call file_parameters%get(section_name=sname, option_name='Dy', val=self%q(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dy)')
         call file_parameters%get(section_name=sname, option_name='Dz', val=self%q(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dz)')
         call file_parameters%get(section_name=sname, option_name='Bx', val=self%q(4,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Bx)')
         call file_parameters%get(section_name=sname, option_name='By', val=self%q(5,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(By)')
         call file_parameters%get(section_name=sname, option_name='Bz', val=self%q(6,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Bz)')

         call file_parameters%get(section_name=sname, option_name='Jx', val=self%q(7,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jx)')
         call file_parameters%get(section_name=sname, option_name='Jy', val=self%q(8,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jy)')
         call file_parameters%get(section_name=sname, option_name='Jz', val=self%q(9,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jz)')

         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_x)')
         call file_parameters%get(section_name=sname, option_name='emin_y', val=self%emin(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_y)')
         call file_parameters%get(section_name=sname, option_name='emin_z', val=self%emin(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_z)')
         call file_parameters%get(section_name=sname, option_name='emax_x', val=self%emax(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_x)')
         call file_parameters%get(section_name=sname, option_name='emax_y', val=self%emax(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_y)')
         call file_parameters%get(section_name=sname, option_name='emax_z', val=self%emax(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_z)')
      enddo
   endif
   if (self%ic_type == IC_TYPE_PLANE_WAVE) then
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='kx', val=self%kx, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(kx)')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ky', val=self%ky, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(ky)')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='kz', val=self%kz, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(kz)')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='lambda', val=self%lambda, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(lambda)')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='B0', val=self%B0, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(B0)')        
   endif 
   if (self%ic_type == IC_TYPE_RMF) then
      call file_parameters%get(section_name='external_fields', option_name='RMF_frequency', val=self%RMF_frequency, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')
      call file_parameters%get(section_name='external_fields', option_name='RMF_B_amplitude', val=self%RMF_B_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')
      call file_parameters%get(section_name='external_fields', option_name='RMF_rotation_axis', &
                           val=buff_char, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_rotation_axis)')
      self%RMF_rotation_axis = trim(buff_char)
      self%RMF_rotation_axis = trim(self%RMF_rotation_axis)
	   select case(self%RMF_rotation_axis)
	   case('X', 'x')
	   	self%alpha = 2_I4P
	   	self%beta  = 3_I4P
	   	self%gamma = 1_I4P
	   case('Y', 'y')
	   	self%alpha = 3_I4P
	   	self%beta  = 1_I4P
	   	self%gamma = 2_I4P
	   case('Z', 'z')
	   	self%alpha = 1_I4P
	   	self%beta  = 2_I4P
	   	self%gamma = 3_I4P
	   endselect
   endif

   if (self%ic_type == IC_TYPE_MAGNETIC_NOZZLE) then
      ! to be added load magnetic nozzle parameters
   endif
   if (self%ic_type == IC_TYPE_RMF_NOZZLE) then
      call file_parameters%get(section_name='external_fields', option_name='RMF_frequency', val=self%RMF_frequency, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_frequency)')
      call file_parameters%get(section_name='external_fields', option_name='RMF_B_amplitude', val=self%RMF_B_amplitude, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(RMF_B_amplitude)')
  endif
   endsubroutine load_from_file

   subroutine set_initial_conditions(self, physics, field, q)
      !< Set initial conditions on PRISM fields.
      class(prism_ic_object),     intent(in)    :: self                 !< IC.
      type(prism_physics_object), intent(in)    :: physics              !< Fluids physiscs.
      type(field_object),         intent(in)    :: field                !< Field object.
      real(R8P),                  intent(inout) :: q(1:,             &
                                                   1-field%grid%ngc:,&
                                                   1-field%grid%ngc:,&
                                                   1-field%grid%ngc:,&
                                                   1:)                  !< Field cell centered variables.
      real(R8P)                                 :: x_cell(1-field%grid%ngc:field%grid%ni+field%grid%ngc), &
                                                   y_cell(1-field%grid%ngc:field%grid%nj+field%grid%ngc), &
                                                   z_cell(1-field%grid%ngc:field%grid%nk+field%grid%ngc) 
                                                                        !< Vettori posizione centro celle del blocco b
      integer(I4P)                              :: b, i, j, k, ri, var  !< Counter.
	   real(R8P) 										   :: B_r, B_theta 			!< Radial and azimuthal components of the rotating magnetic field
	   real(R8P)										   :: theta, alfa, thetaabs!< Angles in cylindrical coordinates
	   real(R8P)										   :: cell_coord(3)	
   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
             nv=>physics%nv, nv_c=>physics%nv_c, nv_cl=>physics%nv_cl)
   select case(self%ic_type)
   case(IC_TYPE_VACUUM) ! vacuum initial conditions
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do var=1, nv
                     q(var,i,j,k,b) = 0.0_R8P
                  enddo
               enddo
            enddo
         enddo
      enddo
   case(IC_TYPE_PLANE_WAVE) !plane wave initial conditions
      do b=1, blocks_number
         call field%grid%cell_xyz(coordinates = field%coordinates(:,b), &
               x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  q(1,i,j,k,b) = self%B0*C0*EPS0*self%kz*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !Dx
                  q(2,i,j,k,b) = self%B0*C0*EPS0*self%kx*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !Dy
                  q(3,i,j,k,b) = self%B0*C0*EPS0*self%ky*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !Dz

                  q(4,i,j,k,b) = self%B0*self%ky*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !Bx
                  q(5,i,j,k,b) = self%B0*self%kz*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !By
                  q(6,i,j,k,b) = self%B0*self%kx*cos(self%kx*2*PI/self%lambda*x_cell(i)+ &
                                 self%ky*2*PI/self%lambda*y_cell(j)+self%kz*2*PI/self%lambda*z_cell(k)) !Bz
                  do var= (nv_c-nv_cl+1), nv
                     q(var,i,j,k,b) = 0.0_R8P
                  enddo
               enddo
            enddo
         enddo
      enddo
   case(IC_TYPE_RMF) !rotating magnetic field initial conditions
   associate(alpha=>self%alpha, beta=>self%beta, gamma=>self%gamma)	
   do b = 1, blocks_number
		call field%grid%cell_xyz(coordinates = field%coordinates(:,b), x_cell = x_cell, y_cell = y_cell, z_cell = z_cell)
         do i = 1, ni
            do j = 1, nj
               do k = 1, nk
					   cell_coord = [x_cell(i), y_cell(j), z_cell(k)]
                  theta = atan(cell_coord(beta)/cell_coord(alpha))
					   thetaabs = abs(atan(cell_coord(beta)/cell_coord(alpha)))
					   alfa = pi/2-thetaabs
					   if (cell_coord(alpha) > 0.0_R8P .and. cell_coord(beta) > 0.0_R8P) then ! 1 quadrante
                     theta = theta
                     B_r = self%RMF_B_amplitude*cos(-theta)
	                  B_theta = self%RMF_B_amplitude*sin(-theta)
					   	q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) + B_r*cos(thetaabs) - B_theta*cos(alfa)
					   	q(beta+3_I4P,i,j,k,b)  = q(beta+3_I4P,i,j,k,b)  + B_r*sin(thetaabs) + B_theta*sin(alfa)
					   	q(gamma,i,j,k,b) = sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*2*PI*self%RMF_frequency* &
					   		self%RMF_B_amplitude*cos(-theta)*EPS0
					   else if (cell_coord(alpha) < 0.0_R8P .and. cell_coord(beta) > 0.0_R8P) then ! 2 quadrante
                     theta = theta+PI
                     B_r = self%RMF_B_amplitude*cos(-theta)
	                  B_theta = self%RMF_B_amplitude*sin(-theta)
					   	q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) - B_r*cos(thetaabs) - B_theta*cos(alfa)
					   	q(beta+3_I4P,i,j,k,b)  = q(beta+3_I4P,i,j,k,b)  + B_r*sin(thetaabs) - B_theta*sin(alfa)
					   	q(gamma,i,j,k,b) = sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*2*PI*self%RMF_frequency* &
					   		self%RMF_B_amplitude*cos(-theta)*EPS0
					   else if (cell_coord(alpha) < 0.0_R8P .and. cell_coord(beta) < 0.0_R8P) then ! 3 quadrante
                     theta = theta+PI
                     B_r = self%RMF_B_amplitude*cos(-theta)
                     B_theta = self%RMF_B_amplitude*sin(-theta)
					   	q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) - B_r*cos(thetaabs) + B_theta*cos(alfa)
					   	q(beta+3_I4P,i,j,k,b)  = q(beta+3_I4P,i,j,k,b)  - B_r*sin(thetaabs) - B_theta*sin(alfa)
					   	q(gamma,i,j,k,b) = sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*2*PI*self%RMF_frequency* &
					   		self%RMF_B_amplitude*cos(-theta)*EPS0
					   else if (cell_coord(alpha) > 0.0_R8P .and. cell_coord(beta) < 0.0_R8P) then ! 4 quadrante
                     theta = theta+2*PI
                     B_r = self%RMF_B_amplitude*cos(-theta)
                     B_theta = self%RMF_B_amplitude*sin(-theta)
					   	q(alpha+3_I4P,i,j,k,b) = q(alpha+3_I4P,i,j,k,b) + B_r*cos(thetaabs) + B_theta*cos(alfa)
					   	q(beta+3_I4P,i,j,k,b)  = q(beta+3_I4P,i,j,k,b)  - B_r*sin(thetaabs) + B_theta*sin(alfa)
					   	q(gamma,i,j,k,b) = sqrt(cell_coord(alpha)**2 + cell_coord(beta)**2)*2*PI*self%RMF_frequency* &
					   		self%RMF_B_amplitude*cos(-theta)*EPS0
					   end if
               enddo
            enddo
         enddo
   enddo
   endassociate
   case default
      ! to be added error print
   endselect
   endassociate
   endsubroutine set_initial_conditions
endmodule adam_prism_ic_object
