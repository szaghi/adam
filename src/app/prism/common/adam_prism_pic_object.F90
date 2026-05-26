!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.
module adam_prism_pic_object
!< ADAM, PRISM Particle-in-Cell class definition, CPU backend.

! ADAM classes, libraries, parameters
use :: adam_field_object, only : field_object
use :: adam_grid_object,  only : grid_object
! ADAM singleton objects
use :: adam_mpih_global,  only : mpih
! PRISM modules
use :: adam_prism_parameters
! third party modules
use :: finer, only : file_ini
use :: penf,  only : I4P, R8P, str

implicit none
private
public :: INI_SECTION_NAME
public :: prism_pic_object
!public :: particle_weighting_interface
!public :: current_weighting_interface
!public :: field_weighting_interface
!public :: particle_weighting
!public :: current_weighting
!public :: field_weighting
public :: PLASMA_TYPE_PROBLEM
public :: SINGLE_PARTICLE_TYPE_PROBLEM
!public :: CIC_WEIGHTING_MODEL
!public :: NGP_WEIGHTING_MODEL
!public :: TSC_WEIGHTING_MODEL
!public :: ZEROD_FIELDS_WEIGHTING_MODEL
!public :: ONED_FIELDS_WEIGHTING_MODEL
public :: NUM_SCHEME_TIME_PIC_LEAPFROG
public :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
!public :: CIC_charge_weighting
!public :: NGP_charge_weighting
!public :: TSC_charge_weighting
!public :: CIC_current_weighting
!public :: NGP_current_weighting
!public :: TSC_current_weighting
!public :: zeroD_field_weighting
!public :: oneD_field_weighting

character(len=3 ), parameter :: INI_SECTION_NAME                 = 'PIC'             !< INI file section name for PIC configuration.
character(len=6 ), parameter :: PLASMA_TYPE_PROBLEM              = 'plasma'
                                                                                     !< Analyzing physical problem involving the
                                                                                     !< presence of plasma
character(len=15), parameter :: SINGLE_PARTICLE_TYPE_PROBLEM     = 'single_particle'
                                                                                     !< Analyzing physical problem involving the
                                                                                     !< presence of a single particle
character(len=3 ), parameter :: CIC_WEIGHTING_MODEL              = 'CIC'             !< CIC weighting model.
character(len=3 ), parameter :: NGP_WEIGHTING_MODEL              = 'NGP'             !< NGP weighting model.
character(len=3 ), parameter :: TSC_WEIGHTING_MODEL              = 'TSC'             !< TSC weighting model.
character(len=2 ), parameter :: ZEROD_FIELDS_WEIGHTING_MODEL     = '0D'              !< 0D field weighting.
character(len=2 ), parameter :: ONED_FIELDS_WEIGHTING_MODEL      = '1D'              !< 1D field weighting.
character(len=8 ), parameter :: NUM_SCHEME_TIME_PIC_LEAPFROG     = 'LEAPFROG'        !< Leapfrog numerical scheme for time operator.
character(len=11), parameter :: NUM_SCHEME_TIME_PIC_RUNGE_KUTTA  = 'RUNGE_KUTTA'
                                                                                     !< Runge-Kutta numerical scheme for time
                                                                                     !< operator.

! PIC variables layout in q_pic array:
!q_pic(1) = x
!q_pic(2) = y
!q_pic(3) = z
!q_pic(4) = vx
!q_pic(5) = vy
!q_pic(6) = vz
!q_pic(7) = charge
!q_pic(8) = mass

type :: prism_pic_object
   real(R8P)                 :: plasma_density             !< Plasma density
   real(R8P)                 :: neutral_fraction = 0.0_R8P !< Neutral fraction
   integer(I4P)              :: particle_number  = 0_I4P   !< Total number of particles.
	integer(I4P)				  :: n_ions = 0_I4P             !< Total ions number
	integer(I4P)				  :: n_electrons = 0_I4P        !< Total electrons number
	integer(I4P)				  :: n_neutrals = 0_I4P         !< Total neutrals number
   integer(I4P), allocatable :: neighbour_list(:,:)        !< Particle grid positions array.
   character(len=99)         :: problem_type               !< Type of problem analyzed
   character(len=99)         :: particle_weighting_model   !< Particle weighting model.
   character(len=99)         :: current_weighting_model    !< Current weighting model.
   character(len=99)         :: field_weighting_model      !< Field weighting model.
   character(len=99)         :: scheme_time                !< Numerical scheme for time operator [runge_kutta, leapfrog,...].
   !< Pointer (abstract) TBP.
   procedure(particle_weighting_interface), pass(self), pointer :: particle_weighting =>null() !< Particle weighting.
   procedure(current_weighting_interface),  pass(self), pointer :: current_weighting  =>null() !< Current weighting.
   procedure(field_weighting_interface),    pass(self), pointer :: field_weighting    =>null() !< field weighting.
contains
   procedure, pass(self) :: description                   !< Return pretty-printed object description.
   procedure, pass(self) :: initialize                    !< Initialize IC.
   procedure, pass(self) :: load_from_file                !< Load config from file.
   procedure, pass(self) :: particle_cartesian_grid_index !< Compute the grid index corresponding to a particle position.
   procedure, pass(self) :: CIC_charge_weighting          !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_charge_weighting          !< Nearest Grid Point weighting of particle quantities to the grid.
   procedure, pass(self) :: TSC_charge_weighting          !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: CIC_current_weighting         !< Cloud-in-Cell weighting of particle quantities to the grid.
   procedure, pass(self) :: NGP_current_weighting         !< Nearest Grid Point weighting of particle quantities to the grid.
   procedure, pass(self) :: TSC_current_weighting         !< Triangular Shaped Cloud weighting of particle quantities to the grid.
   procedure, pass(self) :: zeroD_field_weighting
   procedure, pass(self) :: oneD_field_weighting
endtype prism_pic_object

interface
   subroutine particle_weighting_interface(self, field, grid, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                             !< External fields.
   type(field_object),                  intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(inout) :: q(1:,1-grid%ngc:,&
                                                              1-grid%ngc:,&
                                                              1-grid%ngc:,1:) !< Field variables.
   real(R8P),                           intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),                        intent(in)    :: nv                   !< Number of variables.
   endsubroutine particle_weighting_interface

   subroutine current_weighting_interface(self, field, grid, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                             !< External fields.
   type(field_object),                  intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(inout) :: q(1:,1-grid%ngc:,&
                                                              1-grid%ngc:,&
                                                              1-grid%ngc:,1:) !< Field variables.
   real(R8P),                           intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),                        intent(in)    :: nv                   !< Number of variables.
   endsubroutine current_weighting_interface

   subroutine field_weighting_interface(self, field, grid, pic_fields, q, q_pic, nv)
   import :: prism_pic_object, grid_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                             !< External fields.
   type(field_object),                  intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),                           intent(inout) :: pic_fields(1:,1:)    !< Fields value at particle locations.
   real(R8P),                           intent(in)    :: q(1:,1-grid%ngc:,&
                                                              1-grid%ngc:,&
                                                              1-grid%ngc:,1:) !< Field variables.
   real(R8P),                           intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),                        intent(in)    :: nv                   !< Number of variables.
   endsubroutine field_weighting_interface
endinterface

contains
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_pic_object), intent(in) :: self             !< External fields.
   character(len=:), allocatable       :: desc             !< Description.
   character(len=1), parameter         :: NL=new_line('a') !< New line character.

   desc =            mpih%myrankstr//'PIC object description:'
   desc = desc//NL//mpih%myrankstr//'    Problem type: '//trim(self%problem_type)
   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      desc = desc//NL//mpih%myrankstr//'    Input plasma density [m^(-3)]: '//trim(str(self%plasma_density))
      desc = desc//NL//mpih%myrankstr//'    Neutral fraction: '//trim(str(self%neutral_fraction))
   endif
   desc = desc//NL//mpih%myrankstr//'    Particle weighting model: '//trim(self%particle_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Current weighting model: '//trim(self%current_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Field weighting model: '//trim(self%field_weighting_model)
   desc = desc//NL//mpih%myrankstr//'    Numerical scheme for time operator: '//trim(self%scheme_time)
   endfunction description

   subroutine initialize(self, field, grid, file_parameters)
   !< Initialize PIC.
   class(prism_pic_object), intent(inout) :: self            !< Pic object.
   type(field_object),      intent(in)    :: field           !< Field (sibling realm component, threaded in).
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   type(file_ini),          intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   real(R8P)                              :: domain_volume   !< Total volume of the computational domain where plasma is
                                                             !< present at t0
	character(len=:),        allocatable   :: desc
   character(len=1),        parameter     :: NL=new_line('a')

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize start'

   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, &
               e_min=>grid%domain_emin, e_max=>grid%domain_emax)

   if (self%problem_type == PLASMA_TYPE_PROBLEM) then
      domain_volume = (e_max(1)-e_min(1))*(e_max(2)-e_min(2))*(e_max(3)-e_min(3))
      self%particle_number = nint(self%plasma_density*domain_volume)
      self%n_neutrals = nint(self%neutral_fraction*real(self%particle_number,R8P))
	   self%n_ions = nint(real(self%particle_number-self%n_neutrals, R8P)/2.0_R8P)
	   self%n_electrons = self%n_ions
	   self%n_neutrals = self%particle_number - self%n_ions - self%n_electrons
   elseif (self%problem_type == SINGLE_PARTICLE_TYPE_PROBLEM) then
      self%particle_number = 1_I4P
   endif

   endassociate

   allocate(self%neighbour_list(4, self%particle_number))

   select case(self%particle_weighting_model)
   case(CIC_WEIGHTING_MODEL)
      self%particle_weighting => CIC_charge_weighting
   case(NGP_WEIGHTING_MODEL)
      self%particle_weighting => NGP_charge_weighting
   case(TSC_WEIGHTING_MODEL)
      self%particle_weighting => TSC_charge_weighting
   case default
      call mpih%error_stop(msg=': invalid particle weighting model in prism_cpu_object%initialize')
   endselect

   select case(self%current_weighting_model)
   case(CIC_WEIGHTING_MODEL)
      self%current_weighting => CIC_current_weighting
   case(NGP_WEIGHTING_MODEL)
      self%current_weighting => NGP_current_weighting
   case(TSC_WEIGHTING_MODEL)
      self%current_weighting => TSC_current_weighting
   case default
      call mpih%error_stop(msg=': invalid current weighting model in prism_cpu_object%initialize')
   endselect

   select case(self%field_weighting_model)
   case(ZEROD_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => zeroD_field_weighting
   case(ONED_FIELDS_WEIGHTING_MODEL)
      self%field_weighting => oneD_field_weighting
   !case(TSC_WEIGHTING_MODEL)
   !   current_weighting => TSC_current_weighting
   case default
      call mpih%error_stop(msg=': invalid field weighting model in prism_cpu_object%initialize')
   endselect

   print '(A)', mpih%myrankstr//'prism_pic_object%initialize finish'
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

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='particle_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(particle_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%particle_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%particle_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%particle_weighting_model = TSC_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid particle weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(particle_weighting_model)')
	endselect

	call file_parameters%get(section_name=INI_SECTION_NAME, option_name='current_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(current_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('CIC', 'cic', 'Cic')
      self%current_weighting_model = CIC_WEIGHTING_MODEL
	case('NGP', 'ngp', 'Ngp')
		self%current_weighting_model = NGP_WEIGHTING_MODEL
	case('TSC', 'tsc', 'Tsc')
		self%current_weighting_model = TSC_WEIGHTING_MODEL
	case default
		call mpih%error_stop(msg=': invalid current weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(current_weighting_model)')
	endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='field_weighting_model', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(field_weighting_model) from file')
   select case(trim(adjustl(buff)))
   case('0D', '0d', '0_d', '0_D')
      self%field_weighting_model = ZEROD_FIELDS_WEIGHTING_MODEL
   case('1D', '1d', '1_d', '1_D')
      self%field_weighting_model = ONED_FIELDS_WEIGHTING_MODEL
   case default
      call mpih%error_stop(msg=': invalid field weighting model ['//trim(adjustl(buff))//'] in  &
      ['//INI_SECTION_NAME//'].(field_weighting_model)')
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='scheme_time', val=buff,error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(scheme_time)')
   select case(trim(adjustl(buff)))
   case('LEAPFROG', 'leapfrog', 'Leapfrog')
      self%scheme_time = NUM_SCHEME_TIME_PIC_LEAPFROG
   case('RUNGE_KUTTA', 'runge_kutta', 'Runge_Kutta')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   case default
      call mpih%print_message(msg='warning: numerical scheme "'//trim(adjustl(buff))//'" unknown. Revert back to RK scheme')
      self%scheme_time = NUM_SCHEME_TIME_PIC_RUNGE_KUTTA
   endselect

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='problem_type', &
   val=self%problem_type, error=error)
   if (.not.go_on_fail_.and.error>0) &
   call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(problem_type)')

   if(self%problem_type == PLASMA_TYPE_PROBLEM) then
	   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='plasma_density', &
      val=self%plasma_density, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(plasma_density)')

      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='neutral_fraction', &
      val=self%neutral_fraction, error=error)
      if (.not.go_on_fail_.and.error>0) &
      call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(neutral_fraction)')
   endif
   endsubroutine load_from_file

   subroutine particle_cartesian_grid_index(self, field, grid, q_pic)
   !< Compute the grid index corresponding to a particle position. Good for cartesian grids only.
   class(prism_pic_object), intent(inout) :: self               !< External fields.
   type(field_object),      intent(in)    :: field              !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(in)    :: q_pic(1:,1:)       !< PIC variables.
   real(R8P)                              :: n                  !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p !< Particle grid indices

   associate(blocks_number=>field%blocks_number, ni=>grid%ni, nj=>grid%nj,                                  &
            nk=>grid%nk, ngc=>grid%ngc, dx => field%dxyz(1,:), dy => field%dxyz(2,:), dz => field%dxyz(3,:),&
            np => self%particle_number, e_min => grid%domain_emin, e_max => grid%domain_emax,               &
            neighbour_list => self%neighbour_list)

   !Va completato considerando la presenza di più blocchi, questo funziona per un blocco solo
   do n = 1, np
      i_p = (q_pic(1,n) - e_min(1)) / dx(1)
      j_p = (q_pic(2,n) - e_min(2)) / dy(1)
      k_p = (q_pic(3,n) - e_min(3)) / dz(1)
      b_p = 1 ! Single block only for now

      neighbour_list(1,n) = ceiling(b_p)
      neighbour_list(2,n) = ceiling(i_p)
      neighbour_list(3,n) = ceiling(j_p)
      neighbour_list(4,n) = ceiling(k_p)
   enddo
   endassociate
   endsubroutine particle_cartesian_grid_index

   subroutine NGP_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: dx, dy, dz           !< Cell dimensions
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   !Per iniziare, azzero tutte le cariche altrimenti vado a sommare le cariche del tempo precedente
   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera
      q(nv, i_p, j_p, k_p, b_p) = q(nv, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)
      !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
      !se scritta in questo modo
   enddo
   endsubroutine NGP_charge_weighting

   subroutine CIC_charge_weighting(self, field, grid, q, q_PIC, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_PIC(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le cariche altrimenti vado a sommare le cariche del tempo precedente
   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.0_R8P) then
                  Wx = 1.0_R8P - abs((q_pic(1,n) - cell_coord(1))/dx)
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.0_R8P) then
                  Wy = 1.0_R8P - abs((q_pic(2,n) - cell_coord(2))/dy)
               else
                  Wy = 0.0_R8P
               end if
               if (abs((q_pic(3,n) - cell_coord(3))/dz) <= 1.0_R8P) then
                  Wz = 1.0_R8P - abs((q_pic(3,n) - cell_coord(3))/dz)
               else
                  Wz = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_pic(7,n)/(dx*dy*dz) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine CIC_charge_weighting

   subroutine TSC_charge_weighting(self, field, grid, q, q_pic, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le cariche altrimenti vado a sommare le cariche del tempo precedente
   q(nv,:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 0.5_R8P) then
                  Wx = 0.75_R8P - ((q_pic(1,n) - cell_coord(1))/dx)**2
               elseif (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.5_R8P .and. abs((q_pic(1,n) - cell_coord(1))/dx) > 0.5_R8P) then
                  Wx = 0.5_R8P * (1.5_R8P - abs((q_pic(1,n) - cell_coord(1))/dx))**2
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 0.5_R8P) then
                  Wy = 0.75_R8P - ((q_pic(2,n) - cell_coord(2))/dy)**2
               elseif (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.5_R8P .and. abs((q_pic(2,n) - cell_coord(2))/dy) > 0.5_R8P) then
                  Wy = 0.5_R8P * (1.5_R8P - abs((q_pic(2,n) - cell_coord(2))/dy))**2
               else
                  Wy = 0.0_R8P
               end if
               q(nv, i, j, k, b_p) = q(nv, i, j, k, b_p) + q_pic(7,n)/(dx*dy*dz) * Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine TSC_charge_weighting

   subroutine NGP_current_weighting(self, field, grid, q, q_pic, nv)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera
      q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)
      q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)
      q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)
      !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
      !se scritta in questo modo
   enddo
   endsubroutine NGP_current_weighting

   subroutine CIC_current_weighting(self, field, grid, q, q_pic, nv)
   !< Cloud-in-Cell weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 1.0_R8P) then
                  Wx = 1.0_R8P - abs((q_pic(1,n) - cell_coord(1))/dx)
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.0_R8P) then
                  Wy = 1.0_R8P - abs((q_pic(2,n) - cell_coord(2))/dy)
               else
                  Wy = 0.0_R8P
               end if
               if (abs((q_pic(3,n) - cell_coord(3))/dz) <= 1.0_R8P) then
                  Wz = 1.0_R8P - abs((q_pic(3,n) - cell_coord(3))/dz)
               else
                  Wz = 0.0_R8P
               end if
               q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)* Wx * Wy * Wz
               q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)* Wx * Wy * Wz
               q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)* Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine CIC_current_weighting

   subroutine TSC_current_weighting(self, field, grid, q, q_pic, nv)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   class(prism_pic_object), intent(inout) :: self                 !< External fields.
   type(field_object),      intent(inout) :: field                !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:) !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)         !< PIC variables.
   integer(I4P),            intent(in)    :: nv                   !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b        !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, b_p   !< Particle grid indices
   real(R8P)                              :: dx, dy, dz           !< Grid spacing
   real(R8P)                              :: wx, wy, wz           !< Weighting factors
   real(R8P)                              :: cell_coord(3)        !< Cell coordinates

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)

   !Per iniziare, azzero tutte le correnti altrimenti vado a sommare le cariche del tempo precedente
   q((nv-3):(nv-1),:,:,:,:) = 0.0_R8P

   do n = 1, self%particle_number
      ! Get particle grid indices
      b_p = self%neighbour_list(1,n)
      i_p = self%neighbour_list(2,n)
      j_p = self%neighbour_list(3,n)
      k_p = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli
      dx = field%dxyz(1,b_p)
      dy = field%dxyz(2,b_p)
      dz = field%dxyz(3,b_p)

      !Qua ci va sicuramente un if per le celle di confine, altrimenti darà errore quando arrivo alla frontiera

      do i = i_p-1, i_p+1
         do j = j_p-1, j_p+1
            do k = k_p-1, k_p+1
               cell_coord = [x_cell(i,b_p), y_cell(j,b_p), z_cell(k,b_p)]
               if (abs((q_pic(1,n) - cell_coord(1))/dx) <= 0.5_R8P) then
                  Wx = 0.75_R8P - ((q_pic(1,n) - cell_coord(1))/dx)**2
               elseif (abs((q_PIC(1,n) - cell_coord(1))/dx) <= 1.5_R8P .and. abs((q_pic(1,n) - cell_coord(1))/dx) > 0.5_R8P) then
                  Wx = 0.5_R8P * (1.5_R8P - abs((q_pic(1,n) - cell_coord(1))/dx))**2
               else
                  Wx = 0.0_R8P
               end if
               if (abs((q_pic(2,n) - cell_coord(2))/dy) <= 0.5_R8P) then
                  Wy = 0.75_R8P - ((q_pic(2,n) - cell_coord(2))/dy)**2
               elseif (abs((q_pic(2,n) - cell_coord(2))/dy) <= 1.5_R8P .and. abs((q_pic(2,n) - cell_coord(2))/dy) > 0.5_R8P) then
                  Wy = 0.5_R8P * (1.5_R8P - abs((q_pic(2,n) - cell_coord(2))/dy))**2
               else
                  Wy = 0.0_R8P
               end if
               q(nv-3, i_p, j_p, k_p, b_p) = q(nv-3, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(4,n)* Wx * Wy * Wz
               q(nv-2, i_p, j_p, k_p, b_p) = q(nv-2, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(5,n)* Wx * Wy * Wz
               q(nv-1, i_p, j_p, k_p, b_p) = q(nv-1, i_p, j_p, k_p, b_p) + q_pic(7,n)/(dx*dy*dz)*q_pic(6,n)* Wx * Wy * Wz

               !Ok, ma va normalizzata e la carica nel vettore di stato va necessariamente azzerata a monte di ogni assegnazione
               !se scritta in questo modo
            enddo
         enddo
      enddo
   enddo
   endassociate
   endsubroutine TSC_current_weighting

   subroutine zeroD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   class(prism_pic_object), intent(inout) :: self                                             !< External fields.
   type(field_object),      intent(inout) :: field                                            !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: pic_fields(1:,1:)                                !< Fields value at particle locations
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:)                             !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)                                     !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                               !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                    !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, block_p                           !< Particle grid indices
   real(R8P)                              :: v_p(3), D_p(3), B_p(3), F_p(3), F_l(3), m_p, q_p !< Particle scalars

   do n = 1, self%particle_number
      ! Get particle grid indices
      block_p = self%neighbour_list(1,n)
      i_p     = self%neighbour_list(2,n)
      j_p     = self%neighbour_list(3,n)
      k_p     = self%neighbour_list(4,n)

      D_p(1) = q(1, i_p, j_p, k_p, block_p)
      D_p(2) = q(2, i_p, j_p, k_p, block_p)
      D_p(3) = q(3, i_p, j_p, k_p, block_p)
      B_p(1) = q(4, i_p, j_p, k_p, block_p)
      B_p(2) = q(5, i_p, j_p, k_p, block_p)
      B_p(3) = q(6, i_p, j_p, k_p, block_p)

      !v_p = [q_pic(n,4), q_pic(n,5), q_pic(n,6)]
      !q_p = q_pic(n,7)
      !m_p = q_pic(n,8)
!
      !F_l = crossproduct(a=v_p, b=B_p)
      !F_p(1) = q_p*(D_p(1) + F_l(1))
      !F_p(2) = q_p*(D_p(2) + F_l(2))
      !F_p(3) = q_p*(D_p(3) + F_l(3))

      pic_fields(1,n) = D_p(1)
      pic_fields(2,n) = D_p(2)
      pic_fields(3,n) = D_p(3)
      pic_fields(4,n) = B_p(1)
      pic_fields(5,n) = B_p(2)
      pic_fields(6,n) = B_p(3)

   enddo
   endsubroutine zeroD_field_weighting

   subroutine oneD_field_weighting(self, field, grid, pic_fields, q, q_pic, nv)
   class(prism_pic_object), intent(inout) :: self                                             !< External fields.
   type(field_object),      intent(inout) :: field                                            !< The field.
   type(grid_object),                  intent(in)              :: grid
                                                                                      !< Grid (sibling realm component, threaded
                                                                                      !< in).
   real(R8P),               intent(inout) :: pic_fields(1:,1:)                                !< Fields value at particle locations
   real(R8P),               intent(in)    :: q(1:,1-grid%ngc:,&
                                                  1-grid%ngc:,&
                                                  1-grid%ngc:,1:)                             !< Field variables.
   real(R8P),               intent(in)    :: q_pic(1:,1:)                                     !< PIC variables.
   integer(I4P),            intent(in)    :: nv                                               !< Number of variables.
   real(R8P)                              :: n, i, j, k ,b                                    !< Particle counter
   real(R8P)                              :: i_p, j_p, k_p, block_p                           !< Particle grid indices
   real(R8P)                              :: x_p, y_p, z_p                                    !< Particle position scalar
   real(R8P)                              :: v_p(3), D_p(3), B_p(3), F_p(3), F_l(3), m_p, q_p !< Particle scalars
   real(R8P)                              :: dx, dy, dz                                       !< Grid spacing
   real(R8P)                              :: dDx_dx, dDx_dy, dDx_dz
   real(R8P)                              :: dDy_dx, dDy_dy, dDy_dz
   real(R8P)                              :: dDz_dx, dDz_dy, dDz_dz
   real(R8P)                              :: dBx_dx, dBx_dy, dBx_dz
   real(R8P)                              :: dBy_dx, dBy_dy, dBy_dz
   real(R8P)                              :: dBz_dx, dBz_dy, dBz_dz

   associate(x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   do n = 1, self%particle_number

      ! Get particle grid indices
      block_p = self%neighbour_list(1,n)
      i_p     = self%neighbour_list(2,n)
      j_p     = self%neighbour_list(3,n)
      k_p     = self%neighbour_list(4,n)

      ! Qua va capito come gestire la questione dei blocchi multipli

      !Interpolazione lineare dei campi nella posizione delle particelle
      !x
      if (x_cell(i_p,block_p) >= q_pic(1,n)) then !La particella è nella metà sinistra della cella
         dDx_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(1,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDy_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(2,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDz_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(3,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBx_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(4,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBy_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(5,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBz_dx = lininterp(x2=x_cell(i_p-1,block_p), y2=q(6,i_p-1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
      else
         dDx_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(1,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDy_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(2,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dDz_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(3,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBx_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(4,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBy_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(5,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
         dBz_dx = lininterp(x2=x_cell(i_p+1,block_p), y2=q(6,i_p+1,j_p,k_p,block_p), &
                            x1=x_cell(i_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(1,n))
      endif
      !y
      if (y_cell(j_p,block_p) >= q_pic(2,n)) then !La particella è nella metà sinistra della cella
         dDx_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(1,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDy_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(2,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDz_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(3,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBx_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(4,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBy_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(5,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBz_dy = lininterp(x2=y_cell(j_p-1,block_p), y2=q(6,i_p,j_p-1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
      else
         dDx_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(1,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDy_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(2,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dDz_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(3,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBx_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(4,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBy_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(5,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
         dBz_dy = lininterp(x2=y_cell(j_p+1,block_p), y2=q(6,i_p,j_p+1,k_p,block_p), &
                            x1=y_cell(j_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(2,n))
      endif
      !z
      if (z_cell(k_p,block_p) >= q_pic(3,n)) then !La particella è nella metà sinistra della cella
         dDx_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(1,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDy_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(2,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDz_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(3,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBx_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(4,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBy_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(5,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBz_dz = lininterp(x2=z_cell(k_p-1,block_p), y2=q(6,i_p,j_p,k_p-1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
      else
         dDx_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(1,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(1,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDy_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(2,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(2,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dDz_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(3,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(3,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBx_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(4,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(4,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBy_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(5,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(5,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
         dBz_dz = lininterp(x2=z_cell(k_p+1,block_p), y2=q(6,i_p,j_p,k_p+1,block_p), &
                            x1=z_cell(k_p,block_p), y1=q(6,i_p,j_p,k_p,block_p), xp=q_pic(3,n))
      endif

      D_p(1) = q(1,i_p,j_p,k_p,block_p) + dDx_dx + dDx_dy + dDx_dz
      D_p(2) = q(2,i_p,j_p,k_p,block_p) + dDy_dx + dDy_dy + dDy_dz
      D_p(3) = q(3,i_p,j_p,k_p,block_p) + dDz_dx + dDz_dy + dDz_dz
      B_p(1) = q(4,i_p,j_p,k_p,block_p) + dBx_dx + dBx_dy + dBx_dz
      B_p(2) = q(5,i_p,j_p,k_p,block_p) + dBy_dx + dBy_dy + dBy_dz
      B_p(3) = q(6,i_p,j_p,k_p,block_p) + dBz_dx + dBz_dy + dBz_dz

      !v_p = [q_pic(n,4), q_pic(n,5), q_pic(n,6)]
      !q_p = q_PIC(n,7)
      !m_p = q_PIC(n,8)
      !F_l = crossproduct(a=v_p, b=B_p)
      !F_p(1) = q_p*(D_p(1) + F_l(1))
      !F_p(2) = q_p*(D_p(2) + F_l(2))
      !F_p(3) = q_p*(D_p(3) + F_l(3))

      pic_fields(1,n) = D_p(1)
      pic_fields(2,n) = D_p(2)
      pic_fields(3,n) = D_p(3)
      pic_fields(4,n) = B_p(1)
      pic_fields(5,n) = B_p(2)
      pic_fields(6,n) = B_p(3)

   enddo
   endassociate
   endsubroutine oneD_field_weighting

   function crossproduct(a, b) result(cross)
   real(R8P), intent(in) :: a(3)     !< Left hand side.
   real(R8P), intent(in) :: b(3)     !< Left hand side.
   real(R8P)             :: cross(3) !< Cross product.

   cross(1) = (a(2) * b(3)) - (a(3) * b(2))
   cross(2) = (a(3) * b(1)) - (a(1) * b(3))
   cross(3) = (a(1) * b(2)) - (a(2) * b(1))
   endfunction crossproduct

   function lininterp(x1, y1, x2, y2, xp) result(delta)
   real(R8P), intent(in) :: x1, x2
   real(R8P), intent(in) :: y1, y2
   real(R8P), intent(in) :: xp
   real(R8P)             :: delta
   real(R8P)             :: m

   m = (y2-y1)/(x2-x1)
   delta = m*(xp-x1)
   endfunction lininterp
endmodule adam_prism_pic_object
