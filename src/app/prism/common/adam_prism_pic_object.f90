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
!public :: CIC_weighting
!public :: NGP_weighting
!public :: TSC_weighting

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
   type(mpih_object)     :: mpih                                  !< MPI handler.
   integer(I4P)          :: particle_number = 0.0_R8P             !< Total number of particles.
   character(len=99)     :: particle_weighting_model              !< Particle weighting model.
contains
   procedure, pass(self) :: description      !< Return pretty-printed object description.
   procedure, pass(self) :: initialize       !< Initialize IC.
   procedure, pass(self) :: load_from_file   !< Load config from file.
   !procedure, pass(self) :: CIC_weighting    !< Cloud-in-Cell weighting of particle quantities to the grid.
   !procedure, pass(self) :: NGP_weighting    !< Nearest Grid Point weighting of particle quantities to the grid.
   !procedure, pass(self) :: TSC_weighting    !< Triangular Shaped Cloud weighting of particle quantities to the grid.
endtype prism_pic_object

interface
   subroutine particle_weighting_interface(self, field, q, q_pic)
   import :: prism_pic_object, field_object, I4P, R8P
   class(prism_pic_object), intent(inout) :: self                                                              !< External fields.
   type(field_object),                  intent(inout) :: field                                                             !< The field.
   real(R8P),                           intent(inout) :: q(1:, 1-field%grid%ngc:,1-field%grid%ngc:,1-field%grid%ngc:,1:)   !< Field variables.
   real(R8P),                           intent(inout) :: q_pic(1:,1:)   !< PIC variables.
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

   !subroutine CIC_weighting(self, field, q, q_PIC)
   !!< Cloud-in-Cell weighting of particle quantities to the grid.
   !endsubroutine CIC_weighting

   !subroutine NGP_weighting(self, field, q, q_PIC)
   !!< Nearest Grid Point weighting of particle quantities to the grid.
   !endsubroutine NGP_weighting

   !subroutine TSC_weighting(self, field, q, q_PIC)
   !!< Triangular Shaped Cloud weighting of particle quantities to the grid.
   !endsubroutine TSC_weighting

endmodule adam_prism_pic_object