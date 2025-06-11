!< ADAM, Maxwell IO handler class definition, CPU backend.
module adam_prism_io_object
!< ADAM, Maxwell IO handler class definition, CPU backend.

use adam_mpih_object
use finer
use penf

implicit none
private
public :: prism_io_object

character(len=2), parameter :: INI_SECTION_NAME="IO" !< INI (config) file section name containing IO configs.

type :: prism_io_object
   !< PRISM IO handler class definition, CPU backend.
   type(mpih_object)         :: mpih                       !< MPI handler.
   type(file_ini)            :: file_parameters            !< Prism input file handler.
   integer(I4P)              :: it_save=100_I4P            !< Main output iteration save frequency.
   character(:), allocatable :: output_basename            !< Basename of output files.
   logical                   :: restart=.false.            !< Enable restart from old output data.
   character(:), allocatable :: restart_basename           !< Basename of restart files.
   integer(I4P)              :: restart_save=100_I4P       !< Restart output iteration save frequency.
   logical                   :: save_memory_status=.false. !< Enable save of memory status during allocations.
   integer(I4P)              :: residuals_save=10_I4P      !< Residuals (norm) output iteration save frequency.
   integer(I4P)              :: residuals_unit             !< Residuals file unit.
   contains
      procedure, pass(self) :: description          !< Return pretty-printed object description.
      procedure, pass(self) :: initialize           !< Initialize time handler.
      procedure, pass(self) :: load_from_file       !< Load config from file.
      ! residuals IO
      procedure, pass(self) :: close_file_residuals !< Close file for saving residuals history.
      procedure, pass(self) :: open_file_residuals  !< Open file for saving residuals history.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
endtype prism_io_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_io_object), intent(in) :: self             !< IO handler.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'IO main data'                                               //NL
   desc = desc//self%mpih%myrankstr//'  file parameters:     '//self%file_parameters%filename     //NL
   desc = desc//self%mpih%myrankstr//'  it save:             '//trim(str(self%it_save))           //NL
   desc = desc//self%mpih%myrankstr//'  output basename:     '//self%output_basename              //NL
   desc = desc//self%mpih%myrankstr//'  restart:             '//trim(str(self%restart))           //NL
   desc = desc//self%mpih%myrankstr//'  restart basename:    '//self%restart_basename             //NL
   desc = desc//self%mpih%myrankstr//'  restart (it) save:   '//trim(str(self%restart_save))      //NL
   desc = desc//self%mpih%myrankstr//'  save memory status:  '//trim(str(self%save_memory_status))//NL
   desc = desc//self%mpih%myrankstr//'  residuals (it) save: '//trim(str(self%residuals_save))
   endfunction description

   subroutine initialize(self, filename)
   !< Initialize IO handler.
   class(prism_io_object), intent(inout) :: self     !< IO handler.
   character(*),           intent(in)    :: filename !< File name of parameters file.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_io_object%initialize start'
   call self%file_parameters%initialize(filename=trim(filename))
   call self%file_parameters%load
   call self%load_from_file(file_parameters=self%file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_io_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_io_object), intent(inout)        :: self            !< IO handler.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(99)                                :: buff_c          !< Character buffer.
   integer(I4P)                                 :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='output_basename', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(output_basename)')
   self%output_basename = trim(adjustl(buff_c))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='it_save', val=self%it_save, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(it_save)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='restart', val=self%restart, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(restart)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='restart_basename', val=buff_c, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(restart_basename)')
   self%restart_basename = trim(adjustl(buff_c))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='restart_save', val=self%restart_save, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(restart_save)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='residuals_save', val=self%residuals_save, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(residuals_save)')
   call file_parameters%get(section_name=INI_SECTION_NAME,option_name='save_memory_status',val=self%save_memory_status,error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(save_memory_status)')
   endsubroutine load_from_file

      ! residuals IO
   subroutine  close_file_residuals(self)
   !< Close file for saving residuals history.
   class(prism_io_object), intent(in) :: self !< IO handler.

   close(self%residuals_unit)
   endsubroutine close_file_residuals

   subroutine open_file_residuals(self, nv)
   !< Open file for saving residuals history.
   class(prism_io_object), intent(inout) :: self !< IO handler.
   integer(I4P),           intent(in)    :: nv   !< Number of residuals variables.
   character(:), allocatable             :: rqs  !< String buffer.
   integer(I4P)                          :: v    !< Counter.

   rqs = ''
   do v=1, nv
      rqs = rqs//' "rq'//trim(str(v,.true.))//'"'
   enddo
   open(newunit=self%residuals_unit, file=self%output_basename//'-residuals.dat')
   write(self%residuals_unit, '(A)') 'VARIABLES="it" "time" "blocks_number"'//rqs
   endsubroutine open_file_residuals

   subroutine save_residuals(self, it, time, blocks_number, residuals)
   !< Save residuals history.
   class(prism_io_object), intent(in) :: self          !< IO handler.
   integer(I4P),           intent(in) :: it            !< Current iteration.
   real(R8P),              intent(in) :: time          !< Current time.
   integer(I4P),           intent(in) :: blocks_number !< Current number of blocks.
   real(R8P),              intent(in) :: residuals(1:) !< Residuals (norm) [1:nv].

   write(self%residuals_unit, '(A)') trim(str(it           ))//' '//&
                                     trim(str(time         ))//' '//&
                                     trim(str(blocks_number))//' '//&
                                     trim(str(residuals(1:), separator=' '))
   endsubroutine save_residuals
endmodule adam_prism_io_object
    