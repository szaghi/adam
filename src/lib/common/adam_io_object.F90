!< ADAM, IO class definition.
module adam_io_object
!< ADAM, IO class definition.

! ADAM modules
use :: adam_field_object
use :: adam_grid_object
use :: adam_mpih_object
use :: finer
! third party modules
use :: finer
use :: motion
use :: penf
use :: stringifor

implicit none
private
public :: io_object

character(len=2), parameter :: INI_SECTION_NAME="IO" !< INI (config) file section name containing IO configs.

type :: io_object
   !< ADAM class definition.
   type(mpih_object)           :: mpih                           !< The MPI handler.
   type(file_ini)              :: file_parameters                !< Input file handler.
   integer(I4P)                :: it_save=100_I4P                !< Main output iteration save frequency.
   character(:), allocatable   :: output_basename                !< Basename of output files.
   logical                     :: restart=.false.                !< Enable restart from old output data.
   character(:), allocatable   :: restart_basename               !< Basename of restart files.
   integer(I4P)                :: restart_save=100_I4P           !< Restart output iteration save frequency.
   logical                     :: save_memory_status=.false.     !< Enable save of memory status during allocations.
   integer(I4P)                :: residuals_save=10_I4P          !< Residuals (norm) output iteration save frequency.
   integer(I4P)                :: residuals_unit                 !< Residuals file unit.
   integer(I4P)                :: energy_error_save=10_I4P       !< Energy error output iteration save frequency.
   integer(I4P)                :: energy_error_unit              !< Energy error hystory file unit.
   integer(I4P)                :: energy_history_save=10_I4P     !< Energy history output iteration save frequency.
   integer(I4P)                :: energy_history_unit            !< Energy history file unit.
   integer(I4P)                :: divergence_history_save=10_I4P !< Divergence history output iteration save frequency.
   integer(I4P)                :: divergence_history_unit        !< Divergence history file unit.
   type(grid_object),  pointer :: grid=>null()                   !< The grid.
   type(field_object), pointer :: field=>null()                  !< The field.
   logical                     :: is_initialized=.false.         !< Initialization status.
   ! auxiliary fields saving
   logical :: save_residual_fields  =.false. !< Flag to activate residual fields saving.
   logical :: save_curl_fields      =.false. !< Flag to activate curl fields saving.
   logical :: save_divergence_fields=.false. !< Flag to activate divergence fields saving.
   logical :: save_gradient_fields  =.false. !< Flag to activate gradient fields saving.
   logical :: save_laplacian_fields =.false. !< Flag to activate gradient fields saving.
   contains
      ! public methods
      procedure, pass(self) :: associate_grid_field !< Associate grid and field class.
      procedure, pass(self) :: description          !< Return pretty-printed object description.
      procedure, pass(self) :: initialize           !< Initialize class.
      procedure, pass(self) :: load_from_file       !< Load config from file.
      ! auxiliary fields IO
      procedure, pass(self) :: save_energy_error       !< Save energy error history.
      procedure, pass(self) :: save_energy_history     !< Save energy history.
      procedure, pass(self) :: save_divergence_history !< Save divergence history.
      ! residuals IO
      procedure, pass(self) :: close_file_residuals !< Close file for saving residuals history.
      procedure, pass(self) :: open_file_residuals  !< Open file for saving residuals history.
      procedure, pass(self) :: save_residuals       !< Save residuals history.
      ! XH5F IO
      generic               :: save_field =>           &
                               save_xh5f_field_4D_R8P, &
                               save_xh5f_field_4D_R4P, &
                               save_xh5f_field_4D_I8P, &
                               save_xh5f_field_4D_I4P, &
                               save_xh5f_field_4D_I2P, &
                               save_xh5f_field_4D_I1P, &
                               save_xh5f_field_5D_R8P, &
                               save_xh5f_field_5D_R4P, &
                               save_xh5f_field_5D_I8P, &
                               save_xh5f_field_5D_I4P, &
                               save_xh5f_field_5D_I2P, &
                               save_xh5f_field_5D_I1P !< Save fields by XH5F (XDMF/HDF5) file handler.
      ! private methods
      procedure, pass(self), private :: save_xh5f_field_4D_R8P !< Save fields by XH5F file handler, rank 4D, kind R8P.
      procedure, pass(self), private :: save_xh5f_field_4D_R4P !< Save fields by XH5F file handler, rank 4D, kind R4P.
      procedure, pass(self), private :: save_xh5f_field_4D_I8P !< Save fields by XH5F file handler, rank 4D, kind I8P.
      procedure, pass(self), private :: save_xh5f_field_4D_I4P !< Save fields by XH5F file handler, rank 4D, kind I4P.
      procedure, pass(self), private :: save_xh5f_field_4D_I2P !< Save fields by XH5F file handler, rank 4D, kind I2P.
      procedure, pass(self), private :: save_xh5f_field_4D_I1P !< Save fields by XH5F file handler, rank 4D, kind I1P.
      procedure, pass(self), private :: save_xh5f_field_5D_R8P !< Save fields by XH5F file handler, rank 5D, kind R8P.
      procedure, pass(self), private :: save_xh5f_field_5D_R4P !< Save fields by XH5F file handler, rank 5D, kind R4P.
      procedure, pass(self), private :: save_xh5f_field_5D_I8P !< Save fields by XH5F file handler, rank 5D, kind I8P.
      procedure, pass(self), private :: save_xh5f_field_5D_I4P !< Save fields by XH5F file handler, rank 5D, kind I4P.
      procedure, pass(self), private :: save_xh5f_field_5D_I2P !< Save fields by XH5F file handler, rank 5D, kind I2P.
      procedure, pass(self), private :: save_xh5f_field_5D_I1P !< Save fields by XH5F file handler, rank 5D, kind I1P.
endtype io_object
contains
   ! public methods
   subroutine associate_grid_field(self, grid, field)
   !< Associate grid and field class.
   class(io_object),   intent(inout)      :: self  !< IO handler.
   type(grid_object),  intent(in), target :: grid  !< The grid.
   type(field_object), intent(in), target :: field !< The field.

   self%grid  => grid
   self%field => field
   endsubroutine associate_grid_field

   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(io_object), intent(in)  :: self             !< IO handler.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.

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

   subroutine initialize(self, filename, verbose)
   !< Initialize class.
   class(io_object),   intent(inout)        :: self     !< IO handler.
   character(*),       intent(in)           :: filename !< File name of parameters file.
   logical,            intent(in), optional :: verbose  !< Trigger verbose output.
   logical                                  :: verbose_ !< Trigger verbose output, local variable.

   if (.not.self%is_initialized) then
      verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
      call self%mpih%initialize(verbose=verbose_)
      print '(A)', self%mpih%myrankstr//'io_object%initialize start'
      call self%file_parameters%initialize(filename=trim(filename))
      call self%file_parameters%load
      call self%load_from_file(file_parameters=self%file_parameters)
      self%is_initialized = .true.
      print '(A)', self%description()
      print '(A)', self%mpih%myrankstr//'io_object%initialize finish'
   endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(io_object), intent(inout)        :: self            !< IO handler.
   type(file_ini),   intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,          intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                :: go_on_fail_     !< Go on if load fails.
   character(99)                          :: buff_c          !< Character buffer.
   integer(I4P)                           :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   associate(ISN=>INI_SECTION_NAME)
      call file_parameters%get(section_name=ISN, option_name='output_basename', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(output_basename)')
      self%output_basename = trim(adjustl(buff_c))
      call file_parameters%get(section_name=ISN, option_name='it_save', val=self%it_save, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(it_save)')
      call file_parameters%get(section_name=ISN, option_name='restart', val=self%restart, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(restart)')
      call file_parameters%get(section_name=ISN, option_name='restart_basename', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(restart_basename)')
      self%restart_basename = trim(adjustl(buff_c))
      call file_parameters%get(section_name=ISN, option_name='restart_save', val=self%restart_save, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(restart_save)')
      call file_parameters%get(section_name=ISN, option_name='residuals_save', val=self%residuals_save, error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(residuals_save)')
      call file_parameters%get(section_name=ISN,option_name='save_memory_status',val=self%save_memory_status,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_memory_status)')
      call file_parameters%get(section_name=ISN,option_name='save_residual_fields',val=self%save_residual_fields,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_residual_fields)')
      call file_parameters%get(section_name=ISN,option_name='save_curl_fields',val=self%save_curl_fields,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_curl_fields)')
      call file_parameters%get(section_name=ISN,option_name='save_divergence_fields',val=self%save_divergence_fields,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_divergence_fields)')
      call file_parameters%get(section_name=ISN,option_name='save_gradient_fields',val=self%save_gradient_fields,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_gradient_fields)')
      call file_parameters%get(section_name=ISN,option_name='save_laplacian_fields',val=self%save_laplacian_fields,error=error)
      if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//ISN//'].(save_laplacian_fields)')
   endassociate
   endsubroutine load_from_file

   ! auxiliary fields IO
   subroutine save_energy_error(self,it,time,blocks_number,energy_D,energy_B,rms_energy_error_D,rms_energy_error_B,&
                                is_to_open,is_to_close)
   !< Save energy error history.
   class(io_object), intent(inout)        :: self               !< IO handler.
   integer(I4P),     intent(in)           :: it                 !< Current iteration.
   real(R8P),        intent(in)           :: time               !< Current time.
   integer(I4P),     intent(in)           :: blocks_number      !< Current number of blocks.
   real(R8P),        intent(in)           :: energy_D(1:)       !< Energy history of D.
   real(R8P),        intent(in)           :: energy_B(1:)       !< Energy history of B.
   real(R8P),        intent(in)           :: rms_energy_error_D !< RMS of energy history of D.
   real(R8P),        intent(in)           :: rms_energy_error_B !< RMS of energy history of B.
   logical,          intent(in), optional :: is_to_open         !< Flag to open  file before first saving.
   logical,          intent(in), optional :: is_to_close        !< Flag to close file after last saving.
   logical                                :: is_to_open_        !< Flag to open  file before first saving, local var.
   logical                                :: is_to_close_       !< Flag to close file after last saving, local var.

   if (self%mpih%myrank==0) then
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      if (is_to_open_) then
         open(newunit=self%energy_error_unit, file=self%output_basename//'-energy_error.dat')
         write(self%energy_error_unit,'(A)')'VARIABLES="it" "blocks_number" "time" "error_D" "error_B" "rms_error_D" "rms_error_B"'
      endif
      write(self%energy_error_unit, '(A)') trim(str(it                                                  ))//' '//&
                                           trim(str(blocks_number                                       ))//' '//&
                                           trim(str(time                                                ))//' '//&
                                           trim(str(sqrt(abs(energy_D(it)-energy_D(1))/abs(energy_D(1)))))//' '//&
                                           trim(str(sqrt(abs(energy_B(it)-energy_B(1))/abs(energy_B(1)))))//' '//&
                                           trim(str(rms_energy_error_D                                  ))//' '//&
                                           trim(str(rms_energy_error_B                                  ))
      if (is_to_close_) close(self%energy_error_unit)
   endif
   endsubroutine save_energy_error

   subroutine save_energy_history(self,it,time,blocks_number,energy_D,energy_B,coil_power,Poynting_flux, &
                                  is_to_open,is_to_close)
   !< Save energy history.
   class(io_object), intent(inout)        :: self              !< IO handler.
   integer(I4P),     intent(in)           :: it                !< Current iteration.
   real(R8P),        intent(in)           :: time              !< Current time.
   integer(I4P),     intent(in)           :: blocks_number     !< Current number of blocks.
   real(R8P),        intent(in)           :: energy_D(1:)      !< Energy history of D.
   real(R8P),        intent(in)           :: energy_B(1:)      !< Energy history of B.
   real(R8P),        intent(in)           :: coil_power(1:)    !< Coil power history.
   real(R8P),        intent(in)           :: Poynting_flux(1:) !< Poynting flux history.
   logical,          intent(in), optional :: is_to_open        !< Flag to open  file before first saving.
   logical,          intent(in), optional :: is_to_close       !< Flag to close file after last saving.
   logical                                :: is_to_open_       !< Flag to open  file before first saving, local var.
   logical                                :: is_to_close_      !< Flag to close file after last saving, local var.

   if (self%mpih%myrank==0) then
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      if (is_to_open_) then
         open(newunit=self%energy_history_unit, file=self%output_basename//'-energy_history.dat')
         write(self%energy_history_unit,'(A)')&
               '%VARIABLES="it" "blocks_number" "time" "D_energy [J]" "B_energy [J]" "coil_power [W]" "Poynting_flux [W]"'
      endif
      write(self%energy_history_unit, '(A)') trim(str(it               ))//' '//&
                                             trim(str(blocks_number    ))//' '//&
                                             trim(str(time             ))//' '//&
                                             trim(str(energy_D(it)     ))//' '//&
                                             trim(str(energy_B(it)     ))//' '//&
                                             trim(str(coil_power(it)   ))//' '//&
                                             trim(str(Poynting_flux(it)))
      if (is_to_close_) close(self%energy_history_unit)
   endif
   endsubroutine save_energy_history

   subroutine save_divergence_history(self,it,time,blocks_number,div_D,div_B,div_J,is_to_open,is_to_close)
   !< Save energy history.
   class(io_object), intent(inout)        :: self              !< IO handler.
   integer(I4P),     intent(in)           :: it                !< Current iteration.
   real(R8P),        intent(in)           :: time              !< Current time.
   integer(I4P),     intent(in)           :: blocks_number     !< Current number of blocks.
   real(R8P),        intent(in)           :: div_D             !< Energy history of B.
   real(R8P),        intent(in)           :: div_B             !< Coil power history.
   real(R8P),        intent(in)           :: div_J             !< Poynting flux history.
   logical,          intent(in), optional :: is_to_open        !< Flag to open  file before first saving.
   logical,          intent(in), optional :: is_to_close       !< Flag to close file after last saving.
   logical                                :: is_to_open_       !< Flag to open  file before first saving, local var.
   logical                                :: is_to_close_      !< Flag to close file after last saving, local var.

   if (self%mpih%myrank==0) then
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      if (is_to_open_) then
         open(newunit=self%divergence_history_unit, file=self%output_basename//'-divergence_history.dat')
         write(self%divergence_history_unit,'(A)')&
               '%VARIABLES="it" "blocks_number" "time" "D_divergence" "B_divergence" "J_divergence"'
      endif
      write(self%divergence_history_unit, '(A)') trim(str(it               ))//' '//&
                                                 trim(str(blocks_number    ))//' '//&
                                                 trim(str(time             ))//' '//&
                                                 trim(str(div_D            ))//' '//&
                                                 trim(str(div_B            ))//' '//&
                                                 trim(str(div_J            ))
      if (is_to_close_) close(self%divergence_history_unit)
   endif
   endsubroutine save_divergence_history

   ! residuals IO
   subroutine  close_file_residuals(self)
   !< Close file for saving residuals history.
   class(io_object), intent(in) :: self !< IO handler.

   if (self%mpih%myrank==0) close(self%residuals_unit)
   endsubroutine close_file_residuals

   subroutine open_file_residuals(self, nv)
   !< Open file for saving residuals history.
   class(io_object), intent(inout) :: self !< IO handler.
   integer(I4P),     intent(in)    :: nv   !< Number of residuals variables.
   character(:), allocatable       :: rqs  !< String buffer.
   integer(I4P)                    :: v    !< Counter.

   if (self%mpih%myrank==0) then
      rqs = ''
      do v=1, nv
         rqs = rqs//' "rq'//trim(str(v,.true.))//'"'
      enddo
      open(newunit=self%residuals_unit, file=self%output_basename//'-residuals.dat')
      write(self%residuals_unit, '(A)') 'VARIABLES="it" "time" "blocks_number"'//rqs
   endif
   endsubroutine open_file_residuals

   subroutine save_residuals(self, it, time, blocks_number, residuals)
   !< Save residuals history.
   class(io_object), intent(in) :: self          !< IO handler.
   integer(I4P),     intent(in) :: it            !< Current iteration.
   real(R8P),        intent(in) :: time          !< Current time.
   integer(I4P),     intent(in) :: blocks_number !< Current number of blocks.
   real(R8P),        intent(in) :: residuals(1:) !< Residuals (norm) [1:nv].

   write(self%residuals_unit, '(A)') trim(str(it           ))//' '//&
                                     trim(str(time         ))//' '//&
                                     trim(str(blocks_number))//' '//&
                                     trim(str(residuals(1:), separator=' '))
   endsubroutine save_residuals

   ! XH5F IO
   ! private methods
#define KKP R8P
#define VARTYPE real
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_R8P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_R8P
#include "adam_io_agnostic.INC"

#define KKP R4P
#define VARTYPE real
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_R4P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_R4P
#include "adam_io_agnostic.INC"

#define KKP I8P
#define VARTYPE integer
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I8P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I8P
#include "adam_io_agnostic.INC"

#define KKP I4P
#define VARTYPE integer
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I4P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I4P
#include "adam_io_agnostic.INC"

#define KKP I2P
#define VARTYPE integer
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I2P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I2P
#include "adam_io_agnostic.INC"

#define KKP I1P
#define VARTYPE integer
#define SAVE_XH5F_FIELD_4D_KKP save_xh5f_field_4D_I1P
#define SAVE_XH5F_FIELD_5D_KKP save_xh5f_field_5D_I1P
#include "adam_io_agnostic.INC"
endmodule adam_io_object
