!< ADAM, IO class definition.
module adam_io_object
!< ADAM, IO class definition.

! ADAM singleton objects
use :: adam_ini_library, only : ini_check, ini_get, ini_load, ini_stop
use :: adam_mpih_global, only : mpih
use :: adam_grid_object, only : grid_object
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
   type(file_ini)              :: file_parameters                !< Input file handler.
   integer(I4P)                :: it_save=100_I4P                !< Main output iteration save frequency.
   character(:), allocatable   :: output_basename                !< Basename of output files.
   logical                     :: restart=.false.                !< Enable restart from old output data.
   logical                     :: unused_keys_error=.true.       !< [IO] unused_keys = error (default): an unread key stops the run.
   character(:), allocatable   :: restart_basename               !< Basename of restart files.
   integer(I4P)                :: restart_save=100_I4P           !< Restart output iteration save frequency.
   integer(I4P)                :: residuals_save=10_I4P          !< Residuals (norm) output iteration save frequency.
   integer(I4P)                :: residuals_unit                 !< Residuals file unit.
   integer(I4P)                :: energy_error_save=10_I4P       !< Energy error output iteration save frequency.
   integer(I4P)                :: energy_error_unit              !< Energy error hystory file unit.
   integer(I4P)                :: energy_history_save=10_I4P     !< Energy history output iteration save frequency.
   integer(I4P)                :: energy_history_unit            !< Energy history file unit.
   integer(I4P)                :: divergence_history_unit        !< Divergence history file unit.
   logical                     :: is_initialized=.false.         !< Initialization status.
   integer(I4P) :: charge_conservation_history_unit = 0_I4P
   ! auxiliary fields saving
   logical :: save_residual_fields  =.false. !< Flag to activate residual fields saving.
   contains
      ! public methods
      procedure, pass(self) :: description          !< Return pretty-printed object description.
      procedure, pass(self) :: initialize           !< Initialize class.
      procedure, pass(self) :: load_from_file       !< Load config from file.
      ! auxiliary fields IO
      procedure, pass(self) :: save_energy_error       !< Save energy error history.
      procedure, pass(self) :: save_energy_history     !< Save energy history.
      procedure, pass(self) :: save_divergence_history !< Save divergence history.
      procedure, pass(self) :: save_charge_conservation_history
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
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(io_object), intent(in)  :: self             !< IO handler.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'IO main data'                                               //NL
   desc = desc//mpih%myrankstr//'  file parameters:     '//self%file_parameters%filename     //NL
   desc = desc//mpih%myrankstr//'  it save:             '//trim(str(self%it_save))           //NL
   desc = desc//mpih%myrankstr//'  output basename:     '//self%output_basename              //NL
   desc = desc//mpih%myrankstr//'  restart:             '//trim(str(self%restart))           //NL
   desc = desc//mpih%myrankstr//'  restart basename:    '//self%restart_basename             //NL
   desc = desc//mpih%myrankstr//'  restart (it) save:   '//trim(str(self%restart_save))      //NL
   desc = desc//mpih%myrankstr//'  residuals (it) save: '//trim(str(self%residuals_save))
   endfunction description

   subroutine initialize(self, filename, verbose)
   !< Initialize class.
   class(io_object),   intent(inout)        :: self     !< IO handler.
   character(*),       intent(in)           :: filename !< File name of parameters file.
   logical,            intent(in), optional :: verbose  !< Trigger verbose output.
   logical                                  :: verbose_ !< Trigger verbose output, local variable.

   if (.not.self%is_initialized) then
      verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
      print '(A)', mpih%myrankstr//'io_object%initialize start'
      call ini_load(filename=trim(filename), ini=self%file_parameters) ! resolves "extends = <base file>" (#80 P4)
      call self%load_from_file(file_parameters=self%file_parameters)
      self%is_initialized = .true.
      print '(A)', self%description()
      print '(A)', mpih%myrankstr//'io_object%initialize finish'
   endif
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file.
   class(io_object), intent(inout)        :: self            !< IO handler.
   type(file_ini),   intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   character(99)                          :: buff_c          !< Character buffer.
   character(16)                          :: env_value       !< ADAM_UNUSED_KEYS value.
   integer                                :: env_length      !< ADAM_UNUSED_KEYS length.
   integer                                :: env_status      !< ADAM_UNUSED_KEYS query status.
   integer(I4P)                           :: error           !< Error status.

   associate(ISN=>INI_SECTION_NAME)
      call ini_get(file_parameters, ISN, 'output_basename', buff_c)
      self%output_basename = trim(adjustl(buff_c))
      call ini_get(file_parameters, ISN, 'it_save', self%it_save)
      call ini_get(file_parameters, ISN, 'restart', self%restart, default=.false.)
      call ini_get(file_parameters, ISN, 'restart_basename', buff_c)
      self%restart_basename = trim(adjustl(buff_c))
      call ini_get(file_parameters, ISN, 'restart_save', self%restart_save)
      call ini_get(file_parameters, ISN, 'residuals_save', self%residuals_save)
      call ini_get(file_parameters, ISN, 'save_residual_fields', self%save_residual_fields, default=.false.)
      ! the [IO] options only PRISM consumes (save_memory_status, save_curl_fields, save_divergence_fields, seam_divB_*)
      ! are read by PRISM (prism_io_object, issue #80 P5)
      ! the keys no parser reads (issue #80 P5): an error that stops the run (default) or warnings; the
      ! environment variable ADAM_UNUSED_KEYS, when set, overrides the key (the regression harnesses set it to error)
      call ini_get(file_parameters, ISN, 'unused_keys', buff_c, default='error')
      call get_environment_variable('ADAM_UNUSED_KEYS', env_value, env_length, env_status)
      if (env_status == 0 .and. env_length > 0) buff_c = env_value
      select case(trim(adjustl(buff_c)))
      case('warn')  ; self%unused_keys_error = .false.
      case('error') ; self%unused_keys_error = .true.
      case default  ; call ini_stop(msg='[IO].(unused_keys) (or ADAM_UNUSED_KEYS) = "'//trim(adjustl(buff_c))// &
                                        '"; accepted: warn, error')
      endselect
      call ini_check()
   endassociate
   call ini_check() ! report every option of this parser that failed to load (#80 P2)
   endsubroutine load_from_file

   ! auxiliary fields IO
   subroutine save_energy_error(self,it,time,blocks_number,energy_D,energy_B,rms_energy_error_D,rms_energy_error_B,&
                                is_to_open,is_to_close)
   !< Save energy error history.
   implicit none
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

   if (mpih%myrank==0) then
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      if (is_to_open_) then
         open(newunit=self%energy_error_unit, file=self%output_basename//'-energy_error.dat')
         write(self%energy_error_unit,'(A)')'VARIABLES="it" "blocks_number" "time" "error_D" "error_B" "rms_error_D" "rms_error_B"'
         flush(self%energy_error_unit)
      endif
      write(self%energy_error_unit, '(A)') trim(str(it                                                  ))//' '//&
                                           trim(str(blocks_number                                       ))//' '//&
                                           trim(str(time                                                ))//' '//&
                                           trim(str(sqrt(abs(energy_D(it)-energy_D(1))/abs(energy_D(1)))))//' '//&
                                           trim(str(sqrt(abs(energy_B(it)-energy_B(1))/abs(energy_B(1)))))//' '//&
                                           trim(str(rms_energy_error_D                                  ))//' '//&
                                           trim(str(rms_energy_error_B                                  ))
      flush(self%energy_error_unit)
      if (is_to_close_) close(self%energy_error_unit)
   endif
   endsubroutine save_energy_error

   subroutine save_energy_history(self,it,time,blocks_number,energy_D,energy_B,coil_power,Poynting_flux, &
                                  is_to_open,is_to_close)
   !< Save energy history.
   implicit none
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

   if (mpih%myrank==0.and.it>=0) then !>= aggiunto da FN, era solo >
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      if (is_to_open_) then
         open(newunit=self%energy_history_unit, file=self%output_basename//'-energy_history.dat')
         write(self%energy_history_unit,'(A)')&
               '%VARIABLES="it" "blocks_number" "time" "D_energy [J]" "B_energy [J]" "coil_power [W]" "Poynting_flux [W]"'
         flush(self%energy_history_unit)
      endif
      if (it > 0) then
         write(self%energy_history_unit, '(A)') trim(str(it               ))//' '//&
                                                trim(str(blocks_number    ))//' '//&
                                                trim(str(time             ))//' '//&
                                                trim(str(energy_D(it)     ))//' '//&
                                                trim(str(energy_B(it)     ))//' '//&
                                                trim(str(coil_power(it)   ))//' '//&
                                                trim(str(Poynting_flux(it)))
         flush(self%energy_history_unit)
      endif
      if (is_to_close_) close(self%energy_history_unit)
   endif
   endsubroutine save_energy_history

   subroutine save_charge_conservation_history(self,time,Q_particles,Q_deposited,Q_deposition_error,R_cont_L1,R_cont_L2, &
                                               R_cont_Linf,dQ_num_net,dQ_num_abs,Gauss_L1,Gauss_L2,Gauss_Linf,             &
                                               is_to_open,is_to_close,do_write)
   !< Save PIC charge-conservation diagnostic history.
   implicit none
   class(io_object), intent(inout)        :: self               !< IO handler.
   real(R8P),        intent(in)           :: time               !< Physical time.
   real(R8P),        intent(in)           :: Q_particles        !< Total particle charge.
   real(R8P),        intent(in)           :: Q_deposited        !< Total grid-deposited charge.
   real(R8P),        intent(in)           :: Q_deposition_error !< Deposited-particle charge mismatch.
   real(R8P),        intent(in)           :: R_cont_L1          !< Continuity residual L1 mean.
   real(R8P),        intent(in)           :: R_cont_L2          !< Continuity residual L2 RMS.
   real(R8P),        intent(in)           :: R_cont_Linf        !< Continuity residual Linf.
   real(R8P),        intent(in)           :: dQ_num_net         !< Net numerical charge equivalent.
   real(R8P),        intent(in)           :: dQ_num_abs         !< Absolute numerical charge equivalent.
   real(R8P),        intent(in)           :: Gauss_L1           !< Gauss-law residual L1 mean.
   real(R8P),        intent(in)           :: Gauss_L2           !< Gauss-law residual L2 RMS.
   real(R8P),        intent(in)           :: Gauss_Linf         !< Gauss-law residual Linf.
   logical,          intent(in), optional :: is_to_open         !< Flag to open file before first saving.
   logical,          intent(in), optional :: is_to_close        !< Flag to close file after last saving.
   logical,          intent(in), optional :: do_write           !< Flag to write a data row.
   logical                                :: is_to_open_        !< Open flag, local var.
   logical                                :: is_to_close_       !< Close flag, local var.
   logical                                :: do_write_          !< Write flag, local var.

   if (mpih%myrank==0) then
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      do_write_    = .true.  ; if (present(do_write   )) do_write_    = do_write
      if (is_to_open_) then
         open(newunit=self%charge_conservation_history_unit, file=self%output_basename//'-charge_conservation_history.dat')
         write(self%charge_conservation_history_unit,'(A)') '# time'
         write(self%charge_conservation_history_unit,'(A)') '# Q_particles'
         write(self%charge_conservation_history_unit,'(A)') '# Q_deposited'
         write(self%charge_conservation_history_unit,'(A)') '# Q_deposition_error'
         write(self%charge_conservation_history_unit,'(A)') '# R_cont_L1'
         write(self%charge_conservation_history_unit,'(A)') '# R_cont_L2'
         write(self%charge_conservation_history_unit,'(A)') '# R_cont_Linf'
         write(self%charge_conservation_history_unit,'(A)') '# dQ_num_net'
         write(self%charge_conservation_history_unit,'(A)') '# dQ_num_abs'
         write(self%charge_conservation_history_unit,'(A)') '# Gauss_L1'
         write(self%charge_conservation_history_unit,'(A)') '# Gauss_L2'
         write(self%charge_conservation_history_unit,'(A)') '# Gauss_Linf'
         flush(self%charge_conservation_history_unit)
      endif
      if (do_write_) then
         write(self%charge_conservation_history_unit, '(A)') trim(str(time              ))//' '//&
                                                             trim(str(Q_particles       ))//' '//&
                                                             trim(str(Q_deposited       ))//' '//&
                                                             trim(str(Q_deposition_error))//' '//&
                                                             trim(str(R_cont_L1         ))//' '//&
                                                             trim(str(R_cont_L2         ))//' '//&
                                                             trim(str(R_cont_Linf       ))//' '//&
                                                             trim(str(dQ_num_net        ))//' '//&
                                                             trim(str(dQ_num_abs        ))//' '//&
                                                             trim(str(Gauss_L1          ))//' '//&
                                                             trim(str(Gauss_L2          ))//' '//&
                                                             trim(str(Gauss_Linf        ))
         flush(self%charge_conservation_history_unit)
      endif
      if (is_to_close_) close(self%charge_conservation_history_unit)
   endif
   endsubroutine save_charge_conservation_history

   subroutine save_divergence_history(self,it,time,blocks_number,div_D,div_B,div_J,div_D_name,div_J_name, &
                                      current_solver_residual,is_to_open,is_to_close)
   !< Save energy history.
   implicit none
   class(io_object), intent(inout)        :: self              !< IO handler.
   integer(I4P),     intent(in)           :: it                !< Current iteration.
   real(R8P),        intent(in)           :: time              !< Current time.
   integer(I4P),     intent(in)           :: blocks_number     !< Current number of blocks.
   real(R8P),        intent(in)           :: div_D             !< Energy history of B.
   real(R8P),        intent(in)           :: div_B             !< Coil power history.
   real(R8P),        intent(in)           :: div_J             !< Poynting flux history.
   character(*),     intent(in), optional :: div_D_name        !< Header label for the first divergence monitor.
   character(*),     intent(in), optional :: div_J_name        !< Header label for the third divergence monitor.
   real(R8P),        intent(in), optional :: current_solver_residual !< Current linear-solve residual diagnostic.
   logical,          intent(in), optional :: is_to_open        !< Flag to open  file before first saving.
   logical,          intent(in), optional :: is_to_close       !< Flag to close file after last saving.
   logical                                :: is_to_open_       !< Flag to open  file before first saving, local var.
   logical                                :: is_to_close_      !< Flag to close file after last saving, local var.
   character(len=:), allocatable          :: div_D_name_       !< Resolved header label for the first divergence monitor.
   character(len=:), allocatable          :: div_J_name_       !< Resolved header label for the third divergence monitor.
   real(R8P)                              :: current_residual_ !< Resolved current linear-solve residual diagnostic.

   if (mpih%myrank==0.and.it>=0) then !>= aggiunto da FN, era solo >
      is_to_open_  = .false. ; if (present(is_to_open )) is_to_open_  = is_to_open
      is_to_close_ = .false. ; if (present(is_to_close)) is_to_close_ = is_to_close
      current_residual_ = 0._R8P ; if (present(current_solver_residual)) current_residual_ = current_solver_residual
      div_D_name_ = 'D_divergence'
      if (present(div_D_name)) div_D_name_ = trim(div_D_name)
      div_J_name_ = 'J_divergence'
      if (present(div_J_name)) div_J_name_ = trim(div_J_name)
      if (is_to_open_) then
         open(newunit=self%divergence_history_unit, file=self%output_basename//'-divergence_history.dat')
         write(self%divergence_history_unit,'(A)')&
               '%VARIABLES="it" "blocks_number" "time" "'//div_D_name_//'" "B_divergence" "'//div_J_name_//'" '//&
               '"current_solver_residual"'
         flush(self%divergence_history_unit)
      endif
      if (it>0) then
         write(self%divergence_history_unit, '(A)') trim(str(it               ))//' '//&
                                                    trim(str(blocks_number    ))//' '//&
                                                    trim(str(time             ))//' '//&
                                                    trim(str(div_D            ))//' '//&
                                                    trim(str(div_B            ))//' '//&
                                                    trim(str(div_J            ))//' '//&
                                                    trim(str(current_residual_))
         flush(self%divergence_history_unit)
      endif
      if (is_to_close_) close(self%divergence_history_unit)
   endif
   endsubroutine save_divergence_history

   ! residuals IO
   subroutine  close_file_residuals(self)
   !< Close file for saving residuals history.
   class(io_object), intent(in) :: self !< IO handler.

   if (mpih%myrank==0) close(self%residuals_unit)
   endsubroutine close_file_residuals

   subroutine open_file_residuals(self, nv, is_restart)
   !< Open file for saving residuals history.
   !<
   !< A restarted run (`is_restart = .true.`) appends to the existing history instead of replacing it, so the history of
   !< a run split by restarts is the history of the uninterrupted run.
   implicit none
   class(io_object), intent(inout)        :: self        !< IO handler.
   integer(I4P),     intent(in)           :: nv          !< Number of residuals variables.
   logical,          intent(in), optional :: is_restart  !< Restarted run: append to the existing file.
   logical                                :: is_restart_ !< Restarted run, local variable.
   character(:), allocatable              :: rqs         !< String buffer.
   integer(I4P)                           :: v           !< Counter.

   is_restart_ = .false. ; if (present(is_restart)) is_restart_ = is_restart
   if (mpih%myrank==0) then
      if (is_restart_) then
         open(newunit=self%residuals_unit, file=self%output_basename//'-residuals.dat', status='unknown', position='append')
      else
         rqs = ''
         do v=1, nv
            rqs = rqs//' "rq'//trim(str(v,.true.))//'"'
         enddo
         open(newunit=self%residuals_unit, file=self%output_basename//'-residuals.dat')
         write(self%residuals_unit, '(A)') 'VARIABLES="it" "time" "blocks_number"'//rqs
         flush(self%residuals_unit)
      endif
   endif
   endsubroutine open_file_residuals

   subroutine save_residuals(self, it, time, blocks_number, residuals)
   !< Save residuals history.
   implicit none
   class(io_object), intent(in) :: self          !< IO handler.
   integer(I4P),     intent(in) :: it            !< Current iteration.
   real(R8P),        intent(in) :: time          !< Current time.
   integer(I4P),     intent(in) :: blocks_number !< Current number of blocks.
   real(R8P),        intent(in) :: residuals(1:) !< Residuals (norm) [1:nv].

   write(self%residuals_unit, '(A)') trim(str(it           ))//' '//&
                                     trim(str(time         ))//' '//&
                                     trim(str(blocks_number))//' '//&
                                     trim(str(residuals(1:), separator=' '))
   flush(self%residuals_unit)
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
