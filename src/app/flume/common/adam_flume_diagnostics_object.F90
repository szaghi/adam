!< ADAM, FLUME diagnostics class definition: conservation and div(B) histories.
module adam_flume_diagnostics_object
!< ADAM, FLUME diagnostics class definition: conservation and div(B) histories.
!<
!< Writes one row per save step to `<output_basename>-conservation_history.dat`: iteration, time and the volume
!< integral of every conservative variable. MHD runs also write, on the same cadence,
!< `<output_basename>-divb_history.dat`: iteration, time, `max|div B|`, `sum |div B| dV` and the maximum over the
!< seam-local cells (issue #41, D-10, section 9). Rows are written as they come, nothing is accumulated in memory. The
!< values are computed (and MPI-reduced) by the backends; only rank 0 writes.

! ADAM singleton objects
use :: adam_mpih_global, only : mpih
! third party modules
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P, str
use :: stringifor,       only : string

implicit none
private
public :: flume_diagnostics_object

character(len=11), parameter :: INI_SECTION_NAME="diagnostics" !< INI (config) file section name.

type :: flume_diagnostics_object
   !< FLUME diagnostics class definition: conservation history.
   integer(I4P) :: conservation_history_save=0_I4P !< Conservation history save cadence (<= 0: disabled).
   integer(I4P) :: conservation_unit=0_I4P         !< Conservation history file unit.
   integer(I4P) :: divb_unit=0_I4P                 !< div(B) history file unit (MHD only).
   contains
      ! public methods
      procedure, pass(self) :: close_file            !< Close the history file.
      procedure, pass(self) :: description           !< Return pretty-printed object description.
      procedure, pass(self) :: initialize            !< Initialize diagnostics.
      procedure, pass(self) :: load_from_file        !< Load config from file.
      procedure, pass(self) :: open_file             !< Open the history file and write its header.
      procedure, pass(self) :: save_conservation_row !< Write one conservation history row.
      procedure, pass(self) :: save_divb_row         !< Write one div(B) history row.
endtype flume_diagnostics_object

contains
   ! public methods
   subroutine close_file(self)
   !< Close the history files.
   class(flume_diagnostics_object), intent(inout) :: self    !< Diagnostics.
   logical                                        :: is_open !< Unit status.

   if (mpih%myrank /= 0) return
   if (self%conservation_unit /= 0_I4P) then
      inquire(unit=self%conservation_unit, opened=is_open)
      if (is_open) close(self%conservation_unit)
      self%conservation_unit = 0_I4P
   endif
   if (self%divb_unit /= 0_I4P) then
      inquire(unit=self%divb_unit, opened=is_open)
      if (is_open) close(self%divb_unit)
      self%divb_unit = 0_I4P
   endif
   endsubroutine close_file

   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_diagnostics_object), intent(in) :: self !< Diagnostics.
   character(len=:), allocatable               :: desc !< Description.

   desc = mpih%myrankstr//'Diagnostics main data'//new_line('a')// &
          mpih%myrankstr//'  conservation_history_save: '//trim(str(self%conservation_history_save))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize diagnostics.
   class(flume_diagnostics_object), intent(inout) :: self            !< Diagnostics.
   type(file_ini),                  intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_diagnostics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_diagnostics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key is required.
   class(flume_diagnostics_object), intent(inout) :: self            !< Diagnostics.
   type(file_ini),                  intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P)                                   :: error           !< Error status.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='conservation_history_save', &
                            val=self%conservation_history_save, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(conservation_history_save)')
   endsubroutine load_from_file

   subroutine open_file(self, output_basename, q_name, is_restart, with_divb)
   !< Open the history files and write their headers (append, without header, on restart).
   class(flume_diagnostics_object), intent(inout) :: self            !< Diagnostics.
   character(*),                    intent(in)    :: output_basename !< Output files base name.
   type(string),                    intent(in)    :: q_name(1:)      !< Conservative variables names.
   logical,                         intent(in)    :: is_restart      !< Restarted run flag.
   logical,                         intent(in)    :: with_divb       !< Open the div(B) history too (MHD).
   character(:), allocatable                      :: header          !< File header.
   integer(I4P)                                   :: v               !< Counter.

   if (mpih%myrank /= 0 .or. self%conservation_history_save <= 0_I4P) return
   if (with_divb) then
      if (is_restart) then
         open(newunit=self%divb_unit, file=trim(output_basename)//'-divb_history.dat', &
              status='unknown', position='append', action='write')
      else
         open(newunit=self%divb_unit, file=trim(output_basename)//'-divb_history.dat', status='replace', action='write')
         write(self%divb_unit, '(A)') 'VARIABLES="it" "time" "max_divb" "l1_divb" "seam_max_divb"'
      endif
   endif
   if (is_restart) then
      open(newunit=self%conservation_unit, file=trim(output_basename)//'-conservation_history.dat', &
           status='unknown', position='append', action='write')
   else
      open(newunit=self%conservation_unit, file=trim(output_basename)//'-conservation_history.dat', &
           status='replace', action='write')
      header = 'VARIABLES="it" "time"'
      do v=1, size(q_name, dim=1)
         header = header//' "int_'//q_name(v)%chars()//'"'
      enddo
      write(self%conservation_unit, '(A)') header
   endif
   endsubroutine open_file

   subroutine save_conservation_row(self, it, time, integrals)
   !< Write one conservation history row (rank 0 only; integrals already MPI-reduced).
   class(flume_diagnostics_object), intent(in) :: self          !< Diagnostics.
   integer(I4P),                    intent(in) :: it            !< Time iteration.
   real(R8P),                       intent(in) :: time          !< Time.
   real(R8P),                       intent(in) :: integrals(1:) !< Volume integrals of the conservative variables.

   if (mpih%myrank /= 0 .or. self%conservation_unit == 0_I4P) return
   write(self%conservation_unit, '(I10,*(1X,ES24.16E3))') it, time, integrals
   flush(self%conservation_unit)
   endsubroutine save_conservation_row

   subroutine save_divb_row(self, it, time, norms)
   !< Write one div(B) history row (rank 0 only; norms already MPI-reduced).
   class(flume_diagnostics_object), intent(in) :: self     !< Diagnostics.
   integer(I4P),                    intent(in) :: it       !< Time iteration.
   real(R8P),                       intent(in) :: time     !< Time.
   real(R8P),                       intent(in) :: norms(3) !< max|div B|, sum |div B| dV, seam-local max|div B|.

   if (mpih%myrank /= 0 .or. self%divb_unit == 0_I4P) return
   write(self%divb_unit, '(I10,*(1X,ES24.16E3))') it, time, norms
   flush(self%divb_unit)
   endsubroutine save_divb_row
endmodule adam_flume_diagnostics_object
