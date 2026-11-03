!< ADAM, FLUME mass-flow forcing of a periodic channel (issue #119, T3.1).

module adam_flume_forcing_object
!< ADAM, FLUME mass-flow forcing of a periodic channel (issue #119, T3.1).
!<
!< `[forcing] type = mass-flow`, `axis` (x, y, z), `bulk_momentum` (the target volume mean of the momentum along the axis). Once per step, on the committed
!< state, a uniform momentum `dm = bulk_momentum - mean` is added, so the bulk momentum is held exactly; the total energy gains the change of kinetic energy, the
!< internal energy is untouched. `f = dm / dt` is the driving force per unit volume, written on the conservation history cadence to
!< `<output_basename>-forcing_history.dat` (header `VARIABLES="it" "time" "f"`, the output units of `[reference]`): at a statistically steady state its mean
!< balances the wall friction. Without the section nothing changes, bitwise.

use :: adam_flume_parameters, only : MODEL_EULER
use :: adam_mpih_global,      only : mpih
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P

implicit none
private
public :: flume_forcing_object

character(len=7), parameter :: INI_SECTION_NAME='forcing' !< Section name.

type :: flume_forcing_object
   !< Mass-flow forcing.
   logical      :: is_active=.false.    !< Forcing on.
   integer(I4P) :: axis=1_I4P           !< Axis of the flow.
   real(R8P)    :: bulk_momentum=0._R8P !< Target volume mean of the momentum along the axis.
   real(R8P)    :: f=0._R8P             !< Last driving force per unit volume.
   integer(I4P) :: unit=0_I4P           !< History file unit.
   contains
      procedure, pass(self) :: load_from_file !< Load `[forcing]`.
      procedure, pass(self) :: save_row       !< Write one history row.
endtype flume_forcing_object

contains
   subroutine load_from_file(self, file_parameters, model)
   !< Load `[forcing]`; an unknown value is fatal.
   class(flume_forcing_object), intent(inout) :: self            !< Forcing.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters file handler.
   integer(I4P),                intent(in)    :: model           !< Physical model.
   character(999)                             :: buff            !< Option value buffer.
   integer(I4P)                               :: error           !< Error status.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='type', val=buff, error=error)
   if (error /= 0) return
   select case(trim(adjustl(buff)))
   case('none')
      return
   case('mass-flow')
      self%is_active = .true.
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(type); expected none or mass-flow')
   endselect
   if (model /= MODEL_EULER) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] requires [physics].(physical_model) = euler')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='axis', val=buff, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(axis)')
   select case(trim(adjustl(buff)))
   case('x')
      self%axis = 1_I4P
   case('y')
      self%axis = 2_I4P
   case('z')
      self%axis = 3_I4P
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(axis); expected x, y or z')
   endselect
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='bulk_momentum', val=self%bulk_momentum, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(bulk_momentum)')
   endsubroutine load_from_file

   subroutine save_row(self, output_basename, it, time, force, is_restart)
   !< Write one history row, rank 0: a new file with its header, or appended on a restart, so that the restart continues it. `time` and
   !< `force` are the output factors of the time and of the force per unit volume (1 in code units, issue #119).
   class(flume_forcing_object), intent(inout) :: self            !< Forcing.
   character(*),                intent(in)    :: output_basename !< Output base name.
   integer(I4P),                intent(in)    :: it              !< Time iteration.
   real(R8P),                   intent(in)    :: time            !< Time (output units).
   real(R8P),                   intent(in)    :: force           !< Output factor of the force per unit volume.
   logical,                     intent(in)    :: is_restart      !< Restarted run.

   if (mpih%myrank /= 0) return
   if (self%unit == 0_I4P) then
      if (is_restart) then
         open(newunit=self%unit, file=trim(output_basename)//'-forcing_history.dat', status='unknown', position='append', action='write')
      else
         open(newunit=self%unit, file=trim(output_basename)//'-forcing_history.dat', status='replace', action='write')
         write(self%unit, '(A)') 'VARIABLES="it" "time" "f"'
      endif
   endif
   write(self%unit, '(I10,*(1X,ES24.16E3))') it, time, self%f * force
   flush(self%unit)
   endsubroutine save_row
endmodule adam_flume_forcing_object
