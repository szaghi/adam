!< ADAM, strict reads of INI options with optional defaults (issue #80, P1).
module adam_ini_library
!< ADAM, strict reads of INI options with optional defaults (issue #80, P1).
!<
!< `ini_get(ini, section, option, val[, default])` reads one option of one section:
!<
!< + present and convertible to the type of `val`: `val` takes it;
!< + present but not convertible (e.g. `frequency = abc`): the run stops, as before (FiNeR's own `default=` would silently
!<   take the default here);
!< + absent, with `default`: `val` takes the default;
!< + absent, without `default` (a required option): the run stops with `: failed to load [<section>].(<option>)`, the
!<   message every parser used before.
!<
!< Absence is told apart from a bad value by reading the option first as text, which fails only when the option is not in
!< the section (FiNeR's `has_option` searches every section, not one).
use :: adam_mpih_global, only : mpih
use :: finer, only : file_ini
use :: penf, only : I4P, R8P

implicit none
private
public :: ini_get

integer(I4P), parameter :: TEXT_LEN = 999 !< Length of the text buffer used to detect an option.

interface ini_get
   !< Read an INI option strictly, with an optional default when absent.
   module procedure ini_get_I4P, ini_get_R8P, ini_get_logical, ini_get_character
endinterface ini_get

contains
   logical function is_present(ini, section, option)
   !< True if the option is in the section.
   type(file_ini), intent(in) :: ini     !< INI file.
   character(*),   intent(in) :: section !< Section name.
   character(*),   intent(in) :: option  !< Option name.
   character(len=TEXT_LEN)    :: text    !< Option value as text.
   integer(I4P)               :: error   !< Error status.

   call ini%get(section_name=section, option_name=option, val=text, error=error)
   is_present = error == 0
   endfunction is_present

   subroutine stop_failed(section, option)
   !< Stop the run: the option is missing (and required) or not convertible.
   character(*), intent(in) :: section !< Section name.
   character(*), intent(in) :: option  !< Option name.

   call mpih%error_stop(msg=': failed to load ['//section//'].('//option//')')
   endsubroutine stop_failed

   subroutine ini_get_I4P(ini, section, option, val, default)
   !< Read an integer option.
   type(file_ini), intent(in)           :: ini     !< INI file.
   character(*),   intent(in)           :: section !< Section name.
   character(*),   intent(in)           :: option  !< Option name.
   integer(I4P),   intent(inout)        :: val     !< Value.
   integer(I4P),   intent(in), optional :: default !< Default, when the option is absent.
   integer(I4P)                         :: buff    !< Buffer: val is left untouched on a failed read.
   integer(I4P)                         :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      if (.not.present(default)) call stop_failed(section, option)
      val = default
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) call stop_failed(section, option)
   val = buff
   endsubroutine ini_get_I4P

   subroutine ini_get_R8P(ini, section, option, val, default)
   !< Read a real option.
   type(file_ini), intent(in)           :: ini     !< INI file.
   character(*),   intent(in)           :: section !< Section name.
   character(*),   intent(in)           :: option  !< Option name.
   real(R8P),      intent(inout)        :: val     !< Value.
   real(R8P),      intent(in), optional :: default !< Default, when the option is absent.
   real(R8P)                            :: buff    !< Buffer: val is left untouched on a failed read.
   integer(I4P)                         :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      if (.not.present(default)) call stop_failed(section, option)
      val = default
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) call stop_failed(section, option)
   val = buff
   endsubroutine ini_get_R8P

   subroutine ini_get_logical(ini, section, option, val, default)
   !< Read a logical option.
   type(file_ini), intent(in)           :: ini     !< INI file.
   character(*),   intent(in)           :: section !< Section name.
   character(*),   intent(in)           :: option  !< Option name.
   logical,        intent(inout)        :: val     !< Value.
   logical,        intent(in), optional :: default !< Default, when the option is absent.
   logical                              :: buff    !< Buffer: val is left untouched on a failed read.
   integer(I4P)                         :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      if (.not.present(default)) call stop_failed(section, option)
      val = default
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) call stop_failed(section, option)
   val = buff
   endsubroutine ini_get_logical

   subroutine ini_get_character(ini, section, option, val, default)
   !< Read a text option into a fixed-length buffer (the caller trims it).
   type(file_ini), intent(in)           :: ini     !< INI file.
   character(*),   intent(in)           :: section !< Section name.
   character(*),   intent(in)           :: option  !< Option name.
   character(*),   intent(inout)        :: val     !< Value.
   character(*),   intent(in), optional :: default !< Default, when the option is absent.
   integer(I4P)                         :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      if (.not.present(default)) call stop_failed(section, option)
      val = default
      return
   endif
   call ini%get(section_name=section, option_name=option, val=val, error=error)
   if (error /= 0) call stop_failed(section, option)
   endsubroutine ini_get_character
endmodule adam_ini_library
