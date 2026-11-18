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
!< + absent, without `default` (a required option): the error `failed to load [<section>].(<option>)` is recorded.
!<
!< Recorded errors do not stop the run at once (issue #80, P2, F1): a parser reads all its options, then calls `ini_check`,
!< which stops with every error recorded so far, one line each, so that one run reports every missing or malformed option
!< instead of one per run. A single error stops with the message every parser used before,
!< `: failed to load [<section>].(<option>)`. `ini_stop(msg)` records a validation error (a value outside its accepted set)
!< and stops at once with the whole list; `ini_error(msg)` records one and lets the parser go on. A parser must call
!< `ini_check` before it uses a value it has just read to size, loop or branch, when a failed read would leave that value
!< meaningless (a count key without default, for instance).
!<
!< The recorded errors are module state of this library, kept only between a parser's reads and its `ini_check`.
!< Absence is told apart from a bad value by reading the option first as text, which fails only when the option is not in
!< the section (FiNeR's `has_option` searches every section, not one).
use :: adam_mpih_global, only : mpih
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: ini_check
public :: ini_error
public :: ini_get
public :: ini_stop

integer(I4P), parameter :: TEXT_LEN = 999 !< Length of the text buffer used to detect an option.

character(len=:), allocatable :: pending          !< Errors recorded since the last ini_check, one per line.
integer(I4P)                  :: pending_number=0 !< Number of errors recorded.

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

   subroutine record_failed(section, option)
   !< Record that the option is missing (and required) or not convertible.
   character(*), intent(in) :: section !< Section name.
   character(*), intent(in) :: option  !< Option name.

   call ini_error(msg='failed to load ['//section//'].('//option//')')
   endsubroutine record_failed

   subroutine ini_error(msg)
   !< Record an input error; the run goes on until the next `ini_check` (or `ini_stop`).
   character(*), intent(in) :: msg !< Error message, without the leading ': '.

   if (pending_number == 0) then
      pending = msg
   else
      pending = pending//new_line('a')//msg
   endif
   pending_number = pending_number + 1
   endsubroutine ini_error

   subroutine ini_check()
   !< Stop the run if any input error has been recorded, reporting all of them.
   !<
   !< One error: `: <message>`, as every parser reported it before. More: a count, then one indented message per line.
   character(len=:), allocatable :: msg !< Report.

   if (pending_number == 0) return
   if (pending_number == 1) then
      msg = ': '//pending
   else
      msg = ': '//trim(str(pending_number, .true.))//' input errors:'//new_line('a')//'  '// &
            replace_all(pending, new_line('a'), new_line('a')//'  ')
   endif
   call mpih%error_stop(msg=msg)
   endsubroutine ini_check

   subroutine ini_stop(msg)
   !< Record an input error and stop at once with every error recorded (a value the parser cannot go on with).
   character(*), intent(in) :: msg !< Error message, without the leading ': '.

   call ini_error(msg=msg)
   call ini_check()
   endsubroutine ini_stop

   pure function replace_all(string, old, new) result(replaced)
   !< Return string with every occurrence of old replaced by new.
   character(*), intent(in)      :: string   !< Input.
   character(*), intent(in)      :: old      !< Substring to replace.
   character(*), intent(in)      :: new      !< Replacement.
   character(len=:), allocatable :: replaced !< Output.
   integer(I4P)                  :: i, j     !< Positions.

   replaced = ''
   i = 1
   do
      j = index(string(i:), old)
      if (j == 0) exit
      replaced = replaced//string(i:i+j-2)//new
      i = i + j - 1 + len(old)
   enddo
   replaced = replaced//string(i:)
   endfunction replace_all

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
      if (present(default)) then
         val = default
      else
         call record_failed(section, option)
      endif
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) then
      call record_failed(section, option) ! val left unchanged
      return
   endif
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
      if (present(default)) then
         val = default
      else
         call record_failed(section, option)
      endif
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) then
      call record_failed(section, option) ! val left unchanged
      return
   endif
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
      if (present(default)) then
         val = default
      else
         call record_failed(section, option)
      endif
      return
   endif
   call ini%get(section_name=section, option_name=option, val=buff, error=error)
   if (error /= 0) then
      call record_failed(section, option) ! val left unchanged
      return
   endif
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
      if (present(default)) then
         val = default
      else
         call record_failed(section, option)
      endif
      return
   endif
   call ini%get(section_name=section, option_name=option, val=val, error=error)
   if (error /= 0) call record_failed(section, option)
   endsubroutine ini_get_character
endmodule adam_ini_library
