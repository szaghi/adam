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
!<
!< `ini_load(filename, ini)` loads an input file that may extend another (issue #80, P4): a global option (before the
!< first section) `extends = <file>` names a base file, read first (recursively, its path relative to the including
!< file's directory); the including file's options are then added onto it, replacing the values of the options both
!< define and adding the others, sections included. `extends` itself is consumed: it is not an option of the result.
use :: adam_mpih_global, only : mpih
use :: finer, only : file_ini
use :: penf, only : I4P, R8P, str

implicit none
private
public :: ini_check
public :: ini_env_override
public :: ini_error
public :: ini_get
public :: ini_load
public :: ini_mark_read
public :: ini_stop
public :: ini_warn_unused

integer(I4P), parameter :: TEXT_LEN = 999 !< Length of the text buffer used to detect an option.
integer(I4P), parameter :: EXTENDS_DEPTH_MAX = 8 !< Longest chain of files extending one another.

character(len=:), allocatable :: pending          !< Errors recorded since the last ini_check, one per line.
integer(I4P)                  :: pending_number=0 !< Number of errors recorded.

interface ini_get
   !< Read an INI option strictly, with an optional default when absent.
   module procedure ini_get_I4P, ini_get_R8P, ini_get_logical, ini_get_character, ini_get_I4P_array, ini_get_R8P_array
endinterface ini_get

contains
   recursive subroutine ini_load(filename, ini, depth)
   !< Load an input file into `ini`, resolving `extends = <file>` (a global option) recursively: the base file first, then
   !< the including file's options added onto it (an option both define takes the including file's value).
   !<
   !< A missing file loads empty, as FiNeR does (the required options then fail to load); a missing base file, a chain
   !< longer than `EXTENDS_DEPTH_MAX` (a cycle) and an empty `extends` stop the run.
   character(*),           intent(in)    :: filename  !< Input file name.
   type(file_ini),         intent(inout) :: ini       !< Loaded (and merged) INI.
   integer(I4P), optional, intent(in)    :: depth     !< Extension depth of this file (0 for the input itself).
   type(file_ini)                        :: base      !< The file extended, merged.
   type(file_ini)                        :: own       !< This file alone.
   character(len=TEXT_LEN)               :: extends   !< `extends` value.
   character(len=:), allocatable         :: sections(:) !< Sections of this file.
   character(len=:), allocatable         :: pairs(:)  !< Option name/value pair.
   character(len=:), allocatable         :: path      !< Base file path.
   integer(I4P)                          :: depth_    !< Extension depth, local variable.
   integer(I4P)                          :: error     !< Error status.
   integer(I4P)                          :: s         !< Counter.
   integer(I4P)                          :: slash     !< Last directory separator of filename.
   logical                               :: exists    !< Base file found.

   depth_ = 0_I4P ; if (present(depth)) depth_ = depth
   call own%initialize(filename=trim(filename))
   call own%load
   extends = ''
   call own%get(section_name='', option_name='extends', val=extends, error=error)
   if (error /= 0) then ! no extension: the file as it is
      ini = own
      return
   endif
   if (len_trim(extends) == 0) call mpih%error_stop(msg=': '//trim(filename)//': empty "extends"')
   if (depth_ >= EXTENDS_DEPTH_MAX) &
      call mpih%error_stop(msg=': '//trim(filename)//': "extends" chain longer than '//trim(str(EXTENDS_DEPTH_MAX, .true.))// &
                               ' files (a file extending itself?)')
   path = trim(adjustl(extends))
   if (path(1:1) /= '/') then ! relative to the including file's directory
      slash = index(trim(filename), '/', back=.true.)
      if (slash > 0) path = filename(1:slash)//path
   endif
   inquire(file=path, exist=exists)
   if (.not.exists) call mpih%error_stop(msg=': '//trim(filename)//': the file it extends, "'//path//'", does not exist')
   call ini_load(filename=path, ini=base, depth=depth_ + 1_I4P)
   call own%get_sections_list(sections)
   if (allocated(sections)) then
      do s=1, size(sections, dim=1)
         do while (own%loop(section_name=trim(sections(s)), option_pairs=pairs))
            if (len_trim(sections(s)) == 0 .and. trim(pairs(1)) == 'extends') cycle ! consumed, not an option
            call base%add(section_name=trim(sections(s)), option_name=trim(pairs(1)), val=trim(pairs(2)), error=error)
            if (error /= 0) call mpih%error_stop(msg=': '//trim(filename)//': cannot merge ['//trim(sections(s))//'].('// &
                                                     trim(pairs(1))//') onto "'//path//'"')
         enddo
      enddo
   endif
   ini = base
   ini%filename = trim(filename)
   endsubroutine ini_load

   subroutine ini_env_override(name, val)
   !< Override a real option by the environment variable `name`, when it is set (issue #80 P5). The regression harnesses
   !< set `ADAM_MEMORY_FRACTION` so that a tiny case does not allocate and initialize a block pool that fills the node
   !< (the init of a PRISM FNL case took ~74 s of its ~2 min); a value that is not a number stops the run.
   character(*), intent(in)    :: name  !< Environment variable name.
   real(R8P),    intent(inout) :: val   !< Value, replaced when the variable is set.
   character(64)               :: text  !< Variable value.
   integer                     :: length !< Variable length.
   integer                     :: status !< Query status.
   integer                     :: ios    !< Conversion status.
   real(R8P)                   :: parsed !< Converted value.

   call get_environment_variable(name, text, length, status)
   if (status /= 0 .or. length == 0) return
   read(text, *, iostat=ios) parsed
   if (ios /= 0) call ini_stop(msg='environment variable '//name//' = "'//trim(text)//'" is not a number')
   val = parsed
   endsubroutine ini_env_override

   subroutine ini_mark_read(ini, section, option)
   !< Count an option as read (issue #80 P5). FiNeR counts only the options read through `get`, so a parser that walks a
   !< section with `loop` (to refuse the options it does not know) marks each one it accepts, or `ini_warn_unused` would
   !< report it never read.
   type(file_ini), intent(in) :: ini     !< INI file.
   character(*),   intent(in) :: section !< Section name.
   character(*),   intent(in) :: option  !< Option name.

   if (is_present(ini, section, option)) return
   endsubroutine ini_mark_read

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
   subroutine ini_get_I4P_array(ini, section, option, val)
   !< Read an integer list (values separated by blanks) into val, whose size is the number of values expected.
   type(file_ini), intent(in)    :: ini     !< INI file.
   character(*),   intent(in)    :: section !< Section name.
   character(*),   intent(in)    :: option  !< Option name.
   integer(I4P),   intent(inout) :: val(:)  !< Values.
   integer(I4P)                  :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      call record_failed(section, option)
      return
   endif
   call ini%get(section_name=section, option_name=option, val=val, error=error)
   if (error /= 0) call record_failed(section, option)
   endsubroutine ini_get_I4P_array

   subroutine ini_get_R8P_array(ini, section, option, val)
   !< Read a real list (values separated by blanks) into val, whose size is the number of values expected.
   type(file_ini), intent(in)    :: ini     !< INI file.
   character(*),   intent(in)    :: section !< Section name.
   character(*),   intent(in)    :: option  !< Option name.
   real(R8P),      intent(inout) :: val(:)  !< Values.
   integer(I4P)                  :: error   !< Error status.

   if (.not.is_present(ini, section, option)) then
      call record_failed(section, option)
      return
   endif
   call ini%get(section_name=section, option_name=option, val=val, error=error)
   if (error /= 0) call record_failed(section, option)
   endsubroutine ini_get_R8P_array
   subroutine ini_warn_unused(ini, label, as_error)
   !< Report (rank 0) every option of the file that no parser has read since it was loaded (issue #80, P2, F2): a
   !< misspelled key or a key no longer used is otherwise silently ignored. As warnings, the run goes on; `as_error`
   !< (`[IO].(unused_keys) = error`, P5) prints the same list as errors and stops the run.
   type(file_ini), intent(in)           :: ini      !< INI file, after every parser has read it.
   character(*),   intent(in), optional :: label    !< What the file is (e.g. the realm), for the message.
   logical,        intent(in), optional :: as_error !< Stop the run when an option was never read.
   character(len=:), allocatable        :: items(:,:) !< Unused options, section/option pairs.
   character(len=:), allocatable        :: label_     !< Label, local variable.
   character(len=:), allocatable        :: kind_      !< "warning" or "error".
   logical                              :: as_error_  !< Stop the run, local variable.
   integer(I4P)                         :: o          !< Counter.

   if (mpih%myrank /= 0) return
   call ini%unused_options(items=items)
   if (size(items, dim=1) == 0) return
   as_error_ = .false. ; if (present(as_error)) as_error_ = as_error
   kind_ = 'warning' ; if (as_error_) kind_ = 'error'
   label_ = '' ; if (present(label)) label_ = ' of '//label
   if (allocated(ini%filename)) label_ = label_//' ('//trim(ini%filename)//')'
   if (as_error_) then
      print '(A)', mpih%myrankstr//kind_//': '//trim(str(size(items, dim=1), .true.))//' input option(s)'//label_// &
                   ' never read ([IO] unused_keys = error):'
   else
      print '(A)', mpih%myrankstr//kind_//': '//trim(str(size(items, dim=1), .true.))//' input option(s)'//label_// &
                   ' never read, ignored:'
   endif
   do o=1, size(items, dim=1)
      print '(A)', mpih%myrankstr//kind_//':   ['//trim(items(o, 1))//'].('//trim(items(o, 2))//')'
   enddo
   if (as_error_) call mpih%error_stop(msg=': '//trim(str(size(items, dim=1), .true.))//' input option(s)'//label_// &
                                           ' never read (listed above); remove them, or set [IO] unused_keys = warn')
   endsubroutine ini_warn_unused
endmodule adam_ini_library
