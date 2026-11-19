!< ADAM, test of adam_ini_library%ini_load: an input file extending another (issue #80, P4).
program test_ini_extends
!< ADAM, test of adam_ini_library%ini_load: an input file extending another (issue #80, P4).
!<
!< Without argument: writes a chain top.ini -> mid.ini -> base.ini (and sub/c.ini -> ../base.ini) in the current
!< directory and checks the merged values: the including file's value wins, the base's other options and sections stay,
!< new options and sections are added, `extends` is consumed, a relative path is relative to the including file.
!< With `missing` or `cycle`: loads a file extending a missing file, or a file of a cycle; both must stop the run (the
!< checker, check.sh, expects a non-zero exit and the message). With `dump <file>`: loads the file (extends resolved) and
!< prints every option as `[section].(option) = value`, sorted, to compare two inputs as the applications read them.
use :: adam_ini_library, only : ini_load
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P

implicit none

type(file_ini) :: ini      !< Loaded INI.
character(99)  :: arg      !< Command line argument.
integer(I4P)   :: failed=0 !< Failed checks.

arg = ''
if (command_argument_count() > 0) call get_command_argument(1, arg)
select case(trim(arg))
case('missing')
   call write_file('bad.ini', 'extends = nothere.ini'//new_line('a')//'[grid]'//new_line('a')//'ni = 4')
   call ini_load(filename='bad.ini', ini=ini)
   print '(A)', 'FAIL: a missing base file did not stop the run'
   stop 1
case('dump')
   call get_command_argument(2, arg)
   call ini_load(filename=trim(arg), ini=ini)
   call dump
   stop
case('cycle')
   call write_file('ca.ini', 'extends = cb.ini'//new_line('a')//'[grid]'//new_line('a')//'ni = 4')
   call write_file('cb.ini', 'extends = ca.ini'//new_line('a')//'[grid]'//new_line('a')//'nj = 4')
   call ini_load(filename='ca.ini', ini=ini)
   print '(A)', 'FAIL: a cycle did not stop the run'
   stop 1
endselect
call write_file('base.ini', '[grid]'//new_line('a')//'ni = 16'//new_line('a')//'nj = 8'//new_line('a')// &
                            'emin = 0 0 0'//new_line('a')//'[time]'//new_line('a')//'CFL = 0.5')
call write_file('mid.ini', 'extends = base.ini'//new_line('a')//'[grid]'//new_line('a')//'ni = 32'//new_line('a')// &
                           '[new]'//new_line('a')//'k = 1')
call write_file('top.ini', 'extends = mid.ini'//new_line('a')//'[time]'//new_line('a')//'CFL = 0.25')
call execute_command_line('mkdir -p sub')
call write_file('sub/c.ini', 'extends = ../base.ini'//new_line('a')//'[grid]'//new_line('a')//'nj = 2')
call ini_load(filename='top.ini', ini=ini)
call check_i('top: [grid] ni from mid (overrides base)', 'grid', 'ni', 32)
call check_i('top: [grid] nj from base', 'grid', 'nj', 8)
call check_r('top: [time] CFL from top (overrides base)', 'time', 'CFL', 0.25_R8P)
call check_i('top: [new] k, a section mid adds', 'new', 'k', 1)
call check_a('top: [grid] emin, an array value, from base', 'grid', 'emin', [0._R8P, 0._R8P, 0._R8P])
call check_absent('top: extends is consumed', 'extends')
call check_s('top: the merged INI keeps the input file name', ini%filename, 'top.ini')
call ini%free
call ini_load(filename='sub/c.ini', ini=ini)
call check_i('sub/c: ../base.ini relative to sub/', 'grid', 'ni', 16)
call check_i('sub/c: [grid] nj from c', 'grid', 'nj', 2)
call ini%free
call ini_load(filename='base.ini', ini=ini)
call check_i('base: a file without extends loads as it is', 'grid', 'ni', 16)
if (failed > 0) then
   print '(A,I0,A)', 'test_ini_extends: ', failed, ' check(s) FAILED'
   stop 1
endif
print '(A)', 'test_ini_extends: all checks PASSED'
contains
   subroutine dump()
   !< Print every option as `[section].(option) = value`, one per line, sorted by section then option.
   character(len=:), allocatable :: sections(:) !< Sections.
   character(len=:), allocatable :: pairs(:)    !< Option name/value pair.
   character(len=999), allocatable :: lines(:)  !< Output lines.
   character(len=999)              :: swap      !< Swap buffer.
   integer(I4P)                    :: s, i, j, n !< Counters.

   allocate(lines(0)) ; n = 0
   call ini%get_sections_list(sections)
   if (allocated(sections)) then
      do s=1, size(sections, dim=1)
         do while (ini%loop(section_name=trim(sections(s)), option_pairs=pairs))
            lines = [lines, '['//trim(sections(s))//'].('//trim(pairs(1))//') = '//trim(pairs(2))] ; n = n + 1
         enddo
      enddo
   endif
   do i=2, n ! insertion sort, inputs are small
      swap = lines(i) ; j = i - 1
      do while (j >= 1)
         if (lines(j) <= swap) exit
         lines(j+1) = lines(j) ; j = j - 1
      enddo
      lines(j+1) = swap
   enddo
   do i=1, n
      print '(A)', trim(lines(i))
   enddo
   endsubroutine dump

   subroutine write_file(name, text)
   !< Write a text file.
   character(*), intent(in) :: name !< File name.
   character(*), intent(in) :: text !< Contents.
   integer                  :: u    !< Unit.

   open(newunit=u, file=name, status='replace', action='write')
   write(u, '(A)') text
   close(u)
   endsubroutine write_file

   subroutine verdict(label, ok)
   !< Print one check.
   character(*), intent(in) :: label !< Check.
   logical,      intent(in) :: ok    !< Result.

   if (ok) then
      print '(A)', 'PASS '//label
   else
      print '(A)', 'FAIL '//label ; failed = failed + 1
   endif
   endsubroutine verdict

   subroutine check_i(label, section, option, expected)
   !< Check an integer option.
   character(*), intent(in) :: label, section, option !< Check, section, option.
   integer(I4P), intent(in) :: expected               !< Expected value.
   integer(I4P)             :: val, error             !< Value, error.

   val = -huge(1_I4P) ; call ini%get(section_name=section, option_name=option, val=val, error=error)
   call verdict(label, error == 0 .and. val == expected)
   endsubroutine check_i

   subroutine check_r(label, section, option, expected)
   !< Check a real option.
   character(*), intent(in) :: label, section, option !< Check, section, option.
   real(R8P),    intent(in) :: expected               !< Expected value.
   real(R8P)                :: val                    !< Value.
   integer(I4P)             :: error                  !< Error.

   val = -huge(1._R8P) ; call ini%get(section_name=section, option_name=option, val=val, error=error)
   call verdict(label, error == 0 .and. val == expected)
   endsubroutine check_r

   subroutine check_a(label, section, option, expected)
   !< Check a real array option.
   character(*), intent(in) :: label, section, option !< Check, section, option.
   real(R8P),    intent(in) :: expected(:)            !< Expected values.
   real(R8P)                :: val(size(expected))    !< Values.
   integer(I4P)             :: error                  !< Error.

   val = -huge(1._R8P) ; call ini%get(section_name=section, option_name=option, val=val, error=error)
   call verdict(label, error == 0 .and. all(val == expected))
   endsubroutine check_a

   subroutine check_absent(label, option)
   !< Check that no section holds an option.
   character(*), intent(in) :: label, option !< Check, option.

   call verdict(label, .not.ini%has_option(option_name=option))
   endsubroutine check_absent

   subroutine check_s(label, got, expected)
   !< Check a string.
   character(*), intent(in) :: label, got, expected !< Check, value, expected value.

   call verdict(label, trim(got) == trim(expected))
   endsubroutine check_s
endprogram test_ini_extends
