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
!<
!< It also checks `ini_get_list` (a list takes as many values as are written; absent, its default, possibly empty) and
!< `ini_count_sections` (the enumerated sections `<prefix>N`), F4 and F5 of issue #80. With `gap`: counts sections
!< numbered with a gap, which must stop the run.
use :: adam_ini_library, only : ini_count_sections, ini_get_list, ini_load
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P

implicit none

type(file_ini) :: ini      !< Loaded INI.
character(99)  :: arg      !< Command line argument.
integer(I4P)   :: failed=0 !< Failed checks.
real(R8P),    allocatable :: rlist(:) !< Real list.
integer(I4P), allocatable :: ilist(:) !< Integer list.

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
case('gap')
   call write_file('gap.ini', '[slice_1]'//new_line('a')//'k = 1'//new_line('a')//'[slice_3]'//new_line('a')//'k = 3')
   call ini_load(filename='gap.ini', ini=ini)
   failed = ini_count_sections(ini, 'slice_') ! not inside the print: the function prints its error (recursive I/O)
   print '(A,I0)', 'FAIL: sections numbered with a gap counted ', failed
   stop 1
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
call ini%free
call write_file('lists.ini', '[cooling]'//new_line('a')//'temperatures = 10 100 1000'//new_line('a')// &
                             'ids = 4 5'//new_line('a')//'[slice_1]'//new_line('a')//'k = 1'//new_line('a')// &
                             '[slice_2]'//new_line('a')//'k = 2'//new_line('a')//'[slice_x]'//new_line('a')//'k = 0'// &
                             new_line('a')//'[slices]'//new_line('a')//'k = 0')
call ini_load(filename='lists.ini', ini=ini)
call ini_get_list(ini, 'cooling', 'temperatures', rlist)
call verdict('lists: a real list takes the 3 values written', size(rlist) == 3)
if (size(rlist) == 3) call verdict('lists: the real values', all(rlist == [10._R8P, 100._R8P, 1000._R8P]))
call ini_get_list(ini, 'cooling', 'ids', ilist)
call verdict('lists: an integer list takes the 2 values written', size(ilist) == 2)
call ini_get_list(ini, 'cooling', 'lambdas', rlist, default=[real(R8P)::])
call verdict('lists: an absent list with an empty default is empty', size(rlist) == 0)
call ini_get_list(ini, 'cooling', 'lambdas', rlist, default=[1._R8P, 2._R8P])
call verdict('lists: an absent list takes its default', size(rlist) == 2)
call verdict('sections: [slice_1], [slice_2] counted (not [slice_x], [slices])', ini_count_sections(ini, 'slice_') == 2)
call verdict('sections: no [solid_N] counts 0', ini_count_sections(ini, 'solid_') == 0)
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
