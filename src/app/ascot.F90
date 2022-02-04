!< ADAM Slices Converter, ASCot.
program ascot
!< ADAM Slices Converter, ASCot.
!< Convert ADAM slice binary files into ascii format.
use PENF

implicit none

integer(I4P)                 :: file_unit_input     !< File unit input.
character(:),    allocatable :: file_name_input     !< File name input.
integer(I4P)                 :: file_unit_output    !< File unit output.
character(:),    allocatable :: file_name_output    !< File name output.
character(:),    allocatable :: vars_name           !< Variables name.
logical                      :: is_vars_name_passed !< Flag to check if variables name is passed.
character(1000), allocatable :: args(:)             !< Command line args.
integer(I4P)                 :: na                  !< Command line args number.
integer(I4P)                 :: a, i, j, k, v       !< Counter.
integer(I4P)                 :: nijkv(4)            !< IJKV dimensions.
real(R8P),       allocatable :: vars(:)             !< Variables.

na = command_argument_count()
if (na<2) then
   call print_usage
   stop
endif

allocate(args(na))
! get command line args
do a=1, na
   call get_command_argument(a,args(a))
enddo
! parse command line args
is_vars_name_passed = .false.
file_name_output = 'slice.dat'
do a=1, na
   select case(trim(args(a)))
   case('-i')
      file_name_input = trim(args(a+1))
   case('-v')
      vars_name = trim(args(a+1))
      is_vars_name_passed = .true.
   case('-o')
      file_name_output = trim(args(a+1))
   endselect
enddo

! open files
open(newunit=file_unit_input,  file=file_name_input,  form='unformatted', access='stream')
open(newunit=file_unit_output, file=file_name_output)

! load header
read(unit=file_unit_input) nijkv
allocate(vars(nijkv(4)))
print '(A)', 'slice dimensions [ni, nj, nk, nv]: '//trim(str(nijkv))
if (.not.is_vars_name_passed) then
   vars_name = 'VARIABLES = "x" "y" "Z"'
   do v=4, nijkv(4)
      vars_name = vars_name//' "v'//trim(str(v,.true.))//'"'
   enddo
endif

! save header
write(unit=file_unit_output,'(A)') vars_name
write(unit=file_unit_output,'(A)') 'ZONE'
write(unit=file_unit_output,'(A)') 'I = '//trim(str(nijkv(1)))//', '//'J = '//trim(str(nijkv(2)))//', '//'K = '//trim(str(nijkv(3)))

! load/save data
do k=1, nijkv(3)
   do j=1, nijkv(2)
      do i=1, nijkv(1)
         read(unit=file_unit_input) vars
         write(unit=file_unit_output,'(A)') trim(str(vars, separator=' '))
      enddo
   enddo
enddo

! close files
close(unit=file_unit_input)
close(unit=file_unit_output)

contains
   subroutine print_usage
   !< Print app usage.

   print '(A)', 'Examples of ASCot usage:'
   print '(A)', '  1) ascot -i slice.mat -v "rho rhou rhov rhow rhoe" -o slice.dat'
   print '(A)', '  2) ascot -i slice.mat -o slice.dat'
   print '(A)', '  3) ascot -i slice.mat'
   endsubroutine print_usage
endprogram ascot
