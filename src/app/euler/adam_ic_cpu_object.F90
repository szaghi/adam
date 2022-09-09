!< ADAM, Initial Conditions class definition, CPU backend.
module adam_ic_cpu_object
!< ADAM, Initial Conditions class definition, CPU backend.

use adam_grid_object, only : grid_object
use adam_mpih_object, only : mpih_object
use finer
use penf

implicit none
private
public :: ic_cpu_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions" !< INI (config) file section name containing IC configs.

type :: ic_cpu_object
   !< Initial Conditions class definition, CPU backend.
   type(mpih_object)      :: mpih                 !< MPI handler.
   integer(I4P)           :: regions_number=1     !< Number of IC regions.
   real(R8P), allocatable :: q(:,:)               !< Primitive variables (r,u,v,w,p,s with s being the speccie index) at IC.
   real(R8P), allocatable :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize IC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype ic_cpu_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(ic_cpu_object), intent(in) :: self             !< IC.
   character(len=:), allocatable    :: desc             !< Description.
   character(len=1), parameter      :: NL=new_line('a') !< New line character.
   integer(I4P)                     :: r, v             !< Counter.

   desc =       self%mpih%myrankstr//'IC main data'//NL
   desc = desc//self%mpih%myrankstr//'  regions number: '//trim(str(self%regions_number))
   do r=1, self%regions_number
      desc = desc//NL//self%mpih%myrankstr//'  region('//trim(str(r,.true.))//')'
      do v=1, size(self%q, dim=1)
         desc = desc//NL//self%mpih%myrankstr//'    q('//trim(str(v,.true.))//'): '//trim(str(self%q(v,r)))
      enddo
         desc = desc//NL//self%mpih%myrankstr//'    emin: '//trim(str(self%emin(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    emax: '//trim(str(self%emax(:,r)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(ic_cpu_object), intent(inout) :: self            !< IC.
   type(file_ini),       intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'ic_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'ic_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(ic_cpu_object), intent(inout)        :: self            !< IC.
   type(file_ini),       intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,              intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                    :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                  :: sname           !< Section name.
   integer(I4P)                               :: i               !< Counter.
   integer(I4P)                               :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')

   if (self%regions_number>=1) then
      allocate(   self%q(1:6, 1:self%regions_number))
      allocate(self%emin(1:3, 1:self%regions_number))
      allocate(self%emax(1:3, 1:self%regions_number))
      do i=1, self%regions_number
         sname = INI_SECTION_NAME//'_region_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='r', val=self%q(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r)')
         call file_parameters%get(section_name=sname, option_name='u', val=self%q(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(u)')
         call file_parameters%get(section_name=sname, option_name='v', val=self%q(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(v)')
         call file_parameters%get(section_name=sname, option_name='w', val=self%q(4,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(w)')
         call file_parameters%get(section_name=sname, option_name='p', val=self%q(5,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(p)')
         call file_parameters%get(section_name=sname, option_name='s', val=self%q(6,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(s)')
         call file_parameters%get(section_name=sname, option_name='emin_x', val=self%emin(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_x)')
         call file_parameters%get(section_name=sname, option_name='emin_y', val=self%emin(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_y)')
         call file_parameters%get(section_name=sname, option_name='emin_z', val=self%emin(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emin_z)')
         call file_parameters%get(section_name=sname, option_name='emax_x', val=self%emax(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_x)')
         call file_parameters%get(section_name=sname, option_name='emax_y', val=self%emax(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_y)')
         call file_parameters%get(section_name=sname, option_name='emax_z', val=self%emax(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(emax_z)')
      enddo
   endif
   endsubroutine load_from_file
endmodule adam_ic_cpu_object
