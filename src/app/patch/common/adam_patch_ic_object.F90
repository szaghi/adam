!< ADAM, PATCH Initial Conditions class definition, common CPU backend.
module adam_patch_ic_object
!< ADAM, PATCH Initial Conditions class definition, common CPU backend.

! ADAM modules
use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
! third party modules
use finer
use penf

implicit none
private
public :: IC_TYPE_GAUSS
public :: patch_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions" !< INI (config) file section name containing IC configs.
character(len=5),  parameter :: IC_TYPE_GAUSS="gauss"                 !< Gaussian distribution of rho IC TYPE parameter.

type :: patch_ic_object
   !< Initial Conditions class definition, CPU backend.
   type(mpih_object)         :: mpih                 !< MPI handler.
   integer(I4P)              :: amr_iterations=1_I4P !< Number of AMR iterations imposing IC.
   character(:), allocatable :: ic_type              !< IC type.
   integer(I4P)              :: regions_number=1_I4P !< Number of IC regions.
   real(R8P), allocatable    :: r(:,:)               !< Rho value.
   real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize IC.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_initial_conditions !< Set initial conditions.
endtype patch_ic_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(patch_ic_object), intent(in) :: self             !< IC.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: r, v             !< Counter.

   desc =       self%mpih%myrankstr//'IC main data'//NL
   desc = desc//self%mpih%myrankstr//'  regions number: '//trim(str(self%regions_number))
   do r=1, self%regions_number
      desc = desc//NL//self%mpih%myrankstr//'  region('//trim(str(r,.true.))//')'
      do v=1, size(self%r, dim=1)
         desc = desc//NL//self%mpih%myrankstr//'    r('//trim(str(v,.true.))//'): '//trim(str(self%r(v,r)))
      enddo
         desc = desc//NL//self%mpih%myrankstr//'    emin: '//trim(str(self%emin(:,r)))
         desc = desc//NL//self%mpih%myrankstr//'    emax: '//trim(str(self%emax(:,r)))
   enddo
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(patch_ic_object), intent(inout) :: self            !< IC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'patch_ic_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'patch_ic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(patch_ic_object), intent(inout)        :: self            !< IC.
   type(file_ini),         intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                      :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable                    :: sname           !< Section name.
   integer(I4P)                                 :: i               !< Counter.
   integer(I4P)                                 :: error           !< Error status.
   character(99)                                :: buff_char       !< Option character buffer.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='amr_iterations', val=self%amr_iterations, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(amr_iterations)')
   self%amr_iterations = max(0_I4P, self%amr_iterations)
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='type', val=buff_char, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(type)')
   self%ic_type = trim(adjustl(buff_char))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')

   if (self%regions_number>=1) then
      allocate(   self%r(1:1, 1:self%regions_number))
      allocate(self%emin(1:3, 1:self%regions_number))
      allocate(self%emax(1:3, 1:self%regions_number))
      do i=1, self%regions_number
         sname = INI_SECTION_NAME//'_region_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='r', val=self%r(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(r)')
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

   subroutine set_initial_conditions(self, field, r)
   !< Set initial conditions.
   class(patch_ic_object),     intent(in)    :: self       !< IC.
   type(field_object),         intent(in)    :: field      !< Field object.
   real(R8P),                  intent(inout) :: r(1:,               &
                                                  1-field%grid%ngc:,&
                                                  1-field%grid%ngc:,&
                                                  1-field%grid%ngc:,&
                                                  1:)      !< Conservative variables.
   integer(I4P)                              :: b, i, j, k !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, &
             x=>field%x_cell, y=>field%y_cell, z=>field%z_cell)
   select case(self%ic_type)
   case(IC_TYPE_GAUSS) ! Gauss distribution, only one region (s=1); r(1,1) is the base level
      do b=1, blocks_number
      do k=1, nk
      do j=1, nj
      do i=1, ni
         r(1,i,j,k,b) = exp(-((x(i,b)-0.5_R8P)**2 + (y(j,b)-0.5_R8P)**2 + (z(k,b)-0.5_R8P)**2) / 0.01_R8P)
      enddo
      enddo
      enddo
      enddo
   endselect
   endassociate
   endsubroutine set_initial_conditions
endmodule adam_patch_ic_object
