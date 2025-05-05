!< ADAM, PRISM Initial Conditions class definition, CPU backend.
module adam_prism_ic_object 
    !< ADAM, PRISM Initial Conditions class definition, CPU backend.

use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
use adam_prism_physics_object, only : prism_physics_object
use finer
use penf

implicit none
private
public :: IC_TYPE_RP
public :: IC_TYPE_VACUUM
public :: prism_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions" !< INI (config) file section name containing IC configs.
character(len=6),  parameter :: IC_TYPE_VACUUM="vacuum"               !< Vacuum IC TYPE parameter.
character(len=15), parameter :: IC_TYPE_RP="riemann-problem"          !< Riemann Problem IC TYPE parameter.

type :: prism_ic_object
   !< Initial Conditions class definition, CPU backend.
   type(mpih_object)         :: mpih                 !< MPI handler.
   integer(I4P)              :: amr_iterations=1_I4P !< Number of AMR iterations imposing IC.
   character(:), allocatable :: ic_type              !< IC type.
   integer(I4P)              :: regions_number=1_I4P !< Number of IC regions.
   real(R8P), allocatable    :: q(:,:)               !< Primitive variables (Dx,Dy,Dz,Bx,By,Bz,Jx,Jy,Jz).
   real(R8P), allocatable    :: emin(:,:), emax(:,:) !< IC regions bounding box.
   contains
      ! public methods
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize IC.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_initial_conditions !< Set initial conditions on PRISM fields.
endtype prism_ic_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(prism_ic_object), intent(in) :: self             !< IC.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.
   integer(I4P)                       :: r, v             !< Counter.

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
   class(prism_ic_object), intent(inout) :: self            !< IC.
   type(file_ini),         intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'prism_ic_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'prism_ic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(prism_ic_object), intent(inout)        :: self            !< IC.
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
      allocate(   self%q(1:9, 1:self%regions_number))
      allocate(self%emin(1:3, 1:self%regions_number))
      allocate(self%emax(1:3, 1:self%regions_number))
      do i=1, self%regions_number
         sname = INI_SECTION_NAME//'_region_'//trim(str(i,.true.))
         call file_parameters%get(section_name=sname, option_name='Dx', val=self%q(1,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dx)')
         call file_parameters%get(section_name=sname, option_name='Dy', val=self%q(2,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dy)')
         call file_parameters%get(section_name=sname, option_name='Dz', val=self%q(3,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Dz)')
         call file_parameters%get(section_name=sname, option_name='Bx', val=self%q(4,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Bx)')
         call file_parameters%get(section_name=sname, option_name='By', val=self%q(5,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(By)')
         call file_parameters%get(section_name=sname, option_name='Bz', val=self%q(6,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Bz)')

         call file_parameters%get(section_name=sname, option_name='Jx', val=self%q(7,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jx)')
         call file_parameters%get(section_name=sname, option_name='Jy', val=self%q(8,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jy)')
         call file_parameters%get(section_name=sname, option_name='Jz', val=self%q(9,i), error=error)
         if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//sname//'].(Jz)')

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

   subroutine set_initial_conditions(self, physics, field)  
      !< Set initial conditions on PRISM fields.
      class(prism_ic_object),     intent(in)    :: self                 !< IC.
      type(prism_physics_object), intent(in)    :: physics              !< Fluids physiscs.
      type(field_object),         intent(inout) :: field                !< Field object.
      integer(I4P)                              :: b, i, j, k, ri       !< Counter.

   associate(blocks_number=>field%blocks_number, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, ngc=>field%grid%ngc, &
             q=>field%q, x_cell=>field%x_cell, y_cell=>field%y_cell, z_cell=>field%z_cell)
   select case(self%ic_type)
   case(IC_TYPE_VACUUM) ! vacuum initial conditions
      do b=1, blocks_number
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  q(1,i,j,k,b) = 0.0_R8P
                  q(2,i,j,k,b) = 0.0_R8P
                  q(3,i,j,k,b) = 0.0_R8P
                  q(4,i,j,k,b) = 0.0_R8P
                  q(5,i,j,k,b) = 0.0_R8P
                  q(6,i,j,k,b) = 0.0_R8P
                  q(7,i,j,k,b) = 0.0_R8P
                  q(8,i,j,k,b) = 0.0_R8P
                  q(9,i,j,k,b) = 0.0_R8P
               enddo
            enddo
         enddo
      enddo
   !case(IC_TYPE_RP)      ! Riemann Problem like (uniform regions)
   !   do b=1, blocks_number
   !      do k=1, nk
   !         do j=1, nj
   !            do i=1, ni
   !               do ri=1, self%regions_number
   !                  if ((x_cell(i,b) > self%emin(1,ri).and.x_cell(i,b) <= self%emax(1,ri)).and. &
   !                      (y_cell(j,b) > self%emin(2,ri).and.y_cell(j,b) <= self%emax(2,ri)).and. &
   !                      (z_cell(k,b) > self%emin(3,ri).and.z_cell(k,b) <= self%emax(3,ri))) then
   !                     q(1:5,i,j,k,b) = physics%eos(int(self%q(6,ri)))%primitive2conservative(primitive=self%q(1:5,ri))
   !                  endif
   !               enddo
   !            enddo
   !         enddo
   !      enddo
   !   enddo
   endselect
   endassociate
   endsubroutine set_initial_conditions
endmodule adam_prism_ic_object