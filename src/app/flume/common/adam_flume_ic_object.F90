!< ADAM, FLUME initial conditions class definition.
module adam_flume_ic_object
!< ADAM, FLUME initial conditions class definition.
!<
!< Accepted `[initial_conditions].(type)`: `uniform` (the state of region 1 everywhere) and `riemann-problem`
!< (piecewise-constant axis-aligned regions, a cell belongs to a region iff `emin < center <= emax` on every axis).
!< A cell covered by no region is fatal: leaving it at zero density would divide by zero downstream.

! ADAM classes, libraries, parameters
use :: adam_field_object,         only : field_object
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_euler_library,  only : primitive_to_conservative
use :: adam_flume_parameters,     only : NV_EULER, strip_control
use :: adam_flume_physics_object, only : flume_physics_object
! third party modules
use :: finer,                     only : file_ini
use :: penf,                      only : I4P, R8P, str

implicit none
private
public :: flume_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions"   !< INI section name.
character(len=7),  parameter :: IC_UNIFORM_STR="uniform"                 !< Uniform state.
character(len=15), parameter :: IC_RIEMANN_PROBLEM_STR="riemann-problem" !< Piecewise-constant regions.
character(len=1),  parameter :: PRIM_KEY(5)=['r', 'u', 'v', 'w', 'p']    !< Region primitive state keys.
character(len=6),  parameter :: EXTENT_KEY(6)=['emin_x', 'emin_y', &
                                               'emin_z', 'emax_x', &
                                               'emax_y', 'emax_z']       !< Region extents keys.

type :: flume_ic_object
   !< FLUME initial conditions class definition.
   integer(I4P)              :: amr_iterations=0_I4P !< AMR iterations performed while imposing the initial conditions.
   character(:), allocatable :: ic_type              !< Initial conditions type.
   integer(I4P)              :: regions_number=0_I4P !< Regions number.
   real(R8P),    allocatable :: q_region(:,:)        !< Conservative state of each region [NV_EULER, regions_number].
   real(R8P),    allocatable :: emin(:,:)            !< Minimum corner of each region [3, regions_number].
   real(R8P),    allocatable :: emax(:,:)            !< Maximum corner of each region [3, regions_number].
   contains
      ! public methods
      procedure, pass(self) :: description            !< Return pretty-printed object description.
      procedure, pass(self) :: initialize             !< Initialize initial conditions.
      procedure, pass(self) :: load_from_file         !< Load config from file.
      procedure, pass(self) :: set_initial_conditions !< Set initial conditions on the blocks interior.
endtype flume_ic_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_ic_object), intent(in) :: self             !< Initial conditions.
   character(len=:), allocatable      :: desc             !< Description.
   character(len=1), parameter        :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Initial conditions main data'//NL
   desc = desc//mpih%myrankstr//'  type:           '//self%ic_type//NL
   desc = desc//mpih%myrankstr//'  amr_iterations: '//trim(str(self%amr_iterations))//NL
   desc = desc//mpih%myrankstr//'  regions_number: '//trim(str(self%regions_number))
   endfunction description

   subroutine initialize(self, file_parameters, physics)
   !< Initialize initial conditions.
   class(flume_ic_object),     intent(inout) :: self            !< Initial conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the regions state conversion).

   print '(A)', mpih%myrankstr//'flume_ic_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters, physics=physics)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_ic_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, physics)
   !< Load config from file; every key used by the selected type is required and an unknown type is fatal.
   class(flume_ic_object),     intent(inout) :: self            !< Initial conditions.
   type(file_ini),             intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   type(flume_physics_object), intent(in)    :: physics         !< Physics (for the regions state conversion).
   character(999)                            :: buff            !< Option value buffer.
   character(:), allocatable                 :: sname           !< Region section name.
   real(R8P)                                 :: prim(5)         !< Region primitive state (r, u, v, w, p).
   real(R8P)                                 :: extent(6)       !< Region extents.
   integer(I4P)                              :: error           !< Error status.
   integer(I4P)                              :: r, k            !< Counters.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='amr_iterations', val=self%amr_iterations, &
                            error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(amr_iterations)')
   self%amr_iterations = max(0_I4P, self%amr_iterations)
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='type', val=buff, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(type)')
   self%ic_type = trim(adjustl(strip_control(buff)))
   select case(self%ic_type)
   case(IC_UNIFORM_STR)
      self%regions_number = 1_I4P
   case(IC_RIEMANN_PROBLEM_STR)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, &
                               error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')
      if (self%regions_number < 1_I4P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(regions_number) must be positive')
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(type) "'//self%ic_type//'"; expected one of '// &
                               IC_UNIFORM_STR//', '//IC_RIEMANN_PROBLEM_STR)
   endselect

   if (allocated(self%q_region)) deallocate(self%q_region)
   if (allocated(self%emin)) deallocate(self%emin)
   if (allocated(self%emax)) deallocate(self%emax)
   allocate(self%q_region(NV_EULER,self%regions_number), self%emin(3,self%regions_number), &
            self%emax(3,self%regions_number))
   self%emin = -huge(1._R8P)
   self%emax =  huge(1._R8P)
   do r=1, self%regions_number
      sname = INI_SECTION_NAME//'_region_'//trim(str(r, .true.))
      do k=1, 5
         call file_parameters%get(section_name=sname, option_name=PRIM_KEY(k), val=prim(k), error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//sname//'].('//PRIM_KEY(k)//')')
      enddo
      call primitive_to_conservative(gamma=physics%gamma, r=prim(1), u=prim(2), v=prim(3), w=prim(4), p=prim(5), &
                                     q=self%q_region(:,r))
      if (self%ic_type == IC_RIEMANN_PROBLEM_STR) then
         do k=1, 6
            call file_parameters%get(section_name=sname, option_name=EXTENT_KEY(k), val=extent(k), error=error)
            if (error > 0) call mpih%error_stop(msg=': failed to load ['//sname//'].('//EXTENT_KEY(k)//')')
         enddo
         self%emin(:,r) = extent(1:3)
         self%emax(:,r) = extent(4:6)
      endif
   enddo
   endsubroutine load_from_file

   subroutine set_initial_conditions(self, field, q)
   !< Set initial conditions on the blocks interior; ghost cells are filled by the following ghost update.
   class(flume_ic_object), intent(in)    :: self          !< Initial conditions.
   type(field_object),     intent(in)    :: field         !< Field (realm component, threaded in).
   real(R8P),              intent(inout) :: q(1:,           &
                                              1-field%ngc:, &
                                              1-field%ngc:, &
                                              1-field%ngc:, &
                                              1:)           !< Conservative variables.
   real(R8P)                             :: center(3)     !< Cell center.
   logical                               :: is_set        !< Flag: cell covered by a region.
   integer(I4P)                          :: b, i, j, k, r !< Counters.

   do b=1, field%blocks_number
      do k=1, field%nk
         do j=1, field%nj
            do i=1, field%ni
               center = [field%x_cell(i,b), field%y_cell(j,b), field%z_cell(k,b)]
               is_set = .false.
               do r=1, self%regions_number
                  if (all(center > self%emin(:,r)) .and. all(center <= self%emax(:,r))) then
                     q(:,i,j,k,b) = self%q_region(:,r)
                     is_set = .true.
                     exit
                  endif
               enddo
               if (.not.is_set) &
                  call mpih%error_stop(msg=': cell center ('//trim(str(center(1)))//', '//trim(str(center(2)))//', '// &
                                           trim(str(center(3)))//') is covered by no ['//INI_SECTION_NAME// &
                                           '_region_*]')
            enddo
         enddo
      enddo
   enddo
   endsubroutine set_initial_conditions
endmodule adam_flume_ic_object
