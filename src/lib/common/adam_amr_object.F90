!< ADAM, AMR markers class definition, CPU backend.
module adam_amr_object
!< ADAM, AMR markers class definition, CPU backend.

! ADAM singleton objects
use :: adam_mpih_global, only : mpih
use :: adam_grid_object, only : grid_object
! third party modules
use :: finer
use :: penf

implicit none
private
public :: amr_object
public :: amr_marker_object
public :: AMR_GEO
public :: AMR_GRAD
public :: AMR_TV
public :: AMR_GEO_SOLID
public :: AMR_GEO_PRIMITIVE_BOX
public :: AMR_GEO_STL
public :: AMR_DELTA_T_X
public :: AMR_DELTA_T_Y
public :: AMR_DELTA_T_Z
public :: AMR_DELTA_T_MAX

character(len=3), parameter :: INI_SECTION_NAME="amr" !< INI (config) file section name containing AMR markers configs.

! AMR_GEO geometry kinds. Offset to 100+ so they never numerically collide with the marker
! modes (AMR_GEO/AMR_GRAD/AMR_TV = 1/2/3) — the two enums occupy distinct integer ranges.
integer(I4P), parameter :: AMR_GEO               = 1_I4P   !< Geometrical marker.
integer(I4P), parameter :: AMR_GRAD              = 2_I4P   !< Field gradient marker.
integer(I4P), parameter :: AMR_TV                = 3_I4P   !< Field total variation marker.
integer(I4P), parameter :: AMR_GEO_SOLID         = 101_I4P !< AMR_GEO geometry kind: solid/IB index (legacy default).
integer(I4P), parameter :: AMR_GEO_PRIMITIVE_BOX = 102_I4P !< AMR_GEO geometry kind: primitive axis-aligned box.
integer(I4P), parameter :: AMR_GEO_STL           = 103_I4P !< AMR_GEO geometry kind: STL surface file.
character(1), parameter :: AMR_DELTA_T_X         = 'x'     !< Delta criterion type, dx.
character(1), parameter :: AMR_DELTA_T_Y         = 'y'     !< Delta criterion type, dy.
character(1), parameter :: AMR_DELTA_T_Z         = 'z'     !< Delta criterion type, dz.
character(3), parameter :: AMR_DELTA_T_MAX       = 'max'   !< Delta criterion type, max(dx,dy,dz).

type :: amr_marker_object
   !< AMR marker object.
   integer(I4P)              :: mode=AMR_GRAD          !< Marker mode.
   character(:), allocatable :: delta_type             !< Type of delta criterion: 'x','y','z','max' (max=> max(x,y,z)).
   real(R8P)                 :: delta_fine             !< Fine cell space step.
   real(R8P)                 :: delta_coarse           !< Coarse cell space step.
   integer(I4P)              :: field=1_I4P            !< Field array containing the marker variable, 1=q, 2=q_aux.
   integer(I4P)              :: ivar=1_I4P             !< Array index of variable used by refinement criterion.
   real(R8P)                 :: tol=0.5_R8P            !< Tolerance.
   integer(I4P)              :: solid=1_I4P            !< Solid number (AMR_GEO + AMR_GEO_SOLID).
   integer(I4P)              :: geo_type=AMR_GEO_SOLID !< AMR_GEO geometry kind: solid|primitive|stl (default solid => legacy).
   real(R8P)                 :: box_emin(3)=0._R8P     !< Primitive box minimum corner (AMR_GEO + AMR_GEO_PRIMITIVE_BOX).
   real(R8P)                 :: box_emax(3)=0._R8P     !< Primitive box maximum corner (AMR_GEO + AMR_GEO_PRIMITIVE_BOX).
   integer(I4P)              :: target_level=1_I4P     !< Refine blocks below this level (AMR_GEO + AMR_GEO_PRIMITIVE_BOX).
   character(:), allocatable :: stl_filename           !< STL surface file (AMR_GEO + AMR_GEO_STL).
endtype amr_marker_object

type :: amr_object
   !< AMR markers class definition, CPU backend.
   integer(I4P)                         :: iters=5_I4P          !< AMR updates iterations number.
   integer(I4P)                         :: frequency=100_I4P    !< AMR update time step frequency.
   integer(I4P)                         :: markers_number=0_I4P !< AMR number of markers.
   type(amr_marker_object), allocatable :: markers(:)           !< AMR array of marker objects.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize IC.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype amr_object

contains
   ! public methods
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(amr_object), intent(in) :: self             !< The tree.
   character(len=:), allocatable :: desc             !< Description.
   character(len=1), parameter   :: NL=new_line('a') !< New line character.
   integer(I4P)                  :: m                !< Counter.

   desc =       mpih%myrankstr//'amr main data'                                    //NL
   desc = desc//mpih%myrankstr//'  iters:          '//trim(str(self%iters         ))//NL
   desc = desc//mpih%myrankstr//'  frequency:      '//trim(str(self%frequency     ))//NL
   desc = desc//mpih%myrankstr//'  markers number: '//trim(str(self%markers_number))
   if (self%markers_number>0) then
   do m=1, self%markers_number
   desc = desc//NL//mpih%myrankstr//'    delta type:      '//trim(    self%markers(m)%delta_type   )
   desc = desc//NL//mpih%myrankstr//'    delta fine:      '//trim(str(self%markers(m)%delta_fine)  )
   desc = desc//NL//mpih%myrankstr//'    delta coarse:    '//trim(str(self%markers(m)%delta_coarse))
   desc = desc//NL//mpih%myrankstr//'    ivar:            '//trim(str(self%markers(m)%ivar        ))
   desc = desc//NL//mpih%myrankstr//'    tol:             '//trim(str(self%markers(m)%tol         ))
   desc = desc//NL//mpih%myrankstr//'    solid:           '//trim(str(self%markers(m)%solid       ))
   desc = desc//NL//mpih%myrankstr//'    geo_type:        '//trim(str(self%markers(m)%geo_type    ))
   enddo
   endif
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(amr_object), intent(inout) :: self            !< AMR.
   type(file_ini),    intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call mpih%print_message('amr_object%initialize start')
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   call mpih%print_message('amr_object%initialize finish')
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(amr_object), intent(inout)        :: self            !< AMR.
   type(file_ini),    intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,           intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                 :: go_on_fail_     !< Go on if load fails.
   character(:), allocatable               :: sname           !< Section name.
   character(99)                           :: buff_c          !< Character buffer.
   integer(I4P)                            :: i_marker        !< Counter.
   integer(I4P)                            :: error           !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='frequency', val=self%frequency, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(frequency)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='iters', val=self%iters, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(iters)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='markers_number', val=self%markers_number, error=error)
   if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(markers_number)')

   allocate(self%markers(self%markers_number))
   do i_marker=1, self%markers_number
      sname = INI_SECTION_NAME//'_marker_'//trim(str(i_marker,.true.))
      call file_parameters%get(section_name=sname, option_name='mode', val=self%markers(i_marker)%mode, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(mode)')
      call file_parameters%get(section_name=sname, option_name='delta_type', val=buff_c, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(delta_type)')
      self%markers(i_marker)%delta_type = trim(adjustl(buff_c))
      call file_parameters%get(section_name=sname, option_name='delta_fine', val=self%markers(i_marker)%delta_fine, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(delta_fine)')
      call file_parameters%get(section_name=sname, option_name='delta_coarse', val=self%markers(i_marker)%delta_coarse, error=error)
      if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(delta_coarse)')
      select case(self%markers(i_marker)%mode)
      case(AMR_GEO)
         call file_parameters%get(section_name=sname, option_name='geo_type', val=buff_c, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(geo_type)')
         select case(trim(adjustl(buff_c)))
         case('solid')         ; self%markers(i_marker)%geo_type = AMR_GEO_SOLID
         case('primitive-box') ; self%markers(i_marker)%geo_type = AMR_GEO_PRIMITIVE_BOX
         case('stl')           ; self%markers(i_marker)%geo_type = AMR_GEO_STL
         case default
            call mpih%error_stop(msg=': unknown AMR_GEO geo_type "'//trim(adjustl(buff_c))// &
                                     '" in ['//sname//'] (expected: solid | primitive-box | stl)')
         endselect
         select case(self%markers(i_marker)%geo_type)
         case(AMR_GEO_SOLID)
            call file_parameters%get(section_name=sname, option_name='solid', val=self%markers(i_marker)%solid, error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(solid)')
         case(AMR_GEO_PRIMITIVE_BOX)
            call file_parameters%get(section_name=sname, option_name='box_xmin', &
                                     val=self%markers(i_marker)%box_emin(1), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_xmin)')
            call file_parameters%get(section_name=sname, option_name='box_ymin', &
                                     val=self%markers(i_marker)%box_emin(2), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_ymin)')
            call file_parameters%get(section_name=sname, option_name='box_zmin', &
                                     val=self%markers(i_marker)%box_emin(3), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_zmin)')
            call file_parameters%get(section_name=sname, option_name='box_xmax', &
                                     val=self%markers(i_marker)%box_emax(1), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_xmax)')
            call file_parameters%get(section_name=sname, option_name='box_ymax', &
                                     val=self%markers(i_marker)%box_emax(2), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_ymax)')
            call file_parameters%get(section_name=sname, option_name='box_zmax', &
                                     val=self%markers(i_marker)%box_emax(3), error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(box_zmax)')
            call file_parameters%get(section_name=sname, option_name='target_level', &
                                     val=self%markers(i_marker)%target_level, error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(target_level)')
         case(AMR_GEO_STL)
            call file_parameters%get(section_name=sname, option_name='stl_filename', val=buff_c, error=error)
            if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(stl_filename)')
            self%markers(i_marker)%stl_filename = trim(adjustl(buff_c))
         endselect
      case(AMR_GRAD,AMR_TV)
         call file_parameters%get(section_name=sname, option_name='field', val=self%markers(i_marker)%field, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(field)')
         call file_parameters%get(section_name=sname, option_name='var', val=self%markers(i_marker)%ivar, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(var)')
         call file_parameters%get(section_name=sname, option_name='tol', val=self%markers(i_marker)%tol, error=error)
         if (.not.go_on_fail_.and.error>0) call mpih%error_stop(msg=': failed to load ['//sname//'].(tol)')
      endselect
   enddo
   endsubroutine load_from_file
endmodule adam_amr_object
