!< ADAM, FLUME external gravity: the `[gravity]` input (issue #73, M7).
module adam_flume_gravity_object
!< ADAM, FLUME external gravity: the `[gravity]` input (issue #73, M7).
!<
!< A constant acceleration `g`, given either as its components or as the number it stands for (issue #49, section 4.5):
!<
!< | form | keys | g (code units) |
!< |---|---|---|
!< | acceleration | `g_x`, `g_y`, `g_z` (each optional, default 0) | `(g_x, g_y, g_z)` |
!< | number | `froude` (> 0) and `axis` (`+x`, `-x`, `+y`, `-y`, `+z`, `-z`) | `1 / froude^2` along the signed axis |
!<
!< Giving both forms is fatal, as is a `[gravity]` section with neither. In code units the references are 1, so the
!< number is `1 / Fr^2 = g L0 / u0^2`; the `[reference]` layer converts the dimensional components (length / time^2)
!< and leaves `froude` as it is. `equilibrium` (the well-balanced mode, M7-G2) accepts only `none` and `potential` only
!< `uniform` until their phases land: any other value is fatal, never ignored.
!<
!< Without the section there is no gravity (`is_active` false) and nothing changes. With it, `g` sets the atmosphere
!< initial conditions and the `hydrostatic` boundary condition (M7-G0) and drives the momentum and energy sources
!< (M7-G1): `rho g`, and the energy from the face mass fluxes, the fourth-order
!< `sum_d g_d (-F_{i-3/2} + 13 F_{i-1/2} + 13 F_{i+1/2} - F_{i+3/2}) / 24` (M7-G2a; the two-face average alone is second
!< order), which reads the face beyond each block end: gravity needs one more ghost layer than the WENO stencil.

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : strip_control
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P, str

implicit none
private
public :: flume_gravity_object

character(len=7), parameter :: INI_SECTION_NAME="gravity"                  !< INI section name.
character(len=3), parameter :: G_KEY(3)=['g_x', 'g_y', 'g_z']              !< Acceleration component keys.
character(len=2), parameter :: AXIS_STR(6)=['+x', '-x', '+y', '-y', '+z', '-z'] !< Accepted `axis` spellings.
character(len=4), parameter :: EQUILIBRIUM_NONE="none"                     !< `equilibrium`: no well-balanced mode.
character(len=7), parameter :: POTENTIAL_UNIFORM="uniform"                 !< `potential`: constant acceleration.

type :: flume_gravity_object
   !< FLUME external gravity, code units.
   real(R8P)                 :: g(3)=0._R8P   !< Constant acceleration.
   real(R8P)                 :: froude=0._R8P !< Froude number, 0 when g is given by its components.
   character(:), allocatable :: equilibrium   !< Well-balanced equilibrium (M7-G2): none.
   character(:), allocatable :: potential     !< Potential kind: uniform.
   logical                   :: is_active=.false. !< `[gravity]` given.
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: load_from_file !< Load the gravity from file.
endtype flume_gravity_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_gravity_object), intent(in) :: self             !< Gravity.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   if (.not.self%is_active) then
      desc = mpih%myrankstr//'  gravity:         none'
      return
   endif
   desc = mpih%myrankstr//'  gravity g:       '//trim(str(self%g))
   if (self%froude > 0._R8P) desc = desc//' (froude '//trim(str(self%froude))//')'
   desc = desc//NL//mpih%myrankstr//'  gravity:         equilibrium '//self%equilibrium//', potential '// &
                    self%potential
   endfunction description

   subroutine load_from_file(self, file_parameters)
   !< Load the gravity from file: the components or `froude` + `axis` (both fatal, neither fatal), `equilibrium` and
   !< `potential` (optional, only `none` and `uniform` until M7-G2).
   class(flume_gravity_object), intent(inout) :: self            !< Gravity.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(999)                             :: buff            !< Option value buffer.
   character(:), allocatable                  :: axis            !< Signed axis of the number form.
   logical                                    :: has_component   !< Any `g_*` given.
   integer(I4P)                               :: error           !< Error status.
   integer(I4P)                               :: d               !< Counter.

   self%g = 0._R8P
   self%froude = 0._R8P
   self%equilibrium = EQUILIBRIUM_NONE
   self%potential = POTENTIAL_UNIFORM
   self%is_active = file_parameters%has_section(section_name=INI_SECTION_NAME)
   if (.not.self%is_active) return
   has_component = .false.
   do d=1, 3
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=G_KEY(d), val=self%g(d), error=error)
      if (error > 0) then
         self%g(d) = 0._R8P
      else
         has_component = .true.
      endif
   enddo
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='froude', val=self%froude, error=error)
   if (error > 0) self%froude = 0._R8P
   if (error <= 0) then
      if (has_component) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] takes either g_x, g_y, g_z or froude '// &
                                                  'and axis, not both')
      if (.not.(self%froude > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(froude) must be positive')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='axis', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(froude) needs ['//INI_SECTION_NAME// &
                                              '].(axis), one of +x, -x, +y, -y, +z, -z')
      axis = trim(adjustl(strip_control(buff)))
      d = findloc(AXIS_STR, axis, dim=1)
      if (d == 0) call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(axis) "'//axis//'"; expected one '// &
                                           'of +x, -x, +y, -y, +z, -z')
      self%g((d+1)/2) = merge(1._R8P, -1._R8P, mod(d, 2) == 1) / self%froude**2
   elseif (.not.has_component) then
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs g_x, g_y, g_z or froude and axis')
   endif
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='equilibrium', val=buff, error=error)
   if (error <= 0) self%equilibrium = trim(adjustl(strip_control(buff)))
   if (self%equilibrium /= EQUILIBRIUM_NONE) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(equilibrium) "'//self%equilibrium//'" is not available '// &
                               '(the well-balanced mode lands in issue #73, M7-G2); expected none')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='potential', val=buff, error=error)
   if (error <= 0) self%potential = trim(adjustl(strip_control(buff)))
   if (self%potential /= POTENTIAL_UNIFORM) &
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(potential) "'//self%potential//'"; expected '// &
                               POTENTIAL_UNIFORM)
   endsubroutine load_from_file
endmodule adam_flume_gravity_object
