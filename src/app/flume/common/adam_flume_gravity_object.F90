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
!< and leaves `froude` as it is. `equilibrium` selects the well-balanced mode (M7-G2b): `none` (default) or `atmosphere`
!< (the atmosphere of the initial conditions, at rest, is kept to round-off: the residual of that state is subtracted
!< from every residual). `potential` selects how g reaches the sources (issue #81, M8-P1): `uniform` (default), the
!< constant `g` in every kernel, with the fourth-order energy source; or `linear`, the same field written as the cell
!< potential `Phi = -g . x` and applied through the potential path of self-gravity (face gravity `-dPhi/dx`, cell gravity
!< the mean of its two faces, second order): a test of that path against the constant one. Any other value is fatal,
!< never ignored.
!<
!< **Self-gravity** (issue #81, M8-P3): `self_gravity = .true.` solves `L_h Phi = 4 pi G (rho - mean rho)` every
!< Runge-Kutta stage (the library FFT solver, `poisson_solver = fft`; `poisson_bc = periodic`: every active direction
!< periodic, the grid uniform) and applies it through the potential path. `poisson_bc = isolated` (M8-P4) solves
!< `L_h Phi = 4 pi G rho` with `Phi -> 0` at infinity instead (zero padding and the lattice Green's function: 3-D,
!< every direction non-periodic, cubic cells; the solver refuses anything else). The constant is `G` (dimensional: the
!< `[reference]` layer converts it, dimension `u0^2 / (rho0 L0^2)`) or `G_star` (the number `G rho0 L0^2 / u0^2`, never
!< converted); exactly one of them is required with `self_gravity`, and either without it is fatal. The constant `g` is
!< optional with `self_gravity` (default 0) and is added through its own kernels (`potential = uniform`); `potential =
!< linear` and `equilibrium = atmosphere` are refused with it.
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
public :: EQUILIBRIUM_ATMOSPHERE, EQUILIBRIUM_NONE

character(len=7), parameter :: INI_SECTION_NAME="gravity"                  !< INI section name.
character(len=3), parameter :: G_KEY(3)=['g_x', 'g_y', 'g_z']              !< Acceleration component keys.
character(len=2), parameter :: AXIS_STR(6)=['+x', '-x', '+y', '-y', '+z', '-z'] !< Accepted `axis` spellings.
character(len=4),  parameter :: EQUILIBRIUM_NONE="none"                    !< `equilibrium`: no well-balanced mode.
character(len=10), parameter :: EQUILIBRIUM_ATMOSPHERE="atmosphere"        !< `equilibrium`: the initial atmosphere.
character(len=7), parameter :: POTENTIAL_UNIFORM="uniform"                 !< `potential`: constant acceleration.
character(len=6), parameter :: POTENTIAL_LINEAR="linear"                   !< `potential`: `Phi = -g . x` (M8-P1).
character(len=3), parameter :: POISSON_SOLVER_FFT="fft"                    !< `poisson_solver`: the FFT solver (M8-P2).
character(len=8), parameter :: POISSON_BC_PERIODIC="periodic"              !< `poisson_bc`: periodic (mean removed).
character(len=8), parameter :: POISSON_BC_ISOLATED="isolated"              !< `poisson_bc`: isolated (issue #81, P4).
real(R8P),        parameter :: PI=acos(-1._R8P)                            !< Pi.

type :: flume_gravity_object
   !< FLUME external gravity, code units.
   real(R8P)                 :: g(3)=0._R8P   !< Constant acceleration (the field the atmospheres and `linear` use).
   real(R8P)                 :: g_uniform(3)=0._R8P !< The part of g the constant-g kernels apply: g, or 0 with `linear`.
   real(R8P)                 :: froude=0._R8P !< Froude number, 0 when g is given by its components.
   character(:), allocatable :: equilibrium   !< Well-balanced equilibrium (M7-G2b): none or atmosphere.
   character(:), allocatable :: potential     !< Potential kind: uniform.
   logical                   :: is_active=.false. !< `[gravity]` given.
   logical                   :: is_balanced=.false. !< Well-balanced mode on (`equilibrium` not none).
   logical                   :: has_potential=.false. !< The sources read a cell potential (`linear`, self-gravity).
   logical                   :: self_gravity=.false.  !< Self-gravity on (M8-P3).
   real(R8P)                 :: grav_const=0._R8P     !< Gravitational constant G, code units (self-gravity).
   real(R8P)                 :: four_pi_g=0._R8P      !< 4 pi G, the source factor of the Poisson equation.
   character(:), allocatable :: poisson_solver        !< Poisson solver: fft.
   character(:), allocatable :: poisson_bc            !< Poisson boundary conditions: periodic or isolated.
   logical                   :: isolated=.false.      !< `poisson_bc = isolated` (M8-P4).
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
   if (self%self_gravity) desc = desc//NL//mpih%myrankstr//'  self-gravity:    G '//trim(str(self%grav_const))// &
                                 ', 4 pi G '//trim(str(self%four_pi_g))//', solver '//self%poisson_solver//', bc '// &
                                 self%poisson_bc
   endfunction description

   subroutine load_from_file(self, file_parameters)
   !< Load the gravity from file: the components or `froude` + `axis` (both fatal, neither fatal unless
   !< `self_gravity`), `equilibrium` and `potential` (optional: `none` or `atmosphere`, and `uniform` or `linear`), and
   !< the self-gravity keys `self_gravity`, `G` or `G_star`, `poisson_solver`, `poisson_bc`.
   class(flume_gravity_object), intent(inout) :: self            !< Gravity.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   character(999)                             :: buff            !< Option value buffer.
   character(:), allocatable                  :: axis            !< Signed axis of the number form.
   logical                                    :: has_component   !< Any `g_*` given.
   logical                                    :: has_g, has_g_star !< `G`, `G_star` given.
   real(R8P)                                  :: g_star          !< `G_star`.
   integer(I4P)                               :: error           !< Error status.
   integer(I4P)                               :: d               !< Counter.

   self%g = 0._R8P
   self%froude = 0._R8P
   self%equilibrium = EQUILIBRIUM_NONE
   self%potential = POTENTIAL_UNIFORM
   self%is_balanced = .false.
   self%has_potential = .false.
   self%g_uniform = 0._R8P
   self%self_gravity = .false.
   self%grav_const = 0._R8P
   self%four_pi_g = 0._R8P
   self%poisson_solver = POISSON_SOLVER_FFT
   self%poisson_bc = POISSON_BC_PERIODIC
   self%isolated = .false.
   self%is_active = file_parameters%has_section(section_name=INI_SECTION_NAME)
   if (.not.self%is_active) return
   ! the string keys first: under nvfortran (26.1) a character get that follows the failed get of an absent numeric key
   ! (an unused g_* or froude) returns an empty value (FiNeR, class(*) dummy); found by M7-G2b on FNL
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='equilibrium', val=buff, error=error)
   if (error <= 0) self%equilibrium = trim(adjustl(strip_control(buff)))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='potential', val=buff, error=error)
   if (error <= 0) self%potential = trim(adjustl(strip_control(buff)))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='poisson_solver', val=buff, error=error)
   if (error <= 0) self%poisson_solver = trim(adjustl(strip_control(buff)))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='poisson_bc', val=buff, error=error)
   if (error <= 0) self%poisson_bc = trim(adjustl(strip_control(buff)))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='self_gravity', val=self%self_gravity, error=error)
   if (error > 0) self%self_gravity = .false.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='G', val=self%grav_const, error=error)
   has_g = error <= 0
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='G_star', val=g_star, error=error)
   has_g_star = error <= 0
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
   elseif (.not.has_component .and. .not.self%self_gravity) then
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs g_x, g_y, g_z or froude and axis (or self_gravity)')
   endif
   if (self%equilibrium /= EQUILIBRIUM_NONE .and. self%equilibrium /= EQUILIBRIUM_ATMOSPHERE) &
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(equilibrium) "'//self%equilibrium//'"; expected '// &
                               EQUILIBRIUM_NONE//' or '//EQUILIBRIUM_ATMOSPHERE)
   self%is_balanced = self%equilibrium == EQUILIBRIUM_ATMOSPHERE
   if (self%potential /= POTENTIAL_UNIFORM .and. self%potential /= POTENTIAL_LINEAR) &
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(potential) "'//self%potential//'"; expected '// &
                               POTENTIAL_UNIFORM//' or '//POTENTIAL_LINEAR)
   if (self%potential == POTENTIAL_LINEAR) then
      self%g_uniform = 0._R8P
   else
      self%g_uniform = self%g
   endif
   ! self-gravity (issue #81, M8-P3)
   if (self%self_gravity) then
      if (has_g .and. has_g_star) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] takes either G or G_star, '// &
                                                           'not both')
      if (.not.(has_g .or. has_g_star)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(self_gravity) needs '// &
                                                                 'G or G_star')
      if (has_g_star) self%grav_const = g_star
      if (.not.(self%grav_const > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(G, G_star) must be '// &
                                                                    'positive')
      self%four_pi_g = 4._R8P * PI * self%grav_const
      if (self%poisson_solver /= POISSON_SOLVER_FFT) &
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(poisson_solver) "'//self%poisson_solver// &
                                  '"; expected '//POISSON_SOLVER_FFT)
      if (self%poisson_bc /= POISSON_BC_PERIODIC .and. self%poisson_bc /= POISSON_BC_ISOLATED) &
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(poisson_bc) "'//self%poisson_bc// &
                                  '"; expected '//POISSON_BC_PERIODIC//' or '//POISSON_BC_ISOLATED)
      self%isolated = self%poisson_bc == POISSON_BC_ISOLATED
      ! a linear potential is not periodic: the FFT solver's periodic wrap cannot hold it
      if (self%potential == POTENTIAL_LINEAR) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(potential) = '//POTENTIAL_LINEAR//' is refused with '// &
                                  '['//INI_SECTION_NAME//'].(self_gravity) (issue #81)')
   elseif (has_g .or. has_g_star) then
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(G, G_star) need ['//INI_SECTION_NAME//'].(self_gravity) '// &
                               '= .true.')
   endif
   self%has_potential = self%potential == POTENTIAL_LINEAR .or. self%self_gravity
   ! the well-balanced mode subtracts the residual of a static operator built on the constant g (M7-G2b); with the
   ! potential path it is not defined (issue #81, D-81-5)
   if (self%has_potential .and. self%is_balanced) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(equilibrium) = atmosphere is refused with the potential '// &
                               'path (['//INI_SECTION_NAME//'].(potential) = '//POTENTIAL_LINEAR//' or self_gravity, '// &
                               'issue #81)')
   endsubroutine load_from_file
endmodule adam_flume_gravity_object
