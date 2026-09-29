!< ADAM, FLUME initial conditions class definition.
module adam_flume_ic_object
!< ADAM, FLUME initial conditions class definition.
!<
!< Accepted `[initial_conditions].(type)`:
!<
!< * `uniform`: the state of region 1 everywhere, density and pressure perturbed by the signed relative amplitude `s`,
!<   `rho (1 + s h1)`, `p (1 + s h2)`, with `h1, h2` in [-1, 1) hashed from the block Morton code and the cell indexes:
!<   the same on every rank decomposition and on every backend (no random generator state);
!< * `isentropic-vortex`: the free stream of region 1 plus an isentropic vortex of axis z centred at `(x0, y0)`, of
!<   radius `radius` and velocity strength `strength` (Shu 1998, ICASE 97-65, section 5.1, written for any free stream):
!<   `du = -strength/(2 pi) (y-y0)/radius exp((1-r^2)/2)`, `dv = strength/(2 pi) (x-x0)/radius exp((1-r^2)/2)`,
!<   `d(p/rho) = -(gamma-1)/gamma strength^2/(8 pi^2) exp(1-r^2)`, isentropic, `r = |x-x0| / radius`; an exact steady
!<   solution convected by the free stream;
!< * `riemann-problem`: piecewise-constant axis-aligned regions, a cell belongs to a region iff `emin < center <= emax`
!<   on every axis. A cell covered by no region is fatal: leaving it at zero density would divide by zero downstream;
!< * `glm-pulse` (MHD only, the GLM d'Alembert test, issue #41, MV-3): the state of region 1 plus a Gaussian pulse in the
!<   magnetic field along `pulse_axis` (x, y or z), `B_axis += pulse_amplitude exp(-((s - pulse_center) / pulse_width)^2)`,
!<   `s` the cell coordinate along the axis; pressure and velocity unchanged (`psi` zero);
!< * `divb-peak` (MHD only, the Dedner et al. 2002 peak in `B_x`, issue #41, MV-10): the state of region 1 plus
!<   `B_x += peak_amplitude (1 - (s / peak_radius)^2)^2` for `s < peak_radius`, `s` the distance in the x-y plane from
!<   `(peak_x0, peak_y0)` (Dedner: `peak_amplitude = 1 / sqrt(4 pi)`, `peak_radius = 1/8`); a non-zero initial div(B);
!< * `mhd-linear-wave` (MHD only, issue #41, MV-5; Stone et al. 2008, section 8.2): the state `q0` of region 1 (global
!<   frame) plus `wave_amplitude R sin(2 pi (x cos a + y sin a) / wavelength)`, `a = wave_angle` (degrees, in the x-y
!<   plane), `R` the right eigenvector of the right-going `wave` (`fast`, `alfven`, `slow` or `entropy`) of the
!<   conservative system at `q0` in the frame of the wave normal (the library Roe-Balsara core, `mhd_eigenvectors`),
!<   rotated back to x, y, z; after one period `wavelength / |u_n + c|` the exact solution is the initial state;
!< * `mhd-cpaw` (MHD only, issue #41, MV-6; Toth 2000, J. Comput. Phys. 161): the circularly polarised Alfven wave, an exact
!<   nonlinear solution of any amplitude. Density and pressure of region 1 (whose velocity and field must be zero), the
!<   field `b_par n + wave_amplitude (sin(phi) t1 + h cos(phi) t2)` and the velocity `-sign(b_par) B_perp / sqrt(rho)`,
!<   `phi = 2 pi (x cos a + y sin a) / wavelength`, `a = wave_angle` (degrees), `n = (cos a, sin a, 0)`,
!<   `t1 = (-sin a, cos a, 0)`, `t2 = z`, `h = +1` (`polarisation = right`) or `-1` (`left`, the mirror image under
!<   `z -> -z`); `|B_perp|` is uniform, so the total pressure is too, and the wave travels along `+n` unchanged at
!<   `|b_par| / sqrt(rho)`: after one period `wavelength sqrt(rho) / |b_par|` the exact solution is the initial state;
!< * `mhd-vortex` (MHD only, issue #41, MV-7; Balsara 2004, ApJS 151): the free stream of region 1 (`bx = by = 0`,
!<   fatal otherwise) plus a magnetised vortex of axis z centred at `(x0, y0)`, of radius `radius`:
!<   `dv = kappa/(2 pi) e (-(y-y0), x-x0) / radius`, `dB = mu/(2 pi) e (-(y-y0), x-x0) / radius`,
!<   `dp = (mu^2 (1-r^2) - rho kappa^2) / (8 pi^2) e^2`, `e = exp((1-r^2)/2)`, `r = |x-x0| / radius`, density uniform:
!<   the radial balance `dp/dr = rho v^2/r - B^2/r - d(B^2/2)/dr` (rationalised units) holds exactly, so the vortex is a
!<   steady solution convected by the free stream, with a divergence-free field;
!< * `orszag-tang` (MHD only, issue #41, MV-12; Stone et al. 2008, section 8.4): the Orszag-Tang vortex on the unit
!<   period, keyless and with no region sections: `rho = 25/(36 pi)`, `p = 5/(12 pi)`, `u = -sin(2 pi y)`,
!<   `v = sin(2 pi x)`, `bx = -sin(2 pi y) / sqrt(4 pi)`, `by = sin(4 pi x) / sqrt(4 pi)`, `w = bz = 0` (divergence-free
!<   pointwise, symmetric under the 180 degrees rotation about `(1/2, 1/2)`). Evaluated as odd functions of
!<   `s = 2 x - 1`, `t = 2 y - 1` (`u = sin(pi t)`, `v = -sin(pi s)`, `by = sin(2 pi s) / sqrt(4 pi)`): on cell centres that
!<   are exact binary fractions the rotated cell has exactly `-s, -t`, so the initial state is bitwise symmetric;
!< * `mhd-rotor` (MHD only, issue #41, MV-13; Balsara and Spicer 1999, Toth 2000 first rotor): the ambient state of
!<   region 1 plus a dense disc centred at `(x0, y0)` in rigid rotation: for `r < r0` density `rho_in` and velocity
!<   `v0 (-(y-y0), x-x0) / r0` (speed `v0` at `r0`); for `r0 <= r < r1` the linear taper `f = (r1 - r) / (r1 - r0)`,
!<   density `rho + (rho_in - rho) f`, velocity `f v0 (-(y-y0), x-x0) / r`; pressure and field of region 1 everywhere;
!<   the velocities add to the ambient one;
!< * `field-loop` (MHD only, issue #41, MV-9; Gardiner and Stone 2005, Mignone and Tzeferacos 2010 section 4.4.1): the
!<   state of region 1 plus the field of the vector potential `A_z = loop_amplitude (loop_radius - r)` for
!<   `r < loop_radius`, `r = |x - (x0, y0)|` in the x-y plane: `dB = loop_amplitude (-(y-y0), x-x0) / r` inside the loop,
!<   zero outside (divergence-free pointwise, discontinuous at `r = loop_radius`); with a uniform `w` the discrete
!<   div(B) feeds `B_z` (`d B_z / dt = w div(B)`);
!< * `rotated-riemann` (issue #41, MV-8; Toth 2000, J. Comput. Phys. 161, section 6.3.2): a one-dimensional Riemann problem
!<   rotated in the x-y plane, periodic along its normal. The normal is `(normal_x, normal_y)`, the coordinate
!<   `s = normal_x x + normal_y y` (unnormalised), and with `f = modulo(s - interface_1, period)`,
!<   `d = interface_2 - interface_1` and `w = interface_2_width` a cell is in region 2 iff `0 < f <= d`, in the linear
!<   ramp of the primitive variables from region 2 to region 1 iff `d < f <= d + w` (at `f = d + w` region 1), in
!<   region 1 otherwise: a jump at `interface_1`, and at `interface_2` a jump (`w = 0`) or a ramp. On a periodic line
!<   the velocity jump of a shock tube must be undone somewhere: a second jump is the mirror problem, an expansion that
!<   may approach vacuum (for Ryu-Jones 1a: density ~2e-4, Alfven speed ~100); the ramp spreads it. The region states
!<   are given in the frame of the normal (`u, v` = normal and tangential velocity, `bx, by` = normal and tangential
!<   field, `w, bz` along z) and rotated to x, y. With integer normal components on a square periodic domain of side
!<   `period`, the state is periodic in x and in y (Toth's shifted-periodic strip is not needed).
!<
!< The primitive keys of a region follow the physical model: `r, u, v, w, p` (Euler), plus `bx, by, bz` (MHD; `psi` is
!< zero). `isentropic-vortex` is Euler only (issue #41, section 3.7).

! ADAM classes, libraries, parameters
use :: adam_field_object,         only : field_object
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_euler_library,  only : primitive_to_conservative
use :: adam_flume_mhd_library,    only : mhd_conservative_to_auxiliary, mhd_eigenvectors, mhd_primitive_to_conservative
use :: adam_flume_parameters,     only : IQ_BX, IQ_BY, IQ_RU, IQ_RV, MODEL_EULER, MODEL_MHD, MODEL_MHD_GLM, NV_AUX_MHD, &
                                        NV_EULER, NV_MHD, &
                                        strip_control
use :: adam_flume_physics_object, only : flume_physics_object, primitive_state_to_conservative
! third party modules
use :: finer,                     only : file_ini
use :: penf,                      only : I4P, I8P, R8P, str

implicit none
private
public :: flume_ic_object

character(len=18), parameter :: INI_SECTION_NAME="initial_conditions"   !< INI section name.
character(len=7),  parameter :: IC_UNIFORM_STR="uniform"                 !< Uniform state.
character(len=17), parameter :: IC_ISENTROPIC_VORTEX_STR="isentropic-vortex" !< Isentropic vortex in a free stream.
character(len=15), parameter :: IC_RIEMANN_PROBLEM_STR="riemann-problem" !< Piecewise-constant regions.
character(len=9),  parameter :: IC_GLM_PULSE_STR="glm-pulse"             !< Gaussian pulse in B along an axis (MHD).
character(len=15), parameter :: PULSE_KEY(3)=['pulse_center   ', &
                                              'pulse_width    ', &
                                              'pulse_amplitude']         !< GLM pulse keys.
character(len=9),  parameter :: IC_DIVB_PEAK_STR="divb-peak"             !< Dedner peak in B_x (MHD).
character(len=14), parameter :: PEAK_KEY(4)=['peak_x0       ', &
                                             'peak_y0       ', &
                                             'peak_radius   ', &
                                             'peak_amplitude']           !< div(B) peak keys.
character(len=15), parameter :: IC_MHD_LINEAR_WAVE_STR="mhd-linear-wave" !< MHD linear wave (MHD).
character(len=14), parameter :: WAVE_KEY(3)=['wave_angle    ', &
                                             'wave_amplitude', &
                                             'wavelength    ']           !< Linear wave and CPAW keys.
character(len=8),  parameter :: IC_MHD_CPAW_STR="mhd-cpaw"               !< Circularly polarised Alfven wave (MHD).
character(len=10), parameter :: IC_MHD_VORTEX_STR="mhd-vortex"           !< Magnetised vortex in a free stream (MHD).
character(len=11), parameter :: IC_ORSZAG_TANG_STR="orszag-tang"         !< Orszag-Tang vortex (MHD).
character(len=9),  parameter :: IC_MHD_ROTOR_STR="mhd-rotor"             !< MHD rotor (MHD).
character(len=10), parameter :: IC_FIELD_LOOP_STR="field-loop"           !< Advected field loop (MHD).
character(len=14), parameter :: LOOP_KEY(4)=['x0            ', 'y0            ', &
                                             'loop_radius   ', 'loop_amplitude'] !< Field loop keys.
character(len=15), parameter :: IC_ROTATED_RIEMANN_STR="rotated-riemann" !< Rotated periodic Riemann problem.
character(len=17), parameter :: ROTATED_KEY(6)=['normal_x         ', 'normal_y         ', &
                                                'interface_1      ', 'interface_2      ', &
                                                'period           ', 'interface_2_width'] !< Rotated Riemann keys.
character(len=6),  parameter :: ROTOR_KEY(6)=['x0    ', 'y0    ', &
                                              'r0    ', 'r1    ', &
                                              'rho_in', 'v0    ']       !< MHD rotor keys.
character(len=6),  parameter :: MHD_VORTEX_KEY(5)=['x0    ', 'y0    ', &
                                                   'radius', 'kappa ', &
                                                   'mu    ']             !< Magnetised vortex keys.
character(len=8),  parameter :: VORTEX_KEY(4)=['x0      ', 'y0      ', &
                                               'radius  ', 'strength']   !< Vortex keys.
real(R8P),         parameter :: PI=acos(-1._R8P)                        !< Pi greek.
character(len=2),  parameter :: PRIM_KEY(8)=['r ', 'u ', 'v ', 'w ', &
                                             'p ', 'bx', 'by', 'bz']    !< Region primitive state keys (MHD: all 8).
character(len=6),  parameter :: EXTENT_KEY(6)=['emin_x', 'emin_y', &
                                               'emin_z', 'emax_x', &
                                               'emax_y', 'emax_z']       !< Region extents keys.

type :: flume_ic_object
   !< FLUME initial conditions class definition.
   integer(I4P)              :: amr_iterations=0_I4P !< AMR iterations performed while imposing the initial conditions.
   character(:), allocatable :: ic_type              !< Initial conditions type.
   integer(I4P)              :: regions_number=0_I4P !< Regions number.
   integer(I4P)              :: model=0_I4P          !< Physical model id.
   integer(I4P)              :: nprim=0_I4P          !< Primitive keys number of a region: 5 (Euler) or 8 (MHD).
   real(R8P),    allocatable :: q_region(:,:)        !< Conservative state of each region [nv, regions_number].
   real(R8P),    allocatable :: emin(:,:)            !< Minimum corner of each region [3, regions_number].
   real(R8P),    allocatable :: emax(:,:)            !< Maximum corner of each region [3, regions_number].
   real(R8P)                 :: gamma=0._R8P         !< Specific heats ratio.
   real(R8P)                 :: prim_1(8)=0._R8P     !< Primitive state of region 1, the free stream (first nprim used).
   real(R8P)                 :: s=0._R8P             !< Uniform: relative amplitude of the seeded perturbation.
   real(R8P)                 :: vortex(4)=0._R8P     !< Isentropic vortex: x0, y0, radius, strength.
   integer(I4P)              :: pulse_axis=0_I4P     !< GLM pulse: axis, 1=x, 2=y, 3=z.
   real(R8P)                 :: pulse(3)=0._R8P      !< GLM pulse: center, width, amplitude.
   real(R8P)                 :: peak(4)=0._R8P       !< div(B) peak: x0, y0, radius, amplitude.
   character(:), allocatable :: wave                 !< Linear wave family: fast, alfven, slow, entropy.
   real(R8P)                 :: wave_par(3)=0._R8P   !< Linear wave, CPAW: angle (degrees), amplitude, wavelength.
   real(R8P)                 :: wave_r(NV_MHD)=0._R8P !< Linear wave: right eigenvector, global frame.
   real(R8P)                 :: b_par=0._R8P         !< CPAW: field along the wave normal.
   real(R8P)                 :: polarisation=0._R8P  !< CPAW: +1 right, -1 left.
   real(R8P)                 :: mvortex(5)=0._R8P    !< Magnetised vortex: x0, y0, radius, kappa, mu.
   real(R8P)                 :: rotor(6)=0._R8P      !< MHD rotor: x0, y0, r0, r1, rho_in, v0.
   real(R8P)                 :: loop(4)=0._R8P       !< Field loop: x0, y0, radius, amplitude.
   real(R8P)                 :: rotated(6)=0._R8P    !< Rotated Riemann: normal x, y, interface 1, 2, period, width 2.
   real(R8P)                 :: prim_2(8)=0._R8P     !< Primitive state of region 2 (global frame, first nprim used).
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
   if (self%ic_type == IC_UNIFORM_STR) &
   desc = desc//NL//mpih%myrankstr//'  s:              '//trim(str(self%s))
   if (self%ic_type == IC_ISENTROPIC_VORTEX_STR) &
   desc = desc//NL//mpih%myrankstr//'  vortex:         '//trim(str(self%vortex))
   if (self%ic_type == IC_GLM_PULSE_STR) &
   desc = desc//NL//mpih%myrankstr//'  pulse:          axis '//trim(str(self%pulse_axis))//', '//trim(str(self%pulse))
   if (self%ic_type == IC_DIVB_PEAK_STR) &
   desc = desc//NL//mpih%myrankstr//'  peak:           '//trim(str(self%peak))
   if (self%ic_type == IC_MHD_LINEAR_WAVE_STR) &
   desc = desc//NL//mpih%myrankstr//'  wave:           '//self%wave//', '//trim(str(self%wave_par))
   if (self%ic_type == IC_MHD_CPAW_STR) &
   desc = desc//NL//mpih%myrankstr//'  cpaw:           b_par '//trim(str(self%b_par))//', polarisation '// &
                                    trim(str(self%polarisation))//', '//trim(str(self%wave_par))
   if (self%ic_type == IC_MHD_VORTEX_STR) &
   desc = desc//NL//mpih%myrankstr//'  mhd vortex:     '//trim(str(self%mvortex))
   if (self%ic_type == IC_MHD_ROTOR_STR) &
   desc = desc//NL//mpih%myrankstr//'  rotor:          '//trim(str(self%rotor))
   if (self%ic_type == IC_FIELD_LOOP_STR) &
   desc = desc//NL//mpih%myrankstr//'  field loop:     '//trim(str(self%loop))
   if (self%ic_type == IC_ROTATED_RIEMANN_STR) &
   desc = desc//NL//mpih%myrankstr//'  rotated:        '//trim(str(self%rotated))
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
   real(R8P)                                 :: prim(8)         !< Region primitive state (first nprim used).
   real(R8P)                                 :: extent(6)       !< Region extents.
   integer(I4P)                              :: error           !< Error status.
   integer(I4P)                              :: r, k            !< Counters.

   self%gamma = physics%gamma
   self%model = physics%model
   select case(self%model)
   case(MODEL_EULER)
      self%nprim = 5_I4P
   case(MODEL_MHD, MODEL_MHD_GLM)
      self%nprim = 8_I4P
   case default
      call mpih%error_stop(msg=': no initial conditions for physical model "'//physics%physical_model//'"')
   endselect
   prim = 0._R8P

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
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='s', val=self%s, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(s)')
   case(IC_ISENTROPIC_VORTEX_STR)
      if (self%model /= MODEL_EULER) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_ISENTROPIC_VORTEX_STR//' requires '// &
                                  '[physics].(physical_model) = euler')
      self%regions_number = 1_I4P
      do k=1, 4
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(VORTEX_KEY(k)), val=self%vortex(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(VORTEX_KEY(k))//')')
      enddo
      if (self%vortex(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(radius) must be positive')
   case(IC_GLM_PULSE_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_GLM_PULSE_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='pulse_axis', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(pulse_axis)')
      select case(trim(adjustl(strip_control(buff))))
      case('x')
         self%pulse_axis = 1_I4P
      case('y')
         self%pulse_axis = 2_I4P
      case('z')
         self%pulse_axis = 3_I4P
      case default
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(pulse_axis) "'//trim(adjustl(buff))// &
                                  '"; expected one of x, y, z')
      endselect
      do k=1, 3
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(PULSE_KEY(k)), val=self%pulse(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(PULSE_KEY(k))//')')
      enddo
      if (self%pulse(2) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(pulse_width) must be positive')
   case(IC_DIVB_PEAK_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_DIVB_PEAK_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      do k=1, 4
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(PEAK_KEY(k)), val=self%peak(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(PEAK_KEY(k))//')')
      enddo
      if (self%peak(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(peak_radius) must be positive')
   case(IC_MHD_LINEAR_WAVE_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_LINEAR_WAVE_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='wave', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(wave)')
      self%wave = trim(adjustl(strip_control(buff)))
      select case(self%wave)
      case('fast', 'alfven', 'slow', 'entropy')
      case default
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(wave) "'//self%wave// &
                                  '"; expected one of fast, alfven, slow, entropy')
      endselect
      do k=1, 3
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(WAVE_KEY(k)), val=self%wave_par(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(WAVE_KEY(k))//')')
      enddo
      if (self%wave_par(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(wavelength) must be positive')
   case(IC_MHD_CPAW_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_CPAW_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='polarisation', val=buff, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(polarisation)')
      select case(trim(adjustl(strip_control(buff))))
      case('right')
         self%polarisation =  1._R8P
      case('left')
         self%polarisation = -1._R8P
      case default
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(polarisation) "'//trim(adjustl(buff))// &
                                  '"; expected one of right, left')
      endselect
      do k=1, 3
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(WAVE_KEY(k)), val=self%wave_par(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(WAVE_KEY(k))//')')
      enddo
      if (self%wave_par(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(wavelength) must be positive')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='b_par', val=self%b_par, error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(b_par)')
      if (self%b_par == 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(b_par) must be non-zero')
   case(IC_MHD_VORTEX_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_VORTEX_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      do k=1, 5
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(MHD_VORTEX_KEY(k)), &
                                  val=self%mvortex(k), error=error)
         if (error > 0) &
            call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(MHD_VORTEX_KEY(k))//')')
      enddo
      if (self%mvortex(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(radius) must be positive')
   case(IC_ORSZAG_TANG_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_ORSZAG_TANG_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 0_I4P
   case(IC_MHD_ROTOR_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_ROTOR_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      do k=1, 6
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(ROTOR_KEY(k)), val=self%rotor(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(ROTOR_KEY(k))//')')
      enddo
      if (self%rotor(3) <= 0._R8P .or. self%rotor(4) <= self%rotor(3)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs 0 < r0 < r1')
      if (self%rotor(5) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(rho_in) must be positive')
   case(IC_FIELD_LOOP_STR)
      if (self%model /= MODEL_MHD .and. self%model /= MODEL_MHD_GLM) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_FIELD_LOOP_STR//' requires '// &
                                  '[physics].(physical_model) = mhd-ideal')
      self%regions_number = 1_I4P
      do k=1, 4
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(LOOP_KEY(k)), val=self%loop(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(LOOP_KEY(k))//')')
      enddo
      if (self%loop(3) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(loop_radius) must be positive')
   case(IC_ROTATED_RIEMANN_STR)
      self%regions_number = 2_I4P
      do k=1, 6
         call file_parameters%get(section_name=INI_SECTION_NAME, option_name=trim(ROTATED_KEY(k)), val=self%rotated(k), &
                                  error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//trim(ROTATED_KEY(k))//')')
      enddo
      if (self%rotated(1) == 0._R8P .and. self%rotated(2) == 0._R8P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(normal_x, normal_y) must not both be 0')
      if (self%rotated(5) <= 0._R8P) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(period) must be positive')
      if (self%rotated(4) <= self%rotated(3) .or. self%rotated(6) < 0._R8P .or. &
          self%rotated(4) + self%rotated(6) >= self%rotated(3) + self%rotated(5)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs interface_1 < interface_2, interface_2_width >= 0 '// &
                                  'and interface_2 + interface_2_width < interface_1 + period')
   case(IC_RIEMANN_PROBLEM_STR)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, &
                               error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')
      if (self%regions_number < 1_I4P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(regions_number) must be positive')
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(type) "'//self%ic_type//'"; expected one of '// &
                               IC_UNIFORM_STR//', '//IC_ISENTROPIC_VORTEX_STR//', '//IC_RIEMANN_PROBLEM_STR//', '// &
                               IC_GLM_PULSE_STR//', '//IC_DIVB_PEAK_STR//', '//IC_MHD_LINEAR_WAVE_STR//', '// &
                               IC_MHD_CPAW_STR//', '//IC_MHD_VORTEX_STR//', '//IC_ORSZAG_TANG_STR//', '// &
                               IC_MHD_ROTOR_STR//', '//IC_FIELD_LOOP_STR//', '//IC_ROTATED_RIEMANN_STR)
   endselect

   if (allocated(self%q_region)) deallocate(self%q_region)
   if (allocated(self%emin)) deallocate(self%emin)
   if (allocated(self%emax)) deallocate(self%emax)
   allocate(self%q_region(physics%nv,self%regions_number), self%emin(3,self%regions_number), &
            self%emax(3,self%regions_number))
   self%emin = -huge(1._R8P)
   self%emax =  huge(1._R8P)
   do r=1, self%regions_number
      sname = INI_SECTION_NAME//'_region_'//trim(str(r, .true.))
      do k=1, self%nprim
         call file_parameters%get(section_name=sname, option_name=trim(PRIM_KEY(k)), val=prim(k), error=error)
         if (error > 0) call mpih%error_stop(msg=': failed to load ['//sname//'].('//trim(PRIM_KEY(k))//')')
      enddo
      if (self%ic_type == IC_ROTATED_RIEMANN_STR) call rotate_to_global(normal=self%rotated(1:2), nprim=self%nprim, prim=prim)
      call primitive_state_to_conservative(model=self%model, gamma=physics%gamma, prim=prim, q=self%q_region(:,r))
      if (r == 1_I4P) self%prim_1 = prim
      if (r == 2_I4P) self%prim_2 = prim
      if (r == 1_I4P .and. self%ic_type == IC_MHD_LINEAR_WAVE_STR) &
         call linear_wave_eigenvector(gamma=physics%gamma, R=physics%R, prim=prim, wave=self%wave, &
                                      angle=self%wave_par(1), r_global=self%wave_r)
      if (self%ic_type == IC_RIEMANN_PROBLEM_STR) then
         do k=1, 6
            call file_parameters%get(section_name=sname, option_name=EXTENT_KEY(k), val=extent(k), error=error)
            if (error > 0) call mpih%error_stop(msg=': failed to load ['//sname//'].('//EXTENT_KEY(k)//')')
         enddo
         self%emin(:,r) = extent(1:3)
         self%emax(:,r) = extent(4:6)
      endif
   enddo
   if (self%ic_type == IC_MHD_CPAW_STR) then
      if (any(self%prim_1(2:4) /= 0._R8P) .or. any(self%prim_1(6:8) /= 0._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_CPAW_STR//' builds the velocity and '// &
                                  'the field from b_par and wave_amplitude: ['//INI_SECTION_NAME//'_region_1].(u, v, w, '// &
                                  'bx, by, bz) must be 0')
   endif
   if (self%ic_type == IC_MHD_VORTEX_STR) then
      if (any(self%prim_1(6:7) /= 0._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(type) = '//IC_MHD_VORTEX_STR//' is an equilibrium only '// &
                                  'without an in-plane free-stream field: ['//INI_SECTION_NAME//'_region_1].(bx, by) must '// &
                                  'be 0')
   endif
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
   real(R8P)                             :: h(2)          !< Seeded perturbations, in [-1, 1).
   real(R8P)                             :: prim(8)       !< Perturbed primitive state of one cell.
   real(R8P)                             :: s_            !< Pulse axis coordinate, scaled peak distance, or phase.
   real(R8P)                             :: ca, sa        !< cos, sin of the CPAW normal angle.
   real(R8P)                             :: bt(2)         !< CPAW transverse field, t1 and t2 components.
   logical                               :: is_set        !< Flag: cell covered by a region.
   integer(I4P)                          :: b, i, j, k, r !< Counters.

   select case(self%ic_type)
   case(IC_UNIFORM_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  call hash_cell(code=field%code(b), i=i, j=j, k=k, h=h)
                  prim    = self%prim_1
                  prim(1) = self%prim_1(1) * (1._R8P + self%s * h(1))
                  prim(5) = self%prim_1(5) * (1._R8P + self%s * h(2))
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_ISENTROPIC_VORTEX_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  call isentropic_vortex(gamma=self%gamma, prim=self%prim_1, vortex=self%vortex, x=field%x_cell(i,b), &
                                         y=field%y_cell(j,b), q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_GLM_PULSE_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  center = [field%x_cell(i,b), field%y_cell(j,b), field%z_cell(k,b)]
                  s_     = center(self%pulse_axis)
                  prim   = self%prim_1
                  prim(5+self%pulse_axis) = prim(5+self%pulse_axis) + &
                                            self%pulse(3) * exp(-((s_ - self%pulse(1)) / self%pulse(2))**2)
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_MHD_LINEAR_WAVE_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  s_ = sin(2._R8P * PI * (field%x_cell(i,b) * cos(self%wave_par(1) * PI / 180._R8P) + &
                                          field%y_cell(j,b) * sin(self%wave_par(1) * PI / 180._R8P)) / self%wave_par(3))
                  q(:,i,j,k,b) = self%q_region(:,1)
                  q(1:NV_MHD,i,j,k,b) = q(1:NV_MHD,i,j,k,b) + self%wave_par(2) * s_ * self%wave_r
               enddo
            enddo
         enddo
      enddo
   case(IC_MHD_CPAW_STR)
      ca = cos(self%wave_par(1) * PI / 180._R8P)
      sa = sin(self%wave_par(1) * PI / 180._R8P)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  s_   = 2._R8P * PI * (field%x_cell(i,b) * ca + field%y_cell(j,b) * sa) / self%wave_par(3)
                  bt   = self%wave_par(2) * [sin(s_), self%polarisation * cos(s_)]
                  prim = self%prim_1
                  prim(2:4) = -sign(1._R8P, self%b_par) / sqrt(self%prim_1(1)) * [-sa * bt(1), ca * bt(1), bt(2)]
                  prim(6:8) = [ca * self%b_par - sa * bt(1), sa * self%b_par + ca * bt(1), bt(2)]
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_MHD_VORTEX_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  call mhd_vortex(prim0=self%prim_1, mvortex=self%mvortex, x=field%x_cell(i,b), y=field%y_cell(j,b), &
                                  prim=prim)
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_ORSZAG_TANG_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  prim(1) = 25._R8P / (36._R8P * PI)
                  prim(2) =  sin(PI * (2._R8P * field%y_cell(j,b) - 1._R8P))
                  prim(3) = -sin(PI * (2._R8P * field%x_cell(i,b) - 1._R8P))
                  prim(4) = 0._R8P
                  prim(5) = 5._R8P / (12._R8P * PI)
                  prim(6) = prim(2) / sqrt(4._R8P * PI)
                  prim(7) = sin(2._R8P * PI * (2._R8P * field%x_cell(i,b) - 1._R8P)) / sqrt(4._R8P * PI)
                  prim(8) = 0._R8P
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_MHD_ROTOR_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  call mhd_rotor(prim0=self%prim_1, rotor=self%rotor, x=field%x_cell(i,b), y=field%y_cell(j,b), &
                                 prim=prim)
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_FIELD_LOOP_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  prim = self%prim_1
                  center(1:2) = [field%x_cell(i,b) - self%loop(1), field%y_cell(j,b) - self%loop(2)]
                  s_ = sqrt(center(1) * center(1) + center(2) * center(2))
                  if (s_ > 0._R8P .and. s_ < self%loop(3)) then
                     prim(6) = prim(6) - self%loop(4) * center(2) / s_
                     prim(7) = prim(7) + self%loop(4) * center(1) / s_
                  endif
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_DIVB_PEAK_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  s_   = sqrt((field%x_cell(i,b) - self%peak(1))**2 + (field%y_cell(j,b) - self%peak(2))**2) / self%peak(3)
                  prim = self%prim_1
                  if (s_ < 1._R8P) prim(6) = prim(6) + self%peak(4) * (1._R8P - s_ * s_)**2
                  call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
               enddo
            enddo
         enddo
      enddo
   case(IC_ROTATED_RIEMANN_STR)
      do b=1, field%blocks_number
         do k=1, field%nk
            do j=1, field%nj
               do i=1, field%ni
                  s_ = modulo(self%rotated(1) * field%x_cell(i,b) + self%rotated(2) * field%y_cell(j,b) - self%rotated(3), &
                              self%rotated(5)) - (self%rotated(4) - self%rotated(3))
                  if (s_ > -(self%rotated(4) - self%rotated(3)) .and. s_ <= 0._R8P) then
                     q(:,i,j,k,b) = self%q_region(:,2)
                  elseif (s_ > 0._R8P .and. s_ < self%rotated(6)) then
                     prim = self%prim_2 + (self%prim_1 - self%prim_2) * (s_ / self%rotated(6))
                     call primitive_state_to_conservative(model=self%model, gamma=self%gamma, prim=prim, q=q(:,i,j,k,b))
                  else
                     q(:,i,j,k,b) = self%q_region(:,1)
                  endif
               enddo
            enddo
         enddo
      enddo
   case(IC_RIEMANN_PROBLEM_STR)
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
                     call mpih%error_stop(msg=': cell center ('//trim(str(center(1)))//', '//trim(str(center(2)))// &
                                              ', '//trim(str(center(3)))//') is covered by no ['//INI_SECTION_NAME// &
                                              '_region_*]')
               enddo
            enddo
         enddo
      enddo
   endselect
   endsubroutine set_initial_conditions

   ! private procedures
   pure subroutine hash_cell(code, i, j, k, h)
   !< Return two pseudo-random numbers in [-1, 1) that depend only on the block Morton code and the cell indexes.
   !<
   !< Integer mixing with no integer overflow, so the same bits on every compiler; the cell key never depends on the rank
   !< that owns the block. The key is chained (mix the code, then fold in i, j, k one at a time and mix again), so keys
   !< whose bits overlap, such as small codes and small indexes, do not collide. Xorshift alone is linear over GF(2), which
   !< would make the numbers an affine function of the key bits (exactly equidistributed, spatially structured): each mix
   !< also multiplies by an odd constant modulo 2^52 (xorshift* idea, Vigna 2016), in 26-bit limbs to stay in range.
   integer(I8P), intent(in)  :: code             !< Block Morton code.
   integer(I4P), intent(in)  :: i, j, k          !< Cell indexes.
   real(R8P),    intent(out) :: h(2)             !< Pseudo-random numbers in [-1, 1).
   integer(I8P), parameter   :: MASK=2_I8P**52-1 !< Mask of the 52 low bits.
   integer(I8P)              :: x                !< Hash state.
   integer(I4P)              :: n                !< Counter.

   x = mix(code)
   x = mix(ieor(x, int(i, I8P)))
   x = mix(ieor(x, int(j, I8P)))
   x = mix(ieor(x, int(k, I8P)))
   do n=1, 2
      x = mix(x)
      h(n) = 2._R8P * real(iand(x, MASK), R8P) / real(MASK + 1_I8P, R8P) - 1._R8P
   enddo
   contains
      pure function mix(x_in) result(x_out)
      !< Three xorshift64 steps (13, 7, 17), a multiplication by an odd constant modulo 2^52, a final xorshift of the
      !< high bits onto the low ones; the zero state is mapped to a non-zero one.
      integer(I8P), intent(in) :: x_in                 !< State.
      integer(I8P)             :: x_out                !< Mixed state.
      integer(I8P), parameter  :: M26=2_I8P**26-1      !< Mask of the 26 low bits.
      integer(I8P), parameter  :: C=39083855_I8P       !< Odd multiplier, < 2^26 (so each limb product is < 2^52).
      integer(I4P)             :: m                    !< Counter.

      x_out = x_in
      if (x_out == 0_I8P) x_out = 88172645463325252_I8P
      do m=1, 3
         x_out = ieor(x_out, ishft(x_out, 13))
         x_out = ieor(x_out, ishft(x_out, -7))
         x_out = ieor(x_out, ishft(x_out, 17))
      enddo
      x_out = iand(x_out, MASK)
      x_out = iand(iand(x_out, M26) * C + ishft(iand(ishft(x_out, -26) * C, M26), 26), MASK)
      x_out = ieor(x_out, ishft(x_out, -29))
      endfunction mix
   endsubroutine hash_cell

   pure subroutine linear_wave_eigenvector(gamma, R, prim, wave, angle, r_global)
   !< Return the right eigenvector of the right-going `wave` of the conservative MHD system at the primitive state
   !< `prim` (global frame), for the wave normal at `angle` degrees in the x-y plane: the state is rotated into the frame
   !< of the normal, `R` is the column of `mhd_eigenvectors` (x direction of that frame; waves ordered
   !< `u_n - c_f, u_n - c_a, u_n - c_s, u_n, u_n + c_s, u_n + c_a, u_n + c_f, B_n`), and its momentum and field are
   !< rotated back to x, y, z.
   real(R8P),    intent(in)  :: gamma            !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                !< Gas constant.
   real(R8P),    intent(in)  :: prim(8)          !< Primitive state, global frame (r, u, v, w, p, bx, by, bz).
   character(*), intent(in)  :: wave             !< Wave family: fast, alfven, slow, entropy.
   real(R8P),    intent(in)  :: angle            !< Wave normal angle in the x-y plane (degrees).
   real(R8P),    intent(out) :: r_global(NV_MHD) !< Right eigenvector, global frame.
   real(R8P)                 :: c, s             !< cos, sin of the angle.
   real(R8P)                 :: q(NV_MHD)        !< Conservative state, wave frame.
   real(R8P)                 :: qa(NV_AUX_MHD)   !< Auxiliary state, wave frame.
   real(R8P)                 :: el(NV_MHD,NV_MHD) !< Left eigenvectors.
   real(R8P)                 :: er(NV_MHD,NV_MHD) !< Right eigenvectors.
   integer(I4P)              :: k                !< Wave index.

   c = cos(angle * PI / 180._R8P)
   s = sin(angle * PI / 180._R8P)
   call mhd_primitive_to_conservative(gamma=gamma, r=prim(1), u=c * prim(2) + s * prim(3), v=-s * prim(2) + c * prim(3), &
                                      w=prim(4), p=prim(5), bx=c * prim(6) + s * prim(7), by=-s * prim(6) + c * prim(7),  &
                                      bz=prim(8), q=q)
   call mhd_conservative_to_auxiliary(gamma=gamma, R=R, q=q, qa=qa)
   call mhd_eigenvectors(gamma=gamma, d=1_I4P, qa=qa, el=el, er=er)
   select case(wave)
   case('fast')
      k = 7_I4P
   case('alfven')
      k = 6_I4P
   case('slow')
      k = 5_I4P
   case default ! entropy
      k = 4_I4P
   endselect
   r_global    = er(:,k)
   r_global(IQ_RU) = c * er(IQ_RU,k) - s * er(IQ_RV,k)
   r_global(IQ_RV) = s * er(IQ_RU,k) + c * er(IQ_RV,k)
   r_global(IQ_BX) = c * er(IQ_BX,k) - s * er(IQ_BY,k)
   r_global(IQ_BY) = s * er(IQ_BX,k) + c * er(IQ_BY,k)
   endsubroutine linear_wave_eigenvector

   pure subroutine mhd_vortex(prim0, mvortex, x, y, prim)
   !< Return the primitive state of the magnetised vortex at `(x, y)` (formulas in the module documentation).
   real(R8P), intent(in)  :: prim0(8)   !< Free stream primitive state (r, u, v, w, p, bx, by, bz), bx = by = 0.
   real(R8P), intent(in)  :: mvortex(5) !< Vortex x0, y0, radius, kappa, mu.
   real(R8P), intent(in)  :: x, y       !< Point coordinates.
   real(R8P), intent(out) :: prim(8)    !< Primitive state.
   real(R8P)              :: dx, dy     !< Scaled distances from the vortex centre.
   real(R8P)              :: rr         !< Scaled squared distance.
   real(R8P)              :: e          !< Gaussian factor, exp((1 - r^2) / 2).

   dx = (x - mvortex(1)) / mvortex(3)
   dy = (y - mvortex(2)) / mvortex(3)
   rr = dx * dx + dy * dy
   e  = exp(0.5_R8P * (1._R8P - rr))
   prim    = prim0
   prim(2) = prim0(2) - mvortex(4) / (2._R8P * PI) * dy * e
   prim(3) = prim0(3) + mvortex(4) / (2._R8P * PI) * dx * e
   prim(5) = prim0(5) + (mvortex(5)**2 * (1._R8P - rr) - prim0(1) * mvortex(4)**2) / (8._R8P * PI**2) * e * e
   prim(6) = -mvortex(5) / (2._R8P * PI) * dy * e
   prim(7) =  mvortex(5) / (2._R8P * PI) * dx * e
   endsubroutine mhd_vortex

   pure subroutine mhd_rotor(prim0, rotor, x, y, prim)
   !< Return the primitive state of the MHD rotor at `(x, y)` (formulas in the module documentation).
   real(R8P), intent(in)  :: prim0(8)  !< Ambient primitive state (r, u, v, w, p, bx, by, bz).
   real(R8P), intent(in)  :: rotor(6)  !< Rotor x0, y0, r0, r1, rho_in, v0.
   real(R8P), intent(in)  :: x, y      !< Point coordinates.
   real(R8P), intent(out) :: prim(8)   !< Primitive state.
   real(R8P)              :: dx, dy, r !< Distances from the rotor centre.
   real(R8P)              :: f         !< Taper factor.

   dx = x - rotor(1)
   dy = y - rotor(2)
   r  = sqrt(dx * dx + dy * dy)
   prim = prim0
   if (r < rotor(3)) then
      prim(1) = rotor(5)
      prim(2) = prim0(2) - rotor(6) * dy / rotor(3)
      prim(3) = prim0(3) + rotor(6) * dx / rotor(3)
   elseif (r < rotor(4)) then
      f = (rotor(4) - r) / (rotor(4) - rotor(3))
      prim(1) = prim0(1) + (rotor(5) - prim0(1)) * f
      prim(2) = prim0(2) - f * rotor(6) * dy / r
      prim(3) = prim0(3) + f * rotor(6) * dx / r
   endif
   endsubroutine mhd_rotor

   pure subroutine rotate_to_global(normal, nprim, prim)
   !< Rotate the velocity and (MHD) the field of a primitive state from the frame of `normal` (normal, tangential, z) to
   !< x, y, z: `v = v_n n + v_t t`, `n = normal / |normal|`, `t = (-n_y, n_x)`.
   real(R8P),    intent(in)    :: normal(2) !< Normal, (x, y), not normalised.
   integer(I4P), intent(in)    :: nprim     !< Primitive keys number: 5 (Euler) or 8 (MHD).
   real(R8P),    intent(inout) :: prim(8)   !< Primitive state (r, u, v, w, p, bx, by, bz).
   real(R8P)                   :: n(2)      !< Unit normal.
   real(R8P)                   :: a, b      !< Normal and tangential components.

   n = normal / sqrt(normal(1) * normal(1) + normal(2) * normal(2))
   a = prim(2)
   b = prim(3)
   prim(2) = a * n(1) - b * n(2)
   prim(3) = a * n(2) + b * n(1)
   if (nprim == 8_I4P) then
      a = prim(6)
      b = prim(7)
      prim(6) = a * n(1) - b * n(2)
      prim(7) = a * n(2) + b * n(1)
   endif
   endsubroutine rotate_to_global

   pure subroutine isentropic_vortex(gamma, prim, vortex, x, y, q)
   !< Return the conservative state of the isentropic vortex at `(x, y)` (formulas in the module documentation).
   real(R8P), intent(in)  :: gamma       !< Specific heats ratio.
   real(R8P), intent(in)  :: prim(5)     !< Free stream primitive state (r, u, v, w, p).
   real(R8P), intent(in)  :: vortex(4)   !< Vortex x0, y0, radius, strength.
   real(R8P), intent(in)  :: x, y        !< Point coordinates.
   real(R8P), intent(out) :: q(NV_EULER) !< Conservative variables.
   real(R8P)              :: dx, dy      !< Scaled distances from the vortex centre.
   real(R8P)              :: e           !< Gaussian factor, exp((1 - r^2) / 2).
   real(R8P)              :: T, T0       !< p / rho, perturbed and free stream.
   real(R8P)              :: rho         !< Density.

   dx = (x - vortex(1)) / vortex(3)
   dy = (y - vortex(2)) / vortex(3)
   e  = exp(0.5_R8P * (1._R8P - dx * dx - dy * dy))
   T0 = prim(5) / prim(1)
   T  = T0 - (gamma - 1._R8P) / gamma * vortex(4)**2 / (8._R8P * PI**2) * e * e
   rho = prim(1) * (T / T0)**(1._R8P / (gamma - 1._R8P))
   call primitive_to_conservative(gamma=gamma, r=rho, u=prim(2) - vortex(4) / (2._R8P * PI) * dy * e, &
                                  v=prim(3) + vortex(4) / (2._R8P * PI) * dx * e, w=prim(4), p=rho * T, q=q)
   endsubroutine isentropic_vortex
endmodule adam_flume_ic_object
