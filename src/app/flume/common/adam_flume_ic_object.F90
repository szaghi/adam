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
!<   `s` the cell coordinate along the axis; pressure and velocity unchanged (`psi` zero).
!<
!< The primitive keys of a region follow the physical model: `r, u, v, w, p` (Euler), plus `bx, by, bz` (MHD; `psi` is
!< zero). `isentropic-vortex` is Euler only (issue #41, section 3.7).

! ADAM classes, libraries, parameters
use :: adam_field_object,         only : field_object
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! FLUME modules
use :: adam_flume_euler_library,  only : primitive_to_conservative
use :: adam_flume_parameters,     only : MODEL_EULER, MODEL_MHD, MODEL_MHD_GLM, NV_EULER, strip_control
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
   case(IC_RIEMANN_PROBLEM_STR)
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='regions_number', val=self%regions_number, &
                               error=error)
      if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(regions_number)')
      if (self%regions_number < 1_I4P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(regions_number) must be positive')
   case default
      call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].(type) "'//self%ic_type//'"; expected one of '// &
                               IC_UNIFORM_STR//', '//IC_ISENTROPIC_VORTEX_STR//', '//IC_RIEMANN_PROBLEM_STR//', '// &
                               IC_GLM_PULSE_STR)
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
      call primitive_state_to_conservative(model=self%model, gamma=physics%gamma, prim=prim, q=self%q_region(:,r))
      if (r == 1_I4P) self%prim_1 = prim
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
   real(R8P)                             :: h(2)          !< Seeded perturbations, in [-1, 1).
   real(R8P)                             :: prim(8)       !< Perturbed primitive state of one cell.
   real(R8P)                             :: s_            !< Cell coordinate along the pulse axis.
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
