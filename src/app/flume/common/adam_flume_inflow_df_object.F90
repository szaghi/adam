!< ADAM, FLUME digital-filter turbulent inflow of a boundary layer (issue #119, T3.7).

module adam_flume_inflow_df_object
!< ADAM, FLUME digital-filter turbulent inflow of a boundary layer (issue #119, T3.7).
!<
!< One realm face of kind `inflow-digital-filter` feeds a turbulent boundary layer: the mean profile of a table plus
!< fluctuations with prescribed Reynolds stresses and length scales, built by the digital filter of Klein, Sadiki & Janicka
!< (2003) in the form of Touber & Sandham (2009) and Xie & Castro (2008):
!<
!< 1. **Plane.** A uniform grid of the face, `df_spacing` apart along the wall-normal and the spanwise axes (cell-centred
!<    nodes, the spacing adjusted to divide the domain), periodic along the span when that axis is periodic.
!< 2. **Filter.** Uniform random numbers of zero mean and unit variance (MINSTD, one stream, so every rank builds the same
!<    plane and no communication is needed) are filtered along both plane axes with `b_k = exp(-pi |k| / n)` normalised
!<    to unit sum of squares, `n = L / spacing` the integral length in nodes, the stencil half-width `2 n`: the filtered
!<    field has unit variance and the exponential correlation `exp(-pi r / (2 L))`.
!< 3. **Time.** The streamwise correlation becomes a temporal one, `psi = psi exp(-pi dt / (2 T)) + psi_new sqrt(1 -
!<    exp(-pi dt / T))`, `T = L_x / U_inf` (Xie & Castro 2008): unit variance kept, correlation `exp(-pi t / (2 T))`.
!< 4. **Stresses.** The three unit-variance signals become velocity fluctuations by the Cholesky factor of the Reynolds
!<    stress tensor of the table (Lund, Wu & Squires 1998), with `R_13 = R_23 = 0`: `u' = a11 psi1`, `v' = a21 psi1 +
!<    a22 psi2`, `w' = a33 psi3`.
!< 5. **Thermodynamics.** The strong Reynolds analogy (Morkovin; Touber & Sandham 2009), `T' / T = -(gamma - 1) M^2 u' / U`
!<    with the local mean Mach number, written `T' = -(gamma - 1) U u' rho / (gamma p)` so it holds at U = 0; the
!<    pressure is the free-stream one, `rho' / rho = -T' / T`.
!<
!< The table (`df_profile`, ASCII, `#` comments) has 7 columns: `y / delta`, `U / U_inf`, `T / T_inf`, `R_uu`, `R_vv`,
!< `R_ww`, `R_uv` in units of `U_inf^2`, u streamwise, v wall-normal, w spanwise; the mean density is `rho_inf T_inf / T`
!< (no pressure gradient across the layer). Beyond its last row the free stream, linearly interpolated within.
!<
!< Face keys: `type = inflow-digital-filter`, `r`, `u`, `p` (the free stream: density, speed along the inward normal,
!< pressure), `df_profile`, `df_delta` (the layer thickness), `df_wall_axis` (`x`, `y` or `z`, the wall on its minimum
!< side), `df_spacing`, `df_length_streamwise`, `df_length_normal`, `df_length_spanwise` (integral lengths in units of
!< `df_delta`), `df_seed` (default 1). Only `[physics] physical_model = euler`, one such face per realm.
!<
!< **Ghosts.** Once per step (`advance`) the plane moves on by the step; `build_rows` then lists the ghost cells beyond the
!< face (crown rows beyond that face only) and their states, interpolated bilinearly in the plane; the backends write them
!< after the face pass of the boundary fill. Rows beyond the face and another realm face (edges, corners) are extrapolated.
!< The plane state is saved with the restart (`<restart_basename>.dff`), so a restart continues bit for bit.

use :: adam_flume_parameters, only : MODEL_EULER, strip_control
use :: adam_mpih_global,      only : mpih
use :: adam_parameters,       only : FEC_1_6_ARRAY, FEC_TO_DELTA
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, I8P, R8P, str

implicit none
private
public :: flume_inflow_df_object
public :: load_profile
public :: profile_mean

real(R8P), parameter :: PI=acos(-1._R8P) !< Pi.

type :: flume_inflow_df_object
   !< Digital-filter inflow.
   logical                   :: is_active=.false.        !< A face of this kind exists.
   integer(I4P)              :: face=0_I4P               !< The face (1 to 6).
   integer(I4P)              :: axis_n=1_I4P             !< Streamwise axis (the face normal).
   integer(I4P)              :: axis_w=2_I4P             !< Wall-normal axis.
   integer(I4P)              :: axis_s=3_I4P             !< Spanwise axis.
   real(R8P)                 :: sgn=1._R8P               !< Sign of the inward normal along axis_n.
   real(R8P)                 :: r_inf=1._R8P             !< Free-stream density.
   real(R8P)                 :: u_inf=1._R8P             !< Free-stream speed.
   real(R8P)                 :: p_inf=1._R8P             !< Free-stream pressure.
   real(R8P)                 :: gamma=1.4_R8P            !< Specific heats ratio.
   real(R8P)                 :: delta=1._R8P             !< Layer thickness.
   real(R8P)                 :: spacing=0._R8P           !< Requested plane spacing.
   real(R8P)                 :: length(3)=0._R8P         !< Integral lengths (streamwise, normal, spanwise) / delta.
   integer(I8P)              :: seed=1_I8P               !< Random stream state.
   character(999)            :: profile_file=''          !< Profile table file.
   real(R8P),    allocatable :: prof(:,:)                !< Profile table [7, rows].
   ! plane
   integer(I4P)              :: ny=0_I4P, nz=0_I4P       !< Plane nodes along the wall-normal and spanwise axes.
   real(R8P)                 :: y0=0._R8P, z0=0._R8P     !< Plane origin (the domain minimum along each axis).
   real(R8P)                 :: dy=0._R8P, dz=0._R8P     !< Plane spacing along each axis.
   logical                   :: span_periodic=.false.    !< The spanwise axis is periodic.
   integer(I4P)              :: hy=0_I4P, hz=0_I4P       !< Filter half-widths in nodes.
   real(R8P)                 :: t_scale=0._R8P           !< Time scale T = L_x / U_inf.
   real(R8P),    allocatable :: by(:), bz(:)             !< Filter coefficients [-h:h].
   real(R8P),    allocatable :: psi(:,:,:)               !< Unit-variance signals [3, ny, nz].
   real(R8P),    allocatable :: node_mean(:,:)           !< Per wall-normal node: U, rho, a11, a21, a22, a33 [6, ny].
   real(R8P),    allocatable :: q_plane(:,:,:)           !< Conservative state at the nodes [5, ny, nz].
   ! ghost rows
   integer(I4P)              :: nrows=0_I4P              !< Ghost rows of the face on this rank.
   integer(I4P), allocatable :: rows(:,:)                !< Their block and cell [4, nrows].
   real(R8P),    allocatable :: q_rows(:,:)              !< Their conservative states [5, nrows].
   contains
      procedure, pass(self) :: load            !< Load the face keys.
      procedure, pass(self) :: setup           !< Build the plane, the filter and the first signals.
      procedure, pass(self) :: advance         !< Move the signals on by one step.
      procedure, pass(self) :: build_rows      !< List the ghost rows of the face and their states.
      procedure, pass(self) :: fluctuation     !< Velocity and temperature fluctuations at a plane node.
      procedure, pass(self) :: save_state      !< Save the plane state.
      procedure, pass(self) :: load_state      !< Load the plane state.
endtype flume_inflow_df_object

contains
   subroutine load(self, file_parameters, section, face, model, gamma)
   !< Load the keys of the face `face` (section `section`).
   class(flume_inflow_df_object), intent(inout) :: self            !< Inflow.
   type(file_ini),                intent(in)    :: file_parameters !< Parameters file.
   character(*),                  intent(in)    :: section         !< Face section.
   integer(I4P),                  intent(in)    :: face            !< Face (1 to 6).
   integer(I4P),                  intent(in)    :: model           !< Physical model.
   real(R8P),                     intent(in)    :: gamma           !< Specific heats ratio.
   character(999)                               :: buff            !< Buffer.
   character(20), parameter                     :: LKEY(3)=['df_length_streamwise', 'df_length_normal    ', &
                                                            'df_length_spanwise  ']
   integer(I4P)                                 :: error, d        !< Error status, counter.
   real(R8P)                                    :: seed            !< Seed as read.

   if (self%is_active) call mpih%error_stop(msg=': ['//section//'] a second inflow-digital-filter face: one per realm')
   if (model /= MODEL_EULER) call mpih%error_stop(msg=': ['//section//'] inflow-digital-filter requires [physics].(physical_model) = euler')
   self%is_active = .true.
   self%face = face
   self%axis_n = (face + 1_I4P) / 2_I4P
   self%sgn = merge(1._R8P, -1._R8P, mod(face, 2_I4P) == 1_I4P)
   self%gamma = gamma
   call get_real('r', self%r_inf)
   call get_real('u', self%u_inf)
   call get_real('p', self%p_inf)
   call get_real('df_delta', self%delta)
   call get_real('df_spacing', self%spacing)
   do d=1, 3
      call get_real(trim(LKEY(d)), self%length(d))
   enddo
   call file_parameters%get(section_name=section, option_name='df_wall_axis', val=buff, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//section//'].(df_wall_axis)')
   select case(trim(adjustl(strip_control(buff))))
   case('x')
      self%axis_w = 1_I4P
   case('y')
      self%axis_w = 2_I4P
   case('z')
      self%axis_w = 3_I4P
   case default
      call mpih%error_stop(msg=': ['//section//'].(df_wall_axis) must be x, y or z')
   endselect
   if (self%axis_w == self%axis_n) call mpih%error_stop(msg=': ['//section//'].(df_wall_axis) is the face normal')
   self%axis_s = 6_I4P - self%axis_n - self%axis_w
   call file_parameters%get(section_name=section, option_name='df_seed', val=seed, error=error)
   if (error /= 0) seed = 1._R8P
   self%seed = max(1_I8P, int(seed, I8P))
   call file_parameters%get(section_name=section, option_name='df_profile', val=self%profile_file, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//section//'].(df_profile)')
   self%profile_file = trim(adjustl(strip_control(self%profile_file)))
   call load_profile(filename=trim(self%profile_file), prof=self%prof)
   contains
      subroutine get_real(key, val)
      !< Load a required positive real key.
      character(*), intent(in)  :: key !< Key.
      real(R8P),    intent(out) :: val !< Value.

      call file_parameters%get(section_name=section, option_name=key, val=val, error=error)
      if (error /= 0) call mpih%error_stop(msg=': failed to load ['//section//'].('//key//')')
      if (.not.(val > 0._R8P)) call mpih%error_stop(msg=': ['//section//'].('//key//') must be positive')
      endsubroutine get_real
   endsubroutine load

   subroutine setup(self, emin, emax, periodic)
   !< Build the plane over the domain `[emin, emax]`, the filter coefficients, the mean state of each wall-normal node, and
   !< the first signals (an independent filtered field, so the inflow is turbulent from the first step).
   class(flume_inflow_df_object), intent(inout) :: self        !< Inflow.
   real(R8P),                     intent(in)    :: emin(3)     !< Domain minimum.
   real(R8P),                     intent(in)    :: emax(3)     !< Domain maximum.
   logical,                       intent(in)    :: periodic(3) !< Periodic axes.
   real(R8P)                                    :: n           !< Integral length in nodes.
   real(R8P)                                    :: um, tm, r(4) !< Mean speed and temperature, stresses.
   real(R8P)                                    :: a11, a21     !< Cholesky factors.
   integer(I4P)                                 :: j           !< Counter.

   if (.not.self%is_active) return
   self%y0 = emin(self%axis_w)
   self%z0 = emin(self%axis_s)
   self%ny = max(2_I4P, nint((emax(self%axis_w) - emin(self%axis_w)) / self%spacing, I4P))
   self%nz = max(2_I4P, nint((emax(self%axis_s) - emin(self%axis_s)) / self%spacing, I4P))
   self%dy = (emax(self%axis_w) - emin(self%axis_w)) / real(self%ny, R8P)
   self%dz = (emax(self%axis_s) - emin(self%axis_s)) / real(self%nz, R8P)
   self%span_periodic = periodic(self%axis_s)
   n = self%length(2) * self%delta / self%dy
   call filter_coefficients(n=n, h=self%hy, b=self%by)
   n = self%length(3) * self%delta / self%dz
   call filter_coefficients(n=n, h=self%hz, b=self%bz)
   self%t_scale = self%length(1) * self%delta / self%u_inf
   if (allocated(self%node_mean)) deallocate(self%node_mean)
   allocate(self%node_mean(6, self%ny))
   do j=1, self%ny
      call profile_mean(prof=self%prof, eta=(real(j, R8P) - 0.5_R8P) * self%dy / self%delta, u=um, t=tm, r=r)
      a11 = sqrt(max(r(1), 0._R8P))
      a21 = 0._R8P
      if (a11 > 0._R8P) a21 = r(4) / a11
      self%node_mean(:,j) = [um * self%u_inf, self%r_inf / tm, a11 * self%u_inf, a21 * self%u_inf, &
                             sqrt(max(r(2) - a21**2, 0._R8P)) * self%u_inf, sqrt(max(r(3), 0._R8P)) * self%u_inf]
   enddo
   if (allocated(self%psi)) deallocate(self%psi)
   if (allocated(self%q_plane)) deallocate(self%q_plane)
   allocate(self%psi(3, self%ny, self%nz), self%q_plane(5, self%ny, self%nz))
   call filtered_field(self=self, psi=self%psi)
   call plane_states(self)
   print '(A)', mpih%myrankstr//'inflow-digital-filter: face '//trim(str(self%face, .true.))//', plane '//         &
                trim(str(self%ny, .true.))//' x '//trim(str(self%nz, .true.))//' nodes, filter half-widths '//     &
                trim(str(self%hy, .true.))//' x '//trim(str(self%hz, .true.))//', T = '//trim(str(self%t_scale))
   endsubroutine setup

   subroutine advance(self, dt)
   !< Move the signals on by the step `dt` (Xie & Castro 2008) and refresh the node states.
   class(flume_inflow_df_object), intent(inout) :: self  !< Inflow.
   real(R8P),                     intent(in)    :: dt    !< Time step.
   real(R8P), allocatable                       :: psi_new(:,:,:) !< Fresh filtered field.

   if (.not.self%is_active) return
   allocate(psi_new(3, self%ny, self%nz))
   call filtered_field(self=self, psi=psi_new)
   self%psi = self%psi * exp(-PI * dt / (2._R8P * self%t_scale)) + psi_new * sqrt(1._R8P - exp(-PI * dt / self%t_scale))
   call plane_states(self)
   endsubroutine advance

   subroutine build_rows(self, crown_map, x_cell, y_cell, z_cell)
   !< List the crown rows beyond the face only (block > 0) and their states, bilinear in the plane.
   class(flume_inflow_df_object), intent(inout) :: self              !< Inflow.
   integer(I8P),                  intent(in)    :: crown_map(:,:,:)  !< BC crown map [rows, 9, ngc].
   real(R8P),                     intent(in)    :: x_cell(:,:)       !< X cell centres [1-ngc:ni+ngc, nb] (passed 1-based).
   real(R8P),                     intent(in)    :: y_cell(:,:)       !< Y cell centres.
   real(R8P),                     intent(in)    :: z_cell(:,:)       !< Z cell centres.
   integer(I4P)                                 :: ngc               !< Ghost cells (the offset of the cell arrays).
   integer(I4P)                                 :: crown, c, n, b, fec, idx(3) !< Counters, block, fec, cell.
   real(R8P)                                    :: xyz(3)            !< Cell centre.

   if (.not.self%is_active) return
   ngc = size(crown_map, dim=3)
   n = 0_I4P
   do crown=1, ngc
      do c=1, size(crown_map, dim=1)
         if (row_is_mine(c, crown)) n = n + 1_I4P
      enddo
   enddo
   if (allocated(self%rows)) then
      if (size(self%rows, dim=2) /= n) deallocate(self%rows, self%q_rows)
   endif
   if (.not.allocated(self%rows)) allocate(self%rows(4, n), self%q_rows(5, n))
   self%nrows = n
   n = 0_I4P
   do crown=1, ngc
      do c=1, size(crown_map, dim=1)
         if (.not.row_is_mine(c, crown)) cycle
         n = n + 1_I4P
         b = int(crown_map(c,1,crown), I4P)
         idx = int(crown_map(c,2:4,crown), I4P)
         self%rows(:,n) = [b, idx]
         xyz = [x_cell(idx(1)+ngc,b), y_cell(idx(2)+ngc,b), z_cell(idx(3)+ngc,b)]
         self%q_rows(:,n) = plane_value(self, y=xyz(self%axis_w), z=xyz(self%axis_s))
      enddo
   enddo
   contains
      logical function row_is_mine(c, crown)
      !< True for a live row beyond this face only.
      integer(I4P), intent(in) :: c, crown !< Row and crown.

      row_is_mine = .false.
      if (crown_map(c,1,crown) <= 0_I8P) return
      fec = int(crown_map(c,9,crown), I4P)
      if (count(FEC_TO_DELTA(:,fec) /= 0_I4P) /= 1_I4P) return
      row_is_mine = FEC_1_6_ARRAY(fec) == self%face
      endfunction row_is_mine
   endsubroutine build_rows

   pure function fluctuation(self, j, k) result(f)
   !< Velocity fluctuations (streamwise, wall-normal, spanwise) and temperature fluctuation over the mean temperature at
   !< the node (j, k).
   class(flume_inflow_df_object), intent(in) :: self !< Inflow.
   integer(I4P),                  intent(in) :: j, k !< Node.
   real(R8P)                                 :: f(4) !< u', v', w', T' / T.

   associate(m=>self%node_mean(:,j), s=>self%psi(:,j,k))
   f(1) = m(3) * s(1)
   f(2) = m(4) * s(1) + m(5) * s(2)
   f(3) = m(6) * s(3)
   ! strong Reynolds analogy: T'/T = -(gamma-1) M^2 u'/U = -(gamma-1) U u' rho / (gamma p)
   f(4) = -(self%gamma - 1._R8P) * m(1) * f(1) * m(2) / (self%gamma * self%p_inf)
   endassociate
   endfunction fluctuation

   subroutine save_state(self, basename)
   !< Save the plane state to `<basename>.dff` (rank 0; every rank holds the same state).
   class(flume_inflow_df_object), intent(in) :: self     !< Inflow.
   character(*),                  intent(in) :: basename !< Restart files basename.
   integer(I4P)                              :: unit     !< File unit.

   if (.not.self%is_active) return
   if (mpih%myrank /= 0_I4P) return
   open(newunit=unit, file=trim(adjustl(basename))//'.dff', form='unformatted', access='stream', status='replace', &
        action='write')
   write(unit) self%ny, self%nz, self%seed
   write(unit) self%psi
   close(unit)
   endsubroutine save_state

   subroutine load_state(self, basename)
   !< Load the plane state from `<basename>.dff`; without the file the signals restart from a fresh field (warned).
   class(flume_inflow_df_object), intent(inout) :: self     !< Inflow.
   character(*),                  intent(in)    :: basename !< Restart files basename.
   character(:), allocatable                    :: filename !< File name.
   integer(I4P)                                 :: unit, ios, n(2) !< Unit, status, plane size.
   logical                                      :: is_present      !< File presence.

   if (.not.self%is_active) return
   filename = trim(adjustl(basename))//'.dff'
   inquire(file=filename, exist=is_present)
   if (.not.is_present) then
      call mpih%print_message('warning: '//filename//' is missing: the digital-filter inflow restarts from a fresh field')
      return
   endif
   open(newunit=unit, file=filename, form='unformatted', access='stream', action='read', iostat=ios)
   if (ios == 0) read(unit, iostat=ios) n, self%seed
   if (ios == 0 .and. any(n /= [self%ny, self%nz])) &
      call mpih%error_stop(msg=': '//filename//' holds a plane of '//trim(str(n, .true.))//' nodes, the input '// &
                               trim(str([self%ny, self%nz], .true.)))
   if (ios == 0) read(unit, iostat=ios) self%psi
   if (ios /= 0) call mpih%error_stop(msg=': cannot read '//filename)
   close(unit)
   call plane_states(self)
   endsubroutine load_state

   ! public procedures
   subroutine load_profile(filename, prof)
   !< Load a profile table: 7 columns (y / delta, U / U_inf, T / T_inf, R_uu, R_vv, R_ww, R_uv / U_inf^2), `#` comments,
   !< rows in increasing y; fatal on a malformed table.
   character(*),           intent(in)  :: filename   !< File name.
   real(R8P), allocatable, intent(out) :: prof(:,:)  !< Table [7, rows].
   character(999)                      :: line       !< Line buffer.
   real(R8P)                           :: v(7)       !< Row values.
   integer(I4P)                        :: unit, ios, n, pass !< Unit, status, rows, pass.
   logical                             :: is_present !< File presence.

   inquire(file=filename, exist=is_present)
   if (.not.is_present) call mpih%error_stop(msg=': profile table "'//filename//'" not found')
   do pass=1, 2
      open(newunit=unit, file=filename, action='read', iostat=ios)
      if (ios /= 0) call mpih%error_stop(msg=': cannot open profile table "'//filename//'"')
      n = 0_I4P
      do
         read(unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         read(line, *, iostat=ios) v
         if (ios /= 0) call mpih%error_stop(msg=': profile table "'//filename//'": a row without 7 numbers')
         n = n + 1_I4P
         if (pass == 2) prof(:,n) = v
      enddo
      close(unit)
      if (pass == 1) then
         if (n < 2_I4P) call mpih%error_stop(msg=': profile table "'//filename//'" has fewer than 2 rows')
         allocate(prof(7, n))
      endif
   enddo
   if (any(prof(1,2:) <= prof(1,:n-1))) call mpih%error_stop(msg=': profile table "'//filename//'": y not increasing')
   if (any(prof(3,:) <= 0._R8P)) call mpih%error_stop(msg=': profile table "'//filename//'": T / T_inf must be positive')
   endsubroutine load_profile

   pure subroutine profile_mean(prof, eta, u, t, r)
   !< Return the table's U / U_inf, T / T_inf and stresses (R_uu, R_vv, R_ww, R_uv) at `eta = y / delta`: linear within,
   !< the first row below it, the free stream (1, 1, 0) beyond the last row.
   real(R8P), intent(in)  :: prof(:,:) !< Table [7, rows].
   real(R8P), intent(in)  :: eta       !< y / delta.
   real(R8P), intent(out) :: u, t      !< Mean speed and temperature.
   real(R8P), intent(out) :: r(4)      !< Stresses.
   integer(I4P)           :: lo, hi, mid, n !< Bisection.
   real(R8P)              :: w         !< Weight.

   n = size(prof, dim=2)
   if (eta >= prof(1,n)) then
      u = 1._R8P ; t = 1._R8P ; r = 0._R8P
      return
   endif
   if (eta <= prof(1,1)) then
      u = prof(2,1) ; t = prof(3,1) ; r = prof(4:7,1)
      return
   endif
   lo = 1_I4P ; hi = n
   do while (hi - lo > 1_I4P)
      mid = (lo + hi) / 2_I4P
      if (prof(1,mid) <= eta) then
         lo = mid
      else
         hi = mid
      endif
   enddo
   w = (eta - prof(1,lo)) / (prof(1,hi) - prof(1,lo))
   u = (1._R8P - w) * prof(2,lo) + w * prof(2,hi)
   t = (1._R8P - w) * prof(3,lo) + w * prof(3,hi)
   r = (1._R8P - w) * prof(4:7,lo) + w * prof(4:7,hi)
   endsubroutine profile_mean

   ! private procedures
   subroutine filter_coefficients(n, h, b)
   !< Coefficients `b_k = exp(-pi |k| / n)` normalised to unit sum of squares, half-width `h = ceiling(2 n)` (at least 1).
   real(R8P),              intent(in)  :: n    !< Integral length in nodes.
   integer(I4P),           intent(out) :: h    !< Half-width.
   real(R8P), allocatable, intent(out) :: b(:) !< Coefficients [-h:h].
   integer(I4P)                        :: k    !< Counter.

   h = max(1_I4P, ceiling(2._R8P * n, I4P))
   allocate(b(-h:h))
   do k=-h, h
      b(k) = exp(-PI * real(abs(k), R8P) / max(n, 1.e-12_R8P))
   enddo
   b = b / sqrt(sum(b**2))
   endsubroutine filter_coefficients

   subroutine filtered_field(self, psi)
   !< A fresh filtered field of unit variance: random numbers on the plane padded by the half-widths (wrapped along a
   !< periodic span), filtered along the wall-normal then the spanwise axis.
   type(flume_inflow_df_object), intent(inout) :: self          !< Inflow.
   real(R8P),                    intent(out)   :: psi(:,:,:)    !< Field [3, ny, nz].
   real(R8P), allocatable                      :: r(:,:,:)      !< Random numbers [3, 1-hy:ny+hy, 1-hz:nz+hz].
   real(R8P), allocatable                      :: ry(:,:,:)     !< Filtered along y [3, ny, 1-hz:nz+hz].
   integer(I4P)                                :: c, j, k, a, kk, kz0, kz1 !< Counters, span range drawn.

   kz0 = 1_I4P - self%hz ; kz1 = self%nz + self%hz
   if (self%span_periodic) then
      kz0 = 1_I4P ; kz1 = self%nz
   endif
   allocate(r(3, 1-self%hy:self%ny+self%hy, 1-self%hz:self%nz+self%hz))
   do k=kz0, kz1
      do j=1-self%hy, self%ny+self%hy
         do c=1, 3
            r(c,j,k) = sqrt(12._R8P) * (next_uniform(self%seed) - 0.5_R8P)
         enddo
      enddo
   enddo
   if (self%span_periodic) then
      do k=1-self%hz, 0
         r(:,:,k) = r(:,:,modulo(k-1_I4P, self%nz)+1_I4P)
      enddo
      do k=self%nz+1, self%nz+self%hz
         r(:,:,k) = r(:,:,modulo(k-1_I4P, self%nz)+1_I4P)
      enddo
   endif
   allocate(ry(3, self%ny, 1-self%hz:self%nz+self%hz), source=0._R8P)
   do k=1-self%hz, self%nz+self%hz
      do j=1, self%ny
         do a=-self%hy, self%hy
            ry(:,j,k) = ry(:,j,k) + self%by(a) * r(:,j+a,k)
         enddo
      enddo
   enddo
   psi = 0._R8P
   do k=1, self%nz
      do a=-self%hz, self%hz
         kk = k + a
         psi(:,:,k) = psi(:,:,k) + self%bz(a) * ry(:,:,kk)
      enddo
   enddo
   endsubroutine filtered_field

   subroutine plane_states(self)
   !< The conservative state at every node from the mean and the fluctuations.
   type(flume_inflow_df_object), intent(inout) :: self   !< Inflow.
   real(R8P)                                   :: f(4)   !< Fluctuations.
   real(R8P)                                   :: rho, vel(3), tr !< Density, velocity, T / T_mean.
   integer(I4P)                                :: j, k   !< Counters.

   do k=1, self%nz
      do j=1, self%ny
         f = self%fluctuation(j, k)
         tr = 1._R8P + f(4)
         rho = self%node_mean(2,j) / max(tr, 0.1_R8P)
         vel(self%axis_n) = self%sgn * (self%node_mean(1,j) + f(1))
         vel(self%axis_w) = f(2)
         vel(self%axis_s) = f(3)
         self%q_plane(1,j,k) = rho
         self%q_plane(2:4,j,k) = rho * vel
         self%q_plane(5,j,k) = self%p_inf / (self%gamma - 1._R8P) + 0.5_R8P * rho * sum(vel**2)
      enddo
   enddo
   endsubroutine plane_states

   pure function plane_value(self, y, z) result(q)
   !< The state at (y, z), bilinear between the nodes: clamped along the wall normal, wrapped (periodic span) or clamped
   !< along the span.
   type(flume_inflow_df_object), intent(in) :: self !< Inflow.
   real(R8P),                    intent(in) :: y, z !< Coordinates along the wall-normal and spanwise axes.
   real(R8P)                                :: q(5) !< State.
   real(R8P)                                :: sy, sz, wy, wz !< Node coordinates, weights.
   integer(I4P)                             :: j0, j1, k0, k1 !< Nodes.

   sy = (y - self%y0) / self%dy + 0.5_R8P
   sy = min(max(sy, 1._R8P), real(self%ny, R8P))
   j0 = min(int(sy, I4P), self%ny - 1_I4P) ; j1 = j0 + 1_I4P ; wy = sy - real(j0, R8P)
   sz = (z - self%z0) / self%dz + 0.5_R8P
   if (self%span_periodic) then
      k0 = floor(sz, I4P) ; wz = sz - real(k0, R8P)
      k1 = modulo(k0, self%nz) + 1_I4P ; k0 = modulo(k0 - 1_I4P, self%nz) + 1_I4P
   else
      sz = min(max(sz, 1._R8P), real(self%nz, R8P))
      k0 = min(int(sz, I4P), self%nz - 1_I4P) ; k1 = k0 + 1_I4P ; wz = sz - real(k0, R8P)
   endif
   q = (1._R8P - wy) * (1._R8P - wz) * self%q_plane(:,j0,k0) + wy * (1._R8P - wz) * self%q_plane(:,j1,k0) + &
       (1._R8P - wy) * wz * self%q_plane(:,j0,k1) + wy * wz * self%q_plane(:,j1,k1)
   endfunction plane_value

   function next_uniform(r) result(v)
   !< MINSTD (Park & Miller 1993): the next value in (0, 1), portable, no overflow in 64 bits.
   integer(I8P), intent(inout) :: r !< State.
   real(R8P)                   :: v !< Value.

   r = mod(48271_I8P * r, 2147483647_I8P)
   v = real(r, R8P) / 2147483647._R8P
   endfunction next_uniform
endmodule adam_flume_inflow_df_object
