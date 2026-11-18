!< ADAM, FLUME recycling-rescaling turbulent inflow of a boundary layer (issue #121, S5).

module adam_flume_inflow_rr_object
!< ADAM, FLUME recycling-rescaling turbulent inflow of a boundary layer (issue #121, S5).
!<
!< One realm face of kind `inflow-recycling` feeds a turbulent boundary layer with the mean profile of a table (the
!< digital filter's, `adam_flume_inflow_df_object`) plus the fluctuations of the layer itself, taken from a plane downstream
!< and rescaled, in the variant of the STREAmS solver (Pirozzoli, Bernardini & Grasso; STREAmS-2, `recyc_prepare` and
!< `bcrecyc`) of the method of Lund, Wu & Squires (1998):
!<
!< 1. **Plane.** The `ngc` cell layers ending at the station `rr_x` downstream of the face, a pure translation of the ghost
!<    layers by `L = nint(rr_x / dx) dx` (ghost layer g, g cells beyond the face, takes the cell L downstream of it): every rank
!<    writes the cells it owns into one global plane [6, ngc, ny, nz] (the 6th slot counts the writes) and an allreduce sums it,
!<    so the plane is the same on every rank and independent of the decomposition (every entry is written once).
!< 2. **Fluctuations.** Per layer and wall-normal row, the spanwise Favre mean `u~ = sum(rho u) / sum(rho)` and the mean
!<    density; the fluctuations `u'' = u - u~` (three components) and `rho' = rho - rho_mean`.
!< 3. **Rescaling.** A ghost at height y above the wall takes the fluctuations at `beta y` (inner, the viscous length grows by
!<    `beta = delta_v(x_r) / delta_v(0)`) and at `delta_ratio y` (outer, `delta_ratio = delta(x_r) / delta(0)`), linear in the
!<    plane, blended by `W(eta) = (1 + tanh(4 (eta - 0.2) / (0.2 + 0.6 eta)) / tanh(4)) / 2`, `eta = y / delta`, damped by
!<    `(1 - tanh(4 (eta - 2))) / 2`; the velocity fluctuations are scaled by `beta` (the friction velocity ratio), the density
!<    one clipped to [-0.75, 4] times the mean. A source above the plane's last row centre gives no fluctuation; below its first row
!<    the first row's (STREAmS zeroes the first cell, a DNS cell at y+ < 1; here the first cell is the wall model's log layer).
!< 4. **Span shift.** The plane is read `rr_shift` of the span away (default half), which breaks the streamwise locking of the
!<    recycled structures; a shift needs a periodic span.
!< 5. **State.** Mean (the table) plus fluctuations, the free-stream pressure: `rho = rho_mean + rho'`, `u = U + u''`, and the
!<    energy of the pressure `p_inf`.
!<
!< `beta` and `delta_ratio` are input keys: the case generator takes them from a correlation (make_tbl.py, the momentum
!< integral with the van Driest II skin friction), as STREAmS takes them from its mean-profile model, never from the run.
!<
!< Face keys: `type = inflow-recycling`, `r`, `u`, `p` (the free stream), `rr_profile` (the 7-column table of the digital
!< filter), `rr_delta` (the layer thickness at the face), `rr_wall_axis` (`x`, `y` or `z`, the wall on its minimum side),
!< `rr_spacing` (the cell spacing of the plane, the cells at the face and at the station must have it), `rr_x` (the station,
!< distance from the face), `rr_delta_ratio`, `rr_beta`, `rr_shift` (fraction of the span, default 0.5). Only `[physics]
!< physical_model = euler`, one such face per realm.
!<
!< **Ghosts.** Once per step, before the first stage, from the committed state (`build_rows`): the backends write them after
!< the face pass of the boundary fill, as the digital filter's; rows beyond the face and another realm face (edges, corners)
!< are extrapolated. No state of its own: a restart continues bit for bit.

use :: adam_flume_inflow_df_object, only : load_profile, profile_mean
use :: adam_flume_parameters,       only : MODEL_EULER
use :: adam_ini_library, only : ini_check, ini_get
use :: adam_mpih_global,            only : mpih
use :: adam_parameters,             only : FEC_1_6_ARRAY, FEC_TO_DELTA
use :: finer,                       only : file_ini
use :: penf,                        only : I4P, I8P, R8P, str
use :: mpi

implicit none
private
public :: flume_inflow_rr_object
public :: rr_blend
public :: rr_damping

type :: flume_inflow_rr_object
   !< Recycling-rescaling inflow.
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
   real(R8P)                 :: delta=1._R8P             !< Layer thickness at the face.
   real(R8P)                 :: spacing=0._R8P           !< Cell spacing of the plane.
   real(R8P)                 :: x_r=0._R8P               !< Station, distance from the face (as requested).
   real(R8P)                 :: delta_ratio=1._R8P       !< delta(x_r) / delta(0).
   real(R8P)                 :: beta=1._R8P              !< delta_v(x_r) / delta_v(0).
   real(R8P)                 :: shift=0.5_R8P            !< Span shift, fraction of the span.
   character(999)            :: profile_file=''          !< Profile table file.
   real(R8P),    allocatable :: prof(:,:)                !< Profile table [7, rows].
   ! plane
   real(R8P)                 :: face_x=0._R8P            !< Face coordinate along axis_n.
   real(R8P)                 :: y0=0._R8P, z0=0._R8P     !< Plane origin (the domain minimum along the wall-normal and spanwise axes).
   integer(I4P)              :: ny=0_I4P, nz=0_I4P       !< Plane cells along the wall-normal and spanwise axes.
   integer(I4P)              :: ng=0_I4P                 !< Plane layers (the ghost cells).
   integer(I4P)              :: nr=0_I4P                 !< Translation in cells, L = nr spacing.
   integer(I4P)              :: nshift=0_I4P             !< Span shift in cells.
   logical                   :: has_plane=.false.        !< The plane has been built (the ghost rows can be listed).
   real(R8P),    allocatable :: fluct(:,:,:,:)           !< Plane fluctuations rho', u'', v'', w'' (x, y, z) [4, ng, ny, nz].
   ! recycle cells (this rank's)
   integer(I4P)              :: nrec=0_I4P               !< Recycle cells on this rank.
   integer(I4P), allocatable :: rec_cells(:,:)           !< Their cell and block (i, j, k, b) [4, nrec].
   integer(I4P), allocatable :: rec_idx(:,:)             !< Their plane layer and cell (g, jy, kz) [3, nrec].
   ! ghost rows
   integer(I4P)              :: nrows=0_I4P              !< Ghost rows of the face on this rank.
   integer(I4P), allocatable :: rows(:,:)                !< Their block and cell [4, nrows].
   real(R8P),    allocatable :: q_rows(:,:)              !< Their conservative states [5, nrows].
   contains
      procedure, pass(self) :: load               !< Load the face keys.
      procedure, pass(self) :: setup              !< Size the plane.
      procedure, pass(self) :: list_recycle_cells !< List this rank's cells of the plane.
      procedure, pass(self) :: build_plane        !< Sum the plane over the ranks and take its fluctuations.
      procedure, pass(self) :: build_rows         !< List the ghost rows of the face and their states.
      procedure, pass(self) :: ghost_state        !< Conservative state of a ghost.
endtype flume_inflow_rr_object

contains
   subroutine load(self, file_parameters, section, face, model, gamma)
   !< Load the keys of the face `face` (section `section`).
   class(flume_inflow_rr_object), intent(inout) :: self            !< Inflow.
   type(file_ini),                intent(in)    :: file_parameters !< Parameters file.
   character(*),                  intent(in)    :: section         !< Face section.
   integer(I4P),                  intent(in)    :: face            !< Face (1 to 6).
   integer(I4P),                  intent(in)    :: model           !< Physical model.
   real(R8P),                     intent(in)    :: gamma           !< Specific heats ratio.
   character(999)                               :: buff            !< Buffer.
   integer(I4P)                                 :: error           !< Error status.

   if (self%is_active) call mpih%error_stop(msg=': ['//section//'] a second inflow-recycling face: one per realm')
   if (model /= MODEL_EULER) call mpih%error_stop(msg=': ['//section//'] inflow-recycling requires [physics].(physical_model) = euler')
   self%is_active = .true.
   self%face = face
   self%axis_n = (face + 1_I4P) / 2_I4P
   self%sgn = merge(1._R8P, -1._R8P, mod(face, 2_I4P) == 1_I4P)
   self%gamma = gamma
   call get_real('r', self%r_inf)
   call get_real('u', self%u_inf)
   call get_real('p', self%p_inf)
   call get_real('rr_delta', self%delta)
   call get_real('rr_spacing', self%spacing)
   call get_real('rr_x', self%x_r)
   call get_real('rr_delta_ratio', self%delta_ratio)
   call get_real('rr_beta', self%beta)
   call file_parameters%get(section_name=section, option_name='rr_shift', val=self%shift, error=error)
   if (error /= 0) self%shift = 0.5_R8P
   if (self%shift < 0._R8P .or. self%shift >= 1._R8P) call mpih%error_stop(msg=': ['//section//'].(rr_shift) must be in [0, 1)')
   call ini_get(file_parameters, section, 'rr_wall_axis', buff)
   call ini_check()
   select case(trim(adjustl(buff)))
   case('x')
      self%axis_w = 1_I4P
   case('y')
      self%axis_w = 2_I4P
   case('z')
      self%axis_w = 3_I4P
   case default
      call mpih%error_stop(msg=': ['//section//'].(rr_wall_axis) must be x, y or z')
   endselect
   if (self%axis_w == self%axis_n) call mpih%error_stop(msg=': ['//section//'].(rr_wall_axis) is the face normal')
   self%axis_s = 6_I4P - self%axis_n - self%axis_w
   call ini_get(file_parameters, section, 'rr_profile', self%profile_file)
   self%profile_file = trim(adjustl(self%profile_file))
   call ini_check()
   call load_profile(filename=trim(self%profile_file), prof=self%prof)
   contains
      subroutine get_real(key, val)
      !< Load a required positive real key.
      character(*), intent(in)  :: key !< Key.
      real(R8P),    intent(out) :: val !< Value.

      call ini_get(file_parameters, section, key, val)
      call ini_check()
      if (.not.(val > 0._R8P)) call mpih%error_stop(msg=': ['//section//'].('//key//') must be positive')
      endsubroutine get_real
   endsubroutine load

   subroutine setup(self, emin, emax, periodic, ngc)
   !< Size the plane over the domain `[emin, emax]`: cells of `rr_spacing` along the wall-normal and spanwise axes, `ngc`
   !< layers, the translation and the span shift in cells.
   class(flume_inflow_rr_object), intent(inout) :: self        !< Inflow.
   real(R8P),                     intent(in)    :: emin(3)     !< Domain minimum.
   real(R8P),                     intent(in)    :: emax(3)     !< Domain maximum.
   logical,                       intent(in)    :: periodic(3) !< Periodic axes.
   integer(I4P),                  intent(in)    :: ngc         !< Ghost cells.

   if (.not.self%is_active) return
   self%face_x = merge(emin(self%axis_n), emax(self%axis_n), self%sgn > 0._R8P)
   self%y0 = emin(self%axis_w)
   self%z0 = emin(self%axis_s)
   self%ny = nint((emax(self%axis_w) - emin(self%axis_w)) / self%spacing, I4P)
   self%nz = nint((emax(self%axis_s) - emin(self%axis_s)) / self%spacing, I4P)
   if (abs(real(self%ny, R8P) * self%spacing - (emax(self%axis_w) - emin(self%axis_w))) > 1.e-9_R8P * self%spacing .or. &
       abs(real(self%nz, R8P) * self%spacing - (emax(self%axis_s) - emin(self%axis_s))) > 1.e-9_R8P * self%spacing)     &
      call mpih%error_stop(msg=': inflow-recycling: rr_spacing does not divide the face')
   self%ng = ngc
   self%nr = nint(self%x_r / self%spacing, I4P)
   if (self%nr < ngc) call mpih%error_stop(msg=': inflow-recycling: rr_x is closer to the face than the ghost cells')
   if (real(self%nr, R8P) * self%spacing >= emax(self%axis_n) - emin(self%axis_n)) &
      call mpih%error_stop(msg=': inflow-recycling: rr_x is beyond the domain')
   self%nshift = modulo(nint(self%shift * real(self%nz, R8P), I4P), self%nz)
   if (self%nshift /= 0_I4P .and. .not.periodic(self%axis_s)) &
      call mpih%error_stop(msg=': inflow-recycling: a span shift needs a periodic spanwise axis (rr_shift = 0 otherwise)')
   if (allocated(self%fluct)) deallocate(self%fluct)
   allocate(self%fluct(4, self%ng, self%ny, self%nz), source=0._R8P)
   print '(A)', mpih%myrankstr//'inflow-recycling: face '//trim(str(self%face, .true.))//', plane '//                    &
                trim(str(self%ny, .true.))//' x '//trim(str(self%nz, .true.))//' cells, '//trim(str(self%ng, .true.))// &
                ' layers, station '//trim(str(real(self%nr, R8P) * self%spacing))//' ('//trim(str(self%nr, .true.))// &
                ' cells), shift '//trim(str(self%nshift, .true.))//' cells, delta ratio '//trim(str(self%delta_ratio))// &
                ', beta '//trim(str(self%beta))
   endsubroutine setup

   subroutine list_recycle_cells(self, nb, nc, ngc, x_cell, y_cell, z_cell, dxyz)
   !< List this rank's cells of the plane: the interior cells whose centre is a plane layer's, `ngc` layers ending at the
   !< station, their plane layer and cell; fatal on a cell of the layers at another spacing than the plane's.
   class(flume_inflow_rr_object), intent(inout) :: self              !< Inflow.
   integer(I4P),                  intent(in)    :: nb                !< Blocks in use.
   integer(I4P),                  intent(in)    :: nc(3)             !< Interior cells per block (ni, nj, nk).
   integer(I4P),                  intent(in)    :: ngc               !< Ghost cells (the offset of the cell arrays).
   real(R8P),                     intent(in)    :: x_cell(:,:)       !< X cell centres [1-ngc:ni+ngc, nb] (passed 1-based).
   real(R8P),                     intent(in)    :: y_cell(:,:)       !< Y cell centres.
   real(R8P),                     intent(in)    :: z_cell(:,:)       !< Z cell centres.
   real(R8P),                     intent(in)    :: dxyz(:,:)         !< Cell spacing [3, nb].
   integer(I4P)                                 :: pass, n, b, c(3), g, jy, kz, i, j, k !< Counters.
   integer(I4P)                                 :: a1, a2            !< The axes across the normal.
   real(R8P)                                    :: xn, s, tol        !< Normal coordinate, its layer, tolerance.

   if (.not.self%is_active) return
   tol = 1.e-6_R8P
   a1 = merge(2_I4P, 1_I4P, self%axis_n == 1_I4P) ; a2 = 6_I4P - self%axis_n - a1
   do pass=1, 2
      n = 0_I4P
      do b=1, nb
         do i=1, nc(self%axis_n)
            c = 1_I4P ; c(self%axis_n) = i
            xn = cell_centre(self%axis_n, c, b)
            ! layer g: the centre at face_x + sgn (L - (g - 1/2) spacing)
            s = real(self%nr, R8P) + 0.5_R8P - self%sgn * (xn - self%face_x) / self%spacing
            g = nint(s, I4P)
            if (g < 1_I4P .or. g > self%ng .or. abs(s - real(g, R8P)) > tol) cycle
            if (any(abs(dxyz(:,b) - self%spacing) > tol * self%spacing)) &
               call mpih%error_stop(msg=': inflow-recycling: the station lies on cells of another spacing than rr_spacing')
            do k=1, nc(a2)
               do j=1, nc(a1)
                  c(a1) = j ; c(a2) = k
                  n = n + 1_I4P
                  if (pass == 2_I4P) then
                     jy = nint((cell_centre(self%axis_w, c, b) - self%y0) / self%spacing + 0.5_R8P, I4P)
                     kz = nint((cell_centre(self%axis_s, c, b) - self%z0) / self%spacing + 0.5_R8P, I4P)
                     self%rec_cells(:,n) = [c, b]
                     self%rec_idx(:,n) = [g, jy, kz]
                  endif
               enddo
            enddo
         enddo
      enddo
      if (pass == 1_I4P) then
         if (allocated(self%rec_cells)) deallocate(self%rec_cells, self%rec_idx)
         allocate(self%rec_cells(4, n), self%rec_idx(3, n))
         self%nrec = n
      endif
   enddo
   contains
      pure function cell_centre(a, cc, bb) result(x)
      !< The centre of cell `cc` of block `bb` along axis `a`.
      integer(I4P), intent(in) :: a, cc(3), bb !< Axis, cell, block.
      real(R8P)                :: x            !< Centre.

      select case(a)
      case(1)
         x = x_cell(cc(1)+ngc,bb)
      case(2)
         x = y_cell(cc(2)+ngc,bb)
      case default
         x = z_cell(cc(3)+ngc,bb)
      endselect
      endfunction cell_centre
   endsubroutine list_recycle_cells

   subroutine build_plane(self, q_rec)
   !< Sum this rank's plane cells `q_rec` (the conservative state of `rec_cells`) over the ranks into the global plane and take its
   !< fluctuations: per layer and wall-normal row the spanwise Favre mean velocity and the mean density.
   class(flume_inflow_rr_object), intent(inout) :: self         !< Inflow.
   real(R8P),                     intent(in)    :: q_rec(:,:)   !< Conservative state of this rank's plane cells [5, nrec].
   real(R8P), allocatable                       :: plane(:,:,:,:) !< Plane: rho, rho u, rho v, rho w, writes [5, ng, ny, nz].
   real(R8P)                                    :: srho, sm(3)  !< Spanwise sums.
   integer(I4P)                                 :: m, g, j, k    !< Counters.

   if (.not.self%is_active) return
   allocate(plane(5, self%ng, self%ny, self%nz), source=0._R8P)
   do m=1, self%nrec
      associate(ix=>self%rec_idx(:,m))
      plane(1:4,ix(1),ix(2),ix(3)) = q_rec(1:4,m)
      plane(5,ix(1),ix(2),ix(3)) = plane(5,ix(1),ix(2),ix(3)) + 1._R8P
      endassociate
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, plane, size(plane), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   if (any(plane(5,:,:,:) /= 1._R8P)) &
      call mpih%error_stop(msg=': inflow-recycling: the plane cells are not all covered once (a station on a coarser or a finer level?)')
   do j=1, self%ny
      do g=1, self%ng
         srho = sum(plane(1,g,j,:))
         sm = [sum(plane(2,g,j,:)), sum(plane(3,g,j,:)), sum(plane(4,g,j,:))]
         do k=1, self%nz
            self%fluct(1,g,j,k) = plane(1,g,j,k) - srho / real(self%nz, R8P)
            self%fluct(2:4,g,j,k) = plane(2:4,g,j,k) / plane(1,g,j,k) - sm / srho
         enddo
      enddo
   enddo
   self%has_plane = .true.
   endsubroutine build_plane

   subroutine build_rows(self, crown_map, x_cell, y_cell, z_cell)
   !< List the crown rows beyond the face only (block > 0) and their states (`ghost_state`, from the plane of `build_plane`).
   class(flume_inflow_rr_object), intent(inout) :: self              !< Inflow.
   integer(I8P),                  intent(in)    :: crown_map(:,:,:)  !< BC crown map [rows, 9, ngc].
   real(R8P),                     intent(in)    :: x_cell(:,:)       !< X cell centres [1-ngc:ni+ngc, nb] (passed 1-based).
   real(R8P),                     intent(in)    :: y_cell(:,:)       !< Y cell centres.
   real(R8P),                     intent(in)    :: z_cell(:,:)       !< Z cell centres.
   integer(I4P)                                 :: ngc               !< Ghost cells (the offset of the cell arrays).
   integer(I4P)                                 :: crown, c, n, b, fec, idx(3), g, kz !< Counters, block, fec, cell, layer, span cell.
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
         g = max(1_I4P, min(self%ng, nint(self%sgn * (self%face_x - xyz(self%axis_n)) / self%spacing + 0.5_R8P, I4P)))
         kz = nint((xyz(self%axis_s) - self%z0) / self%spacing + 0.5_R8P, I4P)
         kz = modulo(kz - 1_I4P + self%nshift, self%nz) + 1_I4P
         self%q_rows(:,n) = self%ghost_state(g=g, y=xyz(self%axis_w) - self%y0, kz=kz)
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

   pure function ghost_state(self, g, y, kz) result(q)
   !< Conservative state of a ghost of layer `g`, height `y` above the wall, reading the plane's span cell `kz` (shifted).
   class(flume_inflow_rr_object), intent(in) :: self !< Inflow.
   integer(I4P),                  intent(in) :: g    !< Plane layer.
   real(R8P),                     intent(in) :: y    !< Height above the wall.
   integer(I4P),                  intent(in) :: kz   !< Plane span cell (shifted).
   real(R8P)                                 :: q(5) !< Conservative state.
   real(R8P)                                 :: eta, w, d, um, tm, r(4), rm, rho, vel(3) !< Work.
   real(R8P)                                 :: fi(4), fo(4), f(4) !< Inner, outer, blended fluctuations.
   logical                                   :: oki, oko !< Sources within the plane.

   eta = y / self%delta
   call profile_mean(prof=self%prof, eta=eta, u=um, t=tm, r=r)
   rm = self%r_inf / tm
   call plane_at(y * self%beta, fi, oki)
   call plane_at(y * self%delta_ratio, fo, oko)
   f = 0._R8P
   if (oki .and. oko) then
      w = rr_blend(eta)
      d = rr_damping(eta)
      f = ((1._R8P - w) * fi + w * fo) * d
      f(2:4) = f(2:4) * self%beta
   endif
   rho = rm + max(-0.75_R8P * rm, min(4._R8P * rm, f(1)))
   vel = f(2:4)
   vel(self%axis_n) = vel(self%axis_n) + self%sgn * um * self%u_inf
   q(1) = rho
   q(2:4) = rho * vel
   q(5) = self%p_inf / (self%gamma - 1._R8P) + 0.5_R8P * rho * sum(vel * vel)
   contains
      pure subroutine plane_at(ys, fv, ok)
      !< The fluctuations at height `ys` of layer g, span cell kz, linear between the plane's rows; none above the last row's centre.
      real(R8P), intent(in)  :: ys    !< Height.
      real(R8P), intent(out) :: fv(4) !< Fluctuations.
      logical,   intent(out) :: ok    !< Within the plane.
      real(R8P)              :: s     !< Position in rows.
      integer(I4P)           :: jl    !< Row below.

      s = ys / self%spacing + 0.5_R8P
      ok = s <= real(self%ny, R8P) + 1.e-9_R8P ! up to the last row's centre
      fv = 0._R8P
      if (.not.ok) return
      jl = floor(s, I4P)
      if (jl < 1_I4P) then
         fv = self%fluct(:,g,1,kz)
      elseif (jl >= self%ny) then
         fv = self%fluct(:,g,self%ny,kz)
      else
         s = s - real(jl, R8P)
         fv = (1._R8P - s) * self%fluct(:,g,jl,kz) + s * self%fluct(:,g,jl+1,kz)
      endif
      endsubroutine plane_at
   endfunction ghost_state

   ! public procedures
   pure function rr_blend(eta) result(w)
   !< The blend of the inner and outer fluctuations, `W(eta) = (1 + tanh(4 (eta - 0.2) / (0.2 + 0.6 eta)) / tanh(4)) / 2`
   !< (STREAmS-2 `recyc_prepare`, alfa 4, beta 0.2): 0 at the wall, 1/2 at eta = 0.2, 1 at eta = 1.
   real(R8P), intent(in) :: eta !< y / delta.
   real(R8P)             :: w   !< Weight of the outer fluctuations.

   w = 0.5_R8P * (1._R8P + tanh(4._R8P * (eta - 0.2_R8P) / (0.2_R8P + 0.6_R8P * eta)) / tanh(4._R8P))
   endfunction rr_blend

   pure function rr_damping(eta) result(d)
   !< The damping of the fluctuations beyond the layer, `(1 - tanh(4 (eta - 2))) / 2` (STREAmS-2 `bcrecyc`).
   real(R8P), intent(in) :: eta !< y / delta.
   real(R8P)             :: d   !< Damping.

   d = 0.5_R8P * (1._R8P - tanh(4._R8P * (eta - 2._R8P)))
   endfunction rr_damping
endmodule adam_flume_inflow_rr_object
