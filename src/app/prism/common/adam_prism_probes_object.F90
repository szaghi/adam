!< ADAM, PRISM field probes: time histories of the fields at fixed points and on planes of points.
module adam_prism_probes_object
!< ADAM, PRISM field probes: time histories of the fields at fixed points and on planes of points.
!<
!< The input is optional: without a `[probes]` section, or with neither points nor planes, there are no probes.
!<```ini
!<[probes]
!<save_every = 5                    ; cadence in iterations (> 0), the last step is always saved
!<point_1    = -0.04 0.0 0.0        ; a probe at x y z
!<point_2    = 0.0 0.01 0.0
!<plane_1    = x 0.4364 48 48       ; a plane of probes normal to x at x = 0.4364: the centres of a uniform 48 x 48
!<                                  ; lattice over the domain extent along the other two axes (y first, then z)
!<```
!< Coordinates are in the units of the [grid] extents (divided by L0 in an adimensional model, as the grid). The keys
!< `point_N` and `plane_N` are numbered from 1 without gaps (the first missing number ends the list). The plane lattice
!< spans the whole domain along its two axes: with n equal to the cells of the grid along an axis its points are the
!< cell centres of the uniform grid; a plane normal to `y` runs over (x, z), normal to `z` over (x, y).
!< Each probe samples the whole `q` by trilinear interpolation, written as tent weights on the cell centres,
!< `w = prod_d max(0, 1 - |x_d - xc_d| / dx_d)` (null axes excluded), the weights of the domain-center B diagnostic.
!< Every rank sums the cells of its own blocks and a reduction to rank 0 recombines the interpolation across blocks and
!< ranks: no ghost cell is read. The weights are normalised by their global sum, so within half a cell of the domain
!< boundary, periodic faces included, the interpolation uses the interior cells only.
!< The stencils (local cells and weights) are built once on the current grid by `build_stencils`; a regrid must
!< rebuild them.
!<
!< Output, rank 0 only: one file per point, `<basename>-probe_NN.dat` (one row per save: it, time, the variables), and
!< one file per plane, `<basename>-plane_NN.dat` (per save, one row per probe of the plane: it, time, i, j, the
!< variables, i running fastest), NN with at least two digits. The rows are not flushed: a flush per file and per save
!< costs a file-system round trip each (on a shared file system with thousands of files, seconds to tens of seconds
!< per save, with every rank waiting for rank 0). The runtime buffers reach the files when full and at `close_files`;
!< after a crash the rows still buffered are lost.
!< The former input, `probes_number` and one `[probe_N]` section per probe, is refused.

use :: adam_common_library, only : field_object, grid_object
use :: adam_mpih_global, only : mpih
use :: finer
use :: mpi
use :: penf
use :: stringifor

implicit none
private
public :: prism_probes_object

character(len=6), parameter :: PROBES_SECTION_NAME = 'probes' !< Main section.
character(len=*), parameter :: ROW_FORMAT = '(*(1X,ES23.15E3))' !< Format of the values of a history row.

type :: prism_probes_object
   !< PRISM field probes: configuration, local stencils, history writer.
   integer(I4P)              :: probes_number = 0_I4P !< Number of probes, points and plane probes.
   integer(I4P)              :: points_number = 0_I4P !< Number of single points (the first probes).
   integer(I4P)              :: planes_number = 0_I4P !< Number of planes (their probes follow the points).
   integer(I4P)              :: save_every = 0_I4P    !< Saving cadence in iterations.
   integer(I4P)              :: nv = 0_I4P            !< Number of sampled variables.
   logical                   :: is_built = .false.    !< Stencils built on the current grid.
   character(:), allocatable :: output_basename       !< Basename of output files.
   real(R8P),    allocatable :: point(:,:)            !< Probe points, grid units [3,probes_number].
   real(R8P),    allocatable :: point_input(:,:)      !< Probe points, input units [3,probes_number].
   character(1), allocatable :: plane_axis(:)         !< Normal axis of each plane, x/y/z [planes_number].
   real(R8P),    allocatable :: plane_abscissa(:)     !< Abscissa of each plane along its normal, input units.
   integer(I4P), allocatable :: plane_n(:,:)          !< Lattice of each plane, (n1,n2) [2,planes_number].
   integer(I4P), allocatable :: plane_first(:)        !< First probe of each plane [planes_number+1].
   integer(I4P), allocatable :: unit(:)               !< History file units (rank 0), points then planes.
   ! local stencils, compressed by probe: the cells of probe p are first(p):first(p+1)-1
   integer(I4P), allocatable :: first(:)              !< First local cell of each probe [probes_number+1].
   integer(I4P), allocatable :: cell(:,:)             !< Local cells, (b,i,j,k) [4,entries].
   real(R8P),    allocatable :: weight(:)             !< Local cell weights, normalised [entries].
   real(R8P),    allocatable :: value(:,:)            !< Sampled values, local sums then global [nv,probes_number].
   contains
      ! public methods
      procedure, pass(self) :: build_stencils !< Build the local stencils on the current grid.
      procedure, pass(self) :: close_files    !< Close the history files.
      procedure, pass(self) :: initialize     !< Initialize the probes from the input file.
      procedure, pass(self) :: label          !< Input name of a probe, for the messages.
      procedure, pass(self) :: open_files     !< Open the history files and write their headers.
      procedure, pass(self) :: sample_host    !< Sample the host fields (local sums).
      procedure, pass(self) :: save_history   !< Reduce the local sums and write the history rows.
endtype prism_probes_object

contains
   ! public methods
   subroutine build_stencils(self, grid, field, ni, nj, nk, blocks_number)
   !< Build the local stencils on the current grid; collective (normalises the weights by their global sum).
   class(prism_probes_object), intent(inout) :: self                 !< Probes.
   type(grid_object),          intent(in)    :: grid                 !< Grid.
   type(field_object),         intent(in)    :: field                !< Field (cell centres, steps).
   integer(I4P),               intent(in)    :: ni, nj, nk           !< Interior cells per block.
   integer(I4P),               intent(in)    :: blocks_number        !< Local blocks.
   real(R8P),    allocatable                 :: wx(:), wy(:), wz(:)  !< 1-D tent weights of a block.
   real(R8P),    allocatable                 :: weight_sum(:)        !< Global weight sum of each probe.
   integer(I4P), allocatable                 :: cell(:,:)            !< Local cells, growing.
   real(R8P),    allocatable                 :: weight(:)            !< Local weights, growing.
   integer(I4P)                              :: entries              !< Local cells counter.
   integer(I4P)                              :: p, b, i, j, k, s     !< Counters.

   if (self%probes_number == 0_I4P) return
   allocate(wx(ni), wy(nj), wz(nk))
   allocate(cell(4,8*self%probes_number), weight(8*self%probes_number))
   if (allocated(self%first)) deallocate(self%first)
   allocate(self%first(self%probes_number+1))
   entries = 0_I4P
   do p=1, self%probes_number
      self%first(p) = entries + 1_I4P
      do b=1, blocks_number
         call tent_weights(x=self%point(1,p), xc=field%x_cell(1:ni,b), dx=field%dxyz(1,b), is_null=grid%null_xyz(1), w=wx)
         call tent_weights(x=self%point(2,p), xc=field%y_cell(1:nj,b), dx=field%dxyz(2,b), is_null=grid%null_xyz(2), w=wy)
         call tent_weights(x=self%point(3,p), xc=field%z_cell(1:nk,b), dx=field%dxyz(3,b), is_null=grid%null_xyz(3), w=wz)
         if (all(wx <= 0._R8P) .or. all(wy <= 0._R8P) .or. all(wz <= 0._R8P)) cycle
         do k=1, nk
            if (wz(k) <= 0._R8P) cycle
            do j=1, nj
               if (wy(j) <= 0._R8P) cycle
               do i=1, ni
                  if (wx(i) <= 0._R8P) cycle
                  entries = entries + 1_I4P
                  if (entries > size(weight)) call grow(cell=cell, weight=weight)
                  cell(:,entries) = [b, i, j, k]
                  weight(entries) = wx(i) * wy(j) * wz(k)
               enddo
            enddo
         enddo
      enddo
   enddo
   self%first(self%probes_number+1) = entries + 1_I4P
   allocate(weight_sum(self%probes_number))
   do p=1, self%probes_number
      weight_sum(p) = sum(weight(self%first(p):self%first(p+1)-1))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, weight_sum, self%probes_number, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   do p=1, self%probes_number
      if (weight_sum(p) <= 0._R8P) call mpih%error_stop(msg=': [probes] '//self%label(p)//                          &
                                                            ' has no cell centre within one cell of its point')
      do s=self%first(p), self%first(p+1)-1
         weight(s) = weight(s) / weight_sum(p)
      enddo
   enddo
   if (allocated(self%cell)) deallocate(self%cell)
   if (allocated(self%weight)) deallocate(self%weight)
   allocate(self%cell(4,entries), source=cell(:,1:entries))
   allocate(self%weight(entries), source=weight(1:entries))
   self%is_built = .true.
   contains
      subroutine grow(cell, weight)
      !< Double the capacity of the local stencil buffers.
      !<
      !< The buffers are explicit dummies, not host-associated: with host association nvfortran 26.9 (also at -O0)
      !< keeps the old `size(weight)` in the caller after `move_alloc`, so every later entry doubles them again.
      integer(I4P), allocatable, intent(inout) :: cell(:,:) !< Local cells.
      real(R8P),    allocatable, intent(inout) :: weight(:) !< Local weights.
      integer(I4P), allocatable                :: cell_(:,:)  !< Cells, enlarged.
      real(R8P),    allocatable                :: weight_(:)  !< Weights, enlarged.

      allocate(cell_(4,2*size(weight)), weight_(2*size(weight)))
      cell_(:,1:size(weight)) = cell
      weight_(1:size(weight)) = weight
      call move_alloc(cell_, cell)
      call move_alloc(weight_, weight)
      endsubroutine grow
   endsubroutine build_stencils

   subroutine close_files(self)
   !< Close the history files.
   class(prism_probes_object), intent(inout) :: self    !< Probes.
   logical                                   :: is_open !< Inquiring flag.
   integer(I4P)                              :: p       !< Counter.

   if (.not.allocated(self%unit)) return
   do p=1, size(self%unit)
      if (self%unit(p) > 0_I4P) then
         inquire(unit=self%unit(p), opened=is_open)
         if (is_open) close(self%unit(p))
      endif
   enddo
   endsubroutine close_files

   subroutine initialize(self, file_parameters, output_basename, grid, L0)
   !< Initialize the probes from the input file; the `[probes]` section is optional (absent = no probes).
   class(prism_probes_object), intent(inout)        :: self            !< Probes.
   type(file_ini),             intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   character(*),               intent(in)           :: output_basename !< Basename of output files.
   type(grid_object),          intent(in)           :: grid            !< Grid, for the domain check.
   real(R8P),                  intent(in), optional :: L0              !< Adimensionalization length, as the grid's.
   character(:), allocatable                        :: key             !< Option name.
   character(:), allocatable                        :: val             !< Option value.
   real(R8P),    allocatable                        :: points(:,:)     !< Single points, input units, growing.
   real(R8P)                                        :: scale           !< Grid units per input unit.
   real(R8P)                                        :: emin(3), emax(3)!< Domain extents, input units.
   real(R8P)                                        :: xyz(3)          !< A point.
   character(1)                                     :: axis            !< Plane normal.
   real(R8P)                                        :: abscissa        !< Plane abscissa.
   integer(I4P)                                     :: n(2)            !< Plane lattice.
   integer(I4P)                                     :: a(3)            !< Plane axes: normal, first, second.
   integer(I4P)                                     :: error           !< Error status.
   integer(I4P)                                     :: p, l, i, j      !< Counters.
   character(1), parameter                          :: AXIS_NAME(3) = ['x', 'y', 'z'] !< Coordinate names.

   self%probes_number = 0_I4P
   self%points_number = 0_I4P
   self%planes_number = 0_I4P
   self%is_built = .false.
   self%output_basename = trim(output_basename)
   if (file_parameters%has_section(section_name='probe_1')) &
      call mpih%error_stop(msg=': [probe_N] sections are no longer read: write each probe as ['//                     &
                               PROBES_SECTION_NAME//'] point_N = x y z, or a plane of probes as plane_N = axis abscissa n1 n2')
   if (.not.file_parameters%has_section(section_name=PROBES_SECTION_NAME)) return
   call file_parameters%get_string(section_name=PROBES_SECTION_NAME, option_name='probes_number', val=val, error=error)
   if (error == 0) &
      call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].(probes_number) is no longer read: the probes are the '//   &
                               'keys point_N = x y z and plane_N = axis abscissa n1 n2')
   scale = 1._R8P
   if (present(L0)) scale = 1._R8P / L0
   emin = grid%domain_emin / scale
   emax = grid%domain_emax / scale
   ! single points
   allocate(points(3,8))
   do
      key = 'point_'//trim(str(self%points_number+1_I4P, .true.))
      call file_parameters%get_string(section_name=PROBES_SECTION_NAME, option_name=key, val=val, error=error)
      if (error /= 0) exit
      read(val, *, iostat=error) xyz
      if (error /= 0) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') = "'//val//                  &
                                               '" is not x y z')
      if (any(xyz < emin) .or. any(xyz > emax)) &
         call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') point '//trim(str(xyz))//                 &
                                  ' is outside the domain')
      self%points_number = self%points_number + 1_I4P
      if (self%points_number > size(points, dim=2)) call grow(points)
      points(:,self%points_number) = xyz
   enddo
   ! planes, first pass: count and check
   do
      key = 'plane_'//trim(str(self%planes_number+1_I4P, .true.))
      call file_parameters%get_string(section_name=PROBES_SECTION_NAME, option_name=key, val=val, error=error)
      if (error /= 0) exit
      call parse_plane(key, val, axis, abscissa, n, a)
      self%planes_number = self%planes_number + 1_I4P
   enddo
   if (allocated(self%plane_axis)) deallocate(self%plane_axis)
   if (allocated(self%plane_abscissa)) deallocate(self%plane_abscissa)
   if (allocated(self%plane_n)) deallocate(self%plane_n)
   if (allocated(self%plane_first)) deallocate(self%plane_first)
   allocate(self%plane_axis(self%planes_number), self%plane_abscissa(self%planes_number))
   allocate(self%plane_n(2,self%planes_number), self%plane_first(self%planes_number+1))
   self%probes_number = self%points_number
   do l=1, self%planes_number
      key = 'plane_'//trim(str(l, .true.))
      call file_parameters%get_string(section_name=PROBES_SECTION_NAME, option_name=key, val=val, error=error)
      call parse_plane(key, val, self%plane_axis(l), self%plane_abscissa(l), self%plane_n(:,l), a)
      self%plane_first(l) = self%probes_number + 1_I4P
      self%probes_number = self%probes_number + self%plane_n(1,l) * self%plane_n(2,l)
   enddo
   self%plane_first(self%planes_number+1) = self%probes_number + 1_I4P
   if (self%probes_number == 0_I4P) return
   call file_parameters%get(section_name=PROBES_SECTION_NAME, option_name='save_every', val=self%save_every, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//PROBES_SECTION_NAME//'].(save_every)')
   if (self%save_every <= 0_I4P) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].(save_every) must be > 0')
   if (allocated(self%point)) deallocate(self%point)
   if (allocated(self%point_input)) deallocate(self%point_input)
   allocate(self%point(3,self%probes_number), self%point_input(3,self%probes_number))
   self%point_input(:,1:self%points_number) = points(:,1:self%points_number)
   do l=1, self%planes_number
      a = plane_axes(self%plane_axis(l))
      p = self%plane_first(l)
      do j=1, self%plane_n(2,l)
         do i=1, self%plane_n(1,l)
            self%point_input(a(1),p) = self%plane_abscissa(l)
            self%point_input(a(2),p) = emin(a(2)) + (i - 0.5_R8P) * (emax(a(2)) - emin(a(2))) / self%plane_n(1,l)
            self%point_input(a(3),p) = emin(a(3)) + (j - 0.5_R8P) * (emax(a(3)) - emin(a(3))) / self%plane_n(2,l)
            p = p + 1_I4P
         enddo
      enddo
   enddo
   self%point = self%point_input
   if (present(L0)) self%point = self%point / L0 ! as before the compact input: bitwise the same stencil weights
   call mpih%print_message('prism_probes_object initialized: '//trim(str(self%points_number, .true.))//' points and '// &
                           trim(str(self%planes_number, .true.))//' planes, '//trim(str(self%probes_number, .true.))//  &
                           ' probes, saved every '//trim(str(self%save_every, .true.))//' iterations')
   contains
      subroutine parse_plane(key, val, axis, abscissa, n, a)
      !< Parse and check `plane_N = axis abscissa n1 n2`.
      character(*), intent(in)  :: key      !< Option name.
      character(*), intent(in)  :: val      !< Option value.
      character(1), intent(out) :: axis     !< Normal axis.
      real(R8P),    intent(out) :: abscissa !< Abscissa along the normal.
      integer(I4P), intent(out) :: n(2)     !< Lattice.
      integer(I4P), intent(out) :: a(3)     !< Axes: normal, first, second.
      integer(I4P)              :: ios      !< I/O status.

      read(val, *, iostat=ios) axis, abscissa, n
      if (ios /= 0) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') = "'//val//                    &
                                             '" is not "axis abscissa n1 n2"')
      if (all(axis /= AXIS_NAME)) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') axis "'//axis// &
                                                           '" is not x, y or z')
      a = plane_axes(axis)
      if (abscissa < emin(a(1)) .or. abscissa > emax(a(1))) &
         call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') abscissa '//trim(str(abscissa))//         &
                                  ' is outside the domain')
      if (any(n < 1_I4P)) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].('//key//') n1 n2 must be >= 1')
      endsubroutine parse_plane

      subroutine grow(points)
      !< Double the capacity of the points buffer (explicit dummy: see `build_stencils`).
      real(R8P), allocatable, intent(inout) :: points(:,:) !< Points.
      real(R8P), allocatable                :: points_(:,:) !< Points, enlarged.

      allocate(points_(3,2*size(points, dim=2)))
      points_(:,1:size(points, dim=2)) = points
      call move_alloc(points_, points)
      endsubroutine grow
   endsubroutine initialize

   pure function label(self, p) result(name)
   !< Input name of probe `p`: `point_N`, or `plane_N (i,j)`.
   class(prism_probes_object), intent(in) :: self !< Probes.
   integer(I4P),               intent(in) :: p    !< Probe.
   character(:), allocatable              :: name !< Name.
   integer(I4P)                           :: l, q !< Plane, its probe offset.

   if (p <= self%points_number) then
      name = 'point_'//trim(str(p, .true.))
      return
   endif
   do l=1, self%planes_number
      if (p < self%plane_first(l+1)) exit
   enddo
   q = p - self%plane_first(l)
   name = 'plane_'//trim(str(l, .true.))//' ('//trim(str(mod(q, self%plane_n(1,l)) + 1_I4P, .true.))//','//            &
          trim(str(q / self%plane_n(1,l) + 1_I4P, .true.))//')'
   endfunction label

   subroutine open_files(self, q_name)
   !< Open the history files (rank 0) and write their headers.
   class(prism_probes_object), intent(inout) :: self      !< Probes.
   type(string),               intent(in)    :: q_name(:) !< Names of the sampled variables.
   character(:), allocatable                 :: header    !< Variables header.
   integer(I4P)                              :: p, l, v   !< Counters.
   integer(I4P)                              :: nz        !< Digits of the file numbers.

   self%nv = size(q_name, kind=I4P)
   if (allocated(self%value)) deallocate(self%value)
   allocate(self%value(self%nv,self%probes_number))
   self%value = 0._R8P
   if (mpih%myrank /= 0) return
   call self%close_files
   if (allocated(self%unit)) deallocate(self%unit)
   allocate(self%unit(self%points_number+self%planes_number))
   header = ''
   do v=1, self%nv
      header = header//' "'//q_name(v)%chars()//'"'
   enddo
   ! strz keeps the last nz digits: at least 2, as many as the largest number
   nz = max(2_I4P, len_trim(str(self%points_number, .true.), kind=I4P))
   do p=1, self%points_number
      open(newunit=self%unit(p), file=self%output_basename//'-probe_'//trim(strz(p, nz))//'.dat')
      write(self%unit(p), '(A)') '% probe '//trim(str(p, .true.))//' at x y z = '//trim(str(self%point_input(:,p)))
      write(self%unit(p), '(A)') '%VARIABLES="it" "time"'//header
   enddo
   nz = max(2_I4P, len_trim(str(self%planes_number, .true.), kind=I4P))
   do l=1, self%planes_number
      p = self%points_number + l
      open(newunit=self%unit(p), file=self%output_basename//'-plane_'//trim(strz(l, nz))//'.dat')
      write(self%unit(p), '(A)') '% plane '//trim(str(l, .true.))//' normal to '//self%plane_axis(l)//' at '//         &
                                 trim(str(self%plane_abscissa(l)))//', lattice '//trim(str(self%plane_n(1,l), .true.))// &
                                 ' x '//trim(str(self%plane_n(2,l), .true.))//', i fastest'
      write(self%unit(p), '(A)') '%VARIABLES="it" "time" "i" "j"'//header
   enddo
   endsubroutine open_files

   subroutine sample_host(self, ngc, q)
   !< Sample the host fields: local weighted sums of the stencil cells.
   class(prism_probes_object), intent(inout) :: self        !< Probes.
   integer(I4P),               intent(in)    :: ngc         !< Number of ghost cells.
   real(R8P),                  intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fields [nv,ni,nj,nk,nb], with ghosts.
   integer(I4P)                              :: p, s        !< Counters.

   do p=1, self%probes_number
      self%value(:,p) = 0._R8P
      do s=self%first(p), self%first(p+1)-1
         self%value(:,p) = self%value(:,p) + self%weight(s) * &
                           q(1:self%nv,self%cell(2,s),self%cell(3,s),self%cell(4,s),self%cell(1,s))
      enddo
   enddo
   endsubroutine sample_host

   subroutine save_history(self, it, time)
   !< Reduce the local sums to rank 0 and write the history rows (one per point, one per plane probe); collective.
   class(prism_probes_object), intent(inout) :: self !< Probes.
   integer(I4P),               intent(in)    :: it   !< Current iteration.
   real(R8P),                  intent(in)    :: time !< Current time.
   real(R8P)                                 :: dummy(1) !< Receive buffer of the non-root ranks (not significant).
   integer(I4P)                              :: p, l, q  !< Counters.
   integer(I4P)                              :: u        !< Unit.

   if (mpih%myrank == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, self%value, self%nv*self%probes_number, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, &
                      mpih%error)
   else
      call MPI_REDUCE(self%value, dummy, self%nv*self%probes_number, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,        &
                      mpih%error)
      return
   endif
   do p=1, self%points_number
      write(self%unit(p), '(I0)', advance='no') it
      write(self%unit(p), ROW_FORMAT) time, self%value(:,p)
   enddo
   do l=1, self%planes_number
      u = self%unit(self%points_number+l)
      do p=self%plane_first(l), self%plane_first(l+1)-1
         q = p - self%plane_first(l)
         write(u, '(I0,1X)', advance='no') it
         write(u, '(ES23.15E3,1X,I0,1X,I0)', advance='no') time, mod(q, self%plane_n(1,l)) + 1_I4P, &
                                                            q / self%plane_n(1,l) + 1_I4P
         write(u, ROW_FORMAT) self%value(:,p)
      enddo
   enddo
   endsubroutine save_history

   ! private non TBP
   pure function plane_axes(axis) result(a)
   !< Axes of a plane: its normal, then the first and the second lattice axes.
   character(1), intent(in) :: axis !< Normal axis.
   integer(I4P)             :: a(3) !< Normal, first, second.

   select case(axis)
   case('x')
      a = [1_I4P, 2_I4P, 3_I4P]
   case('y')
      a = [2_I4P, 1_I4P, 3_I4P]
   case default
      a = [3_I4P, 1_I4P, 2_I4P]
   endselect
   endfunction plane_axes

   pure subroutine tent_weights(x, xc, dx, is_null, w)
   !< 1-D tent weights of the cell centres `xc` at `x`; a null axis weighs every cell 1.
   real(R8P),    intent(in)  :: x       !< Point abscissa.
   real(R8P),    intent(in)  :: xc(:)   !< Cell centres.
   real(R8P),    intent(in)  :: dx      !< Cell step.
   logical,      intent(in)  :: is_null !< Null axis.
   real(R8P),    intent(out) :: w(:)    !< Weights.

   if (is_null) then
      w = 1._R8P
   else
      w = max(0._R8P, 1._R8P - abs(x - xc) / dx)
   endif
   endsubroutine tent_weights
endmodule adam_prism_probes_object
