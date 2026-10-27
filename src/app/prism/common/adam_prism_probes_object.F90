!< ADAM, PRISM field probes: time histories of the fields at fixed points.
module adam_prism_probes_object
!< ADAM, PRISM field probes: time histories of the fields at fixed points.
!<
!< The input is optional: without a `[probes]` section there are no probes.
!<```ini
!<[probes]
!<probes_number = 3 ; 0 = no probes
!<save_every    = 5 ; cadence in iterations (> 0), the last step is always saved
!<[probe_1]
!<x = -0.04         ; in the units of the [grid] extents (divided by L0 in an adimensional model, as the grid)
!<y = 0.0
!<z = 0.0
!<```
!< Each probe samples the whole `q` by trilinear interpolation, written as tent weights on the cell centres,
!< `w = prod_d max(0, 1 - |x_d - xc_d| / dx_d)` (null axes excluded), the weights of the domain-center B diagnostic.
!< Every rank sums the cells of its own blocks and a reduction to rank 0 recombines the interpolation across blocks and
!< ranks: no ghost cell is read. The weights are normalised by their global sum, so within half a cell of the domain
!< boundary, periodic faces included, the interpolation uses the interior cells only.
!< The stencils (local cells and weights) are built once on the current grid by `build_stencils`; a regrid must
!< rebuild them. Output: one file per probe, `<basename>-probe_NN.dat` (more digits beyond 99 probes), rank 0 only.

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

type :: prism_probes_object
   !< PRISM field probes: configuration, local stencils, history writer.
   integer(I4P)              :: probes_number = 0_I4P !< Number of probes.
   integer(I4P)              :: save_every = 0_I4P    !< Saving cadence in iterations.
   integer(I4P)              :: nv = 0_I4P            !< Number of sampled variables.
   logical                   :: is_built = .false.    !< Stencils built on the current grid.
   character(:), allocatable :: output_basename       !< Basename of output files.
   real(R8P),    allocatable :: point(:,:)            !< Probe points, grid units [3,probes_number].
   real(R8P),    allocatable :: point_input(:,:)      !< Probe points, input units [3,probes_number].
   integer(I4P), allocatable :: unit(:)               !< History file units (rank 0) [probes_number].
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
      procedure, pass(self) :: open_files     !< Open the history files and write their headers.
      procedure, pass(self) :: sample_host    !< Sample the host fields (local sums).
      procedure, pass(self) :: save_history   !< Reduce the local sums and write one history row per probe.
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
      if (weight_sum(p) <= 0._R8P) call mpih%error_stop(msg=': [probe_'//trim(str(p, .true.))//                       &
                                                            '] no cell centre within one cell of the point')
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
   character(:), allocatable                        :: section         !< Probe section name.
   integer(I4P)                                     :: error           !< Error status.
   integer(I4P)                                     :: p, d            !< Counters.
   character(1), parameter                          :: AXIS(3) = ['x', 'y', 'z'] !< Coordinate keys.

   self%probes_number = 0_I4P
   self%is_built = .false.
   self%output_basename = trim(output_basename)
   if (.not.file_parameters%has_section(section_name=PROBES_SECTION_NAME)) return
   call file_parameters%get(section_name=PROBES_SECTION_NAME, option_name='probes_number', val=self%probes_number, &
                            error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//PROBES_SECTION_NAME//'].(probes_number)')
   if (self%probes_number < 0_I4P) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].(probes_number) < 0')
   if (self%probes_number == 0_I4P) return
   call file_parameters%get(section_name=PROBES_SECTION_NAME, option_name='save_every', val=self%save_every, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//PROBES_SECTION_NAME//'].(save_every)')
   if (self%save_every <= 0_I4P) call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].(save_every) must be > 0')
   if (allocated(self%point)) deallocate(self%point)
   if (allocated(self%point_input)) deallocate(self%point_input)
   allocate(self%point(3,self%probes_number), self%point_input(3,self%probes_number))
   do p=1, self%probes_number
      section = 'probe_'//trim(str(p, .true.))
      if (.not.file_parameters%has_section(section_name=section)) &
         call mpih%error_stop(msg=': ['//PROBES_SECTION_NAME//'].(probes_number) = '//                                &
                                  trim(str(self%probes_number, .true.))//' but ['//section//'] is missing')
      do d=1, 3
         call file_parameters%get(section_name=section, option_name=AXIS(d), val=self%point_input(d,p), error=error)
         if (error /= 0) call mpih%error_stop(msg=': failed to load ['//section//'].('//AXIS(d)//')')
      enddo
      self%point(:,p) = self%point_input(:,p)
      if (present(L0)) self%point(:,p) = self%point(:,p) / L0
      if (any(self%point(:,p) < grid%domain_emin) .or. any(self%point(:,p) > grid%domain_emax)) &
         call mpih%error_stop(msg=': ['//section//'] point '//trim(str(self%point_input(:,p)))//                       &
                                  ' is outside the domain')
   enddo
   call mpih%print_message('prism_probes_object initialized: '//trim(str(self%probes_number, .true.))//                 &
                           ' probes, saved every '//trim(str(self%save_every, .true.))//' iterations')
   endsubroutine initialize

   subroutine open_files(self, q_name)
   !< Open the history files (rank 0) and write their headers.
   class(prism_probes_object), intent(inout) :: self      !< Probes.
   type(string),               intent(in)    :: q_name(:) !< Names of the sampled variables.
   character(:), allocatable                 :: header    !< Variables header.
   integer(I4P)                              :: p, v      !< Counters.
   integer(I4P)                              :: nz        !< Digits of the file numbers.

   self%nv = size(q_name, kind=I4P)
   if (allocated(self%value)) deallocate(self%value)
   allocate(self%value(self%nv,self%probes_number))
   self%value = 0._R8P
   if (mpih%myrank /= 0) return
   call self%close_files
   if (allocated(self%unit)) deallocate(self%unit)
   allocate(self%unit(self%probes_number))
   ! strz keeps the last nz digits: at least 2, as many as the largest probe number
   nz = max(2_I4P, len_trim(str(self%probes_number, .true.), kind=I4P))
   header = '%VARIABLES="it" "time"'
   do v=1, self%nv
      header = header//' "'//q_name(v)%chars()//'"'
   enddo
   do p=1, self%probes_number
      open(newunit=self%unit(p), file=self%output_basename//'-probe_'//trim(strz(p, nz))//'.dat')
      write(self%unit(p), '(A)') '% probe '//trim(str(p, .true.))//' at x y z = '//trim(str(self%point_input(:,p)))
      write(self%unit(p), '(A)') header
      flush(self%unit(p))
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
   !< Reduce the local sums to rank 0 and write one history row per probe; collective.
   class(prism_probes_object), intent(inout) :: self !< Probes.
   integer(I4P),               intent(in)    :: it   !< Current iteration.
   real(R8P),                  intent(in)    :: time !< Current time.
   character(:), allocatable                 :: row      !< History row.
   real(R8P)                                 :: dummy(1) !< Receive buffer of the non-root ranks (not significant).
   integer(I4P)                              :: p, v     !< Counters.

   if (mpih%myrank == 0) then
      call MPI_REDUCE(MPI_IN_PLACE, self%value, self%nv*self%probes_number, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, &
                      mpih%error)
   else
      call MPI_REDUCE(self%value, dummy, self%nv*self%probes_number, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,        &
                      mpih%error)
      return
   endif
   do p=1, self%probes_number
      row = trim(str(it, .true.))//' '//trim(str(time))
      do v=1, self%nv
         row = row//' '//trim(str(self%value(v,p)))
      enddo
      write(self%unit(p), '(A)') row
      flush(self%unit(p))
   enddo
   endsubroutine save_history

   ! private non TBP
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
