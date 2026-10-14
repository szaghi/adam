!< ADAM, FFT Poisson solver on uniform periodic grids.
module adam_poisson_fft_object
!< ADAM, FFT Poisson solver on uniform periodic grids.
!<
!< Solves `L_h phi = factor (f - mean(f))` for a cell field `f` (one variable of a conservative array `q`) on a grid whose
!< blocks all sit at one level and whose active directions are all periodic. `L_h` is the second-order 7-point Laplacian
!< (5-point in 2-D; a null direction contributes nothing). The solve divides each Fourier mode by the eigenvalue of
!< `L_h` itself, `lambda_k = sum_d (2 cos(2 pi k_d / N_d) - 2) / dx_d^2`, so the discrete residual of the returned `phi`
!< is at round-off; the modes with `lambda_k = 0` (the mean, and every mode along a null direction) are set to zero.
!<
!< **Data path (issue #81, D-81-2, first version).** Every rank packs the interior `f` of its blocks with their integer
!< block coordinates; `MPI_Gatherv` brings them to rank 0, which places them by the global cell index, transforms
!< (`transform`, real-to-complex and back), and `MPI_Scatterv`s every block's `phi` with all its ghost layers already
!< filled by the periodic wrap of the global solution, so no ghost exchange follows. Correct for any rank count, ranks
!< without blocks included; not scalable (rank 0 holds the whole grid).
!<
!< **Backends.** This object transforms with FFTW 3 (`-D_FFTW`, `-lfftw3`; FFTW_ESTIMATE plans, deterministic). The FNL
!< twin `adam_fnl_poisson_fft_object` extends it and overrides `make_plans`, `transform` and `destroy_plans` with cuFFT
!< on the GPU of rank 0. Built without either, `initialize` stops: "built without an FFT library".
!<
!< **Entry points.** `solve` packs one variable of a host conservative array; `solve_packed` takes the interior source
!< already packed (block-major, `i` fastest: what a device kernel produces on FNL). Both can return
!< `energy = sum_cells f Phi` (computed on rank 0 from the global arrays, broadcast): `W = energy dV / 2` is the
!< gravitational energy of self-gravity on a uniform grid.
!<
!< **Scope (fatal otherwise).** One level for every block, the blocks covering the domain, every active direction
!< periodic. Isolated boundaries (zero padding) and AMR (multigrid) are later phases of #81.
use :: adam_grid_object, only : grid_object
use :: adam_field_object, only : field_object
use :: adam_mpih_global, only : mpih
use :: iso_c_binding
use :: mpi
use :: penf, only : I4P, R8P, str

implicit none
private
public :: poisson_fft_object

real(R8P), parameter :: PI=acos(-1._R8P) !< Pi.

type :: poisson_fft_object
   !< FFT Poisson solver on uniform periodic grids.
   logical                   :: is_initialized=.false. !< Initialization status.
   integer(I4P)              :: n(3)=1_I4P             !< Global cells per direction.
   integer(I4P)              :: nijk(3)=1_I4P          !< Cells per block per direction.
   integer(I4P)              :: ngc=0_I4P              !< Ghost cells per side.
   integer(I4P)              :: level=0_I4P            !< Level of every block.
   real(R8P)                 :: dxyz(3)=0._R8P         !< Cell spacing per direction.
   real(R8P),    allocatable :: lam1(:)                !< Eigenvalue part along x [0:n(1)/2].
   real(R8P),    allocatable :: lam2(:)                !< Eigenvalue part along y [0:n(2)-1].
   real(R8P),    allocatable :: lam3(:)                !< Eigenvalue part along z [0:n(3)-1].
   real(R8P),    allocatable :: f_global(:,:,:)        !< Rank 0: the source, then the potential [n(1),n(2),n(3)].
   real(R8P),    allocatable :: f_source(:,:,:)        !< Rank 0: a copy of the source, for the energy.
   complex(R8P), allocatable :: f_hat(:,:,:)           !< Rank 0: the transform [n(1)/2+1,n(2),n(3)].
   type(c_ptr)               :: plan_r2c=c_null_ptr    !< FFTW forward plan (rank 0).
   type(c_ptr)               :: plan_c2r=c_null_ptr    !< FFTW backward plan (rank 0).
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy the solver.
      procedure, pass(self) :: initialize !< Check the grid, size the solver, make the plans.
      procedure, pass(self) :: solve        !< Solve for phi from one variable of q.
      procedure, pass(self) :: solve_packed !< Solve for phi from the packed interior source.
      ! backend hooks
      procedure, pass(self) :: destroy_plans !< Destroy the transform plans.
      procedure, pass(self) :: make_plans    !< Make the transform plans (rank 0).
      procedure, pass(self) :: transform     !< f_global <- factor * L_h^-1 f_global, mean removed (rank 0).
endtype poisson_fft_object

#ifdef _FFTW
integer(c_int), parameter :: FFTW_ESTIMATE=64_c_int !< FFTW planner flag, `fftw3.f03`.

interface
   type(c_ptr) function fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, flags) bind(C, name='fftw_plan_dft_r2c_3d')
   !< FFTW: plan a 3-D real-to-complex transform (row-major sizes).
   import :: c_ptr, c_int, c_double, c_double_complex
   integer(c_int), value     :: n0, n1, n2 !< Sizes, slowest first.
   real(c_double)            :: in(*)      !< Input.
   complex(c_double_complex) :: out(*)     !< Output.
   integer(c_int), value     :: flags      !< Planner flags.
   endfunction fftw_plan_dft_r2c_3d

   type(c_ptr) function fftw_plan_dft_c2r_3d(n0, n1, n2, in, out, flags) bind(C, name='fftw_plan_dft_c2r_3d')
   !< FFTW: plan a 3-D complex-to-real transform (row-major sizes).
   import :: c_ptr, c_int, c_double, c_double_complex
   integer(c_int), value     :: n0, n1, n2 !< Sizes, slowest first.
   complex(c_double_complex) :: in(*)      !< Input.
   real(c_double)            :: out(*)     !< Output.
   integer(c_int), value     :: flags      !< Planner flags.
   endfunction fftw_plan_dft_c2r_3d

   subroutine fftw_execute_dft_r2c(p, in, out) bind(C, name='fftw_execute_dft_r2c')
   !< FFTW: execute a real-to-complex plan on the given arrays.
   import :: c_ptr, c_double, c_double_complex
   type(c_ptr), value        :: p      !< Plan.
   real(c_double)            :: in(*)  !< Input.
   complex(c_double_complex) :: out(*) !< Output.
   endsubroutine fftw_execute_dft_r2c

   subroutine fftw_execute_dft_c2r(p, in, out) bind(C, name='fftw_execute_dft_c2r')
   !< FFTW: execute a complex-to-real plan on the given arrays.
   import :: c_ptr, c_double, c_double_complex
   type(c_ptr), value        :: p      !< Plan.
   complex(c_double_complex) :: in(*)  !< Input.
   real(c_double)            :: out(*) !< Output.
   endsubroutine fftw_execute_dft_c2r

   subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
   !< FFTW: destroy a plan.
   import :: c_ptr
   type(c_ptr), value :: p !< Plan.
   endsubroutine fftw_destroy_plan
endinterface
#endif

contains
   ! public methods
   subroutine destroy(self)
   !< Destroy the solver.
   class(poisson_fft_object), intent(inout) :: self !< The solver.

   if (self%is_initialized) call self%destroy_plans
   if (allocated(self%lam1))     deallocate(self%lam1)
   if (allocated(self%lam2))     deallocate(self%lam2)
   if (allocated(self%lam3))     deallocate(self%lam3)
   if (allocated(self%f_global)) deallocate(self%f_global)
   if (allocated(self%f_source)) deallocate(self%f_source)
   if (allocated(self%f_hat))    deallocate(self%f_hat)
   self%n = 1_I4P ; self%nijk = 1_I4P ; self%ngc = 0_I4P ; self%level = 0_I4P ; self%dxyz = 0._R8P
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine initialize(self, grid, field)
   !< Check the grid (one level, the domain covered, active directions periodic), size the solver, make the plans.
   class(poisson_fft_object), intent(inout) :: self         !< The solver.
   type(grid_object),         intent(in)    :: grid         !< The grid.
   type(field_object),        intent(in)    :: field        !< The field (block coordinates).
   integer(I4P)                             :: lmin, lmax   !< Min and max block level over the ranks.
   integer(I4P)                             :: nblk(3)      !< Blocks per direction at the level.
   integer(I4P)                             :: blocks_total !< Blocks over the ranks.
   integer(I4P)                             :: d, k         !< Counters.
   integer(I4P)                             :: ierr         !< MPI error status.

   call self%destroy
   self%nijk = [grid%ni, grid%nj, grid%nk]
   self%ngc  = grid%ngc
   lmin = huge(1_I4P) ; lmax = -huge(1_I4P)
   if (field%blocks_number > 0) then
      lmin = minval(field%coordinates(4,1:field%blocks_number))
      lmax = maxval(field%coordinates(4,1:field%blocks_number))
   endif
   call MPI_ALLREDUCE(MPI_IN_PLACE, lmin, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, lmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
   if (lmin /= lmax) call mpih%error_stop(msg='poisson_fft_object%initialize: the blocks sit at levels '// &
      trim(str(lmin, .true.))//' to '//trim(str(lmax, .true.))//'; the FFT solver needs a uniform grid (one level)')
   self%level = lmin
   self%dxyz  = grid%cell_dxyz(:,self%level)
   nblk = nint((grid%domain_emax - grid%domain_emin) / grid%block_dxyz(:,self%level))
   blocks_total = field%blocks_number
   call MPI_ALLREDUCE(MPI_IN_PLACE, blocks_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   if (blocks_total /= product(nblk)) call mpih%error_stop(msg='poisson_fft_object%initialize: '// &
      trim(str(blocks_total, .true.))//' blocks do not cover the domain ('//trim(str(product(nblk), .true.))// &
      ' at level '//trim(str(self%level, .true.))//')')
   do d=1, 3
      if (.not.grid%null_xyz(d) .and. .not.grid%is_ijk_periodic(d)) call mpih%error_stop(msg= &
         'poisson_fft_object%initialize: direction '//trim(str(d, .true.))//' is not periodic; the FFT solver '// &
         'supports periodic boundaries only')
   enddo
   self%n = nblk * self%nijk
   ! eigenvalue parts of L_h per direction; zero along a null direction (no derivative there)
   allocate(self%lam1(0:self%n(1)/2), self%lam2(0:self%n(2)-1), self%lam3(0:self%n(3)-1))
   self%lam1 = 0._R8P ; self%lam2 = 0._R8P ; self%lam3 = 0._R8P
   if (.not.grid%null_xyz(1)) then
      do k=0, self%n(1)/2
         self%lam1(k) = (2._R8P * cos(2._R8P * PI * k / self%n(1)) - 2._R8P) / self%dxyz(1)**2
      enddo
   endif
   if (.not.grid%null_xyz(2)) then
      do k=0, self%n(2)-1
         self%lam2(k) = (2._R8P * cos(2._R8P * PI * k / self%n(2)) - 2._R8P) / self%dxyz(2)**2
      enddo
   endif
   if (.not.grid%null_xyz(3)) then
      do k=0, self%n(3)-1
         self%lam3(k) = (2._R8P * cos(2._R8P * PI * k / self%n(3)) - 2._R8P) / self%dxyz(3)**2
      enddo
   endif
   if (mpih%myrank == 0) then
      allocate(self%f_global(self%n(1),self%n(2),self%n(3)))
      allocate(self%f_source(self%n(1),self%n(2),self%n(3)))
      allocate(self%f_hat(self%n(1)/2+1,self%n(2),self%n(3)))
      self%f_global = 0._R8P
      self%f_source = 0._R8P
      self%f_hat = (0._R8P, 0._R8P)
      call self%make_plans
   endif
   self%is_initialized = .true.
   endsubroutine initialize

   subroutine solve(self, field, q, iv, factor, phi, energy)
   !< Solve `L_h phi = factor (q(iv) - mean(q(iv)))`; every ghost layer of phi is filled by the periodic wrap.
   class(poisson_fft_object), intent(inout)         :: self                       !< The solver.
   type(field_object),        intent(in)            :: field                      !< The field (coordinates).
   real(R8P),    intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)    !< Conservative [nv,i,j,k,nb].
   integer(I4P), intent(in)    :: iv                                              !< Variable of q, the source.
   real(R8P),    intent(in)    :: factor                                          !< Source factor (4 pi G).
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)     !< Potential [i,j,k,nb].
   real(R8P),    intent(out), optional :: energy                                  !< sum_cells f phi (all ranks).
   real(R8P),    allocatable                        :: buf(:)                     !< Packed interior source.
   integer(I4P)                                     :: b, i, j, k, c              !< Counters.

   allocate(buf(max(field%blocks_number*product(self%nijk),1)))
   c = 0
   do b=1, field%blocks_number
      do k=1, self%nijk(3)
         do j=1, self%nijk(2)
            do i=1, self%nijk(1)
               c = c + 1
               buf(c) = q(iv,i,j,k,b)
            enddo
         enddo
      enddo
   enddo
   call self%solve_packed(field=field, buf=buf, factor=factor, phi=phi, energy=energy)
   endsubroutine solve

   subroutine solve_packed(self, field, buf, factor, phi, energy)
   !< Solve `L_h phi = factor (f - mean(f))` from the packed interior source of this rank's blocks (block-major, `i`
   !< fastest, `field%blocks_number * ni * nj * nk` values); every ghost layer of phi is filled by the periodic wrap.
   class(poisson_fft_object), intent(inout)         :: self            !< The solver.
   type(field_object),        intent(in)            :: field           !< The field (block coordinates).
   real(R8P),                 intent(in)            :: buf(1:)         !< Packed interior source of this rank.
   real(R8P),                 intent(in)            :: factor          !< Source factor (4 pi G for gravity).
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Potential [i,j,k,nb].
   real(R8P),                 intent(out), optional :: energy          !< sum_cells f phi (all ranks).
   integer(I4P)                                     :: nbl             !< Blocks of this rank.
   integer(I4P)                                     :: ncell           !< Interior cells per block.
   integer(I4P)                                     :: nghost          !< Cells per block, ghosts included.
   integer(I4P)                                     :: nr              !< Ranks.
   integer(I4P)                                     :: blocks_total    !< Rank 0: blocks over the ranks.
   integer(I4P), allocatable                        :: counts(:)       !< Rank 0: blocks per rank.
   integer(I4P), allocatable                        :: displs(:)       !< Rank 0: displacement per rank (blocks).
   integer(I4P), allocatable                        :: rcounts(:)      !< Rank 0: counts per rank (values).
   integer(I4P), allocatable                        :: rdispls(:)      !< Rank 0: displacements per rank (values).
   integer(I4P), allocatable                        :: coords(:,:)     !< Block coordinates of this rank [3,nbl].
   integer(I4P), allocatable                        :: coords_all(:,:) !< Rank 0: block coordinates of every rank.
   real(R8P),    allocatable                        :: rbuf(:)         !< Receive buffer of this rank.
   real(R8P),    allocatable                        :: buf_all(:)      !< Rank 0: buffer of every rank.
   real(R8P)                                        :: e               !< sum_cells f phi.
   integer(I4P)                                     :: ub(3)           !< Upper bounds of phi per direction.
   integer(I4P)                                     :: g(3)            !< Global cell index, wrapped.
   integer(I4P)                                     :: b, i, j, k, c, r !< Counters.
   integer(I4P)                                     :: ierr            !< MPI error status.

   if (.not.self%is_initialized) call mpih%error_stop(msg='poisson_fft_object%solve: called before initialize')
   ub = self%nijk + self%ngc
   if (ubound(phi, dim=1) /= ub(1) .or. ubound(phi, dim=2) /= ub(2) .or. ubound(phi, dim=3) /= ub(3) .or. &
       size(phi, dim=4) < field%blocks_number) &
      call mpih%error_stop(msg='poisson_fft_object%solve: phi is not shaped [1-ngc:n+ngc per direction, blocks]')
   nbl    = field%blocks_number
   ncell  = product(self%nijk)
   nghost = product(self%nijk + 2 * self%ngc)
   nr     = mpih%procs_number
   if (size(buf) < nbl*ncell) call mpih%error_stop(msg='poisson_fft_object%solve: the packed source is too short')
   allocate(coords(3,max(nbl,1)))
   do b=1, nbl
      if (field%coordinates(4,b) /= self%level) call mpih%error_stop(msg='poisson_fft_object%solve: a block left '// &
         'level '//trim(str(self%level, .true.))//' after initialize (regrids are not supported)')
      coords(:,b) = field%coordinates(1:3,b)
   enddo
   ! gather to rank 0
   allocate(counts(nr), displs(nr), rcounts(nr), rdispls(nr))
   counts = 0 ; displs = 0 ; rcounts = 0 ; rdispls = 0
   call MPI_GATHER(nbl, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   blocks_total = 0
   if (mpih%myrank == 0) then
      do r=2, nr
         displs(r) = displs(r-1) + counts(r-1)
      enddo
      blocks_total = sum(counts)
      rcounts = 3 * counts ; rdispls = 3 * displs
   endif
   allocate(coords_all(3,max(blocks_total,1)), buf_all(max(blocks_total*nghost,1)))
   call MPI_GATHERV(coords, 3*nbl, MPI_INTEGER, coords_all, rcounts, rdispls, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (mpih%myrank == 0) then
      rcounts = ncell * counts ; rdispls = ncell * displs
   endif
   call MPI_GATHERV(buf, nbl*ncell, MPI_REAL8, buf_all, rcounts, rdispls, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   ! rank 0: place, transform, pack every block with its ghosts by the periodic wrap
   e = 0._R8P
   if (mpih%myrank == 0) then
      c = 0
      do b=1, blocks_total
         do k=1, self%nijk(3)
            do j=1, self%nijk(2)
               do i=1, self%nijk(1)
                  c = c + 1
                  self%f_global(coords_all(1,b)*self%nijk(1)+i, coords_all(2,b)*self%nijk(2)+j, &
                                coords_all(3,b)*self%nijk(3)+k) = buf_all(c)
               enddo
            enddo
         enddo
      enddo
      if (present(energy)) self%f_source = self%f_global
      call self%transform(factor=factor)
      if (present(energy)) e = sum(self%f_source * self%f_global)
      c = 0
      do b=1, blocks_total
         do k=1-self%ngc, ub(3)
            g(3) = modulo(coords_all(3,b)*self%nijk(3) + k - 1, self%n(3)) + 1
            do j=1-self%ngc, ub(2)
               g(2) = modulo(coords_all(2,b)*self%nijk(2) + j - 1, self%n(2)) + 1
               do i=1-self%ngc, ub(1)
                  g(1) = modulo(coords_all(1,b)*self%nijk(1) + i - 1, self%n(1)) + 1
                  c = c + 1
                  buf_all(c) = self%f_global(g(1),g(2),g(3))
               enddo
            enddo
         enddo
      enddo
      rcounts = nghost * counts ; rdispls = nghost * displs
   endif
   ! scatter back, ghosts included
   allocate(rbuf(max(nbl*nghost,1)))
   call MPI_SCATTERV(buf_all, rcounts, rdispls, MPI_REAL8, rbuf, nbl*nghost, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   c = 0
   do b=1, nbl
      do k=1-self%ngc, ub(3)
         do j=1-self%ngc, ub(2)
            do i=1-self%ngc, ub(1)
               c = c + 1
               phi(i,j,k,b) = rbuf(c)
            enddo
         enddo
      enddo
   enddo
   if (present(energy)) then
      call MPI_BCAST(e, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      energy = e
   endif
   endsubroutine solve_packed

   ! backend hooks
   subroutine destroy_plans(self)
   !< Destroy the FFTW plans.
   class(poisson_fft_object), intent(inout) :: self !< The solver.

#ifdef _FFTW
   if (c_associated(self%plan_r2c)) call fftw_destroy_plan(self%plan_r2c)
   if (c_associated(self%plan_c2r)) call fftw_destroy_plan(self%plan_c2r)
#endif
   self%plan_r2c = c_null_ptr
   self%plan_c2r = c_null_ptr
   endsubroutine destroy_plans

   subroutine make_plans(self)
   !< Make the FFTW plans on rank 0's global arrays (FFTW_ESTIMATE: the arrays are not touched, the plans deterministic).
   class(poisson_fft_object), intent(inout) :: self !< The solver.

#ifdef _FFTW
   ! FFTW is row-major: the slowest dimension first
   self%plan_r2c = fftw_plan_dft_r2c_3d(int(self%n(3), c_int), int(self%n(2), c_int), int(self%n(1), c_int), &
                                        self%f_global, self%f_hat, FFTW_ESTIMATE)
   self%plan_c2r = fftw_plan_dft_c2r_3d(int(self%n(3), c_int), int(self%n(2), c_int), int(self%n(1), c_int), &
                                        self%f_hat, self%f_global, FFTW_ESTIMATE)
   if (.not.c_associated(self%plan_r2c) .or. .not.c_associated(self%plan_c2r)) &
      call mpih%error_stop(msg='poisson_fft_object%make_plans: FFTW could not make the plans')
#else
   call mpih%error_stop(msg='poisson_fft_object: built without an FFT library (build a *-fft mode: -D_FFTW on the '// &
                            'CPU, -D_CUFFT on FNL)')
#endif
   endsubroutine make_plans

   subroutine transform(self, factor)
   !< Rank 0: `f_global <- factor * L_h^-1 f_global`, the modes of zero eigenvalue (the mean) set to zero.
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   real(R8P),                 intent(in)    :: factor  !< Source factor.
   real(R8P)                                :: scale   !< Factor over the unnormalised transform pair.
   real(R8P)                                :: lam     !< Eigenvalue of the mode.
   integer(I4P)                             :: i, j, k !< Counters.

#ifdef _FFTW
   call fftw_execute_dft_r2c(self%plan_r2c, self%f_global, self%f_hat)
   scale = factor / real(product(self%n), R8P)
   do k=1, self%n(3)
      do j=1, self%n(2)
         do i=1, self%n(1)/2+1
            lam = self%lam1(i-1) + self%lam2(j-1) + self%lam3(k-1)
            if (lam == 0._R8P) then
               self%f_hat(i,j,k) = (0._R8P, 0._R8P)
            else
               self%f_hat(i,j,k) = self%f_hat(i,j,k) * (scale / lam)
            endif
         enddo
      enddo
   enddo
   call fftw_execute_dft_c2r(self%plan_c2r, self%f_hat, self%f_global)
#else
   call mpih%error_stop(msg='poisson_fft_object%transform: built without an FFT library')
#endif
   endsubroutine transform
endmodule adam_poisson_fft_object
