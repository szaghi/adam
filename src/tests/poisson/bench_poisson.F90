!< Weak- and strong-scaling benchmark of FLAIL's solves (issue #84, E5).
program bench_poisson
!< Weak- and strong-scaling benchmark of FLAIL's solves (issue #84, E5).
!<
!< **Why this program exists.** Both FLAIL solvers gather to rank 0: the FFT solver gathers the whole uniform grid there
!< and transforms it on one core (or one GPU), and the multigrid gathers its coarsest level there in every V-cycle, and
!< the ghost charges of the James step in every isolated solve. E5 replaces those gathers; this program measures what
!< they cost, before and after, phase by phase.
!<
!< **What it times.** A Gaussian source on the cube [-1,1]^3 (octree, `ni`^3 cells per block), uniform at `--level`, or
!< with `--amr` the blocks inside the box [-1/2, 1/2]^3 refined once more (the multigrid path; level 2 or more, so that
!< some block lies inside). FLAIL is set up once (timed), solved once to warm up, then `--reps` times with the source
!< moved by `SHIFT` each time (the multigrid warm-starts from its last solution, as in a time step; `iterations` sums
!< the timed solves'); each solve is bracketed by barriers and timed on every rank, the maximum over the ranks kept.
!< The FFT solver's phase timers (`poisson_fft_object%wall`: gather, rank-0 solve, scatter) are reported per solve as
!< rank 0 measures them (the other ranks wait inside the collectives: their phases mix in rank 0's work), for the
!< uniform FFT path, the multigrid's coarsest level (`mg%fft`) and its James step (`mg%gfft`). Rank 0 prints one
!< `BENCH` line of `key=value` pairs (`bench.sh` collects them).
!<
!< Usage: `mpirun -np N exe/bench_poisson [--level L] [--ni N] [--reps R] [--isolated] [--amr]` (or
!< `exe/bench_poisson_fnl`: the FFT transforms on cuFFT, the multigrid on the device, #84 E5d and E6a).
#ifdef _FNL
#include "fundal.H"
#endif
use :: adam_common_library,         only : realm_object
use :: adam_flail_object,           only : flail_object, FLAIL_BACKEND_AUTO
use :: adam_parameters,             only : BC_PERIODIC, TO_BE_REFINED
use :: adam_tree_node_object,       only : tree_node_object
use :: adam_mpih_global,            only : mpih
#ifdef _FNL
use :: adam_fnl_mpih_global,        only : mpih_fnl
use :: adam_fnl_poisson_fft_object, only : poisson_fft_fnl_object
use :: adam_fnl_poisson_mg_object,  only : poisson_mg_fnl_object
#endif
use :: adam_poisson_fft_object,     only : poisson_fft_object
use :: mpi
use :: penf,                        only : I4P, I8P, R8P, str

implicit none

real(R8P), parameter            :: FACTOR=4._R8P * acos(-1._R8P) !< Source factor, 4 pi G with G = 1.
real(R8P), parameter            :: SIGMA2=0.05_R8P               !< Gaussian width squared.
real(R8P), parameter            :: SHIFT=0.01_R8P                !< Source shift along x per timed solve.
type(realm_object)              :: realm                         !< Realm.
type(flail_object)              :: flail                         !< FLAIL.
type(tree_node_object), pointer :: node                          !< A tree node.
real(R8P), allocatable          :: q(:,:,:,:,:)                  !< Source field (1, i, j, k, b).
real(R8P), allocatable          :: phi(:,:,:,:)                  !< Potential (i, j, k, b).
character(256)                  :: arg                           !< Command argument.
character(:), allocatable       :: filename                      !< Input file name.
character(:), allocatable       :: path                          !< Solver path label.
integer(I4P)                    :: level, ni, reps               !< Base level, cells per block side, timed solves.
logical                         :: isolated, amr                 !< Free-space boundaries; refine the central box.
real(R8P)                       :: t0, t_setup, t_solve          !< Wall times.
real(R8P)                       :: t_min, t_sum                  !< Solve statistics.
real(R8P)                       :: w(3,3)                        !< Phase walls per solve [phase, solver].
real(R8P)                       :: c(3)                          !< Block centre.
integer(I8P)                    :: cells, blocks                 !< Leaf cells and blocks over the ranks.
integer(I4P)                    :: its                           !< Iterations over the timed solves.
integer(I4P)                    :: ngc, nk, a, b, d, ierr        !< Counters.

level = 3 ; ni = 8 ; reps = 5 ; isolated = .false. ; amr = .false.
do a=1, command_argument_count()
   call get_command_argument(a, arg)
   select case(trim(arg))
   case('--level')
      call get_command_argument(a+1, arg) ; read(arg, *) level
   case('--ni')
      call get_command_argument(a+1, arg) ; read(arg, *) ni
   case('--reps')
      call get_command_argument(a+1, arg) ; read(arg, *) reps
   case('--isolated')
      isolated = .true.
   case('--amr')
      amr = .true.
   endselect
enddo
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
if (amr .and. level < 2) call mpih%error_stop(msg='bench_poisson: --amr needs --level 2 or more')
#ifdef _FNL
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
#endif
filename = 'bench_poisson-l'//trim(str(level, .true.))//'-n'//trim(str(ni, .true.))//'.ini'
if (mpih%myrank == 0) call write_input(filename=filename, ni=ni, max_level=level+1)
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
! room for every block on one rank, the refinement of the central box included (8^L uniform, 15/8 of it refined)
call realm%initialize(filename=filename, nv=1_I4P, fields_number=1_I4P, &
                      memory_avail=2.5_R8P * 8._R8P**level * real(ni + 6, R8P)**3 * 8._R8P / 0.85_R8P / 1.e9_R8P)
if (.not.isolated) call realm%adam%grid%set_bc_type(bc_type=[(BC_PERIODIC, d=1, 6)])
ngc = realm%adam%grid%ngc ; nk = realm%adam%grid%nk
allocate(q(1:1,1-ngc:ni+ngc,1-ngc:ni+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
q = 0._R8P
if (level > 0) call realm%adam%refine_uniform(refinement_levels=level, q=q, do_mpi_redistribute=.true., &
                                              do_blocks_reorder=.false.)
if (amr) then
   do b=1, realm%adam%field%blocks_number
      c = 0.5_R8P * (realm%adam%field%emin(:,b) + realm%adam%field%emax(:,b))
      if (maxval(abs(c)) < 0.5_R8P) then
         node => realm%adam%tree%node(code=realm%adam%field%code(b))
         node%refinement_needed = TO_BE_REFINED
      endif
   enddo
   call realm%adam%amr_update(q=q, is_marked_by_tree=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
endif
call set_source(shift=0._R8P)
allocate(phi(1-ngc:ni+ngc,1-ngc:ni+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
phi = 0._R8P
blocks = int(realm%adam%field%blocks_number, I8P)
call MPI_ALLREDUCE(MPI_IN_PLACE, blocks, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
cells = blocks * int(ni, I8P)**2 * int(nk, I8P)

#ifdef _FNL
allocate(poisson_fft_fnl_object :: flail%fft)
allocate(poisson_mg_fnl_object :: flail%mg_mold)
#endif
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
t0 = MPI_WTIME()
call flail%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=isolated, backend=FLAIL_BACKEND_AUTO)
t_setup = MPI_WTIME() - t0
call MPI_ALLREDUCE(MPI_IN_PLACE, t_setup, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
call flail%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi) ! warm-up
call reset_walls
t_min = huge(1._R8P) ; t_sum = 0._R8P ; its = 0_I4P
do a=1, reps
   call set_source(shift=SHIFT * a) ! as a time step moves the source: the multigrid warm-starts from the last solve
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   t0 = MPI_WTIME()
   call flail%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   t_solve = MPI_WTIME() - t0
   call MPI_ALLREDUCE(MPI_IN_PLACE, t_solve, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   t_min = min(t_min, t_solve) ; t_sum = t_sum + t_solve
   its = its + flail%iterations()
enddo
w = 0._R8P
if (flail%is_multigrid()) then
   path = 'mg'
   w(:,2) = flail%mg%fft%wall / reps
   if (allocated(flail%mg%gfft)) w(:,3) = flail%mg%gfft%wall / reps
else
   path = 'fft'
   w(:,1) = flail%fft%wall / reps
endif
if (mpih%myrank == 0) print '(A)', 'BENCH path='//path//' isolated='//merge('T', 'F', isolated)//' amr='//merge('T', 'F', amr)// &
   ' np='//trim(str(mpih%procs_number, .true.))//' level='//trim(str(level, .true.))//' ni='//trim(str(ni, .true.))// &
   ' blocks='//trim(str(blocks, .true.))//' cells='//trim(str(cells, .true.))// &
   ' iterations='//trim(str(its, .true.))//' setup='//trim(str(t_setup, .true.))// &
   ' solve_min='//trim(str(t_min, .true.))//' solve_mean='//trim(str(t_sum / reps, .true.))// &
   ' fft_gather='//trim(str(w(1,1), .true.))//' fft_rank0='//trim(str(w(2,1), .true.))// &
   ' fft_scatter='//trim(str(w(3,1), .true.))//' coarse_gather='//trim(str(w(1,2), .true.))// &
   ' coarse_rank0='//trim(str(w(2,2), .true.))//' coarse_scatter='//trim(str(w(3,2), .true.))// &
   ' james_gather='//trim(str(w(1,3), .true.))//' james_rank0='//trim(str(w(2,3), .true.))// &
   ' james_scatter='//trim(str(w(3,3), .true.))
call flail%destroy
call mpih%finalize

contains
   subroutine set_source(shift)
   !< The Gaussian source centred at (shift, 0, 0).
   real(R8P), intent(in) :: shift   !< Centre along x.
   real(R8P)             :: r2      !< Radius squared.
   integer(I4P)          :: b, i, j, k !< Counters.

   do b=1, realm%adam%field%blocks_number
      do k=1, nk
         do j=1, ni
            do i=1, ni
               r2 = (realm%adam%field%x_cell(i,b) - shift)**2 + realm%adam%field%y_cell(j,b)**2 + &
                    realm%adam%field%z_cell(k,b)**2
               q(1,i,j,k,b) = exp(-r2 / SIGMA2)
            enddo
         enddo
      enddo
   enddo
   endsubroutine set_source

   subroutine reset_walls
   !< Zero the FFT solvers' phase timers (after the warm-up solve).

   if (allocated(flail%fft)) flail%fft%wall = 0._R8P
   if (flail%is_multigrid()) then
      flail%mg%fft%wall = 0._R8P
      if (allocated(flail%mg%gfft)) flail%mg%gfft%wall = 0._R8P
   endif
   endsubroutine reset_walls

   subroutine write_input(filename, ni, max_level)
   !< Write the sections realm_object%initialize reads: the cube [-1,1]^3, one root block of `ni`^3 cells.
   character(*), intent(in) :: filename  !< Input file name.
   integer(I4P), intent(in) :: ni        !< Cells per block side.
   integer(I4P), intent(in) :: max_level !< Deepest level.
   integer(I4P)             :: u         !< File unit.

   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = bench_poisson'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = bench_poisson-restart'
   write(u, '(A)') 'restart_save           = 999999'
   write(u, '(A)') 'residuals_save         = 999999'
   write(u, '(A)') 'save_memory_status     = .false.'
   write(u, '(A)') 'save_residual_fields   = .false.'
   write(u, '(A)') 'save_curl_fields       = .false.'
   write(u, '(A)') 'save_divergence_fields = .false.'
   write(u, '(A)') 'save_gradient_fields   = .false.'
   write(u, '(A)') 'save_laplacian_fields  = .false.'
   write(u, '(A)') '[grid]'
   write(u, '(A)') 'ni     = '//trim(str(ni, .true.))
   write(u, '(A)') 'nj     = '//trim(str(ni, .true.))
   write(u, '(A)') 'nk     = '//trim(str(ni, .true.))
   write(u, '(A)') 'ngc    = 3'
   write(u, '(A)') 'emin_x = -1.0'
   write(u, '(A)') 'emin_y = -1.0'
   write(u, '(A)') 'emin_z = -1.0'
   write(u, '(A)') 'emax_x = 1.0'
   write(u, '(A)') 'emax_y = 1.0'
   write(u, '(A)') 'emax_z = 1.0'
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = .false.'
   write(u, '(A)') 'null_z = .false.'
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level      = '//trim(str(max_level, .true.))
   write(u, '(A)') 'ratio          = 8'
   write(u, '(A)') 'iu_ref_levels  = 0'
   write(u, '(A)') 'i_prune        = 0'
   write(u, '(A)') 'j_prune        = 0'
   write(u, '(A)') 'k_prune        = 0'
   write(u, '(A)') 'l_prune        = -1'
   write(u, '(A)') 'frequency      = 999999'
   write(u, '(A)') 'iters          = 1'
   write(u, '(A)') 'markers_number = 0'
   write(u, '(A)') '[field]'
   write(u, '(A)') 'nv = 1'
   write(u, '(A)') '[runge_kutta]'
   write(u, '(A)') 'scheme = runge-kutta-ssp-33'
   write(u, '(A)') '[weno]'
   write(u, '(A)') 'scheme              = weno-u-5'
   write(u, '(A)') 'ror_number          = 0'
   write(u, '(A)') 'ror_threshold       = 0.9'
   write(u, '(A)') 'ror_vars_number     = 0'
   write(u, '(A)') 'enable_ror_stats    = .false.'
   write(u, '(A)') 'ib_reduction_extent = 0'
   write(u, '(A)') 'ib_reduced_order    = 2'
   write(u, '(A)') '[fdv]'
   write(u, '(A)') 'fdv_scheme = fd'
   write(u, '(A)') 'fdv_order  = 2'
   write(u, '(A)') '[solids]'
   write(u, '(A)') 'number    = 0'
   write(u, '(A)') 'n_eikonal = 0'
   write(u, '(A)') '[slices]'
   write(u, '(A)') 'slices_number = 0'
   close(u)
   endsubroutine write_input
endprogram bench_poisson
