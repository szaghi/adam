!< Unit test of the FFT Poisson solver on uniform periodic grids (issue #81, P2, SV-2).
program test_poisson_fft
!< Unit test of the FFT Poisson solver on uniform periodic grids (issue #81, P2, SV-2).
!<
!< **Why this test exists.** Self-gravity conserves E + rho Phi and momentum exactly only if Phi solves the discrete
!< Poisson equation `L_h Phi = 4 pi G (rho - mean rho)` to round-off (#81 P0, Mullen et al. 2021). The FFT solver
!< divides each Fourier mode by the eigenvalue of `L_h`; a wrong eigenvalue (a missing factor, a swapped axis, a wrong
!< spacing), a misplaced block in the gather, or a wrong ghost wrap in the scatter each break that.
!<
!< **What it pins.** The source is `rho = rho0 + A1 sin(k1 x) cos(k1 y) cos(k1 z) + A2 cos(k2 . x)`, two oblique modes
!< on the anisotropic periodic box [0,1] x [0,2] x [0,3]. Each mode is an eigenfunction of `L_h`, so the discrete
!< solution is known in closed form (each mode divided by its discrete eigenvalue, the mean dropped), and so is the
!< continuum one (divided by `-|k|^2`). Asserted, per configuration and level:
!<   - discrete: `max |Phi - Phi_h| / max |Phi_h|` over every cell, ghost layers included, at round-off (`TOL_DISCRETE`):
!<     the eigenvalues, the gather placement and the ghost wrap;
!<   - residual: `max |L_h Phi - 4 pi G (rho - rho0)| / max |4 pi G (rho - rho0)|` over the interior, from the
!<     scattered ghosts (`TOL_RESIDUAL`);
!<   - continuum: the max error against the exact Phi, second order on the finest pair (`MIN_ORDER`).
!< Configurations: an octree in 3-D (levels 1-3, 16^3 to 64^3 cells), a quadtree in 2-D (`null_z`, levels 1-3) and a
!< single block (level 0: the wrap inside one block; on more than one rank, ranks without blocks).
!<
!< The FNL build (`-D_FNL -D_CUFFT`) runs the spectral solve with cuFFT on rank 0's GPU and is held to the same bounds,
!< so CPU and FNL agree to round-off through the discrete solution.
!<
!< `--no-periodic` builds one configuration without periodic boundaries: the solver must stop (exit status 1).
!< `--dump DIR` writes rank 0's global potential of every level to `DIR/phi-r<ratio>-nk<nk>-l<level>.bin` (3 int32
!< sizes, then float64 values, x fastest): `selfgravity/poisson_oracle.py` compares two builds and plots them.
!<
!< Usage: `mpirun -np N exe/test_poisson_fft [--no-periodic] [--dump DIR]` (or `exe/test_poisson_fft_fnl`).

#ifdef _FNL
#include "fundal.H"
#endif
use :: adam_common_library,         only : realm_object
use :: adam_parameters,             only : BC_PERIODIC
use :: adam_mpih_global,            only : mpih
#ifdef _FNL
use :: adam_fnl_mpih_global,        only : mpih_fnl
use :: adam_fnl_poisson_fft_object, only : poisson_fft_fnl_object
#else
use :: adam_poisson_fft_object,     only : poisson_fft_object
#endif
use :: mpi
use :: penf,                        only : I4P, R8P, str

implicit none

real(R8P), parameter :: PI=acos(-1._R8P)              !< Pi.
real(R8P), parameter :: FACTOR=4._R8P * PI            !< Source factor, 4 pi G with G = 1.
real(R8P), parameter :: RHO0=1.5_R8P                  !< Mean density (dropped by the solve).
real(R8P), parameter :: A1=1._R8P                     !< Amplitude of the first mode.
real(R8P), parameter :: A2=0.5_R8P                    !< Amplitude of the second mode.
real(R8P), parameter :: L(3)=[1._R8P, 2._R8P, 3._R8P] !< Domain lengths.
real(R8P), parameter :: SENTINEL=-1.e30_R8P           !< Value of a phi cell the scatter did not fill.
real(R8P), parameter :: TOL_DISCRETE=1.e-12_R8P       !< Bound on the error against the discrete solution, relative.
real(R8P), parameter :: TOL_RESIDUAL=1.e-10_R8P       !< Bound on the discrete residual, relative.
real(R8P), parameter :: MIN_ORDER=1.9_R8P             !< Bound on the order against the continuum solution.
logical              :: test_passed                   !< Aggregate pass flag.
logical              :: no_periodic                   !< Run the non-periodic refusal only.
character(256)       :: arg                           !< Command argument.
character(256)       :: dump_dir                      !< Directory of the potential dumps ('': no dump).
integer(I4P)         :: a                             !< Counter.

no_periodic = .false.
dump_dir = ''
do a=1, command_argument_count()
   call get_command_argument(a, arg)
   if (trim(arg) == '--no-periodic') no_periodic = .true.
   if (trim(arg) == '--dump') call get_command_argument(a+1, dump_dir)
enddo
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
#ifdef _FNL
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
#endif
test_passed = .true.
if (no_periodic) then
   call check_ladder(label='3-D octree, not periodic', ratio=8_I4P, nk=8_I4P, levels=[1_I4P], periodic=.false.)
   if (mpih%myrank == 0) print '(A)', 'test_poisson_fft: FAILED (the solver accepted a non-periodic grid)'
   call mpih%finalize
   error stop 1
endif
call check_ladder(label='3-D octree,    ', ratio=8_I4P, nk=8_I4P, levels=[1_I4P, 2_I4P, 3_I4P], periodic=.true.)
call check_ladder(label='2-D quadtree,  ', ratio=4_I4P, nk=1_I4P, levels=[1_I4P, 2_I4P, 3_I4P], periodic=.true.)
call check_ladder(label='3-D one block, ', ratio=8_I4P, nk=8_I4P, levels=[0_I4P], periodic=.true.)
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'test_poisson_fft: PASSED'
   else
      print '(A)', 'test_poisson_fft: FAILED'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   subroutine check_ladder(label, ratio, nk, levels, periodic)
   !< Solve on each level of a ladder, assert the discrete error, the residual and the continuum order.
   character(*), intent(in) :: label                !< Configuration label.
   integer(I4P), intent(in) :: ratio                !< Tree ratio.
   integer(I4P), intent(in) :: nk                   !< Cells per block along z (1: 2-D, null z).
   integer(I4P), intent(in) :: levels(:)            !< Uniform refinement levels.
   logical,      intent(in) :: periodic             !< Periodic boundaries.
   real(R8P)                :: e_disc(size(levels)) !< Error against the discrete solution, relative.
   real(R8P)                :: e_res(size(levels))  !< Discrete residual, relative.
   real(R8P)                :: e_cont(size(levels)) !< Max error against the continuum solution.
   integer(I4P)             :: n(size(levels))      !< Cells along x.
   integer(I4P)             :: unfilled             !< Phi cells the scatter left.
   real(R8P)                :: order                !< Observed order of the finest pair.
   logical                  :: ok                   !< Ladder pass flag.
   integer(I4P)             :: m                    !< Counter.

   ok = .true.
   do m=1, size(levels)
      call solve_level(ratio=ratio, nk=nk, level=levels(m), periodic=periodic, n=n(m), e_disc=e_disc(m), &
                       e_res=e_res(m), e_cont=e_cont(m), unfilled=unfilled)
      ok = ok .and. e_disc(m) <= TOL_DISCRETE .and. e_res(m) <= TOL_RESIDUAL .and. unfilled == 0_I4P
      if (mpih%myrank == 0) print '(A)', '   '//label//' N = '//trim(str(n(m), .true.))//                  &
         ': |Phi - Phi_h| '//trim(str(e_disc(m), .true.))//', residual '//trim(str(e_res(m), .true.))//      &
         ', |Phi - Phi_exact| '//trim(str(e_cont(m), .true.))//', unfilled '//trim(str(unfilled, .true.))
   enddo
   if (size(levels) > 1) then
      m = size(levels)
      order = log(e_cont(m-1) / e_cont(m)) / log(real(n(m), R8P) / n(m-1))
      ok = ok .and. order >= MIN_ORDER
      if (mpih%myrank == 0) print '(A)', '   '//label//' order against the exact Phi, N = '//                     &
         trim(str(n(m-1), .true.))//' -> '//trim(str(n(m), .true.))//': '//trim(str(order, .true.))//            &
         ' (min '//trim(str(MIN_ORDER, .true.))//')'
   endif
   if (mpih%myrank == 0) then
      if (ok) then
         print '(A)', '   '//label//' PASS'
      else
         print '(A)', '   '//label//' FAIL'
      endif
   endif
   test_passed = test_passed .and. ok
   endsubroutine check_ladder

   subroutine solve_level(ratio, nk, level, periodic, n, e_disc, e_res, e_cont, unfilled)
   !< Build the uniform grid of one level, solve, measure.
   integer(I4P), intent(in)     :: ratio         !< Tree ratio.
   integer(I4P), intent(in)     :: nk            !< Cells per block along z.
   integer(I4P), intent(in)     :: level         !< Uniform refinement level.
   logical,      intent(in)     :: periodic      !< Periodic boundaries.
   integer(I4P), intent(out)    :: n             !< Cells along x.
   real(R8P),    intent(out)    :: e_disc        !< Error against the discrete solution, relative.
   real(R8P),    intent(out)    :: e_res         !< Discrete residual, relative.
   real(R8P),    intent(out)    :: e_cont        !< Max error against the continuum solution.
   integer(I4P), intent(out)    :: unfilled      !< Phi cells the scatter left.
   type(realm_object)           :: realm         !< Realm.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver        !< Solver.
#else
   type(poisson_fft_object)     :: solver        !< Solver.
#endif
   character(:), allocatable    :: filename      !< Input file name.
   real(R8P),    allocatable    :: q(:,:,:,:,:)  !< Field (1, i, j, k, b).
   real(R8P),    allocatable    :: phi(:,:,:,:)  !< Potential (i, j, k, b).
   real(R8P)                    :: xyz(3)        !< Cell centre.
   real(R8P)                    :: h(3)          !< Cell spacing.
   real(R8P)                    :: lap           !< L_h phi.
   real(R8P)                    :: scale(2)      !< Max |Phi_h|, max |4 pi G (rho - rho0)|.
   logical                      :: active(3)     !< Active directions.
   integer(I4P)                 :: ngc, ni, nj   !< Grid dimensions.
   integer(I4P)                 :: b, i, j, k, d !< Counters.
   integer(I4P)                 :: o(3)          !< Unit offset along a direction.
   integer(I4P)                 :: ierr          !< MPI error status.

   filename = 'test_poisson_fft-r'//trim(str(ratio, .true.))//'-nk'//trim(str(nk, .true.))//'-l'// &
              trim(str(level, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, nk=nk, ratio=ratio, level=level)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=4._R8P, nv=1_I4P) ! 64^3 on one rank: 512 blocks
   if (periodic) call realm%adam%grid%set_bc_type(bc_type=[(BC_PERIODIC, d=1, 6)])
   ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   q = 0._R8P
   if (level > 0) call realm%adam%refine_uniform(refinement_levels=level, q=q, do_mpi_redistribute=.true., &
                                                 do_blocks_reorder=.false.)
   active = .not.realm%adam%grid%null_xyz
   do b=1, realm%adam%field%blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(1,i,j,k,b) = rho([realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), &
                                   realm%adam%field%z_cell(k,b)], active)
            enddo
         enddo
      enddo
   enddo
   allocate(phi(1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   phi = SENTINEL
   call solver%initialize(grid=realm%adam%grid, field=realm%adam%field)
   call solver%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi)
   n = solver%n(1)
   h = solver%dxyz

   e_disc = 0._R8P ; e_res = 0._R8P ; e_cont = 0._R8P ; scale = 0._R8P ; unfilled = 0_I4P
   do b=1, realm%adam%field%blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               if (phi(i,j,k,b) == SENTINEL) then
                  unfilled = unfilled + 1_I4P
                  cycle
               endif
               xyz = [realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), realm%adam%field%z_cell(k,b)]
               e_disc   = max(e_disc, abs(phi(i,j,k,b) - phi_exact(xyz, active, h)))
               scale(1) = max(scale(1), abs(phi_exact(xyz, active, h)))
               if (i < 1 .or. i > ni .or. j < 1 .or. j > nj .or. k < 1 .or. k > nk) cycle
               e_cont = max(e_cont, abs(phi(i,j,k,b) - phi_exact(xyz, active)))
               lap = 0._R8P
               do d=1, 3
                  if (.not.active(d)) cycle
                  o = 0 ; o(d) = 1
                  lap = lap + (phi(i+o(1),j+o(2),k+o(3),b) - 2._R8P * phi(i,j,k,b) + &
                               phi(i-o(1),j-o(2),k-o(3),b)) / h(d)**2
               enddo
               e_res    = max(e_res, abs(lap - FACTOR * (q(1,i,j,k,b) - RHO0)))
               scale(2) = max(scale(2), abs(FACTOR * (q(1,i,j,k,b) - RHO0)))
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_disc,   1, MPI_REAL8,   MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_res,    1, MPI_REAL8,   MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_cont,   1, MPI_REAL8,   MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, scale,    2, MPI_REAL8,   MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   e_disc = e_disc / scale(1)
   e_res  = e_res / scale(2)
   if (mpih%myrank == 0 .and. len_trim(dump_dir) > 0) call dump(filename=trim(dump_dir)//'/phi-r'//         &
      trim(str(ratio, .true.))//'-nk'//trim(str(nk, .true.))//'-l'//trim(str(level, .true.))//'.bin', &
      n=solver%n, f=solver%f_global)
   call solver%destroy
   endsubroutine solve_level

   subroutine dump(filename, n, f)
   !< Write the global potential of rank 0: the 3 sizes (int32), then the values (float64, x fastest).
   character(*), intent(in) :: filename !< File name.
   integer(I4P), intent(in) :: n(3)     !< Sizes.
   real(R8P),    intent(in) :: f(:,:,:) !< Global potential.
   integer(I4P)             :: u        !< File unit.

   open(newunit=u, file=filename, access='stream', form='unformatted', action='write', status='replace')
   write(u) n
   write(u) f
   close(u)
   endsubroutine dump

   pure function wave(m, active) result(kv)
   !< Wave vector of mode m (1 or 2); no wave along an inactive direction.
   integer(I4P), intent(in) :: m         !< Mode.
   logical,      intent(in) :: active(3) !< Active directions.
   real(R8P)                :: kv(3)     !< Wave vector.

   if (m == 1) then
      kv = 2._R8P * PI * [1._R8P, 2._R8P, 1._R8P] / L
   else
      kv = 2._R8P * PI * [2._R8P, 1._R8P, 1._R8P] / L
   endif
   where (.not.active) kv = 0._R8P
   endfunction wave

   pure function rho(xyz, active) result(r)
   !< The source density at a point.
   real(R8P), intent(in) :: xyz(3)       !< Point.
   logical,   intent(in) :: active(3)    !< Active directions.
   real(R8P)             :: r            !< Density.
   real(R8P)             :: k1(3), k2(3) !< Wave vectors.

   k1 = wave(1_I4P, active) ; k2 = wave(2_I4P, active)
   r = RHO0 + A1 * sin(k1(1) * xyz(1)) * cos(k1(2) * xyz(2)) * cos(k1(3) * xyz(3)) + A2 * cos(dot_product(k2, xyz))
   endfunction rho

   pure function phi_exact(xyz, active, h) result(p)
   !< The exact potential: continuum (eigenvalue -|k|^2), or discrete with the spacing h (eigenvalue of L_h).
   real(R8P), intent(in)           :: xyz(3)       !< Point.
   logical,   intent(in)           :: active(3)    !< Active directions.
   real(R8P), intent(in), optional :: h(3)         !< Cell spacing: the discrete solution.
   real(R8P)                       :: p            !< Potential.
   real(R8P)                       :: k1(3), k2(3) !< Wave vectors.
   real(R8P)                       :: lam1, lam2   !< Eigenvalues.

   k1 = wave(1_I4P, active) ; k2 = wave(2_I4P, active)
   if (present(h)) then
      lam1 = eigenvalue(k1, active, h) ; lam2 = eigenvalue(k2, active, h)
   else
      lam1 = -dot_product(k1, k1) ; lam2 = -dot_product(k2, k2)
   endif
   p = FACTOR * (A1 * sin(k1(1) * xyz(1)) * cos(k1(2) * xyz(2)) * cos(k1(3) * xyz(3)) / lam1 + &
                 A2 * cos(dot_product(k2, xyz)) / lam2)
   endfunction phi_exact

   pure function eigenvalue(kv, active, h) result(lam)
   !< Eigenvalue of L_h for a mode of wave vector kv.
   real(R8P), intent(in) :: kv(3)     !< Wave vector.
   logical,   intent(in) :: active(3) !< Active directions.
   real(R8P), intent(in) :: h(3)      !< Cell spacing.
   real(R8P)             :: lam       !< Eigenvalue.
   integer(I4P)          :: d         !< Counter.

   lam = 0._R8P
   do d=1, 3
      if (active(d)) lam = lam + (2._R8P * cos(kv(d) * h(d)) - 2._R8P) / h(d)**2
   enddo
   endfunction eigenvalue

   subroutine write_input(filename, nk, ratio, level)
   !< Write the input of one configuration: the sections realm_object%initialize reads, domain [0,1] x [0,2] x [0,3].
   character(*), intent(in) :: filename !< Input file name.
   integer(I4P), intent(in) :: nk       !< Cells per block along z.
   integer(I4P), intent(in) :: ratio    !< Tree ratio.
   integer(I4P), intent(in) :: level    !< Uniform refinement level.
   integer(I4P)             :: u        !< File unit.

   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = test_poisson_fft'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = test_poisson_fft-restart'
   write(u, '(A)') 'restart_save           = 999999'
   write(u, '(A)') 'residuals_save         = 999999'
   write(u, '(A)') 'save_memory_status     = .false.'
   write(u, '(A)') 'save_residual_fields   = .false.'
   write(u, '(A)') 'save_curl_fields       = .false.'
   write(u, '(A)') 'save_divergence_fields = .false.'
   write(u, '(A)') 'save_gradient_fields   = .false.'
   write(u, '(A)') 'save_laplacian_fields  = .false.'
   write(u, '(A)') '[grid]'
   write(u, '(A)') 'ni     = 8'
   write(u, '(A)') 'nj     = 8'
   write(u, '(A)') 'nk     = '//trim(str(nk, .true.))
   write(u, '(A)') 'ngc    = 3'
   write(u, '(A)') 'emin_x = 0.0'
   write(u, '(A)') 'emin_y = 0.0'
   write(u, '(A)') 'emin_z = 0.0'
   write(u, '(A)') 'emax_x = 1.0'
   write(u, '(A)') 'emax_y = 2.0'
   write(u, '(A)') 'emax_z = 3.0'
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = .false.'
   if (nk == 1_I4P) then
      write(u, '(A)') 'null_z = .true.'
   else
      write(u, '(A)') 'null_z = .false.'
   endif
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level      = '//trim(str(max(level, 1_I4P), .true.))
   write(u, '(A)') 'ratio          = '//trim(str(ratio, .true.))
   write(u, '(A)') 'iu_ref_levels  = '//trim(str(level, .true.))
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
   write(u, '(A)') '[linear-algebra]'
   write(u, '(A)') 'smoothing         = gauss-seidel'
   write(u, '(A)') 'iterations_init   = 3'
   write(u, '(A)') 'iterations_coarse = 10'
   write(u, '(A)') 'iterations_fine   = 3'
   write(u, '(A)') 'iterations        = 10'
   write(u, '(A)') 'tolerance         = 1.e-20'
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
endprogram test_poisson_fft
