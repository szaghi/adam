!< Unit test of the FFT solver's boundary conditions per face (issue #84, E1a).
program test_flail_bc
!< Unit test of the FFT solver's boundary conditions per face (issue #84, E1a).
!<
!< **Why this test exists.** FLAIL's clients need more than periodic and isolated boundaries: PRISM's elliptic solves
!< take Dirichlet (PEC) and Neumann faces, often different at the two ends of one axis. The FFT solver handles them by
!< extending the source per axis (odd about a Dirichlet face, even about a Neumann face, over 2n cells, or 4n when the
!< two faces differ; odd about the ghost centres, 2(n+1) cells, for the #87 ghost-Dirichlet kind) and transforming
!< the extension periodically. A wrong extension, length or eigenvalue, a wrong ghost read in the scatter, or images
!< leaking into the energy each break the assertions below.
!<
!< **What it pins.** Per axis, a discrete eigenfunction of the 7-point operator that satisfies the axis's boundary
!< conditions, at the cell centres `x_i = (i - 1/2) dx`, i = 1..n (`k` its mode):
!<   - P/P: `cos(2 pi k x / L)`;  D/D: `sin(pi k x / L)`;  N/N: `cos(pi k x / L)`;
!<   - D/N: `sin(pi (k + 1/2) x / L)`;  N/D: `cos(pi (k + 1/2) x / L)`;
!<   - G/G (#87 ghost-Dirichlet): `sin(pi k i / (n + 1))`.
!< The source is their product (plus a constant `RHO0` when no axis is Dirichlet, which the solve must drop), so the
!< discrete solution is the product times `factor / sum_d mu_d`, `mu_d` the axis's discrete eigenvalue, known in closed
!< form at every cell including the ghost layer (the formulas extend there and carry the boundary condition). Asserted,
!< per configuration, on 1, 2 or 3 ranks:
!<   - `max |phi - phi_h| / max |phi_h|` over the interior and the first ghost layer at round-off (`TOL_DISCRETE`);
!<   - the energy `sum_interior f phi` against the exact sum, relative (`TOL_DISCRETE`): no images in it;
!<   - for D/D, the order against the continuum solution (the same sine, divided by `-|k|^2`) on a ladder, at least
!<     `MIN_ORDER`: the face placement of the wall is second order.
!< `--refuse-mixed-periodic` (a periodic face against a Dirichlet one) and `--refuse-mixed-ghost` (a ghost-Dirichlet
!< face against a face kind) must stop (exit status 1).
!<
!< Usage: `mpirun -np N exe/test_flail_bc [--refuse-mixed-periodic | --refuse-mixed-ghost]` (or `exe/test_flail_bc_fnl`).
#ifdef _FNL
#include "fundal.H"
#endif
use :: adam_flail_parameters,       only : FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, &
                                           FLAIL_BC_DIRICHLET_GHOST
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
real(R8P), parameter :: FACTOR=4._R8P * PI            !< Source factor.
real(R8P), parameter :: RHO0=1.5_R8P                  !< Constant added without a Dirichlet axis (dropped).
real(R8P), parameter :: L(3)=[1._R8P, 2._R8P, 3._R8P] !< Domain lengths.
real(R8P), parameter :: TOL_DISCRETE=1.e-12_R8P       !< Bound on the discrete error and the energy, relative.
real(R8P), parameter :: MIN_ORDER=1.9_R8P             !< Bound on the D/D order against the continuum.
integer(I4P), parameter :: P=FLAIL_BC_PERIODIC        !< Shorthand.
integer(I4P), parameter :: D=FLAIL_BC_DIRICHLET       !< Shorthand.
integer(I4P), parameter :: N=FLAIL_BC_NEUMANN         !< Shorthand.
integer(I4P), parameter :: G=FLAIL_BC_DIRICHLET_GHOST !< Shorthand.
logical        :: test_passed !< Aggregate pass flag.
character(256) :: arg         !< Command argument.
real(R8P)      :: e_disc, e_cont(3), e_energy !< Errors.
integer(I4P)   :: m           !< Counter.
integer(I4P)   :: cfg_bc(2,3)    !< Current case: kinds per face.
integer(I4P)   :: cfg_kmode(3)   !< Current case: mode per axis.
integer(I4P)   :: cfg_ncell(3)   !< Current case: global cells per axis.
logical        :: cfg_null(3)    !< Current case: null directions.
logical        :: cfg_dirichlet  !< Current case: some axis has a Dirichlet face.
real(R8P)      :: order       !< Observed order.

call get_command_argument(1, arg)
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
#ifdef _FNL
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
#endif
test_passed = .true.
select case(trim(arg))
case('--refuse-mixed-periodic')
   call solve_case(bc=reshape([P,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[8,8,8], null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(1), e_energy=e_energy)
   call refusal_failed
case('--refuse-mixed-ghost')
   call solve_case(bc=reshape([G,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[8,8,8], null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(1), e_energy=e_energy)
   call refusal_failed
endselect
! 3-D, one combination per axis kind
call check('P/P D/D N/N', reshape([P,P, D,D, N,N], [2,3]), null_z=.false.)
call check('D/N N/D P/P', reshape([D,N, N,D, P,P], [2,3]), null_z=.false.)
call check('N/N N/N N/N', reshape([N,N, N,N, N,N], [2,3]), null_z=.false.)
call check('D/D D/D D/D', reshape([D,D, D,D, D,D], [2,3]), null_z=.false.)
call check('N/D D/N N/N', reshape([N,D, D,N, N,N], [2,3]), null_z=.false.)
call check('G/G G/G G/G', reshape([G,G, G,G, G,G], [2,3]), null_z=.false.)
! 2-D (null z)
call check('2-D D/N N/N ', reshape([D,N, N,N, P,P], [2,3]), null_z=.true.)
call check('2-D P/P D/D ', reshape([P,P, D,D, P,P], [2,3]), null_z=.true.)
! the D/D order against the continuum
do m=1, 3
   call solve_case(bc=reshape([D,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[4,4,4]*2**m, null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(m), e_energy=e_energy)
enddo
order = log(e_cont(2) / e_cont(3)) / log(2._R8P)
if (mpih%myrank == 0) print '(A)', '   D/D D/D D/D order against the continuum, 16^3 -> 32^3 -> 64^3 cells per axis: '// &
   trim(str(log(e_cont(1) / e_cont(2)) / log(2._R8P), .true.))//', '//trim(str(order, .true.))//' (min '//            &
   trim(str(MIN_ORDER, .true.))//')'
test_passed = test_passed .and. order >= MIN_ORDER
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'test_flail_bc: PASSED'
   else
      print '(A)', 'test_flail_bc: FAILED'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   subroutine check(label, bc, null_z)
   !< Solve one configuration on 4^3 blocks of 4^3 cells (2-D: 4^2 blocks of 4^2), assert the discrete error and the
   !< energy.
   character(*), intent(in) :: label   !< Configuration label.
   integer(I4P), intent(in) :: bc(2,3) !< Kinds per face [lo/hi, axis].
   logical,      intent(in) :: null_z  !< 2-D.
   integer(I4P)             :: nblk(3), nijk(3) !< Blocks and cells per block per direction.
   logical                  :: ok      !< Pass flag.

   nblk = [4,4,4] ; nijk = [4,4,4]
   if (null_z) then
      nblk(3) = 1 ; nijk(3) = 1
   endif
   call solve_case(bc=bc, nblk=nblk, nijk=nijk, null_z=null_z, e_disc=e_disc, e_cont=e_cont(1), e_energy=e_energy)
   ok = e_disc <= TOL_DISCRETE .and. e_energy <= TOL_DISCRETE
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) print '(A)', '   '//label//': |phi - phi_h| '//trim(str(e_disc, .true.))//', energy '// &
      trim(str(e_energy, .true.))//merge(' ok  ', ' FAIL', ok)
   endsubroutine check

   subroutine solve_case(bc, nblk, nijk, null_z, e_disc, e_cont, e_energy)
   !< Build the source of a configuration on this rank's blocks (round robin over the ranks), solve, and return the
   !< errors: discrete (interior and first ghost layer), continuum (interior, D/D only: the sine over -|k|^2) and
   !< energy.
   integer(I4P), intent(in)  :: bc(2,3)   !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in)  :: nblk(3)   !< Blocks per direction.
   integer(I4P), intent(in)  :: nijk(3)   !< Cells per block per direction.
   logical,      intent(in)  :: null_z    !< 2-D.
   real(R8P),    intent(out) :: e_disc    !< Discrete error, relative.
   real(R8P),    intent(out) :: e_cont    !< Continuum error, max.
   real(R8P),    intent(out) :: e_energy  !< Energy error, relative.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver !< The solver.
#else
   type(poisson_fft_object)     :: solver !< The solver.
#endif
   integer(I4P), parameter   :: GW=1_I4P   !< Ghost layers read back.
   integer(I4P)              :: ncell(3)   !< Global cells per direction.
   real(R8P)                 :: dxyz(3)    !< Spacing.
   real(R8P)                 :: mu(3)      !< Discrete eigenvalue per axis.
   real(R8P)                 :: kap(3)     !< Continuum wavenumber per axis.
   integer(I4P)              :: kmode(3)   !< Mode per axis.
   logical                   :: null_xyz(3) !< Null directions.
   logical                   :: has_dirichlet !< Some axis has a Dirichlet face (no mean to drop).
   integer(I4P), allocatable :: coords(:,:) !< This rank's block coordinates.
   real(R8P),    allocatable :: buf(:), out(:) !< Packed source and potential.
   real(R8P)                 :: s, ph, ph_c, f_all, e, e_exact, emax, pmax, cmax !< Work.
   integer(I4P)              :: nbl, b, bi, bj, bk, i, j, k, c, d_, nb_total, ierr !< Counters.

   null_xyz = [.false., .false., null_z]
   ncell = nblk * nijk
   dxyz = L / ncell
   kmode = [2, 1, 1]
   has_dirichlet = any(bc == D) .or. any(bc == G)
   cfg_bc = bc ; cfg_kmode = kmode ; cfg_ncell = ncell ; cfg_null = null_xyz ; cfg_dirichlet = has_dirichlet
   do d_=1, 3
      call axis_eigen(bc(:,d_), kmode(d_), ncell(d_), dxyz(d_), L(d_), null_xyz(d_), mu(d_), kap(d_))
   enddo
   ! this rank's blocks, round robin
   nb_total = product(nblk)
   nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) == mpih%myrank) nbl = nbl + 1
   enddo
   allocate(coords(3,max(nbl,1)), buf(max(nbl*product(nijk),1)), out(max(nbl*product(nijk+2*GW),1)))
   c = 0 ; nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) /= mpih%myrank) cycle
      nbl = nbl + 1
      bi = mod(b, nblk(1)) ; bj = mod(b / nblk(1), nblk(2)) ; bk = b / (nblk(1) * nblk(2))
      coords(:,nbl) = [bi, bj, bk]
      do k=1, nijk(3)
         do j=1, nijk(2)
            do i=1, nijk(1)
               c = c + 1
               buf(c) = source(bi*nijk(1)+i, bj*nijk(2)+j, bk*nijk(3)+k)
            enddo
         enddo
      enddo
   enddo
   call solver%initialize_blocks(nblk=nblk, nijk=nijk, ngc=GW, dxyz=dxyz, level=0_I4P, null_xyz=null_xyz, &
                                 isolated=.false., bc=bc)
   out = 0._R8P
   call solver%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=GW, out=out, energy=e)
   ! errors over this rank's blocks, then reduced
   emax = 0._R8P ; pmax = 0._R8P ; cmax = 0._R8P ; e_exact = 0._R8P
   c = 0
   do b=1, nbl
      do k=1-GW, nijk(3)+GW ! along a null z the scatter stores the ghost layers too; the mode is 1 there
         do j=1-GW, nijk(2)+GW
            do i=1-GW, nijk(1)+GW
               c = c + 1
               ph = FACTOR * mode(coords(1,b)*nijk(1)+i, coords(2,b)*nijk(2)+j, coords(3,b)*nijk(3)+k) / sum(mu)
               emax = max(emax, abs(out(c) - ph))
               pmax = max(pmax, abs(ph))
               if (i>=1 .and. i<=nijk(1) .and. j>=1 .and. j<=nijk(2) .and. k>=1 .and. k<=nijk(3)) then
                  e_exact = e_exact + source(coords(1,b)*nijk(1)+i, coords(2,b)*nijk(2)+j, coords(3,b)*nijk(3)+k) * ph
                  ph_c = -FACTOR * mode(coords(1,b)*nijk(1)+i, coords(2,b)*nijk(2)+j, coords(3,b)*nijk(3)+k) / &
                         sum(kap**2)
                  cmax = max(cmax, abs(out(c) - ph_c))
               endif
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, emax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, pmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, cmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_exact, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   e_disc = emax / pmax
   e_cont = cmax
   e_energy = abs(e - e_exact) / abs(e_exact)
   call solver%destroy
   endsubroutine solve_case

   function mode(gi, gj, gk) result(v)
   !< The product of the axis eigenfunctions at the global cell (gi, gj, gk), ghosts included (indices 0, n + 1).
   integer(I4P), intent(in) :: gi, gj, gk !< Global cell indices.
   real(R8P)          :: v       !< Value.

   v = axis_fun(cfg_bc(:,1), cfg_kmode(1), cfg_ncell(1), gi, cfg_null(1)) * &
       axis_fun(cfg_bc(:,2), cfg_kmode(2), cfg_ncell(2), gj, cfg_null(2)) * &
       axis_fun(cfg_bc(:,3), cfg_kmode(3), cfg_ncell(3), gk, cfg_null(3))
   endfunction mode

   function source(gi, gj, gk) result(v)
   !< The source at an interior cell: the mode, plus RHO0 without a Dirichlet axis (the solve must drop it).
   integer(I4P), intent(in) :: gi, gj, gk !< Global cell indices.
   real(R8P)          :: v       !< Value.

   v = mode(gi, gj, gk)
   if (.not.cfg_dirichlet) v = v + RHO0
   endfunction source

   pure subroutine axis_eigen(kinds, k, n, dx, len, is_null, mu, kap)
   !< The discrete eigenvalue of the 7-point operator and the continuum wavenumber of an axis's eigenfunction.
   integer(I4P), intent(in)  :: kinds(2) !< Kinds at the two faces.
   integer(I4P), intent(in)  :: k        !< Mode.
   integer(I4P), intent(in)  :: n        !< Cells.
   real(R8P),    intent(in)  :: dx, len  !< Spacing and length.
   logical,      intent(in)  :: is_null  !< Null direction.
   real(R8P),    intent(out) :: mu, kap  !< Discrete eigenvalue, continuum wavenumber.
   real(R8P)                 :: th       !< Phase step per cell.

   mu = 0._R8P ; kap = 0._R8P
   if (is_null) return
   th = phase_step(kinds, k, n)
   mu = (2._R8P * cos(th) - 2._R8P) / dx**2
   kap = th / dx
   endsubroutine axis_eigen

   pure function phase_step(kinds, k, n) result(th)
   !< The eigenfunction's phase step per cell of an axis.
   integer(I4P), intent(in) :: kinds(2) !< Kinds at the two faces.
   integer(I4P), intent(in) :: k, n     !< Mode, cells.
   real(R8P)                :: th       !< Phase step.

   if (all(kinds == P)) then
      th = 2._R8P * PI * k / n
   elseif (all(kinds == G)) then
      th = PI * k / (n + 1)
   elseif (kinds(1) == kinds(2)) then
      th = PI * k / n
   else
      th = PI * (k + 0.5_R8P) / n
   endif
   endfunction phase_step

   pure function axis_fun(kinds, k, n, gi, is_null) result(v)
   !< The axis eigenfunction at the global cell index gi (1..n interior, 0 and n + 1 the ghosts).
   integer(I4P), intent(in) :: kinds(2) !< Kinds at the two faces.
   integer(I4P), intent(in) :: k, n, gi !< Mode, cells, index.
   logical,      intent(in) :: is_null  !< Null direction.
   real(R8P)                :: v        !< Value.
   real(R8P)                :: th, x    !< Phase step, cell phase.

   v = 1._R8P
   if (is_null) return
   th = phase_step(kinds, k, n)
   if (all(kinds == G)) then
      v = sin(th * gi)
      return
   endif
   x = th * (gi - 0.5_R8P)
   if (all(kinds == P)) then
      v = cos(x)
   elseif (kinds(1) == D) then
      v = sin(x)
   else
      v = cos(x)
   endif
   endfunction axis_fun

   subroutine refusal_failed
   !< A refusal configuration ran: the test fails.
   if (mpih%myrank == 0) print '(A)', 'test_flail_bc: FAILED (the solver accepted '//trim(arg)//')'
   call mpih%finalize
   error stop 1
   endsubroutine refusal_failed
endprogram test_flail_bc
