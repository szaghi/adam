!< Unit test of the FFT solver's boundary conditions per face (issue #84, E1a) and its wide operator D(G) (E2a).
program test_flail_bc
!< Unit test of the FFT solver's boundary conditions per face (issue #84, E1a) and its wide operator D(G) (E2a).
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
!< **The wide operator (E2a).** The same configurations with `dg_order = p` (2, 4, 6, 8): PRISM's D(G), the order-p
!< centred first derivative composed with itself. The axis eigenfunctions above are eigenfunctions of it too (they are
!< reflections about the faces, which the wide stencil commutes with); their eigenvalue `mu_d` is computed here in
!< real space, independently of the solver's symbol: `D1` applied twice to the eigenfunction, extended over `p` ghost
!< cells by its own formula, divided by it. Asserted as above, the ghost layers read back to depth `p` (the depth
!< PRISM's D and G read), and the D/D order of p = 4 against the continuum (`MIN_ORDER_WIDE`).
!<
!< **Zero-ghost faces (E2b).** PRISM's Dirichlet face writes zero in every ghost layer; with the wide operator it is no
!< reflection, and the solver runs preconditioned CG. No closed form: the solution is manufactured, a smooth field `u`
!< whose source `D(G) u / factor` the test computes itself (the ghosts by the faces' rules, `D1` applied twice). Per
!< configuration (a zero-ghost axis next to periodic, reflected and zero-ghost axes) and order: `u` recovered to
!< `TOL_PCG`, the ghost layers returned to depth `p` equal to the rules applied to the interior returned, the CG
!< iterations at most `MAX_PCG_IT`; at p = 4 the iterations on 32 cells per axis within one of those on 16 (bounded
!< in the grid: the preconditioner is spectrally equivalent); the compact operator with zero ghosts equal to it with
!< ghost-Dirichlet faces, bitwise.
!<
!< **PRISM's operators (E3a).** The identity PRISM's elliptic solves rest on: PRISM computes the correction as the
!< divergence of the gradient with its own library operators (`adam_fdv_operators_library`, centred FD or FV of order
!< 2s), the gradient read from the solution's ghost layers. With a manufactured source `f = div(grad u)` (the test's
!< ghost rules for u), FLAIL's solution phi, its returned ghosts to depth 2s, gives `div(grad phi) = f` to `TOL_PRISM`
!< (zero-ghost faces, solved by CG: `TOL_PCG`) for FD and FV, s = 1..4, periodic, PEC on a scalar (odd), PEC on the
!< normal component (even), mixed, zero-ghost and zero-ghost next to PEC faces: the solver inverts exactly the operator
!< PRISM applies.
!<
!< **Free-space ghosts (E3c).** PRISM's exact/open face: every ghost layer holds the continuum free-space potential of
!< the source, `-sum_y factor f(y) dV / (4 pi |x - y|)`, computed by the solver on rank 0 from the gathered source. A
!< smooth source on 2^3 blocks of 4^3 cells (round robin over the ranks), per configuration (free faces everywhere,
!< next to Dirichlet and periodic axes, mixed with a zero-ghost face on one axis) and order p = 2, 4, 6: the ghost slabs
!< returned equal to the direct sum computed here (`TOL_DISCRETE`, relative), D(G) of the solution with the ghosts
!< returned equal to `factor f` (`TOL_FREE`), the CG iterations at most `MAX_PCG_IT`; the same on 1, 2 and 3 ranks, so
!< the values do not depend on the decomposition (the legacy sum ran over the rank's own blocks).
!<
!< `--refuse-mixed-periodic` (a periodic face against a Dirichlet one), `--refuse-mixed-ghost` (a ghost-Dirichlet
!< face against a face kind), `--refuse-wide-isolated` (the wide operator with isolated boundaries) and
!< `--refuse-zero-mixed` (a zero-ghost face against a Dirichlet one) and `--refuse-free-compact` (free-space ghosts with
!< the compact operator) must stop (exit status 1).
!<
!< Usage: `mpirun -np N exe/test_flail_bc [--refuse-mixed-periodic | --refuse-mixed-ghost | --refuse-wide-isolated |
!< --refuse-zero-mixed | --refuse-free-compact]` (or `exe/test_flail_bc_fnl`).
#ifdef _FNL
#include "fundal.H"
#endif
use :: adam_fdv_operators_library,  only : compute_divergence_fd_centered, compute_divergence_fv_centered, &
                                           compute_gradient_fd_centered, compute_gradient_fv_centered
use :: adam_flail_parameters,       only : FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, &
                                           FLAIL_BC_DIRICHLET_GHOST, FLAIL_BC_ZERO_GHOSTS, FLAIL_BC_FREE_GHOSTS
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
real(R8P), parameter :: MIN_ORDER_WIDE=3.8_R8P        !< Bound on the D/D order of the wide p = 4 operator.
! centred first-derivative coefficients a_j of order p: D1 q_i = sum_j a_j (q_{i+j} - q_{i-j}) / h (column p/2)
real(R8P), parameter :: A_D1(4,4)=reshape([1._R8P/2._R8P, 0._R8P, 0._R8P, 0._R8P,                         &
                                           2._R8P/3._R8P, -1._R8P/12._R8P, 0._R8P, 0._R8P,                 &
                                           3._R8P/4._R8P, -3._R8P/20._R8P, 1._R8P/60._R8P, 0._R8P,         &
                                           4._R8P/5._R8P, -1._R8P/5._R8P, 4._R8P/105._R8P, -1._R8P/280._R8P], &
                                          [4, 4]) !< a_j [j, p/2].
integer(I4P), parameter :: P=FLAIL_BC_PERIODIC        !< Shorthand.
integer(I4P), parameter :: D=FLAIL_BC_DIRICHLET       !< Shorthand.
integer(I4P), parameter :: N=FLAIL_BC_NEUMANN         !< Shorthand.
integer(I4P), parameter :: G=FLAIL_BC_DIRICHLET_GHOST !< Shorthand.
integer(I4P), parameter :: Z=FLAIL_BC_ZERO_GHOSTS     !< Shorthand.
integer(I4P), parameter :: F=FLAIL_BC_FREE_GHOSTS     !< Shorthand.
real(R8P),    parameter :: TOL_PCG=1.e-9_R8P          !< Bound on the zero-ghost solution's error, relative.
integer(I4P), parameter :: MAX_PCG_IT=40_I4P          !< Bound on the zero-ghost CG iterations.
real(R8P),    parameter :: TOL_PRISM=1.e-11_R8P        !< Bound on |div(grad phi) - f| / max |f|, direct solves.
real(R8P),    parameter :: TOL_FREE=1.e-10_R8P         !< Bound on |D(G) phi - factor f| / max |factor f|, free ghosts.
logical        :: test_passed !< Aggregate pass flag.
character(256) :: arg         !< Command argument.
real(R8P)      :: e_disc, e_cont(3), e_energy !< Errors.
integer(I4P)   :: m, po       !< Counters.
integer(I4P)   :: cfg_bc(2,3)    !< Current case: kinds per face.
integer(I4P)   :: cfg_kmode(3)   !< Current case: mode per axis.
integer(I4P)   :: cfg_ncell(3)   !< Current case: global cells per axis.
integer(I4P)   :: cfg_order      !< Current case: 0 compact, p the wide D(G).
logical        :: cfg_null(3)    !< Current case: null directions.
logical        :: cfg_dirichlet  !< Current case: some axis has a Dirichlet face.
real(R8P)      :: order       !< Observed order.

call get_command_argument(1, arg)
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
#ifdef _FNL
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
#endif
test_passed = .true.
cfg_order = 0 ! the refusals below run before the order loop sets it
select case(trim(arg))
case('--refuse-mixed-periodic')
   call solve_case(bc=reshape([P,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[8,8,8], null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(1), e_energy=e_energy)
   call refusal_failed
case('--refuse-mixed-ghost')
   call solve_case(bc=reshape([G,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[8,8,8], null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(1), e_energy=e_energy)
   call refusal_failed
case('--refuse-wide-isolated')
   call refuse_wide_isolated
   call refusal_failed
case('--refuse-zero-mixed')
   call check_zero('refusal    ', reshape([Z,D, P,P, P,P], [2,3]), 4_I4P, [2,2,2], [8,8,8], .false., m)
   call refusal_failed
case('--refuse-free-compact')
   call check_free('refusal    ', reshape([F,F, F,F, F,F], [2,3]), 0_I4P)
   call refusal_failed
endselect
! the compact operator (E1a) and the wide D(G) of order 2, 4, 6, 8 (E2a)
do po=0, 8, 2
   if (po == 2 .and. mpih%myrank == 0) print '(A)', '   the wide operator D(G) (E2a):'
   cfg_order = po
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
enddo
! zero-ghost faces (E2b): preconditioned CG, a manufactured solution
if (mpih%myrank == 0) print '(A)', '   zero-ghost faces (E2b), preconditioned CG:'
do po=2, 8, 2
   call check_zero('Z/Z P/P P/P', reshape([Z,Z, P,P, P,P], [2,3]), po, [4,4,4], [4,4,4], .false., m)
   call check_zero('Z/Z D/D N/N', reshape([Z,Z, D,D, N,N], [2,3]), po, [4,4,4], [4,4,4], .false., m)
   call check_zero('Z/Z Z/Z Z/Z', reshape([Z,Z, Z,Z, Z,Z], [2,3]), po, [4,4,4], [4,4,4], .false., m)
   call check_zero('2-D Z/Z N/N', reshape([Z,Z, N,N, P,P], [2,3]), po, [4,4,1], [4,4,1], .true.,  m)
enddo
block
integer(I4P) :: it16, it32 !< Iterations on 16 and 32 cells per axis.
call check_zero('Z/Z D/D N/N', reshape([Z,Z, D,D, N,N], [2,3]), 4_I4P, [4,4,4], [4,4,4], .false., it16)
call check_zero('Z/Z D/D N/N', reshape([Z,Z, D,D, N,N], [2,3]), 4_I4P, [4,4,4], [8,8,8], .false., it32)
if (mpih%myrank == 0) print '(A)', '   p = 4, Z/Z D/D N/N, CG iterations on 16 and 32 cells per axis: '// &
   trim(str(it16, .true.))//', '//trim(str(it32, .true.))//merge(' ok  ', ' FAIL', abs(it32 - it16) <= 1)
test_passed = test_passed .and. abs(it32 - it16) <= 1
endblock
call check_zero_compact
! PRISM's own operators (E3a): div(grad phi) = f with the library's FD and FV centred derivatives
if (mpih%myrank == 0) print '(A)', '   PRISM''s operators (E3a), div(grad phi) - f with adam_fdv_operators_library:'
do po=1, 4
   call check_prism_ops('P/P P/P P/P', reshape([P,P, P,P, P,P], [2,3]), po)
   call check_prism_ops('D/D D/D D/D', reshape([D,D, D,D, D,D], [2,3]), po)
   call check_prism_ops('N/N D/D D/D', reshape([N,N, D,D, D,D], [2,3]), po)
   call check_prism_ops('D/N P/P N/D', reshape([D,N, P,P, N,D], [2,3]), po)
   call check_prism_ops('Z/Z Z/Z Z/Z', reshape([Z,Z, Z,Z, Z,Z], [2,3]), po)
   call check_prism_ops('Z/Z D/D N/N', reshape([Z,Z, D,D, N,N], [2,3]), po)
enddo
! free-space ghosts (E3c): PRISM's exact/open face, the ghosts summed by the solver over the gathered source
if (mpih%myrank == 0) print '(A)', '   free-space ghosts (E3c), preconditioned CG:'
do po=2, 6, 2
   call check_free('F/F F/F F/F', reshape([F,F, F,F, F,F], [2,3]), po)
   call check_free('F/F D/D P/P', reshape([F,F, D,D, P,P], [2,3]), po)
   call check_free('F/Z N/N F/F', reshape([F,Z, N,N, F,F], [2,3]), po)
enddo
! the wide p = 4 D/D order against the continuum
cfg_order = 4
do m=1, 3
   call solve_case(bc=reshape([D,D, D,D, D,D], [2,3]), nblk=[2,2,2], nijk=[4,4,4]*2**m, null_z=.false., &
                   e_disc=e_disc, e_cont=e_cont(m), e_energy=e_energy)
enddo
order = log(e_cont(2) / e_cont(3)) / log(2._R8P)
if (mpih%myrank == 0) print '(A)', '   wide p = 4, D/D D/D D/D order against the continuum, 16^3 -> 32^3 -> 64^3: '// &
   trim(str(log(e_cont(1) / e_cont(2)) / log(2._R8P), .true.))//', '//trim(str(order, .true.))//' (min '//       &
   trim(str(MIN_ORDER_WIDE, .true.))//')'
test_passed = test_passed .and. order >= MIN_ORDER_WIDE
! the compact D/D order against the continuum
cfg_order = 0
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
   if (mpih%myrank == 0) print '(A)', '   '//op_label()//' '//label//': |phi - phi_h| '//trim(str(e_disc, .true.))// &
      ', energy '//trim(str(e_energy, .true.))//merge(' ok  ', ' FAIL', ok)
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
   integer(I4P)              :: gw         !< Ghost layers read back: 1, or p for the wide operator.
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

   gw = max(1_I4P, cfg_order)
   null_xyz = [.false., .false., null_z]
   ncell = nblk * nijk
   dxyz = L / ncell
   kmode = [2, 1, 1]
   has_dirichlet = any(bc == D) .or. any(bc == G)
   cfg_bc = bc ; cfg_kmode = kmode ; cfg_ncell = ncell ; cfg_null = null_xyz ; cfg_dirichlet = has_dirichlet
   do d_=1, 3
      call axis_eigen(bc(:,d_), kmode(d_), ncell(d_), dxyz(d_), null_xyz(d_), cfg_order, mu(d_), kap(d_))
   enddo
   ! this rank's blocks, round robin
   nb_total = product(nblk)
   nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) == mpih%myrank) nbl = nbl + 1
   enddo
   allocate(coords(3,max(nbl,1)), buf(max(nbl*product(nijk),1)), out(max(nbl*product(nijk+2*gw),1)))
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
   call solver%initialize_blocks(nblk=nblk, nijk=nijk, ngc=gw, dxyz=dxyz, level=0_I4P, null_xyz=null_xyz, &
                                 isolated=.false., bc=bc, dg_order=cfg_order)
   out = 0._R8P
   call solver%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=gw, out=out, energy=e)
   ! errors over this rank's blocks, then reduced
   emax = 0._R8P ; pmax = 0._R8P ; cmax = 0._R8P ; e_exact = 0._R8P
   c = 0
   do b=1, nbl
      do k=1-gw, nijk(3)+gw ! along a null z the scatter stores the ghost layers too; the mode is 1 there
         do j=1-gw, nijk(2)+gw
            do i=1-gw, nijk(1)+gw
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

   pure subroutine axis_eigen(kinds, k, n, dx, is_null, dg_order, mu, kap)
   !< The discrete eigenvalue and the continuum wavenumber of an axis's eigenfunction: the 7-point operator's in
   !< closed form, or the wide D(G)'s in real space (D1 applied twice to the eigenfunction, extended by its own formula
   !< over the ghost cells, divided by it at the cell where it is largest), independently of the solver's symbol.
   integer(I4P), intent(in)  :: kinds(2) !< Kinds at the two faces.
   integer(I4P), intent(in)  :: k        !< Mode.
   integer(I4P), intent(in)  :: n        !< Cells.
   real(R8P),    intent(in)  :: dx       !< Spacing.
   logical,      intent(in)  :: is_null  !< Null direction.
   integer(I4P), intent(in)  :: dg_order !< 0 compact, p the wide D(G).
   real(R8P),    intent(out) :: mu, kap  !< Discrete eigenvalue, continuum wavenumber.
   real(R8P)                 :: th       !< Phase step per cell.
   real(R8P)                 :: g(-4:4)  !< D1 of the eigenfunction around the cell.
   real(R8P)                 :: dg, fmax !< D1 D1 at the cell, the eigenfunction there.
   integer(I4P)              :: i0, i, j, o !< Cell, counters.

   mu = 0._R8P ; kap = 0._R8P
   if (is_null) return
   th = phase_step(kinds, k, n)
   kap = th / dx
   if (dg_order == 0_I4P) then
      mu = (2._R8P * cos(th) - 2._R8P) / dx**2
      return
   endif
   i0 = 1 ; fmax = 0._R8P
   do i=1, n
      if (abs(axis_fun(kinds, k, n, i, .false.)) > fmax) then
         fmax = abs(axis_fun(kinds, k, n, i, .false.)) ; i0 = i
      endif
   enddo
   g = 0._R8P
   do o=-dg_order/2, dg_order/2
      do j=1, dg_order/2
         g(o) = g(o) + A_D1(j,dg_order/2) * (axis_fun(kinds, k, n, i0+o+j, .false.) - &
                                             axis_fun(kinds, k, n, i0+o-j, .false.))
      enddo
   enddo
   dg = 0._R8P
   do j=1, dg_order/2
      dg = dg + A_D1(j,dg_order/2) * (g(j) - g(-j))
   enddo
   mu = dg / axis_fun(kinds, k, n, i0, .false.) / dx**2
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

   subroutine check_zero(label, bc, p_, nblk, nijk, null_z, iters)
   !< Zero-ghost faces with the wide D(G) of order p_ (E2b): the manufactured `u`, its source computed here, solved; the
   !< error, the ghosts returned against the rules, the CG iterations asserted.
   character(*), intent(in)  :: label   !< Configuration label.
   integer(I4P), intent(in)  :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in)  :: p_      !< Order.
   integer(I4P), intent(in)  :: nblk(3), nijk(3) !< Blocks and cells per block per direction.
   logical,      intent(in)  :: null_z  !< 2-D.
   integer(I4P), intent(out) :: iters   !< CG iterations.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver !< The solver.
#else
   type(poisson_fft_object)     :: solver !< The solver.
#endif
   integer(I4P)              :: ncell(3), nbl, nb_total, b, bi, bj, bk, i, j, k, c, ierr !< Sizes, counters.
   real(R8P)                 :: dxyz(3)  !< Spacing.
   logical                   :: null_xyz(3) !< Null directions.
   real(R8P),    allocatable :: u(:,:,:), f(:,:,:), uo(:,:,:) !< Exact solution, source, solution returned (padded).
   integer(I4P), allocatable :: coords(:,:) !< This rank's block coordinates.
   real(R8P),    allocatable :: buf(:), out(:) !< Packed source and potential.
   real(R8P)                 :: emax, umax, gmax !< Errors.
   logical                   :: ok  !< Pass flag.

   null_xyz = [.false., .false., null_z]
   ncell = nblk * nijk
   dxyz = L / ncell
   allocate(u(1-p_:ncell(1)+p_,1-p_:ncell(2)+p_,1-p_:ncell(3)+p_), f(ncell(1),ncell(2),ncell(3)))
   allocate(uo, mold=u)
   u = 0._R8P
   do k=1, ncell(3)
      do j=1, ncell(2)
         do i=1, ncell(1)
            u(i,j,k) = sin(0.7_R8P * i + 0.3_R8P) * cos(0.5_R8P * j - 0.2_R8P) * cos(0.4_R8P * k + 0.1_R8P) + &
                       0.25_R8P * cos(1.1_R8P * i - 0.9_R8P * j + 0.3_R8P * k)
            if (null_z) u(i,j,k) = sin(0.7_R8P * i + 0.3_R8P) * cos(0.5_R8P * j - 0.2_R8P) + &
                                   0.25_R8P * cos(1.1_R8P * i - 0.9_R8P * j)
         enddo
      enddo
   enddo
   call rule_ghosts(u, ncell, bc, null_xyz, p_)
   call d1d1(u, ncell, dxyz, null_xyz, p_, f)
   f = f / FACTOR
   ! this rank's blocks, round robin
   nb_total = product(nblk)
   nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) == mpih%myrank) nbl = nbl + 1
   enddo
   allocate(coords(3,max(nbl,1)), buf(max(nbl*product(nijk),1)), out(max(nbl*product(nijk+2*p_),1)))
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
               buf(c) = f(bi*nijk(1)+i, bj*nijk(2)+j, bk*nijk(3)+k)
            enddo
         enddo
      enddo
   enddo
   call solver%initialize_blocks(nblk=nblk, nijk=nijk, ngc=p_, dxyz=dxyz, level=0_I4P, null_xyz=null_xyz, &
                                 isolated=.false., bc=bc, dg_order=p_)
   out = 0._R8P
   call solver%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=p_, out=out)
   iters = solver%iterations
   call MPI_BCAST(iters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   ! the interior returned, then the rules applied to it, against the ghosts returned
   uo = 0._R8P
   c = 0
   do b=1, nbl
      do k=1-p_, nijk(3)+p_
         do j=1-p_, nijk(2)+p_
            do i=1-p_, nijk(1)+p_
               c = c + 1
               if (i>=1 .and. i<=nijk(1) .and. j>=1 .and. j<=nijk(2) .and. k>=1 .and. k<=nijk(3)) &
                  uo(coords(1,b)*nijk(1)+i, coords(2,b)*nijk(2)+j, coords(3,b)*nijk(3)+k) = out(c)
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, uo, size(uo), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr) ! each cell from one rank
   emax = maxval(abs(uo(1:ncell(1),1:ncell(2),1:ncell(3)) - u(1:ncell(1),1:ncell(2),1:ncell(3))))
   umax = maxval(abs(u(1:ncell(1),1:ncell(2),1:ncell(3))))
   call rule_ghosts(uo, ncell, bc, null_xyz, p_)
   gmax = 0._R8P
   c = 0
   do b=1, nbl
      do k=1-p_, nijk(3)+p_
         do j=1-p_, nijk(2)+p_
            do i=1-p_, nijk(1)+p_
               c = c + 1
               if (null_z .and. (k < 1 .or. k > nijk(3))) cycle ! a null z: the ghost layers copy the cell
               gmax = max(gmax, abs(out(c) - uo(coords(1,b)*nijk(1)+i, coords(2,b)*nijk(2)+j, coords(3,b)*nijk(3)+k)))
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, gmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ok = emax / umax <= TOL_PCG .and. gmax <= TOL_DISCRETE * umax .and. iters <= MAX_PCG_IT
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) print '(A)', '   D(G) p = '//trim(str(p_, .true.))//' '//label//' ('// &
      trim(str(ncell(1), .true.))//' cells per axis): |u - u_exact| '//trim(str(emax / umax, .true.))// &
      ', ghosts against the rules '//trim(str(gmax, .true.))//', CG iterations '//trim(str(iters, .true.))// &
      merge(' ok  ', ' FAIL', ok)
   call solver%destroy
   endsubroutine check_zero

   subroutine check_free(label, bc, p_)
   !< Free-space ghosts with the wide D(G) of order p_ (E3c): a smooth source on 2^3 blocks of 4^3 cells, solved; the
   !< ghost slabs of the free faces against the direct sum `-FACTOR sum_y f(y) dV / (4 pi |x - y|)` computed here, and
   !< D(G) of the solution (its ghosts: the slabs returned on the free faces, the test's rules on the others) against
   !< `FACTOR f`. p_ = 0 is the refusal (the compact operator).
   character(*), intent(in)  :: label   !< Configuration label.
   integer(I4P), intent(in)  :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in)  :: p_      !< Order.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver !< The solver.
#else
   type(poisson_fft_object)     :: solver !< The solver.
#endif
   integer(I4P), parameter   :: NBLK(3)=[2,2,2], NIJK(3)=[4,4,4] !< Blocks and cells per block per direction.
   integer(I4P)              :: ncell(3), nbl, nb_total, gw, b, bi, bj, bk, i, j, k, ii, jj, kk, c, ax, side, ierr
   integer(I4P)              :: lo(3), hi(3), x(3), iters !< A slab, a global cell, CG iterations.
   integer(I4P)              :: bc_rule(2,3) !< The kinds with the free faces zeroed (the test's rules).
   real(R8P)                 :: dxyz(3), s, gref, gmax, gscale, res !< Spacing, sums, errors.
   real(R8P),    allocatable :: src(:,:,:), uo(:,:,:), y(:,:,:) !< Source, solution returned (padded), D(G) of it.
   integer(I4P), allocatable :: coords(:,:) !< This rank's block coordinates.
   real(R8P),    allocatable :: buf(:), out(:) !< Packed source and potential.
   logical                   :: ok  !< Pass flag.

   gw = max(p_, 1_I4P)
   ncell = NBLK * NIJK
   dxyz = L / ncell
   allocate(src(ncell(1),ncell(2),ncell(3)), uo(1-gw:ncell(1)+gw,1-gw:ncell(2)+gw,1-gw:ncell(3)+gw), y(ncell(1),ncell(2),ncell(3)))
   do k=1, ncell(3)
      do j=1, ncell(2)
         do i=1, ncell(1)
            src(i,j,k) = sin(0.7_R8P * i + 0.3_R8P) * cos(0.5_R8P * j - 0.2_R8P) * cos(0.4_R8P * k + 0.1_R8P) + 0.3_R8P
         enddo
      enddo
   enddo
   nb_total = product(NBLK)
   nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) == mpih%myrank) nbl = nbl + 1
   enddo
   allocate(coords(3,max(nbl,1)), buf(max(nbl*product(NIJK),1)), out(max(nbl*product(NIJK+2*gw),1)))
   c = 0 ; nbl = 0
   do b=0, nb_total-1
      if (mod(b, mpih%procs_number) /= mpih%myrank) cycle
      nbl = nbl + 1
      bi = mod(b, NBLK(1)) ; bj = mod(b / NBLK(1), NBLK(2)) ; bk = b / (NBLK(1) * NBLK(2))
      coords(:,nbl) = [bi, bj, bk]
      do k=1, NIJK(3)
         do j=1, NIJK(2)
            do i=1, NIJK(1)
               c = c + 1
               buf(c) = src(bi*NIJK(1)+i, bj*NIJK(2)+j, bk*NIJK(3)+k)
            enddo
         enddo
      enddo
   enddo
   call solver%initialize_blocks(nblk=NBLK, nijk=NIJK, ngc=gw, dxyz=dxyz, level=0_I4P, null_xyz=[.false.,.false.,.false.], &
                                 isolated=.false., bc=bc, dg_order=p_)
   out = 0._R8P
   call solver%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=gw, out=out)
   iters = solver%iterations
   call MPI_BCAST(iters, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   ! the interior and the face ghost slabs returned (each cell from one block), on every rank
   uo = 0._R8P
   c = 0
   do b=1, nbl
      do k=1-gw, NIJK(3)+gw
         do j=1-gw, NIJK(2)+gw
            do i=1-gw, NIJK(1)+gw
               c = c + 1
               x = coords(:,b) * NIJK + [i, j, k]
               ! a block's own cells, and the face-slab cells outside the domain along one axis and inside this block
               ! along the other two, from the boundary block of that axis only (the ghost layers may be deeper than a
               ! block): edges, corners and inter-block ghosts are taken from no block
               if (count(x < 1 .or. x > ncell) > 1) cycle
               if (count([i, j, k] < 1 .or. [i, j, k] > NIJK) /= count(x < 1 .or. x > ncell)) cycle
               if (any((x < 1 .and. coords(:,b) /= 0) .or. (x > ncell .and. coords(:,b) /= NBLK - 1))) cycle
               uo(x(1),x(2),x(3)) = out(c)
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, uo, size(uo), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   ! the free slabs against the direct sum
   gmax = 0._R8P ; gscale = 0._R8P
   do ax=1, 3
      do side=1, 2
         if (bc(side,ax) /= F) cycle
         lo = 1 ; hi = ncell
         if (side == 1) then
            lo(ax) = 1 - gw ; hi(ax) = 0
         else
            lo(ax) = ncell(ax) + 1 ; hi(ax) = ncell(ax) + gw
         endif
         do k=lo(3), hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  s = 0._R8P
                  do kk=1, ncell(3)
                     do jj=1, ncell(2)
                        do ii=1, ncell(1)
                           s = s + src(ii,jj,kk) / norm2([real(i-ii, R8P), real(j-jj, R8P), real(k-kk, R8P)] * dxyz)
                        enddo
                     enddo
                  enddo
                  gref = -FACTOR * product(dxyz) * s / (4._R8P * PI)
                  gmax = max(gmax, abs(uo(i,j,k) - gref)) ; gscale = max(gscale, abs(gref))
               enddo
            enddo
         enddo
      enddo
   enddo
   gmax = gmax / gscale
   ! D(G) of the solution: the free slabs as returned, the test's rules on the other faces
   bc_rule = merge(Z, bc, bc == F)
   block
   real(R8P), allocatable :: ua(:,:,:) !< The solution with its ghosts for D(G).
   ua = uo
   call rule_ghosts(ua, ncell, bc_rule, [.false.,.false.,.false.], gw)
   do ax=1, 3
      do side=1, 2
         if (bc(side,ax) /= F) cycle
         lo = 1 ; hi = ncell
         if (side == 1) then
            lo(ax) = 1 - gw ; hi(ax) = 0
         else
            lo(ax) = ncell(ax) + 1 ; hi(ax) = ncell(ax) + gw
         endif
         ua(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = uo(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      enddo
   enddo
   call d1d1(ua, ncell, dxyz, [.false.,.false.,.false.], p_, y)
   endblock
   res = maxval(abs(y - FACTOR * src)) / maxval(abs(FACTOR * src))
   ok = gmax <= TOL_DISCRETE .and. res <= TOL_FREE .and. iters <= MAX_PCG_IT
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) print '(A)', '   D(G) p = '//trim(str(p_, .true.))//' '//label//': free ghosts against the '// &
      'direct sum '//trim(str(gmax, .true.))//', |D(G) phi - factor f| '//trim(str(res, .true.))//', CG iterations '// &
      trim(str(iters, .true.))//merge(' ok  ', ' FAIL', ok)
   call solver%destroy
   endsubroutine check_free

   subroutine check_prism_ops(label, bc, s_)
   !< E3a: the manufactured `f = div(grad u)` by the library's centred FD and FV operators of order 2 s_ (u's ghosts by
   !< the test's rules), solved by FLAIL (one 16^3 block, on rank 0), then `div(grad phi)` by the same operators from
   !< FLAIL's phi and its returned ghosts: equal to f.
   character(*), intent(in) :: label   !< Configuration label.
   integer(I4P), intent(in) :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in) :: s_      !< Half stencil (order 2 s_).
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver !< The solver.
#else
   type(poisson_fft_object)     :: solver !< The solver.
#endif
   integer(I4P), parameter   :: NC=16_I4P   !< Cells per axis.
   integer(I4P)              :: p_, sch, nbl, i, j, k, c, ierr !< Order, scheme, counters.
   real(R8P)                 :: dxyz(3), res(2), fmax !< Spacing, residual per scheme, max |f|.
   real(R8P),    allocatable :: u(:,:,:), f(:,:,:), ph(:,:,:), buf(:), out(:) !< Fields.
   integer(I4P)              :: coords(3,1) !< The one block.
   logical                   :: ok          !< Pass flag.
   character(2), parameter   :: SCHEME(2)=['FD', 'FV'] !< Schemes.

   p_ = 2 * s_
   dxyz = L / NC
   allocate(u(1-p_:NC+p_,1-p_:NC+p_,1-p_:NC+p_), f(NC,NC,NC), ph(1-p_:NC+p_,1-p_:NC+p_,1-p_:NC+p_))
   coords(:,1) = 0
   nbl = merge(1, 0, mpih%myrank == 0)
   allocate(buf(max(nbl*NC**3,1)), out(max(nbl*(NC+2*p_)**3,1)))
   res = 0._R8P
   do sch=1, 2
      u = 0._R8P
      do k=1, NC
         do j=1, NC
            do i=1, NC
               u(i,j,k) = sin(0.7_R8P * i + 0.3_R8P) * cos(0.5_R8P * j - 0.2_R8P) * cos(0.4_R8P * k + 0.1_R8P) + &
                          0.25_R8P * cos(1.1_R8P * i - 0.9_R8P * j + 0.3_R8P * k)
            enddo
         enddo
      enddo
      call rule_ghosts(u, [NC, NC, NC], bc, [.false., .false., .false.], p_)
      call div_grad(u, s_, NC, dxyz, sch, f)
      c = 0
      do k=1, NC
         do j=1, NC
            do i=1, NC
               c = c + 1
               buf(min(c,size(buf))) = f(i,j,k) / FACTOR
            enddo
         enddo
      enddo
      call solver%initialize_blocks(nblk=[1,1,1], nijk=[NC,NC,NC], ngc=p_, dxyz=dxyz, level=0_I4P, &
                                    null_xyz=[.false., .false., .false.], isolated=.false., bc=bc, dg_order=p_)
      out = 0._R8P
      call solver%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=p_, out=out)
      call solver%destroy
      if (mpih%myrank == 0) then
         c = 0
         do k=1-p_, NC+p_
            do j=1-p_, NC+p_
               do i=1-p_, NC+p_
                  c = c + 1
                  ph(i,j,k) = out(c)
               enddo
            enddo
         enddo
         call div_grad(ph, s_, NC, dxyz, sch, u(1:NC,1:NC,1:NC))
         fmax = maxval(abs(f))
         res(sch) = maxval(abs(u(1:NC,1:NC,1:NC) - f)) / fmax
      endif
   enddo
   call MPI_BCAST(res, 2, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   ok = all(res <= merge(TOL_PCG, TOL_PRISM, any(bc == Z)))
   test_passed = test_passed .and. ok
   if (mpih%myrank == 0) print '(A)', '   order '//trim(str(p_, .true.))//' '//label//': |div(grad phi) - f| / max|f| '// &
      SCHEME(1)//' '//trim(str(res(1), .true.))//', '//SCHEME(2)//' '//trim(str(res(2), .true.))//merge(' ok  ', ' FAIL', ok)
   endsubroutine check_prism_ops

   subroutine div_grad(q, s_, nc, dxyz, sch_, y)
   !< `y = div(grad q)` on the interior by the library's centred operators (FD: sch_ = 1, FV: 2) of order 2 s_: the
   !< gradient on the interior and s_ ghost layers (from q's 2 s_ ghost layers), then the divergence (E3a).
   integer(I4P), intent(in)  :: s_, nc               !< Half stencil, cells per axis.
   real(R8P),    intent(in)  :: dxyz(3)              !< Spacing.
   real(R8P),    intent(in)  :: q(1-2*s_:,1-2*s_:,1-2*s_:) !< Padded scalar.
   integer(I4P), intent(in)  :: sch_                 !< Scheme.
   real(R8P),    intent(out) :: y(1:,1:,1:)          !< Result [nc^3].
   real(R8P), allocatable    :: gr(:,:,:,:)          !< Gradient, padded by s_.
   integer(I4P)              :: ii, jj, kk           !< Counters.

   allocate(gr(3,1-s_:nc+s_,1-s_:nc+s_,1-s_:nc+s_))
   do kk=1-s_, nc+s_
      do jj=1-s_, nc+s_
         do ii=1-s_, nc+s_
            if (sch_ == 1) then
               call compute_gradient_fd_centered(s=s_, dxyz=dxyz, q=q(ii-s_:ii+s_,jj-s_:jj+s_,kk-s_:kk+s_), &
                                                 gradient=gr(:,ii,jj,kk))
            else
               call compute_gradient_fv_centered(s=s_, dxyz=dxyz, q=q(ii-s_:ii+s_,jj-s_:jj+s_,kk-s_:kk+s_), &
                                                 gradient=gr(:,ii,jj,kk))
            endif
         enddo
      enddo
   enddo
   do kk=1, nc
      do jj=1, nc
         do ii=1, nc
            if (sch_ == 1) then
               call compute_divergence_fd_centered(s=s_, dxyz=dxyz, q=gr(:,ii-s_:ii+s_,jj-s_:jj+s_,kk-s_:kk+s_), &
                                                   divergence=y(ii,jj,kk))
            else
               call compute_divergence_fv_centered(s=s_, dxyz=dxyz, q=gr(:,ii-s_:ii+s_,jj-s_:jj+s_,kk-s_:kk+s_), &
                                                   divergence=y(ii,jj,kk))
            endif
         enddo
      enddo
   enddo
   endsubroutine div_grad

   subroutine check_zero_compact
   !< The compact operator with zero-ghost faces is the one with ghost-Dirichlet faces (one ghost layer read): equal
   !< solutions, every cell and ghost, bitwise.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: s1, s2 !< The solvers.
#else
   type(poisson_fft_object)     :: s1, s2 !< The solvers.
#endif
   real(R8P), allocatable :: buf(:), o1(:), o2(:) !< Source, the two potentials.
   integer(I4P)           :: coords(3,8), nbl, b, i, ierr !< Blocks, counters.
   real(R8P)              :: dmax !< Max difference.

   nbl = 0
   do b=0, 7
      if (mod(b, mpih%procs_number) /= mpih%myrank) cycle
      nbl = nbl + 1
      coords(:,nbl) = [mod(b, 2), mod(b / 2, 2), b / 4]
   enddo
   allocate(buf(max(nbl*512,1)), o1(max(nbl*1000,1)), o2(max(nbl*1000,1)))
   do i=1, size(buf)
      buf(i) = sin(0.37_R8P * i + mpih%myrank)
   enddo
   o1 = 0._R8P ; o2 = 0._R8P
   call s1%initialize_blocks(nblk=[2,2,2], nijk=[8,8,8], ngc=1_I4P, dxyz=[0.1_R8P, 0.1_R8P, 0.1_R8P], level=0_I4P, &
                             null_xyz=[.false., .false., .false.], isolated=.false., bc=reshape([Z,Z, D,N, P,P], [2,3]))
   call s2%initialize_blocks(nblk=[2,2,2], nijk=[8,8,8], ngc=1_I4P, dxyz=[0.1_R8P, 0.1_R8P, 0.1_R8P], level=0_I4P, &
                             null_xyz=[.false., .false., .false.], isolated=.false., bc=reshape([G,G, D,N, P,P], [2,3]))
   call s1%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=1_I4P, out=o1)
   call s2%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=FACTOR, gw=1_I4P, out=o2)
   dmax = 0._R8P
   if (nbl > 0) dmax = maxval(abs(o1(1:nbl*1000) - o2(1:nbl*1000)))
   call MPI_ALLREDUCE(MPI_IN_PLACE, dmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   test_passed = test_passed .and. dmax == 0._R8P
   if (mpih%myrank == 0) print '(A)', '   compact operator, Z/Z D/N P/P against G/G D/N P/P (every cell and ghost): '// &
      trim(str(dmax, .true.))//merge(' ok  ', ' FAIL', dmax == 0._R8P)
   call s1%destroy
   call s2%destroy
   endsubroutine check_zero_compact

   subroutine rule_ghosts(u, nc, bc, null_xyz, ng)
   !< The `ng` ghost layers of the padded global array `u` by the faces' rules (the test's own, for the manufactured
   !< source and the check of the ghosts returned): periodic (and null) the wrap, Dirichlet and Neumann at the face
   !< the odd and even reflections about it, ghost-Dirichlet the odd reflection about the ghost centre, zero ghosts 0.
   integer(I4P), intent(in)    :: nc(3), bc(2,3), ng !< Cells, kinds, ghost layers.
   logical,      intent(in)    :: null_xyz(3)       !< Null directions.
   real(R8P),    intent(inout) :: u(1-ng:,1-ng:,1-ng:) !< Padded array.
   integer(I4P)                :: ax, j, lo, hi      !< Counters.

   do ax=1, 3
      do j=1, ng
         lo = 1 - j ; hi = nc(ax) + j
         if (null_xyz(ax) .or. all(bc(:,ax) == P)) then
            call put(u, ng, ax, lo, modulo(lo - 1, nc(ax)) + 1, 1._R8P)
            call put(u, ng, ax, hi, modulo(hi - 1, nc(ax)) + 1, 1._R8P)
            cycle
         endif
         select case(bc(1,ax))
         case(D) ; call put(u, ng, ax, lo, j, -1._R8P)
         case(N) ; call put(u, ng, ax, lo, j, 1._R8P)
         case(G)
            if (j == 1) then
               call put(u, ng, ax, lo, 1, 0._R8P)
            else
               call put(u, ng, ax, lo, j-1, -1._R8P)
            endif
         case default ; call put(u, ng, ax, lo, 1, 0._R8P)
         endselect
         select case(bc(2,ax))
         case(D) ; call put(u, ng, ax, hi, nc(ax)+1-j, -1._R8P)
         case(N) ; call put(u, ng, ax, hi, nc(ax)+1-j, 1._R8P)
         case(G)
            if (j == 1) then
               call put(u, ng, ax, hi, nc(ax), 0._R8P)
            else
               call put(u, ng, ax, hi, nc(ax)+2-j, -1._R8P)
            endif
         case default ; call put(u, ng, ax, hi, nc(ax), 0._R8P)
         endselect
      enddo
   enddo
   endsubroutine rule_ghosts

   subroutine put(u, ng, d_, t, src, sg)
   !< Plane t of axis d_ of the padded array u <- sg times plane src.
   integer(I4P), intent(in)    :: ng         !< Ghost layers.
   real(R8P),    intent(inout) :: u(1-ng:,1-ng:,1-ng:) !< Padded array.
   integer(I4P), intent(in)    :: d_, t, src !< Axis, indices.
   real(R8P),    intent(in)    :: sg         !< Sign.

   select case(d_)
   case(1) ; u(t,:,:) = sg * u(src,:,:)
   case(2) ; u(:,t,:) = sg * u(:,src,:)
   case(3) ; u(:,:,t) = sg * u(:,:,src)
   endselect
   endsubroutine put

   subroutine d1d1(u, nc, dxyz, null_xyz, p_, y)
   !< `y = D(G) u` on the interior: per axis D1 applied to the padded u (its ghosts set), then D1 again.
   integer(I4P), intent(in)  :: nc(3), p_           !< Cells, order.
   real(R8P),    intent(in)  :: u(1-p_:,1-p_:,1-p_:) !< Padded array.
   real(R8P),    intent(in)  :: dxyz(3)            !< Spacing.
   logical,      intent(in)  :: null_xyz(3)        !< Null directions.
   real(R8P),    intent(out) :: y(1:,1:,1:)        !< Result [nc].
   real(R8P), allocatable    :: g(:)               !< D1 along a line.
   integer(I4P)              :: ax, i, j, k, o, a, idx(3) !< Counters.

   y = 0._R8P
   allocate(g(-p_/2:p_/2))
   do ax=1, 3
      if (null_xyz(ax)) cycle
      do k=1, nc(3)
         do j=1, nc(2)
            do i=1, nc(1)
               do o=-p_/2, p_/2
                  g(o) = 0._R8P
                  do a=1, p_/2
                     idx = [i, j, k] ; idx(ax) = idx(ax) + o + a
                     g(o) = g(o) + A_D1(a,p_/2) * u(idx(1),idx(2),idx(3))
                     idx(ax) = idx(ax) - 2 * a
                     g(o) = g(o) - A_D1(a,p_/2) * u(idx(1),idx(2),idx(3))
                  enddo
               enddo
               do a=1, p_/2
                  y(i,j,k) = y(i,j,k) + A_D1(a,p_/2) * (g(a) - g(-a)) / dxyz(ax)**2
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine d1d1

   function op_label() result(lab)
   !< The operator of the current case, for the lines.
   character(:), allocatable :: lab !< Label.

   if (cfg_order == 0) then
      lab = 'compact   '
   else
      lab = 'D(G) p = '//trim(str(cfg_order, .true.))
   endif
   endfunction op_label

   subroutine refuse_wide_isolated
   !< The wide operator with isolated boundaries: initialize_blocks must stop the run.
#ifdef _FNL
   type(poisson_fft_fnl_object) :: solver !< The solver.
#else
   type(poisson_fft_object)     :: solver !< The solver.
#endif

   call solver%initialize_blocks(nblk=[2,2,2], nijk=[8,8,8], ngc=4_I4P, dxyz=[0.1_R8P, 0.1_R8P, 0.1_R8P], &
                                 level=0_I4P, null_xyz=[.false., .false., .false.], isolated=.true., dg_order=4_I4P)
   endsubroutine refuse_wide_isolated

   subroutine refusal_failed
   !< A refusal configuration ran: the test fails.
   if (mpih%myrank == 0) print '(A)', 'test_flail_bc: FAILED (the solver accepted '//trim(arg)//')'
   call mpih%finalize
   error stop 1
   endsubroutine refusal_failed
endprogram test_flail_bc
