!< Unit test of the multigrid-preconditioned CG Poisson solver on AMR grids (issue #81, P5a, SV-9).
program test_poisson_mg
!< Unit test of the multigrid-preconditioned CG Poisson solver on AMR grids (issue #81, P5a, SV-9).
!<
!< **Why this test exists.** On AMR grids self-gravity needs a Poisson solver for a composite operator (the mimetic
!< mean flux across 2:1 faces, `adam_poisson_mg_object`); exact energy conservation needs that operator symmetric,
!< accuracy needs it consistent: both are pinned here. The source and the box are SV-2's: `rho = rho0 +
!< A1 sin(k1 x) cos(k1 y) cos(k1 z) + A2 cos(k2 . x)` on the anisotropic periodic box [0,1] x [0,2] x [0,3], with the
!< exact potential in closed form.
!<
!< **What it pins.**
!<   - SV-9a, a uniform grid (one level): the multigrid against the FFT solver (the same operator there): the max
!<     difference over max |Phi| (`TOL_FFT`), and at most one FCG iteration (the V-cycle is the exact FFT solve).
!<   - SV-9b, static AMR (an octree, and a quadtree with `null_z`): the central box [L/4, 3L/4] refined
!<     one level more than the base, at two base levels. Asserted: the FCG residual below the tolerance; the symmetry
!<     defect of the assembled rows (not symmetrised by the assembly) for two pseudo-random vectors
!<     (`TOL_SYMMETRY`); every face ghost filled; the max error of Phi against the exact potential and its order
!<     between the two grids (`MIN_ORDER`); the max truncation error `|A Phi_exact - f|` bounded under refinement
!<     (growth at most `MAX_TRUNCATION_GROWTH`: an inconsistent seam coupling grows it as 1/h). Reported: the iterations, the
!<     mean convergence factor, and the face-gravity error on the fine side of the 2:1 faces with its order.
!<   - SV-12 (issue #87, P1), Dirichlet boundaries: the same source on the box without periodic boundaries, the ghost
!<     values the exact potential at the ghost centres, so the discrete problem is `A Phi = factor V (rho - rho0) - B g`.
!<     SV-12a, a uniform octree: the coarsest-level solve (the FFT solver's Dirichlet mode, an odd extension) is the
!<     whole solve, so at most one FCG iteration, and the error against the exact potential at second order between
!<     two grids. SV-12b, static AMR as SV-9b: the residual, the symmetry defect of the Dirichlet rows, every face
!<     ghost filled (the domain faces with the given values), the order of Phi, the truncation error `|A Phi_exact +
!<     B g - f|` bounded under refinement.
!< The printed numbers agree across rank counts to round-off (the reductions sum in another order) and the iteration
!< counts exactly (the driver runs np 1, 2, 3).
!<
!<   - SV-13 (issue #87, P2, `--isolated`), isolated boundaries by the symmetric James step, on the cube [-1,1]^3 with
!<     SV-6's stars (radius 1/2): the n = 1 polytrope and the uniform sphere (cell averages), against their exact
!<     potentials. SV-13a, uniform grids (16^3, 32^3): the multigrid equals the FFT solver's isolated mode (the
!<     Hockney-Eastwood padding, the same discrete free-space solution) to `TOL_FFT_ISOLATED`, each of the two solves
!<     in at most one iteration. SV-13b, base levels 2 and 3 with the star's box [-1/2, 1/2]^3 refined once: the order
!<     of Phi (`MIN_ORDER`; the uniform sphere `MIN_ORDER_SPHERE`, its surface cuts cells, as SV-6). Both: every face
!<     ghost filled, the boundary ghosts against the exact potential (reported), and, on the first grid, the symmetry
!<     of the solution map in situ, `|x . V K y - y . V K x| / |x . V K y|` for two pseudo-random sources
!<     (`TOL_SYMMETRY_ISOLATED`, set by the FCG tolerance `TOL_ISOLATED` of the four solves).
!<
!<   - FLAIL boundary kinds per face (issue #84, E1b, `--bc`): periodic, Dirichlet or Neumann at the face, and
!<     ghost-Dirichlet, mixed across axes and within an axis, on SV-2's box. Uniform grids: the multigrid against the
!<     FFT solver with the same kinds (issue #84, E1a; the same operator: `TOL_FFT`, at most one iteration), the
!<     homogeneous problem with SV-2's source. Static AMR, the corner octant [0, L/2] refined (leaves of both levels on
!<     the low faces, so the boundary rows above the coarsest level and the preconditioner's are exercised): the
!<     manufactured potential
!<     `Phi = prod_d cos(2 pi x_d / L_d + theta_d)` (periodic along every axis, generic at the faces), its Laplacian the
!<     source and its boundary values the ghost values of each kind (the value at the ghost centre or at the face, the
!<     outward derivative at the face); asserted as SV-12b (residual, symmetry, every face ghost filled, order of Phi,
!<     truncation growth) and the order of the domain-face ghosts against the manufactured potential at the ghost
!<     centres (each kind's ghost formula). Without a Dirichlet face the solution is known up to a constant: Phi is
!<     compared after removing both means.
!<
!<   - The levels inside the blocks (issue #84, E5c, `--sub`): below the coarsest tree level the V-cycle coarsens
!<     every block inside down to one cell per block (`initialize(..., root_cells=1)`) and solves that root on every
!<     rank. Each configuration is solved with the levels forced and without them (`root_cells=huge`: the coarsest
!<     tree level solved exactly), and asserted: the levels created (`lmin - lbot`, log2 of the cells per block), both
!<     solves below their stop tolerance, `max |Phi_sub - Phi_exact| / max |Phi_exact|` within `TOL_SUB` (the two
!<     stop at the same tolerance, from different preconditioners), and the iterations at most `MAX_SUB_RATIO` times
!<     the exact ones or `MAX_SUB_UNIFORM` (FCG converges with a wrong preconditioner too: the bound is what sees a
!<     broken level). Uniform and AMR, octree
!<     and quadtree, periodic, every bounded kind (BC_B), mixed periodicity (BC_A), isolated.
!<
!<   - The FLAIL facade with boundary kinds (issue #84, E1c, `--flail`): `flail%setup(..., bc)` and `solve` against the
!<     solver it chooses, called directly, bitwise (every cell and ghost of phi): the FFT solver on a uniform grid, the
!<     multigrid on AMR with inhomogeneous boundary values (`ghost`). A vector whose components take different kinds (a
!<     PEC face normal to z: Dirichlet for the tangential components, Neumann for the normal one) takes one FLAIL per
!<     component: the three instances, solved interleaved, equal three solves each by a FLAIL of its own, bitwise.
!<     `--refuse-flail-ghost-fft`: boundary values on the FFT path; `--refuse-flail-isolated-bc`: isolated with `bc`.
!<     The wide operator D(G) (E2a): `setup(..., dg_order=4)` on a uniform grid against the FFT solver with it, bitwise;
!<     `--refuse-flail-wide-amr`: the wide operator on AMR (not defined across 2:1 faces yet, E7).
!<
!< `--mixed` builds a grid periodic along x only and gives no `bc` (every active direction periodic or none then);
!< `--refuse-bc-mismatch` gives the same grid a periodic y in `bc`; `--refuse-bc-mixed-ghost` an axis with a
!< ghost-Dirichlet face and a Dirichlet one: the solver must stop (exit status 1). `--refuse-boundary-refined` builds
!< the isolated cube at level 1 with one octant refined, a leaf above the coarsest level on the boundary: the isolated
!< solver must stop.
!<
!< Usage: `mpirun -np N exe/test_poisson_mg [--isolated | --bc | --sub | --flail | --mixed | --refuse-bc-mismatch |
!< --refuse-bc-mixed-ghost | --refuse-boundary-refined | --refuse-flail-ghost-fft | --refuse-flail-isolated-bc |
!< --refuse-flail-wide-amr]`.

use :: adam_common_library,     only : realm_object
use :: adam_flail_object,       only : flail_object, FLAIL_BACKEND_AUTO
use :: adam_flail_parameters,   only : FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, FLAIL_BC_DIRICHLET_GHOST
use :: adam_parameters,         only : BC_PERIODIC, TO_BE_REFINED
use :: adam_tree_node_object,   only : tree_node_object
use :: adam_mpih_global,        only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
use :: adam_poisson_mg_object,  only : poisson_mg_object
#ifdef _FNL
! the FNL build (issue #84, E6a): the device multigrid on cuFFT, and the facades on the same backend
use :: adam_fnl_mpih_global,        only : mpih_fnl
use :: adam_fnl_poisson_fft_object, only : poisson_fft_fnl_object
use :: adam_fnl_poisson_mg_object,  only : poisson_mg_fnl_object
#define MG_T poisson_mg_fnl_object
#define FFT_T poisson_fft_fnl_object
#else
#define MG_T poisson_mg_object
#define FFT_T poisson_fft_object
#endif
use :: mpi
use :: penf,                    only : I4P, R8P, str

implicit none

real(R8P),    parameter :: PI=acos(-1._R8P)              !< Pi.
real(R8P),    parameter :: FACTOR=4._R8P * PI            !< Source factor, 4 pi G with G = 1.
real(R8P),    parameter :: RHO0=1.5_R8P                  !< Mean density.
real(R8P),    parameter :: A1=1._R8P                     !< Amplitude of the first mode.
real(R8P),    parameter :: A2=0.5_R8P                    !< Amplitude of the second mode.
real(R8P),    parameter :: L(3)=[1._R8P, 2._R8P, 3._R8P] !< Domain lengths.
real(R8P),    parameter :: SENTINEL=-1.e30_R8P           !< Unfilled phi.
real(R8P),    parameter :: TOL_FFT=1.e-8_R8P             !< Bound on |Phi_mg - Phi_fft| / max |Phi_fft|.
real(R8P),    parameter :: TOL_SYMMETRY=1.e-12_R8P       !< Bound on the operator's symmetry defect.
real(R8P),    parameter :: MIN_ORDER=1.7_R8P             !< Bound on the order of Phi against the exact one.
integer(I4P), parameter :: MAX_ITERATIONS=60_I4P         !< FCG iterations allowed.
real(R8P),    parameter :: MAX_TRUNCATION_GROWTH=1.5_R8P !< Bound on max truncation (fine grid) / (coarse grid).
real(R8P),    parameter :: A_SPHERE=0.5_R8P              !< Radius of the isolated stars (SV-13).
integer(I4P), parameter :: SUB=8_I4P                     !< Samples per direction of a cell average (uniform sphere).
integer(I4P), parameter :: PROFILE_POLYTROPE=1_I4P       !< The n = 1 polytrope.
integer(I4P), parameter :: PROFILE_UNIFORM=2_I4P         !< The uniform sphere.
real(R8P),    parameter :: TOL_ISOLATED=1.e-12_R8P       !< FCG tolerance of the isolated solves.
real(R8P),    parameter :: TOL_FFT_ISOLATED=1.e-9_R8P    !< Bound on |Phi_mg - Phi_fft| / max |Phi_fft|, isolated.
real(R8P),    parameter :: TOL_SYMMETRY_ISOLATED=1.e-9_R8P !< Bound on the solution map's symmetry defect, isolated.
real(R8P),    parameter :: MIN_ORDER_SPHERE=1.2_R8P      !< Bound on the uniform sphere's order.
real(R8P),    parameter :: TOL_SUB=1.e-7_R8P             !< Bound on |Phi_sub - Phi_exact| / max |Phi_exact| (E5c).
real(R8P),    parameter :: MAX_SUB_RATIO=1.5_R8P         !< Bound on the iterations with the levels inside the blocks over
                                                         !< those without (E5c; FCG converges with any preconditioner).
integer(I4P), parameter :: MAX_SUB_UNIFORM=20_I4P        !< The same bound where the exact coarsest needs 1 iteration.
integer(I4P), parameter :: P_=FLAIL_BC_PERIODIC          !< Periodic (the configurations below).
integer(I4P), parameter :: D_=FLAIL_BC_DIRICHLET         !< Dirichlet at the face.
integer(I4P), parameter :: N_=FLAIL_BC_NEUMANN           !< Neumann at the face.
integer(I4P), parameter :: G_=FLAIL_BC_DIRICHLET_GHOST   !< Dirichlet at the ghost centre.
integer(I4P), parameter :: BC_A(2,3)=reshape([P_, P_, D_, N_, N_, D_], [2, 3]) !< Mixed periodicity, mixed axes.
integer(I4P), parameter :: BC_B(2,3)=reshape([D_, D_, N_, N_, G_, G_], [2, 3]) !< Every bounded kind, one per axis.
integer(I4P), parameter :: BC_C(2,3)=reshape([P_, P_, N_, N_, N_, N_], [2, 3]) !< Singular: periodic and Neumann.
integer(I4P), parameter :: BC_D(2,3)=reshape([N_, D_, P_, P_, P_, P_], [2, 3]) !< 2-D: N/D along x, periodic y.
integer(I4P), parameter :: BC_PEC_XY(2,3)=reshape([P_, P_, P_, P_, D_, D_], [2, 3]) !< PEC z faces, tangential comps.
real(R8P),    parameter :: THETA(3)=[0.3_R8P, 0.7_R8P, 1.1_R8P] !< Phases of the manufactured potential.
logical                 :: test_passed                   !< Aggregate pass flag.
logical                 :: bc_suite                      !< Run the boundary-kind suite (issue #84, E1b).
logical                 :: refuse_bc_mismatch            !< Run the bc-against-grid refusal only.
logical                 :: refuse_bc_ghost               !< Run the mixed ghost-Dirichlet refusal only.
logical                 :: mixed                         !< Run the mixed-periodicity-without-bc refusal only.
logical                 :: flail_suite                   !< Run the facade suite (issue #84, E1c).
logical                 :: sub_suite                     !< Run the levels-inside-the-blocks suite (issue #84, E5c).
logical                 :: refuse_flail_ghost            !< Run the facade's ghost-on-FFT refusal only.
logical                 :: refuse_flail_isolated         !< Run the facade's isolated-with-bc refusal only.
logical                 :: refuse_flail_wide             !< Run the facade's wide-operator-on-AMR refusal only.
logical                 :: isolated                      !< Run the isolated suite (SV-13).
logical                 :: refuse_boundary               !< Run the boundary-refined refusal only.
character(256)          :: arg                           !< Command argument.
integer(I4P)            :: a                             !< Counter.

bc_suite = .false. ; refuse_bc_mismatch = .false. ; refuse_bc_ghost = .false. ; mixed = .false.
flail_suite = .false. ; refuse_flail_ghost = .false. ; refuse_flail_isolated = .false. ; refuse_flail_wide = .false.
sub_suite = .false.
isolated = .false. ; refuse_boundary = .false.
do a=1, command_argument_count()
   call get_command_argument(a, arg)
   if (trim(arg) == '--bc') bc_suite = .true.
   if (trim(arg) == '--mixed') mixed = .true.
   if (trim(arg) == '--flail') flail_suite = .true.
   if (trim(arg) == '--sub') sub_suite = .true.
   if (trim(arg) == '--refuse-flail-ghost-fft') refuse_flail_ghost = .true.
   if (trim(arg) == '--refuse-flail-isolated-bc') refuse_flail_isolated = .true.
   if (trim(arg) == '--refuse-flail-wide-amr') refuse_flail_wide = .true.
   if (trim(arg) == '--refuse-bc-mismatch') refuse_bc_mismatch = .true.
   if (trim(arg) == '--refuse-bc-mixed-ghost') refuse_bc_ghost = .true.
   if (trim(arg) == '--isolated') isolated = .true.
   if (trim(arg) == '--refuse-boundary-refined') refuse_boundary = .true.
enddo
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)
#ifdef _FNL
call mpih_fnl%initialize(do_mpi_init=.false., do_device_init=.true.)
#endif
test_passed = .true.
if (refuse_flail_wide) then
   call refuse_flail_wide_amr
   if (mpih%myrank == 0) print '(A)', 'test_poisson_mg: FAILED (the facade accepted the wide operator on AMR)'
   call mpih%finalize
   error stop 1
endif
if (refuse_flail_ghost .or. refuse_flail_isolated) then
   call refuse_flail(ghost_on_fft=refuse_flail_ghost)
   if (mpih%myrank == 0) print '(A)', 'test_poisson_mg: FAILED (the facade accepted an inconsistent call)'
   call mpih%finalize
   error stop 1
endif
if (sub_suite) then
   call check_sub(label='E5c, uniform 3-D octree, periodic,       ', ratio=8_I4P, nk=8_I4P, level=2_I4P, amr=.false.)
   call check_sub(label='E5c, AMR 3-D octree, periodic,           ', ratio=8_I4P, nk=8_I4P, level=2_I4P, amr=.true.)
   call check_sub(label='E5c, AMR 2-D quadtree, periodic,         ', ratio=4_I4P, nk=1_I4P, level=2_I4P, amr=.true.)
   call check_sub(label='E5c, AMR 3-D octree, D/D N/N G/G,        ', ratio=8_I4P, nk=8_I4P, level=2_I4P, amr=.true., &
                  bc=BC_B)
   call check_sub(label='E5c, AMR 3-D octree, P/P D/N N/D,        ', ratio=8_I4P, nk=8_I4P, level=2_I4P, amr=.true., &
                  bc=BC_A)
   call check_sub(label='E5c, AMR 3-D octree, isolated (James),   ', ratio=8_I4P, nk=8_I4P, level=2_I4P, amr=.true., &
                  isolated=.true.)
   if (mpih%myrank == 0) then
      if (test_passed) then
         print '(A)', 'test_poisson_mg: PASSED'
      else
         print '(A)', 'test_poisson_mg: FAILED'
      endif
   endif
   call mpih%finalize
   if (.not.test_passed) error stop 1
   stop
endif
if (flail_suite) then
   call check_flail_forward(label='E1c, facade, uniform 3-D, P/P D/N N/D (FFT),      ', bc=BC_A, ratio=8_I4P, &
                            nk=8_I4P, level=1_I4P, amr=.false.)
   call check_flail_forward(label='E1c, facade, AMR 3-D, D/D N/N G/G (multigrid),    ', bc=BC_B, ratio=8_I4P, &
                            nk=8_I4P, level=2_I4P, amr=.true.)
   call check_flail_forward(label='E1c, facade, AMR 2-D, N/D P/P (multigrid),        ', bc=BC_D, ratio=4_I4P, &
                            nk=1_I4P, level=2_I4P, amr=.true.)
   call check_flail_forward(label='E2a, facade, uniform 3-D, P/P D/N N/D, D(G) p = 4,  ', bc=BC_A, ratio=8_I4P, &
                            nk=8_I4P, level=1_I4P, amr=.false., dg_order=4_I4P)
   call check_flail_vector(label='E1c, facade, PEC z faces, 3 instances, uniform (FFT),', amr=.false.)
   call check_flail_vector(label='E1c, facade, PEC z faces, 3 instances, AMR (multigrid),', amr=.true.)
   if (mpih%myrank == 0) then
      if (test_passed) then
         print '(A)', 'test_poisson_mg: PASSED'
      else
         print '(A)', 'test_poisson_mg: FAILED'
      endif
   endif
   call mpih%finalize
   if (.not.test_passed) error stop 1
   stop
endif
if (mixed) then
   call refuse_mixed
   if (mpih%myrank == 0) print '(A)', 'test_poisson_mg: FAILED (the solver accepted mixed boundaries without bc)'
   call mpih%finalize
   error stop 1
endif
if (refuse_bc_mismatch .or. refuse_bc_ghost) then
   call refuse_bc(mismatch=refuse_bc_mismatch)
   if (mpih%myrank == 0) print '(A)', 'test_poisson_mg: FAILED (the solver accepted inconsistent boundary kinds)'
   call mpih%finalize
   error stop 1
endif
if (refuse_boundary) then
   call refuse_boundary_refined
   if (mpih%myrank == 0) print '(A)', 'test_poisson_mg: FAILED (the isolated solver accepted a refined boundary leaf)'
   call mpih%finalize
   error stop 1
endif
if (bc_suite) then
   call check_bc_uniform(label='E1b, uniform 3-D, P/P D/N N/D,', bc=BC_A, ratio=8_I4P, nk=8_I4P, level=1_I4P)
   call check_bc_uniform(label='E1b, uniform 3-D, D/D N/N G/G,', bc=BC_B, ratio=8_I4P, nk=8_I4P, level=1_I4P)
   call check_bc_uniform(label='E1b, uniform 3-D, P/P N/N N/N,', bc=BC_C, ratio=8_I4P, nk=8_I4P, level=1_I4P)
   call check_bc_uniform(label='E1b, uniform 2-D, N/D P/P,    ', bc=BC_D, ratio=4_I4P, nk=1_I4P, level=2_I4P)
   call check_bc_amr(label='E1b, AMR 3-D, P/P D/N N/D,', bc=BC_A, ratio=8_I4P, nk=8_I4P, levels=[2_I4P, 3_I4P])
   call check_bc_amr(label='E1b, AMR 3-D, D/D N/N G/G,', bc=BC_B, ratio=8_I4P, nk=8_I4P, levels=[2_I4P, 3_I4P])
   call check_bc_amr(label='E1b, AMR 3-D, P/P N/N N/N,', bc=BC_C, ratio=8_I4P, nk=8_I4P, levels=[2_I4P, 3_I4P])
   call check_bc_amr(label='E1b, AMR 2-D, N/D P/P,    ', bc=BC_D, ratio=4_I4P, nk=1_I4P, levels=[2_I4P, 3_I4P])
   if (mpih%myrank == 0) then
      if (test_passed) then
         print '(A)', 'test_poisson_mg: PASSED'
      else
         print '(A)', 'test_poisson_mg: FAILED'
      endif
   endif
   call mpih%finalize
   if (.not.test_passed) error stop 1
   stop
endif
if (isolated) then
   call check_isolated(label='SV-13a, isolated, uniform, polytrope n=1,', prof=PROFILE_POLYTROPE, &
                       levels=[1_I4P, 2_I4P], amr=.false.)
   call check_isolated(label='SV-13b, isolated, AMR, polytrope n=1,    ', prof=PROFILE_POLYTROPE, &
                       levels=[2_I4P, 3_I4P], amr=.true.)
   call check_isolated(label='SV-13b, isolated, AMR, uniform sphere,   ', prof=PROFILE_UNIFORM, &
                       levels=[2_I4P, 3_I4P], amr=.true.)
   if (mpih%myrank == 0) then
      if (test_passed) then
         print '(A)', 'test_poisson_mg: PASSED'
      else
         print '(A)', 'test_poisson_mg: FAILED'
      endif
   endif
   call mpih%finalize
   if (.not.test_passed) error stop 1
   stop
endif
call check_uniform(label='SV-9a, uniform 3-D octree,  ', ratio=8_I4P, nk=8_I4P, level=2_I4P)
call check_uniform(label='SV-9a, uniform 2-D quadtree,', ratio=4_I4P, nk=1_I4P, level=2_I4P)
call check_amr(label='SV-9b, AMR 3-D octree,      ', ratio=8_I4P, nk=8_I4P, levels=[2_I4P, 3_I4P])
call check_amr(label='SV-9b, AMR 2-D quadtree,    ', ratio=4_I4P, nk=1_I4P, levels=[2_I4P, 3_I4P])
call check_dirichlet(label='SV-12a, Dirichlet, uniform 3-D octree,', ratio=8_I4P, nk=8_I4P, levels=[1_I4P, 2_I4P], &
                     sphere=.false.)
call check_dirichlet(label='SV-12b, Dirichlet, AMR 3-D octree,    ', ratio=8_I4P, nk=8_I4P, levels=[2_I4P, 3_I4P], &
                     sphere=.true.)
call check_dirichlet(label='SV-12b, Dirichlet, AMR 2-D quadtree,  ', ratio=4_I4P, nk=1_I4P, levels=[2_I4P, 3_I4P], &
                     sphere=.true.)
if (mpih%myrank == 0) then
   if (test_passed) then
      print '(A)', 'test_poisson_mg: PASSED'
   else
      print '(A)', 'test_poisson_mg: FAILED'
   endif
endif
call mpih%finalize
if (.not.test_passed) error stop 1

contains
   subroutine use_backend(fl)
   !< The facade on the build's backend: the FNL build gives it the cuFFT and device multigrid molds (#84 E6a).
   type(flail_object), intent(inout) :: fl !< FLAIL.

#ifdef _FNL
   allocate(poisson_fft_fnl_object :: fl%fft)
   allocate(poisson_mg_fnl_object :: fl%mg_mold)
#else
   if (allocated(fl%mg_mold)) deallocate(fl%mg_mold)
#endif
   endsubroutine use_backend

   subroutine refuse_mixed
   !< A uniform 3-D grid periodic along x only and no bc: `initialize` must stop the run.
   type(realm_object)      :: realm !< Realm.
   type(MG_T) :: mg    !< Multigrid solver.
   real(R8P), allocatable  :: q(:,:,:,:,:) !< Source.

   call build(realm, 8_I4P, 8_I4P, 1_I4P, sphere=.false., q=q, periodic=.false.)
   call realm%adam%grid%set_bc_type(bc_type=[BC_PERIODIC, BC_PERIODIC, 0_I4P, 0_I4P, 0_I4P, 0_I4P])
   call mg%initialize(grid=realm%adam%grid, field=realm%adam%field)
   endsubroutine refuse_mixed

   subroutine refuse_bc(mismatch)
   !< A uniform 3-D grid periodic along x only, with inconsistent kinds: `initialize` must stop the run. `mismatch`:
   !< bc periodic along y, which the grid bounds; else a ghost-Dirichlet face and a Dirichlet one along y.
   logical, intent(in)     :: mismatch !< The bc-against-grid case, else the mixed ghost one.
   type(realm_object)      :: realm !< Realm.
   type(MG_T) :: mg    !< Multigrid solver.
   real(R8P), allocatable  :: q(:,:,:,:,:) !< Source.

   call build(realm, 8_I4P, 8_I4P, 1_I4P, sphere=.false., q=q, periodic=.false.)
   call realm%adam%grid%set_bc_type(bc_type=[BC_PERIODIC, BC_PERIODIC, 0_I4P, 0_I4P, 0_I4P, 0_I4P])
   if (mismatch) then
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, bc=reshape([P_, P_, P_, P_, D_, D_], [2, 3]))
   else
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, bc=reshape([P_, P_, G_, D_, D_, D_], [2, 3]))
   endif
   endsubroutine refuse_bc

   subroutine refuse_boundary_refined
   !< The isolated cube at level 1 with the octant x, y, z > 0 refined: a leaf at level 2 on the boundary, `initialize`
   !< must stop the run.
   type(realm_object)      :: realm !< Realm.
   type(MG_T) :: mg    !< Multigrid solver.
   real(R8P), allocatable  :: q(:,:,:,:,:) !< Source.

   call build_star(realm, 1_I4P, refine=.true., boundary=.true., prof=PROFILE_POLYTROPE, q=q)
   call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, isolated=.true.)
   endsubroutine refuse_boundary_refined

   subroutine check_isolated(label, prof, levels, amr)
   !< SV-13: isolated boundaries (the James step) on two grids, uniform (against the FFT solver) or with the star's box
   !< refined once; the order against the exact potential, and the symmetry of the solution map on the first grid.
   character(*), intent(in) :: label     !< Label.
   integer(I4P), intent(in) :: prof      !< Source.
   integer(I4P), intent(in) :: levels(2) !< Base levels.
   logical,      intent(in) :: amr       !< Refine the star's box.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi(:,:,:,:), phi_fft(:,:,:,:) !< Fields.
   real(R8P), allocatable   :: x(:), y(:), ux(:), uy(:) !< Sources and solutions of the symmetry check.
   real(R8P)                :: e_phi(2), h(2), e_ghost, sym, order, diff, scale_, dots(2), min_ord !< Measures.
   integer(I4P)             :: unfilled, m, i, nl, nleaf, n, ni, nj, nk, ierr !< Counters.
   logical                  :: ok        !< Pass flag.
   character(:), allocatable :: note     !< The extra measures of a grid.

   ok = .true.
   do m=1, 2
      block
      type(realm_object)       :: realm !< Realm (a fresh one per grid).
      type(MG_T)  :: mg    !< Solver.
      type(FFT_T) :: fft   !< FFT solver (uniform grids).
      call build_star(realm, levels(m), refine=amr, boundary=.false., prof=prof, q=q)
      allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
      phi = SENTINEL
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, tolerance=TOL_ISOLATED, &
                         max_iterations=MAX_ITERATIONS, isolated=.true.)
      call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi)
      ok = ok .and. mg%residual <= mg%tolerance_used
      call measure_star(realm, phi, prof, e_phi(m), e_ghost, unfilled)
      ok = ok .and. unfilled == 0
      h(m) = realm%adam%grid%cell_dxyz(1,levels(m))
      nleaf = realm%adam%field%blocks_number
      call MPI_ALLREDUCE(MPI_IN_PLACE, nleaf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      note = ', leaves '//trim(str(nleaf, .true.))//': iterations '//trim(str(mg%iterations0, .true.))//' + '// &
             trim(str(mg%iterations, .true.))
      n = realm%adam%field%blocks_number ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
      nk = realm%adam%grid%nk
      if (.not.amr) then ! the multigrid against the FFT solver's isolated mode
         allocate(phi_fft, mold=phi)
         phi_fft = SENTINEL
         call fft%initialize(grid=realm%adam%grid, field=realm%adam%field, isolated=.true.)
         call fft%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_fft)
         diff = 0._R8P ; scale_ = 0._R8P
         if (n > 0) then
            diff = maxval(abs(phi(1:ni,1:nj,1:nk,1:n) - phi_fft(1:ni,1:nj,1:nk,1:n)))
            scale_ = maxval(abs(phi_fft(1:ni,1:nj,1:nk,1:n)))
         endif
         call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, scale_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         ok = ok .and. diff / scale_ <= TOL_FFT_ISOLATED .and. mg%iterations0 <= 1 .and. mg%iterations <= 1
         note = note//', |Phi_mg - Phi_fft| / max|Phi| '//trim(str(diff / scale_, .true.))
         call fft%destroy
         deallocate(phi_fft)
      endif
      if (m == 1) then ! the symmetry of the solution map, in situ
         nl = n * ni * nj * nk
         allocate(x(max(nl,1)), y(max(nl,1)), ux(max(nl,1)), uy(max(nl,1)))
         do i=1, nl
            x(i) = sin(1.3_R8P * i + 0.7_R8P * mpih%myrank) ; y(i) = cos(2.1_R8P * i - 0.4_R8P * mpih%myrank)
         enddo
         call mg%solve_packed(field=realm%adam%field, buf=x, factor=1._R8P, phi=phi)
         ux = mg%u(1:max(nl,1))
         call mg%solve_packed(field=realm%adam%field, buf=y, factor=1._R8P, phi=phi)
         uy = mg%u(1:max(nl,1))
         dots = 0._R8P
         do i=1, nl
            dots(1) = dots(1) + mg%vol(i) * x(i) * uy(i)
            dots(2) = dots(2) + mg%vol(i) * y(i) * ux(i)
         enddo
         call MPI_ALLREDUCE(MPI_IN_PLACE, dots, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
         sym = abs(dots(1) - dots(2)) / abs(dots(1))
         ok = ok .and. sym <= TOL_SYMMETRY_ISOLATED
         note = note//', symmetry of V K '//trim(str(sym, .true.))
         deallocate(x, y, ux, uy)
      endif
      if (mpih%myrank == 0) print '(A)', '   '//label//' base level '//trim(str(levels(m), .true.))//note// &
         ', |Phi - Phi_exact| '//trim(str(e_phi(m), .true.))//', boundary ghosts '//trim(str(e_ghost, .true.))// &
         ', unfilled '//trim(str(unfilled, .true.))
      call mg%destroy
      deallocate(phi)
      endblock
   enddo
   order = log(e_phi(1) / e_phi(2)) / log(h(1) / h(2))
   min_ord = MIN_ORDER
   if (prof == PROFILE_UNIFORM) min_ord = MIN_ORDER_SPHERE
   ok = ok .and. order >= min_ord
   if (mpih%myrank == 0) print '(A)', '   '//label//' order of Phi '//trim(str(order, .true.))//' (min '// &
      trim(str(min_ord, .true.))//')'//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   endsubroutine check_isolated

   subroutine measure_star(realm, phi, prof, e_phi, e_ghost, unfilled)
   !< Max error of Phi against the star's exact potential on the leaves; on the first ghost layer beyond a domain face
   !< (the James ghost values); unfilled face ghosts.
   type(realm_object), intent(in)  :: realm    !< Realm.
   real(R8P),          intent(in)  :: phi(1-realm%adam%grid%ngc:,1-realm%adam%grid%ngc:,1-realm%adam%grid%ngc:,1:) !< Phi.
   integer(I4P),       intent(in)  :: prof     !< Source.
   real(R8P),          intent(out) :: e_phi    !< Max |Phi - Phi_exact| on the leaves.
   real(R8P),          intent(out) :: e_ghost  !< Max |Phi - Phi_exact| on the boundary ghosts.
   integer(I4P),       intent(out) :: unfilled !< Unfilled face ghosts.
   real(R8P)                       :: xyz(3), p(3) !< Points.
   integer(I4P)                    :: b, i, j, k, d, s, ijk(3), nb(3), n(3), ierr !< Counters.

   associate(f=>realm%adam%field, g=>realm%adam%grid)
   n = [g%ni, g%nj, g%nk]
   e_phi = 0._R8P ; e_ghost = 0._R8P ; unfilled = 0
   do b=1, f%blocks_number
      do k=1, g%nk
         do j=1, g%nj
            do i=1, g%ni
               xyz = [f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)]
               e_phi = max(e_phi, abs(phi(i,j,k,b) - phi_sphere(xyz, prof)))
               do d=1, 3
                  do s=-1, 1, 2
                     ijk = [i, j, k] ; nb = ijk ; nb(d) = nb(d) + s
                     if (nb(d) >= 1 .and. nb(d) <= n(d)) cycle
                     if (phi(nb(1),nb(2),nb(3),b) == SENTINEL) unfilled = unfilled + 1
                     p = xyz ; p(d) = p(d) + s * f%dxyz(d,b)
                     if (p(d) > g%domain_emin(d) .and. p(d) < g%domain_emax(d)) cycle
                     e_ghost = max(e_ghost, abs(phi(nb(1),nb(2),nb(3),b) - phi_sphere(p, prof)))
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_phi, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_ghost, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   endsubroutine measure_star

   subroutine build_star(realm, level, refine, boundary, prof, q)
   !< The isolated cube [-1,1]^3 (octree, 8^3 cells per block) uniform at `level`; with `refine` the star's box
   !< [-1/2, 1/2]^3 refined once more (block-aligned from level 2 on), or with `boundary` the octant x, y, z > 0.
   type(realm_object),     intent(inout) :: realm    !< Realm.
   integer(I4P),           intent(in)    :: level    !< Base level.
   logical,                intent(in)    :: refine   !< Refine.
   logical,                intent(in)    :: boundary !< Refine the octant (on the boundary), not the star's box.
   integer(I4P),           intent(in)    :: prof     !< Source.
   real(R8P), allocatable, intent(inout) :: q(:,:,:,:,:) !< The source field.
   character(:), allocatable             :: filename !< Input file.
   integer(I4P)                          :: ngc, ni, nj, nk, b, i, j, k, ierr !< Counters.
   logical                               :: mark     !< Refine the block.
   real(R8P)                             :: c(3)     !< Block centre.
   type(tree_node_object), pointer       :: node     !< A tree node.

   filename = 'test_poisson_mg-cube-l'//trim(str(level, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, nk=8_I4P, ratio=8_I4P, level=level+1, cube=.true.)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=8._R8P, nv=1_I4P)
   ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj ; nk = realm%adam%grid%nk
   if (allocated(q)) deallocate(q)
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   q = 0._R8P
   if (level > 0) call realm%adam%refine_uniform(refinement_levels=level, q=q, do_mpi_redistribute=.true., &
                                                 do_blocks_reorder=.false.)
   if (refine) then
      do b=1, realm%adam%field%blocks_number
         c = 0.5_R8P * (realm%adam%field%emin(:,b) + realm%adam%field%emax(:,b))
         if (boundary) then
            mark = all(c > 0._R8P)
         else
            mark = maxval(abs(c)) < A_SPHERE
         endif
         if (mark) then
            node => realm%adam%tree%node(code=realm%adam%field%code(b))
            node%refinement_needed = TO_BE_REFINED
         endif
      enddo
      call realm%adam%amr_update(q=q, is_marked_by_tree=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
   endif
   do b=1, realm%adam%field%blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(1,i,j,k,b) = rho_sphere([realm%adam%field%x_cell(i,b), realm%adam%field%y_cell(j,b), &
                                          realm%adam%field%z_cell(k,b)], prof, realm%adam%field%dxyz(:,b))
            enddo
         enddo
      enddo
   enddo
   endsubroutine build_star

   pure function rho_sphere(xyz, prof, dxyz) result(r)
   !< The isolated source: the n = 1 polytrope at the point, or the uniform sphere's cell average (SUB^3 samples), as
   !< test_poisson_fft (SV-6).
   real(R8P),    intent(in) :: xyz(3)  !< Cell centre.
   integer(I4P), intent(in) :: prof    !< Source.
   real(R8P),    intent(in) :: dxyz(3) !< Cell spacing.
   real(R8P)                :: r       !< Density.
   real(R8P)                :: kr      !< k r.
   real(R8P)                :: p(3)    !< Sample point.
   integer(I4P)             :: i, j, k !< Counters.

   if (prof == PROFILE_POLYTROPE) then
      kr = PI / A_SPHERE * norm2(xyz)
      r = 0._R8P
      if (kr < PI) then
         r = 1._R8P
         if (kr > 0._R8P) r = sin(kr) / kr
      endif
   else
      r = 0._R8P
      do k=1, SUB
         do j=1, SUB
            do i=1, SUB
               p = xyz + dxyz * ([real(i, R8P), real(j, R8P), real(k, R8P)] - 0.5_R8P) / SUB - 0.5_R8P * dxyz
               if (norm2(p) < A_SPHERE) r = r + 1._R8P
            enddo
         enddo
      enddo
      r = r / SUB**3
   endif
   endfunction rho_sphere

   pure function phi_sphere(xyz, prof) result(p)
   !< The exact potential of the isolated source (G = 1).
   real(R8P),    intent(in) :: xyz(3) !< Point.
   integer(I4P), intent(in) :: prof   !< Source.
   real(R8P)                :: p      !< Potential.
   real(R8P)                :: rr, kk !< Radius, k.
   real(R8P)                :: mass   !< Mass.

   rr = norm2(xyz)
   if (prof == PROFILE_POLYTROPE) then
      kk = PI / A_SPHERE
      mass = 4._R8P * A_SPHERE**3 / PI
      if (rr < A_SPHERE) then
         p = -mass / A_SPHERE
         if (rr > 0._R8P) then
            p = p - FACTOR / kk**2 * sin(kk * rr) / (kk * rr)
         else
            p = p - FACTOR / kk**2
         endif
      else
         p = -mass / rr
      endif
   else
      mass = 4._R8P * PI * A_SPHERE**3 / 3._R8P
      if (rr < A_SPHERE) then
         p = -2._R8P * PI * (A_SPHERE**2 - rr**2 / 3._R8P)
      else
         p = -mass / rr
      endif
   endif
   endfunction phi_sphere

   subroutine refuse_flail(ghost_on_fft)
   !< The facade, inconsistent calls: `ghost_on_fft`, boundary values on the FFT path (a uniform grid, bounded along z);
   !< else isolated boundaries with bc. Each must stop the run.
   logical, intent(in)     :: ghost_on_fft !< The ghost case, else the isolated one.
   type(realm_object)      :: realm !< Realm.
   type(flail_object)      :: flail !< FLAIL.
   real(R8P), allocatable  :: q(:,:,:,:,:), phi(:,:,:,:), ghost(:) !< Fields, boundary values.

   call build(realm, 8_I4P, 8_I4P, 1_I4P, sphere=.false., q=q, periodic=.false.)
   call use_backend(flail)
   if (ghost_on_fft) then
      call set_grid_bc(realm, BC_PEC_XY)
      call flail%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.false., backend=FLAIL_BACKEND_AUTO, &
                       bc=BC_PEC_XY)
      allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)), ghost(1))
      ghost = 0._R8P
      call flail%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi, ghost=ghost)
   else
      call flail%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.true., backend=FLAIL_BACKEND_AUTO, &
                       bc=BC_B)
   endif
   endsubroutine refuse_flail

   subroutine refuse_flail_wide_amr
   !< The facade with the wide operator on AMR (the corner refined): setup must stop the run.
   type(realm_object)      :: realm !< Realm.
   type(flail_object)      :: flail !< FLAIL.
   real(R8P), allocatable  :: q(:,:,:,:,:) !< Source.

   call build(realm, 8_I4P, 8_I4P, 2_I4P, sphere=.true., q=q, periodic=.false., corner=.true.)
   call set_grid_bc(realm, BC_B)
   call use_backend(flail)
   call flail%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.false., backend=FLAIL_BACKEND_AUTO, &
                    bc=BC_B, dg_order=4_I4P)
   endsubroutine refuse_flail_wide_amr

   subroutine check_sub(label, ratio, nk, level, amr, bc, isolated)
   !< Issue #84, E5c: one configuration solved with the levels inside the blocks forced (`root_cells=1`) and without
   !< them (`root_cells=huge`); see the module.
   character(*), intent(in)           :: label    !< Label.
   integer(I4P), intent(in)           :: ratio    !< Tree ratio.
   integer(I4P), intent(in)           :: nk       !< Cells per block along z.
   integer(I4P), intent(in)           :: level    !< Base level.
   logical,      intent(in)           :: amr      !< Refine the central box (the star's, isolated).
   integer(I4P), intent(in), optional :: bc(2,3)  !< Kinds per face (default periodic).
   logical,      intent(in), optional :: isolated !< Isolated boundaries (the star's cube).
   type(realm_object)                 :: realm    !< Realm.
   type(MG_T)            :: mg(2)    !< Solvers: the levels forced, none.
   real(R8P), allocatable             :: q(:,:,:,:,:), phi(:,:,:,:,:) !< Source, the two potentials.
   real(R8P)                          :: diff(2)  !< max |Phi_sub - Phi_exact|, max |Phi_exact|.
   integer(I4P)                       :: its(2), nsub_expected, s, b, i, j, k, ierr !< Iterations, levels, counters.
   logical                            :: iso, ok  !< Isolated, pass flag.

   iso = .false. ; if (present(isolated)) iso = isolated
   if (iso) then
      call build_star(realm, level, refine=amr, boundary=.false., prof=PROFILE_POLYTROPE, q=q)
   elseif (present(bc)) then
      call build(realm, ratio, nk, level, sphere=amr, q=q, periodic=.false.)
      call set_grid_bc(realm, bc)
   else
      call build(realm, ratio, nk, level, sphere=amr, q=q)
   endif
   associate(f=>realm%adam%field, g=>realm%adam%grid)
   allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5),2))
   phi = 0._R8P
   do s=1, 2
      if (iso) then
         call mg(s)%initialize(grid=g, field=f, tolerance=TOL_ISOLATED, max_iterations=MAX_ITERATIONS, isolated=.true., &
                               root_cells=merge(1_I4P, huge(1_I4P), s == 1))
      elseif (present(bc)) then
         call mg(s)%initialize(grid=g, field=f, max_iterations=MAX_ITERATIONS, bc=bc, &
                               root_cells=merge(1_I4P, huge(1_I4P), s == 1))
      else
         call mg(s)%initialize(grid=g, field=f, max_iterations=MAX_ITERATIONS, root_cells=merge(1_I4P, huge(1_I4P), s == 1))
      endif
      call mg(s)%solve(field=f, q=q, iv=1_I4P, factor=FACTOR, phi=phi(:,:,:,:,s))
      its(s) = mg(s)%iterations + mg(s)%iterations0
   enddo
   diff = 0._R8P
   do b=1, f%blocks_number
      do k=1, nk
         do j=1, g%nj
            do i=1, g%ni
               diff(1) = max(diff(1), abs(phi(i,j,k,b,1) - phi(i,j,k,b,2)))
               diff(2) = max(diff(2), abs(phi(i,j,k,b,2)))
            enddo
         enddo
      enddo
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 2, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   nsub_expected = nint(log(real(g%ni, R8P)) / log(2._R8P))
   ok = mg(1)%lmin - mg(1)%lbot == nsub_expected .and. mg(2)%lbot == mg(2)%lmin
   ok = ok .and. mg(1)%residual <= mg(1)%tolerance_used .and. mg(2)%residual <= mg(2)%tolerance_used
   ok = ok .and. diff(1) <= TOL_SUB * diff(2)
   ok = ok .and. its(1) <= max(MAX_SUB_UNIFORM, nint(MAX_SUB_RATIO * its(2)))
   if (mpih%myrank == 0) print '(A)', '   '//label//' levels inside the blocks '// &
      trim(str(mg(1)%lmin - mg(1)%lbot, .true.))//', iterations '//trim(str(its(1), .true.))//' (exact coarsest '// &
      trim(str(its(2), .true.))//'), |Phi_sub - Phi_exact| / max |Phi| '//trim(str(diff(1) / diff(2)))//'  '// &
      merge('PASS', 'FAIL', ok)
   test_passed = test_passed .and. ok
   endassociate
   call mg(1)%destroy
   call mg(2)%destroy
   endsubroutine check_sub

   subroutine check_flail_forward(label, bc, ratio, nk, level, amr, dg_order)
   !< Issue #84, E1c: `flail%setup(..., bc)` and `solve` against the solver it chooses called directly, bitwise on every
   !< cell and ghost of phi; on AMR (the corner refined) with the manufactured potential's boundary values as `ghost`.
   character(*), intent(in) :: label   !< Label.
   integer(I4P), intent(in) :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in) :: ratio   !< Tree ratio.
   integer(I4P), intent(in) :: nk      !< Cells per block along z.
   integer(I4P), intent(in) :: level   !< Base level.
   logical,      intent(in) :: amr     !< Refine the corner (the multigrid path), else uniform (the FFT path).
   integer(I4P), intent(in), optional :: dg_order !< The wide D(G) of this order (FFT path), else compact.
   type(realm_object)       :: realm   !< Realm.
   type(flail_object)       :: flail   !< FLAIL.
   type(MG_T)  :: mg      !< Multigrid solver, direct.
   type(FFT_T) :: fft     !< FFT solver, direct.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi_f(:,:,:,:), phi_d(:,:,:,:), ghost(:) !< Fields, boundary values.
   real(R8P)                :: diff    !< Max |phi_flail - phi_direct|, ghosts included.
   logical                  :: ok, active(3) !< Pass flag, active directions.
   integer(I4P)             :: i, ierr !< Counter, MPI status.

   call build(realm, ratio, nk, level, sphere=amr, q=q, periodic=.false., corner=.true.)
   call set_grid_bc(realm, bc)
   active = .not.realm%adam%grid%null_xyz
   allocate(phi_f(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
   allocate(phi_d, mold=phi_f)
   phi_f = SENTINEL ; phi_d = SENTINEL
   call use_backend(flail)
   call flail%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.false., backend=FLAIL_BACKEND_AUTO, bc=bc, &
                    dg_order=dg_order)
   ok = flail%is_multigrid() .eqv. amr
   if (amr) then
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, bc=bc)
      allocate(ghost(max(mg%nbnd,1)))
      do i=1, mg%nbnd
         if (mg%bnd_kind(i) == FLAIL_BC_NEUMANN) then
            ghost(i) = sign(1_I4P, mg%bnd_dir(i)) * grad_m(mg%bnd_x(:,i), active, abs(mg%bnd_dir(i)))
         else
            ghost(i) = phi_m(mg%bnd_x(:,i), active)
         endif
      enddo
      call flail%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_f, ghost=ghost)
      call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_d, ghost=ghost)
      call mg%destroy
   else
      call fft%initialize(grid=realm%adam%grid, field=realm%adam%field, bc=bc, dg_order=dg_order)
      call flail%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_f)
      call fft%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_d)
      call fft%destroy
   endif
   diff = 0._R8P
   if (size(phi_f) > 0) diff = maxval(abs(phi_f - phi_d))
   call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ok = ok .and. diff == 0._R8P
   if (mpih%myrank == 0) print '(A)', '   '//label//' path '//trim(merge('multigrid', 'FFT      ', &
      flail%is_multigrid()))//', |Phi_flail - Phi_direct| (every cell and ghost) '//trim(str(diff, .true.))// &
      merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   call flail%destroy
   endsubroutine check_flail_forward

   subroutine check_flail_vector(label, amr)
   !< Issue #84, E1c: a vector with per-component kinds, a PEC face normal to z (x, y periodic): Dirichlet for the
   !< tangential components, Neumann for the normal one. One FLAIL per component; the three set up together and solved
   !< interleaved (z, x, y) equal three solves each by a FLAIL of its own, set up and destroyed in turn, bitwise.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: amr   !< Refine the corner (the multigrid path), else uniform (the FFT path).
   type(realm_object)       :: realm !< Realm.
   type(flail_object)       :: fl(3) !< One FLAIL per component.
   type(flail_object)       :: one   !< A FLAIL of its own.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi_v(:,:,:,:,:), phi_1(:,:,:,:) !< Fields.
   integer(I4P)             :: bcv(2,3,3) !< Kinds per component.
   real(R8P)                :: diff  !< Max difference, ghosts included.
   logical                  :: ok    !< Pass flag.
   integer(I4P)             :: c, ierr !< Counter, MPI status.
   integer(I4P), parameter  :: ORDER(3)=[3_I4P, 1_I4P, 2_I4P] !< The interleaved order.

   bcv(:,:,1) = BC_PEC_XY ; bcv(:,:,2) = BC_PEC_XY
   bcv(:,:,3) = reshape([P_, P_, P_, P_, N_, N_], [2, 3])
   call build(realm, 8_I4P, 8_I4P, merge(2_I4P, 1_I4P, amr), sphere=amr, q=q, periodic=.false., corner=.true.)
   call set_grid_bc(realm, BC_PEC_XY)
   allocate(phi_v(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5),3))
   allocate(phi_1(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
   phi_v = SENTINEL
   ok = .true.
   do c=1, 3
      call use_backend(fl(c))
      call fl(c)%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.false., backend=FLAIL_BACKEND_AUTO, &
                       bc=bcv(:,:,c))
      ok = ok .and. (fl(c)%is_multigrid() .eqv. amr)
   enddo
   do c=1, 3 ! the factor differs per component: three distinct problems from one source
      call fl(ORDER(c))%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR * ORDER(c), &
                              phi=phi_v(:,:,:,:,ORDER(c)))
   enddo
   diff = 0._R8P
   do c=1, 3
      phi_1 = SENTINEL
      call use_backend(one) ! destroy releases the molds
      call one%setup(grid=realm%adam%grid, field=realm%adam%field, isolated=.false., backend=FLAIL_BACKEND_AUTO, &
                     bc=bcv(:,:,c))
      call one%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR * c, phi=phi_1)
      call one%destroy
      if (size(phi_1) > 0) diff = max(diff, maxval(abs(phi_v(:,:,:,:,c) - phi_1)))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ok = ok .and. diff == 0._R8P
   if (mpih%myrank == 0) print '(A)', '   '//label//' |Phi_interleaved - Phi_alone| (every cell and ghost) '// &
      trim(str(diff, .true.))//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   do c=1, 3
      call fl(c)%destroy
   enddo
   endsubroutine check_flail_vector

   subroutine check_bc_uniform(label, bc, ratio, nk, level)
   !< Issue #84, E1b: the multigrid with face kinds `bc` on a uniform grid against the FFT solver with the same kinds
   !< (the same operator), the homogeneous problem with SV-2's source.
   character(*), intent(in) :: label   !< Label.
   integer(I4P), intent(in) :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in) :: ratio   !< Tree ratio.
   integer(I4P), intent(in) :: nk      !< Cells per block along z.
   integer(I4P), intent(in) :: level   !< Uniform level.
   type(realm_object)       :: realm   !< Realm.
   type(MG_T)  :: mg      !< Multigrid solver.
   type(FFT_T) :: fft     !< FFT solver.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi_mg(:,:,:,:), phi_fft(:,:,:,:) !< Fields.
   real(R8P)                :: diff, scale_ !< Max difference, max |Phi_fft|.
   logical                  :: ok      !< Pass flag.
   integer(I4P)             :: n, ni, nj, ierr !< Sizes, MPI status.

   call build(realm, ratio, nk, level, sphere=.false., q=q, periodic=.false.)
   call set_grid_bc(realm, bc)
   allocate(phi_mg(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
   allocate(phi_fft, mold=phi_mg)
   phi_mg = SENTINEL ; phi_fft = SENTINEL
   call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, max_iterations=MAX_ITERATIONS, bc=bc)
   call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_mg)
   call fft%initialize(grid=realm%adam%grid, field=realm%adam%field, bc=bc)
   call fft%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_fft)
   n = realm%adam%field%blocks_number ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   diff = 0._R8P ; scale_ = 0._R8P
   if (n > 0) then
      diff = maxval(abs(phi_mg(1:ni,1:nj,1:nk,1:n) - phi_fft(1:ni,1:nj,1:nk,1:n)))
      scale_ = maxval(abs(phi_fft(1:ni,1:nj,1:nk,1:n)))
   endif
   call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, scale_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ok = diff / scale_ <= TOL_FFT .and. mg%iterations <= 1 .and. &
        (mg%singular .eqv. .not.any(bc == FLAIL_BC_DIRICHLET .or. bc == FLAIL_BC_DIRICHLET_GHOST))
   if (mpih%myrank == 0) print '(A)', '   '//label//' |Phi_mg - Phi_fft| / max|Phi| '// &
      trim(str(diff / scale_, .true.))//', iterations '//trim(str(mg%iterations, .true.))//', singular '// &
      trim(merge('T', 'F', mg%singular))//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   call mg%destroy
   call fft%destroy
   endsubroutine check_bc_uniform

   subroutine check_bc_amr(label, bc, ratio, nk, levels)
   !< Issue #84, E1b: the face kinds `bc` on static AMR (the central box refined once over two base levels), against the
   !< manufactured potential (`phi_m`), its boundary values the ghost values of each kind.
   character(*), intent(in) :: label     !< Label.
   integer(I4P), intent(in) :: bc(2,3)   !< Kinds per face [lo/hi, axis].
   integer(I4P), intent(in) :: ratio     !< Tree ratio.
   integer(I4P), intent(in) :: nk        !< Cells per block along z.
   integer(I4P), intent(in) :: levels(2) !< Base levels.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi(:,:,:,:), x(:), y(:), ghost(:), pe(:), ye(:), yb(:) !< Fields, vectors.
   real(R8P)                :: e_phi(2), h(2), trunc(2), sym, order, vsum, pmean, umean, v, xyz(3) !< Measures.
   real(R8P)                :: e_gh(2), order_gh, pg(3), hd !< Domain-face ghost error, its order, ghost centre.
   integer(I4P)             :: unfilled, m, i, j, k, b, c, nl, ijk(3), nb(3), dd, s, ierr !< Counters.
   integer(I4P)             :: nbnd, nkind(4) !< Boundary faces over the ranks, per kind.
   integer(I4P)             :: n3(3)     !< Cells per block.
   logical                  :: ok        !< Pass flag.
   logical                  :: active(3) !< Active directions.

   ok = .true.
   do m=1, 2
      block
      type(realm_object)      :: realm !< Realm (a fresh one per grid).
      type(MG_T) :: mg    !< Solver.
      call build(realm, ratio, nk, levels(m), sphere=.true., q=q, periodic=.false., corner=.true.)
      call set_grid_bc(realm, bc)
      active = .not.realm%adam%grid%null_xyz
      associate(f=>realm%adam%field, g=>realm%adam%grid)
      n3 = [g%ni, g%nj, nk]
      do b=1, f%blocks_number ! the source: the manufactured Laplacian over the factor
         do k=1, nk
            do j=1, g%nj
               do i=1, g%ni
                  q(1,i,j,k,b) = lap_m([f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)], active) / FACTOR
               enddo
            enddo
         enddo
      enddo
      allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
      phi = SENTINEL
      call mg%initialize(grid=g, field=f, max_iterations=MAX_ITERATIONS, bc=bc)
      allocate(ghost(max(mg%nbnd,1)))
      nkind = 0
      do i=1, mg%nbnd
         if (mg%bnd_kind(i) == FLAIL_BC_NEUMANN) then
            ghost(i) = sign(1_I4P, mg%bnd_dir(i)) * grad_m(mg%bnd_x(:,i), active, abs(mg%bnd_dir(i)))
         else
            ghost(i) = phi_m(mg%bnd_x(:,i), active)
         endif
         nkind(mg%bnd_kind(i)+1) = nkind(mg%bnd_kind(i)+1) + 1
      enddo
      call mg%solve(field=f, q=q, iv=1_I4P, factor=FACTOR, phi=phi, ghost=ghost)
      nl = f%blocks_number * g%ni * g%nj * nk
      allocate(x(max(nl,1)), y(max(nl,1)), pe(max(nl,1)), ye(max(nl,1)), yb(max(nl,1)))
      do i=1, nl
         x(i) = sin(1.3_R8P * i + 0.7_R8P * mpih%myrank) ; y(i) = cos(2.1_R8P * i - 0.4_R8P * mpih%myrank)
      enddo
      sym = mg%symmetry_defect(x, y)
      ! the manufactured potential at the centres, its volume mean (singular: Phi is compared up to a constant)
      vsum = 0._R8P ; pmean = 0._R8P ; umean = 0._R8P
      c = 0
      do b=1, f%blocks_number
         v = product(merge(g%cell_dxyz(:,f%coordinates(4,b)), 1._R8P, active))
         do k=1, nk
            do j=1, g%nj
               do i=1, g%ni
                  c = c + 1
                  pe(c) = phi_m([f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)], active)
                  vsum = vsum + v ; pmean = pmean + v * pe(c) ; umean = umean + v * phi(i,j,k,b)
               enddo
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, vsum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pmean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, umean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      pmean = pmean / vsum ; umean = umean / vsum
      if (.not.mg%singular) then
         pmean = 0._R8P ; umean = 0._R8P
      endif
      ! the truncation error |(A Phi_m + B g) / V - lap Phi_m|, the error of Phi, the unfilled face ghosts
      call mg%apply_a(pe, ye)
      call mg%apply_b(ghost, yb)
      trunc(m) = 0._R8P ; e_phi(m) = 0._R8P ; e_gh(m) = 0._R8P ; unfilled = 0
      c = 0
      do b=1, f%blocks_number
         v = product(merge(g%cell_dxyz(:,f%coordinates(4,b)), 1._R8P, active))
         do k=1, nk
            do j=1, g%nj
               do i=1, g%ni
                  c = c + 1
                  xyz = [f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)]
                  trunc(m) = max(trunc(m), abs((ye(c) + yb(c)) / v - lap_m(xyz, active)))
                  e_phi(m) = max(e_phi(m), abs((phi(i,j,k,b) - umean) - (pe(c) - pmean)))
                  do dd=1, 3
                     if (.not.active(dd)) cycle
                     do s=-1, 1, 2
                        ijk = [i, j, k] ; nb = ijk ; nb(dd) = nb(dd) + s
                        if (nb(dd) >= 1 .and. nb(dd) <= n3(dd)) cycle
                        if (phi(nb(1),nb(2),nb(3),b) == SENTINEL) unfilled = unfilled + 1
                        if (all(bc(:,dd) == FLAIL_BC_PERIODIC)) cycle ! a domain face of a bounded axis only
                        hd = g%cell_dxyz(dd,f%coordinates(4,b))
                        pg = xyz ; pg(dd) = pg(dd) + s * hd
                        if (pg(dd) > g%domain_emin(dd) .and. pg(dd) < g%domain_emax(dd)) cycle
                        e_gh(m) = max(e_gh(m), abs((phi(nb(1),nb(2),nb(3),b) - umean) - (phi_m(pg, active) - pmean)))
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, trunc(m), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, e_phi(m), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, e_gh(m), 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      nbnd = mg%nbnd
      call MPI_ALLREDUCE(MPI_IN_PLACE, nbnd, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, nkind, 4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      h(m) = g%cell_dxyz(1,levels(m))
      ok = ok .and. mg%residual <= mg%tolerance_used .and. sym <= TOL_SYMMETRY .and. unfilled == 0
      if (mpih%myrank == 0) print '(A)', '   '//label//' base level '//trim(str(levels(m), .true.))// &
         ', boundary faces '//trim(str(nbnd, .true.))//' (D '//trim(str(nkind(FLAIL_BC_DIRICHLET+1), .true.))// &
         ', N '//trim(str(nkind(FLAIL_BC_NEUMANN+1), .true.))//', G '// &
         trim(str(nkind(FLAIL_BC_DIRICHLET_GHOST+1), .true.))//'): iterations '//trim(str(mg%iterations, .true.))// &
         ', symmetry '//trim(str(sym, .true.))//', truncation '//trim(str(trunc(m), .true.))// &
         ', |Phi - Phi_m| '//trim(str(e_phi(m), .true.))//', domain-face ghosts '//trim(str(e_gh(m), .true.))// &
         ', unfilled '//trim(str(unfilled, .true.))
      endassociate
      call mg%destroy
      deallocate(phi, x, y, ghost, pe, ye, yb)
      endblock
   enddo
   order = log(e_phi(1) / e_phi(2)) / log(h(1) / h(2))
   order_gh = log(e_gh(1) / e_gh(2)) / log(h(1) / h(2))
   ok = ok .and. order >= MIN_ORDER .and. order_gh >= MIN_ORDER .and. trunc(2) / trunc(1) <= MAX_TRUNCATION_GROWTH
   if (mpih%myrank == 0) print '(A)', '   '//label//' order of Phi '//trim(str(order, .true.))// &
      ', of the domain-face ghosts '//trim(str(order_gh, .true.))//' (min '// &
      trim(str(MIN_ORDER, .true.))//'), truncation growth '//trim(str(trunc(2) / trunc(1), .true.))//' (max '// &
      trim(str(MAX_TRUNCATION_GROWTH, .true.))//')'//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   endsubroutine check_bc_amr

   subroutine set_grid_bc(realm, bc)
   !< The grid's periodicity from the kinds: an axis periodic in `bc` is periodic in the grid, the others bounded.
   type(realm_object), intent(inout) :: realm   !< Realm.
   integer(I4P),       intent(in)    :: bc(2,3) !< Kinds per face [lo/hi, axis].
   integer(I4P)                      :: t(6)    !< Grid boundary types.
   integer(I4P)                      :: d       !< Counter.

   t = 0_I4P
   do d=1, 3
      if (all(bc(:,d) == FLAIL_BC_PERIODIC) .and. .not.realm%adam%grid%null_xyz(d)) t(2*d-1:2*d) = BC_PERIODIC
   enddo
   call realm%adam%grid%set_bc_type(bc_type=t)
   endsubroutine set_grid_bc

   pure function phi_m(xyz, active) result(p)
   !< The manufactured potential, prod_d cos(2 pi x_d / L_d + theta_d) over the active directions.
   real(R8P), intent(in) :: xyz(3)    !< Point.
   logical,   intent(in) :: active(3) !< Active directions.
   real(R8P)             :: p         !< Potential.
   integer(I4P)          :: d         !< Counter.

   p = 1._R8P
   do d=1, 3
      if (active(d)) p = p * cos(2._R8P * PI * xyz(d) / L(d) + THETA(d))
   enddo
   endfunction phi_m

   pure function lap_m(xyz, active) result(lp)
   !< The Laplacian of the manufactured potential, -sum_d (2 pi / L_d)^2 phi_m.
   real(R8P), intent(in) :: xyz(3)    !< Point.
   logical,   intent(in) :: active(3) !< Active directions.
   real(R8P)             :: lp        !< Laplacian.

   lp = -sum(merge((2._R8P * PI / L)**2, 0._R8P, active)) * phi_m(xyz, active)
   endfunction lap_m

   pure function grad_m(xyz, active, d) result(gr)
   !< d phi_m / d x_d.
   real(R8P),    intent(in) :: xyz(3)    !< Point.
   logical,      intent(in) :: active(3) !< Active directions.
   integer(I4P), intent(in) :: d         !< Direction.
   real(R8P)                :: gr        !< Derivative.
   integer(I4P)             :: e         !< Counter.

   gr = -2._R8P * PI / L(d) * sin(2._R8P * PI * xyz(d) / L(d) + THETA(d))
   do e=1, 3
      if (active(e) .and. e /= d) gr = gr * cos(2._R8P * PI * xyz(e) / L(e) + THETA(e))
   enddo
   endfunction grad_m

   subroutine check_dirichlet(label, ratio, nk, levels, sphere)
   !< SV-12: Dirichlet boundaries, the exact potential as the ghost values; two base levels, uniform or with the box
   !< refined one level more.
   character(*), intent(in) :: label     !< Label.
   integer(I4P), intent(in) :: ratio     !< Tree ratio.
   integer(I4P), intent(in) :: nk        !< Cells per block along z.
   integer(I4P), intent(in) :: levels(2) !< Base levels.
   logical,      intent(in) :: sphere    !< Refine the central box.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi(:,:,:,:), x(:), y(:), ghost(:) !< Fields, vectors, ghost values.
   real(R8P)                :: e_phi(2), e_grav(2), h(2), trunc(2), sym, order !< Measures.
   integer(I4P)             :: unfilled, m, i, nl, iters(2), nbnd, ierr !< Counters.
   logical                  :: ok        !< Pass flag.
   logical                  :: active(3) !< Active directions.
   character(:), allocatable :: note     !< The uniform grid's extra assertion.

   ok = .true.
   do m=1, 2
      block
      type(realm_object)      :: realm !< Realm (a fresh one per grid).
      type(MG_T) :: mg    !< Solver.
      call build(realm, ratio, nk, levels(m), sphere=sphere, q=q, periodic=.false.)
      active = .not.realm%adam%grid%null_xyz
      q = q - RHO0 ! the exact potential solves lap Phi = factor (rho - rho0)
      allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
      phi = SENTINEL
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, max_iterations=MAX_ITERATIONS)
      allocate(ghost(max(mg%nbnd,1)))
      do i=1, mg%nbnd
         ghost(i) = phi_exact(mg%bnd_x(:,i), active)
      enddo
      call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi, ghost=ghost)
      nl = realm%adam%field%blocks_number * realm%adam%grid%ni * realm%adam%grid%nj * nk
      allocate(x(max(nl,1)), y(max(nl,1)))
      do i=1, nl
         x(i) = sin(1.3_R8P * i + 0.7_R8P * mpih%myrank) ; y(i) = cos(2.1_R8P * i - 0.4_R8P * mpih%myrank)
      enddo
      sym = mg%symmetry_defect(x, y)
      call truncation(realm, mg, active, nk, trunc(m), ghost=ghost, source_mean=0._R8P)
      call measure(realm, phi, active, e_phi(m), e_grav(m), unfilled, periodic=.false.)
      h(m) = realm%adam%grid%cell_dxyz(1,levels(m))
      iters(m) = mg%iterations
      nbnd = mg%nbnd
      call MPI_ALLREDUCE(MPI_IN_PLACE, nbnd, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      ok = ok .and. mg%residual <= mg%tolerance .and. sym <= TOL_SYMMETRY .and. unfilled == 0
      if (.not.sphere) ok = ok .and. mg%iterations <= 1
      if (mpih%myrank == 0) print '(A)', '   '//label//' base level '//trim(str(levels(m), .true.))// &
         ', boundary faces '//trim(str(nbnd, .true.))//': iterations '//trim(str(mg%iterations, .true.))// &
         ', symmetry '//trim(str(sym, .true.))//', truncation '//trim(str(trunc(m), .true.))// &
         ', |Phi - Phi_exact| '//trim(str(e_phi(m), .true.))//', unfilled '//trim(str(unfilled, .true.))
      call mg%destroy
      deallocate(phi, x, y, ghost)
      endblock
   enddo
   order = log(e_phi(1) / e_phi(2)) / log(h(1) / h(2))
   ok = ok .and. order >= MIN_ORDER .and. trunc(2) / trunc(1) <= MAX_TRUNCATION_GROWTH
   note = ''
   if (.not.sphere) note = ', iterations at most 1'
   if (mpih%myrank == 0) print '(A)', '   '//label//' order of Phi '//trim(str(order, .true.))//' (min '// &
      trim(str(MIN_ORDER, .true.))//'), truncation growth '//trim(str(trunc(2) / trunc(1), .true.))//' (max '// &
      trim(str(MAX_TRUNCATION_GROWTH, .true.))//')'//note//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   endsubroutine check_dirichlet

   subroutine check_uniform(label, ratio, nk, level)
   !< SV-9a: the multigrid on a uniform grid against the FFT solver.
   character(*), intent(in) :: label !< Label.
   integer(I4P), intent(in) :: ratio !< Tree ratio.
   integer(I4P), intent(in) :: nk    !< Cells per block along z.
   integer(I4P), intent(in) :: level !< Uniform level.
   type(realm_object)       :: realm !< Realm.
   type(MG_T)  :: mg    !< Multigrid solver.
   type(FFT_T) :: fft   !< FFT solver.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi_mg(:,:,:,:), phi_fft(:,:,:,:) !< Fields.
   real(R8P)                :: diff, scale_ !< Max difference, max |Phi_fft|.
   logical                  :: ok    !< Pass flag.
   integer(I4P)             :: n, ni, nj, ierr !< Sizes, MPI status.

   call build(realm, ratio, nk, level, sphere=.false., q=q)
   allocate(phi_mg(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
   allocate(phi_fft, mold=phi_mg)
   phi_mg = SENTINEL ; phi_fft = SENTINEL
   call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, max_iterations=MAX_ITERATIONS)
   call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_mg)
   call fft%initialize(grid=realm%adam%grid, field=realm%adam%field)
   call fft%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi_fft)
   n = realm%adam%field%blocks_number ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   diff = 0._R8P ; scale_ = 0._R8P
   if (n > 0) then
      diff = maxval(abs(phi_mg(1:ni,1:nj,1:nk,1:n) - phi_fft(1:ni,1:nj,1:nk,1:n)))
      scale_ = maxval(abs(phi_fft(1:ni,1:nj,1:nk,1:n)))
   endif
   call MPI_ALLREDUCE(MPI_IN_PLACE, diff, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, scale_, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ok = diff / scale_ <= TOL_FFT .and. mg%iterations <= 1
   if (mpih%myrank == 0) print '(A)', '   '//label//' |Phi_mg - Phi_fft| / max|Phi| '// &
      trim(str(diff / scale_, .true.))//', iterations '//trim(str(mg%iterations, .true.))//', residual '// &
      trim(str(mg%residual, .true.))//merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   call mg%destroy
   call fft%destroy
   endsubroutine check_uniform

   subroutine check_amr(label, ratio, nk, levels)
   !< SV-9b: static AMR, two base levels, each with the blocks near the sphere refined one level more.
   character(*), intent(in) :: label     !< Label.
   integer(I4P), intent(in) :: ratio     !< Tree ratio.
   integer(I4P), intent(in) :: nk        !< Cells per block along z.
   integer(I4P), intent(in) :: levels(2) !< Base levels.
   real(R8P), allocatable   :: q(:,:,:,:,:), phi(:,:,:,:), x(:), y(:) !< Fields, vectors.
   real(R8P)                :: e_phi(2), e_grav(2), h(2), trunc(2), sym, order, order_g !< Measures.
   integer(I4P)             :: unfilled, m, b, i, nl, ierr, nleaf, lmax !< Counters.
   logical                  :: ok        !< Pass flag.
   logical                  :: active(3) !< Active directions.

   ok = .true.
   do m=1, 2
      block
      type(realm_object)      :: realm !< Realm (a fresh one per grid).
      type(MG_T) :: mg    !< Solver.
      call build(realm, ratio, nk, levels(m), sphere=.true., q=q)
      active = .not.realm%adam%grid%null_xyz
      allocate(phi(lbound(q,2):ubound(q,2),lbound(q,3):ubound(q,3),lbound(q,4):ubound(q,4),ubound(q,5)))
      phi = SENTINEL
      call mg%initialize(grid=realm%adam%grid, field=realm%adam%field, max_iterations=MAX_ITERATIONS)
      call mg%solve(field=realm%adam%field, q=q, iv=1_I4P, factor=FACTOR, phi=phi)
      nl = realm%adam%field%blocks_number * realm%adam%grid%ni * realm%adam%grid%nj * nk
      allocate(x(max(nl,1)), y(max(nl,1)))
      do i=1, nl
         x(i) = sin(1.3_R8P * i + 0.7_R8P * mpih%myrank) ; y(i) = cos(2.1_R8P * i - 0.4_R8P * mpih%myrank)
      enddo
      sym = mg%symmetry_defect(x, y)
      call truncation(realm, mg, active, nk, trunc(m))
      call measure(realm, phi, active, e_phi(m), e_grav(m), unfilled)
      h(m) = realm%adam%grid%cell_dxyz(1,levels(m))
      nleaf = realm%adam%field%blocks_number
      lmax = 0
      do b=1, realm%adam%field%blocks_number
         lmax = max(lmax, realm%adam%field%coordinates(4,b))
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, nleaf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, lmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      ok = ok .and. mg%residual <= mg%tolerance .and. sym <= TOL_SYMMETRY .and. unfilled == 0
      if (mpih%myrank == 0) print '(A)', '   '//label//' base level '//trim(str(levels(m), .true.))//', leaves '// &
         trim(str(nleaf, .true.))//' (levels '//trim(str(levels(m), .true.))//'-'//trim(str(lmax, .true.))// &
         '): iterations '//trim(str(mg%iterations, .true.))//', factor '//trim(str(mg%factor_mean, .true.))// &
         ', symmetry '//trim(str(sym, .true.))//', truncation '//trim(str(trunc(m), .true.))// &
         ', |Phi - Phi_exact| '//trim(str(e_phi(m), .true.))// &
         ', 2:1 face gravity error '//trim(str(e_grav(m), .true.))//', unfilled '//trim(str(unfilled, .true.))
      call mg%destroy
      deallocate(phi, x, y)
      endblock
   enddo
   order = log(e_phi(1) / e_phi(2)) / log(h(1) / h(2))
   order_g = log(e_grav(1) / e_grav(2)) / log(h(1) / h(2))
   ok = ok .and. order >= MIN_ORDER .and. trunc(2) / trunc(1) <= MAX_TRUNCATION_GROWTH
   if (mpih%myrank == 0) print '(A)', '   '//label//' order of Phi '//trim(str(order, .true.))//' (min '// &
      trim(str(MIN_ORDER, .true.))//'), truncation growth '//trim(str(trunc(2) / trunc(1), .true.))//' (max '// &
      trim(str(MAX_TRUNCATION_GROWTH, .true.))//'), order of the 2:1 face gravity '//trim(str(order_g, .true.))// &
      merge('  PASS', '  FAIL', ok)
   test_passed = test_passed .and. ok
   endsubroutine check_amr

   subroutine truncation(realm, mg, active, nk, tmax, ghost, source_mean)
   !< The operator's truncation error: S applied to the exact potential against V 4 pi G (rho - rho0), the max over
   !< the leaves per unit volume. A consistent seam keeps it bounded under refinement; an inconsistent one (O(1/h))
   !< doubles it per level. Dirichlet: S Phi_exact + B g with the ghost values g, against V 4 pi G (rho - mean), the
   !< mean `source_mean` (the source the test passes is rho - rho0 already).
   type(realm_object),      intent(in)    :: realm     !< Realm.
   type(MG_T), intent(inout) :: mg        !< Solver.
   logical,                 intent(in)    :: active(3) !< Active directions.
   integer(I4P),            intent(in)    :: nk        !< Cells along z.
   real(R8P),               intent(out)   :: tmax      !< Max truncation.
   real(R8P), optional,     intent(in)    :: ghost(1:)   !< Dirichlet: the ghost values.
   real(R8P), optional,     intent(in)    :: source_mean !< The mean subtracted from rho (default rho0).
   real(R8P), allocatable :: pe(:), ye(:), yb(:)        !< Exact potential, the operator on it, B g.
   real(R8P)              :: rmean                      !< The mean subtracted from rho.
   real(R8P)              :: xyz(3), v, t              !< Measures.
   integer(I4P)           :: b, i, j, k, c, n, ierr    !< Counters.

   associate(f=>realm%adam%field, g=>realm%adam%grid)
   n = f%blocks_number * g%ni * g%nj * nk
   allocate(pe(max(n,1)), ye(max(n,1)))
   c = 0
   do b=1, f%blocks_number
      do k=1, nk
         do j=1, g%nj
            do i=1, g%ni
               c = c + 1
               pe(c) = phi_exact([f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)], active)
            enddo
         enddo
      enddo
   enddo
   call mg%apply_a(pe, ye)
   if (present(ghost)) then
      allocate(yb(max(n,1)))
      call mg%apply_b(ghost, yb)
      if (n > 0) ye(1:n) = ye(1:n) + yb(1:n)
   endif
   rmean = RHO0
   if (present(source_mean)) rmean = RHO0 - source_mean
   tmax = 0._R8P
   c = 0
   do b=1, f%blocks_number
      v = product(merge(g%cell_dxyz(:,f%coordinates(4,b)), 1._R8P, active))
      do k=1, nk
         do j=1, g%nj
            do i=1, g%ni
               c = c + 1
               xyz = [f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)]
               t = abs(ye(c) / v - FACTOR * (rho(xyz, active) - rmean))
               tmax = max(tmax, t)
            enddo
         enddo
      enddo
   enddo
   endassociate
   call MPI_ALLREDUCE(MPI_IN_PLACE, tmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   endsubroutine truncation

   subroutine measure(realm, phi, active, e_phi, e_grav, unfilled, periodic)
   !< Max error of Phi against the exact one (leaves); max error of the face gravity on the faces whose neighbour is
   !< coarser (the fine side of 2:1 faces, from the filled ghost); unfilled face ghosts. The leaves of every rank are
   !< gathered first (extents and levels) so that the neighbour level is found locally.
   type(realm_object), intent(in)  :: realm     !< Realm.
   real(R8P),          intent(in)  :: phi(1-realm%adam%grid%ngc:,1-realm%adam%grid%ngc:,1-realm%adam%grid%ngc:,1:) !< Phi.
   logical,            intent(in)  :: active(3) !< Active directions.
   real(R8P),          intent(out) :: e_phi     !< Max |Phi - Phi_exact|.
   real(R8P),          intent(out) :: e_grav    !< Max face-gravity error on 2:1 faces.
   integer(I4P),       intent(out) :: unfilled  !< Unfilled face ghosts.
   logical, optional,  intent(in)  :: periodic  !< Periodic boundaries (default); else a domain face is no 2:1 face.
   real(R8P),    allocatable :: ext(:,:), ext_all(:,:) !< Leaf extents and levels [7, n].
   integer(I4P), allocatable :: counts(:), displs(:)   !< Allgather counts and displacements.
   real(R8P)                 :: xyz(3), xf(3), p(3), g_num, g_ex, hd !< Measures.
   integer(I4P)              :: b, i, j, k, d, s, lv, nbl, ijk(3), nb(3), n(3), r, ierr !< Counters.

   associate(f=>realm%adam%field, g=>realm%adam%grid)
   nbl = f%blocks_number
   allocate(ext(7,max(nbl,1)), counts(mpih%procs_number), displs(mpih%procs_number))
   do b=1, nbl
      ext(1:3,b) = f%emin(:,b) ; ext(4:6,b) = f%emax(:,b) ; ext(7,b) = real(f%coordinates(4,b), R8P)
   enddo
   call MPI_ALLGATHER(7*nbl, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   displs(1) = 0
   do r=2, mpih%procs_number
      displs(r) = displs(r-1) + counts(r-1)
   enddo
   allocate(ext_all(7,sum(counts)/7))
   call MPI_ALLGATHERV(ext, 7*nbl, MPI_REAL8, ext_all, counts, displs, MPI_REAL8, MPI_COMM_WORLD, ierr)
   n = [g%ni, g%nj, g%nk]
   e_phi = 0._R8P ; e_grav = 0._R8P ; unfilled = 0
   do b=1, nbl
      lv = f%coordinates(4,b)
      do k=1, g%nk
         do j=1, g%nj
            do i=1, g%ni
               xyz = [f%x_cell(i,b), f%y_cell(j,b), f%z_cell(k,b)]
               e_phi = max(e_phi, abs(phi(i,j,k,b) - phi_exact(xyz, active)))
               do d=1, 3
                  if (.not.active(d)) cycle
                  do s=-1, 1, 2
                     ijk = [i, j, k] ; nb = ijk ; nb(d) = nb(d) + s
                     if (nb(d) >= 1 .and. nb(d) <= n(d)) cycle
                     if (phi(nb(1),nb(2),nb(3),b) == SENTINEL) unfilled = unfilled + 1
                     hd = g%cell_dxyz(d,lv)
                     p = xyz ; p(d) = p(d) + s * hd
                     if (present(periodic)) then
                        if (.not.periodic .and. (p(d) < g%domain_emin(d) .or. p(d) > g%domain_emax(d))) cycle
                     endif
                     p(d) = g%domain_emin(d) + modulo(p(d) - g%domain_emin(d), g%domain_emax(d) - g%domain_emin(d))
                     if (level_at(ext_all, p) >= lv) cycle
                     g_num = s * (phi(nb(1),nb(2),nb(3),b) - phi(i,j,k,b)) / hd
                     xf = xyz ; xf(d) = xf(d) + s * hd / 2._R8P
                     g_ex = grad_exact(xf, active, d)
                     e_grav = max(e_grav, abs(g_num - g_ex))
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   endassociate
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_phi, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, e_grav, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(MPI_IN_PLACE, unfilled, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   endsubroutine measure

   pure function level_at(ext_all, p) result(lv)
   !< The level of the leaf containing point p (a null direction always matches), from the gathered extents.
   real(R8P), intent(in) :: ext_all(1:,1:) !< Leaf extents and levels [7, n].
   real(R8P), intent(in) :: p(3)           !< Point.
   integer(I4P)          :: lv             !< Level (-1: none).
   integer(I4P)          :: b              !< Counter.

   lv = -1
   do b=1, size(ext_all, dim=2)
      if (all(p >= ext_all(1:3,b) .and. p <= ext_all(4:6,b))) then
         lv = nint(ext_all(7,b))
         return
      endif
   enddo
   endfunction level_at

   subroutine build(realm, ratio, nk, level, sphere, q, periodic, corner)
   !< The realm: uniform at `level`, and with `sphere` the central box [L/4, 3L/4] refined once more (with `corner`, the
   !< corner octant [0, L/2] instead, its leaves on the low domain faces).
   type(realm_object),     intent(inout) :: realm  !< Realm.
   integer(I4P),           intent(in)    :: ratio  !< Tree ratio.
   integer(I4P),           intent(in)    :: nk     !< Cells per block along z.
   integer(I4P),           intent(in)    :: level  !< Base level.
   logical,                intent(in)    :: sphere !< Refine near the sphere.
   real(R8P), allocatable, intent(inout) :: q(:,:,:,:,:) !< The source field.
   logical, optional,      intent(in)    :: periodic !< Periodic boundaries (default true).
   logical, optional,      intent(in)    :: corner   !< Refine the corner octant instead of the central box.
   character(:), allocatable             :: filename !< Input file.
   integer(I4P)                          :: ngc, ni, nj, b, i, j, k, d, ierr !< Counters.
   logical                               :: active(3) !< Active directions.
   type(tree_node_object), pointer       :: node      !< A tree node.
   logical                               :: corner_   !< Refine the corner octant.
   real(R8P)                             :: c(3)      !< A block centre.

   filename = 'test_poisson_mg-r'//trim(str(ratio, .true.))//'-l'//trim(str(level, .true.))//'.ini'
   if (mpih%myrank == 0) call write_input(filename=filename, nk=nk, ratio=ratio, level=level+1)
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
   call realm%initialize(filename=filename, memory_avail=8._R8P, nv=1_I4P) ! 3-D base level 3 + box: 960 blocks
   if (.not.present(periodic)) then
      call realm%adam%grid%set_bc_type(bc_type=[(BC_PERIODIC, d=1, 6)])
   elseif (periodic) then
      call realm%adam%grid%set_bc_type(bc_type=[(BC_PERIODIC, d=1, 6)])
   endif
   ngc = realm%adam%grid%ngc ; ni = realm%adam%grid%ni ; nj = realm%adam%grid%nj
   if (allocated(q)) deallocate(q)
   allocate(q(1:1,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc,1:realm%adam%field%nb))
   q = 0._R8P
   if (level > 0) call realm%adam%refine_uniform(refinement_levels=level, q=q, do_mpi_redistribute=.true., &
                                                 do_blocks_reorder=.false.)
   corner_ = .false.
   if (present(corner)) corner_ = corner
   if (sphere) then ! refine the central box [L/4, 3L/4] (or the corner [0, L/2]): block-aligned from base level 2 on
      do b=1, realm%adam%field%blocks_number
         c = 0.5_R8P * (realm%adam%field%emin(:,b) + realm%adam%field%emax(:,b))
         if ((.not.corner_ .and. all(abs(c - 0.5_R8P * L) < 0.25_R8P * L .or. realm%adam%grid%null_xyz)) .or. &
             (corner_ .and. all(c < 0.5_R8P * L .or. realm%adam%grid%null_xyz))) then
            node => realm%adam%tree%node(code=realm%adam%field%code(b))
            node%refinement_needed = TO_BE_REFINED
         endif
      enddo
      call realm%adam%amr_update(q=q, is_marked_by_tree=.true., do_mpi_redistribute=.true., do_blocks_reorder=.false.)
   endif
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
   endsubroutine build

   pure function wave(m, active) result(kv)
   !< Wave vector of mode m; none along an inactive direction.
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
   !< The source density.
   real(R8P), intent(in) :: xyz(3)       !< Point.
   logical,   intent(in) :: active(3)    !< Active directions.
   real(R8P)             :: r            !< Density.
   real(R8P)             :: k1(3), k2(3) !< Wave vectors.

   k1 = wave(1_I4P, active) ; k2 = wave(2_I4P, active)
   r = RHO0 + A1 * sin(k1(1) * xyz(1)) * cos(k1(2) * xyz(2)) * cos(k1(3) * xyz(3)) + A2 * cos(dot_product(k2, xyz))
   endfunction rho

   pure function phi_exact(xyz, active) result(p)
   !< The exact (continuum) potential, zero mean.
   real(R8P), intent(in) :: xyz(3)       !< Point.
   logical,   intent(in) :: active(3)    !< Active directions.
   real(R8P)             :: p            !< Potential.
   real(R8P)             :: k1(3), k2(3) !< Wave vectors.

   k1 = wave(1_I4P, active) ; k2 = wave(2_I4P, active)
   p = FACTOR * (A1 * sin(k1(1) * xyz(1)) * cos(k1(2) * xyz(2)) * cos(k1(3) * xyz(3)) / (-dot_product(k1, k1)) + &
                 A2 * cos(dot_product(k2, xyz)) / (-dot_product(k2, k2)))
   endfunction phi_exact

   pure function grad_exact(xyz, active, d) result(gr)
   !< d Phi / d x_d of the exact potential.
   real(R8P),    intent(in) :: xyz(3)              !< Point.
   logical,      intent(in) :: active(3)           !< Active directions.
   integer(I4P), intent(in) :: d                   !< Direction.
   real(R8P)                :: gr                  !< Derivative.
   real(R8P)                :: k1(3), k2(3), t1(3) !< Wave vectors, factors.

   k1 = wave(1_I4P, active) ; k2 = wave(2_I4P, active)
   t1 = [sin(k1(1) * xyz(1)), cos(k1(2) * xyz(2)), cos(k1(3) * xyz(3))]
   select case(d)
   case(1)
      gr = k1(1) * cos(k1(1) * xyz(1)) * t1(2) * t1(3)
   case(2)
      gr = -k1(2) * t1(1) * sin(k1(2) * xyz(2)) * t1(3)
   case default
      gr = -k1(3) * t1(1) * t1(2) * sin(k1(3) * xyz(3))
   endselect
   gr = FACTOR * (A1 * gr / (-dot_product(k1, k1)) - A2 * k2(d) * sin(dot_product(k2, xyz)) / (-dot_product(k2, k2)))
   endfunction grad_exact

   subroutine write_input(filename, nk, ratio, level, cube)
   !< Write the input of one configuration (the sections realm_object%initialize reads), box [0,1] x [0,2] x [0,3], or
   !< the cube [-1,1]^3 (isolated).
   character(*),      intent(in) :: filename !< Input file name.
   integer(I4P),      intent(in) :: nk       !< Cells per block along z.
   integer(I4P),      intent(in) :: ratio    !< Tree ratio.
   integer(I4P),      intent(in) :: level    !< Maximum level.
   logical, optional, intent(in) :: cube     !< The cube [-1,1]^3 (default false).
   integer(I4P)                  :: u        !< File unit.
   logical                       :: cube_    !< The cube.

   cube_ = .false.
   if (present(cube)) cube_ = cube
   open(newunit=u, file=filename, action='write', status='replace')
   write(u, '(A)') '[IO]'
   write(u, '(A)') 'output_basename        = test_poisson_mg'
   write(u, '(A)') 'it_save                = 999999'
   write(u, '(A)') 'restart                = .false.'
   write(u, '(A)') 'restart_basename       = test_poisson_mg-restart'
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
   if (cube_) then
      write(u, '(A)') 'emin_x = -1.0'
      write(u, '(A)') 'emin_y = -1.0'
      write(u, '(A)') 'emin_z = -1.0'
      write(u, '(A)') 'emax_x = 1.0'
      write(u, '(A)') 'emax_y = 1.0'
      write(u, '(A)') 'emax_z = 1.0'
   else
      write(u, '(A)') 'emin_x = 0.0'
      write(u, '(A)') 'emin_y = 0.0'
      write(u, '(A)') 'emin_z = 0.0'
      write(u, '(A)') 'emax_x = 1.0'
      write(u, '(A)') 'emax_y = 2.0'
      write(u, '(A)') 'emax_z = 3.0'
   endif
   write(u, '(A)') 'null_x = .false.'
   write(u, '(A)') 'null_y = .false.'
   if (nk == 1_I4P) then
      write(u, '(A)') 'null_z = .true.'
   else
      write(u, '(A)') 'null_z = .false.'
   endif
   write(u, '(A)') '[amr]'
   write(u, '(A)') 'max_level      = '//trim(str(level, .true.))
   write(u, '(A)') 'ratio          = '//trim(str(ratio, .true.))
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
endprogram test_poisson_mg
