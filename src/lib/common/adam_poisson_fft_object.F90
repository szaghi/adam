!< ADAM, FFT Poisson solver on uniform grids.
module adam_poisson_fft_object
!< ADAM, FFT Poisson solver on uniform grids.
!<
!< Solves `L_h phi = factor (f - mean(f))` for a cell field `f` (one variable of a conservative array `q`) on a grid whose
!< blocks all sit at one level and whose active directions are all periodic. `L_h` is the second-order 7-point Laplacian
!< (5-point in 2-D; a null direction contributes nothing). The solve divides each Fourier mode by the eigenvalue of
!< `L_h` itself, `lambda_k = sum_d (2 cos(2 pi k_d / N_d) - 2) / dx_d^2`, so the discrete residual of the returned `phi`
!< is at round-off; the modes with `lambda_k = 0` (the mean, and every mode along a null direction) are set to zero.
!<
!< **Data path: the distributed slab transform (issue #84, E5d).** Rank r owns the z-planes [zr(r), zr(r+1)) of the
!< transform array (`m` cells per direction: the grid, its extensions, its padding). Every block plane goes to the
!< owners of the target planes that map to it, and each rank builds its planes by the per-axis maps (`tmap`, `tsgn`:
!< a target cell's source cell and sign, the images of `extend_global` and the zero padding). The planes are transformed
!< in 2-D (r2c), transposed to z-lines (`MPI_ALLTOALLV`, the y-rows [yr(r), yr(r+1)) on rank r) and transformed in 1-D,
!< divided by the eigenvalues (or multiplied by the kernel's spectrum, kept in the pencil layout), and back. Every
!< block then fetches its potential with all its ghost layers from the owners of its planes (wrapped modulo `m`), so
!< no ghost exchange follows; the energy sums on the owners. One plan for one plane and one for one row of lines,
!< executed on every plane and row in turn through fixed buffers (`slab_planes`, `slab_lines`, the backend hooks):
!< the arithmetic of a value does not depend on how many planes or rows a rank holds, so the potential is bitwise the
!< same on every rank count (`poisson.sh` leg `slab`). At most `m(3)` ranks hold planes (a pencil layout would lift
!< it). PRISM's CG solves (`pcg`, below), the multigrid's redundant root and the CG's sub-problems keep the rank-0
!< path of the first version (issue #81, D-81-2): rank 0 gathers the source (`MPI_Gatherv`), places it, transforms
!< (`transform`) and scatters every block's `phi` (`MPI_Scatterv`); `global_potential` returns the transform array on
!< rank 0 in either path (diagnostics).
!<
!< **Backends.** This object transforms with FFTW 3 (`-D_FFTW`, `-lfftw3`; FFTW_ESTIMATE plans, deterministic). The FNL
!< twin `adam_fnl_poisson_fft_object` extends it and overrides the slab hooks (`slab_make_plans`, `slab_planes`,
!< `slab_lines`, `slab_destroy_plans`: each rank's GPU) and the rank-0 ones (`make_plans`, `transform`,
!< `destroy_plans`) with cuFFT. Built without either, `initialize` stops: "built without an FFT library".
!<
!< **Entry points.** `solve` packs one variable of a host conservative array; `solve_packed` takes the interior source
!< already packed (block-major, `i` fastest: what a device kernel produces on FNL). Both can return
!< `energy = sum_cells f Phi` (summed on the slab owners, every rank gets it): `W = energy dV / 2` is the
!< gravitational energy of self-gravity on a uniform grid.
!<
!< **Isolated boundaries (issue #81, M8-P4), Hockney and Eastwood.** `initialize(..., isolated=.true.)` solves
!< `L_h phi = factor f` with `phi -> 0` at infinity (no mean removed): the source sits in a zero-padded array of
!< `m >= 2 (n + ngc)` cells per direction, and the circular convolution with the lattice Green's function of `L_h`
!< (`adam_lattice_green_library`, `h^2 G(n)`) equals the free-space one on every cell the scatter reads, the ghost
!< layers included, since `m` exceeds the largest offset, `n - 1 + ngc`, twice. `m` is the smallest even size at or
!< above `2 (n + ngc)` with no prime factor above 13 (#84 E5b): the transform cost no longer depends on the factors of
!< `n + ngc` (medooza, n = 128: 262 = 2 x 131 cells cost FFTW 3.3 times as long as 264; the kernel entries beyond
!< the offsets the read cells need meet only the zero padding). The kernel's spectrum (`khat`, real: the kernel is
!< even) is computed once. Inside the domain the discrete residual is that of the kernel, about 1e-11 of the source
!< (`scripts/lattice_green.py`).
!<
!< **Homogeneous Dirichlet boundaries (issue #87, P1)** are for the multigrid's coarsest level only
!< (`initialize_blocks(..., dirichlet=.true.)`): `phi = 0` at the centres of the first ghost layer. The source is
!< extended odd about those centres to `m = 2 (n + 1)` cells per active direction, and the periodic transform of `m`
!< cells solves it: the odd extension is the type-I sine transform, exact for the 7-point operator with zero ghosts.
!<
!< **Boundary kinds per face (issue #84, E1a).** `bc(2,3)` (`adam_flail_parameters`) gives each face a homogeneous
!< kind: an axis is periodic, or each of its faces is Dirichlet (`phi = 0` at the face) or Neumann (zero normal
!< derivative at the face), possibly different at the two ends. The source is extended odd about a Dirichlet face and
!< even about a Neumann one, over 2n cells when the two faces agree and 4n when they differ (a quarter-wave
!< extension), and the periodic transform of the extension solves the 7-point operator with those boundary
!< conditions exactly; the scatter reads the ghost layers from the extension, so they carry the conditions. The modes
!< of zero eigenvalue are dropped as on a periodic grid: without a Dirichlet face the source mean is removed (the
!< operator is singular). The energy sums the interior only. Inhomogeneous values are the clients' (FLAIL's) work.
!<
!< **The wide operator D(G) (issue #84, E2a).** `dg_order = p` (2, 4, 6, 8; 0, the default, is the compact 7-point
!< Laplacian) solves PRISM's elliptic operator instead: the order-p centred first derivative `D1` composed with itself
!< per axis (the stencil reaches `p` cells, as PRISM's D and G read it), eigenvalue
!< `-sum_d s_p(theta_d)^2 / dx_d^2` with `s_p(theta) = sum_j 2 a_j sin(j theta)` (the `a_j` of `D1`). Every centred
!< `D1` annihilates the checkerboard (`s_p(pi) = 0`, set exactly: `sin(pi)` is 1e-16, not 0), so the modes with every
!< `theta_d` in {0, pi} are in the kernel and dropped, Dirichlet faces or not: G of them is zero, so the gradient a
!< client takes is unique. The odd and even extensions are reflections, which the wide stencil commutes with, so the
!< kinds per face hold as for the compact operator (the scatter fills `p` ghost layers from the extension); PRISM's
!< zero-ghost Dirichlet face is not a reflection (below). Measured on the prototype (`wide_dg_proto.py`): a compact
!< preconditioner needs O(N) iterations for D(G) (the checkerboard branch), so the direct solve here is the method on
!< uniform grids.
!<
!< **Zero-ghost faces (issue #84, E2b).** `FLAIL_BC_ZERO_GHOSTS` (both faces of an axis) is PRISM's elliptic Dirichlet
!< face: every ghost layer zero. With the compact operator it is the ghost-Dirichlet kind (the 7-point stencil reads one
!< layer). With the wide D(G) it is no reflection: the operator is the stencil truncated, symmetric and nonsingular, and
!< no transform diagonalises it. Rank 0 then solves it, on the gathered array, by conjugate gradients preconditioned by
!< the wide p = 2 operator with the same faces, inverted exactly: along a zero-ghost axis that operator splits into the
!< two parity sub-lattices, each the compact Laplacian at spacing 2h with zero at the ghost centres (the ghost-Dirichlet
!< kind), along the other axes it is the p = 2 operator with their own kinds; the 2^z sub-problems (z zero-ghost axes)
!< are FFT solves by sub-solvers of this solver's own type (`sub`). Spectrally equivalent (kappa <= 2.8 at p = 4,
!< 7.1 at p = 8): 1, 17, 24, 30 iterations at p = 2, 4, 6, 8, independent of the grid (`wide_dg_zero_ghost.py`,
!< `wide_dg_e2b.py`, with periodic and with reflected companion axes). The ghost layers returned carry the faces' rules
!< (zero, reflection, periodic wrap).
!<
!< **Free-space ghosts (issue #84, E3c).** `FLAIL_BC_FREE_GHOSTS` is PRISM's exact/open face: every ghost layer holds
!< the continuum free-space potential of the source, `g(x) = -sum_y factor f(y) dV / (4 pi |x - y|)`. The operator is
!< the zero-ghost one (the two kinds mix on an axis, and share the CG path and its preconditioner); the ghost values move
!< to the right-hand side, `D(G)_0 u = factor f - D(G) u_g` (u_g: g in the ghost slabs, zero inside), and the returned
!< ghost layers carry g. Rank 0 sums over the gathered source (`free_ghosts`), so the values do not depend on the
!< decomposition; the cost is O(ghost-slab cells x cells), a direct sum (a lattice free-space solve is #84 E7). The
!< wide D(G) only. Edges and corners of the ghost region stay zero: D(G) reads along the axis lines only.
!<
!< **Scope (fatal otherwise).** One level for every block, the blocks covering the domain. Without `bc`, every active
!< direction periodic. Isolated: 3-D, every direction non-periodic, cubic cells (the kernel is the unit lattice's), the
!< compact operator only.
!< AMR is the multigrid's (`adam_poisson_mg_object`), which calls this solver on its coarsest level.
use :: adam_grid_object, only : grid_object
use :: adam_lattice_green_library, only : lattice_green
use :: adam_field_object, only : field_object
use :: adam_flail_parameters, only : FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, FLAIL_BC_DIRICHLET_GHOST, &
                                      FLAIL_BC_ZERO_GHOSTS, FLAIL_BC_FREE_GHOSTS
use :: adam_mpih_global, only : mpih
use :: iso_c_binding
use :: mpi
use :: penf, only : I4P, I8P, R8P, str

implicit none
private
public :: poisson_fft_object

real(R8P), parameter :: PI=acos(-1._R8P) !< Pi.
! centred first-derivative coefficients a_j of order p: D1 q_i = sum_j a_j (q_{i+j} - q_{i-j}) / h (column p/2)
real(R8P), parameter :: D1_COEF(4,4)=reshape([1._R8P/2._R8P, 0._R8P, 0._R8P, 0._R8P,                         &
                                              2._R8P/3._R8P, -1._R8P/12._R8P, 0._R8P, 0._R8P,                 &
                                              3._R8P/4._R8P, -3._R8P/20._R8P, 1._R8P/60._R8P, 0._R8P,         &
                                              4._R8P/5._R8P, -1._R8P/5._R8P, 4._R8P/105._R8P, -1._R8P/280._R8P], &
                                             [4, 4]) !< a_j [j, p/2].

type :: fft_slot
   !< One sub-solver of the E2b preconditioner (issue #84). A scalar polymorphic pointer in a slot, not a polymorphic
   !< array component: nvfortran 26.9 rejects a recursive allocatable component (NVFORTRAN-S-0155) and crashed
   !< deallocating a polymorphic pointer array allocated with `mold` (pgf90_dealloc_mbr03a, on the FNL type). Owned:
   !< `destroy` frees it.
   class(poisson_fft_object), pointer :: s=>null() !< The sub-solver.
endtype fft_slot

type :: poisson_fft_object
   !< FFT Poisson solver on uniform periodic grids.
   logical                   :: is_initialized=.false. !< Initialization status.
   integer(I4P)              :: n(3)=1_I4P             !< Global cells per direction.
   integer(I4P)              :: m(3)=1_I4P             !< Transform cells per direction: n, or >= 2 (n + ngc) isolated.
   logical                   :: isolated=.false.       !< Isolated boundaries (zero padding), else periodic.
   logical                   :: dirichlet=.false.      !< Homogeneous Dirichlet boundaries (odd extension, #87).
   integer(I4P)              :: bc(2,3)=FLAIL_BC_PERIODIC !< Kinds per face [lo/hi, axis] (issue #84, E1a).
   logical                   :: extended=.false.       !< Some axis is extended (Dirichlet or Neumann faces).
   integer(I4P)              :: dg_order=0_I4P         !< 0: the compact Laplacian; p: the wide D(G) of order p (E2a).
   integer(I4P)              :: axis_order(3)=0_I4P    !< The operator per axis (dg_order, but in E2b's sub-problems).
   ! zero-ghost faces with the wide operator (issue #84, E2b): preconditioned CG on rank 0
   logical                   :: pcg=.false.            !< Some axis has zero-ghost faces and the operator is wide.
   logical                   :: free=.false.           !< Some face has free-space ghosts (E3c).
   real(R8P), allocatable    :: g_free(:,:,:)          !< Rank 0, E3c: the free-space ghost slabs, zero elsewhere.
   logical                   :: null_xyz(3)=.false.    !< Null directions.
   real(R8P)                 :: stencil(-8:8)=0._R8P   !< The 1-D D1 D1 stencil (times h^2), offsets -p..p.
   real(R8P)                 :: tolerance=1.e-12_R8P   !< Relative max-norm residual to stop.
   integer(I4P)              :: max_iterations=200_I4P !< CG iterations at most.
   integer(I4P)              :: iterations=0_I4P       !< Iterations of the last solve.
   real(R8P)                 :: residual=0._R8P        !< Relative residual of the last solve.
   logical                   :: redundant=.false.      !< Every rank holds the global arrays and solves (#84 E5c).
   real(R8P)                 :: wall(3)=0._R8P         !< Wall time (s) of `solve_blocks`, accumulated over calls on
                                                       !< this rank: gather, rank-0 solve, scatter (#84 E5; the
                                                       !< caller resets it).
   integer(I4P)              :: nsub=0_I4P             !< Preconditioner sub-problems, 2^z.
   integer(I4P),              allocatable :: sub_par(:,:) !< Parity per axis of each sub-problem [3,nsub] (-1: no split).
   type(fft_slot),            allocatable :: sub(:)       !< The preconditioner's sub-problems [nsub] (see fft_slot).
   integer(I4P)              :: nt(3)=1_I4P            !< Cells of the periodic transform the eigenvalues refer to: n,
                                                       !< or m with dirichlet.
   integer(I4P)              :: nijk(3)=1_I4P          !< Cells per block per direction.
   integer(I4P)              :: ngc=0_I4P              !< Ghost cells per side.
   integer(I4P)              :: level=0_I4P            !< Level of every block.
   real(R8P)                 :: dxyz(3)=0._R8P         !< Cell spacing per direction.
   real(R8P),    allocatable :: lam1(:)                !< Eigenvalue part along x [0:n(1)/2].
   real(R8P),    allocatable :: lam2(:)                !< Eigenvalue part along y [0:n(2)-1].
   real(R8P),    allocatable :: lam3(:)                !< Eigenvalue part along z [0:n(3)-1].
   real(R8P),    allocatable :: f_global(:,:,:)        !< Rank 0: the source, then the potential [m(1),m(2),m(3)].
   real(R8P),    allocatable :: f_source(:,:,:)        !< Rank 0: a copy of the source, for the energy.
   complex(R8P), allocatable :: f_hat(:,:,:)           !< Rank 0: the transform [m(1)/2+1,m(2),m(3)].
   real(R8P),    allocatable :: khat(:,:,:)            !< Rank 0, isolated: the kernel's spectrum over prod(m).
   type(c_ptr)               :: plan_r2c=c_null_ptr    !< FFTW forward plan (rank 0).
   type(c_ptr)               :: plan_c2r=c_null_ptr    !< FFTW backward plan (rank 0).
   ! the distributed slab transform (issue #84, E5d)
   logical                   :: slab=.false.           !< Every rank transforms its z-planes (see the module).
   integer(I4P)              :: h1=1_I4P               !< Complex cells along x after r2c, m(1)/2 + 1.
   integer(I4P), allocatable :: zr(:)                  !< First z-plane (0-based) of each rank's slab [0:ranks].
   integer(I4P), allocatable :: yr(:)                  !< First y-row (0-based) of each rank's pencils [0:ranks].
   integer(I4P), allocatable :: zown(:)                !< Owner rank of each z-plane [0:m(3)-1].
   integer(I4P), allocatable :: tmap(:,:)              !< Per axis, target cell -> source cell, -1 for zero [0:,3].
   real(R8P),    allocatable :: tsgn(:,:)              !< Per axis, the sign of the image [0:,3].
   real(R8P),    allocatable :: slab_r(:,:,:)          !< This rank's z-planes [m(1),m(2),planes].
   complex(R8P), allocatable :: slab_c(:,:,:)          !< Their 2-D transforms [h1,m(2),planes].
   complex(R8P), allocatable :: pen(:,:,:)             !< This rank's z-lines [m(3),h1,rows].
   real(R8P),    allocatable :: khat_pen(:,:,:)        !< Isolated: the kernel's spectrum in pencils [m(3),h1,rows].
   real(R8P),    allocatable :: w_r(:,:)               !< FFTW: the plane buffer [m(1),m(2)].
   complex(R8P), allocatable :: w_c(:,:)               !< FFTW: the transformed plane buffer [h1,m(2)].
   complex(R8P), allocatable :: w_l(:,:)               !< FFTW: the row-of-lines buffer [m(3),h1].
   type(c_ptr)               :: plan_pr2c=c_null_ptr   !< FFTW: one plane, r2c.
   type(c_ptr)               :: plan_pc2r=c_null_ptr   !< FFTW: one plane, c2r.
   type(c_ptr)               :: plan_lfwd=c_null_ptr   !< FFTW: one row of z-lines, forward.
   type(c_ptr)               :: plan_lbwd=c_null_ptr   !< FFTW: one row of z-lines, backward.
   contains
      ! public methods
      procedure, pass(self) :: destroy    !< Destroy the solver.
      procedure, pass(self) :: initialize !< Check the grid, size the solver, make the plans.
      procedure, pass(self) :: solve        !< Solve for phi from one variable of q.
      procedure, pass(self) :: solve_packed !< Solve for phi from the packed interior source.
      procedure, pass(self) :: initialize_blocks !< Size the solver from block counts and spacing (no grid, no field).
      procedure, pass(self) :: solve_blocks      !< Solve from an explicit list of blocks (coordinates and source).
      ! private methods
      procedure, pass(self), private :: make_kernel !< Rank 0, isolated: the kernel's spectrum.
      procedure, pass(self), private :: extend_global  !< Rank 0: the source extended over the transform cells.
      procedure, pass(self), private :: apply_wide     !< Rank 0: y = D(G) x with the faces' ghost rules (E2b).
      procedure, pass(self), private :: fill_ghosts    !< Rank 0: the ghost layers of a padded array by the faces.
      procedure, pass(self), private :: precondition   !< Rank 0: z = M^-1 r, the sub-lattice solves (E2b).
      procedure, pass(self), private :: solve_pcg      !< Rank 0: f_global <- factor D(G)^-1 f_global by PCG (E2b).
      procedure, pass(self), private :: free_ghosts    !< Rank 0: g_free, the free-space ghost slabs (E3c).
      procedure, pass(self), private :: subtract_free_ghosts !< Rank 0: f <- f - D(G) u_g (E3c).
      procedure, pass(self), private :: put_free_ghosts      !< Rank 0: the free-space ghost slabs into a padded array.
      procedure, pass(self), private :: solve_everywhere     !< `solve_blocks` with `redundant`: every rank solves.
      procedure, pass(self), private :: slab_setup           !< The slab layout, maps, arrays and plans (E5d).
      procedure, pass(self), private :: solve_slabs          !< `solve_blocks` with `slab`: the distributed solve.
      procedure, pass(self), private :: slab_transform       !< slab_r <- factor * L^-1 slab_r, distributed.
      procedure, pass(self)          :: global_potential     !< Gather the slab solution to rank 0 (diagnostics).
      ! backend hooks
      procedure, pass(self) :: destroy_plans    !< Destroy the transform plans.
      procedure, pass(self) :: make_plans       !< Make the transform plans (rank 0).
      procedure, pass(self) :: transform        !< f_global <- factor * L_h^-1 f_global (rank 0).
      procedure, pass(self) :: transform_kernel !< khat <- the spectrum of f_global over prod(m) (rank 0, isolated).
      procedure, nopass     :: root_cells_max   !< Most cells of a multigrid root solved by this backend (#84 E5c).
      procedure, pass(self) :: slab_make_plans    !< Make the slab plans: one plane, one row of lines (E5d).
      procedure, pass(self) :: slab_destroy_plans !< Destroy the slab plans.
      procedure, pass(self) :: slab_planes        !< 2-D transform of every local plane, slab_r <-> slab_c.
      procedure, pass(self) :: slab_lines         !< 1-D transform of every local z-line, in place in pen.
      procedure, pass(self) :: slab_transpose     !< slab_c <-> pen (MPI_ALLTOALLV).
      procedure, pass(self) :: slab_spectral      !< pen <- the spectral solve (eigenvalues or kernel) of pen.
      procedure, pass(self) :: slab_keep_kernel   !< khat_pen <- pen over prod(m) (isolated, at setup).
      procedure, pass(self) :: slab_update        !< slab_r between the host and where the stages run.
endtype poisson_fft_object

#ifdef _FFTW
integer(c_int), parameter :: FFTW_ESTIMATE=64_c_int !< FFTW planner flag, `fftw3.f03`.
integer(c_int), parameter :: FFTW_FORWARD=-1_c_int  !< FFTW sign, `fftw3.f03`.
integer(c_int), parameter :: FFTW_BACKWARD=1_c_int  !< FFTW sign, `fftw3.f03`.

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

   type(c_ptr) function fftw_plan_dft_r2c_2d(n0, n1, in, out, flags) bind(C, name='fftw_plan_dft_r2c_2d')
   !< FFTW: plan a 2-D real-to-complex transform (row-major sizes).
   import :: c_ptr, c_int, c_double, c_double_complex
   integer(c_int), value     :: n0, n1 !< Sizes, slowest first.
   real(c_double)            :: in(*)  !< Input.
   complex(c_double_complex) :: out(*) !< Output.
   integer(c_int), value     :: flags  !< Planner flags.
   endfunction fftw_plan_dft_r2c_2d

   type(c_ptr) function fftw_plan_dft_c2r_2d(n0, n1, in, out, flags) bind(C, name='fftw_plan_dft_c2r_2d')
   !< FFTW: plan a 2-D complex-to-real transform (row-major sizes).
   import :: c_ptr, c_int, c_double, c_double_complex
   integer(c_int), value     :: n0, n1 !< Sizes, slowest first.
   complex(c_double_complex) :: in(*)  !< Input.
   real(c_double)            :: out(*) !< Output.
   integer(c_int), value     :: flags  !< Planner flags.
   endfunction fftw_plan_dft_c2r_2d

   type(c_ptr) function fftw_plan_many_dft(rnk, n, howmany, in, inembed, istride, idist, out, onembed, ostride, &
                                           odist, sgn, flags) bind(C, name='fftw_plan_many_dft')
   !< FFTW: plan a batch of complex transforms.
   import :: c_ptr, c_int, c_double_complex
   integer(c_int), value     :: rnk             !< Rank of each transform.
   integer(c_int)            :: n(*)            !< Sizes.
   integer(c_int), value     :: howmany         !< Transforms.
   complex(c_double_complex) :: in(*)           !< Input.
   type(c_ptr), value        :: inembed         !< NULL: n.
   integer(c_int), value     :: istride, idist  !< Input stride and distance.
   complex(c_double_complex) :: out(*)          !< Output.
   type(c_ptr), value        :: onembed         !< NULL: n.
   integer(c_int), value     :: ostride, odist  !< Output stride and distance.
   integer(c_int), value     :: sgn             !< FFTW_FORWARD or FFTW_BACKWARD.
   integer(c_int), value     :: flags           !< Planner flags.
   endfunction fftw_plan_many_dft

   subroutine fftw_execute_dft(p, in, out) bind(C, name='fftw_execute_dft')
   !< FFTW: execute a complex plan on the given arrays.
   import :: c_ptr, c_double_complex
   type(c_ptr), value        :: p      !< Plan.
   complex(c_double_complex) :: in(*)  !< Input.
   complex(c_double_complex) :: out(*) !< Output.
   endsubroutine fftw_execute_dft

   subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
   !< FFTW: destroy a plan.
   import :: c_ptr
   type(c_ptr), value :: p !< Plan.
   endsubroutine fftw_destroy_plan
endinterface
#endif

contains
   ! public methods
   recursive subroutine destroy(self)
   !< Destroy the solver (and its sub-solvers).
   class(poisson_fft_object), intent(inout) :: self !< The solver.
   integer(I4P)                             :: i    !< Counter.

   if (self%is_initialized) call self%destroy_plans
   if (self%is_initialized .and. self%slab) call self%slab_destroy_plans
   if (allocated(self%zr))       deallocate(self%zr)
   if (allocated(self%yr))       deallocate(self%yr)
   if (allocated(self%zown))     deallocate(self%zown)
   if (allocated(self%tmap))     deallocate(self%tmap)
   if (allocated(self%tsgn))     deallocate(self%tsgn)
   if (allocated(self%slab_r))   deallocate(self%slab_r)
   if (allocated(self%slab_c))   deallocate(self%slab_c)
   if (allocated(self%pen))      deallocate(self%pen)
   if (allocated(self%khat_pen)) deallocate(self%khat_pen)
   if (allocated(self%w_r))      deallocate(self%w_r)
   if (allocated(self%w_c))      deallocate(self%w_c)
   if (allocated(self%w_l))      deallocate(self%w_l)
   self%slab = .false. ; self%h1 = 1_I4P
   if (allocated(self%lam1))     deallocate(self%lam1)
   if (allocated(self%lam2))     deallocate(self%lam2)
   if (allocated(self%lam3))     deallocate(self%lam3)
   if (allocated(self%f_global)) deallocate(self%f_global)
   if (allocated(self%f_source)) deallocate(self%f_source)
   if (allocated(self%f_hat))    deallocate(self%f_hat)
   if (allocated(self%khat))     deallocate(self%khat)
   if (allocated(self%sub)) then
      do i=1, size(self%sub)
         if (associated(self%sub(i)%s)) then
            call self%sub(i)%s%destroy
            deallocate(self%sub(i)%s)
         endif
      enddo
      deallocate(self%sub)
   endif
   if (allocated(self%sub_par))  deallocate(self%sub_par)
   if (allocated(self%g_free))   deallocate(self%g_free)
   self%free = .false.
   self%pcg = .false. ; self%nsub = 0_I4P ; self%iterations = 0_I4P ; self%residual = 0._R8P
   self%dg_order = 0_I4P ; self%axis_order = 0_I4P ; self%stencil = 0._R8P ; self%null_xyz = .false.
   self%n = 1_I4P ; self%m = 1_I4P ; self%nijk = 1_I4P ; self%ngc = 0_I4P ; self%level = 0_I4P ; self%dxyz = 0._R8P
   self%isolated = .false.
   self%redundant = .false.
   self%dirichlet = .false.
   self%bc = FLAIL_BC_PERIODIC
   self%extended = .false.
   self%nt = 1_I4P
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine initialize(self, grid, field, isolated, bc, dg_order)
   !< Check the grid (one level, the domain covered, the boundaries), size the solver, make the plans and, isolated,
   !< the kernel. `bc` (issue #84, E1a): the kind of every face; a periodic axis of the grid must be periodic in `bc`
   !< and the other way round. Without it every active axis must be periodic (or the solve isolated).
   class(poisson_fft_object), intent(inout) :: self         !< The solver.
   type(grid_object),         intent(in)    :: grid         !< The grid.
   type(field_object),        intent(in)    :: field        !< The field (block coordinates).
   logical, optional,         intent(in)    :: isolated     !< Isolated boundaries (default: periodic).
   integer(I4P), optional,    intent(in)    :: bc(2,3)      !< Kinds per face [lo/hi, axis] (default periodic).
   integer(I4P), optional,    intent(in)    :: dg_order     !< 0 compact (default); 2, 4, 6, 8: the wide D(G) (E2a).
   integer(I4P)                             :: lmin, lmax   !< Min and max block level over the ranks.
   integer(I4P)                             :: nblk(3)      !< Blocks per direction at the level.
   integer(I4P)                             :: blocks_total !< Blocks over the ranks.
   integer(I4P)                             :: d            !< Counter.
   integer(I4P)                             :: ierr         !< MPI error status.

   call self%destroy
   self%isolated = .false. ; if (present(isolated)) self%isolated = isolated
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
   if (self%isolated) then
      do d=1, 3
         if (grid%null_xyz(d)) call mpih%error_stop(msg='poisson_fft_object%initialize: isolated boundaries are '// &
            '3-D only, direction '//trim(str(d, .true.))//' is null')
         if (grid%is_ijk_periodic(d)) call mpih%error_stop(msg='poisson_fft_object%initialize: isolated '// &
            'boundaries need non-periodic directions, direction '//trim(str(d, .true.))//' is periodic')
      enddo
      if (abs(self%dxyz(2) - self%dxyz(1)) > 1e-12_R8P * self%dxyz(1) .or. &
          abs(self%dxyz(3) - self%dxyz(1)) > 1e-12_R8P * self%dxyz(1)) call mpih%error_stop(msg= &
         'poisson_fft_object%initialize: isolated boundaries need cubic cells (the lattice Green''s function is the '// &
         'unit lattice''s), cells '//trim(str(self%dxyz(1)))//' x '//trim(str(self%dxyz(2)))//' x '// &
         trim(str(self%dxyz(3))))
   elseif (present(bc)) then
      do d=1, 3
         if (grid%null_xyz(d)) cycle
         if (grid%is_ijk_periodic(d) .neqv. all(bc(:,d) == FLAIL_BC_PERIODIC)) call mpih%error_stop(msg= &
            'poisson_fft_object%initialize: direction '//trim(str(d, .true.))//' is '// &
            merge('    periodic', 'not periodic', grid%is_ijk_periodic(d))//' in the grid and '// &
            merge('    periodic', 'not periodic', all(bc(:,d) == FLAIL_BC_PERIODIC))//' in bc')
      enddo
   else
      do d=1, 3
         if (.not.grid%null_xyz(d) .and. .not.grid%is_ijk_periodic(d)) call mpih%error_stop(msg= &
            'poisson_fft_object%initialize: direction '//trim(str(d, .true.))//' is not periodic; the FFT solver '// &
            'supports periodic boundaries (or isolated ones, 3-D)')
      enddo
   endif
   if (present(bc) .and. .not.self%isolated) then
      call self%initialize_blocks(nblk=nblk, nijk=self%nijk, ngc=self%ngc, dxyz=self%dxyz, level=self%level, &
                                  null_xyz=grid%null_xyz, isolated=self%isolated, bc=bc, dg_order=dg_order)
   else
      call self%initialize_blocks(nblk=nblk, nijk=self%nijk, ngc=self%ngc, dxyz=self%dxyz, level=self%level, &
                                  null_xyz=grid%null_xyz, isolated=self%isolated, dg_order=dg_order)
   endif
   endsubroutine initialize

   recursive subroutine initialize_blocks(self, nblk, nijk, ngc, dxyz, level, null_xyz, isolated, dirichlet, bc, &
                                          dg_order, axis_dg_order, redundant)
   !< Size the solver for a uniform grid of `nblk` blocks per direction, `nijk` cells each, at `level` with spacing
   !< `dxyz`, and make the plans (and, isolated, the kernel) on rank 0. No check of a grid: the caller owns the scope
   !< (the multigrid's coarsest level calls it directly, issue #81 P5; with `dirichlet`, issue #87).
   !<
   !< `bc` (issue #84, E1a) gives the boundary kind of every face (`adam_flail_parameters`), homogeneous: an axis is
   !< periodic (both faces), or its faces are Dirichlet/Neumann at the face (the source extended odd/even about each
   !< face: 2n cells when the two faces agree, 4n when they differ), or ghost-Dirichlet (both faces, odd about the ghost
   !< centres, 2(n+1) cells: what `dirichlet` sets on every active axis). Periodic against non-periodic, and
   !< ghost-Dirichlet against a face kind, on one axis are fatal; a null axis takes no extension.
   !<
   !< `redundant` (issue #84, E5c): every rank holds the global arrays and the plans, and `solve_blocks` gathers the
   !< source to every rank (`MPI_ALLGATHERV`), transforms on every rank and keeps each rank's own blocks: no scatter,
   !< and no rank waits on another's transform. For a small grid (the multigrid's root); the compact operator only, not
   !< isolated, no energy.
   class(poisson_fft_object), intent(inout)        :: self        !< The solver.
   integer(I4P),              intent(in)           :: nblk(3)     !< Blocks per direction.
   integer(I4P),              intent(in)           :: nijk(3)     !< Cells per block per direction.
   integer(I4P),              intent(in)           :: ngc         !< Ghost cells per side (isolated: the padding margin).
   real(R8P),                 intent(in)           :: dxyz(3)     !< Cell spacing.
   integer(I4P),              intent(in)           :: level       !< Level of every block.
   logical,                   intent(in)           :: null_xyz(3) !< Null directions.
   logical,                   intent(in)           :: isolated    !< Isolated boundaries, else periodic.
   logical,                   intent(in), optional :: dirichlet   !< Homogeneous ghost-Dirichlet everywhere (default false).
   integer(I4P),              intent(in), optional :: bc(2,3)     !< Kinds per face [lo/hi, axis] (default periodic).
   integer(I4P),              intent(in), optional :: dg_order    !< 0 compact (default); 2, 4, 6, 8: the wide D(G).
   integer(I4P),              intent(in), optional :: axis_dg_order(3) !< The operator per axis (E2b's sub-problems).
   logical,                   intent(in), optional :: redundant   !< Every rank solves (default false).
   integer(I4P)                                    :: nblk_(3), nijk_(3), ngc_, level_ !< Copies of the arguments.
   integer(I4P)                                    :: dg_order_   !< Copy of the argument.
   integer(I4P)                                    :: axis_order_(3) !< Copy of the argument.
   real(R8P)                                       :: dxyz_(3)    !< Copy of the argument.
   logical                                         :: null_(3), isolated_, dirichlet_ !< Copies of the arguments.
   integer(I4P)                                    :: bc_(2,3)    !< Copy of the argument.
   integer(I4P)                                    :: k, d        !< Counters.

   ! copy first: the arguments may be components of self (initialize passes them), which destroy resets
   nblk_ = nblk ; nijk_ = nijk ; ngc_ = ngc ; level_ = level ; dxyz_ = dxyz ; null_ = null_xyz ; isolated_ = isolated
   dirichlet_ = .false.
   if (present(dirichlet)) dirichlet_ = dirichlet
   bc_ = FLAIL_BC_PERIODIC
   if (present(bc)) bc_ = bc
   if (dirichlet_) bc_ = FLAIL_BC_DIRICHLET_GHOST
   if (isolated_ .and. (dirichlet_ .or. any(bc_ /= FLAIL_BC_PERIODIC))) call mpih%error_stop( &
      msg='poisson_fft_object%initialize_blocks: isolated excludes dirichlet and bc')
   dg_order_ = 0_I4P
   if (present(dg_order)) dg_order_ = dg_order
   if (all(dg_order_ /= [0_I4P, 2_I4P, 4_I4P, 6_I4P, 8_I4P])) call mpih%error_stop(msg= &
      'poisson_fft_object%initialize_blocks: dg_order '//trim(str(dg_order_, .true.))//' is not 0 (compact), 2, 4, 6 or 8')
   if (isolated_ .and. dg_order_ /= 0_I4P) call mpih%error_stop(msg='poisson_fft_object%initialize_blocks: isolated '// &
      'boundaries take the compact operator only (the lattice Green''s function is its), not the wide D(G)')
   axis_order_ = dg_order_
   if (present(axis_dg_order)) axis_order_ = axis_dg_order
   do d=1, 3
      if (null_(d)) cycle
      if (count(bc_(:,d) == FLAIL_BC_PERIODIC) == 1) call mpih%error_stop(msg='poisson_fft_object%initialize_blocks: '//&
         'axis '//trim(str(d, .true.))//' has one periodic face and one non-periodic face')
      if (count(bc_(:,d) == FLAIL_BC_DIRICHLET_GHOST) == 1) call mpih%error_stop(msg= &
         'poisson_fft_object%initialize_blocks: axis '//trim(str(d, .true.))//' mixes a ghost-Dirichlet face with a '// &
         'face kind')
      if (count(bc_(:,d) == FLAIL_BC_ZERO_GHOSTS .or. bc_(:,d) == FLAIL_BC_FREE_GHOSTS) == 1) call mpih%error_stop(msg=&
         'poisson_fft_object%initialize_blocks: axis '//trim(str(d, .true.))//' mixes a zero-ghost face with another '// &
         'kind (both faces or neither: a reflection at the other end couples the parity sub-lattices, issue #84 E2b)')
      if (any(bc_(:,d) == FLAIL_BC_FREE_GHOSTS) .and. axis_order_(d) == 0_I4P) call mpih%error_stop(msg= &
         'poisson_fft_object%initialize_blocks: free-space ghosts on axis '//trim(str(d, .true.))//' need the wide '// &
         'D(G) (dg_order 2, 4, 6 or 8; issue #84 E3c)')
      do k=1, 2
         if (all(bc_(k,d) /= [FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, FLAIL_BC_DIRICHLET_GHOST,  &
                              FLAIL_BC_ZERO_GHOSTS, FLAIL_BC_FREE_GHOSTS])) &
            call mpih%error_stop(msg='poisson_fft_object%initialize_blocks: unknown boundary kind '// &
                                     trim(str(bc_(k,d), .true.))//' on axis '//trim(str(d, .true.)))
      enddo
   enddo
   ! zero ghosts: with the compact operator, the ghost-Dirichlet kind (one ghost layer read); wide, preconditioned CG
   if (all(axis_order_ == 0_I4P)) where (bc_ == FLAIL_BC_ZERO_GHOSTS) bc_ = FLAIL_BC_DIRICHLET_GHOST
   call self%destroy
   self%nijk = nijk_
   self%ngc = ngc_
   self%dxyz = dxyz_
   self%level = level_
   self%isolated = isolated_
   self%dirichlet = dirichlet_
   self%bc = bc_
   self%dg_order = dg_order_
   self%axis_order = axis_order_
   self%null_xyz = null_
   self%pcg = .false.
   do d=1, 3 ! a null axis has no ghosts the operator reads
      if (.not.null_(d)) self%pcg = self%pcg .or. any(self%bc(:,d) == FLAIL_BC_ZERO_GHOSTS .or. &
                                                    self%bc(:,d) == FLAIL_BC_FREE_GHOSTS)
   enddo
   self%free = any(self%bc == FLAIL_BC_FREE_GHOSTS)
   if (self%pcg) self%stencil = dg_stencil(self%dg_order)
   self%n = nblk_ * self%nijk
   self%m = self%n
   if (self%pcg) then
      do d=1, 3
         if (null_(d)) cycle
         if (self%n(d) < self%dg_order) call mpih%error_stop(msg='poisson_fft_object%initialize_blocks: '// &
            trim(str(self%n(d), .true.))//' cells along axis '//trim(str(d, .true.))//' are fewer than the stencil''s '// &
            'reach, '//trim(str(self%dg_order, .true.)))
      enddo
   elseif (self%isolated) then
      do d=1, 3 ! at least 2 (n + ngc); the next 13-smooth size (#84 E5b: FFTW slows on large prime factors)
         self%m(d) = smooth_size(2 * (self%n(d) + self%ngc))
      enddo
   else
      do d=1, 3
         if (null_(d)) cycle
         if (all(self%bc(:,d) == FLAIL_BC_DIRICHLET_GHOST)) then
            self%m(d) = 2 * (self%n(d) + 1) ! odd about the ghost centres 0 and n + 1
         elseif (all(self%bc(:,d) /= FLAIL_BC_PERIODIC)) then
            self%m(d) = merge(2, 4, self%bc(1,d) == self%bc(2,d)) * self%n(d) ! odd/even about each face
         endif
      enddo
   endif
   self%extended = (.not.self%isolated) .and. any(self%m /= self%n) ! never with pcg (m = n)
   self%nt = merge(self%n, self%m, self%isolated)
   ! eigenvalue parts of L_h per direction on the nt-cell periodic transform; zero along a null direction
   allocate(self%lam1(0:self%nt(1)/2), self%lam2(0:self%nt(2)-1), self%lam3(0:self%nt(3)-1))
   self%lam1 = 0._R8P ; self%lam2 = 0._R8P ; self%lam3 = 0._R8P
   if (.not.null_(1)) then
      do k=0, self%nt(1)/2
         self%lam1(k) = axis_eigenvalue(k, self%nt(1), self%dxyz(1), self%axis_order(1))
      enddo
   endif
   if (.not.null_(2)) then
      do k=0, self%nt(2)-1
         self%lam2(k) = axis_eigenvalue(k, self%nt(2), self%dxyz(2), self%axis_order(2))
      enddo
   endif
   if (.not.null_(3)) then
      do k=0, self%nt(3)-1
         self%lam3(k) = axis_eigenvalue(k, self%nt(3), self%dxyz(3), self%axis_order(3))
      enddo
   endif
   if (present(redundant)) self%redundant = redundant
   if (self%redundant .and. (self%isolated .or. self%pcg)) call mpih%error_stop(msg= &
      'poisson_fft_object%initialize_blocks: redundant is for the compact operator, not isolated')
   ! the distributed slab transform (#84 E5d), for every solve but PRISM's CG ones (rank 0, #84 decision 3), the
   ! redundant root and the CG's own sub-problems
   self%slab = .not.(self%pcg .or. self%redundant .or. present(axis_dg_order))
   if (self%slab) then
      call self%slab_setup
   elseif (mpih%myrank == 0 .or. self%redundant) then
      allocate(self%f_global(self%m(1),self%m(2),self%m(3)))
      allocate(self%f_source(self%m(1),self%m(2),self%m(3)))
      allocate(self%f_hat(self%m(1)/2+1,self%m(2),self%m(3)))
      self%f_global = 0._R8P
      self%f_source = 0._R8P
      self%f_hat = (0._R8P, 0._R8P)
      if (.not.self%pcg) call self%make_plans ! the CG solve transforms its sub-problems only
      if (self%isolated) call self%make_kernel
   endif
   if (self%pcg) call make_preconditioner
   self%is_initialized = .true.
   contains
      subroutine make_preconditioner
      !< The 2^z sub-problems of the wide p = 2 operator (z zero-ghost axes): along a zero-ghost axis one parity class
      !< (ceil or floor of n/2 cells), the compact operator at spacing 2h with zero at the ghost centres; along the other
      !< axes the p = 2 operator with their kinds. Sub-solvers of this solver's own type (cuFFT on FNL), one block each.
      integer(I4P) :: z, c, d_, bit       !< Zero-ghost axes, counters.
      integer(I4P) :: n_sub(3), ord_sub(3), bc_sub(2,3) !< A sub-problem's cells, operators, kinds.
      real(R8P)    :: dx_sub(3)           !< A sub-problem's spacing.

      z = count([(zero_class(d_), d_=1, 3)])
      self%nsub = 2_I4P**z
      allocate(self%sub_par(3,self%nsub))
      allocate(self%sub(self%nsub))
      do c=1, self%nsub
         bit = 0
         do d_=1, 3
            self%sub_par(d_,c) = -1
            n_sub(d_) = self%n(d_) ; dx_sub(d_) = self%dxyz(d_) ; ord_sub(d_) = 2_I4P ; bc_sub(:,d_) = self%bc(:,d_)
            if (zero_class(d_)) then
               self%sub_par(d_,c) = ibits(c-1, bit, 1)
               bit = bit + 1
               n_sub(d_) = (self%n(d_) - self%sub_par(d_,c) + 1) / 2 ! cells i (1-based) with mod(i-1, 2) = parity
               dx_sub(d_) = 2._R8P * self%dxyz(d_) ; ord_sub(d_) = 0_I4P ; bc_sub(:,d_) = FLAIL_BC_DIRICHLET_GHOST
            endif
            if (null_(d_)) ord_sub(d_) = 0_I4P
         enddo
         allocate(self%sub(c)%s, mold=self)
         call self%sub(c)%s%initialize_blocks(nblk=[1_I4P, 1_I4P, 1_I4P], nijk=n_sub, ngc=0_I4P, dxyz=dx_sub,     &
                                            level=self%level, null_xyz=null_, isolated=.false., bc=bc_sub, &
                                            dg_order=2_I4P, axis_dg_order=ord_sub)
      enddo
      endsubroutine make_preconditioner

      pure function zero_class(d_)
      !< The axis has zero-ghost (or free-space ghost) faces and is not null: its parity classes decouple.
      integer(I4P), intent(in) :: d_         !< Axis.
      logical                  :: zero_class !< Result.

      zero_class = .not.self%null_xyz(d_) .and. any(self%bc(:,d_) == FLAIL_BC_ZERO_GHOSTS .or. &
                                                    self%bc(:,d_) == FLAIL_BC_FREE_GHOSTS)
      endfunction zero_class
   endsubroutine initialize_blocks

   subroutine solve(self, field, q, iv, factor, phi, energy)
   !< Solve `L_h phi = factor (q(iv) - mean(q(iv)))` (periodic) or `L_h phi = factor q(iv)` with phi -> 0 at infinity
   !< (isolated); every ghost layer of phi is filled from the global solution.
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
   !< Solve `L_h phi = factor (f - mean(f))` (periodic) or `L_h phi = factor f` (isolated) from the packed interior source
   !< of this rank's blocks (block-major, `i` fastest, `field%blocks_number * ni * nj * nk` values); every ghost layer of
   !< phi is filled from the global solution (the periodic wrap, or the padded free-space solution).
   class(poisson_fft_object), intent(inout)         :: self        !< The solver.
   type(field_object),        intent(in)            :: field       !< The field (block coordinates).
   real(R8P),                 intent(in)            :: buf(1:)     !< Packed interior source of this rank.
   real(R8P),                 intent(in)            :: factor      !< Source factor (4 pi G for gravity).
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Potential [i,j,k,nb].
   real(R8P),                 intent(out), optional :: energy      !< sum_cells f phi (all ranks).
   integer(I4P)                                     :: nbl         !< Blocks of this rank.
   integer(I4P), allocatable                        :: coords(:,:) !< Block coordinates of this rank [3,nbl].
   real(R8P),    allocatable                        :: rbuf(:)     !< The potential of this rank's blocks, ghosts included.
   integer(I4P)                                     :: ub(3)       !< Upper bounds of phi per direction.
   integer(I4P)                                     :: b, i, j, k, c !< Counters.

   if (.not.self%is_initialized) call mpih%error_stop(msg='poisson_fft_object%solve: called before initialize')
   ub = self%nijk + self%ngc
   if (ubound(phi, dim=1) /= ub(1) .or. ubound(phi, dim=2) /= ub(2) .or. ubound(phi, dim=3) /= ub(3) .or. &
       size(phi, dim=4) < field%blocks_number) &
      call mpih%error_stop(msg='poisson_fft_object%solve: phi is not shaped [1-ngc:n+ngc per direction, blocks]')
   nbl = field%blocks_number
   allocate(coords(3,max(nbl,1)))
   do b=1, nbl
      if (field%coordinates(4,b) /= self%level) call mpih%error_stop(msg='poisson_fft_object%solve: a block left '// &
         'level '//trim(str(self%level, .true.))//' after initialize (regrids are not supported)')
      coords(:,b) = field%coordinates(1:3,b)
   enddo
   allocate(rbuf(max(nbl*product(self%nijk + 2 * self%ngc),1)))
   call self%solve_blocks(nbl=nbl, coords=coords, buf=buf, factor=factor, gw=self%ngc, out=rbuf, energy=energy)
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
   endsubroutine solve_packed

   subroutine solve_blocks(self, nbl, coords, buf, factor, gw, out, energy)
   !< The solve on an explicit list of this rank's blocks: their integer coordinates at the solver's level and their
   !< packed interior source (block-major, `i` fastest); `out` returns each block's potential with `gw` ghost layers
   !< (block-major, `i` fastest, from `1-gw` to `n+gw`), filled from the global solution. Optionally the energy
   !< `sum_cells f phi` (rank 0 computes it, every rank gets it).
   class(poisson_fft_object), intent(inout)         :: self            !< The solver.
   integer(I4P),              intent(in)            :: nbl             !< Blocks of this rank.
   integer(I4P),              intent(in)            :: coords(1:,1:)   !< Block coordinates [3,nbl].
   real(R8P),                 intent(in)            :: buf(1:)         !< Packed interior source of this rank.
   real(R8P),                 intent(in)            :: factor          !< Source factor (4 pi G for gravity).
   integer(I4P),              intent(in)            :: gw              !< Ghost layers of `out`.
   real(R8P),                 intent(inout)         :: out(1:)         !< Packed potential, ghosts included.
   real(R8P),                 intent(out), optional :: energy          !< sum_cells f phi (all ranks).
   integer(I4P)                                     :: ncell           !< Interior cells per block.
   integer(I4P)                                     :: nghost          !< Cells per block, ghosts included.
   integer(I4P)                                     :: nr              !< Ranks.
   integer(I4P)                                     :: blocks_total    !< Rank 0: blocks over the ranks.
   integer(I4P), allocatable                        :: counts(:)       !< Rank 0: blocks per rank.
   integer(I4P), allocatable                        :: displs(:)       !< Rank 0: displacement per rank (blocks).
   integer(I4P), allocatable                        :: rcounts(:)      !< Rank 0: counts per rank (values).
   integer(I4P), allocatable                        :: rdispls(:)      !< Rank 0: displacements per rank (values).
   integer(I4P), allocatable                        :: coords_all(:,:) !< Rank 0: block coordinates of every rank.
   real(R8P),    allocatable                        :: buf_all(:)      !< Rank 0: buffer of every rank.
   real(R8P),    allocatable                        :: u_pad(:,:,:)    !< Rank 0, E2b: the solution with its ghosts.
   real(R8P)                                        :: e               !< sum_cells f phi.
   integer(I4P)                                     :: ub(3)           !< Upper bounds per direction.
   integer(I4P)                                     :: g(3)            !< Global cell index, wrapped.
   integer(I4P)                                     :: b, i, j, k, c, r !< Counters.
   integer(I4P)                                     :: ierr            !< MPI error status.
   real(R8P)                                        :: t(0:3)          !< Phase marks (MPI_WTIME).

   if (.not.self%is_initialized) call mpih%error_stop(msg='poisson_fft_object%solve: called before initialize')
   if (self%redundant) then
      if (present(energy)) call mpih%error_stop(msg='poisson_fft_object%solve: no energy with redundant')
      call self%solve_everywhere(nbl=nbl, coords=coords, buf=buf, factor=factor, gw=gw, out=out)
      return
   endif
   if (self%slab) then
      call self%solve_slabs(nbl=nbl, coords=coords, buf=buf, factor=factor, gw=gw, out=out, energy=energy)
      return
   endif
   t(0) = MPI_WTIME()
   ub     = self%nijk + gw
   ncell  = product(self%nijk)
   nghost = product(self%nijk + 2 * gw)
   nr     = mpih%procs_number
   if (size(buf) < nbl*ncell) call mpih%error_stop(msg='poisson_fft_object%solve: the packed source is too short')
   if (size(out) < nbl*nghost) call mpih%error_stop(msg='poisson_fft_object%solve: the packed output is too short')
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
   allocate(coords_all(3,max(blocks_total,1)), buf_all(max(blocks_total*nghost,blocks_total*ncell,1)))
   call MPI_GATHERV(coords, 3*nbl, MPI_INTEGER, coords_all, rcounts, rdispls, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   if (mpih%myrank == 0) then
      rcounts = ncell * counts ; rdispls = ncell * displs
   endif
   call MPI_GATHERV(buf, nbl*ncell, MPI_REAL8, buf_all, rcounts, rdispls, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   t(1) = MPI_WTIME()
   ! rank 0: place, transform, pack every block with its ghosts by the periodic wrap
   e = 0._R8P
   if (mpih%myrank == 0) then
      if (self%isolated .or. self%extended) self%f_global = 0._R8P ! the padding (the images) must hold no source
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
      if (self%extended) call self%extend_global
      if (present(energy)) self%f_source = self%f_global
      if (self%free) call self%free_ghosts(factor=factor, gd=max(gw, self%dg_order))
      if (self%pcg) then
         call self%solve_pcg(factor=factor)
      else
         call self%transform(factor=factor)
      endif
      if (present(energy)) then
         if (self%extended) then ! the interior only: the images are not cells
            e = sum(self%f_source(1:self%n(1),1:self%n(2),1:self%n(3)) * self%f_global(1:self%n(1),1:self%n(2),1:self%n(3)))
         else
            e = sum(self%f_source * self%f_global)
         endif
      endif
      if (self%pcg) then ! the ghosts by the faces' rules (no extension holds them)
         allocate(u_pad(1-gw:self%n(1)+gw,1-gw:self%n(2)+gw,1-gw:self%n(3)+gw))
         u_pad = 0._R8P
         u_pad(1:self%n(1),1:self%n(2),1:self%n(3)) = self%f_global
         call self%fill_ghosts(u_pad, gw)
         if (self%free) call self%put_free_ghosts(u_pad, gw)
      endif
      c = 0
      do b=1, blocks_total
         do k=1-gw, ub(3)
            g(3) = modulo(coords_all(3,b)*self%nijk(3) + k - 1, self%m(3)) + 1
            do j=1-gw, ub(2)
               g(2) = modulo(coords_all(2,b)*self%nijk(2) + j - 1, self%m(2)) + 1
               do i=1-gw, ub(1)
                  g(1) = modulo(coords_all(1,b)*self%nijk(1) + i - 1, self%m(1)) + 1
                  c = c + 1
                  if (self%pcg) then
                     buf_all(c) = u_pad(coords_all(1,b)*self%nijk(1)+i, coords_all(2,b)*self%nijk(2)+j, &
                                        coords_all(3,b)*self%nijk(3)+k)
                  else
                     buf_all(c) = self%f_global(g(1),g(2),g(3))
                  endif
               enddo
            enddo
         enddo
      enddo
      rcounts = nghost * counts ; rdispls = nghost * displs
   endif
   t(2) = MPI_WTIME()
   ! scatter back, ghosts included
   call MPI_SCATTERV(buf_all, rcounts, rdispls, MPI_REAL8, out, nbl*nghost, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   if (present(energy)) then
      call MPI_BCAST(e, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      energy = e
   endif
   t(3) = MPI_WTIME()
   self%wall = self%wall + (t(1:3) - t(0:2))
   endsubroutine solve_blocks

   subroutine solve_everywhere(self, nbl, coords, buf, factor, gw, out)
   !< `solve_blocks` with `redundant` (issue #84, E5c): the source of every rank gathered to every rank, placed and
   !< transformed on every rank (the same arrays, the same plans: the same values as rank 0 computes in
   !< `solve_blocks`), and each rank's blocks packed from its own global solution. `wall`: gather, solve, 0.
   class(poisson_fft_object), intent(inout) :: self          !< The solver.
   integer(I4P),              intent(in)    :: nbl           !< Blocks of this rank.
   integer(I4P),              intent(in)    :: coords(1:,1:) !< Block coordinates [3,nbl].
   real(R8P),                 intent(in)    :: buf(1:)       !< Packed interior source of this rank.
   real(R8P),                 intent(in)    :: factor        !< Source factor.
   integer(I4P),              intent(in)    :: gw            !< Ghost layers of `out`.
   real(R8P),                 intent(inout) :: out(1:)       !< Packed potential, ghosts included.
   integer(I4P)                             :: ncell, nr     !< Interior cells per block, ranks.
   integer(I4P), allocatable                :: counts(:), displs(:), rcounts(:), rdispls(:) !< Allgather layout.
   integer(I4P), allocatable                :: coords_all(:,:) !< Block coordinates of every rank.
   real(R8P),    allocatable                :: buf_all(:)    !< Source of every rank.
   integer(I4P)                             :: ub(3), g(3)   !< Upper bounds, wrapped global index.
   integer(I4P)                             :: b, i, j, k, c, r, blocks_total !< Counters.
   integer(I4P)                             :: ierr          !< MPI error status.
   real(R8P)                                :: t(0:2)        !< Phase marks (MPI_WTIME).

   t(0) = MPI_WTIME()
   ub    = self%nijk + gw
   ncell = product(self%nijk)
   nr    = mpih%procs_number
   allocate(counts(nr), displs(nr), rcounts(nr), rdispls(nr))
   call MPI_ALLGATHER(nbl, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   displs(1) = 0
   do r=2, nr
      displs(r) = displs(r-1) + counts(r-1)
   enddo
   blocks_total = sum(counts)
   allocate(coords_all(3,max(blocks_total,1)), buf_all(max(blocks_total*ncell,1)))
   rcounts = 3 * counts ; rdispls = 3 * displs
   call MPI_ALLGATHERV(coords, 3*nbl, MPI_INTEGER, coords_all, rcounts, rdispls, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   rcounts = ncell * counts ; rdispls = ncell * displs
   call MPI_ALLGATHERV(buf, nbl*ncell, MPI_REAL8, buf_all, rcounts, rdispls, MPI_REAL8, MPI_COMM_WORLD, ierr)
   t(1) = MPI_WTIME()
   if (self%extended) self%f_global = 0._R8P
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
   if (self%extended) call self%extend_global
   call self%transform(factor=factor)
   c = 0
   do b=displs(mpih%myrank+1)+1, displs(mpih%myrank+1)+nbl
      do k=1-gw, ub(3)
         g(3) = modulo(coords_all(3,b)*self%nijk(3) + k - 1, self%m(3)) + 1
         do j=1-gw, ub(2)
            g(2) = modulo(coords_all(2,b)*self%nijk(2) + j - 1, self%m(2)) + 1
            do i=1-gw, ub(1)
               g(1) = modulo(coords_all(1,b)*self%nijk(1) + i - 1, self%m(1)) + 1
               c = c + 1
               out(c) = self%f_global(g(1),g(2),g(3))
            enddo
         enddo
      enddo
   enddo
   t(2) = MPI_WTIME()
   self%wall(1:2) = self%wall(1:2) + (t(1:2) - t(0:1))
   endsubroutine solve_everywhere

   subroutine slab_setup(self)
   !< The distributed slab transform (issue #84, E5d): rank r owns the z-planes [zr(r), zr(r+1)) of the transform array
   !< and, after the transpose, the y-rows [yr(r), yr(r+1)) of its z-lines; the per-axis maps from a target cell of the
   !< transform array to its source cell and sign (the extensions of `extend_global`, the isolated padding); the
   !< arrays, the plans and, isolated, the kernel's spectrum in the pencil layout.
   class(poisson_fft_object), intent(inout) :: self !< The solver.
   integer(I4P)                             :: nr, me, r, d, t, mm, nzl, nyl !< Counters, sizes.

   nr = mpih%procs_number ; me = mpih%myrank
   self%h1 = self%m(1) / 2 + 1
   allocate(self%zr(0:nr), self%yr(0:nr), self%zown(0:self%m(3)-1))
   do r=0, nr
      self%zr(r) = int((int(r, I8P) * self%m(3)) / nr, I4P)
      self%yr(r) = int((int(r, I8P) * self%m(2)) / nr, I4P)
   enddo
   do r=0, nr-1
      if (self%zr(r+1) > self%zr(r)) self%zown(self%zr(r):self%zr(r+1)-1) = r
   enddo
   mm = maxval(self%m)
   allocate(self%tmap(0:mm-1,3), self%tsgn(0:mm-1,3))
   self%tmap = -1 ; self%tsgn = 0._R8P
   do d=1, 3
      do t=0, self%m(d)-1
         call axis_map(d, t, self%tmap(t,d), self%tsgn(t,d))
      enddo
   enddo
   nzl = self%zr(me+1) - self%zr(me) ; nyl = self%yr(me+1) - self%yr(me)
   allocate(self%slab_r(self%m(1),self%m(2),max(nzl,1)), self%slab_c(self%h1,self%m(2),max(nzl,1)), &
            self%pen(self%m(3),self%h1,max(nyl,1)))
   self%slab_r = 0._R8P ; self%slab_c = (0._R8P, 0._R8P) ; self%pen = (0._R8P, 0._R8P)
   call self%slab_make_plans
   if (self%isolated) call slab_kernel
   contains
      recursive subroutine axis_map(d_, t_, src, sg)
      !< Target cell t_ (0-based) of axis d_: its source cell (0-based, -1 for zero) and sign, as `extend_global`.
      integer(I4P), intent(in)  :: d_, t_ !< Axis, target cell.
      integer(I4P), intent(out) :: src    !< Source cell.
      real(R8P),    intent(out) :: sg     !< Sign.
      integer(I4P)              :: n_     !< Source cells.

      n_ = self%n(d_)
      src = t_ ; sg = 1._R8P
      if (t_ < n_) return ! the domain itself
      src = -1 ; sg = 0._R8P
      if (self%isolated) return ! the padding
      if (all(self%bc(:,d_) == FLAIL_BC_DIRICHLET_GHOST)) then ! odd about the ghost centres n and m - 1
         if (t_ == n_ .or. t_ == self%m(d_) - 1) return
         src = 2 * n_ - t_ ; sg = -1._R8P
      elseif (t_ < 2 * n_) then ! the reflection about the high face
         src = 2 * n_ - 1 - t_ ; sg = merge(-1._R8P, 1._R8P, self%bc(2,d_) == FLAIL_BC_DIRICHLET)
      else ! the reflection of [0, 2n) about the low face
         call axis_map(d_, 4 * n_ - 1 - t_, src, sg)
         sg = sg * merge(-1._R8P, 1._R8P, self%bc(1,d_) == FLAIL_BC_DIRICHLET)
      endif
      endsubroutine axis_map

      subroutine slab_kernel
      !< Isolated: the kernel `h^2 G(o)` on this rank's planes (offsets wrapped into (-m/2, m/2]), transformed as the
      !< source is, its spectrum over prod(m) kept in the pencil layout (real: the kernel is even).
      integer(I4P) :: o(3), i, j, k !< Offset, counters.

      do k=1, nzl
         o(3) = self%zr(me) + k - 1 ; if (o(3) > self%m(3) / 2) o(3) = o(3) - self%m(3)
         do j=1, self%m(2)
            o(2) = j - 1 ; if (o(2) > self%m(2) / 2) o(2) = o(2) - self%m(2)
            do i=1, self%m(1)
               o(1) = i - 1 ; if (o(1) > self%m(1) / 2) o(1) = o(1) - self%m(1)
               self%slab_r(i,j,k) = self%dxyz(1)**2 * lattice_green(o)
            enddo
         enddo
      enddo
      call self%slab_update(to_device=.true.)
      call self%slab_planes(forward=.true.)
      call self%slab_transpose(forward=.true.)
      call self%slab_lines(forward=.true.)
      allocate(self%khat_pen(self%m(3),self%h1,max(nyl,1)))
      call self%slab_keep_kernel
      self%slab_r = 0._R8P
      endsubroutine slab_kernel
   endsubroutine slab_setup

   subroutine solve_slabs(self, nbl, coords, buf, factor, gw, out, energy)
   !< `solve_blocks` with `slab` (issue #84, E5d): every block plane sent to the owners of the target planes that map to
   !< it, each rank's planes built (the extensions and the padding by the maps), transformed (`slab_transform`), and
   !< every block's potential with its ghost layers fetched back from the owners of its planes (wrapped modulo m, as
   !< `solve_blocks` reads `f_global`). The arithmetic of every value does not depend on the rank count (one plan per
   !< plane and per row of lines). `wall`: redistribution, transform, return.
   class(poisson_fft_object), intent(inout)         :: self          !< The solver.
   integer(I4P),              intent(in)            :: nbl           !< Blocks of this rank.
   integer(I4P),              intent(in)            :: coords(1:,1:) !< Block coordinates [3,nbl].
   real(R8P),                 intent(in)            :: buf(1:)       !< Packed interior source of this rank.
   real(R8P),                 intent(in)            :: factor        !< Source factor.
   integer(I4P),              intent(in)            :: gw            !< Ghost layers of `out`.
   real(R8P),                 intent(inout)         :: out(1:)       !< Packed potential, ghosts included.
   real(R8P),                 intent(out), optional :: energy        !< sum_cells f phi (all ranks).
   integer(I4P), allocatable :: tgo(:,:), ntgo(:)          !< Owner ranks of the targets of each source z-plane.
   integer(I4P), allocatable :: srcidx(:)                  !< Source z-plane -> local index in src (0: not needed).
   real(R8P),    allocatable :: src(:,:,:)                 !< This rank's source planes [n(1),n(2),needed].
   integer(I4P), allocatable :: sc(:), sd(:), rc(:), rd(:) !< Alltoallv counts and displacements (items).
   integer(I4P), allocatable :: shdr(:), rhdr(:), fill(:)  !< Headers, cursors.
   real(R8P),    allocatable :: sval(:), rval(:)           !< Values.
   integer(I4P)              :: nr, me, nb(3), np_, w, nzl, nsrc, b, i, j, k, t, t1, t2, t3, s1, s2, s3, r, c, e, q
   real(R8P)                 :: sg3, sg23, el, tm(0:3)     !< Signs, energy, phase marks.
   integer(I4P)              :: ierr                       !< MPI error status.

   tm(0) = MPI_WTIME()
   nr = mpih%procs_number ; me = mpih%myrank
   nb = self%nijk ; np_ = nb(1) * nb(2)
   nzl = self%zr(me+1) - self%zr(me)
   allocate(sc(0:nr-1), sd(0:nr-1), rc(0:nr-1), rd(0:nr-1), fill(0:nr-1))
   ! the owner ranks of the targets of every source plane (1, 2 or 4 targets, deduplicated)
   allocate(tgo(4,0:self%n(3)-1), ntgo(0:self%n(3)-1))
   ntgo = 0
   do t=0, self%m(3)-1
      s3 = self%tmap(t,3)
      if (s3 < 0) cycle
      r = self%zown(t)
      if (any(tgo(1:ntgo(s3),s3) == r)) cycle
      ntgo(s3) = ntgo(s3) + 1 ; tgo(ntgo(s3),s3) = r
   enddo
   ! this rank's source planes (those its target planes map to)
   allocate(srcidx(0:self%n(3)-1))
   srcidx = 0 ; nsrc = 0
   do t=self%zr(me), self%zr(me+1)-1
      s3 = self%tmap(t,3)
      if (s3 < 0) cycle
      if (srcidx(s3) > 0) cycle
      nsrc = nsrc + 1 ; srcidx(s3) = nsrc
   enddo
   allocate(src(self%n(1),self%n(2),max(nsrc,1)))
   ! forward: every block plane (header c1, c2, s3 and its np_ values) to each other owner that needs it, this rank's
   ! own straight into src
   sc = 0
   do b=1, nbl
      do k=1, nb(3)
         s3 = coords(3,b) * nb(3) + k - 1
         do q=1, ntgo(s3)
            if (tgo(q,s3) /= me) sc(tgo(q,s3)) = sc(tgo(q,s3)) + 1
         enddo
      enddo
   enddo
   call exchange_counts
   allocate(shdr(max(3*sum(sc),1)), sval(max(np_*sum(sc),1)), rhdr(max(3*sum(rc),1)), rval(max(np_*sum(rc),1)))
   fill = sd
   do b=1, nbl
      do k=1, nb(3)
         s3 = coords(3,b) * nb(3) + k - 1
         e = (b-1)*np_*nb(3) + (k-1)*np_
         do q=1, ntgo(s3)
            r = tgo(q,s3)
            if (r == me) then
               call put_plane(coords(1,b), coords(2,b), srcidx(s3), buf(e+1:e+np_))
            else
               c = fill(r) ; fill(r) = fill(r) + 1
               shdr(3*c+1:3*c+3) = [coords(1,b), coords(2,b), s3]
               sval(np_*c+1:np_*c+np_) = buf(e+1:e+np_)
            endif
         enddo
      enddo
   enddo
   if (nr > 1) then
      call MPI_ALLTOALLV(shdr, 3*sc, 3*sd, MPI_INTEGER, rhdr, 3*rc, 3*rd, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call MPI_ALLTOALLV(sval, np_*sc, np_*sd, MPI_REAL8, rval, np_*rc, np_*rd, MPI_REAL8, MPI_COMM_WORLD, ierr)
   endif
   do c=0, sum(rc)-1
      call put_plane(rhdr(3*c+1), rhdr(3*c+2), srcidx(rhdr(3*c+3)), rval(np_*c+1:np_*c+np_))
   enddo
   ! this rank's target planes by the maps
   do k=1, nzl
      t3 = self%zr(me) + k - 1
      s3 = self%tmap(t3,3) ; sg3 = self%tsgn(t3,3)
      if (s3 < 0) then
         self%slab_r(:,:,k) = 0._R8P
         cycle
      endif
      do t2=0, self%m(2)-1
         s2 = self%tmap(t2,2) ; sg23 = sg3 * self%tsgn(t2,2)
         do t1=0, self%m(1)-1
            s1 = self%tmap(t1,1)
            if (s1 < 0 .or. s2 < 0) then
               self%slab_r(t1+1,t2+1,k) = 0._R8P
            else
               self%slab_r(t1+1,t2+1,k) = sg23 * self%tsgn(t1,1) * src(s1+1,s2+1,srcidx(s3))
            endif
         enddo
      enddo
   enddo
   tm(1) = MPI_WTIME()
   call self%slab_transform(factor=factor)
   tm(2) = MPI_WTIME()
   if (present(energy)) then ! the domain's cells only: the images and the padding are not cells
      el = 0._R8P
      do k=1, nzl
         t3 = self%zr(me) + k - 1
         if (t3 >= self%n(3)) cycle
         do j=1, self%n(2)
            do i=1, self%n(1)
               el = el + src(i,j,srcidx(t3)) * self%slab_r(i,j,k)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, el, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      energy = el
   endif
   ! return: each block's planes, ghosts included, from their other owners (requests c1, c2, t3; windows wrapped
   ! modulo m), this rank's own straight from slab_r
   w = (nb(1) + 2 * gw) * (nb(2) + 2 * gw)
   sc = 0
   do b=1, nbl
      do k=1-gw, nb(3)+gw
         t3 = modulo(coords(3,b) * nb(3) + k - 1, self%m(3))
         if (self%zown(t3) /= me) sc(self%zown(t3)) = sc(self%zown(t3)) + 1
      enddo
   enddo
   call exchange_counts
   deallocate(shdr, rhdr, sval, rval)
   allocate(shdr(max(3*sum(sc),1)), rhdr(max(3*sum(rc),1)), sval(max(w*sum(rc),1)), rval(max(w*sum(sc),1)))
   fill = sd
   do b=1, nbl
      do k=1-gw, nb(3)+gw
         t3 = modulo(coords(3,b) * nb(3) + k - 1, self%m(3))
         r = self%zown(t3)
         if (r == me) cycle
         c = fill(r) ; fill(r) = fill(r) + 1
         shdr(3*c+1:3*c+3) = [coords(1,b), coords(2,b), t3]
      enddo
   enddo
   if (nr > 1) call MPI_ALLTOALLV(shdr, 3*sc, 3*sd, MPI_INTEGER, rhdr, 3*rc, 3*rd, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   do c=0, sum(rc)-1
      call get_plane(rhdr(3*c+1), rhdr(3*c+2), rhdr(3*c+3) - self%zr(me) + 1, sval(w*c+1:w*c+w))
   enddo
   if (nr > 1) call MPI_ALLTOALLV(sval, w*rc, w*rd, MPI_REAL8, rval, w*sc, w*sd, MPI_REAL8, MPI_COMM_WORLD, ierr)
   fill = sd
   e = 0
   do b=1, nbl
      do k=1-gw, nb(3)+gw
         t3 = modulo(coords(3,b) * nb(3) + k - 1, self%m(3))
         r = self%zown(t3)
         if (r == me) then
            call get_plane(coords(1,b), coords(2,b), t3 - self%zr(me) + 1, out(e+1:e+w))
         else
            c = fill(r) ; fill(r) = fill(r) + 1
            out(e+1:e+w) = rval(w*c+1:w*c+w)
         endif
         e = e + w
      enddo
   enddo
   tm(3) = MPI_WTIME()
   self%wall = self%wall + (tm(1:3) - tm(0:2))
   contains
      subroutine put_plane(c1, c2, l, v)
      !< The interior plane v of the block (c1, c2) into the local source plane l.
      integer(I4P), intent(in) :: c1, c2, l !< Block coordinates, local source plane.
      real(R8P),    intent(in) :: v(:)      !< Values, i fastest.
      integer(I4P)             :: jj        !< Counter.

      do jj=1, nb(2)
         src(c1*nb(1)+1:c1*nb(1)+nb(1),c2*nb(2)+jj,l) = v((jj-1)*nb(1)+1:jj*nb(1))
      enddo
      endsubroutine put_plane

      subroutine get_plane(c1, c2, l, v)
      !< The window of the block (c1, c2), ghosts included and wrapped modulo m, from the local plane l of slab_r.
      integer(I4P), intent(in)  :: c1, c2, l      !< Block coordinates, local plane.
      real(R8P),    intent(out) :: v(:)           !< Values, i fastest.
      integer(I4P)              :: ii, jj, o, wx, i0, y !< Counters, row width, first column (0-based), row.

      wx = nb(1) + 2 * gw ; i0 = c1 * nb(1) - gw
      o = 0
      do jj=1-gw, nb(2)+gw
         y = modulo(c2 * nb(2) + jj - 1, self%m(2)) + 1
         if (i0 >= 0 .and. i0 + wx <= self%m(1)) then
            v(o+1:o+wx) = self%slab_r(i0+1:i0+wx,y,l)
         else
            do ii=1, wx
               v(o+ii) = self%slab_r(modulo(i0 + ii - 1, self%m(1)) + 1,y,l)
            enddo
         endif
         o = o + wx
      enddo
      endsubroutine get_plane

      subroutine exchange_counts
      !< The receive counts from the send counts, and both displacements (items).
      call MPI_ALLTOALL(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      sd(0) = 0 ; rd(0) = 0
      do r=1, nr-1
         sd(r) = sd(r-1) + sc(r-1) ; rd(r) = rd(r-1) + rc(r-1)
      enddo
      endsubroutine exchange_counts
   endsubroutine solve_slabs

   subroutine slab_transform(self, factor)
   !< `slab_r <- factor L^-1 slab_r` on the distributed slabs: the planes' 2-D transforms, the transpose to z-lines,
   !< their 1-D transforms, the spectral solve, and back; the stages run where the backend keeps the slab arrays
   !< (`slab_update` moves `slab_r` there and back).
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   real(R8P),                 intent(in)    :: factor  !< Source factor.

   call self%slab_update(to_device=.true.)
   call self%slab_planes(forward=.true.)
   call self%slab_transpose(forward=.true.)
   call self%slab_lines(forward=.true.)
   call self%slab_spectral(factor=factor)
   call self%slab_lines(forward=.false.)
   call self%slab_transpose(forward=.false.)
   call self%slab_planes(forward=.false.)
   call self%slab_update(to_device=.false.)
   endsubroutine slab_transform

   subroutine global_potential(self, g)
   !< Rank 0: the transform array after the last solve, [m(1),m(2),m(3)] (the slabs gathered; `f_global` otherwise).
   !< For diagnostics (`test_poisson_fft --dump`); other ranks get an empty array.
   class(poisson_fft_object), intent(in)    :: self     !< The solver.
   real(R8P), allocatable,    intent(inout) :: g(:,:,:) !< The array.
   integer(I4P), allocatable :: rc(:), rd(:)             !< Counts and displacements.
   integer(I4P)              :: nr, r, nzl               !< Counters.
   integer(I4P)              :: ierr                     !< MPI error status.

   if (allocated(g)) deallocate(g)
   if (.not.self%slab) then
      if (mpih%myrank == 0) then
         g = self%f_global
      else
         allocate(g(0,0,0))
      endif
      return
   endif
   nr = mpih%procs_number
   allocate(rc(0:nr-1), rd(0:nr-1))
   do r=0, nr-1
      rc(r) = self%m(1) * self%m(2) * (self%zr(r+1) - self%zr(r))
      rd(r) = self%m(1) * self%m(2) * self%zr(r)
   enddo
   if (mpih%myrank == 0) then
      allocate(g(self%m(1),self%m(2),self%m(3)))
   else
      allocate(g(1,1,1))
   endif
   nzl = self%zr(mpih%myrank+1) - self%zr(mpih%myrank)
   call MPI_GATHERV(self%slab_r, self%m(1)*self%m(2)*nzl, MPI_REAL8, g, rc, rd, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   endsubroutine global_potential

   ! backend hooks
   ! private methods
   subroutine make_kernel(self)
   !< Rank 0, isolated: place the kernel `h^2 G(o)` at every offset `o`, wrapped into `(-m/2, m/2]` per direction, and
   !< keep its spectrum over `prod(m)` (`transform_kernel`).
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   integer(I4P)                             :: o(3)    !< Offset.
   integer(I4P)                             :: i, j, k !< Counters.

   do k=1, self%m(3)
      o(3) = k - 1 ; if (o(3) > self%m(3) / 2) o(3) = o(3) - self%m(3)
      do j=1, self%m(2)
         o(2) = j - 1 ; if (o(2) > self%m(2) / 2) o(2) = o(2) - self%m(2)
         do i=1, self%m(1)
            o(1) = i - 1 ; if (o(1) > self%m(1) / 2) o(1) = o(1) - self%m(1)
            self%f_global(i,j,k) = self%dxyz(1)**2 * lattice_green(o)
         enddo
      enddo
   enddo
   allocate(self%khat(self%m(1)/2+1,self%m(2),self%m(3)))
   call self%transform_kernel
   self%f_global = 0._R8P
   endsubroutine make_kernel

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
   self%plan_r2c = fftw_plan_dft_r2c_3d(int(self%m(3), c_int), int(self%m(2), c_int), int(self%m(1), c_int), &
                                        self%f_global, self%f_hat, FFTW_ESTIMATE)
   self%plan_c2r = fftw_plan_dft_c2r_3d(int(self%m(3), c_int), int(self%m(2), c_int), int(self%m(1), c_int), &
                                        self%f_hat, self%f_global, FFTW_ESTIMATE)
   if (.not.c_associated(self%plan_r2c) .or. .not.c_associated(self%plan_c2r)) &
      call mpih%error_stop(msg='poisson_fft_object%make_plans: FFTW could not make the plans')
#else
   call mpih%error_stop(msg='poisson_fft_object: built without an FFT library (build a *-fft mode: -D_FFTW on the '// &
                            'CPU, -D_CUFFT on FNL)')
#endif
   endsubroutine make_plans

   subroutine extend_global(self)
   !< Rank 0: the source extended over the m cells of every extended axis (cell index t = 0..m-1, array index t + 1),
   !< one axis after the other (the extensions are separable):
   !<   - ghost-Dirichlet (m = 2 (n + 1)): f(t) for t < n, 0 at t = n and t = m - 1 (the ghost centres), -f(2 n - t)
   !<     for n < t < m - 1;
   !<   - face kinds: f(t) for t < n; the reflection about the high face, s_hi f(2 n - 1 - t), for n <= t < 2 n; with
   !<     m = 4 n (the two faces differ) the reflection of [0, 2 n) about the low face, s_lo g(4 n - 1 - t), for
   !<     2 n <= t < 4 n. The sign s is -1 at a Dirichlet face, +1 at a Neumann one; the periodic extension then has
   !<     both reflection symmetries, and so does its periodic solution: the boundary conditions on [0, n).
   class(poisson_fft_object), intent(inout) :: self !< The solver.
   integer(I4P) :: t, d_, src   !< Cell index, axis, source cell index.
   real(R8P)    :: sg           !< Reflection sign.

   do d_=1, 3
      if (self%m(d_) == self%n(d_)) cycle ! a null or periodic axis
      do t=self%n(d_), self%m(d_)-1
         if (all(self%bc(:,d_) == FLAIL_BC_DIRICHLET_GHOST)) then
            if (t == self%n(d_) .or. t == self%m(d_)-1) then
               call set_plane(self%f_global, d_, t, 0, 0._R8P)
            else
               call set_plane(self%f_global, d_, t, 2*self%n(d_)-t, -1._R8P)
            endif
         elseif (t < 2*self%n(d_)) then
            sg = merge(-1._R8P, 1._R8P, self%bc(2,d_) == FLAIL_BC_DIRICHLET)
            call set_plane(self%f_global, d_, t, 2*self%n(d_)-1-t, sg)
         else
            sg = merge(-1._R8P, 1._R8P, self%bc(1,d_) == FLAIL_BC_DIRICHLET)
            src = 4*self%n(d_) - 1 - t
            call set_plane(self%f_global, d_, t, src, sg)
         endif
      enddo
   enddo
   endsubroutine extend_global

   pure subroutine set_plane(f, d_, t, src, sg)
   !< The plane t of axis d_ set to sg times the plane src (cell indices; sg = 0 zeroes it).
   real(R8P),    intent(inout) :: f(1:,1:,1:) !< The array.
   integer(I4P), intent(in)    :: d_, t, src  !< Axis, target and source cell indices.
   real(R8P),    intent(in)    :: sg          !< Sign (0: zero).

   select case(d_)
   case(1)
      if (sg == 0._R8P) then
         f(t+1,:,:) = 0._R8P
      else
         f(t+1,:,:) = sg * f(src+1,:,:)
      endif
   case(2)
      if (sg == 0._R8P) then
         f(:,t+1,:) = 0._R8P
      else
         f(:,t+1,:) = sg * f(:,src+1,:)
      endif
   case(3)
      if (sg == 0._R8P) then
         f(:,:,t+1) = 0._R8P
      else
         f(:,:,t+1) = sg * f(:,:,src+1)
      endif
   endselect
   endsubroutine set_plane

   subroutine fill_ghosts(self, u, ng)
   !< Rank 0: the `ng` ghost layers of the padded global array `u` by the faces' rules, one axis after the other:
   !< periodic (and null) the wrap, Dirichlet and Neumann at the face the odd and even reflections about it,
   !< ghost-Dirichlet the odd reflection about the ghost centre (zero there), zero-ghost zero.
   class(poisson_fft_object), intent(in)    :: self !< The solver.
   integer(I4P),              intent(in)    :: ng   !< Ghost layers.
   real(R8P),                 intent(inout) :: u(1-ng:,1-ng:,1-ng:) !< Padded array [1-ng:n+ng per axis].
   integer(I4P)                             :: d, j, lo, hi, n !< Axis, layer, ghost indices, cells.

   do d=1, 3
      n = self%n(d)
      do j=1, ng
         lo = 1 - j ; hi = n + j
         if (all(self%bc(:,d) == FLAIL_BC_PERIODIC)) then
            call copy_plane(d, lo, modulo(lo - 1, n) + 1, 1._R8P)
            call copy_plane(d, hi, modulo(hi - 1, n) + 1, 1._R8P)
            cycle
         endif
         select case(self%bc(1,d))
         case(FLAIL_BC_DIRICHLET)
            call copy_plane(d, lo, j, -1._R8P)
         case(FLAIL_BC_NEUMANN)
            call copy_plane(d, lo, j, 1._R8P)
         case(FLAIL_BC_DIRICHLET_GHOST)
            if (j == 1) then
               call copy_plane(d, lo, 1, 0._R8P)
            else
               call copy_plane(d, lo, j-1, -1._R8P)
            endif
         case default ! zero ghosts
            call copy_plane(d, lo, 1, 0._R8P)
         endselect
         select case(self%bc(2,d))
         case(FLAIL_BC_DIRICHLET)
            call copy_plane(d, hi, n+1-j, -1._R8P)
         case(FLAIL_BC_NEUMANN)
            call copy_plane(d, hi, n+1-j, 1._R8P)
         case(FLAIL_BC_DIRICHLET_GHOST)
            if (j == 1) then
               call copy_plane(d, hi, n, 0._R8P)
            else
               call copy_plane(d, hi, n+2-j, -1._R8P)
            endif
         case default
            call copy_plane(d, hi, n, 0._R8P)
         endselect
      enddo
   enddo
   contains
      subroutine copy_plane(d_, t, src, sg)
      !< The plane t of axis d_ set to sg times the plane src (sg = 0 zeroes it).
      integer(I4P), intent(in) :: d_, t, src !< Axis, target and source indices.
      real(R8P),    intent(in) :: sg         !< Sign.

      select case(d_)
      case(1)
         u(t,:,:) = sg * u(src,:,:)
      case(2)
         u(:,t,:) = sg * u(:,src,:)
      case(3)
         u(:,:,t) = sg * u(:,:,src)
      endselect
      endsubroutine copy_plane
   endsubroutine fill_ghosts

   subroutine apply_wide(self, x, y)
   !< Rank 0: `y = D(G) x` on the global interior, the ghosts by the faces' rules (`fill_ghosts`).
   class(poisson_fft_object), intent(in)  :: self          !< The solver.
   real(R8P),                 intent(in)  :: x(1:,1:,1:)   !< Interior values [n].
   real(R8P),                 intent(out) :: y(1:,1:,1:)   !< Result [n].
   real(R8P), allocatable                 :: u(:,:,:)      !< Padded copy.
   integer(I4P)                           :: p, i, j, k, o !< Reach, counters.
   real(R8P)                              :: w(3)          !< 1 / dx^2 per axis (0 along a null one).

   p = self%dg_order
   allocate(u(1-p:self%n(1)+p,1-p:self%n(2)+p,1-p:self%n(3)+p))
   u = 0._R8P
   u(1:self%n(1),1:self%n(2),1:self%n(3)) = x
   call self%fill_ghosts(u, p)
   w = 0._R8P
   where (.not.self%null_xyz) w = 1._R8P / self%dxyz**2
   do k=1, self%n(3)
      do j=1, self%n(2)
         do i=1, self%n(1)
            y(i,j,k) = 0._R8P
            do o=-p, p
               y(i,j,k) = y(i,j,k) + self%stencil(o) * (w(1) * u(i+o,j,k) + w(2) * u(i,j+o,k) + w(3) * u(i,j,k+o))
            enddo
         enddo
      enddo
   enddo
   endsubroutine apply_wide

   subroutine precondition(self, r, z)
   !< Rank 0: `z = M^-1 r`, M the wide p = 2 operator with the same faces: the 2^z parity sub-problems, each solved by
   !< its sub-solver's transform (the extension of its kinds, then the division by its eigenvalues).
   class(poisson_fft_object), intent(inout) :: self        !< The solver.
   real(R8P),                 intent(in)    :: r(1:,1:,1:) !< Residual [n].
   real(R8P),                 intent(inout) :: z(1:,1:,1:) !< Correction [n].
   integer(I4P)                             :: c, i, j, k, ii, jj, kk !< Counters, global indices.

   do c=1, self%nsub
      associate(sb=>self%sub(c)%s, par=>self%sub_par(:,c))
      sb%f_global = 0._R8P
      do k=1, sb%n(3)
         kk = merge(k, 2*k-1+par(3), par(3) < 0)
         do j=1, sb%n(2)
            jj = merge(j, 2*j-1+par(2), par(2) < 0)
            do i=1, sb%n(1)
               ii = merge(i, 2*i-1+par(1), par(1) < 0)
               sb%f_global(i,j,k) = r(ii,jj,kk)
            enddo
         enddo
      enddo
      if (sb%extended) call sb%extend_global
      call sb%transform(factor=1._R8P)
      do k=1, sb%n(3)
         kk = merge(k, 2*k-1+par(3), par(3) < 0)
         do j=1, sb%n(2)
            jj = merge(j, 2*j-1+par(2), par(2) < 0)
            do i=1, sb%n(1)
               ii = merge(i, 2*i-1+par(1), par(1) < 0)
               z(ii,jj,kk) = sb%f_global(i,j,k)
            enddo
         enddo
      enddo
      endassociate
   enddo
   endsubroutine precondition

   subroutine solve_pcg(self, factor)
   !< Rank 0: `f_global <- u`, `D(G) u = factor f_global` with zero-ghost faces (E2b), by conjugate gradients
   !< preconditioned by `precondition`, from zero, to `tolerance` (the max-norm residual relative to the source).
   !< D(G) and M are negative definite here (the zero ghosts remove every null mode), which the iteration takes as is.
   class(poisson_fft_object), intent(inout) :: self !< The solver.
   real(R8P),                 intent(in)    :: factor !< Source factor.
   real(R8P), allocatable :: f(:,:,:), u(:,:,:), r(:,:,:), z(:,:,:), d(:,:,:), q(:,:,:) !< Vectors.
   real(R8P)              :: fmax, rz, rz_new, alpha !< Scalars.
   integer(I4P)           :: it !< Counter.

   allocate(f, mold=self%f_global)
   f = factor * self%f_global
   if (self%free) call self%subtract_free_ghosts(f) ! D(G)_0 u = factor f - D(G) u_g (E3c)
   allocate(u, r, z, d, q, mold=f)
   u = 0._R8P ; r = f ; z = 0._R8P
   fmax = maxval(abs(f))
   self%iterations = 0 ; self%residual = 0._R8P
   if (fmax > 0._R8P) then
      call self%precondition(r, z)
      d = z
      rz = sum(r * z)
      self%residual = 1._R8P
      do it=1, self%max_iterations
         call self%apply_wide(d, q)
         alpha = rz / sum(d * q)
         u = u + alpha * d
         r = r - alpha * q
         self%iterations = it
         self%residual = maxval(abs(r)) / fmax
         if (self%residual <= self%tolerance) exit
         call self%precondition(r, z)
         rz_new = sum(r * z)
         d = z + (rz_new / rz) * d
         rz = rz_new
      enddo
      if (self%residual > self%tolerance) print '(A)', mpih%myrankstr//'warning: poisson_fft_object: '// &
         trim(str(self%max_iterations, .true.))//' CG iterations left the residual at '//trim(str(self%residual))// &
         ' (tolerance '//trim(str(self%tolerance))//')'
   endif
   self%f_global = u
   endsubroutine solve_pcg

   subroutine free_ghosts(self, factor, gd)
   !< Rank 0, E3c: `g_free`, to depth `gd`, the ghost slab of every free-space face (the cells outside along that face's
   !< axis, inside along the others): the continuum free-space potential of the gathered source,
   !< `g(x) = -sum_y factor f(y) dV / (4 pi |x - y|)` over the interior cell centres, as PRISM's exact/open face.
   class(poisson_fft_object), intent(inout) :: self   !< The solver.
   real(R8P),                 intent(in)    :: factor !< Source factor.
   integer(I4P),              intent(in)    :: gd     !< Ghost layers.
   real(R8P)                                :: c      !< -factor dV / (4 pi).
   real(R8P)                                :: s, dz2, dyz2 !< Partial sum, squared distances.
   integer(I4P)                             :: lo(3), hi(3) !< The slab.
   integer(I4P)                             :: d, side, i, j, k, ii, jj, kk !< Counters.

   if (allocated(self%g_free)) deallocate(self%g_free)
   allocate(self%g_free(1-gd:self%n(1)+gd,1-gd:self%n(2)+gd,1-gd:self%n(3)+gd))
   self%g_free = 0._R8P
   c = -factor * product(self%dxyz) / (4._R8P * acos(-1._R8P))
   do d=1, 3
      if (self%null_xyz(d)) cycle
      do side=1, 2
         if (self%bc(side,d) /= FLAIL_BC_FREE_GHOSTS) cycle
         lo = 1 ; hi = self%n
         if (side == 1) then
            lo(d) = 1 - gd ; hi(d) = 0
         else
            lo(d) = self%n(d) + 1 ; hi(d) = self%n(d) + gd
         endif
         do k=lo(3), hi(3)
            do j=lo(2), hi(2)
               do i=lo(1), hi(1)
                  s = 0._R8P
                  do kk=1, self%n(3)
                     dz2 = (real(k - kk, R8P) * self%dxyz(3))**2
                     do jj=1, self%n(2)
                        dyz2 = dz2 + (real(j - jj, R8P) * self%dxyz(2))**2
                        do ii=1, self%n(1)
                           s = s + self%f_global(ii,jj,kk) / sqrt(dyz2 + (real(i - ii, R8P) * self%dxyz(1))**2)
                        enddo
                     enddo
                  enddo
                  self%g_free(i,j,k) = c * s
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine free_ghosts

   subroutine subtract_free_ghosts(self, f)
   !< Rank 0, E3c: `f <- f - D(G) u_g` on the interior, u_g the free-space ghost slabs (`g_free`) and zero inside.
   class(poisson_fft_object), intent(in)    :: self        !< The solver.
   real(R8P),                 intent(inout) :: f(1:,1:,1:) !< Right-hand side [n].
   integer(I4P)                             :: p, i, j, k, o !< Reach, counters.
   real(R8P)                                :: w(3), y   !< 1 / dx^2 per axis (0 along a null one), D(G) u_g.

   p = self%dg_order
   w = 0._R8P
   where (.not.self%null_xyz) w = 1._R8P / self%dxyz**2
   associate(g=>self%g_free)
   do k=1, self%n(3)
      do j=1, self%n(2)
         do i=1, self%n(1)
            y = 0._R8P
            do o=-p, p
               y = y + self%stencil(o) * (w(1) * g(i+o,j,k) + w(2) * g(i,j+o,k) + w(3) * g(i,j,k+o))
            enddo
            f(i,j,k) = f(i,j,k) - y
         enddo
      enddo
   enddo
   endassociate
   endsubroutine subtract_free_ghosts

   subroutine put_free_ghosts(self, u, ng)
   !< Rank 0, E3c: the `ng` ghost layers of every free-space face of the padded array `u` set to `g_free`.
   class(poisson_fft_object), intent(in)    :: self !< The solver.
   integer(I4P),              intent(in)    :: ng   !< Ghost layers.
   real(R8P),                 intent(inout) :: u(1-ng:,1-ng:,1-ng:) !< Padded array [1-ng:n+ng per axis].
   integer(I4P)                             :: lo(3), hi(3), d, side !< The slab, counters.

   do d=1, 3
      if (self%null_xyz(d)) cycle
      do side=1, 2
         if (self%bc(side,d) /= FLAIL_BC_FREE_GHOSTS) cycle
         lo = 1 ; hi = self%n
         if (side == 1) then
            lo(d) = 1 - ng ; hi(d) = 0
         else
            lo(d) = self%n(d) + 1 ; hi(d) = self%n(d) + ng
         endif
         u(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = self%g_free(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      enddo
   enddo
   endsubroutine put_free_ghosts

   subroutine transform(self, factor)
   !< Rank 0: `f_global <- factor * L_h^-1 f_global`. Periodic: the modes of zero eigenvalue (the mean) set to zero;
   !< isolated: the convolution with the kernel (`khat`).
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   real(R8P),                 intent(in)    :: factor  !< Source factor.
   real(R8P)                                :: scale   !< Factor over the unnormalised transform pair.
   real(R8P)                                :: lam     !< Eigenvalue of the mode.
   integer(I4P)                             :: i, j, k !< Counters.

#ifdef _FFTW
   call fftw_execute_dft_r2c(self%plan_r2c, self%f_global, self%f_hat)
   if (self%isolated) then
      self%f_hat = self%f_hat * (factor * self%khat)
      call fftw_execute_dft_c2r(self%plan_c2r, self%f_hat, self%f_global)
      return
   endif
   scale = factor / real(product(self%nt), R8P)
   do k=1, self%nt(3)
      do j=1, self%nt(2)
         do i=1, self%nt(1)/2+1
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

   subroutine transform_kernel(self)
   !< Rank 0, isolated: `khat <- FFT(f_global) / prod(m)`, real (the kernel is even in every direction).
   class(poisson_fft_object), intent(inout) :: self !< The solver.

#ifdef _FFTW
   call fftw_execute_dft_r2c(self%plan_r2c, self%f_global, self%f_hat)
   self%khat = real(self%f_hat, R8P) / real(product(self%m), R8P)
#else
   call mpih%error_stop(msg='poisson_fft_object%transform_kernel: built without an FFT library')
#endif
   endsubroutine transform_kernel

   pure function root_cells_max() result(n)
   !< The most cells the multigrid leaves to this backend's exact root solve, below which it does not coarsen inside the
   !< blocks (#84 E5c): 32^3, the redundant root's limit. On the host a 128^3 FFT on one core costs more per V-cycle
   !< than smoothing the level over 8 or more ranks (medooza: the periodic multigrid on 3.9 M cells, 9.2 -> 12.7 s at
   !< np 1, 2.6 -> 2.3 s at np 16). On cuFFT too since the multigrid runs on the device (#84 E6): 3.9 M AMR cells,
   !< 1 / 2 / 4 GPUs, isolated 4.45 / 2.55 / 2.92 s with a 128^3 root, 0.85 / 0.52 / 0.43 s with 32^3; periodic
   !< 0.80 / 0.45 / 0.37 -> 0.33 / 0.19 / 0.12 s. A backend may override it.
   integer(I4P) :: n !< Cells.

   n = 32768_I4P
   endfunction root_cells_max

   subroutine slab_make_plans(self)
   !< The slab plans (#84 E5d): one 2-D plane (r2c, c2r) and one row of z-lines (h1 lines of m(3) cells, c2c), each
   !< made once on its own buffer and executed on every plane and row in turn (FFTW_ESTIMATE: deterministic, and the
   !< arithmetic of a value does not depend on how many planes or rows a rank holds).
   class(poisson_fft_object), intent(inout) :: self !< The solver.
   integer(c_int)                           :: nl(1) !< Line length.

   allocate(self%w_r(self%m(1),self%m(2)), self%w_c(self%h1,self%m(2)), self%w_l(self%m(3),self%h1))
   self%w_r = 0._R8P ; self%w_c = (0._R8P, 0._R8P) ; self%w_l = (0._R8P, 0._R8P)
#ifdef _FFTW
   self%plan_pr2c = fftw_plan_dft_r2c_2d(int(self%m(2), c_int), int(self%m(1), c_int), self%w_r, self%w_c, &
                                         FFTW_ESTIMATE)
   self%plan_pc2r = fftw_plan_dft_c2r_2d(int(self%m(2), c_int), int(self%m(1), c_int), self%w_c, self%w_r, &
                                         FFTW_ESTIMATE)
   nl = int(self%m(3), c_int)
   self%plan_lfwd = fftw_plan_many_dft(1_c_int, nl, int(self%h1, c_int), self%w_l, c_null_ptr, 1_c_int, nl(1), &
                                       self%w_l, c_null_ptr, 1_c_int, nl(1), FFTW_FORWARD, FFTW_ESTIMATE)
   self%plan_lbwd = fftw_plan_many_dft(1_c_int, nl, int(self%h1, c_int), self%w_l, c_null_ptr, 1_c_int, nl(1), &
                                       self%w_l, c_null_ptr, 1_c_int, nl(1), FFTW_BACKWARD, FFTW_ESTIMATE)
   if (.not.(c_associated(self%plan_pr2c) .and. c_associated(self%plan_pc2r) .and. c_associated(self%plan_lfwd) &
       .and. c_associated(self%plan_lbwd))) call mpih%error_stop(msg='poisson_fft_object%slab_make_plans: FFTW '// &
       'could not make the plans')
#else
   nl = 0_c_int
   call mpih%error_stop(msg='poisson_fft_object: built without an FFT library (build a *-fft mode: -D_FFTW on the '// &
                            'CPU, -D_CUFFT on FNL)')
#endif
   endsubroutine slab_make_plans

   subroutine slab_destroy_plans(self)
   !< Destroy the slab plans.
   class(poisson_fft_object), intent(inout) :: self !< The solver.

#ifdef _FFTW
   if (c_associated(self%plan_pr2c)) call fftw_destroy_plan(self%plan_pr2c)
   if (c_associated(self%plan_pc2r)) call fftw_destroy_plan(self%plan_pc2r)
   if (c_associated(self%plan_lfwd)) call fftw_destroy_plan(self%plan_lfwd)
   if (c_associated(self%plan_lbwd)) call fftw_destroy_plan(self%plan_lbwd)
#endif
   self%plan_pr2c = c_null_ptr ; self%plan_pc2r = c_null_ptr
   self%plan_lfwd = c_null_ptr ; self%plan_lbwd = c_null_ptr
   endsubroutine slab_destroy_plans

   subroutine slab_planes(self, forward)
   !< The 2-D transform of every local plane: forward slab_r -> slab_c (r2c), backward slab_c -> slab_r (c2r), one plane
   !< at a time through the plane buffers.
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   logical,                   intent(in)    :: forward !< Direction.
   integer(I4P)                             :: k       !< Counter.

#ifdef _FFTW
   do k=1, self%zr(mpih%myrank+1) - self%zr(mpih%myrank)
      if (forward) then
         self%w_r = self%slab_r(:,:,k)
         call fftw_execute_dft_r2c(self%plan_pr2c, self%w_r, self%w_c)
         self%slab_c(:,:,k) = self%w_c
      else
         self%w_c = self%slab_c(:,:,k)
         call fftw_execute_dft_c2r(self%plan_pc2r, self%w_c, self%w_r)
         self%slab_r(:,:,k) = self%w_r
      endif
   enddo
#else
   k = merge(1, 0, forward)
   call mpih%error_stop(msg='poisson_fft_object%slab_planes: built without an FFT library')
#endif
   endsubroutine slab_planes

   subroutine slab_lines(self, forward)
   !< The 1-D transform of every local z-line, in place in pen, one row of h1 lines at a time through the row buffer.
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   logical,                   intent(in)    :: forward !< Direction.
   integer(I4P)                             :: j       !< Counter.

#ifdef _FFTW
   do j=1, self%yr(mpih%myrank+1) - self%yr(mpih%myrank)
      self%w_l = self%pen(:,:,j)
      if (forward) then
         call fftw_execute_dft(self%plan_lfwd, self%w_l, self%w_l)
      else
         call fftw_execute_dft(self%plan_lbwd, self%w_l, self%w_l)
      endif
      self%pen(:,:,j) = self%w_l
   enddo
#else
   j = merge(1, 0, forward)
   call mpih%error_stop(msg='poisson_fft_object%slab_lines: built without an FFT library')
#endif
   endsubroutine slab_lines

   subroutine slab_transpose(self, forward)
   !< forward: slab_c(x, y, z local) -> pen(z, x, y local); backward the reverse (MPI_ALLTOALLV, complex). This rank's
   !< own block is copied directly.
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   logical,                   intent(in)    :: forward !< Direction.
   integer(I4P), allocatable :: sc(:), sd(:), rc(:), rd(:) !< Counts and displacements (values).
   complex(R8P), allocatable :: sb(:), rb(:)               !< Buffers.
   integer(I4P)              :: nr, me, r, c, x, y, z, nzl, nyl, z0, y0 !< Counters, sizes, first plane and row.
   integer(I4P)              :: ierr                       !< MPI error status.

   nr = mpih%procs_number ; me = mpih%myrank
   nzl = self%zr(me+1) - self%zr(me) ; nyl = self%yr(me+1) - self%yr(me)
   z0 = self%zr(me) ; y0 = self%yr(me)
   allocate(sc(0:nr-1), sd(0:nr-1), rc(0:nr-1), rd(0:nr-1))
   do r=0, nr-1
      if (forward) then
         sc(r) = self%h1 * (self%yr(r+1) - self%yr(r)) * nzl
         rc(r) = self%h1 * nyl * (self%zr(r+1) - self%zr(r))
      else
         sc(r) = self%h1 * nyl * (self%zr(r+1) - self%zr(r))
         rc(r) = self%h1 * (self%yr(r+1) - self%yr(r)) * nzl
      endif
   enddo
   sc(me) = 0 ; rc(me) = 0
   sd(0) = 0 ; rd(0) = 0
   do r=1, nr-1
      sd(r) = sd(r-1) + sc(r-1) ; rd(r) = rd(r-1) + rc(r-1)
   enddo
   allocate(sb(max(sum(sc),1)), rb(max(sum(rc),1)))
   c = 0
   do r=0, nr-1
      if (r == me) cycle
      if (forward) then
         do z=1, nzl
            do y=self%yr(r)+1, self%yr(r+1)
               do x=1, self%h1
                  c = c + 1 ; sb(c) = self%slab_c(x,y,z)
               enddo
            enddo
         enddo
      else
         do z=self%zr(r)+1, self%zr(r+1)
            do y=1, nyl
               do x=1, self%h1
                  c = c + 1 ; sb(c) = self%pen(z,x,y)
               enddo
            enddo
         enddo
      endif
   enddo
   if (nr > 1) call MPI_ALLTOALLV(sb, sc, sd, MPI_DOUBLE_COMPLEX, rb, rc, rd, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
   ! this rank's own block
   do y=1, nyl
      do x=1, self%h1
         do z=1, nzl
            if (forward) then
               self%pen(z0+z,x,y) = self%slab_c(x,y0+y,z)
            else
               self%slab_c(x,y0+y,z) = self%pen(z0+z,x,y)
            endif
         enddo
      enddo
   enddo
   c = 0
   do r=0, nr-1
      if (r == me) cycle
      if (forward) then
         do z=self%zr(r)+1, self%zr(r+1)
            do y=1, nyl
               do x=1, self%h1
                  c = c + 1 ; self%pen(z,x,y) = rb(c)
               enddo
            enddo
         enddo
      else
         do z=1, nzl
            do y=self%yr(r)+1, self%yr(r+1)
               do x=1, self%h1
                  c = c + 1 ; self%slab_c(x,y,z) = rb(c)
               enddo
            enddo
         enddo
      endif
   enddo
   endsubroutine slab_transpose

   subroutine slab_spectral(self, factor)
   !< The spectral solve on this rank's z-lines: the division by the eigenvalues (the zero modes dropped) or, isolated,
   !< the product with the kernel's spectrum.
   class(poisson_fft_object), intent(inout) :: self    !< The solver.
   real(R8P),                 intent(in)    :: factor  !< Source factor.
   real(R8P)                                :: scale   !< Factor over the unnormalised transform pair.
   real(R8P)                                :: lam     !< Eigenvalue of the mode.
   integer(I4P)                             :: i, j, k, nyl, y0 !< Counters, local rows, first row.

   nyl = self%yr(mpih%myrank+1) - self%yr(mpih%myrank) ; y0 = self%yr(mpih%myrank)
   if (self%isolated) then
      do j=1, nyl
         do i=1, self%h1
            do k=1, self%m(3)
               self%pen(k,i,j) = self%pen(k,i,j) * (factor * self%khat_pen(k,i,j))
            enddo
         enddo
      enddo
   else
      scale = factor / real(product(int(self%nt, I8P)), R8P)
      do j=1, nyl
         do i=1, self%h1
            do k=1, self%m(3)
               lam = self%lam1(i-1) + self%lam2(y0+j-1) + self%lam3(k-1)
               if (lam == 0._R8P) then
                  self%pen(k,i,j) = (0._R8P, 0._R8P)
               else
                  self%pen(k,i,j) = self%pen(k,i,j) * (scale / lam)
               endif
            enddo
         enddo
      enddo
   endif
   endsubroutine slab_spectral

   subroutine slab_keep_kernel(self)
   !< Isolated, at setup: `khat_pen <- real(pen) / prod(m)`, the kernel's spectrum (real: the kernel is even).
   class(poisson_fft_object), intent(inout) :: self !< The solver.

   self%khat_pen = real(self%pen, R8P) / real(product(int(self%m, I8P)), R8P)
   endsubroutine slab_keep_kernel

   subroutine slab_update(self, to_device)
   !< Move `slab_r` to (`to_device`) or from where the stages run: nothing on the host.
   class(poisson_fft_object), intent(inout) :: self      !< The solver.
   logical,                   intent(in)    :: to_device !< Direction.

   if (to_device .or. .not.self%slab) return
   endsubroutine slab_update

   ! module helpers
   pure function smooth_size(k) result(m)
   !< The smallest even size `m >= k` with no prime factor above 13: FFTW has hard-coded transforms for these factors
   !< (#84 E5b, medooza: a 7-smooth rule moved the James step's 260 = 2^2 5 13 to 270 and cost it 15%; cuFFT took the
   !< same time on 262 and 270).
   integer(I4P), intent(in) :: k !< Least size.
   integer(I4P)             :: m !< Size.
   integer(I4P), parameter  :: PRIMES(6)=[2_I4P, 3_I4P, 5_I4P, 7_I4P, 11_I4P, 13_I4P] !< Allowed factors.
   integer(I4P)             :: r !< Remainder after the allowed factors.
   integer(I4P)             :: p !< Counter.

   m = k + mod(k, 2_I4P)
   do
      r = m
      do p=1, size(PRIMES)
         do while (mod(r, PRIMES(p)) == 0)
            r = r / PRIMES(p)
         enddo
      enddo
      if (r == 1) return
      m = m + 2_I4P
   enddo
   endfunction smooth_size

   pure function axis_eigenvalue(k, nt, dx, dg_order) result(lam)
   !< The eigenvalue part of mode k of an axis on the nt-cell periodic transform: the compact operator,
   !< `(2 cos(theta) - 2) / dx^2`, or the wide D(G) of order p, `-s_p(theta)^2 / dx^2`, exactly 0 at theta = 0 and pi.
   integer(I4P), intent(in) :: k, nt    !< Mode, transform cells.
   real(R8P),    intent(in) :: dx       !< Spacing.
   integer(I4P), intent(in) :: dg_order !< 0 compact, else the order of D(G).
   real(R8P)                :: lam      !< Eigenvalue part.
   real(R8P)                :: th, s    !< Phase, symbol of dx D1.
   integer(I4P)             :: j        !< Counter.

   if (dg_order == 0_I4P) then
      lam = (2._R8P * cos(2._R8P * PI * k / nt) - 2._R8P) / dx**2
      return
   endif
   lam = 0._R8P
   if (k == 0_I4P .or. 2_I4P * k == nt) return ! theta = 0 or pi: s_p vanishes exactly
   th = 2._R8P * PI * k / nt
   s = 0._R8P
   do j=1, dg_order/2
      s = s + 2._R8P * D1_COEF(j,dg_order/2) * sin(j * th)
   enddo
   lam = -(s * s) / dx**2
   endfunction axis_eigenvalue

   pure function dg_stencil(p) result(c)
   !< The 1-D stencil of D1 D1 (times h^2), offsets -p..p, the centred D1 of order p (`D1_COEF`) composed with itself.
   integer(I4P), intent(in) :: p          !< Order.
   real(R8P)                :: c(-8:8)    !< Stencil.
   real(R8P)                :: a(-4:4)    !< D1 (times h).
   integer(I4P)             :: j, o1, o2  !< Counters.

   a = 0._R8P
   do j=1, p/2
      a(j) = D1_COEF(j,p/2) ; a(-j) = -D1_COEF(j,p/2)
   enddo
   c = 0._R8P
   do o1=-p/2, p/2
      do o2=-p/2, p/2
         c(o1+o2) = c(o1+o2) + a(o1) * a(o2)
      enddo
   enddo
   endfunction dg_stencil
endmodule adam_poisson_fft_object
