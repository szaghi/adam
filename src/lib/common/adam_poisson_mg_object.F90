!< ADAM, multigrid-preconditioned conjugate-gradient Poisson solver on AMR grids (issue #81, M8-P5).
module adam_poisson_mg_object
!< ADAM, multigrid-preconditioned conjugate-gradient Poisson solver on AMR grids (issue #81, M8-P5).
!<
!< Solves `L_h phi = factor (f - mean(f))` on the leaf cells of a 2:1-balanced octree (quadtree in 2-D, `null_z`),
!< periodic in every active direction, with `L_h` a composite second-order Laplacian that is **symmetric and
!< consistent** (issue #81, P5a):
!<   - between leaves of one level, the 7-point (5-point) stencil;
!<   - across a 2:1 face, the **mimetic mean flux**: every fine cell `f_k` (k = 1..n, n = 2^(d-1)) on a coarse cell
!<     C's face gets the same flux `(2/3) w (mean_k phi_fk - phi_C)`, `w` the fine face weight (area over the fine
!<     spacing), and C gets their sum. It is the hybrid mimetic finite-difference flux of the coarse cell seen as a
!<     polygon with hanging nodes, its face unknowns eliminated: consistency (every sub-face flux exact for linear phi)
!<     forces the coarse cell's sub-face block to rank one (it sees only the mean of the face unknowns), and the
!<     elimination gives the mean of the fine cells. Exact for linear phi face by face, so the truncation error stays
!<     bounded at the seams (phi second order, the seam face gradient first order); symmetric in the V-weighted form
!<     (rows times the cell volume), so the solution map `phi = K f` has `V K` symmetric, which exact energy
!<     conservation of self-gravity needs.
!<   Rejected on the way (#81, P5-0 and P5a): Tomida & Stone's flux matching (consistent, not symmetric: the fallback);
!<   its symmetrisation (B + B^T) / 2 (symmetric, inconsistent: O(1/h) truncation, a seam gravity that does not
!<   converge); the Galerkin operator P^T A_fine P (symmetric, first order).
!<
!< **The solve** is flexible preconditioned conjugate gradients (Notay's FCG(1)) on the composite operator: correctness
!< lives in its matrix-vector product only. The preconditioner is one V-cycle over the level hierarchy: level `l` holds the leaves
!< at `l` and a **shadow** block for every tree node at `l` that has children (ADAM keeps leaves only); red-black
!< Gauss-Seidel sweeps down and up on the plain level Laplacian, a face to a cell outside the level reading the
!< coarser correction (piecewise constant) on the way up; restriction the sum of the children's residuals, prolongation
!< piecewise constant; the coarsest tree level (every block at the minimum leaf level: the whole domain) coarsened
!< further **inside the blocks** (issue #84, E5c): each level below it halves the cells per block along every active
!< direction, the same blocks on the same ranks, the plain level Laplacian with the face kinds, the same smoother and
!< transfers, while the cells per block are even and the domain has more than `root_cells` cells at the level
!< (default: the root solver's `root_cells_max`, 32^3 on the host FFT, 128^3 on cuFFT, which outruns the host
!< smoother); the last of them, the root (`lbot`), is solved exactly by the FFT solver
!< (`adam_poisson_fft_object%solve_blocks`; of the caller's backend type with `initialize(coarse=)`): on **every
!< rank** (`redundant`), from the source gathered by `MPI_ALLGATHERV`, when it has at most `REDUNDANT_ROOT_MAX` (32^3)
!< cells, so no rank-0 gather and scatter per V-cycle; a larger root (cuFFT keeps up to 128^3) on rank 0 as before,
!< since every rank transforming it contends for the node's memory (FNL, 4 GPUs, the isolated 128^3 coarsest level:
!< 4.7 -> 6.4 s per solve). A coarsest tree level of at most `root_cells` cells is its own root, solved as before. Measured
!< (`test_poisson_mg --sub`, every level down to one cell per block): 1 -> 14 FCG iterations on a uniform 1:2:3-cell
!< box (point smoothing on anisotropic cells), +1 to +13 on AMR (24-47 -> 25-53). Medooza CPU, 3.9 M AMR cells per
!< solve: isolated 50.7 -> 31.2 s at np 1 and 39.9 -> 6.5 s at np 16 (its 258^3 coarsest level and 42 exact solves
!< were 70% of it); periodic 9.2 -> 12.7 s at np 1, even at np 8, 2.6 -> 2.3 s at np 16 (the generic CSR smoother
!< on the base level costs more than one FFT of it on few ranks; a stencil smoother is #84 E6). It
!< stops at `tolerance` or at the round-off floor of the source, `100 eps max|factor f| / max|factor (f - mean f)|`,
!< whichever is larger (a small perturbation of a uniform density is known only to that floor, #81 P5b).
!<
!< **Distributed and assembled.** Every rank gathers the leaf list once (block coordinates, owner, local index),
!< derives the shadow blocks and their owners (the owner of the first leaf below, in the gathered order: the same on
!< every rank), and assembles, from the geometry only, the rows of its own cells: the operator, the ghost fill, the level
!< operators and transfers. Columns are global cell keys; every operator owns a halo-exchange plan (the remote columns
!< it reads, fetched from their owners by `MPI_ALLTOALLV`). Correct for any rank count. Cell keys carry the level
!< offset by `KEY_LEVEL_OFFSET` (the levels inside the blocks may be negative).
!<
!< **Domain boundaries (issue #87, P1; per face, issue #84, E1b).** Every active axis is periodic (both faces) or
!< bounded, and each face of a bounded axis has a kind (`adam_flail_parameters`; `initialize(..., bc=)`, by default
!< ghost-Dirichlet on every bounded axis, the #87 behaviour). A boundary face of a leaf cell C (weight `w`, the face
!< area over the cell spacing `h`) adds to C's row, and maps its boundary value `g` onto the right-hand side:
!<   - `FLAIL_BC_DIRICHLET_GHOST`, `g` the value at the centre of the first ghost cell: `-w phi_C`, `w g`;
!<   - `FLAIL_BC_DIRICHLET`, `g` the value at the face (`phi_ghost = 2 g - phi_C`): `-2 w phi_C`, `2 w g`;
!<   - `FLAIL_BC_NEUMANN`, `g` the outward normal derivative at the face (`phi_ghost = phi_C + h g`): nothing, `w h g`.
!< So the operator `A` keeps the homogeneous rows, symmetric, and `solve(..., ghost=g)` solves `A phi = factor V f -
!< B g`, `B` the map of the boundary values onto the boundary rows (`apply_b`; its transpose `apply_bt`, the boundary
!< rows onto the faces). The faces are listed per rank (`nbnd`, `bnd_row`, `bnd_w` the weight in B, `bnd_kind`,
!< `bnd_dir` the outward axis signed, `bnd_x` where g lives: the ghost centre or the face centre). Without a
!< Dirichlet face the operator is singular, as on a periodic grid (`singular`): the solve removes the mean of the
!< right-hand side `factor V f - B g` (the compatibility condition, exact for compatible data) and of the solution. The
!< V-cycle's levels take the same face kinds, homogeneous, and its coarsest level the FFT solver with them (odd or even
!< extensions; issue #84, E1a); on a uniform grid the two operators coincide and one V-cycle is the exact solve. The
!< two faces of an axis are both periodic or neither, and both ghost-Dirichlet or neither (the FFT solver's scope).
!<
!< **Isolated boundaries (issue #87, P2), a symmetric James step.** `initialize(..., isolated=.true.)`: each solve is
!< the zero-ghost solve `phi0 = A^-1 factor V f`; the screening charge `s = B^T phi0` on the ghost cells (extended by
!< zero, `phi0` satisfies the infinite operator with the extra source `-s` there); the ghost values of the free-space
!< solution `g = -G s`, `G` the inverse of the V-weighted operator on the infinite coarsest lattice from ghost cell to
!< ghost cell (the lattice Green's function over the coarsest spacing; `ghost_green`, a zero-padded FFT of the coarsest
!< box widened by the ghost layer, `gfft`); and the solve with them. The map `(A^-1 + A^-1 B G B^T A^-1) factor V` is
!< symmetric in the V-weighted form (A and G are), so self-gravity keeps its exact `E + W` exchange; on a uniform grid
!< it is the discrete free-space solution (the FFT solver's isolated mode), on AMR second order (G misses the refined
!< interior, a truncation-order perturbation). `G` is exact on the boundary only if every boundary leaf is at the
!< coarsest level (each face then has its own ghost cell of the coarse lattice): fatal otherwise (D-87-2). The textbook
!< James method (the screening charge on the boundary surface, the potential at the ghost centres) is not symmetric
!< (`selfgravity/james_proto.py`, the P0 prototype).
!<
!< **Scope (fatal otherwise).** The face kinds above: without `bc` every active direction periodic or none (a mixed
!< grid would get a silent zero on its bounded faces), with it the grid's periodicity and `bc` agree; isolated: 3-D, no
!< periodic direction, no `bc`, cubic cells, every boundary leaf at the coarsest level; 2:1 balance across faces;
!< fewer than 2^19 cells per direction at the finest level (the cell keys).
!< On output phi holds the leaf solution and its first ghost layer per face: same level, the neighbour; at a domain
!< face, the ghost of its kind, `g`, `2 g - phi_C` or `phi_C + h g` (g = 0 when none is given); across a 2:1
!< face, the value that makes the one-sided face difference equal the operator's own face gradient,
!< `(mean_f - phi_C) / (1.5 h)` (fine side: `phi_f - (mean_f - phi_C) / 1.5`; coarse side:
!< `phi_C + (mean_f - phi_C) / 0.75`); deeper ghost layers copy the first.
use :: adam_field_object, only : field_object
use :: adam_flail_parameters, only : FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, FLAIL_BC_DIRICHLET_GHOST, &
                                      FLAIL_BC_ZERO_GHOSTS, FLAIL_BC_FREE_GHOSTS
use :: adam_grid_object, only : grid_object
use :: adam_mpih_global, only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
use :: mpi
use :: penf, only : I4P, I8P, R8P, str

implicit none
private
public :: poisson_mg_object
public :: MG_SEAM_GHOST_SCALE

real(R8P),    parameter :: MG_SEAM_GHOST_SCALE=1.5_R8P !< 2:1 fine-side ghost: phi_g = phi_f - (mean_f - phi_C) / 1.5, the
                                                     !< composite operator's face gradient (clients that invert it, as
                                                     !< FLUME's seam gravity work, use this constant).
integer(I4P), parameter :: CLS_SAME=1_I4P    !< Neighbour position: a leaf of the same level.
integer(I4P), parameter :: CLS_FINER=2_I4P   !< Neighbour position: refined (finer leaves).
integer(I4P), parameter :: CLS_COARSER=3_I4P !< Neighbour position: inside a leaf one level coarser.
integer(I4P), parameter :: SPACE_COMP=-100_I4P !< Vector space: the composite leaf cells (levels are the other spaces;
                                              !< below lmin they may be negative, #84 E5c).
integer(I4P), parameter :: KEY_LEVEL_OFFSET=16_I4P !< Level offset in the cell keys (the levels below lmin, E5c).
integer(I4P), parameter :: REDUNDANT_ROOT_MAX=32768_I4P !< A root of at most this many cells is solved on every rank;
                                                        !< a larger one on rank 0 (E5c: four ranks transforming a
                                                        !< 128^3 root each contend for the node's memory).

type :: hash_table
   !< Open-addressing hash, int64 key (>= 0) -> int32 value (0: absent).
   integer(I8P), allocatable :: keys(:) !< Keys (-1: empty).
   integer(I4P), allocatable :: vals(:) !< Values.
   integer(I4P)              :: cap=0   !< Capacity (a prime).
endtype hash_table

type :: dcsr
   !< Distributed CSR operator: local rows; columns local (1..nloc) or halo (nloc+1..nloc+nhalo).
   integer(I4P)              :: nrow=0, nloc=0, nhalo=0 !< Rows, local columns, halo columns.
   integer(I4P), allocatable :: rowptr(:)               !< Row pointers [nrow+1].
   integer(I4P), allocatable :: col(:)                  !< Column slots.
   real(R8P),    allocatable :: val(:)                  !< Values.
   integer(I4P), allocatable :: scount(:), sdispl(:)    !< Send counts and displacements per rank.
   integer(I4P), allocatable :: sidx(:)                 !< Local column indices to send.
   integer(I4P), allocatable :: rcount(:), rdispl(:)    !< Receive counts and displacements per rank (halo order).
endtype dcsr

type :: rvec
   !< A real vector (one per level).
   real(R8P), allocatable :: v(:) !< Values.
endtype rvec

type :: triplets
   !< Assembly buffer: (row, column key, value).
   integer(I4P)              :: n=0  !< Entries.
   integer(I4P), allocatable :: r(:) !< Rows.
   integer(I8P), allocatable :: k(:) !< Column keys.
   real(R8P),    allocatable :: v(:) !< Values.
endtype triplets

type :: poisson_mg_object
   !< Multigrid-preconditioned CG Poisson solver on AMR grids.
   logical                   :: is_initialized=.false. !< Initialization status.
   integer(I4P)              :: nijk(3)=1_I4P          !< Cells per block per direction.
   integer(I4P)              :: ncb=1_I4P              !< Cells per block.
   integer(I4P)              :: ngc=0_I4P              !< Ghost cells per side.
   integer(I4P)              :: nd=3_I4P               !< Active directions.
   integer(I4P)              :: act(3)=[1,2,3]         !< The active directions.
   logical                   :: null(3)=.false.        !< Null directions.
   integer(I4P)              :: nblk0(3)=1_I4P         !< Blocks per direction at level 0.
   real(R8P)                 :: dx0(3)=0._R8P          !< Cell spacing at level 0.
   integer(I4P)              :: lmin=0_I4P, lmax=0_I4P !< Leaf levels range.
   integer(I4P)              :: lbot=0_I4P             !< The root level: lmin less the levels inside the blocks (E5c).
   integer(I4P)              :: root_max=0_I4P         !< Most cells of the root (E5c; the root solver's default).
   integer(I4P)              :: nr=1_I4P, me=0_I4P     !< Ranks, this rank.
   ! the global leaf list
   integer(I4P)              :: nleaf=0_I4P            !< Leaves over the ranks.
   integer(I4P), allocatable :: leaf_c(:,:)            !< Leaf block coordinates (i,j,k,level) [4,nleaf].
   integer(I4P), allocatable :: leaf_owner(:)          !< Owner rank.
   integer(I4P), allocatable :: leaf_loc(:)            !< Local block index (the field's) on its owner.
   type(hash_table)          :: leaf_hash              !< Block key -> leaf id.
   ! the nodes of the hierarchy (leaves and shadows)
   integer(I4P)              :: nnode=0_I4P            !< Nodes.
   integer(I4P), allocatable :: node_c(:,:)            !< Node block coordinates (i,j,k,level) [4,nnode].
   integer(I4P), allocatable :: node_owner(:)          !< Owner rank.
   integer(I4P), allocatable :: node_loc(:)            !< Local block index in its level space on its owner.
   logical,      allocatable :: node_leaf(:)           !< The node is a leaf.
   type(hash_table)          :: node_hash              !< Block key -> node id.
   integer(I4P), allocatable :: nloc_lev(:)            !< This rank's blocks per level [lbot:lmax].
   integer(I4P), allocatable :: lev_node(:,:)          !< This rank's level blocks: node id [max blocks, lbot:lmax]
                                                       !< (below lmin, the lmin blocks).
   integer(I4P)              :: nleaf_loc=0_I4P        !< This rank's leaves.
   integer(I4P), allocatable :: leaf_lev(:)            !< This rank's leaf b: its level.
   integer(I4P), allocatable :: leaf_lev_loc(:)        !< This rank's leaf b: its local index in its level space.
   ! operators
   type(dcsr)                :: a                      !< The composite operator on the leaves.
   type(dcsr)                :: g                      !< Ghost fill: the first ghost layer of the leaves.
   integer(I4P), allocatable :: g_dst(:,:)             !< Ghost fill destinations (i,j,k,b) [4,rows].
   type(dcsr),   allocatable :: al(:)                  !< Level operators [lbot+1:lmax].
   type(dcsr),   allocatable :: el(:)                  !< Level ghost terms from level l-1 [lbot+1:lmax].
   type(dcsr),   allocatable :: rl(:)                  !< Restriction l -> l-1 [lbot+1:lmax].
   type(dcsr),   allocatable :: pl(:)                  !< Prolongation l-1 -> l [lbot+1:lmax].
   type(rvec),   allocatable :: diag(:)                !< Level operator diagonals [lbot+1:lmax].
   type(rvec),   allocatable :: colour(:)              !< Red-black colour (0, 1) of the level cells [lbot+1:lmax].
   real(R8P),    allocatable :: vol(:)                 !< Leaf cell volumes.
   ! domain boundaries (issue #87; per face, issue #84 E1b)
   integer(I4P)              :: bc(2,3)=FLAIL_BC_PERIODIC !< Kind per face [lo/hi, axis] (adam_flail_parameters).
   logical                   :: periodic(3)=.true.     !< Periodic axes (a null axis counts as periodic).
   logical                   :: singular=.true.        !< No Dirichlet face: the operator's kernel is the constants.
   real(R8P)                 :: emin(3)=0._R8P         !< Domain lower corner.
   integer(I4P)              :: nbnd=0_I4P             !< This rank's boundary faces.
   integer(I4P), allocatable :: bnd_row(:)             !< Boundary face: the local leaf cell (row) it bounds.
   real(R8P),    allocatable :: bnd_w(:)               !< Boundary face: its weight in B (see the module).
   real(R8P),    allocatable :: bnd_x(:,:)             !< Boundary face: where its value lives, ghost or face centre.
   integer(I4P), allocatable :: bnd_kind(:)            !< Boundary face: its kind.
   integer(I4P), allocatable :: bnd_dir(:)             !< Boundary face: the outward axis, signed (-3..3).
   real(R8P),    allocatable :: bnd_gs(:)              !< Boundary face: its value's weight in the ghost fill.
   integer(I4P), allocatable :: g_bnd(:)               !< Ghost-fill row: its boundary face, 0 for an interior one.
   ! isolated boundaries, the James step (issue #87, P2)
   logical                   :: isolated=.false.       !< Free-space boundaries: the ghost values by the James step.
   integer(I4P)              :: nc(3)=1_I4P            !< Cells per direction at the coarsest level.
   integer(I4P), allocatable :: bnd_gi(:)              !< Boundary face: its ghost cell in the widened coarsest box.
   real(R8P),    allocatable :: u0(:)                  !< The last zero-ghost solution (warm start of the first solve).
   integer(I4P)              :: iterations0=0_I4P      !< Iterations of the last zero-ghost solve.
   class(poisson_fft_object), allocatable :: fft       !< The root solver, on every rank (FFTW, or the caller's backend).
   class(poisson_fft_object), allocatable :: gfft      !< Isolated: the ghost-to-ghost Green's function (padded FFT).
   ! controls and diagnostics
   real(R8P)                 :: tolerance=1.e-10_R8P   !< Relative max-norm residual to stop.
   integer(I4P)              :: max_iterations=200_I4P !< FCG iterations at most.
   integer(I4P)              :: sweeps=2_I4P           !< Red-black sweeps per level, down and up.
   integer(I4P)              :: iterations=0_I4P       !< Iterations of the last solve.
   real(R8P)                 :: residual=0._R8P        !< Relative residual of the last solve.
   real(R8P)                 :: tolerance_used=0._R8P  !< Stop tolerance of the last solve (the round-off floor of the
                                                       !< source when above `tolerance`).
   real(R8P)                 :: factor_mean=0._R8P     !< Mean convergence factor per iteration of the last solve.
   real(R8P),    allocatable :: u(:)                   !< The last solution (warm start).
   logical                   :: has_guess=.false.      !< u holds a solution of the current grid.
   contains
      ! public methods
      procedure, pass(self) :: destroy         !< Destroy the solver.
      procedure, pass(self) :: initialize      !< Gather the leaves, assemble the operators.
      procedure, pass(self) :: solve           !< Solve for phi from one variable of q.
      procedure, pass(self) :: solve_packed    !< Solve for phi from the packed interior source.
      procedure, pass(self) :: symmetry_defect !< |x . A y - y . A x| / (|x| |A y|) (verification).
      procedure, pass(self) :: apply_a        !< y = A x on the leaves (verification).
      procedure, pass(self) :: apply_b        !< y = B g: the boundary values onto the boundary rows.
      procedure, pass(self) :: apply_bt       !< s = B^T x: the boundary rows onto the boundary faces.
      ! private methods
      procedure, pass(self), private :: apply       !< y = M x for a distributed operator.
      procedure, pass(self), private :: assemble    !< Finish a distributed operator from triplets.
      procedure, pass(self), private :: a_row       !< The row of the composite operator at a cell.
      procedure, pass(self), private :: cell_volume !< Cell volume at a level.
      procedure, pass(self), private :: classify    !< Classify a neighbour position.
      procedure, pass(self), private :: comp_value  !< The level-l value of a cell as leaf-cell weights.
      procedure, pass(self), private :: dx_of       !< Cell spacing at a level.
      procedure, pass(self), private :: face_cells  !< The fine cells of a coarse face.
      procedure, pass(self), private :: face_weight !< V-weighted face weight (area over distance) at a level.
      procedure, pass(self), private :: ghost_green !< Isolated: g = -G s, the ghost-to-ghost Green's function.
      procedure, pass(self), private :: solve_fcg   !< The FCG solve with given ghost values (or none).
      procedure, pass(self), private :: lin         !< Linear index of a cell in a leaf block.
      procedure, pass(self), private :: lin_of      !< Linear index of a cell in a block at a level.
      procedure, pass(self), private :: ncb_of      !< Cells per block at a level.
      procedure, pass(self), private :: ncell       !< Cells along a direction over the domain at a level.
      procedure, pass(self), private :: nijk_of     !< Cells per block per direction at a level.
      procedure, pass(self), private :: outside     !< A cell position beyond a domain face.
      procedure, pass(self), private :: locate      !< Owner and local index of a cell key in a space.
      procedure, pass(self), private :: parent3     !< Parent cell coordinates.
      procedure, pass(self), private :: smooth      !< Red-black Gauss-Seidel on a level.
      procedure, pass(self), private :: vcycle      !< The preconditioner.
      procedure, pass(self), private :: wrap        !< Periodic wrap of cell coordinates at a level.
endtype poisson_mg_object

contains
   ! public methods
   subroutine destroy(self)
   !< Destroy the solver.
   class(poisson_mg_object), intent(inout) :: self !< The solver.

   if (allocated(self%fft)) then
      call self%fft%destroy
      deallocate(self%fft)
   endif
   if (allocated(self%gfft)) then
      call self%gfft%destroy
      deallocate(self%gfft)
   endif
   if (allocated(self%bnd_gi))       deallocate(self%bnd_gi)
   if (allocated(self%u0))           deallocate(self%u0)
   self%isolated = .false.
   self%iterations0 = 0
   if (allocated(self%leaf_c))       deallocate(self%leaf_c, self%leaf_owner, self%leaf_loc)
   if (allocated(self%node_c))       deallocate(self%node_c, self%node_owner, self%node_loc, self%node_leaf)
   if (allocated(self%nloc_lev))     deallocate(self%nloc_lev, self%lev_node)
   if (allocated(self%leaf_lev))     deallocate(self%leaf_lev, self%leaf_lev_loc)
   if (allocated(self%al))           deallocate(self%al, self%el, self%rl, self%pl, self%diag, self%colour)
   if (allocated(self%vol))          deallocate(self%vol)
   if (allocated(self%g_dst))        deallocate(self%g_dst)
   if (allocated(self%u))            deallocate(self%u)
   if (allocated(self%bnd_row))      deallocate(self%bnd_row, self%bnd_w, self%bnd_x, self%bnd_kind, self%bnd_dir, &
                                                self%bnd_gs)
   if (allocated(self%g_bnd))        deallocate(self%g_bnd)
   self%nbnd = 0
   self%bc = FLAIL_BC_PERIODIC
   self%periodic = .true.
   self%singular = .true.
   call hash_free(self%leaf_hash)
   call hash_free(self%node_hash)
   self%a = dcsr()
   self%g = dcsr()
   self%lbot = 0_I4P
   self%root_max = 0_I4P
   self%has_guess = .false.
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine initialize(self, grid, field, tolerance, max_iterations, coarse, isolated, bc, root_cells)
   !< Gather the leaves, derive the hierarchy, assemble the operator, the ghost fill and the preconditioner, set up the
   !< coarsest-level FFT solver: of the dynamic type of `coarse` when given (the FNL backend passes its cuFFT solver),
   !< the FFTW one otherwise. With `isolated` (issue #87), also the ghost-to-ghost Green's function of the James step.
   !< `bc` (issue #84, E1b): the kind of every face (see the module); a periodic axis of the grid must be periodic in
   !< `bc` and the other way round. Without it, every bounded axis takes ghost-Dirichlet faces. `root_cells` (#84 E5c):
   !< the levels inside the blocks below lmin stop once the root has at most this many cells (default: the root
   !< solver's `root_cells_max`, 32^3 on the host, 128^3 on cuFFT; tests force the levels with a small one).
   class(poisson_mg_object), intent(inout)        :: self           !< The solver.
   type(grid_object),        intent(in)           :: grid           !< The grid.
   type(field_object),       intent(in)           :: field          !< The field (block coordinates).
   real(R8P),                intent(in), optional :: tolerance      !< Relative residual to stop.
   integer(I4P),             intent(in), optional :: max_iterations !< FCG iterations at most.
   class(poisson_fft_object), intent(in), optional :: coarse        !< Mold of the coarsest-level solver.
   logical,                  intent(in), optional :: isolated       !< Free-space boundaries (James), else as the grid.
   integer(I4P),             intent(in), optional :: bc(2,3)        !< Kinds per face [lo/hi, axis].
   integer(I4P),             intent(in), optional :: root_cells     !< Most cells of the root (default: the backend's).
   type(poisson_fft_object)  :: host_fft             !< The host FFT solver (its root_cells_max, without coarse).
   integer(I4P), allocatable :: counts(:), displs(:) !< Allgather counts and displacements.
   integer(I4P), allocatable :: mine(:,:)            !< This rank's block coordinates and index [5,nb].
   integer(I4P), allocatable :: all_(:,:)            !< Every rank's [5,nleaf].
   integer(I4P), allocatable :: cnt(:,:)             !< Level blocks per rank and level.
   integer(I4P)              :: d, n, id, l, b, p, r !< Counters.
   integer(I4P)              :: c(4)                 !< Block coordinates.
   integer(I4P)              :: ierr                 !< MPI error status.

   call self%destroy
   self%root_max = host_fft%root_cells_max()
   if (present(coarse)) self%root_max = coarse%root_cells_max()
   if (present(root_cells)) self%root_max = root_cells
   if (present(tolerance)) self%tolerance = tolerance
   if (present(max_iterations)) self%max_iterations = max_iterations
   self%nr = mpih%procs_number
   self%me = mpih%myrank
   self%nijk = [grid%ni, grid%nj, grid%nk]
   self%ncb = product(self%nijk)
   self%ngc = grid%ngc
   self%null = grid%null_xyz
   self%nd = 0
   do d=1, 3
      if (self%null(d)) cycle
      self%nd = self%nd + 1
      self%act(self%nd) = d
   enddo
   ! the face kinds (issue #87; per face, issue #84 E1b): periodic as the grid, ghost-Dirichlet by default elsewhere
   if (present(isolated)) self%isolated = isolated
   if (self%isolated .and. present(bc)) call mpih%error_stop(msg='poisson_mg_object%initialize: isolated boundaries '//&
      'take no bc (the James step computes the ghost values)')
   do p=1, self%nd
      d = self%act(p)
      self%periodic(d) = grid%is_ijk_periodic(d)
      if (present(bc)) then
         if (self%periodic(d) .neqv. all(bc(:,d) == FLAIL_BC_PERIODIC)) call mpih%error_stop(msg= &
            'poisson_mg_object%initialize: direction '//trim(str(d, .true.))//' is '// &
            merge('    periodic', 'not periodic', self%periodic(d))//' in the grid and '// &
            merge('    periodic', 'not periodic', all(bc(:,d) == FLAIL_BC_PERIODIC))//' in bc')
         if (count(bc(:,d) == FLAIL_BC_DIRICHLET_GHOST .or. bc(:,d) == FLAIL_BC_ZERO_GHOSTS) == 1) &
            call mpih%error_stop(msg='poisson_mg_object%initialize: direction '//trim(str(d, .true.))// &
                                     ' mixes a ghost-Dirichlet face (or a zero-ghost one) with a face kind')
         if (any(bc(:,d) == FLAIL_BC_FREE_GHOSTS)) call mpih%error_stop(msg='poisson_mg_object%initialize: '// &
            'free-space ghosts on direction '//trim(str(d, .true.))//' need the FFT solver with the wide D(G) (issue '// &
            '#84 E3c)')
         do r=1, 2
            if (all(bc(r,d) /= [FLAIL_BC_PERIODIC, FLAIL_BC_DIRICHLET, FLAIL_BC_NEUMANN, FLAIL_BC_DIRICHLET_GHOST, &
                                FLAIL_BC_ZERO_GHOSTS])) &
               call mpih%error_stop(msg='poisson_mg_object%initialize: unknown boundary kind '// &
                                        trim(str(bc(r,d), .true.))//' on direction '//trim(str(d, .true.)))
         enddo
         ! zero ghosts (PRISM's Dirichlet, #84 E2b): the compact operator reads one ghost layer, ghost-Dirichlet at 0
         self%bc(:,d) = merge(FLAIL_BC_DIRICHLET_GHOST, bc(:,d), bc(:,d) == FLAIL_BC_ZERO_GHOSTS)
      elseif (.not.self%periodic(d)) then
         self%bc(:,d) = FLAIL_BC_DIRICHLET_GHOST
      endif
   enddo
   ! without bc, every active direction periodic or none: a mixed grid would get a silent zero on its bounded faces
   if (.not.present(bc) .and. any(self%periodic(self%act(1:self%nd))) .and. &
       .not.all(self%periodic(self%act(1:self%nd)))) call mpih%error_stop(msg='poisson_mg_object%initialize: mixed '// &
      'periodic and non-periodic directions without bc; give the kind of every face (issue #84, E1b)')
   self%singular = .not.any(self%bc == FLAIL_BC_DIRICHLET .or. self%bc == FLAIL_BC_DIRICHLET_GHOST)
   if (self%isolated) then
      if (any(self%periodic(self%act(1:self%nd)))) call mpih%error_stop(msg='poisson_mg_object%initialize: '// &
         'isolated boundaries on a periodic grid')
      if (self%nd /= 3) call mpih%error_stop(msg='poisson_mg_object%initialize: isolated boundaries are 3-D only '// &
         '(the kernel is the 3-D lattice Green function)')
   endif
   self%emin = grid%domain_emin
   self%nblk0 = nint((grid%domain_emax - grid%domain_emin) / grid%block_dxyz(:,0))
   self%dx0 = grid%cell_dxyz(:,0)
   ! the global leaf list
   allocate(counts(self%nr), displs(self%nr))
   n = field%blocks_number
   call MPI_ALLGATHER(n, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   displs(1) = 0
   do r=2, self%nr
      displs(r) = displs(r-1) + counts(r-1)
   enddo
   self%nleaf = sum(counts)
   allocate(mine(5,max(n,1)), all_(5,self%nleaf))
   do b=1, n
      mine(1:4,b) = field%coordinates(1:4,b)
      mine(5,b) = b
   enddo
   call MPI_ALLGATHERV(mine, 5*n, MPI_INTEGER, all_, 5*counts, 5*displs, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   allocate(self%leaf_c(4,self%nleaf), self%leaf_owner(self%nleaf), self%leaf_loc(self%nleaf))
   call hash_init(self%leaf_hash, self%nleaf)
   do r=1, self%nr
      do b=1, counts(r)
         id = displs(r) + b
         self%leaf_c(:,id) = all_(1:4,id)
         self%leaf_owner(id) = r - 1
         self%leaf_loc(id) = all_(5,id)
         call hash_put(self%leaf_hash, bkey(self%leaf_c(:,id)), id)
      enddo
   enddo
   self%lmin = minval(self%leaf_c(4,:))
   self%lmax = maxval(self%leaf_c(4,:))
   if (2_I8P**self%lmax * maxval(self%nblk0 * self%nijk) >= 2_I8P**19) call mpih%error_stop(msg= &
      'poisson_mg_object%initialize: the finest level has 2^19 or more cells along a direction (the key limit)')
   ! the levels inside the blocks below lmin (#84 E5c): halve the cells per block along every active direction while
   ! they are even and the root has more than root_max cells
   self%lbot = self%lmin
   do
      c(1:3) = self%nijk_of(self%lbot)
      if (self%lmin - self%lbot >= KEY_LEVEL_OFFSET) exit
      if (product(int(self%nblk0, I8P) * merge(1_I8P, 2_I8P**self%lmin, self%null) * int(c(1:3), I8P)) <= &
          int(self%root_max, I8P)) exit
      if (any([(mod(c(self%act(p)), 2) /= 0, p=1, self%nd)])) exit
      self%lbot = self%lbot - 1
   enddo
   ! the nodes: every leaf and its ancestors down to lmin; a shadow's owner is its first leaf's
   allocate(self%node_c(4,2*self%nleaf+8), self%node_owner(2*self%nleaf+8), self%node_loc(2*self%nleaf+8), &
            self%node_leaf(2*self%nleaf+8))
   call hash_init(self%node_hash, 2*self%nleaf+8)
   self%nnode = 0
   do id=1, self%nleaf
      call add_node(self%leaf_c(:,id), self%leaf_owner(id), .true.)
   enddo
   do id=1, self%nleaf
      c = self%leaf_c(:,id)
      do l=c(4)-1, self%lmin, -1
         do p=1, self%nd
            c(self%act(p)) = c(self%act(p)) / 2
         enddo
         c(4) = l
         if (hash_get(self%node_hash, bkey(c)) > 0) exit
         call add_node(c, self%leaf_owner(id), .false.)
      enddo
   enddo
   ! level spaces: local numbering per owner and level, in node order (the same on every rank)
   allocate(cnt(0:self%nr-1,self%lmin:self%lmax))
   cnt = 0
   do id=1, self%nnode
      l = self%node_c(4,id)
      cnt(self%node_owner(id),l) = cnt(self%node_owner(id),l) + 1
      self%node_loc(id) = cnt(self%node_owner(id),l)
   enddo
   allocate(self%nloc_lev(self%lbot:self%lmax), self%lev_node(max(maxval(cnt(self%me,:)),1),self%lbot:self%lmax))
   self%nloc_lev(self%lmin:self%lmax) = cnt(self%me,:)
   do id=1, self%nnode
      if (self%node_owner(id) == self%me) self%lev_node(self%node_loc(id),self%node_c(4,id)) = id
   enddo
   do l=self%lbot, self%lmin-1 ! the levels inside the blocks: the lmin blocks, the same owners
      self%nloc_lev(l) = self%nloc_lev(self%lmin)
      self%lev_node(:,l) = self%lev_node(:,self%lmin)
   enddo
   self%nleaf_loc = n
   allocate(self%leaf_lev(max(n,1)), self%leaf_lev_loc(max(n,1)))
   do b=1, n
      id = hash_get(self%node_hash, bkey(field%coordinates(1:4,b)))
      self%leaf_lev(b) = field%coordinates(4,b)
      self%leaf_lev_loc(b) = self%node_loc(id)
   enddo
   allocate(self%vol(max(n*self%ncb,1)))
   self%vol = 1._R8P
   do b=1, n
      self%vol((b-1)*self%ncb+1:b*self%ncb) = self%cell_volume(field%coordinates(4,b))
   enddo
   if (self%isolated) then ! the unit lattice's Green's function needs cubic cells
      if (any(abs(self%dx0 - self%dx0(1)) > 1.e-12_R8P * self%dx0(1))) &
         call mpih%error_stop(msg='poisson_mg_object%initialize: isolated boundaries need cubic cells')
   endif
   self%nc = self%nblk0 * merge(1, 2**self%lmin, self%null) * self%nijk
   call assemble_a
   call assemble_ghost_fill
   call assemble_levels
   ! the root: every block at lmin, the whole domain, nijk_of(lbot) cells per block; solved on every rank when small
   ! (#84 E5c), else on rank 0 as before E5c
   if (present(coarse)) then
      allocate(self%fft, mold=coarse)
   else
      allocate(poisson_fft_object :: self%fft)
   endif
   call self%fft%initialize_blocks(nblk=self%nblk0 * merge(1, 2**self%lmin, self%null), nijk=self%nijk_of(self%lbot), &
                                   ngc=0_I4P, dxyz=self%dx_of(self%lbot), level=self%lmin, null_xyz=self%null,        &
                                   isolated=.false., bc=self%bc, redundant=product(int(self%nblk0, I8P) *          &
                                   merge(1_I8P, 2_I8P**self%lmin, self%null) * int(self%nijk_of(self%lbot), I8P)) <= &
                                   int(REDUNDANT_ROOT_MAX, I8P))
   if (self%isolated) then ! the coarsest box widened by the ghost layer, zero-padded: sources and values on its shell
      if (present(coarse)) then
         allocate(self%gfft, mold=coarse)
      else
         allocate(poisson_fft_object :: self%gfft)
      endif
      call self%gfft%initialize_blocks(nblk=[1_I4P, 1_I4P, 1_I4P], nijk=self%nc + 2, ngc=0_I4P, &
                                       dxyz=self%dx_of(self%lmin), level=self%lmin, null_xyz=self%null, isolated=.true.)
   endif
   allocate(self%u(max(n*self%ncb,1)), self%u0(max(n*self%ncb,1)))
   self%u = 0._R8P
   self%u0 = 0._R8P
   self%has_guess = .false.
   self%is_initialized = .true.
   contains
      subroutine add_node(cc, owner, leaf)
      !< Append a node.
      integer(I4P), intent(in) :: cc(4) !< Block coordinates.
      integer(I4P), intent(in) :: owner !< Owner rank.
      logical,      intent(in) :: leaf  !< A leaf.

      self%nnode = self%nnode + 1
      if (self%nnode > size(self%node_owner)) call mpih%error_stop(msg='poisson_mg_object: node storage overflow')
      self%node_c(:,self%nnode) = cc
      self%node_owner(self%nnode) = owner
      self%node_leaf(self%nnode) = leaf
      call hash_put(self%node_hash, bkey(cc), self%nnode)
      endsubroutine add_node

      subroutine assemble_a
      !< The composite operator on this rank's leaf cells, row by row (`a_row`). Not symmetrised: M is symmetric by
      !< construction, so `symmetry_defect` verifies the rows, not an assembly that would force it.
      type(triplets)            :: t        !< Assembly buffer.
      integer(I8P), allocatable :: keys(:)  !< A row: column keys.
      real(R8P),    allocatable :: vals(:)  !< A row: values.
      integer(I4P)              :: i, j, kk, e, lv !< Counters, level.
      integer(I4P)              :: cc(3)    !< Cell coordinates.

      do b=1, field%blocks_number
         lv = field%coordinates(4,b)
         do kk=1, self%nijk(3)
            do j=1, self%nijk(2)
               do i=1, self%nijk(1)
                  cc = field%coordinates(1:3,b) * self%nijk + [i-1, j-1, kk-1]
                  call self%a_row(lv, cc, keys, vals)
                  do e=1, size(keys)
                     call trip_add(t, (b - 1) * self%ncb + self%lin([i, j, kk]), keys(e), vals(e))
                  enddo
               enddo
            enddo
         enddo
      enddo
      call self%assemble(t, field%blocks_number * self%ncb, SPACE_COMP, self%a)
      endsubroutine assemble_a

      subroutine assemble_ghost_fill
      !< The first ghost layer of each face of every local leaf block, as leaf-cell weights.
      type(triplets)            :: t        !< Assembly buffer.
      integer(I8P), allocatable :: keys(:)  !< Value keys.
      real(R8P),    allocatable :: wts(:)   !< Value weights.
      integer(I4P), allocatable :: dst(:,:) !< Destinations.
      integer(I4P)              :: nrow, pp, dd, sgn, i, j, kk, e, a_, lv, cls !< Counters.
      integer(I4P)              :: cc(3), nb3(3), cp(3), ijk(3), f3(3) !< Coordinates.
      integer(I8P), allocatable :: keys_f(:) !< The fine cells of a coarse face.
      real(R8P)                 :: hxyz(3)   !< Cell spacing at the level.

      allocate(dst(4,max(2*self%nd*self%ncb*field%blocks_number,1)))
      allocate(self%g_bnd(max(2*self%nd*self%ncb*field%blocks_number,1)))
      allocate(self%bnd_row(max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_w(max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_x(3,max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_kind(max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_dir(max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_gs(max(2*self%nd*self%ncb*field%blocks_number,1)), &
               self%bnd_gi(max(2*self%nd*self%ncb*field%blocks_number,1)))
      self%g_bnd = 0
      self%bnd_gi = 0
      self%nbnd = 0
      nrow = 0
      do b=1, field%blocks_number
         lv = field%coordinates(4,b)
         do pp=1, self%nd
            dd = self%act(pp)
            do sgn=-1, 1, 2
               do kk=1, self%nijk(3)
                  do j=1, self%nijk(2)
                     do i=1, self%nijk(1)
                        ijk = [i, j, kk]
                        if (sgn < 0 .and. ijk(dd) /= 1) cycle
                        if (sgn > 0 .and. ijk(dd) /= self%nijk(dd)) cycle
                        cc = field%coordinates(1:3,b) * self%nijk + ijk - 1
                        nb3 = cc ; nb3(dd) = nb3(dd) + sgn
                        nb3 = self%wrap(lv, nb3)
                        if (allocated(keys)) deallocate(keys, wts)
                        allocate(keys(0), wts(0))
                        cls = 0
                        if (self%outside(lv, nb3)) then ! a domain face: its kind's ghost, the value set by solve
                           if (self%isolated .and. lv /= self%lmin) call mpih%error_stop(msg= &
                              'poisson_mg_object%initialize: a leaf at level '//trim(str(lv, .true.))//          &
                              ' touches the domain boundary; isolated boundaries need every leaf on the boundary '// &
                              'at the coarsest level, '//trim(str(self%lmin, .true.))//' (issue #87, D-87-2)')
                           self%nbnd = self%nbnd + 1
                           self%bnd_row(self%nbnd) = (b - 1) * self%ncb + self%lin(ijk)
                           self%bnd_kind(self%nbnd) = self%bc(merge(1, 2, sgn < 0),dd)
                           hxyz = self%dx_of(lv)
                           self%bnd_dir(self%nbnd) = sgn * dd
                           self%bnd_x(:,self%nbnd) = self%emin + (real(nb3, R8P) + 0.5_R8P) * self%dx_of(lv)
                           do a_=1, 3
                              if (self%null(a_)) self%bnd_x(a_,self%nbnd) = field%emin(a_,b) + &
                                 0.5_R8P * (field%emax(a_,b) - field%emin(a_,b))
                           enddo
                           select case(self%bnd_kind(self%nbnd))
                           case(FLAIL_BC_DIRICHLET) ! phi_ghost = 2 g - phi_C, g at the face centre
                              self%bnd_w(self%nbnd) = 2._R8P * self%face_weight(lv, dd)
                              self%bnd_gs(self%nbnd) = 2._R8P
                              self%bnd_x(dd,self%nbnd) = self%bnd_x(dd,self%nbnd) - sgn * 0.5_R8P * hxyz(dd)
                              keys = [ckey(lv, cc)] ; wts = [-1._R8P]
                           case(FLAIL_BC_NEUMANN) ! phi_ghost = phi_C + h g, g the outward derivative at the face
                              self%bnd_w(self%nbnd) = self%face_weight(lv, dd) * hxyz(dd)
                              self%bnd_gs(self%nbnd) = hxyz(dd)
                              self%bnd_x(dd,self%nbnd) = self%bnd_x(dd,self%nbnd) - sgn * 0.5_R8P * hxyz(dd)
                              keys = [ckey(lv, cc)] ; wts = [1._R8P]
                           case default ! ghost-Dirichlet: phi_ghost = g, g at the ghost centre
                              self%bnd_w(self%nbnd) = self%face_weight(lv, dd)
                              self%bnd_gs(self%nbnd) = 1._R8P
                           endselect
                           self%g_bnd(nrow+1) = self%nbnd
                           ! its ghost cell in the coarsest box widened by one cell per side, i fastest
                           self%bnd_gi(self%nbnd) = (nb3(1) + 2) + (self%nc(1) + 2) * ((nb3(2) + 1) + &
                                                    (self%nc(2) + 2) * (nb3(3) + 1))
                        else
                           cls = self%classify(lv, nb3)
                        endif
                        select case(cls)
                        case(CLS_SAME)
                           keys = [ckey(lv, nb3)] ; wts = [1._R8P]
                        case(CLS_FINER) ! coarse side: phi_C + (mean_f - phi_C) / 0.75, M's face gradient over H
                           f3 = 2 * cc
                           do a_=1, self%nd
                              f3(self%act(a_)) = 2 * cc(self%act(a_))
                           enddo
                           f3(dd) = 2 * nb3(dd) + merge(0, 1, sgn > 0)
                           f3 = self%wrap(lv+1, f3)
                           call self%face_cells(lv+1, f3, dd, keys_f)
                           keys = [ckey(lv, cc)] ; wts = [1._R8P - 1._R8P / 0.75_R8P]
                           do e=1, size(keys_f)
                              keys = [keys, keys_f(e)] ; wts = [wts, 1._R8P / 0.75_R8P / size(keys_f)]
                           enddo
                        case(CLS_COARSER) ! fine side: phi_f - (mean_f - phi_C) / 1.5, M's face gradient over h
                           cp = self%parent3(nb3)
                           call self%face_cells(lv, cc, dd, keys_f)
                           keys = [ckey(lv, cc), ckey(lv-1, cp)] ; wts = [1._R8P, 1._R8P / MG_SEAM_GHOST_SCALE]
                           do e=1, size(keys_f)
                              keys = [keys, keys_f(e)] ; wts = [wts, -1._R8P / MG_SEAM_GHOST_SCALE / size(keys_f)]
                           enddo
                        endselect
                        nrow = nrow + 1
                        ijk(dd) = ijk(dd) + sgn
                        dst(:,nrow) = [ijk, b]
                        do e=1, size(keys)
                           call trip_add(t, nrow, keys(e), wts(e))
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      allocate(self%g_dst(4,max(nrow,1)))
      if (nrow > 0) self%g_dst(:,1:nrow) = dst(:,1:nrow)
      call self%assemble(t, nrow, SPACE_COMP, self%g)
      endsubroutine assemble_ghost_fill

      subroutine assemble_levels
      !< The preconditioner's level operators, ghost terms and transfers, for every level above the root (lbot): the
      !< tree levels above lmin, lmin itself and the levels inside the blocks when there are some (#84 E5c).
      integer(I4P) :: lv, lt !< Level, top.

      lt = max(self%lmax, self%lbot+1)
      allocate(self%al(self%lbot+1:lt), self%el(self%lbot+1:lt), self%rl(self%lbot+1:lt), self%pl(self%lbot+1:lt), &
               self%diag(self%lbot+1:lt), self%colour(self%lbot+1:lt))
      do lv=self%lbot+1, self%lmax
         call level_ops(lv)
      enddo
      endsubroutine assemble_levels

      subroutine level_ops(lv)
      !< Level lv: A_l (the plain Laplacian on the level's cells, a face to a cell outside the level reading E_l: the
      !< level-(lv-1) cell containing it), R_l (rows: this rank's shadow cells at lv-1, the sum of their children at
      !< lv) and P_l (rows: the level's cells, their parent at lv-1). At lmin and below every cell of the domain is in
      !< the level (no E_l), and below lmin every block has its cells at lv-1 (E5c: the blocks coarsen inside).
      integer(I4P), intent(in) :: lv !< Level.
      type(triplets)           :: ta, te, tr, tp !< Assembly buffers.
      integer(I4P)             :: q, nloc, i, j, kk, pp, dd, sgn, row, o, id_, ch !< Counters.
      integer(I4P)             :: cc(3), nb3(3), pc(3), cch(3) !< Coordinates.
      integer(I4P)             :: nl(3), nlc(3) !< Cells per block at lv and lv-1.
      real(R8P)                :: w !< Face weight.

      nl = self%nijk_of(lv)
      nlc = self%nijk_of(lv-1)
      nloc = self%nloc_lev(lv) * self%ncb_of(lv)
      allocate(self%colour(lv)%v(max(nloc,1)), self%diag(lv)%v(max(nloc,1)))
      self%diag(lv)%v = 1._R8P
      if (nloc > 0) self%diag(lv)%v(1:nloc) = 0._R8P
      do q=1, self%nloc_lev(lv)
         id_ = self%lev_node(q,lv)
         do kk=1, nl(3)
            do j=1, nl(2)
               do i=1, nl(1)
                  cc = self%node_c(1:3,id_) * nl + [i-1, j-1, kk-1]
                  row = (q - 1) * self%ncb_of(lv) + self%lin_of(lv, [i, j, kk])
                  self%colour(lv)%v(row) = real(mod(sum(cc), 2), R8P)
                  do pp=1, self%nd
                     dd = self%act(pp)
                     w = self%face_weight(lv, dd)
                     do sgn=-1, 1, 2
                        nb3 = cc ; nb3(dd) = nb3(dd) + sgn
                        nb3 = self%wrap(lv, nb3)
                        if (self%outside(lv, nb3)) then ! a domain face, homogeneous: the diagonal of its kind
                           select case(self%bc(merge(1, 2, sgn < 0),dd))
                           case(FLAIL_BC_DIRICHLET)
                              call trip_add(ta, row, ckey(lv, cc), -2._R8P * w)
                              self%diag(lv)%v(row) = self%diag(lv)%v(row) - 2._R8P * w
                           case(FLAIL_BC_DIRICHLET_GHOST)
                              call trip_add(ta, row, ckey(lv, cc), -w)
                              self%diag(lv)%v(row) = self%diag(lv)%v(row) - w
                           endselect ! Neumann: no flux
                           cycle
                        endif
                        call trip_add(ta, row, ckey(lv, cc), -w)
                        self%diag(lv)%v(row) = self%diag(lv)%v(row) - w
                        if (lv <= self%lmin) then ! every cell of the domain is in the level
                           call trip_add(ta, row, ckey(lv, nb3), w)
                        elseif (hash_get(self%node_hash, bkey([nb3 / self%nijk, lv])) > 0) then
                           call trip_add(ta, row, ckey(lv, nb3), w)
                        else
                           call trip_add(te, row, ckey(lv-1, self%parent3(nb3)), w)
                        endif
                     enddo
                  enddo
                  call trip_add(tp, row, ckey(lv-1, self%parent3(cc)), 1._R8P)
               enddo
            enddo
         enddo
      enddo
      do q=1, self%nloc_lev(lv-1)
         id_ = self%lev_node(q,lv-1)
         if (lv-1 >= self%lmin .and. self%node_leaf(id_)) cycle ! below lmin every block has its finer cells
         do kk=1, nlc(3)
            do j=1, nlc(2)
               do i=1, nlc(1)
                  pc = self%node_c(1:3,id_) * nlc + [i-1, j-1, kk-1]
                  row = (q - 1) * self%ncb_of(lv-1) + self%lin_of(lv-1, [i, j, kk])
                  do o=0, 2**self%nd-1
                     cch = pc
                     do ch=1, self%nd
                        cch(self%act(ch)) = 2 * pc(self%act(ch)) + ibits(o, ch-1, 1)
                     enddo
                     call trip_add(tr, row, ckey(lv, cch), 1._R8P)
                  enddo
               enddo
            enddo
         enddo
      enddo
      call self%assemble(ta, nloc, lv, self%al(lv))
      call self%assemble(te, nloc, lv-1, self%el(lv))
      call self%assemble(tr, self%nloc_lev(lv-1) * self%ncb_of(lv-1), lv, self%rl(lv))
      call self%assemble(tp, nloc, lv-1, self%pl(lv))
      endsubroutine level_ops

   endsubroutine initialize

   subroutine solve(self, field, q, iv, factor, phi, energy, ghost)
   !< Solve `A phi = factor V q(iv) - B ghost` on the leaves, the mean of the right-hand side removed when the operator
   !< is singular (periodic, or no Dirichlet face); phi's ghosts filled (see the module).
   class(poisson_mg_object), intent(inout)         :: self                     !< The solver.
   type(field_object),       intent(in)            :: field                    !< The field (coordinates).
   real(R8P),    intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Conservative [nv,i,j,k,nb].
   integer(I4P), intent(in)    :: iv                                           !< Variable of q, the source.
   real(R8P),    intent(in)    :: factor                                       !< Source factor (4 pi G).
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Potential [i,j,k,nb].
   real(R8P),    intent(out), optional :: energy                               !< sum_cells V f phi (all ranks).
   real(R8P),    intent(in),  optional :: ghost(1:)                            !< The boundary value per face.
   real(R8P),    allocatable                       :: buf(:)                   !< Packed interior source.
   integer(I4P)                                    :: b, i, j, k, c            !< Counters.

   allocate(buf(max(field%blocks_number*self%ncb,1)))
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
   call self%solve_packed(field=field, buf=buf, factor=factor, phi=phi, energy=energy, ghost=ghost)
   endsubroutine solve

   subroutine solve_packed(self, field, buf, factor, phi, energy, ghost)
   !< Solve from the packed interior source of this rank's leaves (block-major, `i` fastest); then the ghost fill.
   !< Isolated (issue #87, P2), the James step: the zero-ghost solve `phi0 = A^-1 factor V f`, the screening charge
   !< `s = B^T phi0` on the ghost cells, the ghost values of the free-space solution `g = -G s` (`ghost_green`), and
   !< the solve with them, `phi = A^-1 (factor V f - B g)`: the map is `(A^-1 + A^-1 B G B^T A^-1) factor V`,
   !< symmetric in the V-weighted form since A and G are. Each of the two solves warm-starts from its own last solution.
   class(poisson_mg_object), intent(inout)         :: self    !< The solver.
   type(field_object),       intent(in)            :: field   !< The field.
   real(R8P),                intent(in)            :: buf(1:) !< Packed interior source.
   real(R8P),                intent(in)            :: factor  !< Source factor.
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Potential [i,j,k,nb].
   real(R8P),                intent(out), optional :: energy  !< sum_cells V f phi (all ranks).
   real(R8P),                intent(in),  optional :: ghost(1:) !< The boundary value per boundary face.
   real(R8P), allocatable :: g(:), s(:), u_(:) !< Ghost values, screening charge, the last solution aside.

   if (.not.self%isolated) then
      call self%solve_fcg(field=field, buf=buf, factor=factor, phi=phi, energy=energy, ghost=ghost)
      return
   endif
   if (present(ghost)) call mpih%error_stop(msg='poisson_mg_object%solve: ghost values given with isolated '// &
      'boundaries (the James step computes them)')
   allocate(g(max(self%nbnd,1)), s(max(self%nbnd,1)))
   g = 0._R8P ; s = 0._R8P
   call move_alloc(self%u, u_) ; call move_alloc(self%u0, self%u)
   call self%solve_fcg(field=field, buf=buf, factor=factor, phi=phi, ghost=g)
   self%iterations0 = self%iterations
   call self%apply_bt(self%u, s)
   call move_alloc(self%u, self%u0) ; call move_alloc(u_, self%u)
   call self%ghost_green(s, g)
   call self%solve_fcg(field=field, buf=buf, factor=factor, phi=phi, energy=energy, ghost=g)
   endsubroutine solve_packed

   subroutine solve_fcg(self, field, buf, factor, phi, energy, ghost)
   !< Solve from the packed interior source of this rank's leaves (block-major, `i` fastest): FCG on the operator, preconditioned
   !< by the V-cycle, from the previous solution when there is one; then the ghost fill.
   class(poisson_mg_object), intent(inout)         :: self    !< The solver.
   type(field_object),       intent(in)            :: field   !< The field.
   real(R8P),                intent(in)            :: buf(1:) !< Packed interior source.
   real(R8P),                intent(in)            :: factor  !< Source factor.
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Potential [i,j,k,nb].
   real(R8P),                intent(out), optional :: energy  !< sum_cells V f phi (all ranks).
   real(R8P),                intent(in),  optional :: ghost(1:) !< The boundary value per boundary face.
   real(R8P), allocatable :: f(:), r(:), z(:), p(:), qv(:), pp(:), qp(:), gv(:) !< Vectors.
   real(R8P)              :: fmax, alpha, beta, vsum, mean, r0, e_ !< Scalars.
   real(R8P)              :: fraw(1), tol !< Max |factor f| before the mean is removed, the stop tolerance.
   integer(I4P)           :: n, it, b, i, j, k, c, row, l, nn !< Counters.
   logical                :: has_prev                          !< A previous direction exists.
   integer(I4P)           :: ierr                              !< MPI error status.

   if (.not.self%is_initialized) call mpih%error_stop(msg='poisson_mg_object%solve: called before initialize')
   n = field%blocks_number * self%ncb
   nn = max(n,1)
   allocate(f(nn), r(nn), z(nn), p(nn), qv(nn), pp(nn), qp(nn))
   f = 0._R8P ; r = 0._R8P ; z = 0._R8P ; p = 0._R8P ; qv = 0._R8P ; pp = 0._R8P ; qp = 0._R8P
   if (n > 0) f(1:n) = factor * self%vol(1:n) * buf(1:n)
   if (present(ghost) .and. all(self%periodic)) call mpih%error_stop(msg='poisson_mg_object%solve: ghost values '// &
      'given on a periodic grid')
   if (present(ghost)) then
      if (size(ghost) < self%nbnd) call mpih%error_stop(msg='poisson_mg_object%solve: fewer ghost values than faces')
   endif
   ! the boundary values to the right-hand side, f - B g; singular (no Dirichlet face): zero mean, f - V sum(f) / sum(V)
   vsum = gdot(self%vol, [(1._R8P, i=1, nn)], n)
   if (present(ghost)) then
      do i=1, self%nbnd
         f(self%bnd_row(i)) = f(self%bnd_row(i)) - self%bnd_w(i) * ghost(i)
      enddo
   endif
   if (self%singular) then
      mean = gdot(f, [(1._R8P, i=1, nn)], n) / vsum
      if (n > 0) f(1:n) = f(1:n) - mean * self%vol(1:n)
   endif
   fmax = 0._R8P
   do i=1, n
      fmax = max(fmax, abs(f(i) / self%vol(i)))
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, fmax, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   ! the stop: the tolerance, or the round-off floor of the source when its mean dominates (a 1e-6 perturbation of a
   ! uniform density is known to 1e-16 / 1e-6 = 1e-10 of itself, issue #81 P5b: the iteration would never stop); the
   ! recursive residual reaches about 2 epsilon max|factor f| / max|f - mean| (measured 1.7), so the floor is 100 of it
   fraw = 0._R8P
   if (n > 0) fraw(1) = abs(factor) * maxval(abs(buf(1:n)))
   call MPI_ALLREDUCE(MPI_IN_PLACE, fraw, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
   tol = self%tolerance
   if (fmax > 0._R8P) tol = max(tol, 100._R8P * epsilon(1._R8P) * fraw(1) / fmax)
   self%tolerance_used = tol
   if (.not.self%has_guess) self%u = 0._R8P
   call self%apply(self%a, self%u, qv)
   if (n > 0) r(1:n) = f(1:n) - qv(1:n)
   r0 = rel_max(r)
   self%iterations = 0
   self%residual = r0
   has_prev = .false.
   do it=1, self%max_iterations
      if (self%residual <= tol .or. fmax == 0._R8P) exit
      call self%vcycle(r, z)
      if (self%singular) then
         mean = gdot(z, self%vol, n) / vsum
         if (n > 0) z(1:n) = z(1:n) - mean
      endif
      p = z
      if (has_prev) then
         beta = gdot(z, qp, n) / gdot(pp, qp, n)
         if (n > 0) p(1:n) = z(1:n) - beta * pp(1:n)
      endif
      call self%apply(self%a, p, qv)
      alpha = gdot(p, r, n) / gdot(p, qv, n)
      if (n > 0) then
         self%u(1:n) = self%u(1:n) + alpha * p(1:n)
         r(1:n) = r(1:n) - alpha * qv(1:n)
      endif
      pp = p ; qp = qv ; has_prev = .true.
      self%iterations = it
      self%residual = rel_max(r)
   enddo
   if (self%residual > tol .and. fmax > 0._R8P .and. self%me == 0) &
      print '(A)', mpih%myrankstr//'warning: poisson_mg_object: '//trim(str(self%max_iterations, .true.))// &
                   ' FCG iterations left the residual at '//trim(str(self%residual))//' (tolerance '//trim(str(tol))//')'
   self%factor_mean = 0._R8P
   if (self%iterations > 0 .and. r0 > 0._R8P .and. self%residual > 0._R8P) &
      self%factor_mean = (self%residual / r0)**(1._R8P / self%iterations)
   ! remove the solution's mean, keep it for the next solve, fill phi
   if (self%singular) then
      mean = gdot(self%u, self%vol, n) / vsum
      if (n > 0) self%u(1:n) = self%u(1:n) - mean
   endif
   self%has_guess = .true.
   c = 0
   do b=1, field%blocks_number
      do k=1, self%nijk(3)
         do j=1, self%nijk(2)
            do i=1, self%nijk(1)
               c = c + 1
               phi(i,j,k,b) = self%u(c)
            enddo
         enddo
      enddo
   enddo
   allocate(gv(max(self%g%nrow,1)))
   call self%apply(self%g, self%u, gv)
   do row=1, self%g%nrow
      b = self%g_dst(4,row)
      i = self%g_dst(1,row) ; j = self%g_dst(2,row) ; k = self%g_dst(3,row)
      if (self%g_bnd(row) > 0) then ! a domain face: the ghost of its kind (the boundary value 0 without one)
         if (self%bnd_kind(self%g_bnd(row)) == FLAIL_BC_DIRICHLET_GHOST) then
            gv(row) = 0._R8P
            if (present(ghost)) gv(row) = ghost(self%g_bnd(row))
         elseif (present(ghost)) then ! gv holds -phi_C (Dirichlet) or phi_C (Neumann)
            gv(row) = gv(row) + self%bnd_gs(self%g_bnd(row)) * ghost(self%g_bnd(row))
         endif
      endif
      phi(i,j,k,b) = gv(row)
      do l=2, self%ngc ! deeper layers copy the first along the normal
         if (i < 1)            phi(1-l,j,k,b) = gv(row)
         if (i > self%nijk(1)) phi(self%nijk(1)+l,j,k,b) = gv(row)
         if (j < 1)            phi(i,1-l,k,b) = gv(row)
         if (j > self%nijk(2)) phi(i,self%nijk(2)+l,k,b) = gv(row)
         if (k < 1)            phi(i,j,1-l,b) = gv(row)
         if (k > self%nijk(3)) phi(i,j,self%nijk(3)+l,b) = gv(row)
      enddo
   enddo
   if (present(energy)) then
      e_ = 0._R8P
      do c=1, n
         e_ = e_ + self%vol(c) * buf(c) * self%u(c)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, e_, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      energy = e_
   endif
   contains
      function rel_max(rr) result(rm)
      !< max |r / V| over max |f / V|, over the ranks.
      real(R8P), intent(in) :: rr(1:) !< Residual.
      real(R8P)             :: rm     !< Result.
      integer(I4P)          :: ii     !< Counter.

      rm = 0._R8P
      do ii=1, n
         rm = max(rm, abs(rr(ii) / self%vol(ii)))
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, rm, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (fmax > 0._R8P) rm = rm / fmax
      endfunction rel_max
   endsubroutine solve_fcg

   function symmetry_defect(self, x, y) result(defect)
   !< |x . A y - y . A x| / (|x| |A y|) over the ranks, for two leaf vectors: the operator's symmetry, verified.
   class(poisson_mg_object), intent(inout) :: self         !< The solver.
   real(R8P),                intent(in)    :: x(1:)        !< First vector.
   real(R8P),                intent(in)    :: y(1:)        !< Second vector.
   real(R8P)                               :: defect       !< Result.
   real(R8P), allocatable                  :: ax(:), ay(:) !< Products.
   integer(I4P)                            :: n            !< Local length.

   n = self%a%nrow
   allocate(ax(max(n,1)), ay(max(n,1)))
   call self%apply(self%a, x, ax)
   call self%apply(self%a, y, ay)
   defect = abs(gdot(x, ay, n) - gdot(y, ax, n)) / (sqrt(gdot(x, x, n)) * sqrt(gdot(ay, ay, n)))
   endfunction symmetry_defect

   subroutine apply_a(self, x, y)
   !< y = A x on this rank's leaf cells (V-weighted).
   class(poisson_mg_object), intent(in)    :: self  !< The solver.
   real(R8P),                intent(in)    :: x(1:) !< Leaf vector.
   real(R8P),                intent(inout) :: y(1:) !< Result.

   call self%apply(self%a, x, y)
   endsubroutine apply_a

   subroutine apply_b(self, g, y)
   !< y = B g on this rank's leaf cells: every boundary face adds its weight times its boundary value to its row
   !< (`A phi + B g` is the operator with the boundary values g).
   class(poisson_mg_object), intent(in)    :: self  !< The solver.
   real(R8P),                intent(in)    :: g(1:) !< Boundary value per boundary face.
   real(R8P),                intent(inout) :: y(1:) !< Result, a leaf vector.
   integer(I4P)                            :: i     !< Counter.

   y(1:self%a%nrow) = 0._R8P
   do i=1, self%nbnd
      y(self%bnd_row(i)) = y(self%bnd_row(i)) + self%bnd_w(i) * g(i)
   enddo
   endsubroutine apply_b

   subroutine apply_bt(self, x, s)
   !< s = B^T x: every boundary face takes its weight times the value of its row.
   class(poisson_mg_object), intent(in)    :: self  !< The solver.
   real(R8P),                intent(in)    :: x(1:) !< Leaf vector.
   real(R8P),                intent(inout) :: s(1:) !< Result, one value per boundary face.
   integer(I4P)                            :: i     !< Counter.

   do i=1, self%nbnd
      s(i) = self%bnd_w(i) * x(self%bnd_row(i))
   enddo
   endsubroutine apply_bt

   ! private methods
   subroutine apply(self, m, x, y)
   !< y = M x: fetch the halo of x from the owners, then the local rows.
   class(poisson_mg_object), intent(in)    :: self   !< The solver.
   type(dcsr),               intent(in)    :: m      !< Operator.
   real(R8P),                intent(in)    :: x(1:)  !< Local column-space vector.
   real(R8P),                intent(inout) :: y(1:)  !< Result (rows).
   real(R8P), allocatable                  :: sbuf(:), halo(:) !< Send buffer, halo values.
   integer(I4P)                            :: i, k   !< Counters.
   integer(I4P)                            :: ierr   !< MPI error status.
   real(R8P)                               :: acc    !< Accumulator.

   allocate(sbuf(max(size(m%sidx),1)), halo(max(m%nhalo,1)))
   do i=1, size(m%sidx)
      sbuf(i) = x(m%sidx(i))
   enddo
   call MPI_ALLTOALLV(sbuf, m%scount, m%sdispl, MPI_REAL8, halo, m%rcount, m%rdispl, MPI_REAL8, MPI_COMM_WORLD, ierr)
   do i=1, m%nrow
      acc = 0._R8P
      do k=m%rowptr(i), m%rowptr(i+1)-1
         if (m%col(k) <= m%nloc) then
            acc = acc + m%val(k) * x(m%col(k))
         else
            acc = acc + m%val(k) * halo(m%col(k) - m%nloc)
         endif
      enddo
      y(i) = acc
   enddo
   if (.false.) y(1) = y(1) + self%dx0(1)
   endsubroutine apply

   subroutine assemble(self, t, nrow, space, m)
   !< Finish a distributed operator: rows sorted; column keys resolved to local slots, or to halo slots grouped by
   !< owner; the exchange plan agreed with the owners (they receive the keys and return their local indices).
   class(poisson_mg_object), intent(in)    :: self  !< The solver.
   type(triplets),           intent(inout) :: t     !< Triplets (emptied).
   integer(I4P),             intent(in)    :: nrow  !< Rows.
   integer(I4P),             intent(in)    :: space !< Column space.
   type(dcsr),               intent(inout) :: m     !< Operator.
   type(hash_table)          :: hh                  !< Halo key -> provisional halo id.
   integer(I8P), allocatable :: hkey(:), ordered(:), skey(:) !< Halo keys: provisional, by slot, requested.
   integer(I4P), allocatable :: hown(:), slot(:), prov(:), fill(:), next(:) !< Owners, slots, column ids, cursors.
   integer(I4P)              :: e, own, lidx, nh, i, r, k, nsend !< Counters.
   integer(I4P)              :: ierr                !< MPI error status.

   m%nrow = nrow
   if (space == SPACE_COMP) then
      m%nloc = self%nleaf_loc * self%ncb
   else
      m%nloc = self%nloc_lev(space) * self%ncb_of(space)
   endif
   allocate(prov(max(t%n,1)), hkey(max(t%n,1)), hown(max(t%n,1)))
   call hash_init(hh, max(t%n,1))
   nh = 0
   do e=1, t%n
      call self%locate(space, t%k(e), own, lidx)
      if (own == self%me) then
         prov(e) = lidx
      else
         i = hash_get(hh, t%k(e))
         if (i == 0) then
            nh = nh + 1
            hkey(nh) = t%k(e)
            hown(nh) = own
            call hash_put(hh, t%k(e), nh)
            i = nh
         endif
         prov(e) = -i
      endif
   enddo
   call hash_free(hh)
   ! halo slots grouped by owner
   m%nhalo = nh
   allocate(m%rcount(0:self%nr-1), m%rdispl(0:self%nr-1), fill(0:self%nr-1), slot(max(nh,1)), ordered(max(nh,1)))
   m%rcount = 0
   do i=1, nh
      m%rcount(hown(i)) = m%rcount(hown(i)) + 1
   enddo
   m%rdispl(0) = 0
   do r=1, self%nr-1
      m%rdispl(r) = m%rdispl(r-1) + m%rcount(r-1)
   enddo
   fill = m%rdispl
   do i=1, nh
      fill(hown(i)) = fill(hown(i)) + 1
      slot(i) = fill(hown(i))
      ordered(slot(i)) = hkey(i)
   enddo
   ! the requests: the owners translate the keys into their local indices
   allocate(m%scount(0:self%nr-1), m%sdispl(0:self%nr-1))
   call MPI_ALLTOALL(m%rcount, 1, MPI_INTEGER, m%scount, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   m%sdispl(0) = 0
   do r=1, self%nr-1
      m%sdispl(r) = m%sdispl(r-1) + m%scount(r-1)
   enddo
   nsend = sum(m%scount)
   allocate(skey(max(nsend,1)), m%sidx(nsend))
   call MPI_ALLTOALLV(ordered, m%rcount, m%rdispl, MPI_INTEGER8, skey, m%scount, m%sdispl, MPI_INTEGER8, &
                      MPI_COMM_WORLD, ierr)
   do i=1, nsend
      call self%locate(space, skey(i), own, lidx)
      if (own /= self%me) call mpih%error_stop(msg='poisson_mg_object: a halo request reached the wrong owner')
      m%sidx(i) = lidx
   enddo
   ! rows: counting sort of the triplets
   allocate(m%rowptr(nrow+1), m%col(max(t%n,1)), m%val(max(t%n,1)), next(max(nrow,1)))
   m%rowptr = 0
   do e=1, t%n
      m%rowptr(t%r(e)+1) = m%rowptr(t%r(e)+1) + 1
   enddo
   m%rowptr(1) = 1
   do i=1, nrow
      m%rowptr(i+1) = m%rowptr(i+1) + m%rowptr(i)
   enddo
   if (nrow > 0) next(1:nrow) = m%rowptr(1:nrow)
   do e=1, t%n
      k = next(t%r(e))
      next(t%r(e)) = k + 1
      if (prov(e) > 0) then
         m%col(k) = prov(e)
      else
         m%col(k) = m%nloc + slot(-prov(e))
      endif
      m%val(k) = t%v(e)
   enddo
   t%n = 0
   endsubroutine assemble

   subroutine a_row(self, lv, cc, keys, vals)
   !< The row of the composite operator at cell cc of level lv (V-weighted): same-level faces, and the fine and the coarse
   !< side of 2:1 faces (see the module).
   class(poisson_mg_object),  intent(in)    :: self    !< The solver.
   integer(I4P),              intent(in)    :: lv      !< Level.
   integer(I4P),              intent(in)    :: cc(3)   !< Cell coordinates at lv.
   integer(I8P), allocatable, intent(inout) :: keys(:) !< Column keys.
   real(R8P),    allocatable, intent(inout) :: vals(:) !< Values.
   integer(I4P) :: p, d, sgn, cls, o, ch, nf          !< Counters.
   integer(I4P) :: nb3(3), cp(3), f3(3)               !< Coordinates.
   integer(I8P) :: ky                                  !< This cell's key.
   integer(I8P), allocatable :: keys_f(:)              !< The fine cells of a coarse face.
   real(R8P)    :: w, wf                               !< Weights.

   if (allocated(keys)) deallocate(keys, vals)
   allocate(keys(0), vals(0))
   ky = ckey(lv, cc)
   do p=1, self%nd
      d = self%act(p)
      do sgn=-1, 1, 2
         nb3 = cc ; nb3(d) = nb3(d) + sgn
         nb3 = self%wrap(lv, nb3)
         if (self%outside(lv, nb3)) then ! a domain face: the homogeneous row of its kind, the value to the right side
            w = self%face_weight(lv, d)
            select case(self%bc(merge(1, 2, sgn < 0),d))
            case(FLAIL_BC_DIRICHLET)
               keys = [keys, ky] ; vals = [vals, -2._R8P * w]
            case(FLAIL_BC_DIRICHLET_GHOST)
               keys = [keys, ky] ; vals = [vals, -w]
            endselect ! Neumann: no flux
            cycle
         endif
         cls = self%classify(lv, nb3)
         select case(cls)
         case(CLS_SAME)
            w = self%face_weight(lv, d)
            keys = [keys, ky, ckey(lv, nb3)] ; vals = [vals, -w, w]
         case(CLS_COARSER) ! the fine side: the mimetic mean flux (2/3) w (mean of the face's fine cells - phi_C)
            w = self%face_weight(lv, d)
            cp = self%parent3(nb3)
            call self%face_cells(lv, cc, d, keys_f)
            nf = size(keys_f)
            do o=1, nf
               keys = [keys, keys_f(o)] ; vals = [vals, -2._R8P / 3._R8P * w / nf]
            enddo
            keys = [keys, ckey(lv-1, cp)] ; vals = [vals, 2._R8P / 3._R8P * w]
         case(CLS_FINER) ! the coarse side: the adjoint, (2/3) w_f (sum of the face's fine cells - n phi_C)
            wf = self%face_weight(lv+1, d)
            f3 = 2 * cc
            do ch=1, self%nd
               f3(self%act(ch)) = 2 * cc(self%act(ch))
            enddo
            f3(d) = 2 * nb3(d) + merge(0, 1, sgn > 0) ! the children of the neighbour touching this face
            f3 = self%wrap(lv+1, f3)
            call self%face_cells(lv+1, f3, d, keys_f)
            nf = size(keys_f)
            do o=1, nf
               keys = [keys, keys_f(o), ky] ; vals = [vals, 2._R8P / 3._R8P * wf, -2._R8P / 3._R8P * wf]
            enddo
         endselect
      enddo
   enddo
   endsubroutine a_row

   subroutine face_cells(self, lf, fc, d, keys_f)
   !< The keys of the 2^(nd-1) fine cells (level lf) on the coarse face that contains the fine cell fc's face normal
   !< to d: fc's normal index, every pair of tangential indices of its coarse parent.
   class(poisson_mg_object),  intent(in)  :: self      !< The solver.
   integer(I4P),              intent(in)  :: lf        !< Fine level.
   integer(I4P),              intent(in)  :: fc(3)     !< A fine cell of the face.
   integer(I4P),              intent(in)  :: d         !< Normal direction.
   integer(I8P), allocatable, intent(out) :: keys_f(:) !< Keys.
   integer(I4P)                           :: o, ch, a_, nt !< Counters.
   integer(I4P)                           :: c2(3), tang(2) !< Coordinates, tangential directions.

   nt = 0
   do ch=1, self%nd
      if (self%act(ch) == d) cycle
      nt = nt + 1
      tang(nt) = self%act(ch)
   enddo
   allocate(keys_f(2**nt))
   do o=0, 2**nt-1
      c2 = fc
      do a_=1, nt
         c2(tang(a_)) = 2 * (fc(tang(a_)) / 2) + ibits(o, a_-1, 1)
      enddo
      keys_f(o+1) = ckey(lf, c2)
   enddo
   endsubroutine face_cells

   subroutine ghost_green(self, s, g)
   !< Isolated: `g = -G s` on this rank's boundary faces, `G` the inverse of the V-weighted coarsest-level operator on
   !< the infinite lattice, from ghost cell to ghost cell (`G(n) / H`, the unit lattice's Green's function over the
   !< coarsest spacing). Rank 0 gathers the charges into the coarsest box widened by the ghost layer (every face has
   !< its own ghost cell: a boundary leaf is at the coarsest level), convolves with the padded FFT (`gfft`, which solves
   !< `L_H phi = factor f` by the kernel `H^2 G`, so `factor = -1 / H^3`), and scatters the values back.
   class(poisson_mg_object), intent(inout) :: self  !< The solver.
   real(R8P),                intent(in)    :: s(1:) !< Screening charge per boundary face.
   real(R8P),                intent(inout) :: g(1:) !< Ghost value per boundary face.
   integer(I4P), allocatable :: counts(:), displs(:) !< Gather counts and displacements.
   integer(I4P), allocatable :: gi_all(:)            !< Rank 0: the ghost cells of every face.
   real(R8P),    allocatable :: v_all(:)             !< Rank 0: the charges, then the values, of every face.
   real(R8P),    allocatable :: box(:), out(:)       !< Rank 0: the widened box, source and potential.
   real(R8P)                 :: h(3)                 !< Coarsest spacing.
   integer(I4P)              :: ntot, nbl, c, r, ierr !< Faces, blocks of the box solve, counters, MPI status.

   allocate(counts(self%nr), displs(self%nr))
   call MPI_GATHER(self%nbnd, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   displs = 0 ; ntot = 0
   if (self%me == 0) then
      do r=2, self%nr
         displs(r) = displs(r-1) + counts(r-1)
      enddo
      ntot = sum(counts)
   endif
   allocate(gi_all(max(ntot,1)), v_all(max(ntot,1)))
   call MPI_GATHERV(self%bnd_gi, self%nbnd, MPI_INTEGER, gi_all, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_GATHERV(s, self%nbnd, MPI_REAL8, v_all, counts, displs, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   nbl = 0
   if (self%me == 0) then
      nbl = 1
      allocate(box(product(self%nc + 2)), out(product(self%nc + 2)))
      box = 0._R8P
      do c=1, ntot
         box(gi_all(c)) = box(gi_all(c)) + v_all(c)
      enddo
   else
      allocate(box(1), out(1))
   endif
   h = self%dx_of(self%lmin)
   call self%gfft%solve_blocks(nbl=nbl, coords=reshape([0_I4P, 0_I4P, 0_I4P], [3, 1]), buf=box, &
                               factor=-1._R8P / h(1)**3, gw=0_I4P, out=out)
   if (self%me == 0) then
      do c=1, ntot
         v_all(c) = out(gi_all(c))
      enddo
   endif
   call MPI_SCATTERV(v_all, counts, displs, MPI_REAL8, g, self%nbnd, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   endsubroutine ghost_green

   pure function cell_volume(self, lv) result(v)
   !< Cell volume at level lv (the product over the active directions).
   class(poisson_mg_object), intent(in) :: self !< The solver.
   integer(I4P),             intent(in) :: lv   !< Level.
   real(R8P)                            :: v    !< Volume.
   real(R8P)                            :: dx(3) !< Spacing.
   integer(I4P)                         :: p    !< Counter.

   dx = self%dx_of(lv)
   v = 1._R8P
   do p=1, self%nd
      v = v * dx(self%act(p))
   enddo
   endfunction cell_volume

   function classify(self, lv, nb3) result(cls)
   !< Classify the level-lv cell position nb3: a leaf of lv, refined, or inside a leaf of lv-1 (fatal otherwise: the
   !< grid is not 2:1 balanced across faces).
   class(poisson_mg_object), intent(in) :: self   !< The solver.
   integer(I4P),             intent(in) :: lv     !< Level.
   integer(I4P),             intent(in) :: nb3(3) !< Cell coordinates.
   integer(I4P)                         :: cls    !< Class.
   integer(I4P)                         :: id     !< Node id.

   id = hash_get(self%node_hash, bkey([nb3 / self%nijk, lv]))
   if (id > 0) then
      cls = merge(CLS_SAME, CLS_FINER, self%node_leaf(id))
      return
   endif
   if (lv > self%lmin) then
      if (hash_get(self%leaf_hash, bkey([self%parent3(nb3) / self%nijk, lv-1])) > 0) then
         cls = CLS_COARSER
         return
      endif
   endif
   cls = 0
   call mpih%error_stop(msg='poisson_mg_object: a cell neighbours a leaf two levels coarser (the grid is not 2:1 '// &
                            'balanced across faces)')
   endfunction classify

   recursive subroutine comp_value(self, lv, cc, w, keys, wts)
   !< Append w times the level-lv value of cell cc as leaf-cell weights: the leaf itself; the mean of the leaves
   !< inside, if refined; the containing coarser leaf, if inside one.
   class(poisson_mg_object),  intent(in)    :: self    !< The solver.
   integer(I4P),              intent(in)    :: lv      !< Level.
   integer(I4P),              intent(in)    :: cc(3)   !< Cell coordinates.
   real(R8P),                 intent(in)    :: w       !< Weight.
   integer(I8P), allocatable, intent(inout) :: keys(:) !< Keys.
   real(R8P),    allocatable, intent(inout) :: wts(:)  !< Weights.
   integer(I4P)                             :: id, o, ch, l2 !< Counters.
   integer(I4P)                             :: c2(3)   !< Coordinates.

   id = hash_get(self%node_hash, bkey([cc / self%nijk, lv]))
   if (id > 0) then
      if (self%node_leaf(id)) then
         keys = [keys, ckey(lv, cc)] ; wts = [wts, w]
      else
         do o=0, 2**self%nd-1
            c2 = cc
            do ch=1, self%nd
               c2(self%act(ch)) = 2 * cc(self%act(ch)) + ibits(o, ch-1, 1)
            enddo
            call self%comp_value(lv+1, c2, w / 2**self%nd, keys, wts)
         enddo
      endif
      return
   endif
   c2 = cc
   do l2=lv-1, self%lmin, -1
      c2 = self%parent3(c2)
      if (hash_get(self%leaf_hash, bkey([c2 / self%nijk, l2])) > 0) then
         keys = [keys, ckey(l2, c2)] ; wts = [wts, w]
         return
      endif
   enddo
   call mpih%error_stop(msg='poisson_mg_object: a cell is covered by no leaf')
   endsubroutine comp_value

   pure function dx_of(self, lv) result(dx)
   !< Cell spacing at level lv (null directions keep the level-0 one).
   class(poisson_mg_object), intent(in) :: self  !< The solver.
   integer(I4P),             intent(in) :: lv    !< Level.
   real(R8P)                            :: dx(3) !< Spacing.

   dx = merge(self%dx0, self%dx0 / 2._R8P**lv, self%null)
   endfunction dx_of

   pure function face_weight(self, lv, d) result(w)
   !< The V-weighted face weight at level lv along d: the face area over the centre distance, V / dx_d^2.
   class(poisson_mg_object), intent(in) :: self !< The solver.
   integer(I4P),             intent(in) :: lv   !< Level.
   integer(I4P),             intent(in) :: d    !< Direction.
   real(R8P)                            :: w    !< Weight.
   real(R8P)                            :: dx(3) !< Spacing.

   dx = self%dx_of(lv)
   w = self%cell_volume(lv) / dx(d)**2
   endfunction face_weight

   pure function lin(self, ijk) result(l)
   !< Linear index of cell (i,j,k) in a block, i fastest.
   class(poisson_mg_object), intent(in) :: self   !< The solver.
   integer(I4P),             intent(in) :: ijk(3) !< Cell (1-based).
   integer(I4P)                         :: l      !< Index.

   l = ijk(1) + (ijk(2) - 1) * self%nijk(1) + (ijk(3) - 1) * self%nijk(1) * self%nijk(2)
   endfunction lin

   pure function lin_of(self, lv, ijk) result(l)
   !< Linear index of a cell (1-based) in a block at level lv (`nijk_of(lv)` cells).
   class(poisson_mg_object), intent(in) :: self   !< The solver.
   integer(I4P),             intent(in) :: lv     !< Level.
   integer(I4P),             intent(in) :: ijk(3) !< Cell (1-based).
   integer(I4P)                         :: l      !< Index.
   integer(I4P)                         :: nl(3)  !< Cells per block.

   nl = self%nijk_of(lv)
   l = ijk(1) + (ijk(2) - 1) * nl(1) + (ijk(3) - 1) * nl(1) * nl(2)
   endfunction lin_of

   pure function ncb_of(self, lv) result(n)
   !< Cells per block at level lv.
   class(poisson_mg_object), intent(in) :: self !< The solver.
   integer(I4P),             intent(in) :: lv   !< Level.
   integer(I4P)                         :: n    !< Cells.

   n = product(self%nijk_of(lv))
   endfunction ncb_of

   pure function ncell(self, lv, d) result(n)
   !< Cells along direction d over the domain at level lv: the blocks at max(lv, lmin) times their cells.
   class(poisson_mg_object), intent(in) :: self !< The solver.
   integer(I4P),             intent(in) :: lv   !< Level.
   integer(I4P),             intent(in) :: d    !< Direction.
   integer(I4P)                         :: n    !< Cells.
   integer(I4P)                         :: nl(3) !< Cells per block.

   nl = self%nijk_of(lv)
   n = self%nblk0(d) * 2**max(lv, self%lmin) * nl(d)
   endfunction ncell

   pure function nijk_of(self, lv) result(nl)
   !< Cells per block per direction at level lv: nijk at lmin and above, halved per level below it along the active
   !< directions (the levels inside the blocks, #84 E5c).
   class(poisson_mg_object), intent(in) :: self  !< The solver.
   integer(I4P),             intent(in) :: lv    !< Level.
   integer(I4P)                         :: nl(3) !< Cells per block.

   nl = self%nijk
   if (lv < self%lmin) nl = merge(self%nijk, self%nijk / 2**(self%lmin - lv), self%null)
   endfunction nijk_of

   subroutine locate(self, space, key, owner, lidx)
   !< Owner rank and local index of a cell key in a space (the composite leaves, or a level).
   class(poisson_mg_object), intent(in)  :: self  !< The solver.
   integer(I4P),             intent(in)  :: space !< Space.
   integer(I8P),             intent(in)  :: key   !< Cell key.
   integer(I4P),             intent(out) :: owner !< Owner rank.
   integer(I4P),             intent(out) :: lidx  !< Local index on the owner.
   integer(I4P)                          :: lv, cc(3), id, loc !< Decoded key, ids.
   integer(I4P)                          :: nl(3)              !< Cells per block at the level.

   call kdecode(key, lv, cc)
   if (space == SPACE_COMP) then
      id = hash_get(self%leaf_hash, bkey([cc / self%nijk, lv]))
      if (id == 0) call mpih%error_stop(msg='poisson_mg_object: a composite key names no leaf')
      owner = self%leaf_owner(id)
      loc = self%leaf_loc(id)
   elseif (lv < self%lmin) then ! a level inside the blocks (#84 E5c): the lmin block holding the cell
      if (lv /= space) call mpih%error_stop(msg='poisson_mg_object: a level key names another level')
      nl = self%nijk_of(lv)
      id = hash_get(self%node_hash, bkey([cc / nl, self%lmin]))
      if (id == 0) call mpih%error_stop(msg='poisson_mg_object: a level key names no node')
      owner = self%node_owner(id)
      lidx = (self%node_loc(id) - 1) * self%ncb_of(lv) + self%lin_of(lv, cc - (cc / nl) * nl + 1)
      return
   else
      if (lv /= space) call mpih%error_stop(msg='poisson_mg_object: a level key names another level')
      id = hash_get(self%node_hash, bkey([cc / self%nijk, lv]))
      if (id == 0) call mpih%error_stop(msg='poisson_mg_object: a level key names no node')
      owner = self%node_owner(id)
      loc = self%node_loc(id)
   endif
   lidx = (loc - 1) * self%ncb + self%lin(cc - (cc / self%nijk) * self%nijk + 1)
   endsubroutine locate

   pure function parent3(self, cc) result(pc)
   !< Parent cell coordinates (one level coarser) along the active directions.
   class(poisson_mg_object), intent(in) :: self  !< The solver.
   integer(I4P),             intent(in) :: cc(3) !< Cell coordinates.
   integer(I4P)                         :: pc(3) !< Parent.

   pc = merge(cc, cc / 2, self%null)
   endfunction parent3

   subroutine smooth(self, lv, e, rhs, order)
   !< Red-black Gauss-Seidel sweeps on level lv: A_l e = rhs (the ghost term already in rhs).
   class(poisson_mg_object), intent(in)    :: self    !< The solver.
   integer(I4P),             intent(in)    :: lv      !< Level.
   real(R8P),                intent(inout) :: e(1:)   !< Correction.
   real(R8P),                intent(in)    :: rhs(1:) !< Right-hand side.
   integer(I4P),             intent(in)    :: order   !< 0: colour 0 first, 1: colour 1 first.
   real(R8P), allocatable                  :: sbuf(:), halo(:) !< Exchange buffers.
   integer(I4P)                            :: sw, cc_, colr, i, k, ierr !< Counters.
   real(R8P)                               :: acc     !< Off-diagonal sum.

   associate(m=>self%al(lv))
   allocate(sbuf(max(size(m%sidx),1)), halo(max(m%nhalo,1)))
   do sw=1, self%sweeps
      do cc_=0, 1
         colr = mod(cc_ + order, 2)
         do i=1, size(m%sidx)
            sbuf(i) = e(m%sidx(i))
         enddo
         call MPI_ALLTOALLV(sbuf, m%scount, m%sdispl, MPI_REAL8, halo, m%rcount, m%rdispl, MPI_REAL8, &
                            MPI_COMM_WORLD, ierr)
         do i=1, m%nrow
            if (nint(self%colour(lv)%v(i)) /= colr) cycle
            acc = 0._R8P
            do k=m%rowptr(i), m%rowptr(i+1)-1
               if (m%col(k) == i) cycle
               if (m%col(k) <= m%nloc) then
                  acc = acc + m%val(k) * e(m%col(k))
               else
                  acc = acc + m%val(k) * halo(m%col(k) - m%nloc)
               endif
            enddo
            e(i) = (rhs(i) - acc) / self%diag(lv)%v(i)
         enddo
      enddo
   enddo
   endassociate
   endsubroutine smooth

   subroutine vcycle(self, r_leaf, z_leaf)
   !< The preconditioner: z = M r, one V-cycle over the levels (see the module), down to the root (lbot), solved
   !< exactly by FFT on every rank (#84 E5c: below lmin the levels coarsen inside the blocks).
   class(poisson_mg_object), intent(inout) :: self       !< The solver.
   real(R8P),                intent(in)    :: r_leaf(1:) !< Leaf residual (V-weighted).
   real(R8P),                intent(inout) :: z_leaf(1:) !< Leaf correction.
   type(rvec),   allocatable :: r(:), e(:)                !< Level residuals and corrections.
   real(R8P),    allocatable :: tmp(:), res(:), buf(:), out(:) !< Work vectors.
   integer(I4P), allocatable :: coords(:,:)               !< Root blocks.
   integer(I4P)              :: lv, b, q, n, nc           !< Counters.
   real(R8P)                 :: v0                        !< Root cell volume.

   allocate(r(self%lbot:self%lmax), e(self%lbot:self%lmax))
   do lv=self%lbot, self%lmax
      n = max(self%nloc_lev(lv) * self%ncb_of(lv), 1)
      allocate(r(lv)%v(n), e(lv)%v(n))
      r(lv)%v = 0._R8P ; e(lv)%v = 0._R8P
   enddo
   ! the leaf residual into the level spaces
   do b=1, self%nleaf_loc
      lv = self%leaf_lev(b)
      q = self%leaf_lev_loc(b)
      r(lv)%v((q-1)*self%ncb+1:q*self%ncb) = r_leaf((b-1)*self%ncb+1:b*self%ncb)
   enddo
   ! down
   do lv=self%lmax, self%lbot+1, -1
      n = self%nloc_lev(lv) * self%ncb_of(lv)
      nc = self%nloc_lev(lv-1) * self%ncb_of(lv-1)
      call self%smooth(lv, e(lv)%v, r(lv)%v, 0_I4P)
      allocate(tmp(max(n,1)), res(max(n,1)))
      call self%apply(self%al(lv), e(lv)%v, tmp)
      if (n > 0) res(1:n) = r(lv)%v(1:n) - tmp(1:n)
      deallocate(tmp)
      allocate(tmp(max(nc,1)))
      call self%apply(self%rl(lv), res, tmp)
      if (nc > 0) r(lv-1)%v(1:nc) = r(lv-1)%v(1:nc) + tmp(1:nc)
      deallocate(tmp, res)
   enddo
   ! the root: exact, by FFT on every rank (L e = r / V, the mean removed)
   n = self%nloc_lev(self%lbot)
   nc = n * self%ncb_of(self%lbot)
   v0 = self%cell_volume(self%lbot)
   allocate(coords(3,max(n,1)), buf(max(nc,1)), out(max(nc,1)))
   coords = 0
   do q=1, n
      coords(:,q) = self%node_c(1:3,self%lev_node(q,self%lbot))
   enddo
   buf = 0._R8P
   if (nc > 0) buf(1:nc) = r(self%lbot)%v(1:nc) / v0
   call self%fft%solve_blocks(nbl=n, coords=coords, buf=buf, factor=1._R8P, gw=0_I4P, out=out)
   if (nc > 0) e(self%lbot)%v(1:nc) = out(1:nc)
   ! up
   do lv=self%lbot+1, self%lmax
      n = self%nloc_lev(lv) * self%ncb_of(lv)
      allocate(tmp(max(n,1)), res(max(n,1)))
      call self%apply(self%pl(lv), e(lv-1)%v, tmp)
      if (n > 0) e(lv)%v(1:n) = e(lv)%v(1:n) + tmp(1:n)
      if (lv > self%lmin) then ! a face to a cell outside the level reads the coarser correction
         call self%apply(self%el(lv), e(lv-1)%v, res)
         if (n > 0) tmp(1:n) = r(lv)%v(1:n) - res(1:n)
      elseif (n > 0) then ! at lmin and below every cell is in the level
         tmp(1:n) = r(lv)%v(1:n)
      endif
      call self%smooth(lv, e(lv)%v, tmp, 1_I4P)
      deallocate(tmp, res)
   enddo
   ! the leaf correction
   do b=1, self%nleaf_loc
      lv = self%leaf_lev(b)
      q = self%leaf_lev_loc(b)
      z_leaf((b-1)*self%ncb+1:b*self%ncb) = e(lv)%v((q-1)*self%ncb+1:q*self%ncb)
   enddo
   endsubroutine vcycle

   pure function outside(self, lv, cc) result(beyond)
   !< A cell position at level lv beyond a domain face along a bounded active direction (never along a periodic one).
   class(poisson_mg_object), intent(in) :: self   !< The solver.
   integer(I4P),             intent(in) :: lv     !< Level.
   integer(I4P),             intent(in) :: cc(3)  !< Cell coordinates.
   logical                              :: beyond !< Result.
   integer(I4P)                         :: p, d   !< Counters.

   beyond = .false.
   do p=1, self%nd
      d = self%act(p)
      if (self%periodic(d)) cycle
      if (cc(d) < 0 .or. cc(d) >= self%ncell(lv, d)) beyond = .true.
   enddo
   endfunction outside

   pure function wrap(self, lv, cc) result(cw)
   !< Periodic wrap of cell coordinates at level lv along the periodic active directions (bounded ones unchanged).
   class(poisson_mg_object), intent(in) :: self  !< The solver.
   integer(I4P),             intent(in) :: lv    !< Level.
   integer(I4P),             intent(in) :: cc(3) !< Cell coordinates.
   integer(I4P)                         :: cw(3) !< Wrapped.
   integer(I4P)                         :: p, d  !< Counters.

   cw = cc
   do p=1, self%nd
      d = self%act(p)
      if (.not.self%periodic(d)) cycle
      cw(d) = modulo(cc(d), self%ncell(lv, d))
   enddo
   endfunction wrap

   ! module helpers
   function gdot(x, y, n) result(s)
   !< sum_i x_i y_i over the local n entries and the ranks.
   real(R8P),    intent(in) :: x(1:) !< First vector.
   real(R8P),    intent(in) :: y(1:) !< Second vector.
   integer(I4P), intent(in) :: n     !< Local entries.
   real(R8P)                :: s     !< Result.
   integer(I4P)             :: i     !< Counter.
   integer(I4P)             :: ierr  !< MPI error status.

   s = 0._R8P
   do i=1, n
      s = s + x(i) * y(i)
   enddo
   call MPI_ALLREDUCE(MPI_IN_PLACE, s, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
   endfunction gdot

   pure function bkey(c4) result(key)
   !< Key of a block or cell: (i, j, k, level) packed in 19 + 19 + 19 + 7 bits.
   integer(I4P), intent(in) :: c4(4) !< Coordinates (i, j, k, level).
   integer(I8P)             :: key   !< Key.

   key = int(c4(4), I8P) * 2_I8P**57 + int(c4(1), I8P) * 2_I8P**38 + int(c4(2), I8P) * 2_I8P**19 + int(c4(3), I8P)
   endfunction bkey

   pure function ckey(lv, cc) result(key)
   !< Key of a cell at level lv.
   integer(I4P), intent(in) :: lv    !< Level.
   integer(I4P), intent(in) :: cc(3) !< Cell coordinates.
   integer(I8P)             :: key   !< Key.

   key = bkey([cc, lv + KEY_LEVEL_OFFSET]) ! the levels inside the blocks may be negative (#84 E5c)
   endfunction ckey

   pure subroutine kdecode(key, lv, cc)
   !< Decode a cell key.
   integer(I8P), intent(in)  :: key   !< Key.
   integer(I4P), intent(out) :: lv    !< Level.
   integer(I4P), intent(out) :: cc(3) !< Cell coordinates.

   lv = int(key / 2_I8P**57, I4P) - KEY_LEVEL_OFFSET
   cc(1) = int(mod(key / 2_I8P**38, 2_I8P**19), I4P)
   cc(2) = int(mod(key / 2_I8P**19, 2_I8P**19), I4P)
   cc(3) = int(mod(key, 2_I8P**19), I4P)
   endsubroutine kdecode

   subroutine hash_init(h, n)
   !< A table for about n keys (capacity: the first prime above 2 n).
   type(hash_table), intent(inout) :: h        !< Table.
   integer(I4P),     intent(in)    :: n        !< Expected keys.
   integer(I4P)                    :: c, f     !< Capacity candidate, factor.
   logical                         :: is_prime !< Primality.

   c = max(2 * n + 1, 7)
   do
      is_prime = .true.
      do f=2, int(sqrt(real(c))) + 1
         if (f < c .and. mod(c, f) == 0) then
            is_prime = .false.
            exit
         endif
      enddo
      if (is_prime) exit
      c = c + 1
   enddo
   call hash_free(h)
   h%cap = c
   allocate(h%keys(0:c-1), h%vals(0:c-1))
   h%keys = -1_I8P
   h%vals = 0
   endsubroutine hash_init

   subroutine hash_free(h)
   !< Free a table.
   type(hash_table), intent(inout) :: h !< Table.

   if (allocated(h%keys)) deallocate(h%keys, h%vals)
   h%cap = 0
   endsubroutine hash_free

   subroutine hash_put(h, key, val)
   !< Insert or overwrite (linear probing).
   type(hash_table), intent(inout) :: h   !< Table.
   integer(I8P),     intent(in)    :: key !< Key (>= 0).
   integer(I4P),     intent(in)    :: val !< Value (> 0).
   integer(I4P)                    :: s   !< Slot.

   s = int(modulo(key, int(h%cap, I8P)), I4P)
   do
      if (h%keys(s) == -1_I8P .or. h%keys(s) == key) exit
      s = mod(s + 1, h%cap)
   enddo
   h%keys(s) = key
   h%vals(s) = val
   endsubroutine hash_put

   pure function hash_get(h, key) result(val)
   !< The value of a key, 0 when absent.
   type(hash_table), intent(in) :: h    !< Table.
   integer(I8P),     intent(in) :: key  !< Key.
   integer(I4P)                 :: val  !< Value.
   integer(I4P)                 :: s, n !< Slot, probes.

   val = 0
   if (h%cap == 0) return
   s = int(modulo(key, int(h%cap, I8P)), I4P)
   do n=1, h%cap
      if (h%keys(s) == -1_I8P) return
      if (h%keys(s) == key) then
         val = h%vals(s)
         return
      endif
      s = mod(s + 1, h%cap)
   enddo
   endfunction hash_get

   subroutine trip_add(t, row, key, val)
   !< Append a triplet.
   type(triplets), intent(inout) :: t     !< Buffer.
   integer(I4P),   intent(in)    :: row   !< Row.
   integer(I8P),   intent(in)    :: key   !< Column key.
   real(R8P),      intent(in)    :: val   !< Value.
   integer(I4P),   allocatable   :: r2(:) !< Grown rows.
   integer(I8P),   allocatable   :: k2(:) !< Grown keys.
   real(R8P),      allocatable   :: v2(:) !< Grown values.

   if (.not.allocated(t%r)) then
      allocate(t%r(1024), t%k(1024), t%v(1024))
   elseif (t%n == size(t%r)) then
      allocate(r2(2*t%n), k2(2*t%n), v2(2*t%n))
      r2(1:t%n) = t%r ; k2(1:t%n) = t%k ; v2(1:t%n) = t%v
      call move_alloc(r2, t%r) ; call move_alloc(k2, t%k) ; call move_alloc(v2, t%v)
   endif
   t%n = t%n + 1
   t%r(t%n) = row
   t%k(t%n) = key
   t%v(t%n) = val
   endsubroutine trip_add
endmodule adam_poisson_mg_object
