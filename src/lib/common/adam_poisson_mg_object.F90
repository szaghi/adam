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
!< piecewise constant; the coarsest level (every block at the minimum leaf level: the whole domain) solved exactly by
!< the FFT solver (`adam_poisson_fft_object%solve_blocks`; of the caller's backend type with `initialize(coarse=)`). It
!< stops at `tolerance` or at the round-off floor of the source, `100 eps max|factor f| / max|factor (f - mean f)|`,
!< whichever is larger (a small perturbation of a uniform density is known only to that floor, #81 P5b).
!<
!< **Distributed and assembled.** Every rank gathers the leaf list once (block coordinates, owner, local index),
!< derives the shadow blocks and their owners (the owner of the first leaf below, in the gathered order: the same on
!< every rank), and assembles, from the geometry only, the rows of its own cells: the operator, the ghost fill, the level
!< operators and transfers. Columns are global cell keys; every operator owns a halo-exchange plan (the remote columns
!< it reads, fetched from their owners by `MPI_ALLTOALLV`). Correct for any rank count; not tuned for speed.
!<
!< **Scope (fatal otherwise).** Every active direction periodic (isolated boundaries on AMR: #81 P5b); 2:1 balance
!< across faces; fewer than 2^19 cells per direction at the finest level (the cell keys).
!< On output phi holds the leaf solution and its first ghost layer per face: same level, the neighbour; across a 2:1
!< face, the value that makes the one-sided face difference equal the operator's own face gradient,
!< `(mean_f - phi_C) / (1.5 h)` (fine side: `phi_f - (mean_f - phi_C) / 1.5`; coarse side:
!< `phi_C + (mean_f - phi_C) / 0.75`); deeper ghost layers copy the first.
use :: adam_field_object, only : field_object
use :: adam_grid_object, only : grid_object
use :: adam_mpih_global, only : mpih
use :: adam_poisson_fft_object, only : poisson_fft_object
use :: mpi
use :: penf, only : I4P, I8P, R8P, str

implicit none
private
public :: poisson_mg_object

integer(I4P), parameter :: CLS_SAME=1_I4P    !< Neighbour position: a leaf of the same level.
integer(I4P), parameter :: CLS_FINER=2_I4P   !< Neighbour position: refined (finer leaves).
integer(I4P), parameter :: CLS_COARSER=3_I4P !< Neighbour position: inside a leaf one level coarser.
integer(I4P), parameter :: SPACE_COMP=-1_I4P !< Vector space: the composite leaf cells (levels are spaces >= 0).

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
   integer(I4P), allocatable :: nloc_lev(:)            !< This rank's blocks per level [lmin:lmax].
   integer(I4P), allocatable :: lev_node(:,:)          !< This rank's level blocks: node id [max blocks, lmin:lmax].
   integer(I4P)              :: nleaf_loc=0_I4P        !< This rank's leaves.
   integer(I4P), allocatable :: leaf_lev(:)            !< This rank's leaf b: its level.
   integer(I4P), allocatable :: leaf_lev_loc(:)        !< This rank's leaf b: its local index in its level space.
   ! operators
   type(dcsr)                :: a                      !< The composite operator on the leaves.
   type(dcsr)                :: g                      !< Ghost fill: the first ghost layer of the leaves.
   integer(I4P), allocatable :: g_dst(:,:)             !< Ghost fill destinations (i,j,k,b) [4,rows].
   type(dcsr),   allocatable :: al(:)                  !< Level operators [lmin+1:lmax].
   type(dcsr),   allocatable :: el(:)                  !< Level ghost terms from level l-1 [lmin+1:lmax].
   type(dcsr),   allocatable :: rl(:)                  !< Restriction l -> l-1 [lmin+1:lmax].
   type(dcsr),   allocatable :: pl(:)                  !< Prolongation l-1 -> l [lmin+1:lmax].
   type(rvec),   allocatable :: diag(:)                !< Level operator diagonals [lmin+1:lmax].
   type(rvec),   allocatable :: colour(:)              !< Red-black colour (0, 1) of the level cells [lmin+1:lmax].
   real(R8P),    allocatable :: vol(:)                 !< Leaf cell volumes.
   class(poisson_fft_object), allocatable :: fft       !< The coarsest-level solver (FFTW, or the caller's backend).
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
      procedure, pass(self), private :: lin         !< Linear index of a cell in a block.
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
   if (allocated(self%leaf_c))       deallocate(self%leaf_c, self%leaf_owner, self%leaf_loc)
   if (allocated(self%node_c))       deallocate(self%node_c, self%node_owner, self%node_loc, self%node_leaf)
   if (allocated(self%nloc_lev))     deallocate(self%nloc_lev, self%lev_node)
   if (allocated(self%leaf_lev))     deallocate(self%leaf_lev, self%leaf_lev_loc)
   if (allocated(self%al))           deallocate(self%al, self%el, self%rl, self%pl, self%diag, self%colour)
   if (allocated(self%vol))          deallocate(self%vol)
   if (allocated(self%g_dst))        deallocate(self%g_dst)
   if (allocated(self%u))            deallocate(self%u)
   call hash_free(self%leaf_hash)
   call hash_free(self%node_hash)
   self%a = dcsr()
   self%g = dcsr()
   self%has_guess = .false.
   self%is_initialized = .false.
   endsubroutine destroy

   subroutine initialize(self, grid, field, tolerance, max_iterations, coarse)
   !< Gather the leaves, derive the hierarchy, assemble the operator, the ghost fill and the preconditioner, set up the
   !< coarsest-level FFT solver: of the dynamic type of `coarse` when given (the FNL backend passes its cuFFT solver),
   !< the FFTW one otherwise.
   class(poisson_mg_object), intent(inout)        :: self           !< The solver.
   type(grid_object),        intent(in)           :: grid           !< The grid.
   type(field_object),       intent(in)           :: field          !< The field (block coordinates).
   real(R8P),                intent(in), optional :: tolerance      !< Relative residual to stop.
   integer(I4P),             intent(in), optional :: max_iterations !< FCG iterations at most.
   class(poisson_fft_object), intent(in), optional :: coarse        !< Mold of the coarsest-level solver.
   integer(I4P), allocatable :: counts(:), displs(:) !< Allgather counts and displacements.
   integer(I4P), allocatable :: mine(:,:)            !< This rank's block coordinates and index [5,nb].
   integer(I4P), allocatable :: all_(:,:)            !< Every rank's [5,nleaf].
   integer(I4P), allocatable :: cnt(:,:)             !< Level blocks per rank and level.
   integer(I4P)              :: d, n, id, l, b, p, r !< Counters.
   integer(I4P)              :: c(4)                 !< Block coordinates.
   integer(I4P)              :: ierr                 !< MPI error status.

   call self%destroy
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
      if (.not.grid%is_ijk_periodic(d)) call mpih%error_stop(msg='poisson_mg_object%initialize: direction '// &
         trim(str(d, .true.))//' is not periodic; the multigrid solver supports periodic boundaries (isolated: #81 P5b)')
      self%nd = self%nd + 1
      self%act(self%nd) = d
   enddo
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
   allocate(self%nloc_lev(self%lmin:self%lmax), self%lev_node(max(maxval(cnt(self%me,:)),1),self%lmin:self%lmax))
   self%nloc_lev = cnt(self%me,:)
   do id=1, self%nnode
      if (self%node_owner(id) == self%me) self%lev_node(self%node_loc(id),self%node_c(4,id)) = id
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
   call assemble_a
   call assemble_ghost_fill
   call assemble_levels
   ! the coarsest level: every block at lmin, the whole domain
   if (present(coarse)) then
      allocate(self%fft, mold=coarse)
   else
      allocate(poisson_fft_object :: self%fft)
   endif
   call self%fft%initialize_blocks(nblk=self%nblk0 * merge(1, 2**self%lmin, self%null), nijk=self%nijk, ngc=0_I4P, &
                                   dxyz=self%dx_of(self%lmin), level=self%lmin, null_xyz=self%null, isolated=.false.)
   allocate(self%u(max(n*self%ncb,1)))
   self%u = 0._R8P
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

      allocate(dst(4,max(2*self%nd*self%ncb*field%blocks_number,1)))
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
                        cls = self%classify(lv, nb3)
                        if (allocated(keys)) deallocate(keys, wts)
                        allocate(keys(0), wts(0))
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
                           keys = [ckey(lv, cc), ckey(lv-1, cp)] ; wts = [1._R8P, 1._R8P / 1.5_R8P]
                           do e=1, size(keys_f)
                              keys = [keys, keys_f(e)] ; wts = [wts, -1._R8P / 1.5_R8P / size(keys_f)]
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
      !< The preconditioner's level operators, ghost terms and transfers, for every level above lmin.
      integer(I4P) :: lv, lt !< Level, top.

      lt = max(self%lmax, self%lmin+1)
      allocate(self%al(self%lmin+1:lt), self%el(self%lmin+1:lt), self%rl(self%lmin+1:lt), self%pl(self%lmin+1:lt), &
               self%diag(self%lmin+1:lt), self%colour(self%lmin+1:lt))
      do lv=self%lmin+1, self%lmax
         call level_ops(lv)
      enddo
      endsubroutine assemble_levels

      subroutine level_ops(lv)
      !< Level lv: A_l (the plain Laplacian on the level's cells, a face to a cell outside the level reading E_l: the
      !< level-(lv-1) cell containing it), R_l (rows: this rank's shadow cells at lv-1, the sum of their children at
      !< lv) and P_l (rows: the level's cells, their parent at lv-1).
      integer(I4P), intent(in) :: lv !< Level.
      type(triplets)           :: ta, te, tr, tp !< Assembly buffers.
      integer(I4P)             :: q, nloc, i, j, kk, pp, dd, sgn, row, o, id_, ch !< Counters.
      integer(I4P)             :: cc(3), nb3(3), pc(3), cch(3) !< Coordinates.
      real(R8P)                :: w !< Face weight.

      nloc = self%nloc_lev(lv) * self%ncb
      allocate(self%colour(lv)%v(max(nloc,1)), self%diag(lv)%v(max(nloc,1)))
      self%diag(lv)%v = 1._R8P
      if (nloc > 0) self%diag(lv)%v(1:nloc) = 0._R8P
      do q=1, self%nloc_lev(lv)
         id_ = self%lev_node(q,lv)
         do kk=1, self%nijk(3)
            do j=1, self%nijk(2)
               do i=1, self%nijk(1)
                  cc = self%node_c(1:3,id_) * self%nijk + [i-1, j-1, kk-1]
                  row = (q - 1) * self%ncb + self%lin([i, j, kk])
                  self%colour(lv)%v(row) = real(mod(sum(cc), 2), R8P)
                  do pp=1, self%nd
                     dd = self%act(pp)
                     w = self%face_weight(lv, dd)
                     do sgn=-1, 1, 2
                        nb3 = cc ; nb3(dd) = nb3(dd) + sgn
                        nb3 = self%wrap(lv, nb3)
                        call trip_add(ta, row, ckey(lv, cc), -w)
                        self%diag(lv)%v(row) = self%diag(lv)%v(row) - w
                        if (hash_get(self%node_hash, bkey([nb3 / self%nijk, lv])) > 0) then
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
         if (self%node_leaf(id_)) cycle
         do kk=1, self%nijk(3)
            do j=1, self%nijk(2)
               do i=1, self%nijk(1)
                  pc = self%node_c(1:3,id_) * self%nijk + [i-1, j-1, kk-1]
                  row = (q - 1) * self%ncb + self%lin([i, j, kk])
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
      call self%assemble(tr, self%nloc_lev(lv-1) * self%ncb, lv, self%rl(lv))
      call self%assemble(tp, nloc, lv-1, self%pl(lv))
      endsubroutine level_ops

   endsubroutine initialize

   subroutine solve(self, field, q, iv, factor, phi, energy)
   !< Solve `L_h phi = factor (q(iv) - mean(q(iv)))` on the leaves; phi's ghosts filled (see the module).
   class(poisson_mg_object), intent(inout)         :: self                     !< The solver.
   type(field_object),       intent(in)            :: field                    !< The field (coordinates).
   real(R8P),    intent(in)    :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Conservative [nv,i,j,k,nb].
   integer(I4P), intent(in)    :: iv                                           !< Variable of q, the source.
   real(R8P),    intent(in)    :: factor                                       !< Source factor (4 pi G).
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:)  !< Potential [i,j,k,nb].
   real(R8P),    intent(out), optional :: energy                               !< sum_cells V f phi (all ranks).
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
   call self%solve_packed(field=field, buf=buf, factor=factor, phi=phi, energy=energy)
   endsubroutine solve

   subroutine solve_packed(self, field, buf, factor, phi, energy)
   !< Solve from the packed interior source of this rank's leaves (block-major, `i` fastest): FCG on the operator, preconditioned
   !< by the V-cycle, from the previous solution when there is one; then the ghost fill.
   class(poisson_mg_object), intent(inout)         :: self    !< The solver.
   type(field_object),       intent(in)            :: field   !< The field.
   real(R8P),                intent(in)            :: buf(1:) !< Packed interior source.
   real(R8P),                intent(in)            :: factor  !< Source factor.
   real(R8P),    intent(inout) :: phi(1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< Potential [i,j,k,nb].
   real(R8P),                intent(out), optional :: energy  !< sum_cells V f phi (all ranks).
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
   ! zero mean (periodic): f - V sum(f) / sum(V)
   vsum = gdot(self%vol, [(1._R8P, i=1, nn)], n)
   mean = gdot(f, [(1._R8P, i=1, nn)], n) / vsum
   if (n > 0) f(1:n) = f(1:n) - mean * self%vol(1:n)
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
      mean = gdot(z, self%vol, n) / vsum
      if (n > 0) z(1:n) = z(1:n) - mean
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
   mean = gdot(self%u, self%vol, n) / vsum
   if (n > 0) self%u(1:n) = self%u(1:n) - mean
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
   endsubroutine solve_packed

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
      m%nloc = self%nloc_lev(space) * self%ncb
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

   subroutine locate(self, space, key, owner, lidx)
   !< Owner rank and local index of a cell key in a space (the composite leaves, or a level).
   class(poisson_mg_object), intent(in)  :: self  !< The solver.
   integer(I4P),             intent(in)  :: space !< Space.
   integer(I8P),             intent(in)  :: key   !< Cell key.
   integer(I4P),             intent(out) :: owner !< Owner rank.
   integer(I4P),             intent(out) :: lidx  !< Local index on the owner.
   integer(I4P)                          :: lv, cc(3), id, loc !< Decoded key, ids.

   call kdecode(key, lv, cc)
   if (space == SPACE_COMP) then
      id = hash_get(self%leaf_hash, bkey([cc / self%nijk, lv]))
      if (id == 0) call mpih%error_stop(msg='poisson_mg_object: a composite key names no leaf')
      owner = self%leaf_owner(id)
      loc = self%leaf_loc(id)
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
   !< The preconditioner: z = M r, one V-cycle over the levels (see the module).
   class(poisson_mg_object), intent(inout) :: self       !< The solver.
   real(R8P),                intent(in)    :: r_leaf(1:) !< Leaf residual (V-weighted).
   real(R8P),                intent(inout) :: z_leaf(1:) !< Leaf correction.
   type(rvec),   allocatable :: r(:), e(:)                !< Level residuals and corrections.
   real(R8P),    allocatable :: tmp(:), res(:), buf(:), out(:) !< Work vectors.
   integer(I4P), allocatable :: coords(:,:)               !< Coarsest blocks.
   integer(I4P)              :: lv, b, q, n, nc           !< Counters.
   real(R8P)                 :: v0                        !< Coarsest cell volume.

   allocate(r(self%lmin:self%lmax), e(self%lmin:self%lmax))
   do lv=self%lmin, self%lmax
      n = max(self%nloc_lev(lv) * self%ncb, 1)
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
   do lv=self%lmax, self%lmin+1, -1
      n = self%nloc_lev(lv) * self%ncb
      nc = self%nloc_lev(lv-1) * self%ncb
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
   ! coarsest: exact, by FFT (L e = r / V, the mean removed)
   n = self%nloc_lev(self%lmin)
   v0 = self%cell_volume(self%lmin)
   allocate(coords(3,max(n,1)), buf(max(n*self%ncb,1)), out(max(n*self%ncb,1)))
   coords = 0
   do q=1, n
      coords(:,q) = self%node_c(1:3,self%lev_node(q,self%lmin))
   enddo
   buf = 0._R8P
   if (n > 0) buf(1:n*self%ncb) = r(self%lmin)%v(1:n*self%ncb) / v0
   call self%fft%solve_blocks(nbl=n, coords=coords, buf=buf, factor=1._R8P, gw=0_I4P, out=out)
   if (n > 0) e(self%lmin)%v(1:n*self%ncb) = out(1:n*self%ncb)
   ! up
   do lv=self%lmin+1, self%lmax
      n = self%nloc_lev(lv) * self%ncb
      allocate(tmp(max(n,1)), res(max(n,1)))
      call self%apply(self%pl(lv), e(lv-1)%v, tmp)
      if (n > 0) e(lv)%v(1:n) = e(lv)%v(1:n) + tmp(1:n)
      call self%apply(self%el(lv), e(lv-1)%v, res)
      if (n > 0) tmp(1:n) = r(lv)%v(1:n) - res(1:n)
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

   pure function wrap(self, lv, cc) result(cw)
   !< Periodic wrap of cell coordinates at level lv along the active directions.
   class(poisson_mg_object), intent(in) :: self  !< The solver.
   integer(I4P),             intent(in) :: lv    !< Level.
   integer(I4P),             intent(in) :: cc(3) !< Cell coordinates.
   integer(I4P)                         :: cw(3) !< Wrapped.
   integer(I4P)                         :: p, d  !< Counters.

   cw = cc
   do p=1, self%nd
      d = self%act(p)
      cw(d) = modulo(cc(d), self%nblk0(d) * 2**lv * self%nijk(d))
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

   key = bkey([cc, lv])
   endfunction ckey

   pure subroutine kdecode(key, lv, cc)
   !< Decode a cell key.
   integer(I8P), intent(in)  :: key   !< Key.
   integer(I4P), intent(out) :: lv    !< Level.
   integer(I4P), intent(out) :: cc(3) !< Cell coordinates.

   lv = int(key / 2_I8P**57, I4P)
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
