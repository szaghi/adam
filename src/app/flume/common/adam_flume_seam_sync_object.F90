!< ADAM, FLUME per-stage seam flux synchronisation data of the positivity limiter (issue #50, D1).
module adam_flume_seam_sync_object
!< ADAM, FLUME per-stage seam flux synchronisation data of the positivity limiter (issue #50, D1).
!<
!< **Why.** The positivity limiter keeps every forward-Euler stage admissible cell by cell, but at a 2:1 AMR seam the
!< end-of-step Berger-Colella reflux then replaces the coarse face flux by the restricted fine one, a correction no cell
!< factor bounds: on the Balsara-Spicer blast across a refined box it drives coarse seam cells to `rho e < 0` (issue
!< #50). The limited runs therefore synchronise the seam flux at every stage, before the update:
!<
!< - the fine faces of the seam take as backbone `F^LF(q_C, q_f)`, with the coarse donor cell `q_C` as outer state;
!< - the coarse face takes the means of the fine high-order and backbone fluxes, `F_H` and `F_LF` (2x2 fine faces per
!<   coarse face cell, 2x1 on a quadtree face normal to x or y): the coarse backbone update with the mean of
!<   `F^LF(q_C, q_f)` is a convex combination of standard Lax-Friedrichs partial updates of `q_C`, hence admissible, so
!<   the coarse cell factor `Lambda_C` is computed with this pair;
!< - both sides blend with one factor `theta_s = min(Lambda_C, min Lambda_f)`, so each side stays in its corner box and
!<   the coarse flux is exactly the mean of the fine ones: conservative at every stage, the reflux reduces to round-off.
!<
!< **What this object holds.** The register's intra-realm AMR faces are replicated on every rank (`flux_register_object`,
!< issue #28), so per-face skins indexed like the register are complete after one MPI_ALLREDUCE each: the coarse donor
!< states `qc` (SUM: one owner writes), the fine means `flo`, `fhi` (SUM over the fine quadrants), the fine factors
!< `lmin` (MIN) and the seam factors `th` (SUM: the coarse owner writes). The fine seam faces of this rank keep their
!< backbone and high-order fluxes in a compact store (`fine_lo`, `fine_hi`): the blend overwrites the face fluxes before
!< the seam faces are set. The skin cell `c` runs over the two tangential axes, inner fastest, the register order.

! ADAM classes, libraries, parameters
use :: adam_flux_register_object, only : flux_register_object, SEAM_KIND_INTRA_REALM_AMR
! ADAM singleton objects
use :: adam_mpih_global,          only : mpih
! third party modules
use :: mpi
use :: penf,                      only : I4P, R8P

implicit none
private
public :: flume_seam_sync_object
public :: seam_face_cells
public :: seam_fine_to_coarse
public :: seam_skin_index

type :: flume_seam_sync_object
   !< Per-stage seam flux synchronisation data (host).
   integer(I4P)              :: nv = 0_I4P          !< Variables number.
   integer(I4P)              :: ncells = 0_I4P      !< Skin cells of all intra-realm AMR register faces.
   integer(I4P)              :: nfine = 0_I4P       !< Face cells of this rank's fine seam faces.
   integer(I4P), allocatable :: off(:)              !< Skin offset of each register face (-1: not intra-realm AMR).
   integer(I4P), allocatable :: fine_off(:,:)       !< Compact offset of each (block, face) fine seam face (-1: none).
   real(R8P),    allocatable :: qc(:,:)             !< Coarse donor states [nv, ncells].
   real(R8P),    allocatable :: flo(:,:)            !< Mean fine backbone fluxes [nv, ncells].
   real(R8P),    allocatable :: fhi(:,:)            !< Mean fine high-order fluxes [nv, ncells].
   real(R8P),    allocatable :: lmin(:)             !< Smallest fine factor [ncells].
   real(R8P),    allocatable :: th(:)               !< Seam factors [ncells].
   real(R8P),    allocatable :: fine_lo(:,:)        !< Fine seam face backbone fluxes [nv, nfine].
   real(R8P),    allocatable :: fine_hi(:,:)        !< Fine seam face high-order fluxes [nv, nfine].
   contains
      ! public methods
      procedure, pass(self) :: build          !< Size the skins from the register, reset them; .true. if any seam.
      procedure, pass(self) :: reduce_states  !< Complete the coarse donor states over the ranks.
      procedure, pass(self) :: reduce_fluxes  !< Complete the fine means over the ranks.
      procedure, pass(self) :: reduce_factors !< Complete the smallest fine factors over the ranks.
      procedure, pass(self) :: reduce_theta   !< Complete the seam factors over the ranks.
endtype flume_seam_sync_object

contains
   ! public methods
   function build(self, flux_register, register_index, blocks_number, ni, nj, nk, nv) result(active)
   !< Size the skins of the register's intra-realm AMR faces and the compact store of this rank's fine seam faces, and
   !< reset them (`qc`, `flo`, `fhi`, `th` to 0, `lmin` to huge). The result depends on the replicated register only,
   !< so every rank takes the same branch and issues the same collectives.
   class(flume_seam_sync_object), intent(inout) :: self                  !< The synchronisation data.
   type(flux_register_object),    intent(in)    :: flux_register         !< Forest's flux register.
   integer(I4P),                  intent(in)    :: register_index(1:,1:) !< (block, face) signed register index.
   integer(I4P),                  intent(in)    :: blocks_number         !< Actual blocks number.
   integer(I4P),                  intent(in)    :: ni, nj, nk            !< Grid dimensions.
   integer(I4P),                  intent(in)    :: nv                    !< Variables number.
   logical                                      :: active                !< The register has intra-realm AMR faces.
   integer(I4P)                                 :: f, b, fec, s, n       !< Counters.

   active = .false.
   self%ncells = 0_I4P
   self%nfine = 0_I4P
   if (.not.flux_register%is_initialized_) return
   if (flux_register%nfaces <= 0_I4P .or. .not.allocated(flux_register%face)) return
   if (allocated(self%off)) deallocate(self%off)
   allocate(self%off(flux_register%nfaces))
   n = 0_I4P
   do f=1, flux_register%nfaces
      self%off(f) = -1_I4P
      if (flux_register%face(f)%seam_kind /= SEAM_KIND_INTRA_REALM_AMR) cycle
      if (.not.allocated(flux_register%face(f)%F_coarse)) cycle
      self%off(f) = n
      n = n + flux_register%face(f)%nface_cells
   enddo
   if (n == 0_I4P) return
   active = .true.
   self%nv = nv
   self%ncells = n
   call resize2(self%qc, nv, n) ; call resize2(self%flo, nv, n) ; call resize2(self%fhi, nv, n)
   if (allocated(self%lmin)) then
      if (size(self%lmin) /= n) deallocate(self%lmin, self%th)
   endif
   if (.not.allocated(self%lmin)) allocate(self%lmin(n), self%th(n))
   self%qc = 0._R8P ; self%flo = 0._R8P ; self%fhi = 0._R8P ; self%th = 0._R8P ; self%lmin = huge(1._R8P)
   if (allocated(self%fine_off)) deallocate(self%fine_off)
   allocate(self%fine_off(blocks_number, 6))
   self%fine_off = -1_I4P
   do b=1, min(blocks_number, size(register_index, dim=1))
      do fec=1, 6
         s = register_index(b, fec)
         if (s >= 0_I4P .or. -s > flux_register%nfaces) cycle
         if (self%off(-s) < 0_I4P) cycle
         self%fine_off(b, fec) = self%nfine
         self%nfine = self%nfine + seam_face_cells(fec=fec, ni=ni, nj=nj, nk=nk)
      enddo
   enddo
   call resize2(self%fine_lo, nv, max(1_I4P, self%nfine)) ; call resize2(self%fine_hi, nv, max(1_I4P, self%nfine))
   contains
      subroutine resize2(a, n1, n2)
      !< (Re)allocate a rank-2 array to [n1, n2] if its shape differs.
      real(R8P), allocatable, intent(inout) :: a(:,:) !< Array.
      integer(I4P),           intent(in)    :: n1, n2 !< Extents.

      if (allocated(a)) then
         if (size(a, dim=1) == n1 .and. size(a, dim=2) == n2) return
         deallocate(a)
      endif
      allocate(a(n1, n2))
      endsubroutine resize2
   endfunction build

   subroutine reduce_states(self)
   !< Complete the coarse donor states over the ranks (each skin cell written by its coarse owner only).
   class(flume_seam_sync_object), intent(inout) :: self !< The synchronisation data.

   if (mpih%procs_number <= 1_I4P) return
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%qc, size(self%qc), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   endsubroutine reduce_states

   subroutine reduce_fluxes(self)
   !< Complete the fine means over the ranks (each quadrant written by its fine owner only).
   class(flume_seam_sync_object), intent(inout) :: self !< The synchronisation data.

   if (mpih%procs_number <= 1_I4P) return
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%flo, size(self%flo), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%fhi, size(self%fhi), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   endsubroutine reduce_fluxes

   subroutine reduce_factors(self)
   !< Complete the smallest fine factors over the ranks.
   class(flume_seam_sync_object), intent(inout) :: self !< The synchronisation data.

   if (mpih%procs_number <= 1_I4P) return
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%lmin, size(self%lmin), MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, mpih%error)
   endsubroutine reduce_factors

   subroutine reduce_theta(self)
   !< Complete the seam factors over the ranks (each skin cell written by its coarse owner only).
   class(flume_seam_sync_object), intent(inout) :: self !< The synchronisation data.

   if (mpih%procs_number <= 1_I4P) return
   call MPI_ALLREDUCE(MPI_IN_PLACE, self%th, size(self%th), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpih%error)
   endsubroutine reduce_theta

   ! public procedures
   pure function seam_face_cells(fec, ni, nj, nk) result(n)
   !< Cells of a block face `fec` (1..6: -x, +x, -y, +y, -z, +z).
   integer(I4P), intent(in) :: fec        !< Face.
   integer(I4P), intent(in) :: ni, nj, nk !< Grid dimensions.
   integer(I4P)             :: n          !< Face cells.
   !$acc routine seq
   !$omp declare target

   select case((fec + 1_I4P) / 2_I4P)
   case(1_I4P)
      n = nj * nk
   case(2_I4P)
      n = ni * nk
   case default
      n = ni * nj
   endselect
   endfunction seam_face_cells

   pure function seam_skin_index(fec, ni, nj, i, j, k) result(c)
   !< Skin index of the interior cell `(i, j, k)` beside face `fec` (the inverse of `seam_skin_cell`).
   integer(I4P), intent(in) :: fec     !< Face.
   integer(I4P), intent(in) :: ni, nj  !< Grid dimensions.
   integer(I4P), intent(in) :: i, j, k !< Cell.
   integer(I4P)             :: c       !< Skin index.
   !$acc routine seq
   !$omp declare target

   select case((fec + 1_I4P) / 2_I4P)
   case(1_I4P)
      c = (k - 1_I4P) * nj + j
   case(2_I4P)
      c = (k - 1_I4P) * ni + i
   case default
      c = (j - 1_I4P) * ni + i
   endselect
   endfunction seam_skin_index

   pure function seam_fine_to_coarse(fec, ni, nj, nk, ioff, joff, ri, ro, c) result(cc)
   !< Coarse skin cell covering the fine skin cell `c` of a fine block whose quadrant on the coarse face is `(ioff,
   !< joff)` (the convention of `restrict_fine_face_to_quadrant`: `ri x ro` fine face cells per coarse face cell, 2 per
   !< refined tangential axis, 1 with a zero offset along an axis the tree does not refine, z of a quadtree, issue #46).
   integer(I4P), intent(in) :: fec        !< Fine block face.
   integer(I4P), intent(in) :: ni, nj, nk !< Grid dimensions.
   integer(I4P), intent(in) :: ioff, joff !< Quadrant offsets (inner, outer).
   integer(I4P), intent(in) :: ri, ro     !< Refinement ratios along the inner and outer tangential axes.
   integer(I4P), intent(in) :: c          !< Fine skin cell.
   integer(I4P)             :: cc         !< Coarse skin cell.
   integer(I4P)             :: inner_n, outer_n, fi, fo !< Tangential extents, fine tangential indexes.
   !$acc routine seq
   !$omp declare target

   select case((fec + 1_I4P) / 2_I4P)
   case(1_I4P)
      inner_n = nj ; outer_n = nk
   case(2_I4P)
      inner_n = ni ; outer_n = nk
   case default
      inner_n = ni ; outer_n = nj
   endselect
   fi = 1_I4P + mod(c - 1_I4P, inner_n)
   fo = 1_I4P + (c - 1_I4P) / inner_n
   cc = (joff * outer_n / 2_I4P + (fo + ro - 1_I4P) / ro - 1_I4P) * inner_n + ioff * inner_n / 2_I4P + (fi + ri - 1_I4P) / ri
   endfunction seam_fine_to_coarse
endmodule adam_flume_seam_sync_object
