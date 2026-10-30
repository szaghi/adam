!< ADAM, inter-realm seam ghost fill across ranks (issue #40).
module adam_seam_exchange
!< ADAM, inter-realm seam ghost fill across ranks (issue #40).
!<
!< A seam ghost of realm `is` is filled from the interior cell of its peer realm that contains the ghost centre. The
!< realms are partitioned over the ranks independently, so that cell may be owned by another rank. The forest splits
!< the seam rows of every realm at topology time:
!<
!<   * local rows (`maps%seam_local_map_ghost_cell`): ghost and peer cell on this rank, copied by the realm's own
!<     `fill_seam_from_peer_forest`;
!<   * cross-rank rows (`maps%seam_mpi_recv_cell` on the realm owning the ghosts, `maps%seam_mpi_send_cell` on the
!<     realm owning the cells), enumerated in one canonical order on every rank, so the rows rank A sends to rank B
!<     for a seam and the rows B receives from A match one to one.
!<
!< Across a 2:1 seam (`coupling = refined`, issue #52) a send row carries a kind: copy, interpolate (the coarse anchor
!< cell and the packed octant and anchor shifts of `seam_meta_pack`, evaluated by `interp_seam_ghost`, the intra-realm
!< coarse->fine fill) or restrict (the base fine cell, mean of its 2x2x2 block, the intra-realm fine->coarse fill), so
!< the owner of the cells sends the ghost values. Same-rank copy rows keep the local path; same-rank interpolate and
!< restrict rows travel as messages to self (rank = this rank in the send/receive rows), without MPI.
!<
!< `seam_fill` moves both kinds for one (realm, peer slot) pair: the local copy, then the owner of the cells packs
!< (`pack_seam_cells_forest`), the buffers travel point to point, the owner of the ghosts unpacks
!< (`unpack_seam_cells_forest`). Each rank calls it for the same pairs in the same order (the peer slots come from the
!< manifest, not from the rows a rank holds); a rank with no rows of a pair returns at once, so the exchange is not
!< collective, but a rank that skipped a pair it has rows for would leave its partners waiting.
use :: adam_field_object, only : interp_seam_ghost, restrict_seam_high_order, seam_guard_t, seam_guard_prepare, &
                                 seam_guard_value
use :: adam_maps_object,  only : maps_object
use :: adam_mpih_global,  only : mpih
use :: adam_realm_object, only : realm_object, realm_pointer
use :: adam_seam_interpolation_library, only : SEAM_FILL_INJECTION
use :: mpi
use :: penf

implicit none
private
public :: seam_fill
public :: seam_fill_all
public :: seam_fill_batch
public :: seam_ghosts_to_stages
public :: seam_peer_slot
public :: pack_seam_rows
public :: unpack_seam_rows
public :: SEAM_ROW_COPY, SEAM_ROW_INTERPOLATE, SEAM_ROW_RESTRICT

integer(I4P), parameter :: SEAM_TAG = 4040_I4P !< MPI tag of the cross-rank seam messages.
! Kinds of seam row (issue #52), the values of the intra-realm ghost-map flag (`adam_field_object%update_ghost_local`):
integer(I4P), parameter :: SEAM_ROW_COPY        = 1_I4P !< Same cell size: copy the cell.
integer(I4P), parameter :: SEAM_ROW_INTERPOLATE = 4_I4P !< Coarse cell to a fine ghost: interpolate (`interp_seam_ghost`).
integer(I4P), parameter :: SEAM_ROW_RESTRICT    = 8_I4P !< Fine cells to a coarse ghost: mean of the 2x2x2 fine cells.

type :: seam_post_object
   !< A posted seam fill of one (realm, peer slot) pair (issue #53, F3): its buffers and requests live until the wait.
   integer(I4P)              :: is = 0_I4P       !< Realm owning the ghosts.
   integer(I4P)              :: p = 0_I4P        !< Peer slot.
   integer(I4P)              :: n_recv = 0_I4P   !< Cross-rank receive rows.
   integer(I4P)              :: n_req = 0_I4P    !< Posted requests.
   logical                   :: active = .false. !< Cross-rank (or to-self) rows were posted.
   logical                   :: timing = .false. !< Accumulate the phase timers.
   integer(I4P), allocatable :: req(:)           !< MPI requests.
   real(R8P),    allocatable :: send_buf(:)      !< Packed values sent (alive until the wait).
   real(R8P),    allocatable :: recv_buf(:)      !< Values received.
endtype seam_post_object

contains
   subroutine seam_fill(self, realm, p)
   !< Fill the seam ghosts of `self` (the forest realm `self%realm_index`) from its peer slot `p`: local rows, then
   !< cross-rank rows. Only `self` is written; `realm(:)` (which contains `self`) is read for the peer.
   class(realm_object), intent(inout), target :: self   !< Realm owning the ghosts.
   type(realm_pointer), intent(in),    target :: realm(:) !< Forest realms.
   integer(I4P),        intent(in)            :: p      !< Peer slot of `self`.
   type(seam_post_object)                     :: post   !< The posted exchange of the pair.
   real(R8P)                                  :: t0     !< Wait start (MPI_Wtime).

   call seam_fill_post(self=self, realm=realm, p=p, post=post)
   if (.not.post%active) return
   if (post%timing) t0 = MPI_Wtime()
   call MPI_WAITALL(post%n_req, post%req, MPI_STATUSES_IGNORE, mpih%error)
   if (post%timing) self%adam%maps%seam_wtime(4) = self%adam%maps%seam_wtime(4) + (MPI_Wtime() - t0)
   call seam_fill_complete(self=self, post=post)
   endsubroutine seam_fill

   subroutine seam_fill_batch(realm, cadence, wtime)
   !< Fill the seam ghosts of every (realm, peer slot) pair of seam cadence `cadence` with one wait (issue #53, F3):
   !< post every pair (local rows, pack, receives and sends), wait once on all the requests, then unpack every pair in
   !< pair order. Every rank reaches the same pairs in the same order (the slots come from the manifest), and each pair
   !< has its own tag. Batching moves every pack before every unpack: the packs read real (interior) cells only, never
   !< a seam ghost, so the values sent are the ones the pair-by-pair exchange sends. `wtime` (timing only) accumulates
   !< the wait at a barrier put before the posts (the partners' arrival skew, load imbalance) and the single wait.
   type(realm_pointer), intent(inout), target :: realm(:) !< Forest realms.
   integer(I4P),        intent(in)            :: cadence  !< Seam cadence of the pairs to fill.
   real(R8P),           intent(inout)         :: wtime(2) !< Timing: imbalance before the posts, wait.
   type(seam_post_object), allocatable        :: post(:)  !< The posted exchanges, in pair order.
   integer(I4P), allocatable                  :: req(:)   !< All the requests.
   integer(I4P)                               :: is, p    !< Realm and slot.
   integer(I4P)                               :: n, k     !< Pairs counter.
   integer(I4P)                               :: n_req    !< Requests counter.
   logical                                    :: timing   !< Timing on.
   real(R8P)                                  :: t0       !< Phase start (MPI_Wtime).

   n = 0_I4P
   do is=1_I4P, int(size(realm), I4P)
      if (.not.allocated(realm(is)%ptr%adam%maps%seam_local_peer_realm)) cycle
      n = n + count(realm(is)%ptr%adam%maps%seam_local_cadence == cadence)
   enddo
   if (n == 0_I4P) return
   timing = realm(1)%ptr%adam%maps%seam_timing
   if (timing) then
      t0 = MPI_Wtime()
      call MPI_BARRIER(MPI_COMM_WORLD, mpih%error)
      wtime(1) = wtime(1) + (MPI_Wtime() - t0)
   endif
   allocate(post(n)) ; k = 0_I4P
   do is=1_I4P, int(size(realm), I4P)
      if (.not.allocated(realm(is)%ptr%adam%maps%seam_local_peer_realm)) cycle
      do p=1_I4P, int(size(realm(is)%ptr%adam%maps%seam_local_peer_realm), I4P)
         if (realm(is)%ptr%adam%maps%seam_local_cadence(p) /= cadence) cycle
         k = k + 1_I4P
         associate(r=>realm(is)%ptr)
            call seam_fill_post(self=r, realm=realm, p=p, post=post(k))
         endassociate
      enddo
   enddo
   n_req = sum([(post(k)%n_req, k=1_I4P, n)])
   allocate(req(max(n_req, 1_I4P))) ; n_req = 0_I4P
   do k=1_I4P, n
      if (post(k)%n_req == 0_I4P) cycle
      req(n_req+1_I4P:n_req+post(k)%n_req) = post(k)%req(1:post(k)%n_req)
      n_req = n_req + post(k)%n_req
   enddo
   if (timing) t0 = MPI_Wtime()
   call MPI_WAITALL(n_req, req, MPI_STATUSES_IGNORE, mpih%error)
   if (timing) wtime(2) = wtime(2) + (MPI_Wtime() - t0)
   do k=1_I4P, n
      if (.not.post(k)%active) cycle
      associate(r=>realm(post(k)%is)%ptr)
         call seam_fill_complete(self=r, post=post(k))
      endassociate
   enddo
   endsubroutine seam_fill_batch

   subroutine seam_fill_post(self, realm, p, post)
   !< Post the fill of the seam ghosts of `self` from its peer slot `p`: the local rows (copied now), the pack of the
   !< rows `self`'s ghosts receive from the cells this rank owns, the receives and the sends (one message per partner
   !< rank, tagged by the pair) and the message to self. `seam_fill_complete` unpacks once the requests are complete.
   class(realm_object),    intent(inout), target :: self                       !< Realm owning the ghosts.
   type(realm_pointer),    intent(in),    target :: realm(:)                   !< Forest realms.
   integer(I4P),           intent(in)            :: p                          !< Peer slot of `self`.
   type(seam_post_object), intent(out)           :: post                       !< The posted exchange.
   integer(I4P)                                  :: is                         !< Forest index of `self`.
   integer(I4P)                                  :: ip                         !< Peer realm owning the cells.
   integer(I4P)                                  :: q                          !< Slot of `is` in the peer's lists.
   integer(I4P)                                  :: nv                         !< Values per cell.
   integer(I4P)                                  :: n_send                     !< Cross-rank send rows of this rank.
   integer(I4P)                                  :: s0, r0                     !< First send/receive row of the slot.
   integer(I4P)                                  :: r, r_end, rank             !< Row segment bounds and partner rank.
   integer(I4P)                                  :: tag                        !< MPI tag of the pair.
   integer(I4P)                                  :: self_recv(2), self_send(2) !< Same-rank row segments (0 if none).
   real(R8P)                                     :: tw(0:4)                    !< Phase boundaries (MPI_Wtime).

   is = self%realm_index
   ip = self%adam%maps%seam_local_peer_realm(p)
   post%is = is ; post%p = p
   ! The phase boundaries need no device synchronisation: the FNL seam kernels are synchronous and the packed buffers
   ! cross the host through blocking copies (issue #53).
   post%timing = self%adam%maps%seam_timing
   if (post%timing) tw(0) = MPI_Wtime()
   ! local rows: the receiver's own copy kernel
   if (self%adam%maps%seam_local_peer_row_count(p) > 0_I4P) &
      call self%fill_seam_from_peer_forest(peer=realm(ip)%ptr, p_idx=p)
   if (post%timing) tw(1) = MPI_Wtime()
   ! cross-rank rows
   q = seam_peer_slot(realm(ip)%ptr, is)
   post%n_recv = 0_I4P ; n_send = 0_I4P
   if (allocated(self%adam%maps%seam_mpi_recv_row_count)) post%n_recv = self%adam%maps%seam_mpi_recv_row_count(p)
   if (allocated(realm(ip)%ptr%adam%maps%seam_mpi_send_row_count)) n_send = realm(ip)%ptr%adam%maps%seam_mpi_send_row_count(q)
   if (post%n_recv + n_send == 0_I4P) then
      if (post%timing) self%adam%maps%seam_wtime(1) = self%adam%maps%seam_wtime(1) + (tw(1) - tw(0))
      return
   endif
   post%active = .true.
   ! one tag per (realm, slot) pair, so the batched messages of different pairs between two ranks never match each other
   if (p > 100_I4P .or. SEAM_TAG + 100_I4P * is > 32767_I4P) &
      call mpih%error_stop(msg='adam_seam_exchange%seam_fill_post: realm '//trim(str(is, .true.))//', slot '// &
                           trim(str(p, .true.))//' exceeds the seam tag space (100 slots per realm, tag <= 32767)')
   tag = SEAM_TAG + 100_I4P * (is - 1_I4P) + (p - 1_I4P)
   self_recv = 0_I4P ; self_send = 0_I4P
   nv = self%nv
   allocate(post%req(post%n_recv + n_send)) ; post%n_req = 0_I4P
   allocate(post%recv_buf(nv * post%n_recv), post%send_buf(nv * n_send))
   if (post%timing) tw(2) = MPI_Wtime()
   if (n_send > 0_I4P) then
      associate(r_ip=>realm(ip)%ptr)
         call r_ip%pack_seam_cells_forest(p_idx=q, buf=post%send_buf)
      endassociate
   endif
   if (post%timing) tw(3) = MPI_Wtime()
   ! one message per partner rank: the rows of a slot are grouped by rank
   if (post%n_recv > 0_I4P) then
      associate(rows=>self%adam%maps%seam_mpi_recv_cell)
      r0 = self%adam%maps%seam_mpi_recv_row_start(p)
      r = 1_I4P
      do while (r <= post%n_recv)
         rank = rows(r0 + r - 1_I4P, 1) ; r_end = r
         do while (r_end < post%n_recv)
            if (rows(r0 + r_end, 1) /= rank) exit
            r_end = r_end + 1_I4P
         enddo
         if (rank == mpih%myrank) then
            self_recv = [r, r_end]
         else
            post%n_req = post%n_req + 1_I4P
            call MPI_IRECV(post%recv_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, tag, &
                           MPI_COMM_WORLD, post%req(post%n_req), mpih%error)
         endif
         r = r_end + 1_I4P
      enddo
      endassociate
   endif
   if (n_send > 0_I4P) then
      associate(rows=>realm(ip)%ptr%adam%maps%seam_mpi_send_cell)
      s0 = realm(ip)%ptr%adam%maps%seam_mpi_send_row_start(q)
      r = 1_I4P
      do while (r <= n_send)
         rank = rows(s0 + r - 1_I4P, 1) ; r_end = r
         do while (r_end < n_send)
            if (rows(s0 + r_end, 1) /= rank) exit
            r_end = r_end + 1_I4P
         enddo
         if (rank == mpih%myrank) then
            self_send = [r, r_end]
         else
            post%n_req = post%n_req + 1_I4P
            call MPI_ISEND(post%send_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, tag, &
                           MPI_COMM_WORLD, post%req(post%n_req), mpih%error)
         endif
         r = r_end + 1_I4P
      enddo
      endassociate
   endif
   ! message to self: the same-rank interpolate/restrict rows, in the same canonical order on both lists
   if (self_recv(1) > 0_I4P) then
      if (self_send(1) == 0_I4P .or. self_send(2) - self_send(1) /= self_recv(2) - self_recv(1)) &
         call mpih%error_stop(msg='adam_seam_exchange%seam_fill: same-rank seam rows do not match (realm '// &
                              trim(str(is, .true.))//', slot '//trim(str(p, .true.))//')')
      post%recv_buf(nv*(self_recv(1)-1_I4P)+1_I4P:nv*self_recv(2)) = &
         post%send_buf(nv*(self_send(1)-1_I4P)+1_I4P:nv*self_send(2))
   endif
   if (post%timing) then
      tw(4) = MPI_Wtime()
      self%adam%maps%seam_wtime(1:4) = self%adam%maps%seam_wtime(1:4) + (tw(1:4) - tw(0:3))
      self%adam%maps%seam_fills = self%adam%maps%seam_fills + 1_I8P
      self%adam%maps%seam_rows = self%adam%maps%seam_rows + int([post%n_recv, n_send], I8P)
   endif
   endsubroutine seam_fill_post

   subroutine seam_fill_complete(self, post)
   !< Complete a posted fill (`seam_fill_post`), its requests complete: unpack the received rows into `self`'s ghosts.
   class(realm_object),    intent(inout) :: self !< Realm owning the ghosts.
   type(seam_post_object), intent(in)    :: post !< The posted exchange.
   real(R8P)                             :: t0   !< Phase start (MPI_Wtime).

   if (post%timing) t0 = MPI_Wtime()
   if (post%n_recv > 0_I4P) call self%unpack_seam_cells_forest(p_idx=post%p, buf=post%recv_buf)
   if (post%timing) self%adam%maps%seam_wtime(5) = self%adam%maps%seam_wtime(5) + (MPI_Wtime() - t0)
   endsubroutine seam_fill_complete

   subroutine seam_fill_all(self, realm)
   !< Fill the seam ghosts of `self` from every peer slot, whatever the cadence (diagnostics on the committed state).
   class(realm_object), intent(inout), target :: self     !< Realm owning the ghosts.
   type(realm_pointer), intent(in),    target :: realm(:) !< Forest realms.
   integer(I4P)                               :: p        !< Peer slot.

   if (.not.allocated(self%adam%maps%seam_local_peer_realm)) return
   do p=1_I4P, int(size(self%adam%maps%seam_local_peer_realm), I4P)
      call seam_fill(self=self, realm=realm, p=p)
   enddo
   endsubroutine seam_fill_all

   subroutine seam_ghosts_to_stages(maps, cadence, ngc, q, q_rk)
   !< Copy the seam ghosts of `q` into every stage buffer of `q_rk`, for the peer slots of the given cadence (issue
   !< #113): the local rows and the receive rows of each slot, the cells the seam fill writes in this realm.
   !<
   !< The `end_of_step` (alpha) fill writes the committed `q` only, while the stages compute on `q_rk`, whose ghosts
   !< the stage loop never refills from a peer (`rk_object%initialize_stages` copies the interior). Called when a step
   !< opens, right after the forest's alpha fill, it gives every stage of the step the peer state at the start of the
   !< step, the lagged coupling alpha means.
   type(maps_object), intent(in)    :: maps                                !< Maps of the realm owning the ghosts.
   integer(I4P),      intent(in)    :: cadence                             !< Cadence of the slots copied.
   integer(I4P),      intent(in)    :: ngc                                 !< Ghost cells number.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Committed state (v, i, j, k, b).
   real(R8P),         intent(inout) :: q_rk(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Stage buffers (v, i, j, k, b, s).
   integer(I4P)                     :: p, c, row, s, b, i, j, k            !< Counters and cell.

   if (.not.allocated(maps%seam_local_cadence)) return
   do p=1_I4P, int(size(maps%seam_local_cadence), I4P)
      if (maps%seam_local_cadence(p) /= cadence) cycle
      do c=1_I4P, maps%seam_local_peer_row_count(p)
         row = maps%seam_local_peer_row_start(p) + c - 1_I4P
         b = maps%seam_local_map_ghost_cell(row,3) ; i = maps%seam_local_map_ghost_cell(row,7)
         j = maps%seam_local_map_ghost_cell(row,8) ; k = maps%seam_local_map_ghost_cell(row,9)
         do s=1_I4P, int(size(q_rk, dim=6), I4P)
            q_rk(:,i,j,k,b,s) = q(:,i,j,k,b)
         enddo
      enddo
      do c=1_I4P, maps%seam_mpi_recv_row_count(p)
         row = maps%seam_mpi_recv_row_start(p) + c - 1_I4P
         b = maps%seam_mpi_recv_cell(row,2) ; i = maps%seam_mpi_recv_cell(row,3)
         j = maps%seam_mpi_recv_cell(row,4) ; k = maps%seam_mpi_recv_cell(row,5)
         do s=1_I4P, int(size(q_rk, dim=6), I4P)
            q_rk(:,i,j,k,b,s) = q(:,i,j,k,b)
         enddo
      enddo
   enddo
   endsubroutine seam_ghosts_to_stages

   subroutine pack_seam_rows(rows, row_start, row_count, regime, ngc, q, buf, density_var, ratio)
   !< Pack the ghost values of the send rows `[rank, b, i, j, k, kind, meta]` `row_start..row_start+row_count-1` from
   !< the field `q` (the cells' owner active buffer), `nv` values per row in row order (issue #40, #52): a copy of the
   !< cell, the coarse->fine interpolant of `interp_seam_ghost` (the intra-realm fill; injection copies the anchor) or the
   !< restriction from the base cell: the mean of the 2x2x2 fine cells (summed in the intra-realm order, then scaled) when
   !< the meta is 0, otherwise the 6-point Lagrange or the WENO restriction of its packed codes and regime
   !< (`[amr] seam_restriction`, #68 C2, C3). Along an axis the tree does not refine (`ratio` 1, #71) the mean covers one
   !< cell: 2x2x1 on a quadtree.
   integer(I4P), intent(in)  :: rows(1:,1:)                     !< Send rows.
   integer(I4P), intent(in)  :: row_start                       !< First row.
   integer(I4P), intent(in)  :: row_count                       !< Number of rows.
   integer(I4P), intent(in)  :: regime                          !< Coarse->fine fill regime (`maps%seam_ghost_fill`).
   integer(I4P), intent(in)  :: ngc                             !< Ghost cells number.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field (v, i, j, k, b).
   real(R8P),    intent(out) :: buf(1:)                         !< Packed values.
   integer(I4P), intent(in)  :: density_var                     !< Seam fill density guard (`maps%seam_density_var`).
   integer(I4P), intent(in)  :: ratio(1:3)                      !< Refinement ratio per axis (`maps%refine_ratio`).
   integer(I4P)              :: c, row, v, nv                   !< Counters.
   integer(I4P)              :: b, i, j, k                      !< Cell.
   integer(I4P)              :: ic, jc, kc                      !< 2x2x2 counters.
   real(R8P)                 :: total                           !< Restriction sum.
   type(seam_guard_t)        :: guard                           !< The ghost's part of the seam fill guard.
   real(R8P)                 :: value_                          !< A seam ghost interpolant.
   real(R8P)                 :: mean_factor                     !< 1 / product(ratio), exact.

   nv = int(size(q, dim=1), I4P)
   mean_factor = 1._R8P / real(product(ratio), R8P)
   do c=1_I4P, row_count
      row = row_start + c - 1_I4P
      b = rows(row,2) ; i = rows(row,3) ; j = rows(row,4) ; k = rows(row,5)
      select case(rows(row,6))
      case(SEAM_ROW_INTERPOLATE)
         ! the seam fill guard (#101): the ghost's part once, before its variables
         if (regime /= SEAM_FILL_INJECTION .and. density_var > 0_I4P) &
            call seam_guard_prepare(regime=regime, meta=rows(row,7), ngc=ngc, q=q, i_send=i, j_send=j, k_send=k, &
                                    b_send=b, density_var=density_var, guard=guard)
         do v=1_I4P, nv
            if (regime == SEAM_FILL_INJECTION) then
               buf(nv*(c-1)+v) = q(v,i,j,k,b)
            elseif (density_var > 0_I4P) then
               value_ = 0._R8P ! the guarded density row takes the prepared value: its interpolant is not needed
               if (v /= density_var) value_ = interp_seam_ghost(regime=regime, meta=rows(row,7), ngc=ngc, q=q, v=v, &
                                                                i_send=i, j_send=j, k_send=k, b_send=b, density_var=0_I4P)
               buf(nv*(c-1)+v) = seam_guard_value(guard=guard, ngc=ngc, q=q, v=v, i_send=i, j_send=j, k_send=k, &
                                                  b_send=b, density_var=density_var, value_in=value_)
            else
               buf(nv*(c-1)+v) = interp_seam_ghost(regime=regime, meta=rows(row,7), ngc=ngc, q=q, v=v, &
                                                   i_send=i, j_send=j, k_send=k, b_send=b,       &
                                                   density_var=density_var)
            endif
         enddo
      case(SEAM_ROW_RESTRICT)
         if (rows(row,7) /= 0_I4P) then
            ! high-order restriction (#68 C2, C3): the fine owner packed its codes and regime in the meta; 0 is the mean
            do v=1_I4P, nv
               buf(nv*(c-1)+v) = restrict_seam_high_order(meta=rows(row,7), ngc=ngc, q=q, v=v, &
                                                          i_send=i, j_send=j, k_send=k, b_send=b)
            enddo
            cycle
         endif
         do v=1_I4P, nv
            total = 0._R8P
            do kc=0,ratio(3)-1 ; do jc=0,ratio(2)-1 ; do ic=0,ratio(1)-1
               total = total + q(v,i+ic,j+jc,k+kc,b)
            enddo ; enddo ; enddo
            buf(nv*(c-1)+v) = total * mean_factor
         enddo
      case default
         buf(nv*(c-1)+1:nv*c) = q(:,i,j,k,b)
      endselect
   enddo
   endsubroutine pack_seam_rows

   subroutine unpack_seam_rows(rows, row_start, row_count, ngc, buf, q)
   !< Unpack `buf` (`nv` values per row in row order) into the ghost cells of the receive rows `[rank, b, i, j, k]`
   !< `row_start..row_start+row_count-1` of the field `q` (issue #40).
   integer(I4P), intent(in)    :: rows(1:,1:)                   !< Receive rows.
   integer(I4P), intent(in)    :: row_start                     !< First row.
   integer(I4P), intent(in)    :: row_count                     !< Number of rows.
   integer(I4P), intent(in)    :: ngc                           !< Ghost cells number.
   real(R8P),    intent(in)    :: buf(1:)                       !< Packed values.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field (v, i, j, k, b).
   integer(I4P)                :: c, row, nv                    !< Counters.

   nv = int(size(q, dim=1), I4P)
   do c=1_I4P, row_count
      row = row_start + c - 1_I4P
      q(:,rows(row,3),rows(row,4),rows(row,5),rows(row,2)) = buf(nv*(c-1)+1:nv*c)
   enddo
   endsubroutine unpack_seam_rows

   pure function seam_peer_slot(peer, is) result(q)
   !< Return the slot of realm `is` in the peer slots of `peer` (0 if `is` is not a peer of it).
   class(realm_object), intent(in) :: peer !< Realm whose slots are searched.
   integer(I4P),        intent(in) :: is   !< Realm index looked up.
   integer(I4P)                    :: q    !< Slot.

   q = 0_I4P
   if (.not.allocated(peer%adam%maps%seam_local_peer_realm)) return
   do q=1_I4P, int(size(peer%adam%maps%seam_local_peer_realm), I4P)
      if (peer%adam%maps%seam_local_peer_realm(q) == is) return
   enddo
   q = 0_I4P
   endfunction seam_peer_slot
endmodule adam_seam_exchange
