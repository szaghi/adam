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
use :: adam_field_object, only : interp_seam_ghost
use :: adam_mpih_global,  only : mpih
use :: adam_realm_object, only : realm_object
use :: adam_seam_interpolation_library, only : SEAM_FILL_INJECTION
use :: mpi
use :: penf

implicit none
private
public :: seam_fill
public :: seam_fill_all
public :: seam_peer_slot
public :: pack_seam_rows
public :: unpack_seam_rows
public :: SEAM_ROW_COPY, SEAM_ROW_INTERPOLATE, SEAM_ROW_RESTRICT

integer(I4P), parameter :: SEAM_TAG = 4040_I4P !< MPI tag of the cross-rank seam messages.
! Kinds of seam row (issue #52), the values of the intra-realm ghost-map flag (`adam_field_object%update_ghost_local`):
integer(I4P), parameter :: SEAM_ROW_COPY        = 1_I4P !< Same cell size: copy the cell.
integer(I4P), parameter :: SEAM_ROW_INTERPOLATE = 4_I4P !< Coarse cell to a fine ghost: interpolate (`interp_seam_ghost`).
integer(I4P), parameter :: SEAM_ROW_RESTRICT    = 8_I4P !< Fine cells to a coarse ghost: mean of the 2x2x2 fine cells.

contains
   subroutine seam_fill(self, realm, p)
   !< Fill the seam ghosts of `self` (the forest realm `self%realm_index`) from its peer slot `p`: local rows, then
   !< cross-rank rows. Only `self` is written; `realm(:)` (which contains `self`) is read for the peer.
   class(realm_object), intent(inout), target :: self                       !< Realm owning the ghosts.
   class(realm_object), intent(in),    target :: realm(:)                   !< Forest realms.
   integer(I4P),        intent(in)            :: p                          !< Peer slot of `self`.
   integer(I4P)                               :: is                         !< Forest index of `self`.
   integer(I4P)                               :: ip                         !< Peer realm owning the cells.
   integer(I4P)                               :: q                          !< Slot of realm `is` in the peer's lists.
   integer(I4P)                               :: nv                         !< Values per cell.
   integer(I4P)                               :: n_send, n_recv             !< Cross-rank rows of this rank.
   integer(I4P)                               :: s0, r0                     !< First send/receive row of the slot.
   real(R8P), allocatable                     :: send_buf(:), recv_buf(:)   !< Packed values.
   integer(I4P), allocatable                  :: req(:)                     !< MPI requests.
   integer(I4P)                               :: n_req                      !< Posted requests.
   integer(I4P)                               :: r, r_end, rank             !< Row segment bounds and partner rank.
   integer(I4P)                               :: self_recv(2), self_send(2) !< Same-rank row segments (0 if none).
   logical                                    :: timing                     !< Accumulate the phase timers (issue #53).
   real(R8P)                                  :: tw(0:5)                    !< Phase boundaries (MPI_Wtime).

   is = self%realm_index
   ip = self%adam%maps%seam_local_peer_realm(p)
   ! The phase boundaries need no device synchronisation: the FNL seam kernels are synchronous and the packed buffers
   ! cross the host through blocking copies (issue #53).
   timing = self%adam%maps%seam_timing
   if (timing) tw(0) = MPI_Wtime()
   ! local rows: the receiver's own copy kernel
   if (self%adam%maps%seam_local_peer_row_count(p) > 0_I4P) call self%fill_seam_from_peer_forest(peer=realm(ip), p_idx=p)
   if (timing) tw(1) = MPI_Wtime()
   ! cross-rank rows
   q = seam_peer_slot(realm(ip), is)
   n_recv = 0_I4P ; n_send = 0_I4P
   if (allocated(self%adam%maps%seam_mpi_recv_row_count)) n_recv = self%adam%maps%seam_mpi_recv_row_count(p)
   if (allocated(realm(ip)%adam%maps%seam_mpi_send_row_count)) n_send = realm(ip)%adam%maps%seam_mpi_send_row_count(q)
   if (n_recv + n_send == 0_I4P) then
      if (timing) self%adam%maps%seam_wtime(1) = self%adam%maps%seam_wtime(1) + (tw(1) - tw(0))
      return
   endif
   self_recv = 0_I4P ; self_send = 0_I4P
   nv = self%nv
   allocate(req(n_recv + n_send)) ; n_req = 0_I4P
   allocate(recv_buf(nv * n_recv), send_buf(nv * n_send))
   if (timing) tw(2) = MPI_Wtime()
   if (n_send > 0_I4P) then
      associate(r_ip=>realm(ip)) ! nvfortran 26.1 polymorphic array element dispatch workaround (0062a237)
         call r_ip%pack_seam_cells_forest(p_idx=q, buf=send_buf)
      endassociate
   endif
   if (timing) tw(3) = MPI_Wtime()
   ! one message per partner rank: the rows of a slot are grouped by rank
   if (n_recv > 0_I4P) then
      associate(rows=>self%adam%maps%seam_mpi_recv_cell)
      r0 = self%adam%maps%seam_mpi_recv_row_start(p)
      r = 1_I4P
      do while (r <= n_recv)
         rank = rows(r0 + r - 1_I4P, 1) ; r_end = r
         do while (r_end < n_recv)
            if (rows(r0 + r_end, 1) /= rank) exit
            r_end = r_end + 1_I4P
         enddo
         if (rank == mpih%myrank) then
            self_recv = [r, r_end]
         else
            n_req = n_req + 1_I4P
            call MPI_IRECV(recv_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, SEAM_TAG, &
                           MPI_COMM_WORLD, req(n_req), mpih%error)
         endif
         r = r_end + 1_I4P
      enddo
      endassociate
   endif
   if (n_send > 0_I4P) then
      associate(rows=>realm(ip)%adam%maps%seam_mpi_send_cell)
      s0 = realm(ip)%adam%maps%seam_mpi_send_row_start(q)
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
            n_req = n_req + 1_I4P
            call MPI_ISEND(send_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, SEAM_TAG, &
                           MPI_COMM_WORLD, req(n_req), mpih%error)
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
      recv_buf(nv*(self_recv(1)-1_I4P)+1_I4P:nv*self_recv(2)) = send_buf(nv*(self_send(1)-1_I4P)+1_I4P:nv*self_send(2))
   endif
   call MPI_WAITALL(n_req, req, MPI_STATUSES_IGNORE, mpih%error)
   if (timing) tw(4) = MPI_Wtime()
   if (n_recv > 0_I4P) call self%unpack_seam_cells_forest(p_idx=p, buf=recv_buf)
   if (timing) then
      tw(5) = MPI_Wtime()
      self%adam%maps%seam_wtime = self%adam%maps%seam_wtime + (tw(1:5) - tw(0:4))
      self%adam%maps%seam_fills = self%adam%maps%seam_fills + 1_I8P
      self%adam%maps%seam_rows = self%adam%maps%seam_rows + int([n_recv, n_send], I8P)
   endif
   endsubroutine seam_fill

   subroutine seam_fill_all(self, realm)
   !< Fill the seam ghosts of `self` from every peer slot, whatever the cadence (diagnostics on the committed state).
   class(realm_object), intent(inout), target :: self     !< Realm owning the ghosts.
   class(realm_object), intent(in),    target :: realm(:) !< Forest realms.
   integer(I4P)                               :: p        !< Peer slot.

   if (.not.allocated(self%adam%maps%seam_local_peer_realm)) return
   do p=1_I4P, int(size(self%adam%maps%seam_local_peer_realm), I4P)
      call seam_fill(self=self, realm=realm, p=p)
   enddo
   endsubroutine seam_fill_all

   subroutine pack_seam_rows(rows, row_start, row_count, regime, ngc, q, buf)
   !< Pack the ghost values of the send rows `[rank, b, i, j, k, kind, meta]` `row_start..row_start+row_count-1` from
   !< the field `q` (the cells' owner active buffer), `nv` values per row in row order (issue #40, #52): a copy of the
   !< cell, the coarse->fine interpolant of `interp_seam_ghost` (the intra-realm fill; injection copies the anchor) or the
   !< mean of the 2x2x2 fine cells from the base cell (summed in the intra-realm order, then scaled).
   integer(I4P), intent(in)  :: rows(1:,1:)                     !< Send rows.
   integer(I4P), intent(in)  :: row_start                       !< First row.
   integer(I4P), intent(in)  :: row_count                       !< Number of rows.
   integer(I4P), intent(in)  :: regime                          !< Coarse->fine fill regime (`maps%seam_ghost_fill`).
   integer(I4P), intent(in)  :: ngc                             !< Ghost cells number.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field (v, i, j, k, b).
   real(R8P),    intent(out) :: buf(1:)                         !< Packed values.
   integer(I4P)              :: c, row, v, nv                   !< Counters.
   integer(I4P)              :: b, i, j, k                      !< Cell.
   integer(I4P)              :: ic, jc, kc                      !< 2x2x2 counters.
   real(R8P)                 :: total                           !< Restriction sum.

   nv = int(size(q, dim=1), I4P)
   do c=1_I4P, row_count
      row = row_start + c - 1_I4P
      b = rows(row,2) ; i = rows(row,3) ; j = rows(row,4) ; k = rows(row,5)
      select case(rows(row,6))
      case(SEAM_ROW_INTERPOLATE)
         do v=1_I4P, nv
            if (regime == SEAM_FILL_INJECTION) then
               buf(nv*(c-1)+v) = q(v,i,j,k,b)
            else
               buf(nv*(c-1)+v) = interp_seam_ghost(regime=regime, meta=rows(row,7), ngc=ngc, q=q, v=v, &
                                                   i_send=i, j_send=j, k_send=k, b_send=b)
            endif
         enddo
      case(SEAM_ROW_RESTRICT)
         do v=1_I4P, nv
            total = 0._R8P
            do kc=0,1 ; do jc=0,1 ; do ic=0,1
               total = total + q(v,i+ic,j+jc,k+kc,b)
            enddo ; enddo ; enddo
            buf(nv*(c-1)+v) = total * 0.125_R8P
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
