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
!< `seam_fill` moves both kinds for one (realm, peer slot) pair: the local copy, then the owner of the cells packs
!< (`pack_seam_cells_forest`), the buffers travel point to point, the owner of the ghosts unpacks
!< (`unpack_seam_cells_forest`). Each rank calls it for the same pairs in the same order (the peer slots come from the
!< manifest, not from the rows a rank holds); a rank with no rows of a pair returns at once, so the exchange is not
!< collective, but a rank that skipped a pair it has rows for would leave its partners waiting.
use :: adam_mpih_global,  only : mpih
use :: adam_realm_object, only : realm_object
use :: mpi
use :: penf

implicit none
private
public :: seam_fill
public :: seam_fill_all
public :: seam_peer_slot

integer(I4P), parameter :: SEAM_TAG = 4040_I4P !< MPI tag of the cross-rank seam messages.

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

   is = self%realm_index
   ip = self%adam%maps%seam_local_peer_realm(p)
   ! local rows: the receiver's own copy kernel
   if (self%adam%maps%seam_local_peer_row_count(p) > 0_I4P) call self%fill_seam_from_peer_forest(peer=realm(ip), p_idx=p)
   ! cross-rank rows
   q = seam_peer_slot(realm(ip), is)
   n_recv = 0_I4P ; n_send = 0_I4P
   if (allocated(self%adam%maps%seam_mpi_recv_row_count)) n_recv = self%adam%maps%seam_mpi_recv_row_count(p)
   if (allocated(realm(ip)%adam%maps%seam_mpi_send_row_count)) n_send = realm(ip)%adam%maps%seam_mpi_send_row_count(q)
   if (n_recv + n_send == 0_I4P) return
   nv = self%nv
   allocate(req(n_recv + n_send)) ; n_req = 0_I4P
   allocate(recv_buf(nv * n_recv), send_buf(nv * n_send))
   if (n_send > 0_I4P) then
      associate(r_ip=>realm(ip)) ! nvfortran 26.1 polymorphic array element dispatch workaround (0062a237)
         call r_ip%pack_seam_cells_forest(p_idx=q, buf=send_buf)
      endassociate
   endif
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
         n_req = n_req + 1_I4P
         call MPI_IRECV(recv_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, SEAM_TAG, &
                        MPI_COMM_WORLD, req(n_req), mpih%error)
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
         n_req = n_req + 1_I4P
         call MPI_ISEND(send_buf(nv * (r - 1_I4P) + 1_I4P), nv * (r_end - r + 1_I4P), MPI_REAL8, rank, SEAM_TAG, &
                        MPI_COMM_WORLD, req(n_req), mpih%error)
         r = r_end + 1_I4P
      enddo
      endassociate
   endif
   call MPI_WAITALL(n_req, req, MPI_STATUSES_IGNORE, mpih%error)
   if (n_recv > 0_I4P) call self%unpack_seam_cells_forest(p_idx=p, buf=recv_buf)
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
