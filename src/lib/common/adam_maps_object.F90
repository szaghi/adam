!< ADAM, maps class definition
module adam_maps_object
!< ADAM, maps class definition

! ADAM classes, libraries, parameters
use :: adam_tree_object, only : tree_object, tree_iterator_object, NODE_BOUNDARY_CONDITION, &
                                 NODE_LESS_REFINED, NODE_MORE_REFINED, NODE_STANDARD
use :: adam_tree_node_object
use :: adam_parameters
! ADAM singleton objects
use :: adam_mpih_global, only : mpih
use :: adam_grid_object, only : grid_object
! third party modules
use :: penf
! sdk modules
use :: mpi

implicit none
private
public :: maps_object
public :: inter_realm_neighbor_t
public :: FACE_X_MAX, FACE_X_MIN, FACE_Y_MAX, FACE_Y_MIN, FACE_Z_MAX, FACE_Z_MIN
public :: COUPLING_MIRROR, COUPLING_PERIODIC, COUPLING_INTERPOLATE
public :: face_axis_sign

!< Face direction codes (see [[inter_realm_neighbor_t]]).
!<
!< Encoded as small integers rather than a Fortran enum because the values
!< travel through the manifest INI parser (D.3 task of [issue #13]) and need
!< stable, INI-friendly representations.
integer(I4P), parameter :: FACE_X_MAX = 1_I4P  !< +x face (max-x boundary of a block).
integer(I4P), parameter :: FACE_X_MIN = 2_I4P  !< -x face (min-x boundary of a block).
integer(I4P), parameter :: FACE_Y_MAX = 3_I4P  !< +y face.
integer(I4P), parameter :: FACE_Y_MIN = 4_I4P  !< -y face.
integer(I4P), parameter :: FACE_Z_MAX = 5_I4P  !< +z face.
integer(I4P), parameter :: FACE_Z_MIN = 6_I4P  !< -z face.

!< Inter-realm coupling kinds.
!<
!< COUPLING_MIRROR is the only kind exercised by the first Phase D use case
!< (two-realm x-split rmf, [issue #13]). PERIODIC and INTERPOLATE are
!< schema-reserved for follow-up use cases and are not implemented yet.
integer(I4P), parameter :: COUPLING_MIRROR      = 1_I4P  !< Peer cells copied as-is (the trivial pass-through coupling).
integer(I4P), parameter :: COUPLING_PERIODIC    = 2_I4P
                                                         !< Periodic identification across the inter-realm face (reserved, not
                                                         !< implemented).
integer(I4P), parameter :: COUPLING_INTERPOLATE = 3_I4P  !< Interpolate across mismatched grids (reserved, not implemented).

type :: inter_realm_neighbor_t
   !< Description of a single inter-realm face coupling.
   !<
   !< One entry per (my realm, my block, my face) that touches a face on
   !< another realm. Populated by the forest during initialization from the
   !< `[forest.topology]` section of the manifest (D.3 task of [issue #13]).
   !< Consumed by the realm-side `exchange_inter_realm_halos_forest` TBP.
   !<
   !< Architectural note: this type carries only intrinsic-typed components
   !< (R1 of [issue #10] — no derived-type pointer components reachable from
   !< kernels). The `peer_realm` index is used by the consumer to look up the
   !< peer realm in the `realm(:)` array the forest passes through; the
   !< lookup happens on the host side outside any device region (R2).
   integer(I4P) :: my_realm   = 0_I4P
                                      !< Realm index that owns this entry (1-based; redundant with the enclosing maps' realm, kept
                                      !< for self-describing dumps).
   integer(I4P) :: my_block   = 0_I4P !< Block index in MY realm whose face touches the peer.
   integer(I4P) :: my_face    = 0_I4P !< Face code on MY block (FACE_X_MAX..FACE_Z_MIN).
   integer(I4P) :: peer_realm = 0_I4P !< Realm index of the peer.
   integer(I4P) :: peer_block = 0_I4P !< Block index in the peer realm whose face touches mine.
   integer(I4P) :: peer_face  = 0_I4P !< Face code on the peer block (typically opposite to my_face).
   integer(I4P) :: coupling   = COUPLING_MIRROR !< Coupling kind (COUPLING_MIRROR | COUPLING_PERIODIC | COUPLING_INTERPOLATE).
endtype inter_realm_neighbor_t

type :: maps_object
   !< Maps class definition
   logical :: is_initialized_=.false. !< Flag: maps have been initialized.
   ! local maps
   integer(I8P), allocatable :: local_map(:,:)            !< Local map, list block index changes of my nodes.
   integer(I8P), allocatable :: local_map_ghost(:,:)      !< Local map for ghost cells updating [fec_number, 4].
   integer(I8P), allocatable :: local_map_ghost_cell(:,:) !< Local map ghost cells update, cells order.
   integer(I8P), allocatable :: local_map_bc_face(:,:)    !< Local map for face BC ghost cells.
   integer(I8P), allocatable :: local_map_bc_edge(:,:)    !< Local map for edge BC ghost cells.
   integer(I8P), allocatable :: local_map_bc_corner(:,:)  !< Local map for corner BC ghost cells.
   integer(I8P), allocatable :: local_map_bc_crown(:,:,:) !< Local map for face BC ghost cells, crown order.
   ! MPI send/recv maps
   integer(I4P)              :: my_nodes_number=0_I4P     !< Number of my nodes, keep_nodes_number + recv_nodes_number.
   integer(I4P)              :: send_nodes_number=0_I4P   !< Number of nodes to be sent.
   integer(I4P)              :: recv_nodes_number=0_I4P   !< Number of nodes to be received.
   integer(I4P)              :: keep_nodes_number=0_I4P   !< Number of nodes to be keept.
   integer(I4P)              :: inner_blocks_number=0_I4P !< Number of inner blocks where I need fecs.
   integer(I4P), allocatable :: inner_outer_block_map(:)  !< Inner/outer blocks map.
   integer(I4P), allocatable :: comm_map_n_send(:)        !< Communication map, number of blocks to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv(:)        !< Communication map, number of blocks to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr(:)      !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr(:)      !< Communication map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable :: comm_map_send(:)          !< Communication map, blocks to send [sum(comm_map_n_send)].
   integer(I8P), allocatable :: comm_map_recv(:)          !< Communication map, blocks to receive [sum(comm_map_n_recv)].
   ! MPI send/recv ghost cells maps
   integer(I4P), allocatable :: comm_map_n_send_ghost(:)      !< Communication map, number of ghost celss to send [procs_number].
   integer(I4P), allocatable :: comm_map_n_recv_ghost(:)      !< Communication map, number of ghost celss to recv [procs_number].
   integer(I4P), allocatable :: comm_map_send_ptr_ghost(:)    !< Communication map, pointers in list to send [procs_number+1].
   integer(I4P), allocatable :: comm_map_recv_ptr_ghost(:)    !< Communication map, pointers in list to recv [procs_number+1].
   integer(I8P), allocatable :: comm_map_send_ghost(:,:)      !< Communication map, send ghost cells, fec order.
   integer(I8P), allocatable :: comm_map_send_ghost_cell(:,:) !< Communication map, send ghost cells, cells order.
   integer(I8P), allocatable :: comm_map_recv_ghost(:,:)      !< Communication map, recv ghost cells, fec order.
   integer(I8P), allocatable :: comm_map_recv_ghost_cell(:,:) !< Communication map, recv ghost cells, cells order.
   ! MPI send/recv ghost buffers
   real(R8P), allocatable :: send_buffer_ghost(:) !< Send buffer of ghost cells.
   real(R8P), allocatable :: recv_buffer_ghost(:) !< Receive buffer of ghost cells.
   ! Inter-realm coupling (Phase D of issue #10, design in issue #13).
   ! Unallocated for single-realm runs (the no-op case). Populated by the
   ! forest during initialization from the manifest's [forest.topology]
   ! section. Consumed by realm_object%exchange_inter_realm_halos_forest.
   type(inter_realm_neighbor_t), allocatable :: inter_realm_neighbors(:)
                                                                         !< Face couplings to other realms; unallocated when this
                                                                         !< realm has no inter-realm neighbors.
   !
   ! Inter-realm ghost-cell map (Phase A of issue #13 — seam comm-map).
   !
   ! One row per ghost cell that lives in THIS realm and pulls its value
   ! from another realm. Mirrors `local_map_ghost_cell` (the intra-realm
   ! ghost map) with `peer_realm` prepended as column 1.
   !
   ! Layout: `inter_realm_ghost_cell(c, 1:10)` =
   !   [peer_realm, b_send, b_recv, i_send, j_send, k_send,
   !    i_recv, j_recv, k_recv, one_or_eight]
   !
   ! where:
   !   * `peer_realm`   — 1-based index of the realm to read from;
   !   * `b_send`       — block index in the peer realm;
   !   * `b_recv`       — block index in THIS realm;
   !   * `i/j/k_send`   — INTERIOR cell coordinates in the peer block
   !                      (1..ni_peer / 1..nj_peer / 1..nk_peer);
   !   * `i/j/k_recv`   — GHOST cell coordinates in self's block
   !                      (in [1-ngc..ni+ngc] etc.);
   !   * `one_or_eight` — 1 for same-resolution mirror copy, reserved 8
   !                      for the 2:1-coarser case (Phase A v1 only
   !                      populates same-resolution = 1; AMR jumps are
   !                      a follow-up).
   !
   ! Unallocated when there are no inter-realm seams. Populated by
   ! `forest_object%populate_inter_realm_topology` (D.4a + Phase A) by
   ! enumerating every ghost cell in every self-seam-block's ghost
   ! region and resolving the matching peer (realm, block, interior
   ! cell) — including the corner / edge ghosts that the previous
   ! geometric `exchange_inter_realm_halos_forest` missed (defect B of
   ! issue #13).
   !
   ! Consumer: a flat-loop exchange in
   ! `realm_object%exchange_inter_realm_halos_forest` overrides
   ! (currently `prism_cpu_object`), replacing the per-call geometric
   ! `find_peer_block` + `copy_peer_face` pair with an indexed read.
   integer(I4P), allocatable :: inter_realm_ghost_cell(:,:) !< Per-cell inter-realm ghost map; layout per the comment above.
   !
   ! Inter-realm seam → flux register lookup (Phase A step 4 of issue #13).
   !
   ! Maps a (block, face_1_6) pair on THIS realm to the index of the
   ! corresponding entry in `forest_flux_register%face(:)`. The
   ! face_1_6 axis is the BC-routine numbering (FEC_1_6_ARRAY space:
   ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z) so PRISM's residual routines —
   ! which already compute `fec_1_6` in their BC loops — can look up the
   ! register index directly without a second translation.
   !
   ! Layout: `inter_realm_face_register_index(1:nb, 1:6)` — encodes both
   ! the register index AND this realm's side (coarse vs fine) in one
   ! signed integer:
   !
   !   0       : not a seam face for this (block, face) pair
   !   +idx    : seam face; this realm is the COARSE side of register entry idx
   !   -idx    : seam face; this realm is the FINE   side of register entry idx
   !
   ! The consumer recovers the register index via `abs()` and the
   ! coarse/fine role via the sign — no realm-identity field is needed
   ! on the realm object. The sign mirrors the manifest convention
   ! `pair%realm_a = coarse, pair%realm_b = fine` (currently same
   ! resolution; the labels become load-bearing when true AMR coarse-fine
   ! adjacency lands).
   !
   ! Populated by `forest_object%populate_inter_realm_topology` (Phase A
   ! step 4) at the same time as the per-cell ghost map and the BC_SEAM
   ! override; ALL three derive from the same manifest face-pair walk.
   !
   ! Consumer: PRISM's `compute_residuals_fv_centered` consults the table
   ! after reconstructing each face's flux. If non-zero, the routine
   ! calls `forest_flux_register%accumulate_coarse_flux` (sign > 0) or
   ! `accumulate_fine_flux` (sign < 0) on the indexed entry.
   integer(I4P), allocatable :: inter_realm_face_register_index(:,:)
                                                                     !< (block, face_1_6) → signed flux register face index; 0 = not
                                                                     !< a seam face, +idx = coarse side, -idx = fine side.
   !
   ! Seam ghost-cell maps — agnostic-dummy redesign (issue #13 follow-up).
   !
   ! These supersede `inter_realm_ghost_cell` (which remains allocated during
   ! the migration). The redesign separates the same-rank fast path from the
   ! cross-rank MPI path — mirroring the intra-realm split between
   ! `local_map_ghost_cell` and `comm_map_send/recv_ghost_cell` — so the
   ! forest's Phase 2 seam fill can dispatch identically to intra-realm ghost
   ! fill: pack a buffer from the peer's active `q`, unpack into THIS realm's
   ! ghost slots, with no knowledge of the peer's integrator buffer layout.
   !
   ! Layout: `seam_local_map_ghost_cell(c, 1:9)` =
   !   [peer_realm, b_send, b_recv, i_send, j_send, k_send, i_recv, j_recv, k_recv]
   !
   ! Same columns as `inter_realm_ghost_cell` minus `one_or_eight` (Phase A
   ! is same-resolution only; restore when AMR coarse-fine seams land). Rows
   ! are SORTED BY peer_realm so per-peer row ranges can be extracted in O(1)
   ! via the seam_local_peer_* index arrays below.
   integer(I4P), allocatable :: seam_local_map_ghost_cell(:,:)
                                                                     !< Per-cell same-rank inter-realm ghost map; layout per the comment
                                                                     !< above; rows sorted by peer_realm.
   integer(I4P), allocatable :: seam_local_peer_realm(:)     !< Peer realm index for each row range (one entry per distinct peer).
   integer(I4P), allocatable :: seam_local_peer_row_start(:) !< First row index in seam_local_map_ghost_cell for each peer.
   integer(I4P), allocatable :: seam_local_peer_row_count(:) !< Row count in seam_local_map_ghost_cell for each peer.
   !
   ! Per-peer pack/unpack buffers. Shape: `(nv * max_rows_per_peer, n_peers)`.
   ! Each peer's pack-then-unpack roundtrip uses an independent column — the
   ! forest may process peers serially OR in parallel without aliasing.
   ! Allocated at topology init; reused every stage (no per-step allocation).
   real(R8P), allocatable :: seam_local_send_buf(:,:) !< Per-peer pack buffer for same-rank seams.
   real(R8P), allocatable :: seam_local_recv_buf(:,:) !< Per-peer unpack buffer for same-rank seams.
   !
   ! Seam MPI maps and buffers — Phase A: DECLARED but NOT POPULATED.
   ! Reserved for the cross-rank seam exchange (Phase B). If any of these
   ! become allocated the forest Phase 2 loop must dispatch to
   ! `update_ghost_seam_mpi`; until that path lands the forest error_stops
   ! on any allocated `seam_comm_map_send_ghost_cell`.
   integer(I4P), allocatable :: seam_comm_map_send_ghost_cell(:,:)
   integer(I4P), allocatable :: seam_comm_map_recv_ghost_cell(:,:)
   integer(I4P), allocatable :: seam_comm_map_send_ptr_ghost(:)
   integer(I4P), allocatable :: seam_comm_map_recv_ptr_ghost(:)
   real(R8P),    allocatable :: seam_mpi_send_buf(:)
   real(R8P),    allocatable :: seam_mpi_recv_buf(:)
   contains
      ! public methods
      procedure, pass(self) :: blocks_reorder             !< Reorder blocks indexes in field.
      procedure, pass(self) :: get_block_layout          !< Return block Morton codes and/or coordinates.
      procedure, pass(self) :: initialize                 !< Initialize maps.
      procedure, pass(self) :: make_comm_local_maps       !< Make communication/local maps.
      procedure, pass(self) :: make_comm_local_maps_ghost !< Make communication/local maps of ghost cells.
      procedure, pass(self) :: make_local_maps_bc         !< Make local maps of boundary conditions.
      procedure, pass(self) :: mpi_gather_nodes_data      !< Gather nodes data between MPI processes.
      procedure, pass(self) :: save_local_map             !< Save local map.
      ! private methods
      procedure, pass(self), private :: alloc_comm_local_maps_ghost    !< Allocate communication/local maps of ghost cells.
      procedure, pass(self), private :: count_bc_numbers               !< Count BC numbers.
      procedure, pass(self), private :: ijk_mmd                        !< Return IJK min/max/delta.
      procedure, pass(self), private :: ijk_mmd_ghost                  !< Return IJK min/max/delta, ghost maps.
      procedure, pass(self), private :: make_comm_map_recv_ghost_cell  !< Make communication recv map of ghost cells in cell order.
      procedure, pass(self), private :: make_comm_map_send_ghost_cell  !< Make communication send map of ghost cells in cell order.
      procedure, pass(self), private :: make_local_map_ghost_cell      !< Make local map of ghost cells in cell order.
      procedure, pass(self), private :: populate_comm_local_maps_ghost !< Populate communication/local maps of ghost cells.
endtype maps_object

contains
   ! public module-scope helpers
   pure subroutine face_axis_sign(face_code, axis, sgn)
   !< Translate FACE_X_MAX / FACE_X_MIN / ... into (axis 1..3, sign ±1).
   !<
   !< Canonical decoder for the FACE_* face-direction constants. Used by
   !< the forest's seam-topology builder, by inter-realm ghost-fill paths,
   !< and by realm-side reflux TBPs that need to map a `face_code` onto
   !< the corresponding axis and sign.
   integer(I4P), intent(in)  :: face_code !< Face code (FACE_X_MAX..FACE_Z_MIN).
   integer(I4P), intent(out) :: axis      !< 1=x, 2=y, 3=z. 0 on malformed input.
   integer(I4P), intent(out) :: sgn       !< +1 if MAX, -1 if MIN. 0 on malformed input.

   select case (face_code)
   case (FACE_X_MAX); axis = 1_I4P; sgn = +1_I4P
   case (FACE_X_MIN); axis = 1_I4P; sgn = -1_I4P
   case (FACE_Y_MAX); axis = 2_I4P; sgn = +1_I4P
   case (FACE_Y_MIN); axis = 2_I4P; sgn = -1_I4P
   case (FACE_Z_MAX); axis = 3_I4P; sgn = +1_I4P
   case (FACE_Z_MIN); axis = 3_I4P; sgn = -1_I4P
   case default;      axis = 0_I4P; sgn = 0_I4P
   end select
   endsubroutine face_axis_sign

   ! public methods
   subroutine blocks_reorder(self, tree)
   !< Reorder blocks indexes in field.
   class(maps_object), intent(inout) :: self                !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   type(tree_node_object), pointer   :: node_ptr            !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I4P)                      :: outer_blocks_number !< Number of outer blocks where I need fecs.
   integer(I4P)                      :: inner_blocks_number !< Number of inner blocks where I need fecs.
   logical                           :: is_inner_block      !< Flag to check if a neighbor block is inner or not.
   integer(I8P), allocatable         :: neighbor(:)         !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh               !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type       !< Neighbors type.
   integer(I4P), allocatable         :: inner_block_map(:)  !< Inner blocks map.
   integer(I4P), allocatable         :: outer_block_map(:)  !< Outer blocks map.
   integer(I4P)                      :: fec, n              !< Counter.

   outer_blocks_number = 0
   inner_blocks_number = 0
   allocate(inner_block_map(self%my_nodes_number))
   allocate(outer_block_map(self%my_nodes_number))
   if (allocated(self%inner_outer_block_map)) deallocate(self%inner_outer_block_map)
   allocate(self%inner_outer_block_map(self%my_nodes_number))
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      if (mpih%myrank == node_ptr%myrank) then
         is_inner_block = .true.
         do fec=1, 26
            neighbor      = node_ptr%neighbor(fec)%codes
            neighbor_type = node_ptr%neighbor(fec)%ntype
            if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
               do n=1, size(neighbor, dim=1)
                  neigh => tree%node(code=neighbor(n))
                  if (mpih%myrank /= neigh%myrank) is_inner_block = .false.
               enddo
            endif
         enddo
         if (is_inner_block) then
            inner_blocks_number = inner_blocks_number + 1
            inner_block_map(inner_blocks_number) = node_ptr%block_index
            node_ptr%block_index_new = inner_blocks_number
         else
            outer_blocks_number = outer_blocks_number + 1
            outer_block_map(outer_blocks_number) = node_ptr%block_index
            node_ptr%block_index_new = -outer_blocks_number
         endif
      endif
   enddo
   self%inner_blocks_number = inner_blocks_number
   self%inner_outer_block_map(1:inner_blocks_number ) = inner_block_map(1:inner_blocks_number)
   self%inner_outer_block_map(inner_blocks_number+1:) = outer_block_map(1:outer_blocks_number)
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      if (mpih%myrank == node_ptr%myrank) then
         if (node_ptr%block_index_new < 0) then
            node_ptr%block_index = -node_ptr%block_index_new + inner_blocks_number
         else
            node_ptr%block_index =  node_ptr%block_index_new
         endif
      endif
   enddo
   call self%mpi_gather_nodes_data(tree=tree, node_member='block_index')
   endsubroutine blocks_reorder

   subroutine get_block_layout(self, tree, code, coordinates)
   !< Return block Morton codes and/or coordinates from the underlying tree.
   !< Provides a single-level accessor so callers do not traverse tree directly.
   class(maps_object), intent(in)            :: self           !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   integer(I8P),       intent(out), optional :: code(:)        !< Block Morton codes [blocks_number].
   integer(I4P),       intent(out), optional :: coordinates(:,:) !< Block coordinates IJKL [4,blocks_number].

   if (present(code))        code        = tree%block_code
   if (present(coordinates)) coordinates = tree%block_coordinates
   endsubroutine get_block_layout

   subroutine initialize(self, tree, verbose)
   !< Initialize maps.
   class(maps_object), intent(inout)        :: self     !< The maps.
   type(tree_object),  intent(inout) :: tree !< Tree (sibling realm component, threaded in).
   logical,            intent(in), optional :: verbose  !< Trigger verbose output.
   logical                                  :: verbose_ !< Trigger verbose output, local variable.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih%print_message('maps_object%initialize start')
   if (.not.tree%is_initialized_) &
      call mpih%error_stop(': maps_object%initialize: tree is not initialized')
   allocate(self%comm_map_n_send(0:mpih%procs_number-1))
   allocate(self%comm_map_n_recv(0:mpih%procs_number-1))
   allocate(self%comm_map_send_ptr(0:mpih%procs_number))
   allocate(self%comm_map_recv_ptr(0:mpih%procs_number))
   self%comm_map_n_send = 0_I4P
   self%comm_map_n_recv = 0_I4P
   self%comm_map_send_ptr = 0_I4P
   self%comm_map_recv_ptr = 0_I4P
   allocate(self%comm_map_n_send_ghost(0:mpih%procs_number-1))
   allocate(self%comm_map_n_recv_ghost(0:mpih%procs_number-1))
   allocate(self%comm_map_send_ptr_ghost(0:mpih%procs_number))
   allocate(self%comm_map_recv_ptr_ghost(0:mpih%procs_number))
   call self%make_comm_local_maps(tree=tree)
   self%is_initialized_ = .true.
   if (verbose_) call mpih%print_message('maps_object%initialize finish')
   endsubroutine initialize

   subroutine make_comm_local_maps(self, tree)
   !< Make communication/local maps.
   !<```
   !<   comm_map_send     = [ 17, 511, 92, 3, 54, 56, 11, 12...] (block index).
   !<                          |   |       |  ||       |
   !<   comm_map_send_ptr = [  0,  1,  3,  4,  4, 6, (8)]        (pointer to comm_map_send)
   !<   comm_map_recv     = [ 23, 4, 51, 69, 145, 2, 72, 16, 6]  (block index).
   !<                         |       |  ||          |       |
   !<   comm_map_recv_prt = [ 0,  2,  3,  3,  6, 8, (9)]         (pointer to comm_map_recv)
   !<```
   class(maps_object), intent(inout) :: self                 !< The maps.
   type(tree_object),  intent(inout) :: tree !< Tree (sibling realm component, threaded in).
   type(tree_node_object), pointer   :: node_ptr             !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I8P), allocatable         :: codes_sorted(:)      !< List of (sorted) codes.
   integer(I4P), allocatable         :: comm_map_send_ctr(:) !< Communication map, counters in list to send [procs_number+1].
   integer(I4P), allocatable         :: comm_map_recv_ctr(:) !< Communication map, counters in list to recv [procs_number+1].
   integer(I8P)                      :: my_nodes_number      !< Number of my nodes.
   integer(I4P)                      :: n_send               !< Number of nodes that I have to send.
   integer(I4P)                      :: n_recv               !< Number of nodes that I have to receive.
   integer(I4P)                      :: n_keep               !< Number of nodes that I have to keep.
   integer(I8P)                      :: c, l                 !< Counter.
   integer(I4P)                      :: p                    !< Counter.

   call mpih%print_message('maps_object%make_comm_local_maps start')
   ! initialize communication maps
   self%comm_map_n_send = 0_I4P
   self%comm_map_n_recv = 0_I4P
   if (allocated(self%comm_map_send)) deallocate(self%comm_map_send)
   if (allocated(self%comm_map_recv)) deallocate(self%comm_map_recv)
   if (allocated(self%local_map    )) deallocate(self%local_map    )

   ! compute the number of blocks to send/receive
   my_nodes_number = 0_I8P
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      if (node_ptr%myrank==mpih%myrank) my_nodes_number = my_nodes_number + 1_I8P
      if     (is_node_to_send(n=node_ptr)) then
         ! I have this node that must be sent to node_ptr%myrank_new
         self%comm_map_n_send(node_ptr%myrank_new) = self%comm_map_n_send(node_ptr%myrank_new) + 1
      elseif (is_node_to_receive(n=node_ptr)) then
         ! node_ptr%rank has this node that must be sent to me
         self%comm_map_n_recv(node_ptr%myrank) = self%comm_map_n_recv(node_ptr%myrank) + 1
      endif
   enddo
   n_send = sum(self%comm_map_n_send, dim=1)
   n_recv = sum(self%comm_map_n_recv, dim=1)
   n_keep = my_nodes_number - n_send
   self%my_nodes_number = n_keep + n_recv
   self%recv_nodes_number = n_recv
   self%send_nodes_number = n_send
   if (n_keep > 0_I4P) allocate(self%local_map(n_keep,2))

   ! allocate communications maps
   if (n_send>0_I4P) then
      allocate(self%comm_map_send(n_send))
      self%comm_map_send = 0_I8P
   endif
   if (n_recv>0_I4P) then
      allocate(self%comm_map_recv(n_recv))
      self%comm_map_recv = 0_I8P
   endif

   ! compute communication maps pointers/counters
   self%comm_map_send_ptr = 0_I4P
   self%comm_map_recv_ptr = 0_I4P
   do p=1, mpih%procs_number
      self%comm_map_send_ptr(p) = self%comm_map_send_ptr(p-1) + self%comm_map_n_send(p-1)
      self%comm_map_recv_ptr(p) = self%comm_map_recv_ptr(p-1) + self%comm_map_n_recv(p-1)
   enddo
   comm_map_send_ctr = self%comm_map_send_ptr
   comm_map_recv_ctr = self%comm_map_recv_ptr

   ! populate communication maps
   codes_sorted = tree%codes() ! sorted list of Morton codes
   ! send map
   if (n_send>0_I4P) then
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => tree%node(code=codes_sorted(c))
         if (is_node_to_send(n=node_ptr)) then
            self%comm_map_send(comm_map_send_ctr(node_ptr%myrank_new)+1) = node_ptr%block_index
            comm_map_send_ctr(node_ptr%myrank_new) = comm_map_send_ctr(node_ptr%myrank_new) + 1
         endif
      enddo
   endif
   ! receive map
   if (n_recv>0_I4P) then
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => tree%node(code=codes_sorted(c))
         if (is_node_to_receive(n=node_ptr)) then
            self%comm_map_recv(comm_map_recv_ctr(node_ptr%myrank)+1) = node_ptr%block_index_new
            comm_map_recv_ctr(node_ptr%myrank) = comm_map_recv_ctr(node_ptr%myrank) + 1
         endif
      enddo
   endif
   ! local maps
   if (n_keep > 0_I4P) then
      l = 0
      do c=1, size(codes_sorted, dim=1) ! sort blocks in Morton order
         node_ptr => tree%node(code=codes_sorted(c))
         if (is_node_to_keep(n=node_ptr)) then
            l = l + 1
            self%local_map(l,1) = node_ptr%block_index_new
            self%local_map(l,2) = node_ptr%block_index
         endif
      enddo
   endif

   ! update tree blocks coordinates
   call tree%update_blocks_coordinates(my_nodes_number=self%my_nodes_number)
   call mpih%print_message('maps_object%make_comm_local_maps finish')
   contains
      function is_node_to_keep(n)
      !< Check if node `n` must be kept.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_keep !< Check result.

      is_node_to_keep = ((mpih%myrank == n%myrank).and.(n%myrank == n%myrank_new))
      endfunction is_node_to_keep

      function is_node_to_send(n)
      !< Check if node `n` must be sent.
      type(tree_node_object), intent(in), pointer :: n               !< Pointer to current node.
      logical                                     :: is_node_to_send !< Check result.

      is_node_to_send = ((mpih%myrank == n%myrank).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_send

      function is_node_to_receive(n)
      !< Check if node `n` must be received.
      type(tree_node_object), intent(in), pointer :: n                  !< Pointer to current node.
      logical                                     :: is_node_to_receive !< Check result.

      is_node_to_receive = ((mpih%myrank == n%myrank_new).and.(n%myrank /= n%myrank_new))
      endfunction is_node_to_receive
   endsubroutine make_comm_local_maps

   subroutine make_comm_local_maps_ghost(self, tree, grid, nv)
   !< Make communication/local maps of ghost cells.
   class(maps_object), intent(inout) :: self !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P),       intent(in)    :: nv   !< Number of field variables.

   call mpih%print_message('maps_object%make_comm_local_maps_ghost start')
   call self%alloc_comm_local_maps_ghost(tree=tree, grid=grid, nv=nv)
   call self%populate_comm_local_maps_ghost(tree=tree, grid=grid)
   call self%make_local_map_ghost_cell
   call self%make_comm_map_recv_ghost_cell(nv=nv)
   call self%make_comm_map_send_ghost_cell(nv=nv)
   call mpih%print_message('maps_object%make_comm_local_maps_ghost finish')
   endsubroutine make_comm_local_maps_ghost

   subroutine mpi_gather_nodes_data(self, tree, node_member)
   !< Gather nodes data status between MPI processes.
   class(maps_object), intent(inout) :: self           !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   character(*),       intent(in)    :: node_member    !< Node member to be shared.
   type(tree_node_object), pointer   :: node_ptr       !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I8P), allocatable         :: send_buffer(:) !< Send buffer nodes data.
   integer(I8P), allocatable         :: recv_buffer(:) !< Recv buffer nodes data.
   integer(I4P), allocatable         :: recv_count(:)  !< Number of nodes that are received from each process.
   integer(I4P), allocatable         :: disp_count(:)  !< Displacement of nodes that are received from each process.
   integer(I8P)                      :: n, p           !< Counter.

   allocate(send_buffer(self%my_nodes_number   * 2)) ! [Morton code, refinement_needed]
   allocate(recv_buffer(tree%nodes_number * 2)) ! [Morton code, refinement_needed]
   allocate(recv_count(0:mpih%procs_number - 1))
   allocate(disp_count(0:mpih%procs_number - 1))
   ! populating receive counters and send buffer
   recv_count = 0_I4P
   n = 0_I8P
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      recv_count(node_ptr%myrank) = recv_count(node_ptr%myrank) + 2
      if (mpih%myrank == node_ptr%myrank) then
         n = n + 1
         send_buffer(n) = node_ptr%code
         n = n + 1
         select case(trim(node_member))
         case('refinement_needed')
            send_buffer(n) = node_ptr%refinement_needed
         case('block_index')
            send_buffer(n) = node_ptr%block_index
         endselect
      endif
   enddo
   ! computing displacement counts
   disp_count = 0_I4P
   do p=1, mpih%procs_number - 1
      disp_count(p) = disp_count(p-1) + recv_count(p-1)
   enddo

   call MPI_ALLGATHERV(send_buffer, self%my_nodes_number * 2, MPI_INTEGER8, &
                       recv_buffer, recv_count, disp_count, MPI_INTEGER8, MPI_COMM_WORLD, mpih%error)

   ! update nodes data
   do n=1, tree%nodes_number*2, 2
      node_ptr => tree%node(code=recv_buffer(n))
      select case(trim(node_member))
      case('refinement_needed')
         node_ptr%refinement_needed = int(recv_buffer(n+1), I4P)
      case('block_index')
         node_ptr%block_index = recv_buffer(n+1)
      endselect
   enddo
   endsubroutine mpi_gather_nodes_data

   subroutine save_local_map(self, file_name)
   !< Save local map.
   class(maps_object), intent(in)  :: self      !< The maps.
   character(*),       intent(in)  :: file_name !< Output file name.
   integer(I4P)                    :: file_unit !< Output file unit.
   integer(I4P)                    :: i1        !< Counter.

   if (allocated(self%local_map_ghost)) then
      open(newunit=file_unit, file=trim(adjustl(file_name)), form='FORMATTED')
      do i1=lbound(self%local_map_ghost, dim=1), ubound(self%local_map_ghost, dim=1)
         write(unit=file_unit, fmt='(13(I6,1X))') self%local_map_ghost(i1, :)
      enddo
      close(file_unit)
   endif
   endsubroutine save_local_map

   subroutine make_local_maps_bc(self, tree, grid, verbose)
   !< Make local maps of boundary conditions.
   class(maps_object), intent(inout)        :: self                   !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   logical,            intent(in), optional :: verbose                !< Flag to activate verbose mode.
   logical                                  :: verbose_               !< Flag to activate verbose mode, local var.
   type(tree_node_object), pointer          :: node_ptr               !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I4P)                             :: neighbor_type          !< Neighbors type.
   integer(I4P)                             :: neighbor_bc_fec        !< Neighbors fec for BC.
   integer(I4P)                             :: fec                    !< Counter.
   integer(I8P)                             :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
   integer(I4P)                             :: fec_bc_faces_number    !< BC faces number.
   integer(I4P)                             :: fec_bc_edges_number    !< BC edges number.
   integer(I4P)                             :: fec_bc_corners_number  !< BC corners number.
   integer(I4P)                             :: fec_bc_type            !< BC type.
   integer(I4P), allocatable                :: c_crown(:), c          !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih%print_message('maps_object%make_local_maps_bc start')

   ! find faces/edges/cornes BC
   if (allocated(self%local_map_bc_face))   deallocate(self%local_map_bc_face)
   if (allocated(self%local_map_bc_edge))   deallocate(self%local_map_bc_edge)
   if (allocated(self%local_map_bc_corner)) deallocate(self%local_map_bc_corner)
   call self%count_bc_numbers(tree=tree, fec_bc_faces_number=fec_bc_faces_number, &
                              fec_bc_edges_number=fec_bc_edges_number, &
                              fec_bc_corners_number=fec_bc_corners_number)
   ! allocate BC arrays
   if (fec_bc_faces_number > 0.or.fec_bc_edges_number > 0.or.fec_bc_corners_number > 0) then
      if (fec_bc_faces_number > 0)   allocate(self%local_map_bc_face(fec_bc_faces_number, 12))
      if (fec_bc_edges_number > 0)   allocate(self%local_map_bc_edge(fec_bc_edges_number, 12))
      if (fec_bc_corners_number > 0) allocate(self%local_map_bc_corner(fec_bc_corners_number, 12))
      fec_bc_faces_number = 0_I4P
      fec_bc_edges_number = 0_I4P
      fec_bc_corners_number = 0_I4P
      iter%b = 1_I4P ; iter%p => null()
      do while(tree%loop(iter, node_ptr=node_ptr))
         if (node_ptr%myrank==mpih%myrank) then
            do fec=1, 26
               neighbor_type   = node_ptr%neighbor(fec)%ntype
               neighbor_bc_fec = node_ptr%neighbor(fec)%bc_fec
               if (neighbor_type == NODE_BOUNDARY_CONDITION) then
                  fec_bc_type = grid%fec_bc_type(fec=neighbor_bc_fec)
                  ijk_min_max_delta = self%ijk_mmd(grid=grid, fec=fec, neighbor_bc_fec=neighbor_bc_fec)
                  if (fec<=6) then
                     fec_bc_faces_number = fec_bc_faces_number + 1_I4P
                     call set_local_map_bc(fec=fec_bc_faces_number,                  &
                                           block_index=node_ptr%block_index,         &
                                           neighbor_bc_fec=int(neighbor_bc_fec,I8P), &
                                           ijk_min_max_delta=ijk_min_max_delta,      &
                                           fec_bc_type=int(fec_bc_type,I8P),         &
                                           local_map_bc=self%local_map_bc_face)
                  endif
                  if (fec>=7.and.fec<=18) then
                     fec_bc_edges_number = fec_bc_edges_number + 1_I4P
                     call set_local_map_bc(fec=fec_bc_edges_number,                  &
                                           block_index=node_ptr%block_index,         &
                                           neighbor_bc_fec=int(neighbor_bc_fec,I8P), &
                                           ijk_min_max_delta=ijk_min_max_delta,      &
                                           fec_bc_type=int(fec_bc_type,I8P),         &
                                           local_map_bc=self%local_map_bc_edge)
                  endif
                  if (fec>=19) then
                     fec_bc_corners_number = fec_bc_corners_number + 1_I4P
                     call set_local_map_bc(fec=fec_bc_corners_number,                &
                                           block_index=node_ptr%block_index,         &
                                           neighbor_bc_fec=int(neighbor_bc_fec,I8P), &
                                           ijk_min_max_delta=ijk_min_max_delta,      &
                                           fec_bc_type=int(fec_bc_type,I8P),         &
                                           local_map_bc=self%local_map_bc_corner)
                  endif
               endif
            enddo
         endif
      enddo
   endif
   ! put faces/edges/corners arrays into crowns one
   if (allocated(self%local_map_bc_crown)) deallocate(self%local_map_bc_crown)
   if (allocated(self%local_map_bc_face).or.allocated(self%local_map_bc_edge).or.allocated(self%local_map_bc_corner)) then
      c = 0
      if (allocated(self%local_map_bc_face  )) c = c + bc_cells_number(local_map_bc=self%local_map_bc_face  )
      if (allocated(self%local_map_bc_edge  )) c = c + bc_cells_number(local_map_bc=self%local_map_bc_edge  )
      if (allocated(self%local_map_bc_corner)) c = c + bc_cells_number(local_map_bc=self%local_map_bc_corner)
      allocate(self%local_map_bc_crown(1:c,1:9,1:grid%ngc))
      allocate(c_crown(1:grid%ngc))
      self%local_map_bc_crown = -1
      c_crown = 1
      if (allocated(self%local_map_bc_face  )) call populate_local_map_bc_crown(local_map_bc=self%local_map_bc_face,        &
                                                                                local_map_bc_crown=self%local_map_bc_crown, &
                                                                                c_crown=c_crown)
      if (allocated(self%local_map_bc_edge  )) call populate_local_map_bc_crown(local_map_bc=self%local_map_bc_edge,        &
                                                                                local_map_bc_crown=self%local_map_bc_crown, &
                                                                                c_crown=c_crown)
      if (allocated(self%local_map_bc_corner)) call populate_local_map_bc_crown(local_map_bc=self%local_map_bc_corner,      &
                                                                                local_map_bc_crown=self%local_map_bc_crown, &
                                                                                c_crown=c_crown)
      deallocate(c_crown)
   endif
   if (verbose_) call mpih%print_message('tree_object%make_local_maps_bc finish')
   endsubroutine make_local_maps_bc

   ! private methods
   subroutine alloc_comm_local_maps_ghost(self, tree, grid, nv)
   !< Allocate communication/local maps of ghost cells, `local_map_ghost, comm_map_send_ghost, comm_map_recv_ghost`.
   class(maps_object), intent(inout) :: self             !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P),       intent(in)    :: nv               !< Number of field variables.
   type(tree_node_object), pointer   :: node_ptr         !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I8P), allocatable         :: neighbor(:)      !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh            !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type    !< Neighbors type.
   integer(I4P)                      :: my_fec_number    !< Number of faces/edges/corners ghost cells locally exchanged.
   integer(I4P)                      :: recv_fec_number  !< Number of faces/edges/corners ghost cells externally received.
   integer(I4P)                      :: send_fec_number  !< Number of faces/edges/corners ghost cells externally sent.
   integer(I4P)                      :: fec, n, p        !< Counter.
   integer(I4P)                      :: weight_reduction !< Neighbor weight reduction for send/recv at different level.

   if (allocated(self%local_map_ghost    )) deallocate(self%local_map_ghost    )
   if (allocated(self%comm_map_send_ghost)) deallocate(self%comm_map_send_ghost)
   if (allocated(self%comm_map_recv_ghost)) deallocate(self%comm_map_recv_ghost)
   my_fec_number = 0
   recv_fec_number = 0
   send_fec_number = 0
   self%comm_map_n_recv_ghost = 0
   self%comm_map_n_send_ghost = 0
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      do fec=1, 26
         weight_reduction = 2 ** count(FEC_TO_DELTA(:, fec)==0_I4P, dim=1)
         if (allocated(node_ptr%neighbor(fec)%codes)) then
            neighbor      = node_ptr%neighbor(fec)%codes
            neighbor_type = node_ptr%neighbor(fec)%ntype
            if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
               do n=1, size(neighbor, dim=1)
                  neigh => tree%node(code=neighbor(n))
                  if (.not.associated(neigh)) then
                     call mpih%print_message(msg=': error neighbor n='//trim(str(n))//&
                                                      ' of neighbors='//trim(str(neighbor))//' not associated')
                  endif
                  if     ((mpih%myrank == neigh%myrank).and.(mpih%myrank == node_ptr%myrank)) then
                     my_fec_number = my_fec_number + 1
                  elseif ((mpih%myrank /= neigh%myrank).and.(mpih%myrank == node_ptr%myrank)) then
                     ! when receiving from same or less refined than me the size of the message is full, when receiving from
                     ! more refined the message is an averaged portion (reduced size)
                     recv_fec_number = recv_fec_number + 1
                     if (neighbor_type==NODE_MORE_REFINED) then
                        self%comm_map_n_recv_ghost(neigh%myrank) = self%comm_map_n_recv_ghost(neigh%myrank) + &
                                                                   grid%weight_neighbor(fec)/ weight_reduction
                     else
                        self%comm_map_n_recv_ghost(neigh%myrank) = self%comm_map_n_recv_ghost(neigh%myrank) + &
                                                                   grid%weight_neighbor(fec)
                     endif
                  elseif ((mpih%myrank == neigh%myrank).and.(mpih%myrank /= node_ptr%myrank)) then
                     ! when sending to same or more refined than me the size of the message is full, when sending to less
                     ! refined the message is an averaged portion (reduced size)
                     send_fec_number = send_fec_number + 1
                     if (neighbor_type==NODE_MORE_REFINED) then
                        self%comm_map_n_send_ghost(node_ptr%myrank) = self%comm_map_n_send_ghost(node_ptr%myrank) + &
                                                                      grid%weight_neighbor(fec)/ weight_reduction
                     else
                        self%comm_map_n_send_ghost(node_ptr%myrank) = self%comm_map_n_send_ghost(node_ptr%myrank) + &
                                                                      grid%weight_neighbor(fec)
                     endif
                  endif
               enddo
            endif
         endif
      enddo
   enddo
   ! create pointer for each process in the send/recv buffers
   self%comm_map_send_ptr_ghost = 0_I4P
   self%comm_map_recv_ptr_ghost = 0_I4P
   do p=1, mpih%procs_number
      self%comm_map_send_ptr_ghost(p) = self%comm_map_send_ptr_ghost(p-1) + self%comm_map_n_send_ghost(p-1)
      self%comm_map_recv_ptr_ghost(p) = self%comm_map_recv_ptr_ghost(p-1) + self%comm_map_n_recv_ghost(p-1)
   enddo
   ! allocate maps
   if (my_fec_number>0  ) allocate(self%local_map_ghost(    1:my_fec_number,  1:13))
   if (send_fec_number>0) allocate(self%comm_map_send_ghost(1:send_fec_number,1:15))
   if (recv_fec_number>0) allocate(self%comm_map_recv_ghost(1:recv_fec_number,1:15))
   ! mpi buffers
   ! Allocate only when there is something to exchange: on a single rank
   ! (or any rank with no neighbours) the ghost counts sum to zero, and a
   ! zero-size buffer must not be allocated — it would later be offloaded
   ! to the device and fault. Leaving it unallocated is the correct
   ! "nothing to exchange" state, which the FNL copy layer skips.
   if (allocated(self%send_buffer_ghost)) deallocate(self%send_buffer_ghost)
   if (sum(self%comm_map_n_send_ghost, dim=1) > 0) &
      allocate(self%send_buffer_ghost(1:sum(self%comm_map_n_send_ghost, dim=1)*nv))
   if (allocated(self%recv_buffer_ghost)) deallocate(self%recv_buffer_ghost)
   if (sum(self%comm_map_n_recv_ghost, dim=1) > 0) &
      allocate(self%recv_buffer_ghost(1:sum(self%comm_map_n_recv_ghost, dim=1)*nv))
   endsubroutine alloc_comm_local_maps_ghost

   subroutine count_bc_numbers(self, tree, fec_bc_faces_number, fec_bc_edges_number, fec_bc_corners_number)
   !< Count faces/edges/corners BC numbers.
   class(maps_object), intent(inout) :: self                   !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   integer(I4P),       intent(out)   :: fec_bc_faces_number    !< BC faces number.
   integer(I4P),       intent(out)   :: fec_bc_edges_number    !< BC edges number.
   integer(I4P),       intent(out)   :: fec_bc_corners_number  !< BC corners number.
   type(tree_node_object), pointer   :: node_ptr               !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I4P)                      :: fec                    !< Counter.

   fec_bc_faces_number = 0_I4P
   fec_bc_edges_number = 0_I4P
   fec_bc_corners_number = 0_I4P
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      if (node_ptr%myrank==mpih%myrank) then
         do fec=1, 26
            if (node_ptr%neighbor(fec)%ntype == NODE_BOUNDARY_CONDITION) then
               if (fec<=6)             fec_bc_faces_number   = fec_bc_faces_number   + 1_I4P
               if (fec>=7.and.fec<=18) fec_bc_edges_number   = fec_bc_edges_number   + 1_I4P
               if (fec>=19)            fec_bc_corners_number = fec_bc_corners_number + 1_I4P
            endif
         enddo
      endif
   enddo
   endsubroutine count_bc_numbers

   subroutine make_comm_map_recv_ghost_cell(self, nv)
   !< Make communication receive map of ghost cells in cell order.
   class(maps_object), intent(inout) :: self               !< The maps.
   integer(I4P),       intent(in)    :: nv                 !< Number of field variables.
   integer(I4P)                      :: i, j, k, c, f, v   !< Counter.
   integer(I4P)                      :: iii, jjj, kkk      !< Counter.
   integer(I4P)                      :: ic, jc, kc         !< Counter.
   integer(I4P)                      :: portion            !< Portion of fec updated (0=>whole fec).
   integer(I4P)                      :: b_recv             !< Index of receiving block.
   integer(I4P)                      :: imin               !< Lower limit of i indexes.
   integer(I4P)                      :: jmin               !< Lower limit of j indexes.
   integer(I4P)                      :: kmin               !< Lower limit of j indexes.
   integer(I4P)                      :: imax               !< Upper limit of i indexes.
   integer(I4P)                      :: jmax               !< Upper limit of j indexes.
   integer(I4P)                      :: kmax               !< Upper limit of k indexes.
   integer(I4P)                      :: idelta             !< Delta offset for ghost-inner cells of i.
   integer(I4P)                      :: jdelta             !< Delta offset for ghost-inner cells of j.
   integer(I4P)                      :: kdelta             !< Delta offset for ghost-inner cells of k.
   integer(I4P)                      :: recv_ptr, recv_ctr !< Counter.

   if (allocated(self%comm_map_recv_ghost)) self%comm_map_recv_ghost(:,15) = self%comm_map_recv_ghost(:,15) * nv
   if (allocated(self%comm_map_recv_ptr_ghost)) self%comm_map_recv_ptr_ghost = self%comm_map_recv_ptr_ghost * nv
   if (allocated(self%comm_map_recv_ghost_cell)) deallocate(self%comm_map_recv_ghost_cell)
   if (allocated(self%comm_map_recv_ghost)) then
      c = 0
      do f=1, size(self%comm_map_recv_ghost, dim=1)
         portion = self%comm_map_recv_ghost(f, 5 )
         imin    = self%comm_map_recv_ghost(f, 6 )
         jmin    = self%comm_map_recv_ghost(f, 7 )
         kmin    = self%comm_map_recv_ghost(f, 8 )
         imax    = self%comm_map_recv_ghost(f, 9 )
         jmax    = self%comm_map_recv_ghost(f, 10)
         kmax    = self%comm_map_recv_ghost(f, 11)
         if (portion>=0_I4P) then
            ! receiving from a block at the same level or finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo

      ! Nv variables map
      allocate(self%comm_map_recv_ghost_cell(1:c*nv,1:6))
      c = 1
      do f=1, size(self%comm_map_recv_ghost, dim=1)
         b_recv   = self%comm_map_recv_ghost(f, 1 )
         portion  = self%comm_map_recv_ghost(f, 5 )
         imin     = self%comm_map_recv_ghost(f, 6 )
         jmin     = self%comm_map_recv_ghost(f, 7 )
         kmin     = self%comm_map_recv_ghost(f, 8 )
         imax     = self%comm_map_recv_ghost(f, 9 )
         jmax     = self%comm_map_recv_ghost(f, 10)
         kmax     = self%comm_map_recv_ghost(f, 11)
         idelta   = self%comm_map_recv_ghost(f, 12)
         jdelta   = self%comm_map_recv_ghost(f, 13)
         kdelta   = self%comm_map_recv_ghost(f, 14)
         recv_ptr = self%comm_map_recv_ghost(f, 15)
         if (portion>=0_I4P) then
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1, nv
                        self%comm_map_recv_ghost_cell(c, 1 ) = recv_ptr + recv_ctr
                        self%comm_map_recv_ghost_cell(c,2:6) = [b_recv,i,j,k,v]
                        recv_ctr = recv_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        do v=1, nv
                           self%comm_map_recv_ghost_cell(c, 1 ) = recv_ptr + recv_ctr
                           self%comm_map_recv_ghost_cell(c,2:6) = [b_recv,iii+ic,jjj+jc,kkk+kc,v]
                           recv_ctr = recv_ctr + 1
                           c = c + 1
                        enddo
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
   endif
   endsubroutine make_comm_map_recv_ghost_cell

   subroutine make_comm_map_send_ghost_cell(self, nv)
   !< Make communication send map of ghost cells in cell order.
   class(maps_object), intent(inout) :: self                !< The maps.
   integer(I4P),       intent(in)    :: nv                  !< Number of field variables.
   integer(I4P)                      :: i, j, k, c, f, v, n !< Counter.
   integer(I4P)                      :: iii, jjj, kkk       !< Counter.
   integer(I4P)                      :: ic, jc, kc          !< Counter.
   integer(I4P)                      :: portion             !< Portion of fec updated (0=>whole fec).
   integer(I4P)                      :: b_send              !< Index of sending block.
   integer(I4P)                      :: imin                !< Lower limit of i indexes.
   integer(I4P)                      :: jmin                !< Lower limit of j indexes.
   integer(I4P)                      :: kmin                !< Lower limit of j indexes.
   integer(I4P)                      :: imax                !< Upper limit of i indexes.
   integer(I4P)                      :: jmax                !< Upper limit of j indexes.
   integer(I4P)                      :: kmax                !< Upper limit of k indexes.
   integer(I4P)                      :: idelta              !< Delta offset for ghost-inner cells of i.
   integer(I4P)                      :: jdelta              !< Delta offset for ghost-inner cells of j.
   integer(I4P)                      :: kdelta              !< Delta offset for ghost-inner cells of k.
   integer(I4P)                      :: send_ptr, send_ctr  !< Counter.

   if (allocated(self%comm_map_send_ghost)) self%comm_map_send_ghost(:,15) = self%comm_map_send_ghost(:,15) * nv
   if (allocated(self%comm_map_send_ptr_ghost)) self%comm_map_send_ptr_ghost = self%comm_map_send_ptr_ghost * nv
   if (allocated(self%comm_map_send_ghost_cell)) deallocate(self%comm_map_send_ghost_cell)
   if (allocated(self%comm_map_send_ghost)) then
      c = 0
      do f=1, size(self%comm_map_send_ghost, dim=1)
         portion = self%comm_map_send_ghost(f, 5 )
         imin    = self%comm_map_send_ghost(f, 6 )
         jmin    = self%comm_map_send_ghost(f, 7 )
         kmin    = self%comm_map_send_ghost(f, 8 )
         imax    = self%comm_map_send_ghost(f, 9 )
         jmax    = self%comm_map_send_ghost(f, 10)
         kmax    = self%comm_map_send_ghost(f, 11)
         if (portion>=0_I4P) then
            ! sending to a block at my level or to a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! sending to a block finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      ! Nv variables map
      allocate(self%comm_map_send_ghost_cell(1:c*nv,1:7))
      c = 1
      do f=1, size(self%comm_map_send_ghost, dim=1)
         b_send    = self%comm_map_send_ghost(f, 2 )
         portion   = self%comm_map_send_ghost(f, 5 )
         imin      = self%comm_map_send_ghost(f, 6 )
         jmin      = self%comm_map_send_ghost(f, 7 )
         kmin      = self%comm_map_send_ghost(f, 8 )
         imax      = self%comm_map_send_ghost(f, 9 )
         jmax      = self%comm_map_send_ghost(f, 10)
         kmax      = self%comm_map_send_ghost(f, 11)
         idelta    = self%comm_map_send_ghost(f, 12)
         jdelta    = self%comm_map_send_ghost(f, 13)
         kdelta    = self%comm_map_send_ghost(f, 14)
         send_ptr  = self%comm_map_send_ghost(f, 15)
         if (portion==0_I4P) then
            ! sending to a block at my level
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1, nv
                        self%comm_map_send_ghost_cell(c, 1:4) = [b_send,i+idelta,j+jdelta,k+kdelta]
                        self%comm_map_send_ghost_cell(c,  5 ) = v
                        self%comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                        self%comm_map_send_ghost_cell(c,  7 ) = 1
                        send_ctr = send_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
            ! sending to a block finer than me
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do n=1,8
                        do v=1, nv
                           self%comm_map_send_ghost_cell(c, 1:4) = [b_send,i,j,k]
                           self%comm_map_send_ghost_cell(c,  5 ) = v
                           self%comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                           self%comm_map_send_ghost_cell(c,  7 ) = 1
                           send_ctr = send_ctr + 1
                           c = c + 1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         else
            ! sending to a block coarser than me, loop is over the coarser grid
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     do v=1, nv
                        self%comm_map_send_ghost_cell(c, 1:4) = [b_send,iii,jjj,kkk]
                        self%comm_map_send_ghost_cell(c,  5 ) = v
                        self%comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                        self%comm_map_send_ghost_cell(c,  7 ) = 8
                        send_ctr = send_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
   endif
   endsubroutine make_comm_map_send_ghost_cell

   subroutine make_local_map_ghost_cell(self)
   !< Make local map of ghost cells in cell order.
   class(maps_object), intent(inout) :: self          !< The maps.
   integer(I4P)                      :: i, j, k, c, f !< Counter.
   integer(I4P)                      :: iii, jjj, kkk !< Counter.
   integer(I4P)                      :: ic, jc, kc    !< Counter.
   integer(I4P)                      :: portion       !< Portion of fec updated (0=>whole fec).
   integer(I4P)                      :: b_recv        !< Index of receiving block.
   integer(I4P)                      :: b_send        !< Index of sending block.
   integer(I4P)                      :: imin          !< Lower limit of i indexes.
   integer(I4P)                      :: jmin          !< Lower limit of j indexes.
   integer(I4P)                      :: kmin          !< Lower limit of j indexes.
   integer(I4P)                      :: imax          !< Upper limit of i indexes.
   integer(I4P)                      :: jmax          !< Upper limit of j indexes.
   integer(I4P)                      :: kmax          !< Upper limit of k indexes.
   integer(I4P)                      :: idelta        !< Delta offset for ghost-inner cells of i.
   integer(I4P)                      :: jdelta        !< Delta offset for ghost-inner cells of j.
   integer(I4P)                      :: kdelta        !< Delta offset for ghost-inner cells of k.

   if (allocated(self%local_map_ghost_cell)) deallocate(self%local_map_ghost_cell)
   if (allocated(self%local_map_ghost)) then
      c = 0
      do f=1, size(self%local_map_ghost, dim=1)
         portion = self%local_map_ghost(f, 4 )
         imin    = self%local_map_ghost(f, 5 )
         jmin    = self%local_map_ghost(f, 6 )
         kmin    = self%local_map_ghost(f, 7 )
         imax    = self%local_map_ghost(f, 8 )
         jmax    = self%local_map_ghost(f, 9 )
         kmax    = self%local_map_ghost(f, 10)
         if (portion>=0) then
            ! receiving from a block with the same refinement or finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      allocate(self%local_map_ghost_cell(1:c,1:9))
      c = 1
      do f=1, size(self%local_map_ghost, dim=1)
         b_recv  = self%local_map_ghost(f, 1 )
         b_send  = self%local_map_ghost(f, 2 )
         portion = self%local_map_ghost(f, 4 )
         imin    = self%local_map_ghost(f, 5 )
         jmin    = self%local_map_ghost(f, 6 )
         kmin    = self%local_map_ghost(f, 7 )
         imax    = self%local_map_ghost(f, 8 )
         jmax    = self%local_map_ghost(f, 9 )
         kmax    = self%local_map_ghost(f, 10)
         idelta  = self%local_map_ghost(f, 11)
         jdelta  = self%local_map_ghost(f, 12)
         kdelta  = self%local_map_ghost(f, 13)
         if     (portion==0) then
            ! receiving from a block with the same refinement
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     self%local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                     self%local_map_ghost_cell(c,3:5) = [i+idelta, j+jdelta, k+kdelta]
                     self%local_map_ghost_cell(c,6:8) = [i, j, k]
                     self%local_map_ghost_cell(c, 9 ) = 1
                     c = c + 1
                  enddo
               enddo
            enddo
         elseif (portion>0) then
            ! receiving from a block finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     self%local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                     self%local_map_ghost_cell(c,3:5) = [iii, jjj, kkk]
                     self%local_map_ghost_cell(c,6:8) = [i, j, k]
                     self%local_map_ghost_cell(c, 9 ) = 8
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        self%local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                        self%local_map_ghost_cell(c,3:5) = [i, j, k]
                        self%local_map_ghost_cell(c,6:8) = [iii+ic,  jjj+jc, kkk+kc]
                        self%local_map_ghost_cell(c, 9 ) = 1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
   endif
   endsubroutine make_local_map_ghost_cell

   function ijk_mmd(self, grid, fec, neighbor_bc_fec) result(ijk_min_max_delta)
   !< Return IJK min/max/delta.
   class(maps_object), intent(in) :: self                   !< The maps.
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P),       intent(in) :: fec                    !< Current fec number.
   integer(I4P),       intent(in) :: neighbor_bc_fec        !< Neighbors fec for BC.
   integer(I8P)                   :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
   integer(I4P)                   :: ijkmin(3)              !< Lower limit of ijk indexes.
   integer(I4P)                   :: ijkmax(3)              !< Upper limit of ijk indexes.
   integer(I4P)                   :: ijkdelta(3)            !< Delta offset for ghost-inner cells mapping same refinement.
   integer(I4P)                   :: ijkdelta_global(3)     !< Delta offset for ghost-inner cells mapping same refinement.
   integer(I4P)                   :: nijk(3)                !< Ni, Nj , Nk stored in array.
   integer(I4P)                   :: i                      !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc)
   nijk = [ni, nj, nk]

   ijkdelta = FEC_TO_DELTA(1:3, fec)
   ijkdelta_global = FEC_TO_DELTA(1:3, neighbor_bc_fec)

   do i=1, 3
      if     (ijkdelta(i)==1) then
         ijkmin(i) = nijk(i) + 1
         ijkmax(i) = nijk(i) + ngc
      elseif (ijkdelta(i)==-1) then
         ijkmin(i) = 0
         ijkmax(i) = 1 - ngc
      elseif (ijkdelta(i)==0) then
         ijkmin(i) = 1
         ijkmax(i) = nijk(i)
      endif
   enddo
   ijk_min_max_delta = [ijkmin, ijkmax, ijkdelta_global]
   endassociate
   endfunction ijk_mmd

   function ijk_mmd_ghost(self, grid, fec, portion) result(ijk_min_max_delta)
   !< Return IJK min/max/delta, ghost maps.
   class(maps_object), intent(in) :: self                   !< The maps.
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   integer(I4P),       intent(in) :: fec                    !< Current fec number.
   integer(I8P),       intent(in) :: portion                !< Current portion.
   integer(I8P)                   :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
   integer(I4P)                   :: abs_portion            !< ABS of current portion.
   integer(I4P)                   :: portion_cur            !< Current portion of fec updated (0=>whole fec).
   integer(I4P)                   :: ijkmin(3)              !< Lower limit of ijk indexes.
   integer(I4P)                   :: ijkmax(3)              !< Upper limit of ijk indexes.
   integer(I4P)                   :: ijkdelta(3)            !< Delta offset for ghost-inner cells mapping same refinement.
   integer(I4P)                   :: delta(3)               !< Neighbor delta of current fec.
   integer(I4P)                   :: nijk(3)                !< Ni, Nj , Nk stored in array.
   integer(I4P)                   :: portion_array(3)       !< Portion array binary useful for less refined case.
   integer(I4P)                   :: i                      !< Counter.

   associate(ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, ngc=>grid%ngc)
   nijk = [ni, nj, nk]

   abs_portion = abs(portion)
   delta = FEC_TO_DELTA(1:3, fec)

   if     (portion==0) then
      do i=1, 3
         if     (delta(i)==1) then
            ijkmin(i) = nijk(i) + 1
            ijkmax(i) = nijk(i) + ngc
            ijkdelta(i) = -nijk(i)
         elseif (delta(i)==-1) then
            ijkmin(i) = 1 - ngc
            ijkmax(i) = 0
            ijkdelta(i) = nijk(i)
         elseif (delta(i)==0) then
            ijkmin(i) = 1
            ijkmax(i) = nijk(i)
            ijkdelta(i) = 0
         endif
      enddo
   elseif (portion>0) then
      do i=1, 3
         if     (delta(i)==1) then
            ijkmin(i) = nijk(i) + 1
            ijkmax(i) = nijk(i) + ngc
            ijkdelta(i) = -1 - 2 * nijk(i)
         elseif (delta(i)==-1) then
            ijkmin(i) = 1 - ngc
            ijkmax(i) = 0
            ijkdelta(i) = -1 + nijk(i)
         elseif (delta(i)==0) then
            portion_cur = mod(abs_portion-1, 2)
            abs_portion = 2 * portion_cur + (abs_portion - 1) / 2 + 1
            ijkmin(i) = portion_cur * nijk(i) / 2 + 1
            ijkmax(i) = ijkmin(i) + nijk(i) / 2 - 1
            ijkdelta(i) = -2 * ijkmin(i) + 1
         endif
      enddo
   else
      portion_array(1) = mod(abs_portion-1,2)
      portion_array(2) = mod((abs_portion-1)/2,2)
      portion_array(3) = mod((abs_portion-1)/4,2)
      do i=1, 3
         if     (delta(i)==1) then
            ijkmin(i) = nijk(i) / 2 * portion_array(i) + 1
            ijkmax(i) = ijkmin(i) + ngc / 2 - 1
            ijkdelta(i) = - nijk(i) * portion_array(i) + nijk(i) - 1
         elseif (delta(i)==-1) then
            ijkmin(i) = nijk(i) / 2 + nijk(i) / 2 * portion_array(i) + 1 - ngc / 2
            ijkmax(i) = ijkmin(i) + ngc / 2 - 1
            ijkdelta(i) = - nijk(i) -  nijk(i) * portion_array(i) - 1
         elseif (delta(i)==0) then
            ijkmin(i) = nijk(i) / 2 * portion_array(i) + 1
            ijkmax(i) = ijkmin(i) + nijk(i) / 2 - 1
            ijkdelta(i) = - nijk(i) * portion_array(i) - 1
         endif
      enddo
   endif
   ijk_min_max_delta = [ijkmin, ijkmax, ijkdelta]
   endassociate
   endfunction ijk_mmd_ghost

   subroutine populate_comm_local_maps_ghost(self, tree, grid)
   !< Populate communication/local maps of ghost cells, `local_map_ghost, comm_map_send_ghost, comm_map_recv_ghost`.
   class(maps_object), intent(inout) :: self                       !< The maps.
   type(tree_object),  intent(in) :: tree !< Tree (sibling realm component, threaded in).
   type(grid_object),     intent(in)              :: grid !< Grid (sibling realm component, threaded in).
   type(tree_node_object), pointer   :: node_ptr                   !< Pointer to current node.
   type(tree_iterator_object)                      :: iter                 !< Tree traversal cursor (issue #13: re-entrant loop).
   integer(I8P), allocatable         :: neighbor(:)                !< List of code neighbors.
   type(tree_node_object), pointer   :: neigh                      !< Pointer to node neighbor.
   integer(I4P)                      :: neighbor_type              !< Neighbors type.
   integer(I4P)                      :: neighbor_portion           !< Neighbors portion.
   integer(I4P)                      :: fec, mf, rf, sf, n         !< Counter.
   integer(I4P)                      :: weight_reduction           !< Neighbor weight reduction for send/recv at different level.
   integer(I4P), allocatable         :: comm_map_send_ctr_ghost(:) !< Communication map, counters to send.
   integer(I4P), allocatable         :: comm_map_recv_ctr_ghost(:) !< Communication map, counters to recv.

   comm_map_send_ctr_ghost = self%comm_map_send_ptr_ghost
   comm_map_recv_ctr_ghost = self%comm_map_recv_ptr_ghost
   sf = 0
   rf = 0
   mf = 0
   iter%b = 1_I4P ; iter%p => null()
   do while(tree%loop(iter, node_ptr=node_ptr))
      do fec=1, 26
         weight_reduction = 2 ** count(FEC_TO_DELTA(:, fec)==0_I4P, dim=1)
         if (allocated(node_ptr%neighbor(fec)%codes)) then
            neighbor         = node_ptr%neighbor(fec)%codes
            neighbor_type    = node_ptr%neighbor(fec)%ntype
            neighbor_portion = node_ptr%neighbor(fec)%portion
            if (neighbor_type /= NODE_BOUNDARY_CONDITION) then
               do n=1, size(neighbor, dim=1)
                  neigh => tree%node(code=neighbor(n))
                  if     ((mpih%myrank == neigh%myrank).and.(mpih%myrank == node_ptr%myrank)) then
                     mf = mf + 1
                     self%local_map_ghost(mf, 1) = node_ptr%block_index
                     self%local_map_ghost(mf, 2) = neigh%block_index
                     self%local_map_ghost(mf, 3) = fec
                     if     (neighbor_type==NODE_STANDARD) then
                        self%local_map_ghost(mf, 4) = 0
                     elseif (neighbor_type==NODE_MORE_REFINED) then
                        self%local_map_ghost(mf, 4) = n
                     elseif (neighbor_type==NODE_LESS_REFINED) then
                        self%local_map_ghost(mf, 4) = -neighbor_portion
                     endif
                     self%local_map_ghost(mf, 5:13)=self%ijk_mmd_ghost(grid=grid, fec=fec, portion=self%local_map_ghost(mf, 4))
                  elseif ((mpih%myrank /= neigh%myrank).and.(mpih%myrank == node_ptr%myrank)) then
                     rf = rf + 1
                     self%comm_map_recv_ghost(rf, 15) = comm_map_recv_ctr_ghost(neigh%myrank)
                     self%comm_map_recv_ghost(rf, 1) = node_ptr%block_index
                     self%comm_map_recv_ghost(rf, 2) = neigh%block_index
                     self%comm_map_recv_ghost(rf, 3) = neigh%myrank
                     self%comm_map_recv_ghost(rf, 4) = fec
                     if     (neighbor_type==NODE_STANDARD) then
                        self%comm_map_recv_ghost(rf, 5) = 0
                        comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                                grid%weight_neighbor(fec)
                     elseif (neighbor_type==NODE_MORE_REFINED) then
                        self%comm_map_recv_ghost(rf, 5) = n
                        comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                                grid%weight_neighbor(fec) / weight_reduction
                     elseif (neighbor_type==NODE_LESS_REFINED) then
                        self%comm_map_recv_ghost(rf, 5) = -neighbor_portion
                        comm_map_recv_ctr_ghost(neigh%myrank) = comm_map_recv_ctr_ghost(neigh%myrank) + &
                                                                grid%weight_neighbor(fec)
                     endif
                     self%comm_map_recv_ghost(rf, 6:14)=self%ijk_mmd_ghost(grid=grid, fec=fec, &
                                                         portion=self%comm_map_recv_ghost(rf, 5))
                  elseif ((mpih%myrank == neigh%myrank).and.(mpih%myrank /= node_ptr%myrank)) then
                     sf = sf + 1
                     self%comm_map_send_ghost(sf, 15) = comm_map_send_ctr_ghost(node_ptr%myrank)
                     self%comm_map_send_ghost(sf, 1) = node_ptr%block_index
                     self%comm_map_send_ghost(sf, 2) = neigh%block_index
                     self%comm_map_send_ghost(sf, 3) = node_ptr%myrank
                     self%comm_map_send_ghost(sf, 4) = fec
                     if     (neighbor_type==NODE_STANDARD) then
                        self%comm_map_send_ghost(sf, 5) = 0
                        comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                   grid%weight_neighbor(fec)
                     elseif (neighbor_type==NODE_MORE_REFINED) then
                        self%comm_map_send_ghost(sf, 5) = n
                        comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                   grid%weight_neighbor(fec) / weight_reduction
                     elseif (neighbor_type==NODE_LESS_REFINED) then
                        self%comm_map_send_ghost(sf, 5) = -neighbor_portion
                        comm_map_send_ctr_ghost(node_ptr%myrank) = comm_map_send_ctr_ghost(node_ptr%myrank) + &
                                                                   grid%weight_neighbor(fec)
                     endif
                     self%comm_map_send_ghost(sf, 6:14)=self%ijk_mmd_ghost(grid=grid, fec=fec, &
                                                         portion=self%comm_map_send_ghost(sf, 5))
                  endif
               enddo
            endif
         endif
      enddo
   enddo
   endsubroutine populate_comm_local_maps_ghost

   ! non TBP
   function bc_cells_number(local_map_bc) result(cells_number)
   !< Return BC cells number.
   integer(I8P), intent(in) :: local_map_bc(:,:) !< Local map for BC ghost cells.
   integer(I8P)             :: f, i, j, k        !< Counter.
   integer(I8P)             :: imin              !< Lower limit of ijk indexes.
   integer(I8P)             :: jmin              !< Lower limit of ijk indexes.
   integer(I8P)             :: kmin              !< Lower limit of ijk indexes.
   integer(I8P)             :: imax              !< Upper limit of ijk indexes.
   integer(I8P)             :: jmax              !< Upper limit of ijk indexes.
   integer(I8P)             :: kmax              !< Upper limit of ijk indexes.
   integer(I4P)             :: cells_number      !< Number of BC cells.

   cells_number = 0
   do f=1, size(local_map_bc, dim=1)
      imin    = local_map_bc(f, 3 )
      jmin    = local_map_bc(f, 4 )
      kmin    = local_map_bc(f, 5 )
      imax    = local_map_bc(f, 6 )
      jmax    = local_map_bc(f, 7 )
      kmax    = local_map_bc(f, 8 )
      do k=kmin, kmax, sign(1_I8P, kmax-kmin)
         do j=jmin, jmax, sign(1_I8P, jmax-jmin)
            do i=imin, imax, sign(1_I8P, imax-imin)
               cells_number = cells_number + 1
            enddo
         enddo
      enddo
   enddo
   endfunction bc_cells_number

   subroutine populate_local_map_bc_crown(local_map_bc, local_map_bc_crown, c_crown)
   !< Populate map of BC cells in crown order.
   integer(I8P), intent(in)    :: local_map_bc(1:,1:)          !< Local map for BC ghost cells.
   integer(I8P), intent(inout) :: local_map_bc_crown(1:,1:,1:) !< Local map for face BC ghost cells, crown order.
   integer(I4P), intent(inout) :: c_crown(1:)                  !< Counter.
   integer(I8P)                :: f, i, j, k, b                !< Counter.
   integer(I8P)                :: imin                         !< Lower limit of ijk indexes.
   integer(I8P)                :: jmin                         !< Lower limit of ijk indexes.
   integer(I8P)                :: kmin                         !< Lower limit of ijk indexes.
   integer(I8P)                :: imax                         !< Upper limit of ijk indexes.
   integer(I8P)                :: jmax                         !< Upper limit of ijk indexes.
   integer(I8P)                :: kmax                         !< Upper limit of ijk indexes.
   integer(I8P)                :: idelta                       !< IJK delta step for extrapolation.
   integer(I8P)                :: jdelta                       !< IJK delta step for extrapolation.
   integer(I8P)                :: kdelta                       !< IJK delta step for extrapolation.
   integer(I8P)                :: bc_type                      !< Boundary condition type.
   integer(I8P)                :: crown                        !< Crown counter.
   integer(I8P)                :: fec                          !< Boundary condition fec.

   do f=1, size(local_map_bc, dim=1)
      b       = local_map_bc(f, 1 )
      fec     = local_map_bc(f, 2 )
      imin    = local_map_bc(f, 3 )
      jmin    = local_map_bc(f, 4 )
      kmin    = local_map_bc(f, 5 )
      imax    = local_map_bc(f, 6 )
      jmax    = local_map_bc(f, 7 )
      kmax    = local_map_bc(f, 8 )
      idelta  = local_map_bc(f, 9 )
      jdelta  = local_map_bc(f, 10)
      kdelta  = local_map_bc(f, 11)
      bc_type = local_map_bc(f, 12)
      do k=kmin, kmax, sign(1_I8P, kmax-kmin)
         do j=jmin, jmax, sign(1_I8P, jmax-jmin)
            do i=imin, imax, sign(1_I8P, imax-imin)
               crown = maxval([abs(i-imin), abs(j-jmin), abs(k-kmin)], mask=[imin/=1, jmin/=1, kmin/=1]) + 1
               local_map_bc_crown(c_crown(crown), 1, crown) = b
               local_map_bc_crown(c_crown(crown), 2, crown) = i
               local_map_bc_crown(c_crown(crown), 3, crown) = j
               local_map_bc_crown(c_crown(crown), 4, crown) = k
               local_map_bc_crown(c_crown(crown), 5, crown) = idelta
               local_map_bc_crown(c_crown(crown), 6, crown) = jdelta
               local_map_bc_crown(c_crown(crown), 7, crown) = kdelta
               local_map_bc_crown(c_crown(crown), 8, crown) = bc_type
               local_map_bc_crown(c_crown(crown), 9, crown) = fec
               c_crown(crown) = c_crown(crown) + 1
            enddo
         enddo
      enddo
   enddo
   endsubroutine populate_local_map_bc_crown

   subroutine set_local_map_bc(fec, block_index, neighbor_bc_fec, ijk_min_max_delta, fec_bc_type, local_map_bc)
   !< Set data into faces/edges/corners local map BC.
   integer(I4P), intent(in)    :: fec                    !< Current fec.
   integer(I8P), intent(in)    :: block_index            !< Block index in field blocks array.
   integer(I8P), intent(in)    :: neighbor_bc_fec        !< Neighbors fec for BC.
   integer(I8P), intent(in)    :: ijk_min_max_delta(1:9) !< IJK min/max/delta.
   integer(I8P), intent(in)    :: fec_bc_type            !< BC type.
   integer(I8P), intent(inout) :: local_map_bc(1:,1:)    !< Local map BC data.

   local_map_bc(fec,1)    = block_index
   local_map_bc(fec,2)    = neighbor_bc_fec
   local_map_bc(fec,3:11) = ijk_min_max_delta
   local_map_bc(fec,12)   = fec_bc_type
   endsubroutine set_local_map_bc
endmodule adam_maps_object
