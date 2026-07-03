!< ADAM, maps class definition, FNL backend.
module adam_fnl_maps_object
!< ADAM, maps class definition, FNL backend.

! ADAM objects
use :: adam_maps_object,              only : maps_object
use :: adam_seam_interpolation_library, only : SEAM_FILL_INJECTION
! ADAM FNL singleton objects
use :: adam_fnl_mpih_global, only : mpih_fnl
! third party modules
use :: fundal
use :: penf
! sdk modules
use :: mpi

implicit none
private
public :: maps_fnl_object

type :: maps_fnl_object
   !< Maps class, FNL backend.
   ! GPU data.
   ! The `=> null()` initializers are mandatory, not cosmetic: dev_assign_to_device
   ! does `if (associated(dst)) call dev_free(dst)` before allocating, and
   ! `associated()` on a pointer that was never nullified is undefined behaviour
   ! (nvfortran traps it as "Null pointer for tmp$r"). A default-initialised
   ! pointer component starts every instance in a defined, disassociated state.
   integer(I8P), pointer :: local_map_ghost_cell_gpu(:,:)     => null() !< Local map for ghost cells updating, cells order.
   integer(I8P), pointer :: comm_map_recv_ghost_cell_gpu(:,:) => null() !< Communication map, `fec` information, cell order.
   integer(I8P), pointer :: comm_map_send_ghost_cell_gpu(:,:) => null() !< Communication map, `fec` information, cell order.
   real(R8P),    pointer :: send_buffer_ghost_gpu(:)          => null() !< Send buffer of ghost cells.
   real(R8P),    pointer :: recv_buffer_ghost_gpu(:)          => null() !< Receive buffer of ghost cells.
   integer(I8P), pointer :: local_map_bc_crown_gpu(:,:,:)     => null() !< Local map for face BC ghost cells, "crown" order.
   ! Inter-realm seam ghost-fill — device-resident counterparts of the host
   ! seam maps/buffers on `maps_object` (see adam_maps_object.F90). The host
   ! arrays drive layout; these are sized/populated 1:1 by `copy_cpu_gpu`.
   ! Per-peer buffers shaped `(buf_len, n_peers)` so concurrent peers do not
   ! alias. Cross-rank seam path (`seam_*_mpi_*`) uses GPU-direct MPI: device
   ! pointers handed straight to MPI_Isend/MPI_Irecv (same model as
   ! update_ghost_mpi_gpu).
   integer(I4P), pointer :: seam_local_map_ghost_cell_gpu(:,:)  => null() !< Per-cell seam ghost map (sorted by peer_realm).
   real(R8P),    pointer :: seam_local_send_buf_gpu(:,:)        => null() !< Per-peer pack buffer, device-resident.
   real(R8P),    pointer :: seam_local_recv_buf_gpu(:,:)        => null() !< Per-peer unpack buffer, device-resident.
   ! Cross-rank seam — declared but not yet populated (same-rank fast path only).
   integer(I4P), pointer :: seam_comm_map_send_ghost_cell_gpu(:,:) => null()
   integer(I4P), pointer :: seam_comm_map_recv_ghost_cell_gpu(:,:) => null()
   real(R8P),    pointer :: seam_mpi_send_buf_gpu(:)               => null()
   real(R8P),    pointer :: seam_mpi_recv_buf_gpu(:)               => null()
   ! Intra-realm coarse->fine seam ghost fill (issue #21 N2 / #22 F3): host-side
   ! mirror of `maps_object%seam_ghost_fill`, the active fill regime consumed by
   ! the flag-4 branches of the ghost kernels (passed as a scalar kernel
   ! argument, not device-resident). Refreshed on every `copy_cpu_gpu`.
   integer(I4P)          :: seam_ghost_fill = SEAM_FILL_INJECTION             !< Active seam ghost-fill regime (host mirror).
   contains
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from (maps global singleton) CPU to (maps_fnl_object) GPU.
      procedure, pass(self) :: initialize   !< Initialize MPI handler data.
endtype maps_fnl_object

contains
   subroutine copy_cpu_gpu(self, maps, verbose)
   !< Copy data from the (realm-local) CPU `maps` to (maps_fnl_object) GPU.
   class(maps_fnl_object), intent(inout)        :: self     !< The maps.
   type(maps_object),      intent(in)           :: maps     !< Realm-local CPU maps; was the `maps` singleton.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu start')
   ! Each source map is an `allocatable` component of the (realm-local) CPU
   ! `maps`, populated conditionally by maps_object's make_* routines.
   ! On a single MPI rank the inter-rank communication maps and buffers
   ! have nothing to build and are legitimately left unallocated. Passing
   ! an unallocated allocatable as the non-optional assumed-shape
   ! `intent(in)` `src` of dev_assign_to_device is illegal Fortran. Guard
   ! every copy with allocated(): an empty CPU map leaves the GPU pointer
   ! null, the correct device-side state for "nothing to exchange".
   ! Intra-realm seam ghost fill (issue #22 F3): flag-4 map rows (coarse->fine
   ! seam interpolation, metadata column included in the verbatim map copy) are
   ! handled on device by the flag-4 branches of adam_fnl_field_kernels; the
   ! active regime is mirrored here for the kernel callers.
   self%seam_ghost_fill = maps%seam_ghost_fill
   if (allocated(maps%local_map_ghost_cell)) then
      call dev_assign_to_device(dst=self%local_map_ghost_cell_gpu, src=maps%local_map_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy local_map_ghost_cell_gpu done ('// &
         trim(str(count(maps%local_map_ghost_cell(:,9) == 4_I8P)))//' seam flag-4 rows)')
   else if (verbose_) then
      call mpih_fnl%print_message('skip local_map_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%comm_map_send_ghost_cell)) then
      call dev_assign_to_device(dst=self%comm_map_send_ghost_cell_gpu, src=maps%comm_map_send_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy comm_map_send_ghost_cell done ('// &
         trim(str(count(maps%comm_map_send_ghost_cell(:,7) == 4_I8P)))//' seam flag-4 rows)')
   else if (verbose_) then
      call mpih_fnl%print_message('skip comm_map_send_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%comm_map_recv_ghost_cell)) then
      call dev_assign_to_device(dst=self%comm_map_recv_ghost_cell_gpu, src=maps%comm_map_recv_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy comm_map_recv_ghost_cell done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip comm_map_recv_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%send_buffer_ghost)) then
      call dev_assign_to_device(dst=self%send_buffer_ghost_gpu, src=maps%send_buffer_ghost)
      if (verbose_) call mpih_fnl%print_message('copy send_buffer_ghost done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip send_buffer_ghost_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%recv_buffer_ghost)) then
      call dev_assign_to_device(dst=self%recv_buffer_ghost_gpu, src=maps%recv_buffer_ghost)
      if (verbose_) call mpih_fnl%print_message('copy recv_buffer_ghost done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip recv_buffer_ghost_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%local_map_bc_crown)) then
      call dev_assign_to_device(dst=self%local_map_bc_crown_gpu, src=maps%local_map_bc_crown)
      if (verbose_) call mpih_fnl%print_message('copy local_map_bc_crown done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip local_map_bc_crown_gpu (CPU map not allocated)')
   endif
   ! Inter-realm seam ghost-fill — device-resident counterparts.
   if (allocated(maps%seam_local_map_ghost_cell)) then
      call dev_assign_to_device(dst=self%seam_local_map_ghost_cell_gpu, src=maps%seam_local_map_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_map_ghost_cell done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip seam_local_map_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%seam_local_send_buf)) then
      call dev_assign_to_device(dst=self%seam_local_send_buf_gpu, src=maps%seam_local_send_buf)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_send_buf done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip seam_local_send_buf_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%seam_local_recv_buf)) then
      call dev_assign_to_device(dst=self%seam_local_recv_buf_gpu, src=maps%seam_local_recv_buf)
      if (verbose_) call mpih_fnl%print_message('copy seam_local_recv_buf done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip seam_local_recv_buf_gpu (CPU map not allocated)')
   endif
   ! Cross-rank seam — host counterparts not yet populated; skip silently when absent.
   if (allocated(maps%seam_comm_map_send_ghost_cell)) then
      call dev_assign_to_device(dst=self%seam_comm_map_send_ghost_cell_gpu, src=maps%seam_comm_map_send_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_comm_map_send_ghost_cell done')
   endif
   if (allocated(maps%seam_comm_map_recv_ghost_cell)) then
      call dev_assign_to_device(dst=self%seam_comm_map_recv_ghost_cell_gpu, src=maps%seam_comm_map_recv_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy seam_comm_map_recv_ghost_cell done')
   endif
   if (allocated(maps%seam_mpi_send_buf)) then
      call dev_assign_to_device(dst=self%seam_mpi_send_buf_gpu, src=maps%seam_mpi_send_buf)
   endif
   if (allocated(maps%seam_mpi_recv_buf)) then
      call dev_assign_to_device(dst=self%seam_mpi_recv_buf_gpu, src=maps%seam_mpi_recv_buf)
   endif
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine initialize(self, maps)
   !< Initialize maps from the (realm-local) CPU `maps`.
   !< Requires `mpih_fnl` (adam_fnl_mpih_global) initialized and `maps` populated before calling.
   class(maps_fnl_object), intent(inout) :: self !< The maps, FNL backend.
   type(maps_object),      intent(in)    :: maps !< Realm-local CPU maps; was the `maps` singleton.

   call mpih_fnl%print_message('maps_fnl_object%initialize start')
   call self%copy_cpu_gpu(maps=maps, verbose=.true.)
   call mpih_fnl%print_message('maps_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_fnl_maps_object
