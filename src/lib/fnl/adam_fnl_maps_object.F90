!< ADAM, maps class definition, FNL backend.
module adam_fnl_maps_object
!< ADAM, maps class definition, FNL backend.

! ADAM singleton objects
use :: adam_maps_global,     only : maps
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
   contains
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from (maps global singleton) CPU to (maps_fnl_object) GPU.
      procedure, pass(self) :: initialize   !< Initialize MPI handler data.
endtype maps_fnl_object

contains
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (maps global singleton) CPU to (maps_fnl_object) GPU.
   class(maps_fnl_object), intent(inout)        :: self     !< The maps.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu start')
   ! Each source map is an `allocatable` component of the CPU `maps`
   ! singleton, populated conditionally by maps_object's make_* routines.
   ! On a single MPI rank the inter-rank communication maps and buffers
   ! have nothing to build and are legitimately left unallocated. Passing
   ! an unallocated allocatable as the non-optional assumed-shape
   ! `intent(in)` `src` of dev_assign_to_device is illegal Fortran. Guard
   ! every copy with allocated(): an empty CPU map leaves the GPU pointer
   ! null, the correct device-side state for "nothing to exchange".
   if (allocated(maps%local_map_ghost_cell)) then
      call dev_assign_to_device(dst=self%local_map_ghost_cell_gpu, src=maps%local_map_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy local_map_ghost_cell_gpu done')
   else if (verbose_) then
      call mpih_fnl%print_message('skip local_map_ghost_cell_gpu (CPU map not allocated)')
   endif
   if (allocated(maps%comm_map_send_ghost_cell)) then
      call dev_assign_to_device(dst=self%comm_map_send_ghost_cell_gpu, src=maps%comm_map_send_ghost_cell)
      if (verbose_) call mpih_fnl%print_message('copy comm_map_send_ghost_cell done')
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
   if (verbose_) call mpih_fnl%print_message('maps_fnl_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine initialize(self)
   !< Initialize maps from the program-scope `maps` singleton (adam_maps_global).
   !< Requires `mpih_fnl` (adam_fnl_mpih_global) and `maps` (adam_maps_global) to be
   !< initialized before calling.
   class(maps_fnl_object), intent(inout) :: self !< The maps, FNL backend.

   call mpih_fnl%print_message('maps_fnl_object%initialize start')
   call self%copy_cpu_gpu(verbose=.true.)
   call mpih_fnl%print_message('maps_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_fnl_maps_object
