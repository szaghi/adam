!< ADAM, maps class definition, NVF backend.
module adam_maps_nvf_object
!< ADAM, maps class definition, NVF backend.

use adam_maps_object, only : maps_object
use adam_memory_nvf_library, only : assign_allocatable_gpu
use adam_mpih_nvf_object, only : mpih_nvf_object
use penf
use mpi

implicit none
private
public :: maps_nvf_object

type :: maps_nvf_object
   !< Maps class, NVF backend.
   ! ADAM library objects
   type(mpih_nvf_object)      :: mpih         !< MPI handler, NVF backend.
   type(maps_object), pointer :: maps=>null() !< The maps.
   ! GPU data
   integer(I8P), allocatable, device :: local_map_ghost_cell_gpu(:,:)     !< Local map for ghost cells updating, cells order.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_cell_gpu(:,:) !< Communication map, `fec` information, cell order.
   integer(I8P), allocatable, device :: comm_map_send_ghost_cell_gpu(:,:) !< Communication map, `fec` information, cell order.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)          !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)          !< Receive buffer of ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_crown_gpu(:,:,:)     !< Local map for face BC ghost cells, "crown" order.
   contains
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from (maps_object) CPU to (maps_nvf_object) GPU.
      procedure, pass(self) :: initialize   !< Initialize MPI handler data.
endtype maps_nvf_object

contains
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (maps_object) CPU to (maps_nvf_object) GPU.
   class(maps_nvf_object), intent(inout)        :: self     !< The maps.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   character(:), allocatable                    :: r        !< My rank stringified.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('maps_nvf_object%copy_cpu_gpu start')
   r = self%mpih%myrankstr
   call assign_allocatable_gpu(lhs=self%local_map_ghost_cell_gpu, rhs=self%maps%local_map_ghost_cell, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(local_map_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%comm_map_send_ghost_cell_gpu, rhs=self%maps%comm_map_send_ghost_cell, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(comm_map_send_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%comm_map_recv_ghost_cell_gpu, rhs=self%maps%comm_map_recv_ghost_cell, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(comm_map_recv_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%send_buffer_ghost_gpu, rhs=self%maps%send_buffer_ghost, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(send_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%recv_buffer_ghost_gpu, rhs=self%maps%recv_buffer_ghost, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(recv_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%local_map_bc_crown_gpu, rhs=self%maps%local_map_bc_crown, &
                               msg=r//'maps_nvf_object%copy_cpu_gpu(local_map_bc_crown_gpu) ', verbose=verbose)
   if (verbose_) call self%mpih%print_message('maps_nvf_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine initialize(self, maps)
   !< Initialize maps.
   class(maps_nvf_object), intent(inout)      :: self !< The maps, NVF backend.
   type(maps_object),      intent(in), target :: maps !< The maps.

   call self%mpih%initialize
   call self%mpih%print_message('maps_nvf_object%initialize start')
   self%maps => maps
   call self%copy_cpu_gpu(verbose=.true.)
   call self%mpih%print_message('maps_nvf_object%initialize finish')
   endsubroutine initialize
endmodule adam_maps_nvf_object
