!< ADAM, maps class definition, GMP backend.
module adam_maps_gmp_object
!< ADAM, maps class definition, GMP backend.

use adam_maps_object, only : maps_object
use adam_memory_gmp_library, only : assign_allocatable_gpu
use adam_mpih_gmp_object, only : mpih_gmp_object
use penf
use mpi

implicit none
private
public :: maps_gmp_object

type :: maps_gmp_object
   !< Maps class, GMP backend.
   ! ADAM library objects
   type(mpih_gmp_object), pointer :: mpih=>null() !< MPI handler, GMP backend.
   type(maps_object),     pointer :: maps=>null() !< The maps.
   ! GPU data
   integer(I8P), pointer, contiguous :: local_map_ghost_cell_gpu(:,:)=>null()     !< Local map for ghost cells updating, cells order.
   integer(I8P), pointer, contiguous :: comm_map_recv_ghost_cell_gpu(:,:)=>null() !< Communication map, `fec` information, cell order.
   integer(I8P), pointer, contiguous :: comm_map_send_ghost_cell_gpu(:,:)=>null() !< Communication map, `fec` information, cell order.
   real(R8P),    pointer, contiguous :: send_buffer_ghost_gpu(:)=>null()          !< Send buffer of ghost cells.
   real(R8P),    pointer, contiguous :: recv_buffer_ghost_gpu(:)=>null()          !< Receive buffer of ghost cells.
   integer(I8P), pointer, contiguous :: local_map_bc_crown_gpu(:,:,:)=>null()     !< Local map for face BC ghost cells, "crown" order.
   contains
      procedure, pass(self) :: copy_cpu_gpu !< Copy data from (maps_object) CPU to (maps_gmp_object) GPU.
      procedure, pass(self) :: initialize   !< Initialize MPI handler data.
endtype maps_gmp_object

contains
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (maps_object) CPU to (maps_gmp_object) GPU.
   class(maps_gmp_object), intent(inout)        :: self     !< The maps.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.
   character(:), allocatable                    :: r        !< My rank stringified.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('maps_gmp_object%copy_cpu_gpu start')
   r = self%mpih%myrankstr//'maps_gmp_object%copy_cpu_gpu '
   call assign_allocatable_gpu(lhs=self%local_map_ghost_cell_gpu, rhs=self%maps%local_map_ghost_cell, &
                               omp_dev=self%mpih%mydev, msg=r//' local_map_ghost_cell_gpu ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%comm_map_send_ghost_cell_gpu, rhs=self%maps%comm_map_send_ghost_cell, &
                               omp_dev=self%mpih%mydev, msg=r//' comm_map_send_ghost_cell_gpu ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%comm_map_recv_ghost_cell_gpu, rhs=self%maps%comm_map_recv_ghost_cell, &
                               omp_dev=self%mpih%mydev, msg=r//' comm_map_recv_ghost_cell_gpu ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%send_buffer_ghost_gpu, rhs=self%maps%send_buffer_ghost, &
                               omp_dev=self%mpih%mydev, msg=r//' send_buffer_ghost_gpu ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%recv_buffer_ghost_gpu, rhs=self%maps%recv_buffer_ghost, &
                               omp_dev=self%mpih%mydev, msg=r//' recv_buffer_ghost_gpu ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%local_map_bc_crown_gpu, rhs=self%maps%local_map_bc_crown, &
                               omp_dev=self%mpih%mydev, msg=r//' local_map_bc_crown_gpu ', verbose=verbose)
   if (verbose_) call self%mpih%print_message('maps_gmp_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine initialize(self, mpih, maps)
   !< Initialize maps.
   class(maps_gmp_object), intent(inout)      :: self !< The maps, gmp backend.
   type(mpih_gmp_object),  intent(in), target :: mpih !< MPI handler, GMP backend.
   type(maps_object),      intent(in), target :: maps !< The maps.

   self%mpih => mpih
   call self%mpih%print_message('maps_gmp_object%initialize start')
   self%maps => maps
   call self%copy_cpu_gpu(verbose=.true.)
   call self%mpih%print_message('maps_gmp_object%initialize finish')
   endsubroutine initialize
endmodule adam_maps_gmp_object
